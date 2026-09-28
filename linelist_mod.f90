!=======================================================================================
! linelist_mod -- structure-of-arrays (SoA) bulk linelist handling.
!
! Every line-list line of a composition particle lives in its linelist_type, one per
! isotope. Metals and molecules go through the kernels here; the H I lines keep their
! Stark-profile kernel (hprofl_mod%add_h_opacities, driven by composition_mod) and add
! only their principal quantum numbers to the store (append_h_line).
! The plain linelist_type data
! container lives in the leaf module line_types_mod, so that particle_type (ion_mod) can
! embed it without depending on this heavy module. That inversion lets linelist_mod, in
! turn, USE ion_mod and drive the opacity kernels directly from a particle_type (see
! add_particle_opacity), while the low-level kernels still take plain patch work arrays.
!
! The metal/molecule kernel:
!   * The opacity kernel (add_opacity) loops the SPATIAL domain in the outer loop, then
!     lines (sorted ascending in wavelength), then frequencies. A monotonic frequency
!     cursor ("march") locates each line centre.
!   * Works in log(lambda): the thermal Doppler width is constant per cell and the
!     Doppler shift is additive, so the profile abscissa is vs = idf*(llam(f) - llc).
!   * loggf is folded into the Boltzmann exponential, exp(loggf*ln10 - ex), instead of
!     pre-converting to Einstein Bij.
!
! The orchestration that walks elements/ions/molecules stays in composition_mod
! (add_linelist_opacities / init_linelist), which owns the patch thermodynamic arrays
! and the departure-coefficient lookup; per particle it now calls add_particle_opacity.
!=======================================================================================
MODULE linelist_mod
  USE prec_mod
  USE const_mod
  USE hop_env_mod, only: trace
  USE handyfuncs_mod, only: air2vac, vac2air
  USE op_math_mod, only: find_closest
  USE bbtrans_mod, only: disable_ga, disable_stark, disable_vdw, &
                         f_conv, tol => tolerance, line_reject_frac
  USE line_types_mod, only: linelist_type
  USE ion_mod, only: particle_type
  ! The one Voigt kernel. qlg_bot_sweep is the window form, used by add_opacity_equi;
  ! qlg_bot_phi is the same profile at a single point, in the convention of the retired
  ! voigt_function2, for the two routines here whose abscissa is not a uniform ramp.
  USE qlg, only: qlg_init, qlg_bot_sweep, qlg_bot_phi
  implicit none
  private

  public :: build_linelist, finalize_linelist, add_opacity, add_opacity_equi, voigt_peak
  public :: add_particle_opacity
  public :: append_line, append_h_line, write_ll_header, write_ll_record, bij_const
  public :: reserve_linelist, print_ll_particle, print_ll_atom, merge_ll_diag
  public :: voigt_kernel_init, merge_intervals, in_any_interval

  ! Interval type: represents a wavelength range [blue_vac, red_vac] in vacuum Angstrom
  type, public :: interval_t
    real(dp) :: blue_vac, red_vac
  end type interval_t

  ! Bij = bij_const * (10**loggf) / nu0  (see init_einstein_coefficients in bbtrans_mod)
  real(dp), parameter :: bij_const = c2_2h * f_conv * h_4pi
  real(dp), parameter :: over_h4pi = 1.0_dp / h_4pi
  real(4),  parameter :: over_pi   = 0.3183098862
  real(4),  parameter :: sqrt_pi   = 1.7724538509
  real(4),  parameter :: ex_cutoff = -85.0     ! skip lines whose log(Rij) underflows
  ! 1/sqrt(pi): the peak of the profile in qlg_bot_phi's convention, at a = v = 0. The
  ! tolerance pre-tests below rely on it being an exact upper bound.
  real(4), parameter  :: voigt_peak = 0.5641895835

  ! --- The Voigt kernel ---------------------------------------------------------------
  ! QLG-bot (HOP/qlg_mod.f90), unconditionally, everywhere. It is the tabulated member of
  ! the QLG family of the Voigt-QLG letter: three basis functions {q.G, q.Z(0.2), q.L^2}
  ! at x-degree 1, six fitted constants, G and Z(0.2)/L^2 read from two rank-1 tables at
  ! one shared snapped index. Synspec's cost for a much smaller flux error -- 5.7e-3 max
  ! relative on the validation domain, which the kernel-vs-abundance study bounded at
  ! 0.006 dex per line and 0.0007 dex per element, i.e. below anything this code resolves.
  !
  ! The comparison against the ~22 rival and QLG-family kernels this used to dispatch
  ! between at runtime -- haw, mclean, moog, siu, synspec, w4, voigt21, the fp64
  ! reference, and the QLG tiers and ablations -- lives on the `voigt-kernel-zoo` branch,
  ! together with the voigt_kernel dispatch chain that selected them.
  !
  ! add_opacity_equi calls qlg_bot_sweep on the whole line window at once, which is where
  ! QLG-bot's advantage lies: the per-line setup is amortised over the window and the
  ! point loop vectorises. add_opacity walks a non-uniform abscissa and calls the
  ! one-point qlg_bot_phi instead.
  !
  ! voigt_kernel IS READ AND IGNORED. It stays only so that the ~900 namelists and the
  ! standalone library's hop_set_option that still set it keep loading; there is nothing
  ! left to select. Delete it once those have been swept.
  integer, public :: voigt_kernel = 0

  ! --- profiling work counters (add_opacity); cheap integer increments ----------------
  logical,    public :: ll_profile = .false.   ! enable counting
  integer(8), public :: ll_nvisit = 0          ! active line-visits (cell x line)
  integer(8), public :: ll_npass  = 0          ! lines passing the ex/strength cutoff
  integer(8), public :: ll_nvoigt = 0          ! Voigt evaluations (centre + wings)
  integer(8), public :: ll_nmarch = 0          ! cursor advance steps

contains

  !=====================================================================================
  ! Parse + set up one per-species block of linelist records into the SoA arrays.
  ! The block's records live in the byte buffer buf; lstart(i)/lend(i) give the first
  ! and last byte of record i (so no per-record string copy is made). Records are
  ! assumed pre-sorted ascending in wavelength. If do_mask, only records with lam_air
  ! in [lam_min, lam_max] are kept (contiguous slice, located by bisection).
  !-------------------------------------------------------------------------------------
  !=====================================================================================
  ! Build QLG-bot's two lookup tables. MUST be called once before any opacity or rate is
  ! accumulated: qlg_bot_sweep and qlg_bot_phi are PURE and cannot check, so an omitted
  ! call is not an error but a silently zero profile. Idempotent, and a few milliseconds.
  !
  ! Entry points that must make it: SRC/composition_patch_mod.f90 after reading
  ! &linelist_params, and the standalone library's hop_init.
  !
  ! The name is kept from when it built the tables of whichever kernel voigt_kernel had
  ! selected; it is now a one-line forward to qlg_init.
  !-------------------------------------------------------------------------------------
  SUBROUTINE voigt_kernel_init
    call qlg_init
  END SUBROUTINE voigt_kernel_init

  SUBROUTINE build_linelist(ll, buf, lstart, lend, n, awgt, ion_j, xion, is_mol, &
                            vbar2e4, do_mask, lam_min, lam_max, intervals)
    type(linelist_type), intent(inout) :: ll
    character(len=*),    intent(in)    :: buf
    integer(8),          intent(in)    :: lstart(:), lend(:)
    integer,             intent(in)    :: n, ion_j
    real(dp),            intent(in)    :: awgt, xion, vbar2e4
    logical,             intent(in)    :: is_mol, do_mask
    real(dp),            intent(in)    :: lam_min, lam_max
    type(interval_t), intent(in), optional :: intervals(:)
    real(dp) :: lam_air, lam_vac, lower_ev, loggf, gv_in, ga_in, gq_in, r_i, r_j
    real(dp) :: nu0, ga_t, gw1_t, gw2_t, gq_t
    real(dp) :: fields(10)
    integer  :: lo, hi, nkeep, i, k, iostat, nvalue
    integer  :: ib, nblk
    integer, allocatable :: blk_lo(:), blk_hi(:)
    integer, save :: it_mask=0, it_parse=0
    !-----------------------------------------------------------------------------------
    if (n <= 0) return

    ! select records based on intervals (multi-window) or single lam_min/lam_max range
    call trace%begin ('ll_build_mask', itimer=it_mask)
    if (do_mask) then
      if (present(intervals)) then
        ! multi-window selection. The file is wavelength-sorted, so bisect out one record
        ! block per interval rather than parsing every record: the intervals are sorted and
        ! disjoint (merge_intervals), so the blocks are too. Each block is widened by one
        ! record on either side and membership is still decided by in_any_interval in
        ! vacuum, so the air-space bisection only bounds the work and never the selection.
        nblk = size(intervals)
        allocate(blk_lo(nblk), blk_hi(nblk))
        nkeep = 0
        do ib = 1, nblk
          lo = bisect_block(buf, lstart, lend, n, vac2air(intervals(ib)%blue_vac))
          hi = bisect_block(buf, lstart, lend, n, vac2air(intervals(ib)%red_vac)) - 1
          if (hi < n) then
            if (parse_lam(buf(lstart(hi+1):lend(hi+1))) <= vac2air(intervals(ib)%red_vac)) &
              hi = hi + 1
          end if
          lo = max(1, lo - 1)
          hi = min(n, hi + 1)
          if (ib > 1) lo = max(lo, blk_hi(ib-1) + 1)   ! never visit a record twice
          blk_lo(ib) = lo
          blk_hi(ib) = hi
          do i = lo, hi
            lam_vac = air2vac(parse_lam(buf(lstart(i):lend(i))))
            if (in_any_interval(lam_vac, intervals)) nkeep = nkeep + 1
          end do
        end do
      else
        ! single-window selection: use contiguous bisection
        ! note: lam_min and lam_max are in vacuum Angstrom
        lo = bisect_block(buf, lstart, lend, n, vac2air(lam_min))
        hi = bisect_block(buf, lstart, lend, n, vac2air(lam_max)) - 1
        ! bisect_block returns the first index with lam >= target; refine the upper edge
        if (hi < n) then
          if (parse_lam(buf(lstart(hi+1):lend(hi+1))) <= vac2air(lam_max)) hi = hi + 1
        end if
        nkeep = hi - lo + 1
      end if
    else
      lo = 1
      hi = n
      nkeep = hi - lo + 1
    end if
    call trace%end(it_mask)
    if (nkeep <= 0) return

    call grow(ll, nkeep)

    call trace%begin ('ll_build_parse', itimer=it_parse)
    if (present(intervals)) then
      ! multi-window: iterate the record blocks found above and check intervals
      k = ll%nline
      do ib = 1, nblk
      do i = blk_lo(ib), blk_hi(ib)
        call parse_record(buf(lstart(i):lend(i)), fields, iostat, nvalue)
        lam_air = fields(1); lower_ev = fields(2); loggf = fields(3)
        gv_in = fields(4); ga_in = fields(5); gq_in = fields(6)
        r_i = fields(7); r_j = fields(8)
        lam_vac = air2vac(lam_air)
        if (.not. in_any_interval(lam_vac, intervals)) cycle

        k = k + 1
        ll%nline = k
        nu0   = HZ_TO_AA / lam_vac
        ll%llam0(k)    = log(lam_vac)
        ll%nu0(k)      = nu0
        ll%lower_ev(k) = lower_ev
        ll%gpref(k)    = loggf*ln10 + log(bij_const / nu0)
        ll%lev_i(k)    = nint(r_i)
        ll%lev_j(k)    = nint(r_j)
        ll%active(k)   = .true.
        ll%used(k)     = .false.
        ll%strength(k) = 0
        ll%reach(k)    = 0
        ll%reach_in(k) = merge(real(fields(10) / lam_vac, 4), -1.0, nvalue >= 10)
        ll%loggf(k)    = loggf
        ll%gv_in(k)    = gv_in
        ll%ga_in(k)    = ga_in
        ll%gq_in(k)    = gq_in

        call prepare_radiative   (ga_t,  ga_in, nu0, loggf)
        call prepare_stark       (gq_t,  gq_in, nu0, lower_ev, ion_j, xion)
        call prepare_vanderwaals (gw1_t, gw2_t, ll%unsold(k), gv_in, nu0, &
                                  lower_ev, ion_j, xion, vbar2e4, is_mol)
        ll%ga(k)  = ga_t
        ll%gq(k)  = gq_t
        ll%gw1(k) = gw1_t
        ll%gw2(k) = gw2_t
      end do
      end do
    else
      ! single-window or no mask: use the contiguous slice [lo:hi]
      do i = lo, hi
        call parse_record(buf(lstart(i):lend(i)), fields, iostat, nvalue)
        lam_air = fields(1); lower_ev = fields(2); loggf = fields(3)
        gv_in = fields(4); ga_in = fields(5); gq_in = fields(6)
        r_i = fields(7); r_j = fields(8)

        k = ll%nline + 1
        ll%nline = k
        lam_vac = air2vac(lam_air)
        nu0   = HZ_TO_AA / lam_vac
        ll%llam0(k)    = log(lam_vac)
        ll%nu0(k)      = nu0
        ll%lower_ev(k) = lower_ev
        ll%gpref(k)    = loggf*ln10 + log(bij_const / nu0)
        ll%lev_i(k)    = nint(r_i)
        ll%lev_j(k)    = nint(r_j)
        ll%active(k)   = .true.
        ll%used(k)     = .false.
        ll%strength(k) = 0
        ll%reach(k)    = 0
        ll%reach_in(k) = merge(real(fields(10) / lam_vac, 4), -1.0, nvalue >= 10)
        ll%loggf(k)    = loggf
        ll%gv_in(k)    = gv_in
        ll%ga_in(k)    = ga_in
        ll%gq_in(k)    = gq_in

        call prepare_radiative   (ga_t,  ga_in, nu0, loggf)
        call prepare_stark       (gq_t,  gq_in, nu0, lower_ev, ion_j, xion)
        call prepare_vanderwaals (gw1_t, gw2_t, ll%unsold(k), gv_in, nu0, &
                                  lower_ev, ion_j, xion, vbar2e4, is_mol)
        ll%ga(k)  = ga_t
        ll%gq(k)  = gq_t
        ll%gw1(k) = gw1_t
        ll%gw2(k) = gw2_t
      end do
    end if
    call trace%end(it_parse)
  END SUBROUTINE build_linelist

  !=====================================================================================
  ! Append a single line with explicit (already prepared) parameters; used to convert
  ! model-atom bound-bound transitions (atom_mod%build_ll) into a linelist_type.
  ! Broadening parameters carry the same over_4pi/nu0 scaling as the prepare_* routines.
  !-------------------------------------------------------------------------------------
  SUBROUTINE append_line(ll, llam0, nu0, lower_ev, gpref, ga, gw1, gw2, gq, unsold, &
                         lev_i, lev_j, itrans, ipart, wA, cwt, dej, iso_out, &
                         loggf, gv_in, ga_in, gq_in)
    type(linelist_type), intent(inout) :: ll
    real(8), intent(in) :: llam0, gpref
    real(8), intent(in) :: nu0, lower_ev, ga, gw1, gw2, gq, wA, cwt, dej
    real(8), intent(in) :: loggf, gv_in, ga_in, gq_in
    logical, intent(in) :: unsold
    integer, intent(in) :: lev_i, lev_j, itrans, ipart, iso_out
    integer :: k
    !-----------------------------------------------------------------------------------
    call grow(ll, 1)
    call grow_ext(ll, 1)
    k = ll%nline + 1
    ll%nline = k
    ll%llam0(k)    = llam0
    ll%nu0(k)      = nu0
    ll%lower_ev(k) = lower_ev
    ll%gpref(k)    = gpref
    ll%ga(k)       = ga
    ll%gw1(k)      = gw1
    ll%gw2(k)      = gw2
    ll%gq(k)       = gq
    ll%unsold(k)   = unsold
    ll%lev_i(k)    = lev_i
    ll%lev_j(k)    = lev_j
    ll%active(k)   = .true.
    ll%used(k)     = .false.
    ll%strength(k) = 0
    ll%reach(k)    = 0
    ll%reach_in(k) = -1
    ll%loggf(k)    = loggf
    ll%gv_in(k)    = gv_in
    ll%ga_in(k)    = ga_in
    ll%gq_in(k)    = gq_in
    ll%itrans(k)   = itrans
    ll%ipart(k)    = ipart
    ll%wA(k)       = wA
    ll%cwt(k)      = cwt
    ll%dej(k)      = dej
    ll%iso_out(k)  = iso_out
    ll%iws(k)      = 0
    ll%iwe(k)      = 0
    ll%icen(k)     = 0
  END SUBROUTINE append_line

  !=====================================================================================
  ! Append one H I line-list line with principal quantum numbers n_lo -> n_up. Only the
  ! identification is stored: add_h_opacities rebuilds the Stark, resonance and van der
  ! Waals widths from (n_lo, n_up, nu0) with hydrogen_line, so ga/gw1/gw2/gq stay 0.
  ! lev_i/lev_j start as n_lo/n_up.
  !-------------------------------------------------------------------------------------
  SUBROUTINE append_h_line(ll, llam0, nu0, lower_ev, loggf, gv_in, ga_in, gq_in, &
                           n_lo, n_up)
    type(linelist_type), intent(inout) :: ll
    real(8), intent(in) :: llam0, nu0, lower_ev, loggf, gv_in, ga_in, gq_in
    integer, intent(in) :: n_lo, n_up
    integer :: k
    !-----------------------------------------------------------------------------------
    call grow(ll, 1)
    call growi(ll%n_lo, ll%nline, ll%nline + 1)
    call growi(ll%n_up, ll%nline, ll%nline + 1)
    k = ll%nline + 1
    ll%nline = k
    ll%llam0(k)    = llam0
    ll%nu0(k)      = nu0
    ll%lower_ev(k) = lower_ev
    ll%gpref(k)    = loggf*ln10 + log(bij_const / nu0)
    ll%ga(k)       = 0
    ll%gw1(k)      = 0
    ll%gw2(k)      = 0
    ll%gq(k)       = 0
    ll%unsold(k)   = .false.
    ll%lev_i(k)    = n_lo
    ll%lev_j(k)    = n_up
    ll%n_lo(k)     = n_lo
    ll%n_up(k)     = n_up
    ll%active(k)   = .true.
    ll%used(k)     = .false.
    ll%strength(k) = 0
    ll%reach(k)    = 0
    ll%reach_in(k) = -1
    ll%loggf(k)    = loggf
    ll%gv_in(k)    = gv_in
    ll%ga_in(k)    = ga_in
    ll%gq_in(k)    = gq_in
  END SUBROUTINE append_h_line

  !=====================================================================================
  ! After all blocks for a particle have been appended, sort every parallel array by
  ! llam0 ascending so the kernel's monotonic march is valid across isotope blocks.
  !-------------------------------------------------------------------------------------
  SUBROUTINE finalize_linelist(ll)
    type(linelist_type), intent(inout) :: ll
    integer, allocatable :: idx(:)
    integer :: i
    integer, save :: it_argsort=0, it_permute=0
    logical :: sorted
    !-----------------------------------------------------------------------------------
    ! lines select by reach only when every one of them carries it
    ll%has_reach = .false.
    if (ll%nline > 0) then
      ll%has_reach = all(ll%reach_in(1:ll%nline) >= 0)
      if (ll%has_reach) ll%reach_in_max = maxval(ll%reach_in(1:ll%nline))
    end if
    if (ll%nline <= 1) return
    ! the common case (one block per particle) is already wavelength-sorted on read,
    ! so an O(n) check lets us skip the sort entirely
    sorted = .true.
    do i = 2, ll%nline
      if (ll%llam0(i) < ll%llam0(i-1)) then
        sorted = .false.
        exit
      end if
    end do
    if (sorted) return

    call trace%begin ('ll_argsort', itimer=it_argsort)
    allocate(idx(ll%nline))
    idx = [(i, i=1, ll%nline)]
    call argsort8(ll%llam0(1:ll%nline), idx)
    call trace%end(it_argsort)

    call trace%begin ('ll_permute', itimer=it_permute)
    call permute8(ll%llam0,    idx)
    call permute4(ll%nu0,      idx)
    call permute4(ll%lower_ev, idx)
    call permute4(ll%gpref,    idx)
    call permute4(ll%loggf,    idx)
    call permute4(ll%gv_in,    idx)
    call permute4(ll%ga_in,    idx)
    call permute4(ll%gq_in,    idx)
    call permute4(ll%ga,       idx)
    call permute4(ll%gw1,      idx)
    call permute4(ll%gw2,      idx)
    call permute4(ll%gq,       idx)
    call permute4(ll%strength, idx)
    call permute4(ll%reach,    idx)
    call permute4(ll%reach_in, idx)
    call permutei(ll%lev_i,    idx)
    call permutei(ll%lev_j,    idx)
    call permutel(ll%unsold,   idx)
    call permutel(ll%active,   idx)
    call permutel(ll%used,     idx)
    if (allocated(ll%n_lo)) then
      call permutei(ll%n_lo,    idx)
      call permutei(ll%n_up,    idx)
    end if
    ! model-atom extensions (bulk linelists leave these unallocated)
    if (allocated(ll%ipart)) then
      call permutei(ll%ipart,   idx)
      call permutei(ll%itrans,  idx)
      call permutei(ll%iso_out, idx)
      call permutei(ll%iws,     idx)
      call permutei(ll%iwe,     idx)
      call permutei(ll%icen,    idx)
      call permute4(ll%wA,      idx)
      call permute4(ll%cwt,     idx)
      call permute4(ll%dej,     idx)
    end if
    call trace%end(it_permute)
  END SUBROUTINE finalize_linelist

  !=====================================================================================
  ! Opacity kernel: spatial-outer, line, frequency-inner; log-lambda march.
  !
  !   chi, eta  : (nf, m1, m2, m3)  accumulated into (the "local" work arrays)
  !   cont      : (nf, m1, m2, m3)  continuum reference for the tolerance window
  !   logref    : (nf, m1, m2, m3)  log(tol*cont/voigt_peak); log-space tolerance pre-test
  !   llam      : (nf)              log(lambda) of the wavelength grid (globally ascending
  !                                 unless ll is windowed, then ascending per window)
  !   inv_df    : (m1,m2,m3)       inverse Doppler width c/v_th (this ll's isotope)
  !   vbar2     : (m1,m2,m3)       mean square velocity (this ll's isotope)
  !   lnbase    : (m1,m2,m3, np)    log(n) - log_part + dxij*e_kT + log(idf) per
  !                                 population group (idf = this ll's isotope Doppler
  !                                 factor; bulk linelist: np=1; model atoms:
  !                                 one group per ion stage, selected by ll%ipart)
  !   per-cell  : vproj, e_kT, ne, tothi, stark_pe
  !   dep_mem   : (m1,m2,m3, nlevel) NLTE departure coefficients (or size 0 if LTE)
  !   do_nlte   : whether dep_mem applies to this particle
  !
  ! The wing cutoff is evaluated on the LTE line opacity (phi*Rij), so the profile
  ! extent does not depend on the departure coefficients (no flip-flop between
  ! NLTE iterations); the departures only scale the added chi/eta.
  !
  ! LINES WITH A REACH (ll%has_reach, a list written by hop_prune.x or a run's dump).
  ! A patch is a chunk of the run's grid. Without reaches every chunk visits every line
  ! of the slice, and a line centred outside the chunk takes the nearest chunk end as its
  ! centre index f0: its tolerance reference is the continuum and accumulated opacity
  ! there, and add_opacity_equi centres its evaluation window on that end point, so the
  ! chunk receives wing beyond the line's extent. A spectrum then depends on where the
  ! chunks fall. With reaches the kernel instead
  !   * visits only lines whose [llam0 - reach, llam0 + reach], widened by the largest
  !     |vproj| and half a grid step, overlaps the chunk (select_lines);
  !   * evaluates a line whose shifted centre lies outside the chunk about its true
  !     centre. Its tolerance reference is tol times the continuum at the nearer chunk
  !     end, which differs from the one at the centre only by the continuum slope over
  !     the distance; the peak test and the analytic extent use it as for a line inside,
  !     and the extent is capped by the reach and clipped to the chunk. The cap only
  !     cuts wing below the tolerance level in the states the reach was recorded from.
  !     The accumulated-opacity floor is not applied: it depends on the lines
  !     processed before in the same chunk, and leaving it out can only widen the window;
  !   * treats a line centred inside the chunk exactly as without reaches.
  ! A list without reaches keeps the whole-slice visit and the nearest-end treatment.
  !-------------------------------------------------------------------------------------
  SUBROUTINE add_opacity(ll, llam, chi, eta, cont, logref, vproj, inv_df, vbar2, e_kT, &
                         lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
                         one_m_eps, chi_scat, do_scat, i0, i1, llam_lines)
    type(linelist_type), intent(inout) :: ll
    real(dp), dimension(:),       intent(in)    :: llam
    real(4),  dimension(:,:,:,:), intent(inout) :: chi, eta
    real(4),  dimension(:,:,:,:), intent(in)    :: cont, logref
    real(4),  dimension(:,:,:),   intent(in)    :: vproj, e_kT, ne, tothi, stark_pe
    real(4),  dimension(:,:,:),   intent(in)    :: inv_df, vbar2
    real(4),  dimension(:,:,:,:), intent(in)    :: lnbase
    real(4),  dimension(:,:,:,:), intent(in)    :: dep_mem
    logical,  intent(in) :: do_nlte
    ! line_scattering: precomputed stage-appropriate scattering fraction 1-eps
    ! on the frequency grid, and the scattering-opacity accumulator (both
    ! size-0 dummies with do_scat=.false. when disabled)
    real(4),  dimension(:,:,:,:), intent(in)    :: one_m_eps
    real(4),  dimension(:,:,:,:), intent(inout) :: chi_scat
    logical,  intent(in) :: do_scat
    integer,  intent(in), optional :: i0, i1
    ! llam_lines: log(lambda) range of the line centres this patch visits. The caller
    ! sets it to the neighbourhood the selection applied to the patch's window when the
    ! list was built. A disjoint multi-window list otherwise lets a foreign line clamp
    ! onto this patch's edge through the forward march, so a window would not reproduce
    ! its own run.
    real(dp), intent(in), optional :: llam_lines(2)
    integer  :: l0, l1, nsel, isel
    integer, allocatable :: sel(:)
    integer  :: f0_prev
    real(4)  :: refl_c, lref_c, flo
    integer  :: i, j, k, l, f, f0, fc, nf, m(3), nl, np, ws, we, fa, fb, lo, hi, io, mem_offset
    real(8)  :: llc, dl, hb(2)
    integer  :: fe
    real(4)  :: idf, unsold_c, ekt, ne_c, tothi_c, starkpe, vproj_c
    real(4)  :: a, ex, Rij, bi, bj, ref, vs, phi, lte, lvb2, vm2, edge, lref_o
    real(4)  :: arg
    real(4)  :: base_s(size(lnbase,4))
    logical  :: nlte, has_ipart, outside
    integer(8) :: c_visit, c_pass, c_voigt, c_march
    !-----------------------------------------------------------------------------------
    m  = shape(vproj)
    nf = size(llam)
    nl = ll%nline
    np = size(lnbase, 4)
    ! patch-local line slice (llam0 is ascending, see finalize_linelist)
    l0 = 1
    l1 = nl
    if (present(llam_lines) .and. nl > 0) then
      l0 = lower_bound8(ll%llam0, nl, llam_lines(1))
      l1 = upper_bound8(ll%llam0, nl, llam_lines(2))
    end if
    if (nl == 0) return
    call select_lines(ll, llam, vproj, l0, l1, sel, nsel, hb)
    has_ipart = allocated(ll%ipart)
    c_visit = 0; c_pass = 0; c_voigt = 0; c_march = 0
    lo = 1;  hi = m(1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    mem_offset = lo - 1
    if (hi - lo + 1 /= m(1)) &
      error stop 'linelist add_opacity: per-cell arrays are not sized to the cell range'
    if (size(chi, 2) /= m(1)) &
      error stop 'linelist add_opacity: chi/eta extent does not match the cell range'
    ! the spectrum machinery always delivers lam (hence llam) ascending; the monotonic
    ! march relies on this, so assert it once rather than branching per line.
    if (nf > 1 .and. llam(nf) < llam(1)) &
      error stop 'linelist add_opacity: log(lambda) grid must be ascending'

    do k = 1, m(3)
      do j = 1, m(2)
        do i = lo, hi
          io       = i - mem_offset      ! chunk-local index; cont/chi_scat/dep_mem stay at i
          f0_prev  = -1                  ! running-floor cache is per cell
          idf      = inv_df(io,j,k)
          unsold_c = vbar2(io,j,k)**0.3 * tothi(io,j,k)
          lvb2     = log(vbar2(io,j,k))                   ! vb2**gw2 as exp(gw2*log(vb2))
          ekt      = e_kT(io,j,k)
          ne_c     = ne(io,j,k)
          tothi_c  = tothi(io,j,k)
          starkpe  = stark_pe(io,j,k)
          vproj_c  = vproj(io,j,k)
          ! line-independent part of log(Rij) per population group (idf folded in)
          base_s(1:np) = lnbase(io,j,k,1:np)
          fc = 1
          ws = 1
          we = nf

          do isel = 1, nsel
            l = sel(isel)
            c_visit = c_visit + 1

            ! shifted line centre in log(lambda)
            llc = ll%llam0(l) - vproj_c
            ! monotonic in l (vproj_c fixed), so the cursor fc only advances
            ! forward through the (ascending) grid.
            do while (fc < nf)
              if (llam(fc) >= llc) exit
              fc = fc + 1
              c_march = c_march + 1
            end do
            if (fc > 1) then
              if (abs(llam(fc) - llc) < abs(llam(fc-1) - llc)) then
                f0 = fc
              else
                f0 = fc - 1
              end if
            else
              f0 = 1
            end if

            ! log(Rij) for the main isotope = base + gpref - ex.
            ex  = ll%lower_ev(l) * ekt
            if (has_ipart) then
              arg = base_s(ll%ipart(l)) + ll%gpref(l) - ex
            else
              arg = base_s(1) + ll%gpref(l) - ex
            end if
            if (arg < ex_cutoff) cycle

            outside = ll%has_reach .and. (llc < llam(1) .or. llc > llam(nf))
            if (outside) then
              ! centred beyond this chunk: the tolerance reference is the continuum at
              ! the nearer chunk end, without the accumulated-opacity floor (see above)
              fe = merge(1, nf, llc < llam(1))
              lref_o = logref(fe,io,j,k)
              if (arg <= lref_o) cycle
              Rij = exp(arg)
              ref = tol * cont(fe,i,j,k)
            else
              ! log-space tolerance pre-test: voigt_peak*Rij <= ref  <=>  arg <= lref_c.
              ! This rejects the vast majority of lines without evaluating exp().
              ! The reference is the continuum times tol, floored by the ALREADY
              ! ACCUMULATED line opacity: tol*tol*chi, so that narrow overlapping lines
              ! cannot cancel broader weaker ones. line_reject_frac (default 0) raises
              ! that floor further. Cached on f0 -- but chi is written by this very loop,
              ! so any line that actually deposits opacity invalidates the cache below.
              if (f0 /= f0_prev) then
                f0_prev = f0
                refl_c  = tol * cont(f0,i,j,k)
                flo     = max(tol, line_reject_frac) * tol * chi(f0,io,j,k)
                if (flo > refl_c) then
                  refl_c = flo
                  lref_c = log(flo * sqrt_pi)      ! == log(ref/voigt_peak)
                else
                  lref_c = logref(f0,io,j,k)
                end if
              end if

              if (arg <= lref_c) cycle

              ! this ll's isotope Doppler factor idf is folded into lnbase (hence arg)
              Rij = exp(arg)
              ref = refl_c

              ! exact pre-filter: the profile <= 1/sqrt(pi) everywhere, so if even the
              ! peak cannot reach the tolerance the whole line is negligible -- skip
              ! the damping and Voigt evaluations entirely (most lines drop out here).
              if (voigt_peak * Rij <= ref) cycle
            end if
            c_pass = c_pass + 1

            nlte = do_nlte .and. ll%lev_i(l) > 0 .and. ll%lev_j(l) > 0
            if (nlte) then
              bi = dep_mem(i,j,k, ll%lev_i(l))
              bj = dep_mem(i,j,k, ll%lev_j(l))
            else
              bi = 1.0
              bj = 1.0
            end if

            ! damping parameter (physical, then scaled by idf)
            a = ll%ga(l)
            if (ll%unsold(l)) then
              a = a + ll%gw1(l) * unsold_c
            else
              a = a + ll%gw1(l) * exp(ll%gw2(l) * lvb2) * tothi_c
            end if
            if (ll%gq(l) >= 0.0) then
              a = a + ll%gq(l) * ne_c
            else
              a = a - ll%gq(l) * starkpe
            end if
            a = a * idf

            if (outside) then
              vm2 = max(arg - lref_o, a*Rij*over_pi/ref)
              dl  = sqrt(vm2) / idf
              ll%reach(l) = max(ll%reach(l), real(dl, 4))
              dl  = min(dl, real(ll%reach_in(l), 8))
              if (llc + dl < llam(1) - hb(1) .or. llc - dl > llam(nf) + hb(2)) cycle
              ll%used(l) = .true.
              fa = bisect_nearest(llam, ws, we, llc - dl)
              fb = bisect_nearest(llam, ws, we, llc + dl)
            else
              ll%used(l) = .true.
              ll%strength(l) = max(ll%strength(l), voigt_peak*Rij / cont(f0,i,j,k))

              ! analytic profile extent where phi*Rij drops to ref: Gaussian core
              ! v^2 = ln(peak/ref) = arg - logref, Lorentz wing v^2 = a*Rij/(pi*ref);
              ! a straight (data-independent) loop over [fa,fb], centre included,
              ! pipelines much better than a walk whose exit condition depends on
              ! the previous Voigt value. The unclipped half-width is the line's reach.
              vm2  = max(arg - lref_c, a*Rij*over_pi/ref)
              dl   = sqrt(vm2) / idf
              ll%reach(l) = max(ll%reach(l), real(dl, 4))
              edge = idf * real(max(llam(we) - llc, llc - llam(ws)))
              if (vm2 >= edge*edge) then
                fa = ws         ! line fills the window; skip the extent bisections
                fb = we
              else
                fa = bisect_nearest(llam, ws, f0, llc - dl)
                fb = bisect_nearest(llam, f0, we, llc + dl)
              end if
            end if
            c_voigt = c_voigt + fb - fa + 1

            if (do_scat .and. .not. nlte) then
              ! LTE line (bi=bj=1) with the scattering fraction accumulated
              do f = fa, fb
                vs   = idf * real(llam(f) - llc)
                phi  = qlg_bot_phi(a, vs)
                lte  = phi * Rij
                chi(f,io,j,k) = chi(f,io,j,k) + lte
                eta(f,io,j,k) = eta(f,io,j,k) + lte
                chi_scat(f,i,j,k) = chi_scat(f,i,j,k) + one_m_eps(f,i,j,k)*lte
              end do
            else
              do f = fa, fb
                vs   = idf * real(llam(f) - llc)
                phi  = qlg_bot_phi(a, vs)
                lte  = phi * Rij
                chi(f,io,j,k) = chi(f,io,j,k) + lte*bi
                eta(f,io,j,k) = eta(f,io,j,k) + lte*bj
              end do
            end if
            f0_prev = -1        ! chi just changed; the running floor is stale
          end do
        end do
      end do
    end do
    if (ll_profile) then
      !$omp atomic
      ll_nvisit = ll_nvisit + c_visit
      !$omp atomic
      ll_npass  = ll_npass  + c_pass
      !$omp atomic
      ll_nvoigt = ll_nvoigt + c_voigt
      !$omp atomic
      ll_nmarch = ll_nmarch + c_march
    endif
  END SUBROUTINE add_opacity

  !=====================================================================================
  ! add_opacity variant for a grid that is EQUIDISTANT in log(lambda) (constant-R grid;
  ! caller dispatches on spectrum%equidistant). Same physics as add_opacity, but the
  ! nearest grid index and both wing limits are pure arithmetic -- no forward march and
  ! no bisection are needed, because the grid spacing dllam is constant. The Voigt window
  ! limits come from inverting the far-wing approximation (Gaussian core vs. Lorentz wing,
  ! v = sqrt(vm2)/idf in log(lambda)) and converting that half-width straight to an index
  ! offset via 1/dllam. Windowed (per-line ll%iws/iwe) mode is not supported here.
  !-------------------------------------------------------------------------------------
  SUBROUTINE add_opacity_equi(ll, llam, chi, eta, cont, logref, vproj, inv_df, vbar2, &
                              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
                              one_m_eps, chi_scat, do_scat, i0, i1, llam_lines)
    type(linelist_type), intent(inout) :: ll
    real(dp), dimension(:),       intent(in)    :: llam
    real(4),  dimension(:,:,:,:), intent(inout) :: chi, eta
    real(4),  dimension(:,:,:,:), intent(in)    :: cont, logref
    real(4),  dimension(:,:,:),   intent(in)    :: vproj, e_kT, ne, tothi, stark_pe
    real(4),  dimension(:,:,:),   intent(in)    :: inv_df, vbar2
    real(4),  dimension(:,:,:,:), intent(in)    :: lnbase
    real(4),  dimension(:,:,:,:), intent(in)    :: dep_mem
    logical,  intent(in) :: do_nlte
    ! line_scattering: see add_opacity
    real(4),  dimension(:,:,:,:), intent(in)    :: one_m_eps
    real(4),  dimension(:,:,:,:), intent(inout) :: chi_scat
    logical,  intent(in) :: do_scat
    integer,  intent(in), optional :: i0, i1
    ! llam_lines: log(lambda) range of the line centres this patch visits. The caller
    ! sets it to the neighbourhood the selection applied to the patch's window when the
    ! list was built. A disjoint multi-window list otherwise lets a foreign line clamp
    ! onto this patch's edge through the forward march, so a window would not reproduce
    ! its own run.
    real(dp), intent(in), optional :: llam_lines(2)
    integer  :: l0, l1, nsel, isel
    integer, allocatable :: sel(:)
    integer  :: f0_prev
    real(4)  :: refl_c, lref_c, flo
    integer  :: i, j, k, l, f, f0, nf, m(3), nl, np, fa, fb, noff, nw, lo, hi, io, mem_offset
    real(8)  :: llc, llam1, dllam, inv_dllam, hb(2), dl
    integer  :: fe
    real(4)  :: idf, unsold_c, ekt, ne_c, tothi_c, starkpe, vproj_c
    real(4)  :: a, ex, Rij, bi, bj, ref, vs, phi, lte, lvb2, vm2, lref_o
    real(4)  :: arg
    real(4)  :: base_s(size(lnbase,4))
    logical  :: nlte, has_ipart, outside
    integer(8) :: c_visit, c_pass, c_voigt, c_fill
    ! kernel scratch: the window's K values and the fp32 abscissa ramp. Kbuf is
    ! an automatic array (stack) but only [1:nw] is ever touched, and windows average
    ! ~8 points, so it stays L1-resident.
    real(4)  :: Kbuf(size(llam)), v0, dv, Rn, rnoff, nfac
    !-----------------------------------------------------------------------------------
    m  = shape(vproj)
    nf = size(llam)
    nl = ll%nline
    np = size(lnbase, 4)
    ! patch-local line slice (llam0 is ascending, see finalize_linelist)
    l0 = 1
    l1 = nl
    if (present(llam_lines) .and. nl > 0) then
      l0 = lower_bound8(ll%llam0, nl, llam_lines(1))
      l1 = upper_bound8(ll%llam0, nl, llam_lines(2))
    end if
    if (nl == 0) return
    if (nf < 2) return
    call select_lines(ll, llam, vproj, l0, l1, sel, nsel, hb)
    has_ipart = allocated(ll%ipart)
    c_visit = 0; c_pass = 0; c_voigt = 0; c_fill = 0
    lo = 1;  hi = m(1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    mem_offset = lo - 1
    if (hi - lo + 1 /= m(1)) &
      error stop 'linelist add_opacity_equi: per-cell arrays are not sized to the cell range'
    if (size(chi, 2) /= m(1)) &
      error stop 'linelist add_opacity_equi: chi/eta extent does not match the cell range'
    if (llam(nf) < llam(1)) &
      error stop 'linelist add_opacity_equi: log(lambda) grid must be ascending'

    ! constant grid spacing in log(lambda); index = nint((llc - llam1)/dllam) + 1
    llam1     = llam(1)
    dllam     = (llam(nf) - llam(1)) / real(nf - 1, 8)
    inv_dllam = 1.0_dp / dllam

    do k = 1, m(3)
      do j = 1, m(2)
        do i = lo, hi
          io       = i - mem_offset      ! chunk-local index; cont/chi_scat/dep_mem stay at i
          f0_prev  = -1                  ! running-floor cache is per cell
          idf      = inv_df(io,j,k)
          ! CELL-INVARIANT PARTS OF THE PER-LINE WINDOW SETUP, hoisted out of the line
          ! loop. dv is only a multiply the compiler is free to lift anyway, but nfac
          ! replaces a DIVISION by idf that it is not: gfortran will not turn x/idf into
          ! x*(1/idf) without -ffast-math, so `sqrt(vm2)/idf*inv_dllam` was a
          ! reciprocal-latency chain on every line that survived the tolerance filters.
          ! Every kernel pays it, so this moves the floor rather than the ranking.
          nfac     = real(inv_dllam) / idf
          dv       = idf * real(dllam)
          unsold_c = vbar2(io,j,k)**0.3 * tothi(io,j,k)
          lvb2     = log(vbar2(io,j,k))                   ! vb2**gw2 as exp(gw2*log(vb2))
          ekt      = e_kT(io,j,k)
          ne_c     = ne(io,j,k)
          tothi_c  = tothi(io,j,k)
          starkpe  = stark_pe(io,j,k)
          vproj_c  = vproj(io,j,k)
          ! line-independent part of log(Rij) per population group (idf folded in)
          base_s(1:np) = lnbase(io,j,k,1:np)

          do isel = 1, nsel
            l = sel(isel)
            c_visit = c_visit + 1

            ! shifted line centre in log(lambda); nearest index by arithmetic
            llc = ll%llam0(l) - vproj_c
            f0  = nint((llc - llam1) * inv_dllam) + 1
            f0  = min(nf, max(1, f0))

            ! log(Rij) for the main isotope = base + gpref - ex.
            ex  = ll%lower_ev(l) * ekt
            if (has_ipart) then
              arg = base_s(ll%ipart(l)) + ll%gpref(l) - ex
            else
              arg = base_s(1) + ll%gpref(l) - ex
            end if
            if (arg < ex_cutoff) cycle

            outside = ll%has_reach .and. (llc < llam(1) .or. llc > llam(nf))
            if (outside) then
              ! centred beyond this chunk: reference at the nearer chunk end (add_opacity)
              fe = merge(1, nf, llc < llam(1))
              lref_o = logref(fe,io,j,k)
              if (arg <= lref_o) cycle
              Rij = exp(arg)
              ref = tol * cont(fe,i,j,k)
            else
              ! log-space tolerance pre-test: voigt_peak*Rij <= ref  <=>  arg <= lref_c.
              ! This rejects the vast majority of lines without evaluating exp().
              ! The reference is the continuum times tol, floored by the ALREADY
              ! ACCUMULATED line opacity: tol*tol*chi, so that narrow overlapping lines
              ! cannot cancel broader weaker ones. line_reject_frac (default 0) raises
              ! that floor further. Cached on f0 -- but chi is written by this very loop,
              ! so any line that actually deposits opacity invalidates the cache below.
              if (f0 /= f0_prev) then
                f0_prev = f0
                refl_c  = tol * cont(f0,i,j,k)
                flo     = max(tol, line_reject_frac) * tol * chi(f0,io,j,k)
                if (flo > refl_c) then
                  refl_c = flo
                  lref_c = log(flo * sqrt_pi)      ! == log(ref/voigt_peak)
                else
                  lref_c = logref(f0,io,j,k)
                end if
              end if

              if (arg <= lref_c) cycle

              ! this ll's isotope Doppler factor idf is folded into lnbase (hence arg)
              Rij = exp(arg)
              ref = refl_c

              ! exact pre-filter: the profile <= 1/sqrt(pi) everywhere, so if even the
              ! peak cannot reach the tolerance the whole line is negligible -- skip
              ! the damping and Voigt evaluations entirely (most lines drop out here).
              if (voigt_peak * Rij <= ref) cycle
            end if
            c_pass = c_pass + 1

            nlte = do_nlte .and. ll%lev_i(l) > 0 .and. ll%lev_j(l) > 0
            if (nlte) then
              bi = dep_mem(i,j,k, ll%lev_i(l))
              bj = dep_mem(i,j,k, ll%lev_j(l))
            else
              bi = 1.0
              bj = 1.0
            end if

            ! damping parameter (physical, then scaled by idf)
            a = ll%ga(l)
            if (ll%unsold(l)) then
              a = a + ll%gw1(l) * unsold_c
            else
              a = a + ll%gw1(l) * exp(ll%gw2(l) * lvb2) * tothi_c
            end if
            if (ll%gq(l) >= 0.0) then
              a = a + ll%gq(l) * ne_c
            else
              a = a - ll%gq(l) * starkpe
            end if
            a = a * idf

            ! WINDOW EXTENT. Analytic, where phi*Rij drops to ref, in v-units: Gaussian
            ! core v^2 = arg - logref, Lorentz far wing v^2 = a*Rij/(pi*ref). On a
            ! uniform log(lambda) grid the half-width in llam is sqrt(vm2)/idf, so the
            ! index half-width is that over dllam -- no bisection required. nfac is
            ! inv_dllam/idf, hoisted to the cell loop above.
            !
            ! THE ROUNDING IS nint, NOT int()+1, AND nw = 1 IS REACHABLE. f0 is the
            ! NEAREST bin to the centre (nint, above), so the centre sits at some offset
            ! d in [-0.5,0.5] bins from f0 and bin f0+n is |n+d| bins out. A symmetric
            ! window covers every bin within W = sqrt(vm2)*nfac iff noff >= W + |d|, and
            ! nint(W) = floor(W + 0.5) is the smallest integer that satisfies that for
            ! every d -- tight and still conservative. The old floor(W)+1 exceeded it by
            ! a bin per side whenever frac(W) < 0.5, and could never go below nw = 3: a
            ! line whose profile dies within a twentieth of a bin was evaluated at three
            ! points, two of them provably under tolerance. Sub-half-bin lines now cost
            ! the one centre evaluation they are worth. Every kernel sweep is `do f=1,n`
            ! and handles n = 1.
            !
            ! One thing this gives up: ref is sampled at f0 only (cont/chi vary across
            ! the window), and the old +1 was an incidental one-bin cushion against that.
            ! It is dropped exactly where it matters least -- narrow windows, over which
            ! cont is flat.
            !
            ! The rnoff guard is not cosmetic: nint() of a half-width that overruns the
            ! grid overflows the integer, so the "line fills the window" case is taken
            ! arithmetically rather than through the clamps.
            if (outside) then
              ! the window is the extent at the chunk-end reference, capped by the reach
              ! and rounded to the nearest grid index like f0
              vm2 = max(arg - lref_o, a*Rij*over_pi/ref)
              dl  = sqrt(vm2) / idf
              ll%reach(l) = max(ll%reach(l), real(dl, 4))
              dl  = min(dl, real(ll%reach_in(l), 8))
              fa  = max(1,  nint((llc - dl - llam1) * inv_dllam) + 1)
              fb  = min(nf, nint((llc + dl - llam1) * inv_dllam) + 1)
              if (fa > fb) cycle
              ll%used(l) = .true.
            else
              ll%used(l) = .true.
              ll%strength(l) = max(ll%strength(l), voigt_peak*Rij / cont(f0,i,j,k))
              vm2  = max(arg - lref_c, a*Rij*over_pi/ref)
              rnoff = sqrt(vm2) * nfac
              ! the unclipped half-width in log(lambda) is the line's reach
              ll%reach(l) = max(ll%reach(l), rnoff * real(dllam, 4))
              if (rnoff >= real(nf, 4)) then
                fa = 1
                fb = nf
                c_fill = c_fill + 1
              else
                noff = nint(rnoff)
                fa   = max(1,  f0 - noff)
                fb   = min(nf, f0 + noff)
              end if
            end if
            c_voigt = c_voigt + fb - fa + 1
            ! Branchless fp32 kernel. The whole window is evaluated first, so the
            ! profile loop and the accumulation loop each vectorise; K stays in L1.
            ! The abscissa is a pure fp32 ramp -- on a uniform log(lambda) grid
            ! llam(f) - llc = (llam(fa) - llc) + (f-fa)*dllam exactly, so the real(8)
            ! llam stream drops out of the inner loop altogether. dv is cell-invariant
            ! and hoisted with nfac.
            nw = fb - fa + 1
            v0 = idf * real(llam(fa) - llc)
            call qlg_bot_sweep (a, v0, dv, nw, Kbuf)
            ! qlg_bot_sweep returns K = H(a,v) in the K(0,0) = 1 normalisation, so
            ! voigt_peak = 1/sqrt(pi) converts it to the per-Doppler-width profile the
            ! accumulation below wants -- the historical Rij/sqrt(pi).
            Rn = Rij * voigt_peak
            if (do_scat .and. .not. nlte) then
              do f = fa, fb
                lte = Kbuf(f-fa+1) * Rn
                chi(f,io,j,k) = chi(f,io,j,k) + lte
                eta(f,io,j,k) = eta(f,io,j,k) + lte
                chi_scat(f,i,j,k) = chi_scat(f,i,j,k) + one_m_eps(f,i,j,k)*lte
              end do
            else
              do f = fa, fb
                lte = Kbuf(f-fa+1) * Rn
                chi(f,io,j,k) = chi(f,io,j,k) + lte*bi
                eta(f,io,j,k) = eta(f,io,j,k) + lte*bj
              end do
            end if
            f0_prev = -1        ! chi just changed; the running floor is stale
          end do
        end do
      end do
    end do
    if (ll_profile) then
      !$omp atomic
      ll_nvisit = ll_nvisit + c_visit
      !$omp atomic
      ll_npass  = ll_npass  + c_pass
      !$omp atomic
      ll_nvoigt = ll_nvoigt + c_voigt
    endif
  END SUBROUTINE add_opacity_equi

  !=====================================================================================
  ! Per-particle SoA opacity: relocated here (from composition_mod) is the isotope loop
  ! that turns a particle's thermodynamic state into the kernel's population base and
  ! Doppler width, then dispatches to add_opacity[_equi]. This is what lets linelist_mod
  ! consume ion_mod's particle_type directly. The caller (composition_mod) still owns the
  ! patch thermodynamic arrays and the departure-coefficient lookup, and passes the
  ! resolved dep_ev/dep_mem in (a size-1 zero dep_mem and do_nlte=.false. in the LTE case).
  !
  !   ptc          : the ion or molecule (polymorphic); its ll(0:niso) are the per-isotope
  !                  SoA line lists (index 0 = main isotope)
  !   temp, vmic   : (m1,m2,m3) temperature and microturbulence (Doppler width inputs)
  !   dxi          : (m1,m2,m3) continuum lowering per ion stage
  !   e_kT,ne,tothi,stark_pe,vproj : (m1,m2,m3) per-cell kernel inputs
  !   chi, eta     : (nf,m1,m2,m3) accumulated line opacity/emissivity
  !   cont, logref : (nf,m1,m2,m3) continuum + log-space tolerance reference
  !   llam         : (nf) log(lambda) grid
  !   dep_ev       : level energies of the departure set (used only when do_nlte)
  !   dep_mem      : (m1,m2,m3,nlevel) departure coefficients (dummy (.,.,.,1) when LTE)
  !-------------------------------------------------------------------------------------
  SUBROUTINE add_particle_opacity(ptc, temp, vmic, dxi, llam, chi, eta, cont, logref, &
                                  vproj, e_kT, ne, tothi, stark_pe, dep_ev, dep_mem, &
                                  do_nlte, equidistant, one_m_eps, chi_scat, i0, i1, &
                                  population, population_scale, llam_lines)
    class(particle_type), intent(inout) :: ptc
    real(4),  dimension(:,:,:),   intent(in)    :: temp, vmic, vproj, e_kT, ne, tothi, stark_pe
    real(dp), dimension(:,:,:),   intent(in)    :: dxi
    real(dp), dimension(:),       intent(in)    :: llam
    real(4),  dimension(:,:,:,:), intent(inout) :: chi, eta
    real(4),  dimension(:,:,:,:), intent(in)    :: cont, logref
    real(8),  dimension(:),       intent(in)    :: dep_ev
    real(4),  dimension(:,:,:,:), intent(in)    :: dep_mem
    logical,  intent(in) :: do_nlte, equidistant
    ! line_scattering: stage-appropriate scattering fraction (already the slice
    ! for this particle's stage) and the shared scattering-opacity accumulator
    real(4), dimension(:,:,:,:), intent(in),    optional :: one_m_eps
    real(4), dimension(:,:,:,:), intent(inout), optional :: chi_scat
    integer, intent(in), optional :: i0, i1
    ! passed straight to the kernels: this patch's line neighbourhood in log(lambda)
    real(dp), intent(in), optional :: llam_lines(2)
    real(dp), dimension(:,:,:), pointer, intent(in), optional :: population
    real(dp), intent(in), optional :: population_scale
    real(4),  allocatable :: idf3(:,:,:), vb3(:,:,:), lnbase(:,:,:,:)
    real(4)  :: eps_dummy(0,0,0,0), scat_dummy(0,0,0,0)
    real(dp) :: iso_frac
    integer  :: m(3), kiso, ll_n, l, lo, hi
    logical  :: do_scat, use_population
    !-----------------------------------------------------------------------------------
    m = shape(temp)
    lo = 1;  hi = size(ptc%n, 1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    if (hi - lo + 1 /= m(1)) error stop 'add_particle_opacity: cell range does not match the per-cell arrays'
    allocate(idf3(m(1),m(2),m(3)), vb3(m(1),m(2),m(3)), lnbase(m(1),m(2),m(3),1))
    do_scat = present(one_m_eps) .and. present(chi_scat)
    use_population = .false.
    if (present(population)) use_population = associated(population)

    ! one kernel call per isotope: its own line list ll(kiso), Doppler width, mean
    ! velocity and isotope fraction (the fraction is applied here, per isotope, rather
    ! than folded into each line's gpref -- see build_linelist)
    do kiso = 0, ubound(ptc%ll, 1)
      associate(sll => ptc%ll(kiso))
        ll_n = sll%nline
        if (ll_n == 0) cycle

        if (kiso == 0) then
          idf3 = ptc%inverse_doppler_factor(temp, vmic)
          vb3  = ptc%compute_vbar2(temp)
          iso_frac = 1
        else
          idf3 = ptc%isotope(kiso)%inverse_doppler_factor(temp, vmic)
          vb3  = ptc%isotope(kiso)%compute_vbar2(temp)
          iso_frac = ptc%isotope(kiso)%frac
        end if

        ! ptc%part is the log partition function; lnbase carries the density, the continuum
        ! lowering, this isotope's Doppler factor and its fraction. Floor the density before
        ! the log: molecules are legitimately absent (n==0) where they do not form, and
        ! log(0)=-Inf raises the divide-by-zero FP trap under full_debug. The floored value
        ! gives a large negative lnbase -> exp(lnbase)=0 downstream, i.e. no opacity.
        if (use_population) then
          lnbase(:,:,:,1) = log(max(population(lo:hi,:,:), tiny(1.0_dp)))
        else
          lnbase(:,:,:,1) = log(max(ptc%n(lo:hi,:,:), tiny(1.0_dp)))
          if (present(population_scale)) &
            lnbase(:,:,:,1) = lnbase(:,:,:,1) + real(log(population_scale),4)
        endif
        lnbase(:,:,:,1) = lnbase(:,:,:,1) - ptc%part(lo:hi,:,:) &
                        + real(dxi) * (ptc%j-1) * e_kT + log(idf3) &
                        + real(log(iso_frac), 4)

        if (do_nlte) then
          do l = 1, ll_n
            sll%lev_i(l) = find_closest(dep_ev, 100*ptc%j + real(sll%lower_ev(l), 8))
            sll%lev_j(l) = find_closest(dep_ev, 100*ptc%j + real(sll%lower_ev(l), 8) &
                                                + H_E*sll%nu0(l))
          end do
        end if

        ! arithmetic kernel on a uniform log(lambda) grid, else the general march
        if (equidistant) then
          if (do_scat) then
            call add_opacity_equi(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              one_m_eps, chi_scat, .true., i0, i1, llam_lines)
          else
            call add_opacity_equi(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              eps_dummy, scat_dummy, .false., i0, i1, llam_lines)
          end if
        else
          if (do_scat) then
            call add_opacity(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              one_m_eps, chi_scat, .true., i0, i1, llam_lines)
          else
            call add_opacity(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              eps_dummy, scat_dummy, .false., i0, i1, llam_lines)
          end if
        end if
      end associate
    end do
  END SUBROUTINE add_particle_opacity

  !=====================================================================================
  ! index of the grid point in llam(ws:we) (ascending) closest to llc
  !-------------------------------------------------------------------------------------
  PURE FUNCTION bisect_nearest(llam, ws, we, llc) result(f0)
    real(dp), dimension(:), intent(in) :: llam
    integer,  intent(in) :: ws, we
    real(8),  intent(in) :: llc
    integer :: f0, lo, hi, mid
    !-----------------------------------------------------------------------------------
    lo = ws
    hi = we
    do while (lo < hi)                 ! first index with llam >= llc
      mid = (lo + hi) / 2
      if (llam(mid) < llc) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    f0 = lo
    if (f0 > ws) then
      if (abs(llam(f0-1) - llc) < abs(llam(f0) - llc)) f0 = f0 - 1
    end if
  END FUNCTION bisect_nearest

  !=====================================================================================
  ! The lines of ll(l0:l1) a kernel call visits, ascending, in sel(1:nsel). Without
  ! reaches that is every active line of the slice. With reaches (see add_opacity) it is
  ! the active lines whose reach window can touch the chunk llam(1:nf): the slice is
  ! first narrowed by bisection on the centres with the list's largest reach as margin,
  ! then each line is compared once with its own reach, so the depth loop never visits
  ! a line that cannot contribute. Both margins include the largest |vproj| and half a
  ! grid step at each end, hb, which the kernels' nearest-point rounding can add.
  !
  ! A per-chunk list is rebuilt on every call rather than stored: it costs one pass over
  ! the narrowed slice, a small fraction of one depth point's visit, and needs no state
  ! tied to how a run's patches divide the grid.
  !-------------------------------------------------------------------------------------
  SUBROUTINE select_lines(ll, llam, vproj, l0, l1, sel, nsel, hb)
    type(linelist_type), intent(in) :: ll
    real(dp), intent(in) :: llam(:)
    real(4),  intent(in) :: vproj(:,:,:)
    integer,  intent(inout) :: l0, l1
    integer, allocatable, intent(out) :: sel(:)
    integer,  intent(out) :: nsel
    real(8),  intent(out) :: hb(2)
    real(8)  :: vmax, blue, red
    integer  :: l, nf
    !-----------------------------------------------------------------------------------
    nf = size(llam)
    hb = 0
    if (nf > 1) hb = 0.5d0 * [llam(2) - llam(1), llam(nf) - llam(nf-1)]
    if (ll%has_reach) then
      vmax = 0
      if (size(vproj) > 0) vmax = maxval(abs(vproj))
      blue = llam(1)  - hb(1) - vmax
      red  = llam(nf) + hb(2) + vmax
      l0 = max(l0, lower_bound8(ll%llam0, ll%nline, blue - ll%reach_in_max))
      l1 = min(l1, upper_bound8(ll%llam0, ll%nline, red  + ll%reach_in_max))
    end if
    allocate(sel(max(0, l1 - l0 + 1)))
    nsel = 0
    do l = l0, l1
      if (.not. ll%active(l)) cycle
      if (ll%has_reach) then
        if (ll%llam0(l) + ll%reach_in(l) < blue .or. ll%llam0(l) - ll%reach_in(l) > red) cycle
      end if
      nsel = nsel + 1
      sel(nsel) = l
    end do
  END SUBROUTINE select_lines

  !=====================================================================================
  ! Ported broadening setup (scalar versions of the bbtrans_mod prepare_* routines).
  ! Results are scaled by over_4pi/nu0, as there.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE prepare_radiative(ga, ga_in, nu0, loggf)
    real(dp), intent(out) :: ga
    real(dp), intent(in)  :: ga_in, nu0, loggf
    !-----------------------------------------------------------------------------------
    if (ga_in == 0) then
      ! as in Korg, g_u A_ul; see bbtrans_mod prepare_Radiative for the other codes
      ga = f_conv * 10**loggf * nu0**2
    else
      ga = ga_in                      ! radiative damping rate, 1/s
    end if
    ga = ga * over_4pi / nu0
    if (disable_ga) ga = 0
  END SUBROUTINE prepare_radiative

  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE prepare_stark(gq, gq_in, nu0, lower_ev, j, cnt_ev)
    real(dp), intent(out) :: gq
    real(dp), intent(in)  :: gq_in, nu0, lower_ev, cnt_ev
    integer,  intent(in)  :: j
    real(dp), parameter :: scale = 1e4**(-one_sixth)   ! (T/1e4 K)**(1/6), with T**(1/6) in stark_pe
    real(dp), parameter :: c19 = log10(39. * (8. / (PI*EM*BK**5))**one_sixth) ! ~19.4
    real(dp) :: upper_ev, c1
    !-----------------------------------------------------------------------------------
    gq = 0
    if (gq_in < -10) then
      gq = - BK * 10**(c19 + two_third*gq_in)   ! recipe c4 SIU (Gray 2005, eqs. 11.26-11.27); sign marks recipe
    elseif (gq_in < 0) then
      gq = - 10**gq_in * scale                  ! sign marks recipe
    elseif (gq_in == 0) then
      upper_ev = lower_ev + H_E * nu0
      if (lower_ev < cnt_ev) then
        ! Cowley (1971) eqs. 5 and 6 at 1e4 K, which he recommends using at all T;
        ! c1 is n*^4 of the upper level (his eq. 4)
        c1 = (j**2 * RydV / (cnt_ev - upper_ev))**2
        if (j == 1) then
          gq = 2.26e-7 * c1
        else
          gq = 5.42e-7 * c1 / (j+1)**2
        end if
      end if
    else
      gq = gq_in
    end if
    gq = gq * over_4pi / nu0
    if (disable_stark) gq = 0
  END SUBROUTINE prepare_stark

  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE prepare_vanderwaals(gw1, gw2, unsold, gv, nu0, lower_ev, lower_ion, &
                                      cnt_ev, vbar2e4, is_mol)
    real(dp), intent(out) :: gw1, gw2
    logical,  intent(out) :: unsold
    real(dp), intent(in)  :: gv, nu0, lower_ev, cnt_ev, vbar2e4
    integer,  intent(in)  :: lower_ion
    logical,  intent(in)  :: is_mol
    real(dp) :: upper_ev, sigma, alpha, gx, gammaf, balfac, c1, c625
    ! uns = 17*C_U**0.4 (CGS, energies in eV). Unsold (1955, eqs. 82.48 and 82.55) gives
    ! gamma_6 = 17 C6**0.4 v**0.6 N with C6 = 1.61e-33*(13.5*Z)**2*[...], i.e. C_U = 2.93e-31
    ! and uns = 1.041e-11; Gray (2005, eq. 11.30) uses C_U = 3.0e-31. The factors below
    ! correspond to C_U = 3.22e-31, which makes these widths 3.7% larger than Unsold's.
    real(dp), parameter :: uns = 8.411 * 1.283984e-12
    real(dp), parameter :: barklem = (1e-2 * Rbohr)**2
    !-----------------------------------------------------------------------------------
    gw1 = 0
    gw2 = 0
    unsold = .false.
    if (gv >= 20.0) then
      ! Barklem cross section + velocity exponent
      sigma  = barklem * int(gv)
      alpha  = gv - int(gv)
      gx     = 1.0 - 0.5*alpha
      gammaf = 1 + (-.5748646 + (.9512363 + (-.6998588 + (.4245549 - .1010678*gx)*gx)*gx)*gx)*gx
      balfac = (4.0/pi)**(0.5*alpha) * gammaf * 2.e10
      gw2 = 0.5*(1 - alpha)
      gw1 = balfac * sigma * 10**(-12*gw2)
    elseif (.not. is_mol) then    ! Unsold is not used for molecules
      unsold = .true.
      if (gv < -20.0) then
        c625 = 10**(0.4 * gv)       ! recipe c6 SIU
        gw1  = 17. * c625
      elseif (gv < 0.0) then
        gw1 = 10**gv / (vbar2e4)**0.3  ! recipe 'W' Turbospectrum
      elseif (gv /= 0.0) then
        upper_ev = lower_ev + H_E * nu0
        if (lower_ev < cnt_ev) then    ! Unsold approximation
          c1 = 1./(cnt_ev - upper_ev)**2 - 1./(cnt_ev - lower_ev)**2
          if (c1 > 0) then
            c625 = (lower_ion**2 * c1)**0.4
            gw1  = gv * uns * c625
          end if
        end if
      end if
    end if
    gw1 = gw1 * over_4pi / nu0
    if (disable_vdw) then
      unsold = .true.
      gw1 = 0
    end if
  END SUBROUTINE prepare_vanderwaals

  !=====================================================================================
  ! linelist_out.txt writers -- one block per (species, isotope), in the same format
  ! the linelist reader (init_linelist) and tumult's LineList.from_tumult_format
  ! parse:  "'Fe I'" / "'iso' nlines" / records "lam_air evi loggf gv ga gq i j
  ! strength reach". strength is a diagnostic both readers ignore. reach [Angstrom] is
  ! the half-width of the line's evaluation window, where its opacity falls to the
  ! tolerance times the continuum, maximised over the states the list was written from
  ! (ll%reach times the vacuum wavelength). init_linelist reads it back so a run can
  ! select lines by it (add_opacity); hydrogen lines record none and write 0.
  !-------------------------------------------------------------------------------------
  SUBROUTINE write_ll_header(unit, id, iso, n)
    integer, intent(in) :: unit, iso, n
    character(len=*), intent(in) :: id
    !-----------------------------------------------------------------------------------
    write(unit, "(A)") "'"//trim(id)//"'"
    write(unit, "(A,I0,A,I0)") "'", iso, "' ", n
  END SUBROUTINE write_ll_header

  !-------------------------------------------------------------------------------------
  SUBROUTINE write_ll_record(unit, lam_air, evi, loggf, gv, ga, gq, i, j, strength, reach)
    integer, intent(in) :: unit, i, j
    real(8), intent(in) :: lam_air, evi, loggf, gv, ga, gq, strength, reach
    !-----------------------------------------------------------------------------------
    write(unit, "(F12.4, F10.4, F9.3, F14.6, 2(ES12.3), 2(I6), ES12.3, ES13.5)") &
      lam_air, evi, loggf, gv, ga, gq, i, j, strength, reach
  END SUBROUTINE write_ll_record

  !=====================================================================================
  ! Accumulate the used/strength/reach diagnostics of src into dst, two copies of the
  ! same linelist (a composition and one of its per-patch copies), so dst ends up with
  ! the maximum over every copy merged into it.
  !-------------------------------------------------------------------------------------
  SUBROUTINE merge_ll_diag(dst, src)
    type(linelist_type), intent(inout) :: dst
    type(linelist_type), intent(in) :: src
    integer :: n
    !-----------------------------------------------------------------------------------
    n = dst%nline
    if (n == 0) return
    if (src%nline /= n) error stop 'merge_ll_diag: the linelists are not copies'
    dst%used(1:n)     = dst%used(1:n) .or. src%used(1:n)
    dst%strength(1:n) = max(dst%strength(1:n), src%strength(1:n))
    dst%reach(1:n)    = max(dst%reach(1:n), src%reach(1:n))
  END SUBROUTINE merge_ll_diag

  !=====================================================================================
  ! All linelist_out.txt blocks of one composition particle: for each isotope (main
  ! first), a header plus the records that were used. ll is the per-isotope SoA array
  ! (0 = main, k = isotope iso_values(k)).
  !-------------------------------------------------------------------------------------
  SUBROUTINE print_ll_particle(unit, id, iso_values, ll)
    integer,             intent(in) :: unit
    character(len=*),    intent(in) :: id
    integer,             intent(in) :: iso_values(:)
    type(linelist_type), intent(in) :: ll(0:)
    integer :: l, kiso, isoval, n
    !-----------------------------------------------------------------------------------
    do kiso=0, ubound(ll,1)
      isoval = 0
      if (kiso > 0) isoval = iso_values(kiso)
      ! first pass: block size (the header carries the record count)
      n = 0
      do l=1, ll(kiso)%nline
        if (ll(kiso)%active(l) .and. ll(kiso)%used(l)) n = n + 1
      end do
      if (n == 0) cycle

      call write_ll_header(unit, id, isoval, n)
      do l=1, ll(kiso)%nline
        associate(sll => ll(kiso))
          if (.not. (sll%active(l) .and. sll%used(l))) cycle
          call write_ll_record(unit, vac2air(exp(sll%llam0(l))), &
            dble(sll%lower_ev(l)), dble(sll%loggf(l)), dble(sll%gv_in(l)), &
            dble(sll%ga_in(l)), dble(sll%gq_in(l)), sll%lev_i(l), sll%lev_j(l), &
            dble(sll%strength(l)), dble(sll%reach(l)) * exp(sll%llam0(l)))
        end associate
      end do
    end do
  END SUBROUTINE print_ll_particle

  !=====================================================================================
  ! All linelist_out.txt blocks of one model-atom SoA linelist (atom_trans_t%ll), one
  ! per (ion stage, isotope). ids(s-s0+1) is the species id of stage s; the isotope
  ! numbers are taken from iso_out directly (ascending, so the main block leads). evi
  ! is relative to the stage ground level and i/j are the atom level indices, so a
  ! block is a self-consistent input linelist for that species.
  !-------------------------------------------------------------------------------------
  SUBROUTINE print_ll_atom(unit, ll, ids, s0)
    integer,             intent(in) :: unit, s0
    character(len=*),    intent(in) :: ids(:)
    type(linelist_type), intent(in) :: ll
    integer :: l, s, n, k, nv, iv, isoval
    integer :: ivals(64)
    !-----------------------------------------------------------------------------------
    if (ll%nline == 0) return
    ! collect and sort the isotope numbers present among the exported lines
    nv = 0
    do l=1, ll%nline
      if (.not. (ll%active(l) .and. ll%used(l))) cycle
      if (any(ivals(1:nv) == ll%iso_out(l))) cycle
      nv = min(nv + 1, size(ivals))
      ivals(nv) = ll%iso_out(l)
    end do
    do k=2, nv                     ! insertion sort; nv is tiny
      iv = ivals(k)
      l  = k - 1
      do while (l >= 1)
        if (ivals(l) <= iv) exit
        ivals(l+1) = ivals(l)
        l = l - 1
      end do
      ivals(l+1) = iv
    end do

    do s = minval(ll%ipart(1:ll%nline)), maxval(ll%ipart(1:ll%nline))
      do k=1, nv
        isoval = ivals(k)
        n = 0
        do l=1, ll%nline
          if (ll%active(l) .and. ll%used(l) .and. ll%ipart(l) == s &
              .and. ll%iso_out(l) == isoval) n = n + 1
        end do
        if (n == 0) cycle
        call write_ll_header(unit, ids(s-s0+1), isoval, n)
        do l=1, ll%nline
          if (.not. (ll%active(l) .and. ll%used(l) .and. ll%ipart(l) == s &
              .and. ll%iso_out(l) == isoval)) cycle
          call write_ll_record(unit, vac2air(exp(ll%llam0(l))), &
            dble(ll%lower_ev(l)), dble(ll%loggf(l)), dble(ll%gv_in(l)), &
            dble(ll%ga_in(l)), dble(ll%gq_in(l)), ll%lev_i(l), ll%lev_j(l), &
            dble(ll%strength(l)), dble(ll%reach(l)) * exp(ll%llam0(l)))
        end do
      end do
    end do
  END SUBROUTINE print_ll_atom

  !=====================================================================================
  ! --- small private helpers ---------------------------------------------------------
  !=====================================================================================

  ! first real number on a record
  PURE FUNCTION parse_lam(str) result(lam)
    character(len=*), intent(in) :: str
    real(dp) :: lam
    integer  :: iostat
    read(str, *, iostat=iostat) lam
    if (iostat /= 0) lam = huge(1.0_dp)
  END FUNCTION parse_lam

  ! Decode ordinary numeric records without entering the formatted-I/O runtime once per
  ! transition. Unusual syntax falls back to the list-directed reader, preserving the
  ! input format's existing compatibility. Three columns are required; absent optional
  ! columns retain their historical zero defaults. values(9:10) are the trailing
  ! strength and reach a written line list carries, and nvalue counts the numbers read,
  ! so nvalue >= 10 marks a record with a reach. A non-numeric token after the eight
  ! input columns ends the ten-column parse, so the eight are then parsed on their own.
  SUBROUTINE parse_record(str, values, iostat, nvalue)
    character(len=*), intent(in) :: str
    real(dp), intent(out) :: values(10)
    integer, intent(out) :: iostat, nvalue
    logical :: ok
    values = 0
    iostat = 0
    call parse_fields_fast(str, values, nvalue, ok)
    if (ok .and. nvalue >= 3) return
    values = 0
    call parse_fields_fast(str, values(1:8), nvalue, ok)
    if (ok .and. nvalue >= 3) return
    values = 0
    nvalue = 0
    read(str, *, iostat=iostat) values(1:8)
  END SUBROUTINE parse_record

  PURE SUBROUTINE parse_fields_fast(str, values, nvalue, ok)
    character(len=*), intent(in) :: str
    real(dp), intent(out) :: values(:)
    integer, intent(out) :: nvalue
    logical, intent(out) :: ok
    integer :: p, last, digit, exponent, exponent_sign, fractional_digits, ndigit
    integer(8) :: significand
    real(dp) :: sign
    logical :: have_digit
    character :: c

    values = 0
    nvalue = 0
    ok = .true.
    p = 1
    last = len_trim(str)
    do while (p <= last .and. nvalue < size(values))
      do while (p <= last)
        c = str(p:p)
        if (c /= ' ' .and. c /= achar(9)) exit
        p = p + 1
      end do
      if (p > last) exit

      sign = 1.0_dp
      if (str(p:p) == '-' .or. str(p:p) == '+') then
        if (str(p:p) == '-') sign = -1.0_dp
        p = p + 1
      end if
      significand = 0_8
      ndigit = 0
      fractional_digits = 0
      have_digit = .false.
      do while (p <= last)
        digit = iachar(str(p:p)) - iachar('0')
        if (digit < 0 .or. digit > 9) exit
        if (ndigit >= 18) then
          ok = .false.
          return
        end if
        significand = 10_8*significand + digit
        ndigit = ndigit + 1
        have_digit = .true.
        p = p + 1
      end do
      if (p <= last) then
        if (str(p:p) == '.') then
          p = p + 1
          do while (p <= last)
            digit = iachar(str(p:p)) - iachar('0')
            if (digit < 0 .or. digit > 9) exit
            if (ndigit >= 18) then
              ok = .false.
              return
            end if
            significand = 10_8*significand + digit
            ndigit = ndigit + 1
            fractional_digits = fractional_digits + 1
            have_digit = .true.
            p = p + 1
          end do
        end if
      end if
      if (.not. have_digit) then
        ok = .false.
        return
      end if

      exponent = 0
      if (p <= last) then
        c = str(p:p)
        if (c == 'e' .or. c == 'E' .or. c == 'd' .or. c == 'D') then
          p = p + 1
          exponent_sign = 1
          if (p <= last) then
            if (str(p:p) == '-' .or. str(p:p) == '+') then
              if (str(p:p) == '-') exponent_sign = -1
              p = p + 1
            end if
          end if
          have_digit = .false.
          do while (p <= last)
            digit = iachar(str(p:p)) - iachar('0')
            if (digit < 0 .or. digit > 9) exit
            exponent = 10*exponent + digit
            have_digit = .true.
            p = p + 1
          end do
          if (.not. have_digit) then
            ok = .false.
            return
          end if
          exponent = exponent_sign*exponent
        end if
      end if
      if (p <= last) then
        c = str(p:p)
        if (c /= ' ' .and. c /= achar(9)) then
          ok = .false.
          return
        end if
      end if

      nvalue = nvalue + 1
      values(nvalue) = sign*real(significand, dp) * &
                       10.0_dp**(exponent - fractional_digits)
    end do
  END SUBROUTINE parse_fields_fast

  ! first record index whose leading wavelength is >= target (records sorted ascending)
  FUNCTION bisect_block(buf, lstart, lend, n, target) result(idx)
    character(len=*), intent(in) :: buf
    integer(8),       intent(in) :: lstart(:), lend(:)
    integer,          intent(in) :: n
    real(dp),         intent(in) :: target
    integer :: idx, lo, hi, mid
    !-----------------------------------------------------------------------------------
    lo = 1
    hi = n + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (parse_lam(buf(lstart(mid):lend(mid))) < target) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    idx = lo
  END FUNCTION bisect_block

  ! ensure every SoA array can hold (nline + extra) entries
  SUBROUTINE grow(ll, extra)
    type(linelist_type), intent(inout) :: ll
    integer, intent(in) :: extra
    integer :: need
    need = ll%nline + extra
    call grow8(ll%llam0,    ll%nline, need)
    call grow4(ll%nu0,      ll%nline, need)
    call grow4(ll%lower_ev, ll%nline, need)
    call grow4(ll%gpref,    ll%nline, need)
    call grow4(ll%loggf,    ll%nline, need)
    call grow4(ll%gv_in,    ll%nline, need)
    call grow4(ll%ga_in,    ll%nline, need)
    call grow4(ll%gq_in,    ll%nline, need)
    call grow4(ll%ga,       ll%nline, need)
    call grow4(ll%gw1,      ll%nline, need)
    call grow4(ll%gw2,      ll%nline, need)
    call grow4(ll%gq,       ll%nline, need)
    call grow4(ll%strength, ll%nline, need)
    call grow4(ll%reach,    ll%nline, need)
    call grow4(ll%reach_in, ll%nline, need)
    call growi(ll%lev_i,    ll%nline, need)
    call growi(ll%lev_j,    ll%nline, need)
    call growl(ll%unsold,   ll%nline, need)
    call growl(ll%active,   ll%nline, need)
    call growl(ll%used,     ll%nline, need)
  END SUBROUTINE grow

  ! pre-size all arrays (incl. the model-atom extensions) for nline+extra entries,
  ! so per-line append_line calls do not reallocate
  SUBROUTINE reserve_linelist(ll, extra)
    type(linelist_type), intent(inout) :: ll
    integer, intent(in) :: extra
    call grow(ll, extra)
    call grow_ext(ll, extra)
  END SUBROUTINE reserve_linelist

  ! model-atom extension arrays (see append_line)
  SUBROUTINE grow_ext(ll, extra)
    type(linelist_type), intent(inout) :: ll
    integer, intent(in) :: extra
    integer :: need
    need = ll%nline + extra
    call growi(ll%ipart,   ll%nline, need)
    call growi(ll%itrans,  ll%nline, need)
    call growi(ll%iso_out, ll%nline, need)
    call growi(ll%iws,     ll%nline, need)
    call growi(ll%iwe,     ll%nline, need)
    call growi(ll%icen,    ll%nline, need)
    call grow4(ll%wA,      ll%nline, need)
    call grow4(ll%cwt,     ll%nline, need)
    call grow4(ll%dej,     ll%nline, need)
  END SUBROUTINE grow_ext

  SUBROUTINE grow8(a, used, need)
    real(dp), allocatable, intent(inout) :: a(:)
    integer, intent(in) :: used, need
    real(dp), allocatable :: tmp(:)
    if (.not. allocated(a)) then
      allocate(a(need)); return
    end if
    if (size(a) >= need) return
    allocate(tmp(need))
    tmp(1:used) = a(1:used)
    call move_alloc(tmp, a)
  END SUBROUTINE grow8

  SUBROUTINE grow4(a, used, need)
    real(4), allocatable, intent(inout) :: a(:)
    integer, intent(in) :: used, need
    real(4), allocatable :: tmp(:)
    if (.not. allocated(a)) then
      allocate(a(need)); return
    end if
    if (size(a) >= need) return
    allocate(tmp(need))
    tmp(1:used) = a(1:used)
    call move_alloc(tmp, a)
  END SUBROUTINE grow4

  SUBROUTINE growi(a, used, need)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: used, need
    integer, allocatable :: tmp(:)
    if (.not. allocated(a)) then
      allocate(a(need)); return
    end if
    if (size(a) >= need) return
    allocate(tmp(need))
    tmp(1:used) = a(1:used)
    call move_alloc(tmp, a)
  END SUBROUTINE growi

  SUBROUTINE growl(a, used, need)
    logical, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: used, need
    logical, allocatable :: tmp(:)
    if (.not. allocated(a)) then
      allocate(a(need)); return
    end if
    if (size(a) >= need) return
    allocate(tmp(need))
    tmp(1:used) = a(1:used)
    call move_alloc(tmp, a)
  END SUBROUTINE growl

  ! reorder a(1:n) <- a(idx(1:n))
  SUBROUTINE permute8(a, idx)
    real(dp), allocatable, intent(inout) :: a(:)
    integer, intent(in) :: idx(:)
    real(dp), allocatable :: tmp(:)
    integer :: n, i
    n = size(idx)
    allocate(tmp(n))
    do i = 1, n
      tmp(i) = a(idx(i))
    end do
    a(1:n) = tmp
  END SUBROUTINE permute8

  SUBROUTINE permute4(a, idx)
    real(4), allocatable, intent(inout) :: a(:)
    integer, intent(in) :: idx(:)
    real(4), allocatable :: tmp(:)
    integer :: n, i
    n = size(idx)
    allocate(tmp(n))
    do i = 1, n
      tmp(i) = a(idx(i))
    end do
    a(1:n) = tmp
  END SUBROUTINE permute4

  SUBROUTINE permutei(a, idx)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: idx(:)
    integer, allocatable :: tmp(:)
    integer :: n, i
    n = size(idx)
    allocate(tmp(n))
    do i = 1, n
      tmp(i) = a(idx(i))
    end do
    a(1:n) = tmp
  END SUBROUTINE permutei

  SUBROUTINE permutel(a, idx)
    logical, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: idx(:)
    logical, allocatable :: tmp(:)
    integer :: n, i
    n = size(idx)
    allocate(tmp(n))
    do i = 1, n
      tmp(i) = a(idx(i))
    end do
    a(1:n) = tmp
  END SUBROUTINE permutel

  ! argsort: reorder idx so that key(idx) is ascending (recursive quicksort on indices)
  RECURSIVE SUBROUTINE argsort8(key, idx)
    real(dp), intent(in)    :: key(:)
    integer,  intent(inout) :: idx(:)
    call qs(1, size(idx))
  contains
    recursive subroutine qs(lo, hi)
      integer, intent(in) :: lo, hi
      integer  :: i, j, p
      real(dp) :: pivot
      if (lo >= hi) return
      pivot = key(idx((lo+hi)/2))
      i = lo; j = hi
      do
        do while (key(idx(i)) < pivot); i = i + 1; end do
        do while (key(idx(j)) > pivot); j = j - 1; end do
        if (i <= j) then
          p = idx(i); idx(i) = idx(j); idx(j) = p
          i = i + 1; j = j - 1
        end if
        if (i > j) exit
      end do
      if (lo < j) call qs(lo, j)
      if (i < hi) call qs(i, hi)
    end subroutine qs
  END SUBROUTINE argsort8

  !=====================================================================================
  ! Check if a wavelength (in vacuum) falls within any of the given intervals.
  ! An interval includes its boundaries (closed interval [blue, red]).
  !-------------------------------------------------------------------------------------
  ! First index with a(i) >= x, or n+1 when every entry is below x. Together with
  ! upper_bound8 this brackets a sorted array; an empty bracket gives l0 > l1 and so an
  ! empty line loop, which is the correct result for a patch no line reaches.
  PURE INTEGER FUNCTION lower_bound8(a, n, x) result(idx)
    real(8), intent(in) :: a(:)
    integer, intent(in) :: n
    real(8), intent(in) :: x
    integer :: lo, hi, mid
    lo = 1
    hi = n + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (a(mid) >= x) then
        hi = mid
      else
        lo = mid + 1
      end if
    end do
    idx = lo
  END FUNCTION lower_bound8

  ! Last index with a(i) <= x, or 0 when every entry is above x.
  PURE INTEGER FUNCTION upper_bound8(a, n, x) result(idx)
    real(8), intent(in) :: a(:)
    integer, intent(in) :: n
    real(8), intent(in) :: x
    integer :: lo, hi, mid
    lo = 1
    hi = n + 1
    do while (lo < hi)
      mid = (lo + hi) / 2
      if (a(mid) > x) then
        hi = mid
      else
        lo = mid + 1
      end if
    end do
    idx = lo - 1
  END FUNCTION upper_bound8

  PURE FUNCTION in_any_interval(lam_vac, intervals) result(in_interval)
    real(dp), intent(in) :: lam_vac
    type(interval_t), intent(in) :: intervals(:)
    logical :: in_interval
    integer :: i
    !-----------------------------------------------------------------------------------
    in_interval = .false.
    do i = 1, size(intervals)
      if (lam_vac >= intervals(i)%blue_vac .and. lam_vac <= intervals(i)%red_vac) then
        in_interval = .true.
        return
      end if
    end do
  END FUNCTION in_any_interval

  !=====================================================================================
  ! Build a union of intervals from output segments and return them sorted and merged.
  ! Input:
  !   windows(:,1) = blue edges (in vacuum Angstrom)
  !   windows(:,2) = red edges (in vacuum Angstrom)
  !   dlam = expansion margin (in Angstrom)
  ! Output:
  !   merged_intervals = sorted, non-overlapping intervals
  !   n_merged = number of merged intervals
  !
  ! The routine:
  ! 1. Expands each window by dlam on both sides
  ! 2. Clamps the lower edge to a physically valid wavelength
  ! 3. Sorts the expanded intervals
  ! 4. Merges overlapping or touching intervals
  !-------------------------------------------------------------------------------------
  SUBROUTINE merge_intervals(windows, n_windows, dlam, merged_intervals, n_merged)
    real(dp), intent(in) :: windows(:,:)
    integer, intent(in) :: n_windows
    real(dp), intent(in) :: dlam
    type(interval_t), intent(out), allocatable :: merged_intervals(:)
    integer, intent(out) :: n_merged
    type(interval_t), allocatable :: expanded(:)
    integer :: i, j
    real(dp) :: min_lam = 1.0_dp  ! physically valid minimum wavelength
    !-----------------------------------------------------------------------------------

    if (n_windows <= 0) then
      allocate(merged_intervals(0))
      n_merged = 0
      return
    end if

    ! Step 1: Expand windows by dlam
    allocate(expanded(n_windows))
    do i = 1, n_windows
      expanded(i)%blue_vac = max(windows(i,1) - dlam, min_lam)
      expanded(i)%red_vac  = windows(i,2) + dlam
    end do

    ! Step 2: Sort by blue edge
    call sort_intervals(expanded, n_windows)

    ! Step 3: Merge overlapping/touching intervals
    allocate(merged_intervals(n_windows))  ! worst case: no merging
    n_merged = 1
    merged_intervals(1) = expanded(1)

    do i = 2, n_windows
      ! if current interval starts before or at the end of the last merged one, merge them
      if (expanded(i)%blue_vac <= merged_intervals(n_merged)%red_vac) then
        ! merge: extend the red edge if needed
        merged_intervals(n_merged)%red_vac = max(merged_intervals(n_merged)%red_vac, &
                                                  expanded(i)%red_vac)
      else
        ! gap: start a new merged interval
        n_merged = n_merged + 1
        merged_intervals(n_merged) = expanded(i)
      end if
    end do

    ! trim to actual size
    if (n_merged < n_windows) then
      merged_intervals = merged_intervals(1:n_merged)
    end if

  END SUBROUTINE merge_intervals

  !=====================================================================================
  ! Sort an array of intervals by blue edge (ascending).
  !-------------------------------------------------------------------------------------
  SUBROUTINE sort_intervals(intervals, n)
    type(interval_t), intent(inout) :: intervals(:)
    integer, intent(in) :: n
    integer :: i, j
    type(interval_t) :: temp
    !-----------------------------------------------------------------------------------
    ! simple bubble sort for small arrays
    do i = 1, n-1
      do j = i+1, n
        if (intervals(j)%blue_vac < intervals(i)%blue_vac) then
          temp = intervals(i)
          intervals(i) = intervals(j)
          intervals(j) = temp
        end if
      end do
    end do
  END SUBROUTINE sort_intervals

END MODULE linelist_mod
