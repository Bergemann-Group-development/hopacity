!=======================================================================================
! linelist_mod -- structure-of-arrays (SoA) bulk linelist handling.
!
! This is a fresh, performance-oriented reimplementation of the bulk metal/molecular
! linelist opacity that lived, in array-of-structures (AoS) form, as particle%line(:)
! (bbtrans_type); that AoS bb form now survives only for hydrogen (Stark profiles).
! The plain linelist_type data
! container lives in the leaf module line_types_mod, so that particle_type (ion_mod) can
! embed it without depending on this heavy module. That inversion lets linelist_mod, in
! turn, USE ion_mod and drive the opacity kernels directly from a particle_type (see
! add_particle_opacity), while the low-level kernels still take plain patch work arrays.
!
! Key differences from the AoS path:
!   * SoA storage: one linelist_type per particle, parallel arrays over lines.
!   * The opacity kernel (add_opacity) loops the SPATIAL domain in the outer loop, then
!     lines (sorted ascending in wavelength), then frequencies. A monotonic frequency
!     cursor ("march") replaces the per-line get_center search.
!   * Works in log(lambda): the thermal Doppler width is constant per cell and the
!     Doppler shift is additive, so the profile abscissa is vs = idf*(llam(f) - llc).
!   * loggf is folded into the Boltzmann exponential, exp(loggf*ln10 - ex), instead of
!     pre-converting to Einstein Bij.
!
! The orchestration that walks elements/ions/molecules stays in composition_mod
! (add_linelist_opacities2 / init_linelist2), which owns the patch thermodynamic arrays
! and the departure-coefficient lookup; per particle it now calls add_particle_opacity.
!=======================================================================================
MODULE linelist_mod
  USE prec_mod
  USE const_mod
  USE hop_env_mod, only: trace
  USE handyfuncs_mod, only: air2vac, vac2air
  USE op_math_mod, only: find_closest
  USE bbtrans_mod, only: voigt_function2, disable_ga, disable_stark, disable_vdw, &
                         f_conv, tol => tolerance, bbtrans_type, vdw_shift
  USE line_types_mod, only: linelist_type
  USE ion_mod, only: particle_type
  USE voigt_av_mod, only: voigt_av_t, voigt_av_setup, voigt_av_block
  implicit none
  private

  public :: build_linelist, finalize_linelist, add_opacity, add_opacity_equi, voigt_peak
  public :: add_particle_opacity
  public :: append_line, write_ll_header, write_ll_record, bij_const
  public :: reserve_linelist, print_ll_particle, print_ll_atom

  ! Bij = bij_const * (10**loggf) / nu0  (see init_einstein_coefficients in bbtrans_mod)
  real(dp), parameter :: bij_const = c2_2h * f_conv * h_4pi
  real(dp), parameter :: over_h4pi = 1.0_dp / h_4pi
  real(4),  parameter :: over_pi   = 0.3183098862
  real(dp), parameter :: ex_cutoff = -85.0_dp  ! skip lines whose log(Rij) underflows
  real(4), parameter  :: voigt_peak = 0.5641895835 ! 1/sqrt(pi): max of voigt_function2

  ! --- Voigt kernel selection (add_opacity_equi only) ---------------------------------
  ! 0 = voigt_function2 (bbtrans_mod, Traving/Eriksson; the historical default)
  ! 1 = voigt_av        (voigt_av_mod, branchless fp32; vectorises, and is the more
  !                      accurate of the two for damped lines -- see voigt_av_mod header)
  ! Selecting 1 shifts spectra by ~2e-4 in the weak-line regime and by up to a few
  ! percent for a > 0.2 lines, where voigt_function2 is the one that was wrong.
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
  SUBROUTINE build_linelist(ll, buf, lstart, lend, n, awgt, ion_j, xion, is_mol, &
                            vbar2e4, do_mask, lam_min, lam_max)
    type(linelist_type), intent(inout) :: ll
    character(len=*),    intent(in)    :: buf
    integer(8),          intent(in)    :: lstart(:), lend(:)
    integer,             intent(in)    :: n, ion_j
    real(dp),            intent(in)    :: awgt, xion, vbar2e4
    logical,             intent(in)    :: is_mol, do_mask
    real(dp),            intent(in)    :: lam_min, lam_max
    real(dp) :: lam_air, lower_ev, loggf, gv_in, ga_in, gq_in, r_i, r_j
    real(dp) :: nu0, lam0, ga_t, gw1_t, gw2_t, gq_t
    integer  :: lo, hi, nkeep, i, k, iostat
    integer, save :: it_mask=0, it_parse=0
    !-----------------------------------------------------------------------------------
    if (n <= 0) return

    ! locate the contiguous in-range slice [lo:hi]
    call trace%begin ('ll_build_mask', itimer=it_mask)
    if (do_mask) then
      lo = bisect_block(buf, lstart, lend, n, lam_min)     ! first record with lam >= lam_min
      hi = bisect_block(buf, lstart, lend, n, lam_max) - 1 ! last  record with lam <= lam_max
      ! bisect_block returns the first index with lam >= target; refine the upper edge
      if (hi < n) then
        if (parse_lam(buf(lstart(hi+1):lend(hi+1))) <= lam_max) hi = hi + 1
      end if
    else
      lo = 1
      hi = n
    end if
    call trace%end(it_mask)
    nkeep = hi - lo + 1
    if (nkeep <= 0) return

    call grow(ll, nkeep)

    call trace%begin ('ll_build_parse', itimer=it_parse)
    do i = lo, hi
      gv_in = 0; ga_in = 0; gq_in = 0; r_i = 0; r_j = 0
      read(buf(lstart(i):lend(i)), *, iostat=iostat) lam_air, lower_ev, loggf, gv_in, &
                                            ga_in, gq_in, r_i, r_j
      ! only lam_air, lower_ev, loggf are required; the rest default to 0

      k = ll%nline + 1
      ll%nline = k
      lam0  = air2vac(lam_air)
      nu0   = HZ_TO_AA / lam0
      ll%llam0(k)    = log(lam0)
      ll%nu0(k)      = nu0
      ll%lower_ev(k) = lower_ev
      ll%gpref(k)    = loggf*ln10 + log(bij_const / nu0)
      ll%lev_i(k)    = nint(r_i)
      ll%lev_j(k)    = nint(r_j)
      ll%active(k)   = .true.
      ll%used(k)     = .false.
      ll%strength(k) = 0
      ll%loggf(k)    = loggf
      ll%gv_in(k)    = gv_in
      ll%ga_in(k)    = ga_in
      ll%gq_in(k)    = gq_in

      call prepare_radiative   (ga_t,  ga_in, nu0, loggf)
      call prepare_stark       (gq_t,  gq_in, nu0, lower_ev, ion_j, xion)
      call prepare_vanderwaals (gw1_t, gw2_t, ll%unsold(k), gv_in, nu0, &
                                lower_ev, ion_j, xion, vbar2e4, is_mol, ll%vdw_ratio(k))
      ll%ga(k)  = ga_t
      ll%gq(k)  = gq_t
      ll%gw1(k) = gw1_t
      ll%gw2(k) = gw2_t
    end do
    call trace%end(it_parse)
  END SUBROUTINE build_linelist

  !=====================================================================================
  ! Append a single line with explicit (already prepared) parameters; used to convert
  ! model-atom bound-bound transitions (atom_mod%build_ll) into a linelist_type.
  ! Broadening parameters carry the same over_4pi/nu0 scaling as the prepare_* routines.
  !-------------------------------------------------------------------------------------
  SUBROUTINE append_line(ll, llam0, nu0, lower_ev, gpref, ga, gw1, gw2, gq, unsold, &
                         vdw_ratio, lev_i, lev_j, itrans, ipart, wA, cwt, dej, iso_out, &
                         loggf, gv_in, ga_in, gq_in)
    type(linelist_type), intent(inout) :: ll
    real(8), intent(in) :: llam0, gpref
    real(8), intent(in) :: nu0, lower_ev, ga, gw1, gw2, gq, wA, cwt, dej
    real(8), intent(in) :: loggf, gv_in, ga_in, gq_in
    logical, intent(in) :: unsold
    real(4), intent(in) :: vdw_ratio
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
    ll%vdw_ratio(k) = vdw_ratio
    ll%lev_i(k)    = lev_i
    ll%lev_j(k)    = lev_j
    ll%active(k)   = .true.
    ll%used(k)     = .false.
    ll%strength(k) = 0
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
    call permutei(ll%lev_i,    idx)
    call permutei(ll%lev_j,    idx)
    call permutel(ll%unsold,   idx)
    call permutel(ll%active,   idx)
    call permutel(ll%used,     idx)
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
  !-------------------------------------------------------------------------------------
  SUBROUTINE add_opacity(ll, llam, chi, eta, cont, logref, vproj, inv_df, vbar2, e_kT, &
                         lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
                         one_m_eps, chi_scat, do_scat)
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
    integer  :: i, j, k, l, f, f0, fc, nf, m(3), nl, np, ws, we, fa, fb
    real(8)  :: llc, llc_s, dl
    real(4)  :: idf, unsold_c, ekt, ne_c, tothi_c, starkpe, vproj_c
    real(4)  :: a, avdw, ex, Rij, bi, bj, ref, vs, phi, lte, lvb2, vm2, edge
    real(8)  :: arg
    real(4)  :: base_s(size(lnbase,4))
    logical  :: nlte, has_ipart
    integer(8) :: c_visit, c_pass, c_voigt, c_march
    !-----------------------------------------------------------------------------------
    m  = shape(vproj)
    nf = size(llam)
    nl = ll%nline
    np = size(lnbase, 4)
    if (nl == 0) return
    has_ipart = allocated(ll%ipart)
    c_visit = 0; c_pass = 0; c_voigt = 0; c_march = 0
    ! the spectrum machinery always delivers lam (hence llam) ascending; the monotonic
    ! march relies on this, so assert it once rather than branching per line.
    if (nf > 1 .and. llam(nf) < llam(1)) &
      error stop 'linelist add_opacity: log(lambda) grid must be ascending'

    do k = 1, m(3)
      do j = 1, m(2)
        do i = 1, m(1)
          idf      = inv_df(i,j,k)
          unsold_c = vbar2(i,j,k)**0.3 * tothi(i,j,k)
          lvb2     = log(vbar2(i,j,k))                    ! vb2**gw2 as exp(gw2*log(vb2))
          ekt      = e_kT(i,j,k)
          ne_c     = ne(i,j,k)
          tothi_c  = tothi(i,j,k)
          starkpe  = stark_pe(i,j,k)
          vproj_c  = vproj(i,j,k)
          ! line-independent part of log(Rij) per population group (idf folded in)
          base_s(1:np) = lnbase(i,j,k,1:np)
          fc = 1
          ws = 1
          we = nf

          do l = 1, nl
            if (.not. ll%active(l)) cycle
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

            ! log-space tolerance pre-test: voigt_peak*Rij <= ref  <=>  arg <= logref(f0).
            ! This rejects the vast majority of lines without evaluating exp().
            if (arg <= logref(f0,i,j,k)) cycle

            ! this ll's isotope Doppler factor idf is folded into lnbase (hence arg)
            Rij = exp(arg)

            nlte = do_nlte .and. ll%lev_i(l) > 0 .and. ll%lev_j(l) > 0
            if (nlte) then
              bi = dep_mem(i,j,k, ll%lev_i(l))
              bj = dep_mem(i,j,k, ll%lev_j(l))
            else
              bi = 1.0
              bj = 1.0
            end if

            ! folding b_i into ref makes every phi*Rij test an actual-opacity test
            ! (phi*Rij >= ref <=> phi*Rij*bi >= tol*cont) at no per-frequency cost
            ref = tol * cont(f0,i,j,k) / bi

            ! exact pre-filter: voigt_function2 <= 1/sqrt(pi) everywhere, so if even the
            ! peak cannot reach the tolerance the whole line is negligible -- skip
            ! the damping and Voigt evaluations entirely (most lines drop out here).
            if (voigt_peak * Rij <= ref) cycle
            c_pass = c_pass + 1

            ! damping parameter (physical, then scaled by idf as in the AoS kernel).
            ! avdw = the van der Waals contribution (pre-idf) is captured separately so
            ! it can drive the Lindholm pressure shift below.
            a = ll%ga(l)
            if (ll%unsold(l)) then
              avdw = ll%gw1(l) * unsold_c
            else
              avdw = ll%gw1(l) * exp(ll%gw2(l) * lvb2) * tothi_c
            end if
            a = a + avdw
            if (ll%gq(l) >= 0.0) then
              a = a + ll%gq(l) * ne_c
            else
              a = a - ll%gq(l) * starkpe
            end if
            a = a * idf

            ! Lindholm van der Waals pressure shift (red): the pre-idf vdW damping term
            ! is the vdW HWHM in ln(lambda) units, so shifting the profile centre red by
            ! vdw_shift*ratio*HWHM is a simple offset in llc. == llc when vdw_shift==0.
            ! ratio is the per-line impact-theory d/w (ll%vdw_ratio, precomputed at build
            ! time). The window [fa,fb] and f0 stay on the unshifted centre (shift << window).
            llc_s = llc + vdw_shift * ll%vdw_ratio(l) * avdw

            ll%used(l)     = .true.
            ll%strength(l) = max(ll%strength(l), voigt_peak*Rij / cont(f0,i,j,k))

            ! analytic profile extent where phi*Rij drops to ref: Gaussian core
            ! v^2 = ln(peak/ref) = arg - logref, Lorentz wing v^2 = a*Rij/(pi*ref);
            ! a straight (data-independent) loop over [fa,fb], centre included,
            ! pipelines much better than a walk whose exit condition depends on
            ! the previous Voigt value
            vm2  = max(real(arg - logref(f0,i,j,k)), a*Rij*over_pi/ref)
            edge = idf * real(max(llam(we) - llc, llc - llam(ws)))
            if (vm2 >= edge*edge) then
              fa = ws         ! line fills the window; skip the extent bisections
              fb = we
            else
              dl = sqrt(vm2) / idf
              fa = bisect_nearest(llam, ws, f0, llc - dl)
              fb = bisect_nearest(llam, f0, we, llc + dl)
            end if
            c_voigt = c_voigt + fb - fa + 1

            if (do_scat .and. .not. nlte) then
              ! LTE line (bi=bj=1) with the scattering fraction accumulated
              do f = fa, fb
                vs   = idf * real(llam(f) - llc_s)
                phi  = voigt_function2(a, vs)
                lte  = phi * Rij
                chi(f,i,j,k) = chi(f,i,j,k) + lte
                eta(f,i,j,k) = eta(f,i,j,k) + lte
                chi_scat(f,i,j,k) = chi_scat(f,i,j,k) + one_m_eps(f,i,j,k)*lte
              end do
            else
              do f = fa, fb
                vs   = idf * real(llam(f) - llc_s)
                phi  = voigt_function2(a, vs)
                lte  = phi * Rij
                chi(f,i,j,k) = chi(f,i,j,k) + lte*bi
                eta(f,i,j,k) = eta(f,i,j,k) + lte*bj
              end do
            end if
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
  SUBROUTINE add_opacity_equi(ll, llam, chi, eta, cont, logref, vproj, inv_df, vbar2, e_kT, &
                              lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
                              one_m_eps, chi_scat, do_scat)
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
    integer  :: i, j, k, l, f, f0, nf, m(3), nl, np, fa, fb, noff, nw
    real(8)  :: llc, llc_s, llam1, dllam, inv_dllam
    real(4)  :: idf, unsold_c, ekt, ne_c, tothi_c, starkpe, vproj_c
    real(4)  :: a, avdw, ex, Rij, bi, bj, ref, vs, phi, lte, lvb2, vm2
    real(8)  :: arg
    real(4)  :: base_s(size(lnbase,4))
    logical  :: nlte, has_ipart
    integer(8) :: c_visit, c_pass, c_voigt
    ! voigt_kernel==1 scratch: per-line profile state, the window's K values, and the
    ! fp32 abscissa ramp. Kbuf is an automatic array (stack) but only [1:nw] is ever
    ! touched, and windows average ~8 points, so it stays L1-resident.
    type(voigt_av_t) :: vst
    real(4)  :: Kbuf(size(llam)), v0, dv, Rn
    !-----------------------------------------------------------------------------------
    m  = shape(vproj)
    nf = size(llam)
    nl = ll%nline
    np = size(lnbase, 4)
    if (nl == 0) return
    if (nf < 2) return
    has_ipart = allocated(ll%ipart)
    c_visit = 0; c_pass = 0; c_voigt = 0
    if (llam(nf) < llam(1)) &
      error stop 'linelist add_opacity_equi: log(lambda) grid must be ascending'

    ! constant grid spacing in log(lambda); index = nint((llc - llam1)/dllam) + 1
    llam1     = llam(1)
    dllam     = (llam(nf) - llam(1)) / real(nf - 1, 8)
    inv_dllam = 1.0_dp / dllam

    do k = 1, m(3)
      do j = 1, m(2)
        do i = 1, m(1)
          idf      = inv_df(i,j,k)
          unsold_c = vbar2(i,j,k)**0.3 * tothi(i,j,k)
          lvb2     = log(vbar2(i,j,k))                    ! vb2**gw2 as exp(gw2*log(vb2))
          ekt      = e_kT(i,j,k)
          ne_c     = ne(i,j,k)
          tothi_c  = tothi(i,j,k)
          starkpe  = stark_pe(i,j,k)
          vproj_c  = vproj(i,j,k)
          ! line-independent part of log(Rij) per population group (idf folded in)
          base_s(1:np) = lnbase(i,j,k,1:np)

          do l = 1, nl
            if (.not. ll%active(l)) cycle
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

            ! log-space tolerance pre-test: voigt_peak*Rij <= ref  <=>  arg <= logref(f0).
            ! This rejects the vast majority of lines without evaluating exp().
            if (arg <= logref(f0,i,j,k)) cycle

            ! this ll's isotope Doppler factor idf is folded into lnbase (hence arg)
            Rij = exp(arg)

            nlte = do_nlte .and. ll%lev_i(l) > 0 .and. ll%lev_j(l) > 0
            if (nlte) then
              bi = dep_mem(i,j,k, ll%lev_i(l))
              bj = dep_mem(i,j,k, ll%lev_j(l))
            else
              bi = 1.0
              bj = 1.0
            end if

            ! folding b_i into ref makes every phi*Rij test an actual-opacity test
            ! (phi*Rij >= ref <=> phi*Rij*bi >= tol*cont) at no per-frequency cost
            ref = tol * cont(f0,i,j,k) / bi

            ! exact pre-filter: voigt_function2 <= 1/sqrt(pi) everywhere, so if even the
            ! peak cannot reach the tolerance the whole line is negligible -- skip
            ! the damping and Voigt evaluations entirely (most lines drop out here).
            if (voigt_peak * Rij <= ref) cycle
            c_pass = c_pass + 1

            ! damping parameter (physical, then scaled by idf as in the AoS kernel).
            ! avdw = the van der Waals contribution (pre-idf) is captured separately so
            ! it can drive the Lindholm pressure shift below.
            a = ll%ga(l)
            if (ll%unsold(l)) then
              avdw = ll%gw1(l) * unsold_c
            else
              avdw = ll%gw1(l) * exp(ll%gw2(l) * lvb2) * tothi_c
            end if
            a = a + avdw
            if (ll%gq(l) >= 0.0) then
              a = a + ll%gq(l) * ne_c
            else
              a = a - ll%gq(l) * starkpe
            end if
            a = a * idf

            ! Lindholm van der Waals pressure shift (red): the pre-idf vdW damping term
            ! is the vdW HWHM in ln(lambda) units, so shifting the profile centre red by
            ! vdw_shift*ratio*HWHM is a simple offset in llc. == llc when vdw_shift==0.
            ! ratio is the per-line impact-theory d/w (ll%vdw_ratio, precomputed at build
            ! time). The window [fa,fb] and f0 stay on the unshifted centre (shift << window).
            llc_s = llc + vdw_shift * ll%vdw_ratio(l) * avdw

            ll%used(l)     = .true.
            ll%strength(l) = max(ll%strength(l), voigt_peak*Rij / cont(f0,i,j,k))

            ! analytic profile extent where phi*Rij drops to ref, in v-units: Gaussian
            ! core v^2 = arg - logref, Lorentz far wing v^2 = a*Rij/(pi*ref). On a
            ! uniform log(lambda) grid the half-width in llam is sqrt(vm2)/idf, so the
            ! index half-width is that over dllam -- no bisection required.
            vm2  = max(real(arg - logref(f0,i,j,k)), a*Rij*over_pi/ref)
            noff = int(sqrt(vm2) / idf * inv_dllam) + 1
            fa   = max(1,  f0 - noff)
            fb   = min(nf, f0 + noff)
            c_voigt = c_voigt + fb - fa + 1

            if (voigt_kernel == 1) then
              ! Branchless fp32 kernel. The whole window is evaluated first, so the
              ! profile loop and the accumulation loop each vectorise; K stays in L1.
              ! The abscissa is a pure fp32 ramp -- on a uniform log(lambda) grid
              ! llam(f) - llc = (llam(fa) - llc) + (f-fa)*dllam exactly, so the real(8)
              ! llam stream drops out of the inner loop altogether.
              nw = fb - fa + 1
              v0 = idf * real(llam(fa) - llc)
              dv = idf * real(dllam)
              call voigt_av_setup (a, vst)
              call voigt_av_block (vst, v0, dv, nw, Kbuf)
              ! voigt_av returns K = H(a,v); voigt_function2 returns H(a,v)/sqrt(pi)
              Rn = Rij * voigt_peak
              if (do_scat .and. .not. nlte) then
                do f = fa, fb
                  lte = Kbuf(f-fa+1) * Rn
                  chi(f,i,j,k) = chi(f,i,j,k) + lte
                  eta(f,i,j,k) = eta(f,i,j,k) + lte
                  chi_scat(f,i,j,k) = chi_scat(f,i,j,k) + one_m_eps(f,i,j,k)*lte
                end do
              else
                do f = fa, fb
                  lte = Kbuf(f-fa+1) * Rn
                  chi(f,i,j,k) = chi(f,i,j,k) + lte*bi
                  eta(f,i,j,k) = eta(f,i,j,k) + lte*bj
                end do
              end if

            else if (do_scat .and. .not. nlte) then
              ! LTE line (bi=bj=1) with the scattering fraction accumulated
              do f = fa, fb
                vs   = idf * real(llam(f) - llc_s)
                phi  = voigt_function2(a, vs)
                lte  = phi * Rij
                chi(f,i,j,k) = chi(f,i,j,k) + lte
                eta(f,i,j,k) = eta(f,i,j,k) + lte
                chi_scat(f,i,j,k) = chi_scat(f,i,j,k) + one_m_eps(f,i,j,k)*lte
              end do
            else
              do f = fa, fb
                vs   = idf * real(llam(f) - llc_s)
                phi  = voigt_function2(a, vs)
                lte  = phi * Rij
                chi(f,i,j,k) = chi(f,i,j,k) + lte*bi
                eta(f,i,j,k) = eta(f,i,j,k) + lte*bj
              end do
            end if
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
                                  do_nlte, equidistant, one_m_eps, chi_scat)
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
    real(4),  allocatable :: idf3(:,:,:), vb3(:,:,:), lnbase(:,:,:,:)
    real(4)  :: eps_dummy(0,0,0,0), scat_dummy(0,0,0,0)
    real(dp) :: iso_frac
    integer  :: m(3), kiso, ll_n, l
    logical  :: do_scat
    !-----------------------------------------------------------------------------------
    m = shape(temp)
    allocate(idf3(m(1),m(2),m(3)), vb3(m(1),m(2),m(3)), lnbase(m(1),m(2),m(3),1))
    do_scat = present(one_m_eps) .and. present(chi_scat)

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
        lnbase(:,:,:,1) = log(max(ptc%n, tiny(1.0_dp))) - ptc%part &
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
              one_m_eps, chi_scat, .true.)
          else
            call add_opacity_equi(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              eps_dummy, scat_dummy, .false.)
          end if
        else
          if (do_scat) then
            call add_opacity(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              one_m_eps, chi_scat, .true.)
          else
            call add_opacity(sll, llam, chi, eta, cont, logref, vproj, idf3, vb3, &
              e_kT, lnbase, ne, tothi, stark_pe, dep_mem, do_nlte, &
              eps_dummy, scat_dummy, .false.)
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
  ! Ported broadening setup (scalar versions of the bbtrans_mod prepare_* routines).
  ! Results are scaled by over_4pi/nu0, identical to the AoS path.
  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE prepare_radiative(ga, ga_in, nu0, loggf)
    real(dp), intent(out) :: ga
    real(dp), intent(in)  :: ga_in, nu0, loggf
    !-----------------------------------------------------------------------------------
    if (ga_in == 0) then
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
    real(dp), parameter :: scale = 1e-4 * one_sixth
    real(dp), parameter :: c19 = log10(39. * (8. / (PI*EM*BK**5))**one_sixth) ! ~19.4
    real(dp) :: upper_ev, c1
    !-----------------------------------------------------------------------------------
    gq = 0
    if (gq_in < -10) then
      gq = - BK * 10**(c19 + two_third*gq_in)   ! recipe c4 SIU (Gray p.237); sign marks recipe
    elseif (gq_in < 0) then
      gq = - 10**gq_in * scale                  ! sign marks recipe
    elseif (gq_in == 0) then
      upper_ev = lower_ev + H_E * nu0
      if (lower_ev < cnt_ev) then
        c1 = (j / (cnt_ev - upper_ev))**2
        if (j == 1) then
          c1 = 2.26e-7 * c1
        else
          c1 = 5.42e-7 / (j+1)**2
        end if
        gq = RydV * c1
      end if
    else
      gq = gq_in
    end if
    gq = gq * over_4pi / nu0
    if (disable_stark) gq = 0
  END SUBROUTINE prepare_stark

  !-------------------------------------------------------------------------------------
  PURE SUBROUTINE prepare_vanderwaals(gw1, gw2, unsold, gv, nu0, lower_ev, lower_ion, &
                                      cnt_ev, vbar2e4, is_mol, vdw_ratio)
    real(dp), intent(out) :: gw1, gw2
    logical,  intent(out) :: unsold
    real(dp), intent(in)  :: gv, nu0, lower_ev, cnt_ev, vbar2e4
    integer,  intent(in)  :: lower_ion
    logical,  intent(in)  :: is_mol
    real(4),  intent(out) :: vdw_ratio
    real(dp) :: upper_ev, sigma, alpha, gx, gammaf, balfac, c1, c625
    real(dp), parameter :: uns = 8.411 * 1.283984e-12
    real(dp), parameter :: barklem = (1e-2 * Rbohr)**2
    !-----------------------------------------------------------------------------------
    gw1 = 0
    gw2 = 0
    unsold = .false.
    vdw_ratio = 0
    if (gv >= 20.0) then
      ! Barklem cross section + velocity exponent
      sigma  = barklem * int(gv)
      alpha  = gv - int(gv)
      gx     = 1.0 - 0.5*alpha
      gammaf = 1 + (-.5748646 + (.9512363 + (-.6998588 + (.4245549 - .1010678*gx)*gx)*gx)*gx)*gx
      balfac = (4.0/pi)**(0.5*alpha) * gammaf * 2.e10
      gw2 = 0.5*(1 - alpha)
      gw1 = balfac * sigma * 10**(-12*gw2)
      ! Lindholm shift/HWHM ratio for an effective inverse-power interaction:
      ! d/w = tan(pi*alpha/2). The clamp guards against pathological alpha (near 1,
      ! where the tangent blows up).
      vdw_ratio = min(tan(0.5*pi * min(max(alpha, 0.02_dp), 0.96_dp)), 3.0_dp)
    elseif (.not. is_mol) then    ! Unsold is not used for molecules
      unsold = .true.
      vdw_ratio = 0.7265  ! pure van der Waals (p=6, alpha=0.4) -> tan(0.2*pi)
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
      gw1 = 0                       ! vdw_ratio is then irrelevant: it multiplies avdw = 0
    end if
  END SUBROUTINE prepare_vanderwaals

  !=====================================================================================
  ! linelist_out.txt writers -- one block per (species, isotope), in the same format
  ! the linelist reader (init_linelist2) and tumult's LineList.from_tumult_format
  ! parse:  "'Fe I'" / "'iso' nlines" / records "lam_air evi loggf gv ga gq i j
  ! strength" (the trailing strength column is a diagnostic; both readers ignore it).
  !-------------------------------------------------------------------------------------
  SUBROUTINE write_ll_header(unit, id, iso, n)
    integer, intent(in) :: unit, iso, n
    character(len=*), intent(in) :: id
    !-----------------------------------------------------------------------------------
    write(unit, "(A)") "'"//trim(id)//"'"
    write(unit, "(A,I0,A,I0)") "'", iso, "' ", n
  END SUBROUTINE write_ll_header

  !-------------------------------------------------------------------------------------
  SUBROUTINE write_ll_record(unit, lam_air, evi, loggf, gv, ga, gq, i, j, strength)
    integer, intent(in) :: unit, i, j
    real(8), intent(in) :: lam_air, evi, loggf, gv, ga, gq, strength
    !-----------------------------------------------------------------------------------
    write(unit, "(F12.4, F10.4, F9.3, F14.6, 2(ES12.3), 2(I6), ES12.3)") &
      lam_air, evi, loggf, gv, ga, gq, i, j, strength
  END SUBROUTINE write_ll_record

  !=====================================================================================
  ! All linelist_out.txt blocks of one composition particle: for each isotope (main
  ! first), a header plus the detailed AoS records (incl. hydrogen) and the bulk SoA
  ! records that were used. ll is the per-isotope SoA array (0 = main, k = isotope
  ! iso_values(k)); the AoS lines still carry their own iso_idx into iso_values.
  !-------------------------------------------------------------------------------------
  SUBROUTINE print_ll_particle(unit, id, iso_values, ll, lines)
    integer,             intent(in) :: unit
    character(len=*),    intent(in) :: id
    integer,             intent(in) :: iso_values(:)
    type(linelist_type), intent(in) :: ll(0:)
    type(bbtrans_type),  intent(in), optional :: lines(:)
    integer :: l, kiso, isoval, n
    !-----------------------------------------------------------------------------------
    do kiso=0, ubound(ll,1)
      isoval = 0
      if (kiso > 0) isoval = iso_values(kiso)
      ! first pass: block size (the header carries the record count)
      n = 0
      if (present(lines)) then
        do l=1, size(lines)
          if (lines(l)%active .and. lines(l)%used .and. lines(l)%iso_idx == kiso) n = n+1
        end do
      endif
      do l=1, ll(kiso)%nline
        if (ll(kiso)%active(l) .and. ll(kiso)%used(l)) n = n + 1
      end do
      if (n == 0) cycle

      call write_ll_header(unit, id, isoval, n)
      ! detailed (AoS) transitions, incl. hydrogen
      if (present(lines)) then
        do l=1, size(lines)
          associate(line => lines(l))
            if (line%active .and. line%used .and. line%iso_idx == kiso) then
              call write_ll_record(unit, line%lam_air, line%lower_ev, line%loggf, &
                line%gv_in, line%ga_in, line%gq_in, line%i, line%j, line%strength)
            endif
          end associate
        end do
      endif
      ! bulk (SoA) linelist lines
      do l=1, ll(kiso)%nline
        associate(sll => ll(kiso))
          if (.not. (sll%active(l) .and. sll%used(l))) cycle
          call write_ll_record(unit, vac2air(exp(sll%llam0(l))), &
            dble(sll%lower_ev(l)), dble(sll%loggf(l)), dble(sll%gv_in(l)), &
            dble(sll%ga_in(l)), dble(sll%gq_in(l)), sll%lev_i(l), sll%lev_j(l), &
            dble(sll%strength(l)))
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
            dble(ll%strength(l)))
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
    call grow4(ll%vdw_ratio, ll%nline, need)
    call grow4(ll%strength, ll%nline, need)
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

END MODULE linelist_mod
