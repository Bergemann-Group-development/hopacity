!=======================================================================================
! linelist_mod -- structure-of-arrays (SoA) bulk linelist handling.
!
! This is a fresh, performance-oriented reimplementation of the bulk metal/molecular
! linelist opacity that lives, in array-of-structures (AoS) form, as particle%line(:)
! (bbtrans_type) + composition%add_linelist_opacities. It is intentionally standalone:
! it only references low-level modules (prec/const/handyfuncs/op_math/bbtrans), never
! particle_type or composition_type, so that particle_type may hold a linelist_type
! component without creating a circular module dependency.
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
! The orchestration that walks elements/ions/molecules and gathers the per-particle
! thermodynamic arrays stays in composition_mod (add_linelist_opacities2 / init_linelist2);
! it must, because it needs the composition types.
!=======================================================================================
MODULE linelist_mod
  USE prec_mod
  USE const_mod
  USE trace_mod
  USE handyfuncs_mod, only: air2vac
  USE bbtrans_mod, only: voigt_function2, disable_ga, disable_stark, disable_vdw, &
                         f_conv, tol => tolerance
  implicit none
  private

  public :: linelist_type, build_linelist, finalize_linelist, add_opacity

  !-------------------------------------------------------------------------------------
  type :: linelist_type
    integer :: nline = 0
    ! line centres
    real(dp), allocatable :: nu0(:)      ! rest frequency [Hz]
    real(dp), allocatable :: lam0(:)     ! vacuum wavelength [AA]
    real(dp), allocatable :: llam0(:)    ! log(lam0) (natural log) -- sort key for the march
    real(dp), allocatable :: lower_ev(:) ! lower level energy [eV]
    ! oscillator strength + Einstein prefactor folded into one log term, so
    ! Rij = inv_part*idf*exp(gpref - ex), with gpref = loggf*ln10 + log(bij_const*frac/nu0)
    real(dp), allocatable :: gpref(:)
    ! broadening (already scaled by over_4pi/nu0, exactly as the AoS prepare_* routines)
    real(dp), allocatable :: ga(:)       ! radiative damping
    real(dp), allocatable :: gw1(:)      ! van der Waals coefficient
    real(dp), allocatable :: gw2(:)      ! van der Waals exponent
    real(dp), allocatable :: gq(:)       ! Stark constant (sign encodes recipe)
    logical,  allocatable :: unsold(:)   ! .true. -> Unsold vdW, .false. -> Barklem
    ! bookkeeping
    integer,  allocatable :: iso_idx(:)  ! isotope index (0 = main isotope)
    integer,  allocatable :: lev_i(:)    ! NLTE lower level index (filled at opacity time)
    integer,  allocatable :: lev_j(:)    ! NLTE upper level index
    logical,  allocatable :: active(:)   ! line is active
    logical,  allocatable :: used(:)     ! diagnostic: line contributed somewhere
    real(dp), allocatable :: strength(:) ! diagnostic: max opacity/continuum at centre
  end type linelist_type

  ! Bij = bij_const * (10**loggf) / nu0  (see init_einstein_coefficients in bbtrans_mod)
  real(dp), parameter :: bij_const = c2_2h * f_conv * h_4pi
  real(dp), parameter :: ex_cutoff = -85.0_dp  ! skip cells where exp(-ex) underflows

contains

  !=====================================================================================
  ! Parse + set up one per-species block of linelist records into the SoA arrays.
  ! The block's records live in the byte buffer buf; lstart(i)/lend(i) give the first
  ! and last byte of record i (so no per-record string copy is made). Records are
  ! assumed pre-sorted ascending in wavelength. If do_mask, only records with lam_air
  ! in [lam_min, lam_max] are kept (contiguous slice, located by bisection).
  !-------------------------------------------------------------------------------------
  SUBROUTINE build_linelist(ll, buf, lstart, lend, n, awgt, ion_j, xion, is_mol, &
                            iso_idx, iso_frac, vbar2e4, do_mask, lam_min, lam_max)
    type(linelist_type), intent(inout) :: ll
    character(len=*),    intent(in)    :: buf
    integer(8),          intent(in)    :: lstart(:), lend(:)
    integer,             intent(in)    :: n, ion_j, iso_idx
    real(dp),            intent(in)    :: awgt, xion, iso_frac, vbar2e4
    logical,             intent(in)    :: is_mol, do_mask
    real(dp),            intent(in)    :: lam_min, lam_max
    real(dp) :: lam_air, lower_ev, loggf, gv_in, ga_in, gq_in, r_i, r_j
    real(dp) :: nu0, lam0
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
      ll%lam0(k)     = lam0
      ll%llam0(k)    = log(lam0)
      ll%nu0(k)      = nu0
      ll%lower_ev(k) = lower_ev
      ll%gpref(k)    = loggf*ln10 + log(bij_const * iso_frac / nu0)
      ll%iso_idx(k)  = iso_idx
      ll%lev_i(k)    = nint(r_i)
      ll%lev_j(k)    = nint(r_j)
      ll%active(k)   = .true.
      ll%used(k)     = .false.
      ll%strength(k) = 0

      call prepare_radiative   (ll%ga(k),  ga_in, nu0, loggf)
      call prepare_stark       (ll%gq(k),  gq_in, nu0, lower_ev, ion_j, xion)
      call prepare_vanderwaals (ll%gw1(k), ll%gw2(k), ll%unsold(k), gv_in, nu0, &
                                lower_ev, ion_j, xion, vbar2e4, is_mol)
    end do
    call trace%end(it_parse)
  END SUBROUTINE build_linelist

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
    call permute8(ll%nu0,      idx)
    call permute8(ll%lam0,     idx)
    call permute8(ll%llam0,    idx)
    call permute8(ll%lower_ev, idx)
    call permute8(ll%gpref,    idx)
    call permute8(ll%ga,       idx)
    call permute8(ll%gw1,      idx)
    call permute8(ll%gw2,      idx)
    call permute8(ll%gq,       idx)
    call permute8(ll%strength, idx)
    call permutei(ll%iso_idx,  idx)
    call permutei(ll%lev_i,    idx)
    call permutei(ll%lev_j,    idx)
    call permutel(ll%unsold,   idx)
    call permutel(ll%active,   idx)
    call permutel(ll%used,     idx)
    call trace%end(it_permute)
  END SUBROUTINE finalize_linelist

  !=====================================================================================
  ! Opacity kernel: spatial-outer, line, frequency-inner; log-lambda march.
  !
  !   chi, eta  : (nf, m1, m2, m3)  accumulated into (the "local" work arrays)
  !   cont      : (nf, m1, m2, m3)  continuum reference for the tolerance window
  !   llam      : (nf)              log(lambda) of the (ascending) wavelength grid
  !   inv_df    : (m1,m2,m3, 0:niso) inverse Doppler width c/v_th, indexed by isotope
  !   vbar2     : (m1,m2,m3, 0:niso) mean square velocity, indexed by isotope
  !   per-cell  : vproj, e_kT, dxij, inv_part, ne, tothi, stark_pe
  !   dep_mem   : (m1,m2,m3, nlevel) NLTE departure coefficients (or size 0 if LTE)
  !   do_nlte   : whether dep_mem applies to this particle
  !-------------------------------------------------------------------------------------
  SUBROUTINE add_opacity(ll, llam, chi, eta, cont, vproj, inv_df, vbar2, e_kT, dxij, &
                         inv_part, ne, tothi, stark_pe, dep_mem, do_nlte)
    type(linelist_type), intent(inout) :: ll
    real(dp), dimension(:),       intent(in)    :: llam
    real(4),  dimension(:,:,:,:), intent(inout) :: chi, eta
    real(4),  dimension(:,:,:,:), intent(in)    :: cont
    real(4),  dimension(:,:,:),   intent(in)    :: vproj, e_kT, dxij, inv_part, ne, tothi, stark_pe
    real(4),  dimension(:,:,:,0:), intent(in)   :: inv_df, vbar2
    real(4),  dimension(:,:,:,:), intent(in)    :: dep_mem
    logical,  intent(in) :: do_nlte
    integer  :: i, j, k, l, f, f0, fc, nf, m(3), iso, nl
    real(8)  :: llc
    real(4)  :: idf0, idf, vb2, unsold_c, ekt, dxij_c, invpart, ne_c, tothi_c, starkpe, vproj_c
    real(4)  :: a, ex, Rij, cchi, ceta, ref, vs, phi, lchi, leta
    real(8)  :: arg
    logical  :: nlte, asc
    !-----------------------------------------------------------------------------------
    m  = shape(vproj)
    nf = size(llam)
    nl = ll%nline
    if (nl == 0) return
    ! the march is a monotonic cursor; it works for either grid direction (lam is
    ! normally ascending, i.e. nu descending) as long as the grid is monotonic.
    asc = llam(nf) >= llam(1)

    do k = 1, m(3)
      do j = 1, m(2)
        do i = 1, m(1)
          idf0     = inv_df(i,j,k,0)
          unsold_c = vbar2(i,j,k,0)**0.3 * tothi(i,j,k)   ! Unsold uses the main isotope
          ekt      = e_kT(i,j,k)
          dxij_c   = dxij(i,j,k)
          invpart  = inv_part(i,j,k)
          ne_c     = ne(i,j,k)
          tothi_c  = tothi(i,j,k)
          starkpe  = stark_pe(i,j,k)
          vproj_c  = vproj(i,j,k)
          if (asc) then; fc = 1; else; fc = nf; end if

          do l = 1, nl
            if (.not. ll%active(l)) cycle

            ! shifted line centre in log(lambda); monotonic in l (vproj_c fixed),
            ! so the cursor fc only moves one way through the grid.
            llc = ll%llam0(l) - vproj_c
            if (asc) then
              do while (fc < nf)
                if (llam(fc) >= llc) exit
                fc = fc + 1
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
            else
              do while (fc > 1)
                if (llam(fc) >= llc) exit
                fc = fc - 1
              end do
              if (fc < nf) then
                if (abs(llam(fc) - llc) < abs(llam(fc+1) - llc)) then
                  f0 = fc
                else
                  f0 = fc + 1
                end if
              else
                f0 = nf
              end if
            end if

            ! exponent of the folded Boltzmann * gf * Einstein prefactor; the cutoff
            ! now accounts for line strength (gpref), skipping negligible lines before
            ! the damping and Voigt work (the tolerance test below still backstops).
            ex  = (ll%lower_ev(l) - dxij_c) * ekt
            arg = invpart + ll%gpref(l) - ex
            if (arg < ex_cutoff) cycle

            iso = ll%iso_idx(l)
            if (iso == 0) then
              idf = idf0
              vb2 = vbar2(i,j,k,0)
            else
              idf = inv_df(i,j,k,iso)
              vb2 = vbar2(i,j,k,iso)
            end if

            ! damping parameter (physical, then scaled by idf as in the AoS kernel)
            a = ll%ga(l)
            if (ll%unsold(l)) then
              a = a + ll%gw1(l) * unsold_c
            else
              a = a + ll%gw1(l) * vb2**ll%gw2(l) * tothi_c
            end if
            if (ll%gq(l) >= 0.0) then
              a = a + ll%gq(l) * ne_c
            else
              a = a - ll%gq(l) * starkpe
            end if
            a = a * idf

            Rij = idf * exp(arg)

            nlte = do_nlte .and. ll%lev_i(l) > 0 .and. ll%lev_j(l) > 0
            if (nlte) then
              cchi = Rij * dep_mem(i,j,k, ll%lev_i(l))
              ceta = Rij * dep_mem(i,j,k, ll%lev_j(l))
            else
              cchi = Rij
              ceta = Rij
            end if

            ref = tol * cont(f0,i,j,k)

            ! centre
            vs   = idf * real(llam(f0) - llc)
            phi  = voigt_function2(a, vs)
            lchi = phi * cchi
            leta = phi * ceta
            if (lchi <= ref) cycle

            ll%used(l)     = .true.
            ll%strength(l) = max(ll%strength(l), real(lchi / cont(f0,i,j,k), dp))
            chi(f0,i,j,k)  = chi(f0,i,j,k) + lchi
            eta(f0,i,j,k)  = eta(f0,i,j,k) + leta

            ! red wing
            f = f0 + 1
            do while (lchi >= ref .and. f <= nf)
              vs   = idf * real(llam(f) - llc)
              phi  = voigt_function2(a, vs)
              lchi = phi * cchi
              leta = phi * ceta
              chi(f,i,j,k) = chi(f,i,j,k) + lchi
              eta(f,i,j,k) = eta(f,i,j,k) + leta
              f = f + 1
            end do

            ! blue wing
            f    = f0 - 1
            lchi = ref
            do while (lchi >= ref .and. f >= 1)
              vs   = idf * real(llam(f) - llc)
              phi  = voigt_function2(a, vs)
              lchi = phi * cchi
              leta = phi * ceta
              chi(f,i,j,k) = chi(f,i,j,k) + lchi
              eta(f,i,j,k) = eta(f,i,j,k) + leta
              f = f - 1
            end do
          end do
        end do
      end do
    end do
  END SUBROUTINE add_opacity

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
    elseif (ga_in < 100) then
      ga = ga_in * 1e8
    else
      ga = ga_in
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
                                      cnt_ev, vbar2e4, is_mol)
    real(dp), intent(out) :: gw1, gw2
    logical,  intent(out) :: unsold
    real(dp), intent(in)  :: gv, nu0, lower_ev, cnt_ev, vbar2e4
    integer,  intent(in)  :: lower_ion
    logical,  intent(in)  :: is_mol
    real(dp) :: upper_ev, sigma, alpha, gx, gammaf, balfac, c1, c625
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
    call grow8(ll%nu0,      ll%nline, need)
    call grow8(ll%lam0,     ll%nline, need)
    call grow8(ll%llam0,    ll%nline, need)
    call grow8(ll%lower_ev, ll%nline, need)
    call grow8(ll%gpref,    ll%nline, need)
    call grow8(ll%ga,       ll%nline, need)
    call grow8(ll%gw1,      ll%nline, need)
    call grow8(ll%gw2,      ll%nline, need)
    call grow8(ll%gq,       ll%nline, need)
    call grow8(ll%strength, ll%nline, need)
    call growi(ll%iso_idx,  ll%nline, need)
    call growi(ll%lev_i,    ll%nline, need)
    call growi(ll%lev_j,    ll%nline, need)
    call growl(ll%unsold,   ll%nline, need)
    call growl(ll%active,   ll%nline, need)
    call growl(ll%used,     ll%nline, need)
  END SUBROUTINE grow

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
