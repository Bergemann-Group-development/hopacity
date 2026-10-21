MODULE hprofl_mod
  USE bbtrans_mod, only: bbtrans_type, disable_rad, disable_stk, disable_vdw
  USE qlg, only: qlg_bot_coef_t, qlg_bot_coef, qlg_bot_eval   ! the one Voigt kernel; see linelist_mod's kernel note
  USE prec_mod, only: RYD, dp
  USE const_mod
  USE hvcs_mod, only: vcs_supported, vcs_resonance, vcs_radiative, ensure_vcs_loaded, &
                      vcs_profile_t, vcs_prepare, vcs_eval, vcs_tail
  implicit none

  real(4), dimension(4,3), parameter ::  xknmtb = reshape([&
    0.0001716, 0.0090190, 0.1001000, 0.5820000,  &
    0.0005235, 0.0177200, 0.1710000, 0.8660000,  &
    0.0008912, 0.0250700, 0.2230000, 1.0200000], [4,3])


  !  Einstein A-value sums for H lines
  real(4), dimension(100), parameter :: asum = [&
    0.000E+00, 4.696E+08, 9.980E+07, 3.017E+07, 1.155E+07, 5.189E+06,  &
    2.616E+06, 1.437E+06, 8.444E+05, 5.234E+05, 3.389E+05, 2.275E+05,  &
    1.575E+05, 1.120E+05, 8.142E+04, 6.040E+04, 4.560E+04, 3.496E+04,  &
    2.719E+04, 2.141E+04, 1.711E+04, 1.377E+04, 1.119E+04, 9.166E+03,  &
    7.572E+03, 6.341E+03, 5.338E+03, 4.523E+03, 3.854E+03, 3.302E+03,  &
    2.844E+03, 2.460E+03, 2.138E+03, 1.866E+03, 1.635E+03, 1.438E+03,  &
    1.269E+03, 1.124E+03, 9.983E+02, 8.894E+02, 7.947E+02, 7.120E+02,  &
    6.396E+02, 5.759E+02, 5.198E+02, 4.703E+02, 4.263E+02, 3.873E+02,  &
    3.526E+02, 3.215E+02, 2.938E+02, 2.689E+02, 2.465E+02, 2.264E+02,  &
    2.082E+02, 1.918E+02, 1.769E+02, 1.634E+02, 1.512E+02, 1.400E+02,  &
    1.298E+02, 1.206E+02, 1.121E+02, 1.043E+02, 9.720E+01, 9.066E+01,  &
    8.465E+01, 7.912E+01, 7.403E+01, 6.933E+01, 6.498E+01, 6.097E+01,  &
    5.725E+01, 5.381E+01, 5.061E+01, 4.765E+01, 4.489E+01, 4.232E+01,  &
    3.994E+01, 3.771E+01, 3.563E+01, 3.369E+01, 3.188E+01, 3.019E+01,  &
    2.860E+01, 2.712E+01, 2.572E+01, 2.442E+01, 2.319E+01, 2.204E+01,  &
    2.096E+01, 1.994E+01, 1.898E+01, 1.808E+01, 1.722E+01, 1.642E+01,  &
    1.566E+01, 1.495E+01, 1.427E+01, 1.363E+01]

  !  For Lyman lines only the s-p transition is allowed.
  real(4), dimension(100), parameter :: asumlyman = [&
    0.000E+00, 6.265E+08, 1.897E+08, 8.126E+07, 4.203E+07, 2.450E+07,  &
    1.236E+07, 8.249E+06, 5.782E+06, 4.208E+06, 3.158E+06, 2.430E+06,  &
    1.910E+06, 1.567E+06, 1.274E+06, 1.050E+06, 8.752E+05, 7.373E+05,  &
    6.269E+05, 5.375E+05, 4.643E+05, 4.038E+05, 3.534E+05, 3.111E+05,  &
    2.752E+05, 2.447E+05, 2.185E+05, 1.959E+05, 1.763E+05, 1.593E+05,  &
    1.443E+05, 1.312E+05, 1.197E+05, 1.094E+05, 1.003E+05, 9.216E+04,  &
    8.489E+04, 7.836E+04, 7.249E+04, 6.719E+04, 6.239E+04, 5.804E+04,  &
    5.408E+04, 5.048E+04, 4.719E+04, 4.418E+04, 4.142E+04, 3.888E+04,  &
    3.655E+04, 3.440E+04, 3.242E+04, 3.058E+04, 2.888E+04, 2.731E+04,  &
    2.585E+04, 2.449E+04, 2.322E+04, 2.204E+04, 2.094E+04, 1.991E+04,  &
    1.894E+04, 1.804E+04, 1.720E+04, 1.640E+04, 1.566E+04, 1.496E+04,  &
    1.430E+04, 1.368E+04, 1.309E+04, 1.254E+04, 1.201E+04, 1.152E+04,  &
    1.105E+04, 1.061E+04, 1.019E+04, 9.796E+03, 9.419E+03, 9.061E+03,  &
    8.721E+03, 8.398E+03, 8.091E+03, 7.799E+03, 7.520E+03, 7.255E+03,  &
    7.002E+03, 6.760E+03, 6.530E+03, 6.310E+03, 6.100E+03, 5.898E+03,  &
    5.706E+03, 5.522E+03, 5.346E+03, 5.177E+03, 5.015E+03, 4.860E+03,  &
    4.711E+03, 4.569E+03, 4.432E+03, 4.300E+03]

  !  Fine structure components for alpha lines in nu*10**-7
  real(4), dimension(34), parameter :: stalph = [&
    -730.,  370.,  188.,  515.,  327.,  619., -772.,  &
    -473., -369.,  120.,  256.,  162.,  285., -161.,  &
    -38.3,  6.82, -174., -147., -101., -77.5,   55.,  &
    126.,   75.,  139.,  -60.,   3.7,   27.,  -69.,  &
    -42.,  -18.,  -5.5,  -9.1,  -33.,  -24.]*1e7

  !  Alpha component weights
  real(4), dimension(34), parameter :: stwtal = [&
    1.,2.,1.,2.,1.,2.,1.,2.,3.,1.,2.,1.,2.,1.,4.,6.,1.,  &
    2.,3.,4.,1.,2.,1.,2.,1.,4.,6.,1.,7.,6.,4.,4.,4.,5.]

  !  Fine structure for M.EQ.INFINITY IN nu*10**-7
  real(4), dimension(5,4), parameter :: stcomp = reshape([&
    0.,    0.,    0.,    0.,    0.,  &
    468.,  576., -522.,    0.,    0.,  &
    260.,  290.,  -33., -140.,    0.,  &
    140.,  150.,   18.,  -27.,  -51.], [5,4])*1e7

  !  Weights for fine structure components
  real(4), dimension(5,4), parameter :: stcpwt = reshape([&
    1., 0., 0., 0., 0., &
    1., 1., 2., 0., 0., &
    1., 1., 4., 3., 0., &
    1., 1., 4., 6., 4.], [5,4])



  !  Lyman alpha red wing profiles following Allard et al (1998, A&A 335,
  !  1124). Taken from ATLAS12 code (see Castelli & Kurucz 2001, A&A 372,
  !  260)

  !  Store log10(I(d omega)) for N(H) = 1e14 cm^-3.

  !     LYMAN ALPHA QUASI H2+ PROFILE
  !     DELTA WAVENO =  -15000+100*(N-1) N=1,111   UP TO -4000
  real(8), dimension(111), parameter :: lymanh2plus = [&
    -15.14, -15.06, -14.97, -14.88, -14.80, -14.71, -14.62, -14.53,  &
    -14.44, -14.36, -14.27, -14.18, -14.09, -14.01, -13.92, -13.83,  &
    -13.74, -13.65, -13.57, -13.48, -13.39, -13.30, -13.21, -13.13,  &
    -13.04, -12.95, -12.86, -12.77, -12.69, -12.60, -12.51, -12.40,  &
    -12.29, -12.15, -12.02, -11.90, -11.76, -11.63, -11.53, -11.41,  &
    -11.30, -11.22, -11.15, -11.09, -11.07, -11.06, -11.07, -11.09,  &
    -11.12, -11.16, -11.19, -11.21, -11.24, -11.27, -11.30, -11.33,  &
    -11.36, -11.39, -11.42, -11.45, -11.48, -11.48, -11.48, -11.48,  &
    -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48,  &
    -11.48, -11.48, -11.48, -11.48, -11.41, -11.40, -11.39, -11.38,  &
    -11.37, -11.36, -11.35, -11.34, -11.33, -11.32, -11.30, -11.29,  &
    -11.28, -11.27, -11.27, -11.27, -11.26, -11.25, -11.24, -11.23,  &
    -11.22, -11.21, -11.20, -11.19, -11.18, -11.17, -11.15, -11.14,  &
    -11.13, -11.12, -11.11, -11.10, -11.09, -11.08, -11.07] - 14.

  !     LYMAN ALPHA QUASI H2 PROFILE
  !     DELTA WAVENO = -22000+200*(N-1)  N=1,91  -4000
  real(8), dimension(91), parameter :: lymanh2 = [&
    -13.43, -13.32, -13.21, -13.10, -12.98, -12.86, -12.79, -12.72,  &
    -12.65, -12.58, -12.51, -12.47, -12.45, -12.45, -12.48, -12.51,  &
    -12.53, -12.56, -12.59, -12.62, -12.65, -12.69, -12.73, -12.77,  &
    -12.81, -12.85, -12.87, -12.89, -12.90, -12.90, -12.90, -12.90,  &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,  &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,  &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,  &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,  &
    -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.89, -12.88,  &
    -12.87, -12.86, -12.85, -12.84, -12.83, -12.81, -12.80, -12.79,  &
    -12.78, -12.76, -12.74, -12.72, -12.70, -12.68, -12.65, -12.62,  &
    -12.59, -12.56, -12.53] - 14.

  !  Data for self-broadening from calculations of Barklem, Piskunov and
  !  O'Mara (2000, A&A 363, 1091).
  !  SIGMA is in m^2.
  !  BALFAC = (4/Pi)**(alf/2) Gamma(2-2/alf) the leading factor for
  !  computing the width. 1E6 factor is due to SI -> CGS.
  real(4), dimension(3), parameter :: sigma  = [3.304E-18, 6.497E-18, 1.178E-17]*1.e4
  real(4), dimension(3), parameter :: alf    = [    0.677,     0.455,     0.380]
  real(4), dimension(3), parameter :: balfac = [  0.97875,   0.97659,   0.97794]*1.e6

  integer, dimension(4), parameter :: istal  = [1, 3, 10, 21]
  integer, dimension(4), parameter :: lnghal = [2, 7, 11, 14]
  integer, dimension(4), parameter :: lncomp = [1, 3,  4,  5]

  logical, public :: hydrogen_BPO=.true.

  ! Profile parameters of one hydrogen line (see hydrogen_line). The widths carry the
  ! over_4pi/nu0 scaling of the bbtrans_type broadening constants; dq(1:io), wt(1:io)
  ! are the fine-structure offsets [nu/nu0] and weights, dqbar their weighted mean and
  ! dqmax the largest |dq|. bpo marks the lines whose self-broadening follows Barklem,
  ! Piskunov & O'Mara (2000) instead of Ali & Griem.
  type, public :: hline_t
    integer  :: n=0, m=0          ! lower and upper principal quantum numbers
    real(dp) :: nu0=0             ! rest frequency [Hz]
    real(dp) :: grad=0, gstk=0, gw1=0, gw2=0, gw3=0
    logical  :: bpo=.false.
    integer  :: io=0
    real(dp) :: dq(14)=0, wt(14)=0
    real(dp) :: dqbar=0, dqmax=0
  end type hline_t

  ! The q-independent part of stark_eval for one line in one cell (stark_prepare)
  type :: stark_cell_t
    logical  :: off = .true.      ! no Stark broadening (hwstk = 0)
    integer  :: n = 0, m = 0
    real(8)  :: hw5, pp, cne, ex, c1coef, g1, gcon1, gcon2, gam0
  end type stark_cell_t

  ! One hydrogen line in one cell, prepared by hcell_prepare for hprofile
  type, public :: hcell_t
    real(4) :: idf, ne, tothi, resont, vdw, hwstk, hwlor, hfwid
    real(4) :: hwadd = 0                  ! table path: Lorentzian half width added to the table
    logical :: use_vcs = .false., lymanalf = .false.
    type(stark_cell_t)   :: st
    type(qlg_bot_coef_t) :: vcore, vwing   ! Voigt coefficients inside / outside hfwid
    type(vcs_profile_t)  :: vcs
    real(8) :: xholt4000 = 0              ! Lyman alpha: half the Stark profile at q4000
  end type hcell_t

  ! Lyman alpha wing regions, in relative detuning q
  real(8), parameter :: LyA0 = 82259.105 ! wavenumber
  real(8), parameter :: rayqcut = (1e8/1400.d0 - LyA0)/LyA0 !Changed from 1240 to 1400 by BPz 04/10-2007
  real(8), parameter :: redthresh = 20000./LyA0 - 1d0
  real(8), parameter :: q4000 = -4000./LyA0

contains

  !====================================================================================
  SUBROUTINE add_h_opacities(h, Bij, f0, nu, vproj, inv_df, temp, ni, nj, &
    ne, tothi, vbar2h, fo, t3nhe, bi, bj, chi, eta, cont, do_nlte, line_used, &
    line_strength, i0, i1)
    type(hline_t), intent(in) :: h
    real(dp), intent(in) :: Bij                 ! Einstein Bij, incl. the h/4pi factor
    ! diagnostics: set when the line is evaluated in some cell; the maximum centre
    ! line/continuum opacity ratio
    logical, intent(inout) :: line_used
    real(4), intent(inout) :: line_strength
    integer, intent(in) :: f0
    real(8), dimension(:), intent(in), contiguous :: nu
    real(4), dimension(:,:,:), intent(in), contiguous :: vproj, inv_df, ni, nj
    real(4), dimension(:,:,:), intent(in), pointer :: bi, bj
    real(4), dimension(:,:,:), intent(in), contiguous :: temp, ne, tothi, vbar2h, fo, t3nhe
    real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
    real(4), dimension(:,:,:,:), intent(in) :: cont
    logical, intent(in) :: do_nlte
    integer, intent(in), optional :: i0, i1
    real(4) :: lchi, leta, phi, idf, cchi, ceta, v, ref, strength, nup, nlo
    real(4) :: inv_nu0, Rij
    real(8) :: q, qt, vp
    real(4), parameter :: tol=1e-7
    integer :: i,j,k,f, nf, m(3), lo, hi, io, mem_offset
    logical used, tail
    type(hcell_t) :: hc                  ! the line in the current cell
    !  ln|d| and 1/d for d = (nu-nu0)/nu0 - dqbar, the tail detuning at rest. They are
    !  the same in every cell, so they are filled once per line, lazily over the
    !  contiguous range f_lo:f_hi that the frequency loops have reached.
    real(8), allocatable :: lnd(:), invd(:)
    integer :: f_lo, f_hi
    !----------------------------------------------------------------------------------
    m = shape(vproj)
    nf = size(nu)
    inv_nu0 = 1.0 / h%nu0
    used = .false.
    strength = 0.0
    lo = 1;  hi = m(1)
    if (present(i0)) lo = i0
    if (present(i1)) hi = i1
    mem_offset = lo - 1
    tail = vcs_supported(h%n, h%m)
    if (tail) allocate(lnd(nf), invd(nf))
    f_lo = f0;  f_hi = f0 - 1
    if (hi - lo + 1 /= m(1)) &
      error stop 'hprofl add_h_opacities: per-cell arrays are not sized to the cell range'
    if (size(chi, 2) /= m(1)) &
      error stop 'hprofl add_h_opacities: chi/eta extent does not match the cell range'

    do k=1, m(3)
      do j=1, m(2)
        do i=lo, hi
          io = i - mem_offset            ! chunk-local index; cont and bi/bj stay whole-tile at i
          ! fetch variables
          idf = inv_df(io,j,k)
          v   = (vproj(io,j,k) + 1)*h%nu0
          nlo = ni(io,j,k)
          nup = nj(io,j,k)

          ! the idf normalisation is applied to the Voigt terms inside hprofile, since
          ! the Stark profile is normalised differently; Rij does not carry it
          Rij = nlo * Bij
          if (do_nlte) then
            ceta = Rij * bj(i,j,k)
            cchi = Rij * bi(i,j,k)
          else
            ceta = Rij
            cchi = Rij
          endif

          call hcell_prepare(h, temp(io,j,k), ne(io,j,k), fo(io,j,k), vbar2h(io,j,k), &
                             tothi(io,j,k), t3nhe(io,j,k), idf, hc)
          qt = huge(qt)
          if (tail .and. hc%use_vcs) qt = hc%vcs%q_tail + h%dqmax
          vp = vproj(io,j,k)
          ref = tol * cont(f0,i,j,k)

          f = f0
          lchi = ref
          do while (lchi >= ref .and. f <= nf) ! if line opacity below tolerance exit
            q = (nu(f) - v) * inv_nu0
            phi = profile_at(f, q)
            lchi = phi * cchi
            leta = phi * ceta
            chi(f,io,j,k) = chi(f,io,j,k) + lchi
            eta(f,io,j,k) = eta(f,io,j,k) + leta
            if (f==f0) then
              used = .true.
              strength = max(strength, lchi / cont(f0,i,j,k))
            endif
            f = f + 1
          end do

          f = f0-1
          lchi = ref
          do while (lchi >= ref .and. f >= 1) ! if line opacity below tolerance exit
            q = (nu(f) - v) * inv_nu0
            phi = profile_at(f, q)
            lchi = phi * cchi
            leta = phi * ceta
            chi(f,io,j,k) = chi(f,io,j,k) + lchi
            eta(f,io,j,k) = eta(f,io,j,k) + leta
            f = f - 1
          end do
        end do
      end do
    end do
    if (used) then
      line_used = .true.
      line_strength = max(line_strength, strength)
    endif

  contains

    !  The profile at frequency index f, detuning q. Past qt every fine-structure
    !  component is on the VCS power-law tail (see vcs_fine), evaluated at
    !  ln|q - dqbar| = ln|d| + ln(1 - vp/d). |vp/d| < 1e-2 there, so the series to
    !  second order leaves an error below 1e-6 in ln|q|.
    FUNCTION profile_at (f, q) result(phi)
      integer, intent(in) :: f
      real(8), intent(in) :: q
      real(4) :: phi
      real(8) :: r
      !--------------------------------------------------------------------------------
      if (abs(q) > qt) then
        if (f > f_hi .or. f < f_lo) call extend(f)
        r = vp*invd(f)
        phi = real(vcs_tail(hc%vcs, lnd(f) - r*(1.d0 + 0.5d0*r)), 4)
        if (hc%hwadd > 0) phi = phi + added_wing(hc, q - h%dqbar)
      else
        phi = hprofile(h, hc, q)
      endif
    END FUNCTION profile_at

    SUBROUTINE extend (f)
      integer, intent(in) :: f
      integer :: g
      !--------------------------------------------------------------------------------
      do g=f_hi+1, f
        call fill(g)
      end do
      do g=f, f_lo-1
        call fill(g)
      end do
      f_lo = min(f, f_lo)
      f_hi = max(f, f_hi)
    END SUBROUTINE extend

    SUBROUTINE fill (g)
      integer, intent(in) :: g
      real(8) :: d
      !--------------------------------------------------------------------------------
      d = (nu(g) - h%nu0)/h%nu0 - h%dqbar
      d = sign(max(abs(d), tiny(d)), d)
      lnd(g) = log(abs(d))
      invd(g) = 1.d0/d
    END SUBROUTINE fill
  END SUBROUTINE add_h_opacities

  !===============================================================================
  ! The per-cell part of a hydrogen line's profile: widths, the Stark setup, the Voigt
  ! coefficients of the core and wing branches, and the tabulated VCS profile. The
  ! caller keeps c across cells, so an unchanged VCS stamp is not rebuilt.
  !===============================================================================
  SUBROUTINE hcell_prepare (h, temp, ne, fo, vbar2h, tothi, t3nhe, idf, c)
    type(hline_t), intent(in) :: h
    real(4), intent(in) :: temp, ne, fo, vbar2h, tothi, t3nhe, idf
    type(hcell_t), intent(inout) :: c
    real(4) :: nh_table
    !---------------------------------------------------------------------------
    c%idf = idf
    c%ne = ne
    c%tothi = tothi
    c%resont = tothi * h%gw1
    if (h%gw2 /= 0) c%resont = c%resont * vbar2h**h%gw2
    c%hwstk = h%gstk * fo
    c%vdw = h%gw3 * t3nhe
    c%hwlor = h%grad + c%resont + c%vdw
    c%hfwid = max(1/idf, c%hwlor, c%hwstk) ! defines where we switch from Voigt to adding

    ! Tabulated profile (hvcs_mod): whole-profile replacement, with the same
    ! fine-structure superposition as the analytic path. Broadening the table
    ! lacks is added as one Lorentzian of half width hwadd outside
    ! hfwid = max(Doppler width, hwadd): radiative damping for series tabulated
    ! without it (Stehle & Hutcheon), broadening by neutral helium, which no table
    ! has, and the broadening by neutral hydrogen of a BPO line or of a line
    ! without an N(HI) axis. Tables with that axis carry Ali &
    ! Griem resonance broadening; a BPO line reads them at the lowest node
    ! (N(HI) = 1e12 cm^-3, where that width is below 1e-7 AA). Without broadening
    ! by neutral hydrogen (disable_vdw) resont is zero and the lowest node is used.
    c%use_vcs = vcs_supported(h%n, h%m)
    if (c%use_vcs) then
      c%hwadd = c%vdw
      if (.not. vcs_radiative(h%n)) c%hwadd = c%hwadd + real(h%grad, 4)
      if (h%bpo .or. .not. vcs_resonance(h%n, h%m)) c%hwadd = c%hwadd + c%resont
      nh_table = tothi
      if (h%bpo .or. disable_vdw) nh_table = 0
      c%hfwid = max(1/idf, c%hwadd)
      call vcs_prepare(h%n, h%m, h%nu0, temp, ne, nh_table, c%vcs)
      return
    endif

    c%lymanalf = (h%n == 1) .and. (h%m == 2)
    call stark_prepare(h%n, h%m, h%nu0, temp, ne, c%hwstk, c%st)
    c%vcore = qlg_bot_coef((c%hwlor + c%hwstk) * idf)
    c%vwing = qlg_bot_coef(c%hwlor * idf)
    if (c%lymanalf) c%xholt4000 = 0.5*stark_eval(c%st, q4000)
  END SUBROUTINE hcell_prepare

  !===============================================================================
  !  The profile at q = (nu - nu0)/nu0 in a cell prepared by hcell_prepare.
  !  Based on code by Deane Peterson and Bob Kurucz
  !===============================================================================
  FUNCTION hprofile (h, c, q) result(phi)
    type(hline_t), intent(in) :: h
    type(hcell_t), intent(in) :: c
    real(8), intent(in) :: q
    real(4) :: phi
    !---------------------------------------------------------------------------
    if (c%use_vcs) then
      phi = vcs_fine(h, c, q)
    else if (c%lymanalf) then
      call lyman_voigt(h, c, q, phi)
    else
      phi = hvoigt(h, c, q)
    endif
  END FUNCTION hprofile

  !====================================================================================
  FUNCTION hvoigt (h, c, q) result(phi)
    type(hline_t), intent(in) :: h
    type(hcell_t), intent(in) :: c
    real(8), intent(in) :: q
    real(4) :: vs, xv, phi
    integer :: i
    !----------------------------------------------------------------------------------
    if (abs(q) > c%hfwid) then
      !  Treat stark broadening separately
      phi = stark_eval(c%st, q)
      !  Including fine structure
      do i=1, h%io
        vs = (q - h%dq(i)) * c%idf
        xv = qlg_bot_eval(c%vwing, vs) * c%idf
        phi = phi + h%wt(i)*xv
      end do
    else
      phi = 0.d0
      do i=1, h%io
        vs = (q - h%dq(i)) * c%idf
        xv = qlg_bot_eval(c%vcore, vs) * c%idf
        phi = phi + h%wt(i)*xv
      end do
    endif
  END FUNCTION

  !====================================================================================
  !  Tabulated VCS profile with the same fine-structure superposition as hvoigt:
  !  the prepared (symmetric, unit-area) table profile is shifted to each
  !  component dq(i) and weighted by wt(i). The weights sum to one, so the
  !  normalisation is preserved; lines with a single component reduce to
  !  vcs_eval at q. Mirrors SIU's HYDL_AG_CON_FS, with HOP's component pattern.
  !  Broadening the table lacks is added to each component as a Lorentzian
  !  (added_wing), as SIU's HYDL_BPO_CON_FS and the adding method of Barklem's
  !  HLINPROF do for BPO.
  !
  !  Once every component lies beyond the table's outer edge, the profile is the
  !  power law A*|q|**s, and the weighted sum equals one evaluation at the mean
  !  shift dqbar up to a relative s(s-1)/2 * var(dq)/q**2, below 1e-5 there. Most
  !  evaluations of the strong lines fall in this regime; add_h_opacities
  !  evaluates it directly, sharing ln|q| between cells.
  !====================================================================================
  FUNCTION vcs_fine (h, c, q) result(phi)
    type(hline_t), intent(in) :: h
    type(hcell_t), intent(in) :: c
    real(8), intent(in) :: q
    real(4) :: phi
    integer :: i
    !----------------------------------------------------------------------------------
    if (abs(q) - h%dqmax > c%vcs%q_tail) then
      phi = real(vcs_tail(c%vcs, log(abs(q - h%dqbar))), 4)
      if (c%hwadd > 0) phi = phi + added_wing(c, q - h%dqbar)
      return
    endif
    phi = 0.0
    do i=1, h%io
      phi = phi + h%wt(i)*real(vcs_eval(c%vcs, q - h%dq(i)), 4)
    end do
    if (c%hwadd > 0) then
      do i=1, h%io
        phi = phi + h%wt(i)*added_wing(c, q - h%dq(i))
      end do
    endif
  END FUNCTION

  !====================================================================================
  !  The added Lorentzian (unit area in q) at detuning d from a component, added
  !  to the table profile beyond hfwid = max(Doppler width, hwadd), the same
  !  switch between core and adding as in hvoigt. Inside it the Doppler core
  !  dominates, and leaving the Lorentzian out there keeps the summed area close
  !  to one; the excess, about 2/pi * hwadd/hfwid, is what the analytic path's
  !  adding also carries.
  !====================================================================================
  PURE FUNCTION added_wing (c, d) result(phi)
    type(hcell_t), intent(in) :: c
    real(8), intent(in) :: d
    real(4) :: phi
    !----------------------------------------------------------------------------------
    if (abs(d) > c%hfwid) then
      phi = real(c%hwadd / (PI * (d*d + c%hwadd*c%hwadd)), 4)
    else
      phi = 0.0
    endif
  END FUNCTION added_wing

  !====================================================================================
  SUBROUTINE lyman_voigt (h, c, q, phi)
    type(hline_t), intent(in) :: h
    type(hcell_t), intent(in) :: c
    real(8), intent(in) :: q
    real(4), intent(out) :: phi
    real(4) :: a, vs, xv
    real(4) :: xfac, xlymanh2
    real(8) :: xstark, sat4000
    type(qlg_bot_coef_t) :: vc
    integer :: i
    !----------------------------------------------------------------------------------
    if (abs(q) > c%hfwid) then
      a = 0
      !  Lyman alpha:
      if (q > rayqcut) then
        !  Natural broadening
        !  This is cut at 1400 AA in the red.
        !  Note that this is Rayleigh scattering and therefore the rayleigh scattering
        !  in composition_mod has to be cut below 1400 AA to avoid double counting.
        a = h%grad*(q+1d0)**2
      endif
      if (q > q4000) then
        !  Self broadening of Lyman alpha:
        !  Red wing of Lyman Alpha before the satellite region and blue wing
        !  (ie. redward of ~1190 A).
        !  The impact approximation breaks down quickly here validity is ~ 1 A
        !  at temp~5000 K (Lortet & Roueff 1969, A&A 3, 462) but we use Ali & Griem
        !  anyway (~ Barklem et al for 2p).
        !  A factor 1.13 is needed to match this region to the far red wing
        !  (below) at -4000 cm^-1 (62 A) detuning. A linear transition is used.
        xfac = 1.0
        if (q < 0) xfac = xfac + 0.13*q/q4000
        a = a + xfac*c%resont + c%vdw
      endif
      ! the wing damping depends on q, so its coefficients are built per point
      vc = qlg_bot_coef(a * c%idf)
    else
      vc = c%vcore
    endif

    phi = 0.d0

    !  Doppler core including fine structure
    do i=1, h%io
      vs = (q - h%dq(i)) * c%idf
      !  if INSIDE HFWID WE USE VOIGT FUNCTION WITH DOP, STARK AND LOR
      !  if OUTSIDE JUST ADD THE DOPPLER PART - REST IS ADDED LATER IN DETAIL
      xv = qlg_bot_eval(vc, vs) * c%idf
      phi = phi + h%wt(i)*xv
    end do

    if (abs(q) > c%hfwid) then
      if (q > redthresh) then
        !  Far red wing of Lyman Alpha in satellite region:
        !  Self-broadening following Kurucz's implementation of Allard et al
        !  (A&A 335, 1124, 1998), see Castelli & Kurucz (A&A 372, 260 2001) for
        !  more details.
        !  Tables LYMANH2() are on a delta wavenumber grid spacing of 200
        !  waves/cm -22000, -21800, ...., -4000 waves/cm (91 points). The
        !  profiles store log10(I(d omega)) for N(H) = 1e14 cm^-3. Assumed
        !  insensitive to temp, and linear scaling with N(H).  Note I(d nu) =
        !  I(d omega)/c
        if (q <= q4000) then
          !  If redward of last point (1660 -> 5000 Angstrom) extrapolate,
          !  otherwise (1278 -> 1660 Angstrom) linear (in log10 profiles)
          !  interpolation
          xlymanh2 = lyman_interp(lymanh2, 200., -22000., q*LyA0)
          phi = phi + LyA0 * c%tothi * 10.**xlymanh2

          !  Stark section

          !  Stark wings due to protons in the satellite region:
          !  Following Kurucz's implementation of Allard et al (A&A 335, 1124,
          !  1998), see Castelli & Kurucz (A&A 372, 260 2001) for more details.
          !  Tables LYMANH2PLUS() are on a delta wavenumber grid spacing of 100
          !  waves/cm -15000, -14900, ...., -4000 waves/cm (111 points).
          !  The profiles store log10(I(d omega)) for N(H+) = 1e14 cm^-3. Assumed
          !  insensitive to temp, and linear scaling with N(H+). Note I(d nu) =
          !  I(d omega)/c and we assume N(H+) = N(e-).

          !  If redward of last point (1487 -> 5000 Angstroem) extrapolate,
          !  otherwise (1278 -> 1487 Angstroem) interpolation.
          xlymanh2 = lyman_interp(lymanh2plus, 100., -15000., q*LyA0)
          phi = phi + LyA0 * c%ne * 10.**xlymanh2

          !  There may be a contribution due to electrons also at long range.
          !  The static Holtsmark profile (1/2 the total STARK1) is used noting
          !  that the correction for quantum effects suggested by Stehle (1994,
          !  A&AS 104, 509 eqn 7) has been applied already.
          phi = phi + 0.5*stark_eval(c%st, q)

        else
          !  Stark wings due to ions and electrons.  If Lyman alpha this is for
          !  before (blueward of) the satellite region.
          xstark = 0.5*stark_eval(c%st, q)

          !  If Lyman alpha we match the static ion part to the Allard et al
          !  value at -4000 cm^-1 detuning on red wing.
          if (q < 0 .and. xstark > 0) then
            sat4000 = LyA0 * c%ne * 10.**lymanh2plus(111)
            xfac = 1. + (sat4000/c%xholt4000 - 1.)*q/q4000
            xstark = xstark * (1. + xfac)
          end if
          phi = phi + xstark
        end if
      endif
      !  End Stark section
    end if
  END SUBROUTINE

  !====================================================================================
  ! Profile parameters of the hydrogen line n -> m with rest frequency nu0 [Hz]:
  ! Stark, natural, resonance and van der Waals widths, and the fine-structure
  ! pattern. Depends only on its arguments and the disable_*/hydrogen_BPO switches,
  ! so callers may rebuild it wherever it is needed instead of storing it.
  !------------------------------------------------------------------------------------
  FUNCTION hydrogen_line (n, m, nu0) result(h)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: nu0
    type(hline_t) :: h
    integer :: mmn, n2, m2, ipos, i, k, ifins
    real(4) :: mn2, ignm
    real(4) :: xknm, radamp, resont
    real(8), dimension(size(h%dq)) :: finest, finswt
    !----------------------------------------------------------------------------------
    h%n = n
    h%m = m
    h%nu0 = nu0
    mmn = m-n
    n2 = n*n
    m2 = m*m
    mn2 = m2*n2
    ignm = mn2/real(m2-n2)

    if ((mmn <= 3) .and. (n <= 4)) then
      xknm = xknmtb(n,mmn)
    else
      xknm = 5.5E-5*ignm * mn2/(1. + .13/real(mmn))
    end if
    h%gstk = 5.*xknm*h%nu0/(cc*1e8)

    !  Natural damping taken from tables
    !  dnu/ nu_0 = 1/2 * 1/nu_0 * 1/2 pi * Gamma
    if (n == 1) then
      radamp = asumlyman(m)
    else
      radamp = asum(n)+asum(m)
    end if
    h%grad = radamp * over_4pi / h%nu0

    !  Resonance broadening following Ali & Griem (1966, Phys Rev 144, 366).
    !  RESONT is dnu/nu per unit H density. Only the lower state is included
    !  (p-d approx for Balmer lines).  For N > 3 p-states play a very small
    !  role, and van der Waals might even dominate, there is however no
    !  available theory for this.  The lower state resonance broadening is
    !  used as a guess.
    if (n /= 1) then
      resont = hfnm(1,n)/(1. - 1./n2)
    else
      resont = hfnm(1,m)/(1. - 1./m2)
    end if
    h%gw1 = resont * 2.07E-24*ignm
    h%gw2 = 0

    if (hydrogen_BPO .and. (n == 2) .and. (m <= 5)) then
      !  For low Balmer lines use Barklem, Piskunov & O'Mara (2000, A&A 363, 1091)
      !  instead of Ali & Griem for self-broadening.
      resont = balfac(mmn)*sigma(mmn) ! already contains factors 1e4 and 1e6
      h%gw2 = 0.5*(1. - alf(mmn)) ! 0.5 is needed because we will be using VBARH squared

      ! The factor 1e-12 reconciles two differences from the original hbop.f:
      ! 1. The velocity constant here is 10,000x larger (4.2e8 vs 42008).
      ! 2. The original code normalized velocity by 1e4 (1e8 when squared).
      ! To compensate, we need (1/1e4) * (1/1e8) = 1e-12.
      h%gw1 = resont * 1e-12**h%gw2 / (h%nu0 * 2*PI)
      h%bpo = .true.
    endif

    h%gw3 = 4.45E-26*ignm*(m2*(7.*m2 + 5.))**0.4

    !  Approximately include fine structure.  Exact pattern for alpha lines,
    !  M infinite pattern used for others.
    if (n > 4 .or. m > 10) then
      ifins = 1
      finest(1) = 0.
      finswt(1) = 1.
    else if (m-n > 1) then
      ifins = lncomp(n)
      do i=1, ifins
        finest(i) = stcomp(i, n)
        finswt(i) = stcpwt(i, n)/(n*n)
      end do
    else
      !  eg: Ly alpha IFINS=2, IPOS=1, FINEST=-7.3E9,3.7E9, FINSWT=1/3, 2/3
      ifins = lnghal(n)
      ipos = istal(n)
      do i=1, ifins
        k = ipos-1+i
        finest(i) = stalph(k)
        finswt(i) = stwtal(k)/(3*n*n)
      end do
    end if

    h%io = ifins
    h%dq(1:ifins) = finest(1:ifins) / h%nu0
    h%wt(1:ifins) = finswt(1:ifins)
    h%dqbar = sum(h%wt(1:ifins)*h%dq(1:ifins))
    h%dqmax = maxval(abs(h%dq(1:ifins)))

    if (disable_rad) h%grad = 0
    if (disable_stk) h%gstk = 0
    if (disable_vdw) then
      h%gw1 = 0
      h%gw2 = 0
      h%gw3 = 0
      h%bpo = .false.
    endif
  END FUNCTION hydrogen_line

  !====================================================================================
  ! Fill a model-atom hydrogen transition (bbtrans_type, n and m in self%i, self%j)
  ! with the hydrogen_line parameters, and load the VCS tables if they are enabled.
  !------------------------------------------------------------------------------------
  SUBROUTINE prepare_hydrogen_damping (self)
    class(bbtrans_type), intent(inout) :: self
    type(hline_t) :: h
    !----------------------------------------------------------------------------------
    call ensure_vcs_loaded()
    h = hydrogen_line(self%i, self%j, self%nu0)
    self%grad  = h%grad
    self%gstk  = h%gstk
    self%gw1 = h%gw1
    self%gw2 = h%gw2
    self%gw3 = h%gw3
  END SUBROUTINE prepare_hydrogen_damping

  !===============================================================================
  FUNCTION lyman_interp(arr, dx, x0, newx) result(out)
    real(8), intent(in), dimension(:) :: arr
    real(4), intent(in) :: dx, x0
    real(8), intent(in) :: newx
    real(8) :: out, x
    integer :: i
    !---------------------------------------------------------------------------
    x = (newx - x0)/dx
    i = max(0.0, x)
    out = (arr(i+2) - arr(i+1))*(x - dble(i)) + arr(i+1)
  END FUNCTION lyman_interp

  !===============================================================================
  !  HFNM calculates hydrogen oscillator strengths
  !===============================================================================
  FUNCTION hfnm(n,m)
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(8)             :: hfnm
    real(8) :: wt, xmn12, fk, xmn
    real(8), save :: fnm, wtc, fkn, gca, ginf
    integer, save :: xn=0, xm=0
    !$omp threadprivate(fnm, wtc, fkn, gca, ginf, xn, xm)
    !---------------------------------------------------------------------------
    hfnm=0.
    if (m <= n) return
    if (n /= xn) then
      xn=n
      xm=m
      ginf=0.2027/xn**0.71
      gca=0.124/xn
      fkn=xn*1.9603
      wtc=0.45-2.4/xn**3*(xn-1.)
      xmn=m-n
      fk=fkn*(xm/(xmn*(xm+xn)))**3
      xmn12=xmn**1.2
      wt=(xmn12-1.)/(xmn12+wtc)
      fnm=fk*(1.-wt*ginf-(0.222+gca/xm)*(1.-wt))
    else if (m /= xm) then
      xm=m
      xmn=m-n
      fk=fkn*(xm/(xmn*(xm+xn)))**3
      xmn12=xmn**1.2
      wt=(xmn12-1.)/(xmn12+wtc)
      fnm=fk*(1.-wt*ginf-(0.222+gca/xm)*(1.-wt))
    end if
    hfnm=fnm
  END FUNCTION hfnm


  !===============================================================================
  !  Returns the Stark broadened line profile.  The code follows Griem's
  !  theories (mostly 1960 ApJ, 132, 883) with corrections to approximate
  !  the Vidal, Cooper & Smith (1973, ApJS 25, 37) profiles.

  !  Area normalised to unity with frequency.

  !  by Deane Peterson & Bob Kurucz.
  !  (adapted, corrected and comments added by PB)
  !===============================================================================
  !
  !  stark_prepare holds everything that does not depend on the detuning q, which is
  !  the bulk of the work; stark_eval is the profile at one q.
  !===============================================================================
  SUBROUTINE stark_prepare(n, m, nu0, temp, ne, hwstk, s)
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(8), intent(in) :: nu0
    real(4), intent(in) :: temp
    real(4), intent(in) :: ne
    real(4), intent(in) :: hwstk
    type(stark_cell_t), intent(out) :: s
    real(8) :: t4, xne16, y1b, y1s, gcon, xm2mn2, wty1, y1scal, c1
    real(8) :: y1num, y1wht
    integer :: mmn
    real(4), dimension(2,2), parameter :: y1wtm = reshape([1.e18, 1.e17, 1.e16, 1.e14], [2,2])
    !---------------------------------------------------------------------------
    s%off = hwstk == 0
    if (s%off) return
    s%n = n
    s%m = m

    ! hw5 = xknm * nu0 * fo / c (fundamental Stark shift scale)
    s%hw5 = hwstk * 0.2d0

    ! Normalized temperature and density
    t4 = temp * 1.d-4
    xne16 = ne ** 0.1666667d0

    ! Plasma parameters
    s%pp = xne16 * 0.08989d0 / sqrt(dble(temp))  ! shielding parameter
    s%cne = sqrt(5.96d-23 * ne)  ! characteristic electron detuning scale

    ! Y1 scaling factors
    y1b = 2.d0 / (1.d0 + 0.012d0 / temp * sqrt(ne / temp))
    y1s = t4 ** 0.3d0 / xne16

    ! Gamma correction factors
    s%gcon1 = 0.2d0 + 0.09d0 * sqrt(t4) / (1.d0 + ne * 1.d-13)
    s%gcon2 = 0.2d0 / (1.d0 + ne * 1.d-15)
    gcon = 1.d0 - s%gcon1 - s%gcon2

    ! Line-dependent variables
    mmn = m - n
    xm2mn2 = dble(m * m - n * n)

    ! Y1 weighting factors
    if (m == 2) then
      y1num = 550.d0
    else if (m == 3) then
      y1num = 380.d0
    else
      y1num = 320.d0
    end if

    if (mmn <= 2 .and. n <= 2) then
      y1wht = y1wtm(n, mmn)
    else if (mmn <= 3) then
      y1wht = 1.d14
    else
      y1wht = 1.d13
    end if

    ! Y1 scale factor
    wty1 = 1.d0 / (1.d0 + ne / y1wht)
    y1scal = y1num * y1s * wty1 + y1b * (1.d0 - wty1)

    s%ex = h_k * nu0 / temp
    s%c1coef = 0.5d0 * y1scal * s%ex * xm2mn2

    ! c1 and g1 (depend on hw5)
    c1 = s%hw5 * s%c1coef
    s%g1 = 6.77d0 * sqrt(c1)

    ! electron impact width in the small-detuning limit
    s%gam0 = s%g1 * max(0.d0, 0.2114d0 - log(c1 * s%cne)) * gcon
  END SUBROUTINE stark_prepare

  !===============================================================================
  FUNCTION stark_eval(s, q) result(stark)
    type(stark_cell_t), intent(in) :: s
    real(8), intent(in) :: q
    real(8) :: stark
    real(8) :: f, p1, q1, y1, y2, beta, gam, prqs, fns
    !---------------------------------------------------------------------------
    if (s%off) then
      stark = 0
      return
    endif

    ! Absolute relative detuning
    q1 = abs(q)

    y1 = s%c1coef * q1
    y2 = (q1 / s%cne) ** 2

    ! Reduced detuning
    beta = q1 / s%hw5

    ! Compute gamma (electron impact width)
    if ((y2 <= 1.d-4) .and. (y1 <= 1.d-5)) then
      gam = s%gam0
    else
      gam = s%g1 * (0.5d0 * exp(-min(80.d0, y1)) + vcse1f(y1) - 0.5d0 * vcse1f(y2)) &
        * (1.d0 - s%gcon1 / (1.d0 + (90.d0 * y1) ** 3) - s%gcon2 / (1.d0 + 2000.d0 * y1))
      if (gam <= 1.d-20) gam = 0.d0
    end if

    ! Quasistatic ion profile
    prqs = holtsmark(beta, s%pp, s%n, s%m)

    ! Impact (Lorentzian) electron profile
    if (gam > 0.d0) then
      f = gam / pi / (gam * gam + beta * beta)
    else
      f = 0.d0
    end if

    ! Quasistatic electron fraction
    p1 = (0.9d0 * y1) ** 2
    fns = (p1 + 0.03d0 * sqrt(y1)) / (p1 + 1.d0)

    ! Combined profile
    stark = (prqs * (1.d0 + fns) + f) / s%hw5

    ! Boltzmann factor for red wing (q < 0 means lower frequency)
    if (q < 0.d0) stark = stark * exp(s%ex * q)
  END FUNCTION stark_eval


  !===============================================================================
  !  Calculates S(BETA,P) for Hydrogen lines ie. the Holtsmark profile for
  !  quasistatic charged particles.  The alpha and beta lines of the first
  !  three series are explicitly included. All other cases use the H18
  !  profile. Profiles are normalised to full oscillator strength. Method
  !  is based on Griem (1960, ApJ 132, 883).

  !  By Deane Peterson and Bob Kurucz.

  !  STORAGE FOR CORRECTIONS (P,BETA,IND),(P,IND),(P,IND)
  !===============================================================================
  FUNCTION holtsmark(b,p,n,m)
    real(8), intent(in) :: b
    real(8), intent(in) :: p
    integer, intent(in) :: n
    integer, intent(in) :: m
    real                :: holtsmark
    real :: propbm(5,15,7),c(5,7),d(5,7), pp(5),beta(15)

    real(8) :: wtbm,wtbp,wtpp,sb,dd,cc,wt,pr2,cbm,cbp,wtpm,corr,&
      b2,pr1
    integer :: jm,jp,ip,im,mmn,indx

    real, dimension(75) :: prob1,prob2,prob3,prob4,prob5,prob6,prob7
    real, dimension(5)  :: c1,c2,c3,c4,c5,c6,c7,d1,d2,d3,d4,d5,d6,d7
    equivalence (propbm(1,1,1),prob1(1)),(propbm(1,1,2),prob2(1))
    equivalence (propbm(1,1,3),prob3(1)),(propbm(1,1,4),prob4(1))
    equivalence (propbm(1,1,5),prob5(1)),(propbm(1,1,6),prob6(1))
    equivalence (propbm(1,1,7),prob7(1))
    equivalence (c(1,1),c1(1)),(c(1,2),c2(1)),(c(1,3),c3(1))
    equivalence (c(1,4),c4(1)),(c(1,5),c5(1)),(c(1,6),c6(1))
    equivalence (c(1,7),c7(1))
    equivalence (d(1,1),d1(1)),(d(1,2),d2(1)),(d(1,3),d3(1))
    equivalence (d(1,4),d4(1)),(d(1,5),d5(1)),(d(1,6),d6(1))
    equivalence (d(1,7),d7(1))
    save propbm,c,d,pp,beta

    !  Lyman alpha

    DATA prob1/ -.980,-.967,-.948,-.918,-.873,-.968,-.949,-.921,-.879,-.821,  &
      -.950,-.922,-.883,-.830,-.764,-.922,-.881,-.830,-.770,-.706,  &
      -.877,-.823,-.763,-.706,-.660,-.806,-.741,-.682,-.640,-.625,  &
      -.691,-.628,-.588,-.577,-.599,-.511,-.482,-.484,-.514,-.568,  &
      -.265,-.318,-.382,-.455,-.531,-.013,-.167,-.292,-.394,-.478,  &
      .166,-.056,-.216,-.332,-.415, .251, .035,-.122,-.237,-.320,  &
      .221, .059,-.068,-.168,-.247, .160, .055,-.037,-.118,-.189,  &
      .110, .043,-.022,-.085,-.147/
    DATA c1 /-18.396, 84.674,-96.273,  3.927, 55.191/
    DATA d1 / 11.801,  9.079, -0.651,-11.071,-26.545/

    !  Lyman beta

    DATA prob2/ -.242, .060, .379, .671, .894, .022, .314, .569, .746, .818,  &
      .273, .473, .605, .651, .607, .432, .484, .489, .442, .343,  &
      .434, .366, .294, .204, .091, .304, .184, .079,-.025,-.135,  &
      .167, .035,-.082,-.189,-.290, .085,-.061,-.183,-.287,-.374,  &
      .032,-.127,-.249,-.344,-.418,-.024,-.167,-.275,-.357,-.420,  &
      -.061,-.170,-.257,-.327,-.384,-.047,-.124,-.192,-.252,-.306,  &
      -.043,-.092,-.142,-.190,-.238,-.038,-.070,-.107,-.146,-.187,  &
      -.030,-.049,-.075,-.106,-.140/
    DATA c2 / 95.740, 18.489, 14.902, 24.466, 42.456/
    DATA d2 / -6.665, -7.136,-10.605,-15.882,-23.632/

    !  Balmer alpha

    DATA prob3/ -.484,-.336,-.206,-.111,-.058,-.364,-.264,-.192,-.154,-.144,  &
      -.299,-.268,-.250,-.244,-.246,-.319,-.333,-.337,-.336,-.337,  &
      -.397,-.414,-.415,-.413,-.420,-.456,-.455,-.451,-.456,-.478,  &
      -.446,-.441,-.446,-.469,-.512,-.358,-.381,-.415,-.463,-.522,  &
      -.214,-.288,-.360,-.432,-.503,-.063,-.196,-.304,-.394,-.468,  &
      .063,-.108,-.237,-.334,-.409, .151,-.019,-.148,-.245,-.319,  &
      .149, .016,-.091,-.177,-.246, .115, .023,-.056,-.126,-.189,  &
      .078, .021,-.036,-.091,-.145/
    DATA c3 /-25.088,145.882,-50.165,  7.902, 51.003/
    DATA d3 /  7.872,  5.592, -2.716,-12.180,-25.661/

    !  Balmer beta

    DATA prob4/ -.082, .163, .417, .649, .829, .096, .316, .515, .660, .729,  &
      .242, .393, .505, .556, .534, .320, .373, .394, .369, .290,  &
      .308, .274, .226, .152, .048, .232, .141, .052,-.046,-.154,  &
      .148, .020,-.094,-.200,-.299, .083,-.070,-.195,-.299,-.385,  &
      .031,-.130,-.253,-.348,-.422,-.023,-.167,-.276,-.359,-.423,  &
      -.053,-.165,-.254,-.326,-.384,-.038,-.119,-.190,-.251,-.306,  &
      -.034,-.088,-.140,-.190,-.239,-.032,-.066,-.103,-.144,-.186,  &
      -.027,-.048,-.075,-.106,-.142/
    DATA c4 / 93.783, 10.066,  9.224, 20.685, 40.136/
    DATA d4 / -5.918, -6.501,-10.130,-15.588,-23.570/

    !  Paschen alpha

    DATA prob5/ -.819,-.759,-.689,-.612,-.529,-.770,-.707,-.638,-.567,-.498,  &
      -.721,-.659,-.595,-.537,-.488,-.671,-.617,-.566,-.524,-.497,  &
      -.622,-.582,-.547,-.523,-.516,-.570,-.545,-.526,-.521,-.537,  &
      -.503,-.495,-.496,-.514,-.551,-.397,-.418,-.448,-.492,-.547,  &
      -.246,-.315,-.384,-.453,-.522,-.080,-.210,-.316,-.406,-.481,  &
      .068,-.107,-.239,-.340,-.418, .177,-.006,-.143,-.246,-.324,  &
      .184, .035,-.082,-.174,-.249, .146, .042,-.046,-.123,-.190,  &
      .103, .036,-.027,-.088,-.146/
    DATA c5 /-19.819, 94.981,-79.606,  3.159, 52.106/
    DATA d5 / 10.938,  8.028, -1.267,-11.375,-26.047/

    !  Paschen beta

    DATA prob6/ -.073, .169, .415, .636, .809, .102, .311, .499, .639, .710,  &
      .232, .372, .479, .531, .514, .294, .349, .374, .354, .279,  &
      .278, .253, .212, .142, .040, .215, .130, .044,-.051,-.158,  &
      .141, .015,-.097,-.202,-.300, .080,-.072,-.196,-.299,-.385,  &
      .029,-.130,-.252,-.347,-.421,-.022,-.166,-.275,-.359,-.423,  &
      -.050,-.164,-.253,-.325,-.384,-.035,-.118,-.189,-.252,-.306,  &
      -.032,-.087,-.139,-.190,-.240,-.029,-.064,-.102,-.143,-.185,  &
      -.025,-.046,-.074,-.106,-.142/
    DATA c6 /111.107, 11.910,  9.857, 21.371, 41.006/
    DATA d6 / -5.899, -6.381,-10.044,-15.574,-23.644/

    !  Balmer 18

    DATA prob7/ .005, .128, .260, .389, .504, .004, .109, .220, .318, .389,  &
      -.007, .079, .162, .222, .244,-.018, .041, .089, .106, .080,  &
      -.026,-.003, .003,-.023,-.086,-.025,-.048,-.087,-.148,-.234,  &
      -.008,-.085,-.165,-.251,-.343, .018,-.111,-.223,-.321,-.407,  &
      .032,-.130,-.255,-.354,-.431, .014,-.148,-.269,-.359,-.427,  &
      -.005,-.140,-.243,-.323,-.386, .005,-.095,-.178,-.248,-.307,  &
      -.002,-.068,-.129,-.187,-.241,-.007,-.049,-.094,-.139,-.186,  &
      -.010,-.036,-.067,-.103,-.143/
    DATA c7 /511.318,  1.532,  4.044, 19.266, 41.812/
    DATA d7 / -6.070, -4.528, -8.759,-14.984,-23.956/
    DATA pp/0.,.2,.4,.6,.8/
    DATA beta/1.,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.310,  &
      7.943,10.,12.59,15.85,19.95,25.12/
    !---------------------------------------------------------------------------
    b2=b*b
    sb=sqrt(b)

    if(b > 500.) then
      !  Very large B
      holtsmark=(1.5/sb + 27./b2)/b2
      return
    end if

    !  Other cases
    corr=1.
    indx=7
    mmn=m-n
    if(n <= 3.and.mmn <= 2) indx=2*(n-1)+mmn

    !  Determine relevant Debye range
    im=min(int(5.*p)+1,4)
    ip=im+1
    wtpp=5.*(p-pp(im))
    wtpm=1.-wtpp
    if(b <= 25.12) then
      jp=2
      do while (b > beta(jp) .and. jp < 15)
        jp = jp + 1
      end do
      jm=jp-1

      wtbp=(b-beta(jm))/(beta(jp)-beta(jm))
      wtbm=1.-wtbp
      cbp=propbm(ip,jp,indx)*wtpp+propbm(im,jp,indx)*wtpm
      cbm=propbm(ip,jm,indx)*wtpp+propbm(im,jm,indx)*wtpm
      corr=1.+cbp*wtbp+cbm*wtbm

      !  Get approximate profile for the inner part
      pr1=0.
      pr2=0.
      wt=max(min(0.5*(10.-b),1.),0.)
      if(b <= 10.) pr1=8./(83.+(2.+0.95*b2)*b)
      if(b >= 8.)  pr2=(1.5/sb+27./b2)/b2
      holtsmark=(pr1*wt+pr2*(1.-wt))*corr

    else
      !  Asymptotic part for medium B's
      cc=c(ip,indx)*wtpp+c(im,indx)*wtpm
      dd=d(ip,indx)*wtpp+d(im,indx)*wtpm
      corr=1.+dd/(cc+b*sb)
      holtsmark=(1.5/sb+27./b2)/b2*corr
    end if
  END FUNCTION holtsmark

  !===============================================================================
  !  E1 function calculator for VCS approximation. It's rough, but
  !  arranged to be fast. X must be >=0.
  !  From Kurucz codes.
  !===============================================================================
  FUNCTION vcse1f(x)
    real(8), intent(in)  :: x
    real(8) :: vcse1f
    !---------------------------------------------------------------------------
    vcse1f = 0.0
    if(x <= 0.0) return
    if(x <= 0.01) then
      vcse1f = -log(x) - euler_const + x
    elseif (x <= 1.0) then
      vcse1f = -log(x) - euler_const + x*(0.99999193+x*(-0.24991055+  &
        x*(0.05519968+x*(-0.00976004+ x*0.00107857))))
    elseif (x <= 30.) then
      vcse1f = (x*(x+2.334733)+0.25062)/(x*(x+3.330657)+ 1.681534)/x*exp(-x)
    end if
  END FUNCTION vcse1f

END MODULE
