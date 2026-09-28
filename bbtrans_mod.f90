MODULE bbtrans_mod
  USE prec_mod
  USE const_mod
  implicit none

  type :: bbtrans_type
    real(dp) :: nu0
    real(dp) :: lam0
    real(dp) :: lam_air
    integer :: nnu ! number of frequencies in the line
    real(dp) :: ga=0, gw1=0, gw2=0, gw3=0, gq=0, gv=0 ! line broadening (damping) constants
    real(dp) :: ga_in=0, gq_in=0, gv_in=0 ! input values for line broadening
    real(dp) :: lower_ev, upper_g=1, loggf
    real(dp) :: Aji, Aij, Bji, Bij ! Einstein coefficients, Bij is the
    ! transition from the lower to the upper
    ! level
    logical :: active=.true., unsold
    ! Added to deal with multi component lines
    ! IO > 1 leads to the following being read from the atom file
    integer :: io = 0 ! number of components in term (lifted from Multi 2.3)
    real(8), dimension(:), allocatable :: dq
    real(8), dimension(:), allocatable :: wt
    real(8), dimension(:), allocatable :: da
    integer :: i=0 ! lower level
    integer :: j=0 ! upper level
  contains
    procedure :: init_einstein_coefficients
    procedure :: damping => line_damping
    procedure :: prepare_Radiative
    procedure :: prepare_Stark
    procedure :: prepare_VanderWaals
    procedure :: deactivate
  end type

  logical, public :: disable_stark=.false., disable_vdw=.false., disable_ga=.false.
  logical, public :: ga_fallback=.false.
  real(4), public, save :: tolerance=1e-4
  real(4), public, save :: line_reject_frac=0.0
  real(4), parameter :: pqs  = 1.128379167  ! 2./pi**0.5
  real(8), parameter :: f_conv = 8.0d0 * pi**2 * ec**2 / em / cc**3

contains

  !====================================================================================
  SUBROUTINE init_einstein_coefficients(self, f, gij)
    class(bbtrans_type), intent(inout) :: self
    real(8), intent(in) :: f, gij
    !----------------------------------------------------------------------------------
    ! Einstein coefficients include an additional factor of h_4pi for convenience
    self%Aij = f_conv * f * self%nu0**2
    self%Aij = self%Aij * h_4pi
    self%Bij = c2_2h / self%nu0**3 * self%Aij
    self%Aji = self%Aij * gij
    self%Bji = self%Bij * gij
  END SUBROUTINE init_einstein_coefficients

  !====================================================================================
  ! Voigt profiles are evaluated in HOP/qlg_mod.f90 (qlg_bot_sweep for a line window,
  ! qlg_bot_phi / qlg_bot_profile for single points). pseudo_voigt, further down, is a
  ! fitted line shape for the convolution kernels.
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE line_damping (self, adamp, ne, tothi, vbar2, stark_pe, unsold)
    class(bbtrans_type), intent(in)  :: self
    real(4), intent(in)  :: ne, tothi, vbar2, stark_pe, unsold
    real(4), intent(inout) :: adamp
    !----------------------------------------------------------------------------------
    adamp = 0
    ! Van der Waals Broadening
    if (self%unsold) then
      adamp = adamp + self%gw1 * unsold
    else
      adamp = adamp + self%gw1 * vbar2**self%gw2 * tothi
    endif

    ! Stark Broadening
    if (self%gq >= 0.0) then
      adamp = adamp + self%gq * ne
    else! if line%gq < 0
      adamp = adamp - self%gq * stark_pe
    end if

    adamp = adamp + self%ga
  END SUBROUTINE line_damping
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE prepare_Stark(self, gq, j, cnt_ev)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: gq, cnt_ev
    integer, intent(in) :: j
    real(dp), parameter :: scale = 1e4**(-one_sixth)   ! (T/1e4 K)**(1/6), with T**(1/6) in stark_pe
    real(dp), parameter :: c19=log10(39.*(8./(PI*EM*BK**5))**one_sixth) ! 19.4
    real(dp) :: nu0_ev, upper_ev, c1
    !----------------------------------------------------------------------------------
    if (gq < -10) then
      ! recipe c4 SIU (Gray 2005, eqs. 11.26-11.27)
      ! make negative to mark this recipe
      self%gq = - BK * 10**(c19 + two_third*gq)
    elseif (gq < 0) then
      ! make negative to mark this recipe
      self%gq = - 10**gq * scale
    elseif (gq == 0) then
      nu0_ev = H_E * self%nu0
      upper_ev = self%lower_ev + nu0_ev
      if (self%lower_ev < cnt_ev) then
        ! Cowley (1971) eqs. 5 and 6 at 1e4 K, which he recommends using at all T;
        ! c1 is n*^4 of the upper level (his eq. 4)
        c1 = (j**2 * RydV / (cnt_ev - upper_ev))**2
        if (j == 1) then
          self%gq = 2.26e-7 * c1
        else
          self%gq = 5.42e-7 * c1 / (j+1)**2
        endif
      endif
    else
      self%gq = gq
    endif
    self%gq = self%gq * over_4pi / self%nu0
    if (disable_stark) self%gq = 0
  END SUBROUTINE prepare_Stark
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE prepare_Radiative(self, ga)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: ga
    !----------------------------------------------------------------------------------
    if (ga == 0) then
      ! Default as in Korg: 8 pi^2 e^2 nu0^2 gf / (m_e c^3) = g_u A_ul. Turbospectrum
      ! divides this by g_u, giving A_ul; MOOG uses the classical damping constant
      ! 8 pi^2 e^2 nu0^2 / (3 m_e c^3), which equals A_ul only for gf/g_u = 1/3.
      ! We cannot reuse Aji at this point, because it may contain isotope fractions.
      self%ga = f_conv * 10**self%loggf * self%nu0**2
    else
      self%ga = ga                    ! radiative damping rate, 1/s
    endif
    self%ga = self%ga * over_4pi / self%nu0
    if (disable_ga) self%ga = 0
  END SUBROUTINE prepare_Radiative
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE SUBROUTINE prepare_VanderWaals(self, gv, lower_ion, cnt_ev, vbar2e4, is_mol)
    class(bbtrans_type), intent(inout)  :: self
    integer, intent(in) :: lower_ion
    real(dp), intent(in) :: gv, cnt_ev, vbar2e4
    logical, intent(in) :: is_mol
    real(dp) :: upper_ev, sigma, alpha, gx, gammaf, balfac, c1, c625, vbar2
    ! uns = 17*C_U**0.4 (CGS, energies in eV). Unsold (1955, eqs. 82.48 and 82.55) gives
    ! gamma_6 = 17 C6**0.4 v**0.6 N with C6 = 1.61e-33*(13.5*Z)**2*[...], i.e. C_U = 2.93e-31
    ! and uns = 1.041e-11; Gray (2005, eq. 11.30) uses C_U = 3.0e-31. The factors below
    ! correspond to C_U = 3.22e-31, which makes these widths 3.7% larger than Unsold's.
    real(dp), parameter :: uns=8.411*1.283984e-12
    real(dp), parameter :: barklem = (1e-2*Rbohr)**2 ! 2.80028e-21 bohr radius in m and squared
    !----------------------------------------------------------------------------------
    ! compute van der waals broadening due to hydrogen
    ! following Barklem
    self%unsold = .false.
    if(gv >= 20.0) then
      sigma = barklem*int(gv)
      alpha = gv - int(gv)
      gx = 1.0 - 0.5*alpha
      gammaf = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*gx)*gx)*gx)*gx)*gx
      balfac = (4.0/pi)**(0.5*alpha) * gammaf * 2.e10

      self%gw2 = 0.5*(1 - alpha)
      self%gw1 = balfac * sigma * 10**(-12*self%gw2)
    elseif (.not. is_mol) then ! cannot use Unsold for molecules
      self%unsold = .true.
      if (gv < -20.0) then
        ! recipe c6 SIU
        c625 = 10**(0.4 * gv)
        self%gw1 = 17. * c625
      elseif (gv < 0.0) then
        ! recipe 'W' Turbospectrum
        self%gw1 = 10**gv / (vbar2e4)**0.3
      elseif (gv.ne.0.0) then
        upper_ev = self%lower_ev + H_E * self%nu0
        ! Unsold approximation
        if (self%lower_ev < cnt_ev) then
          c1 = 1./(cnt_ev - upper_ev)**2 - 1./(cnt_ev - self%lower_ev)**2
          if (c1 > 0) then
            c625 = (lower_ion**2 * c1)**0.4 ! c6^(2/5) impact approximation
            self%gw1 = gv * uns * c625
          else
            self%gw1 = 0
          endif
        endif
      endif
    endif
    self%gw1 = self%gw1 * over_4pi / self%nu0
    if (disable_vdw) then
      self%unsold = .true.
      self%gw1 = 0
    endif
  END SUBROUTINE prepare_VanderWaals
  !------------------------------------------------------------------------------------

  
  !====================================================================================
  PURE SUBROUTINE deactivate(self)
    class(bbtrans_type), intent(inout) :: self
    !----------------------------------------------------------------------------------
    self%active = .false.
  END SUBROUTINE deactivate

  !====================================================================================                                                                                     !!
  SUBROUTINE hum1zpf16  (nx, x,y, w)
    !  complex probability function for complex argument Z=X+iY
    !  real part = voigt function K(x,y)
    !
    !  large x+y     J  Humlicek, JQSRT 27, 437, 1982       asymptotic region I approximation
    !  small x+y:    J  Humlicek, JQSRT 21, 309-313, 1979   16-term rational approximation Eq. (6)
    !====================================================================================
    real(4), parameter ::  sqrtpi    = 1.7724538509055159d0 ! sqrt(pi)
    real(4), parameter ::  recsqrtpi = one/sqrtpi
    integer, intent(in) :: nx
    real(4), intent(in) :: x(nx), y
    complex(8), intent(out) ::  w(nx)
    ! humlicek prbfct region i bounds
    real(4), parameter :: s15=15.d0
    real(4) ::            x12, x21
    ! humlicek zpf16 constants
    complex(8), parameter ::  a(16) = &
      [(41445.0374210222, 0.0),   (0.0, -136631.072925829),  &
      (-191726.143960199, 0.0),  (0.0, 268628.568621291),   &
      (173247.907201704,  0.0),  (0.0, -179862.56759178),   &
      (-63310.0020563537, 0.0),  (0.0, 56893.7798630723),   &
      (11256.4939105413,  0.0),  (0.0, -9362.62673144278),  &
      (-1018.67334277366, 0.0),  (0.0, 810.629101627698),   &
      (44.5707404545965, 0.0),   (0.0, -34.5401929182016),  &
      (-0.740120821385939, 0.0), (0.0, 0.564189583547714)]
    real(4), parameter ::  b(8) =  &
      [ 7918.06640624997, -126689.0625, 295607.8125, &
      -236486.25, 84459.375, -15015.0, 1365.0, -60.0]
    integer :: i
    complex(8) :: t, z,zz, numer, denom
    !----------------------------------------------------------------------------------
    x12 = y - s15 !        left wing -- center
    x21 = -x12    ! 15-y   center -- right wing
    if (y>s15 .or. x(1)>x21 .or. x(nx)<x12) then
      ! all points are in humlicek's region i
      do i=1,nx
        t    = cmplx(y,-x(i))
        w(i) = (recsqrtpi*t) / (half + t*t)
      end do
    else
      do i=1,nx
        if (abs(x(i))>x21) then
          t    = cmplx(y,-x(i))
          w(i) = (recsqrtpi*t) / (half + t*t)
        else
          z     = cmplx(x(i), y+1.31183)
          zz    = z*z
          numer = ((((((a(16)*z+a(15))*z+a(14))*z+a(13))*z+a(12))*z+a(11))*z+a(10))*z+a(9)
          numer = ((((((((numer*z+a(8))*z+a(7))*z+a(6))*z+a(5))*z+a(4))*z+a(3))*z+a(2))*z+a(1))
          denom = b(1)+(b(2)+(b(3)+(b(4)+(b(5)+(b(6)+(b(7)+b(8)*zz)*zz)*zz)*zz)*zz)*zz)*zz
          w(i)  = numer/denom
        end if
      end do
    end if
    w = w * recsqrtpi
  END SUBROUTINE hum1zpf16

  !====================================================================================
  !DIR$ ATTRIBUTES FORCEINLINE :: pseudo_voigt
  PURE FUNCTION pseudo_voigt(x, peak, eta0, eta1, eta2) result(phi)
    real, intent(in) :: x, peak, eta0, eta1, eta2
    real :: x2, L, G, psd, phi, e
    integer :: ix
    real, parameter :: iln2 = 1/log(2.0)
    !----------------------------------------------------------------------------------
    x2 = x*x
    L = 1 / (1 + x2)
    psd = L*(eta1 + L*eta2)

    ! Reinterpret the bits of ref as integer
    ix = transfer(psd, ix)
    ! Shift right to get exponent bits and subtract bias (127)
    ! and subtract extra 10 to check if k is larger than e*2^-10
    e = 137 - ishft(ix, -23)

    if (x2*iln2 < e) then
      G = exp(-x2)
    else
      G = 0
    endif

    psd = psd + eta0*G
    phi = peak*psd
  END FUNCTION pseudo_voigt

END MODULE bbtrans_mod
