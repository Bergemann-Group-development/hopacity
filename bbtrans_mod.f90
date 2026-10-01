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
    procedure :: prepare_Radiative => bb_prepare_Radiative
    procedure :: prepare_Stark => bb_prepare_Stark
    procedure :: prepare_VanderWaals => bb_prepare_VanderWaals
    procedure :: deactivate
  end type

  logical, public :: disable_stark=.false., disable_vdw=.false., disable_ga=.false.
  real(4), public, save :: tolerance=1e-4
  real(4), public, save :: line_reject_frac=0.0
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
  ! qlg_bot_phi / qlg_bot_profile for single points).
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
  ! Line broadening constants for model-atom lines (the type-bound forms) and
  ! line-list lines (linelist_mod%build_linelist, which calls the scalar forms
  ! directly).  All results are scaled by over_4pi/nu0.
  !------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE bb_prepare_Stark(self, gq, j, cnt_ev)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: gq, cnt_ev
    integer, intent(in) :: j
    !----------------------------------------------------------------------------------
    call prepare_stark(self%gq, gq, self%nu0, self%lower_ev, j, cnt_ev)
  END SUBROUTINE bb_prepare_Stark

  PURE ELEMENTAL SUBROUTINE bb_prepare_Radiative(self, ga)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: ga
    !----------------------------------------------------------------------------------
    call prepare_radiative(self%ga, ga, self%nu0, self%loggf)
  END SUBROUTINE bb_prepare_Radiative

  PURE SUBROUTINE bb_prepare_VanderWaals(self, gv, lower_ion, cnt_ev, vbar2e4, is_mol)
    class(bbtrans_type), intent(inout)  :: self
    integer, intent(in) :: lower_ion
    real(dp), intent(in) :: gv, cnt_ev, vbar2e4
    logical, intent(in) :: is_mol
    !----------------------------------------------------------------------------------
    call prepare_vanderwaals(self%gw1, self%gw2, self%unsold, gv, self%nu0, &
                             self%lower_ev, lower_ion, cnt_ev, vbar2e4, is_mol)
  END SUBROUTINE bb_prepare_VanderWaals

  !====================================================================================
  PURE SUBROUTINE prepare_radiative(ga, ga_in, nu0, loggf)
    real(dp), intent(out) :: ga
    real(dp), intent(in)  :: ga_in, nu0, loggf
    !----------------------------------------------------------------------------------
    if (ga_in == 0) then
      ! Default as in Korg: 8 pi^2 e^2 nu0^2 gf / (m_e c^3) = g_u A_ul. Turbospectrum
      ! divides this by g_u, giving A_ul; MOOG uses the classical damping constant
      ! 8 pi^2 e^2 nu0^2 / (3 m_e c^3), which equals A_ul only for gf/g_u = 1/3.
      ! Aji cannot be reused here, because it may contain isotope fractions.
      ga = f_conv * 10**loggf * nu0**2
    else
      ga = ga_in                      ! radiative damping rate, 1/s
    end if
    ga = ga * over_4pi / nu0
    if (disable_ga) ga = 0
  END SUBROUTINE prepare_radiative

  !====================================================================================
  PURE SUBROUTINE prepare_stark(gq, gq_in, nu0, lower_ev, j, cnt_ev)
    real(dp), intent(out) :: gq
    real(dp), intent(in)  :: gq_in, nu0, lower_ev, cnt_ev
    integer,  intent(in)  :: j
    real(dp), parameter :: scale = 1e4**(-one_sixth)   ! (T/1e4 K)**(1/6), with T**(1/6) in stark_pe
    real(dp), parameter :: c19 = log10(39. * (8. / (PI*EM*BK**5))**one_sixth) ! ~19.4
    real(dp) :: upper_ev, c1
    !----------------------------------------------------------------------------------
    gq = 0
    if (gq_in < -10) then
      gq = - BK * 10**(c19 + two_third*gq_in)   ! recipe c4 SIU (Gray 2005, eqs. 11.26-11.27); sign marks recipe
    elseif (gq_in < 0) then
      gq = - 10**gq_in * scale                  ! sign marks recipe
    elseif (gq_in == 0) then
      upper_ev = lower_ev + H_E * nu0
      ! The estimate is hydrogenic in the upper level's distance below the ionisation
      ! limit cnt_ev. It has no meaning for an upper level at or above the limit
      ! (autoionising, converging on an excited parent term, or listed against a
      ! rounded limit), where n*^4 would be taken from |cnt_ev - upper_ev| and grows
      ! without bound towards it, so such a line gets no Cowley estimate.
      if (upper_ev < cnt_ev) then
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

  !====================================================================================
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
    real(dp), parameter :: barklem = (1e-2 * Rbohr)**2 ! 2.80028e-21 bohr radius in m and squared
    !----------------------------------------------------------------------------------
    gw1 = 0
    gw2 = 0
    unsold = .false.
    if (gv >= 20.0) then
      ! van der Waals broadening by hydrogen: Barklem cross section + velocity exponent
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
        ! Unsold approximation, hydrogenic in the levels' distances below the
        ! ionisation limit cnt_ev. As for the Cowley Stark estimate (prepare_stark), a
        ! line whose upper level is at or above the limit gets none.
        if (upper_ev < cnt_ev) then
          c1 = 1./(cnt_ev - upper_ev)**2 - 1./(cnt_ev - lower_ev)**2
          if (c1 > 0) then
            c625 = (lower_ion**2 * c1)**0.4  ! c6^(2/5) impact approximation
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
  !------------------------------------------------------------------------------------

  
  !====================================================================================
  PURE SUBROUTINE deactivate(self)
    class(bbtrans_type), intent(inout) :: self
    !----------------------------------------------------------------------------------
    self%active = .false.
  END SUBROUTINE deactivate

END MODULE bbtrans_mod
