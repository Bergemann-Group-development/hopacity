!===============================================================================
!> Bound-bound transitions: the transition type, its Einstein coefficients and the
!> line broadening constants shared by model-atom and line-list lines.
!===============================================================================
MODULE bbtrans_mod
  USE prec_mod
  USE const_mod
  implicit none

  !===============================================================================
  !> A bound-bound transition of a model atom: its rest frequency and wavelengths,
  !> Einstein coefficients and line broadening constants.
  !===============================================================================
  type :: bbtrans_type
    real(dp) :: nu0
    real(dp) :: lam0
    real(dp) :: lam_air
    integer :: nnu                                          !< number of frequencies in the line
    real(dp) :: grad=0, gw1=0, gw2=0, gw3=0, gstk=0, gvdw=0 !< line broadening (damping) constants
    real(dp) :: grad_in=0, gstk_in=0, gv_in=0               !< input values for line broadening
    real(dp) :: lower_ev, upper_g=1, loggf
    real(dp) :: Aji, Aij, Bji, Bij !< Einstein coefficients, ij from the lower to the upper level
    logical :: active=.true., unsold
    integer :: i=0 !< lower level
    integer :: j=0 !< upper level
  contains
    procedure :: init_einstein_coefficients
    procedure :: damping => line_damping
    procedure :: prepare_Radiative => bb_prepare_Radiative
    procedure :: prepare_Stark => bb_prepare_Stark
    procedure :: prepare_VanderWaals => bb_prepare_VanderWaals
    procedure :: deactivate
  end type

  logical, public :: disable_stk=.false. !< no Stark broadening
  logical, public :: disable_vdw=.false. !< no van der Waals broadening
  logical, public :: disable_rad=.false. !< no radiative damping
  !> Molecular lines may take damping estimates from the molecule's ionisation energy:
  !> the Cowley Stark estimate (prepare_stark) for a line with no tabulated Stark width
  !> and the Unsold van der Waals estimate (prepare_vanderwaals) for a line whose gvdw
  !> is an enhancement factor. The same switch applies the tabulated negative gvdw forms
  !> (log C6, Turbospectrum 'W') to molecular lines, as Korg does. MULTI and Turbospectrum
  !> apply none of these to molecules, and their presets switch them off.
  logical, public :: molecule_damping_estimates=.true.
  real(4), public, save :: tolerance=1e-4 !< line opacity below tolerance times the continuum is left out
  !> lines below max(tolerance, line_reject_frac)*tolerance times the line opacity
  !> already accumulated are left out
  real(4), public, save :: line_reject_frac=0.0
  real(8), parameter :: f_conv = 8.0d0 * pi**2 * ec**2 / em / cc**3 !< 8 pi^2 e^2/(m_e c^3); f_conv gf nu0^2 = g_u A_ul [s^-1]

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
  !> Damping parameter adamp = Gamma/(4 pi nu0) of the line in each cell, the sum of its
  !> van der Waals, Stark and radiative damping.
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
    if (self%gstk >= 0.0) then
      adamp = adamp + self%gstk * ne
    else! if line%gstk < 0
      adamp = adamp - self%gstk * stark_pe
    end if

    adamp = adamp + self%grad
  END SUBROUTINE line_damping
  !------------------------------------------------------------------------------------

  !====================================================================================
  !> Sets the line's Stark constant gstk with prepare_stark.
  !====================================================================================
  PURE ELEMENTAL SUBROUTINE bb_prepare_Stark(self, gstk, j, cnt_ev, is_mol)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: gstk, cnt_ev
    integer, intent(in) :: j
    logical, intent(in) :: is_mol
    !----------------------------------------------------------------------------------
    call prepare_stark(self%gstk, gstk, self%nu0, self%lower_ev, j, cnt_ev, is_mol)
  END SUBROUTINE bb_prepare_Stark

  !====================================================================================
  !> Sets the line's radiative damping grad with prepare_radiative.
  !====================================================================================
  PURE ELEMENTAL SUBROUTINE bb_prepare_Radiative(self, grad)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: grad
    !----------------------------------------------------------------------------------
    call prepare_radiative(self%grad, grad, self%nu0, self%loggf)
  END SUBROUTINE bb_prepare_Radiative

  !====================================================================================
  !> Sets the line's van der Waals constants gw1, gw2 and unsold with prepare_vanderwaals.
  !====================================================================================
  PURE SUBROUTINE bb_prepare_VanderWaals(self, gvdw, lower_ion, cnt_ev, vbar2e4, is_mol)
    class(bbtrans_type), intent(inout)  :: self
    integer, intent(in) :: lower_ion
    real(dp), intent(in) :: gvdw, cnt_ev, vbar2e4
    logical, intent(in) :: is_mol
    !----------------------------------------------------------------------------------
    call prepare_vanderwaals(self%gw1, self%gw2, self%unsold, gvdw, self%nu0, &
                             self%lower_ev, lower_ion, cnt_ev, vbar2e4, is_mol)
  END SUBROUTINE bb_prepare_VanderWaals

  !====================================================================================
  !> Radiative damping of a line, scaled by 1/(4 pi nu0): grad_in [s^-1], or for
  !> grad_in = 0 the value g_u A_ul = 8 pi^2 e^2 nu0^2 gf/(m_e c^3) that Korg uses.
  !====================================================================================
  PURE SUBROUTINE prepare_radiative(grad, grad_in, nu0, loggf)
    real(dp), intent(out) :: grad
    real(dp), intent(in)  :: grad_in, nu0, loggf
    !----------------------------------------------------------------------------------
    if (grad_in == 0) then
      ! Default as in Korg: 8 pi^2 e^2 nu0^2 gf / (m_e c^3) = g_u A_ul. Turbospectrum
      ! divides this by g_u, giving A_ul; MOOG uses the classical damping constant
      ! 8 pi^2 e^2 nu0^2 / (3 m_e c^3), which equals A_ul only for gf/g_u = 1/3.
      ! Aji cannot be reused here, because it may contain isotope fractions.
      grad = f_conv * 10**loggf * nu0**2
    else
      grad = grad_in                      ! radiative damping rate, 1/s
    end if
    grad = grad * over_4pi / nu0
    if (disable_rad) grad = 0
  END SUBROUTINE prepare_radiative

  !====================================================================================
  !> Stark broadening constant of a line, scaled by 1/(4 pi nu0).  gstk_in selects
  !> the recipe, and a negative result marks the T-dependent ones for line_damping:
  !>
  !> | gstk_in  | meaning                                                                         |
  !> |----------|---------------------------------------------------------------------------------|
  !> | < -10    | log10 C4 as in SIU, converted with Gray (2005), eqs. 11.26-11.27                |
  !> | -10 to 0 | log10 gamma_4/n_e at 10^4 K, scaled as T^(1/6)                                  |
  !> | 0        | Cowley (1971) estimate at 10^4 K, none if the upper level is at or above cnt_ev |
  !> | > 0      | gamma_4/n_e as given                                                            |
  !>
  !> [Gray (2005)](https://ui.adsabs.harvard.edu/abs/2005oasp.book.....G),
  !> [Cowley (1971)](https://ui.adsabs.harvard.edu/abs/1971Obs....91..139C)
  !====================================================================================
  PURE SUBROUTINE prepare_stark(gstk, gstk_in, nu0, lower_ev, j, cnt_ev, is_mol)
    real(dp), intent(out) :: gstk
    real(dp), intent(in)  :: gstk_in, nu0, lower_ev, cnt_ev
    integer,  intent(in)  :: j
    logical,  intent(in)  :: is_mol
    real(dp), parameter :: scale = 1e4**(-one_sixth)   ! (T/1e4 K)**(1/6), with T**(1/6) in stark_pe
    real(dp), parameter :: c19 = log10(39. * (8. / (PI*EM*BK**5))**one_sixth) ! ~19.4
    real(dp) :: upper_ev, c1
    !----------------------------------------------------------------------------------
    gstk = 0
    if (gstk_in < -10) then
      gstk = - BK * 10**(c19 + two_third*gstk_in)   ! recipe c4 SIU (Gray 2005, eqs. 11.26-11.27); sign marks recipe
    elseif (gstk_in < 0) then
      gstk = - 10**gstk_in * scale                  ! sign marks recipe
    elseif (gstk_in == 0) then
      upper_ev = lower_ev + H_E * nu0
      ! The estimate is hydrogenic in the upper level's distance below the ionisation
      ! limit cnt_ev. It has no meaning for an upper level at or above the limit
      ! (autoionising, converging on an excited parent term, or listed against a
      ! rounded limit), where n*^4 would be taken from |cnt_ev - upper_ev| and grows
      ! without bound towards it, so such a line gets no Cowley estimate. A molecular
      ! line gets none either unless molecule_damping_estimates is set.
      if (upper_ev < cnt_ev .and. (molecule_damping_estimates .or. .not. is_mol)) then
        ! Cowley (1971) eqs. 5 and 6 at 1e4 K, which he recommends using at all T;
        ! c1 is n*^4 of the upper level (his eq. 4)
        c1 = (j**2 * RydV / (cnt_ev - upper_ev))**2
        if (j == 1) then
          gstk = 2.26e-7 * c1
        else
          gstk = 5.42e-7 * c1 / (j+1)**2
        end if
      end if
    else
      gstk = gstk_in
    end if
    gstk = gstk * over_4pi / nu0
    if (disable_stk) gstk = 0
  END SUBROUTINE prepare_stark

  !====================================================================================
  !> Van der Waals broadening constants of a line: gw1, scaled by 1/(4 pi nu0), and
  !> the velocity exponent gw2; unsold marks the recipes that line_damping evaluates
  !> with the Unsold temperature dependence.  gvdw selects the recipe:
  !>
  !> | gvdw     | meaning                                                                              |
  !> |----------|--------------------------------------------------------------------------------------|
  !> | >= 20    | ABO cross section [a0^2] at 10^4 m/s (integer part) and velocity exponent (fraction) |
  !> | < -20    | log10 C6 as in SIU                                                                   |
  !> | -20 to 0 | log10 gamma_6/n_H at 10^4 K (Turbospectrum 'W')                                      |
  !> | 0        | none                                                                                 |
  !> | 0 to 20  | Unsold estimate times gvdw, none if the upper level is at or above cnt_ev            |
  !>
  !> Molecular lines take the recipes below 20 only with molecule_damping_estimates.
  !====================================================================================
  PURE SUBROUTINE prepare_vanderwaals(gw1, gw2, unsold, gvdw, nu0, lower_ev, lower_ion, &
                                      cnt_ev, vbar2e4, is_mol)
    real(dp), intent(out) :: gw1, gw2
    logical,  intent(out) :: unsold
    real(dp), intent(in)  :: gvdw, nu0, lower_ev, cnt_ev, vbar2e4
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
    if (gvdw >= 20.0) then
      ! van der Waals broadening by hydrogen: Barklem cross section + velocity exponent
      sigma  = barklem * int(gvdw)
      alpha  = gvdw - int(gvdw)
      gx     = 1.0 - 0.5*alpha
      gammaf = 1 + (-.5748646 + (.9512363 + (-.6998588 + (.4245549 - .1010678*gx)*gx)*gx)*gx)*gx
      balfac = (4.0/pi)**(0.5*alpha) * gammaf * 2.e10
      gw2 = 0.5*(1 - alpha)
      gw1 = balfac * sigma * 10**(-12*gw2)
    elseif (.not. is_mol .or. molecule_damping_estimates) then
      ! Atomic lines take every recipe below; molecular lines only with
      ! molecule_damping_estimates, the Unsold estimate from the molecule's ionisation
      ! energy.
      unsold = .true.
      if (gvdw < -20.0) then
        c625 = 10**(0.4 * gvdw)       ! recipe c6 SIU
        gw1  = 17. * c625
      elseif (gvdw < 0.0) then
        gw1 = 10**gvdw / (vbar2e4)**0.3  ! recipe 'W' Turbospectrum
      elseif (gvdw /= 0.0) then
        upper_ev = lower_ev + H_E * nu0
        ! Unsold approximation, hydrogenic in the levels' distances below the
        ! ionisation limit cnt_ev. As for the Cowley Stark estimate (prepare_stark), a
        ! line whose upper level is at or above the limit gets none.
        if (upper_ev < cnt_ev) then
          c1 = 1./(cnt_ev - upper_ev)**2 - 1./(cnt_ev - lower_ev)**2
          if (c1 > 0) then
            c625 = (lower_ion**2 * c1)**0.4  ! c6^(2/5) impact approximation
            gw1  = gvdw * uns * c625
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
  !> Marks the line inactive.
  !====================================================================================
  PURE SUBROUTINE deactivate(self)
    class(bbtrans_type), intent(inout) :: self
    !----------------------------------------------------------------------------------
    self%active = .false.
  END SUBROUTINE deactivate

END MODULE bbtrans_mod
