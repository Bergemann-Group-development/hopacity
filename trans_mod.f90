MODULE trans_mod
  USE const_mod
  USE bbtrans_mod
  USE op_math_mod, only: searchsorted
  USE interpol_mod, only: hermite_interpolation, hermite_derivative
  USE handyfuncs_mod, only: symmetric_log_grid
  ! *******************************************************************
  !
  !     Contains variables connected to the atom(s)
  !     included in the calculations
  !
  ! *******************************************************************
  implicit none

  type, extends(bbtrans_type) :: transition_type
    integer :: ntrans !< transition identifier

    real(8) :: upper_ev

    integer :: iblue = 0 !< start index where line is active in overall wavelength grid
    integer :: ired  = 0 !< end index where line is active in overall wavelength grid
    real(8) :: lam_min, lam_max

    real(8), dimension(:), pointer, contiguous :: nu => null()  !< frequencies [Hz]
    real(8), dimension(:), pointer, contiguous :: lam => null() !< wavelength [AA]
    real(8), dimension(:), pointer, contiguous :: lam_input => null()
    ! The frequency arrays of a transition are pointers, so the copy each rt_patch
    ! keeps of a transition shares them with the atom instead of duplicating them
    ! for every window, direction and depth slab.  The patch indexes them with the
    ! global frequency index of its own range.
    !> frequency quadrature weight of the rate integrals: dnu/nu for a line, whose
    !> rates take Aij as a factor, and dnu/(h/4pi*nu) for a continuum (atom_mod)
    real(8), dimension(:), pointer, contiguous :: wnu => null()

  contains
    procedure :: check_range
  end type transition_type

  type, extends(transition_type) :: line_type

    logical :: hydrogen
    real(4), dimension(:), pointer, contiguous :: q => null()  !< frequencies [Doppler units]
    real(8), dimension(:), pointer, contiguous :: wq => null() !< frequency quadrature weight [Doppler units]

    real(4),dimension(:,:,:,:), pointer :: phi => null()            !< actual line profile data
    real(4),dimension(:,:,:), pointer, contiguous :: wphi => null() !< normalization such that
    ! \Sigma_mu=1^{nrays} \Sigma_nu=1^{nnu} phi(nu,v(mu)) = 1
    ! in case of strong velocities this does not automatically follow from
    ! the independent normalization of the angle and frequency quadrature

    real(8) :: f    !< oscillator strength
    real(8) :: qmax !< width of profile in doppler units
    real(8) :: q0   !< switch from linear to logarithmic frequency grid spacing

    real(8), dimension(:), allocatable :: fterm
    real(8), dimension(:), allocatable :: gwterm
    real(8), dimension(:), allocatable :: gqterm
    integer :: iso=0

  contains
    procedure :: get_lam

  end type line_type

  type, extends(transition_type) ::  cont_type
    character(len=20) :: bf_type !< either 'detail' or 'hydrogenic'

    !real(8) :: nu0 ! ionization edge frequency
    real(4) :: alpha0                                             !< cross section at ionization threshold [cm^2]
    real(4), dimension(:), pointer, contiguous :: alpha => null() !< cross sections [cm^2]
    real(8), dimension(:), pointer, contiguous :: alpha_input => null()

  contains
    procedure :: hermite_interp
  end type cont_type

contains

  !====================================================================================
  PURE SUBROUTINE hermite_interp(self, alpha, nu)
    class(cont_type), intent(in) :: self
    real(8), intent(in), dimension(:) :: nu
    real(4), intent(inout), dimension(:) :: alpha
    real(8), allocatable, dimension(:) :: dalpha_dnu, nu_input, alpha_input
    real(8) :: minl, maxl
    integer :: n
    !----------------------------------------------------------------------------------
    ! The detailed bf cross sections are Hermite-interpolated in frequency (as in the
    ! original frequency code) for backwards compatibility: Hermite interpolation is
    ! nonlinear in the abscissa, so a nu grid and a lam grid give different results.
    ! lam_input is stored increasing, so reverse it to get increasing frequencies.
    n = size(self%lam_input)
    allocate(nu_input(n), alpha_input(n), dalpha_dnu(n))
    nu_input    = HZ_TO_AA / self%lam_input(n:1:-1)
    alpha_input = self%alpha_input(n:1:-1)
    minl = minval(nu_input)
    maxl = maxval(nu_input)
    call hermite_derivative(dalpha_dnu, nu_input, alpha_input)
    call hermite_interpolation(alpha, nu, nu_input, alpha_input, dalpha_dnu)
    where (nu > maxl) alpha = 0
    where (nu < minl) alpha = 0
    deallocate(dalpha_dnu, nu_input, alpha_input)
  END SUBROUTINE hermite_interp

  !====================================================================================
  PURE SUBROUTINE check_range(self, fs_in, fe_in, fs_out, fe_out, skip)
    class(transition_type), intent(in) :: self
    integer, intent(in) :: fs_in, fe_in
    integer, intent(out) :: fs_out, fe_out
    logical, intent(out) :: skip
    !--------------------------------------------------------------------------------
    skip = .false.

    fs_out = max(self%iblue,  fs_in)
    fe_out = min(self%ired,  fe_in)

    if (fs_out > fe_out) skip = .true.
    if (fs_out > fe_in)  skip = .true.
    if (fe_out < fs_in)  skip = .true.
  END SUBROUTINE check_range
  !------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE get_lam(self, qnorm, lam, legacy)
    class(line_type), intent(inout)  :: self
    real(8), intent(in) :: qnorm
    real(8), intent(inout), dimension(:), pointer, contiguous :: lam
    real(8), dimension(:), allocatable :: q
    logical, intent(in) :: legacy
    real(8) :: dq, a, qmax, q0
    integer :: i, nnu
    logical :: linear
    !---------------------------------------------------------------------
    q0   = self%q0
    qmax = self%qmax
    nnu  = self%nnu
    if (nnu < 0) then
      call symmetric_log_grid(q0, qmax, q)
      self%nnu = size(q)
      ! q is a symmetric wavelength offset [AA] about the line center; lam0+q keeps
      ! the grid increasing
      allocate(lam(self%nnu))
      lam = self%lam0 + q
    else
      linear = qmax.le.q0

      if(.not. linear) then
        a = 10.**(q0+.5)
        qmax = log10(a*max(0.5,qmax-q0-0.5))
        a = 1/a
      endif

      allocate(q(nnu))
      dq = 2.*qmax / (nnu-1)
      q = [(i*dq - qmax, i=0, nnu-1)]

      if(.not. linear) then
        if (legacy) then
          q = q + a*sign(1.d0, q)*(10.**abs(q) - 1)
        else
          q = q + a*(10.**q - 10.**(-q))
        endif
      endif

      ! convert to wavelength [AA]
      allocate(lam(nnu))
      q = q * (qnorm * KM_TO_CM) / cc
      lam = self%lam0*(1.0 + q)
    endif
  END SUBROUTINE
  !------------------------------------------------------------------------------------

END MODULE trans_mod
