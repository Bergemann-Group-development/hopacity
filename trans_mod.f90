MODULE trans_mod
  USE const_mod
  USE bbtrans_mod
  USE op_math_mod, only: searchsorted
  USE mmath_mod, only: hermite_interpolation, hermite_derivative
  USE handyfuncs_mod, only: symmetric_log_grid
  ! *******************************************************************
  !
  !     Contains variables connected to the atom(s)
  !     included in the calculations
  !
  ! *******************************************************************
  implicit none

  type anisotropy_type
    real(8),dimension(:,:,:),allocatable :: J00,J20
    complex(8),dimension(:,:,:),allocatable :: J2p1,J2m1,J2p2,J2m2
  end type anisotropy_type

  type, extends(bbtrans_type) :: transition_type
    integer :: ntrans ! transition identifier

    real(8) :: upper_ev

    integer :: ired  = 0 ! start index where line is active in overall frequency grid
    integer :: iblue = 0 ! end index where line is active in overall frequency grid
    real(8) :: nu_min, nu_max

    real(8), dimension(:), pointer, contiguous :: nu => null()  ! frequencies [Hz]
    real(8), dimension(:), pointer, contiguous :: nu_input => null()
    real(8), dimension(:), allocatable :: wnu ! frequency quadrature weight [Hz]

    integer, dimension(:), pointer, contiguous :: if => null()
    !real(8),dimension(:,:,:),  allocatable :: rij,rji
    !real(8),dimension(:,:,:),  allocatable :: qrad ! heating rate
    !real(8),dimension(:,:,:),  allocatable :: gamma_ij,gamma_ji

  contains
    procedure :: overlaps
    procedure :: check_range
  end type transition_type

  type, extends(transition_type) :: line_type

    character(len=72) :: profile_type ! what kind of profile to use,
    ! for now only 'voigt' or
    ! 'hlinop' for hydrogen

    real(8) :: lam0       ! line center vacuum wavelength [cm]
    logical :: hydrogen
    real(4),dimension(:), allocatable :: q  ! frequencies [Doppler units]
    real(8),dimension(:), allocatable :: wq ! frequency quadrature weight [Doppler units]

    real(4),dimension(:,:,:,:), pointer :: phi => null() ! actual line profile data
    real(4),dimension(:,:,:), pointer, contiguous :: wphi => null() ! normalization such that
    ! \Sigma_mu=1^{nrays} \Sigma_nu=1^{nnu} phi(nu,v(mu)) = 1
    ! in case of strong velocities this does not automatically follow from
    ! the independent normalization of the angle and frequency quadrature

    real(8) :: f    ! oscillator strength
    real(8) :: qmax ! width of profile in doppler units
    real(8) :: q0   ! switch from linear to logarithmic frequency grid spacing

    real(8), dimension(:), allocatable :: fterm
    real(8), dimension(:), allocatable :: gwterm
    real(8), dimension(:), allocatable :: gqterm
    integer :: iso=0

  contains
    procedure :: get_nu
    procedure :: freql

  end type line_type

  type, extends(transition_type) ::  cont_type

    character(len=20) :: bf_type ! either 'detail' or 'hydrogenic'

    !real(8) :: nu0 ! ionization edge frequency
    real(4) :: alpha0 ! cross section at ionization threshold [cm^2]
    real(4), dimension(:), pointer, contiguous :: alpha => null() ! cross sections [cm^2]
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
    real(8), allocatable, dimension(:) :: dalpha_dnu
    real(8) :: minl, maxl
    !----------------------------------------------------------------------------------
    allocate(dalpha_dnu, mold=self%nu_input)
    minl = minval(self%nu_input)
    maxl = maxval(self%nu_input)
    call hermite_derivative(dalpha_dnu, self%nu_input, self%alpha_input)
    call hermite_interpolation(alpha, nu, self%nu_input, self%alpha_input, dalpha_dnu)
    where (nu > maxl) alpha = 0
    where (nu < minl) alpha = 0
    deallocate(dalpha_dnu)
  END SUBROUTINE hermite_interp

  !====================================================================================
  PURE SUBROUTINE check_range(self, fs_in, fe_in, fs_out, fe_out, skip)
    class(transition_type), intent(in) :: self
    integer, intent(in) :: fs_in, fe_in
    integer, intent(out) :: fs_out, fe_out
    logical, intent(out) :: skip
    !--------------------------------------------------------------------------------
    skip = .false.

    fs_out = max(self%ired,  fs_in)
    fe_out = min(self%iblue, fe_in)

    if (fs_out > fe_out) skip = .true.
    if (fs_out > fe_in)  skip = .true.
    if (fe_out < fs_in)  skip = .true.
  END SUBROUTINE check_range
  !------------------------------------------------------------------------------------

  !====================================================================================
  logical pure function overlaps(self, other)
    class(transition_type), intent(in) :: self, other
    !--------------------------------------------------------------------------------
    overlaps = .false.

    if (self%ired > other%ired .and. self%ired < other%iblue) overlaps = .true.
    if (other%ired > self%ired .and. other%ired < self%iblue) overlaps = .true.

  end function overlaps
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE get_nu(self, qnorm, nu, legacy)
    class(line_type), intent(inout)  :: self
    real(8), intent(in) :: qnorm
    real(8), intent(inout), dimension(:), pointer, contiguous :: nu
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
      ! convert to [Hz]
      allocate(nu(self%nnu))
      q = (CC / self%nu0 * CM_TO_AA) - q
      nu = CC / (q * AA_TO_CM)
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

      ! convert to [Hz]
      allocate(nu(nnu))
      q = q * qnorm * KM_TO_CM/cc
      nu = self%nu0*(1.0 + q)
    endif
  END SUBROUTINE

  !====================================================================================
  SUBROUTINE freql(self, symmetric)
    !***********************************************************************

    !  gives  frequency quadrature points and corresponding weights.
    !  one- or two-sided cases.  the points and weights are based on
    !  the trapezoidal rule applied to a mapping function  that maps
    !  equidistant  points in  x to  equidistant points in  q  in the
    !  core (q .lt. q0), and equidistant points in log(q) in the wings
    !  (q .gt. q0).

    !  nq     : total number of quadrature points                (in)
    !  qmax   : maximum value of frequency (doppler units)       (in)
    !  ind    : =1 for one-sided and =2 for two-sided cases      (in)
    !  q0     : transition from core to wing (doppler units)     (in)
    !  q      : quadrature points                               (out)
    !  wq     : quadrature weights                              (out)

    !  coded by: ake nordlund (oct-1981).

    !  revised by mats carlsson (jan-1984): new options
    !  qmax.le.q0   q linear in frequency
    !  qmax or q0.lt.0  q given and only weights are calculated

    !:
    !: freql  88-05-03  modifications: (mats carlsson)
    !:        added options gave incorrect weights for ind=1
    !:        this affected cooling rates for lines with frequency
    !:        quadrature given with qmax=q0 or qmax.lt.0 or q0.lt.0
    !:        the corrected weights are twice the old ones
    !:
    !***********************************************************************

    class(line_type), intent(inout)  :: self
    real(8),dimension(:),allocatable :: q,wq

    real(8)  :: dq, x, x10, a, half, xmax, dx,qmax,q0, start
    integer :: nnu, j
    logical :: symmetric
    !---------------------------------------------------------------------
    q0=self%q0
    qmax=self%qmax
    nnu=self%nnu
    allocate(q(nnu))
    allocate(wq(nnu))

    if (qmax.le.q0) then
      ! linear spacing
      dq = qmax / (nnu-1)
      q(1) = 0

      if (symmetric) then
        dq = 2.*dq
        q(1) = -qmax
      endif

      do j=2,nnu
        q(j) = q(j-1) + dq
      end do

    else if(qmax.ge.0.0 .and. q0.ge.0.0) then
      a=10.**(q0+.5)
      half=0.5
      xmax=log10(a*max(half,qmax-q0-half))

      dx=xmax/(nnu-1)
      start=0

      if (symmetric) then
        dx = 2.*dx
        start = -xmax
      endif

      do j=1,nnu
        x=start+(j-1)*dx
        x10=10.**x
        q(j)=x+(x10-1./x10)/a
      end do
    end if

    ! convert to [Hz]
    allocate(self%q(nnu))
    allocate(self%nu(nnu))
    self%q = q*KM_TO_CM/cc
    self%nu = self%nu0*(1.0 + self%q)
    self%nu_min = self%nu(1)
    self%nu_max = self%nu(nnu)
  END SUBROUTINE freql
  !------------------------------------------------------------------------------------

END MODULE trans_mod
