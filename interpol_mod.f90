!===============================================================================
!> All 1-D cubic Hermite and monotone PCHIP interpolation of the code: the cubic
!> basis (p = 3t^2 - 2t^3) and the Fritsch-Carlson node-slope rule live here once.
!>
!> The interior slope is the reciprocal harmonic mean documented by
!> [Ibgui et al. (2013)](https://ui.adsabs.harvard.edu/abs/2013A&A...549A.126I).
!> Boundary slopes use quadratic extrapolation from the neighbouring node slope.
!> This is the one node-slope rule of the code.
!===============================================================================
MODULE interpol_mod
  ! The node slopes (the gradient half, pchip_slopes) are computed once per driver
  ! call, while the per-point evaluation stays inline in the drivers.  Both halves
  ! are in this file, so they inline together under -Ofast / -O3 in the hot loops
  ! (rt_patch and ray_patch ghost zones).
  USE prec_mod                       ! dp; single precision written as real(4)
  USE op_math_mod, only: searchsorted, option
  USE console_mod
  implicit none
  private

  public :: pchip_mono, pchip_mono4, pchip_mono8, pchip_mono4_reg_period, pchip_slopes
  public :: hermite_interpolation, hermite_derivative, hermite_coefs, hermite_value
  public :: interp1d_sorted, interp1d_reg_period
  public :: linterpolate_3d, linterpolate_3d_rotated
  public :: pchip_integral, optical_depth, heights_from_tau

  !===============================================================================
  !> Monotone piecewise cubic Hermite (PCHIP) interpolation of yy(xx) to newx, with
  !> the node slopes of pchip_slopes.  newx must lie within [xx(1), xx(n)].
  !>
  !> Written by Richard Hoppe (2023).
  !===============================================================================
  interface pchip_mono
    module procedure pchip_mono4, pchip_mono8
  end interface pchip_mono

  !===============================================================================
  !> Monotone node slopes dydx of yy(xx), the gradient half: the reciprocal harmonic
  !> mean at interior nodes (Ibgui et al. 2013), zero at non-monotone or flat
  !> segments, quadratic extrapolation at the boundaries.  Only the nodes imin to
  !> imax are computed, and dydx is zero elsewhere.  dx and dy return the node
  !> spacings and differences for reuse by the caller's evaluation loop.
  !===============================================================================
  interface pchip_slopes
    module procedure pchip_slopes4, pchip_slopes8
  end interface pchip_slopes

  !===============================================================================
  !> Cubic Hermite interpolation of the nodes (x, y, dydx) at xx, held constant
  !> outside [x(1), x(n)].
  !===============================================================================
  interface hermite_interpolation
    module procedure hermite_interpolation_single, &
      hermite_interpolation_1d, &
      hermite_interpolation_1d4
  end interface hermite_interpolation

  !===============================================================================
  !> Cubic Hermite interpolant on an interval of width h at the local parameter t in
  !> [0,1], with end values y0, y1 and end slopes m0, m1:
  !> (1-p) y0 + p y1 + h ((t^3-2t^2+t) m0 + (t^3-t^2) m1), with p = 3t^2-2t^3.  It is
  !> the evaluation half shared by the drivers.
  !===============================================================================
  interface hermite_value
    module procedure hermite_value4, hermite_value8
  end interface hermite_value

contains

  !===============================================================================
  !> The cumulative integral of y over x from x(1), exact for the monotone PCHIP
  !> interpolant of y: h (y_i + y_i+1)/2 + h^2 (s_i - s_i+1)/12 per interval of
  !> width h, with s the pchip_slopes8 node slopes.  It has the sign of the x steps.
  !===============================================================================
  PURE SUBROUTINE pchip_integral (x, y, cum)
    real(8), intent(in) :: x(:)    !< abscissae, monotone
    real(8), intent(in) :: y(:)    !< values at x
    real(8), intent(out) :: cum(:) !< integral of y from x(1) to each x
    real(8) :: dx(size(x)-1), dy(size(x)-1), slope(size(x))
    integer :: i
    !-----------------------------------------------------------------------------
    call pchip_slopes8(x, y, 1, size(x), dx, dy, slope)
    cum(1) = 0d0
    do i=2, size(x)
      cum(i) = cum(i-1) + dx(i-1) * (0.5d0 * (y(i-1) + y(i)) &
                                     + dx(i-1) * (slope(i-1) - slope(i)) / 12d0)
    end do
  END SUBROUTINE pchip_integral

  !===============================================================================
  !> Optical depth down one column, index 1 at the top: the PCHIP integral of chi
  !> in the height z (cm, either sign convention).  The top cell gets tau = chi H,
  !> the depth of an atmosphere above the column in which chi falls off with the
  !> scale height H.  H is the smaller of the density and extinction scale heights
  !> between the two top cells: chi falls faster than rho where H- dominates
  !> (chi ~ rho n_e), and in a temperature inversion, where chi can grow outward,
  !> the density scale height keeps H finite.  The floor of 0.01 on the log step
  !> does the same in a density inversion.
  !===============================================================================
  PURE SUBROUTINE optical_depth (chi, lnrho, z, tau)
    real(8), intent(in) :: chi(:)   !< extinction [cm^-1], index 1 at the top
    real(8), intent(in) :: lnrho(:) !< ln of the density
    real(8), intent(in) :: z(:)     !< height [cm]
    real(8), intent(out) :: tau(:)  !< optical depth
    !-----------------------------------------------------------------------------
    call pchip_integral(z, chi, tau)
    tau = chi(1) * abs(z(2) - z(1)) / max(lnrho(2) - lnrho(1), log(chi(2) / chi(1)), 1d-2) &
        + abs(tau)
  END SUBROUTINE optical_depth

  !===============================================================================
  !> The heights (cm, z(1) = 0 at the top, negative below) of a column given on an
  !> optical depth scale, index 1 at the top: dz = (tau/chi) d ln tau, the scale
  !> height tau/chi integrated as a PCHIP in ln tau, which is exact for an
  !> exponential atmosphere.  Only the differences of tau enter, so the depth of
  !> the top cell is not kept; optical_depth assigns it chi H again.
  !===============================================================================
  PURE SUBROUTINE heights_from_tau (tau, chi, z)
    real(8), intent(in) :: tau(:) !< optical depth, index 1 at the top
    real(8), intent(in) :: chi(:) !< extinction [cm^-1]
    real(8), intent(out) :: z(:)  !< height [cm], 0 at the top
    !-----------------------------------------------------------------------------
    call pchip_integral(log(tau), tau / chi, z)
    z = -z
  END SUBROUTINE heights_from_tau

  !===============================================================================
  !> hermite_value for real(4).
  !===============================================================================
  elemental PURE FUNCTION hermite_value4 (t, h, y0, y1, m0, m1) result(val)
    real(4), intent(in) :: t  !< position in the interval, 0 to 1
    real(4), intent(in) :: h  !< interval width
    real(4), intent(in) :: y0 !< value at the start
    real(4), intent(in) :: y1 !< value at the end
    real(4), intent(in) :: m0 !< slope at the start
    real(4), intent(in) :: m1 !< slope at the end
    real(4) :: val, t2, t3, p
    !---------------------------------------------------------------------------
    t2 = t*t
    t3 = t2*t
    p  = 3.0*t2 - 2.0*t3
    val = (1.0 - p)*y0 + p*y1 + h*((t3 - 2.0*t2 + t)*m0 + (t3 - t2)*m1)
  END FUNCTION hermite_value4

  !===============================================================================
  !> hermite_value for real(8).
  !===============================================================================
  elemental PURE FUNCTION hermite_value8 (t, h, y0, y1, m0, m1) result(val)
    real(8), intent(in) :: t  !< position in the interval, 0 to 1
    real(8), intent(in) :: h  !< interval width
    real(8), intent(in) :: y0 !< value at the start
    real(8), intent(in) :: y1 !< value at the end
    real(8), intent(in) :: m0 !< slope at the start
    real(8), intent(in) :: m1 !< slope at the end
    real(8) :: val, t2, t3, p
    !---------------------------------------------------------------------------
    t2 = t*t
    t3 = t2*t
    p  = 3.0d0*t2 - 2.0d0*t3
    val = (1.0d0 - p)*y0 + p*y1 + h*((t3 - 2.0d0*t2 + t)*m0 + (t3 - t2)*m1)
  END FUNCTION hermite_value8

  !===============================================================================
  !> Cubic Hermite basis coefficients for the local parameter t in [0,1] on an
  !> interval of width h.  The interpolated value on the interval is
  !> c0 y0 + c1 y1 + cm0 m0 + cm1 m1 with end values y0, y1 and end slopes m0, m1.
  !> p = 3t^2 - 2t^3 is also returned for carrying auxiliary (cross-axis) slope
  !> fields in the value basis: g = (1-p) g0 + p g1.
  !===============================================================================
  elemental PURE SUBROUTINE hermite_coefs (t, h, p, c0, c1, cm0, cm1)
    real(4), intent(in) :: t    !< position in the interval, 0 to 1
    real(4), intent(in) :: h    !< interval width
    real(4), intent(out) :: p   !< 3t^2 - 2t^3
    real(4), intent(out) :: c0  !< coefficient of y0
    real(4), intent(out) :: c1  !< coefficient of y1
    real(4), intent(out) :: cm0 !< coefficient of m0
    real(4), intent(out) :: cm1 !< coefficient of m1
    real(4) :: t2, t3
    !---------------------------------------------------------------------------
    t2 = t*t
    t3 = t2*t
    p   = 3.0*t2 - 2.0*t3
    c0  = 1.0 - p
    c1  = p
    cm0 = h*(t3 - 2.0*t2 + t)
    cm1 = h*(t3 - t2)
  END SUBROUTINE hermite_coefs

  !===============================================================================
  !> pchip_slopes for real(4).
  !===============================================================================
  PURE SUBROUTINE pchip_slopes4 (xx, yy, imin, imax, dx, dy, dydx)
    real(4), intent(in) :: xx(:)    !< nodes, ascending
    real(4), intent(in) :: yy(:)    !< values at the nodes
    integer, intent(in) :: imin     !< first node whose slope is needed
    integer, intent(in) :: imax     !< last node whose slope is needed
    real(4), intent(out) :: dx(:)   !< node spacings xx(i+1) - xx(i)
    real(4), intent(out) :: dy(:)   !< differences yy(i+1) - yy(i)
    real(4), intent(out) :: dydx(:) !< node slopes
    integer :: n, i
    real(8) :: alpha, d1, d2
    !---------------------------------------------------------------------------
    n = size(xx)

    dydx = 0.
    dx = xx(2:n) - xx(1:n-1)
    dy = yy(2:n) - yy(1:n-1)

    if (n == 2) then
      dydx(1) = dy(1)/dx(1)
      dydx(2) = dydx(1)
    endif

    do i=max(imin, 2), min(imax, n-1)
      if (dy(i-1) == 0 .or. dy(i) == 0) cycle

      d1 = dy(i-1)/dx(i-1)
      d2 = dy(i)/dx(i)

      if (sign(1.d0, d1) == sign(1.d0, d2)) then
        ! reciprocal (harmonic) mean derivative
        alpha = (1. + dx(i-1)/(dx(i) + dx(i-1))) / 3.0
        dydx(i) = real(1.d0 / ((1.d0 - alpha)/d1 + alpha/d2))
      endif
    end do

    ! boundary slopes from quadratic polynomial using 2 points and 1 dydx
    if (imin == 1) then
      dydx(1) = 2.*dy(1)/dx(1) - dydx(2)
    endif

    if (imax == n) then
      dydx(n) = 2.*dy(n-1)/dx(n-1) - dydx(n-1)
    endif
  END SUBROUTINE pchip_slopes4

  !===============================================================================
  !> pchip_slopes for real(8).
  !===============================================================================
  PURE SUBROUTINE pchip_slopes8 (xx, yy, imin, imax, dx, dy, dydx)
    real(8), intent(in) :: xx(:)    !< nodes, ascending
    real(8), intent(in) :: yy(:)    !< values at the nodes
    integer, intent(in) :: imin     !< first node whose slope is needed
    integer, intent(in) :: imax     !< last node whose slope is needed
    real(8), intent(out) :: dx(:)   !< node spacings xx(i+1) - xx(i)
    real(8), intent(out) :: dy(:)   !< differences yy(i+1) - yy(i)
    real(8), intent(out) :: dydx(:) !< node slopes
    integer :: n, i
    real(8) :: alpha, d1, d2
    !---------------------------------------------------------------------------
    n = size(xx)

    dydx = 0.
    dx = xx(2:n) - xx(1:n-1)
    dy = yy(2:n) - yy(1:n-1)

    if (n == 2) then
      dydx(1) = dy(1)/dx(1)
      dydx(2) = dydx(1)
    endif

    do i=max(imin, 2), min(imax, n-1)
      if (dy(i-1) == 0 .or. dy(i) == 0) cycle

      d1 = dy(i-1)/dx(i-1)
      d2 = dy(i)/dx(i)

      if (sign(1.d0, d1) == sign(1.d0, d2)) then
        ! reciprocal (harmonic) mean derivative
        alpha = (1. + dx(i-1)/(dx(i) + dx(i-1))) / 3.0
        dydx(i) = 1. / ((1. - alpha)/d1 + alpha/d2)
      endif
    end do

    ! boundary slopes from quadratic polynomial using 2 points and 1 dydx
    if (imin == 1) then
      dydx(1) = 2.*dy(1)/dx(1) - dydx(2)
    endif

    if (imax == n) then
      dydx(n) = 2.*dy(n-1)/dx(n-1) - dydx(n-1)
    endif
  END SUBROUTINE pchip_slopes8

  !===============================================================================
  !> pchip_mono for real(4), with the slope newdydx on request.
  !===============================================================================
  SUBROUTINE pchip_mono4(xx, yy, newx, newy, sorted, newdydx)
    real(4), intent(in) :: xx(:)                 !< X values of the data points
    real(4), intent(in) :: yy(:)                 !< Y values of the data points
    real(4), intent(in) :: newx(:)               !< X values to interpolate at
    real(4), intent(out) :: newy(:)              !< Interpolated Y values
    logical, intent(in), optional :: sorted      !< newx ascending, scanned instead of bisected (default .false.)
    real(4), intent(out), optional :: newdydx(:) !< Derivative at interpolation points
    real(4), save, allocatable :: dydx(:)   ! derivative
    real(4), save, allocatable :: dx(:)
    real(4), save, allocatable :: dy(:)
    logical :: is_sorted, get_derivative
    integer :: n                 ! Number of data points
    integer :: n_interp          ! Number of interpolation points
    integer :: j                 ! Loop index for interpolation points
    integer :: i                 ! Loop index for data points
    integer :: imin, imax
    real(4) :: t1, t2, t3, p
    !$omp threadprivate(dydx, dx, dy)
    !-----------------------------------------------------------------------------
    n = size(xx)
    n_interp = size(newx)
    is_sorted = option(sorted, .false.)
    get_derivative = present(newdydx)

    if (allocated(dydx)) then
      if (size(dydx) /= n) then
        deallocate(dydx)
        deallocate(dx)
        deallocate(dy)
        allocate(dydx(n))
        allocate(dx(n-1))
        allocate(dy(n-1))
      endif
    else
      allocate(dydx(n))
      allocate(dx(n-1))
      allocate(dy(n-1))
    endif

    imin = searchsorted(xx, minval(newx), outbounds_ok=.true.)
    imax = searchsorted(xx, maxval(newx), outbounds_ok=.true.)+1
    imin = max(imin, 1)
    imax = min(imax, n)

    call pchip_slopes4(xx, yy, imin, imax, dx, dy, dydx)

    i = 1
    do j=1, n_interp
      if (is_sorted) then
        ! Find the interval containing newx(j)
        do i=i, n - 1
          if (newx(j) <= xx(i + 1)) then
            exit ! Found the interval
          end if
        end do
        i = min(i, n-1)
      else
        i = searchsorted(xx, newx(j), outbounds_ok=.true.)
      endif

      ! Calculate the interpolation parameter
      t1 = (newx(j) - xx(i)) / dx(i)

      t2 = t1*t1
      t3 = t2*t1
      p = 3.0*t2 - 2.0*t3

      newy(j) = (1.0 - p)*yy(i) + p*yy(i+1) + dx(i)*((t3 - 2.0*t2 + t1)*dydx(i) + (t3-t2)*dydx(i+1))

      if (get_derivative) then
        newdydx(j) =  6.0*(t1 - t2)          * dy(i)/dx(i) &
                   + (3.0*t2 - 4.0*t1 + 1.0) * dydx(i)     &
                   + (3.0*t2 - 2.0*t1)       * dydx(i+1)
      end if
    end do
  END SUBROUTINE pchip_mono4

  !===============================================================================
  !> pchip_mono for real(8).
  !===============================================================================
  SUBROUTINE pchip_mono8(xx, yy, newx, newy, sorted)
    real(8), intent(in) :: xx(:)            !< X values of the data points
    real(8), intent(in) :: yy(:)            !< Y values of the data points
    real(8), intent(in) :: newx(:)          !< X values to interpolate at
    real(8), intent(out) :: newy(:)         !< Interpolated Y values
    logical, intent(in), optional :: sorted !< newx ascending, scanned instead of bisected (default .false.)
    real(8), save, allocatable :: dydx(:)   ! derivative
    real(8), save, allocatable :: dx(:)
    real(8), save, allocatable :: dy(:)
    logical :: is_sorted
    integer :: n                 ! Number of data points
    integer :: n_interp          ! Number of interpolation points
    integer :: j                 ! Loop index for interpolation points
    integer :: i                 ! Loop index for data points
    integer :: imin, imax
    real(8) :: t1, t2, t3, p
    !$omp threadprivate(dydx, dx, dy)
    !-----------------------------------------------------------------------------
    n = size(xx)
    n_interp = size(newx)
    is_sorted = option(sorted, .false.)

    if (allocated(dydx)) then
      if (size(dydx) /= n) then
        deallocate(dydx)
        deallocate(dx)
        deallocate(dy)
        allocate(dydx(n))
        allocate(dx(n-1))
        allocate(dy(n-1))
      endif
    else
      allocate(dydx(n))
      allocate(dx(n-1))
      allocate(dy(n-1))
    endif

    imin = searchsorted(xx, minval(newx), outbounds_ok=.true.)
    imax = searchsorted(xx, maxval(newx), outbounds_ok=.true.)+1
    imin = max(imin, 1)
    imax = min(imax, n)

    call pchip_slopes8(xx, yy, imin, imax, dx, dy, dydx)

    i = 1
    do j=1, n_interp
      if (is_sorted) then
        ! Find the interval containing newx(j)
        do i=i, n - 1
          if (newx(j) <= xx(i + 1)) then
            exit ! Found the interval
          end if
        end do
        i = min(i, n-1)
      else
        i = searchsorted(xx, newx(j), outbounds_ok=.true.)
      endif

      ! Calculate the interpolation parameter
      t1 = (newx(j) - xx(i)) / dx(i)

      t2 = t1*t1
      t3 = t2*t1
      p = 3.0*t2 - 2.0*t3

      newy(j) = (1.0 - p)*yy(i) + p*yy(i+1) + dx(i)*((t3 - 2.0*t2 + t1)*dydx(i) + (t3-t2)*dydx(i+1))
    end do
  END SUBROUTINE pchip_mono8

  !===============================================================================
  !> Monotone PCHIP interpolation of the periodic data yy, given on the grid points
  !> 0 to n-1, to the points new_x0 + (j-1) new_dx, in grid units.
  !>
  !> Written by Richard Hoppe (2023).
  !===============================================================================
  SUBROUTINE pchip_mono4_reg_period(yy, new_x0, new_dx, newy, newdydx)
    real(4), intent(in) :: yy(:)                 !< Y values of the data points
    real(8), intent(in) :: new_x0                !< Starting x value in grid units
    real(8), intent(in) :: new_dx                !< Spacing of new x values in grid units
    real(4), intent(out) :: newy(:)              !< Interpolated Y values
    real(4), intent(out), optional :: newdydx(:) !< Derivative per grid unit at interpolation points
    real(4), save, allocatable :: dydx(:)   ! derivative
    real(4), save, allocatable :: dy(:)
    integer :: n                 ! Number of data points
    integer :: n_interp          ! Number of interpolation points
    integer :: j                 ! Loop index for interpolation points
    integer :: i, i1, i2         ! Loop index for data points
    real(4) :: t1, t2, t3, p, ratio, newx, x0
    real(8) :: dydx8, d1, d2
    logical :: get_derivative
    !$omp threadprivate(dydx, dy)
    !-----------------------------------------------------------------------------
    n = size(yy)
    n_interp = size(newy)
    get_derivative = present(newdydx)

    ! the grid-aligned copy fast path has no node slopes to hand out, so it is
    ! bypassed when the derivative is requested; exact node points then fall out
    ! of the general evaluation with t1=0, giving newy=yy and newdydx=dydx
    if (.not. get_derivative .and. &
        new_x0 == 0 .and. new_dx - nint(new_dx) < 2.0 * epsilon(new_dx)) then
      ! newx is on the same grid as xx, just copy the values
      x0 = new_x0 - new_dx
      do j=1, n_interp
        newx = x0 + j*new_dx
        i = modulo(nint(newx), n) + 1
        newy(j) = yy(i)
      end do
      return
    endif

    if (allocated(dydx)) then
      if (size(dydx) /= n) then
        deallocate(dydx)
        deallocate(dy)
        allocate(dydx(n))
        allocate(dy(n))
      end if
    else
      allocate(dydx(n))
      allocate(dy(n))
    end if

    dydx = 0.
    dy(1:n-1) = yy(2:n) - yy(1:n-1)
    dy(n) = yy(1) - yy(n)

    do i=1, n
      i1 = modulo(i-2, n) + 1
      i2 = modulo(i-1, n) + 1
      if (dy(i1) == 0 .or. dy(i2) == 0) cycle

      d1 = dy(i1)
      d2 = dy(i2)

      if (sign(1.d0, d1) == sign(1.d0, d2)) then
        ! harmonic mean derivative
        dydx8 = 2*d1*d2 / (d1 + d2)
        dydx(i) = real(dydx8)
      endif
    end do

    x0 = new_x0 - new_dx
    do j=1, n_interp
      ratio = x0 + j*new_dx
      i = floor(ratio)
      i1 = modulo(i,   n) + 1
      i2 = modulo(i+1, n) + 1

      ! Calculate the interpolation parameter
      t1 = ratio - i

      t2 = t1*t1
      t3 = t2*t1
      p = 3.0*t2 - 2.0*t3

      newy(j) = (1.0 - p)*yy(i1) + p*yy(i2) + (t3 - 2.0*t2 + t1)*dydx(i1) + (t3-t2)*dydx(i2)

      if (get_derivative) then
        ! d(newy)/dx per grid unit (interval width h = 1); dy(i1) is the segment
        ! difference across [i1,i2]
        newdydx(j) =  6.0*(t1 - t2)          * dy(i1)   &
                   + (3.0*t2 - 4.0*t1 + 1.0) * dydx(i1) &
                   + (3.0*t2 - 2.0*t1)       * dydx(i2)
      end if
    end do
  END SUBROUTINE pchip_mono4_reg_period

  !===============================================================================
  !>  Node slopes over a whole array via the canonical monotone rule (thin wrapper
  !>  used by the Hermite-interpolation drivers below and by external callers that
  !>  precompute the gradient once, e.g. cont_type%hermite_interp in trans_mod).
  !===============================================================================
  PURE SUBROUTINE hermite_derivative(dydx,x,y)
    real(8),dimension(:),intent(in)  :: x    !< nodes, ascending
    real(8),dimension(:),intent(in)  :: y    !< values at the nodes
    real(8),dimension(:),intent(out) :: dydx !< node slopes
    real(8) :: dx(size(x)-1), dy(size(x)-1)
    integer :: n
    !---------------------------------------------------------------------------
    n = size(x)
    call pchip_slopes8(x, y, 1, n, dx, dy, dydx)
  END SUBROUTINE hermite_derivative

  !===============================================================================
  !> hermite_interpolation at one point.
  !===============================================================================
  PURE SUBROUTINE hermite_interpolation_single(yy,xx,x,y,dydx)
    real(kind=8), dimension(:), intent(in)  :: x    !< nodes, ascending
    real(kind=8), dimension(:), intent(in)  :: y    !< values at the nodes
    real(kind=8), dimension(:), intent(in)  :: dydx !< node slopes (hermite_derivative)
    real(8),intent(in) :: xx                        !< point to interpolate to
    real(8),intent(out):: yy                        !< interpolated value
    integer :: n,i
    real(8) :: dx
    !---------------------------------------------------------------------------
    n=size(x)

    if (xx<=x(1)) then
      yy = y(1)
    else if (xx>=x(n)) then
      yy = y(n)
    else
      i = searchsorted(x, xx, outbounds_ok=.true.)
      dx = x(i+1)-x(i)
      yy = hermite_value8((xx-x(i)) / dx, dx, y(i), y(i+1), dydx(i), dydx(i+1))
    end if
  END SUBROUTINE hermite_interpolation_single

  !===============================================================================
  !> hermite_interpolation at an array of points.
  !===============================================================================
  PURE SUBROUTINE hermite_interpolation_1d(yy,xx,x,y,dydx)
    real(kind=8), dimension(:), intent(in)  :: x    !< nodes, ascending
    real(kind=8), dimension(:), intent(in)  :: y    !< values at the nodes
    real(kind=8), dimension(:), intent(in)  :: dydx !< node slopes (hermite_derivative)
    real(8), dimension(:),intent(in) :: xx          !< points to interpolate to
    real(8), dimension(:),intent(out):: yy          !< interpolated values
    integer :: k
    !---------------------------------------------------------------------------
    do k=1,size(xx)
      call hermite_interpolation_single(yy(k),xx(k),x,y,dydx)
    end do
  END SUBROUTINE hermite_interpolation_1d

  !===============================================================================
  !> hermite_interpolation at an array of points, into real(4).
  !===============================================================================
  PURE SUBROUTINE hermite_interpolation_1d4(yy,xx,x,y,dydx)
    real(kind=8), dimension(:), intent(in)  :: x    !< nodes, ascending
    real(kind=8), dimension(:), intent(in)  :: y    !< values at the nodes
    real(kind=8), dimension(:), intent(in)  :: dydx !< node slopes (hermite_derivative)
    real(8), dimension(:),intent(in) :: xx          !< points to interpolate to
    real(4), dimension(:),intent(out):: yy          !< interpolated values
    real(8) :: y8
    integer :: k
    !---------------------------------------------------------------------------
    do k=1,size(xx)
      call hermite_interpolation_single(y8,xx(k),x,y,dydx)
      yy(k) = y8
    end do
  END SUBROUTINE hermite_interpolation_1d4

  !===============================================================================
  !> Linear interpolation of yy(xx) to the ascending points newx.  Outside
  !> [xx(1), xx(n)] newy is fill, or the end value without it.
  !===============================================================================
  SUBROUTINE interp1d_sorted (xx, yy, newx, newy, fill)
    real(4), dimension(:), intent(in) :: newx  !< points to interpolate to, ascending
    real(4), dimension(:), intent(in) :: xx    !< nodes, ascending
    real(4), dimension(:), intent(in) :: yy    !< values at the nodes
    real(4), dimension(:), intent(out) :: newy !< interpolated values
    real(4) :: m, left, right
    real(4), optional :: fill                  !< value outside the nodes (default: the end values)
    integer :: ii, i, nnew, nold
    !-----------------------------------------------------------------------------
    nnew = size(newx, 1)
    nold = size(xx, 1)

    if (present(fill)) then
      left = fill
      right = fill
    else
      left = yy(1)
      right = yy(nold)
    endif

    i = 1
    do ii=1, nnew
      do while (xx(i+1) < newx(ii))
        i = i + 1
        if (i >= nold) then
          i = nold
          exit
        endif
      end do

      if (i == nold) then
        newy(ii:nnew) = right
        exit
      elseif (xx(1) > newx(ii)) then
        newy(ii) = left
      else
        m = (yy(i+1) - yy(i)) /  (xx(i+1) - xx(i))
        newy(ii) = m * (newx(ii) - xx(i)) + yy(i)
      endif
    end do
  END SUBROUTINE interp1d_sorted

  !===============================================================================
  !> Linear interpolation of the periodic data yy, given on the grid points 0 to n-1,
  !> to the points new_x0 + (j-1) new_dx, in grid units.
  !===============================================================================
  SUBROUTINE interp1d_reg_period (yy, new_x0, new_dx, newy)
    real(8), intent(in) :: new_x0              !< Starting x value in grid units
    real(8), intent(in) :: new_dx              !< Spacing of new x values in grid units
    real(4), dimension(:), intent(in) :: yy    !< values on the periodic grid
    real(4), dimension(:), intent(out) :: newy !< interpolated values
    real(4) :: ratio, delta
    integer :: i, i1, i2, nnew, nold
    real(8) :: x0
    !-----------------------------------------------------------------------------
    nnew = size(newy, 1)
    nold = size(yy, 1)
    x0 = new_x0 - new_dx

    do i=1, nnew
      ratio = x0 + i*new_dx
      i1 = floor(ratio)
      delta = ratio - i1
      i1 = modulo(i1, nold) + 1
      i2 = modulo(i1, nold) + 1
      newy(i) = (1.0 - delta) * yy(i1) + delta * yy(i2)
    end do
  END SUBROUTINE interp1d_reg_period

  !===============================================================================
  !> Interpolates the 3D field f, periodic in x and y, from its grid (spacings dx,
  !> dy and heights oldz) to the grid of fout (horizontal origin lower, spacings ds
  !> and heights newz), axis by axis, with monotone PCHIP or linear interpolation.
  !>
  !> dfdx/dfdy/dfdz optionally return the partial derivatives of the interpolant
  !> at the output points, in units of the respective input coordinate (dx, dy,
  !> oldz).  The horizontal partials are captured in the corresponding pass and
  !> then carried through the remaining passes as plain value fields; they (and
  !> dfdz) require pchip interpolation, since the linear branches have no
  !> continuous derivative to hand out.
  !===============================================================================
  SUBROUTINE linterpolate_3d(f, fout, dx, dy, oldz, lower, ds, newz, pchip, dfdx, dfdy, dfdz)
    real(4), intent(in),  dimension(:,:,:) :: f              !< field on the input grid
    real(4), intent(out), dimension(:,:,:) :: fout           !< field on the output grid
    real(8), intent(in) :: dx                                !< input grid spacing in x
    real(8), intent(in) :: dy                                !< input grid spacing in y
    real(8), intent(in) :: lower(2)                          !< x and y of the first output point
    real(8), intent(in) :: ds(2)                             !< output grid spacings in x and y
    real(4), dimension(:), intent(in) :: oldz                !< heights of the input grid, ascending
    real(4), dimension(:), intent(in) :: newz                !< heights of the output grid, ascending
    logical, optional, intent(in) :: pchip                   !< monotone PCHIP (default), else linear
    real(4), intent(out), dimension(:,:,:), optional :: dfdx !< df/dx at the output points
    real(4), intent(out), dimension(:,:,:), optional :: dfdy !< df/dy at the output points
    real(4), intent(out), dimension(:,:,:), optional :: dfdz !< df/dz at the output points
    real(4), dimension(:,:,:), allocatable :: fx
    real(4), dimension(:,:,:), allocatable :: fxy
    real(4), dimension(:,:,:), allocatable :: gx, gxy, hxy
    real(8) :: df(2), offset(2)
    integer :: i,j,k, nx, ny, nz, ox,oy,oz, zs, ze
    logical :: do_pchip, get_dx, get_dy, get_dz
    !-----------------------------------------------------------------------------
    ox = size(f,1); nx = size(fout,1)
    oy = size(f,2); ny = size(fout,2)
    oz = size(f,3); nz = size(fout,3)
    do_pchip = option(pchip, .true.)
    get_dx = present(dfdx)
    get_dy = present(dfdy)
    get_dz = present(dfdz)
    if ((get_dx .or. get_dy .or. get_dz) .and. .not. do_pchip) &
      error stop 'linterpolate_3d: derivative output requires pchip interpolation'
    allocate(fx(nx,oy,oz))
    allocate(fxy(nx,ny,oz))
    if (get_dx) allocate(gx(nx,oy,oz), gxy(nx,ny,oz))
    if (get_dy) allocate(hxy(nx,ny,oz))
    df = ds / [dx, dy]
    offset = lower / [dx, dy]

    zs = searchsorted(oldz, newz(1), outbounds_ok=.true.)
    ze = searchsorted(oldz, newz(nz), outbounds_ok=.true.)
    zs = max(1, zs-1)
    ze = min(size(oldz), ze+2)

    ! x-dimension interpolation
    do k = zs, ze
      do j = 1, oy
        if (get_dx) then
          call pchip_mono4_reg_period(f(:,j,k), offset(1), df(1), fx(:,j,k), newdydx=gx(:,j,k))
        elseif (do_pchip) then
          call pchip_mono4_reg_period(f(:,j,k), offset(1), df(1), fx(:,j,k))
        else
          call interp1d_reg_period(f(:,j,k), offset(1), df(1), fx(:,j,k))
        endif
      end do
      ! y-dimension interpolation
      do i = 1, nx
        if (do_pchip) then
          if (get_dy) then
            call pchip_mono4_reg_period(fx(i,:,k), offset(2), df(2), fxy(i,:,k), newdydx=hxy(i,:,k))
          else
            call pchip_mono4_reg_period(fx(i,:,k), offset(2), df(2), fxy(i,:,k))
          endif
          if (get_dx) &
            call pchip_mono4_reg_period(gx(i,:,k), offset(2), df(2), gxy(i,:,k))
        else
          call interp1d_reg_period(fx(i,:,k), offset(2), df(2), fxy(i,:,k))
        endif
      end do
    end do
    ! z-dimension interpolation
    do j = 1, ny
      do i = 1, nx
        if (do_pchip) then
          if (get_dz) then
            call pchip_mono4(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:), sorted=.true., &
                             newdydx=dfdz(i,j,:))
          else
            call pchip_mono4(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:), sorted=.true.)
          endif
          if (get_dx) &
            call pchip_mono4(oldz(zs:ze), gxy(i,j,zs:ze), newz, dfdx(i,j,:), sorted=.true.)
          if (get_dy) &
            call pchip_mono4(oldz(zs:ze), hxy(i,j,zs:ze), newz, dfdy(i,j,:), sorted=.true.)
        else
          call interp1d_sorted(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:))
        endif
      end do
    end do
    ! grid-unit horizontal derivatives -> per input-coordinate unit
    if (get_dx) dfdx = dfdx / real(dx)
    if (get_dy) dfdy = dfdy / real(dy)
  END SUBROUTINE linterpolate_3d

  !===============================================================================
  !> Interpolates f as linterpolate_3d, but with input level k shifted horizontally
  !> by (bot2pivot - oldz(k)) vec(1:2)/vec(3), so that the output columns run along
  !> the ray direction vec.
  !>
  !> The horizontal passes shear per snapshot z-level, so the z-columns of the
  !> output run along the ray direction vec: the optional dfdz is therefore the
  !> along-ray derivative of the interpolant, parameterized by z (in oldz units).
  !===============================================================================
  SUBROUTINE linterpolate_3d_rotated(f, fout, dx, dy, oldz, lower, ds, newz, vec, bot2pivot, pchip, dfdz)
    real(4), intent(in),  dimension(:,:,:) :: f    !< field on the input grid
    real(4), intent(out), dimension(:,:,:) :: fout !< field on the output grid
    real(4), dimension(:), intent(in) :: oldz      !< heights of the input grid, ascending
    real(4), dimension(:), intent(in) :: newz      !< heights of the output grid, ascending
    real(4), dimension(:,:,:), allocatable :: fx, fxy
    real(8), intent(in) :: dx                                !< input grid spacing in x
    real(8), intent(in) :: dy                                !< input grid spacing in y
    real(8), intent(in) :: lower(2)                          !< x and y of the first output point
    real(8), intent(in) :: ds(2)                             !< output grid spacings in x and y
    real(8), intent(in) :: vec(3)                            !< ray direction
    real(8), intent(in) :: bot2pivot                         !< height without horizontal shift
    logical, optional, intent(in) :: pchip                   !< monotone PCHIP (default), else linear
    real(4), intent(out), dimension(:,:,:), optional :: dfdz !< along-ray df/dz at the output points
    real(8) :: delta(2), deltak(2), offset(2), df(2)
    integer :: i,j,k, nx, ny, nz, ox,oy,oz, zs, ze
    logical :: do_pchip, get_dz
    !-----------------------------------------------------------------------------
    ox = size(f,1); nx = size(fout,1)
    oy = size(f,2); ny = size(fout,2)
    oz = size(f,3); nz = size(fout,3)
    do_pchip = option(pchip, .true.)
    get_dz = present(dfdz)
    if (get_dz .and. .not. do_pchip) &
      error stop 'linterpolate_3d_rotated: derivative output requires pchip interpolation'
    allocate(fx(nx,oy,oz))
    allocate(fxy(nx,ny,oz))
    delta = vec(1:2) / (vec(3) * [dx, dy])
    offset = lower / [dx, dy]
    df = ds / [dx, dy]

    zs = searchsorted(oldz, newz(1), outbounds_ok=.true.)
    ze = searchsorted(oldz, newz(nz), outbounds_ok=.true.)
    zs = max(1, zs-1)
    ze = min(size(oldz), ze+2)

    do k = zs, ze
      deltak = (bot2pivot - oldz(k))*delta + offset
      ! x-dimension linear interpolation
      do j = 1, oy
        if (do_pchip) then
          call pchip_mono4_reg_period(f(:,j,k), deltak(1), df(1), fx(:,j,k))
        else
          call interp1d_reg_period(f(:,j,k), deltak(1), df(1), fx(:,j,k))
        endif
      end do
      ! y-dimension linear interpolation
      do i = 1, nx
        if (do_pchip) then
          call pchip_mono4_reg_period(fx(i,:,k), deltak(2), df(2), fxy(i,:,k))
        else
          call interp1d_reg_period(fx(i,:,k), deltak(2), df(2), fxy(i,:,k))
        endif
      end do
    end do
    ! z-dimension interpolation
    do j = 1, ny
      do i = 1, nx
        if (do_pchip) then
          if (get_dz) then
            call pchip_mono4(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:), sorted=.true., &
                             newdydx=dfdz(i,j,:))
          else
            call pchip_mono4(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:), sorted=.true.)
          endif
        else
          call interp1d_sorted(oldz(zs:ze), fxy(i,j,zs:ze), newz, fout(i,j,:))
        endif
      end do
    end do
  END SUBROUTINE linterpolate_3d_rotated

END MODULE interpol_mod
