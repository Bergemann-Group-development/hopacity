MODULE level_mod
  USE prec_mod
  implicit none

  type :: state_type
    real(dp) :: xion=999.99
    real(dp), dimension(:,:,:), allocatable :: n, f
  end type

  type, extends(state_type) :: level_type
    integer :: i
    real(dp) :: xl, g=0
    real(dp), dimension(:,:,:), allocatable :: w
  contains
    procedure :: fischel_sparks
  end type

  public :: wcalc
  private :: exp_neg, exp_neg0, exp4, exp4_array   ! from exp_neg.inc, local to each includer

contains

  !===============================================================================
  elemental real(dp) FUNCTION fischel_sparks(self, ion, e_kT, ne)
    class(level_type), intent(in) :: self
    real(dp), intent(in) :: e_kT, ne
    integer, intent(in) :: ion
    real(dp) :: a, a2, r, r2, r3, axl2, p, p2, p3, p4
    real(dp) :: bk_err
    !---------------------------------------------------------------------------
    bk_err = BK / uop_BK ! keep around for consistency

    a = 31.321*ion**2 * log10e * e_kT
    p = ion * 4.2E3 / (ne*bk_err)**0.166666667

    if (p <= self%xl) then
      axl2 = a / self%xl**2
      ! FISCHEL AND SPARKS, EQ. (27)
      fischel_sparks = p**4 / self%xl * (1 + axl2*(0.33333333 + 0.1*axl2))
    else
      r = self%xl - 1
      r2 = r*r
      r3 = r2*r
      p2 = p*p
      p3 = p2*p
      p4 = p3*p
      a2 = a*a

      ! FISCHEL AND SPARKS, EQ. (26)
      fischel_sparks = 1.33333333*p3 + 0.5*p2 + 0.16666667*p + 1.33333333*a*p - 0.4*a2/p &
        - 0.33333333*r3 - 0.5*r2 - 0.16666667*r - a*r + 0.5*a2/self%xl
    endif
  end FUNCTION fischel_sparks

  !===============================================================================
  PURE SUBROUTINE calc_x_and_k(ns, x1, x2, k)
    real(dp), intent(in) :: ns
    real(dp), intent(out) :: x1, x2, k
    real(dp) :: ns2, rih
    !---------------------------------------------------------------------------
    ns2 = ns*ns
    rih = ns * sqrt(2.5*ns2 + 0.5)
    x1 = (Rbohr*(rih + 1.73))**3
    x2 = (Rbohr*(rih + 1.02))**3

    k = 8.3e14 / (ns2*ns2)
    if(ns > 3.) then
      k = k * 5.33333333*ns / (ns + 1.d0)**2
    endif
    k = k*sqrt(k)
  END SUBROUTINE calc_x_and_k

  !===============================================================================
  ! Computes the occupation probability for a H atom in
  ! state with effective principle quantum number NS in a plasma
  ! enviroment with NH, NE, NHE number densities of H, ions,
  ! and He atoms respectively.  This code assumes the perturbing
  ! neutral H and He are in the ground state, (noting the hard
  ! sphere approximation used is quite crude anyway) and that ions
  ! are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
  ! See eqn 4.71 Hummer & Milhalas (1988)
  ! Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
  ! which is a fit to detailed calculation including correlations,
  ! thus the temperature dependence
  ! Sizes of neutral atoms adopted are sqrt(<r^2>)
  ! Coded by Paul Barklem and Kjell Eriksson, Aug 2003
  ! Optimized by Richard Hoppe, Jul 2023
  !===============================================================================
  PURE SUBROUTINE wcalc(nh, nhe, ne, shield, ns, w)
    real(dp), intent(in) , dimension(:,:,:) :: nh, nhe, ne, shield
    real(dp), intent(in) :: ns
    real(dp), intent(out), dimension(:,:,:) :: w
    real(dp) :: x1, x2, neutr, k1, k2, x, frac, charg, res
    real(dp), parameter :: pi43=4.18879d0 ! 4*pi/3
    integer :: i,j,k, m(3)
    !---------------------------------------------------------------------------
    call calc_x_and_k(ns, x1, x2, k1)
    k2 = k1*k1

    m = shape(ne)
    do k=1, m(3)
      do j=1, m(2)
        do i=1, m(1)
          ! neutral perturbers
          neutr = pi43 * (nh(i,j,k)*x1 + nhe(i,j,k)*x2)
          res = exp(-neutr)

          ! charged perturbers
          x = shield(i,j,k)
          if (k1*x < 1e4) then
            frac = 0.1402d0*x*k2 / (ne(i,j,k) * (1. + 0.1285*x*k1))
            charg = frac / (1.d0 + frac)
            res = res * charg
          endif

          w(i,j,k) = res
        end do
      end do
    end do
  END SUBROUTINE wcalc

  !===============================================================================
  SUBROUTINE wcalc2(self, chi, lam, ab, nh, nhe, ne, shield, wocc, nlvl)
    class(level_type), target :: self
    real(4), intent(inout), dimension(:,:,:,:) :: chi
    real(dp), intent(in), dimension(:) :: lam
    real(4), intent(in), dimension(:) :: ab
    real(dp), intent(in) , dimension(:,:,:) :: nh, nhe, ne, shield
    real(dp), intent(in), dimension(:,:,:), optional, target :: wocc, nlvl
    real(dp), dimension(:,:,:), pointer :: wp, np
    real(dp), dimension(:), allocatable :: kd
    real(4), dimension(:), allocatable :: x1, x2, kk, g
    integer, dimension(:), allocatable :: f0, f1
    real(dp) :: ns, redcut, nu, cutfac, x1d, x2d, kkd, xd
    real(4) :: w, omw, n, x, xk, k1, nh43, nhe43, ne1, frac, a, E, omE, C, omC, d
    real(dp), parameter :: pi43=4.18879d0 ! 4*pi/3
    integer :: f,i,j,k,l, m(3), nf, n0, n1
    !---------------------------------------------------------------------------
    if (present(wocc)) then
      wp => wocc
    else
      wp => self%w
    endif
    if (present(nlvl)) then
      np => nlvl
    else
      np => self%n
    endif
    m = shape(ne)
    nf = size(lam)
    if (size(chi,2) /= m(1) .or. size(wp,1) /= m(1) .or. size(np,1) /= m(1) .or. &
        size(nh,1) /= m(1) .or. size(nhe,1) /= m(1) .or. size(shield,1) /= m(1)) &
      error stop 'wcalc2: chi/nh/nhe/shield/w/n span different cells than ne'

    ! wavelength of alpha line
    redcut = hc_e / ((self%i**(-2.d0) - (self%i+1)**(-2.d0))*RYD) * 1d8

    ! Wavelengths within the cutoff split into f0, beyond the series limit, and
    ! f1, below it, where the photon reaches the dissolved upper levels of
    ! effective quantum number ns and their occupation probability wstar enters
    allocate(g(nf), f0(nf), f1(nf), x1(nf), x2(nf), kk(nf), kd(nf))
    n0 = 0
    n1 = 0
    do f=1, nf
      if (lam(f) > redcut) then
        cutfac = exp((redcut - lam(f))*1.d-2)
      else
        cutfac = 1.0
      endif
      if (cutfac > 1d-1) then
        g(f) = real(2.0898E-26 / (2*self%i**5) * cutfac * ab(f), 4)

        nu = cc/lam(f)*1e8
        ns = 1. - (self%xion + h_e*nu)/RYD
        if(ns > 0) then
          ns = 1.d0/sqrt(ns)
          call calc_x_and_k(ns, x1d, x2d, kkd)
          n1 = n1 + 1
          f1(n1) = f
          x1(n1) = real(x1d, 4)
          x2(n1) = real(x2d, 4)
          kk(n1) = real(kkd, 4)
          kd(n1) = kkd
        else
          n0 = n0 + 1
          f0(n0) = f
        endif
      endif
    end do

    do k=1, m(3)
      do j=1, m(2)
        do i=1, m(1)
          w = real(wp(i,j,k), 4)
          omw = real(1.d0 - wp(i,j,k), 4)
          xd = shield(i,j,k)
          x = real(xd, 4)
          n = real(np(i,j,k), 4)
          ne1 = real(ne(i,j,k), 4)
          nh43  = real(pi43 * nh(i,j,k), 4)
          nhe43 = real(pi43 * nhe(i,j,k), 4)

          ! n does not contain the occupation probability, therefore the
          ! (1 - wstar/wi) from the original hbop code becomes (wi - wstar)
          ! which is more efficient as we avoid division for a frequency dependent quantity
          do l=1, n0
            f = f0(l)
            chi(f,i,j,k) = chi(f,i,j,k) + w * n * g(f)
          end do
          do l=1, n1
            f = f1(l)
            ! neutral perturbers, E = exp(-a), with 1 - E formed without
            ! cancellation
            a = nh43*x1(l) + nhe43*x2(l)
            call exp_neg0(a, E, omE)

            ! charged perturbers, C = frac/(1+frac); k1 reaches ~2e22, so k1**2
            ! is formed as xk*k1 with xk = x*k1 < 1e4 to stay in the real(4) range.
            ! C is discontinuous at xk = 1e4, so that test is made on the double
            ! precision product, where a real(4) rounding could flip it
            k1 = kk(l)
            xk = x*k1
            C = 1.
            omC = 0.
            if (xd*kd(l) < 1d4) then
              frac = 0.1402*xk*k1 / (ne1 * (1. + 0.1285*xk))
              C = frac / (1. + frac)
              omC = 1. / (1. + frac)
            endif

            ! wstar = E*C.  Where w is close to 1, so is wstar near the cutoff,
            ! and w - wstar is formed as (1 - wstar) - (1 - w) with
            ! 1 - wstar = (1-E) + E*(1-C)
            if (w > 0.5) then
              d = (omE + E*omC) - omw
            else
              d = w - E*C
            endif
            chi(f,i,j,k) = chi(f,i,j,k) + d * n * g(f)
          end do
        end do
      end do
    end do
  END SUBROUTINE wcalc2

  include 'exp_neg.inc'

END MODULE level_mod
