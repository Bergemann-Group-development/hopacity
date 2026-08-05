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
  SUBROUTINE wcalc2(self, chi, lam, ab, nh, nhe, ne, shield, i0, i1)
    class(level_type) :: self
    real(4), intent(inout), dimension(:,:,:,:) :: chi
    real(dp), intent(in), dimension(:) :: lam, ab
    real(dp), intent(in) , dimension(:,:,:) :: nh, nhe, ne, shield
    integer, intent(in), optional :: i0, i1
    real(dp), dimension(:), allocatable :: x1, x2, kk, lfac, cutfac
    real(dp) :: ns, ns2, rih, neutr, k1, k2, x, charg, res, nh43, nhe43
    real(dp) :: redcut, nu, w, n, frac, wstar
    real(dp), parameter :: pi43=4.18879d0 ! 4*pi/3
    integer :: f,i,j,k, m(3), nf, xlo, xhi
    !---------------------------------------------------------------------------
    m = shape(ne)
    nf = size(lam)
    xlo = 1; xhi = m(1)
    if (present(i0)) xlo = i0
    if (present(i1)) xhi = i1
    allocate(x1, x2, kk, lfac, cutfac, mold=lam)
    kk = 0.

    ! wavelength of alpha line
    redcut = hc_e / ((self%i**(-2.d0) - (self%i+1)**(-2.d0))*RYD) * 1d8

    do f=1, nf
      if (lam(f) > redcut) then
        cutfac(f) = exp((redcut - lam(f))*1.d-2)
      else
        cutfac(f) = 1.0
      endif
      if (cutfac(f) > 1d-1) then
        lfac(f) = 2.0898E-26 / (2*self%i**5) * cutfac(f)

        nu = cc/lam(f)*1e8
        ns = 1. - (self%xion + h_e*nu)/RYD
        if(ns > 0) then
          ns = 1.d0/sqrt(ns)
          call calc_x_and_k(ns, x1(f), x2(f), kk(f))
        endif
      endif
    end do

    do k=1, m(3)
      do j=1, m(2)
        do i=xlo, xhi
          w = self%w(i,j,k)
          x = shield(i,j,k)
          n = self%n(i,j,k)
          nh43  = pi43 * nh(i,j,k)
          nhe43 = pi43 * nhe(i,j,k)
          do f=1, nf
            if (cutfac(f) > 1d-1) then
              if(kk(f) /= 0.d0) then
                ! neutral perturbers
                neutr = nh43*x1(f) + nhe43*x2(f)
                wstar = exp(-neutr)

                ! charged perturbers
                k1 = kk(f)
                if (k1*x < 1e4) then
                  k2 = k1*k1
                  frac = 0.1402d0*x*k2 / (ne(i,j,k) * (1. + 0.1285*x*k1))
                  charg = frac / (1.d0 + frac)
                  wstar = wstar * charg
                endif

                res = w - wstar
              else
                res = w
              endif

              ! n does not contain the occupation probability, therefore the
              ! (1 - wstar/wi) from the original hbop code becomes (wi - wstar)
              ! which is more efficient as we avoid division for a frequency dependent quantity
              res = res * lfac(f) * n
              chi(f,i,j,k) = chi(f,i,j,k) + res*ab(f)
            endif
          end do
        end do
      end do
    end do
  END SUBROUTINE wcalc2

END MODULE level_mod
