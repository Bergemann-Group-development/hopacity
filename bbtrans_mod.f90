MODULE bbtrans_mod
  USE prec_mod
  USE const_mod
  USE op_math_mod, only: searchsorted, find_closest
  USE hop_env_mod, only: trace
  USE handyfuncs_mod, only: find_closest_from_center
  implicit none

  type :: bbtrans_type
    real(dp) :: nu0
    real(dp) :: lam0
    real(dp) :: lam_air
    integer :: nnu ! number of frequencies in the line
    real(dp) :: ga=0, gw1=0, gw2=0, gw3=0, gq=0, gv=0 ! line broadening (damping) constants
    real(4) :: vdw_ratio=0 ! Lindholm vdW shift/HWHM ratio d/w (see prepare_VanderWaals)
    real(dp) :: ga_in=0, gq_in=0, gv_in=0 ! input values for line broadening
    real(dp) :: lower_ev, upper_g=1, loggf
    real(dp) :: Aji, Aij, Bji, Bij ! Einstein coefficients, Bij is the
    ! transition from the lower to the upper
    ! level
    logical :: active=.true., unsold
    logical, pointer :: used => null()
    real(dp), pointer :: strength => null()
    integer :: iso_idx=0 ! isotope index for lines in line list
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
    procedure :: get_center
    procedure :: deactivate
  end type

  integer, save :: count1, count2, count3, count4, count5, count6
  logical, public :: disable_stark=.false., disable_vdw=.false., disable_ga=.false.
  logical, public :: ga_fallback=.false.
  ! Lindholm van der Waals pressure line shift: a global scale on the per-line
  ! impact-theory ratio d/w = tan(pi*alpha/2), built from each line's width velocity
  ! parameter alpha by prepare_VanderWaals and stored per line (vdw_ratio here,
  ! ll%vdw_ratio in linelist_mod). vdw_shift=1.0 applies the full physical (red-shifted)
  ! ratio; 0 = neglected (default, reproduces the pre-shift synthesis).
  real(4), public :: vdw_shift=0.0
  real(4), public, save :: tolerance=1e-4
  real(4), public, save :: line_reject_frac=0.0
  real(4), parameter :: pqs  = 1.128379167  ! 2./pi**0.5
  real(4), parameter :: isqp = 0.5641895835 ! 1/pi**0.5
  real(4), parameter :: ipi  = 0.3183098862 ! 1/pi
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
  SUBROUTINE add_opacities(self, f0, nu, vproj, inv_df, e_kT, dxij, inv_part, &
    ne, tothi, vbar2, stark_pe, unsold, bi, bj, chi, eta, cont, do_nlte)
    type(bbtrans_type), intent(inout) :: self
    integer, intent(in) :: f0
    real(8), dimension(:), intent(in), contiguous :: nu
    real(4), dimension(:,:,:), intent(in), contiguous :: vproj, inv_df, inv_part
    real(4), dimension(:,:,:), intent(in), pointer :: bi, bj
    real(4), dimension(:,:,:), intent(in), contiguous :: e_kT, dxij, ne, tothi, vbar2, stark_pe, unsold
    real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
    real(4), dimension(:,:,:,:), intent(in) :: cont
    logical, intent(in) :: do_nlte
    real(4) :: ex, ni, vs, phi, lchi, leta, ceta, cchi, Rij
    real(4) :: a, idf, idf_nu, v, ref, ekt, inv_nu0, strength
    integer :: i,j,k,f, nf, m(3), itimer=0
    logical :: used
    !----------------------------------------------------------------------------------
    !call trace%begin('add_lte_opacities_1', itimer=itimer)
    m = shape(vproj)
    nf = size(nu)
    inv_nu0 = 1.0 / self%nu0
    used = .false.
    strength = 0.0

    do k=1, m(3)
      do j=1, m(2)
        do i=1, m(1)
          ! Fixed Broadening
          a = self%ga
          ! Van der Waals Broadening
          if (self%unsold) then
            a = a + self%gw1 * unsold(i,j,k)
          else
            a = a + self%gw1 * vbar2(i,j,k)**self%gw2 * tothi(i,j,k)
          endif
          ! Stark Broadening
          if (self%gq >= 0.0) then
            a = a + self%gq * ne(i,j,k)
          else! if line%gq < 0
            a = a - self%gq * stark_pe(i,j,k)
          end if

          ! fetch variables
          idf = inv_df(i,j,k)
          idf_nu = idf * inv_nu0
          v   = (vproj(i,j,k) + 1)*self%nu0
          ekt = e_kT(i,j,k)
          ex  = (self%lower_ev - dxij(i,j,k)) * ekt

          if (ex < 85.) then
            a = a * idf
            ni = inv_part(i,j,k) * exp(-ex)
            ! ni is really ni/gi, the factor gi is contained in Bij
            ! since we have used loggf to compute Aij and hence Bij

            Rij = ni * self%Bij * idf
            if (do_nlte) then
              ceta = Rij * bj(i,j,k)
              cchi = Rij * bi(i,j,k)
            else
              ceta = Rij
              cchi = Rij
            endif

            ! check if line peak is larger than continuum tolerance
            ref = tolerance * cont(f0,i,j,k)

            vs = (nu(f0) - v) * idf_nu
            phi = voigt_function2(a, vs)
            lchi = phi * cchi
            leta = phi * ceta
            if (lchi > ref) then
              used = .true.
              strength = max(strength, lchi / cont(f0,i,j,k))

              chi(f0,i,j,k) = chi(f0,i,j,k) + lchi
              eta(f0,i,j,k) = eta(f0,i,j,k) + leta

              f = f0+1
              do while (lchi >= ref .and. f <= nf) ! if line opacity below tolerance exit
                vs = (nu(f) - v) * idf_nu
                phi = voigt_function2(a, vs)
                lchi = phi * cchi
                leta = phi * ceta
                chi(f,i,j,k) = chi(f,i,j,k) + lchi
                eta(f,i,j,k) = eta(f,i,j,k) + leta
                f = f + 1
              end do

              f = f0-1
              lchi = ref
              do while (lchi >= ref .and. f >= 1) ! if line opacity below tolerance exit
                vs = (nu(f) - v) * idf_nu
                phi = voigt_function2(a, vs)
                lchi = phi * cchi
                leta = phi * ceta
                chi(f,i,j,k) = chi(f,i,j,k) + lchi
                eta(f,i,j,k) = eta(f,i,j,k) + leta
                f = f - 1
              end do
            endif
          endif
        end do
      end do
    end do
    if (used) then
      self%used = used
      self%strength = max(self%strength, strength)
    endif
    !call trace%end(itimer)
  END SUBROUTINE add_opacities
  !------------------------------------------------------------------------------------

  !====================================================================================
  !DIR$ ATTRIBUTES FORCEINLINE :: voigt_function2
  ELEMENTAL FUNCTION voigt_function2(as,vs) result(h)
    !---------------------------------------------------------------------
    !  computes a voigt function  h = h(a,v) / sqrt(pi)
    !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
    !  traving (landolt-b\rnstein, p. 449).
    !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
    !---------------------------------------------------------------------
    real(4), intent(in) :: as,vs
    real(4) :: h, a
    real(4) :: v,u,v2,v4,ex
    real(4) :: a2,u2
    complex(8) :: tav,uav,w4,num,den
    !----------------------------------------------------------------------------------
    v = abs(vs)
    a = as
    u = v + a
    v2 = v*v

    if (a > 0.2) then
      if (a > 1.4 .or. u > 3.2) then
        a2 = a*a
        u = 1 / (a2 + v2)
        u2 = 0.5 * u*u
        ! a gt 1.4  or  a + v gt 3.2
        ! make sure to never compute a2**2 as this may overflow
        h = ipi*a*u * (1. + u2*(7.5 - a2 + 3.*v2 - 12*u2*a2*a2))

      else
        tav = cmplx(a,-v)
        if (a >= 0.195e0*v - 0.176e0) then
          num = 16.4955e0 + tav*(20.20933e0 + tav*(11.96482e0      &
            + tav*(3.778987e0 + tav*0.5642236e0)))
          den = 16.4955e0 + tav*(38.82363e0 + tav*(39.27121e0 &
            + tav*(21.69274e0 + tav*(6.699398e0 + tav))))
          w4 = num/den
        else
          uav = tav*tav
          num = tav*(36183.31e0 - uav*(3321.9905e0 - uav*(1540.787e0 &
            - uav*(219.0313e0 - uav*( 35.76683e0 - uav*(1.320522e0 &
            - uav*0.56419e0))))))
          den = (32066.6e0 - uav*(24322.84e0 - uav*(9022.228e0       &
            - uav*(2186.181e0 - uav*(364.2191e0                    &
            - uav*(61.57037e0 - uav*(1.841439e0 - uav)))))))
          w4 = num/den
          w4 = exp(uav) - w4
        endif
        h = real(w4) * isqp
      endif

    elseif (v >= 5.0) then
      v4 = 4*v2*v2
      h = ipi * a*(15. + 6.*v2 + v4) / (v4*v2)

    else
      if (v >= 2.4) then
        h = 1./(v2 - 1.5) * (0.554153432 + v*(0.278711796 + v*(-0.188325687 + v*(0.042991293 + v*-0.003278278))))
      elseif (v >= 1.3) then
        h = -4.48480194 + v*(9.39456063 + v*(-6.61487486 + v*(1.98919585 + v*-0.22041650)))
      else
        h = -1.12470432 + v*(-0.15516677 + v*(3.28867591 + v*(-2.34357915 + v*0.42139162)))
      endif

      ! a le 0.2  and v lt 5.
      ex = exp(-v2)
      h = (h*a + ex*(1. + a*a*(1. - 2.*v2))) * isqp
    endif
  END FUNCTION voigt_function2
  !------------------------------------------------------------------------------------

  !====================================================================================
  !DIR$ ATTRIBUTES FORCEINLINE :: voigt_profile
  ELEMENTAL FUNCTION voigt_profile(a,q,vpr,idf) result(h)
    !---------------------------------------------------------------------
    !  computes a normalised voigt profile  h = h(a,v) / sqrt(pi) * nu/dnud
    !  a=gamma/(4*pi*dnud)   and  q=(nu-nu0)/nu0.  this  is  done after
    !---------------------------------------------------------------------
    real(4), intent(in) :: a,q,vpr,idf
    real(4) :: h
    real(4) :: v,u,v2,v4,ex
    real(4) :: a2,u2
    complex(8) :: tav,uav,w4,num,den
    !----------------------------------------------------------------------------------
    v = abs(q - vpr) * idf
    u = a + v
    v2 = v*v

    if (a > 0.2) then
      if (a > 1.4 .or. u > 3.2) then
        a2 = a*a
        u = 1 / (a2 + v2)
        u2 = 0.5 * u*u
        ! a gt 1.4  or  a + v gt 3.2
        h = ipi*a*u * (1. + u2*((3.*v2 - a2) + u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2)))

      else
        tav = cmplx(a,-v)
        if (a >= 0.195e0*v - 0.176e0) then
          num = 16.4955e0 + tav*(20.20933e0 + tav*(11.96482e0      &
            + tav*(3.778987e0 + tav*0.5642236e0)))
          den = 16.4955e0 + tav*(38.82363e0 + tav*(39.27121e0 &
            + tav*(21.69274e0 + tav*(6.699398e0 + tav))))
          w4 = num/den
        else
          uav = tav*tav
          num = tav*(36183.31e0 - uav*(3321.9905e0 - uav*(1540.787e0 &
            - uav*(219.0313e0 - uav*( 35.76683e0 - uav*(1.320522e0 &
            - uav*0.56419e0))))))
          den = (32066.6e0 - uav*(24322.84e0 - uav*(9022.228e0       &
            - uav*(2186.181e0 - uav*(364.2191e0                    &
            - uav*(61.57037e0 - uav*(1.841439e0 - uav)))))))
          w4 = num/den
          w4 = exp(uav) - w4
        endif
        h = real(w4) * isqp
      endif

    elseif (v >= 5.0) then
      v4 = 4*v2*v2
      h = ipi * a*(15. + 6.*v2 + v4) / (v4*v2)

    else
      if (v >= 2.4) then
        h = 1./(v2 - 1.5) * (0.554153432 + v*(0.278711796 + v*(-0.188325687 + v*(0.042991293 + v*-0.003278278))))
      elseif (v >= 1.3) then
        h = -4.48480194 + v*(9.39456063 + v*(-6.61487486 + v*(1.98919585 + v*-0.22041650)))
      else
        h = -1.12470432 + v*(-0.15516677 + v*(3.28867591 + v*(-2.34357915 + v*0.42139162)))
      endif

      ! a le 0.2  and v lt 5.
      ex = exp(-v2)
      h = (h*a + ex*(1. + a*a*(1. - 2.*v2))) * isqp
    endif
    h = h * idf * (q + 1)
  END FUNCTION voigt_profile
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION voigt_function(a,vs) result(h)
    !---------------------------------------------------------------------
    !  computes a voigt function  h = h(a,v)
    !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
    !  traving (landolt-b\rnstein, p. 449).
    !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
    !---------------------------------------------------------------------
    real(4), intent(in) :: a,vs
    real(4) :: h
    integer :: k,m,i
    real(8),dimension(5) :: a1
    real(8) :: v,u,v2,ex,quo,h1,pqs,h1p,h2p,h3p,h4p,psi
    real(8) :: a2,u2

    real(8),dimension(19),parameter :: ak = (/ &
      -1.12470432, -0.15516677,  3.28867591, -2.34357915,             &
      0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,&
      -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,&
      -0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278 /)

    real(8),parameter :: sqp = 1.772453851
    real(8),parameter :: sq2 = 1.414213562
    !----------------------------------------------------------------------------------
    v = abs(vs)
    u = a + v
    v2 = v*v

    k = 1

    if (a == 0.0) then
      h=0.
      if(v2.lt.100.)h=exp(-v2)
      return

    elseif (a > 0.2) then
      if (a > 1.4 .or. u > 3.2) then
        a2 = a*a
        u = sq2*(a2 + v2)
        u2 = 1./(u*u)
        ! a gt 1.4  or  a + v gt 3.2
        h = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) + &
          u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2))
        return
      else
        k = 2
      endif
    elseif (v >= 5.0) then
      h = a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp)
      return
    endif

    ex=0.
    if(v2 < 100.) ex = exp(-v2)

    quo = 1.
    if (v >= 2.4) then
      quo = 1./(v2 - 1.5)
      m = 11
    else
      m = 6
      if (v < 1.3) m = 1
    endif

    do i=1,5
      a1(i) = ak(m)
      m = m + 1
    end do

    h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))

    if (k==1) then
      ! a le 0.2  and v lt 5.
      h = h1*a + ex*(1. + a*a*(1. - 2.*v2))
      return
    endif

    pqs = 2./sqp
    h1p = h1 + pqs*ex
    h2p = pqs*h1p - 2.*v2*ex
    h3p = (pqs*(1. - ex*(1. - 2.*v2)) - 2.*v2*h1p)/3. + pqs*h2p
    h4p = (2.*v2*v2*ex - pqs*h1p)/3. + pqs*h3p
    psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
    ! 0.2 < a <= 1.4  and  a + v <= 3.2
    h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
    return
  END FUNCTION voigt_function
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION new_voigt(a,v) result(h)
    ! Champion Voigt routine in speed and accuracy test by BE+NP 97-11-18
    !  VOIGT FUNCTION CALCULATION: HUMLICEK'S APPROXIMATION
    ! The results are normalized to an area of 1/sqrt(pi).
    real(4), intent(in) :: a, v
    complex(8) :: tav,uav,w4,v4
    real(4)    :: sav,aa,vv,h
    integer :: i,j,k
    !----------------------------------------------------------------------------------
    aa = 0.01! a
    vv = abs(v)
    tav = cmplx(aa,-vv)
    sav = vv + aa
    uav = tav*tav
    if (sav >= 15.e0) then
      w4 = tav*0.5641896e0 / (0.5e0 + uav)
    elseif (sav >= 5.5e0) then
      w4 = tav*(1.410474e0 + uav*0.5641896e0) / (0.75e0 + uav*(3.e0 + uav))
    elseif (aa >= 0.195e0*vv - 0.176e0) then
      w4 = (16.4955e0 + tav*(20.20933e0 + tav*(11.96482e0       &
        +  tav*(3.778987e0 + tav*0.5642236e0)))) / (16.4955e0  &
        +  tav*(38.82363e0 + tav*(39.27121e0 + tav*(21.69274e0 &
        +  tav*(6.699398e0 + tav)))))
    else
      w4 = tav*(36183.31e0 - uav*(3321.9905e0 - uav*(1540.787e0 &
        - uav*(219.0313e0 - uav*( 35.76683e0 - uav*(1.320522e0 &
        - uav*0.56419e0))))))
      v4 = (32066.6e0 - uav*(24322.84e0 - uav*(9022.228e0       &
        - uav*(2186.181e0 - uav*(364.2191e0                    &
        - uav*(61.57037e0 - uav*(1.841439e0 - uav)))))))
      w4 = exp(uav) - w4/v4
    endif
    h = real(w4) * isqp
  END FUNCTION new_voigt
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION dirty_voigt(a, x) result(h)
    real(4), intent(in) :: a, x
    real(4) :: h, z
    ! Coefficients for the approximation function (Table A1, Smith+15)
    real(4), parameter :: A0 = 15.75328153963877
    real(4), parameter :: A1 = 286.9341762324778
    real(4), parameter :: A2 = 19.05706700907019
    real(4), parameter :: A3 = 28.22644017233441
    real(4), parameter :: A4 = 9.526399802414186
    real(4), parameter :: A5 = 35.29217026286130
    real(4), parameter :: A6 = 0.8681020834678775
    real(4), parameter :: B0 = 0.0003300469163682737
    real(4), parameter :: B1 = 0.5403095364583999
    real(4), parameter :: B2 = 2.676724102580895
    real(4), parameter :: B3 = 12.82026082606220
    real(4), parameter :: B4 = 3.21166435627278
    real(4), parameter :: B5 = 32.032981933420
    real(4), parameter :: B6 = 9.0328158696
    real(4), parameter :: B7 = 23.7489999060
    real(4), parameter :: B8 = 1.82106170570
    !----------------------------------------------------------------------------------
    z = x * x
    if (z <= 3.0) then
      h = exp(-z) * (1.0 - a * (A0 + A1 / (z - A2 + A3 / (z - A4 + A5 / (z - A6)))))
    else if (z < 25.0) then
      h = exp(-z) + a * (B0 + B1 / (z - B2 + B3 / (z + B4 + B5 / (z - B6 + B7 / (z - B8)))))
    else
      h = 0.5 * sqrt(2.0 / pi) * a / (z - 1.5 - 1.5 / (z - 3.5 - 5.0 / (z - 5.5)))
    end if
    h = h * isqp
  END FUNCTION dirty_voigt

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
    real(dp), parameter :: scale=1e-4*one_sixth
    real(dp), parameter :: c19=log10(39.*(8./(PI*EM*BK**5))**one_sixth) ! 19.4
    real(dp) :: nu0_ev, upper_ev, c1
    !----------------------------------------------------------------------------------
    if (gq < -10) then
      ! recipe c4 SIU (formula from gray, p.237)
      ! make negative to mark this recipe
      self%gq = - BK * 10**(c19 + two_third*gq)
    elseif (gq < 0) then
      ! make negative to mark this recipe
      self%gq = - 10**gq * scale
    elseif (gq == 0) then
      nu0_ev = H_E * self%nu0
      upper_ev = self%lower_ev + nu0_ev
      if (self%lower_ev < cnt_ev) then
        c1 = (j/(cnt_ev - upper_ev))**2
        if (j == 1) then
          c1 = 2.26e-7 * c1
        else
          c1 = 5.42e-7 / (j+1)**2
        endif
        self%gq = RydV*c1
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
      ! Here we assume Aji for the radiative broadening, Turbospectrum devides this by
      ! the upper statistical weight, but KORG does not. MOOg assumes GF/GU=3 instead.
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
    real(dp), parameter :: uns=8.411*1.283984e-12 !8.411~(13.598**2)**0.4, 1.283984e-12~17*(1.61e-33)**0.4
    real(dp), parameter :: barklem = (1e-2*Rbohr)**2 ! 2.80028e-21 bohr radius in m and squared
    !----------------------------------------------------------------------------------
    ! compute van der waals broadening due to hydrogen
    ! following Barklem
    self%unsold = .false.
    self%vdw_ratio = 0
    if(gv >= 20.0) then
      sigma = barklem*int(gv)
      alpha = gv - int(gv)
      gx = 1.0 - 0.5*alpha
      gammaf = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*gx)*gx)*gx)*gx)*gx
      balfac = (4.0/pi)**(0.5*alpha) * gammaf * 2.e10

      self%gw2 = 0.5*(1 - alpha)
      self%gw1 = balfac * sigma * 10**(-12*self%gw2)
      ! Lindholm shift/HWHM ratio for an effective inverse-power interaction:
      ! d/w = tan(pi*alpha/2). The clamp guards against pathological alpha (near 1,
      ! where the tangent blows up).
      self%vdw_ratio = min(tan(0.5*pi * min(max(alpha, 0.02_dp), 0.96_dp)), 3.0_dp)
    elseif (.not. is_mol) then ! cannot use Unsold for molecules
      self%unsold = .true.
      self%vdw_ratio = 0.7265 ! pure van der Waals (p=6, alpha=0.4) -> tan(0.2*pi)
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
        ! Jul  5 2010 jorritl: I'm sceptical about the numerical
        ! constant  1.283984e-12. According to Gray it should be
        ! 0.3d-30^0.4 = 6.1780059e-13 and  8.411 should be 17
        ! I haven't found the original reference to Unsold 1955 yet...
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
      self%gw1 = 0    ! vdw_ratio is then irrelevant: it multiplies a zero vdW width
    endif
  END SUBROUTINE prepare_VanderWaals
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE get_center(self, lam, dlam, i0)
    class(bbtrans_type), intent(in)  :: self
    real(dp), dimension(:), intent(in) :: lam
    real(dp), intent(in) :: dlam
    integer, intent(out) :: i0
    !----------------------------------------------------------------------------------
    if (dlam == 0) then
      i0 = findloc(lam, self%lam0, 1)
    else
      i0 = find_closest(lam, self%lam0)
      if (abs(lam(i0) - self%lam0) > dlam) then
        i0  = 0
      endif
    endif
  END SUBROUTINE get_center
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
