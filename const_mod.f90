! *******************************************************************
!>
!> physical and mathematical constants
!>
! *******************************************************************
MODULE const_mod
  implicit none

  ! math
  real(8), parameter :: pi = 3.14159265359d0 !< Pi
  real(8), parameter :: euler_const=0.5772156649d0 !< Euler's constant
  real(8), parameter :: sqrtpi   =1.772453850906d0
  real(8), parameter :: pi4 = 4.d0*pi
  complex(8), parameter :: ci = (0,1)

  ! numbers
  real(8), parameter :: zero=0.0d0
  real(8), parameter :: one=1.0d0
  real(8), parameter :: half=0.5d0
  real(8), parameter :: one_twelfth=1.0d0/12.0d0
  real(8), parameter :: one_sixth=1.0d0/6.0d0
  real(8), parameter :: one_third=1.0d0/3.0d0
  real(8), parameter :: two_third=2.0d0/3.0d0
  real(8), parameter :: one_eigth=1.0d0/8.0d0
  real(8), parameter :: one_24th =1.0d0/2.4d1
  real(8), parameter :: ln2 = 0.69314718056d0
  real(8), parameter :: ln10 = log(10.d0)
  real(8), parameter :: iln2 = 1.44269504089d0 !< 1/ln(2)

  ! physics
  real(8), parameter :: ee = 1.602189d-12  !< electron volt [erg]
  real(8), parameter :: hh = 6.626176d-27  !< Planck's constant [erg s]
  real(8), parameter :: cc = 2.99792458d10 !< Speed of light [cm/s]
  real(8), parameter :: em = 9.109534d-28  !< Electron mass [g]
  real(8), parameter :: uu = 1.6605655d-24 !< Atomic mass unit [g]
  real(8), parameter :: bk = 1.380662d-16  !< Boltzman's cst. [erg/K]
  real(8), parameter :: ec = 4.80325d-10   !< electron charge [statcoulomb]
  real(8), parameter :: mh = 1.6726D-24    !< hydrogen atom mass [g]
  real(8), parameter :: Eryd = 13.605d0*ee !< Rydberg constant [erg]
  real(8), parameter :: Rbohr= 5.29177349d-9 !< Bohr radius [cm]
  real(8), parameter :: sigma_th= 8.0d0/3.0d0*pi/em**2*(ec/cc)**4 !< Thomson cross-section
  real(8), parameter :: saha_const=0.5*(hh/sqrt(2.*pi*em)/sqrt(bk))**3 !< constants in saha equation
  real(8), parameter :: XMH = 1.67339d-24 !< uu*1.008

  ! unit conversions
  real(8), parameter :: KM_TO_CM = 1d5
  real(8), parameter :: CM_TO_KM = 1d-5
  real(8), parameter :: AA_TO_CM = 1d-8
  real(8), parameter :: CM_TO_AA = 1d8
  real(8), parameter :: HZ_TO_AA = CC / AA_TO_CM
  real(8), parameter :: AA_TO_HZ = AA_TO_CM / CC

  ! convenient compound constants

  real(8), parameter :: hc      = hh * cc
  real(8), parameter :: hc2     = hh * cc**2
  real(8), parameter :: h_e     = hh / ee
  real(8), parameter :: hc_e    = hh * cc / ee
  real(8), parameter :: twoh_c2 = 2.0d0 * hh / cc**2
  real(8), parameter :: c2_2h   = cc**2 / 2.0d0 / hh
  real(8), parameter :: twohc   = 2.0d0 * hh * cc
  real(8), parameter :: hc_k    = hh * cc / bk
  real(8), parameter :: e_k     = ee/bk
  real(8), parameter :: h_k     = hh / bk
  real(8), parameter :: over_4pi= 1.0 / pi4
  real(8), parameter :: h_4pi   = hh * over_4pi
  ! 24.09.21 PE: Include parameter for line opacity
  real(8), parameter :: pie2_mc = pi * ec**2 / em / cc
  real(8), parameter :: RydV = 13.605d0     !< Rydberg constant [eV]
  real(8), parameter :: Rydk = 13.605d0*e_k !< Rydberg constant [K]

contains

  !===============================================================================
  ! PE: 22-09-21: Convert to elemental FUNCTION
  !>  calculates planck FUNCTION bnu at frequency nu, temp t
  !===============================================================================
  real(8) elemental FUNCTION planck(nu,t) result(bnu)
    real(8), intent(in)  :: nu
    real(8), intent(in)  :: t
    real(8) :: x, c1
    !---------------------------------------------------------------------
    c1 = twoh_c2 * nu**3
    x = h_k*nu/t

    if(x.lt.80.) then
      bnu = c1 / (exp(x) - 1.d0)
    else
      bnu = c1 * exp(-x)
    endif
  end FUNCTION planck
END MODULE const_mod

