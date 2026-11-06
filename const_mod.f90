! *******************************************************************
!>
!> physical and mathematical constants
!>
! *******************************************************************
MODULE const_mod
  USE, intrinsic :: iso_fortran_env, only: REAL64
  implicit none

  ! kinds
  integer, parameter :: dp = REAL64 !< double-precision real kind

  ! math
  real(8), parameter :: pi = 3.14159265359d0       !< Pi
  real(8), parameter :: euler_const=0.5772156649d0 !< Euler's constant
  real(8), parameter :: sqrtpi   =1.772453850906d0 !< sqrt(pi)
  real(8), parameter :: pi4 = 4.d0*pi              !< 4 pi
  complex(8), parameter :: ci = (0,1)              !< imaginary unit

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
  real(8), parameter :: ln2 = 0.69314718056d0     !< ln 2
  real(8), parameter :: ln10 = log(10.d0)         !< ln 10
  real(8), parameter :: iln2 = 1.44269504089d0    !< 1/ln(2)
  real(8), parameter :: log10e = log10(exp(1.d0)) !< log10(e)
  real(8), parameter :: inf8 = huge(1.d0)         !< largest real(8)
  real(4), parameter :: inf4 = huge(1.e0)         !< largest real(4)

  ! physics
  real(8), parameter :: ee = 1.602189d-12                              !< 1 eV in erg [erg eV^-1]
  real(8), parameter :: hh = 6.626176d-27                              !< Planck's constant [erg s]
  real(8), parameter :: cc = 2.99792458d10                             !< Speed of light [cm/s]
  real(8), parameter :: em = 9.109534d-28                              !< Electron mass [g]
  real(8), parameter :: uu = 1.6605655d-24                             !< Atomic mass unit [g]
  real(8), parameter :: bk = 1.380662d-16                              !< Boltzman's cst. [erg/K]
  real(8), parameter :: ec = 4.80325d-10                               !< electron charge [statcoulomb]
  real(8), parameter :: mh = 1.6726D-24                                !< hydrogen atom mass [g]
  real(8), parameter :: Eryd = 13.605d0*ee                             !< Rydberg constant [erg]
  real(8), parameter :: Rbohr= 5.29177349d-9                           !< Bohr radius [cm]
  real(8), parameter :: sigma_th= 8.0d0/3.0d0*pi/em**2*(ec/cc)**4      !< Thomson cross-section
  real(8), parameter :: saha_const=0.5*(hh/sqrt(2.*pi*em)/sqrt(bk))**3 !< constants in saha equation
  real(8), parameter :: XMH = 1.67339d-24                              !< uu*1.008

  ! unit conversions
  real(8), parameter :: KM_TO_CM = 1d5
  real(8), parameter :: CM_TO_KM = 1d-5
  real(8), parameter :: AA_TO_CM = 1d-8
  real(8), parameter :: CM_TO_AA = 1d8
  real(8), parameter :: HZ_TO_AA = CC / AA_TO_CM !< c [AA s^-1]: lambda [AA] = HZ_TO_AA / nu [Hz]
  real(8), parameter :: AA_TO_HZ = AA_TO_CM / CC !< 1/c [s AA^-1]: nu [Hz] = 1 / (AA_TO_HZ lambda [AA])

  ! convenient compound constants

  real(8), parameter :: hc      = hh * cc            !< hc [erg cm]
  real(8), parameter :: hc2     = hh * cc**2         !< hc^2 [erg cm^2 s^-1]
  real(8), parameter :: h_e     = hh / ee            !< h/e [eV s]
  real(8), parameter :: hc_e    = hh * cc / ee       !< hc/e [eV cm]
  real(8), parameter :: twoh_c2 = 2.0d0 * hh / cc**2 !< 2h/c^2 [erg s^3 cm^-2]
  real(8), parameter :: c2_2h   = cc**2 / 2.0d0 / hh !< c^2/2h [cm^2 s^-3 erg^-1]
  real(8), parameter :: twohc   = 2.0d0 * hh * cc    !< 2hc [erg cm]
  real(8), parameter :: hc_k    = hh * cc / bk       !< hc/k [cm K]
  real(8), parameter :: e_k     = ee/bk              !< e/k [K eV^-1]
  real(8), parameter :: h_k     = hh / bk            !< h/k [K s]
  real(8), parameter :: over_4pi= 1.0 / pi4          !< 1/(4 pi)
  real(8), parameter :: h_4pi   = hh * over_4pi      !< h/(4 pi) [erg s]
  ! 24.09.21 PE: Include parameter for line opacity
  real(8), parameter :: pie2_mc = pi * ec**2 / em / cc            !< pi e^2/(m_e c) [cm^2 s^-1]
  real(8), parameter :: RydV = 13.605d0                           !< Rydberg constant [eV]
  real(8), parameter :: Rydk = 13.605d0*e_k                       !< Rydberg constant [K]
  real(8), parameter :: RYD = 13.598434d0                         !< hydrogen ionisation energy, Rydberg energy of H [eV]
  real(8), parameter :: temp_to_theta = log10e*e_k                !< theta = temp_to_theta/T, log10(e) e/k [K]
  real(8), parameter :: op_saha_const = (2*pi*em*bk/hh**2)**1.5d0 !< (2 pi m_e k/h^2)^(3/2) [cm^-3 K^-3/2]

contains

  !===============================================================================
  !> Planck function B_nu = 2 h nu^3/c^2 / (exp(h nu/kT) - 1)
  !> [erg s^-1 cm^-2 Hz^-1 sr^-1].  Beyond h nu/kT = 80 the -1 is below rounding
  !> and the Wien form exp(-h nu/kT) is used, which cannot overflow.
  !===============================================================================
  real(8) elemental FUNCTION planck(nu,t) result(bnu)
    real(8), intent(in)  :: nu !< frequency [Hz]
    real(8), intent(in)  :: t  !< temperature [K]
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

