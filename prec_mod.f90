MODULE prec_mod
  USE, intrinsic :: iso_fortran_env
  USE const_mod
  implicit none

  !integer, parameter :: sp = REAL32 ! in conflict with spectrum object
  integer, parameter :: dp = REAL64                !< double-precision real kind
  integer, parameter :: qp = REAL128               !< quadruple-precision real kind
  real(dp), parameter :: log10e = log10(exp(1.d0)) !< log10(e)

  ! physics
  real(dp), parameter :: EK  = EE/BK  !< e/k [K eV^-1]
  !real(dp), parameter :: op_saha_const = (2*PI*EM*BK/HH**2)**1.5
  real(dp), parameter :: op_saha_const = 1.202e9/BK/2/5040**2.5  !< (2 pi m_e k/h^2)^(3/2) [cm^-3 K^-3/2] of the Saha equation
  !real(dp), parameter :: temp_to_theta = log10e*EK
  real(dp), parameter :: temp_to_theta = 5040 !< theta = temp_to_theta/T, the rounded log10(e) e/k [K]
  real(dp), parameter :: uop_BK = 1.38053e-16 !< Boltzmann constant [erg K^-1], the older value fischel_sparks is written for
  real(8),  parameter :: inf8 = huge(1.d0)    !< largest real(8)
  real(4),  parameter :: inf4 = huge(1.e0)    !< largest real(4)
  !real(dp), parameter :: uop_BK = BK

  ! This is public so it can be overwritten by input from absdat
  real(dp), public :: RYD = 13.598408 !< Rydberg constant [eV]

  integer, public, parameter :: pp(4)=[1,1,1,1]

contains

END MODULE prec_mod
