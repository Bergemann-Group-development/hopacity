MODULE prec_mod
  USE, intrinsic :: iso_fortran_env
  USE const_mod
  implicit none

  !integer, parameter :: sp = REAL32 ! in conflict with spectrum object
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128
  real(dp), parameter :: log10e = log10(exp(1.d0))

  ! physics
  real(dp), parameter :: EK  = EE/BK
  !real(dp), parameter :: op_saha_const = (2*PI*EM*BK/HH**2)**1.5
  real(dp), parameter :: op_saha_const = 1.202e9/BK/2/5040**2.5
  !real(dp), parameter :: temp_to_theta = log10e*EK
  real(dp), parameter :: temp_to_theta = 5040
  real(dp), parameter :: uop_BK = 1.38053e-16
  real(8),  parameter :: inf8 = huge(1.d0)
  real(4),  parameter :: inf4 = huge(1.e0)
  !real(dp), parameter :: uop_BK = BK

  ! This is public so it can be overwritten by input from absdat
  real(dp), public :: RYD = 13.598408 !< Rydberg constant [eV]

  integer, public, parameter :: pp(4)=[1,1,1,1]

contains

END MODULE prec_mod
