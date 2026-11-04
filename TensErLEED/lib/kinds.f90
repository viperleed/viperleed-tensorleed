!> \file kinds.f90
!! \brief Kind parameters of ref-calc.
!!
!! lib/rfacsb.f90 holds an identical module kinds for rfactor; the two must
!! stay identical (tests/unit/test_kinds.py).

!> Kind parameters for reals and complexes.
!!
!! sp is the precision of the legacy default REAL (REAL(4)); dp is double
!! precision (REAL(8)). Every real and complex is declared with one of these,
!! so that a change of precision is a change of one line here.
module kinds
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none

  private

  public :: sp, dp

  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
end module kinds
