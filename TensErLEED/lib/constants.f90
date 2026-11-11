!> \file constants.f90
!! \brief Physical and numerical constants of ref-calc.
!!
!! The values are those of the legacy file src/GLOBAL and of the literal
!! 3.14159265 that the legacy routines use for pi
!! (tests/unit/test_constants.py).

!> Physical and numerical constants, in the precision of the legacy code.
module constants
  use kinds, only: sp
  implicit none

  private

  public :: bohr, hartree, memach, pi

  !> Bohr radius in Angstrom (BOHR of the legacy file GLOBAL).
  real(sp), parameter :: bohr = 0.529177_sp
  !> Hartree in eV (HARTREE of the legacy file GLOBAL).
  real(sp), parameter :: hartree = 27.211396_sp
  !> Machine accuracy used in matrix inversion routines (MEMACH of GLOBAL).
  real(sp), parameter :: memach = 1.0e-6_sp
  !> Pi, as the legacy literal 3.14159265.
  real(sp), parameter :: pi = 3.14159265_sp
end module constants
