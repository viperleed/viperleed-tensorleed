! Commit one runtime fault, selected by the first command-line argument:
!   bounds    read one element past the end of an array
!   invalid   square root of a negative number
!   zero      division by zero
!   overflow  product larger than huge()
! Built with the debug flags (harness.build.DEBUG), each must stop the
! program with a runtime error or SIGFPE. The operand z is 0, but it is
! taken from the number of arguments, so the compiler cannot fold the
! faulty operation away; the probe prints its result if it survives.
program trap_probe
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  character(len=16) :: mode
  real(real64) :: a(4), z, r
  integer :: i

  call get_command_argument(1, mode)
  z = real(command_argument_count() - 1, real64)
  a = 1.0_real64
  select case (mode)
  case ('bounds')
    i = size(a) + 1 + int(z)
    r = a(i)
  case ('invalid')
    r = sqrt(z - 1.0_real64)
  case ('zero')
    r = 1.0_real64 / z
  case ('overflow')
    r = huge(z) * (z + 2.0_real64)
  case default
    error stop 'unknown mode'
  end select
  print *, r
end program trap_probe
