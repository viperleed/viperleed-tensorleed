!> \file leed_output.f90
!! \brief Writer of ref-calc: messages of the input reader.
!!
!! Every routine takes the unit to write to; the unit is opened by the
!! caller.

!> Output of ref-calc.
module leed_output
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  private

  public :: stop_with_error
  public :: write_no_bulk_calculation, write_tensor_not_required
  public :: write_tensor_required, write_tensor_file

contains

  !> Write message on stderr and stop with a non-zero exit code.
  subroutine stop_with_error(message)
    character(len=*), intent(in) :: message

    write (error_unit, '(a)') 'ref-calc: ' // message
    error stop 1
  end subroutine stop_with_error

  !> Echo that the input requires no bulk calculation.
  subroutine write_no_bulk_calculation(unit)
    integer, intent(in) :: unit

    write (unit, *) "No bulk calculation by subref required!"
  end subroutine write_no_bulk_calculation

  !> Echo that no tensor output is required for a stacked layer.
  subroutine write_tensor_not_required(unit, istack)
    integer, intent(in) :: unit
    integer, intent(in) :: istack  !< stacked layer number

    write (unit, *) "Tensor output for layer no. ", istack, &
                    " not required."
  end subroutine write_tensor_not_required

  !> Echo that tensor output is required for a stacked layer.
  subroutine write_tensor_required(unit, istack)
    integer, intent(in) :: unit
    integer, intent(in) :: istack  !< stacked layer number

    write (unit, *) "Tensor output for layer no. ", istack, &
                    " required. Filename(s): "
  end subroutine write_tensor_required

  !> Echo the tensor file of one subplane of a stacked layer.
  subroutine write_tensor_file(unit, istack, isub, filename)
    integer, intent(in) :: unit
    integer, intent(in) :: istack  !< stacked layer number
    integer, intent(in) :: isub  !< subplane number
    character(len=*), intent(in) :: filename

    write (unit, *) "Layer ", istack, ", subplane ", isub, ": ", filename
  end subroutine write_tensor_file
end module leed_output
