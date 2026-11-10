!> \file leed_output.f90
!! \brief Writer of ref-calc: messages of the input reader.
!!
!! Every routine takes the unit to write to; the unit is opened by the
!! caller.

!> Output of ref-calc.
module leed_output
  use, intrinsic :: iso_fortran_env, only: error_unit
  use kinds, only: sp
  implicit none

  private

  public :: stop_with_error
  public :: write_no_bulk_calculation, write_tensor_not_required
  public :: write_tensor_required, write_tensor_file
  public :: write_title, write_substrate_vectors, write_superlattice_vectors
  public :: write_inner_potential, write_beam_criterion, write_output_beams
  public :: write_incidence, write_convergence, write_lmax
  public :: write_phase_shift_format, write_iform, write_phase_shift_header
  public :: write_inner_potential_header, write_phase_shift_coefficients
  public :: write_phase_shift_file_title, write_phase_shift_parameters
  public :: write_first_phase_shifts, write_end_of_input

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

  !> Echo the title of the calculation.
  subroutine write_title(unit, title)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: title

    write (unit, *) title
  end subroutine write_title

  !> Echo the substrate lattice vectors, in Angstrom.
  subroutine write_substrate_vectors(unit, a1, a2)
    integer, intent(in) :: unit
    real(sp), intent(in) :: a1(2), a2(2)

    write (unit, '(" SURF VECS",2(5X,2F9.4))') a1
    write (unit, '(15X,2F9.4)') a2
  end subroutine write_substrate_vectors

  !> Echo the header of the surface parameters and the superlattice
  !! vectors, in Angstrom.
  subroutine write_superlattice_vectors(unit, a1, a2)
    integer, intent(in) :: unit
    real(sp), intent(in) :: a1(2), a2(2)

    write (unit, '(" PARAMETERS FOR SURFACE ")')
    write (unit, '(" SURF VECS",2(5X,2F9.4))') a1
    write (unit, '(15X,2F9.4)') a2
  end subroutine write_superlattice_vectors

  !> Echo imaginary part of the inner potential, surface distance (both as
  !! read) and the work function.
  subroutine write_inner_potential(unit, vpi, ase, workfn)
    integer, intent(in) :: unit
    real(sp), intent(in) :: vpi, ase, workfn

    write (unit, '(" VPI =",F7.4,"    ASE = ",F7.4,"    WORKFN = ",F7.4)') &
      vpi, ase, workfn
  end subroutine write_inner_potential

  !> Echo the criterion for the selection of beams.
  subroutine write_beam_criterion(unit, tst)
    integer, intent(in) :: unit
    real(sp), intent(in) :: tst

    write (unit, '(" TST = ",1F11.8)') tst
  end subroutine write_beam_criterion

  !> Echo the beams to be output.
  subroutine write_output_beams(unit, beams)
    integer, intent(in) :: unit
    integer, intent(in) :: beams(:)

    write (unit, '(" PUNCHED BEAMS   ",15I5)') beams
  end subroutine write_output_beams

  !> Echo the direction of incidence, in degrees.
  subroutine write_incidence(unit, theta, fi)
    integer, intent(in) :: unit
    real(sp), intent(in) :: theta, fi

    write (unit, '(" THETA  FI = ",2F7.2)') theta, fi
  end subroutine write_incidence

  !> Echo the convergence criteria of the layer doubling.
  subroutine write_convergence(unit, eps, liter)
    integer, intent(in) :: unit
    real(sp), intent(in) :: eps
    integer, intent(in) :: liter

    write (unit, '(" EPS = ",1F7.4,"  LITER = ",1I3)') eps, liter
  end subroutine write_convergence

  !> Echo the highest angular momentum.
  subroutine write_lmax(unit, lmax)
    integer, intent(in) :: unit
    integer, intent(in) :: lmax

    write (unit, '(" LMAX = ",1I3)') lmax
  end subroutine write_lmax

  !> Echo the header of the phase-shift block.
  subroutine write_phase_shift_header(unit)
    integer, intent(in) :: unit

    write (unit, '(" PHASE SHIFTS")')
  end subroutine write_phase_shift_header

  !> Echo the format of the phase-shift file.
  subroutine write_phase_shift_format(unit, psfmt)
    integer, intent(in) :: unit
    integer, intent(in) :: psfmt

    write (unit, '('' PHASESHIFT FORMAT IS'',I3)') psfmt
  end subroutine write_phase_shift_format

  !> Echo the format of the tensor output.
  subroutine write_iform(unit, iform)
    integer, intent(in) :: unit
    integer, intent(in) :: iform

    write (unit, *) 'IFORM = ', iform
  end subroutine write_iform

  !> Echo the header of the parameters for the inner potential.
  subroutine write_inner_potential_header(unit)
    integer, intent(in) :: unit

    write (unit, '(" PARAMETERS FOR INNER POTENTIAL")')
  end subroutine write_inner_potential_header

  !> Echo the inner-potential coefficients of a phase-shift file of format 1.
  subroutine write_phase_shift_coefficients(unit, c)
    integer, intent(in) :: unit
    real(sp), intent(in) :: c(4)

    write (unit, *) "   C1      C2      C3      C4   "
    write (unit, '(4F8.2)') c
  end subroutine write_phase_shift_coefficients

  !> Echo the title line of a phase-shift file of format 2.
  subroutine write_phase_shift_file_title(unit, title)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: title

    write (unit, '(A70)') title
  end subroutine write_phase_shift_file_title

  !> Echo the inner-potential parameters of a phase-shift file of format 2.
  subroutine write_phase_shift_parameters(unit, em, c)
    integer, intent(in) :: unit
    real(sp), intent(in) :: em, c(8)

    write (unit, '(F6.2,1x,4F7.2,2x,F7.2,3F8.4)') em, c
  end subroutine write_phase_shift_parameters

  !> Echo the phase shifts at the first energy for one element.
  subroutine write_first_phase_shifts(unit, ielement, energy, phase)
    integer, intent(in) :: unit
    integer, intent(in) :: ielement
    real(sp), intent(in) :: energy
    real(sp), intent(in) :: phase(:)

    write (unit, '(" 1st ENERGY FOR ELEMENT/SITE:",I3)') ielement
    write (unit, '(" E = ",1F7.4,3X,20F8.4)') energy, phase
  end subroutine write_first_phase_shifts

  !> Write the empty lines that close the echo of the input.
  subroutine write_end_of_input(unit)
    integer, intent(in) :: unit

    write (unit, '(/)')
  end subroutine write_end_of_input
end module leed_output
