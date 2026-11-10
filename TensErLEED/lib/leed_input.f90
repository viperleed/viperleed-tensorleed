!> \file leed_input.f90
!! \brief Reader of the input of ref-calc, into the types of leed_types.
!!
!! Every routine takes the unit to read from, which the caller has opened.
!! Invalid input stops with a message on stderr and a non-zero exit code.

! Ported from SUBROUTINE READGEO of lib.tleed.f, which reads structural
! (and some calculational) information for the reference surface in
! question.
! VB, 06/98

!> Reader of the input of ref-calc.
module leed_input
  use kinds, only: sp
  use leed_output, only: stop_with_error, write_no_bulk_calculation, &
                         write_tensor_file, write_tensor_not_required, &
                         write_tensor_required
  use leed_types, only: geometry
  implicit none

  private

  public :: read_geometry

  !> Bohr radius in Angstrom (BOHR of the legacy file GLOBAL).
  real(sp), parameter :: bohr = 0.529177_sp

contains

  !> Read the structural information for the reference calculation.
  !!
  !! The counts of the input are compared with the sizes the program was
  !! built with; a mismatch stops the program. The first line of each block
  !! of the input is a header and skipped. Lengths are converted from
  !! Angstrom to atomic units. The stacked layers are read from the lowest
  !! (index n_stack) to the top one (index 1).
  subroutine read_geometry(unit, out_unit, n_elements, n_sites, &
                           n_layer_types, max_subplanes, n_stack, &
                           n_bravais, n_sublattices, geom)
    !> Unit to read from
    integer, intent(in) :: unit
    !> Unit to echo the tensor-output choices to
    integer, intent(in) :: out_unit
    !> Number of elements in the phase shifts; NEL
    integer, intent(in) :: n_elements
    !> Number of site types, must equal the input; NSITE
    integer, intent(in) :: n_sites
    !> Number of layer types, must equal the input; NLTYPE
    integer, intent(in) :: n_layer_types
    !> Maximum number of subplanes of a layer type; NNSUB
    integer, intent(in) :: max_subplanes
    !> Number of layers stacked on the bulk, must equal the input; NSTACK
    integer, intent(in) :: n_stack
    !> Number of Bravais layer types, must equal the input; NBRAV
    integer, intent(in) :: n_bravais
    !> Number of superlattice sublattices in the substrate lattice; NL
    integer, intent(in) :: n_sublattices
    type(geometry), intent(out) :: geom

    integer :: n_sites_input, n_layer_types_input, n_stack_input
    integer :: n_bravais_input, n_composite
    integer :: isite, iel, iltype, isub, ipos, istack, i

    ! First, read information for different atom types (distinguished by
    ! element, vibrational amplitude, and stoichiometry).
    read (unit, *)
    read (unit, *)
    read (unit, *)

    read (unit, '(I3)') n_sites_input
    if (n_sites_input /= n_sites) then
      call stop_with_error("NSITE .ne. MNSITE - please correct!")
    end if

    allocate (geom%sites%concentration(n_sites, n_elements))
    allocate (geom%sites%vibration_amplitude(n_sites, n_elements))
    do isite = 1, n_sites
      read (unit, *)
      do iel = 1, n_elements
        ! AMI: note change F7.4 -> F9.4 in v1.71 by LH
        read (unit, '(2F9.4)') geom%sites%concentration(isite, iel), &
                               geom%sites%vibration_amplitude(isite, iel)
        geom%sites%vibration_amplitude(isite, iel) &
            = geom%sites%vibration_amplitude(isite, iel) / bohr
      end do
    end do

    ! Next, read information on different layer types to be included.
    read (unit, *)
    read (unit, *)
    read (unit, *)
    read (unit, '(I3)') n_layer_types_input
    if (n_layer_types_input /= n_layer_types) then
      call stop_with_error("NLTYPE .ne. MNLTYPE - please correct!")
    end if

    ! n_bravais_input and n_composite only count layer types
    n_bravais_input = 0
    n_composite = 0
    allocate (geom%layer_types(n_layer_types))

    do iltype = 1, n_layer_types
      associate (layer => geom%layer_types(iltype))
        read (unit, *)
        read (unit, '(I3)') layer%lattice
        read (unit, '(I3)') layer%n_subplanes

        if (layer%n_subplanes > max_subplanes) then
          call stop_with_error("Number of subplanes in layer type " &
                               // int_text(iltype) &
                               // " exceeds dimension MNSUB - please correct!")
        else if (layer%n_subplanes <= 0) then
          call stop_with_error("Number of subplanes in layer type " &
                               // int_text(iltype) &
                               // " less than 1 ? - please correct!")
        end if

        if (layer%n_subplanes == 1) then
          n_bravais_input = n_bravais_input + 1
          layer%bravais_index = n_bravais_input
          layer%composite_index = 0
        else
          n_composite = n_composite + 1
          layer%composite_index = n_composite
          layer%bravais_index = 0
        end if

        allocate (layer%site_type(layer%n_subplanes))
        allocate (layer%subplane_position(3, layer%n_subplanes))

        do isub = 1, layer%n_subplanes
          ! used to be 3F7.4 pre v1.7
          read (unit, '(I3,3F9.4)') layer%site_type(isub), &
              (layer%subplane_position(ipos, isub), ipos=1, 3)

          if (layer%site_type(isub) < 1 .or. layer%site_type(isub) > n_sites) then
            call stop_with_error("Site type " // int_text(layer%site_type(isub)) &
                                 // " in sublayer " // int_text(isub) &
                                 // " of layer " // int_text(iltype) &
                                 // ": No such site defined - please correct!")
          end if

          do ipos = 1, 3
            layer%subplane_position(ipos, isub) &
                = layer%subplane_position(ipos, isub) / bohr
          end do
        end do

        ! Dead check, kept as in the legacy code: after the loop above isub
        ! is n_subplanes + 1 >= 2, so it never fires.
        if (isub == 1 .and. any(layer%subplane_position(:, 1) /= 0.0_sp)) then
          call stop_with_error("For a Bravais layer, a subplane position " &
                               // "makes no sense and should be (0,0,0). " &
                               // "To make sure, please correct.")
        end if
      end associate
    end do

    if (n_bravais_input /= n_bravais) then
      call stop_with_error("Number of Bravais layers in actual input does " &
                           // "not match dimension MNBRAV - please correct!")
    end if

    ! Now gather info for layer doubling.
    read (unit, *)
    read (unit, *)
    read (unit, *)

    read (unit, '(I3)') geom%bulk_slab

    if (geom%bulk_slab == 0) then
      ! ASA is the substrate interlayer vector (linking reference points
      ! on successive layers) (Angstrom). Used to be 3F7.4 pre v1.7.
      read (unit, '(3F9.4)') (geom%bulk_unit_vector(i), i=1, 3)
      read (unit, '(I3)') geom%bulk_top_layer
      read (unit, '(I3)') geom%bulk_bottom_layer
      read (unit, '(3F9.4)') (geom%bulk_internal_vector(i), i=1, 3)

      if (geom%bulk_top_layer < 1 .or. geom%bulk_top_layer > n_layer_types) then
        call stop_with_error("Top bulk layer is illegal layer type. " &
                             // "Please correct.")
      end if

      if (geom%bulk_bottom_layer < 1 &
          .or. geom%bulk_bottom_layer > n_layer_types) then
        call stop_with_error("Bottom bulk layer is illegal layer type. " &
                             // "Please correct.")
      end if

      do i = 1, 3
        geom%bulk_unit_vector(i) = geom%bulk_unit_vector(i) / bohr
        geom%bulk_internal_vector(i) = geom%bulk_internal_vector(i) / bohr
      end do

      geom%min_distance_vector(1) = geom%bulk_internal_vector(1)

    else if (geom%bulk_slab == 1) then
      call write_no_bulk_calculation(out_unit)
      read (unit, *)
      read (unit, *)
      read (unit, *)
      read (unit, *)
    else
      call stop_with_error("TSLAB should be 0 or 1 - please correct!")
    end if

    ! Next, read info for layer stacking.
    read (unit, *)
    read (unit, *)
    read (unit, *)

    read (unit, '(I3)') n_stack_input
    if (n_stack_input /= n_stack) then
      call stop_with_error("NSTACK .ne. MNSTACK - please correct!")
    end if

    allocate (geom%stack(n_stack))

    do istack = n_stack, 1, -1
      associate (layer => geom%stack(istack))
        ! used to be 3F7.4 pre v1.7
        read (unit, '(I3,3F9.4)') layer%layer_type, (layer%distance(i), i=1, 3)

        if (layer%layer_type < 1 .or. layer%layer_type > n_layer_types) then
          call stop_with_error("Stacked layer number " // int_text(istack) &
                               // " has illegal layer type. Please correct.")
        end if

        do ipos = 1, 3
          layer%distance(ipos) = layer%distance(ipos) / bohr
        end do

        if (abs(layer%distance(1)) <= geom%min_distance_vector(1)) then
          geom%min_distance_vector(1) = layer%distance(1)
        end if

        read (unit, '(I3)') layer%tensor_output

        if (layer%tensor_output == 0) then
          call write_tensor_not_required(out_unit, istack)
        else if (layer%tensor_output == 1) then
          if (n_sublattices > 1 &
              .and. geom%layer_types(layer%layer_type)%lattice == 2) then
            call stop_with_error("In the presence of a non-trivial " &
                                 // "superlattice, Tensor LEED is currently " &
                                 // "not possible for a substrate type layer.")
          end if

          call write_tensor_required(out_unit, istack)

          allocate (layer%tensor_file( &
                    geom%layer_types(layer%layer_type)%n_subplanes))
          do isub = 1, geom%layer_types(layer%layer_type)%n_subplanes
            read (unit, '(A20)') layer%tensor_file(isub)
            call write_tensor_file(out_unit, istack, isub, &
                                   layer%tensor_file(isub))
          end do
        else
          call stop_with_error("Tensor LEED flag for layer " &
                               // int_text(istack) &
                               // " not 0 or 1 - please correct!")
        end if
      end associate
    end do
  end subroutine read_geometry

  !> The integer as text, without blanks.
  function int_text(value) result(text)
    integer, intent(in) :: value
    character(len=:), allocatable :: text

    character(len=16) :: buffer

    write (buffer, '(I0)') value
    text = trim(buffer)
  end function int_text
end module leed_input
