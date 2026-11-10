!> \file leed_input.f90
!! \brief Reader of the input of ref-calc, into the types of leed_types.
!!
!! Every routine takes the unit to read from, which the caller has opened.
!! Invalid input stops with a message on stderr and a non-zero exit code.

! read_input is a port of SUBROUTINE READIN of lib.tleed.f, which reads in
! most of the input for the LEED programs (lateral lattices, calculational
! features, beam list and phase shifts).
!   TST is read in F11.8 format (07/05/92). Phase shifts are output for the
!   first energy only. The number of elements and lmax + 1 are required for
!   the dimension of the phase shifts (07/05/92).

! Ported from SUBROUTINE READGEO of lib.tleed.f, which reads structural
! (and some calculational) information for the reference surface in
! question.
! VB, 06/98

!> Reader of the input of ref-calc.
module leed_input
  use kinds, only: sp
  use leed_output
  use leed_types, only: beam_list, control, geometry, inner_potential, &
                        lattice, phase_shifts
  implicit none

  private

  public :: read_input, read_geometry

  !> Bohr radius in Angstrom (BOHR of the legacy file GLOBAL).
  real(sp), parameter :: bohr = 0.529177_sp
  !> Hartree in eV (HARTREE of the legacy file GLOBAL).
  real(sp), parameter :: hartree = 27.211396_sp
  real(sp), parameter :: pi = 3.14159265359_sp

contains

  !> Read the lattices, calculational features, beam list and phase shifts.
  !!
  !! The sizes the program was built with are compared with the counts of
  !! the input; a mismatch stops the program. IDEG, NL1, NL2, NPUN and NPSI
  !! are records of the input (IDEG, NL1 and NL2 after the inner-potential
  !! record, NPUN before the output beams, NPSI after IFORM); their values
  !! must equal those the program was built with. Lengths are converted from
  !! Angstrom to atomic units, angles from degrees to radian, the imaginary
  !! part of the inner potential from eV to Hartree (and made negative).
  !! The phase shifts are echoed for the first energy only.
  subroutine read_input(unit, out_unit, n_sublattices, n_beam_sets, &
                        n_beams, n_output, n_energies, n_elements, lmax1, &
                        ideg, nl1, nl2, title, lat, ctrl, pot, vpi, ei, ef, &
                        de, beams, tst, theta, fi, iform, phase)
    !> Unit to read from
    integer, intent(in) :: unit
    !> Unit to echo the input to
    integer, intent(in) :: out_unit
    !> Number of sublattices, must equal n_beam_sets; NL
    integer, intent(in) :: n_sublattices
    !> Number of beam sets; KNBS
    integer, intent(in) :: n_beam_sets
    !> Total number of beams, must equal the input; KNT
    integer, intent(in) :: n_beams
    !> Number of beams to output; NPUN
    integer, intent(in) :: n_output
    !> Number of energies at which phase shifts are read; NPSI
    integer, intent(in) :: n_energies
    !> Number of elements in the phase shifts, must equal the input; NEL
    integer, intent(in) :: n_elements
    !> Highest angular momentum plus one, must be the input lmax + 1;
    !! LMAX1
    integer, intent(in) :: lmax1
    !> Degree of the symmetry used in the layer-matrix routines, must equal
    !! the input; IDEG
    integer, intent(in) :: ideg
    !> Repetitions of the substrate mesh along each superlattice vector,
    !! must equal the input; NL1, NL2
    integer, intent(in) :: nl1, nl2
    character(len=80), intent(out) :: title
    type(lattice), intent(out) :: lat
    type(control), intent(out) :: ctrl
    type(inner_potential), intent(out) :: pot
    !> Imaginary part of the inner potential; VPI
    real(sp), intent(out) :: vpi
    !> Initial energy, final energy and step; EI, EF, DE
    real(sp), intent(out) :: ei, ef, de
    type(beam_list), intent(out) :: beams
    !> Criterion for the selection of beams; TST
    real(sp), intent(out) :: tst
    !> Polar and azimuthal angle of incidence in radian; THETA, FI
    real(sp), intent(out) :: theta, fi
    !> Format of the tensor output; IFORM
    integer, intent(out) :: iform
    type(phase_shifts), intent(out) :: phase

    character(len=80) :: ps_comment, ps_title
    real(sp) :: atv, ase, workfn
    real(sp) :: sub_a(2, 2), super_a(2, 2)
    integer :: psfmt, ps_lmax, nel_input, n_beams_input
    integer :: ideg_input, nl1_input, nl2_input, n_output_input, n_energies_input
    integer :: i, j, k, kk, n_in_set, ienergy, iel, l

    ! Consistency check
    if (n_sublattices /= n_beam_sets) then
      call stop_with_error("Number of beam sets, MKNBS, differs from " &
                           // "number of sublattices, MNL = MNL1*MNL2! " &
                           // "Assuming incorrect input, I stop!")
    end if

    read (unit, '(A80)') title
    call write_title(out_unit, title)

    read (unit, '(3F9.2)') ei, ef, de

    ! The two basis vectors of the substrate layer lattice, in terms of the
    ! planar cartesian Y- and Z-axes (the X-axis is perpendicular to the
    ! surface).
    read (unit, '(3F9.4)') (sub_a(i, 1), i=1, 2)
    read (unit, '(3F9.4)') (sub_a(i, 2), i=1, 2)
    call write_substrate_vectors(out_unit, sub_a(:, 1), sub_a(:, 2))
    sub_a = sub_a / bohr
    lat%substrate_vector = sub_a
    lat%substrate_area = abs(sub_a(1, 1) * sub_a(2, 2) &
                             - sub_a(2, 1) * sub_a(1, 2))
    atv = 2.0_sp * pi / lat%substrate_area
    ! Reciprocal-lattice basis vectors of the substrate
    lat%substrate_reciprocal(1, 1) = sub_a(2, 2) * atv
    lat%substrate_reciprocal(2, 1) = -sub_a(1, 2) * atv
    lat%substrate_reciprocal(1, 2) = -sub_a(2, 1) * atv
    lat%substrate_reciprocal(2, 2) = sub_a(1, 1) * atv

    ! The same for the superlattice of an overlayer
    read (unit, '(3F9.4)') (super_a(i, 1), i=1, 2)
    read (unit, '(3F9.4)') (super_a(i, 2), i=1, 2)
    call write_superlattice_vectors(out_unit, super_a(:, 1), super_a(:, 2))
    super_a = super_a / bohr
    lat%superlattice_vector = super_a
    lat%superlattice_area = abs(super_a(1, 1) * super_a(2, 2) &
                                - super_a(2, 1) * super_a(1, 2))
    atv = 2.0_sp * pi / lat%superlattice_area
    lat%superlattice_reciprocal(1, 1) = super_a(2, 2) * atv
    lat%superlattice_reciprocal(2, 1) = -super_a(1, 2) * atv
    lat%superlattice_reciprocal(1, 2) = -super_a(2, 1) * atv
    lat%superlattice_reciprocal(2, 2) = super_a(1, 1) * atv

    ! Imaginary part of the inner potential, and spacing between the
    ! surface (where muffin-tin constant and damping set in) and the
    ! top-layer nuclei, in Angstrom.
    read (unit, '(3F9.4)') vpi, ase, workfn
    call write_inner_potential(out_unit, vpi, ase, workfn)
    pot%surface_distance = ase / bohr
    pot%work_function = workfn
    vpi = -vpi / hartree
    pot%vpi_substrate = vpi
    pot%vpi_overlayer = vpi

    read (unit, *) ideg_input, nl1_input, nl2_input
    if (ideg_input /= ideg) then
      call stop_with_error("IDEG = " // int_text(ideg_input) &
                           // " in the input differs from MIDEG = " &
                           // int_text(ideg) // " in PARAM")
    end if
    if (nl1_input /= nl1) then
      call stop_with_error("NL1 = " // int_text(nl1_input) &
                           // " in the input differs from MNL1 = " &
                           // int_text(nl1) // " in PARAM")
    end if
    if (nl2_input /= nl2) then
      call stop_with_error("NL2 = " // int_text(nl2_input) &
                           // " in the input differs from MNL2 = " &
                           // int_text(nl2) // " in PARAM")
    end if
    lat%symmetry_degree = ideg_input
    lat%repetitions(1) = nl1_input
    lat%repetitions(2) = nl2_input

    ! The beam sets, each as reciprocal-lattice vectors in terms of the
    ! substrate ones and with their symmetry codes.
    allocate (beams%set_size(n_beam_sets))
    allocate (beams%index(2, n_beams), beams%momentum(2, n_beams))
    allocate (beams%symmetry(2, n_beams))
    n_beams_input = 0
    do j = 1, n_beam_sets
      read (unit, '(15I5)') beams%set_size(j)
      n_in_set = beams%set_size(j)
      do k = 1, n_in_set
        kk = k + n_beams_input
        read (unit, '(2F10.5,4I3)') beams%index(1, kk), beams%index(2, kk), &
                                    (beams%symmetry(i, kk), i=1, 2)
        do i = 1, 2
          beams%momentum(i, kk) &
              = beams%index(1, kk) * lat%substrate_reciprocal(i, 1) &
                + beams%index(2, kk) * lat%substrate_reciprocal(i, 2)
        end do
      end do
      n_beams_input = n_beams_input + beams%set_size(j)
    end do
    if (n_beams_input /= n_beams) then
      call stop_with_error("Number of beams in full dynamic beam list, " &
                           // int_text(n_beams_input) &
                           // " is inconsistent with dimension MKNT = " &
                           // int_text(n_beams) &
                           // " from PARAM! Please correct!")
    end if

    ! Criterion for the selection of beams: a beam may not decay to less
    ! than this fraction of its initial value travelling from any layer to
    ! the next.
    read (unit, '(1F11.8)') tst
    call write_beam_criterion(out_unit, tst)

    ! Beams for which intensities are output, as positions in the input list
    read (unit, *) n_output_input
    if (n_output_input /= n_output) then
      call stop_with_error("NPUN = " // int_text(n_output_input) &
                           // " in the input differs from MNPUN = " &
                           // int_text(n_output) // " in PARAM")
    end if
    allocate (beams%output(n_output_input))
    read (unit, '(15I5)') (beams%output(k), k=1, n_output_input)
    call write_output_beams(out_unit, beams%output)

    ! Direction of incidence (degrees)
    read (unit, '(2F9.4)') theta, fi
    call write_incidence(out_unit, theta, fi)
    theta = theta * pi / 180.0_sp
    fi = fi * pi / 180.0_sp

    ! Convergence criteria for layer doubling
    read (unit, '(F7.4)') ctrl%eps
    read (unit, '(26I3)') ctrl%liter
    call write_convergence(out_unit, ctrl%eps, ctrl%liter)

    read (unit, '(26I3)') ctrl%lmax
    call write_lmax(out_unit, ctrl%lmax)
    if (ctrl%lmax /= lmax1 - 1) then
      call stop_with_error("Please enforce LMAX in input = MLMAX in " &
                           // "PARAM. Inaccurate dimensions may lead to " &
                           // "problems.")
    end if

    call write_phase_shift_header(out_unit)
    read (unit, '(26I3)') psfmt
    call write_phase_shift_format(out_unit, psfmt)
    read (unit, '(26I3)') iform
    call write_iform(out_unit, iform)
    read (unit, *) n_energies_input
    if (n_energies_input /= n_energies) then
      call stop_with_error("NPSI = " // int_text(n_energies_input) &
                           // " in the input differs from MNPSI = " &
                           // int_text(n_energies) // " in PARAM")
    end if
    call write_inner_potential_header(out_unit)

    allocate (phase%energy(n_energies_input))
    allocate (phase%value(n_energies_input, n_elements, lmax1))

    select case (psfmt)
    case (1)
      read (unit, '(I3,1x,4F8.2)') nel_input, pot%c(4), pot%c(1), pot%c(2), &
                                   pot%c(3)
      if (nel_input /= n_elements) then
        call stop_with_error("NEL from phaseshift file not consistent " &
                             // "with MNEL! Please correct!")
      end if
      call write_phase_shift_coefficients(out_unit, pot%c(1:4))

      ! Energies (Hartree) at which phase shifts are input, and the phase
      ! shifts (radian). Linear interpolation occurs for energies between
      ! the tabulated ones (and linear extrapolation above the highest).
      do ienergy = 1, n_energies_input
        read (unit, '(F7.4)') phase%energy(ienergy)
        do iel = 1, n_elements
          read (unit, '(10F7.4)') (phase%value(ienergy, iel, l), l=1, lmax1)
          if (ctrl%lmax < 10) read (unit, *)
        end do
      end do

    case (2)
      read (unit, '(26I3)') nel_input
      if (nel_input /= n_elements) then
        call stop_with_error("NEL from phaseshift file not consistent " &
                             // "with MNEL! Please correct!")
      end if

      read (unit, '(1A8,1A70)') ps_comment, ps_title
      call write_phase_shift_file_title(out_unit, ps_title)
      read (unit, '(3x,I3,2x,F6.2,1x,4F7.2,2x,F7.2,3F8.4)') ps_lmax, &
          pot%em, pot%c
      call write_phase_shift_parameters(out_unit, pot%em, pot%c)

      do iel = 1, n_elements
        read (unit, *)
        read (unit, *)
        do ienergy = 1, n_energies_input
          read (unit, '(F8.3,20F8.4)') phase%energy(ienergy), &
              (phase%value(ienergy, iel, l), l=1, lmax1)
          phase%energy(ienergy) = phase%energy(ienergy) / hartree
        end do
      end do

    case default
      call stop_with_error("Phaseshift format incorrectly assigned!")
    end select

    do iel = 1, n_elements
      call write_first_phase_shifts(out_unit, iel, phase%energy(1), &
                                    phase%value(1, iel, :))
    end do
    call write_end_of_input(out_unit)
  end subroutine read_input

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
