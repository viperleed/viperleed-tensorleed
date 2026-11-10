!> \file leed_types.f90
!! \brief Derived types holding the input of the reference calculation.
!!
!! The types are filled by the reader in leed_input.f90 and are sized from
!! the counts in the input. The legacy name of a quantity is given in a
!! comment where it helps to compare with lib.tleed.f and with the
!! TensErLEED paper (Blum and Heinz, Comput. Phys. Commun. 134, 392 (2001)).

!> Types of the input of ref-calc.
module leed_types
  use kinds, only: sp
  implicit none

  private

  public :: atom_sites, layer_type, stacked_layer, geometry
  public :: lattice, control, inner_potential, beam_list, phase_shifts

  !> The kinds of atoms, distinguished by element, vibrational amplitude
  !! and concentration.
  type :: atom_sites
    !> Concentration of element iel on site type isite; CONC(isite, iel).
    real(sp), allocatable :: concentration(:, :)
    !> Vibrational amplitude of element iel on site type isite in atomic
    !! units; VIB(isite, iel).
    real(sp), allocatable :: vibration_amplitude(:, :)
  end type atom_sites

  !> One layer type: a Bravais layer (one subplane) or a composite layer.
  type :: layer_type
    !> Lattice of the layer type: 1 overlayer, 2 substrate; LATT.
    integer :: lattice = 0
    !> Number of subplanes; NSUB.
    integer :: n_subplanes = 0
    !> Number of this layer among the Bravais layers, 0 for a composite
    !! layer; LBRAV.
    integer :: bravais_index = 0
    !> Number of this layer among the composite layers, 0 for a Bravais
    !! layer; LCOMP.
    integer :: composite_index = 0
    !> Site type of each subplane, (n_subplanes); STYPE(ilayer, :).
    integer, allocatable :: site_type(:)
    !> Position of each subplane in atomic units, (3, n_subplanes);
    !! SUBPOS(ilayer, :, :) with the coordinate index first.
    real(sp), allocatable :: subplane_position(:, :)
  end type layer_type

  !> One layer of the stack on top of the bulk.
  type :: stacked_layer
    !> Layer type of this layer; LTYPE.
    integer :: layer_type = 0
    !> Interlayer vector to the next lower layer in atomic units; LDIST.
    real(sp) :: distance(3) = 0.0_sp
    !> Tensor output wanted for this layer (0 no, 1 yes); TENS.
    integer :: tensor_output = 0
    !> Name of the tensor file of each subplane, (n_subplanes of the layer
    !! type); only allocated if tensor_output is 1; LAYFILE(istack, :).
    character(len=40), allocatable :: tensor_file(:)
  end type stacked_layer

  !> Atoms, layer types and layer stacking of the reference surface.
  type :: geometry
    type(atom_sites) :: sites
    !> The layer types; index as in the input.
    type(layer_type), allocatable :: layer_types(:)
    !> 0 if a bulk calculation is required, 1 if not; TSLAB.
    integer :: bulk_slab = 0
    !> Interlayer vector between bulk units in atomic units; ASA.
    real(sp) :: bulk_unit_vector(3) = 0.0_sp
    !> Layer type of the top layer of a bulk unit; TOPLAYB.
    integer :: bulk_top_layer = 0
    !> Layer type of the bottom layer of a bulk unit; BOTLAYB.
    integer :: bulk_bottom_layer = 0
    !> Interlayer vector between top and bottom layer of a bulk unit in
    !! atomic units; ASBULK.
    real(sp) :: bulk_internal_vector(3) = 0.0_sp
    !> Interlayer vector of smallest distance, for the convergence
    !! criterion; ASB. Only the first component is set, and only if a bulk
    !! calculation is required.
    real(sp) :: min_distance_vector(3)
    !> Layers on top of the bulk; index 1 is the top layer; (LTYPE, LDIST,
    !! TENS, LAYFILE).
    type(stacked_layer), allocatable :: stack(:)
  end type geometry

  !> Two-dimensional lattices of the substrate and of the superlattice
  !! (legacy common block /SL/ and the vectors READIN returned).
  !! Lengths are in atomic units, reciprocal vectors in inverse atomic
  !! units. Each vector is (2): the components in the surface plane.
  type :: lattice
    !> Basis vectors of the substrate lattice; ARA1, ARA2.
    real(sp) :: substrate_vector(2, 2) = 0.0_sp
    !> Reciprocal basis vectors of the substrate lattice; RAR1, RAR2.
    real(sp) :: substrate_reciprocal(2, 2) = 0.0_sp
    !> Basis vectors of the superlattice; ARB1, ARB2.
    real(sp) :: superlattice_vector(2, 2) = 0.0_sp
    !> Reciprocal basis vectors of the superlattice; RBR1, RBR2.
    real(sp) :: superlattice_reciprocal(2, 2) = 0.0_sp
    !> Area of the substrate unit mesh; TVA.
    real(sp) :: substrate_area = 0.0_sp
    !> Area of the superlattice unit mesh; TVB.
    real(sp) :: superlattice_area = 0.0_sp
  end type lattice

  !> Convergence and angular-momentum controls (legacy common block /MS/).
  type :: control
    !> Highest angular momentum considered; LMAX.
    integer :: lmax = 0
    !> Convergence criterion of the layer doubling; EPS.
    real(sp) :: eps = 0.0_sp
    !> Maximum number of iterations of the layer doubling; LITER.
    integer :: liter = 0
  end type control

  !> Parameters of the inner potential (legacy common blocks /ADS/ and
  !! /INPOT/). Components that a phase-shift file format does not provide
  !! stay zero.
  type :: inner_potential
    !> Distance between the topmost atomic nuclei and the onset of the
    !! inner potential; ASE.
    real(sp) :: surface_distance = 0.0_sp
    !> Imaginary part of the inner potential at the substrate and at the
    !! overlayer; VPIS, VPIO.
    real(sp) :: vpi_substrate = 0.0_sp
    real(sp) :: vpi_overlayer = 0.0_sp
    !> Parameters of the phase-shift program for the real part of the
    !! inner potential; EM, C1 ... C8.
    real(sp) :: em = 0.0_sp
    real(sp) :: c(8) = 0.0_sp
    !> Work function; WORKFN.
    real(sp) :: work_function = 0.0_sp
  end type inner_potential

  !> The beams of the calculation, in sets, and those to be output.
  type :: beam_list
    !> Number of beams in each beam set, (n_beam_sets); KNB.
    integer, allocatable :: set_size(:)
    !> Beams in units of the reciprocal substrate vectors, (2, n_beams);
    !! SPQF.
    real(sp), allocatable :: index(:, :)
    !> Beams in atomic units, (2, n_beams); SPQ.
    real(sp), allocatable :: momentum(:, :)
    !> Symmetry codes of the beams for the superlattice and the substrate
    !! region, (2, n_beams); KSYM.
    integer, allocatable :: symmetry(:, :)
    !> Beams of which the intensities are written out, as positions in
    !! the beam list, (n_output); NPU.
    integer, allocatable :: output(:)
  end type beam_list

  !> The tabulated phase shifts.
  type :: phase_shifts
    !> Energies in Hartree at which phase shifts are tabulated,
    !! (n_energies); ES.
    real(sp), allocatable :: energy(:)
    !> Phase shifts in radian, (n_energies, n_elements, lmax + 1); PHSS.
    real(sp), allocatable :: value(:, :, :)
  end type phase_shifts
end module leed_types
