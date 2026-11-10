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
end module leed_types
