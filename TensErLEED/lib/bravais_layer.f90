!> \file bravais_layer.f90
!! \brief Lattice sums and scattering of Bravais layers: sublattice setup.
!!
!! Each routine keeps the precision of every quantity of the legacy routine
!! it follows, and its results are bit-for-bit those of the legacy code
!! (tests/tier1/test_refcalc_routines.py).

!> Setup for the lattice sums of Bravais layers.
module bravais_layer
  use constants, only: pi
  use kinds, only: sp
  use leed_types, only: lattice
  implicit none

  private

  public :: sublattice_equivalences

  !> Small number to compare the phase factors of sublattice vectors
  !! (the legacy argument EPSD).
  real(sp), parameter :: epsd = 2.0e-4_sp

contains

  !> Sets up the sublattice vectors and their equivalences under rotation
  !! (legacy SLIND).
  !!
  !! Finds how the sublattices js of the superlattice are transformed into
  !! one another by rotations through j * 2 pi / ideg: jjs(js, j) is the
  !! sublattice into which sublattice js is rotated by the j-th rotation.
  !! If no sublattice matches, jjs is nl + 1. Author: Pendry.
  !!
  !! The substrate lattice vectors define the sublattices, the reciprocal
  !! vectors of the superlattice identify them. The arrays v, vl and jjs are
  !! explicitly shaped on purpose: ref-calc declares them with maximum sizes,
  !! and the lattice-sum routines take the first nl * 2 (nl * ideg) elements
  !! as arrays of the shape (nl, 2) ((nl, ideg)).
  subroutine sublattice_equivalences(lat, nl, v, vl, jjs)
    !> Lattices; uses substrate_vector, superlattice_reciprocal,
    !! symmetry_degree (ideg) and repetitions (nl1, nl2)
    type(lattice), intent(in) :: lat
    !> Number of sublattices; nl1 * nl2
    integer, intent(in) :: nl
    !> Vectors pointing to the sublattices, in polar coordinates (modulus,
    !! angle)
    real(sp), intent(out) :: v(nl, 2)
    !> Phase factors of the sublattice vectors for the two reciprocal vectors
    !! of the superlattice
    complex(sp), intent(out) :: vl(nl, 2)
    !> Relationship of the sublattices under rotation
    integer, intent(out) :: jjs(nl, lat%symmetry_degree)

    integer :: i, j, k, ideg, jl1, kl2, nl1, nl2
    real(sp) :: s1, s2, adr1, adr2, adr, ang, ainc, cang, sang, a, b, test
    complex(sp) :: vla, vlb, ci
    ! Legacy names: br1, br2 = substrate_vector; rar1, rar2 =
    ! superlattice_reciprocal
    real(sp) :: br1(2), br2(2), rar1(2), rar2(2)

    br1 = lat%substrate_vector(:, 1)
    br2 = lat%substrate_vector(:, 2)
    rar1 = lat%superlattice_reciprocal(:, 1)
    rar2 = lat%superlattice_reciprocal(:, 2)
    ideg = lat%symmetry_degree
    nl1 = lat%repetitions(1)
    nl2 = lat%repetitions(2)

    ci = cmplx(0.0_sp, 1.0_sp, sp)
    ! Set up vectors v defining sublattices and quantities vl for later
    ! reference.
    i = 1
    s1 = 0.0_sp
    do jl1 = 1, nl1
      s2 = 0.0_sp
      do kl2 = 1, nl2
        adr1 = s1 * br1(1) + s2 * br2(1)
        adr2 = s1 * br1(2) + s2 * br2(2)
        v(i, 1) = sqrt(adr1 * adr1 + adr2 * adr2)
        v(i, 2) = 0.0_sp
        if (v(i, 1) > 0.0_sp) v(i, 2) = atan2(adr2, adr1)
        vl(i, 1) = exp(ci * (adr1 * rar1(1) + adr2 * rar1(2)))
        vl(i, 2) = exp(ci * (adr1 * rar2(1) + adr2 * rar2(2)))
        i = i + 1
        s2 = s2 + 1.0_sp
      end do
      s1 = s1 + 1.0_sp
    end do

    ! Rotate each vector v and find to which v it becomes equivalent in
    ! terms of the quantities vl. This equivalence means belonging to the
    ! same sublattice.
    ainc = 2.0_sp * pi / real(ideg, sp)
    do i = 1, nl
      adr = v(i, 1)
      ang = v(i, 2)
      do k = 1, ideg
        ang = ang + ainc
        cang = cos(ang)
        sang = sin(ang)
        a = rar1(1) * cang + rar1(2) * sang
        b = rar2(1) * cang + rar2(2) * sang
        vla = exp(ci * adr * a)
        vlb = exp(ci * adr * b)
        do j = 1, nl
          test = abs(vla - vl(j, 1)) + abs(vlb - vl(j, 2))
          if (test - 5.0_sp * epsd <= 0.0_sp) exit
        end do
        jjs(i, k) = j
      end do
    end do
  end subroutine sublattice_equivalences
end module bravais_layer
