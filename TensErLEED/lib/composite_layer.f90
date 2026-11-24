!> \file composite_layer.f90
!! \brief Composite layers (several subplanes): ordering of the subplanes
!! and the propagators between them.
!!
!! Each routine keeps the precision of every quantity of the legacy routine
!! it follows, and its results are bit-for-bit those of the legacy code
!! (tests/tier1/test_refcalc_routines.py).

!> Subplanes of a composite layer and the propagators between them.
module composite_layer
  use kinds, only: sp
  implicit none

  private

  public :: sort_subplanes

contains

  !> Orders the subplanes of a composite layer by increasing position along
  !! +x (legacy SRTLAY); that is from the vacuum side towards the solid.
  !! Author: Pendry.
  !!
  !! The chemical-element assignment is reordered along. Subplanes at equal
  !! position are not permuted. The topmost subplane becomes the origin of
  !! the layer. Also finds the interplanar vectors and which of the
  !! propagators between subplanes are equal to an earlier one.
  !!
  !! The arrays are explicitly shaped; the driver passes its larger arrays
  !! whole, so that the leading part of their storage holds the packed
  !! (nlay, ...) arrays.
  subroutine sort_subplanes(pos, poss, lps, lpss, mgh, nlay, drl, nlay2, &
                            nugh, ngeq, nord)
    !> Number of subplanes.
    integer, intent(in) :: nlay
    !> nlay * (nlay - 1) / 2.
    integer, intent(in) :: nlay2
    !> Subplane positions (atomic positions in the unit cell).
    real(sp), intent(in) :: pos(nlay, 3)
    !> Reordered pos; the first x is zero.
    real(sp), intent(out) :: poss(nlay, 3)
    !> Chemical element assignment of each subplane.
    integer, intent(in) :: lps(nlay)
    !> Reordered lps.
    integer, intent(out) :: lpss(nlay)
    !> Organization of the propagators: mgh(i, j) = k means that the
    !! propagator from subplane j to subplane i is the k-th block of the
    !! columnar matrix gh (nlay2 blocks for i < j, then nlay2 for i > j).
    integer, intent(out) :: mgh(nlay, nlay)
    !> Interplanar vectors.
    real(sp), intent(out) :: drl(nlay2, 3)
    !> nugh(k) = 1: propagator k must be computed; 0: it is equal to an
    !! earlier one (see ngeq).
    integer, intent(out) :: nugh(nlay2)
    !> ngeq(k) = m: propagator k can be copied from the new propagator m;
    !! 0: no copying.
    integer, intent(out) :: ngeq(nlay2)
    !> Order of the subplanes, to trace the reordering back.
    integer, intent(out) :: nord(nlay)

    integer :: i, j, k, m, in, km, lpsa, norda
    real(sp) :: pm, posa(3)

    do i = 1, nlay
      lpss(i) = lps(i)
      nord(i) = i
      poss(i, :) = pos(i, :)
    end do

    ! Analyse the order of the subplane positions along the x axis and
    ! reorder them in ascending position
    do i = 1, nlay - 1
      km = i
      pm = poss(i, 1)
      do k = i + 1, nlay
        if (poss(k, 1) >= pm) cycle
        pm = poss(k, 1)
        km = k
      end do
      if (km == i) cycle
      posa = poss(km, :)
      lpsa = lpss(km)
      norda = nord(km)
      do k = km, i + 1, -1
        poss(k, :) = poss(k - 1, :)
        nord(k) = nord(k - 1)
        lpss(k) = lpss(k - 1)
      end do
      poss(i, :) = posa
      lpss(i) = lpsa
      nord(i) = norda
    end do

    ! Set the origin of the composite layer to the topmost subplane
    do i = 2, nlay
      poss(i, 1) = poss(i, 1) - poss(1, 1)
    end do
    poss(1, 1) = 0.0_sp

    ! Interplanar vectors, organization of gh, and equal interplanar
    ! vectors
    in = 1
    do i = 1, nlay
      do j = i, nlay
        mgh(i, j) = 0
        mgh(j, i) = 0
        if (i == j) cycle
        drl(in, :) = poss(i, :) - poss(j, :)
        mgh(i, j) = in
        mgh(j, i) = in + nlay2
        nugh(in) = 1
        ngeq(in) = 0
        do m = 1, in - 1
          posa = abs(drl(in, :) - drl(m, :))
          if ((posa(1) + posa(2) + posa(3)) - 0.001_sp <= 0.0_sp) then
            ngeq(in) = m
            nugh(in) = 0
            exit
          end if
        end do
        in = in + 1
      end do
    end do
  end subroutine sort_subplanes
end module composite_layer
