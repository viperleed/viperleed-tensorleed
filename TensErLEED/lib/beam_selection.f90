!> \file beam_selection.f90
!! \brief Selection of the beams in use at an energy.
!!
!! The routine keeps the precision of every quantity of the legacy routine
!! it follows, and its results are bit-for-bit those of the legacy code
!! (tests/tier1/test_refcalc_routines.py).

!> Beams in use at an energy.
module beam_selection
  use kinds, only: sp
  implicit none

  private

  public :: select_beams

contains

  !> Selects the beams of the input beam list that are needed at the current
  !! energy (legacy BEAMS).
  !!
  !! A beam is needed if its out-of-plane wave vector is real or has an
  !! imaginary part not larger than sqrt(tst), i.e., if it is not evanescent
  !! beyond the limit tst, which must already include the interlayer spacing.
  !! The selected beams of each beam set are stored one after the other in
  !! pq and pqf. The positions of the beams to be output are translated into
  !! positions in the selected list (0 for a beam that is not selected).
  !!
  !! The dummies are explicitly shaped like those of the legacy routine; the
  !! actual arrays of the driver are larger than the beam list, and their
  !! leading part is passed.
  subroutine select_beams(knbs, knb, spq, spqf, knt, npu, npun, ak2, ak3, e, &
                          tst, nb, pq, pqf, npuc, mpu, nt, np)
    !> Number of input beam sets.
    integer, intent(in) :: knbs
    !> Number of beams in each input beam set.
    integer, intent(in) :: knb(knbs)
    !> Total number of input beams.
    integer, intent(in) :: knt
    !> Input beams in atomic units (reciprocal lattice vectors G).
    real(sp), intent(in) :: spq(2, knt)
    !> Input beams in units of the reciprocal substrate vectors.
    real(sp), intent(in) :: spqf(2, knt)
    !> Number of beams of which intensities are output.
    integer, intent(in) :: npun
    !> Positions in the input beam list of the beams to be output.
    integer, intent(in) :: npu(npun)
    !> Components of the incident wave vector parallel to the surface.
    real(sp), intent(in) :: ak2, ak3
    !> Current energy.
    real(sp), intent(in) :: e
    !> Criterion for the selection of beams.
    real(sp), intent(in) :: tst
    !> Number of selected beams in each beam set.
    integer, intent(out) :: nb(knbs)
    !> Selected beams, as spq; only the first nt are set.
    real(sp), intent(inout) :: pq(2, knt)
    !> Selected beams, as spqf; only the first nt are set.
    real(sp), intent(inout) :: pqf(2, knt)
    !> Position in the selected list of the beams to be output; 0 for a beam
    !! that is not selected. Only the first mpu are set.
    integer, intent(inout) :: npuc(npun)
    !> Number of beams to be output that were reached in the input list.
    integer, intent(out) :: mpu
    !> Total number of selected beams.
    integer, intent(out) :: nt
    !> Largest number of beams selected in any of the beam sets.
    integer, intent(out) :: np

    integer :: knbj, j, n, k, kk, mpu1, i

    knbj = 0
    nt = 0
    np = 0
    mpu = 0
    do j = 1, knbs
      n = knb(j)
      nb(j) = 0
      do k = 1, n
        kk = k + knbj
        mpu1 = mpu + 1
        ! Fortran conditions are not evaluated lazily: npu(mpu1) is only
        ! read once mpu1 <= npun is known
        if (mpu1 <= npun) then
          if (kk == npu(mpu1)) then
            mpu = mpu + 1
            npuc(mpu) = 0
          end if
        end if
        ! beam propagates
        if (2.0_sp * e - (ak2 + spq(1, kk))**2 - (ak3 + spq(2, kk))**2 &
            >= -tst) then
          nb(j) = nb(j) + 1
          nt = nt + 1
          if (mpu1 <= npun) then
            if (kk == npu(mpu1)) npuc(mpu) = nt
          end if
          do i = 1, 2
            pq(i, nt) = spq(i, kk)
            pqf(i, nt) = spqf(i, kk)
          end do
        end if
      end do
      knbj = knbj + knb(j)
      np = max(np, nb(j))
    end do
  end subroutine select_beams
end module beam_selection
