!> \file atomic_scattering.f90
!! \brief Atomic scattering: tables for the temperature-dependent phase
!! shifts.
!!
!! As in angular_coefficients, every routine keeps the precision of each
!! quantity of the legacy routine it follows, and its results are bit-for-bit
!! those of the legacy code (tests/tier1/test_refcalc_routines.py).

!> Tables for the atomic t-matrices.
module atomic_scattering
  use kinds, only: sp, dp
  implicit none

  private

  public :: legendre_product_integrals

contains

  !> Tabulates ppp(i1, i2, i3), each element containing the integral of the
  !! product of three Legendre functions P(i1), P(i2), P(i3) (legacy CPPP).
  !!
  !! The integrals follow Gaunt's summation scheme, set out by Slater,
  !! Atomic Structure, vol. 1, 309, 310. They are used to compute
  !! temperature-dependent phase shifts. Author: Pendry.
  !!
  !! The dummy ppp is explicitly shaped on purpose: ref-calc declares ppp
  !! with its maximum size and the readers of the table take the first
  !! n1 * n2 * n3 elements as an array of shape (n1, n2, n3).
  !! The factorial table f is a work array, and the same as the one that
  !! matrix_coefficients fills: f(n + 1) is set to n!, in double precision
  !! (REAL*8 on purpose: it avoids a divergence of the sums to infinity).
  !! The largest n1 + n2 + n3 - 1 = 4 * lmax + 2 needs size(f) >= 4 * lmax + 3.
  subroutine legendre_product_integrals(n1, n2, n3, ppp, f)
    integer, intent(in) :: n1, n2, n3
    real(sp), intent(out) :: ppp(n1, n2, n3)
    real(dp), intent(out) :: f(:)

    integer :: i, i1, i2, i3, im1, im2, im3, is, it, j
    integer :: ia1, ia2, ia3, ia4, ia5
    real(sp) :: a, sign_
    real(dp) :: sum_

    f(1) = 1.0_dp
    do i = 1, size(f) - 1
      f(i + 1) = f(i) * real(i, dp)
    end do

    do i3 = 1, n3
      do i2 = 1, n2
        do i1 = 1, n1
          ! Order the three indices: im1 >= im2 >= im3
          im1 = i1
          im2 = i2
          im3 = i3
          if (i1 < i2) then
            im1 = i2
            im2 = i1
          end if
          if (im2 < i3) then
            im3 = im2
            im2 = i3
            if (im1 < im2) then
              j = im1
              im1 = im2
              im2 = j
            end if
          end if

          a = 0.0_sp
          is = i1 + i2 + i3 - 3
          if (mod(is, 2) /= 1 .and. abs(im2 - im1) - im3 + 1 <= 0) then
            sum_ = 0.0_dp
            is = is / 2
            sign_ = 1.0_sp
            do it = 1, im3
              sign_ = -sign_
              ia1 = im1 + it - 1
              ia2 = im1 - it + 1
              ia3 = im3 - it + 1
              ia4 = im2 + im3 - it
              ia5 = im2 - im3 + it
              sum_ = sum_ - sign_ * f(ia1) * f(ia4) &
                            / (f(ia2) * f(ia3) * f(ia5) * f(it))
            end do
            ia1 = 2 + is - im1
            ia2 = 2 + is - im2
            ia3 = 2 + is - im3
            ia4 = 3 + 2 * (is - im3)
            a = -(-1.0_sp)**(is - im2) * f(ia4) * f(is + 1) * f(im3) * sum_ &
                / (f(ia1) * f(ia2) * f(ia3) * f(2 * is + 2))
          end if
          ppp(i1, i2, i3) = a
        end do
      end do
    end do
  end subroutine legendre_product_integrals
end module atomic_scattering
