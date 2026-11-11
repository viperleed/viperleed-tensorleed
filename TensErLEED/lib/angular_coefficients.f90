!> \file angular_coefficients.f90
!! \brief Clebsch-Gordan and Gaunt coefficients, factorials and index tables.
!!
!! The setup routines of ref-calc that depend on LMAX only. Each routine
!! keeps the precision of every quantity of the legacy routine it follows:
!! some work in single precision on purpose, and the results are bit-for-bit
!! those of the legacy code (tests/tier1/test_refcalc_routines.py).

!> Angular-momentum coefficients and (L,M) index tables.
module angular_coefficients
  use kinds, only: sp, dp
  implicit none

  private

  public :: composite_layer_coefficients

contains

  !> Factorial(l) / 10**l, using an asymptotic expansion for l > 4
  !! (legacy FACT). Used by composite_layer_coefficient.
  !!
  !! The expansion works extremely well, much higher than l = 15. As in
  !! the legacy routine, the result is undefined for l < 0.
  function scaled_factorial(l) result(fact)
    integer, intent(in) :: l
    real(dp) :: fact

    real(dp) :: x

    if (l > 4) then
      x = real(l + 1, dp)
      fact = exp(-x) * (10.0_dp**((x - 0.5_dp) * log10(x) - (x - 1.0_dp))) &
             * sqrt(6.283185307179586_dp) &
             * (1.0_dp + (1.0_dp / (12.0_dp * x)) &
                + (1.0_dp / (288.0_dp * (x**2))) &
                - (139.0_dp / (51840.0_dp * (x**3))) &
                - (571.0_dp / (2488320.0_dp * (x**4))))
      return
    end if
    ! The single-precision literals are converted to dp as the legacy
    ! routine does; 6.0 * 0.001 and 24.0 * 0.0001 are products in sp.
    if (l == 0) fact = 1.0_dp
    if (l == 1) fact = real(0.1_sp, dp)
    if (l == 2) fact = real(0.02_sp, dp)
    if (l == 3) fact = real(6.0_sp * 0.001_sp, dp)
    if (l == 4) fact = real(24.0_sp * 0.0001_sp, dp)
  end function scaled_factorial

  !> Clebsch-Gordan type coefficient of the composite-layer propagators for
  !! one combination of (l, m) (legacy CA); the same purpose for
  !! composite_layer_coefficients as harmonics_product_integral has for
  !! matrix_coefficients.
  !!
  !! The legacy routine mixes precisions on purpose: delt, total and yint
  !! are sp, the prefactors dp, and 12.56637 is a single-precision value of
  !! 4 pi that is not the sp pi times four.
  function composite_layer_coefficient(l1, ma1, l2, m2, l3, ma3) result(ca)
    integer, intent(in) :: l1, ma1, l2, m2, l3, ma3
    real(sp) :: ca

    integer :: m1, m3, igt, ilst, num, it, ita
    real(sp) :: delt, total, yint
    real(dp) :: pref, prefa, prefb, preft

    m1 = -ma1
    m3 = -ma3
    if (abs(m1) > l1 .or. abs(m2) > l2 .or. abs(m3) > l3 &
        .or. (m1 + m2) /= -m3 &
        .or. .not. (l3 >= abs(l1 - l2) .and. l3 <= (l1 + l2)) &
        .or. mod(l1 + l2 + l3, 2) /= 0) then
      yint = 0.0_sp
    else if (l1 == 0 .or. l2 == 0 .or. l3 == 0) then
      yint = ((2.0_sp * real(l1, sp) + 1.0_sp) &
              * (2.0_sp * real(l2, sp) + 1.0_sp) &
              * (2.0_sp * real(l3, sp) + 1.0_sp) / 12.56637_sp)**0.5_sp &
             * (-1.0_sp)**((abs(m1) + abs(m2) + abs(m3)) / 2) &
             / (real(abs(l1) + abs(l2) + abs(l3), sp) + 1.0_sp)
    else
      delt = real(scaled_factorial(l1 + l2 - l3) &
                  * scaled_factorial(l2 + l3 - l1) &
                  * scaled_factorial(l3 + l1 - l2) &
                  / scaled_factorial(l1 + l2 + l3 + 1) &
                  * ((-1)**(l1 - l2 - m3 + ((l1 + l2 + l3) / 2))), sp)
      pref = scaled_factorial(l1 + m1) * scaled_factorial(l1 - m1) &
             * scaled_factorial(l2 + m2) * scaled_factorial(l2 - m2) &
             * scaled_factorial(l3 - m3) * scaled_factorial(l3 + m3)
      prefa = sqrt(pref) * (10.0_dp**(l1 + l2 + l3))
      prefb = scaled_factorial((l1 + l2 + l3) / 2) &
              / (scaled_factorial((l1 + l2 - l3) / 2) &
                 * scaled_factorial((l2 + l3 - l1) / 2) &
                 * scaled_factorial((l3 + l1 - l2) / 2)) * 0.1_sp
      preft = delt * prefa * prefb
      igt = min(l2 + m2, l1 - m1, l1 + l2 - l3)
      ilst = max(l1 + m2 - l3, -m1 + l2 - l3, 0)
      num = igt - ilst + 1
      total = 0.0_sp
      do it = 1, num
        ita = ilst + it - 1
        total = total + ((-1.0_sp)**ita) &
                / (scaled_factorial(ita) &
                   * scaled_factorial(l3 - l2 + ita + m1) &
                   * scaled_factorial(l3 - l1 + ita - m2) &
                   * scaled_factorial(l1 + l2 - l3 - ita) &
                   * scaled_factorial(l1 - ita - m1) &
                   * scaled_factorial(l2 - ita + m2)) &
                * (10.0_dp**(-l1 - l2 - l3))
      end do
      yint = ((2.0_sp * real(l1, sp) + 1.0_sp) &
              * (2.0_sp * real(l2, sp) + 1.0_sp) &
              * (2.0_sp * real(l3, sp) + 1.0_sp) / 12.56637_sp)**0.5_sp &
             * preft * total
    end if
    ca = ((-1.0_sp)**(ma1 + ma3)) * yint
  end function composite_layer_coefficient

  !> Clebsch-Gordan coefficients for the composite layer, in the order in
  !! which the composite-layer propagators use them (legacy CAAA).
  !!
  !! caa must hold the coefficients: the legacy size NCA(lmax) of the
  !! driver. At lmax = 1 the legacy routine writes past that size; the
  !! behaviour is kept.
  subroutine composite_layer_coefficients(lmax, caa)
    integer, intent(in) :: lmax  !< largest angular momentum l
    real(sp), intent(out) :: caa(:)

    integer :: lmmax, i, j, ii, l1, l2, m1, m2, m3, il, im, lmin, lsum, ila
    integer :: la

    lmmax = (lmax + 1)**2
    ii = 1
    do i = 1, lmmax
      l1 = int(sqrt(real(i - 1, sp) + 0.00001_sp))
      m1 = i - l1 - l1 * l1 - 1
      do j = 1, lmmax
        l2 = int(sqrt(real(j - 1, sp) + 0.00001_sp))
        m2 = j - l2 - l2 * l2 - 1
        m3 = m2 - m1
        il = abs(l1 - l2)
        im = abs(m3)
        lmin = max(il, im + mod(il + im, 2))
        lsum = l1 + l2
        lmin = lmin + 1
        lsum = lsum + 1
        if (i == j) then
          if (m1 > 0) cycle
        else
          if (l1 < l2 .or. (l1 == l2 .and. abs(m1) < abs(m2))) cycle
        end if
        do ila = lmin, lsum, 2
          la = ila - 1
          caa(ii) = composite_layer_coefficient(l1, m1, l2, m2, la, m3)
          ii = ii + 1
        end do
      end do
    end do
  end subroutine composite_layer_coefficients

end module angular_coefficients
