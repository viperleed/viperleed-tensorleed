!> \file angular_coefficients.f90
!! \brief Clebsch-Gordan and Gaunt coefficients, factorials and index tables.
!!
!! The setup routines of ref-calc that depend on LMAX only. Each routine
!! keeps the precision of every quantity of the legacy routine it follows:
!! some work in single precision on purpose, and the results are bit-for-bit
!! those of the legacy code (tests/tier1/test_refcalc_routines.py).

!> Angular-momentum coefficients and (L,M) index tables.
module angular_coefficients
  use constants, only: pi
  use kinds, only: sp, dp
  use leed_output, only: stop_with_error
  implicit none

  private

  public :: composite_layer_coefficients, matrix_coefficients

  !> The single-precision literal 3.14159265 as a double precision number,
  !! which is what the legacy REAL*8 variables PI of BLM and CELMG held.
  real(dp), parameter :: pi_dp = real(pi, dp)

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

  !> Integral of the product of three spherical harmonics (legacy BLM).
  !!
  !! Each harmonic can be expressed as a prefactor times a Legendre
  !! function. The three prefactors are lumped together as factor 'C', and
  !! the integral of the three Legendre functions follows Gaunt's summation
  !! scheme, set out by Slater, Atomic Structure, vol. 1, 309, 310.
  !! Author: Pendry.
  !!
  !! fac(n + 1) must hold n! (set by matrix_coefficients).
  function harmonics_product_integral(l1, m1, l2, m2, l3, m3, lmax, fac) &
      result(blm)
    integer, intent(in) :: l1, m1, l2, m2, l3, m3
    integer, intent(in) :: lmax  !< largest angular momentum l
    real(dp), intent(in) :: fac(:)
    real(sp) :: blm

    integer :: nl1, nl2, nl3, nm1, nm2, nm3, ic, ntemp
    integer :: is, ia1, ia2, ia3, ia4, ia5, ia6, ia7, ia8, ia9
    integer :: ib1, ib2, ib3, ib4, ib5, it1, it2, it
    integer :: ic1, ic2, ic3, ic4, ic5, ic6
    real(dp) :: a, b, bn, c, cn, alt_sign
    character(len=32) :: arguments

    if (m1 + m2 + m3 /= 0) then
      blm = 0.0_sp
      return
    end if
    if (l1 > lmax + lmax .or. l2 > lmax .or. l3 > lmax &
        .or. l1 < abs(m1) .or. l2 < abs(m2) .or. l3 < abs(m3)) then
      ! The arguments of matrix_coefficients stay within these bounds for
      ! every lmax, so no input reaches this branch.
      write (arguments, '(6(i3, ","))') l1, m1, l2, m2, l3, m3
      call stop_with_error("INVALID ARGUMENTS FOR BLM  " // trim(arguments))
    end if
    if (mod(l1 + l2 + l3, 2) /= 0) then
      blm = 0.0_sp
      return
    end if

    nl1 = l1
    nl2 = l2
    nl3 = l3
    nm1 = abs(m1)
    nm2 = abs(m2)
    nm3 = abs(m3)
    ic = (nm1 + nm2 + nm3) / 2
    if (max(nm1, nm2, nm3) > nm1) then
      if (max(nm2, nm3) <= nm2) then
        nl1 = l2
        nl2 = l1
        nm1 = nm2
        nm2 = abs(m1)
      else
        nl1 = l3
        nl3 = l1
        nm1 = nm3
        nm3 = abs(m1)
      end if
    end if
    if (nl2 < nl3) then
      ntemp = nl2
      nl2 = nl3
      nl3 = ntemp
      ntemp = nm2
      nm2 = nm3
      nm3 = ntemp
    end if
    if (nl3 < abs(nl2 - nl1)) then
      blm = 0.0_sp
      return
    end if

    ! Factor 'A'
    is = (nl1 + nl2 + nl3) / 2
    ia1 = is - nl2 - nm3
    ia2 = nl2 + nm2
    ia3 = nl2 - nm2
    ia4 = nl3 + nm3
    ia5 = nl1 + nl2 - nl3
    ia6 = is - nl1
    ia7 = is - nl2
    ia8 = is - nl3
    ia9 = nl1 + nl2 + nl3 + 1
    a = ((-1.0_dp)**ia1) / fac(ia3 + 1) * fac(ia2 + 1) / fac(ia6 + 1) &
        * fac(ia4 + 1)
    a = a / fac(ia7 + 1) * fac(ia5 + 1) / fac(ia8 + 1) * fac(is + 1) &
        / fac(ia9 + 1)

    ! Sum 'B'
    ib1 = nl1 + nm1
    ib2 = nl2 + nl3 - nm1
    ib3 = nl1 - nm1
    ib4 = nl2 - nl3 + nm1
    ib5 = nl3 - nm3
    it1 = max(0, -ib4) + 1
    it2 = min(ib2, ib3, ib5) + 1
    b = 0.0_dp
    alt_sign = (-1.0_dp)**it1
    ib1 = ib1 + it1 - 2
    ib2 = ib2 - it1 + 2
    ib3 = ib3 - it1 + 2
    ib4 = ib4 + it1 - 2
    ib5 = ib5 - it1 + 2
    do it = it1, it2
      alt_sign = -alt_sign
      ib1 = ib1 + 1
      ib2 = ib2 - 1
      ib3 = ib3 - 1
      ib4 = ib4 + 1
      ib5 = ib5 - 1
      bn = alt_sign / fac(it) * fac(ib1 + 1) / fac(ib3 + 1) * fac(ib2 + 1)
      bn = bn / fac(ib4 + 1) / fac(ib5 + 1)
      b = b + bn
    end do

    ! Factor 'C'
    ic1 = nl1 - nm1
    ic2 = nl1 + nm1
    ic3 = nl2 - nm2
    ic4 = nl2 + nm2
    ic5 = nl3 - nm3
    ic6 = nl3 + nm3
    cn = real((2 * nl1 + 1) * (2 * nl2 + 1) * (2 * nl3 + 1), dp) / pi_dp
    c = cn / fac(ic2 + 1) * fac(ic1 + 1) / fac(ic4 + 1) * fac(ic3 + 1)
    c = c / fac(ic6 + 1) * fac(ic5 + 1)
    c = sqrt(c) / 2.0_dp

    blm = real(((-1.0_dp)**ic) * a * b * c, sp)
  end function harmonics_product_integral

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

  !> Clebsch-Gordan type coefficients clm (CLM and ELM of the legacy code)
  !! for the matrices X and TAU (legacy CELMG).
  !!
  !! The non-zero values are tabulated first for (l2+m2) and (l3+m3) odd,
  !! and then for (l2+m2) and (l3+m3) even, the same scheme by which they
  !! are retrieved in the routines XM and XMT. Author: Pendry.
  !!
  !! On return fac(n + 1) holds n! for n = 0, ..., 4 * lmax + 1; the other
  !! elements are not touched. fac needs at least 4 * lmax + 2 elements.
  subroutine matrix_coefficients(lmax, clm, fac)
    integer, intent(in) :: lmax  !< largest angular momentum l
    real(sp), intent(out) :: clm(:)
    real(dp), intent(inout) :: fac(:)

    integer :: l2max, nf, i, n, ll, lm, ln, lo, l, m, lm2, lm3, lp, lq, k, ii
    integer :: il2, i2, il3, i3, l1, l2, l3, m1, m2, m3, la1, lb1, la11, lb11
    integer :: l11, phase_exponent, ylm_index
    real(dp) :: cl, a, b, asg, cm, alm
    real(sp), allocatable :: ylm(:), fac2(:), fac1(:)

    l2max = lmax + lmax
    nf = 4 * lmax + 1
    n = 2 * lmax + 1
    allocate (ylm(n * n), fac2(n * n), fac1(n))

    fac(1) = 1.0_dp
    do i = 1, nf
      fac(i + 1) = real(i, dp) * fac(i)
    end do

    ! The array ylm is first loaded with spherical harmonics, arguments
    ! theta = pi/2, phi = 0.
    lm = 0
    cl = 0.0_dp
    a = 1.0_dp
    b = 1.0_dp
    asg = 1.0_dp
    ll = l2max + 1
    ! Multiplicative factors required
    do l = 1, ll
      fac1(l) = asg * sqrt((2.0_dp * cl + 1.0_dp) * a &
                           / (4.0_dp * pi_dp * b * b))
      cm = -cl
      ln = l + l - 1
      do m = 1, ln
        lo = lm + m
        fac2(lo) = sqrt((cl + 1.0_dp + cm) * (cl + 1.0_dp - cm) &
                        / ((2.0_dp * cl + 3.0_dp) * (2.0_dp * cl + 1.0_dp)))
        cm = cm + 1.0_dp
      end do
      cl = cl + 1.0_dp
      a = a * 2.0_dp * cl * (2.0_dp * cl - 1.0_dp) / 4.0_dp
      b = b * cl
      asg = -asg
      lm = lm + ln
    end do

    ! First all the ylm for m = +-l and m = +-(l-1) are calculated by
    ! explicit formulae.
    lm = 1
    asg = -1.0_dp
    ylm(1) = fac1(1)
    do l = 1, l2max
      ln = lm + l + l + 1
      ylm(ln) = fac1(l + 1)
      ylm(lm + 1) = asg * fac1(l + 1)
      ylm(ln - 1) = 0.0_sp
      ylm(lm + 2) = 0.0_sp
      asg = -asg
      lm = ln
    end do

    ! Using ylm and yl(m-1) in a recurrence relation, yl(m+1) is calculated.
    lm = 1
    ll = l2max - 1
    do l = 1, ll
      ln = l + l - 1
      lm2 = lm + ln + 4
      lm3 = lm - ln
      do m = 1, ln
        lo = lm2 + m
        lp = lm3 + m
        lq = lm + m + 1
        ylm(lo) = -(fac2(lp) * ylm(lp)) / fac2(lq)
      end do
      lm = lm + l + l + 1
    end do

    k = 1
    do ii = 0, 1
      ll = lmax + ii
      do il2 = 1, ll
        l2 = il2 - ii
        m2 = -l2 + 1 - ii
        do i2 = 1, il2
          do il3 = 1, ll
            l3 = il3 - ii
            m3 = -l3 + 1 - ii
            do i3 = 1, il3
              la1 = max(abs(l2 - l3), abs(m2 - m3))
              lb1 = l2 + l3
              la11 = la1 + 1
              lb11 = lb1 + 1
              m1 = m3 - m2
              do l11 = la11, lb11, 2
                l1 = la11 + lb11 - l11 - 1
                phase_exponent = (l3 - l2 - l1) / 2 + m3
                ylm_index = l1 * (l1 + 1) - m1 + 1
                alm = ((-1.0_dp)**phase_exponent) * 4.0_dp * pi_dp &
                      * harmonics_product_integral(l1, m1, l2, m2, l3, -m3, &
                                                   lmax, fac)
                clm(k) = ylm(ylm_index) * alm
                k = k + 1
              end do
              m3 = m3 + 2
            end do
          end do
          m2 = m2 + 2
        end do
      end do
    end do
  end subroutine matrix_coefficients

end module angular_coefficients
