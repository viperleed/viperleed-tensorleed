!> \file rfacsb.f90
!! \brief Modules of the R-factor program rfactor, in dependency order.

!> Kind parameters for reals and complexes.
!!
!! sp is the precision of the legacy default REAL (REAL(4)); dp is double
!! precision (REAL(8)). Every real and complex is declared with one of these,
!! so that a change of precision is a change of one line here.
module kinds
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none

  private

  public :: sp, dp

  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
end module kinds

!> Interpolation of intensity curves and the energy working grid.
!!
!! yval interpolates a tabulated curve at one abscissa: 3rd-order polynomial
!! interpolation through the four nearest points, or 2nd/1st order if the
!! curve has fewer than four points. The four points of the last interval
!! found are kept in COMMON /DATBLK/ and reused when the next call finds the
!! same interval (itil, itih).
module interpolation
  use kinds, only: sp
  implicit none

  private

  public :: grid, intpol, yval

contains

  !> Finds the interpolation interval [worx(il), worx(ih)] containing x by
  !! successive halving. Outside the table, the first or last interval.
  subroutine binsrx(il, ih, x, worx, length)
    integer, intent(out) :: il, ih
    real(sp), intent(in) :: x
    integer, intent(in) :: length
    real(sp), intent(in) :: worx(length)

    integer :: i, ihi

    i = (length + 1) / 2
    ihi = (i + 1) / 2
    if (x < worx(1) .or. x >= worx(length)) then
      if (x < worx(1)) then
        il = 1
        ih = 2
      else
        il = length - 1
        ih = length
      end if
      return
    end if
    do
      if (i < 1) then
        i = ihi + i
      else if (i >= length) then
        i = i - ihi
      else if (x >= worx(i) .and. x <= worx(i + 1)) then
        il = i
        ih = i + 1
        return
      else if (x < worx(i)) then
        i = i - ihi
      else
        i = ihi + i
      end if
      ihi = (ihi + 1) / 2
    end do
  end subroutine binsrx

  !> Stores the four grid points nearest to the interval (il, ih) and their
  !! ordinates in /DATBLK/, shifted inwards at the ends of the table.
  subroutine stfpts(il, ih, worx, wory, length)
    integer, intent(in) :: il, ih, length
    real(sp), intent(in) :: worx(length), wory(length)

    real(sp) :: x0, x1, x2, x3, y0, y1, y2, y3
    common /datblk/ x0, x1, x2, x3, y0, y1, y2, y3
    integer :: i

    i = il - 1
    if (il <= 1) i = i + 1
    if (ih >= length) i = i - 1
    x0 = worx(i)
    x1 = worx(i + 1)
    x2 = worx(i + 2)
    x3 = worx(i + 3)
    y0 = wory(i)
    y1 = wory(i + 1)
    y2 = wory(i + 2)
    y3 = wory(i + 3)
  end subroutine stfpts

  !> 3rd-order (Newton) polynomial through the four points in /DATBLK/,
  !! evaluated at x.
  real(sp) function xnterp(x)
    real(sp), intent(in) :: x

    real(sp) :: x0, x1, x2, x3, y0, y1, y2, y3
    common /datblk/ x0, x1, x2, x3, y0, y1, y2, y3
    real(sp) :: term, fact1, fact2, temp

    term = y0
    fact1 = x - x0
    fact2 = (y1 - y0) / (x1 - x0)
    term = term + fact1 * fact2
    fact1 = fact1 * (x - x1)
    fact2 = ((y2 - y1)/(x2 - x1) - fact2) / (x2 - x0)
    term = term + fact1 * fact2
    fact1 = fact1 * (x - x2)
    temp = ((y3 - y2)/(x3 - x2) - (y2 - y1)/(x2 - x1))/(x3 - x1)
    fact2 = (temp - fact2) / (x3 - x0)
    term = term + fact1 * fact2
    xnterp = term
  end function xnterp

  !> 2nd-order polynomial through the first three points (xs, ys), or
  !! 1st order through the first two if n <= 2, evaluated at x.
  subroutine xntrp2(x, y, xs, ys, n)
    real(sp), intent(in) :: x
    real(sp), intent(out) :: y
    integer, intent(in) :: n
    real(sp), intent(in) :: xs(n), ys(n)

    real(sp) :: a, b, c

    if (n <= 2) then
      y = (ys(2) - ys(1))*(x - xs(1))/(xs(2) - xs(1)) + ys(1)
      return
    end if
    a = (ys(1) - ys(2))/(xs(1) - xs(2))/(xs(2) - xs(3)) - &
        (ys(1) - ys(3))/(xs(1) - xs(3))/(xs(2) - xs(3))
    b = (ys(1) - ys(2))/(xs(1) - xs(2)) - a*(xs(1) + xs(2))
    c = ys(1) - (a*xs(1) + b)*xs(1)
    y = c + x*(b + a*x)
  end subroutine xntrp2

  !> Interpolates the curve (worx, wory) at x. itil, itih are the interval
  !! of the previous call (0 to force a new one); they are updated.
  real(sp) function yval(x, wory, worx, length, itil, itih)
    real(sp), intent(in) :: x
    integer, intent(in) :: length
    real(sp), intent(in) :: wory(length), worx(length)
    integer, intent(inout) :: itil, itih

    real(sp) :: y
    integer :: il, ih

    if (length < 4) then
      call xntrp2(x, y, worx, wory, length)
      yval = y
      return
    end if
    call binsrx(il, ih, x, worx, length)
    if (il /= itil .or. ih /= itih) call stfpts(il, ih, worx, wory, length)
    y = xnterp(x)
    itih = ih
    itil = il
    yval = y
  end function yval

  !> Energy working grid (steps of eincr) for one beam. Drops the leading
  !! points with zero intensity from (a, e) and shortens ne accordingly;
  !! ne = 0 if at most one point is left. egrid(1:npts) is the grid and
  !! imin the index of the first nonzero intensity on entry.
  subroutine grid(a, ne, e, ngp, eincr, egrid, npts, imin)
    integer, intent(in) :: ngp
    real(sp), intent(inout) :: a(ngp), e(ngp)
    integer, intent(inout) :: ne
    real(sp), intent(in) :: eincr
    real(sp), intent(inout) :: egrid(ngp)
    integer, intent(inout) :: npts, imin

    real(sp) :: xmin
    integer :: lmin, lmax, i, ie

10  format('0 IN THIS BEAM TOO FEW ENERGY VALUES FOR INTERPOLATION')

    do ie = 1, ne
      imin = ie
      if (a(ie) > 1.e-8_sp) exit
    end do

    if (imin == ne) then
      ne = 0
      write(6, 10)
      return
    end if

    lmin = int((e(imin) - 0.0001_sp) / eincr) + 1
    lmin = max(lmin, 0)
    xmin = real(lmin, sp) * eincr
    lmax = int((e(ne) + 0.0001_sp) / eincr)

    do i = imin, ne
      e(i - imin + 1) = e(i)
      a(i - imin + 1) = a(i)
    end do
    ne = ne - imin + 1

    npts = lmax - lmin + 1

    egrid(1) = xmin
    do ie = 2, npts
      egrid(ie) = egrid(ie - 1) + eincr
    end do
  end subroutine grid

  !> Interpolates the intensities a(e) of beam ib onto the working grid
  !! (steps of eincr), replacing (a, e, ne); negative values become zero.
  !! Leading zero intensities are dropped; ne = 0 if at most one point is
  !! left. x and woryt are work space; ipr = 2 prints the result.
  subroutine intpol(a, ne, e, ngp, eincr, ipr, x, woryt, ib)
    integer, intent(in) :: ngp
    real(sp), intent(inout) :: a(ngp), e(ngp)
    integer, intent(inout) :: ne
    real(sp), intent(in) :: eincr
    integer, intent(in) :: ipr, ib
    real(sp), intent(inout) :: x(ngp), woryt(ngp)

    real(sp) :: xmin, xval
    integer :: itil, itih, imin, lmin, lmax, nem, npts, i, ie

10  format('0 IN THIS BEAM TOO FEW ENERGY VALUES FOR INTERPOLATION')
20  format('*INTENSITIES AFTER INTERPOLATION IN BEAM', 1I4, &
           /, 100(5(1F7.2, 1E13.4, 3X), /))

    itil = 0
    itih = 0

    do ie = 1, ne
      imin = ie
      if (a(ie) > 1.e-8_sp) exit
    end do

    if (imin == ne) then
      ne = 0
      write(6, 10)
      return
    end if

    lmin = int((e(imin) - 0.0001_sp) / eincr) + 1
    lmin = max(lmin, 0)
    xmin = real(lmin, sp) * eincr
    lmax = int((e(ne) + 0.0001_sp) / eincr)

    do i = imin, ne
      x(i - imin + 1) = e(i)
      woryt(i - imin + 1) = a(i)
    end do
    nem = ne - imin + 1

    npts = lmax - lmin + 1
    ne = npts

    xval = xmin - eincr

    do i = 1, npts
      xval = xval + eincr
      a(i) = yval(xval, woryt(1:nem), x(1:nem), nem, itil, itih)
      if (a(i) < 0.0_sp) then
        a(i) = 0.0_sp
      end if
    end do

    e(1) = xmin
    do ie = 2, npts
      e(ie) = e(ie - 1) + eincr
    end do

    if (ipr == 2 .and. ne > 0) then
      write(6, 20) ib, (e(ie), a(ie), ie = 1, ne)
    end if
  end subroutine intpol
end module interpolation

!> Derivatives, smoothing and the smooth-R-factor Y function of intensity
!! curves.
module derivatives
  use kinds, only: sp
  implicit none

  private

  public :: der, smooth, ysmooth

contains

  !> First derivative y1 of y (tabulated on ne points with step h), after
  !! Zanazzi and Jona, by ne:
  !! ne >= 23: 7-point central formula, 4-point one-sided at 3 end points;
  !! 3 <= ne <= 22: central differences, 3-point one-sided at end points;
  !! ne = 2: two-point difference; ne = 1: derivative set to zero.
  subroutine der(y, ne, y1, h)
    integer, intent(in) :: ne
    real(sp), intent(in) :: y(ne), h
    real(sp), intent(out) :: y1(ne)

    integer :: j, jf, m

    if (ne >= 23) then
      jf = ne - 3
      do j = 4, jf
        y1(j) = (y(j + 3) - 9._sp*y(j + 2) + 45._sp*y(j + 1) - 45._sp*y(j - 1) &
                 + 9._sp*y(j - 2) - y(j - 3)) / (60._sp*h)
      end do
      do j = 1, 3
        y1(j) = (2._sp*y(j + 3) - 9._sp*y(j + 2) + 18._sp*y(j + 1) - 11._sp*y(j)) &
                / (6._sp*h)
        m = ne - j + 1
        y1(m) = (11._sp*y(m) - 18._sp*y(m - 1) + 9._sp*y(m - 2) - 2._sp*y(m - 3)) &
                / (6._sp*h)
      end do
    else if (ne >= 3) then
      jf = ne - 1
      do j = 2, jf
        y1(j) = (y(j + 1) - y(j - 1)) / (2._sp*h)
      end do
      y1(1) = (-3._sp*y(1) + 4._sp*y(2) - y(3)) / (2._sp*h)
      y1(jf + 1) = (3._sp*y(jf + 1) - 4._sp*y(jf) + y(jf - 1)) / (2._sp*h)
    else if (ne == 1) then
      y1(1) = 0._sp
    else
      y1(1) = (y(2) - y(1)) / h
      y1(2) = y1(1)
    end if
  end subroutine der

  !> Smooths the curve (ee, ae) of n points ismoth times by weighted
  !! three-point averaging; on equal intervals a simplified formula is used.
  !! Curves of fewer than 3 points are left unchanged. ipr = 2 prints the
  !! result for beam ib.
  subroutine smooth(ismoth, ae, ee, n, ipr, ib)
    integer, intent(in) :: ismoth, n, ipr, ib
    real(sp), intent(inout) :: ae(n)
    real(sp), intent(in) :: ee(n)

    real(sp) :: am, af, e21, e32
    integer :: i, ie

10  format('*EXP. ENERG. AND INTENS. AFTER SMOOTHING IN BEAM', 1I3, &
           /, 2000(5(1F7.2, 1E13.4, 3X), /))

    do i = 1, ismoth
      if (n >= 3) then
        am = ae(1)
        do ie = 2, n - 1
          af = ae(ie)
          e21 = ee(ie) - ee(ie - 1)
          e32 = ee(ie + 1) - ee(ie)
          if (abs(e21 - e32) < 0.0001_sp) then
            ae(ie) = 0.25_sp * (2._sp*af + am + ae(ie + 1))
          else
            ae(ie) = 0.5_sp * (af + (e32*am + e21*ae(ie + 1))/(e21 + e32))
          end if
          am = af
        end do
      end if
    end do

    if (ipr == 2 .and. n > 0) then
      write(6, 10) ib, (ee(ie), ae(ie), ie = 1, n)
    end if
  end subroutine smooth

  !> Y function of the smooth R factor R_S (Imre et al., J. Phys.: Condens.
  !! Matter 38, 105001) for beams 1..nb of the curves a(ib, 1:ne(ib)) with
  !! first derivatives ap (per eV). The output y is 2*vi times Y_S of the
  !! paper; the constant factor cancels in R_S (same vi for experiment and
  !! theory). The second derivative is taken from three points of a (grid
  !! step eincr); at the ends from the nearest three, and zero for curves of
  !! fewer than 3 points. Entries beyond ne(ib) are not set. e is used only
  !! for printing (ipr >= 2).
  subroutine ysmooth(a, ap, nbd, ngp, nb, ne, e, y, vi, eincr, ipr)
    integer, intent(in) :: nbd, ngp, nb, ipr
    real(sp), intent(in) :: a(nbd, ngp), ap(nbd, ngp), e(nbd, ngp)
    integer, intent(in) :: ne(nbd)
    real(sp), intent(inout) :: y(nbd, ngp)
    real(sp), intent(in) :: vi, eincr

    real(sp), parameter :: alpha = 4.0_sp, beta = 0.15_sp
    real(sp) :: af, deriv1, deriv2, ydensq, y1, y2
    integer :: ib, ie, n

40  format(' SMOOTH Y FUNCTION IN BEAM', 1I5, /, 1000(5(1F7.2, 1E13.4, 3X), /))

    do ib = 1, nb
      n = ne(ib)
      if (n /= 0) then
        do ie = 1, n
          af = a(ib, ie)
          ! avoids division by 0 in I'/sqrt(I*I+...)
          if (af <= 1.0e-10_sp) then
            af = 1.0e-10_sp
          end if
          deriv1 = ap(ib, ie)*vi
          if (n < 3) then
            deriv2 = 0.0_sp
          else if (ie == 1) then
            deriv2 = a(ib, ie + 2) + a(ib, ie) - 2*a(ib, ie + 1)
          else if (ie < n) then
            deriv2 = a(ib, ie + 1) + a(ib, ie - 1) - 2*af
          else
            deriv2 = a(ib, ie) + a(ib, ie - 2) - 2*a(ib, ie - 1)
          end if
          ydensq = af*af + 4*deriv1*deriv1
          ! threshold 1e-10, not 0, avoids overflow in the division by I''
          if (deriv2 > 1.0e-10_sp) then
            deriv2 = deriv2*(vi/eincr)**2
            y1 = alpha*af/deriv2 - (0.5_sp*alpha)*deriv1*deriv1/(deriv2*deriv2) + beta
            if (y1 > 0) then
              y2 = y1/sqrt(1.0_sp + y1*y1)
              ydensq = ydensq + y2*y2*deriv2*deriv2
            end if
          end if
          y(ib, ie) = 2*deriv1/sqrt(ydensq)
        end do
      end if
    end do
    if (ipr >= 2) then
      do ib = 1, nb
        n = ne(ib)
        if (n /= 0) then
          write(6, 40) ib, (e(ib, ie), y(ib, ie), ie = 1, n)
        end if
      end do
    end if
  end subroutine ysmooth
end module derivatives

!> Pendry Y function, integrals over the energy range, and the Zanazzi-Jona
!! normalisation, as they enter the R factors.
module rfactor_terms
  use kinds, only: sp
  use interpolation, only: yval
  implicit none

  private

  public :: epszj, varsum, ypend

contains

  !> Largest |aep(ib, ie)| for ie = ne1..ne2 (eps of the Zanazzi-Jona
  !! R factor, aep being the derivative of the experimental curve).
  subroutine epszj(aep, nbd, ngp, ib, ne1, ne2, eps)
    integer, intent(in) :: nbd, ngp, ib, ne1, ne2
    real(sp), intent(in) :: aep(nbd, ngp)
    real(sp), intent(out) :: eps

    real(sp) :: a
    integer :: ie

    eps = 0._sp
    do ie = ne1, ne2
      a = abs(aep(ib, ie))
      if (a > eps) eps = a
    end do
  end subroutine epszj

  !> Integral s of y(1:n) (step h) by the trapezoid rule, after Zanazzi and
  !! Jona.
  subroutine trapezoid(y, n, h, s)
    integer, intent(in) :: n
    real(sp), intent(in) :: y(n), h
    real(sp), intent(out) :: s

    real(sp) :: a
    integer :: j

    a = 0._sp
    do j = 1, n
      a = a + y(j)
    end do
    s = a - 0.5_sp*(y(1) + y(n))
    s = s*h
  end subroutine trapezoid

  !> Integral s over energy points ie = ie1..ie2 of the curves a1, b1 of
  !! beam ib1 and, shifted by nv points, a2, b2 of beam ib2 (step eincr).
  !! The integrand depends on nf:
  !!   nf = 1 (and any nf other than 2, 3, 4): a1
  !!   nf = 2: a1**2
  !!   nf = 3: (a1 - c*a2)**2
  !!   nf = 4: |b1 - c*b2| * |a1 - c*a2| / (|a1| + eps), the Zanazzi-Jona
  !!           integrand, interpolated onto a 10 times denser grid first.
  !! For nf = 1..3 the integrand is left in y(1:ie2-ie1+1); for nf = 4,
  !! y(1:ie2-ie1+1) holds the energies relative to the first point.
  subroutine varsum(a1, a2, b1, b2, nbd1, nbd2, ngp, ib1, ib2, ie1, ie2, nv, &
                    eincr, eps, c, nf, s, y)
    integer, intent(in) :: nbd1, nbd2, ngp, ib1, ib2, ie1, ie2, nv, nf
    real(sp), intent(in) :: a1(nbd1, ngp), b1(nbd1, ngp)
    real(sp), intent(in) :: a2(nbd2, ngp), b2(nbd2, ngp)
    real(sp), intent(in) :: eincr, eps, c
    real(sp), intent(out) :: s
    real(sp), intent(inout) :: y(ngp)

    real(sp) :: y1(ngp), y2(ngp), y3(ngp), y4(ngp), yy(10*ngp)
    real(sp) :: aa1, aa2, ab1, ab2, de, x
    integer :: ie, ies, itil, itih, n, nn

    n = 0
    if (nf /= 4) then
      do ie = ie1, ie2
        n = n + 1
        ies = ie + nv
        select case (nf)
        case (2)
          y(n) = a1(ib1, ie)**2
        case (3)
          y(n) = (a1(ib1, ie) - c*a2(ib2, ies))**2
        case default
          y(n) = a1(ib1, ie)
        end select
      end do
      call trapezoid(y(1:n), n, eincr, s)
      return
    end if

    do ie = ie1, ie2
      n = n + 1
      ies = ie + nv
      y(n) = real(n - 1, sp)*eincr
      y1(n) = a1(ib1, ie)
      y2(n) = a2(ib2, ies)
      y3(n) = b1(ib1, ie)
      y4(n) = b2(ib2, ies)
    end do
    de = 0.1_sp*eincr
    nn = 10*(n - 1) + 1
    do ie = 1, nn
      x = real(ie - 1, sp)*de
      itil = 0
      itih = 0
      aa1 = yval(x, y1(1:n), y(1:n), n, itil, itih)
      itil = 0
      itih = 0
      aa2 = yval(x, y2(1:n), y(1:n), n, itil, itih)
      itil = 0
      itih = 0
      ab1 = yval(x, y3(1:n), y(1:n), n, itil, itih)
      itil = 0
      itih = 0
      ab2 = yval(x, y4(1:n), y(1:n), n, itil, itih)
      yy(ie) = abs(ab1 - c*ab2)*abs(aa1 - c*aa2)/(abs(aa1) + eps)
    end do
    call trapezoid(yy(1:nn), nn, de, s)
  end subroutine varsum

  !> Pendry Y function y = L/(1 + vi**2 L**2) of the curves a(ib, 1:ne(ib))
  !! for beams 1..nb, with L = ap/a the logarithmic derivative (ap: first
  !! derivative per eV, vi: imaginary part of the inner potential). Where
  !! |a| < 1e-7 the equivalent form (a/ap)/((a/ap)**2 + vi**2) is used, and
  !! y = 0 where also |ap| < 1e-7. Entries beyond ne(ib) are not set. e is
  !! used only for printing (ipr >= 2).
  subroutine ypend(a, ap, nbd, ngp, nb, ne, e, y, vi, ipr)
    integer, intent(in) :: nbd, ngp, nb, ipr
    real(sp), intent(in) :: a(nbd, ngp), ap(nbd, ngp), e(nbd, ngp)
    integer, intent(in) :: ne(nbd)
    real(sp), intent(inout) :: y(nbd, ngp)
    real(sp), intent(in) :: vi

    real(sp) :: af, apf
    integer :: ib, ie, n

40  format(' PENDRY Y FUNCTION IN BEAM', 1I3, /, 1000(5(1F7.2, 1E13.4, 3X), /))

    do ib = 1, nb
      n = ne(ib)
      if (n /= 0) then
        do ie = 1, n
          af = a(ib, ie)
          if (abs(af) >= 1.e-7_sp) then
            af = ap(ib, ie)/af
            y(ib, ie) = af/(1._sp + vi*vi*af*af)
          else
            apf = ap(ib, ie)
            if (abs(apf) >= 1.e-7_sp) then
              af = af/apf
              y(ib, ie) = af/(af*af + vi*vi)
            else
              y(ib, ie) = 0._sp
            end if
          end if
        end do
      end if
    end do
    if (ipr >= 2) then
      do ib = 1, nb
        n = ne(ib)
        if (n /= 0) then
          write(6, 40) ib, (e(ib, ie), y(ib, ie), ie = 1, n)
        end if
      end do
    end if
  end subroutine ypend
end module rfactor_terms

!> Averaging of theoretical beams over domains and beam groups, and the
!! check of averaged intensities.
module domain_averaging
  use kinds, only: sp
  implicit none

  private

  public :: maxint, rintav

contains

  !> Largest a(is, ib, ie) over all is = 1..ns, ie = 1..ne, and 0 if none
  !! is positive.
  subroutine maxint(a, ns, nbd, ib, ne, am)
    integer, intent(in) :: ns, nbd, ib, ne
    real(sp), intent(in) :: a(ns, nbd, ne)
    real(sp), intent(out) :: am

    integer :: ie, is

    am = 0._sp
    do is = 1, ns
      do ie = 1, ne
        if (a(is, ib, ie) > am) am = a(is, ib, ie)
      end do
    end do
  end subroutine maxint

  !> Averages the intensities at(:, jb, :) of all beams jb with the same
  !! group label kav(jb) into atav(:, lav, :), lav = 1, 2, ... being the
  !! group labels in increasing order (lavm groups found). Each beam is
  !! weighted by the number of symmetry-equivalent beams its symmetry code
  !! sym(jb) stands for in the LEED program (all 1 without symmetry). The
  !! group takes the name pq1(:, lav) of its first beam in pq. Entries of
  !! atav and pq1 for labels without beams are not set.
  subroutine rintav(at, ns, nb, ne, pq, pq1, kav, sym, lavm, atav)
    integer, intent(in) :: ns, nb, ne
    real(sp), intent(in) :: at(ns, nb, ne), pq(2, nb)
    real(sp), intent(inout) :: pq1(2, nb), atav(ns, nb, ne)
    integer, intent(in) :: kav(nb), sym(nb)
    integer, intent(out) :: lavm

    real(sp) :: wght, wt
    integer :: ib, ie, is, jb, lab, lav

    lavm = 0
    do lav = 1, nb
      do ib = 1, nb
        ! the first beam of group lav
        if (lav /= kav(ib)) cycle
        lavm = lavm + 1
        pq1(:, lav) = pq(:, ib)
        wt = 0.0_sp
        do jb = ib, nb
          if (lav /= kav(jb)) cycle
          lab = sym(jb)
          if (lab < 2) then
            wght = 1.0_sp
          else if (lab < 7) then
            wght = 2.0_sp
          else if (lab < 10) then
            wght = 4.0_sp
          else if (lab == 10) then
            wght = 8.0_sp
          else if (lab < 12) then
            wght = 3.0_sp
          else if (lab < 15) then
            wght = 6.0_sp
          else
            wght = 12.0_sp
          end if
          if (wt < 0.001_sp) then
            do is = 1, ns
              do ie = 1, ne
                atav(is, lav, ie) = wght*at(is, jb, ie)
              end do
            end do
          else
            do is = 1, ns
              do ie = 1, ne
                atav(is, lav, ie) = atav(is, lav, ie) + wght*at(is, jb, ie)
              end do
            end do
          end if
          wt = wt + wght
        end do
        do is = 1, ns
          do ie = 1, ne
            atav(is, lav, ie) = atav(is, lav, ie)/wt
          end do
        end do
        exit
      end do
    end do
  end subroutine rintav
end module domain_averaging

!> Reading of the experimental (or pseudo-experimental) beams and averaging
!! of experimental beams with the same group label.
module experiment_input
  use kinds, only: sp
  implicit none

  private

  public :: expav, gapin, reade, readt

contains

  !> Reads nbed experimental beams from unit 8: a format line, then per beam
  !! its name, the number of points nee(ib) with a factor fac, and the pairs
  !! (ee(ib, ie), ae(ib, ie)). The intensities are divided by fac, meant to
  !! bring them to the order of 1 and to match curves averaged together.
  subroutine reade(ae, ee, nbed, ngp, nee, bename, ipr)
    integer, intent(in) :: nbed, ngp, ipr
    real(sp), intent(inout) :: ae(nbed, ngp), ee(nbed, ngp)
    integer, intent(out) :: nee(nbed)
    character(len=4), intent(out) :: bename(5, nbed)

    character(len=4) :: fmt(19)
    real(sp) :: fac
    integer :: i, ib, ie, n

10  format(19A4)
35  format(1I4, 1E14.4)
70  format('*EXP. BEAM ', 1I4, ' (', 5A4, ')')
80  format('* EXP. ENERGIES AND INTENSITIES', /, 100(5(1F7.2, 1E13.4, 3X), /))

    read(8, 10) (fmt(i), i = 1, 19)
    do ib = 1, nbed
      read(8, 10) (bename(i, ib), i = 1, 5)
      write(6, 70) ib, (bename(i, ib), i = 1, 5)
      read(8, 35) nee(ib), fac
      n = nee(ib)
      read(8, fmt) (ee(ib, ie), ae(ib, ie), ie = 1, n)
      if (ipr >= 1) write(6, 80) (ee(ib, ie), ae(ib, ie), ie = 1, n)
      do ie = 1, n
        ae(ib, ie) = ae(ib, ie)/fac
      end do
    end do
  end subroutine reade

  !> Reads theoretical curves from unit 8 as pseudo-experiment: the number
  !! of beams (must equal nbed), the beam names, then per energy the
  !! intensities of all beams. Leading points with intensity <= 1e-20 are
  !! dropped per beam; nee(ib) counts the points above 1e-20.
  subroutine readt(ae, ee, nbed, ngp, nee, bename)
    integer, intent(in) :: nbed, ngp
    real(sp), intent(inout) :: ae(nbed, ngp), ee(nbed, ngp)
    integer, intent(out) :: nee(nbed)
    character(len=4), intent(out) :: bename(5, nbed)

    real(sp) :: eeh(ngp), dummy
    integer :: ioff(nbed)
    integer :: i, ib, inn, j, k, nbeams

101 format(I5, 22A4)
103 format(I5, E13.4)
105 format(F7.2, F7.4, 4E14.5, /, 2000(5E14.5, /))

    read(8, 103) nbeams
    write(6, 103) nbeams
    if (nbed /= nbeams) then
      write(7, '(A33)') '*********************************'
      write(7, '(A33)') 'SOMETHING IS WRONG WITH NUMBER OF'
      write(7, '(A33)') 'PSEUDOEXPERIMENTAL BEAMS!        '
      write(7, '(A33)') 'NO AVERAGING ALLOWED!            '
      write(7, '(A33)') '*********************************'
      stop 'ERROR 1 IN PSEUDOEXPERIMENT'
    end if
    do ib = 1, nbed
      read(8, 101) inn, (bename(i, ib), i = 1, 5)
      write(6, 101) inn, (bename(i, ib), i = 1, 5)
    end do

    do k = 1, nbed
      nee(k) = 0
      ioff(k) = 0
      do j = 1, ngp
        ae(k, j) = 0._sp
      end do
    end do

    do j = 1, ngp
      read(8, 105, end=250) eeh(j), dummy, (ae(k, j), k = 1, nbed)
      do k = 1, nbed
        if (ae(k, j) > 1.e-20_sp) then
          nee(k) = nee(k) + 1
        else if (nee(k) == 0) then
          ioff(k) = ioff(k) + 1
        end if
      end do
    end do

250 continue
    do k = 1, nbed
      do j = 1, nee(k)
        ee(k, j) = eeh(j + ioff(k))
        ae(k, j) = ae(k, j + ioff(k))
      end do
    end do
  end subroutine readt

  !> Reads the gaps of the experimental spectra from unit 8: their number
  !! ngap, then per gap the beam bm(ig) and the energy range
  !! [eg1(ig), eg2(ig)], which is left out of smoothing, interpolation and
  !! the R factor. Gaps must be ordered by beam, and within a beam by
  !! energy; otherwise the program stops. gcount(ib) is incremented for
  !! every gap of beam ib.
  subroutine gapin(ngap, bm, eg1, eg2, gcount, nbed)
    integer, intent(out) :: ngap
    integer, intent(inout) :: bm(:)
    real(sp), intent(inout) :: eg1(:), eg2(:)
    integer, intent(in) :: nbed
    integer, intent(inout) :: gcount(nbed)

    integer :: ig

10  format(I3, 2F7.2)
20  format('*THERE ARE', I3, ' GAPS IN THE EXPERIMENTAL SPECTRA')
30  format(' *** GAPS MUST BE IN ORDER OF EXPERIMENTAL BEAMS ***', /, &
           ' ***               PROGRAM STOPPED               ***')
40  format(' *** GAPS IN BEAM', I3, ' MUST BE ORDERED BY ENERGY ***', /, &
           ' ***              PROGRAM STOPPED              ***')
50  format(' GAP', I3, ' IN BEAM', I3, ' BETWEEN', F7.2, ' AND', F7.2, ' EV')

    read(8, 10) ngap
    write(6, 20) ngap

    do ig = 1, ngap
      read(8, 10) bm(ig), eg1(ig), eg2(ig)
      if (ig > 1) then
        if (bm(ig) < bm(ig - 1)) then
          write(6, 30)
          stop
        else if (bm(ig) == bm(ig - 1) .and. eg1(ig) < eg1(ig - 1)) then
          write(6, 40) bm(ig)
          stop
        end if
      end if
      write(6, 50) ig, bm(ig), eg1(ig), eg2(ig)
      gcount(bm(ig)) = gcount(bm(ig)) + 1
    end do
  end subroutine gapin

  !> Averages the experimental beams with the same group label nbea(ib) onto
  !! a common energy grid in increasing order, leaving out points inside
  !! the gaps (eg1, eg2) of each beam. The nbe averaged beams replace the
  !! first rows of ae, ee, nee and bename (each named after its first beam);
  !! their gaps are returned in eg1, eg2 (energies), ng1, ng2 (indices),
  !! ngap in total and gcount per beam. nbea is set to 0. a and nv are work
  !! space.
  subroutine expav(ae, ee, nbed, ngp, nee, bename, nbea, nbe, ipr, a, nv, &
                   gcount, ngap, eg1, eg2, ng1, ng2)
    integer, intent(in) :: nbed, ngp, ipr
    real(sp), intent(inout) :: ae(nbed, ngp), ee(nbed, ngp)
    integer, intent(inout) :: nee(nbed), nbea(nbed), gcount(nbed)
    character(len=4), intent(inout) :: bename(5, nbed)
    integer, intent(out) :: nbe, ngap
    real(sp), intent(inout) :: a(ngp), eg1(:), eg2(:)
    integer, intent(inout) :: nv(ngp), ng1(:), ng2(:)

    integer :: gc(nbed)
    real(sp) :: de, e, ediff, emax, emin
    integer :: gin, i, ib, ib1, ibmin, ie, iemax, ien, ig, igap, igc, n, &
               nbeat, nemax

10  format('*EXP. ENERG. AND INTENS. AFTER AVERAGING IN BEAM', 1I4, /, &
           100(5(1F7.2, 1E13.4, 3X), /))

    nbe = 1
    igap = 0
    ngap = 0

    do ib = 1, nbed
      if (nbea(ib) /= 0) then
        emin = 1.e6_sp
        emax = 0._sp
        de = 1.e6_sp

        ! energy range of all beams averaged with ib (including ib)
        do ib1 = ib, nbed
          if (nbea(ib1) == nbea(ib)) then
            do ie = 1, nee(ib1)
              if (emin > ee(ib1, ie)) then
                emin = ee(ib1, ie)
                ibmin = ib1
              end if
              if (emax < ee(ib1, ie)) then
                emax = ee(ib1, ie)
              end if
            end do
          end if
        end do

        ! energy step: smallest step above 0.01 from emin in beam ibmin
        do ie = 1, nee(ibmin)
          ediff = ee(ibmin, ie) - emin
          if (ediff < de .and. ediff > 0.01_sp) then
            de = ediff
          end if
        end do
        nemax = int((emax - emin)/de + 0.0001_sp) + 1
        do ie = 1, nemax
          nv(ie) = 0
          a(ie) = 0._sp
        end do

        ! sum over the same beams again, leaving out points inside gaps
        iemax = 0
        nbeat = nbea(ib)
        do ib1 = ib, nbed
          if (nbea(ib1) == nbeat) then
            do ie = 1, nee(ib1)
              e = ee(ib1, ie)
              ien = int((e - emin)/de + 0.0001_sp) + 1
              if (ien > iemax) then
                iemax = ien
              end if
              if (gcount(ib1) == 0) then
                nv(ien) = nv(ien) + 1
                a(ien) = a(ien) + ae(ib1, ie)
              else
                gin = 0
                igc = igap
                do i = 1, gcount(ib1)
                  igc = igc + 1
                  if (e >= eg1(igc) .and. e <= eg2(igc)) then
                    gin = 1
                  end if
                end do
                if (gin /= 1) then
                  nv(ien) = nv(ien) + 1
                  a(ien) = a(ien) + ae(ib1, ie)
                end if
              end if
            end do
            nbea(ib1) = 0
          end if
        end do

        ! averages; grid points without data form the gaps of beam nbe
        gin = 0
        ig = 0
        do ie = 1, iemax
          ee(nbe, ie) = emin + real(ie - 1, sp)*de
          if (nv(ie) /= 0) then
            ae(nbe, ie) = a(ie)/real(nv(ie), sp)
            gin = 0
          else if (gin == 0) then
            ig = ig + 1
            eg1(ngap + ig) = ee(nbe, ie)
            eg2(ngap + ig) = ee(nbe, ie)
            ng1(ngap + ig) = ie
            ng2(ngap + ig) = ie
            gin = 1
          else
            eg2(ngap + ig) = ee(nbe, ie)
            ng2(ngap + ig) = ie
          end if
        end do

        nee(nbe) = iemax
        gc(nbe) = ig
        do i = 1, 5
          bename(i, nbe) = bename(i, ib)
        end do
        ngap = ngap + ig
        nbe = nbe + 1
      end if
      igap = igap + gcount(ib)
    end do

    nbe = nbe - 1
    do ib = 1, nbe
      gcount(ib) = gc(ib)
    end do

    if (ipr == 2) then
      do ib = 1, nbe
        n = nee(ib)
        write(6, 10) ib, (ee(ib, ie), ae(ib, ie), ie = 1, n)
      end do
    end if
  end subroutine expav
end module experiment_input

!> Reading of the theoretical beams written by the LEED program, for all
!! geometries, and extraction of one geometry.
module theory_input
  use kinds, only: sp
  implicit none

  private

  public :: rdpl, strip2

contains

  !> Reads the theoretical beams from unit 5: per beam its label, name
  !! pq(:, j) and symmetry code sym(j); then per energy e(i) and geometry
  !! its value spac and the nb intensities ydata(:, :, i). Energies above
  !! emax end the reading, those below emin are overwritten by the next.
  !! At an end of file, up to limfil files are read, skipping the nb + 2
  !! header records of each further file. Geometries listed in nssk are
  !! skipped; nss is the number left. ne returns the number of energies.
  subroutine rdpl(e, ngp, ydata, nss, nb, mne, ne, spac, pq, kav, sym, &
                  emin, emax, limfil, nssk, ns, mnbtd)
    integer, intent(in) :: ngp, nss, nb, mne, limfil, ns, mnbtd
    real(sp), intent(inout) :: e(ngp), ydata(nss, nb, mne), spac(nss)
    real(sp), intent(out) :: pq(2, nb)
    integer, intent(in) :: kav(nb), nssk(ns)
    integer, intent(out) :: sym(nb), ne
    real(sp), intent(in) :: emin, emax

    integer :: lab(mnbtd)
    real(sp) :: spacj(nss), yd(mnbtd)
    character(len=4) :: idum
    integer :: i, ifile, j, js, kk

20  format(I5, 2F10.5, I3)
25  format('DIMENSION MNET TOO SMALL, PLEASE CORRECT!', /, &
           '      ****   PROGRAM STOPPED   ****')
27  format('*THEOR. BEAM GROUPINGS')
29  format(' ', 2F10.5, I5)
40  format(F7.2, F7.4, 4E14.5, /, 2000(5E14.5, /))
57  format(A4)
65  format('*SURFACE STRUCTURE NO. AND CHARACTERIZATION')
75  format(18X, I3, F15.4)

    do j = 1, nb
      read(5, 20) lab(j), (pq(i, j), i = 1, 2), sym(j)
    end do
    write(6, 27)
    do j = 1, nb
      write(6, 29) (pq(i, j), i = 1, 2), kav(j)
    end do

    ifile = 1
    i = 1
    energies: do
      js = 0
      if (i > mne + 1) then
        write(6, 25)
        stop
      end if
      do j = 1, ns
        read(5, 40, end=55) e(i), spacj(j), yd(1:nb)
        if (any(nssk(1:ns) == j)) cycle
        js = js + 1
        spac(js) = spacj(j)
        ydata(js, 1:nb, i) = yd(1:nb)
      end do
      if (e(i) > emax) exit energies
      if (e(i) >= emin) i = i + 1
      cycle energies
55    ifile = ifile + 1
      if (ifile > limfil) exit energies
      do kk = 1, nb + 2
        read(5, 57) idum
      end do
    end do energies
    ne = i - 1

    write(6, 65)
    do j = 1, nss
      write(6, 75) j, spac(j)
    end do
  end subroutine rdpl

  !> Copies the intensities ys(is, :, 1:ne) of geometry is to y and the
  !! energies es to every beam of et.
  subroutine strip2(ys, ns, nbd, ne, is, y, es, et, ngp)
    integer, intent(in) :: ns, nbd, ne, is, ngp
    real(sp), intent(in) :: ys(ns, nbd, ne), es(ngp)
    real(sp), intent(inout) :: y(nbd, ngp), et(nbd, ngp)

    integer :: ib, ie

    do ib = 1, nbd
      do ie = 1, ne
        et(ib, ie) = es(ie)
        y(ib, ie) = ys(is, ib, ie)
      end do
    end do
  end subroutine strip2
end module theory_input

!> Energy interval common to an experimental and a theoretical beam.
module common_interval
  use kinds, only: sp
  implicit none

  private

  public :: comnei

contains

  !> Finds the energy interval common to experimental beam ibe (grid ee)
  !! and theoretical beam ibt (grid et, shifted by the inner potential v0),
  !! both on steps of eincr: it is ee(ibe, ne1:ne2) and et(ibt, nt1:nt2).
  subroutine comnei(ee, nbed, nee, et, nbtd, ngp, net, ibe, ibt, v0, eincr, &
                    ne1, ne2, nt1, nt2)
    integer, intent(in) :: nbed, nbtd, ngp, ibe, ibt
    real(sp), intent(in) :: ee(nbed, ngp), et(nbtd, ngp), v0, eincr
    integer, intent(in) :: nee(nbed), net(nbtd)
    integer, intent(out) :: ne1, ne2, nt1, nt2

    real(sp) :: de1, de2
    integer :: ne, nt

    ne = nee(ibe)
    nt = net(ibt)
    de1 = et(ibt, 1) + v0 - ee(ibe, 1)
    de2 = et(ibt, nt) + v0 - ee(ibe, ne)
    if (de1 < 0._sp) then
      ne1 = 1
      nt1 = int((-de1/eincr) + 0.0001_sp) + 1
    else
      ne1 = int((de1/eincr) + 0.0001_sp) + 1
      nt1 = 1
    end if
    if (de2 < 0._sp) then
      ne2 = ne - int((-de2/eincr) + 0.0001_sp)
      nt2 = nt
    else
      ne2 = ne
      nt2 = nt - int((de2/eincr) + 0.0001_sp)
    end if
  end subroutine comnei
end module common_interval

!> Plot output: the experimental and best-fit theoretical spectra in column
!! format, and gnuplot and LaTeX scripts that plot them.
module plot_output
  use kinds, only: sp
  implicit none

  private

  public :: column, script

contains

  !> Writes the experimental spectra to exp.column and the best-fit
  !! theoretical spectra to theo.column, one line per energy step with an
  !! (energy, intensity) pair per beam. Energies outside a beam's range or
  !! inside one of its gcount(ibe) gaps (indices nge1:nge2 resp. nget1:nget2,
  !! numbered over all beams) are written as zero intensity. Intensities are
  !! scaled to intmax: with norm = 1 after division by the integrals tse
  !! resp. tst, to the largest value of both curves; with norm = 2 each
  !! curve to its own maximum. ngap is the size of the gap arrays.
  subroutine column(neminb, nemaxb, ntminb, ntmaxb, ee, ae, bet, bat, nbed, &
                    nbtd, ngp, nbe, ibp, gc, nge1, nge2, ngt1, ngt2, ngap, &
                    tse, tst, v0, intmax, norm)
    integer, intent(in) :: nbed, nbtd, ngp, nbe, ngap, norm
    integer, intent(in) :: neminb(nbed), nemaxb(nbed)
    integer, intent(in) :: ntminb(nbtd), ntmaxb(nbtd)
    real(sp), intent(in) :: ee(nbed, ngp), ae(nbed, ngp)
    real(sp), intent(in) :: bet(nbtd, ngp), bat(nbtd, ngp)
    integer, intent(in) :: ibp(nbed), gc(nbed)
    integer, intent(in) :: nge1(:), nge2(:), ngt1(:), ngt2(:)
    real(sp), intent(in) :: tse(nbed), tst(nbtd), v0, intmax

    real(sp) :: aei(nbed), eei(nbed), bati(nbed), beti(nbed)
    real(sp) :: batint(nbed, ngp), betint(nbed, ngp), tstint(nbed)
    real(sp) :: amax(nbed), amaxt(nbed)
    integer :: ntmini(nbed), ntmaxi(nbed)
    integer :: igc(nbed), nge(nbed, ngap, 2), ngt(nbed, ngap, 2)
    integer :: neemin, neemax, netmin, netmax
    integer :: i, ibe, iee, iet, igap, igdum, igp
    character(len=10) :: fmt

10  format('(', I4, 'F8.2)')
11  format(A10)

    open(10, file='exp.column')
    open(11, file='theo.column')

    ! spectra output format for exp.column and theo.column
    write(10, 10) 2*nbe
    rewind(10)
    read(10, 11) fmt
    rewind(10)

    ! lowest lower and highest upper energy indices of all beams
    neemin = ngp
    neemax = 1
    netmin = ngp
    netmax = 1

    do ibe = 1, nbe
      if (ibp(ibe) == 0) then
        ntmini(ibe) = ngp
        ntmaxi(ibe) = 1
        tstint(ibe) = 1._sp
        do igp = 1, ngp
          betint(ibe, igp) = 0._sp
          batint(ibe, igp) = 0._sp
        end do
      else
        ntmini(ibe) = ntminb(ibp(ibe))
        ntmaxi(ibe) = ntmaxb(ibp(ibe))
        tstint(ibe) = tst(ibp(ibe))
        do igp = 1, ngp
          betint(ibe, igp) = bet(ibp(ibe), igp)
          batint(ibe, igp) = bat(ibp(ibe), igp)
        end do
      end if
    end do

    do ibe = 1, nbe
      amax(ibe) = 0._sp
      if (neminb(ibe) < neemin) then
        neemin = neminb(ibe)
      end if
      if (nemaxb(ibe) > neemax) then
        neemax = nemaxb(ibe)
      end if
      amaxt(ibe) = 0._sp
      if (ntmini(ibe) < netmin) then
        netmin = ntmini(ibe)
      end if
      if (ntmaxi(ibe) > netmax) then
        netmax = ntmaxi(ibe)
      end if
    end do

    ! gap indices per beam
    do ibe = 1, nbe
      do i = 1, ngap
        nge(ibe, i, 1) = 0
        nge(ibe, i, 2) = 0
        ngt(ibe, i, 1) = 0
        ngt(ibe, i, 2) = 0
      end do
    end do

    igap = 0
    do ibe = 1, nbe
      do i = 1, gc(ibe)
        nge(ibe, i, 1) = nge1(i + igap)
        nge(ibe, i, 2) = nge2(i + igap)
        ngt(ibe, i, 1) = ngt1(i + igap)
        ngt(ibe, i, 2) = ngt2(i + igap)
      end do
      igap = igap + gc(ibe)
    end do

    ! highest intensity amax(ibe) of each experimental beam
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iee = neemin, neemax
      do ibe = 1, nbe
        if (iee >= neminb(ibe) .and. iee <= nemaxb(ibe)) then
          igdum = igc(ibe) + 1
          if (gc(ibe) /= 0 .and. iee >= nge(ibe, igdum, 1) .and. &
              iee <= nge(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
            aei(ibe) = 0._sp
            if (iee == nge(ibe, igdum, 2)) then
              igc(ibe) = igdum
            end if
          else if (norm == 1) then
            aei(ibe) = ae(ibe, iee)/tse(ibe)
          else if (norm == 2) then
            aei(ibe) = ae(ibe, iee)
          end if

          if (aei(ibe) > amax(ibe)) then
            amax(ibe) = aei(ibe)
          end if
        end if
      end do
    end do

    ! highest intensity of each theoretical beam: into amax (norm = 1) or
    ! amaxt (norm = 2)
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iet = netmin, netmax
      do ibe = 1, nbe
        if (iet >= ntmini(ibe) .and. iet <= ntmaxi(ibe)) then
          igdum = igc(ibe) + 1
          if (gc(ibe) /= 0 .and. iet >= ngt(ibe, igdum, 1) .and. &
              iet <= ngt(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
            bati(ibe) = 0._sp
            if (iet == ngt(ibe, igdum, 2)) then
              igc(ibe) = igdum
            end if
          else if (norm == 1) then
            bati(ibe) = batint(ibe, iet)/tstint(ibe)
            if (bati(ibe) > amax(ibe)) then
              amax(ibe) = bati(ibe)
            end if
          else if (norm == 2) then
            bati(ibe) = batint(ibe, iet)
            if (bati(ibe) > amaxt(ibe)) then
              amaxt(ibe) = bati(ibe)
            end if
          end if
        end if
      end do
    end do

    ! normalised experimental beams to exp.column
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iee = neemin, neemax
      do ibe = 1, nbe
        igdum = igc(ibe) + 1
        if (iee < neminb(ibe) .or. iee > nemaxb(ibe)) then
          eei(ibe) = 0._sp
          aei(ibe) = 0._sp
        else if (gc(ibe) /= 0 .and. iee >= nge(ibe, igdum, 1) .and. &
                 iee <= nge(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
          eei(ibe) = 0._sp
          aei(ibe) = 0._sp
          if (iee == nge(ibe, igdum, 2)) then
            igc(ibe) = igdum
          end if
        else
          eei(ibe) = ee(ibe, iee)
          if (norm == 1) then
            aei(ibe) = ae(ibe, iee)*intmax/(tse(ibe)*amax(ibe))
          else if (norm == 2) then
            aei(ibe) = ae(ibe, iee)*intmax/amax(ibe)
          end if
        end if
      end do

      write(10, fmt) (eei(ibe), aei(ibe), ibe = 1, nbe)
    end do

    ! normalised theoretical beams to theo.column
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iet = netmin, netmax
      do ibe = 1, nbe
        igdum = igc(ibe) + 1
        if (iet < ntmini(ibe) .or. iet > ntmaxi(ibe)) then
          beti(ibe) = -v0
          bati(ibe) = 0._sp
        else if (gc(ibe) /= 0 .and. iet >= ngt(ibe, igdum, 1) .and. &
                 iet <= ngt(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
          beti(ibe) = -v0
          bati(ibe) = 0._sp
          if (iet == ngt(ibe, igdum, 2)) then
            igc(ibe) = igdum
          end if
        else
          beti(ibe) = betint(ibe, iet)
          if (norm == 1) then
            bati(ibe) = batint(ibe, iet)*intmax/(tstint(ibe)*amax(ibe))
          else if (norm == 2) then
            bati(ibe) = batint(ibe, iet)*intmax/amaxt(ibe)
          end if
        end if
      end do

      write(11, fmt) (beti(ibe) + v0, bati(ibe), ibe = 1, nbe)
    end do
  end subroutine column

  !> Writes plot.gnu, a gnuplot script that plots the spectra of
  !! exp.column and theo.column to one PDF per beam (energy range emin to
  !! emax, labelled with the beam name and its R factor rb), and plot.tex,
  !! a LaTeX document with one page per plot under the title and the
  !! average R factor brav.
  subroutine script(emin, emax, intmax, plsize, xtics, nbe, bename, rb, &
                    nbed, title, brav, ibp)
    integer, intent(in) :: nbe, nbed, xtics
    real(sp), intent(in) :: emin, emax, intmax, plsize(2), rb(nbed), brav
    character(len=4), intent(in) :: bename(5, nbed)
    character(len=80), intent(in) :: title
    integer, intent(in) :: ibp(nbed)

    integer :: i, ibe

10  format('set terminal pdfcairo enhanced font "Arial,12"', /, &
           'set encoding iso_8859_1')
11  format('set xrange [', I3, ':', I3, ']', /, 'set yrange [0:', I4, ']')
12  format('set size ', F5.3, ',', F5.3, /, 'set style data lines', /, &
           'set xtics ', I3, /, 'set noytics', /, &
           'set xlabel "Energy [eV]"', /, 'set noxzeroaxis', /)
15  format('set output "pic_', I3.3, '.pdf"', /, 'set nolabel')
16  format('set label "', 5A4, '\nR =', F7.4, '" at ', I3, &
           ',', I4, ' center font "Helvetica,12"')
17  format('plot ''exp.column'' using ($', I3, '>0?$', I3, ':1/0) :', &
           '($', I3, '>0?$', I3, ':1/0) title "(exp.)" ', &
           'lt rgb "red" lw 2,\', / &
           ' ''theo.column'' using ($', I3, '>0?$', I3, ':1/0) :', &
           '($', I3, '>0?$', I3, ':1/0) title "(theory)" ', &
           'lt rgb "blue" lw 2', /)
20  format('\documentclass[12pt,a4paper]{article}', /, &
           '\usepackage{graphicx}', /, '\usepackage{float}', /)
21  format('\begin{document}', /, '\pagestyle{empty}', /)
22  format('\newcommand{\bigtitle}{', /, '\begin{center}', /, '\large', &
           /, '\textbf{', A80, '}\\', /, '$\mathbf{R_{P} =', F7.4, '}$', &
           /, '\end{center}', /, '}', /)
25  format('\newpage', /, '\bigtitle', /, '\begin{figure}[H]')
26  format('\includegraphics[angle=0,width=14cm]{pic_', I3.3, &
           '.pdf}', /, '\end{figure}', /)
27  format('\end{document}')

    open(12, file='plot.gnu')
    open(13, file='plot.tex')

    ! general settings of the plots
    write(12, 10)
    write(12, 11) int(emin), int(emax), int(1.5_sp*intmax)
    write(12, 12) plsize(1), plsize(2), xtics

    ! one plot per beam
    do ibe = 1, nbe
      if (ibp(ibe) /= 0) then
        write(12, 15) ibe
        write(12, 16) (bename(i, ibe), i = 1, 5), rb(ibe), &
          int(emin + (emax - emin)/2), int(1.2_sp*intmax)
        write(12, 17) (2*ibe) - 1, (2*ibe) - 1, (2*ibe) - 1, 2*ibe, &
          (2*ibe) - 1, (2*ibe) - 1, (2*ibe) - 1, 2*ibe
      end if
    end do

    ! document settings and title
    write(13, 20)
    write(13, 21)
    write(13, 22) title, brav

    ! one page per plot
    do ibe = 1, nbe
      if (ibp(ibe) /= 0) then
        write(13, 25)
        write(13, 26) ibe
      end if
    end do

    write(13, 27)
  end subroutine script
end module plot_output
