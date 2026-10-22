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
