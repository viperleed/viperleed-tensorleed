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
