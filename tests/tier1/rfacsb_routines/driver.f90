!> Calls one routine of module derivatives (lib/rfacsb.f90) on the
!! inputs read from stdin and prints its results.
!!
!! Every real is read and written as the 8 hex digits of its bit pattern,
!! one value per line, so that inputs and results are exact. stdin:
!!   DER:     name, ne, h, y(1:ne)                 -> y1(1:ne)
!!   SMOOTH:  name, ismoth, n, ee(1:n), ae(1:n)    -> ae(1:n)
!!   YSMOOTH: name, n, vi, eincr, a(1:n), ap(1:n)  -> y(1:n) (one beam)
!! generate/rfacsb_routines_legacy_driver.f is the same driver for the
!! legacy routines.
program rfacsb_routines_driver
  use, intrinsic :: iso_fortran_env, only: int32
  use kinds, only: sp
  use derivatives, only: der, smooth, ysmooth
  implicit none

  character(len=16) :: name
  integer :: n, ismoth
  real(sp) :: h, vi, eincr
  real(sp), allocatable :: x(:), y(:), ap(:)
  real(sp), allocatable :: a2(:, :), ap2(:, :), e2(:, :), y2(:, :)

  read(5, '(A)') name
  select case (trim(name))
  case ('DER')
    read(5, *) n
    h = read_real()
    allocate(y(n), x(n))
    call read_reals(y)
    call der(y, n, x, h)
    call write_reals(x)
  case ('SMOOTH')
    read(5, *) ismoth
    read(5, *) n
    allocate(x(n), y(n))
    call read_reals(x)
    call read_reals(y)
    call smooth(ismoth, y, x, n, 0, 1)
    call write_reals(y)
  case ('YSMOOTH')
    read(5, *) n
    vi = read_real()
    eincr = read_real()
    allocate(x(n), ap(n))
    call read_reals(x)
    call read_reals(ap)
    allocate(a2(1, n), ap2(1, n), e2(1, n), y2(1, n))
    a2(1, :) = x
    ap2(1, :) = ap
    e2 = 0.0_sp
    y2 = 0.0_sp
    call ysmooth(a2, ap2, 1, n, 1, [n], e2, y2, vi, eincr, 0)
    call write_reals(y2(1, :))
  case default
    error stop 'unknown routine'
  end select

contains

  real(sp) function read_real()
    integer(int32) :: bits
    read(5, '(Z8)') bits
    read_real = transfer(bits, read_real)
  end function read_real

  subroutine read_reals(v)
    real(sp), intent(out) :: v(:)
    integer :: i
    do i = 1, size(v)
      v(i) = read_real()
    end do
  end subroutine read_reals

  subroutine write_reals(v)
    real(sp), intent(in) :: v(:)
    integer :: i
    do i = 1, size(v)
      write(6, '(Z8.8)') transfer(v(i), 0_int32)
    end do
  end subroutine write_reals
end program rfacsb_routines_driver
