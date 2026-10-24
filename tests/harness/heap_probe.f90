! Print one never-written byte of a small and of a large (8 MiB) allocation.
! With harness.runner.POISONED_HEAP_ENV both are 127 (0x7f); without it the
! large block comes from fresh, zeroed pages (mmap) and reads 0.
program heap_probe
  use, intrinsic :: iso_fortran_env, only: int8
  implicit none
  integer(int8), allocatable :: small(:), large(:)
  allocate(small(64), large(8 * 1024 * 1024))
  print '(i0, 1x, i0)', small(33), large(4 * 1024 * 1024)
end program heap_probe
