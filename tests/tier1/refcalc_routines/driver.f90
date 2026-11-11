!> Calls the routines that set up the angular coefficients (CAAA, CELMG,
!! LXGENT, CPPP) and the sublattices (SLIND) of ref-calc on the inputs read
!! from stdin and prints their results.
!!
!! Every group of routines is called from one wrapper below (caa_table,
!! clm_table, lm_permutations, gaunt_table, sublattice_table), so that a
!! group can be switched to its module version without touching the others.
!! Output, one block per array: a header line "name kind count", then the
!! values: kind r4 / r8 as the hex digits of the bit patterns (8 / 16 per
!! value), kind i4 in decimal. stdin:
!!   ANGULAR: name, lmax
!!   SLIND:   name, ideg nl1 nl2, then br1(2) br2(2) ar1(2) ar2(2)
!!            rar1(2) rar2(2) as hex bit patterns, one per line
!! generate/refcalc_routines_legacy_driver.f is the same driver for the
!! legacy routines of the baseline.
program refcalc_routines_driver
  use, intrinsic :: iso_fortran_env, only: int32, int64
  use angular_coefficients, only: composite_layer_coefficients, matrix_coefficients
  use kinds, only: sp, dp
  implicit none

  ! dimensions of the legacy CAA and CLM as a function of LMAX
  ! (the DATA tables NCA and NLMS of ref-calc.f)
  integer, parameter :: nca(18) = [1, 70, 264, 759, 1820, 3836, 7344, 13053, &
      21868, 34914, 53560, 79443, 114492, 160952, 221408, 298809, 396492, &
      518206]
  integer, parameter :: nlms(18) = [0, 76, 284, 809, 1925, 4032, 7680, 13593, &
      22693, 36124, 55276, 81809, 117677, 165152, 226848, 305745, 405213, &
      529036]
  ! size of FAC of ref-calc.f for MLMAX = 18
  integer, parameter :: nfac = 4 * 18 + 3

  character(len=16) :: name

  read(5, '(A)') name
  select case (trim(name))
  case ('ANGULAR')
    call angular_coefficients_case()
  case ('SLIND')
    call sublattice_case()
  case default
    error stop 'unknown routine'
  end select

contains

  subroutine angular_coefficients_case()
    integer :: lmax, lmmax
    real(sp), allocatable :: caa(:), clm(:), ppp(:, :, :)
    real(dp) :: fac(nfac)
    integer, allocatable :: lx(:), lxi(:), lt(:), lxm(:)

    read(5, *) lmax
    lmmax = (lmax + 1)**2

    allocate(caa(nca(lmax)))
    call caa_table(lmax, caa)
    call write_real('CAA', caa)

    allocate(clm(nlms(lmax)))
    fac = 0.0_dp
    call clm_table(lmax, clm, fac)
    call write_real('CLM', clm)
    call write_double('FAC_CELMG', fac)

    allocate(lx(lmmax), lxi(lmmax), lt(lmmax), lxm(lmmax))
    call lm_permutations(lmax, lx, lxi, lt, lxm)
    call write_int('LX', lx)
    call write_int('LXI', lxi)
    call write_int('LT', lt)
    call write_int('LXM', lxm)

    allocate(ppp(2 * lmax + 1, lmax + 1, lmax + 1))
    call gaunt_table(lmax, ppp, fac)
    call write_real('PPP', reshape(ppp, [size(ppp)]))
    call write_double('FAC_CPPP', fac)
  end subroutine angular_coefficients_case

  !> CAAA: Clebsch-Gordan coefficients of the composite layer
  subroutine caa_table(lmax, caa)
    integer, intent(in) :: lmax
    real(sp), intent(out) :: caa(:)

    caa = 0.0_sp
    call composite_layer_coefficients(lmax, caa)
  end subroutine caa_table

  !> CELMG: Clebsch-Gordan coefficients of the matrices X and TAU
  subroutine clm_table(lmax, clm, fac)
    integer, intent(in) :: lmax
    real(sp), intent(out) :: clm(:)
    real(dp), intent(inout) :: fac(:)

    clm = 0.0_sp
    call matrix_coefficients(lmax, clm, fac)
  end subroutine clm_table

  !> LXGENT: permutations of the (L,M) sequence
  subroutine lm_permutations(lmax, lx, lxi, lt, lxm)
    integer, intent(in) :: lmax
    integer, intent(out) :: lx(:), lxi(:), lt(:), lxm(:)

    lx = 0
    lxi = 0
    lt = 0
    lxm = 0
    call lxgent(lx, lxi, lt, lxm, lmax, (lmax + 1)**2)
  end subroutine lm_permutations

  !> CPPP: Gaunt coefficients for the temperature-dependent phase shifts
  subroutine gaunt_table(lmax, ppp, fac)
    integer, intent(in) :: lmax
    real(sp), intent(out) :: ppp(:, :, :)
    real(dp), intent(inout) :: fac(:)

    ppp = 0.0_sp
    call cppp(ppp, size(ppp, 1), size(ppp, 2), size(ppp, 3), size(fac), fac)
  end subroutine gaunt_table

  subroutine sublattice_case()
    integer :: ideg, nl1, nl2, nl
    real(sp) :: t(12)
    real(sp), allocatable :: v(:, :)
    complex(sp), allocatable :: vl(:, :)
    integer, allocatable :: jjs(:, :)

    read(5, *) ideg, nl1, nl2
    call read_reals(t)
    nl = nl1 * nl2
    allocate(v(nl, 2), vl(nl, 2), jjs(nl, ideg))
    call sublattice_table(ideg, nl1, nl2, t, v, vl, jjs)
    call write_real('V', reshape(v, [size(v)]))
    call write_real('VL', transfer(vl, [0.0_sp]))
    call write_int('JJS', reshape(jjs, [size(jjs)]))
  end subroutine sublattice_case

  !> SLIND, with the lattice (br1, br2, ar1, ar2, rar1, rar2) in COMMON /SL/
  subroutine sublattice_table(ideg, nl1, nl2, t, v, vl, jjs)
    integer, intent(in) :: ideg, nl1, nl2
    real(sp), intent(in) :: t(12)
    real(sp), intent(out) :: v(:, :)
    complex(sp), intent(out) :: vl(:, :)
    integer, intent(out) :: jjs(:, :)
    real(sp) :: br1(2), br2(2), ar1(2), ar2(2), rar1(2), rar2(2)
    integer :: nl1_common, nl2_common
    common /sl/ br1, br2, ar1, ar2, rar1, rar2, nl1_common, nl2_common

    br1 = t(1:2)
    br2 = t(3:4)
    ar1 = t(5:6)
    ar2 = t(7:8)
    rar1 = t(9:10)
    rar2 = t(11:12)
    nl1_common = nl1
    nl2_common = nl2
    v = 0.0_sp
    vl = (0.0_sp, 0.0_sp)
    jjs = 0
    call slind(v, vl, jjs, nl1 * nl2, ideg, 2.0e-4_sp)
  end subroutine sublattice_table

  subroutine read_reals(v)
    real(sp), intent(out) :: v(:)
    integer(int32) :: bits
    integer :: i
    do i = 1, size(v)
      read(5, '(Z8)') bits
      v(i) = transfer(bits, v(i))
    end do
  end subroutine read_reals

  subroutine write_real(label, v)
    character(len=*), intent(in) :: label
    real(sp), intent(in) :: v(:)
    integer :: i
    write(6, '(A,1X,A,1X,I0)') label, 'r4', size(v)
    write(6, '(8Z8.8)') (transfer(v(i), 0_int32), i = 1, size(v))
  end subroutine write_real

  subroutine write_double(label, v)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: v(:)
    integer :: i
    write(6, '(A,1X,A,1X,I0)') label, 'r8', size(v)
    write(6, '(4Z16.16)') (transfer(v(i), 0_int64), i = 1, size(v))
  end subroutine write_double

  subroutine write_int(label, v)
    character(len=*), intent(in) :: label
    integer, intent(in) :: v(:)
    write(6, '(A,1X,A,1X,I0)') label, 'i4', size(v)
    write(6, '(8(I0,1X))') v
  end subroutine write_int
end program refcalc_routines_driver
