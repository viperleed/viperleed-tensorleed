!> Calls the routines that set up the angular coefficients (CAAA, CELMG,
!! LXGENT, CPPP) and the sublattices (SLIND) of ref-calc on the inputs read
!! from stdin and prints their results.
!!
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
  use angular_coefficients, only: composite_layer_coefficients, lm_permutations, &
                                  matrix_coefficients
  use atomic_scattering, only: legendre_product_integrals
  use bravais_layer, only: sublattice_equivalences
  use kinds, only: sp, dp
  use leed_types, only: lattice
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

    ! CAAA: Clebsch-Gordan coefficients of the composite layer
    allocate(caa(nca(lmax)))
    call composite_layer_coefficients(lmax, caa)
    call write_real('CAA', caa)

    ! CELMG: Clebsch-Gordan coefficients of the matrices X and TAU; it sets
    ! only fac(1:4*lmax+2), the rest of the printed fac stays zero
    allocate(clm(nlms(lmax)))
    fac = 0.0_dp
    call matrix_coefficients(lmax, clm, fac)
    call write_real('CLM', clm)
    call write_double('FAC_CELMG', fac)

    ! LXGENT: permutations of the (L,M) sequence
    allocate(lx(lmmax), lxi(lmmax), lt(lmmax), lxm(lmmax))
    call lm_permutations(lmax, lx, lxi, lt, lxm)
    call write_int('LX', lx)
    call write_int('LXI', lxi)
    call write_int('LT', lt)
    call write_int('LXM', lxm)

    ! CPPP: Gaunt coefficients for the temperature-dependent phase shifts
    allocate(ppp(2 * lmax + 1, lmax + 1, lmax + 1))
    call legendre_product_integrals(size(ppp, 1), size(ppp, 2), size(ppp, 3), &
                                    ppp, fac)
    call write_real('PPP', reshape(ppp, [size(ppp)]))
    call write_double('FAC_CPPP', fac)
  end subroutine angular_coefficients_case

  !> SLIND: the lattice is filled from the stdin values br1, br2, ar1, ar2,
  !! rar1, rar2
  subroutine sublattice_case()
    integer :: ideg, nl1, nl2, nl
    real(sp) :: t(12)
    real(sp), allocatable :: v(:, :)
    complex(sp), allocatable :: vl(:, :)
    integer, allocatable :: jjs(:, :)
    type(lattice) :: lat

    read(5, *) ideg, nl1, nl2
    call read_reals(t)
    nl = nl1 * nl2
    lat%substrate_vector(:, 1) = t(1:2)
    lat%substrate_vector(:, 2) = t(3:4)
    lat%superlattice_vector(:, 1) = t(5:6)
    lat%superlattice_vector(:, 2) = t(7:8)
    lat%superlattice_reciprocal(:, 1) = t(9:10)
    lat%superlattice_reciprocal(:, 2) = t(11:12)
    lat%symmetry_degree = ideg
    lat%repetitions = [nl1, nl2]
    allocate(v(nl, 2), vl(nl, 2), jjs(nl, ideg))
    call sublattice_equivalences(lat, nl, v, vl, jjs)
    call write_real('V', reshape(v, [size(v)]))
    call write_real('VL', transfer(vl, [0.0_sp]))
    call write_int('JJS', reshape(jjs, [size(jjs)]))
  end subroutine sublattice_case

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
