!> Calls the routines that set up the angular coefficients (CAAA, CELMG,
!! LXGENT, CPPP) and the sublattices (SLIND) of ref-calc, and those of its
!! energy loop up to the layer diffraction matrices of a Bravais layer
!! (TSCATF, BEAMS, SPHRM_MOD, BEAMS_YLM, FMAT, RSMF_SIMPLE), and the
!! composite layer (SRTLAY and GHMAT; TAUMAT; RTINV_SIMPLE), on the inputs read
!! from stdin and prints their results.
!!
!! Output, one block per array: a header line "name kind count", then the
!! values: kind r4 / r8 as the hex digits of the bit patterns (8 / 16 per
!! value), kind i4 in decimal. stdin:
!!   ANGULAR: name, lmax
!!   SLIND:   name, ideg nl1 nl2, then br1(2) br2(2) ar1(2) ar2(2)
!!            rar1(2) rar2(2) as hex bit patterns, one per line
!!   TSCATF:  name, lmax nel npsi, es(npsi), phss(npsi, nel, lmax + 1), the
!!            energy, the number of calls, then for each call the element,
!!            its concentration and vibration amplitudes perpendicular and
!!            parallel to the surface (the site of
!!            ref-calc.f: calls with concentration > 0, averaged by it)
!!   BEAMS:   name, knbs, the set sizes, npun, npu, tst, spq(2, knt),
!!            spqf(2, knt), the number of energies, then for each energy
!!            e, ak2, ak3
!!   SPHRM:   name, number of points, ct, st, cf of all points (real and
!!            imaginary part of each); runs LMAX = 2 to 18
!!   BYLM:    name, e, vpi, ak2, ak3, lmax, then the beam list as for
!!            BEAMS (npun = 1, npu = 1, one energy): knbs, set sizes, tst,
!!            spq, spqf
!!   FMAT:    name, the lattice as for SLIND, lmax, the number of energies,
!!            then for each e, vpi, ak2, ak3
!!   RSMF:    name, the lattice as for SLIND, lmax, e, vpi, ak2, ak3, tv,
!!            mode (O: overlayer layer, one call, one sublattice; B:
!!            substrate layer, one call per beam set, all sublattices),
!!            the atomic t-matrix (lmax + 1 complex values), then the beam
!!            list as for BEAMS (npun = 1, npu = 1, one energy): knbs,
!!            set sizes, tst, spq, spqf
!!   GHMAT, TAUMAT, RTINV: name, the lattice as for SLIND, lmax, then
!!            mlmax mnsub nlay lay nltype iltype (the sizes MLMAX and MNSUB
!!            of the arrays of ref-calc, which set L2M = 2 * mlmax + 1 and the
!!            leading dimension of the t-matrices; the number of sublayers
!!            of the layer iltype, its lateral periodicity (1 overlayer, 2
!!            bulk), the number of layer types and the layer type), mode (O:
!!            one call for all beams; B: one call per beam set, with the
!!            offset of the set in the beam list), e, vpi, ak2, ak3, tv (the
!!            unit cell area of the layer), tst (the convergence criterion of
!!            GHMAT), the sublayer positions of all layer types (nltype,
!!            mnsub, 3), the number of different t-matrices in use, their
!!            number for each sublayer, the t-matrices (mnsub, mlmax + 1),
!!            then the beam list as for BEAMS (npun = 1, npu = 1, one energy)
!! generate/refcalc_routines_legacy_driver.f is the same driver for the
!! legacy routines of the baseline.
program refcalc_routines_driver
  use, intrinsic :: iso_fortran_env, only: int32, int64
  use angular_coefficients, only: composite_layer_coefficients, lm_permutations, &
                                  matrix_coefficients
  use atomic_scattering, only: legendre_product_integrals, atomic_tmatrix, &
                               accumulate_atomic_tmatrix
  use beam_selection, only: select_beams
  use bravais_layer, only: sublattice_equivalences, lattice_sums, &
                           spherical_harmonics, bravais_layer_matrices, &
                           beams_ylm
  use composite_layer, only: sort_subplanes, interplanar_propagators, &
                             subplane_scattering, copy_subplane_positions, &
                             composite_layer_matrices
  use kinds, only: sp, dp
  use leed_types, only: lattice, energy_point
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
  case ('TSCATF')
    call atomic_tmatrix_case()
  case ('BEAMS')
    call beam_selection_case()
  case ('SPHRM')
    call spherical_harmonics_case()
  case ('BYLM')
    call beams_ylm_case()
  case ('FMAT')
    call lattice_sums_case()
  case ('RSMF')
    call bravais_layer_case()
  case ('GHMAT', 'TAUMAT', 'RTINV')
    call composite_layer_case(trim(name))
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
    integer :: ideg, nl
    real(sp), allocatable :: v(:, :)
    complex(sp), allocatable :: vl(:, :)
    integer, allocatable :: jjs(:, :)
    type(lattice) :: lat

    call read_lattice(lat, nl)
    ideg = lat%symmetry_degree
    allocate(v(nl, 2), vl(nl, 2), jjs(nl, ideg))
    call sublattice_equivalences(lat, nl, v, vl, jjs)
    call write_real('V', reshape(v, [size(v)]))
    call write_real('VL', transfer(vl, [0.0_sp]))
    call write_int('JJS', reshape(jjs, [size(jjs)]))
  end subroutine sublattice_case

  !> TSCATF for the elements of one site, and their concentration-weighted
  !! sum (the ATA average of the atomic t-matrix), as the energy loop of
  !! ref-calc.f does it: T0 = T = 100, DR0 = VSITE = 0.
  subroutine atomic_tmatrix_case()
    integer :: lmax, l1, nel, npsi, nn1, nn2, nn3, ncall, k, iel
    real(sp) :: eb, dr0, vsite, t0, t, drper, drpar, conc, vib
    real(sp), allocatable :: es(:), phss(:, :, :), phs(:), ppp(:, :, :)
    complex(sp), allocatable :: af(:), caf(:), del(:), ctab(:), sum_(:), tmat(:)
    complex(dp), allocatable :: bj(:)

    read(5, *) lmax, nel, npsi
    l1 = lmax + 1
    nn3 = l1
    nn2 = l1
    nn1 = nn2 + nn3 - 1
    allocate(es(npsi), phss(npsi, nel, l1), phs(l1), ppp(nn1, nn2, nn3))
    allocate(af(l1), caf(l1), del(l1), ctab(nn3), sum_(nn2), tmat(l1))
    allocate(bj(nn1))
    call read_reals(es)
    call read_array(phss)
    eb = read_real()
    call gaunt_table(lmax, ppp)
    dr0 = 0.0_sp
    vsite = 0.0_sp
    t0 = 100.0_sp
    t = 100.0_sp
    tmat = (0.0_sp, 0.0_sp)
    read(5, *) ncall
    do k = 1, ncall
      read(5, *) iel
      conc = read_real()
      vib = read_real()
      drpar = read_real()
      if (conc > 0.0_sp) then
        drper = vib
        call atomic_tmatrix(iel, l1, npsi, nel, nn1, nn2, nn3, es, phss, eb, vsite, &
                            ppp, dr0, drper, drpar, t0, t, af, caf, phs, del, &
                            ctab, sum_, bj)
        call write_real(numbered('PHS', k), phs)
        call write_complex(numbered('AF', k), af)
        call write_complex(numbered('DEL', k), del)
        call write_complex(numbered('CAF', k), caf)
        call accumulate_atomic_tmatrix(conc, caf, tmat)
      end if
    end do
    call write_complex('TMAT', tmat)
  end subroutine atomic_tmatrix_case

  !> BEAMS on one beam list at several energies.
  subroutine beam_selection_case()
    integer :: knbs, knt, npun, ne, ie, nt, np, mpu
    integer, allocatable :: knb(:), npu(:), nb(:), npuc(:)
    real(sp), allocatable :: spq(:, :), spqf(:, :), pq(:, :), pqf(:, :)
    real(sp) :: tst, e, ak2, ak3

    call read_beam_list(knbs, knb, knt, tst, spq, spqf, npu)
    npun = size(npu)
    allocate(nb(knbs), npuc(npun), pq(2, knt), pqf(2, knt))
    read(5, *) ne
    ! npuc is not reset between energies, as in ref-calc: the routine must
    ! overwrite every entry it reports
    npuc = -7
    do ie = 1, ne
      e = read_real()
      ak2 = read_real()
      ak3 = read_real()
      nb = 0
      pq = 0.0_sp
      pqf = 0.0_sp
      call select_beams(knbs, knb, spq, spqf, knt, npu, npun, ak2, ak3, e, tst, &
                        nb, pq, pqf, npuc, mpu, nt, np)
      call write_int(numbered('NB', ie), nb)
      call write_int(numbered('NT_NP_MPU', ie), [nt, np, mpu])
      call write_int(numbered('NPUC', ie), npuc)
      call write_real(numbered('PQ', ie), reshape(pq(:, 1:nt), [2 * nt]))
      call write_real(numbered('PQF', ie), reshape(pqf(:, 1:nt), [2 * nt]))
    end do
  end subroutine beam_selection_case

  !> SPHRM_MOD on a grid of CT, ST, CF for every LMAX from 2 to 18.
  subroutine spherical_harmonics_case()
    integer :: npts, lmax, lmmax, i
    complex(sp), allocatable :: ct(:), st(:), cf(:), ylm(:), all_ylm(:)

    read(5, *) npts
    allocate(ct(npts), st(npts), cf(npts))
    call read_complexes(ct)
    call read_complexes(st)
    call read_complexes(cf)
    do lmax = 2, 18
      lmmax = (lmax + 1)**2
      allocate(ylm(lmmax), all_ylm(lmmax * npts))
      do i = 1, npts
        ylm = (0.0_sp, 0.0_sp)
        call spherical_harmonics(lmax, ylm, lmmax, ct(i), st(i), cf(i))
        all_ylm((i - 1) * lmmax + 1:i * lmmax) = ylm
      end do
      call write_complex(numbered('YLM_L', lmax), all_ylm)
      deallocate(ylm, all_ylm)
    end do
  end subroutine spherical_harmonics_case

  !> Read the beam list as for BEAMS: set count and sizes, the sets of
  !! incident-beam indices (only if npu is present), the criterion tst, and
  !! the beam momenta and index pairs.
  subroutine read_beam_list(knbs, knb, knt, tst, spq, spqf, npu)
    integer, intent(out) :: knbs, knt
    integer, allocatable, intent(out) :: knb(:)
    real(sp), intent(out) :: tst
    real(sp), allocatable, intent(out) :: spq(:, :), spqf(:, :)
    integer, allocatable, intent(out), optional :: npu(:)
    integer :: npun

    read(5, *) knbs
    allocate(knb(knbs))
    read(5, *) knb
    knt = sum(knb)
    if (present(npu)) then
      read(5, *) npun
      allocate(npu(npun))
      read(5, *) npu
    end if
    tst = read_real()
    allocate(spq(2, knt), spqf(2, knt))
    call read_array(spq)
    call read_array(spqf)
  end subroutine read_beam_list

  !> Read a lattice as for SLIND and lmax, and set up the sublattices.
  subroutine set_up_lattice(lat, nl, lmax, klm, v, jjs, flms)
    type(lattice), intent(out) :: lat
    integer, intent(out) :: nl, lmax, klm
    real(sp), allocatable, intent(out) :: v(:, :)
    integer, allocatable, intent(out) :: jjs(:, :)
    complex(sp), allocatable, intent(out) :: flms(:, :)
    complex(sp), allocatable :: vl(:, :)

    call read_lattice(lat, nl)
    read(5, *) lmax
    klm = (2 * lmax + 1) * (2 * lmax + 2) / 2
    allocate(v(nl, 2), vl(nl, 2), jjs(nl, lat%symmetry_degree), flms(nl, klm))
    call sublattice_equivalences(lat, nl, v, vl, jjs)
  end subroutine set_up_lattice

  !> FMAT at one energy point, with the cutoff distance of ref-calc.f.
  subroutine run_lattice_sums(point, lat, nl, lmax, klm, v, jjs, flms)
    type(energy_point), intent(in) :: point
    type(lattice), intent(in) :: lat
    integer, intent(in) :: nl, lmax, klm, jjs(:, :)
    real(sp), intent(in) :: v(:, :)
    complex(sp), intent(out) :: flms(:, :)
    real(sp) :: dcut

    dcut = -5.0_sp * sqrt(2.0_sp * point%energy) &
           / (min(point%imaginary_potential, -0.05_sp))
    flms = (0.0_sp, 0.0_sp)
    call lattice_sums(flms, v, jjs, nl, nl, dcut, lat, lmax, klm, point)
  end subroutine run_lattice_sums

  !> Read a lattice as for SLIND, and return the number of sublattices.
  subroutine read_lattice(lat, nl)
    type(lattice), intent(out) :: lat
    integer, intent(out) :: nl
    real(sp) :: t(12)

    read(5, *) lat%symmetry_degree, lat%repetitions
    call read_reals(t)
    lat%substrate_vector(:, 1) = t(1:2)
    lat%substrate_vector(:, 2) = t(3:4)
    lat%superlattice_vector(:, 1) = t(5:6)
    lat%superlattice_vector(:, 2) = t(7:8)
    lat%superlattice_reciprocal(:, 1) = t(9:10)
    lat%superlattice_reciprocal(:, 2) = t(11:12)
    nl = lat%repetitions(1) * lat%repetitions(2)
  end subroutine read_lattice

  !> BEAMS_YLM for each beam set of a beam list (with its offset in the
  !! list), and for a window of the list that starts inside the first set
  !! and ends before the last beam. CYLM is reset to a sentinel before every
  !! call and written in full after it, so that rows that the routine must
  !! leave alone show.
  subroutine beams_ylm_case()
    type(energy_point) :: point
    integer :: lmax, lmmax, knbs, knt, npun, nt, np, mpu, ibs, na, n
    real(sp) :: tst
    integer, allocatable :: knb(:), npu(:), nb(:), npuc(:)
    real(sp), allocatable :: spq(:, :), spqf(:, :), pq(:, :), pqf(:, :)
    complex(sp), allocatable :: xa(:), cylm(:, :), ylm(:)
    complex(sp), parameter :: sentinel = (7.0_sp, -3.0_sp)

    point%energy = read_real()
    point%imaginary_potential = read_real()
    point%k_parallel(1) = read_real()
    point%k_parallel(2) = read_real()
    read(5, *) lmax
    lmmax = (lmax + 1)**2
    call read_beam_list(knbs, knb, knt, tst, spq, spqf)
    npun = 1
    allocate(nb(knbs), npu(npun), npuc(npun), pq(2, knt), pqf(2, knt))
    npu = 1
    pq = 0.0_sp
    pqf = 0.0_sp
    nb = 0
    npuc = -7
    call select_beams(knbs, knb, spq, spqf, knt, npu, npun, point%k_parallel(1), &
                      point%k_parallel(2), point%energy, tst, nb, pq, pqf, &
                      npuc, mpu, nt, np)
    allocate(cylm(nt, lmmax), ylm(lmmax))
    na = 0
    do ibs = 1, knbs
      allocate(xa(nb(ibs)))
      cylm = sentinel
      call beams_ylm(cylm, xa, pq, nt, nb(ibs), na, point, lmax, lmmax, ylm)
      call write_complex(numbered('XA', ibs), xa)
      call write_complex(numbered('CYLM', ibs), reshape(cylm, [nt * lmmax]))
      deallocate(xa)
      na = na + nb(ibs)
    end do
    call write_int('NB', nb)
    n = nt - 3
    allocate(xa(n))
    cylm = sentinel
    call beams_ylm(cylm, xa, pq, nt, n, 1, point, lmax, lmmax, ylm)
    call write_complex('XA_WINDOW', xa)
    call write_complex('CYLM_WINDOW', reshape(cylm, [nt * lmmax]))
  end subroutine beams_ylm_case

  !> FMAT on the sublattices of the lattice at several (E, VPI, k-parallel).
  subroutine lattice_sums_case()
    type(lattice) :: lat
    integer :: nl, lmax, klm, ne, ie
    type(energy_point) :: point
    real(sp), allocatable :: v(:, :)
    complex(sp), allocatable :: flms(:, :)
    integer, allocatable :: jjs(:, :)

    call set_up_lattice(lat, nl, lmax, klm, v, jjs, flms)
    read(5, *) ne
    do ie = 1, ne
      point%energy = read_real()
      point%imaginary_potential = read_real()
      point%k_parallel(1) = read_real()
      point%k_parallel(2) = read_real()
      call run_lattice_sums(point, lat, nl, lmax, klm, v, jjs, flms)
      call write_complex(numbered('FLMS', ie), reshape(flms, [nl * klm]))
    end do
  end subroutine lattice_sums_case

  !> BEAMS, FMAT and RSMF_SIMPLE chained as the energy loop of ref-calc.f
  !! does: for an overlayer-type layer (mode O: one call for all beams, one
  !! sublattice) or a substrate-type layer (mode B: one call per beam set
  !! with its offset in the beam list, all sublattices).
  subroutine bravais_layer_case()
    type(lattice) :: lat
    type(energy_point) :: point
    integer :: nl, lmax, l1, lmmax, klm, lev, lod, nlm
    integer :: knbs, knt, npun, nt, np, mpu, ibs, na, naa, nll
    real(sp) :: e, ak2, ak3, vpi, tv, tst
    character(len=1) :: mode
    integer, allocatable :: knb(:), npu(:), nb(:), npuc(:), jjs(:, :)
    integer, allocatable :: lx(:), lxi(:), lt(:), lxm(:), iple(:), iplo(:)
    real(sp), allocatable :: spq(:, :), spqf(:, :), pq(:, :), pqf(:, :)
    real(sp), allocatable :: v(:, :), clm(:)
    complex(sp), allocatable :: flms(:, :)
    complex(sp), allocatable :: tmat(:), ra(:, :), ta(:, :), amult(:)
    complex(sp), allocatable :: cylm(:, :), flm(:), xev(:, :), xod(:, :)
    complex(sp), allocatable :: ylm(:), ylme(:), ylmo(:)
    complex(sp), allocatable :: xevst(:, :), xodst(:, :)

    call set_up_lattice(lat, nl, lmax, klm, v, jjs, flms)
    l1 = lmax + 1
    lmmax = l1 * l1
    lev = l1 * (l1 + 1) / 2
    lod = lmmax - lev
    nlm = nlms(lmax)
    e = read_real()
    vpi = read_real()
    ak2 = read_real()
    ak3 = read_real()
    tv = read_real()
    read(5, *) mode
    allocate(tmat(l1))
    call read_complexes(tmat)
    call read_beam_list(knbs, knb, knt, tst, spq, spqf)
    npun = 1
    allocate(nb(knbs), npu(npun), npuc(npun), pq(2, knt), pqf(2, knt))
    npu = 1

    allocate(clm(nlm), lx(lmmax), lxi(lmmax), lt(lmmax), lxm(lmmax))
    call x_tables(lmax, clm, lx, lxi, lt, lxm)

    pq = 0.0_sp
    pqf = 0.0_sp
    nb = 0
    npuc = -7
    call select_beams(knbs, knb, spq, spqf, knt, npu, npun, ak2, ak3, e, tst, &
                      nb, pq, pqf, npuc, mpu, nt, np)
    point = energy_point(e, [ak2, ak3], vpi)
    call run_lattice_sums(point, lat, nl, lmax, klm, v, jjs, flms)

    allocate(ra(nt, nt), ta(nt, nt), amult(nt), cylm(nt, lmmax))
    allocate(flm(klm), xev(lev, lev), xod(lod, lod), ylm(lmmax))
    allocate(ylme(lev), ylmo(lod), iple(lev), iplo(lod))
    allocate(xevst(lev, lev), xodst(lod, lod))
    ra = (0.0_sp, 0.0_sp)
    ta = (0.0_sp, 0.0_sp)
    amult = (0.0_sp, 0.0_sp)
    cylm = (0.0_sp, 0.0_sp)
    xevst = (0.0_sp, 0.0_sp)
    xodst = (0.0_sp, 0.0_sp)
    if (mode == 'O') then
      nll = 1
      naa = nt
      call bravais_layer_matrices(ra, ta, naa, amult, cylm, pq, nt, flms, flm, &
                                  v, nl, 0, nll, tmat, tmat, lmax, lx, lxi, &
                                  lmmax, klm, xev, xod, lev, lod, ylm, ylme, &
                                  ylmo, iple, iplo, clm, nlm, xevst, xodst, &
                                  point, tv)
    else
      nll = nl
      na = 0
      do ibs = 1, knbs
        naa = nb(ibs)
        call bravais_layer_matrices(ra, ta, naa, amult, cylm, pq, nt, flms, &
                                    flm, v, nl, na, nll, tmat, tmat, lmax, lx, &
                                    lxi, lmmax, klm, xev, xod, lev, lod, ylm, &
                                    ylme, ylmo, iple, iplo, clm, nlm, xevst, &
                                    xodst, point, tv)
        na = na + nb(ibs)
      end do
    end if
    call write_int('NB', nb)
    call write_complex('RA', reshape(ra, [nt * nt]))
    call write_complex('TA', reshape(ta, [nt * nt]))
    call write_complex('AMULT', amult)
    call write_complex('CYLM', reshape(cylm, [nt * lmmax]))
    call write_complex('XEVST', reshape(xevst, [lev * lev]))
    call write_complex('XODST', reshape(xodst, [lod * lod]))
  end subroutine bravais_layer_case

  !> SRTLAY and GHMAT, TAUMAT, or RTINV_SIMPLE on a composite layer, with the
  !! arguments, array sizes and COMMON blocks of ref-calc.f. Mode O: one call
  !! for all beams (offsets na = ns = 0); mode B: one call per beam set, with
  !! the offset of the set in the beam list as na and the sum of the earlier
  !! set sizes as ns.
  subroutine composite_layer_case(routine)
    character(len=*), intent(in) :: routine
    type(lattice) :: lat
    type(energy_point) :: point
    integer :: nl, lmax, klm, l1, lmmax, lev, lod, lev2, nlm, ncaa
    integer :: mlmax, mnsub, nlay, layt, nltype, iltype, nnsub, lmax1
    integer :: nlay2, lmg, lmn, lm2n, lmni, l2m, lms, lmt, ntauv
    integer :: knbs, knt, npun, nt, np, mpu, ncalls, k, naa, mlmmax
    real(sp) :: tstb, tsts, tv1, dcuts
    integer :: lay, lm, ntau, na, ns
    character(len=1) :: mode
    character(len=16) :: label
    integer, allocatable :: knb(:), npu(:), nb(:), npuc(:), jjs(:, :)
    integer, allocatable :: lx(:), lxi(:), lt(:), lxm(:), lps(:), lpss(:)
    integer, allocatable :: nord(:), mgh(:, :), nugh(:), ngeq(:), ipl(:)
    real(sp), allocatable :: spq(:, :), spqf(:, :), pq(:, :), pqf(:, :)
    real(sp), allocatable :: v(:, :), clm(:), caa(:), subpos(:, :, :)
    real(sp), allocatable :: pos(:, :), poss(:, :), drl(:, :), test(:)
    complex(sp), allocatable :: flms(:, :), tsf(:, :), tau(:, :), gh(:, :)
    complex(sp), allocatable :: th(:, :), flm(:), xev(:, :), xh(:), hghd(:)
    complex(sp), allocatable :: amult(:), cylm(:, :), rtop(:, :), ttop(:, :)
    complex(sp), allocatable :: rbot(:, :), tbot(:, :), rg(:, :, :), rgp(:)
    complex(sp), allocatable :: ts(:), tg(:, :), vt(:), taug(:), taugm(:)
    complex(sp), allocatable :: rg1(:), rg2(:), tstore(:, :, :), buffer(:)
    ! the COMMON blocks of the legacy routines: blank (E, AK2, AK3, VPI, TV,
    ! EMACH), /SL/ (the lattice) and /MS/
    real(sp) :: e_c, ak2_c, ak3_c, vpi_c, tv_c, emach_c
    real(sp) :: sl_vectors(2, 6), eps_c
    integer :: sl_repetitions(2), lmax_c, liter_c
    common e_c, ak2_c, ak3_c, vpi_c, tv_c, emach_c
    common /sl/ sl_vectors, sl_repetitions
    common /ms/ lmax_c, eps_c, liter_c

    call set_up_lattice(lat, nl, lmax, klm, v, jjs, flms)
    sl_vectors(:, 1:2) = lat%substrate_vector
    sl_vectors(:, 3:4) = lat%superlattice_vector
    sl_vectors(:, 5:6) = lat%superlattice_reciprocal
    sl_repetitions = lat%repetitions
    read(5, *) mlmax, mnsub, nlay, layt, nltype, iltype
    read(5, *) mode
    e_c = read_real()
    vpi_c = read_real()
    ak2_c = read_real()
    ak3_c = read_real()
    tv1 = read_real()
    tsts = read_real()
    tv_c = tv1
    emach_c = 1.0e-6_sp
    lmax_c = lmax
    eps_c = 0.001_sp
    liter_c = 10
    l1 = lmax + 1
    lmmax = l1 * l1
    lev = l1 * (l1 + 1) / 2
    lod = lmmax - lev
    lev2 = 2 * lev
    nlm = nlms(lmax)
    ncaa = nca(lmax)
    nnsub = mnsub
    lmax1 = mlmax + 1
    mlmmax = lmax1 * lmax1
    l2m = 2 * mlmax + 1
    lms = l2m * l2m
    nlay2 = nlay * (nlay - 1) / 2
    lmg = nlay2 * lmmax * 2
    lmn = nlay * lmmax
    lm2n = 2 * lmn
    lmni = lmn
    point = energy_point(e_c, [ak2_c, ak3_c], vpi_c)
    dcuts = -5.0_sp * sqrt(2.0_sp * e_c) / (min(vpi_c, -0.05_sp))

    allocate(subpos(nltype, mnsub, 3), pos(nlay, 3), poss(nlay, 3))
    call read_array(subpos)
    read(5, *) ntauv
    allocate(lps(nlay), tsf(mnsub, lmax1), buffer(mnsub * lmax1))
    read(5, *) lps
    call read_complexes(buffer)
    tsf = reshape(buffer, [mnsub, lmax1])
    call read_beam_list(knbs, knb, knt, tstb, spq, spqf)
    npun = 1
    allocate(nb(knbs), npu(npun), npuc(npun), pq(2, knt), pqf(2, knt))
    npu = 1
    pq = 0.0_sp
    pqf = 0.0_sp
    nb = 0
    npuc = -7
    call select_beams(knbs, knb, spq, spqf, knt, npu, npun, ak2_c, ak3_c, e_c, &
                      tstb, nb, pq, pqf, npuc, mpu, nt, np)
    call run_lattice_sums(point, lat, nl, lmax, klm, v, jjs, flms)

    allocate(clm(nlm), lx(lmmax), lxi(lmmax), lt(lmmax), lxm(lmmax))
    call x_tables(lmax, clm, lx, lxi, lt, lxm)
    allocate(caa(ncaa))
    call composite_layer_coefficients(lmax, caa)
    call copy_subplane_positions(pos, subpos, nltype, nnsub, nlay, iltype)

    lay = layt
    lm = l1
    ntau = ntauv
    lmt = ntauv * lmmax
    allocate(mgh(nlay, nlay), drl(nlay2, 3))
    allocate(nugh(nlay2), ngeq(nlay2), test(nlay2))
    allocate(lpss(nlay), nord(nlay), gh(lmg, lmmax), hghd(l2m))
    allocate(tg(2, 2 * mnsub * mlmmax), vt(2 * mnsub * mlmmax))
    allocate(tau(lmt, lev), xev(lev, lev2), xh(lev), flm(klm))
    mgh = 0
    drl = 0.0_sp
    nugh = 0
    ngeq = 0
    test = 0.0_sp
    lpss = 0
    nord = 0
    gh = (0.0_sp, 0.0_sp)
    hghd = (0.0_sp, 0.0_sp)
    tg = (0.0_sp, 0.0_sp)
    vt = (0.0_sp, 0.0_sp)
    tau = (0.0_sp, 0.0_sp)
    xev = (0.0_sp, 0.0_sp)
    xh = (0.0_sp, 0.0_sp)
    flm = (0.0_sp, 0.0_sp)
    if (routine == 'RTINV') then
      allocate(th(lmni, lmni), ipl(lmni), amult(nt), cylm(nt, lmmax))
      allocate(rtop(nt, nt), ttop(nt, nt), rbot(nt, nt), tbot(nt, nt))
      allocate(rg(2, nlay, nt), rgp(nt), ts(lmn), taug(lmt))
      allocate(taugm(lmt), rg1(nlay), rg2(nlay), tstore(2, lmn, nt))
      th = (0.0_sp, 0.0_sp)
      ipl = 0
      amult = (0.0_sp, 0.0_sp)
      cylm = (0.0_sp, 0.0_sp)
      rtop = (0.0_sp, 0.0_sp)
      ttop = (0.0_sp, 0.0_sp)
      rbot = (0.0_sp, 0.0_sp)
      tbot = (0.0_sp, 0.0_sp)
      rg = (0.0_sp, 0.0_sp)
      rgp = (0.0_sp, 0.0_sp)
      ts = (0.0_sp, 0.0_sp)
      taug = (0.0_sp, 0.0_sp)
      taugm = (0.0_sp, 0.0_sp)
      rg1 = (0.0_sp, 0.0_sp)
      rg2 = (0.0_sp, 0.0_sp)
      tstore = (0.0_sp, 0.0_sp)
    end if

    ncalls = merge(1, knbs, mode == 'O')
    na = 0
    ns = 0
    do k = 1, ncalls
      naa = merge(nt, nb(k), mode == 'O')
      select case (routine)
      case ('GHMAT')
        call sort_subplanes(pos, poss, lps, lpss, mgh, nlay, drl, nlay2, nugh, &
                            ngeq, nord)
        call interplanar_propagators(gh, lmg, lmmax, nugh, ngeq, nlay2, tsts, &
                                     test, vt, l2m, lm, tg, lms, drl, tv1, &
                                     lxm, dcuts, caa, ncaa, lay, hghd, &
                                     pq(:, 1 + na), lat, point)
        call write_complex(numbered('GH', k), reshape(gh, [lmg * lmmax]))
      case ('TAUMAT')
        call subplane_scattering(tau, lmt, ntau, xev, lev, lev2, lod, tsf, &
                                 lmmax, flms, nl, klm, lm, clm, nlm, lxi, &
                                 pq(:, 1 + na), nnsub, lmax1, xh, flm, &
                                 lay, lat, point)
        call write_complex(numbered('TAU', k), reshape(tau, [lmt * lev]))
        call write_complex(numbered('FLM', k), flm)
      case default
        call composite_layer_matrices(rtop, ttop, rbot, tbot, naa, amult, &
            cylm, pq, nt, flms, flm, nl, lxi, lt, lxm, lmmax, klm, xev, lev, &
            lev2, tau, lmt, taug, taugm, clm, nlm, pos, poss, mgh, nlay, drl, &
            nlay2, nugh, ngeq, test, gh, lmg, rg, rgp, ts, lmn, tg, lm2n, vt, &
            caa, ncaa, th, lmni, ipl, tsf, lps, lpss, nord, nnsub, lmax1, xh, &
            hghd, l2m, rg1, rg2, tstore, na, ns, lay, ntau, tsts, tv1, dcuts, &
            lmax, lat, point)
      end select
      if (mode == 'B') then
        na = na + naa
        ns = ns + naa
      end if
    end do

    call write_int('NB', nb)
    select case (routine)
    case ('GHMAT')
      call write_int('MGH', reshape(mgh, [nlay * nlay]))
      call write_real('DRL', reshape(drl, [nlay2 * 3]))
      call write_real('POSS', reshape(poss, [nlay * 3]))
      call write_int('LPSS', lpss)
      call write_int('NORD', nord)
      call write_int('NUGH', nugh)
      call write_int('NGEQ', ngeq)
      call write_real('TEST', test)
    case ('RTINV')
      call write_complex('R_TOP', reshape(rtop, [nt * nt]))
      call write_complex('T_TOP', reshape(ttop, [nt * nt]))
      call write_complex('R_BOT', reshape(rbot, [nt * nt]))
      call write_complex('T_BOT', reshape(tbot, [nt * nt]))
      call write_complex('TSTORE', reshape(tstore, [2 * lmn * nt]))
      call write_complex('TAU', reshape(tau, [lmt * lev]))
      call write_complex('GH', reshape(gh, [lmg * lmmax]))
      call write_real('POSS', reshape(poss, [nlay * 3]))
      call write_int('NORD', nord)
      call write_complex('AMULT', amult)
      call write_complex('CYLM', reshape(cylm, [nt * lmmax]))
    end select
  end subroutine composite_layer_case

  !> The Gaunt coefficients PPP as ref-calc.f computes them.
  subroutine gaunt_table(lmax, ppp)
    integer, intent(in) :: lmax
    real(sp), intent(out) :: ppp(:, :, :)
    real(sp) :: clm(nlms(lmax))
    real(dp) :: fac(nfac)
    fac = 0.0_dp
    call matrix_coefficients(lmax, clm, fac)
    call legendre_product_integrals(size(ppp, 1), size(ppp, 2), size(ppp, 3), &
                                    ppp, fac)
  end subroutine gaunt_table

  !> CLM and the (L,M) permutations LX, LXI, LT and LXM as ref-calc.f computes them.
  subroutine x_tables(lmax, clm, lx, lxi, lt, lxm)
    integer, intent(in) :: lmax
    real(sp), intent(out) :: clm(:)
    integer, intent(out) :: lx(:), lxi(:), lt(:), lxm(:)
    real(dp) :: fac(nfac)
    fac = 0.0_dp
    clm = 0.0_sp
    call matrix_coefficients(lmax, clm, fac)
    call lm_permutations(lmax, lx, lxi, lt, lxm)
  end subroutine x_tables

  !> label with the two-digit number k appended after an underscore
  function numbered(label, k) result(text)
    character(len=*), intent(in) :: label
    integer, intent(in) :: k
    character(len=16) :: text
    write(text, '(A,A,I2.2)') label, '_', k
  end function numbered

  function read_real() result(x)
    real(sp) :: x
    real(sp) :: buffer(1)
    call read_reals(buffer)
    x = buffer(1)
  end function read_real

  !> Read the values of a (rank 2 or 3) real array in storage order.
  subroutine read_array(a)
    real(sp), intent(out) :: a(..)
    real(sp), allocatable :: buffer(:)
    allocate(buffer(size(a)))
    call read_reals(buffer)
    select rank (a)
    rank (2)
      a = reshape(buffer, shape(a))
    rank (3)
      a = reshape(buffer, shape(a))
    end select
  end subroutine read_array

  subroutine read_complexes(c)
    complex(sp), intent(out) :: c(:)
    real(sp) :: t(2 * size(c))
    call read_reals(t)
    c = cmplx(t(1::2), t(2::2), kind=sp)
  end subroutine read_complexes

  subroutine write_complex(label, c)
    character(len=*), intent(in) :: label
    complex(sp), intent(in) :: c(:)
    call write_real(label, transfer(c, [0.0_sp]))
  end subroutine write_complex

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
