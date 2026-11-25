!> \file bravais_layer.f90
!! \brief Lattice sums and scattering of Bravais layers: sublattice setup,
!! lattice sums, spherical harmonics.
!!
!! Each routine keeps the precision of every quantity of the legacy routine
!! it follows, and its results are bit-for-bit those of the legacy code
!! (tests/tier1/test_refcalc_routines.py).

!> Setup and lattice sums of Bravais layers, and spherical harmonics.
module bravais_layer
  use constants, only: pi
  use kinds, only: sp
  use leed_output, only: stop_with_error
  use leed_types, only: lattice, energy_point
  implicit none

  private

  public :: sublattice_equivalences, lattice_sums, spherical_harmonics
  public :: bravais_layer_matrices, beams_ylm

  !> Small number to compare the phase factors of sublattice vectors
  !! (the legacy argument EPSD).
  real(sp), parameter :: epsd = 2.0e-4_sp

contains

  !> Sets up the sublattice vectors and their equivalences under rotation
  !! (legacy SLIND).
  !!
  !! Finds how the sublattices js of the superlattice are transformed into
  !! one another by rotations through j * 2 pi / ideg: jjs(js, j) is the
  !! sublattice into which sublattice js is rotated by the j-th rotation.
  !! If no sublattice matches (ideg does not fit the substrate lattice and
  !! the superlattice; possible only for nl > 1), the routine stops with an
  !! error. Author: Pendry.
  !!
  !! The substrate lattice vectors define the sublattices, the reciprocal
  !! vectors of the superlattice identify them. The arrays v, vl and jjs are
  !! explicitly shaped: their shapes (nl, 2) and (nl, ideg) equal the PARAM
  !! shapes of the ref-calc actuals, because the input reader requires NL1,
  !! NL2 and IDEG to equal their PARAM values (sizing from the input comes
  !! later).
  subroutine sublattice_equivalences(lat, nl, v, vl, jjs)
    !> Lattices; uses substrate_vector, superlattice_reciprocal,
    !! symmetry_degree (ideg) and repetitions (nl1, nl2)
    type(lattice), intent(in) :: lat
    !> Number of sublattices; nl1 * nl2
    integer, intent(in) :: nl
    !> Vectors pointing to the sublattices, in polar coordinates (modulus,
    !! angle)
    real(sp), intent(out) :: v(nl, 2)
    !> Phase factors of the sublattice vectors for the two reciprocal vectors
    !! of the superlattice
    complex(sp), intent(out) :: vl(nl, 2)
    !> Relationship of the sublattices under rotation
    integer, intent(out) :: jjs(nl, lat%symmetry_degree)

    integer :: i, j, k, ideg, jl1, kl2, nl1, nl2
    real(sp) :: s1, s2, adr1, adr2, adr, ang, ainc, cang, sang, a, b, test
    complex(sp) :: vla, vlb, ci
    character(len=200) :: message
    ! Legacy names: br1, br2 = substrate_vector; rar1, rar2 =
    ! superlattice_reciprocal
    real(sp) :: br1(2), br2(2), rar1(2), rar2(2)

    br1 = lat%substrate_vector(:, 1)
    br2 = lat%substrate_vector(:, 2)
    rar1 = lat%superlattice_reciprocal(:, 1)
    rar2 = lat%superlattice_reciprocal(:, 2)
    ideg = lat%symmetry_degree
    nl1 = lat%repetitions(1)
    nl2 = lat%repetitions(2)

    ci = cmplx(0.0_sp, 1.0_sp, sp)
    ! Set up vectors v defining sublattices and quantities vl for later
    ! reference.
    i = 1
    s1 = 0.0_sp
    do jl1 = 1, nl1
      s2 = 0.0_sp
      do kl2 = 1, nl2
        adr1 = s1 * br1(1) + s2 * br2(1)
        adr2 = s1 * br1(2) + s2 * br2(2)
        v(i, 1) = sqrt(adr1 * adr1 + adr2 * adr2)
        v(i, 2) = 0.0_sp
        if (v(i, 1) > 0.0_sp) v(i, 2) = atan2(adr2, adr1)
        vl(i, 1) = exp(ci * (adr1 * rar1(1) + adr2 * rar1(2)))
        vl(i, 2) = exp(ci * (adr1 * rar2(1) + adr2 * rar2(2)))
        i = i + 1
        s2 = s2 + 1.0_sp
      end do
      s1 = s1 + 1.0_sp
    end do

    ! Rotate each vector v and find to which v it becomes equivalent in
    ! terms of the quantities vl. This equivalence means belonging to the
    ! same sublattice.
    ainc = 2.0_sp * pi / real(ideg, sp)
    do i = 1, nl
      adr = v(i, 1)
      ang = v(i, 2)
      do k = 1, ideg
        ang = ang + ainc
        cang = cos(ang)
        sang = sin(ang)
        a = rar1(1) * cang + rar1(2) * sang
        b = rar2(1) * cang + rar2(2) * sang
        vla = exp(ci * adr * a)
        vlb = exp(ci * adr * b)
        do j = 1, nl
          test = abs(vla - vl(j, 1)) + abs(vlb - vl(j, 2))
          if (test - 5.0_sp * epsd <= 0.0_sp) exit
        end do
        if (j > nl) then
          write (message, '(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'IDEG = ', ideg, &
            ' does not fit the substrate lattice and the superlattice: ' &
            // 'sublattice ', i, ' rotated by ', k, ' x 360/', ideg, &
            ' degrees matches none of the ', nl, ' sublattices'
          call stop_with_error(trim(message))
        end if
        jjs(i, k) = j
      end do
    end do
  end subroutine sublattice_equivalences

  !> Lattice sums over each sublattice of the superlattice (legacy FMAT).
  !! Author: Pendry.
  !!
  !! flms(js, lm) is the sum over the lattice points of sublattice js, with
  !! lm = (0,0), (1,-1), (1,1), (2,-2), (2,0), (2,2), ...; it is zero for
  !! odd l + m. The sum is split into one over a single sector of the
  !! lattice and the other ideg - 1 sectors, which are related by symmetry
  !! except for factors that depend on the direction of the lattice vector.
  !! ideg is limited to 6 (arrays anc, ans). Do not use ideg = 1; ideg = 3
  !! is preferable over ideg = 6.
  !!
  !! The arrays flms, v and jjs are explicitly shaped; the driver passes the
  !! contiguous leading part of its larger arrays (its nl equals the leading
  !! extent, see sublattice_equivalences). The work arrays scc and sa are
  !! small (ideg <= 6) and automatic.
  subroutine lattice_sums(flms, v, jjs, nl, nls, dcut, lat, lmax, klm, point)
    !> Number of sublattices.
    integer, intent(in) :: nl
    !> Number of sublattice sums desired (1 or nl).
    integer, intent(in) :: nls
    !> Largest angular momentum.
    integer, intent(in) :: lmax
    !> (2 * lmax + 1) * (2 * lmax + 2) / 2.
    integer, intent(in) :: klm
    !> Lattice sums. Set to zero before the sum.
    complex(sp), intent(out) :: flms(nl, klm)
    !> Lattice; uses superlattice_vector (AR1, AR2) and symmetry_degree
    !! (ideg)
    type(lattice), intent(in) :: lat
    !> Sublattice vectors in polar coordinates, from
    !! sublattice_equivalences.
    real(sp), intent(in) :: v(nl, 2)
    !> Relationship of the sublattices under rotation, from
    !! sublattice_equivalences.
    integer, intent(in) :: jjs(nl, lat%symmetry_degree)
    !> Cutoff distance of the lattice sum.
    real(sp), intent(in) :: dcut
    !> Energy, parallel momentum and imaginary potential; uses all three.
    type(energy_point), intent(in) :: point

    integer :: ideg
    integer :: asst  ! Flag: when ideg <= 2, the sum is run twice
    complex(sp) :: scc(lat%symmetry_degree, 4 * lmax + 1)
    complex(sp) :: sa(2 * lmax * lat%symmetry_degree)
    complex(sp) :: rtab(4), czero, ci, kappa, sc, sd, se, z, acs, acc, rf
    real(sp) :: anc(6), ans(6), r(2), ak(2)
    real(sp) :: ang, d, ag, dcut2, det, add, ant, ast, adr1, adr2
    real(sp) :: b1, b2, bp1, bp2, an1, an2, ar, abc, abb, ad, al
    real(sp) :: ar1(2), ar2(2), e, vpi
    integer :: l2max, lim, liml, li1, li2, i1, i2, j, m, js, jsp
    integer :: mj, mp, mm, jf, jl, jm, km

    ar1 = lat%superlattice_vector(:, 1)
    ar2 = lat%superlattice_vector(:, 2)
    ideg = lat%symmetry_degree
    e = point%energy
    ak = point%k_parallel
    vpi = point%imaginary_potential

    czero = cmplx(0.0_sp, 0.0_sp, sp)
    ci = cmplx(0.0_sp, 1.0_sp, sp)
    kappa = cmplx(2.0_sp * e, -2.0_sp * vpi + 0.000001_sp, sp)
    kappa = sqrt(kappa)
    ag = sqrt(ak(1) * ak(1) + ak(2) * ak(2))

    ! anc, ans and sa are prepared to be used in the sum over symmetrically
    ! related sectors of the lattice
    l2max = lmax + lmax
    lim = l2max * ideg
    liml = l2max + 1
    ang = 2.0_sp * pi / real(ideg, sp)
    d = 1.0_sp
    do j = 1, ideg
      anc(j) = cos(d * ang)
      ans(j) = sin(d * ang)
      d = d + 1.0_sp
    end do
    d = 1.0_sp
    do j = 1, lim
      sa(j) = exp(-ci * d * ang)
      d = d + 1.0_sp
    end do

    flms = czero

    ! The range of summation is limited by dcut
    d = sqrt(ar1(1) * ar1(1) + ar1(2) * ar1(2))
    li1 = int(dcut / d) + 2
    d = sqrt(ar2(1) * ar2(1) + ar2(2) * ar2(2))
    li2 = int(dcut / d) + 2
    dcut2 = dcut * dcut
    det = ar1(1) * ar2(2) - ar1(2) * ar2(1)

    ! One sublattice at a time is treated in the first sector
    do js = 1, nls
      asst = 0
      add = 1.0_sp
      ant = -1.0_sp
      adr1 = v(js, 1) * cos(v(js, 2))
      adr2 = v(js, 1) * sin(v(js, 2))
      ! Shift point adr1,2 by multiples of ar1 and ar2 into the superlattice
      ! unit cell near the origin, defined by the limits (-0.001 to
      ! 0.999) * ar1 and (0.001 to 1.001) * ar2
      b1 = (adr1 * ar2(2) - adr2 * ar2(1)) / det
      b2 = (ar1(1) * adr2 - ar1(2) * adr1) / det
      bp1 = mod(b1, 1.0_sp)
      if (bp1 < -0.001_sp) then
        bp1 = bp1 + 1.0_sp
      else if (bp1 > 0.999_sp) then
        bp1 = bp1 - 1.0_sp
      end if
      bp2 = mod(b2, 1.0_sp)
      if (bp2 < 0.001_sp) then
        bp2 = bp2 + 1.0_sp
      end if
      adr1 = bp1 * ar1(1) + bp2 * ar2(1)
      adr2 = bp1 * ar1(2) + bp2 * ar2(2)
      ast = -1.0_sp
      ! Restart from here once more if ideg <= 2
      do
        an1 = ant
        do i1 = 1, li1  ! First in-plane direction
          an1 = an1 + add
          an2 = ast
          do i2 = 1, li2  ! Second in-plane direction
            an2 = an2 + 1.0_sp
            ! r is the current lattice vector in the sum; ar = |r|;
            ! rtab(1) = -exp(i * phi(r))
            r(1) = an1 * ar1(1) + an2 * ar2(1) + adr1
            r(2) = an1 * ar1(2) + an2 * ar2(2) + adr2
            ar = r(1) * r(1) + r(2) * r(2)
            if (ar > dcut2) cycle  ! Beyond the limit of the lattice sum
            ar = sqrt(ar)
            rtab(1) = -cmplx(r(1) / ar, r(2) / ar, sp)
            abc = 1.0_sp
            abb = 0.0_sp
            if (ag > 1.0e-4_sp) then  ! off-normal beam
              abc = (ak(1) * r(1) + ak(2) * r(2)) / (ag * ar)
              abb = (-ak(2) * r(1) + ak(1) * r(2)) / (ag * ar)
            end if
            sc = ci * ag * ar
            ! scc contains factors in the summation that depend on the
            ! direction of r. Contributions from symmetrically related
            ! sectors can be generated simply and are accumulated for each
            ! sector, indexed by j. The second subscript is ordered
            ! m = -l2max, -l2max + 1, ..., l2max
            do j = 1, ideg
              ad = abc * anc(j) - abb * ans(j)
              sd = exp(sc * ad)
              scc(j, liml) = sd
              mj = 0
              se = rtab(1)
              do m = 1, l2max
                mj = mj + j
                mp = liml + m
                mm = liml - m
                scc(j, mp) = sd * sa(mj) / se
                scc(j, mm) = sd * se / sa(mj)
                se = se * rtab(1)
              end do
            end do
            z = ar * kappa
            acs = sin(z)
            acc = cos(z)
            ! rtab(3), rtab(4) = spherical Hankel functions of the first
            ! kind, l = 0, 1
            rtab(3) = (acs - ci * acc) / z
            rtab(4) = ((acs / z - acc) - ci * (acc / z + acs)) / z
            al = 0.0_sp
            ! The sum over factors independent of the direction of r is
            ! accumulated in flms for each sublattice. The second
            ! subscript orders l and m as (0,0), (1,-1), (1,1), (2,-2),
            ! (2,0), (2,2), ...
            jf = 1
            do jl = 1, liml
              rf = rtab(3) * ci
              jm = l2max + 2 - jl
              do km = 1, jl
                ! The corresponding lattice points in the other sectors
                ! contribute to the appropriate sublattice
                do j = 1, ideg
                  jsp = jjs(js, j)
                  flms(jsp, jf) = flms(jsp, jf) + scc(j, jm) * rf
                end do
                jf = jf + 1
                jm = jm + 2
              end do
              ! Spherical Hankel functions for higher l by recurrence
              acs = (2.0_sp * al + 3.0_sp) * rtab(4) / z - rtab(3)
              rtab(3) = rtab(4)
              rtab(4) = acs
              al = al + 1.0_sp
            end do
          end do
        end do

        ! Special treatment is required if ideg <= 2: two sectors remain to
        ! be summed over
        if (ideg <= 2 .and. asst == 0) then
          asst = 1
          add = -1.0_sp
          ast = -1.0_sp
          if (bp2 > 0.999_sp) ast = -2.0_sp
          ant = 0.0_sp
        else
          exit
        end if
      end do
    end do
  end subroutine lattice_sums

  !> Computes the complex spherical harmonics up to angular momentum lmax for
  !! one direction (legacy SPHRM_MOD). Author: Pendry.
  !!
  !! ylm((l + 1)**2 - l + m) is Y(l, m); lmmax = (lmax + 1)**2. The
  !! argument of the harmonics is given by ct = cos(theta), st = sin(theta)
  !! and cf = exp(i * phi), all complex.
  subroutine spherical_harmonics(lmax, ylm, lmmax, ct, st, cf)
    !> Largest angular momentum.
    integer, intent(in) :: lmax
    !> Number of harmonics, (lmax + 1)**2.
    integer, intent(in) :: lmmax
    !> Spherical harmonics.
    complex(sp), intent(out) :: ylm(lmmax)
    !> cos(theta), sin(theta) and exp(i * phi).
    complex(sp), intent(in) :: ct, st, cf

    ! The auxiliary arrays have fixed sizes: 21 = lmax + 1 and
    ! 441 = (lmax + 1)**2 for lmax = 20.
    integer, parameter :: max_lmax = 20
    real(sp) :: fac1(max_lmax + 1), fac3(max_lmax + 1)
    real(sp) :: fac2((max_lmax + 1)**2)
    complex(sp) :: sf, sa
    real(sp) :: cl, a, b, asg, cm
    integer :: l, ll, lm, lm2, lm3, ln, lo, lp, lq, m

    ! Not reachable from ref-calc: its input check keeps lmax at or below 18
    if (lmax > max_lmax) then
      call stop_with_error('spherical_harmonics: the auxiliary arrays are too ' &
                           // 'small for lmax > 20')
    end if

    lm = 0
    cl = 0.0_sp
    a = 1.0_sp
    b = 1.0_sp
    asg = 1.0_sp
    ll = lmax + 1
    do l = 1, ll
      fac1(l) = asg * sqrt((2.0_sp * cl + 1.0_sp) * a / (4.0_sp * pi * b * b))
      fac3(l) = sqrt(2.0_sp * cl)
      cm = -cl
      ln = l + l - 1
      do m = 1, ln
        lo = lm + m
        fac2(lo) = sqrt((cl + 1.0_sp + cm) * (cl + 1.0_sp - cm) &
                        / ((2.0_sp * cl + 3.0_sp) * (2.0_sp * cl + 1.0_sp)))
        cm = cm + 1.0_sp
      end do
      cl = cl + 1.0_sp
      a = a * 2.0_sp * cl * (2.0_sp * cl - 1.0_sp) / 4.0_sp
      b = b * cl
      asg = -asg
      lm = lm + ln
    end do

    lm = 1
    asg = -1.0_sp
    sf = cf
    sa = cmplx(1.0_sp, 0.0_sp, sp)
    ylm(1) = cmplx(fac1(1), 0.0_sp, sp)
    do l = 1, lmax
      ln = lm + l + l + 1
      ylm(ln) = fac1(l + 1) * sa * sf * st
      ylm(lm + 1) = asg * fac1(l + 1) * sa * st / sf
      ylm(ln - 1) = -fac3(l + 1) * fac1(l + 1) * sa * sf * ct / cf
      ylm(lm + 2) = asg * fac3(l + 1) * fac1(l + 1) * sa * ct * cf / sf
      sa = st * sa
      sf = sf * cf
      asg = -asg
      lm = ln
    end do

    lm = 1
    ll = lmax - 1
    do l = 1, ll
      ln = l + l - 1
      lm2 = lm + ln + 4
      lm3 = lm - ln
      do m = 1, ln
        lo = lm2 + m
        lp = lm3 + m
        lq = lm + m + 1
        ylm(lo) = -(fac2(lp) * ylm(lp) - ct * ylm(lq)) / fac2(lq)
      end do
      lm = lm + l + l + 1
    end do
  end subroutine spherical_harmonics

  !> Computes the spherical harmonics of the propagation directions inside
  !! the crystal of the beams n_offset + 1 ... n_offset + n of pq (legacy
  !! BEAMS_YLM) and stores them in rows 1 ... n of cylm.
  !!
  !! Rows n + 1 ... nt of cylm are left unchanged. xa(ig) returns the
  !! complex perpendicular wave-vector component of beam n_offset + ig inside
  !! the crystal. Used by bravais_layer_matrices, by composite_layer_matrices
  !! (one beam set) and by the driver (all beams, before FINAL1).
  subroutine beams_ylm(cylm, xa, pq, nt, n, n_offset, point, lmax, lmmax, ylm)
    !> Total number of beams at the current energy
    integer, intent(in) :: nt
    !> Number of beams to compute
    integer, intent(in) :: n
    !> Offset of the first of them in the list pq
    integer, intent(in) :: n_offset
    !> Highest angular momentum
    integer, intent(in) :: lmax
    !> (lmax + 1)**2
    integer, intent(in) :: lmmax
    !> Beams (reciprocal lattice vectors) at the current energy
    real(sp), intent(in) :: pq(2, nt)
    !> Energy, parallel momentum of the incident beam, imaginary potential
    type(energy_point), intent(in) :: point
    !> Spherical harmonics of the beams, (beam within its set, (L, M)); rows
    !! 1 ... n are set, the others are kept
    complex(sp), intent(inout) :: cylm(nt, lmmax)
    !> Perpendicular wave-vector component of the beams
    complex(sp), intent(out) :: xa(n)
    !> Work space
    complex(sp), intent(out) :: ylm(lmmax)

    integer :: ig, jg, k
    real(sp) :: bk2, bk3, c, b
    complex(sp) :: ya, ct, st, cf
    complex(sp), parameter :: ru = cmplx(1.0_sp, 0.0_sp, sp)

    ya = cmplx(2.0_sp * point%energy, &
               -2.0_sp * point%imaginary_potential + 0.000001_sp, sp)
    ya = sqrt(ya)
    do ig = 1, n
      jg = ig + n_offset
      bk2 = pq(1, jg) + point%k_parallel(1)
      bk3 = pq(2, jg) + point%k_parallel(2)
      c = bk2 * bk2 + bk3 * bk3
      xa(ig) = cmplx(2.0_sp * point%energy - c, &
                     -2.0_sp * point%imaginary_potential + 0.000001_sp, sp)
      xa(ig) = sqrt(xa(ig))
      b = 0.0_sp
      cf = ru
      if (c > 1.0e-7_sp) then  ! off-normal beam
        b = sqrt(c)
        cf = cmplx(bk2 / b, bk3 / b, sp)
      end if
      ct = xa(ig) / ya
      st = b / ya
      call spherical_harmonics(lmax, ylm, lmmax, ct, st, cf)
      do k = 1, lmmax
        cylm(ig, k) = ylm(k)
      end do
    end do
  end subroutine beams_ylm

  !> Reflection and transmission matrices of a single Bravais-lattice layer
  !! for the beams of one beam set (legacy RSMF_SIMPLE), by Pendry's
  !! formula.
  !!
  !! The matrices of the beam set are stored into ra and ta at the rows and
  !! columns n_offset + 1 ... n_offset + n; entries of other beam sets are
  !! not touched. Thermal vibrations enter through the temperature-dependent
  !! atomic t-matrix caf. Origin shifts and beam symmetry codes are not
  !! supported.
  !!
  !! The beam-space matrices are packed with leading dimension nt, the
  !! number of beams at the current energy, into the leading nt * nt
  !! elements of the (larger) arrays of the caller; the dummies are
  !! explicit-shape for that reason. The rows of amult and cylm (computed by
  !! beams_ylm) are indexed by the beam within its set, without n_offset.
  !! With several beam sets each call therefore overwrites rows 1 ... n of
  !! the previous one. After the call, amult and rows 1 ... n of cylm hold
  !! the values of this beam set; rows n + 1 ... nt of cylm are left
  !! unchanged, as the legacy routine does.
  !!
  !! LAPACK status flags are not checked.
  subroutine bravais_layer_matrices(ra, ta, n, amult, cylm, pq, nt, flms, &
                                    flm, v, nl, n_offset, nll, af, caf, &
                                    lmax, lx, lxi, lmmax, klm, xev, xod, &
                                    lev, lod, ylm, ylme, ylmo, iple, iplo, &
                                    clm, nlm, xevst, xodst, point, tv)
    !> Number of beams in the current beam set
    integer, intent(in) :: n
    !> Total number of beams at the current energy
    integer, intent(in) :: nt
    !> Number of sublattices (rows of flms and v)
    integer, intent(in) :: nl
    !> Offset of the current beam set in the list pq
    integer, intent(in) :: n_offset
    !> Number of sublattices to sum over (1 for an overlayer, nl for a
    !! substrate layer); NLL
    integer, intent(in) :: nll
    !> Highest angular momentum
    integer, intent(in) :: lmax
    !> (lmax + 1)**2
    integer, intent(in) :: lmmax
    !> (2 * lmax + 1) * (2 * lmax + 2) / 2
    integer, intent(in) :: klm
    !> (lmax + 1) * (lmax + 2) / 2
    integer, intent(in) :: lev
    !> lmax * (lmax + 1) / 2
    integer, intent(in) :: lod
    !> Dimension of clm
    integer, intent(in) :: nlm
    !> Permutations of the (L, M) sequence
    integer, intent(in) :: lx(lmmax), lxi(lmmax)
    !> Clebsch-Gordan coefficients
    real(sp), intent(in) :: clm(nlm)
    !> Beams (reciprocal lattice vectors) at the current energy
    real(sp), intent(in) :: pq(2, nt)
    !> Sublattice vectors in polar coordinates (modulus, angle)
    real(sp), intent(in) :: v(nl, 2)
    !> Area of the unit cell
    real(sp), intent(in) :: tv
    !> Energy, parallel momentum of the incident beam, imaginary potential
    type(energy_point), intent(in) :: point
    !> Lattice sums
    complex(sp), intent(in) :: flms(nl, klm)
    !> Atomic t-matrix elements, without and with thermal vibrations
    complex(sp), intent(in) :: af(lmax + 1), caf(lmax + 1)
    !> Reflection and transmission matrices; the entries of the beam set
    !! are set, the others are kept
    complex(sp), intent(inout) :: ra(nt, nt), ta(nt, nt)
    !> Prefactors of the matrix elements of the beams of the set
    complex(sp), intent(out) :: amult(n)
    !> Spherical harmonics of the beams of the set, (beam within its set,
    !! (L, M)); rows 1 ... n are set, the others are kept
    complex(sp), intent(inout) :: cylm(nt, lmmax)
    !> Matrices 1 - X for the L + M even and odd blocks, with the sign
    !! reversed (stored for perturbative LEED)
    complex(sp), intent(out) :: xevst(lev, lev), xodst(lod, lod)
    !> Work spaces
    complex(sp), intent(out) :: flm(klm), ylm(lmmax), ylme(lev), ylmo(lod)
    !> Work spaces, on return the LU factors of 1 - X
    complex(sp), intent(out) :: xev(lev, lev), xod(lod, lod)
    !> Pivot indices of the LU factorizations
    integer, intent(out) :: iple(lev), iplo(lod)

    integer :: ig, jg, jgp, jgps, jlm, jj, js, i, j, k, l, ll, jm, lev1, lm, info
    real(sp) :: bk2, bk3, abr
    complex(sp) :: xa, ya, st, am
    complex(sp) :: rt(2)
    complex(sp), parameter :: ci = cmplx(0.0_sp, 1.0_sp, sp)
    complex(sp), parameter :: cz = cmplx(0.0_sp, 0.0_sp, sp)
    complex(sp), parameter :: ru = cmplx(1.0_sp, 0.0_sp, sp)
    external :: cgetrf, cgetrs

    lm = lmax + 1
    lev1 = lev + 1

    ! Spherical harmonics for all beam directions; amult returns the
    ! perpendicular wave-vector component of each beam
    call beams_ylm(cylm, amult, pq, nt, n, n_offset, point, lmax, lmmax, ylm)
    ya = cmplx(2.0_sp * point%energy, &
               -2.0_sp * point%imaginary_potential + 0.000001_sp, sp)
    ya = sqrt(ya)
    do ig = 1, n
      ! Prefactor for reflection and transmission matrix elements
      amult(ig) = -8.0_sp * pi * pi * ci / (ya * tv * amult(ig))
    end do

    ! Sum over sublattices
    do k = 1, klm
      flm(k) = cz
    end do
    bk2 = pq(1, 1 + n_offset)
    bk3 = pq(2, 1 + n_offset)
    do js = 1, nll
      abr = v(js, 1) * (bk2 * cos(v(js, 2)) + bk3 * sin(v(js, 2)))
      xa = exp(abr * ci)
      do i = 1, klm
        flm(i) = flm(i) + flms(js, i) * xa
      end do
    end do

    ! Intralayer multiple-scattering matrix X, split into L + M even (xev)
    ! and L + M odd (xod)
    call intralayer_matrix(flm, xev, xod, lev, lod, caf, lm, lxi, lmmax, klm, &
                           clm, nlm)
    do j = 1, lev
      xev(j, j) = xev(j, j) - ru
    end do
    do j = 1, lod
      xod(j, j) = xod(j, j) - ru
    end do

    ! Store 1 - X (with reversed sign) for perturbative LEED
    do j = 1, lev
      do i = 1, lev
        xevst(i, j) = -xev(i, j)
      end do
    end do
    do j = 1, lod
      do i = 1, lod
        xodst(i, j) = -xod(i, j)
      end do
    end do

    ! LU factorization of 1 - X; the sign is taken care of in amult
    call cgetrf(lev, lev, xev, lev, iple, info)
    call cgetrf(lod, lod, xod, lod, iplo, info)

    ! Loop over scattered beams
    do jgp = 1, n
      jgps = jgp + n_offset

      ! Quantities to be multiplied from the left by the inverse of 1 - X
      jlm = 1
      st = ru
      do l = 1, lm
        ll = l + l - 1
        do jm = 1, ll
          ylm(jlm) = cylm(jgp, jlm) * af(l) / st
          jlm = jlm + 1
        end do
        st = st * ci
      end do
      ! Change ordering of (L, M) pairs
      do jlm = 1, lev
        jj = lx(jlm)
        ylme(jlm) = ylm(jj)
      end do
      do jlm = lev1, lmmax
        jj = lx(jlm)
        ylmo(jlm - lev) = ylm(jj)
      end do

      ! Complete the inversion and multiply by back-substitution
      call cgetrs('N', lev, 1, xev, lev, iple, ylme, lev, info)
      call cgetrs('N', lod, 1, xod, lod, iplo, ylmo, lod, info)
      ! Reorder the (L, M) pairs
      do jlm = 1, lev
        jj = lx(jlm)
        ylm(jj) = ylme(jlm)
      end do
      do jlm = lev1, lmmax
        jj = lx(jlm)
        ylm(jj) = ylmo(jlm - lev)
      end do

      am = amult(jgp)
      ! Loop over incident beams
      do jg = 1, n
        call matrix_elements(jg, lm, ylm, cylm, lmmax, nt, rt)
        ra(jgps, jg + n_offset) = rt(1) * am
        ta(jgps, jg + n_offset) = rt(2) * am
      end do
      ! Add the unscattered beam to the transmission
      ta(jgps, jgps) = ta(jgps, jgps) + ru
    end do
  end subroutine bravais_layer_matrices

  !> Reflection and transmission matrix elements of a Bravais-lattice layer
  !! for one incident beam, from the (1 - X)**(-1) Y T vector and the
  !! spherical harmonics (legacy MFOLD_SIMPLE, simplified: no Debye-Waller
  !! factor, origin shifts, or symmetry folding). The sum runs over (L, M).
  subroutine matrix_elements(jg, lm, ylm, cylm, lmmax, nt, rt)
    !> Index of the current incident beam within its beam set
    integer, intent(in) :: jg
    !> lmax + 1
    integer, intent(in) :: lm
    !> Number of (L, M) pairs, (lmax + 1)**2
    integer, intent(in) :: lmmax
    !> Total number of beams
    integer, intent(in) :: nt
    !> (1 - X)**(-1) Y T
    complex(sp), intent(in) :: ylm(lmmax)
    !> Spherical harmonics of the beams
    complex(sp), intent(in) :: cylm(nt, lmmax)
    !> Reflection (1) and transmission (2) matrix element
    complex(sp), intent(out) :: rt(2)

    integer :: l, ll, l1, m, mm, jlm
    complex(sp) :: r, t, st, sm, sl, cy, ctr, ctt
    complex(sp), parameter :: ci = cmplx(0.0_sp, 1.0_sp, sp)
    complex(sp), parameter :: cz = cmplx(0.0_sp, 0.0_sp, sp)
    complex(sp), parameter :: ru = cmplx(1.0_sp, 0.0_sp, sp)

    r = cz
    t = cz
    st = ru
    sl = ru
    jlm = 1
    l1 = 1
    do l = 1, lm
      sm = sl
      ll = l + l - 1
      l1 = l1 + ll
      do m = 1, ll
        mm = l1 - m
        cy = cylm(jg, mm)
        ctt = ylm(jlm) * st
        ctr = ctt * sl
        ctt = ctt * sm
        r = r + ctr * cy
        t = t + ctt * cy
        jlm = jlm + 1
        sm = -sm
      end do
      sl = -sl
      st = st * ci
    end do

    rt(1) = r
    rt(2) = t
  end subroutine matrix_elements

  !> Intralayer multiple-scattering matrix X of a single Bravais-lattice
  !! layer (legacy XM).
  !!
  !! The (L, M) sequence is the symmetrized one and X is split into two
  !! blocks, L + M even (xev) and L + M odd (xod), by block diagonalization.
  !! The Clebsch-Gordan coefficients are used in the order in which they
  !! are stored.
  subroutine intralayer_matrix(flm, xev, xod, lev, lod, caf, lm, lxi, lmmax, &
                               klm, clm, nlm)
    !> (lmax + 1) * (lmax + 2) / 2
    integer, intent(in) :: lev
    !> lmax * (lmax + 1) / 2
    integer, intent(in) :: lod
    !> lmax + 1
    integer, intent(in) :: lm
    !> (lmax + 1)**2
    integer, intent(in) :: lmmax
    !> (2 * lmax + 1) * (2 * lmax + 2) / 2
    integer, intent(in) :: klm
    !> Dimension of clm
    integer, intent(in) :: nlm
    !> Permutation of the (L, M) sequence
    integer, intent(in) :: lxi(lmmax)
    !> Clebsch-Gordan coefficients
    real(sp), intent(in) :: clm(nlm)
    !> Lattice sums of the layer
    complex(sp), intent(in) :: flm(klm)
    !> Temperature-dependent atomic t-matrix elements
    complex(sp), intent(in) :: caf(lm)
    !> The two blocks of X
    complex(sp), intent(out) :: xev(lev, lev), xod(lod, lod)

    integer :: lmax, l2max, jset, mm, n, j, k, l, m, jl, lpp, mpp, mpa, lpa
    integer :: mp1, lp1, jlm, jx, kx
    complex(sp) :: acc
    complex(sp), parameter :: czero = cmplx(0.0_sp, 0.0_sp, sp)

    lmax = lm - 1
    xev = czero
    xod = czero
    l2max = lmax + lmax
    n = 1
    ! First the odd block (jset = 1, built in xev and copied to xod), then
    ! the even one. The index n runs through clm across both.
    do jset = 1, 0, -1
      mm = 0
      if (jset == 1) mm = lev
      j = 1
      do l = jset, lmax
        jl = l + 1
        do m = -l + jset, l, 2
          k = 1
          do lpp = jset, lmax
            do mpp = -lpp + jset, lpp, 2
              mpa = abs(mpp - m)
              lpa = abs(lpp - l)
              if (lpa > mpa) mpa = lpa
              mp1 = mpp - m + l2max + 1
              lp1 = l + lpp + 1
              acc = czero
              do
                jlm = (lp1 * lp1 + mp1 - l2max) / 2
                acc = acc + clm(n) * flm(jlm)
                n = n + 1
                lp1 = lp1 - 2
                if (lp1 - 1 - mpa < 0) exit
              end do
              jx = lxi(j + mm)
              kx = lxi(k + mm)
              xev(jx, kx) = acc * caf(jl)
              k = k + 1
            end do
          end do
          j = j + 1
        end do
      end do
      if (jset == 1) xod = xev(1:lod, 1:lod)
    end do
  end subroutine intralayer_matrix
end module bravais_layer
