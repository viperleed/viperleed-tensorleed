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
  !!
  !! The legacy copy in lib.tleed.f is still used by RTINV_SIMPLE, which
  !! cannot use this module; this routine is a temporary duplicate of it.
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
end module bravais_layer
