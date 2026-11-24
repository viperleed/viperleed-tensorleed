!> \file composite_layer.f90
!! \brief Composite layers (several subplanes): ordering of the subplanes
!! and the propagators between them.
!!
!! Each routine keeps the precision of every quantity of the legacy routine
!! it follows, and its results are bit-for-bit those of the legacy code
!! (tests/tier1/test_refcalc_routines.py).

!> Subplanes of a composite layer, the propagators between them and the
!! multiple scattering within them.
module composite_layer
  use constants, only: pi
  use kinds, only: sp, dp
  use leed_types, only: lattice, energy_point
  implicit none

  private

  public :: sort_subplanes, interplanar_propagators, subplane_scattering

  !> The literal for pi in the estimate of the number of reciprocal-lattice
  !! points of interplanar_propagators; it differs from constants::pi in the
  !! last digit and is kept as in the legacy code.
  real(sp), parameter :: pi_estimate = 3.1415926_sp

contains

  !> Orders the subplanes of a composite layer by increasing position along
  !! +x (legacy SRTLAY); that is from the vacuum side towards the solid.
  !! Author: Pendry.
  !!
  !! The chemical-element assignment is reordered along. Subplanes at equal
  !! position are not permuted. The topmost subplane becomes the origin of
  !! the layer. Also finds the interplanar vectors and which of the
  !! propagators between subplanes are equal to an earlier one.
  !!
  !! The arrays are explicitly shaped; the driver passes its larger arrays
  !! whole, so that the leading part of their storage holds the packed
  !! (nlay, ...) arrays.
  subroutine sort_subplanes(pos, poss, lps, lpss, mgh, nlay, drl, nlay2, &
                            nugh, ngeq, nord)
    !> Number of subplanes.
    integer, intent(in) :: nlay
    !> nlay * (nlay - 1) / 2.
    integer, intent(in) :: nlay2
    !> Subplane positions (atomic positions in the unit cell).
    real(sp), intent(in) :: pos(nlay, 3)
    !> Reordered pos; the first x is zero.
    real(sp), intent(out) :: poss(nlay, 3)
    !> Chemical element assignment of each subplane.
    integer, intent(in) :: lps(nlay)
    !> Reordered lps.
    integer, intent(out) :: lpss(nlay)
    !> Organization of the propagators: mgh(i, j) = k means that the
    !! propagator from subplane j to subplane i is the k-th block of the
    !! columnar matrix gh (nlay2 blocks for i < j, then nlay2 for i > j).
    integer, intent(out) :: mgh(nlay, nlay)
    !> Interplanar vectors.
    real(sp), intent(out) :: drl(nlay2, 3)
    !> nugh(k) = 1: propagator k must be computed; 0: it is equal to an
    !! earlier one (see ngeq).
    integer, intent(out) :: nugh(nlay2)
    !> ngeq(k) = m: propagator k can be copied from the new propagator m;
    !! 0: no copying.
    integer, intent(out) :: ngeq(nlay2)
    !> Order of the subplanes, to trace the reordering back.
    integer, intent(out) :: nord(nlay)

    integer :: i, j, k, m, in, km, lpsa, norda
    real(sp) :: pm, posa(3)

    do i = 1, nlay
      lpss(i) = lps(i)
      nord(i) = i
      poss(i, :) = pos(i, :)
    end do

    ! Analyse the order of the subplane positions along the x axis and
    ! reorder them in ascending position
    do i = 1, nlay - 1
      km = i
      pm = poss(i, 1)
      do k = i + 1, nlay
        if (poss(k, 1) >= pm) cycle
        pm = poss(k, 1)
        km = k
      end do
      if (km == i) cycle
      posa = poss(km, :)
      lpsa = lpss(km)
      norda = nord(km)
      do k = km, i + 1, -1
        poss(k, :) = poss(k - 1, :)
        nord(k) = nord(k - 1)
        lpss(k) = lpss(k - 1)
      end do
      poss(i, :) = posa
      lpss(i) = lpsa
      nord(i) = norda
    end do

    ! Set the origin of the composite layer to the topmost subplane
    do i = 2, nlay
      poss(i, 1) = poss(i, 1) - poss(1, 1)
    end do
    poss(1, 1) = 0.0_sp

    ! Interplanar vectors, organization of gh, and equal interplanar
    ! vectors
    in = 1
    do i = 1, nlay
      do j = i, nlay
        mgh(i, j) = 0
        mgh(j, i) = 0
        if (i == j) cycle
        drl(in, :) = poss(i, :) - poss(j, :)
        mgh(i, j) = in
        mgh(j, i) = in + nlay2
        nugh(in) = 1
        ngeq(in) = 0
        do m = 1, in - 1
          posa = abs(drl(in, :) - drl(m, :))
          if ((posa(1) + posa(2) + posa(3)) - 0.001_sp <= 0.0_sp) then
            ngeq(in) = m
            nugh(in) = 0
            exit
          end if
        end do
        in = in + 1
      end do
    end do
  end subroutine sort_subplanes

  !> Interplanar propagators between the subplanes of a composite layer, in
  !! (l,m) space (legacy GHMAT). Author: Pendry.
  !!
  !! Depending on the interplanar spacing, the propagators are computed by
  !! a reciprocal-space summation or by a direct-space summation (in
  !! direct_propagator_sum). The (l,m) elements of each propagator are
  !! ordered (0,0), (2,-2), (2,0), (2,2), (4,-4), (4,-2), ..., (1,-1),
  !! (1,1), (3,-3), (3,1), ..., (1,0), (3,-2), (3,0), (3,2), (5,-4), ...,
  !! (2,-1), (2,1), (4,-3), (4,-1), (4,1), ... The propagators gh(i, j)
  !! (i, j subplane indices) are stacked on top of each other in the
  !! columnar matrix gh, as organized by sort_subplanes. The propagators
  !! that are equal to an earlier one are copied from it.
  !!
  !! nugh is modified: it is zero on exit for every propagator that was
  !! computed. The work arrays work and s are overwritten. The array
  !! work(l2m, l2m) is used as the table of spherical harmonics by the
  !! direct sum (whole) and by the reciprocal sum (leading lm by lm part).
  !! The array pq is the parallel momentum of the first beam of the
  !! current beam set.
  subroutine interplanar_propagators(gh, lmg, lmmax, nugh, ngeq, nlay2, tst, &
                                     test, work, l2m, lm, s, lms, drl, tv, &
                                     lxm, dcut, caa, ncaa, lay, hghd, pq, lat, &
                                     point)
    !> 2 * nlay2 * lmmax.
    integer, intent(in) :: lmg
    !> (lmax + 1)**2.
    integer, intent(in) :: lmmax
    !> nlay * (nlay - 1) / 2.
    integer, intent(in) :: nlay2
    !> 2 * lmax + 1, or larger (as legacy).
    integer, intent(in) :: l2m
    !> lmax + 1.
    integer, intent(in) :: lm
    !> l2m * l2m.
    integer, intent(in) :: lms
    !> Number of Clebsch-Gordan coefficients.
    integer, intent(in) :: ncaa
    !> All interplanar propagators, stacked vertically.
    complex(sp), intent(inout) :: gh(lmg, lmmax)
    !> Propagators that must still be computed (see sort_subplanes).
    integer, intent(inout) :: nugh(nlay2)
    !> Propagators that can be copied from newly computed ones.
    integer, intent(in) :: ngeq(nlay2)
    !> Convergence criterion of the reciprocal-lattice sum.
    real(sp), intent(in) :: tst
    !> Cutoff of the reciprocal-lattice sum of each propagator; only the
    !! entries of the propagators that are computed are set.
    real(sp), intent(inout) :: test(nlay2)
    !> Work array: spherical harmonics.
    complex(sp), intent(out) :: work(l2m, l2m)
    !> Work array: lattice sum.
    complex(sp), intent(out) :: s(lms)
    !> Interplanar vectors, from sort_subplanes.
    real(sp), intent(in) :: drl(nlay2, 3)
    !> Area of the unit cell of each subplane.
    real(sp), intent(in) :: tv
    !> Permutation of the (l,m) sequence.
    integer, intent(in) :: lxm(lmmax)
    !> Cutoff radius of the direct-lattice sum.
    real(sp), intent(in) :: dcut
    !> Clebsch-Gordan coefficients.
    real(sp), intent(in) :: caa(ncaa)
    !> 1 for an overlayer, otherwise a substrate layer.
    integer, intent(in) :: lay
    !> Work array of direct_propagator_sum.
    complex(sp), intent(inout) :: hghd(l2m)
    !> Parallel momentum of the first beam of the current beam set.
    real(sp), intent(in) :: pq(2)
    !> Lattices; uses substrate_vector, superlattice_vector and
    !! superlattice_reciprocal.
    type(lattice), intent(in) :: lat
    !> Energy, parallel momentum and imaginary potential.
    type(energy_point), intent(in) :: point

    integer :: i, j, k, iz, ii, in, iin, inn, is, it, kk, nlylm
    integer :: ncount, jj1, jj2, j1, j2, ng1, ng2, nog
    integer :: ip, jp, nii, nni, njm, nnm, im1, imp, jm2, jmp
    integer :: l1, l2, m1, m2
    real(sp) :: e, vpi, ak2, ak3
    real(sp) :: rxr1(2), rxr2(2), atv, tsts, fact1, akp2, gx, gy, gkx, gky, gk2
    real(sp) :: cfy, fy, phase
    complex(sp) :: ci, cz, kprg, k0, z, t1, t1p, t2, t2p, t3, bs, cs, cfac

    e = point%energy
    vpi = point%imaginary_potential
    ak2 = point%k_parallel(1)
    ak3 = point%k_parallel(2)

    cz = cmplx(0.0_sp, 0.0_sp, sp)
    ci = cmplx(0.0_sp, 1.0_sp, sp)
    if (lay == 1) then  ! overlayer
      rxr1 = lat%superlattice_reciprocal(:, 1)
      rxr2 = lat%superlattice_reciprocal(:, 2)
    else  ! substrate layer
      atv = 2.0_sp * pi / tv
      rxr1(1) = lat%substrate_vector(2, 2) * atv
      rxr1(2) = -lat%substrate_vector(1, 2) * atv
      rxr2(1) = -lat%substrate_vector(2, 1) * atv
      rxr2(2) = lat%substrate_vector(1, 1) * atv
    end if
    k0 = sqrt(cmplx(2.0_sp * e, -2.0_sp * vpi + 0.000001_sp, sp))
    cfac = -16.0_sp * pi * pi * ci / tv
    nlylm = nlay2 * lmmax

    ! Initialize the propagators to be computed
    do i = 1, nlay2
      if (nugh(i) == 0) cycle
      is = (i - 1) * lmmax
      it = is + nlylm
      do j = 1, lmmax
        do k = 1, lmmax
          gh(is + k, j) = cz
          gh(it + k, j) = cz
        end do
      end do
    end do

    ! Estimated number of points in the direct-lattice sum
    tsts = dcut * dcut * pi / tv
    do iz = 1, nlay2
      if (nugh(iz) == 0) cycle
      if (abs(drl(iz, 1)) > 0.001_sp) then
        ! Estimated number of points in the reciprocal-lattice sum
        fact1 = log(tst) / drl(iz, 1)
        akp2 = (2.0_sp * e + fact1 * fact1) * tv / (4.0_sp * pi_estimate)
        ! Skip the direct sum if the reciprocal one is faster (but the
        ! number of reciprocal-lattice points must be restricted due to a
        ! convergence problem)
        if ((tsts >= 2.0_sp * akp2) .and. (akp2 < 80.0_sp)) cycle
      end if
      ! gh(i, j) and gh(j, i) with a direct-lattice sum, for the two
      ! propagation directions
      do k = 1, 2
        call direct_propagator_sum(iz, k, gh, lmg, lmmax, s, lms, work, l2m, &
                                   drl, nlay2, k0, dcut, caa, ncaa, lxm, lay, &
                                   pq, hghd, lat, point)
      end do
      nugh(iz) = 0
    end do

    ! Cutoff of the reciprocal-lattice sum for each remaining propagator.
    ! The convergence of the planar sum is not related to the convergence
    ! of the layer-doubling scheme controlled by tst.
    tsts = 0.0_sp
    do i = 1, nlay2
      if ((nugh(i) /= 0) .and. (abs(drl(i, 1)) > 0.001_sp)) then
        test(i) = 8.0_sp * abs(log(0.002_sp) / drl(i, 1))
        tsts = max(test(i), tsts)
      end if
    end do

    if (tsts > 0.00001_sp) then
      ! Two one-dimensional summation loops in one quadrant of reciprocal
      ! space
      ncount = 0
      jj1 = 0
      rows: do
        jj1 = jj1 + 1
        jj2 = 0
        columns: do
          jj2 = jj2 + 1
          nog = 0
          j1 = jj1 - 1
          j2 = jj2 - 1
          quadrants: do kk = 1, 4
            select case (kk)
            case (1)
              ng1 = j1
              ng2 = j2
            case (2)
              if (j1 == 0 .and. j2 == 0) exit quadrants
              if (j2 == 0) cycle quadrants
              ng1 = j1
              ng2 = -j2
            case (3)
              if (j1 == 0) cycle quadrants
              ng1 = -j1
              ng2 = j2
            case default
              if (j1 == 0 .or. j2 == 0) cycle quadrants
              ng1 = -j1
              ng2 = -j2
            end select

            ! Current reciprocal-lattice point
            gx = real(ng1, sp) * rxr1(1) + real(ng2, sp) * rxr2(1)
            gy = real(ng1, sp) * rxr1(2) + real(ng2, sp) * rxr2(2)
            gkx = gx + ak2 + pq(1)
            gky = gy + ak3 + pq(2)
            gk2 = gkx * gkx + gky * gky

            ! Test for the cutoff
            kprg = sqrt(cmplx(2.0_sp * e - gk2, -2.0_sp * vpi + 0.000001_sp, sp))
            akp2 = aimag(kprg)
            if (akp2 > tsts) cycle quadrants
            nog = 1
            z = kprg / k0
            fy = 0.0_sp
            if (gk2 > 1.0e-8_sp) then
              cfy = gkx / sqrt(gk2)
              if (abs(abs(cfy) - 1.0_sp) <= 1.0e-6_sp) then
                if (cfy < 0.0_sp) fy = pi
              else
                fy = acos(cfy)
              end if
              if (gky < 0.0_sp) fy = -fy
            end if

            ! Spherical harmonics
            call harmonics_table(lm, z, fy, work(1:lm, 1:lm))

            ! Loop over the interplanar vectors
            do iz = 1, nlay2
              ! Skip if the propagator is not needed or this contribution
              ! is outside the cutoff
              if (nugh(iz) == 0) cycle
              if (akp2 > test(iz)) cycle
              t2 = gkx * drl(iz, 2) + gky * drl(iz, 3)
              t2p = -t2
              t3 = kprg * abs(drl(iz, 1))
              t1 = (exp(ci * (t2 + t3)) / kprg) * cfac
              t1p = (exp(ci * (t2p + t3)) / kprg) * cfac

              ! A limited set of matrix elements
              do i = 1, lmmax
                ip = lxm(i)
                nii = ip + (iz - 1) * lmmax
                nni = nii + nlylm
                l1 = int(sqrt(real(i - 1, sp) + 0.00001_sp))
                m1 = i - l1 - l1 * l1 - 1
                if (m1 > 0) then
                  bs = work(m1, l1 + 1)
                else
                  bs = (-1.0_sp)**(mod(m1, 2)) * work(l1 + 1, -m1 + 1)
                end if
                do j = 1, lmmax
                  jp = lxm(j)
                  l2 = int(sqrt(real(j - 1, sp) + 0.00001_sp))
                  m2 = j - l2 - l2 * l2 - 1
                  if (i == j) then
                    if (m1 > 0) cycle
                  else if (l1 < l2 .or. (l1 == l2 .and. abs(m1) < abs(m2))) then
                    cycle
                  end if
                  if (m2 >= 0) then
                    cs = work(l2 + 1, m2 + 1)
                  else
                    cs = (-1.0_sp)**(mod(m2, 2)) * work(-m2, l2 + 1)
                  end if
                  ! gh(i, j) and gh(j, i) are produced together
                  gh(nni, jp) = gh(nni, jp) + t1p * bs * cs
                  gh(nii, jp) = gh(nii, jp) &
                      + (-1.0_sp)**(mod(l1 + m1 + l2 + m2, 2)) * t1 * bs * cs
                end do
              end do
            end do
          end do quadrants

          if (nog == 1) cycle columns
          if (jj2 /= 1) cycle rows
          ncount = ncount + 1
          if (ncount > 3) exit rows
        end do columns
      end do rows

      ! Fill in the missing matrix elements from the elements just
      ! produced, using symmetry relations
      do iz = 1, nlay2
        if (nugh(iz) == 0) cycle
        do i = 1, lmmax
          ip = lxm(i)
          nii = ip + (iz - 1) * lmmax
          nni = nii + nlylm
          l1 = int(sqrt(real(i - 1, sp) + 0.00001_sp))
          m1 = i - l1 - l1 * l1 - 1
          im1 = i - 2 * m1
          imp = lxm(im1)
          do j = 1, lmmax
            jp = lxm(j)
            l2 = int(sqrt(real(j - 1, sp) + 0.00001_sp))
            m2 = j - l2 - l2 * l2 - 1
            jm2 = j - 2 * m2
            jmp = lxm(jm2)
            njm = jmp + (iz - 1) * lmmax
            nnm = njm + nlylm
            phase = (-1.0_sp)**mod(m1 + m2, 2)
            if (i == j) then
              if (m1 <= 0) cycle
            else if (l1 >= l2 .and. &
                     (l1 /= l2 .or. abs(m1) >= abs(m2))) then
              cycle
            end if
            gh(nni, jp) = phase * gh(nnm, imp)
            gh(nii, jp) = phase * gh(njm, imp)
          end do
        end do
        nugh(iz) = 0
      end do
    end if

    ! Copy the new propagators onto the others that are identical
    do iz = 1, nlay2
      if (ngeq(iz) == 0) cycle
      ii = (iz - 1) * lmmax
      iin = ii + nlylm
      in = (ngeq(iz) - 1) * lmmax
      inn = in + nlylm
      do j = 1, lmmax
        do i = 1, lmmax
          gh(ii + i, j) = gh(in + i, j)
          gh(iin + i, j) = gh(inn + i, j)
        end do
      end do
    end do
  end subroutine interplanar_propagators

  !> Intra-subplane multiple scattering in (l,m) space for one Bravais
  !! lattice layer, for each of several chemical elements (legacy TAUMAT).
  !! Author: Pendry.
  !!
  !! The (l,m) sequence is the symmetrized one. tau is split into
  !! block-diagonal parts for even and odd l+m, of dimensions lev and lod;
  !! in tau the odd part is left-justified and reduced in width for storage
  !! efficiency. The tau matrices of the different chemical elements are
  !! stacked on top of each other: first the one built from the atomic
  !! t-matrix elements tsf(1, l), then those from tsf(2, l), and so on. Only
  !! the lattice sum of the first sublattice is used.
  !!
  !! The work arrays x, xh and flm are overwritten.
  subroutine subplane_scattering(tau, lmt, ntau, x, lev, lev2, lod, tsf, &
                                 lmmax, flms, nl, klm, lm, clm, nlm, lxi, &
                                 pq, nnsub, lmax1, xh, flm, lay, lat, point)
    !> ntau * lmmax.
    integer, intent(in) :: lmt
    !> Number of chemical elements.
    integer, intent(in) :: ntau
    !> (lmax + 1) * (lmax + 2) / 2: dimension of the even-l+m block.
    integer, intent(in) :: lev
    !> 2 * lev.
    integer, intent(in) :: lev2
    !> lmmax - lev: dimension of the odd-l+m block.
    integer, intent(in) :: lod
    integer, intent(in) :: nnsub, lmax1
    !> (lmax + 1)**2.
    integer, intent(in) :: lmmax
    !> Number of sublattices of the lattice sums.
    integer, intent(in) :: nl
    integer, intent(in) :: klm
    !> lmax + 1.
    integer, intent(in) :: lm
    !> Dimension of clm.
    integer, intent(in) :: nlm
    !> Output: tau matrices.
    complex(sp), intent(inout) :: tau(lmt, lev)
    !> Work array.
    complex(sp), intent(out) :: x(lev, lev2)
    !> Atomic t-matrix elements.
    complex(sp), intent(in) :: tsf(nnsub, lmax1)
    !> Lattice sums, from the layer-matrix routines.
    complex(sp), intent(in) :: flms(nl, klm)
    !> Clebsch-Gordan coefficients.
    real(sp), intent(in) :: clm(nlm)
    !> Permutation of the (l,m) sequence.
    integer, intent(in) :: lxi(lmmax)
    !> Parallel momentum of the first beam of the current beam set.
    real(sp), intent(in) :: pq(2)
    !> Work array of the linear solver.
    complex(sp), intent(out) :: xh(lev)
    !> Work array: total lattice sum.
    complex(sp), intent(out) :: flm(klm)
    !> 1 for an overlayer, otherwise a substrate layer.
    integer, intent(in) :: lay
    !> Lattices; uses substrate_vector and repetitions.
    type(lattice), intent(in) :: lat
    !> Energy and imaginary potential.
    type(energy_point), intent(in) :: point

    integer :: i, j, k, it, il, is, l, ld, ld1, ll, js, clm_index, lmax
    real(sp) :: bk2, bk3, s1, s2, adr1, adr2, abr
    complex(sp) :: ak, cz, ci, xa

    lmax = lm - 1
    cz = cmplx(0.0_sp, 0.0_sp, sp)
    ci = cmplx(0.0_sp, 1.0_sp, sp)
    ak = -0.5_sp / sqrt(cmplx(2.0_sp * point%energy, &
                              -2.0_sp * point%imaginary_potential &
                              + 0.000001_sp, sp))
    flm = cz
    bk2 = pq(1)
    bk3 = pq(2)

    ! Total lattice sum: the sum over the sublattices with phase factors
    ! (only the first term for an overlayer)
    js = 1
    s1 = 0.0_sp
    sublattices: do j = 1, lat%repetitions(1)
      s2 = 0.0_sp
      do k = 1, lat%repetitions(2)
        adr1 = s1 * lat%substrate_vector(1, 1) + s2 * lat%substrate_vector(1, 2)
        adr2 = s1 * lat%substrate_vector(2, 1) + s2 * lat%substrate_vector(2, 2)
        abr = adr1 * bk2 + adr2 * bk3
        xa = exp(abr * ci)
        do i = 1, klm
          flm(i) = flm(i) + flms(js, i) * xa
        end do
        if (lay == 1) exit sublattices
        js = js + 1
        s2 = s2 + 1.0_sp
      end do
      s1 = s1 + 1.0_sp
    end do sublattices

    do it = 1, ntau
      do il = 1, 2
        x = cz
        ll = lod
        if (il == 2) ll = lev
        ! Matrix 1 - X for odd l+m (il = 1), then for even l+m (il = 2).
        ! clm_index runs on from the first call to the second
        call multiple_scattering_matrix(il, flm, x, lev, ll, tsf, it, lm, &
                                        lxi, lmmax, klm, clm, nlm, clm_index, &
                                        nnsub, lmax1)

        ! Quantities into which the inverse of 1 - X is multiplied: the
        ! diagonal blocks, first those of l = 0 (1), 2 (3), ... then those
        ! of l = 1 (2), 3 (4), ... for even (odd) l+m
        ld1 = 0
        do is = 0, 1
          if (il == 2) then
            do l = is, lmax, 2
              ld = ld1 + 1
              ld1 = ld + l
              do i = ld, ld1
                x(i, i + lev) = ak * tsf(it, l + 1)
              end do
            end do
          else
            do l = 1 + is, lmax, 2
              ld = ld1 + 1
              ld1 = ld + l - 1
              do i = ld, ld1
                x(i, i + lod) = ak * tsf(it, l + 1)
              end do
            end do
          end if
        end do

        ! Inversion and multiplication
        call solve_pivoted(x, lev, ll, ll, xh)

        ld = (it - 1) * lmmax
        if (il == 1) ld = ld + lev
        do j = 1, ll
          do i = 1, ll
            tau(ld + i, j) = x(i, j + ll)
          end do
        end do
      end do
    end do
  end subroutine subplane_scattering

  !> Intra-layer multiple-scattering matrix 1 - X for subplane_scattering
  !! (legacy XMT), in the Tong convention, from input in the Pendry
  !! convention. Author: Pendry.
  !!
  !! It must be called twice, first for odd l+m (il = 1), then for even l+m
  !! (il = 2); clm_index, the running index of clm, is set by the first call
  !! and must not be changed in between.
  !! x is the output matrix 1 - X, block-diagonalized, with the second block
  !! left-adjusted.
  subroutine multiple_scattering_matrix(il, flm, x, lev, ll, tsf, it, lm, lxi, &
                                        lmmax, klm, clm, nlm, clm_index, &
                                        nnsub, lmax1)
    !> 1 for odd l+m, 2 for even l+m.
    integer, intent(in) :: il
    integer, intent(in) :: lev, lmmax, klm, nlm, nnsub, lmax1
    !> lev for even l+m, otherwise lod.
    integer, intent(in) :: ll
    !> Index of the current chemical element.
    integer, intent(in) :: it
    !> lmax + 1.
    integer, intent(in) :: lm
    !> Total lattice sum.
    complex(sp), intent(in) :: flm(klm)
    !> The leading ll columns of the work array of subplane_scattering.
    complex(sp), intent(inout) :: x(lev, ll)
    complex(sp), intent(in) :: tsf(nnsub, lmax1)
    integer, intent(in) :: lxi(lmmax)
    real(sp), intent(in) :: clm(nlm)
    !> Running index of clm: set to 1 if il = 1, advanced by every call.
    integer, intent(inout) :: clm_index

    integer :: lmax, l2max, jset, mm, j, k, l, m, lpp, mpp, mpa, lpa
    integer :: mp1, lp1, jlm, jx, kx
    complex(sp) :: acc, ci, st, su, ru

    lmax = lm - 1
    ru = cmplx(1.0_sp, 0.0_sp, sp)
    ci = cmplx(0.0_sp, 1.0_sp, sp)
    l2max = lmax + lmax
    if (il < 2) then  ! odd l+m only
      jset = 1
      mm = lev
      clm_index = 1
    else  ! even l+m only
      jset = 0
      mm = 0
    end if

    j = 1
    l = jset
    l_loop: do
      m = -l + jset
      st = (-ci)**mod(l, 4)
      m_loop: do
        k = 1
        lpp = jset
        lpp_loop: do
          mpp = -lpp + jset
          su = ((-ci)**mod(lpp, 4)) / st
          mpp_loop: do
            mpa = abs(mpp - m)
            lpa = abs(lpp - l)
            if (lpa > mpa) mpa = lpa
            mp1 = mpp - m + l2max + 1
            lp1 = l + lpp + 1
            acc = cmplx(0.0_sp, 0.0_sp, sp)
            do
              jlm = (lp1 * lp1 + mp1 - l2max) / 2
              acc = acc + clm(clm_index) * flm(jlm)
              clm_index = clm_index + 1
              lp1 = lp1 - 2
              if (lp1 - 1 - mpa < 0) exit
            end do
            jx = lxi(j + mm)
            kx = lxi(k + mm)
            x(kx, jx) = -acc * tsf(it, lpp + 1) * su
            if (j == k) x(kx, jx) = x(kx, jx) + ru
            k = k + 1
            mpp = mpp + 2
            if (lpp - mpp < 0) exit mpp_loop
          end do mpp_loop
          lpp = lpp + 1
          if (lmax - lpp < 0) exit lpp_loop
        end do lpp_loop
        j = j + 1
        m = m + 2
        if (l - m < 0) exit m_loop
      end do m_loop
      l = l + 1
      if (lmax - l < 0) exit l_loop
    end do l_loop
  end subroutine multiple_scattering_matrix

  !> Solves a set of linear equations by Gaussian elimination with partial
  !! pivoting (legacy CXMTXT, option "system solution"). Author: Pendry.
  !!
  !! The coefficient matrix is a(1:nr, 1:nr), the nsys constant vectors are
  !! the columns a(1:nr, nr+1:nr+nsys). The matrix is destroyed (it becomes
  !! upper triangular) and the solutions overwrite the constant vectors.
  !! The pivot search and the threshold below which an element is not
  !! eliminated use the real moduli; singularity is not checked.
  subroutine solve_pivoted(a, nc, nr, nsys, work)
    !> First dimension of a, at least nr.
    integer, intent(in) :: nc
    !> Order of the matrix.
    integer, intent(in) :: nr
    !> Number of systems.
    integer, intent(in) :: nsys
    complex(sp), intent(inout) :: a(nc, nr + nsys)
    !> Work array, at least nc long.
    complex(sp), intent(out) :: work(nc)

    integer :: i, j, k, l, kmax, nplsy
    real(sp) :: amax, qz
    complex(sp) :: temp, const

    nplsy = nr + nsys

    ! Triangularize
    do i = 1, nr - 1
      ! Pivot element
      kmax = i
      amax = abs(a(i, i))
      do k = i + 1, nr
        qz = abs(a(k, i))
        if (amax - qz < 0.0_sp) then
          kmax = k
          amax = abs(a(k, i))
        end if
      end do
      if (kmax /= i) then  ! interchange rows
        do l = i, nplsy
          temp = a(i, l)
          a(i, l) = a(kmax, l)
          a(kmax, l) = temp
        end do
      end if
      ! Eliminate a(i+1, i) to a(nr, i)
      do j = i + 1, nr
        temp = a(j, i)
        qz = abs(temp)
        if (qz < 1.0e-10_sp) cycle
        const = -temp / a(i, i)
        do l = i, nplsy
          a(j, l) = a(j, l) + a(i, l) * const
        end do
      end do
    end do

    ! Back substitution
    do i = nr + 1, nplsy
      do k = nr, 1, -1
        temp = a(k, i)
        do j = k + 1, nr
          temp = temp - a(k, j) * work(j)
        end do
        work(k) = temp / a(k, k)
      end do
      do l = 1, nr
        a(l, i) = work(l)
      end do
    end do
  end subroutine solve_pivoted

  !> Direct-lattice sum for the propagators between two subplanes of a
  !! composite layer (legacy GHD). The subplanes may be coplanar.
  !!
  !! Computes the propagator for the interplanar vector number iz, from the
  !! first to the second subplane (is = 1) or from the second to the first
  !! (is = 2), and adds it to gh (see interplanar_propagators).
  subroutine direct_propagator_sum(iz, is, gh, lmg, lmmax, s, lms, y, l2m, &
                                   drl, nlay2, k0, dcut, caa, ncaa, lxm, lay, &
                                   pq, h, lat, point)
    !> Number of the interplanar vector.
    integer, intent(in) :: iz
    !> 1 for the propagation from the first to the second subplane, 2 for
    !! the other direction.
    integer, intent(in) :: is
    integer, intent(in) :: lmg, lmmax, lms, l2m, nlay2, ncaa
    complex(sp), intent(inout) :: gh(lmg, lmmax)
    !> Work array: lattice sum.
    complex(sp), intent(out) :: s(lms)
    !> Work array: spherical harmonics.
    complex(sp), intent(out) :: y(l2m, l2m)
    real(sp), intent(in) :: drl(nlay2, 3)
    !> Complex magnitude of the wave vector.
    complex(sp), intent(in) :: k0
    !> Cutoff radius of the lattice sum.
    real(sp), intent(in) :: dcut
    real(sp), intent(in) :: caa(ncaa)
    integer, intent(in) :: lxm(lmmax)
    !> 1 for an overlayer, otherwise a substrate layer.
    integer, intent(in) :: lay
    !> Parallel momentum of the first beam of the current beam set.
    real(sp), intent(in) :: pq(2)
    !> Work array: spherical Bessel functions.
    complex(sp), intent(inout) :: h(l2m)
    !> Lattices; uses substrate_vector and superlattice_vector.
    type(lattice), intent(in) :: lat
    !> Energy, parallel momentum and imaginary potential; uses k_parallel.
    type(energy_point), intent(in) :: point

    integer :: l, m, l1, kk, jj1, jj2, j1, j2, nr1, nr2, nor
    real(sp) :: dcut2, v(3), px, py, x1, x2, fy, cfy, ak2, ak3, ar1(2), ar2(2)
    complex(sp) :: ru, ci, cz, ff, z, z1, z2, z3, st

    ak2 = point%k_parallel(1)
    ak3 = point%k_parallel(2)
    if (lay == 1) then  ! overlayer
      ar1 = lat%superlattice_vector(:, 1)
      ar2 = lat%superlattice_vector(:, 2)
    else
      ar1 = lat%substrate_vector(:, 1)
      ar2 = lat%substrate_vector(:, 2)
    end if

    ru = cmplx(1.0_sp, 0.0_sp, sp)
    cz = cmplx(0.0_sp, 0.0_sp, sp)
    ci = cmplx(0.0_sp, 1.0_sp, sp)
    dcut2 = dcut * dcut

    s = cz

    v = drl(iz, :)
    ! Turn the interplanar vector around if is = 2
    if (is == 2) v = -v

    ! Two one-dimensional lattice loops for the summation over one quadrant
    jj1 = 0
    rows: do
      jj1 = jj1 + 1
      jj2 = 0
      columns: do
        jj2 = jj2 + 1
        nor = 0
        j1 = jj1 - 1
        j2 = jj2 - 1
        quadrants: do kk = 1, 4
          select case (kk)
          case (1)
            nr1 = j1
            nr2 = j2
          case (2)
            if (j1 == 0 .and. j2 == 0) exit quadrants
            if (j2 == 0) cycle quadrants
            nr1 = j1
            nr2 = -j2
          case (3)
            if (j1 == 0) cycle quadrants
            nr1 = -j1
            nr2 = j2
          case default
            if (j1 == 0 .or. j2 == 0) cycle quadrants
            nr1 = -j1
            nr2 = -j2
          end select

          px = real(nr1, sp) * ar1(1) + real(nr2, sp) * ar2(1)
          py = real(nr1, sp) * ar1(2) + real(nr2, sp) * ar2(2)
          x1 = (px + v(2))**2 + (py + v(3))**2

          ! Cutoff of the lattice summation at the radius dcut
          if (x1 > dcut2) cycle quadrants
          nor = 1
          z1 = exp(-ci * (px * (ak2 + pq(1)) + py * (ak3 + pq(2))))
          x2 = sqrt(x1 + v(1)**2)
          x1 = sqrt(x1)
          z2 = k0 * x2
          z = cmplx(v(1) / x2, 0.0_sp, sp)
          fy = 0.0_sp
          if (abs(x1) >= 1.0e-6_sp) then
            cfy = (px + v(2)) / x1
            if (abs(abs(cfy) - 1.0_sp) > 1.0e-6_sp) then
              fy = acos(cfy)
            else
              if (cfy < 0.0_sp) fy = pi
            end if
            if ((py + v(3)) < 0.0_sp) fy = -fy
          end if

          ! Spherical Bessel functions and spherical harmonics
          call spherical_bessel(z2, h, l2m)
          call harmonics_table(l2m, z, fy, y)

          ! Lattice sum into s
          st = ru
          do l = 1, l2m
            st = st * ci
            z3 = st * h(l) * z1
            l1 = l * l - l
            do m = 1, l
              s(l1 + m) = s(l1 + m) + z3 * y(l, m)
            end do
            do m = 1, l - 1
              s(l1 - m + 1) = s(l1 - m + 1) &
                  + z3 * y(m, l) * (-1.0_sp)**mod(m, 2)
            end do
          end do
        end do quadrants

        if (nor /= 1) exit columns
      end do columns
      if (jj2 == 1) exit rows
    end do rows

    ff = -8.0_sp * pi * k0
    ! Multiply the lattice sum into the Clebsch-Gordan coefficients
    call propagator_from_sum(iz, is, gh, lmg, lmmax, s, lms, caa, ncaa, ff, &
                             lxm, nlay2)
  end subroutine direct_propagator_sum

  !> Propagators between two subplanes from their direct-lattice sum
  !! (legacy GHSC), multiplied into the Clebsch-Gordan coefficients.
  !!
  !! iz is the number of the interplanar vector, is the direction (1: from
  !! the first to the second subplane, 2: the other way); ff is the
  !! prefactor.
  subroutine propagator_from_sum(iz, is, gh, lmg, lmmax, s, lms, caa, ncaa, &
                                 ff, lxm, nlay2)
    integer, intent(in) :: iz, is, lmg, lmmax, lms, ncaa, nlay2
    complex(sp), intent(inout) :: gh(lmg, lmmax)
    complex(sp), intent(in) :: s(lms)
    real(sp), intent(in) :: caa(ncaa)
    complex(sp), intent(in) :: ff
    integer, intent(in) :: lxm(lmmax)

    integer :: i, j, ii, ila, la, ip, jp, nii, njm, nlylm
    integer :: l1, l2, m1, m2, m3, il, im, lmin, lmax_sum, im1, imp, jm2, jmp
    real(sp) :: cc, phase

    ii = 1
    nlylm = lmmax * nlay2
    ! A limited set of matrix elements
    do i = 1, lmmax
      ip = lxm(i)
      nii = ip + (iz - 1) * lmmax
      if (is == 2) nii = nii + nlylm
      l1 = int(sqrt(real(i - 1, sp) + 0.00001_sp))
      m1 = i - l1 - l1 * l1 - 1
      do j = 1, lmmax
        jp = lxm(j)
        l2 = int(sqrt(real(j - 1, sp) + 0.00001_sp))
        m2 = j - l2 - l2 * l2 - 1
        m3 = m2 - m1
        il = abs(l1 - l2)
        im = abs(m3)
        lmin = max(il, im + mod(il + im, 2))
        lmax_sum = l1 + l2
        lmin = lmin + 1
        lmax_sum = lmax_sum + 1
        if (i == j) then
          if (m1 > 0) cycle
        else if (l1 < l2 .or. (l1 == l2 .and. abs(m1) < abs(m2))) then
          cycle
        end if
        do ila = lmin, lmax_sum, 2
          la = ila - 1
          cc = caa(ii)
          ii = ii + 1
          gh(nii, jp) = gh(nii, jp) + cc * s(la * la + la + m3 + 1)
        end do
        gh(nii, jp) = ff * gh(nii, jp)
      end do
    end do

    ! Fill in the missing matrix elements from the elements just produced
    do i = 1, lmmax
      ip = lxm(i)
      nii = ip + (iz - 1) * lmmax
      if (is == 2) nii = nii + nlylm
      l1 = int(sqrt(real(i - 1, sp) + 0.00001_sp))
      m1 = i - l1 - l1 * l1 - 1
      im1 = i - 2 * m1
      imp = lxm(im1)
      do j = 1, lmmax
        jp = lxm(j)
        l2 = int(sqrt(real(j - 1, sp) + 0.00001_sp))
        m2 = j - l2 - l2 * l2 - 1
        jm2 = j - 2 * m2
        jmp = lxm(jm2)
        njm = jmp + (iz - 1) * lmmax
        if (is == 2) njm = njm + nlylm
        phase = (-1.0_sp)**mod(m1 + m2, 2)
        if (i == j) then
          if (m1 <= 0) cycle
        else if (l1 >= l2 .and. (l1 /= l2 .or. abs(m1) >= abs(m2))) then
          cycle
        end if
        gh(nii, jp) = phase * gh(njm, imp)
      end do
    end do
  end subroutine propagator_from_sum

  !> Spherical Bessel functions (legacy SB): hh(1:n3) are the functions of
  !! the complex argument x for l = 0 to n3 - 1. hh is not changed if x is
  !! (nearly) zero.
  subroutine spherical_bessel(x, hh, n3)
    complex(sp), intent(in) :: x
    integer, intent(in) :: n3
    complex(sp), intent(inout) :: hh(n3)

    integer :: j
    real(sp) :: f
    complex(sp) :: a, b, c

    a = cmplx(0.0_sp, 1.0_sp, sp)
    c = x
    b = a * c
    f = abs(x)
    if (real(f, dp) - 1.0e-38_sp <= 0.0_dp) return
    hh(1) = exp(b) * (-a) / c
    hh(2) = exp(b) * (1.0_sp / (-c) - a / c**2)
    do j = 3, n3
      hh(j) = (2.0_sp * real(j - 2, sp) + 1.0_sp) / c * hh(j - 1) - hh(j - 2)
    end do
  end subroutine spherical_bessel

  !> Spherical harmonics of the direction with polar angle cosine z
  !! (complex) and azimuth fi (legacy SH).
  !!
  !! y(l, m) for the angular momenta up to nharm - 1 is stored in the
  !! lower triangle for m >= 0 and, conjugated up to a sign, in the upper
  !! triangle. The recursion runs in double precision; the result and the
  !! quantities that enter it as factors (bn, bm) are single precision.
  subroutine harmonics_table(nharm, z, fi, y)
    integer, intent(in) :: nharm
    complex(sp), intent(in) :: z
    real(sp), intent(in) :: fi
    complex(sp), intent(out) :: y(:, :)

    integer :: l, ll
    real(sp) :: rz, anora, bm, bn
    real(dp) :: am
    complex(sp) :: a, anor, q1a, ystar, yy, znw
    complex(dp) :: bb, bb1, aa

    am = 1.0_dp
    a = cmplx(0.0_sp, 1.0_sp, sp) * fi
    rz = real(z, sp)
    znw = sqrt(cmplx(1.0_sp, 0.0_sp, sp) - z * z)
    anora = 0.2820948_sp
    yy = cmplx(1.0_sp, 0.0_sp, sp)

    do l = 1, nharm
      bm = real(l, sp) - 1.0_sp
      q1a = exp(bm * a)
      anor = anora * q1a * (-1.0_sp)**(l + 1)
      anora = anora * sqrt(1.0_sp + 0.5_sp / (bm + 1.0_sp)) / (2.0_sp * bm + 1.0_sp)
      if ((abs(abs(rz) - 1.0_sp) >= 1.0e-10_sp) .or. l /= 1) then
        yy = znw**int(bm + 0.1_sp)
      end if
      aa = am * yy
      am = am * (2.0_sp * bm + 1.0_sp)
      bb = am * z * yy
      do ll = l, nharm
        bn = real(ll, sp)
        y(ll, l) = cmplx(aa * anor, kind=sp)
        if (l /= 1) then
          ystar = y(ll, l) / (q1a * q1a)
          y(l - 1, ll) = ystar
        end if
        bb1 = bb
        bb = bb1 * z + (bb1 * z - aa) * (bn + bm) / (bn - bm + 1.0_sp)
        aa = bb1
        anor = anor * sqrt((1.0_sp + 1.0_sp / (bn - 0.5_sp)) &
                           * ((bn - bm) / (bn + bm)))
      end do
    end do
  end subroutine harmonics_table
end module composite_layer
