!> \file composite_layer.f90
!! \brief Composite layers (several subplanes): ordering of the subplanes,
!! the propagators between them, and the reflection and transmission
!! matrices by Beeby's matrix inversion.
!!
!! Each routine keeps the precision of every quantity of the legacy routine
!! it follows, and its results are bit-for-bit those of the legacy code
!! (tests/tier1/test_refcalc_routines.py).

!> Subplanes of a composite layer, the propagators between them and the
!! multiple scattering within them; the reflection and transmission matrices
!! of the layer.
module composite_layer
  use constants, only: pi
  use bravais_layer, only: beams_ylm
  use kinds, only: sp, dp
  use leed_types, only: lattice, energy_point
  implicit none

  private

  public :: sort_subplanes, interplanar_propagators, subplane_scattering
  public :: copy_subplane_positions, composite_layer_matrices

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
  !! direct sum (whole) and by the reciprocal sum (its leading lm * lm
  !! elements, as an lm by lm table).
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
    complex(sp) :: ci, cz, kprg, k0, z, t1, t1p, t2, t2p, t3, cfac

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
            call harmonics_table(lm, z, fy, work)

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
              call add_reciprocal_terms(gh, lmg, lmmax, lxm, iz, nlylm, lm, &
                                        work, t1, t1p)
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
        call multiple_scattering_matrix(il, flm, x(:, 1:ll), lev, ll, tsf, it, &
                                        lm, lxi, lmmax, klm, clm, nlm, &
                                        clm_index, nnsub, lmax1)

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
        call solve_pivoted(x(:, 1:2 * ll), lev, ll, ll, xh)

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

  !> Adds the contribution of one reciprocal-lattice point to the limited
  !! set of matrix elements of the propagators between subplanes, for the
  !! interplanar vector iz, from the table y of spherical harmonics.
  subroutine add_reciprocal_terms(gh, lmg, lmmax, lxm, iz, nlylm, lm, y, t1, &
                                  t1p)
    integer, intent(in) :: lmg, lmmax, iz, nlylm, lm
    !> Propagators, stacked vertically.
    complex(sp), intent(inout) :: gh(lmg, lmmax)
    !> Permutation of the (l,m) sequence.
    integer, intent(in) :: lxm(lmmax)
    !> Table of spherical harmonics.
    complex(sp), intent(in) :: y(lm, lm)
    !> Phase factors of the two propagation directions.
    complex(sp), intent(in) :: t1, t1p

    integer :: i, j, ip, jp, nii, nni, l1, l2, m1, m2
    complex(sp) :: bs, cs

    do i = 1, lmmax
      ip = lxm(i)
      nii = ip + (iz - 1) * lmmax
      nni = nii + nlylm
      l1 = int(sqrt(real(i - 1, sp) + 0.00001_sp))
      m1 = i - l1 - l1 * l1 - 1
      if (m1 > 0) then
        bs = y(m1, l1 + 1)
      else
        bs = (-1.0_sp)**(mod(m1, 2)) * y(l1 + 1, -m1 + 1)
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
          cs = y(l2 + 1, m2 + 1)
        else
          cs = (-1.0_sp)**(mod(m2, 2)) * y(-m2, l2 + 1)
        end if
        ! gh(i, j) and gh(j, i) are produced together
        gh(nni, jp) = gh(nni, jp) + t1p * bs * cs
        gh(nii, jp) = gh(nii, jp) &
            + (-1.0_sp)**(mod(l1 + m1 + l2 + m2, 2)) * t1 * bs * cs
      end do
    end do
  end subroutine add_reciprocal_terms

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
    complex(sp), intent(out) :: y(nharm, nharm)

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

  !> Copies the positions of the subplanes of layer type iltype out of the
  !! array of all composite layer types (legacy COPYPOS), so that the
  !! dimensions are those of composite_layer_matrices.
  subroutine copy_subplane_positions(pos, subpos, nltype, nnsub, nlay, iltype)
    integer, intent(in) :: nltype, nnsub, nlay, iltype
    !> Positions of the subplanes of all layer types.
    real(sp), intent(in) :: subpos(nltype, nnsub, 3)
    !> Positions of the subplanes of the layer type.
    real(sp), intent(out) :: pos(nlay, 3)

    integer :: ipos, ilay

    do ipos = 1, 3
      do ilay = 1, nlay
        pos(ilay, ipos) = subpos(iltype, ilay, ipos)
      end do
    end do
  end subroutine copy_subplane_positions

  !> Reflection and transmission matrices of a composite layer of nlay
  !! subplanes for the beams of one beam set (legacy RTINV_SIMPLE), by
  !! Beeby's matrix inversion. Beam symmetry codes and registry shifts are
  !! not supported.
  !!
  !! Sorts the subplanes (sort_subplanes), computes the propagators between
  !! them (interplanar_propagators) and the scattering within each of them
  !! (subplane_scattering), builds and factorizes the matrix to invert
  !! (beeby_matrix), and for each incident beam finds the scattering
  !! amplitudes of the subplanes, the elements of the reflection and
  !! transmission matrices, and the amplitudes stored for tensor LEED.
  !!
  !! r_top, t_top, r_bot, t_bot: reflection and transmission matrices for
  !! beams incident from above (top) or below (bot) the composite layer;
  !! element (i, j) connects incident beam j to beam i. The amplitudes of
  !! both beams are taken right next to the layer: r_top(i, j) beams right
  !! above; t_top(i, j) j right above, i right below; r_bot(i, j) both beams
  !! right below; t_bot(i, j) j right below, i right above. The entries of
  !! the beam set are set at the rows and columns ns + 1 ... ns + n; the
  !! others are kept.
  !!
  !! As in bravais_layer_matrices, the rows of amult and cylm are indexed
  !! by the beam within its set (without n_offset), so with several beam
  !! sets each call overwrites rows 1 ... n of the previous one. The stored
  !! amplitudes tstore(:, :, jgp) are indexed by the beam jgp within its set
  !! as well, without ns: with several beam sets each call overwrites the
  !! columns of the previous one.
  !!
  !! The beam-space matrices are packed with leading dimension nt, the
  !! number of beams at the current energy, into the leading nt * nt
  !! elements of the (larger) arrays of the caller; rg, th, tau, gh and the
  !! work arrays likewise take the leading part of the storage of larger
  !! arrays.
  !!
  !! LAPACK status flags are not checked.
  subroutine composite_layer_matrices(r_top, t_top, r_bot, t_bot, n, amult, &
                                      cylm, pq, nt, flms, flm, nl, lxi, lt, &
                                      lxm, lmmax, klm, xev, lev, lev2, tau, &
                                      lmt, taug, taugm, clm, nlm, pos, poss, &
                                      mgh, nlay, drl, nlay2, nugh, ngeq, test, &
                                      gh, lmg, rg, rg_perp, ts, lmn, tg, lm2n, &
                                      vt, caa, ncaa, th, lmni, ipl, tsf, lps, &
                                      lpss, nord, nnsub, lmax1, xh, hghd, l2m, &
                                      rg1, rg2, tstore, n_offset, ns, lay, &
                                      ntau, tst, tv, dcut, lmax, lat, point)
    !> Number of beams in the current beam set.
    integer, intent(in) :: n
    !> Total number of beams at the current energy.
    integer, intent(in) :: nt
    !> Number of sublattices (rows of flms).
    integer, intent(in) :: nl
    !> (lmax + 1)**2.
    integer, intent(in) :: lmmax
    !> (2 * lmax + 1) * (2 * lmax + 2) / 2.
    integer, intent(in) :: klm
    !> (lmax + 1) * (lmax + 2) / 2.
    integer, intent(in) :: lev
    !> 2 * lev.
    integer, intent(in) :: lev2
    !> ntau * lmmax.
    integer, intent(in) :: lmt
    !> Dimension of clm.
    integer, intent(in) :: nlm
    !> Number of subplanes.
    integer, intent(in) :: nlay
    !> nlay * (nlay - 1) / 2.
    integer, intent(in) :: nlay2
    !> 2 * nlay2 * lmmax.
    integer, intent(in) :: lmg
    !> nlay * lmmax.
    integer, intent(in) :: lmn
    !> 2 * lmn.
    integer, intent(in) :: lm2n
    !> Number of Clebsch-Gordan coefficients.
    integer, intent(in) :: ncaa
    !> lmn.
    integer, intent(in) :: lmni
    !> Dimensions of tsf.
    integer, intent(in) :: nnsub, lmax1
    !> 2 * lmax + 1, or larger.
    integer, intent(in) :: l2m
    !> Offset of the current beam set in the list pq.
    integer, intent(in) :: n_offset
    !> Offset of the current beam set in the matrices r_top ... t_bot: the
    !! sum of the sizes of the earlier beam sets.
    integer, intent(in) :: ns
    !> 1 for an overlayer, otherwise a substrate layer.
    integer, intent(in) :: lay
    !> Number of chemical elements.
    integer, intent(in) :: ntau
    !> Highest angular momentum.
    integer, intent(in) :: lmax
    !> Convergence criterion of the reciprocal-lattice sums.
    real(sp), intent(in) :: tst
    !> Area of the unit cell.
    real(sp), intent(in) :: tv
    !> Cutoff radius of the direct-lattice sums.
    real(sp), intent(in) :: dcut
    !> Permutations of the (l,m) sequence.
    integer, intent(in) :: lxi(lmmax), lt(lmmax), lxm(lmmax)
    !> Chemical element of each subplane.
    integer, intent(in) :: lps(nlay)
    !> Beams (reciprocal lattice vectors) at the current energy.
    real(sp), intent(in) :: pq(2, nt)
    !> Clebsch-Gordan coefficients.
    real(sp), intent(in) :: clm(nlm), caa(ncaa)
    !> Positions of the subplanes.
    real(sp), intent(in) :: pos(nlay, 3)
    !> Atomic t-matrix elements.
    complex(sp), intent(in) :: tsf(nnsub, lmax1)
    !> Lattice sums.
    complex(sp), intent(in) :: flms(nl, klm)
    !> Lattices.
    type(lattice), intent(in) :: lat
    !> Energy, parallel momentum of the incident beam, imaginary potential.
    type(energy_point), intent(in) :: point
    !> Reflection and transmission matrices; the entries of the beam set are
    !! set, the others are kept.
    complex(sp), intent(inout) :: r_top(nt, nt), t_top(nt, nt)
    complex(sp), intent(inout) :: r_bot(nt, nt), t_bot(nt, nt)
    !> Prefactors of the matrix elements of the beams of the set.
    complex(sp), intent(out) :: amult(n)
    !> Spherical harmonics of the beams of the set, (beam within its set,
    !! (l, m)); rows 1 ... n are set, the others are kept.
    complex(sp), intent(inout) :: cylm(nt, lmmax)
    !> Amplitudes of the scattering within the subplanes for each incident
    !! beam, for tensor LEED: (1, :, jgp) for incidence towards the solid,
    !! (2, :, jgp) away from it; the subplanes are in their original order.
    complex(sp), intent(inout) :: tstore(2, lmn, nt)
    !> Subplanes sorted by position (see sort_subplanes).
    real(sp), intent(out) :: poss(nlay, 3)
    integer, intent(out) :: lpss(nlay), nord(nlay)
    integer, intent(out) :: mgh(nlay, nlay)
    real(sp), intent(out) :: drl(nlay2, 3)
    integer, intent(out) :: nugh(nlay2), ngeq(nlay2)
    real(sp), intent(out) :: test(nlay2)
    !> Propagators between the subplanes, in (l,m) space.
    complex(sp), intent(inout) :: gh(lmg, lmmax)
    !> Tau matrices of the chemical elements.
    complex(sp), intent(inout) :: tau(lmt, lev)
    !> Matrix to invert; on return its LU factors.
    complex(sp), intent(out) :: th(lmni, lmni)
    !> Pivot indices of the LU factorization.
    integer, intent(out) :: ipl(lmni)
    !> Plane-wave propagators between the subplanes and the layer limits.
    complex(sp), intent(out) :: rg(2, nlay, n), rg_perp(n)
    !> Work spaces.
    complex(sp), intent(out) :: flm(klm), xev(lev, lev2), xh(lev)
    complex(sp), intent(out) :: taug(lmt), taugm(lmt), ts(lmn)
    complex(sp), intent(out) :: tg(2 * lm2n), vt(lm2n)
    complex(sp), intent(inout) :: hghd(l2m)
    complex(sp), intent(out) :: rg1(nlay), rg2(nlay)

    integer :: ig, jg, jgp, jgps, jgs, i, k, in, lp, insort, inc
    integer :: lod, lee, loe, lm, lms, info
    real(sp) :: bk2, bk3, c, x, e, vpi, ak2, ak3, gp(2)
    complex(sp) :: xa, st, am
    complex(sp) :: rt(2)
    complex(sp), parameter :: ci = cmplx(0.0_sp, 1.0_sp, sp)
    complex(sp), parameter :: ru = cmplx(1.0_sp, 0.0_sp, sp)
    external :: cgetrf, cgetrs

    e = point%energy
    vpi = point%imaginary_potential
    ak2 = point%k_parallel(1)
    ak3 = point%k_parallel(2)
    lm = lmax + 1
    lms = l2m * l2m
    lod = lmmax - lev
    lee = (lmax / 2 + 1)**2
    loe = ((lmax - 1) / 2 + 1)**2

    ! Tau for each chemical element
    call subplane_scattering(tau, lmt, ntau, xev, lev, lev2, lod, tsf, lmmax, &
                             flms, nl, klm, lm, clm, nlm, lxi, &
                             pq(:, 1 + n_offset), nnsub, lmax1, xh, flm, lay, &
                             lat, point)

    ! Sort the subplanes by increasing position along +x, i.e., from the
    ! vacuum side towards the solid; find which propagators can be kept,
    ! which must be recomputed and which are mutually identical
    call sort_subplanes(pos, poss, lps, lpss, mgh, nlay, drl, nlay2, nugh, &
                        ngeq, nord)

    ! Propagators between the subplanes in (l,m) space
    call interplanar_propagators(gh, lmg, lmmax, nugh, ngeq, nlay2, tst, test, &
                                 vt, l2m, lm, tg, lms, drl, tv, lxm, dcut, &
                                 caa, ncaa, lay, hghd, pq(:, 1 + n_offset), &
                                 lat, point)

    ! Matrix to be inverted in Beeby's scheme
    call beeby_matrix(th, lmni, gh, lmg, lmmax, mgh, nlay, tau, lmt, lev, lpss)

    ! Gaussian elimination
    call cgetrf(lmni, lmni, th, lmni, ipl, info)

    ! Spherical harmonics of the beam directions; amult is overwritten below
    call beams_ylm(cylm, amult, pq, nt, n, n_offset, point, lmax, lmmax, vt)
    do ig = 1, n
      ! Plane-wave propagators, used as prefactors for gh and th
      jg = ig + n_offset
      bk2 = pq(1, jg) + ak2
      bk3 = pq(2, jg) + ak3
      c = bk2 * bk2 + bk3 * bk3
      xa = cmplx(2.0_sp * e - c, -2.0_sp * vpi + 0.000001_sp, sp)
      xa = sqrt(xa)
      rg_perp(ig) = exp(ci * xa * (poss(nlay, 1) - poss(1, 1)))
      do i = 1, nlay
        x = bk2 * poss(i, 2) + bk3 * poss(i, 3)
        rg(1, i, ig) = exp(ci * (xa * poss(i, 1) + x))
        rg(2, i, ig) = exp(ci * (xa * poss(i, 1) - x))
      end do
      ! Prefactor of the reflection and transmission matrix elements
      amult(ig) = -16.0_sp * pi * pi * ci / (tv * xa)
    end do

    ! Loop over incident beams
    do jgp = 1, n
      ! Tau times the spherical harmonics, into which the inverse of th is
      ! multiplied
      call tau_ylm(taug, taugm, lmt, tau, lev, cylm, nt, lmmax, lt, ntau, lod, &
                   lee, loe, jgp)

      jgps = jgp + ns
      gp(1) = pq(1, jgp + n_offset)
      gp(2) = pq(2, jgp + n_offset)
      ! First incidence towards +x, i.e., beams propagating towards the
      ! solid; then incidence towards -x
      do inc = 1, -1, -2
        ! Plane-wave propagating factors into taug/taugm
        do i = 1, nlay
          in = (i - 1) * lmmax
          lp = (lpss(i) - 1) * lmmax
          if (inc == 1) then
            st = rg(1, i, jgp)
            do k = 1, lmmax
              ts(in + k) = taug(lp + k) * st
            end do
          else
            st = ru / rg(2, i, jgp)
            do k = 1, lmmax
              ts(in + k) = taugm(lp + k) * st
            end do
          end if
        end do

        ! Left multiplication with the inverse of th
        call cgetrs('N', lmni, 1, th, lmni, ipl, ts, lmni, info)

        ! Further plane-wave propagating factors; and store the
        ! intralayer scattering for the tensors, with the sublayers in
        ! their original order
        do i = 1, nlay
          in = (i - 1) * lmmax
          insort = (nord(i) - 1) * lmmax
          if (inc == 1) then
            st = ru / rg(1, i, jgp)
          else
            st = rg(2, i, jgp)
          end if
          do k = 1, lmmax
            if (inc == 1) then
              tstore(1, insort + k, jgp) = 4.0_sp * pi * ts(in + k)
            else
              tstore(2, insort + k, jgp) = 4.0_sp * pi * ts(in + k) &
                                           * rg_perp(jgp)
            end if
            ts(in + k) = ts(in + k) * st
          end do
        end do

        ! Loop over scattered beams
        do jg = 1, n
          call reflection_transmission_element(jg, n_offset, point%k_parallel, &
                                               gp, lm, cylm, n, lmmax, pq, nt, &
                                               ts, lmn, rg, rg1, rg2, nlay, &
                                               jgp, lxm, rt, inc, poss)
          jgs = jg + ns
          am = amult(jg)
          if (inc == -1) then
            r_bot(jgs, jgps) = rt(1) * am
            t_bot(jgs, jgps) = rt(2) * am
          else
            r_top(jgs, jgps) = rt(1) * am
            t_top(jgs, jgps) = rt(2) * am
          end if
        end do

        ! The unscattered plane wave adds to the transmission
        if (inc == 1) then
          t_top(jgps, jgps) = t_top(jgps, jgps) + rg_perp(jgp)
        else
          t_bot(jgps, jgps) = t_bot(jgps, jgps) + rg_perp(jgp)
        end if
      end do
    end do
  end subroutine composite_layer_matrices

  !> Matrix that is inverted in Beeby's scheme (legacy THMAT): the identity
  !! minus tau times the propagators between the subplanes.
  !!
  !! th: the matrix, of dimension lmni = nlay * lmmax. gh: the propagators,
  !! organized by mgh (see sort_subplanes). lps: the chemical element of
  !! each subplane, sorted.
  subroutine beeby_matrix(th, lmni, gh, lmg, lmmax, mgh, nlay, tau, lmt, lev, &
                          lps)
    integer, intent(in) :: lmni, lmg, lmmax, nlay, lmt, lev
    integer, intent(in) :: mgh(nlay, nlay), lps(nlay)
    complex(sp), intent(in) :: gh(lmg, lmmax), tau(lmt, lev)
    complex(sp), intent(out) :: th(lmni, lmni)

    integer :: i, j, k, in, inp, in1, ic, icl, id, idl, m, m1, lp, lp1
    integer :: icv, idv, nlay1, lod
    complex(sp) :: st, su
    complex(sp), parameter :: ru = cmplx(1.0_sp, 0.0_sp, sp)
    complex(sp), parameter :: cz = cmplx(0.0_sp, 0.0_sp, sp)

    lod = lmmax - lev

    ! Start from the identity matrix
    do j = 1, lmni
      do i = 1, lmni
        th(i, j) = cz
      end do
      th(j, j) = ru
    end do

    ic = 0
    nlay1 = nlay - 1
    do in = 1, nlay1
      ic = ic + 1
      icl = (ic - 1) * lmmax
      id = ic
      in1 = in + 1
      do inp = in1, nlay
        id = id + 1
        idl = (id - 1) * lmmax
        m = (mgh(in, inp) - 1) * lmmax
        m1 = (mgh(inp, in) - 1) * lmmax
        lp = (lps(in) - 1) * lmmax
        lp1 = (lps(inp) - 1) * lmmax

        ! -tau * gh for even l+m (in the first index)
        do i = 1, lev
          do j = 1, lmmax
            st = cz
            su = cz
            do k = 1, lev
              st = st + tau(lp + i, k) * gh(m + k, j)
              su = su + tau(lp1 + i, k) * gh(m1 + k, j)
            end do
            th(icl + i, idl + j) = -st
            th(idl + i, icl + j) = -su
          end do
        end do

        lp = lp + lev
        lp1 = lp1 + lev
        m = m + lev
        m1 = m1 + lev
        icv = icl + lev
        idv = idl + lev

        ! -tau * gh for odd l+m (in the first index)
        do i = 1, lod
          do j = 1, lmmax
            st = cz
            su = cz
            do k = 1, lod
              st = st + tau(lp + i, k) * gh(m + k, j)
              su = su + tau(lp1 + i, k) * gh(m1 + k, j)
            end do
            th(icv + i, idl + j) = -st
            th(idv + i, icl + j) = -su
          end do
        end do
      end do
    end do
  end subroutine beeby_matrix

  !> Tau times the spherical harmonics of the incident beam jgp, for each
  !! chemical element (legacy TAUY_SIMPLE); the input of the Beeby inversion.
  !!
  !! taug and taugm are the vectors for the beam moving towards the solid
  !! and away from it, of dimension lmt = ntau * lmmax. lt is the
  !! permutation of the (l,m) sequence. lee, loe: number of (l,m)
  !! combinations of even l+m up to which, and of odd l+m up to which, the
  !! sign of the contribution is reversed.
  subroutine tau_ylm(taug, taugm, lmt, tau, lev, cylm, nt, lmmax, lt, ntau, &
                     lod, lee, loe, jgp)
    integer, intent(in) :: lmt, lev, nt, lmmax, ntau, lod, lee, loe, jgp
    integer, intent(in) :: lt(lmmax)
    complex(sp), intent(in) :: tau(lmt, lev), cylm(nt, lmmax)
    complex(sp), intent(out) :: taug(lmt), taugm(lmt)

    integer :: i, is, jlm, ilm, jtau, klp
    complex(sp) :: cz, st, su, cf

    cz = cmplx(0.0_sp, 0.0_sp, sp)

    ! Matrix product tau * ylm for each chemical element
    do i = 1, ntau
      is = (i - 1) * lmmax
      do jlm = 1, lev
        st = cz
        jtau = is + jlm
        do ilm = 1, lev
          klp = lt(ilm)
          if (ilm > lee) then
            st = st - tau(jtau, ilm) * cylm(jgp, klp)
          else
            st = st + tau(jtau, ilm) * cylm(jgp, klp)
          end if
        end do
        taug(jtau) = st
        taugm(jtau) = st
      end do

      is = is + lev
      do jlm = 1, lod
        st = cz
        su = cz
        jtau = is + jlm
        do ilm = 1, lod
          klp = lt(ilm + lev)
          cf = tau(jtau, ilm) * cylm(jgp, klp)
          if (ilm > loe) then
            st = st - cf
            su = su + cf
          else
            st = st + cf
            su = su - cf
          end if
        end do
        taug(jtau) = st
        taugm(jtau) = su
      end do
    end do
  end subroutine tau_ylm

  !> Reflection (rt_out(1)) and transmission (rt_out(2)) of incident beam
  !! jgp into the beam jg of the same beam set, for a composite layer
  !! (legacy MFOLT_SIMPLE). Beam symmetry codes and registry shifts are not
  !! supported.
  !!
  !! inc = +1: the incident beam moves towards the solid; -1: away from it.
  !! gp and ak: reciprocal-space part of the incident beam and parallel
  !! momentum of the primary beam. ts: amplitudes scattered by each
  !! (l,m) combination of each subplane, from composite_layer_matrices.
  !! rg(1, ...) is the propagation into the solid, rg(2, ...) towards
  !! vacuum; rg_r and rg_t are work spaces for the propagation to each
  !! subplane relative to incidence or exit at the extremes of the layer,
  !! for reflection and transmission. poss: the positions of the subplanes
  !! (sorted from the surface), z first. lxm: permutation of the (l,m)
  !! sequence. lm: lmax + 1.
  subroutine reflection_transmission_element(jg, na, ak, gp, lm, cylm, n, &
                                             lmmax, pq, nt, ts, lmn, rg, rg_r, &
                                             rg_t, nlay, jgp, lxm, rt_out, &
                                             inc, poss)
    integer, intent(in) :: jg, na, lm, n, lmmax, nt, lmn, nlay, jgp, inc
    integer, intent(in) :: lxm(lmmax)
    real(sp), intent(in) :: ak(2), gp(2), pq(2, nt), poss(nlay, 3)
    complex(sp), intent(in) :: cylm(nt, lmmax), ts(lmn), rg(2, nlay, n)
    complex(sp), intent(out) :: rg_r(nlay), rg_t(nlay), rt_out(2)

    integer :: jga, i1, i2, j, l, ll, m, jlm, klm, l1, l2
    real(sp) :: k_in(2), k_out(2), br, bt
    complex(sp) :: st, sm, cy, ctr, ctt, r, t, cr, ct, cs, edw_r, edw_t
    complex(sp), parameter :: cz = cmplx(0.0_sp, 0.0_sp, sp)
    complex(sp), parameter :: ru = cmplx(1.0_sp, 0.0_sp, sp)
    complex(sp), parameter :: ci = cmplx(0.0_sp, 1.0_sp, sp)

    jga = jg + na
    k_out(1) = pq(1, jga) + ak(1)
    k_out(2) = pq(2, jga) + ak(2)
    k_in(1) = gp(1) + ak(1)
    k_in(2) = gp(2) + ak(2)

    if (inc == 1) then
      i1 = 1
      i2 = nlay
    else
      i1 = nlay
      i2 = 1
    end if

    ! Total in-plane propagation phase shift for reflection at the first
    ! subplane encountered (i1), with both k_in and k_out there
    br = (k_in(1) - k_out(1)) * poss(i1, 2) &
         + (k_in(2) - k_out(2)) * poss(i1, 3)

    ! Total in-plane propagation phase shift for transmission between the
    ! first (i1) and the last (i2) subplanes encountered: k_in at the
    ! first, k_out at the last
    bt = k_in(1) * poss(i1, 2) + k_in(2) * poss(i1, 3) &
         - (k_out(1) * poss(i2, 2) + k_out(2) * poss(i2, 3))

    ! In-plane propagators for reflection and transmission
    edw_r = exp(ci * br)
    edw_t = exp(ci * bt)

    ! Effect of propagation on the scattered amplitudes, to be combined
    ! with the (l,m) contributions below
    do j = 1, nlay
      if (inc == 1) then
        ! Fraction of the incident beam reaching plane j relative to the
        ! first, both moving into the solid
        st = rg(1, j, jgp) / rg(1, 1, jgp)
        ! Fraction of the scattered beam leaving plane j relative to the
        ! one reaching the first, both moving out
        rg_r(j) = st * rg(2, j, jg) / rg(2, 1, jg)
        ! Fraction of the scattered beam leaving the last plane relative to
        ! the one reaching plane j, both moving into the solid
        rg_t(j) = st * rg(1, nlay, jg) / rg(1, j, jg)
      else
        ! Fraction of the incident beam reaching the last plane relative to
        ! the one at plane j, both moving out
        st = rg(2, nlay, jgp) / rg(2, j, jgp)
        ! Fraction of the scattered beam leaving the last plane relative to
        ! the one reaching plane j, both moving out
        rg_r(j) = st * rg(1, nlay, jg) / rg(1, j, jg)
        ! Fraction of the scattered beam leaving plane j relative to the
        ! one leaving the first, both moving into the solid
        rg_t(j) = st * rg(2, j, jg) / rg(2, 1, jg)
      end if
    end do

    ! Combine propagation with the (l,m)-dependent scattering. The factor
    ! sm = (-1)**(l+m) * (-1)*(-1)**l = -(-1)**m replaces the products of
    ! older code.
    r = cz
    t = cz
    jlm = 1
    l1 = 1
    do l = 1, lm
      sm = ru
      ll = l + l - 1
      l1 = l1 + ll
      l2 = l1 - 2 * l
      do m = 1, ll
        klm = lxm(jlm)
        ! Sum over all subplanes
        cr = cz
        ct = cz
        do j = 1, nlay
          cs = ts((j - 1) * lmmax + klm)
          cr = cr + rg_r(j) * cs
          ct = ct + rg_t(j) * cs
        end do

        ! Spherical harmonics and further prefactors
        cy = cylm(jg, l2 + m)
        if (inc == 1) then
          ctr = sm
          ctt = ru
        else
          ctr = ru
          ctt = sm
        end if
        r = r + ctr * cy * cr
        t = t + ctt * cy * ct
        jlm = jlm + 1
        sm = -sm
      end do
    end do

    rt_out(1) = r * edw_r
    rt_out(2) = t * edw_t
  end subroutine reflection_transmission_element
end module composite_layer
