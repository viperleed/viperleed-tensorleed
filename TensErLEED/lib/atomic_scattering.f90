!> \file atomic_scattering.f90
!! \brief Atomic scattering: temperature-dependent phase shifts and the atomic
!! t-matrices.
!!
!! As in angular_coefficients, every routine keeps the precision of each
!! quantity of the legacy routine it follows, and its results are bit-for-bit
!! those of the legacy code (tests/tier1/test_refcalc_routines.py).

!> Atomic t-matrices.
module atomic_scattering
  use kinds, only: sp, dp
  use leed_output, only: stop_with_error
  implicit none

  private

  public :: legendre_product_integrals, atomic_tmatrix
  public :: accumulate_atomic_tmatrix

contains

  !> Tabulates ppp(i1, i2, i3), each element containing the integral of the
  !! product of three Legendre functions P(i1), P(i2), P(i3) (legacy CPPP).
  !!
  !! The integrals follow Gaunt's summation scheme, set out by Slater,
  !! Atomic Structure, vol. 1, 309, 310. They are used to compute
  !! temperature-dependent phase shifts. Author: Pendry.
  !!
  !! The dummy ppp is explicitly shaped: its shape (n1, n2, n3) equals the
  !! PARAM shape of the ref-calc actual, because the input reader requires
  !! LMAX to equal MLMAX (sizing from LMAX comes later).
  !! The factorial table f is workspace: the routine fills it before use and
  !! no caller reads it afterwards. f(n + 1) is set to n!, in double
  !! precision (REAL*8 on purpose: it avoids a divergence of the sums to
  !! infinity). The largest n1 + n2 + n3 - 1 = 4 * lmax + 2 needs
  !! size(f) >= 4 * lmax + 3.
  subroutine legendre_product_integrals(n1, n2, n3, ppp, f)
    integer, intent(in) :: n1, n2, n3
    real(sp), intent(out) :: ppp(n1, n2, n3)
    real(dp), intent(out) :: f(:)

    integer :: i, i1, i2, i3, im1, im2, im3, is, it, j
    integer :: ia1, ia2, ia3, ia4, ia5
    real(sp) :: a, sign_
    real(dp) :: sum_

    f(1) = 1.0_dp
    do i = 1, size(f) - 1
      f(i + 1) = f(i) * real(i, dp)
    end do

    do i3 = 1, n3
      do i2 = 1, n2
        do i1 = 1, n1
          ! Order the three indices: im1 >= im2 >= im3
          im1 = i1
          im2 = i2
          im3 = i3
          if (i1 < i2) then
            im1 = i2
            im2 = i1
          end if
          if (im2 < i3) then
            im3 = im2
            im2 = i3
            if (im1 < im2) then
              j = im1
              im1 = im2
              im2 = j
            end if
          end if

          a = 0.0_sp
          is = i1 + i2 + i3 - 3
          if (mod(is, 2) /= 1 .and. abs(im2 - im1) - im3 + 1 <= 0) then
            sum_ = 0.0_dp
            is = is / 2
            sign_ = 1.0_sp
            do it = 1, im3
              sign_ = -sign_
              ia1 = im1 + it - 1
              ia2 = im1 - it + 1
              ia3 = im3 - it + 1
              ia4 = im2 + im3 - it
              ia5 = im2 - im3 + it
              sum_ = sum_ - sign_ * f(ia1) * f(ia4) &
                            / (f(ia2) * f(ia3) * f(ia5) * f(it))
            end do
            ia1 = 2 + is - im1
            ia2 = 2 + is - im2
            ia3 = 2 + is - im3
            ia4 = 3 + 2 * (is - im3)
            a = -(-1.0_sp)**(is - im2) * f(ia4) * f(is + 1) * f(im3) * sum_ &
                / (f(ia1) * f(ia2) * f(ia3) * f(2 * is + 2))
          end if
          ppp(i1, i2, i3) = a
        end do
      end do
    end do
  end subroutine legendre_product_integrals

  !> Adds the t-matrix of one element, weighted with its concentration, to
  !! the average t-matrix of a site (the average-t-matrix approximation).
  !! Part of the energy loop of the driver, tmat(il) += conc * caf(il).
  subroutine accumulate_atomic_tmatrix(conc, caf, tmat)
    !> Concentration of the element on the site.
    real(sp), intent(in) :: conc
    !> Temperature-dependent atomic t-matrix of the element, (l1).
    complex(sp), intent(in) :: caf(:)
    !> Average t-matrix of the site, (l1).
    complex(sp), intent(inout) :: tmat(:)

    integer :: il

    do il = 1, size(tmat)
      tmat(il) = tmat(il) + conc * caf(il)
    end do
  end subroutine accumulate_atomic_tmatrix

  !> Interpolates tabulated phase shifts to the current energy and produces
  !! the atomic t-matrix elements (legacy TSCATF): af without and caf with
  !! the correction for thermal vibrations. The energy is eb - v (v can
  !! describe local variations of the muffin-tin constant), and the
  !! phase shifts are extrapolated above the tabulated range. Below it the
  !! program stops with an error.
  !!
  !! The dimensions are those of the legacy routine: nn2 = nn3 = l1 is the
  !! number of phase shifts, nn1 = nn2 + nn3 - 1, ppp the Gaunt table of
  !! legendre_product_integrals, and (phs, del, ctab, sum_, bj) are
  !! workspace. phs and del are also output: the interpolated and the
  !! temperature-dependent phase shifts.
  !! The vibration amplitudes drper and drpar (perpendicular and parallel
  !! to the surface) are averaged to an isotropic one, valid at the
  !! temperature t0; t is the current temperature, dr0 the fourth power of
  !! the zero-temperature amplitude.
  subroutine atomic_tmatrix(iel, l1, npsi, nel, nn1, nn2, nn3, es, phss, eb, v, &
                            ppp, dr0, drper, drpar, t0, t, af, caf, phs, del, &
                            ctab, sum_, bj)
    !> Element, in the order of the phase-shift input.
    integer, intent(in) :: iel
    !> lmax + 1.
    integer, intent(in) :: l1
    !> Number of energies of the phase-shift table.
    integer, intent(in) :: npsi
    !> Number of elements of the phase-shift table.
    integer, intent(in) :: nel
    integer, intent(in) :: nn1, nn2, nn3
    !> Energies of the table.
    real(sp), intent(in) :: es(npsi)
    !> Tabulated phase shifts, (npsi, nel, l1).
    real(sp), intent(in) :: phss(npsi, nel, l1)
    real(sp), intent(in) :: eb, v
    real(sp), intent(in) :: ppp(nn1, nn2, nn3)
    real(sp), intent(in) :: dr0, drper, drpar, t0, t
    complex(sp), intent(out) :: af(l1), caf(l1)
    real(sp), intent(out) :: phs(l1)
    complex(sp), intent(out) :: del(l1)
    complex(sp), intent(out) :: ctab(nn3), sum_(nn2)
    complex(dp), intent(out) :: bj(nn1)

    complex(sp), parameter :: ci = (0.0_sp, 1.0_sp)
    integer :: i, l
    real(sp) :: a, dr, e, fac
    complex(sp) :: ca

    e = eb - v
    if (e < es(1)) call stop_with_error('TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS')

    ! Find the interval of the table that contains e, or the last one
    i = 1
    do
      if ((e - es(i)) * (e - es(i + 1)) <= 0.0_sp) exit
      i = i + 1
      if (i >= npsi) then
        i = i - 1
        exit
      end if
    end do
    fac = (e - es(i)) / (es(i + 1) - es(i))
    do l = 1, l1
      phs(l) = phss(i, iel, l) + fac * (phss(i + 1, iel, l) - phss(i, iel, l))
    end do

    ! Temperature-independent t-matrix elements
    do l = 1, l1
      a = phs(l)
      af(l) = a * ci
      af(l) = exp(af(l))
      a = sin(a)
      af(l) = a * af(l)
    end do

    ! Average any anisotropy of the vibration amplitudes
    dr = sqrt((drper * drper + 2.0_sp * drpar * drpar) / 3.0_sp)
    call thermal_phase_shifts(nn1, nn2, nn3, ppp, dr0, dr, t0, t, e, phs, del, &
                              ctab, sum_, bj)

    ! Temperature-dependent t-matrix elements
    do l = 1, l1
      ca = del(l)
      caf(l) = ca * ci
      caf(l) = exp(caf(l))
      ca = sin(ca)
      caf(l) = ca * caf(l)
    end do
  end subroutine atomic_tmatrix

  !> Incorporates the thermal vibrations in the phase shifts through a
  !! Debye-Waller factor, for isotropic vibration amplitudes (legacy PSTEMP).
  !!
  !! n3 is the number of input phase shifts phs, n2 the number of output
  !! (complex) temperature-dependent phase shifts del, n1 = n2 + n3 - 1, and
  !! ppp the Gaunt table. dr0 is the fourth power of the root-mean-square
  !! zero-temperature amplitude, dr the isotropic amplitude at the reference
  !! temperature t0, temp the temperature and e the energy. ctab, sum_ and
  !! bj are workspace. After the sums of two consecutive values of l have
  !! been below the threshold, the remaining elements of del stay zero.
  !! The spherical Bessel functions come from the series of spherical_bessel.
  !! (Uses BESSEL of the TensErLEED package for a better convergence at
  !! larger vibration amplitudes than the original scheme, U. Löffler, 1995.)
  subroutine thermal_phase_shifts(n1, n2, n3, ppp, dr0, dr, t0, temp, e, phs, &
                                  del, ctab, sum_, bj)
    integer, intent(in) :: n1, n2, n3
    real(sp), intent(in) :: ppp(n1, n2, n3)
    real(sp), intent(in) :: dr0, dr, t0, temp, e
    real(sp), intent(in) :: phs(n3)
    complex(sp), intent(out) :: del(n2)
    complex(sp), intent(out) :: ctab(n3), sum_(n2)
    complex(dp), intent(out) :: bj(n1)

    complex(sp), parameter :: ci = (0.0_sp, 1.0_sp)
    integer :: i, itest, l, ll, lll, llmax, llmin
    real(sp) :: absdel, alfa, falfe, fl
    complex(sp) :: cs, z

    del = (0.0_sp, 0.0_sp)
    alfa = dr * dr * temp / t0
    ! 0.166667 is not 1/6 in single precision, and the result depends on it
    alfa = 0.166667_sp * sqrt(alfa * alfa + dr0)
    falfe = -4.0_sp * alfa * e
    if (abs(falfe) - 1.0e-3_sp < 0.0_sp) then
      do i = 1, n3
        del(i) = cmplx(phs(i), 0.0_sp, sp)
      end do
      return
    end if

    ! bj(1:n1) holds the spherical Bessel functions of the first kind
    z = falfe * ci
    call spherical_bessel(bj, z, n1)

    cs = cmplx(1.0_sp, 0.0_sp, sp)
    fl = 1.0_sp
    do i = 1, n1
      bj(i) = exp(falfe) * fl * cs * bj(i)
      fl = fl + 2.0_sp
      cs = cs * ci
    end do

    fl = 1.0_sp
    do i = 1, n3
      ctab(i) = (exp(2.0_sp * phs(i) * ci) - (1.0_sp, 0.0_sp)) * fl
      fl = fl + 2.0_sp
    end do

    itest = 0
    do lll = 1, n2
      sum_(lll) = cmplx(0.0_sp, 0.0_sp, sp)
      do l = 1, n3
        llmin = abs(l - lll) + 1
        llmax = l + lll - 1
        do ll = llmin, llmax
          sum_(lll) = sum_(lll) + ppp(ll, lll, l) * ctab(l) * bj(ll)
        end do
      end do

      ! sum_(lll) is already the temperature-dependent t-matrix element. It
      ! is converted to a phase shift only to be converted back right after
      ! the call in atomic_tmatrix; kept for the sake of compatibility with
      ! the book of Van Hove and Tong.
      del(lll) = -ci * log(sum_(lll) + (1.0_sp, 0.0_sp)) / (2.0_sp, 0.0_sp)
      absdel = abs(del(lll))
      if (absdel - 1.0e-2_sp < 0.0_sp) then
        ! Small for the second time in a row: stop
        if (itest == 1) exit
        itest = 1
      else
        itest = 0
      end if
    end do
  end subroutine thermal_phase_shifts

  !> Spherical Bessel functions j_l(z), l = 0, ..., n1 - 1, of a complex
  !! argument, by direct evaluation of Pendry's expression (A.27) (legacy
  !! BESSEL). Downward recursion is said to be more stable; no deviation
  !! was found in the range of interest. V. Blum.
  !!
  !! n1 is the highest angular momentum, 2 * lmax + 1, for which bj is
  !! computed: for low lmax, the sums required for the temperature-dependent
  !! phase shifts do not converge fast enough with fewer.
  subroutine spherical_bessel(bj, z, n1)
    integer, intent(in) :: n1
    complex(dp), intent(out) :: bj(n1)
    complex(sp), intent(in) :: z

    complex(dp) :: pre, term, sum_, zsq
    integer :: k, l

    zsq = z * z / 2.0_sp

    ! pre is the l-dependent prefactor
    pre = 1.0_sp
    do l = 0, n1 - 1
      ! k indexes the series expansion, which ends when term is really small.
      ! The test is done in single precision.
      term = 1.0_sp
      sum_ = 1.0_sp
      k = 1
      do
        term = -term * (zsq / real(k * (2 * l + 2 * k + 1), sp))
        sum_ = sum_ + term
        k = k + 1
        if (.not. abs(cmplx(term, kind=sp)) > 1.0e-16_sp) exit
      end do
      bj(l + 1) = pre * sum_
      pre = pre * z / (2.0_sp * real(l + 1, sp) + 1.0_sp)
    end do
  end subroutine spherical_bessel
end module atomic_scattering
