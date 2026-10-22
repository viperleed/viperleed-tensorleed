!> \file rfacsb.f90
!! \brief Modules of the R-factor program rfactor, in dependency order.

! Circumvented allocation error in subroutine RINTAV by introduction
! of a second field ATAV and PQ1.        Oct. 2018 L. Hammer
! Changed the intensity limit of GRID and INTPOL to 1.E-8.
! Output for gnuplot adapted to particular computer settings. LH
! INPUT formats expanded to allow for 9999 exp. and theo. beams
!                                                    L.Hammer Mar. 2021
! Note: Y function output in YPEND and YSMOOTH is still limited to 999 beams
! Y_Smooth implementation (and nicer format of Y_Pendry) M. Schmid 2026
! Fortran 90 version: modules, explicit interfaces, arrays sized from the
! input.                                                              2026

!-----------------------------------------------------------------------
! R factors
!-----------------------------------------------------------------------
!
! Each experimental beam Ie(E) is compared with its theory beam It(E) over
! the energy interval common to both, without the gaps of the experiment.
! The theory is normalised to the experiment by c = int(Ie)/int(It); all
! integrals run over that interval (trapezoid rule on the working grid of
! step eincr), ' is d/dE:
!   R2   = 0.5 int((Ie - c*It)**2) / int(Ie**2)
!   R_ZJ = int(|Ie'' - c*It''| |Ie' - c*It'| / (|Ie'| + max|Ie'|))
!          / (0.027 int(Ie))
!          (reduced Zanazzi-Jona R factor; Zanazzi and Jona, Surf. Sci. 62,
!          61 (1977); integrand on a 10 times denser grid)
!   R_P  = int((Ye - Yt)**2) / int(Ye**2 + Yt**2)
!          (Pendry, J. Phys. C 13, 937 (1980)), with the Y function
!          Y = L/(1 + V0i**2 L**2) of the logarithmic derivative L = I'/I,
!          V0i the imaginary part of the inner potential (vi)
!   R_S  = as R_P with the smooth Y function Y_S (Imre et al., J. Phys.:
!          Condens. Matter 38, 105001 (2026)), which differs from Y_P
!          near intensity minima by a term in I'' (ysmooth)
! Averages over the beams: R2 and R_ZJ weighted by the beam weight wb
! times the length of the common interval; R_P and R_S as the ratio of the
! wb-weighted sums of numerators and denominators. The average R factor is
! the wr-weighted mean of the four, also per integer and fractional order
! beams (mittel).
!
! Inner potential: the theory energies are shifted rigidly by
! v0 = v01, v01 + vincr, ..., v02 (the real part of the inner potential
! relative to v0rr), and the common interval is found anew for each shift.
! The best fit is the lowest average R factor over all geometries and
! shifts.
!
! Beam averaging: experimental beams with the same label nbea are averaged
! (expav); theory beams with the same label kav are averaged with weights
! from their symmetry codes (rintav), e.g. for symmetry-equivalent beams
! or domains.
!
! Gaps: energy ranges of an experimental beam that are left out of the
! smoothing, the interpolation and the R-factor integrals; the theory
! integrals over the shifted gap are subtracted as well.
!
!-----------------------------------------------------------------------
! Legacy names (rfactor.f and rfacsb.f, fixed form)
!-----------------------------------------------------------------------
!
! Routines (same name unless given):
!   BINSRX, STFPTS, XNTERP, XNTRP2, YVAL, GRID, INTPOL     -> interpolation
!   COMMON /DATBLK/ (YVAL interval and points)  -> type interpolation_cache
!   DER, SMOOTH, YSMOOTH                                   -> derivatives
!   YPEND, VARSUM, EPSZJ; SUM -> trapezoid                 -> rfactor_terms
!   RINTAV, MAXINT                                      -> domain_averaging
!   READE, READT, GAPIN                                 -> experiment_input
!   EXPAV                                         -> experiment_preparation
!   RDPL -> read_theory_records                             -> theory_input
!        -> select_energies, average_theory (calls RINTAV, MAXINT)
!   STRIP2                                           -> theory_preparation
!   COMNEI                                               -> common_interval
!   COLUMN, SCRIPT; all WRITEs -> write_*                 -> rfactor_output
! Main program rfactor.f -> program rfactor (rfactor.f90) and
!   reading NL1-NL3, KAV, IBP, MITTEL                -> read_control
!   reading the experiment and gaps                  -> read_experiment
!   smoothing, interpolation, derivatives, Y functions and integrals of
!   the experiment / of one geometry -> prepare_experiment / prepare_theory
!   V0 loop -> scan_v0; beam loop -> beam_rfactor; averages over the
!   beams -> average_over_beams; best fit -> keep_best
!
! Variables of the main program (same name unless given):
!   NL1-NL3, KAV, IBP, MITTEL, NBEA, NINT, NFRAC
!                                           -> ctl%  (control_parameters)
!   theory TEXT, NBTD, PQ, SYM -> tr%title, nbtd, pq, sym (theory_records)
!   AE, EE, NEE, BENAME, AEP, AEPP, YPE, YSE, TSE, TSE2, TSEYP2, TSEYS2,
!     NBE, NGAP, GCOUNT, GBM, EG1, EG2, NG1, NG2  -> ex%  (experiment_beams)
!   AT, ET, NET, ATP, ATPP, YPT, YST, TST, TSTYP2, TSTYS2, NBT, PQAV
!                                                 -> th%  (theory_beams)
!   NEMIN, NEMAX, NTMIN, NTMAX, NGT1, NGT2, SEST, STST, RAVB
!                                                 -> res% (beam_results)
!   BRAV, NSB, BV0, BCIFE, BCIFT -> best%rav, geometry, v0, cife, cift
!   NEMINB, NEMAXB, NTMINB, NTMAXB, NG1B, NG2B, NGT1B, NGT2B, TSEB, TSTB,
!     BRAVB, BET, BAT -> best%  without the B (TSEB, TSTB -> sest, stst;
!     BRAVB -> ravb)                                    (best_fit)
!   AR, ARM -> sums%rf, rf_m; RAZP, RANP, RAZS, RANS -> sums%rp_num,
!     rp_den, rs_num, rs_den; RAZMP, RANMP, RAZMS, RANMS -> sums%rp_num_m,
!     rp_den_m, rs_num_m, rs_den_m; ERANG, ERANGM, SEM, STM, RAV, RAVM,
!     CIFE, CIFT -> sums%                               (rfactor_sums)
!   integrals SS, SU, SG, SE, ST, S2, SRZJ, SYP2, SYS2, ... and the beam
!     R factors R2, RRZJ, RPE, RSM: locals of beam_rfactor
!   V0, V0R, NV0, BGRAV: locals of scan_v0

!> Kind parameters for reals and complexes.
!!
!! sp is the precision of the legacy default REAL (REAL(4)); dp is double
!! precision (REAL(8)). Every real and complex is declared with one of these,
!! so that a change of precision is a change of one line here.
module kinds
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none

  private

  public :: sp, dp

  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
end module kinds

!> Enlarging allocatable arrays whose final size is known only after the
!! input has been read.
module array_growth
  use kinds, only: sp
  implicit none

  private

  public :: grow

  !> Makes the last dimension of an allocated array at least n long. The
  !! existing elements are kept; new elements are undefined.
  interface grow
    module procedure grow_integer, grow_real, grow_real_2d
  end interface grow

contains

  subroutine grow_integer(a, n)
    integer, allocatable, intent(inout) :: a(:)
    integer, intent(in) :: n

    integer, allocatable :: b(:)

    if (n <= size(a)) return
    allocate(b(n))
    b(1:size(a)) = a
    call move_alloc(b, a)
  end subroutine grow_integer

  subroutine grow_real(a, n)
    real(sp), allocatable, intent(inout) :: a(:)
    integer, intent(in) :: n

    real(sp), allocatable :: b(:)

    if (n <= size(a)) return
    allocate(b(n))
    b(1:size(a)) = a
    call move_alloc(b, a)
  end subroutine grow_real

  subroutine grow_real_2d(a, n)
    real(sp), allocatable, intent(inout) :: a(:, :)
    integer, intent(in) :: n

    real(sp), allocatable :: b(:, :)

    if (n <= size(a, 2)) return
    allocate(b(size(a, 1), n))
    b(:, 1:size(a, 2)) = a
    call move_alloc(b, a)
  end subroutine grow_real_2d
end module array_growth

!> Output of rfactor. The log on stdout (unit 6): echoes of the input,
!! intermediate curves depending on the print level ipr, and the R factors
!! of every beam, geometry and inner potential shift. ROUT (unit 7): the
!! R factors per beam and their averages; ROUTSHORT (unit 9): the best
!! average R factor of each geometry; both are opened by the main program.
!! The spectra of the best fit go to exp.column and theo.column, the plot
!! scripts to plot.gnu and plot.tex.
module rfactor_output
  use kinds, only: sp
  implicit none

  private

  public :: column, script
  public :: write_averaged_experiment, write_beam_counts, write_beam_list, &
            write_exp_beam_name, write_exp_curve, write_gap, &
            write_gap_count, write_nl1, write_pseudo_beam_count, &
            write_pseudo_beam_name, write_title
  public :: write_interpolated, write_pendry_y, write_smooth_y, &
            write_smoothed, write_too_few_points
  public :: write_rout_header, write_structure, write_suspect_intensity, &
            write_theory_header
  public :: write_averages, write_beam_header, write_beam_rfactors, &
            write_best_fit, write_no_overlap, write_normalization, &
            write_scan_end, write_v0_shift

contains

  !> Echoes the namelist NL1 as read.
  subroutine write_nl1(emin, emax, eincr, limfil, ipr, vi, v0rr, v01, v02, &
                       vincr, ismoth, eot, plot, gap)
    real(sp), intent(in) :: emin, emax, eincr, vi, v0rr, v01, v02, vincr, plot
    integer, intent(in) :: limfil, ipr, ismoth, eot, gap

    namelist /nl1/ emin, emax, eincr, limfil, ipr, vi, v0rr, v01, v02, &
                   vincr, ismoth, eot, plot, gap

    write(6, nl1)
  end subroutine write_nl1

  !> Echoes the number of theory (nbtd) and experimental (nbed) beams.
  subroutine write_beam_counts(nbtd, nbed)
    integer, intent(in) :: nbtd, nbed

3   format('*NBTD,NBED = ', 2I4)

    write(6, 3) nbtd, nbed
  end subroutine write_beam_counts

  !> Echoes a list with one entry per beam under a heading, 25 per line.
  subroutine write_beam_list(heading, list)
    character(len=*), intent(in) :: heading
    integer, intent(in) :: list(:)

10  format(A, /, 400(25I4, /))

    write(6, 10) heading, list
  end subroutine write_beam_list

  !> Echoes a title (of the experiment or the theory).
  subroutine write_title(title)
    character(len=80), intent(in) :: title

23  format('*', A80)

    write(6, 23) title
  end subroutine write_title

  !> Echoes the name of experimental beam ib.
  subroutine write_exp_beam_name(ib, bename)
    integer, intent(in) :: ib
    character(len=4), intent(in) :: bename(5)

70  format('*EXP. BEAM ', 1I4, ' (', 5A4, ')')

    write(6, 70) ib, bename
  end subroutine write_exp_beam_name

  !> Echoes the energies ee and intensities ae of an experimental beam.
  subroutine write_exp_curve(ee, ae)
    real(sp), intent(in) :: ee(:), ae(:)

    integer :: ie

80  format('* EXP. ENERGIES AND INTENSITIES', /, 100(5(1F7.2, 1E13.4, 3X), /))

    write(6, 80) (ee(ie), ae(ie), ie = 1, size(ee))
  end subroutine write_exp_curve

  !> Echoes the number of beams of a pseudo-experiment.
  subroutine write_pseudo_beam_count(nbeams)
    integer, intent(in) :: nbeams

103 format(I5, E13.4)

    write(6, 103) nbeams
  end subroutine write_pseudo_beam_count

  !> Echoes the label inn and name of a beam of a pseudo-experiment.
  subroutine write_pseudo_beam_name(inn, bename)
    integer, intent(in) :: inn
    character(len=4), intent(in) :: bename(5)

101 format(I5, 22A4)

    write(6, 101) inn, bename
  end subroutine write_pseudo_beam_name

  !> Echoes the number of gaps in the experimental spectra.
  subroutine write_gap_count(ngap)
    integer, intent(in) :: ngap

20  format('*THERE ARE', I3, ' GAPS IN THE EXPERIMENTAL SPECTRA')

    write(6, 20) ngap
  end subroutine write_gap_count

  !> Echoes gap ig: beam bm, energy range eg1 to eg2.
  subroutine write_gap(ig, bm, eg1, eg2)
    integer, intent(in) :: ig, bm
    real(sp), intent(in) :: eg1, eg2

50  format(' GAP', I3, ' IN BEAM', I3, ' BETWEEN', F7.2, ' AND', F7.2, ' EV')

    write(6, 50) ig, bm, eg1, eg2
  end subroutine write_gap

  !> Prints the averaged experimental beams 1..nbe (nee points each on
  !! the energies ee, intensities ae).
  subroutine write_averaged_experiment(nbe, nee, ee, ae)
    integer, intent(in) :: nbe, nee(:)
    real(sp), intent(in) :: ee(:, :), ae(:, :)

    integer :: ib, ie, n

10  format('*EXP. ENERG. AND INTENS. AFTER AVERAGING IN BEAM', 1I4, /, &
           100(5(1F7.2, 1E13.4, 3X), /))

    do ib = 1, nbe
      n = nee(ib)
      write(6, 10) ib, (ee(ib, ie), ae(ib, ie), ie = 1, n)
    end do
  end subroutine write_averaged_experiment

  !> Prints experimental beam ib (energies ee, intensities ae) after
  !! smoothing.
  subroutine write_smoothed(ib, ee, ae)
    integer, intent(in) :: ib
    real(sp), intent(in) :: ee(:), ae(:)

    integer :: ie

10  format('*EXP. ENERG. AND INTENS. AFTER SMOOTHING IN BEAM', 1I3, &
           /, 2000(5(1F7.2, 1E13.4, 3X), /))

    write(6, 10) ib, (ee(ie), ae(ie), ie = 1, size(ee))
  end subroutine write_smoothed

  !> Reports a beam with too few points to be interpolated.
  subroutine write_too_few_points()
10  format('0 IN THIS BEAM TOO FEW ENERGY VALUES FOR INTERPOLATION')

    write(6, 10)
  end subroutine write_too_few_points

  !> Prints beam ib (energies e, intensities a) after interpolation onto
  !! the working grid.
  subroutine write_interpolated(ib, e, a)
    integer, intent(in) :: ib
    real(sp), intent(in) :: e(:), a(:)

    integer :: ie

20  format('*INTENSITIES AFTER INTERPOLATION IN BEAM', 1I4, &
           /, 100(5(1F7.2, 1E13.4, 3X), /))

    write(6, 20) ib, (e(ie), a(ie), ie = 1, size(e))
  end subroutine write_interpolated

  !> Prints the Pendry Y function y(ib, 1:ne(ib)) on the energies e of the
  !! beams 1..nb that have points.
  subroutine write_pendry_y(nb, ne, e, y)
    integer, intent(in) :: nb, ne(:)
    real(sp), intent(in) :: e(:, :), y(:, :)

    integer :: ib, ie, n

40  format(' PENDRY Y FUNCTION IN BEAM', 1I3, /, 1000(5(1F7.2, 1E13.4, 3X), /))

    do ib = 1, nb
      n = ne(ib)
      if (n /= 0) then
        write(6, 40) ib, (e(ib, ie), y(ib, ie), ie = 1, n)
      end if
    end do
  end subroutine write_pendry_y

  !> Prints the smooth Y function y(ib, 1:ne(ib)) on the energies e of the
  !! beams 1..nb that have points.
  subroutine write_smooth_y(nb, ne, e, y)
    integer, intent(in) :: nb, ne(:)
    real(sp), intent(in) :: e(:, :), y(:, :)

    integer :: ib, ie, n

40  format(' SMOOTH Y FUNCTION IN BEAM', 1I5, /, 1000(5(1F7.2, 1E13.4, 3X), /))

    do ib = 1, nb
      n = ne(ib)
      if (n /= 0) then
        write(6, 40) ib, (e(ib, ie), y(ib, ie), ie = 1, n)
      end if
    end do
  end subroutine write_smooth_y

  !> Writes the theory title to ROUT and stdout, and echoes the beams
  !! (parallel momenta pq) with their groups kav and the code of every
  !! geometry (values).
  subroutine write_theory_header(title, pq, kav, values)
    character(len=80), intent(in) :: title
    real(sp), intent(in) :: pq(:, :), values(:)
    integer, intent(in) :: kav(:)

    integer :: i, j

22  format(A80)
23  format('*', A80)
27  format('*THEOR. BEAM GROUPINGS')
30  format(' ', 2F10.5, I5)
65  format('*SURFACE STRUCTURE NO. AND CHARACTERIZATION')
75  format(18X, I3, F15.4)

    write(7, 22) title
    write(6, 23) title

    write(6, 27)
    do j = 1, size(kav)
      write(6, 30) (pq(i, j), i = 1, 2), kav(j)
    end do
    write(6, 65)
    do j = 1, size(values)
      write(6, 75) j, values(j)
    end do
  end subroutine write_theory_header

  !> Reports an averaged theory beam ib whose maximum intensity am
  !! exceeds 1.
  subroutine write_suspect_intensity(ib, am)
    integer, intent(in) :: ib
    real(sp), intent(in) :: am

28  format(' ', //, ' MAX INTENSITY IN THEOR BEAM ', 1I3, ' IS SUSPECT(', &
           1E13.5, ')- **** STOP PROGRAM ****')

    write(6, 28) ib, am
  end subroutine write_suspect_intensity

  !> Writes the column headings of ROUT.
  subroutine write_rout_header()
29  format(1X, /, 6X, 'BEAM', 12X, 'IBE', 3X, 'D12', 5X, 'V0R', 4X, 'EMIN', &
           4X, 'EMAX', 5X, 'OVL', 5X, 'RAV')

    write(7, 29)
  end subroutine write_rout_header

  !> Starts the log of geometry is with code value.
  subroutine write_structure(is, value)
    integer, intent(in) :: is
    real(sp), intent(in) :: value

30  format('*STRUCTURE NO.', 1I5, ' CHARACTERIZED BY VALUE', 1F7.4)

    write(6, 30) is, value
  end subroutine write_structure

  !> Starts the log of inner potential shift v0.
  subroutine write_v0_shift(v0)
    real(sp), intent(in) :: v0

60  format('*THEOR. INNER POTENTIAL SHIFTED BY', 1F7.2, ' EV')

    write(6, 60) v0
  end subroutine write_v0_shift

  !> Starts the log of experimental beam ibe (name bename) and its theory
  !! beam ibt (parallel momentum pq).
  subroutine write_beam_header(ibe, bename, ibt, pq)
    integer, intent(in) :: ibe, ibt
    character(len=4), intent(in) :: bename(5)
    real(sp), intent(in) :: pq(2)

70  format('*BEAM NO.', 1I3, ' IN EXP. (', 5A4, '), WHICH IS BEAM NO.', 1I3, &
           ' IN THEORY (', 2F6.3, ')')

    write(6, 70) ibe, bename, ibt, pq
  end subroutine write_beam_header

  !> Reports experimental beam ibe (name bename) without overlap with its
  !! theory beam, and writes its ROUT line for geometry code d12, inner
  !! potential v0r and energy range emin to emax.
  subroutine write_no_overlap(ibe, bename, d12, v0r, emin, emax)
    integer, intent(in) :: ibe
    character(len=4), intent(in) :: bename(5)
    real(sp), intent(in) :: d12, v0r, emin, emax

75  format('  **** NO OVERLAP THEOR./EXP. IN (EXP.) BEAM NO.', 1I5)
80  format(5A4, I5, 1F8.4, 3F8.2, 1F10.2, 1F8.4)

    write(6, 75) ibe
    write(7, 80) bename, ibe, d12, v0r, emin, emax, 0.0_sp, 0.0_sp
  end subroutine write_no_overlap

  !> Logs the normalisation c of the experiment to the theory.
  subroutine write_normalization(c)
    real(sp), intent(in) :: c

112 format(' NORMALIZATION C EXP/THEOR =', 1E13.5)

    write(6, 112) c
  end subroutine write_normalization

  !> Logs the R factors rf (R2, R_ZJ, R_P, R_S) and their weighted average
  !! ravb of experimental beam ibe (name bename) at shift v0 over the common
  !! energy range eet, and writes its ROUT line for geometry code d12,
  !! inner potential v0r and energy range emin to emax.
  subroutine write_beam_rfactors(v0, eet, rf, ravb, ibe, bename, d12, v0r, &
                                 emin, emax)
    real(sp), intent(in) :: v0, eet, rf(:), ravb, d12, v0r, emin, emax
    integer, intent(in) :: ibe
    character(len=4), intent(in) :: bename(5)

80  format(5A4, I5, 1F8.4, 3F8.2, 1F10.2, 1F8.4)
115 format(a)
120 format(' ', 1F7.2, 1F10.4, 5F7.4)

    write(6, 115) '   V0        EET    R_2    R_ZJ   R_P    R_S     RAV'
    write(6, 120) v0, eet, rf, ravb
    write(7, 80) bename, ibe, d12, v0r, emin, emax, eet, ravb
  end subroutine write_beam_rfactors

  !> Logs the R factors at shift v0 averaged over the integer order beams
  !! (if nint > 0), the fractional order beams (if nfrac > 0) and all
  !! beams, with their energy ranges and weighted averages, and the ratios
  !! cife and cift of fractional to integer order intensities (if both
  !! exist); writes the averages to ROUT for geometry code d12, inner
  !! potential v0r and energy range emin to emax. arm, erangm and ravm hold
  !! the integer (1) and fractional (2) order averages.
  subroutine write_averages(v0, d12, v0r, emin, emax, nint, nfrac, erangm, &
                            arm, ravm, erang, ar, rav, cife, cift)
    real(sp), intent(in) :: v0, d12, v0r, emin, emax
    integer, intent(in) :: nint, nfrac
    real(sp), intent(in) :: erangm(2), arm(:, :), ravm(2), erang, ar(:), &
                            rav, cife, cift

120 format(' ', 1F7.2, 1F10.4, 5F7.4)
125 format('AV.-INT             ', I5, F8.4, 3F8.2, 1F10.2, F8.4, '  <---')
126 format('AV.-FRAC ORD        ', I5, F8.4, 3F8.2, 1F10.2, F8.4, '  <---')
127 format('AVERAGE             ', I5, F8.4, 3F8.2, 1F10.2, F8.4, '  <---')
128 format('RELATION FRAC TO INT BEAM INTENSITIES:', F6.3, &
           ' FOR EXP DATA,', F6.3, ' FOR THEOR DATA')
140 format('*AVERAGE OVER ALL BEAMS')

    write(6, 140)
    if (nint /= 0) then
      write(6, 120) v0, erangm(1), arm(:, 1), ravm(1)
      write(7, 125) -1, d12, v0r, emin, emax, erangm(1), ravm(1)
    end if
    if (nfrac /= 0) then
      write(6, 120) v0, erangm(2), arm(:, 2), ravm(2)
      write(7, 126) -1, d12, v0r, emin, emax, erangm(2), ravm(2)
    end if
    write(6, 120) v0, erang, ar, rav
    write(7, 127) 0, d12, v0r, emin, emax, erang, rav
    if (nfrac /= 0 .and. nint /= 0) then
      write(6, 128) cife, cift
    end if
  end subroutine write_averages

  !> Ends the scan of one geometry: a blank line in ROUT and its best
  !! average R factor bgrav in ROUTSHORT.
  subroutine write_scan_end(bgrav)
    real(sp), intent(in) :: bgrav

129 format(' ')
130 format(F8.4)

    write(7, 129)
    write(9, 130) bgrav
  end subroutine write_scan_end

  !> Writes the best fit, geometry nsb (code value) at shift bv0 with the
  !! average R factor brav, to stdout and ROUT, the ratios bcife and bcift
  !! of fractional to integer order intensities to ROUT (if nint and nfrac
  !! are both nonzero), and ends the log.
  subroutine write_best_fit(nsb, value, bv0, brav, nint, nfrac, bcife, bcift)
    integer, intent(in) :: nsb, nint, nfrac
    real(sp), intent(in) :: value, bv0, brav, bcife, bcift

128 format('RELATION FRAC TO INT BEAM INTENSITIES:', F6.3, &
           ' FOR EXP DATA,', F6.3, ' FOR THEOR DATA')
175 format('BEST GEOMETRY ENCOUNTERED, NO.', 1I5, &
           ', CHARACTERIZED BY VALUE ', 1F7.4, /, &
           'WITH INNER POTENTIAL SHIFT ', 1F7.2, &
           ' EV, AVERAGE R-FACTOR = ', 1F7.4)
180 format('*CORRECT TERMINATION')

    write(6, 175) nsb, value, bv0, brav
    write(7, 175) nsb, value, bv0, brav
    if (nfrac /= 0 .and. nint /= 0) then
      write(7, 128) bcife, bcift
    end if
    write(6, 180)
  end subroutine write_best_fit

  !> Writes the experimental spectra to exp.column and the best-fit
  !! theoretical spectra to theo.column, one line per energy step with an
  !! (energy, intensity) pair per beam. Energies outside a beam's range or
  !! inside one of its gcount(ibe) gaps (indices nge1:nge2 resp. nget1:nget2,
  !! numbered over all beams) are written as zero intensity. Intensities are
  !! scaled to intmax: with norm = 1 after division by the integrals tse
  !! resp. tst, to the largest value of both curves; with norm = 2 each
  !! curve to its own maximum. ngap is the size of the gap arrays. The gap
  !! tests look up the gap after the last one of a beam before checking that
  !! there is one; the per-beam gap indices nge and ngt have a zero entry
  !! for it.
  subroutine column(neminb, nemaxb, ntminb, ntmaxb, ee, ae, bet, bat, nbed, &
                    nbtd, ngp, nbe, ibp, gc, nge1, nge2, ngt1, ngt2, ngap, &
                    tse, tst, v0, intmax, norm)
    integer, intent(in) :: nbed, nbtd, ngp, nbe, ngap, norm
    integer, intent(in) :: neminb(nbed), nemaxb(nbed)
    integer, intent(in) :: ntminb(nbtd), ntmaxb(nbtd)
    real(sp), intent(in) :: ee(nbed, ngp), ae(nbed, ngp)
    real(sp), intent(in) :: bet(nbtd, ngp), bat(nbtd, ngp)
    integer, intent(in) :: ibp(nbed), gc(nbed)
    integer, intent(in) :: nge1(:), nge2(:), ngt1(:), ngt2(:)
    real(sp), intent(in) :: tse(nbed), tst(nbtd), v0, intmax

    real(sp) :: aei(nbed), eei(nbed), bati(nbed), beti(nbed)
    real(sp) :: batint(nbed, ngp), betint(nbed, ngp), tstint(nbed)
    real(sp) :: amax(nbed), amaxt(nbed)
    integer :: ntmini(nbed), ntmaxi(nbed)
    integer :: igc(nbed), nge(nbed, ngap + 1, 2), ngt(nbed, ngap + 1, 2)
    integer :: neemin, neemax, netmin, netmax
    integer :: i, ibe, iee, iet, igap, igdum, igp
    character(len=10) :: fmt

10  format('(', I4, 'F8.2)')

    ! spectra output format for exp.column and theo.column
    write(fmt, 10) 2*nbe

    open(10, file='exp.column')
    open(11, file='theo.column')

    ! lowest lower and highest upper energy indices of all beams
    neemin = ngp
    neemax = 1
    netmin = ngp
    netmax = 1

    do ibe = 1, nbe
      if (ibp(ibe) == 0) then
        ntmini(ibe) = ngp
        ntmaxi(ibe) = 1
        tstint(ibe) = 1._sp
        do igp = 1, ngp
          betint(ibe, igp) = 0._sp
          batint(ibe, igp) = 0._sp
        end do
      else
        ntmini(ibe) = ntminb(ibp(ibe))
        ntmaxi(ibe) = ntmaxb(ibp(ibe))
        tstint(ibe) = tst(ibp(ibe))
        do igp = 1, ngp
          betint(ibe, igp) = bet(ibp(ibe), igp)
          batint(ibe, igp) = bat(ibp(ibe), igp)
        end do
      end if
    end do

    do ibe = 1, nbe
      amax(ibe) = 0._sp
      if (neminb(ibe) < neemin) then
        neemin = neminb(ibe)
      end if
      if (nemaxb(ibe) > neemax) then
        neemax = nemaxb(ibe)
      end if
      amaxt(ibe) = 0._sp
      if (ntmini(ibe) < netmin) then
        netmin = ntmini(ibe)
      end if
      if (ntmaxi(ibe) > netmax) then
        netmax = ntmaxi(ibe)
      end if
    end do

    ! gap indices per beam
    do ibe = 1, nbe
      do i = 1, ngap + 1
        nge(ibe, i, 1) = 0
        nge(ibe, i, 2) = 0
        ngt(ibe, i, 1) = 0
        ngt(ibe, i, 2) = 0
      end do
    end do

    igap = 0
    do ibe = 1, nbe
      do i = 1, gc(ibe)
        nge(ibe, i, 1) = nge1(i + igap)
        nge(ibe, i, 2) = nge2(i + igap)
        ngt(ibe, i, 1) = ngt1(i + igap)
        ngt(ibe, i, 2) = ngt2(i + igap)
      end do
      igap = igap + gc(ibe)
    end do

    ! highest intensity amax(ibe) of each experimental beam
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iee = neemin, neemax
      do ibe = 1, nbe
        if (iee >= neminb(ibe) .and. iee <= nemaxb(ibe)) then
          igdum = igc(ibe) + 1
          if (gc(ibe) /= 0 .and. iee >= nge(ibe, igdum, 1) .and. &
              iee <= nge(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
            aei(ibe) = 0._sp
            if (iee == nge(ibe, igdum, 2)) then
              igc(ibe) = igdum
            end if
          else if (norm == 1) then
            aei(ibe) = ae(ibe, iee)/tse(ibe)
          else if (norm == 2) then
            aei(ibe) = ae(ibe, iee)
          end if

          if (aei(ibe) > amax(ibe)) then
            amax(ibe) = aei(ibe)
          end if
        end if
      end do
    end do

    ! highest intensity of each theoretical beam: into amax (norm = 1) or
    ! amaxt (norm = 2)
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iet = netmin, netmax
      do ibe = 1, nbe
        if (iet >= ntmini(ibe) .and. iet <= ntmaxi(ibe)) then
          igdum = igc(ibe) + 1
          if (gc(ibe) /= 0 .and. iet >= ngt(ibe, igdum, 1) .and. &
              iet <= ngt(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
            bati(ibe) = 0._sp
            if (iet == ngt(ibe, igdum, 2)) then
              igc(ibe) = igdum
            end if
          else if (norm == 1) then
            bati(ibe) = batint(ibe, iet)/tstint(ibe)
            if (bati(ibe) > amax(ibe)) then
              amax(ibe) = bati(ibe)
            end if
          else if (norm == 2) then
            bati(ibe) = batint(ibe, iet)
            if (bati(ibe) > amaxt(ibe)) then
              amaxt(ibe) = bati(ibe)
            end if
          end if
        end if
      end do
    end do

    ! normalised experimental beams to exp.column
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iee = neemin, neemax
      do ibe = 1, nbe
        igdum = igc(ibe) + 1
        if (iee < neminb(ibe) .or. iee > nemaxb(ibe)) then
          eei(ibe) = 0._sp
          aei(ibe) = 0._sp
        else if (gc(ibe) /= 0 .and. iee >= nge(ibe, igdum, 1) .and. &
                 iee <= nge(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
          eei(ibe) = 0._sp
          aei(ibe) = 0._sp
          if (iee == nge(ibe, igdum, 2)) then
            igc(ibe) = igdum
          end if
        else
          eei(ibe) = ee(ibe, iee)
          if (norm == 1) then
            aei(ibe) = ae(ibe, iee)*intmax/(tse(ibe)*amax(ibe))
          else if (norm == 2) then
            aei(ibe) = ae(ibe, iee)*intmax/amax(ibe)
          end if
        end if
      end do

      write(10, fmt) (eei(ibe), aei(ibe), ibe = 1, nbe)
    end do

    ! normalised theoretical beams to theo.column
    do ibe = 1, nbe
      igc(ibe) = 0
    end do

    do iet = netmin, netmax
      do ibe = 1, nbe
        igdum = igc(ibe) + 1
        if (iet < ntmini(ibe) .or. iet > ntmaxi(ibe)) then
          beti(ibe) = -v0
          bati(ibe) = 0._sp
        else if (gc(ibe) /= 0 .and. iet >= ngt(ibe, igdum, 1) .and. &
                 iet <= ngt(ibe, igdum, 2) .and. igc(ibe) < gc(ibe)) then
          beti(ibe) = -v0
          bati(ibe) = 0._sp
          if (iet == ngt(ibe, igdum, 2)) then
            igc(ibe) = igdum
          end if
        else
          beti(ibe) = betint(ibe, iet)
          if (norm == 1) then
            bati(ibe) = batint(ibe, iet)*intmax/(tstint(ibe)*amax(ibe))
          else if (norm == 2) then
            bati(ibe) = batint(ibe, iet)*intmax/amaxt(ibe)
          end if
        end if
      end do

      write(11, fmt) (beti(ibe) + v0, bati(ibe), ibe = 1, nbe)
    end do
  end subroutine column

  !> Writes plot.gnu, a gnuplot script that plots the spectra of
  !! exp.column and theo.column to one PDF per beam (energy range emin to
  !! emax, labelled with the beam name and its R factor rb), and plot.tex,
  !! a LaTeX document with one page per plot under the title and the
  !! average R factor brav.
  subroutine script(emin, emax, intmax, plsize, xtics, nbe, bename, rb, &
                    nbed, title, brav, ibp)
    integer, intent(in) :: nbe, nbed, xtics
    real(sp), intent(in) :: emin, emax, intmax, plsize(2), rb(nbed), brav
    character(len=4), intent(in) :: bename(5, nbed)
    character(len=80), intent(in) :: title
    integer, intent(in) :: ibp(nbed)

    integer :: i, ibe

10  format('set terminal pdfcairo enhanced font "Arial,12"', /, &
           'set encoding iso_8859_1')
11  format('set xrange [', I3, ':', I3, ']', /, 'set yrange [0:', I4, ']')
12  format('set size ', F5.3, ',', F5.3, /, 'set style data lines', /, &
           'set xtics ', I3, /, 'set noytics', /, &
           'set xlabel "Energy [eV]"', /, 'set noxzeroaxis', /)
15  format('set output "pic_', I3.3, '.pdf"', /, 'set nolabel')
16  format('set label "', 5A4, '\nR =', F7.4, '" at ', I3, &
           ',', I4, ' center font "Helvetica,12"')
17  format('plot ''exp.column'' using ($', I3, '>0?$', I3, ':1/0) :', &
           '($', I3, '>0?$', I3, ':1/0) title "(exp.)" ', &
           'lt rgb "red" lw 2,\', / &
           ' ''theo.column'' using ($', I3, '>0?$', I3, ':1/0) :', &
           '($', I3, '>0?$', I3, ':1/0) title "(theory)" ', &
           'lt rgb "blue" lw 2', /)
20  format('\documentclass[12pt,a4paper]{article}', /, &
           '\usepackage{graphicx}', /, '\usepackage{float}', /)
21  format('\begin{document}', /, '\pagestyle{empty}', /)
22  format('\newcommand{\bigtitle}{', /, '\begin{center}', /, '\large', &
           /, '\textbf{', A80, '}\\', /, '$\mathbf{R_{P} =', F7.4, '}$', &
           /, '\end{center}', /, '}', /)
25  format('\newpage', /, '\bigtitle', /, '\begin{figure}[H]')
26  format('\includegraphics[angle=0,width=14cm]{pic_', I3.3, &
           '.pdf}', /, '\end{figure}', /)
27  format('\end{document}')

    open(12, file='plot.gnu')
    open(13, file='plot.tex')

    ! general settings of the plots
    write(12, 10)
    write(12, 11) int(emin), int(emax), int(1.5_sp*intmax)
    write(12, 12) plsize(1), plsize(2), xtics

    ! one plot per beam
    do ibe = 1, nbe
      if (ibp(ibe) /= 0) then
        write(12, 15) ibe
        write(12, 16) (bename(i, ibe), i = 1, 5), rb(ibe), &
          int(emin + (emax - emin)/2), int(1.2_sp*intmax)
        write(12, 17) (2*ibe) - 1, (2*ibe) - 1, (2*ibe) - 1, 2*ibe, &
          (2*ibe) - 1, (2*ibe) - 1, (2*ibe) - 1, 2*ibe
      end if
    end do

    ! document settings and title
    write(13, 20)
    write(13, 21)
    write(13, 22) title, brav

    ! one page per plot
    do ibe = 1, nbe
      if (ibp(ibe) /= 0) then
        write(13, 25)
        write(13, 26) ibe
      end if
    end do

    write(13, 27)
  end subroutine script
end module rfactor_output

!> Interpolation of intensity curves and the energy working grid.
!!
!! yval interpolates a tabulated curve at one abscissa: 3rd-order polynomial
!! interpolation through the four nearest points, or 2nd/1st order if the
!! curve has fewer than four points. The four points of the last interval
!! found are kept in an interpolation_cache and reused when the next call
!! with the same cache finds the same interval, even if that call passes a
!! different curve.
module interpolation
  use kinds, only: sp
  implicit none

  private

  public :: grid, grid_points, intpol, interpolation_cache, yval

  !> Interval (il, ih) of the last yval call and the four points (x, y)
  !! stored for it. il = ih = 0 forces yval to store new points.
  type :: interpolation_cache
    integer :: il, ih
    real(sp) :: x(0:3), y(0:3)
  end type interpolation_cache

contains

  !> Finds the interpolation interval [worx(il), worx(ih)] containing x by
  !! successive halving. Outside the table, the first or last interval.
  subroutine binsrx(il, ih, x, worx, length)
    integer, intent(out) :: il, ih
    real(sp), intent(in) :: x
    integer, intent(in) :: length
    real(sp), intent(in) :: worx(length)

    integer :: i, ihi

    i = (length + 1) / 2
    ihi = (i + 1) / 2
    if (x < worx(1) .or. x >= worx(length)) then
      if (x < worx(1)) then
        il = 1
        ih = 2
      else
        il = length - 1
        ih = length
      end if
      return
    end if
    do
      if (i < 1) then
        i = ihi + i
      else if (i >= length) then
        i = i - ihi
      else if (x >= worx(i) .and. x <= worx(i + 1)) then
        il = i
        ih = i + 1
        return
      else if (x < worx(i)) then
        i = i - ihi
      else
        i = ihi + i
      end if
      ihi = (ihi + 1) / 2
    end do
  end subroutine binsrx

  !> Stores the four grid points nearest to the interval (il, ih) and their
  !! ordinates in cache, shifted inwards at the ends of the table.
  subroutine stfpts(il, ih, worx, wory, length, cache)
    integer, intent(in) :: il, ih, length
    real(sp), intent(in) :: worx(length), wory(length)
    type(interpolation_cache), intent(inout) :: cache

    integer :: i

    i = il - 1
    if (il <= 1) i = i + 1
    if (ih >= length) i = i - 1
    cache%x(0:3) = worx(i:i + 3)
    cache%y(0:3) = wory(i:i + 3)
  end subroutine stfpts

  !> 3rd-order (Newton) polynomial through the four points in cache,
  !! evaluated at x.
  real(sp) function xnterp(x, cache)
    real(sp), intent(in) :: x
    type(interpolation_cache), intent(in) :: cache

    real(sp) :: term, fact1, fact2, temp

    associate (x0 => cache%x(0), x1 => cache%x(1), x2 => cache%x(2), &
               x3 => cache%x(3), y0 => cache%y(0), y1 => cache%y(1), &
               y2 => cache%y(2), y3 => cache%y(3))
      term = y0
      fact1 = x - x0
      fact2 = (y1 - y0) / (x1 - x0)
      term = term + fact1 * fact2
      fact1 = fact1 * (x - x1)
      fact2 = ((y2 - y1)/(x2 - x1) - fact2) / (x2 - x0)
      term = term + fact1 * fact2
      fact1 = fact1 * (x - x2)
      temp = ((y3 - y2)/(x3 - x2) - (y2 - y1)/(x2 - x1))/(x3 - x1)
      fact2 = (temp - fact2) / (x3 - x0)
      term = term + fact1 * fact2
    end associate
    xnterp = term
  end function xnterp

  !> 2nd-order polynomial through the first three points (xs, ys), or
  !! 1st order through the first two if n <= 2, evaluated at x.
  subroutine xntrp2(x, y, xs, ys, n)
    real(sp), intent(in) :: x
    real(sp), intent(out) :: y
    integer, intent(in) :: n
    real(sp), intent(in) :: xs(n), ys(n)

    real(sp) :: a, b, c

    if (n <= 2) then
      y = (ys(2) - ys(1))*(x - xs(1))/(xs(2) - xs(1)) + ys(1)
      return
    end if
    a = (ys(1) - ys(2))/(xs(1) - xs(2))/(xs(2) - xs(3)) - &
        (ys(1) - ys(3))/(xs(1) - xs(3))/(xs(2) - xs(3))
    b = (ys(1) - ys(2))/(xs(1) - xs(2)) - a*(xs(1) + xs(2))
    c = ys(1) - (a*xs(1) + b)*xs(1)
    y = c + x*(b + a*x)
  end subroutine xntrp2

  !> Interpolates the curve (worx, wory) at x. cache holds the interval and
  !! points of the previous call (cache%il = cache%ih = 0 to force new ones);
  !! for four or more points it is updated.
  real(sp) function yval(x, wory, worx, length, cache)
    real(sp), intent(in) :: x
    integer, intent(in) :: length
    real(sp), intent(in) :: wory(length), worx(length)
    type(interpolation_cache), intent(inout) :: cache

    real(sp) :: y
    integer :: il, ih

    if (length < 4) then
      call xntrp2(x, y, worx, wory, length)
      yval = y
      return
    end if
    call binsrx(il, ih, x, worx, length)
    if (il /= cache%il .or. ih /= cache%ih) then
      call stfpts(il, ih, worx, wory, length, cache)
    end if
    y = xnterp(x, cache)
    cache%ih = ih
    cache%il = il
    yval = y
  end function yval

  !> Energy working grid (steps of eincr) for one beam. Drops the leading
  !! points with zero intensity from (a, e) and shortens ne accordingly;
  !! ne = 0 if at most one point is left. egrid(1:npts) is the grid and
  !! imin the index of the first nonzero intensity on entry.
  subroutine grid(a, ne, e, ngp, eincr, egrid, npts, imin, too_few)
    integer, intent(in) :: ngp
    real(sp), intent(inout) :: a(ngp), e(ngp)
    integer, intent(inout) :: ne
    real(sp), intent(in) :: eincr
    real(sp), intent(inout) :: egrid(ngp)
    integer, intent(inout) :: npts, imin
    logical, intent(out) :: too_few

    real(sp) :: xmin
    integer :: lmin, lmax, i, ie

    do ie = 1, ne
      imin = ie
      if (a(ie) > 1.e-8_sp) exit
    end do

    too_few = imin == ne
    if (too_few) then
      ne = 0
      return
    end if

    lmin = int((e(imin) - 0.0001_sp) / eincr) + 1
    lmin = max(lmin, 0)
    xmin = real(lmin, sp) * eincr
    lmax = int((e(ne) + 0.0001_sp) / eincr)

    do i = imin, ne
      e(i - imin + 1) = e(i)
      a(i - imin + 1) = a(i)
    end do
    ne = ne - imin + 1

    npts = lmax - lmin + 1

    egrid(1) = xmin
    do ie = 2, npts
      egrid(ie) = egrid(ie - 1) + eincr
    end do
  end subroutine grid

  !> Upper bound of the number of working grid points (steps of eincr)
  !! that grid and intpol make from a curve on the energies e.
  pure integer function grid_points(e, eincr)
    real(sp), intent(in) :: e(:), eincr

    integer :: lmin, lmax

    if (size(e) < 2) then
      grid_points = 0
      return
    end if
    lmin = int((minval(e) - 0.0001_sp) / eincr) + 1
    lmin = max(lmin, 0)
    lmax = int((maxval(e) + 0.0001_sp) / eincr)
    grid_points = lmax - lmin + 1
  end function grid_points

  !> Interpolates the intensities a(e) onto the working grid (steps of
  !! eincr), replacing (a, e, ne); negative values become zero. Leading
  !! zero intensities are dropped; ne = 0 if at most one point is left
  !! (too_few). x and woryt are work space.
  subroutine intpol(a, ne, e, ngp, eincr, x, woryt, too_few)
    integer, intent(in) :: ngp
    real(sp), intent(inout) :: a(ngp), e(ngp)
    integer, intent(inout) :: ne
    real(sp), intent(in) :: eincr
    real(sp), intent(inout) :: x(ngp), woryt(ngp)
    logical, intent(out) :: too_few

    type(interpolation_cache) :: cache
    real(sp) :: xmin, xval
    integer :: imin, lmin, lmax, nem, npts, i, ie

    cache%il = 0
    cache%ih = 0

    do ie = 1, ne
      imin = ie
      if (a(ie) > 1.e-8_sp) exit
    end do

    too_few = imin == ne
    if (too_few) then
      ne = 0
      return
    end if

    lmin = int((e(imin) - 0.0001_sp) / eincr) + 1
    lmin = max(lmin, 0)
    xmin = real(lmin, sp) * eincr
    lmax = int((e(ne) + 0.0001_sp) / eincr)

    do i = imin, ne
      x(i - imin + 1) = e(i)
      woryt(i - imin + 1) = a(i)
    end do
    nem = ne - imin + 1

    npts = lmax - lmin + 1
    ne = npts

    xval = xmin - eincr

    do i = 1, npts
      xval = xval + eincr
      a(i) = yval(xval, woryt(1:nem), x(1:nem), nem, cache)
      if (a(i) < 0.0_sp) then
        a(i) = 0.0_sp
      end if
    end do

    e(1) = xmin
    do ie = 2, npts
      e(ie) = e(ie - 1) + eincr
    end do
  end subroutine intpol
end module interpolation

!> Derivatives, smoothing and the smooth-R-factor Y function of intensity
!! curves.
module derivatives
  use kinds, only: sp
  implicit none

  private

  public :: der, smooth, ysmooth

contains

  !> First derivative y1 of y (tabulated on ne points with step h), after
  !! Zanazzi and Jona, by ne:
  !! ne >= 23: 7-point central formula, 4-point one-sided at 3 end points;
  !! 3 <= ne <= 22: central differences, 3-point one-sided at end points;
  !! ne = 2: two-point difference; ne = 1: derivative set to zero.
  subroutine der(y, ne, y1, h)
    integer, intent(in) :: ne
    real(sp), intent(in) :: y(ne), h
    real(sp), intent(out) :: y1(ne)

    integer :: j, jf, m

    if (ne >= 23) then
      jf = ne - 3
      do j = 4, jf
        y1(j) = (y(j + 3) - 9._sp*y(j + 2) + 45._sp*y(j + 1) - 45._sp*y(j - 1) &
                 + 9._sp*y(j - 2) - y(j - 3)) / (60._sp*h)
      end do
      do j = 1, 3
        y1(j) = (2._sp*y(j + 3) - 9._sp*y(j + 2) + 18._sp*y(j + 1) - 11._sp*y(j)) &
                / (6._sp*h)
        m = ne - j + 1
        y1(m) = (11._sp*y(m) - 18._sp*y(m - 1) + 9._sp*y(m - 2) - 2._sp*y(m - 3)) &
                / (6._sp*h)
      end do
    else if (ne >= 3) then
      jf = ne - 1
      do j = 2, jf
        y1(j) = (y(j + 1) - y(j - 1)) / (2._sp*h)
      end do
      y1(1) = (-3._sp*y(1) + 4._sp*y(2) - y(3)) / (2._sp*h)
      y1(jf + 1) = (3._sp*y(jf + 1) - 4._sp*y(jf) + y(jf - 1)) / (2._sp*h)
    else if (ne == 1) then
      y1(1) = 0._sp
    else
      y1(1) = (y(2) - y(1)) / h
      y1(2) = y1(1)
    end if
  end subroutine der

  !> Smooths the curve (ee, ae) of n points ismoth times by weighted
  !! three-point averaging; on equal intervals a simplified formula is used.
  !! Curves of fewer than 3 points are left unchanged.
  subroutine smooth(ismoth, ae, ee, n)
    integer, intent(in) :: ismoth, n
    real(sp), intent(inout) :: ae(n)
    real(sp), intent(in) :: ee(n)

    real(sp) :: am, af, e21, e32
    integer :: i, ie

    do i = 1, ismoth
      if (n >= 3) then
        am = ae(1)
        do ie = 2, n - 1
          af = ae(ie)
          e21 = ee(ie) - ee(ie - 1)
          e32 = ee(ie + 1) - ee(ie)
          if (abs(e21 - e32) < 0.0001_sp) then
            ae(ie) = 0.25_sp * (2._sp*af + am + ae(ie + 1))
          else
            ae(ie) = 0.5_sp * (af + (e32*am + e21*ae(ie + 1))/(e21 + e32))
          end if
          am = af
        end do
      end if
    end do
  end subroutine smooth

  !> Y function of the smooth R factor R_S (Imre et al., J. Phys.: Condens.
  !! Matter 38, 105001) for beams 1..nb of the curves a(ib, 1:ne(ib)) with
  !! first derivatives ap (per eV). The output y is 2*vi times Y_S of the
  !! paper; the constant factor cancels in R_S (same vi for experiment and
  !! theory). The second derivative is taken from three points of a (grid
  !! step eincr); at the ends from the nearest three, and zero for curves of
  !! fewer than 3 points. Entries beyond ne(ib) are not set.
  subroutine ysmooth(a, ap, nbd, ngp, nb, ne, y, vi, eincr)
    integer, intent(in) :: nbd, ngp, nb
    real(sp), intent(in) :: a(nbd, ngp), ap(nbd, ngp)
    integer, intent(in) :: ne(nbd)
    real(sp), intent(inout) :: y(nbd, ngp)
    real(sp), intent(in) :: vi, eincr

    real(sp), parameter :: alpha = 4.0_sp, beta = 0.15_sp
    real(sp) :: af, deriv1, deriv2, ydensq, y1, y2
    integer :: ib, ie, n

    do ib = 1, nb
      n = ne(ib)
      if (n /= 0) then
        do ie = 1, n
          af = a(ib, ie)
          ! avoids division by 0 in I'/sqrt(I*I+...)
          if (af <= 1.0e-10_sp) then
            af = 1.0e-10_sp
          end if
          deriv1 = ap(ib, ie)*vi
          if (n < 3) then
            deriv2 = 0.0_sp
          else if (ie == 1) then
            deriv2 = a(ib, ie + 2) + a(ib, ie) - 2*a(ib, ie + 1)
          else if (ie < n) then
            deriv2 = a(ib, ie + 1) + a(ib, ie - 1) - 2*af
          else
            deriv2 = a(ib, ie) + a(ib, ie - 2) - 2*a(ib, ie - 1)
          end if
          ydensq = af*af + 4*deriv1*deriv1
          ! threshold 1e-10, not 0, avoids overflow in the division by I''
          if (deriv2 > 1.0e-10_sp) then
            deriv2 = deriv2*(vi/eincr)**2
            y1 = alpha*af/deriv2 - (0.5_sp*alpha)*deriv1*deriv1/(deriv2*deriv2) + beta
            if (y1 > 0) then
              y2 = y1/sqrt(1.0_sp + y1*y1)
              ydensq = ydensq + y2*y2*deriv2*deriv2
            end if
          end if
          y(ib, ie) = 2*deriv1/sqrt(ydensq)
        end do
      end if
    end do
  end subroutine ysmooth
end module derivatives

!> Pendry Y function, integrals over the energy range, and the Zanazzi-Jona
!! normalisation, as they enter the R factors.
module rfactor_terms
  use kinds, only: sp
  use interpolation, only: interpolation_cache, yval
  implicit none

  private

  public :: epszj, varsum, ypend

contains

  !> Largest |aep(ib, ie)| for ie = ne1..ne2 (eps of the Zanazzi-Jona
  !! R factor, aep being the derivative of the experimental curve).
  subroutine epszj(aep, nbd, ngp, ib, ne1, ne2, eps)
    integer, intent(in) :: nbd, ngp, ib, ne1, ne2
    real(sp), intent(in) :: aep(nbd, ngp)
    real(sp), intent(out) :: eps

    real(sp) :: a
    integer :: ie

    eps = 0._sp
    do ie = ne1, ne2
      a = abs(aep(ib, ie))
      if (a > eps) eps = a
    end do
  end subroutine epszj

  !> Integral s of y(1:n) (step h) by the trapezoid rule, after Zanazzi and
  !! Jona.
  subroutine trapezoid(y, n, h, s)
    integer, intent(in) :: n
    real(sp), intent(in) :: y(n), h
    real(sp), intent(out) :: s

    real(sp) :: a
    integer :: j

    a = 0._sp
    do j = 1, n
      a = a + y(j)
    end do
    s = a - 0.5_sp*(y(1) + y(n))
    s = s*h
  end subroutine trapezoid

  !> Integral s over energy points ie = ie1..ie2 of the curves a1, b1 of
  !! beam ib1 and, shifted by nv points, a2, b2 of beam ib2 (step eincr).
  !! The integrand depends on nf:
  !!   nf = 1 (and any nf other than 2, 3, 4): a1
  !!   nf = 2: a1**2
  !!   nf = 3: (a1 - c*a2)**2
  !!   nf = 4: |b1 - c*b2| * |a1 - c*a2| / (|a1| + eps), the Zanazzi-Jona
  !!           integrand, interpolated onto a 10 times denser grid first.
  !! For nf = 1..3 the integrand is left in y(1:ie2-ie1+1); for nf = 4,
  !! y(1:ie2-ie1+1) holds the energies relative to the first point.
  subroutine varsum(a1, a2, b1, b2, nbd1, nbd2, ngp, ib1, ib2, ie1, ie2, nv, &
                    eincr, eps, c, nf, s, y)
    integer, intent(in) :: nbd1, nbd2, ngp, ib1, ib2, ie1, ie2, nv, nf
    real(sp), intent(in) :: a1(nbd1, ngp), b1(nbd1, ngp)
    real(sp), intent(in) :: a2(nbd2, ngp), b2(nbd2, ngp)
    real(sp), intent(in) :: eincr, eps, c
    real(sp), intent(out) :: s
    real(sp), intent(inout) :: y(ngp)

    real(sp) :: y1(ngp), y2(ngp), y3(ngp), y4(ngp), yy(10*ngp)
    type(interpolation_cache) :: cache
    real(sp) :: aa1, aa2, ab1, ab2, de, x
    integer :: ie, ies, n, nn

    n = 0
    if (nf /= 4) then
      do ie = ie1, ie2
        n = n + 1
        ies = ie + nv
        select case (nf)
        case (2)
          y(n) = a1(ib1, ie)**2
        case (3)
          y(n) = (a1(ib1, ie) - c*a2(ib2, ies))**2
        case default
          y(n) = a1(ib1, ie)
        end select
      end do
      call trapezoid(y(1:n), n, eincr, s)
      return
    end if

    do ie = ie1, ie2
      n = n + 1
      ies = ie + nv
      y(n) = real(n - 1, sp)*eincr
      y1(n) = a1(ib1, ie)
      y2(n) = a2(ib2, ies)
      y3(n) = b1(ib1, ie)
      y4(n) = b2(ib2, ies)
    end do
    de = 0.1_sp*eincr
    nn = 10*(n - 1) + 1
    do ie = 1, nn
      x = real(ie - 1, sp)*de
      cache%il = 0
      cache%ih = 0
      aa1 = yval(x, y1(1:n), y(1:n), n, cache)
      cache%il = 0
      cache%ih = 0
      aa2 = yval(x, y2(1:n), y(1:n), n, cache)
      cache%il = 0
      cache%ih = 0
      ab1 = yval(x, y3(1:n), y(1:n), n, cache)
      cache%il = 0
      cache%ih = 0
      ab2 = yval(x, y4(1:n), y(1:n), n, cache)
      yy(ie) = abs(ab1 - c*ab2)*abs(aa1 - c*aa2)/(abs(aa1) + eps)
    end do
    call trapezoid(yy(1:nn), nn, de, s)
  end subroutine varsum

  !> Pendry Y function y = L/(1 + vi**2 L**2) of the curves a(ib, 1:ne(ib))
  !! for beams 1..nb, with L = ap/a the logarithmic derivative (ap: first
  !! derivative per eV, vi: imaginary part of the inner potential). Where
  !! |a| < 1e-7 the equivalent form (a/ap)/((a/ap)**2 + vi**2) is used, and
  !! y = 0 where also |ap| < 1e-7. Entries beyond ne(ib) are not set.
  subroutine ypend(a, ap, nbd, ngp, nb, ne, y, vi)
    integer, intent(in) :: nbd, ngp, nb
    real(sp), intent(in) :: a(nbd, ngp), ap(nbd, ngp)
    integer, intent(in) :: ne(nbd)
    real(sp), intent(inout) :: y(nbd, ngp)
    real(sp), intent(in) :: vi

    real(sp) :: af, apf
    integer :: ib, ie, n

    do ib = 1, nb
      n = ne(ib)
      if (n /= 0) then
        do ie = 1, n
          af = a(ib, ie)
          if (abs(af) >= 1.e-7_sp) then
            af = ap(ib, ie)/af
            y(ib, ie) = af/(1._sp + vi*vi*af*af)
          else
            apf = ap(ib, ie)
            if (abs(apf) >= 1.e-7_sp) then
              af = af/apf
              y(ib, ie) = af/(af*af + vi*vi)
            else
              y(ib, ie) = 0._sp
            end if
          end if
        end do
      end if
    end do
  end subroutine ypend
end module rfactor_terms

!> Averaging of theoretical beams over domains and beam groups, and the
!! check of averaged intensities.
module domain_averaging
  use kinds, only: sp
  implicit none

  private

  public :: maxint, rintav

contains

  !> Largest a(is, ib, ie) over all is = 1..ns, ie = 1..ne, and 0 if none
  !! is positive.
  subroutine maxint(a, ns, nbd, ib, ne, am)
    integer, intent(in) :: ns, nbd, ib, ne
    real(sp), intent(in) :: a(ns, nbd, ne)
    real(sp), intent(out) :: am

    integer :: ie, is

    am = 0._sp
    do is = 1, ns
      do ie = 1, ne
        if (a(is, ib, ie) > am) am = a(is, ib, ie)
      end do
    end do
  end subroutine maxint

  !> Averages the intensities at(:, jb, :) of all beams jb with the same
  !! group label kav(jb) into atav(:, lav, :), lav = 1, 2, ... being the
  !! group labels in increasing order (lavm groups found). Each beam is
  !! weighted by the number of symmetry-equivalent beams its symmetry code
  !! sym(jb) stands for in the LEED program (all 1 without symmetry). The
  !! group takes the name pq1(:, lav) of its first beam in pq. Entries of
  !! atav and pq1 for labels without beams are not set.
  subroutine rintav(at, ns, nb, ne, pq, pq1, kav, sym, lavm, atav)
    integer, intent(in) :: ns, nb, ne
    real(sp), intent(in) :: at(ns, nb, ne), pq(2, nb)
    real(sp), intent(inout) :: pq1(2, nb), atav(ns, nb, ne)
    integer, intent(in) :: kav(nb), sym(nb)
    integer, intent(out) :: lavm

    real(sp) :: wght, wt
    integer :: ib, ie, is, jb, lab, lav

    lavm = 0
    do lav = 1, nb
      do ib = 1, nb
        ! the first beam of group lav
        if (lav /= kav(ib)) cycle
        lavm = lavm + 1
        pq1(:, lav) = pq(:, ib)
        wt = 0.0_sp
        do jb = ib, nb
          if (lav /= kav(jb)) cycle
          lab = sym(jb)
          if (lab < 2) then
            wght = 1.0_sp
          else if (lab < 7) then
            wght = 2.0_sp
          else if (lab < 10) then
            wght = 4.0_sp
          else if (lab == 10) then
            wght = 8.0_sp
          else if (lab < 12) then
            wght = 3.0_sp
          else if (lab < 15) then
            wght = 6.0_sp
          else
            wght = 12.0_sp
          end if
          if (wt < 0.001_sp) then
            do is = 1, ns
              do ie = 1, ne
                atav(is, lav, ie) = wght*at(is, jb, ie)
              end do
            end do
          else
            do is = 1, ns
              do ie = 1, ne
                atav(is, lav, ie) = atav(is, lav, ie) + wght*at(is, jb, ie)
              end do
            end do
          end if
          wt = wt + wght
        end do
        do is = 1, ns
          do ie = 1, ne
            atav(is, lav, ie) = atav(is, lav, ie)/wt
          end do
        end do
        exit
      end do
    end do
  end subroutine rintav
end module domain_averaging

!> Control parameters of an R-factor run, read from the head of WEXPEL
!! (unit 8): namelists NL1 to NL3, the beam averaging, relation and
!! labelling lists, and the title of the experiment.
module control_input
  use kinds, only: sp
  use array_growth, only: grow
  use rfactor_output, only: write_beam_counts, write_beam_list, write_nl1, &
                            write_title
  implicit none

  private

  public :: check_beam_relation, control_parameters, nrfact, read_control

  !> Number of R factors: R2, Zanazzi-Jona, Pendry and Smooth.
  integer, parameter :: nrfact = 4

  !> Input parameters of an R-factor run.
  type :: control_parameters
    ! NL1: theory energy range (emin, emax) and grid step (eincr), all in
    ! eV; number of theory files (limfil); print level (ipr); imaginary
    ! (vi) and real (v0rr) inner potential; range and step of the inner
    ! potential shift (v01, v02, vincr); number of 3-point smoothings of
    ! the experiment (ismoth); experiment (eot = 0) or pseudo-experiment
    ! (eot = 1); plot output if plot = 1; gaps in the experiment if gap = 1
    real(sp) :: emin, emax, eincr, vi, v0rr, v01, v02, vincr, plot
    integer :: limfil, ipr, ismoth, eot, gap
    ! NL2: weights of the experimental beams (wb) and of the R factors (wr;
    ! wrsum is their sum)
    real(sp), allocatable :: wb(:)
    real(sp) :: wr(nrfact), wrsum
    ! NL3: plot normalisation (norm), maximum intensity (intmax), aspect
    ! ratio (plsize) and energy tic interval (xtics)
    integer :: norm, xtics
    real(sp) :: intmax, plsize(2)
    ! theory beam groups to average (kav), experimental beams to average
    ! (nbea), theory beam group of each experimental beam (ibp), label of
    ! each theory beam group (mittel: 1 integer, 2 fractional order) and
    ! the number of integer (nint) and fractional (nfrac) groups
    integer, allocatable :: kav(:), nbea(:), ibp(:), mittel(:)
    integer :: nint, nfrac
  end type control_parameters

contains

  !> Reads the control parameters from unit 8 for nbtd theory beams and
  !! ns geometries, echoes them to stdout together with the title of the
  !! experiment, and returns the number of geometries to use (nss). The
  !! number of experimental beams is the length of their averaging scheme
  !! nbea, which has a nonzero entry for every beam. Stops on unsupported
  !! input.
  subroutine read_control(ctl, nbtd, ns, nss)
    type(control_parameters), intent(out) :: ctl
    integer, intent(in) :: nbtd, ns
    integer, intent(out) :: nss

    real(sp) :: emin, emax, eincr, vi, v0rr, v01, v02, vincr, plot
    integer :: limfil, ipr, ismoth, eot, gap
    integer :: nssk(ns)
    real(sp), allocatable :: wb(:)
    real(sp) :: wr(nrfact), wrsum
    integer :: norm, xtics
    real(sp) :: intmax, plsize(2)
    character(len=80) :: title
    integer :: nbed, nbea(25), i, j, is

    namelist /nl1/ emin, emax, eincr, limfil, ipr, vi, v0rr, v01, v02, &
                   vincr, ismoth, eot, plot, gap
    namelist /nl2/ nssk, wb, wr
    namelist /nl3/ norm, intmax, plsize, xtics

4   format(400(25I4, /))
9   format('MITTEL(,', I4, ') MUST BE EITHER 0, 1 OR 2', /, &
           '   *****   PLEASE CORRECT   *****   ')
22  format(A80)

    emin = 20.0_sp
    emax = 220.0_sp
    eincr = 1.0_sp
    limfil = 1
    ipr = 1
    vi = 4.0_sp
    v01 = -5.0_sp
    v02 = 5.0_sp
    vincr = 2.0_sp
    v0rr = 0.0_sp
    ismoth = 2
    plot = 0.0_sp
    gap = 0

    nssk = 0
    ! default: Pendry R factor only
    wr(1) = 0.0_sp
    wr(2) = 0.0_sp
    wr(3) = 1.0_sp
    wr(4) = 0.0_sp

    norm = 1
    intmax = 999.99_sp
    plsize(1) = 1.0_sp
    plsize(2) = 1.0_sp
    xtics = 50

    ctl%nint = 0
    ctl%nfrac = 0

    read(8, nl1)
    call write_nl1(emin, emax, eincr, limfil, ipr, vi, v0rr, v01, v02, vincr, &
                   ismoth, eot, plot, gap)

    if (limfil /= 1) then
      write(0, *) 'rfactor: LIMFIL must be 1; reading several ', &
                  'theory files is not supported'
      stop 1
    end if
    if (v0rr /= 0.0_sp) then
      write(0, *) 'rfactor: V0RR must be 0; a real inner-potential ', &
                  'offset is not supported'
      stop 1
    end if

    allocate(ctl%kav(nbtd), ctl%mittel(nbtd))

    read(8, 4) (ctl%kav(j), j = 1, nbtd)

    ! nbea: 25 entries per record; it ends with the first record whose
    ! last entry is zero (an empty record after a multiple of 25 beams)
    allocate(ctl%nbea(0))
    do
      read(8, '(25I4)') nbea
      nbed = size(ctl%nbea)
      j = findloc(nbea /= 0, .true., dim=1, back=.true.)
      call grow(ctl%nbea, nbed + j)
      ctl%nbea(nbed + 1:) = nbea(1:j)
      if (j < 25) exit
    end do
    nbed = size(ctl%nbea)
    allocate(ctl%ibp(nbed), wb(nbed))
    wb = 1.0_sp

    call write_beam_counts(nbtd, nbed)
    call write_beam_list('*AVERAGING SCHEME OF THEOR BEAMS', ctl%kav(1:nbtd))
    call write_beam_list('*AVERAGING SCHEME OF EXP. BEAMS', ctl%nbea(1:nbed))
    if (nbed == 0 .or. any(ctl%nbea == 0)) then
      write(0, *) 'rfactor: NBEA = 0 (skipping an experimental ', &
                  'beam) is not supported'
      stop 1
    end if

    read(8, 4) (ctl%ibp(j), j = 1, nbed)
    call write_beam_list('*RELATIONSHIP OF THEOR AND EXP BEAMS', &
                         ctl%ibp(1:nbed))

    read(8, 4) (ctl%mittel(j), j = 1, nbtd)
    call write_beam_list('*LABELLING OF THEOR BEAM GROUPINGS', &
                         ctl%mittel(1:nbtd))

    do i = 1, nbtd
      if (ctl%mittel(i) == 1) then
        ctl%nint = ctl%nint + 1
      else if (ctl%mittel(i) == 2) then
        ctl%nfrac = ctl%nfrac + 1
      else if (ctl%mittel(i) /= 0) then
        write(6, 9) i
        stop
      end if
    end do

    read(8, nl2)
    do is = 1, ns
      if (nssk(is) /= 0) then
        write(0, *) 'rfactor: NSSK must be 0; skipping geometries ', &
                    'is not supported'
        stop 1
      end if
    end do

    ! viperleed writes only the first three weights; if they are all zero,
    ! the smooth R factor R_S is used
    wrsum = 0.0_sp
    do i = 1, nrfact
      wrsum = wrsum + wr(i)
    end do
    if (wrsum == 0.0_sp) then
      wr(4) = 1.0_sp
      wrsum = 1.0_sp
    end if

    nss = 0
    do is = 1, ns
      if (nssk(is) == 0) nss = nss + 1
    end do

    read(8, nl3)

    read(8, 22) title
    call write_title(title)

    ctl%emin = emin
    ctl%emax = emax
    ctl%eincr = eincr
    ctl%vi = vi
    ctl%v0rr = v0rr
    ctl%v01 = v01
    ctl%v02 = v02
    ctl%vincr = vincr
    ctl%plot = plot
    ctl%limfil = limfil
    ctl%ipr = ipr
    ctl%ismoth = ismoth
    ctl%eot = eot
    ctl%gap = gap
    ctl%wb = wb
    ctl%wr = wr
    ctl%wrsum = wrsum
    ctl%norm = norm
    ctl%xtics = xtics
    ctl%intmax = intmax
    ctl%plsize = plsize
  end subroutine read_control

  !> Stops if one of the nbe averaged experimental beams has no theory beam
  !! (ibp = 0).
  subroutine check_beam_relation(ibp, nbe)
    integer, intent(in) :: ibp(:), nbe

    integer :: ib

    do ib = 1, nbe
      if (ibp(ib) == 0) then
        write(0, *) 'rfactor: IBP = 0 (skipping a theoretical beam) ', &
                    'is not supported'
        stop 1
      end if
    end do
  end subroutine check_beam_relation
end module control_input

!> Reading of the experimental (or pseudo-experimental) beams and their
!! gaps.
module experiment_input
  use kinds, only: sp
  use array_growth, only: grow
  use control_input, only: control_parameters
  use rfactor_output, only: write_exp_beam_name, write_exp_curve, write_gap, &
                            write_gap_count, write_pseudo_beam_count, &
                            write_pseudo_beam_name
  implicit none

  private

  public :: experiment_beams, gapin, read_experiment, reade, readt

  !> Experimental beams; after expav the first nbe are the averaged
  !! beams.
  type :: experiment_beams
    integer :: nbe
    ! number of points (nee), names (bename), intensities (ae), energies
    ! (ee), first and second derivatives (aep, aepp), Pendry and smooth
    ! Y functions (ype, yse) of each beam
    integer, allocatable :: nee(:)
    character(len=4), allocatable :: bename(:, :)
    real(sp), allocatable :: ae(:, :), ee(:, :), aep(:, :), aepp(:, :)
    real(sp), allocatable :: ype(:, :), yse(:, :)
    ! integrals over the intensity, the squared intensity and the squared
    ! Pendry and smooth Y functions of each beam
    real(sp), allocatable :: tse(:), tse2(:), tseyp2(:), tseys2(:)
    ! gaps: total number (ngap), number per beam (gcount), beam of each
    ! gap (gbm), its energy range (eg1, eg2) and grid indices (ng1, ng2)
    integer :: ngap
    integer, allocatable :: gcount(:), gbm(:), ng1(:), ng2(:)
    real(sp), allocatable :: eg1(:), eg2(:)
  end type experiment_beams

contains

  !> Reads the experiment (eot = 0) or pseudo-experiment from unit 8 into
  !! ex, one beam per entry of ctl%nbea, and the gaps if ctl%gap = 1.
  subroutine read_experiment(ex, ctl)
    type(experiment_beams), intent(out) :: ex
    type(control_parameters), intent(in) :: ctl

    integer :: nbed

    nbed = size(ctl%nbea)
    allocate(ex%nee(nbed), ex%bename(5, nbed), ex%gcount(nbed))
    allocate(ex%ae(nbed, 0), ex%ee(nbed, 0))
    allocate(ex%tse(nbed), ex%tse2(nbed), ex%tseyp2(nbed), ex%tseys2(nbed))
    allocate(ex%gbm(0), ex%ng1(0), ex%ng2(0), ex%eg1(0), ex%eg2(0))
    ex%bename = repeat(char(0), 4)

    if (ctl%eot <= 0) then
      call reade(ex%ae, ex%ee, nbed, ex%nee, ex%bename, ctl%ipr)
    else
      call readt(ex%ae, ex%ee, nbed, ex%nee, ex%bename)
    end if

    ex%gcount = 0
    if (ctl%gap == 1) then
      call gapin(ex%ngap, ex%gbm, ex%eg1, ex%eg2, ex%gcount, nbed)
    end if
  end subroutine read_experiment

  !> Reads nbed experimental beams from unit 8: a format line, then per beam
  !! its name, the number of points nee(ib) with a factor fac, and the pairs
  !! (ee(ib, ie), ae(ib, ie)). The intensities are divided by fac, meant to
  !! bring them to the order of 1 and to match curves averaged together.
  !! ae and ee grow to the longest beam.
  subroutine reade(ae, ee, nbed, nee, bename, ipr)
    integer, intent(in) :: nbed, ipr
    real(sp), allocatable, intent(inout) :: ae(:, :), ee(:, :)
    integer, intent(out) :: nee(nbed)
    character(len=4), intent(out) :: bename(5, nbed)

    character(len=4) :: fmt(19)
    real(sp) :: fac
    integer :: i, ib, ie, n

10  format(19A4)
35  format(1I4, 1E14.4)

    read(8, 10) (fmt(i), i = 1, 19)
    do ib = 1, nbed
      read(8, 10) (bename(i, ib), i = 1, 5)
      call write_exp_beam_name(ib, bename(:, ib))
      read(8, 35) nee(ib), fac
      n = nee(ib)
      call grow(ae, n)
      call grow(ee, n)
      read(8, fmt) (ee(ib, ie), ae(ib, ie), ie = 1, n)
      if (ipr >= 1) call write_exp_curve(ee(ib, 1:n), ae(ib, 1:n))
      do ie = 1, n
        ae(ib, ie) = ae(ib, ie)/fac
      end do
    end do
  end subroutine reade

  !> Reads theoretical curves from unit 8 as pseudo-experiment: the number
  !! of beams (must equal nbed), the beam names, then per energy the
  !! intensities of all beams, up to the end of the file. Leading points
  !! with intensity <= 1e-20 are dropped per beam; nee(ib) counts the points
  !! above 1e-20. ae and ee grow to the number of energies.
  subroutine readt(ae, ee, nbed, nee, bename)
    integer, intent(in) :: nbed
    real(sp), allocatable, intent(inout) :: ae(:, :), ee(:, :)
    integer, intent(out) :: nee(nbed)
    character(len=4), intent(out) :: bename(5, nbed)

    real(sp), allocatable :: eeh(:)
    real(sp) :: dummy
    integer :: ioff(nbed)
    integer :: i, ib, inn, j, k, nbeams

101 format(I5, 22A4)
103 format(I5, E13.4)
105 format(F7.2, F7.4, 4E14.5, /, 2000(5E14.5, /))

    read(8, 103) nbeams
    call write_pseudo_beam_count(nbeams)
    if (nbed /= nbeams) then
      write(7, '(A33)') '*********************************'
      write(7, '(A33)') 'SOMETHING IS WRONG WITH NUMBER OF'
      write(7, '(A33)') 'PSEUDOEXPERIMENTAL BEAMS!        '
      write(7, '(A33)') 'NO AVERAGING ALLOWED!            '
      write(7, '(A33)') '*********************************'
      stop 'ERROR 1 IN PSEUDOEXPERIMENT'
    end if
    do ib = 1, nbed
      read(8, 101) inn, (bename(i, ib), i = 1, 5)
      call write_pseudo_beam_name(inn, bename(:, ib))
    end do

    do k = 1, nbed
      nee(k) = 0
      ioff(k) = 0
    end do

    allocate(eeh(0))
    j = 0
    do
      j = j + 1
      call grow(eeh, j)
      call grow(ae, j)
      read(8, 105, end=250) eeh(j), dummy, (ae(k, j), k = 1, nbed)
      do k = 1, nbed
        if (ae(k, j) > 1.e-20_sp) then
          nee(k) = nee(k) + 1
        else if (nee(k) == 0) then
          ioff(k) = ioff(k) + 1
        end if
      end do
    end do

250 continue
    call grow(ee, size(ae, 2))
    do k = 1, nbed
      do j = 1, nee(k)
        ee(k, j) = eeh(j + ioff(k))
        ae(k, j) = ae(k, j + ioff(k))
      end do
    end do
  end subroutine readt

  !> Reads the gaps of the experimental spectra from unit 8: their number
  !! ngap, then per gap the beam bm(ig) and the energy range
  !! [eg1(ig), eg2(ig)], which is left out of smoothing, interpolation and
  !! the R factor. Gaps must be ordered by beam, and within a beam by
  !! energy; otherwise the program stops. gcount(ib) is incremented for
  !! every gap of beam ib. bm, eg1 and eg2 are sized for ngap gaps.
  subroutine gapin(ngap, bm, eg1, eg2, gcount, nbed)
    integer, intent(out) :: ngap
    integer, allocatable, intent(inout) :: bm(:)
    real(sp), allocatable, intent(inout) :: eg1(:), eg2(:)
    integer, intent(in) :: nbed
    integer, intent(inout) :: gcount(nbed)

    integer :: ig

10  format(I3, 2F7.2)
30  format(' *** GAPS MUST BE IN ORDER OF EXPERIMENTAL BEAMS ***', /, &
           ' ***               PROGRAM STOPPED               ***')
40  format(' *** GAPS IN BEAM', I3, ' MUST BE ORDERED BY ENERGY ***', /, &
           ' ***              PROGRAM STOPPED              ***')

    read(8, 10) ngap
    call write_gap_count(ngap)
    call grow(bm, ngap)
    call grow(eg1, ngap)
    call grow(eg2, ngap)

    do ig = 1, ngap
      read(8, 10) bm(ig), eg1(ig), eg2(ig)
      if (ig > 1) then
        if (bm(ig) < bm(ig - 1)) then
          write(6, 30)
          stop
        else if (bm(ig) == bm(ig - 1) .and. eg1(ig) < eg1(ig - 1)) then
          write(6, 40) bm(ig)
          stop
        end if
      end if
      call write_gap(ig, bm(ig), eg1(ig), eg2(ig))
      gcount(bm(ig)) = gcount(bm(ig)) + 1
    end do
  end subroutine gapin
end module experiment_input

!> Reading of the theoretical beams written by the LEED program, for all
!! geometries.
module theory_input
  use kinds, only: sp
  use array_growth, only: grow
  implicit none

  private

  public :: read_theory_records, theory_records

  !> The theory as read from stdin: title, number of beams nbtd, their
  !! names pq(:, j) and symmetry codes sym(j), and nrec records of one
  !! geometry at one energy each: energy e, geometry code spac and the
  !! intensities y(:, irec) of all beams. ns geometries follow each other
  !! at every energy.
  type :: theory_records
    character(len=80) :: title
    integer :: nbtd, ns, nrec
    integer, allocatable :: sym(:)
    real(sp), allocatable :: pq(:, :), e(:), spac(:), y(:, :)
  end type theory_records

contains

  !> Reads the theory from unit 5 up to the end of the file into tr: the
  !! title, the number of beams, per beam its label (not used), name and
  !! symmetry code, then the records of all geometries and energies. The
  !! number of geometries is the number of leading records at the first
  !! energy. Stops if there is no record.
  subroutine read_theory_records(tr)
    type(theory_records), intent(out) :: tr

    integer, allocatable :: lab(:)
    integer :: i, j, n

20  format(I5, 2F10.5, I3)
22  format(A80)
26  format(I5)
40  format(F7.2, F7.4, 4E14.5, /, 2000(5E14.5, /))

    read(5, 22) tr%title
    read(5, 26) tr%nbtd

    allocate(lab(tr%nbtd), tr%pq(2, tr%nbtd), tr%sym(tr%nbtd))
    do j = 1, tr%nbtd
      read(5, 20) lab(j), (tr%pq(i, j), i = 1, 2), tr%sym(j)
    end do

    allocate(tr%e(0), tr%spac(0), tr%y(tr%nbtd, 0))
    n = 0
    do
      if (n == size(tr%e)) then
        call grow(tr%e, 2*n + 1)
        call grow(tr%spac, 2*n + 1)
        call grow(tr%y, 2*n + 1)
      end if
      read(5, 40, end=10) tr%e(n + 1), tr%spac(n + 1), tr%y(:, n + 1)
      n = n + 1
    end do
10  tr%nrec = n

    tr%ns = 0
    do while (tr%ns < n)
      if (tr%e(tr%ns + 1) /= tr%e(1)) exit
      tr%ns = tr%ns + 1
    end do
    if (tr%ns == 0) then
      write(0, *) 'rfactor: no theoretical intensities on stdin'
      stop 1
    end if
  end subroutine read_theory_records
end module theory_input

!> Energy interval common to an experimental and a theoretical beam.
module common_interval
  use kinds, only: sp
  implicit none

  private

  public :: comnei

contains

  !> Finds the energy interval common to experimental beam ibe (grid ee)
  !! and theoretical beam ibt (grid et, shifted by the inner potential v0),
  !! both on steps of eincr: it is ee(ibe, ne1:ne2) and et(ibt, nt1:nt2).
  subroutine comnei(ee, nbed, nee, et, nbtd, ngp, net, ibe, ibt, v0, eincr, &
                    ne1, ne2, nt1, nt2)
    integer, intent(in) :: nbed, nbtd, ngp, ibe, ibt
    real(sp), intent(in) :: ee(nbed, ngp), et(nbtd, ngp), v0, eincr
    integer, intent(in) :: nee(nbed), net(nbtd)
    integer, intent(out) :: ne1, ne2, nt1, nt2

    real(sp) :: de1, de2
    integer :: ne, nt

    ne = nee(ibe)
    nt = net(ibt)
    de1 = et(ibt, 1) + v0 - ee(ibe, 1)
    de2 = et(ibt, nt) + v0 - ee(ibe, ne)
    if (de1 < 0._sp) then
      ne1 = 1
      nt1 = int((-de1/eincr) + 0.0001_sp) + 1
    else
      ne1 = int((de1/eincr) + 0.0001_sp) + 1
      nt1 = 1
    end if
    if (de2 < 0._sp) then
      ne2 = ne - int((-de2/eincr) + 0.0001_sp)
      nt2 = nt
    else
      ne2 = ne
      nt2 = nt - int((de2/eincr) + 0.0001_sp)
    end if
  end subroutine comnei
end module common_interval

!> Preparation of the experimental beams for the R factors: averaging,
!! smoothing, interpolation onto the working grid (steps of eincr),
!! derivatives, Y functions and integrals over the whole curves.
module experiment_preparation
  use kinds, only: sp
  use control_input, only: control_parameters
  use array_growth, only: grow
  use interpolation, only: grid, grid_points, interpolation_cache, intpol, &
                           yval
  use derivatives, only: der, smooth, ysmooth
  use rfactor_terms, only: varsum, ypend
  use experiment_input, only: experiment_beams
  use rfactor_output, only: write_interpolated, write_pendry_y, &
                            write_smooth_y, write_smoothed, &
                            write_too_few_points
  implicit none

  private

  public :: expav, experiment_grid_size, prepare_experiment

contains

  !> Averages the experimental beams with the same group label nbea(ib) onto
  !! a common energy grid in increasing order, leaving out points inside
  !! the gaps (eg1, eg2) of each beam. The nbe averaged beams replace the
  !! first rows of ae, ee, nee and bename (each named after its first beam);
  !! their gaps are returned in eg1, eg2 (energies), ng1, ng2 (indices),
  !! ngap in total and gcount per beam. nbea is set to 0. ae and ee grow
  !! to the longest averaged beam, the gap arrays to the number of gaps.
  subroutine expav(ae, ee, nbed, nee, bename, nbea, nbe, gcount, ngap, eg1, &
                   eg2, ng1, ng2)
    integer, intent(in) :: nbed
    real(sp), allocatable, intent(inout) :: ae(:, :), ee(:, :)
    integer, intent(inout) :: nee(nbed), nbea(nbed), gcount(nbed)
    character(len=4), intent(inout) :: bename(5, nbed)
    integer, intent(out) :: nbe, ngap
    real(sp), allocatable, intent(inout) :: eg1(:), eg2(:)
    integer, allocatable, intent(inout) :: ng1(:), ng2(:)

    ! sum of the intensities (a) and number of values (nv) per energy
    real(sp), allocatable :: a(:)
    integer, allocatable :: nv(:)
    integer :: gc(nbed)
    real(sp) :: de, e, ediff, emax, emin
    integer :: gin, i, ib, ib1, ibmin, ie, iemax, ien, ig, igap, igc, &
               nbeat, nemax

    nbe = 1
    igap = 0
    ngap = 0
    allocate(a(0), nv(0))

    do ib = 1, nbed
      if (nbea(ib) /= 0) then
        emin = 1.e6_sp
        emax = 0._sp
        de = 1.e6_sp

        ! energy range of all beams averaged with ib (including ib)
        do ib1 = ib, nbed
          if (nbea(ib1) == nbea(ib)) then
            do ie = 1, nee(ib1)
              if (emin > ee(ib1, ie)) then
                emin = ee(ib1, ie)
                ibmin = ib1
              end if
              if (emax < ee(ib1, ie)) then
                emax = ee(ib1, ie)
              end if
            end do
          end if
        end do

        ! energy step: smallest step above 0.01 from emin in beam ibmin
        do ie = 1, nee(ibmin)
          ediff = ee(ibmin, ie) - emin
          if (ediff < de .and. ediff > 0.01_sp) then
            de = ediff
          end if
        end do
        nemax = int((emax - emin)/de + 0.0001_sp) + 1
        call grow(a, nemax)
        call grow(nv, nemax)
        call grow(ae, nemax)
        call grow(ee, nemax)
        do ie = 1, nemax
          nv(ie) = 0
          a(ie) = 0._sp
        end do

        ! sum over the same beams again, leaving out points inside gaps
        iemax = 0
        nbeat = nbea(ib)
        do ib1 = ib, nbed
          if (nbea(ib1) == nbeat) then
            do ie = 1, nee(ib1)
              e = ee(ib1, ie)
              ien = int((e - emin)/de + 0.0001_sp) + 1
              if (ien > iemax) then
                iemax = ien
              end if
              if (gcount(ib1) == 0) then
                nv(ien) = nv(ien) + 1
                a(ien) = a(ien) + ae(ib1, ie)
              else
                gin = 0
                igc = igap
                do i = 1, gcount(ib1)
                  igc = igc + 1
                  if (e >= eg1(igc) .and. e <= eg2(igc)) then
                    gin = 1
                  end if
                end do
                if (gin /= 1) then
                  nv(ien) = nv(ien) + 1
                  a(ien) = a(ien) + ae(ib1, ie)
                end if
              end if
            end do
            nbea(ib1) = 0
          end if
        end do

        ! averages; grid points without data form the gaps of beam nbe
        gin = 0
        ig = 0
        do ie = 1, iemax
          ee(nbe, ie) = emin + real(ie - 1, sp)*de
          if (nv(ie) /= 0) then
            ae(nbe, ie) = a(ie)/real(nv(ie), sp)
            gin = 0
          else if (gin == 0) then
            ig = ig + 1
            if (ngap + ig > min(size(eg1), size(eg2), size(ng1), size(ng2))) then
              call grow(eg1, 2*(ngap + ig))
              call grow(eg2, 2*(ngap + ig))
              call grow(ng1, 2*(ngap + ig))
              call grow(ng2, 2*(ngap + ig))
            end if
            eg1(ngap + ig) = ee(nbe, ie)
            eg2(ngap + ig) = ee(nbe, ie)
            ng1(ngap + ig) = ie
            ng2(ngap + ig) = ie
            gin = 1
          else
            eg2(ngap + ig) = ee(nbe, ie)
            ng2(ngap + ig) = ie
          end if
        end do

        nee(nbe) = iemax
        gc(nbe) = ig
        do i = 1, 5
          bename(i, nbe) = bename(i, ib)
        end do
        ngap = ngap + ig
        nbe = nbe + 1
      end if
      igap = igap + gcount(ib)
    end do

    nbe = nbe - 1
    do ib = 1, nbe
      gcount(ib) = gc(ib)
    end do
  end subroutine expav

  !> Number of points the arrays of the averaged experimental beams need
  !! for the beams and their interpolation onto the working grid (steps of
  !! eincr).
  integer function experiment_grid_size(ex, eincr) result(ngp)
    type(experiment_beams), intent(in) :: ex
    real(sp), intent(in) :: eincr

    integer :: ib, n

    ngp = size(ex%ae, 2)
    do ib = 1, ex%nbe
      n = ex%nee(ib)
      ngp = max(ngp, grid_points(ex%ee(ib, 1:n), eincr))
    end do
  end function experiment_grid_size

  !> Smooths the averaged experimental beams, interpolates them onto the
  !! working grid, and computes their derivatives, Y functions and
  !! integrals, in arrays of ngp points per beam (at least
  !! experiment_grid_size). Beams with gaps are smoothed and interpolated
  !! section by section; intensity and derivatives are zero inside a gap.
  !! Logs the curves after each step as ctl%ipr asks.
  subroutine prepare_experiment(ex, ctl, ngp)
    type(experiment_beams), intent(inout) :: ex
    type(control_parameters), intent(in) :: ctl
    integer, intent(in) :: ngp

    real(sp), dimension(ngp) :: a, e, ap, app, xpl, ypl, worx
    type(interpolation_cache) :: ycache
    real(sp) :: xval
    logical :: too_few
    integer :: nbed, ib, ie, ie2, ne, ne1, igap, igc, nz, newnz, npts, &
               imin, nlen, ilen, ipts

    nbed = size(ex%ae, 1)
    call grow(ex%ae, ngp)
    call grow(ex%ee, ngp)
    allocate(ex%aep(nbed, ngp), ex%aepp(nbed, ngp), ex%ype(nbed, ngp), &
             ex%yse(nbed, ngp))

    igap = 0
    do ib = 1, ex%nbe
      if (ex%gcount(ib) == 0) then
        ne = ex%nee(ib)
        do ie = 1, ne
          a(ie) = ex%ae(ib, ie)
          e(ie) = ex%ee(ib, ie)
        end do

        if (ctl%ismoth > 0) then
          call smooth(ctl%ismoth, a(1:ne), e(1:ne), ne)
          if (ctl%ipr == 2 .and. ne > 0) then
            call write_smoothed(ib, e(1:ne), a(1:ne))
          end if
        end if

        call intpol(a, ne, e, ngp, ctl%eincr, xpl, ypl, too_few)
        if (too_few) call write_too_few_points()
        if (ctl%ipr == 2 .and. ne > 0) then
          call write_interpolated(ib, e(1:ne), a(1:ne))
        end if

        if (ne > 0) then
          call der(a(1:ne), ne, ap(1:ne), ctl%eincr)
          call der(ap(1:ne), ne, app(1:ne), ctl%eincr)
        end if

        ex%nee(ib) = ne
        do ie = 1, ne
          ex%ae(ib, ie) = a(ie)
          ex%ee(ib, ie) = e(ie)
          ex%aep(ib, ie) = ap(ie)
          ex%aepp(ib, ie) = app(ie)
        end do

      else
        if (ctl%ismoth > 0) then
          ! smooth the sections between the gaps
          nz = 0
          do igc = igap + 1, igap + ex%gcount(ib)
            ne = ex%ng1(igc) - (nz + 1)
            do ie = 1, ne
              a(ie) = ex%ae(ib, ie + nz)
              e(ie) = ex%ee(ib, ie + nz)
            end do
            call smooth(ctl%ismoth, a(1:ne), e(1:ne), ne)
            if (ctl%ipr == 2 .and. ne > 0) then
              call write_smoothed(ib, e(1:ne), a(1:ne))
            end if
            do ie = 1, ne
              ex%ae(ib, ie + nz) = a(ie)
              ex%ee(ib, ie + nz) = e(ie)
            end do
            nz = ex%ng2(igc)
          end do
          ne = ex%nee(ib) - (nz + 1)
          do ie = 1, ne
            a(ie) = ex%ae(ib, ie + nz)
            e(ie) = ex%ee(ib, ie + nz)
          end do
          call smooth(ctl%ismoth, a(1:ne), e(1:ne), ne)
          if (ctl%ipr == 2 .and. ne > 0) then
            call write_smoothed(ib, e(1:ne), a(1:ne))
          end if
          do ie = 1, ne
            ex%ae(ib, ie + nz) = a(ie)
            ex%ee(ib, ie + nz) = e(ie)
          end do
        end if

        ! energy working grid
        ne = ex%nee(ib)
        do ie = 1, ne
          a(ie) = ex%ae(ib, ie)
          e(ie) = ex%ee(ib, ie)
        end do
        call grid(a, ne, e, ngp, ctl%eincr, worx, npts, imin, too_few)
        if (too_few) call write_too_few_points()
        ex%nee(ib) = ne

        ! interpolate the sections between the gaps onto the working grid;
        ! all sections of the beam share one interpolation cache
        nz = 0
        newnz = 0
        ycache%il = 0
        ycache%ih = 0
        do igc = igap + 1, igap + ex%gcount(ib)
          if (ex%ng1(igc) > imin) then
            ne1 = ex%ng1(igc) - imin - nz
            do ie = 1, ne1
              xpl(ie) = ex%ee(ib, ie + nz)
              ypl(ie) = ex%ae(ib, ie + nz)
            end do
            nz = ex%ng2(igc)

            ex%ng1(igc) = 0
            ex%ng2(igc) = npts
            do ipts = 1, npts
              if (worx(ipts) <= ex%eg1(igc)) then
                ex%ng1(igc) = ipts
              end if
              if (worx(npts - ipts + 1) >= ex%eg2(igc)) then
                ex%ng2(igc) = npts - ipts + 1
              end if
            end do
            ex%eg1(igc) = worx(ex%ng1(igc))
            ex%eg2(igc) = worx(ex%ng2(igc))

            nlen = ex%ng1(igc) - newnz - 1
            do ilen = 1, nlen
              xval = worx(ilen + newnz)
              a(ilen) = yval(xval, ypl(1:ne1), xpl(1:ne1), ne1, ycache)
            end do

            if (nlen > 0) then
              call der(a(1:nlen), nlen, ap(1:nlen), ctl%eincr)
              call der(ap(1:nlen), nlen, app(1:nlen), ctl%eincr)
            end if

            do ilen = 1, nlen
              ex%ae(ib, ilen + newnz) = a(ilen)
              ex%aep(ib, ilen + newnz) = ap(ilen)
              ex%aepp(ib, ilen + newnz) = app(ilen)
            end do

            ! intensity and derivatives are zero inside the gap
            do ipts = ex%ng1(igc), ex%ng2(igc)
              ex%ae(ib, ipts) = 0.0_sp
              ex%aep(ib, ipts) = 0.0_sp
              ex%aepp(ib, ipts) = 0.0_sp
            end do
            newnz = ex%ng2(igc)

          else if (imin >= ex%ng2(igc)) then
            ! no nonzero intensity before the end of the gap
            ex%ng1(igc) = 0
            ex%ng2(igc) = 0
            ex%eg1(igc) = 0.0_sp
            ex%eg2(igc) = 0.0_sp
          end if
        end do

        ! section after the last gap
        ne1 = ne - imin - nz
        do ie = 1, ne1
          xpl(ie) = ex%ee(ib, ie + nz)
          ypl(ie) = ex%ae(ib, ie + nz)
        end do

        nlen = npts - newnz
        do ilen = 1, nlen
          xval = worx(ilen + newnz)
          a(ilen) = yval(xval, ypl(1:ne1), xpl(1:ne1), ne1, ycache)
        end do

        if (nlen > 0) then
          call der(a(1:nlen), nlen, ap(1:nlen), ctl%eincr)
          call der(ap(1:nlen), nlen, app(1:nlen), ctl%eincr)
        end if

        do ilen = 1, nlen
          ex%ae(ib, ilen + newnz) = a(ilen)
          ex%aep(ib, ilen + newnz) = ap(ilen)
          ex%aepp(ib, ilen + newnz) = app(ilen)
        end do
        ex%nee(ib) = npts
      end if
      igap = igap + ex%gcount(ib)
    end do

    call ypend(ex%ae, ex%aep, nbed, ngp, ex%nbe, ex%nee, ex%ype, ctl%vi)
    if (ctl%ipr >= 2) call write_pendry_y(ex%nbe, ex%nee, ex%ee, ex%ype)
    call ysmooth(ex%ae, ex%aep, nbed, ngp, ex%nbe, ex%nee, ex%yse, ctl%vi, &
                 ctl%eincr)
    if (ctl%ipr >= 2) call write_smooth_y(ex%nbe, ex%nee, ex%ee, ex%yse)

    ! integrals over the whole curves: intensity, squared intensity (R2
    ! only) and squared Y functions (R_P and R_S only)
    do ib = 1, ex%nbe
      ie2 = ex%nee(ib)
      if (ie2 > 0) then
        ex%tse(ib) = 0.0_sp
        ex%tse2(ib) = 0.0_sp
        ex%tseyp2(ib) = 0.0_sp
        ex%tseys2(ib) = 0.0_sp

        call varsum(ex%ae, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ib, 1, 1, &
                    ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 1, ex%tse(ib), ypl)

        if (ctl%wr(1) > 1.e-6_sp) then
          call varsum(ex%ae, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ib, 1, 1, &
                      ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ex%tse2(ib), ypl)
        end if

        if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
          call varsum(ex%ype, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ib, 1, &
                      1, ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ex%tseyp2(ib), &
                      ypl)
          call varsum(ex%yse, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ib, 1, &
                      1, ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ex%tseys2(ib), &
                      ypl)
        end if
      end if
    end do
  end subroutine prepare_experiment
end module experiment_preparation

!> Theory beams: selection of the energy range, domain averaging of all
!! geometries, and the preparation of one geometry for the R factors
!! (interpolation onto the working grid, derivatives, Y functions and
!! integrals).
module theory_preparation
  use kinds, only: sp
  use control_input, only: control_parameters
  use interpolation, only: grid_points, intpol
  use derivatives, only: der, ysmooth
  use rfactor_terms, only: varsum, ypend
  use domain_averaging, only: maxint, rintav
  use theory_input, only: theory_records
  use rfactor_output, only: write_interpolated, write_pendry_y, &
                            write_smooth_y, write_too_few_points
  implicit none

  private

  public :: average_theory, prepare_theory, select_energies, theory_beams, &
            theory_grid_size

  !> Averaged theory beams of one geometry.
  type :: theory_beams
    ! number of averaged beams (nbt) and their parallel momenta (pqav)
    integer :: nbt
    real(sp), allocatable :: pqav(:, :)
    ! number of points (net), intensities (at), energies (et), first and
    ! second derivatives (atp, atpp), Pendry and smooth Y functions (ypt,
    ! yst) of each beam
    integer, allocatable :: net(:)
    real(sp), allocatable :: at(:, :), et(:, :), atp(:, :), atpp(:, :)
    real(sp), allocatable :: ypt(:, :), yst(:, :)
    ! integrals over the intensity and the squared Pendry and smooth Y
    ! functions of each beam
    real(sp), allocatable :: tst(:), tstyp2(:), tstys2(:)
  end type theory_beams

contains

  !> Number of points the arrays of the theory beams need for the neti
  !! energies es and their interpolation onto the working grid (steps of
  !! eincr).
  pure integer function theory_grid_size(es, neti, eincr) result(ngp)
    real(sp), intent(in) :: es(:), eincr
    integer, intent(in) :: neti

    ngp = max(neti, grid_points(es(1:neti), eincr))
  end function theory_grid_size

  !> Takes the energies and intensities of all geometries from the records
  !! tr, ns records per energy: per energy e(i) and geometry its code
  !! spac and the intensities ydata(:, :, i). Energies above emax end the
  !! selection, those below emin are overwritten by the next. ne returns
  !! the number of energies; e and ydata hold one more, the energy that
  !! ended the selection (or the incomplete last one).
  subroutine select_energies(tr, emin, emax, e, ydata, spac, ne)
    type(theory_records), intent(in) :: tr
    real(sp), intent(in) :: emin, emax
    real(sp), allocatable, intent(out) :: e(:), ydata(:, :, :), spac(:)
    integer, intent(out) :: ne

    integer :: i, irec, j

    allocate(e((tr%nrec + tr%ns - 1)/tr%ns), spac(tr%ns))
    allocate(ydata(tr%ns, tr%nbtd, size(e)))

    irec = 0
    i = 1
    energies: do
      do j = 1, tr%ns
        irec = irec + 1
        if (irec > tr%nrec) exit energies
        e(i) = tr%e(irec)
        spac(j) = tr%spac(irec)
        ydata(j, 1:tr%nbtd, i) = tr%y(1:tr%nbtd, irec)
      end do
      if (e(i) > emax) exit energies
      if (e(i) >= emin) i = i + 1
    end do energies
    ne = i - 1
  end subroutine select_energies

  !> Copies the intensities ys(is, :, 1:ne) of geometry is to y and the
  !! energies es to every beam of et.
  subroutine strip2(ys, ns, nbd, ne, is, y, es, et, ngp)
    integer, intent(in) :: ns, nbd, ne, is, ngp
    real(sp), intent(in) :: ys(ns, nbd, ne), es(ne)
    real(sp), intent(inout) :: y(nbd, ngp), et(nbd, ngp)

    integer :: ib, ie

    do ib = 1, nbd
      do ie = 1, ne
        et(ib, ie) = es(ie)
        y(ib, ie) = ys(is, ib, ie)
      end do
    end do
  end subroutine strip2

  !> Averages the beam groups of the nss geometries and neti energies in
  !! ats into atsav and th%pqav. suspect is the first averaged beam whose
  !! maximum intensity am exceeds 1, or 0 if there is none. th is sized for
  !! the averaged beams of one geometry on ngp points.
  subroutine average_theory(ctl, th, tr, nss, neti, ngp, ats, atsav, suspect, &
                            am)
    type(control_parameters), intent(in) :: ctl
    type(theory_beams), intent(out) :: th
    type(theory_records), intent(in) :: tr
    integer, intent(in) :: nss, neti, ngp
    real(sp), contiguous, intent(in) :: ats(:, :, :)
    real(sp), allocatable, intent(out) :: atsav(:, :, :)
    integer, intent(out) :: suspect
    real(sp), intent(out) :: am

    integer :: nbtd, ib

    nbtd = tr%nbtd
    allocate(th%pqav(2, nbtd), th%net(nbtd), th%tst(nbtd), &
             th%tstyp2(nbtd), th%tstys2(nbtd))
    allocate(th%at(nbtd, ngp), th%et(nbtd, ngp), th%atp(nbtd, ngp), &
             th%atpp(nbtd, ngp), th%ypt(nbtd, ngp), th%yst(nbtd, ngp))
    allocate(atsav(nss, nbtd, neti))

    call rintav(ats, nss, nbtd, neti, tr%pq, th%pqav, ctl%kav, tr%sym, &
                th%nbt, atsav)

    suspect = 0
    do ib = 1, th%nbt
      call maxint(atsav, nss, nbtd, ib, neti, am)
      if (am > 1.0_sp) then
        suspect = ib
        return
      end if
    end do
  end subroutine average_theory

  !> Prepares geometry is of the averaged theory atsav (nss geometries,
  !! nbtd beams, neti points on the energies es): interpolation onto the
  !! working grid, derivatives, Y functions and integrals over the whole
  !! curves. Logs the curves after each step as ctl%ipr asks.
  subroutine prepare_theory(th, ctl, atsav, nss, nbtd, neti, is, es)
    type(theory_beams), intent(inout) :: th
    type(control_parameters), intent(in) :: ctl
    real(sp), contiguous, intent(in) :: atsav(:, :, :), es(:)
    integer, intent(in) :: nss, nbtd, neti, is

    real(sp), dimension(size(th%at, 2)) :: a, e, ap, app, xpl, ypl
    logical :: too_few
    integer :: ngp, ib, ie, ie2, ne

    ngp = size(th%at, 2)

    call strip2(atsav, nss, nbtd, neti, is, th%at, es, th%et, ngp)

    do ib = 1, th%nbt
      do ie = 1, neti
        a(ie) = th%at(ib, ie)
        e(ie) = th%et(ib, ie)
      end do

      ne = neti
      call intpol(a, ne, e, ngp, ctl%eincr, xpl, ypl, too_few)
      if (too_few) call write_too_few_points()
      if (ctl%ipr == 2 .and. ne > 0) then
        call write_interpolated(ib, e(1:ne), a(1:ne))
      end if

      if (ne > 0) then
        call der(a(1:ne), ne, ap(1:ne), ctl%eincr)
        call der(ap(1:ne), ne, app(1:ne), ctl%eincr)
      end if

      th%net(ib) = ne
      do ie = 1, ne
        th%at(ib, ie) = a(ie)
        th%et(ib, ie) = e(ie)
        th%atp(ib, ie) = ap(ie)
        th%atpp(ib, ie) = app(ie)
      end do
    end do

    call ypend(th%at, th%atp, nbtd, ngp, th%nbt, th%net, th%ypt, ctl%vi)
    if (ctl%ipr >= 2) call write_pendry_y(th%nbt, th%net, th%et, th%ypt)
    call ysmooth(th%at, th%atp, nbtd, ngp, th%nbt, th%net, th%yst, ctl%vi, &
                 ctl%eincr)
    if (ctl%ipr >= 2) call write_smooth_y(th%nbt, th%net, th%et, th%yst)

    ! integrals over the whole curves: intensity and squared Y functions
    do ib = 1, th%nbt
      ie2 = th%net(ib)
      if (ie2 > 0) then
        call varsum(th%at, th%at, th%at, th%at, nbtd, nbtd, ngp, ib, 1, 1, &
                    ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 1, th%tst(ib), ypl)
        call varsum(th%ypt, th%at, th%at, th%at, nbtd, nbtd, ngp, ib, 1, 1, &
                    ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, th%tstyp2(ib), ypl)
        call varsum(th%yst, th%at, th%at, th%at, nbtd, nbtd, ngp, ib, 1, 1, &
                    ie2, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, th%tstys2(ib), ypl)
      end if
    end do
  end subroutine prepare_theory
end module theory_preparation

!> Scan of the inner potential shift for one geometry: R factors per beam,
!! their averages over all, integer and fractional order beams, and the
!! best fit over all geometries and shifts.
module rfactor_scan
  use kinds, only: sp
  use control_input, only: control_parameters, nrfact
  use experiment_input, only: experiment_beams
  use theory_preparation, only: theory_beams
  use common_interval, only: comnei
  use rfactor_terms, only: epszj, varsum
  use rfactor_output, only: write_averages, write_beam_header, &
                            write_beam_rfactors, write_no_overlap, &
                            write_normalization, write_scan_end, &
                            write_v0_shift
  implicit none

  private

  public :: beam_results, best_fit, scan_v0, start_scan

  !> Per-beam results of the latest inner potential shift: grid indices of
  !! the common energy interval in experiment (nemin, nemax) and theory
  !! (ntmin, ntmax), theory grid indices of the gaps (ngt1, ngt2),
  !! integrals over experiment (sest) and theory (stst) in the common
  !! interval, and the weighted average R factor of each beam (ravb).
  type :: beam_results
    integer, allocatable :: nemin(:), nemax(:), ntmin(:), ntmax(:)
    integer, allocatable :: ngt1(:), ngt2(:)
    real(sp), allocatable :: sest(:), stst(:), ravb(:)
  end type beam_results

  !> Best fit: average R factor (rav), geometry index, inner potential
  !! shift (v0), ratio of fractional to integer order intensities in
  !! experiment and theory (cife, cift), and copies of the beam_results,
  !! experimental gap indices (ng1, ng2) and theory curves (et, at) of
  !! that fit.
  type :: best_fit
    real(sp) :: rav, v0, cife, cift
    integer :: geometry
    integer, allocatable :: nemin(:), nemax(:), ntmin(:), ntmax(:)
    integer, allocatable :: ng1(:), ng2(:), ngt1(:), ngt2(:)
    real(sp), allocatable :: sest(:), stst(:), ravb(:)
    real(sp), allocatable :: et(:, :), at(:, :)
  end type best_fit

  !> Sums over the beams for one inner potential shift. Index 1:nrfact of
  !! rf is R2, R_ZJ, R_P, R_S averaged over all beams, rf_m(:, kmit) the
  !! same over the integer (kmit = 1) or fractional (kmit = 2) order beams;
  !! rp_num/rp_den and rs_num/rs_den are numerator and denominator of R_P
  !! and R_S (_m: per kmit); erang is the total common energy range; sem
  !! and stm the experimental and theory integrals per kmit; rav and ravm
  !! the weighted averages; cife and cift the fractional to integer order
  !! intensity ratios of experiment and theory.
  type :: rfactor_sums
    real(sp) :: rf(nrfact) = 0.0_sp, rf_m(nrfact, 2) = 0.0_sp
    real(sp) :: rp_num = 0.0_sp, rp_den = 0.0_sp
    real(sp) :: rs_num = 0.0_sp, rs_den = 0.0_sp
    real(sp) :: rp_num_m(2) = 0.0_sp, rp_den_m(2) = 0.0_sp
    real(sp) :: rs_num_m(2) = 0.0_sp, rs_den_m(2) = 0.0_sp
    real(sp) :: erang = 0.0_sp, erangm(2) = 0.0_sp
    real(sp) :: sem(2) = 0.0_sp, stm(2) = 0.0_sp
    real(sp) :: rav = 0.0_sp, ravm(2) = 0.0_sp
    real(sp) :: cife = 0.0_sp, cift = 0.0_sp
  end type rfactor_sums

contains

  !> Sizes the per-beam results and the best fit for nbed experimental
  !! and nbtd theory beams, ngp grid points and ngap gaps.
  subroutine start_scan(res, best, nbed, nbtd, ngp, ngap)
    type(beam_results), intent(out) :: res
    type(best_fit), intent(out) :: best
    integer, intent(in) :: nbed, nbtd, ngp, ngap

    allocate(res%nemin(nbed), res%nemax(nbed), res%sest(nbed), res%ravb(nbed))
    allocate(res%ntmin(nbtd), res%ntmax(nbtd), res%stst(nbtd))
    allocate(res%ngt1(ngap), res%ngt2(ngap))

    allocate(best%nemin(nbed), best%nemax(nbed), best%sest(nbed), &
             best%ravb(nbed))
    allocate(best%ntmin(nbtd), best%ntmax(nbtd), best%stst(nbtd))
    allocate(best%et(nbtd, ngp), best%at(nbtd, ngp))
    allocate(best%ng1(ngap), best%ng2(ngap), best%ngt1(ngap), &
             best%ngt2(ngap))
    best%rav = 1.e2_sp
  end subroutine start_scan

  !> Scans the inner potential shift v0 = v01, v01 + vincr, ... up to v02
  !! for geometry is (code d12) of the prepared theory th (nbtd beams as
  !! read) against the experiment ex, and updates best. Writes the beam
  !! R factors and averages to stdout and ROUT (unit 7), and the best
  !! average R factor of the geometry to ROUTSHORT (unit 9).
  subroutine scan_v0(ctl, ex, th, nbtd, is, d12, res, best)
    type(control_parameters), intent(in) :: ctl
    type(experiment_beams), intent(in) :: ex
    type(theory_beams), intent(in) :: th
    integer, intent(in) :: nbtd, is
    real(sp), intent(in) :: d12
    type(beam_results), intent(inout) :: res
    type(best_fit), intent(inout) :: best

    type(rfactor_sums) :: sums
    real(sp) :: v0, v0r, bgrav, c, eet, rf(nrfact)
    logical :: overlap
    integer :: nv0, iv, ibe, ibt, igap

    bgrav = 1.e2_sp
    nv0 = int((ctl%v02 - ctl%v01)/ctl%vincr + 0.0001_sp) + 1
    v0 = ctl%v01

    do iv = 1, nv0
      call write_v0_shift(v0)
      v0r = -ctl%v0rr + v0

      sums = rfactor_sums()
      igap = 0
      do ibe = 1, ex%nbe
        ibt = ctl%ibp(ibe)
        call write_beam_header(ibe, ex%bename(:, ibe), ibt, th%pqav(:, ibt))
        call beam_rfactor(ctl, ex, th, nbtd, ibe, v0, igap, res, sums, &
                          overlap, c, eet, rf)
        if (overlap) then
          call write_normalization(c)
          call write_beam_rfactors(v0, eet, rf, res%ravb(ibe), ibe, &
                                   ex%bename(:, ibe), d12, v0r, ctl%emin, &
                                   ctl%emax)
        else
          call write_no_overlap(ibe, ex%bename(:, ibe), d12, v0r, ctl%emin, &
                                ctl%emax)
        end if
      end do

      call average_over_beams(ctl, sums)
      call write_averages(v0, d12, v0r, ctl%emin, ctl%emax, ctl%nint, &
                          ctl%nfrac, sums%erangm, sums%rf_m, sums%ravm, &
                          sums%erang, sums%rf, sums%rav, sums%cife, sums%cift)

      if (sums%rav < bgrav) then
        bgrav = sums%rav
      end if
      if (sums%rav < best%rav) then
        call keep_best(ctl, ex, th, res, sums, is, v0, best)
      end if

      v0 = v0 + ctl%vincr
    end do

    call write_scan_end(bgrav)
  end subroutine scan_v0

  !> R factors rf (R2, R_ZJ, R_P, R_S) of experimental beam ibe and its
  !! theory beam at inner potential shift v0, added to sums. Integrals run
  !! over the energy interval common to both, of length eet without the
  !! gaps; c is the normalisation of the experiment to the theory. If the
  !! interval is empty (no overlap), only res is set. igap counts the gaps
  !! of the previous beams and is not advanced for a beam without overlap.
  subroutine beam_rfactor(ctl, ex, th, nbtd, ibe, v0, igap, res, sums, &
                          overlap, c, eet, rf)
    type(control_parameters), intent(in) :: ctl
    type(experiment_beams), intent(in) :: ex
    type(theory_beams), intent(in) :: th
    integer, intent(in) :: nbtd, ibe
    real(sp), intent(in) :: v0
    integer, intent(inout) :: igap
    type(beam_results), intent(inout) :: res
    type(rfactor_sums), intent(inout) :: sums
    logical, intent(out) :: overlap
    real(sp), intent(out) :: c, eet, rf(nrfact)

    ! epse: maximum of the derivative of the experiment; ss*, su*, sg*:
    ! integrals below, above and in the gaps of the common interval; se*,
    ! st*: integrals of experiment and theory over it (2: squared
    ! intensity, yp2/ys2: squared Pendry/smooth Y function); s2, srzj,
    ! syp2, sys2: integrands of R2, R_ZJ, R_P, R_S (sg2, sgzj: in the gaps)
    real(sp) :: epse
    real(sp) :: ss, su, ss2, su2, ssyp2, suyp2, ssys2, suys2, sg, sgyp2, sgys2
    real(sp) :: se, st, se2, seyp2, styp2, seys2, stys2
    real(sp) :: s2, sg2, srzj, sgzj, syp2, sys2
    real(sp) :: r2, rrzj, rpe, rsm
    real(sp) :: ypl(size(ex%ae, 2))
    integer :: nbed, ngp, ibt, kmit, ne1, ne2, nt1, nt2, neet, ne, nt, igc

    nbed = size(ex%ae, 1)
    ngp = size(ex%ae, 2)

    r2 = 0.0_sp
    rrzj = 0.0_sp
    rpe = 0.0_sp
    rsm = 0.0_sp
    res%ravb(ibe) = 0.0_sp
    res%nemin(ibe) = 1
    res%nemax(ibe) = ngp
    res%sest(ibe) = 1.0_sp

    ibt = ctl%ibp(ibe)

    kmit = ctl%mittel(ibt)

    ! energy interval of length eet common to experiment and theory,
    ! bounded by the grid points (ne1, ne2) and (nt1, nt2)
    call comnei(ex%ee, nbed, ex%nee, th%et, nbtd, ngp, th%net, ibe, ibt, v0, &
                ctl%eincr, ne1, ne2, nt1, nt2)

    if (ne2 > ne1) then
      neet = nt1 - ne1
      eet = ex%ee(ibe, ne2) - ex%ee(ibe, ne1)

      ! theory grid points of the gaps; eet is reduced by the gaps inside
      ! the common interval
      if (ex%gcount(ibe) /= 0) then
        do igc = igap + 1, igap + ex%gcount(ibe)
          if (ex%ng1(igc) <= ne1) then
            if (ex%ng2(igc) < ne2) then
              res%ngt1(igc) = 0
              res%ngt2(igc) = 0
              ne1 = max(ne1, ex%ng1(igc))
              nt1 = ne1 + neet
              eet = ex%ee(ibe, ne2) - ex%ee(ibe, ne1)
            else
              eet = 0.0_sp
            end if
          else if (ex%ng1(igc) <= ne2) then
            if (ex%ng2(igc) <= ne2) then
              res%ngt1(igc) = ex%ng1(igc) + neet
              res%ngt2(igc) = ex%ng2(igc) + neet
              eet = eet - (ex%ng2(igc) - ex%ng1(igc)) * ctl%eincr
            else
              res%ngt1(igc) = 0
              res%ngt2(igc) = 0
              ne2 = min(ne2, ex%ng2(igc))
              nt2 = ne2 + neet
              eet = ex%ee(ibe, ne2) - ex%ee(ibe, ne1)
            end if
          else
            res%ngt1(igc) = 0
            res%ngt2(igc) = 0
          end if
        end do
      end if
    else
      eet = 0.0_sp
    end if

    res%ntmin(ibt) = nt1
    res%ntmax(ibt) = nt2
    res%nemin(ibe) = ne1
    res%nemax(ibe) = ne2

    overlap = eet /= 0.0_sp
    if (.not. overlap) return

    ! maximum |derivative| of the experiment, for R_ZJ
    call epszj(ex%aep, nbed, ngp, ibe, ne1, ne2, epse)

    ! integrals of the experiment over the common interval: whole curve
    ! minus the parts below and above
    ne = ex%nee(ibe)
    nt = th%net(ibt)

    ss = 0.0_sp
    ss2 = 0.0_sp
    ssyp2 = 0.0_sp
    ssys2 = 0.0_sp
    su = 0.0_sp
    su2 = 0.0_sp
    suyp2 = 0.0_sp
    suys2 = 0.0_sp

    if (ne1 > 1) then
      call varsum(ex%ae, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, 1, &
                  ne1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 1, ss, ypl)
      if (ctl%wr(1) > 1.e-6_sp) then
        call varsum(ex%ae, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, 1, &
                    ne1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ss2, ypl)
      end if
      if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
        call varsum(ex%ype, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, 1, &
                    ne1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ssyp2, ypl)
        call varsum(ex%yse, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, 1, &
                    ne1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ssys2, ypl)
      end if
    end if

    if (ne2 < ne) then
      call varsum(ex%ae, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, ne2, &
                  ne, 0, ctl%eincr, 0.0_sp, 0.0_sp, 1, su, ypl)
      if (ctl%wr(1) > 1.e-6_sp) then
        call varsum(ex%ae, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, ne2, &
                    ne, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, su2, ypl)
      end if
      if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
        call varsum(ex%ype, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, &
                    ne2, ne, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, suyp2, ypl)
        call varsum(ex%yse, ex%ae, ex%ae, ex%ae, nbed, nbed, ngp, ibe, 1, &
                    ne2, ne, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, suys2, ypl)
      end if
    end if

    se = ex%tse(ibe) - ss - su
    se2 = ex%tse2(ibe) - ss2 - su2
    seyp2 = ex%tseyp2(ibe) - ssyp2 - suyp2
    seys2 = ex%tseys2(ibe) - ssys2 - suys2
    res%sest(ibe) = se

    ! integrals of the theory over the common interval: whole curve minus
    ! the parts below, above and in the gaps
    ss = 0.0_sp
    ssyp2 = 0.0_sp
    ssys2 = 0.0_sp
    sg = 0.0_sp
    sgyp2 = 0.0_sp
    sgys2 = 0.0_sp
    su = 0.0_sp
    suyp2 = 0.0_sp
    suys2 = 0.0_sp

    if (nt1 /= 1) then
      call varsum(th%at, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, 1, &
                  nt1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 1, ss, ypl)
      if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
        call varsum(th%ypt, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, 1, &
                    nt1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ssyp2, ypl)
        call varsum(th%yst, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, 1, &
                    nt1, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, ssys2, ypl)
      end if
    end if

    if (ex%gcount(ibe) /= 0) then
      do igc = igap + 1, igap + ex%gcount(ibe)
        call varsum(th%at, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, &
                    res%ngt1(igc), res%ngt2(igc), 0, ctl%eincr, 0.0_sp, &
                    0.0_sp, 1, sg, ypl)
        if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
          call varsum(th%ypt, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, &
                      res%ngt1(igc), res%ngt2(igc), 0, ctl%eincr, 0.0_sp, &
                      0.0_sp, 2, sgyp2, ypl)
          call varsum(th%yst, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, &
                      res%ngt1(igc), res%ngt2(igc), 0, ctl%eincr, 0.0_sp, &
                      0.0_sp, 2, sgys2, ypl)
        end if
      end do
    end if

    if (nt2 /= nt) then
      call varsum(th%at, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, nt2, &
                  nt, 0, ctl%eincr, 0.0_sp, 0.0_sp, 1, su, ypl)
      if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
        call varsum(th%ypt, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, &
                    nt2, nt, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, suyp2, ypl)
        call varsum(th%yst, th%at, th%at, th%at, nbtd, nbtd, ngp, ibt, 1, &
                    nt2, nt, 0, ctl%eincr, 0.0_sp, 0.0_sp, 2, suys2, ypl)
      end if
    end if

    st = th%tst(ibt) - ss - sg - su
    styp2 = th%tstyp2(ibt) - ssyp2 - sgyp2 - suyp2
    stys2 = th%tstys2(ibt) - ssys2 - sgys2 - suys2
    res%stst(ibt) = st

    c = se / st

    if (kmit /= 0) then
      sums%sem(kmit) = sums%sem(kmit) + se
      sums%stm(kmit) = sums%stm(kmit) + st
    end if

    ! integrals over experiment and theory together, without the gaps
    s2 = 0.0_sp
    sg2 = 0.0_sp
    srzj = 0.0_sp
    sgzj = 0.0_sp
    syp2 = 0.0_sp
    sys2 = 0.0_sp
    sgyp2 = 0.0_sp
    sgys2 = 0.0_sp

    if (ctl%wr(1) > 1.e-6_sp) then
      call varsum(ex%ae, th%at, ex%ae, ex%ae, nbed, nbtd, ngp, ibe, ibt, ne1, &
                  ne2, neet, ctl%eincr, 0.0_sp, c, 3, s2, ypl)
      if (ex%gcount(ibe) /= 0) then
        do igc = igap + 1, igap + ex%gcount(ibe)
          call varsum(ex%ae, th%at, ex%ae, ex%ae, nbed, nbtd, ngp, ibe, ibt, &
                      ex%ng1(igc), ex%ng2(igc), neet, ctl%eincr, 0.0_sp, c, &
                      3, sg2, ypl)
        end do
      end if
    end if

    if (ctl%wr(2) > 1.e-6_sp) then
      call varsum(ex%aep, th%atp, ex%aepp, th%atpp, nbed, nbtd, ngp, ibe, ibt, &
                  ne1, ne2, neet, ctl%eincr, epse, c, 4, srzj, ypl)
      if (ex%gcount(ibe) /= 0) then
        do igc = igap + 1, igap + ex%gcount(ibe)
          call varsum(ex%aep, th%atp, ex%aepp, th%atpp, nbed, nbtd, ngp, ibe, &
                      ibt, ex%ng1(igc), ex%ng2(igc), neet, ctl%eincr, epse, &
                      c, 4, sgzj, ypl)
        end do
      end if
    end if

    if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
      call varsum(ex%ype, th%ypt, ex%ype, ex%ype, nbed, nbtd, ngp, ibe, ibt, &
                  ne1, ne2, neet, ctl%eincr, 0.0_sp, 1.0_sp, 3, syp2, ypl)
      call varsum(ex%yse, th%yst, ex%yse, ex%yse, nbed, nbtd, ngp, ibe, ibt, &
                  ne1, ne2, neet, ctl%eincr, 0.0_sp, 1.0_sp, 3, sys2, ypl)
      if (ex%gcount(ibe) /= 0) then
        do igc = igap + 1, igap + ex%gcount(ibe)
          call varsum(ex%ype, th%ypt, ex%ype, ex%ype, nbed, nbtd, ngp, ibe, &
                      ibt, ex%ng1(igc), ex%ng2(igc), neet, ctl%eincr, 0.0_sp, &
                      1.0_sp, 3, sgyp2, ypl)
          call varsum(ex%yse, th%yst, ex%yse, ex%yse, nbed, nbtd, ngp, ibe, &
                      ibt, ex%ng1(igc), ex%ng2(igc), neet, ctl%eincr, 0.0_sp, &
                      1.0_sp, 3, sgys2, ypl)
        end do
      end if
    end if

    s2 = s2 - sg2
    srzj = srzj - sgzj
    syp2 = syp2 - sgyp2
    sys2 = sys2 - sgys2

    ! R factors, normalised to about 1 for anticorrelated curves (sin**2
    ! against cos**2 over one period); the averages over the beams are
    ! weighted by the beam weight times the common energy range
    eet = ctl%wb(ibe) * eet

    ! R2 from the integral of (experiment - c*theory)**2
    if (ctl%wr(1) > 1.e-6_sp) then
      r2 = 0.5_sp * s2 / se2
      sums%rf(1) = sums%rf(1) + eet*r2
      if (kmit /= 0) then
        sums%rf_m(1, kmit) = sums%rf_m(1, kmit) + eet*r2
      end if
    end if

    ! reduced Zanazzi-Jona R factor (times 0.5)
    if (ctl%wr(2) > 1.e-6_sp) then
      rrzj = srzj / (0.027_sp*se)
      sums%rf(2) = sums%rf(2) + eet*rrzj
      if (kmit /= 0) then
        sums%rf_m(2, kmit) = sums%rf_m(2, kmit) + eet*rrzj
      end if
    end if

    ! Pendry and smooth R factors; the averages over the beams are ratios
    ! of the weighted sums of numerators and denominators
    if (ctl%wr(3) > 1.e-6_sp .or. ctl%wr(4) > 1.e-6_sp) then
      rpe = syp2 / (seyp2 + styp2)
      rsm = sys2 / (seys2 + stys2)
      sums%rp_num = sums%rp_num + ctl%wb(ibe) * syp2
      sums%rp_den = sums%rp_den + ctl%wb(ibe) * (seyp2 + styp2)
      sums%rs_num = sums%rs_num + ctl%wb(ibe) * sys2
      sums%rs_den = sums%rs_den + ctl%wb(ibe) * (seys2 + stys2)
      if (kmit /= 0) then
        sums%rp_num_m(kmit) = sums%rp_num_m(kmit) + ctl%wb(ibe) * syp2
        sums%rp_den_m(kmit) = sums%rp_den_m(kmit) &
                              + ctl%wb(ibe) * (seyp2 + styp2)
        sums%rs_num_m(kmit) = sums%rs_num_m(kmit) + ctl%wb(ibe) * sys2
        sums%rs_den_m(kmit) = sums%rs_den_m(kmit) &
                              + ctl%wb(ibe) * (seys2 + stys2)
      end if
    end if

    ! weighted average of the R factors of this beam (read by viperleed)
    res%ravb(ibe) = (ctl%wr(1)*r2 + ctl%wr(2)*rrzj + &
                     ctl%wr(3)*rpe + ctl%wr(4)*rsm)/ctl%wrsum
    sums%erang = sums%erang + eet
    if (kmit /= 0) then
      sums%erangm(kmit) = sums%erangm(kmit) + eet
    end if

    eet = eet / ctl%wb(ibe)
    rf(1) = r2
    rf(2) = rrzj
    rf(3) = rpe
    rf(4) = rsm
    igap = igap + ex%gcount(ibe)
  end subroutine beam_rfactor

  !> Averages the beam sums of one inner potential shift into R factors
  !! over all, integer and fractional order beams.
  subroutine average_over_beams(ctl, sums)
    type(control_parameters), intent(in) :: ctl
    type(rfactor_sums), intent(inout) :: sums

    integer :: i, imit

    do i = 1, 2
      sums%rf(i) = sums%rf(i)/sums%erang
      do imit = 1, 2
        if (sums%erangm(imit) /= 0.0_sp) then
          sums%rf_m(i, imit) = sums%rf_m(i, imit)/sums%erangm(imit)
        end if
      end do
    end do

    if (ctl%wr(3) > 1e-6_sp .or. ctl%wr(4) > 1e-6_sp) then
      sums%rf(3) = sums%rp_num/sums%rp_den
      sums%rf(4) = sums%rs_num/sums%rs_den
      do imit = 1, 2
        if (sums%rp_den_m(imit) /= 0.0_sp) then
          sums%rf_m(3, imit) = sums%rp_num_m(imit)/sums%rp_den_m(imit)
        end if
        if (sums%rs_den_m(imit) /= 0.0_sp) then
          sums%rf_m(4, imit) = sums%rs_num_m(imit)/sums%rs_den_m(imit)
        end if
      end do
    end if

    sums%rav = (ctl%wr(1)*sums%rf(1) + ctl%wr(2)*sums%rf(2) + &
                ctl%wr(3)*sums%rf(3) + ctl%wr(4)*sums%rf(4)) / ctl%wrsum
    do imit = 1, 2
      sums%ravm(imit) = (ctl%wr(1)*sums%rf_m(1, imit) + &
                         ctl%wr(2)*sums%rf_m(2, imit) + &
                         ctl%wr(3)*sums%rf_m(3, imit) + &
                         ctl%wr(4)*sums%rf_m(4, imit)) / ctl%wrsum
    end do

    ! ratio of fractional to integer order intensities
    if (ctl%nint /= 0 .and. ctl%nfrac /= 0) then
      sums%cife = (sums%sem(2) * real(ctl%nint, sp)) / &
                  (sums%sem(1) * real(ctl%nfrac, sp))
      sums%cift = (sums%stm(2) * real(ctl%nint, sp)) / &
                  (sums%stm(1) * real(ctl%nfrac, sp))
    end if
  end subroutine average_over_beams

  !> Makes the fit of geometry is at shift v0 the best fit.
  subroutine keep_best(ctl, ex, th, res, sums, is, v0, best)
    type(control_parameters), intent(in) :: ctl
    type(experiment_beams), intent(in) :: ex
    type(theory_beams), intent(in) :: th
    type(beam_results), intent(in) :: res
    type(rfactor_sums), intent(in) :: sums
    integer, intent(in) :: is
    real(sp), intent(in) :: v0
    type(best_fit), intent(inout) :: best

    integer :: igap, ibe, ibt, igc, iet

    best%rav = sums%rav
    best%geometry = is
    best%v0 = v0

    igap = 0
    do ibe = 1, ex%nbe
      ibt = ctl%ibp(ibe)

      if (ex%gcount(ibe) /= 0) then
        do igc = igap + 1, igap + ex%gcount(ibe)
          best%ng1(igc) = ex%ng1(igc)
          best%ng2(igc) = ex%ng2(igc)
          best%ngt1(igc) = res%ngt1(igc)
          best%ngt2(igc) = res%ngt2(igc)
        end do
      end if
      igap = igap + ex%gcount(ibe)

      best%nemin(ibe) = res%nemin(ibe)
      best%nemax(ibe) = res%nemax(ibe)
      best%sest(ibe) = res%sest(ibe)
      best%ravb(ibe) = res%ravb(ibe)

      best%ntmin(ibt) = res%ntmin(ibt)
      best%ntmax(ibt) = res%ntmax(ibt)
      do iet = 1, th%net(ibt)
        best%at(ibt, iet) = th%at(ibt, iet)
        best%et(ibt, iet) = th%et(ibt, iet)
      end do
      best%stst(ibt) = res%stst(ibt)
    end do

    if (ctl%nfrac /= 0 .and. ctl%nint /= 0) then
      best%cife = sums%cife
      best%cift = sums%cift
    end if
  end subroutine keep_best
end module rfactor_scan
