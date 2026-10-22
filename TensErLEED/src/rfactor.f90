!  VAN HOVE R-FACTOR PROGRAM  /  MODIFIED BY G. BESOLD (JAN. '82)
!                                ADDITION       (JAN.  '87)              260187
!  MODIFIED BY RD: AVERAGE-R-FACTORS EXCHANGED                           141290
!  Allocation error in subroutine RINTAV removed
!  via introduction of a second field ATAV.       L.Hammer Oct. 2018
!  OUTPUT expanded for error curve calculations   L.Hammer Mar. 2021
!  INPUT formats expanded to allow for 9999 exp. and theo. beams  L.Hammer Mar. 2021
!  R_Smooth added (see J. Phys.: Condens. Matter 38, 105001).     M.Schmid Mar. 2026
!
!##############################################################################
!  R-FACTOR PROGRAM FOR COMPARING LEED EXPERIMENT AND THEORY
!##############################################################################

!  Four different R-factors (the R2-factor, the Zanazzi-Jona R-factor, Pendry's
!  R factor and the Smooth R factor), together with their average can be produced.
!  Experimental IV-curves from different experiments can be averaged together
!  (the same energy values, but no energy limits are required).
!  The input energies may be in any order.
!  It is possible to specify gaps in the experimental IV-curves. These energy
!  ranges will be ignored in the R-factor calculation. The gaps will also be
!  taken into account in the processing of the experimenta data.
!  The experiment can be smoothed before interpolation onto a regular energy
!  grid.
!  The theoretical IV-curves (read in with the format used in the van Hove-Tong
!  LEED programs) for related beams can be averaged together to account for
!  different beams.
!  The averaged theoretical beams can be labelled as integer order or fractio-
!  nal order beam. In this case the average for all integer order and/or all
!  fractional order beams is calculated and written to output (see MITTEL).

!******************************************************************************
!  Input formats (see description of individual quantities further down)
!******************************************************************************

!  - &NL1 EMIN=...,...,&END (some compilers expect a symbol different from &)
!  - KAV (400(25I4,/))
!  - IBP (400(25I4,/))
!  - MITTEL (400(25I4,/))
!  - $NL2 NSSK=...,...,WB=...,..., $END
!  - TEXT (describing experiment) (CHARACTER*80)
!  - NBEA (400(25I4,/))
!  - FORMAT (for exp energies and intensities) (19A4)
!  - BENAME(IB) (5A4)
!  - NEE(IB),FAC (1I4,1E13.4)
!  - EE(IB,IE),AE(IB,IE) (FORMAT)
!  - repeat last step until NEE(IB) pairs of values have been read in
!  - read BENAME(IB+1),NEE(IB+1),... until NBED beams have been read in
!  - TEXT (describing theory) (CHARACTER*80)
!  - NBTD (1I3)
!  - LAB(IB),PQ(IB),SYM(IB) (1I3,2F10.5,1I3) (LAB not used)
!  - repeat last step until NBTD beams have been read in
!  - read energy, geometry and theoretical intensities for all beams
!    (van Hove-output format)
!  - repeat last step until NS geometries have been read in
!  - repeat last two steps until an END-OF-FILE is encountered

!******************************************************************************
!  Limitations
!******************************************************************************

!  Energies to be input in eV, intensities in any unit (may be beam-dependent),
!  cf purpose of quantity FAC in subroutine READE.

program rfactor
  use kinds, only: sp
  use control_input, only: control_parameters, read_control
  use experiment_preparation, only: experiment_beams, experiment_grid_size, &
                                    prepare_experiment, read_experiment
  use theory_input, only: read_theory_records, select_energies, &
                          theory_records
  use theory_preparation, only: average_theory, prepare_theory, &
                                theory_beams, theory_grid_size
  use rfactor_scan, only: beam_results, best_fit, report_best_fit, &
                          scan_v0, start_scan

  implicit none

  ! tr: the theory as read; nss: number of geometries; neti: theory
  ! energies selected; ngp: points of the energy working grid
  ! es, ats, atsav: theory energies and intensities of all geometries as
  !         selected and after averaging the beam groups
  ! values: code of each geometry
  type(control_parameters) :: ctl
  type(theory_records) :: tr
  type(experiment_beams) :: ex
  type(theory_beams) :: th
  type(beam_results) :: res
  type(best_fit) :: best
  integer :: nss, neti, ngp, is
  real(sp), allocatable :: es(:), ats(:, :, :), atsav(:, :, :), values(:)

30 format('*STRUCTURE NO.', 1I5, ' CHARACTERIZED BY VALUE', 1F7.4)

  open(7, file='ROUT')
  open(9, file='ROUTSHORT')
  open(8, file='WEXPEL', status='OLD')
  rewind(8)

  call read_theory_records(tr)
  call read_control(ctl, tr%nbtd, tr%ns, nss)
  call read_experiment(ex, ctl)
  call select_energies(tr, ctl%emin, ctl%emax, es, ats, values, neti)
  deallocate(tr%e, tr%spac, tr%y)

  ngp = max(experiment_grid_size(ex, ctl%eincr), &
            theory_grid_size(es, neti, ctl%eincr))
  call prepare_experiment(ex, ctl, ngp)
  close(8)

  call average_theory(ctl, th, tr, nss, neti, ngp, ats, values, atsav)
  deallocate(ats)

  call start_scan(res, best, size(ctl%nbea), tr%nbtd, ngp, ex%ngap)
  do is = 1, nss
    write(6, 30) is, values(is)
    call prepare_theory(th, ctl, atsav, nss, tr%nbtd, neti, is, es)
    call scan_v0(ctl, ex, th, tr%nbtd, is, values(is), res, best)
  end do

  call report_best_fit(ctl, ex, best, tr%nbtd, tr%title, values)

  stop

end program rfactor
