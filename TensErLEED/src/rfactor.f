C  VAN HOVE R-FACTOR PROGRAM  /  MODIFIED BY G. BESOLD (JAN. '82)
C                                ADDITION       (JAN.  '87)              260187
C  MODIFIED BY RD: AVERAGE-R-FACTORS EXCHANGED                           141290
C  Allocation error in subroutine RINTAV removed
C  via introduction of a second field ATAV.       L.Hammer Oct. 2018
C  OUTPUT expanded for error curve calculations   L.Hammer Mar. 2021
C  INPUT formats expanded to allow for 9999 exp. and theo. beams  L.Hammer Mar. 2021
C  R_Smooth added (see J. Phys.: Condens. Matter 38, 105001).     M.Schmid Mar. 2026
C
C##############################################################################
C  R-FACTOR PROGRAM FOR COMPARING LEED EXPERIMENT AND THEORY
C##############################################################################

C  Four different R-factors (the R2-factor, the Zanazzi-Jona R-factor, Pendry's
C  R factor and the Smooth R factor), together with their average can be produced.
C  Experimental IV-curves from different experiments can be averaged together
C  (the same energy values, but no energy limits are required).
C  The input energies may be in any order.
C  It is possible to specify gaps in the experimental IV-curves. These energy
C  ranges will be ignored in the R-factor calculation. The gaps will also be
C  taken into account in the processing of the experimenta data.
C  The experiment can be smoothed before interpolation onto a regular energy
C  grid.
C  The theoretical IV-curves (read in with the format used in the van Hove-Tong
C  LEED programs) for related beams can be averaged together to account for
C  different beams.
C  The averaged theoretical beams can be labelled as integer order or fractio-
C  nal order beam. In this case the average for all integer order and/or all
C  fractional order beams is calculated and written to output (see MITTEL).

C******************************************************************************
C  Input formats (see description of individual quantities further down)
C******************************************************************************

C  - &NL1 EMIN=...,...,&END (some compilers expect a symbol different from &)
C  - KAV (400(25I4,/))
C  - IBP (400(25I4,/))
C  - MITTEL (400(25I4,/))
C  - $NL2 NSSK=...,...,WB=...,..., $END
C  - TEXT (describing experiment) (CHARACTER*80)
C  - NBEA (400(25I4,/))
C  - FORMAT (for exp energies and intensities) (19A4)
C  - BENAME(IB) (5A4)
C  - NEE(IB),FAC (1I4,1E13.4)
C  - EE(IB,IE),AE(IB,IE) (FORMAT)
C  - repeat last step until NEE(IB) pairs of values have been read in
C  - read BENAME(IB+1),NEE(IB+1),... until NBED beams have been read in
C  - TEXT (describing theory) (CHARACTER*80)
C  - NBTD (1I3)
C  - LAB(IB),PQ(IB),SYM(IB) (1I3,2F10.5,1I3) (LAB not used)
C  - repeat last step until NBTD beams have been read in
C  - read energy, geometry and theoretical intensities for all beams
C    (van Hove-output format)
C  - repeat last step until NS geometries have been read in
C  - repeat last two steps until an END-OF-FILE is encountered

C******************************************************************************
C  Limitations
C******************************************************************************

C  Energies to be input in eV, intensities in any unit (may be beam-dependent),
C  cf purpose of quantity FAC in subroutine READE.


      program rfactor

      use kinds, only: sp
      use control_input, only: control_parameters, read_control
      use experiment_preparation, only: experiment_beams,
     &    prepare_experiment, read_experiment
      use theory_preparation, only: prepare_theory, read_theory,
     &    theory_beams
      use rfactor_scan, only: beam_results, best_fit, report_best_fit,
     &    scan_v0, start_scan

      implicit none

!  array sizes: experimental (mnbed) and theory (mnbtd) beams before
!  averaging, theory energies as read (mnet), points of the energy
!  working grid (mngp), geometries (mns) and gaps (mngap)

      integer mnbed, mnbtd, mnet, mngp, mns, mngap
      include 'PARAM'

!  title: description of the experiment, then of the theory
!  nbtd: number of theory beams (re-read from the theory input)
!  nss: number of geometries; neti: theory energies as read
!  es, ats, atsav: theory energies and intensities of all geometries as
!          read and after averaging the beam groups
!  values: code of each geometry; pq, sym: momenta and symmetry codes of
!          the theory beams
!  ngp: never set; passed to grid as its declared array size (known
!          issue, see tests/README.md)

      type(control_parameters) ctl
      type(experiment_beams) ex
      type(theory_beams) th
      type(beam_results) res
      type(best_fit) best
      character(len=80) title
      integer nbtd, nss, neti, is, ngp
      real(sp) es(mngp), ats(mns,mnbtd,mnet), atsav(mns,mnbtd,mnet)
      real(sp) values(mns), pq(2,mnbtd)
      integer sym(mnbtd)

 30   format('*STRUCTURE NO.',1I5,' CHARACTERIZED BY VALUE',1F7.4)

      open(7,file='ROUT')
      open(9,file='ROUTSHORT')
      open(8,file='WEXPEL',status='OLD')
      rewind(8)

      call read_control(ctl, mnbed, mnbtd, mns, nss, title)
      call read_experiment(ex, ctl, mnbed, mngp, mngap)
      call prepare_experiment(ex, ctl, ngp)
      close(8)

      call read_theory(ctl, th, nss, nbtd, neti, es, ats, atsav, values,
     &                 pq, sym, title)

      call start_scan(res, best, mnbed, mnbtd, mngp, mngap)
      do is = 1, nss
        write(6,30) is, values(is)
        call prepare_theory(th, ctl, atsav, nss, nbtd, neti, is, es)
        call scan_v0(ctl, ex, th, nbtd, is, values(is), res, best)
      end do

      call report_best_fit(ctl, ex, best, nbtd, title, values)

      stop

      end program rfactor
