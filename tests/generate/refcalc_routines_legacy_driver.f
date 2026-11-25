C  Calls the legacy routines of lib/lib.tleed.f that set up the angular
C  coefficients (CAAA, CELMG, LXGENT, CPPP) and the sublattices (SLIND),
C  and those of the energy loop up to the layer diffraction matrices of
C  a Bravais layer (TSCATF, BEAMS, SPHRM_MOD, BEAMS_YLM, FMAT, RSMF_SIMPLE),
C  and the composite layer (SRTLAY and GHMAT; TAUMAT; RTINV_SIMPLE), on the
C  inputs read from stdin and prints their results. Same protocol as
C  tier1/refcalc_routines/driver.f90 (see there).
      MODULE LEGHLP
      IMPLICIT NONE
      CONTAINS

      SUBROUTINE RDBEAMS(KNBS, KNB, KNT, TST, SPQ, SPQF, NPU)
C  The beam list as for BEAMS: set count and sizes, the incident-beam
C  index sets (only if NPU is present), the criterion TST, and the beam
C  momenta and index pairs
      INTEGER KNBS, KNT
      INTEGER, ALLOCATABLE :: KNB(:)
      REAL TST, RDR1
      REAL, ALLOCATABLE :: SPQ(:,:), SPQF(:,:)
      INTEGER, ALLOCATABLE, OPTIONAL :: NPU(:)
      INTEGER NPUN

      READ(5,*) KNBS
      ALLOCATE(KNB(KNBS))
      READ(5,*) KNB
      KNT = SUM(KNB)
      IF (PRESENT(NPU)) THEN
        READ(5,*) NPUN
        ALLOCATE(NPU(NPUN))
        READ(5,*) NPU
      ENDIF
      TST = RDR1()
      ALLOCATE(SPQ(2,KNT), SPQF(2,KNT))
      CALL RDR(SPQ, 2*KNT)
      CALL RDR(SPQF, 2*KNT)
      END SUBROUTINE

      SUBROUTINE SETLAT(IDEG, NL, LMAX, KLM, V, JJS, FLMS, SCC, SA)
C  The lattice of an SLIND case into /SL/, LMAX, and the sublattices
      INTEGER IDEG, NL, LMAX, KLM, NL1, NL2
      REAL, ALLOCATABLE :: V(:,:)
      INTEGER, ALLOCATABLE :: JJS(:,:)
      COMPLEX, ALLOCATABLE :: FLMS(:,:), SCC(:,:), SA(:), VL(:,:)

      CALL RDLAT(IDEG, NL1, NL2)
      NL = NL1*NL2
      READ(5,*) LMAX
      KLM = (2*LMAX+1)*(2*LMAX+2)/2
      ALLOCATE(V(NL,2), VL(NL,2), JJS(NL,IDEG), FLMS(NL,KLM))
      ALLOCATE(SCC(IDEG,4*LMAX+1), SA(2*LMAX*IDEG))
      V = 0.
      VL = (0., 0.)
      JJS = 0
      CALL SLIND(V, VL, JJS, NL, IDEG, 2.0E-4)
      END SUBROUTINE

      SUBROUTINE RUNFMAT(E, VPI, V, JJS, NL, IDEG, LMAX, KLM, FLMS,
     +                   SCC, SA)
C  FMAT with the cutoff distance of ref-calc.f
      INTEGER NL, IDEG, LMAX, KLM, JJS(NL,IDEG)
      REAL E, VPI, V(NL,2), DCUT
      COMPLEX FLMS(NL,KLM), SCC(IDEG,4*LMAX+1), SA(2*LMAX*IDEG)

      DCUT = - 5.0 * SQRT(2.0*E)/(AMIN1(VPI, - 0.05))
      FLMS = (0., 0.)
      CALL FMAT(FLMS, V, JJS, NL, NL, DCUT, IDEG, LMAX, KLM, SCC, SA)
      END SUBROUTINE
      END MODULE

      PROGRAM LEGDRV
      USE LEGHLP
      IMPLICIT NONE
      CHARACTER*16 NAME

      READ(5,'(A)') NAME
      IF (NAME.EQ.'ANGULAR') THEN
        CALL ANGULR
      ELSE IF (NAME.EQ.'SLIND') THEN
        CALL SLDRV
      ELSE IF (NAME.EQ.'TSCATF') THEN
        CALL TSDRV
      ELSE IF (NAME.EQ.'BEAMS') THEN
        CALL BMDRV
      ELSE IF (NAME.EQ.'SPHRM') THEN
        CALL SPDRV
      ELSE IF (NAME.EQ.'BYLM') THEN
        CALL BYDRV
      ELSE IF (NAME.EQ.'FMAT') THEN
        CALL FMDRV
      ELSE IF (NAME.EQ.'RSMF') THEN
        CALL RSDRV
      ELSE IF (NAME.EQ.'GHMAT' .OR. NAME.EQ.'TAUMAT' .OR.
     +         NAME.EQ.'RTINV') THEN
        CALL CPDRV(NAME)
      ELSE
        STOP 'unknown routine'
      ENDIF
      END

      SUBROUTINE ANGULR
C  The calls of ref-calc.f: CAAA, CELMG, LXGENT, CPPP, in this order, with
C  the sizes of ref-calc.f for MLMAX = 18
      IMPLICIT NONE
      INTEGER LMAX, LMMAX, NCAA, N, NN, NLM, NFAC, NN1, NN2, NN3
      INTEGER NCA(18), NLMS(18)
      REAL, ALLOCATABLE :: CAA(:), CLM(:), YLMC(:), FAC2(:), FAC1(:)
      REAL, ALLOCATABLE :: PPP(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: FAC(:)
      INTEGER, ALLOCATABLE :: LX(:), LXI(:), LT(:), LXM(:)

      DATA NCA/1,70,264,759,1820,3836,7344,13053,21868,34914,53560,
     +     79443,114492,160952,221408,298809,396492,518206/
      DATA NLMS/0,76,284,809,1925,4032,7680,13593,22693,36124,55276,
     +     81809,117677,165152,226848,305745,405213,529036/

      READ(5,*) LMAX
      LMMAX = (LMAX+1)**2
      NFAC = 4*18 + 3

      NCAA = NCA(LMAX)
      ALLOCATE(CAA(NCAA))
      CAA = 0.
      CALL CAAA(CAA,NCAA,LMMAX)
      CALL WRR('CAA', CAA, NCAA)

      N = 2*LMAX + 1
      NN = N*N
      NLM = NLMS(LMAX)
      ALLOCATE(CLM(NLM), YLMC(NN), FAC2(NN), FAC1(N), FAC(NFAC))
      CLM = 0.
      FAC = 0.
      CALL CELMG(CLM,NLM,YLMC,FAC2,NN,FAC1,N,LMAX,NFAC,FAC)
      CALL WRR('CLM', CLM, NLM)
      CALL WRD('FAC_CELMG', FAC, NFAC)

      ALLOCATE(LX(LMMAX), LXI(LMMAX), LT(LMMAX), LXM(LMMAX))
      LX = 0
      LXI = 0
      LT = 0
      LXM = 0
      CALL LXGENT(LX,LXI,LT,LXM,LMAX,LMMAX)
      CALL WRI('LX', LX, LMMAX)
      CALL WRI('LXI', LXI, LMMAX)
      CALL WRI('LT', LT, LMMAX)
      CALL WRI('LXM', LXM, LMMAX)

      NN3 = LMAX + 1
      NN2 = LMAX + 1
      NN1 = NN2 + NN3 - 1
      ALLOCATE(PPP(NN1,NN2,NN3))
      PPP = 0.
      CALL CPPP(PPP,NN1,NN2,NN3,NFAC,FAC)
      CALL WRR('PPP', PPP, NN1*NN2*NN3)
      CALL WRD('FAC_CPPP', FAC, NFAC)
      END

      SUBROUTINE SLDRV
C  SLIND on the lattice read from stdin (substrate and superlattice
C  vectors, reciprocal superlattice vectors), EPSD as in ref-calc.f
      IMPLICIT NONE
      INTEGER IDEG, NL, NL1, NL2
      REAL, ALLOCATABLE :: V(:,:), VLR(:)
      COMPLEX, ALLOCATABLE :: VL(:,:)
      INTEGER, ALLOCATABLE :: JJS(:,:)

      CALL RDLAT(IDEG, NL1, NL2)
      NL = NL1*NL2
      ALLOCATE(V(NL,2), VL(NL,2), VLR(4*NL), JJS(NL,IDEG))
      V = 0.
      VL = (0., 0.)
      JJS = 0
      CALL SLIND(V, VL, JJS, NL, IDEG, 2.0E-4)
      VLR = TRANSFER(VL, VLR)
      CALL WRR('V', V, 2*NL)
      CALL WRR('VL', VLR, 4*NL)
      CALL WRI('JJS', JJS, NL*IDEG)
      END

      REAL FUNCTION RDR1()
      INTEGER BITS
      READ(5,'(Z8)') BITS
      RDR1 = TRANSFER(BITS, RDR1)
      END

      SUBROUTINE RDR(V, N)
      INTEGER N, I
      REAL V(N), RDR1
      DO I = 1,N
        V(I) = RDR1()
      ENDDO
      END

      SUBROUTINE WRR(NAME, V, N)
      CHARACTER*(*) NAME
      INTEGER N, I
      REAL V(N)
      WRITE(6,'(A,1X,A,1X,I0)') NAME, 'r4', N
      WRITE(6,'(8Z8.8)') (TRANSFER(V(I), 0), I = 1,N)
      END

      SUBROUTINE WRD(NAME, V, N)
      CHARACTER*(*) NAME
      INTEGER N, I
      DOUBLE PRECISION V(N)
      INTEGER*8 ZERO8
      ZERO8 = 0
      WRITE(6,'(A,1X,A,1X,I0)') NAME, 'r8', N
      WRITE(6,'(4Z16.16)') (TRANSFER(V(I), ZERO8), I = 1,N)
      END

      SUBROUTINE WRI(NAME, V, N)
      CHARACTER*(*) NAME
      INTEGER N, I
      INTEGER V(N)
      WRITE(6,'(A,1X,A,1X,I0)') NAME, 'i4', N
      WRITE(6,'(8(I0,1X))') (V(I), I = 1,N)
      END

      SUBROUTINE RDC(C, N)
C  N complex values from hex bit patterns, real part first
      INTEGER N
      COMPLEX C(N)
      REAL T(2*N)
      CALL RDR(T, 2*N)
      C = CMPLX(T(1:2*N:2), T(2:2*N:2))
      END

      SUBROUTINE WRC(NAME, C, N)
      CHARACTER*(*) NAME
      INTEGER N
      COMPLEX C(N)
      REAL T(2*N)
      T(1:2*N:2) = REAL(C)
      T(2:2*N:2) = AIMAG(C)
      CALL WRR(NAME, T, 2*N)
      END

      SUBROUTINE NUMBERED(NAME, K, LABEL)
C  LABEL = NAME with the two-digit number K appended
      CHARACTER*(*) NAME
      CHARACTER*16 LABEL
      INTEGER K
      WRITE(LABEL,'(A,A,I2.2)') NAME, '_', K
      END

      SUBROUTINE GAUNT(LMAX, PPP, NN1, NN2, NN3)
C  PPP as ref-calc.f computes it: CELMG fills FAC, CPPP uses it
      IMPLICIT NONE
      INTEGER LMAX, NN1, NN2, NN3, NFAC, N, NN, NLM, NLMS(18)
      REAL PPP(NN1,NN2,NN3)
      REAL, ALLOCATABLE :: CLM(:), YLMC(:), FAC2(:), FAC1(:)
      DOUBLE PRECISION, ALLOCATABLE :: FAC(:)
      DATA NLMS/0,76,284,809,1925,4032,7680,13593,22693,36124,55276,
     +     81809,117677,165152,226848,305745,405213,529036/
      NFAC = 4*18 + 3
      N = 2*LMAX + 1
      NN = N*N
      NLM = NLMS(LMAX)
      ALLOCATE(CLM(NLM), YLMC(NN), FAC2(NN), FAC1(N), FAC(NFAC))
      CLM = 0.
      FAC = 0.
      CALL CELMG(CLM,NLM,YLMC,FAC2,NN,FAC1,N,LMAX,NFAC,FAC)
      PPP = 0.
      CALL CPPP(PPP,NN1,NN2,NN3,NFAC,FAC)
      END

      SUBROUTINE XTABLES(LMAX, CLM, NLM, LX, LXI, LT, LXM, LMMAX)
C  CLM and the (L,M) permutations LX, LXI, LT, LXM as ref-calc.f
C  computes them
      IMPLICIT NONE
      INTEGER LMAX, NLM, LMMAX, LX(LMMAX), LXI(LMMAX), N, NN, NFAC
      INTEGER LT(LMMAX), LXM(LMMAX)
      REAL CLM(NLM)
      REAL, ALLOCATABLE :: YLMC(:), FAC2(:), FAC1(:)
      DOUBLE PRECISION, ALLOCATABLE :: FAC(:)
      NFAC = 4*18 + 3
      N = 2*LMAX + 1
      NN = N*N
      ALLOCATE(YLMC(NN), FAC2(NN), FAC1(N), FAC(NFAC))
      CLM = 0.
      FAC = 0.
      CALL CELMG(CLM,NLM,YLMC,FAC2,NN,FAC1,N,LMAX,NFAC,FAC)
      LX = 0
      LXI = 0
      LT = 0
      LXM = 0
      CALL LXGENT(LX,LXI,LT,LXM,LMAX,LMMAX)
      END

      SUBROUTINE TSDRV
C  The atomic t-matrix of one site as ref-calc.f builds it: TSCATF for
C  each element with CONC > 0, summed with the weights CONC (the ATA
C  average); T0 = T = 100, DR0 = 0, VSITE = 0
      IMPLICIT NONE
      REAL RDR1
      INTEGER LMAX, L1, NEL, NPSI, NN1, NN2, NN3, NCALL, K, IEL, IL
      REAL EB, DR0, VSITE, T0, T, DRPER, DRPAR, CONC, VIB
      REAL, ALLOCATABLE :: ES(:), PHSS(:,:,:), PHS(:), PPP(:,:,:)
      COMPLEX, ALLOCATABLE :: AF(:), CAF(:), DEL(:), CTAB(:), SUM(:)
      COMPLEX, ALLOCATABLE :: TMAT(:)
      COMPLEX*16, ALLOCATABLE :: BJ(:)
      CHARACTER*16 LABEL

      READ(5,*) LMAX, NEL, NPSI
      L1 = LMAX + 1
      NN3 = L1
      NN2 = L1
      NN1 = NN2 + NN3 - 1
      ALLOCATE(ES(NPSI), PHSS(NPSI,NEL,L1), PHS(L1), PPP(NN1,NN2,NN3))
      ALLOCATE(AF(L1), CAF(L1), DEL(L1), CTAB(NN3), SUM(NN2), TMAT(L1))
      ALLOCATE(BJ(NN1))
      CALL RDR(ES, NPSI)
      CALL RDR(PHSS, NPSI*NEL*L1)
      EB = RDR1()
      CALL GAUNT(LMAX, PPP, NN1, NN2, NN3)
      DR0 = 0.
      VSITE = 0.
      T0 = 100.
      T = 100.
      TMAT = (0., 0.)
      READ(5,*) NCALL
      DO K = 1, NCALL
        READ(5,*) IEL
        CONC = RDR1()
        VIB = RDR1()
        DRPAR = RDR1()
        IF (CONC.GT.0.) THEN
          DRPER = VIB
          CALL TSCATF(IEL,L1,ES,PHSS,NPSI,EB,VSITE,PPP,NN1,NN2,NN3,
     +                DR0,DRPER,DRPAR,T0,T,AF,CAF,NEL,L1,PHS,DEL,
     +                CTAB,SUM,BJ)
          CALL NUMBERED('PHS', K, LABEL)
          CALL WRR(LABEL, PHS, L1)
          CALL NUMBERED('AF', K, LABEL)
          CALL WRC(LABEL, AF, L1)
          CALL NUMBERED('DEL', K, LABEL)
          CALL WRC(LABEL, DEL, L1)
          CALL NUMBERED('CAF', K, LABEL)
          CALL WRC(LABEL, CAF, L1)
          DO IL = 1, L1
            TMAT(IL) = TMAT(IL) + CONC*CAF(IL)
          ENDDO
        ENDIF
      ENDDO
      CALL WRC('TMAT', TMAT, L1)
      END

      SUBROUTINE BMDRV
C  BEAMS at several energies on one beam list
      USE LEGHLP
      IMPLICIT NONE
      REAL RDR1
      INTEGER KNBS, KNT, NPUN, NE, IE, NT, NP, MPU
      INTEGER, ALLOCATABLE :: KNB(:), NPU(:), NB(:), NPUC(:)
      REAL, ALLOCATABLE :: SPQ(:,:), SPQF(:,:), PQ(:,:), PQF(:,:)
      REAL TST, E, AK2, AK3
      CHARACTER*16 LABEL

      CALL RDBEAMS(KNBS, KNB, KNT, TST, SPQ, SPQF, NPU)
      NPUN = SIZE(NPU)
      ALLOCATE(NB(KNBS), NPUC(NPUN), PQ(2,KNT), PQF(2,KNT))
      READ(5,*) NE
C  NPUC is not reset between energies, as in ref-calc. -7 is a sentinel
C  that no beam number takes, so an entry the routine leaves unset shows
C  in the output.
      NPUC = -7
      DO IE = 1, NE
        E = RDR1()
        AK2 = RDR1()
        AK3 = RDR1()
        NB = 0
        PQ = 0.
        PQF = 0.
        CALL BEAMS(KNBS,KNB,SPQ,SPQF,KNT,NPU,NPUN,AK2,AK3,E,TST,NB,PQ,
     +             PQF,NPUC,MPU,NT,NP)
        CALL NUMBERED('NB', IE, LABEL)
        CALL WRI(LABEL, NB, KNBS)
        CALL NUMBERED('NT_NP_MPU', IE, LABEL)
        CALL WRI(LABEL, (/NT, NP, MPU/), 3)
        CALL NUMBERED('NPUC', IE, LABEL)
        CALL WRI(LABEL, NPUC, NPUN)
        CALL NUMBERED('PQ', IE, LABEL)
        CALL WRR(LABEL, PQ, 2*NT)
        CALL NUMBERED('PQF', IE, LABEL)
        CALL WRR(LABEL, PQF, 2*NT)
      ENDDO
      END

      SUBROUTINE SPDRV
C  SPHRM_MOD on a grid of CT, ST, CF for every LMAX from 2 to 18
      IMPLICIT NONE
      INTEGER NPTS, LMAX, LMMAX, I
      COMPLEX, ALLOCATABLE :: CT(:), ST(:), CF(:), YLM(:), ALL(:)
      CHARACTER*16 LABEL

      READ(5,*) NPTS
      ALLOCATE(CT(NPTS), ST(NPTS), CF(NPTS))
      CALL RDC(CT, NPTS)
      CALL RDC(ST, NPTS)
      CALL RDC(CF, NPTS)
      DO LMAX = 2, 18
        LMMAX = (LMAX+1)**2
        ALLOCATE(YLM(LMMAX), ALL(LMMAX*NPTS))
        DO I = 1, NPTS
          YLM = (0., 0.)
          CALL SPHRM_MOD(LMAX, YLM, LMMAX, CT(I), ST(I), CF(I))
          ALL((I-1)*LMMAX+1:I*LMMAX) = YLM
        ENDDO
        CALL NUMBERED('YLM_L', LMAX, LABEL)
        CALL WRC(LABEL, ALL, LMMAX*NPTS)
        DEALLOCATE(YLM, ALL)
      ENDDO
      END

      SUBROUTINE RDLAT(IDEG, NL1, NL2)
C  The lattice of an SLIND case into /SL/
      IMPLICIT NONE
      INTEGER IDEG, NL1, NL2, I
      REAL BR1(2), BR2(2), AR1(2), AR2(2), RAR1(2), RAR2(2), T(12)
      COMMON /SL/BR1, BR2, AR1, AR2, RAR1, RAR2, NL1X, NL2X
      INTEGER NL1X, NL2X

      READ(5,*) IDEG, NL1, NL2
      CALL RDR(T, 12)
      DO I = 1,2
        BR1(I) = T(I)
        BR2(I) = T(2+I)
        AR1(I) = T(4+I)
        AR2(I) = T(6+I)
        RAR1(I) = T(8+I)
        RAR2(I) = T(10+I)
      ENDDO
      NL1X = NL1
      NL2X = NL2
      END

      SUBROUTINE FMDRV
C  FMAT on the sublattices of SLIND at several (E, VPI, k-parallel)
      USE LEGHLP
      IMPLICIT NONE
      REAL RDR1
      INTEGER IDEG, NL, LMAX, KLM, NE, IE
      REAL E, AK2, AK3, VPI, TV, EMACH
      COMMON E, AK2, AK3, VPI, TV, EMACH
      REAL, ALLOCATABLE :: V(:,:)
      COMPLEX, ALLOCATABLE :: FLMS(:,:), SCC(:,:), SA(:)
      INTEGER, ALLOCATABLE :: JJS(:,:)
      CHARACTER*16 LABEL

      CALL SETLAT(IDEG, NL, LMAX, KLM, V, JJS, FLMS, SCC, SA)
      TV = 0.
      EMACH = 1.0E-6
      READ(5,*) NE
      DO IE = 1, NE
        E = RDR1()
        VPI = RDR1()
        AK2 = RDR1()
        AK3 = RDR1()
        CALL RUNFMAT(E, VPI, V, JJS, NL, IDEG, LMAX, KLM, FLMS, SCC, SA)
        CALL NUMBERED('FLMS', IE, LABEL)
        CALL WRC(LABEL, FLMS, NL*KLM)
      ENDDO
      END

      SUBROUTINE RSDRV
C  BEAMS, FMAT and RSMF_SIMPLE chained as in ref-calc.f: for the layer
C  of the overlayer type (mode O: one call for all beams, one sublattice)
C  or of the substrate type (mode B: one call per beam set with the
C  offset of the set in the beam list, all sublattices)
      USE LEGHLP
      IMPLICIT NONE
      REAL RDR1
      INTEGER IDEG, NL, LMAX, L1, LMMAX, KLM, LEV, LOD, NLM
      INTEGER KNBS, KNT, NPUN, NT, NP, MPU, IBS, NA, NAA, NLL
      INTEGER NLMS(18)
      REAL E, AK2, AK3, VPI, TV, EMACH, TST
      COMMON E, AK2, AK3, VPI, TV, EMACH
      INTEGER LMAXC, LITER
      REAL EPS
      COMMON /MS/ LMAXC, EPS, LITER
      CHARACTER*1 MODE
      INTEGER, ALLOCATABLE :: KNB(:), NPU(:), NB(:), NPUC(:)
      INTEGER, ALLOCATABLE :: JJS(:,:), LX(:), LXI(:), IPLE(:), IPLO(:)
      INTEGER, ALLOCATABLE :: LT(:), LXM(:)
      REAL, ALLOCATABLE :: SPQ(:,:), SPQF(:,:), PQ(:,:), PQF(:,:)
      REAL, ALLOCATABLE :: V(:,:), CLM(:)
      COMPLEX, ALLOCATABLE :: FLMS(:,:), SCC(:,:), SA(:)
      COMPLEX, ALLOCATABLE :: TMAT(:), RA(:,:), TA(:,:), AMULT(:)
      COMPLEX, ALLOCATABLE :: CYLM(:,:), FLM(:), XEV(:,:), XOD(:,:)
      COMPLEX, ALLOCATABLE :: YLM(:), YLME(:), YLMO(:)
      COMPLEX, ALLOCATABLE :: XEVST(:,:), XODST(:,:)
      DATA NLMS/0,76,284,809,1925,4032,7680,13593,22693,36124,55276,
     +     81809,117677,165152,226848,305745,405213,529036/

      CALL SETLAT(IDEG, NL, LMAX, KLM, V, JJS, FLMS, SCC, SA)
      L1 = LMAX + 1
      LMMAX = L1*L1
      LEV = L1*(L1+1)/2
      LOD = LMMAX - LEV
      NLM = NLMS(LMAX)
      LMAXC = LMAX
      EPS = 0.001
      LITER = 10
      E = RDR1()
      VPI = RDR1()
      AK2 = RDR1()
      AK3 = RDR1()
      TV = RDR1()
      EMACH = 1.0E-6
      READ(5,*) MODE
      ALLOCATE(TMAT(L1))
      CALL RDC(TMAT, L1)
      CALL RDBEAMS(KNBS, KNB, KNT, TST, SPQ, SPQF)
      NPUN = 1
      ALLOCATE(NB(KNBS), NPU(NPUN), NPUC(NPUN), PQ(2,KNT), PQF(2,KNT))
      NPU = 1

      ALLOCATE(CLM(NLM), LX(LMMAX), LXI(LMMAX), LT(LMMAX), LXM(LMMAX))
      CALL XTABLES(LMAX, CLM, NLM, LX, LXI, LT, LXM, LMMAX)

      PQ = 0.
      PQF = 0.
      NB = 0
C  Sentinel that no beam number takes: an unset entry shows in the output
      NPUC = -7
      CALL BEAMS(KNBS,KNB,SPQ,SPQF,KNT,NPU,NPUN,AK2,AK3,E,TST,NB,PQ,
     +           PQF,NPUC,MPU,NT,NP)
      CALL RUNFMAT(E, VPI, V, JJS, NL, IDEG, LMAX, KLM, FLMS, SCC, SA)

      ALLOCATE(RA(NT,NT), TA(NT,NT), AMULT(NT), CYLM(NT,LMMAX))
      ALLOCATE(FLM(KLM), XEV(LEV,LEV), XOD(LOD,LOD), YLM(LMMAX))
      ALLOCATE(YLME(LEV), YLMO(LOD), IPLE(LEV), IPLO(LOD))
      ALLOCATE(XEVST(LEV,LEV), XODST(LOD,LOD))
      RA = (0., 0.)
      TA = (0., 0.)
      AMULT = (0., 0.)
      CYLM = (0., 0.)
      XEVST = (0., 0.)
      XODST = (0., 0.)
      IF (MODE.EQ.'O') THEN
        NLL = 1
        NAA = NT
        CALL RSMF_SIMPLE(RA,TA,NAA,AMULT,CYLM,PQ,NT,FLMS,FLM,V,NL,0,
     +                   NLL,TMAT,TMAT,L1,LX,LXI,LMMAX,KLM,XEV,XOD,
     +                   LEV,LOD,YLM,YLME,YLMO,IPLE,IPLO,CLM,NLM,
     +                   XEVST,XODST)
      ELSE
        NLL = NL
        NA = 0
        DO IBS = 1, KNBS
          NAA = NB(IBS)
          CALL RSMF_SIMPLE(RA,TA,NAA,AMULT,CYLM,PQ,NT,FLMS,FLM,V,NL,NA,
     +                     NLL,TMAT,TMAT,L1,LX,LXI,LMMAX,KLM,XEV,XOD,
     +                     LEV,LOD,YLM,YLME,YLMO,IPLE,IPLO,CLM,NLM,
     +                     XEVST,XODST)
          NA = NA + NB(IBS)
        ENDDO
      ENDIF
      CALL WRI('NB', NB, KNBS)
      CALL WRC('RA', RA, NT*NT)
      CALL WRC('TA', TA, NT*NT)
      CALL WRC('AMULT', AMULT, NT)
      CALL WRC('CYLM', CYLM, NT*LMMAX)
      CALL WRC('XEVST', XEVST, LEV*LEV)
      CALL WRC('XODST', XODST, LOD*LOD)
      END

      SUBROUTINE BYDRV
C  BEAMS_YLM for each beam set of a beam list (with its offset in the
C  list), and for a window of the list that starts inside the first set
C  and ends before the last beam. CYLM is reset to a sentinel before every
C  call and written in full after it, so that rows that the routine must
C  leave alone show.
      USE LEGHLP
      IMPLICIT NONE
      REAL RDR1
      INTEGER LMAX, LMMAX, KNBS, KNT, NPUN, NT, NP, MPU, IBS, NA, N
      REAL TST, E, AK2, AK3, VPI
      INTEGER, ALLOCATABLE :: KNB(:), NPU(:), NB(:), NPUC(:)
      REAL, ALLOCATABLE :: SPQ(:,:), SPQF(:,:), PQ(:,:), PQF(:,:)
      COMPLEX, ALLOCATABLE :: XA(:), CYLM(:,:), YLM(:)
      COMPLEX SENT
      CHARACTER*16 LABEL

      SENT = (7.0, -3.0)
      E = RDR1()
      VPI = RDR1()
      AK2 = RDR1()
      AK3 = RDR1()
      READ(5,*) LMAX
      LMMAX = (LMAX+1)**2
      CALL RDBEAMS(KNBS, KNB, KNT, TST, SPQ, SPQF)
      NPUN = 1
      ALLOCATE(NB(KNBS), NPU(NPUN), NPUC(NPUN), PQ(2,KNT), PQF(2,KNT))
      NPU = 1
      PQ = 0.
      PQF = 0.
      NB = 0
C  Sentinel that no beam number takes: an unset entry shows in the output
      NPUC = -7
      CALL BEAMS(KNBS,KNB,SPQ,SPQF,KNT,NPU,NPUN,AK2,AK3,E,TST,NB,PQ,
     +           PQF,NPUC,MPU,NT,NP)
      ALLOCATE(CYLM(NT,LMMAX), YLM(LMMAX))
      NA = 0
      DO IBS = 1, KNBS
        ALLOCATE(XA(NB(IBS)))
        CYLM = SENT
        CALL BEAMS_YLM(CYLM,XA,PQ,NT,NB(IBS),NA,E,AK2,AK3,VPI,LMAX,
     +                 LMMAX,YLM)
        CALL NUMBERED('XA', IBS, LABEL)
        CALL WRC(LABEL, XA, NB(IBS))
        CALL NUMBERED('CYLM', IBS, LABEL)
        CALL WRC(LABEL, CYLM, NT*LMMAX)
        DEALLOCATE(XA)
        NA = NA + NB(IBS)
      ENDDO
      CALL WRI('NB', NB, KNBS)
      N = NT - 3
      ALLOCATE(XA(N))
      CYLM = SENT
      CALL BEAMS_YLM(CYLM,XA,PQ,NT,N,1,E,AK2,AK3,VPI,LMAX,LMMAX,YLM)
      CALL WRC('XA_WINDOW', XA, N)
      CALL WRC('CYLM_WINDOW', CYLM, NT*LMMAX)
      END

      SUBROUTINE CPDRV(NAME)
C  SRTLAY and GHMAT, TAUMAT, or RTINV_SIMPLE on a composite layer, with
C  the arguments, array sizes and COMMON blocks of ref-calc.f. Mode O:
C  one call for all beams (offsets NA = NS = 0); mode B: one call per
C  beam set, with the offset of the set in the beam list as NA and the
C  sum of the earlier set sizes as NS.
      USE LEGHLP
      IMPLICIT NONE
      CHARACTER*16 NAME, LABEL
      CHARACTER*1 MODE
      REAL RDR1
      INTEGER IDEG, NL, LMAX, KLM, L1, LMMAX, LEV, LOD, LEV2, NLM, NCAA
      INTEGER MLMAX, MNSUB, NLAY, LAYT, NLTYPE, ILTYPE, NNSUB, LMAX1
      INTEGER NLAY2, LMG, LMN, LM2N, LMNI, L2M, LMS, LMT, NTAUV, NLL
      INTEGER KNBS, KNT, NPUN, NT, NP, MPU, NCALLS, K, NAA, MLMMAX
      INTEGER NCA(18), NLMS(18), INFO
      REAL E, AK2, AK3, VPI, TV, EMACH, TSTB, TSTS, TV1, DCUTS
      COMMON E, AK2, AK3, VPI, TV, EMACH
      INTEGER LMAXC, LITER
      REAL EPS
      COMMON /MS/ LMAXC, EPS, LITER
      INTEGER NA, NS, ID, LAY, LM, NTAU, NOPT, NEW
      COMMON /MPT/ NA, NS, ID, LAY, LM, NTAU, TSTS, TV1, DCUTS, NOPT,
     +             NEW
      INTEGER, ALLOCATABLE :: KNB(:), NPU(:), NB(:), NPUC(:), JJS(:,:)
      INTEGER, ALLOCATABLE :: LX(:), LXI(:), LT(:), LXM(:), LPS(:)
      INTEGER, ALLOCATABLE :: LPSS(:), NORD(:), MGH(:,:), NUGH(:)
      INTEGER, ALLOCATABLE :: NGEQ(:), NGOL(:), IPL(:)
      REAL, ALLOCATABLE :: SPQ(:,:), SPQF(:,:), PQ(:,:), PQF(:,:)
      REAL, ALLOCATABLE :: V(:,:), CLM(:), CAA(:), SUBPOS(:,:,:)
      REAL, ALLOCATABLE :: POS(:,:), POSS(:,:), DRL(:,:), SDRL(:,:)
      REAL, ALLOCATABLE :: TEST(:)
      COMPLEX, ALLOCATABLE :: FLMS(:,:), SCC(:,:), SA(:), TSF(:,:)
      COMPLEX, ALLOCATABLE :: TAU(:,:), GH(:,:), TH(:,:), FLM(:)
      COMPLEX, ALLOCATABLE :: XEV(:,:), XH(:), HGHD(:), AMULT(:)
      COMPLEX, ALLOCATABLE :: CYLM(:,:), RTOP(:,:), TTOP(:,:)
      COMPLEX, ALLOCATABLE :: RBOT(:,:), TBOT(:,:), RG(:,:,:), RGP(:)
      COMPLEX, ALLOCATABLE :: TS(:), TG(:,:), VT(:), TAUG(:), TAUGM(:)
      COMPLEX, ALLOCATABLE :: RG1(:), RG2(:), TSTORE(:,:,:)
      DATA NCA/1,70,264,759,1820,3836,7344,13053,21868,34914,53560,
     +     79443,114492,160952,221408,298809,396492,518206/
      DATA NLMS/0,76,284,809,1925,4032,7680,13593,22693,36124,55276,
     +     81809,117677,165152,226848,305745,405213,529036/

      CALL SETLAT(IDEG, NL, LMAX, KLM, V, JJS, FLMS, SCC, SA)
      READ(5,*) MLMAX, MNSUB, NLAY, LAYT, NLTYPE, ILTYPE
      READ(5,*) MODE
      E = RDR1()
      VPI = RDR1()
      AK2 = RDR1()
      AK3 = RDR1()
      TV1 = RDR1()
      TSTS = RDR1()
      TV = TV1
      EMACH = 1.0E-6
      LMAXC = LMAX
      EPS = 0.001
      LITER = 10
      L1 = LMAX + 1
      LMMAX = L1*L1
      LEV = L1*(L1+1)/2
      LOD = LMMAX - LEV
      LEV2 = 2*LEV
      NLM = NLMS(LMAX)
      NCAA = NCA(LMAX)
      NNSUB = MNSUB
      LMAX1 = MLMAX + 1
      MLMMAX = LMAX1*LMAX1
      L2M = 2*MLMAX + 1
      LMS = L2M*L2M
      NLAY2 = NLAY*(NLAY-1)/2
      LMG = NLAY2*LMMAX*2
      LMN = NLAY*LMMAX
      LM2N = 2*LMN
      LMNI = LMN
      IF (LAYT.EQ.1) THEN
        NLL = 1
      ELSE
        NLL = NL
      ENDIF
      DCUTS = - 5.0 * SQRT(2.0*E)/(AMIN1(VPI, - 0.05))

      ALLOCATE(SUBPOS(NLTYPE,MNSUB,3), POS(NLAY,3), POSS(NLAY,3))
      CALL RDR(SUBPOS, 3*NLTYPE*MNSUB)
      READ(5,*) NTAUV
      ALLOCATE(LPS(NLAY), TSF(MNSUB,LMAX1))
      READ(5,*) LPS
      CALL RDC(TSF, MNSUB*LMAX1)
      CALL RDBEAMS(KNBS, KNB, KNT, TSTB, SPQ, SPQF)
      NPUN = 1
      ALLOCATE(NB(KNBS), NPU(NPUN), NPUC(NPUN), PQ(2,KNT), PQF(2,KNT))
      NPU = 1
      PQ = 0.
      PQF = 0.
      NB = 0
C  Sentinel that no beam number takes: an unset entry shows in the output
      NPUC = -7
      CALL BEAMS(KNBS,KNB,SPQ,SPQF,KNT,NPU,NPUN,AK2,AK3,E,TSTB,NB,PQ,
     +           PQF,NPUC,MPU,NT,NP)
      CALL RUNFMAT(E, VPI, V, JJS, NL, IDEG, LMAX, KLM, FLMS, SCC, SA)

      ALLOCATE(CLM(NLM), LX(LMMAX), LXI(LMMAX), LT(LMMAX), LXM(LMMAX))
      CALL XTABLES(LMAX, CLM, NLM, LX, LXI, LT, LXM, LMMAX)
      ALLOCATE(CAA(NCAA))
      CAA = 0.
      CALL CAAA(CAA, NCAA, LMMAX)
      CALL COPYPOS(POS, SUBPOS, NLTYPE, NNSUB, NLAY, ILTYPE)

      ID = 1
      NOPT = 1
      NEW = 1
      LAY = LAYT
      LM = L1
      NTAU = NTAUV
      LMT = NTAU*LMMAX
      ALLOCATE(MGH(NLAY,NLAY), DRL(NLAY2,3), SDRL(NLAY2,3))
      ALLOCATE(NUGH(NLAY2), NGEQ(NLAY2), NGOL(NLAY2), TEST(NLAY2))
      ALLOCATE(LPSS(NLAY), NORD(NLAY), GH(LMG,LMMAX), HGHD(L2M))
      ALLOCATE(TG(2,2*MNSUB*MLMMAX), VT(2*MNSUB*MLMMAX))
      ALLOCATE(TAU(LMT,LEV), XEV(LEV,LEV2), XH(LEV), FLM(KLM))
      MGH = 0
      DRL = 0.
      SDRL = 0.
      NUGH = 0
      NGEQ = 0
      NGOL = 0
      TEST = 0.
      LPSS = 0
      NORD = 0
      GH = (0., 0.)
      HGHD = (0., 0.)
      TG = (0., 0.)
      VT = (0., 0.)
      TAU = (0., 0.)
      XEV = (0., 0.)
      XH = (0., 0.)
      FLM = (0., 0.)
      IF (NAME.EQ.'RTINV') THEN
        ALLOCATE(TH(LMNI,LMNI), IPL(LMNI), AMULT(NT), CYLM(NT,LMMAX))
        ALLOCATE(RTOP(NT,NT), TTOP(NT,NT), RBOT(NT,NT), TBOT(NT,NT))
        ALLOCATE(RG(2,NLAY,NT), RGP(NT), TS(LMN), TAUG(LMT))
        ALLOCATE(TAUGM(LMT), RG1(NLAY), RG2(NLAY))
        ALLOCATE(TSTORE(2,LMN,NT))
        TH = (0., 0.)
        IPL = 0
        AMULT = (0., 0.)
        CYLM = (0., 0.)
        RTOP = (0., 0.)
        TTOP = (0., 0.)
        RBOT = (0., 0.)
        TBOT = (0., 0.)
        RG = (0., 0.)
        RGP = (0., 0.)
        TS = (0., 0.)
        TAUG = (0., 0.)
        TAUGM = (0., 0.)
        RG1 = (0., 0.)
        RG2 = (0., 0.)
        TSTORE = (0., 0.)
      ENDIF

      IF (MODE.EQ.'O') THEN
        NCALLS = 1
      ELSE
        NCALLS = KNBS
      ENDIF
      NA = 0
      NS = 0
      DO K = 1, NCALLS
        IF (MODE.EQ.'O') THEN
          NAA = NT
        ELSE
          NAA = NB(K)
        ENDIF
        IF (NAME.EQ.'GHMAT') THEN
          CALL SRTLAY(POS,POSS,LPS,LPSS,MGH,NLAY,DRL,SDRL,NLAY2,NUGH,
     +                NGEQ,NGOL,NEW,NORD)
          CALL GHMAT(GH,LMG,LMMAX,MGH,NLAY,NUGH,NGEQ,NGOL,NLAY2,TSTS,
     +               TEST,VT,L2M,VT,LM,TG,LMS,DRL,TV1,LXM,LEV,DCUTS,
     +               CAA,NCAA,LAY,HGHD,PQ(1,1+NA))
          CALL NUMBERED('GH', K, LABEL)
          CALL WRC(LABEL, GH, LMG*LMMAX)
        ELSE IF (NAME.EQ.'TAUMAT') THEN
          CALL TAUMAT(TAU,LMT,NTAU,XEV,LEV,LEV2,LOD,TSF,LMMAX,LMAX,
     +                FLMS,NL,KLM,LM,CLM,NLM,LXI,NT,PQ,NA,NLL,FLM,
     +                NNSUB,LMAX1,XH)
          CALL NUMBERED('TAU', K, LABEL)
          CALL WRC(LABEL, TAU, LMT*LEV)
          CALL NUMBERED('FLM', K, LABEL)
          CALL WRC(LABEL, FLM, KLM)
        ELSE
          CALL RTINV_SIMPLE(
     +         RTOP, TTOP, RBOT, TBOT, NAA, AMULT, CYLM, PQ, NT, FLMS,
     +         FLM, NL, LXI, LT, LXM, LMMAX, KLM, XEV, LEV, LEV2,
     +         TAU, LMT, TAUG, TAUGM, CLM, NLM, POS, POSS, MGH, NLAY,
     +         DRL, SDRL, NUGH, NGEQ, NGOL, NLAY2, TEST, GH, LMG, RG,
     +         RGP, TS, LMN, TG, LM2N, VT, CAA, NCAA, TH, LMNI,
     +         IPL, TSF, LPS, LPSS, NORD, NNSUB, LMAX1, XH, HGHD,
     +         L2M, RG1, RG2, TSTORE)
          IF (MODE.EQ.'B') THEN
            CALL NUMBERED('R_TOP', K, LABEL)
            CALL WRC(LABEL, RTOP, NT*NT)
            CALL NUMBERED('T_TOP', K, LABEL)
            CALL WRC(LABEL, TTOP, NT*NT)
            CALL NUMBERED('R_BOT', K, LABEL)
            CALL WRC(LABEL, RBOT, NT*NT)
            CALL NUMBERED('T_BOT', K, LABEL)
            CALL WRC(LABEL, TBOT, NT*NT)
            CALL NUMBERED('TSTORE', K, LABEL)
            CALL WRC(LABEL, TSTORE, 2*LMN*NT)
            CALL NUMBERED('AMULT', K, LABEL)
            CALL WRC(LABEL, AMULT, NAA)
            CALL NUMBERED('CYLM', K, LABEL)
            CALL WRC(LABEL, CYLM, NT*LMMAX)
          ENDIF
        ENDIF
        IF (MODE.EQ.'B') THEN
          NA = NA + NAA
          NS = NS + NAA
        ENDIF
      ENDDO

      CALL WRI('NB', NB, KNBS)
      IF (NAME.EQ.'GHMAT') THEN
        CALL WRI('MGH', MGH, NLAY*NLAY)
        CALL WRR('DRL', DRL, NLAY2*3)
        CALL WRR('POSS', POSS, NLAY*3)
        CALL WRI('LPSS', LPSS, NLAY)
        CALL WRI('NORD', NORD, NLAY)
        CALL WRI('NUGH', NUGH, NLAY2)
        CALL WRI('NGEQ', NGEQ, NLAY2)
        CALL WRR('TEST', TEST, NLAY2)
      ELSE IF (NAME.EQ.'RTINV') THEN
        CALL WRC('R_TOP', RTOP, NT*NT)
        CALL WRC('T_TOP', TTOP, NT*NT)
        CALL WRC('R_BOT', RBOT, NT*NT)
        CALL WRC('T_BOT', TBOT, NT*NT)
        CALL WRC('TSTORE', TSTORE, 2*LMN*NT)
        CALL WRC('TAU', TAU, LMT*LEV)
        CALL WRC('GH', GH, LMG*LMMAX)
        CALL WRR('POSS', POSS, NLAY*3)
        CALL WRI('NORD', NORD, NLAY)
        CALL WRC('AMULT', AMULT, NT)
        CALL WRC('CYLM', CYLM, NT*LMMAX)
      ENDIF
      END
