C  Calls the legacy routines of lib/lib.tleed.f that set up the angular
C  coefficients (CAAA, CELMG, LXGENT, CPPP) and the sublattices (SLIND) on
C  the inputs read from stdin and prints their results. Same protocol as
C  tier1/refcalc_routines/driver.f90 (see there).
      PROGRAM LEGDRV
      IMPLICIT NONE
      CHARACTER*16 NAME

      READ(5,'(A)') NAME
      IF (NAME.EQ.'ANGULAR') THEN
        CALL ANGULR
      ELSE IF (NAME.EQ.'SLIND') THEN
        CALL SLDRV
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
      INTEGER IDEG, NL, I, NL1, NL2
      REAL BR1(2), BR2(2), AR1(2), AR2(2), RAR1(2), RAR2(2), T(12)
      COMMON /SL/BR1, BR2, AR1, AR2, RAR1, RAR2, NL1, NL2
      REAL, ALLOCATABLE :: V(:,:), VLR(:)
      COMPLEX, ALLOCATABLE :: VL(:,:)
      INTEGER, ALLOCATABLE :: JJS(:,:)

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
