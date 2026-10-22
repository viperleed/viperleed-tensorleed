C  Calls one legacy routine of lib/rfacsb.f (DER, SMOOTH or YSMOOTH) on
C  the inputs read from stdin and prints its results. Same protocol as
C  tier1/rfacsb_routines/driver.f90: every real is read and written as
C  the 8 hex digits of its bit pattern, one value per line.
      PROGRAM LEGDRV
      IMPLICIT NONE
      CHARACTER*16 NAME
      INTEGER N, ISMOTH, NE(1)
      REAL H, VI, EINCR, RDR1
      REAL, ALLOCATABLE :: X(:), Y(:), AP(:), E(:)

      READ(5,'(A)') NAME
      IF (NAME.EQ.'DER') THEN
        READ(5,*) N
        H = RDR1()
        ALLOCATE(Y(N), X(N))
        CALL RDR(Y, N)
        CALL DER(Y, N, X, H)
        CALL WRR(X, N)
      ELSE IF (NAME.EQ.'SMOOTH') THEN
        READ(5,*) ISMOTH
        READ(5,*) N
        ALLOCATE(X(N), Y(N))
        CALL RDR(X, N)
        CALL RDR(Y, N)
        CALL SMOOTH(ISMOTH, Y, X, N, 0, 1)
        CALL WRR(Y, N)
      ELSE IF (NAME.EQ.'YSMOOTH') THEN
        READ(5,*) N
        VI = RDR1()
        EINCR = RDR1()
        ALLOCATE(X(N), AP(N), E(N), Y(N))
        CALL RDR(X, N)
        CALL RDR(AP, N)
        E = 0.
        Y = 0.
        NE(1) = N
        CALL YSMOOTH(X, AP, 1, 1, N, 1, 1, NE, E, Y, VI, EINCR, 0)
        CALL WRR(Y, N)
      ELSE
        STOP 'unknown routine'
      ENDIF
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

      SUBROUTINE WRR(V, N)
      INTEGER N, I
      REAL V(N)
      DO I = 1,N
        WRITE(6,'(Z8.8)') TRANSFER(V(I), 0)
      ENDDO
      END
