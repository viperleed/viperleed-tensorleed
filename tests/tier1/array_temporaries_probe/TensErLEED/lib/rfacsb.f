C  Probe fixture: the strided section A(1:N:2) is passed on purpose,
C  so gfortran must create an array temporary. No IMPLICIT NONE or
C  INTENT on purpose: the fixture imitates the legacy code.
      SUBROUTINE SHIFT(A, N)
      INTEGER N
      REAL A(N)
      CALL SCALE(A(1:N:2), N/2)
      END
      SUBROUTINE SCALE(B, M)
      INTEGER M
      REAL B(M)
      B = 2.0*B
      END
