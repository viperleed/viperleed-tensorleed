      PROGRAM RFACTOR
      INCLUDE 'PARAM'
      REAL A(NP)
      A = 1.0
      CALL SHIFT(A, NP)
      PRINT *, A(1)
      END
