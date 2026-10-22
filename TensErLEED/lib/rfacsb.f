C  Circumvented allocation error in subroutine RINTAV by introduction 
C  of a second field ATAV and PQ1.        Oct. 2018 L. Hammer
C
C  In line 617: changed limit to 1.E-8.
C
C  Output for gnuplot adapted to particular computer settings. LH
C  INPUT formats expanded to allow for 9999 exp. and theo. beams  L.Hammer Mar. 2021
C  Note: Y function output in YPEND and YSMOOTH is still limited to 999 beams
C  Y_Smooth implementation (and nicer format of Y_Pendry) M. Schmid 2026
C
C-----------------------------------------------------------------------
C  Subroutine COLUMN writes the experimental spectra and the theoretical 110205
C  best-fit spectra in column format into the files exp.column and       110205
C  theo.column, respectively.                                            110205

      SUBROUTINE COLUMN(NEMINB,NEMAXB,NTMINB,NTMAXB,EE,AE,BET,BAT,NBED,
     *                  NBTD,NGP,NBE,IBP,GC,NGE1,NGE2,NGT1,NGT2,NGAP,
     *                  TSE,TST,V0,INTMAX,NORM)

      REAL AE(NBED,NGP), EE(NBED,NGP), AEI(NBED), EEI(NBED)
      INTEGER NEMINB(NBED), NEMAXB(NBED)

      REAL BAT(NBTD,NGP), BET(NBTD,NGP), BATI(NBED), BETI(NBED)
      INTEGER NTMINB(NBTD),NTMAXB(NBTD)
      REAL BATINT(NBED,NGP),BETINT(NBED,NGP),TSTINT(NBED)
      INTEGER NTMINI(NBED),NTMAXI(NBED)

      INTEGER IBP(NBED)
      INTEGER GC(NBED), NGE1(NGAP), NGE2(NGAP), NGT1(NGAP), NGT2(NGAP)
      REAL TSE(NBED), TST(NBTD), V0, INTMAX, AMAX(NBED), AMAXT(NBED)

      INTEGER IGC(NBED), NGE(NBED,NGAP,2), NGT(NBED,NGAP,2)
      INTEGER NETMIN,NETMAX

      CHARACTER*10 FMT

 10   FORMAT('(',I4,'F8.2)')
 11   FORMAT(A10)

      OPEN(10,FILE='exp.column')
      OPEN(11,FILE='theo.column')

C  determine spectra output format for exp.column and theo.column

      WRITE(10,10) 2*NBE
      REWIND(10)
      READ(10,11) FMT
      REWIND(10)

C  get lowest lower and highest upper energy indices of all beams 

      NEEMIN = NGP
      NEEMAX = 1
      NETMIN = NGP
      NETMAX = 1

      DO IBE = 1,NBE
         IF (IBP(IBE).EQ.0) THEN
            NTMINI(IBE) = NGP
            NTMAXI(IBE) = 1
            TSTINT(IBE) = 1.
            DO IGP = 1,NGP
            BETINT(IBE,IGP) = 0.
            BATINT(IBE,IGP) = 0.
            ENDDO
         ELSE
            NTMINI(IBE) = NTMINB(IBP(IBE))
            NTMAXI(IBE) = NTMAXB(IBP(IBE))
            TSTINT(IBE) = TST(IBP(IBE))
            DO IGP = 1,NGP
            BETINT(IBE,IGP) = BET(IBP(IBE),IGP)
            BATINT(IBE,IGP) = BAT(IBP(IBE),IGP)
            ENDDO
         ENDIF
      ENDDO

      DO IBE = 1,NBE
         AMAX(IBE) = 0.
         IF (NEMINB(IBE).LT.NEEMIN) THEN
            NEEMIN = NEMINB(IBE)
         ENDIF
         IF (NEMAXB(IBE).GT.NEEMAX) THEN
            NEEMAX = NEMAXB(IBE)
         ENDIF
         AMAXT(IBE) = 0.
         IF (NTMINI(IBE).LT.NETMIN) THEN
            NETMIN = NTMINI(IBE)
         ENDIF
         IF (NTMAXI(IBE).GT.NETMAX) THEN
            NETMAX = NTMAXI(IBE)
         ENDIF
      ENDDO

C  write gap indices in storage arrays

      DO IBE = 1,NBE
        DO I = 1,NGAP
          NGE(IBE,I,1) = 0
          NGE(IBE,I,2) = 0
          NGT(IBE,I,1) = 0
          NGT(IBE,I,2) = 0          
        ENDDO
      ENDDO

      IGAP = 0
      DO IBE = 1,NBE
        DO I = 1,GC(IBE)
          NGE(IBE,I,1) = NGE1(I+IGAP)
          NGE(IBE,I,2) = NGE2(I+IGAP)
          NGT(IBE,I,1) = NGT1(I+IGAP)
          NGT(IBE,I,2) = NGT2(I+IGAP)
        ENDDO
        IGAP = IGAP + GC(IBE)
      ENDDO

C  determine highest intensity AMAX(IBE) for each experimental beam (gaps are
C  neglected in the scheme) 

      DO IBE = 1,NBE
        IGC(IBE) = 0
      ENDDO

      DO IEE=NEEMIN,NEEMAX

        DO IBE = 1,NBE

          IF (IEE.GE.NEMINB(IBE).AND.IEE.LE.NEMAXB(IBE)) THEN

            IGDUM = IGC(IBE) + 1
            IF (GC(IBE).NE.0.AND.IEE.GE.NGE(IBE,IGDUM,1).AND.
     *         IEE.LE.NGE(IBE,IGDUM,2).AND.IGC(IBE).LT.GC(IBE)) THEN
              AEI(IBE) = 0.
              IF (IEE.EQ.NGE(IBE,IGDUM,2)) THEN
                IGC(IBE) = IGDUM
              ENDIF

            ELSEIF (NORM.EQ.1) THEN
              AEI(IBE) = AE(IBE,IEE)/TSE(IBE)

            ELSEIF (NORM.EQ.2) THEN
              AEI(IBE) = AE(IBE,IEE)
            ENDIF

            IF (AEI(IBE).GT.AMAX(IBE)) THEN
              AMAX(IBE) = AEI(IBE)
            ENDIF

          ENDIF
        ENDDO
      ENDDO

C  determine highest intensity AMAX(IBE) for theoretical beam (gaps are
C  neglected in the scheme) 

      DO IBE = 1,NBE
        IGC(IBE) = 0
      ENDDO

      DO IET=NETMIN,NETMAX

        DO IBE = 1,NBE

          IF (IET.GE.NTMINI(IBE).AND.IET.LE.NTMAXI(IBE)) THEN

            IGDUM =IGC(IBE) + 1
            IF (GC(IBE).NE.0.AND.IET.GE.NGT(IBE,IGDUM,1).AND.
     *         IET.LE.NGT(IBE,IGDUM,2).AND.IGC(IBE).LT.GC(IBE)) THEN
              BATI(IBE) = 0.
              IF (IET.EQ.NGT(IBE,IGDUM,2)) THEN
                IGC(IBE) = IGDUM
              ENDIF

            ELSEIF (NORM.EQ.1) THEN
              BATI(IBE) = BATINT(IBE,IET)/TSTINT(IBE)
              IF (BATI(IBE).GT.AMAX(IBE)) THEN
                AMAX(IBE) = BATI(IBE)
              ENDIF

            ELSEIF (NORM.EQ.2) THEN
              BATI(IBE) = BATINT(IBE,IET)
              IF (BATI(IBE).GT.AMAXT(IBE)) THEN
                AMAXT(IBE) = BATI(IBE)
              ENDIF
            ENDIF

          ENDIF
        ENDDO
      ENDDO

C  normalization for integral beam intensities (NORM = 1) or highest beam
C  intensity (NORM = 2)

      DO IBE = 1,NBE
        IGC(IBE) = 0
      ENDDO

      DO IEE=NEEMIN,NEEMAX

        DO IBE = 1,NBE

C  energies outside the common energy interval are ignored

          IGDUM =IGC(IBE) + 1
          IF (IEE.LT.NEMINB(IBE).OR.IEE.GT.NEMAXB(IBE)) THEN
            EEI(IBE) = 0.
            AEI(IBE) = 0.

          ELSEIF (GC(IBE).NE.0.AND.IEE.GE.NGE(IBE,IGDUM,1).AND.
     *           IEE.LE.NGE(IBE,IGDUM,2).AND.IGC(IBE).LT.GC(IBE)) THEN
            EEI(IBE) = 0.
            AEI(IBE) = 0.
            IF (IEE.EQ.NGE(IBE,IGDUM,2)) THEN
              IGC(IBE) = IGDUM
            ENDIF

C  do normalization

          ELSE
            EEI(IBE) = EE(IBE,IEE)
            IF (NORM.EQ.1) THEN
              AEI(IBE) = AE(IBE,IEE)*INTMAX/(TSE(IBE)*AMAX(IBE))
            ELSEIF (NORM.EQ.2) THEN
              AEI(IBE) = AE(IBE,IEE)*INTMAX/AMAX(IBE)
            ENDIF
          ENDIF

        ENDDO

C  write experimental beam data to exp.column

        WRITE(10,FMT) (EEI(IBE),AEI(IBE),IBE=1,NBE)
      ENDDO

      DO IBE = 1,NBE
        IGC(IBE) = 0
      ENDDO

      DO IET=NETMIN,NETMAX

        DO IBE = 1,NBE

C  energies outside the common energy interval are ignored

          IGDUM =IGC(IBE) + 1
          IF (IET.LT.NTMINI(IBE).OR.IET.GT.NTMAXI(IBE)) THEN
            BETI(IBE) = -V0
            BATI(IBE) = 0.

          ELSEIF (GC(IBE).NE.0.AND.IET.GE.NGT(IBE,IGDUM,1).AND.
     *           IET.LE.NGT(IBE,IGDUM,2).AND.IGC(IBE).LT.GC(IBE)) THEN
            BETI(IBE) = -V0
            BATI(IBE) = 0.
            IF (IET.EQ.NGT(IBE,IGDUM,2)) THEN
              IGC(IBE) = IGDUM
            ENDIF

C  do normalization

          ELSE
            BETI(IBE) = BETINT(IBE,IET)
            IF (NORM.EQ.1) THEN
              BATI(IBE) = BATINT(IBE,IET)*INTMAX/(TSTINT(IBE)*AMAX(IBE))
            ELSEIF (NORM.EQ.2) THEN
              BATI(IBE) = BATINT(IBE,IET)*INTMAX/AMAXT(IBE)
            ENDIF
          ENDIF

        ENDDO

C  write theoretical beam data in theo.column

        WRITE(11,FMT) (BETI(IBE)+V0,BATI(IBE),IBE=1,NBE)
      ENDDO

      END

C-----------------------------------------------------------------------
C  SUBROUTINE COMNEI FINDS ENERGY INTERVAL COMMON TO EXP. AND THEORY     040280
      SUBROUTINE COMNEI(EE,NBED,NEE,ET,NBTD,NGP,NET,IBE,IBT,V0,EINCR,
     1NE1,NE2,NT1,NT2)
C  FOR CDC ONLY  ***************************
C     LEVEL 2, ET,EE
      DIMENSION EE(NBED,NGP),NEE(NBED),ET(NBTD,NGP),NET(NBTD)
      NE=NEE(IBE)
      NT=NET(IBT)
      DE1=ET(IBT,1)+V0-EE(IBE,1)
      DE2=ET(IBT,NT)+V0-EE(IBE,NE)
      IF (DE1.LT.0.) GO TO 10
      NE1=INT((DE1/EINCR)+0.0001)+1
      NT1=1
      GO TO 20
10    NE1=1
      NT1=INT((-DE1/EINCR)+0.0001)+1
20    IF (DE2.LT.0.) GO TO 30
      NE2=NE
      NT2=NT-INT((DE2/EINCR)+0.0001)
      GO TO 40
30    NE2=NE-INT((-DE2/EINCR)+0.0001)
      NT2=NT
40    CONTINUE
      RETURN
      END
C-------------------------------------------------------------------------------
C  Subroutine EXPAV averages intensities from different experiments. It also
C  accounts for gaps in the input data, calculating gap info for the averaged
C  beams. The averaged intensities are output in order of increasing energy.

      SUBROUTINE EXPAV(AE,EE,NBED,NGP,NEE,BENAME,NBEA,NBE,IPR,A,NV,
     *                 GCOUNT,NGAP,EG1,EG2,NG1,NG2)

      REAL AE(NBED,NGP), EE(NBED,NGP)
      INTEGER NBED, NEE(NBED), NBEA(NBED), NBE, IPR
      REAL BENAME(5,NBED), A(NGP)
      INTEGER NV(NGP), GCOUNT(NBED), NGAP, NG1(NGP), NG2(NGP)
      REAL EG1(NGP), EG2(NGP)
      INTEGER IGAP, GIN, GC(NBED)
      REAL EMIN, EMAX, DE

 10   FORMAT(48H*EXP. ENERG. AND INTENS. AFTER AVERAGING IN BEAM,1I4,/,
     *       100(5(1F7.2,1E13.4,3X),/))

      NBE  = 1
      IGAP = 0                                                           201103
      NGAP = 0                                                           201103

      DO IB = 1,NBED
        IF (NBEA(IB).NE.0) THEN                                          040280

          EMIN = 1.E6                                                    040280
          EMAX = 0.                                                      040280
          DE   = 1.E6                                                    100281

C  loop to find all beams to be averaged with IB

          DO IB1 = IB,NBED                                               040280

            IF (NBEA(IB1).EQ.NBEA(IB)) THEN                              040280

C  find minimum and maximum energies of all beams to be averaged with beam IB
C  (including beam IB)

              DO IE = 1,NEE(IB1)                                         100281
                IF (EMIN.GT.EE(IB1,IE)) THEN                             100281
                  EMIN=EE(IB1,IE)                                        100281
                  IBMIN=IB1                                              040280
                ENDIF
                IF (EMAX.LT.EE(IB1,IE)) THEN                             100281
                  EMAX=EE(IB1,IE)                                        100281
                ENDIF
              ENDDO

            ENDIF

          ENDDO

C  find energy increment

          DO IE = 1,NEE(IBMIN)                                           110281
            EDIFF = EE(IBMIN,IE) - EMIN                                  110281
            IF (EDIFF.LT.DE.AND.EDIFF.GT.0.01) THEN                      110281
              DE = EDIFF                                                 110281
            ENDIF
          ENDDO
          NEMAX=INT((EMAX-EMIN)/DE+0.0001)+1                             040280
          DO IE = 1,NEMAX                                                040280
            NV(IE) = 0
            A(IE)  = 0.
          ENDDO

C  loop over same beams again, averaging over these beams, excluding gaps, and
C  reordering energies

          IEMAX = 0
          NBEAT = NBEA(IB)                                               040280

          DO IB1 = IB,NBED                                               040280
            
            IF (NBEA(IB1).EQ.NBEAT) THEN                                 040280
              DO IE = 1,NEE(IB1)

C  find energy index and maximum energy index

                E   = EE(IB1,IE)
                IEN = INT((E-EMIN)/DE + 0.0001) + 1
                IF (IEN.GT.IEMAX) THEN
                  IEMAX=IEN
                ENDIF

                IF (GCOUNT(IB1).EQ.0) THEN                               201103

C  for beams without gaps increase averaged beams counter and add intensities

                  NV(IEN) = NV(IEN) + 1
                  A(IEN) = A(IEN) + AE(IB1,IE)

                ELSE

C  for beams with gaps check whether current energy is located within a gap,
C  if not increase averaged beams counter and add intensities

                  GIN = 0                                                201103
                  IGC = IGAP                                             201103
                  DO I = 1,GCOUNT(IB1)                                   201103
                    IGC = IGC + 1                                        201103
                    IF (E.GE.EG1(IGC).AND.E.LE.EG2(IGC)) THEN            201103
                      GIN = 1                                            201103
                    ENDIF
                  ENDDO
                  IF (GIN.NE.1) THEN                                     201103
                    NV(IEN) = NV(IEN) + 1
                    A(IEN) = A(IEN) + AE(IB1,IE)
                  ENDIF

                ENDIF

              ENDDO
              NBEA(IB1) = 0                                                040280
            ENDIF

          ENDDO

C  loop over energy index to finish averaging scheme

          GIN = 0
          IG  = 0

          DO IE = 1,IEMAX

            EE(NBE,IE) = EMIN + FLOAT(IE-1)*DE

C  when there is no gap for current perform averaging, else update gap data

            IF (NV(IE).NE.0) THEN
              AE(NBE,IE) = A(IE)/FLOAT(NV(IE))
              GIN = 0
            ELSE IF (GIN.EQ.0) THEN
              IG = IG + 1
              EG1(NGAP+IG) = EE(NBE,IE)
              EG2(NGAP+IG) = EE(NBE,IE)
              NG1(NGAP+IG) = IE
              NG2(NGAP+IG) = IE
              GIN = 1
            ELSE
              EG2(NGAP+IG) = EE(NBE,IE)
              NG2(NGAP+IG) = IE
            ENDIF

          ENDDO

C  set properties of averaged beam NBE and take name of first beam encountered 

          NEE(NBE) = IEMAX
          GC(NBE)  = IG
          DO I = 1,5
            BENAME(I,NBE) = BENAME(I,IB)
          ENDDO

          NGAP = NGAP + IG
          NBE = NBE + 1

        ENDIF
        IGAP = IGAP + GCOUNT(IB)
      ENDDO

      NBE = NBE - 1
      DO IB = 1,NBE
        GCOUNT(IB) = GC(IB)
      ENDDO

C  print energies and averaged beam intensities if desired

      IF (IPR.EQ.2) THEN
        DO IB = 1,NBE
          N = NEE(IB)
          WRITE(6,10) IB,(EE(IB,IE),AE(IB,IE),IE=1,N)
        ENDDO
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------
C  Subroutine GAPIN reads informations about energy gaps in the experimental
C  spectra. If such gaps exist, this energy range will be ignored in the
C  smoothing and interpolation operations as well as for the R-factor calcu-
C  lations. The gaps must be in the order of the experimental beams and for 
C  more than one gap in a beam ordered by energy, otherwise the program stops.

      SUBROUTINE GAPIN(NGAP,BM,EG1,EG2,GCOUNT,NBED)

      INTEGER NGAP, BM(NGAP), GCOUNT(NBED)
      REAL EG1(NGAP), EG2(NGAP)

 10   FORMAT(I3,2F7.2)
 20   FORMAT(10H*THERE ARE,I3,33H GAPS IN THE EXPERIMENTAL SPECTRA)
 30   FORMAT(52H *** GAPS MUST BE IN ORDER OF EXPERIMENTAL BEAMS ***,/,
     *       52H ***               PROGRAM STOPPED               ***)
 40   FORMAT(17H *** GAPS IN BEAM,I3,30H MUST BE ORDERED BY ENERGY ***,
     *       /,50H ***              PROGRAM STOPPED              ***)
 50   FORMAT(4H GAP,I3,8H IN BEAM,I3,8H BETWEEN,F7.2,4H AND,F7.2,3H EV)

      READ (8,10) NGAP
      WRITE(6,20)NGAP

      DO IG = 1,NGAP
        READ(8,10) BM(IG),EG1(IG),EG2(IG)
        IF (IG.GT.1) THEN
          IF (BM(IG).LT.BM(IG-1)) THEN
            WRITE(6,30)
            STOP
          ELSE IF (BM(IG).EQ.BM(IG-1).AND.EG1(IG).LT.EG1(IG-1)) THEN
            WRITE(6,40) BM(IG)
            STOP
          ENDIF
        ENDIF
        WRITE(6,50) IG,BM(IG),EG1(IG),EG2(IG)
        GCOUNT(BM(IG)) = GCOUNT(BM(IG)) + 1
      ENDDO

      RETURN
      END
      
C-----------------------------------------------------------------------
C  SUBROUTINE RDPL READS LEED INTENSITIES PUNCHED OUT BY LEED            281279
C  SUBROUTINE RINT. DATA FROM DIFFERENT RUNS OF THE LEED PROGRAMS OVER
C  SUCCESSIVE ENERGY RANGES  MAY BE CONCATENATED (LEAVING SEPARATING EOF
C  MARKS) WITHOUT NEED FOR DELETION OF FIRST PART OF DATA (TITLE,BEAMS,
C  ETC), BUT USER SUPPLIED INPUT SHOULD NOT BE REPEATED. IT IS ASSUMED
C  THAT IF A RUN OF THE LEED PROGRAMS HAS NOT TERMINATED A GIVEN
C  ENERGY (I.E. DONE SOME BUT NOT ALL OF THE GEOMETRIES AT THAT
C  ENERGY), A NEXT RUN WILL HAVE STARTED AT THE UNCOMPLETED ENERGY
C  (REDOING IT ENTIRELY).
C  INDIVIDUAL GEOMETRIES MAY BE SKIPPED.
C   E= READ ENERGY GRID.
C   YDATA(K,J,I)= READ INTENSITIES  K INDEXES GEOMETRIES, J INDEXES
C    BEAMS, I INDEXES ENERGIES.
C   NS= INPUT NO. OF GEOMETRIES USED IN THE LEED COMPUTATION.
C   NB= INPUT NO. OF BEAMS FOR WHICH INTENSITIES WERE PUNCHED OUT.
C   NE= OUTPUT NO. OF ENERGIES ON GRID.
C   SPAC= NUMBERS SPECIFYING THE VARIOUS GEOMETRIES (E.G. INTERLAYER
C    SPACINGS).
C   PQ= READ BEAMS.
C   KAV= INPUT BEAM GROUPINGS (BEAMS WITH THE SAME VALUE OF KAV WILL BE
C    AVERAGED TOGETHER).
C   SYM= READ BEAM SYMMETRY CODES (USED TO DETERMINE THE AVERAGING
C    WEIGHTS).
C   EMIN,EMAX - INPUT AT ENERGIES BELOW EMIN OR ABOVE EMAX IS IGNORED.
C   LIMFIL= NO. OF FILES (SEPARATED BY EOF MARKS) TO BE READ
C   NSSK= LIST OF GEOMETRIES TO BE SKIPPED
C   NSS= NO. OF GEOMETRIES LEFT AFTER SKIPPING.
      SUBROUTINE RDPL(E,NGP,YDATA,NSS,NB,MNE,NE,SPAC,PQ,KAV,SYM,
     *                EMIN,EMAX,LIMFIL,NSSK,NS,MNBTD)
      DIMENSION E(NGP),YDATA(NSS,NB,MNE),SPAC(NSS),LAB(MNBTD),PQ(2,NB)
      DIMENSION KAV(NB),SP(NSS),YD(MNBTD),NSSK(NS)
      INTEGER SYM(NB)
C  FOR CDC ONLY  **************************
C     LEVEL 2, YDATA
C  READ BEAM NAMES (PQ), IDENTIFIERS (LAB=NPU OF LEED PROGRAM;
C  NOT USED HERE) AND SYMMETRY PROPERTIES (SYM), AS SUPPLIED BY
C  THE LEED PROGRAM
      DO 10 J=1,NB
10    READ(5,20) LAB(J),(PQ(I,J),I=1,2),SYM(J)
20    FORMAT(I5,2F10.5,I3)
      WRITE(6,27)
27    FORMAT(22H*THEOR. BEAM GROUPINGS)                                  040481
      DO 28 J=1,NB
      WRITE(6,29)(PQ(I,J),I=1,2),KAV(J)
28                                        CONTINUE
29    FORMAT(1H ,2F10.5,I5)
C  READ ENERGIES (E), GEOMETRICAL CHARACTERIZATION (SPAC) AND
C  INTENSITIES (YDATA), AS SUPPLIED BY THE LEED PROGRAM
      IFILE=1
      I=1
30    JS=0
      IF (I.GT.MNE+1) THEN
        WRITE(6,25)
 25     FORMAT(41HDIMENSION MNET TOO SMALL, PLEASE CORRECT!/
     *         35H      ****   PROGRAM STOPPED   ****)
        STOP
      ENDIF
      DO 50 J=1,NS
      READ (5,40,END=55) E(I),SP(J),(YD(K),K=1,NB)
C      WRITE(6,*) E(I),NB
40    FORMAT(F7.2,F7.4,4E14.5,/,2000(5E14.5,/))
C  CHECK FOR END OF INPUT
      DO 32 IS=1,NS
32    IF (NSSK(IS).EQ.J) GO TO 50
      JS=JS+1
      SPAC(JS)=SP(J)
      DO 35 K=1,NB
35    YDATA(JS,K,I)=YD(K)
50    CONTINUE
      IF (E(I).GT.EMAX) GO TO 60
      IF (E(I).GE.EMIN) I=I+1
      GO TO 30
55    IFILE=IFILE+1
      IF (IFILE.GT.LIMFIL) GO TO 60
C  START READING NEXT FILE, SKIPPING FIRST PART OF DATA
      NSK=NB+2
      DO 56 KK=1,NSK
56    READ(5,57)IDUM
57    FORMAT(A4)
      GO TO 30
60    NE=I-1
      WRITE(6,65)
65    FORMAT(43H*SURFACE STRUCTURE NO. AND CHARACTERIZATION)
      DO 70 J=1,NSS
70    WRITE(6,75)J,SPAC(J)
75    FORMAT(18X,I3,F15.4)
C      DO 90 IB=1,MNBTD
C      WRITE(6,80)(E(IE),YDATA(1,IB,IE),IE=1,NE)
C80    FORMAT(32H* THEO. ENERGIES AND INTENSITIES,
C     1/,83(5(1F7.2,1E13.4,3X),/))
C90    CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE READE INPUTS EXPERIMENTAL IV-CURVES

      SUBROUTINE READE(AE,EE,NBED,NGP,NEE,BENAME,IPR)
      DIMENSION AE(NBED,NGP),EE(NBED,NGP),NEE(NBED),FMT(19),
     *          BENAME(5,NBED)

C  READ INPUT FORMAT OF EXP. INTENSITIES
      READ(8,60)(FMT(I),I=1,19)
60    FORMAT(19A4)
      DO 90 IB=1,NBED
C  READ AND PRINT NAME OF THIS BEAM
      READ(8,10)(BENAME(I,IB),I=1,5)
10    FORMAT(19A4)
      WRITE(6,70)IB,(BENAME(I,IB),I=1,5)
70    FORMAT(11H*EXP. BEAM ,1I4,2H (,5A4,1H))
C  READ IN NO. OF DATA POINTS TO BE INPUT FOR THIS BEAM AND CONSTANT
C  CORRECTION FACTOR FOR INTENSITIES. THIS FACTOR IS MEANT TO ALLOW
C  NORMALIZATION TO INTENSITIES OF THE ORDER OF 1 (NOT NECESSARY, BUT
C  SAFE), AND TO MATCH UP CURVES TO BE AVERAGED TOGETHER WHEN THEIR
C  ENERGY RANGES DIFFER (TO AVOID DISCONTINUITIES AT ENERGIES WHERE
C  THE NUMBER OF CURVES AVERAGED TOGETHER CHANGES)
      READ(8,35)NEE(IB),FAC
35    FORMAT(1I4,1E14.4)
      N=NEE(IB)
C  READ (AND MAYBE PRINT) EXP. INTENSITIES
      READ(8,FMT)(EE(IB,IE),AE(IB,IE),IE=1,N)
      IF (IPR.LT.1) GO TO 82
      WRITE(6,80)(EE(IB,IE),AE(IB,IE),IE=1,N)
80    FORMAT(31H* EXP. ENERGIES AND INTENSITIES,
     1/,100(5(1F7.2,1E13.4,3X),/))
82    DO 86 IE=1,N
86    AE(IB,IE)=AE(IB,IE)/FAC
90    CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE READT READS IN THEORETICAL IV-CURVES TO BE USED AS 
C  PSEUDOEXPERIMENT
C  AUTHOR: R.DOELL, 06.11.92

      SUBROUTINE READT(AE,EE,NBED,NGP,NEE,BENAME,IPR)

      DIMENSION AE(NBED,NGP),EE(NBED,NGP),NEE(NBED),BENAME(5,NBED),
     *          EEH(NGP),IOFF(NBED)

C  IOFF     OFFSET FOR DATA POINT NUMBER

 101  FORMAT(I5,22A4)
 102  FORMAT(A8)
 103  FORMAT(I5,E13.4)
 104  FORMAT(A27)
 105  FORMAT(F7.2,F7.4,4E14.5,/,2000(5E14.5,/))
 106  FORMAT(8F10.3)

C BEGIN READ IN PART

      READ (8,103) NBEAMS
      WRITE (6,103) NBEAMS
      IF (NBED .NE. NBEAMS) THEN
          WRITE(7,'(A33)') '*********************************'
          WRITE(7,'(A33)') 'SOMETHING IS WRONG WITH NUMBER OF'
          WRITE(7,'(A33)') 'PSEUDOEXPERIMENTAL BEAMS!        '
          WRITE(7,'(A33)') 'NO AVERAGING ALLOWED!            '
          WRITE(7,'(A33)') '*********************************'
          STOP 'ERROR 1 IN PSEUDOEXPERIMENT'
      ENDIF
      DO  IB=1,NBED   
         READ(8,101) INN,(BENAME(I,IB),I=1,5)
         WRITE(6,101) INN,(BENAME(I,IB),I=1,5)
      ENDDO

C  PRESET ALL ARRAYS

      DO K=1,NBED
        NEE(K) = 0
        IOFF(K) = 0
        DO J=1,NGP    
            AE(K,J) = 0.
        ENDDO
      ENDDO
C
      DO  J=1,NGP
         READ (8,105,END=250) EEH(J),DUMMY,(AE(K,J),K=1,NBED)
C         WRITE (6,105) EEH(J),DUMMY,(AE(K,J),K=1,NBED)
         DO  K=1,NBED
            IF (AE(K,J) .GT. 1.E-20) THEN
               NEE(K) = NEE(K)+1
            ELSE IF (NEE(K) .EQ. 0) THEN
               IOFF(K) = IOFF(K)+1
            ENDIF
        ENDDO
      ENDDO

C END READ IN PART

  250 CONTINUE

      DO K=1,NBED
        DO J=1,NEE(K)
          EE(K,J)=EEH(J+IOFF(K))
          AE(K,J)=AE(K,J+IOFF(K))
        ENDDO
      ENDDO
C
      RETURN
      END
C------------------------------------------------------------------------------
C  Subroutine SCRIPT creates a gnuplot script that plots the experimen-  110205
C  tal and the best-fit theoretical spectra of each beam provided by     110205
C  exp.column and theo.column to eps files. Furthermore a latex script   110205
C  is created that integrates these eps files and additional data into   110205
C  an postscript file                                                    110205

      SUBROUTINE SCRIPT(EMIN,EMAX,INTMAX,PLSIZE,XTICS,NBE,BENAME,RB,
     *                  NBED,TITLE,BRAV,IBP)

      REAL INTMAX,PLSIZE(2),BENAME(5,NBED),RB(NBED)
      INTEGER XTICS, IBP(NBED)
      CHARACTER*80 TITLE

      OPEN(12,FILE='plot.gnu')
      OPEN(13,FILE='plot.tex')

C  write the general settings of the plots to plot.gnu

 10   FORMAT('set terminal pdfcairo enhanced font "Arial,12"',/,
     *       'set encoding iso_8859_1')
 11   FORMAT('set xrange [',I3,':',I3,']',/,'set yrange [0:',I4,']')
 12   FORMAT('set size ',F5.3,',',F5.3,/,'set style data lines',/,
     *       'set xtics ',I3,/,'set noytics',/,
     *       'set xlabel "Energy [eV]"',/,'set noxzeroaxis',/)

      WRITE(12,10)
      WRITE(12,11) INT(EMIN),INT(EMAX),INT(1.5*INTMAX)
      WRITE(12,12) PLSIZE(1),PLSIZE(2),XTICS

C  create one plot for each beam
 15   FORMAT('set output "pic_',I3.3,'.pdf"',/,'set nolabel')
 16   FORMAT('set label "',5A4,'\nR =',F7.4,'" at ',I3,
     *       ',',I4,' center font "Helvetica,12"')
 17   FORMAT('plot ''exp.column'' using ($',I3,'>0?$',I3,':1/0) :',
     *       '($',I3,'>0?$',I3,':1/0) title "(exp.)" ',
     *       'lt rgb "red" lw 2,\',/ 
     *       ' ''theo.column'' using ($',I3,'>0?$',I3,':1/0) :',
     *       '($',I3,'>0?$',I3,':1/0) title "(theory)" ',
     *       'lt rgb "blue" lw 2',/)

      DO IBE = 1,NBE
        IF (IBP(IBE).NE.0) THEN
           WRITE(12,15) IBE
           WRITE(12,16) (BENAME(I,IBE),I=1,5),RB(IBE),
     *          INT(EMIN+(EMAX-EMIN)/2),INT(1.2*INTMAX)
           WRITE(12,17) (2*IBE)-1,(2*IBE)-1,(2*IBE)-1,2*IBE,
     *          (2*IBE)-1,(2*IBE)-1,(2*IBE)-1,2*IBE
        ENDIF
      ENDDO

C  write the general settings of the output document to plot.tex and create 
C  title

 20   FORMAT('\documentclass[12pt,a4paper]{article}',/,
     *       '\usepackage{graphicx}',/,'\usepackage{float}',/)
 21   FORMAT('\begin{document}',/,'\pagestyle{empty}',/)
 22   FORMAT('\newcommand{\bigtitle}{',/,'\begin{center}',/,'\large',
     *       /,'\textbf{',A80,'}\\',/,'$\mathbf{R_{P} =',F7.4,'}$',
     *       /,'\end{center}',/,'}',/)

      WRITE(13,20)
      WRITE(13,21)
      WRITE(13,22) TITLE,BRAV

C  create a new page for each plot (beam) in the output document

 25   FORMAT('\newpage',/,'\bigtitle',/,'\begin{figure}[H]')
 26   FORMAT('\includegraphics[angle=0,width=14cm]{pic_',I3.3,
     *       '.pdf}',/,'\end{figure}',/)
 27   FORMAT('\end{document}')

      DO IBE = 1,NBE
        IF (IBP(IBE).NE.0) THEN
           WRITE(13,25)
           WRITE(13,26) IBE
        ENDIF
      ENDDO

      WRITE(13,27)

      RETURN
      END
C-----------------------------------------------------------------------
C  Subroutine STRIP2 copies energies and intensities for a given geometry IS

      SUBROUTINE STRIP2(YS,NS,NBD,NE,IS,Y,ES,ET,NGP)

      REAL YS(NS,NBD,NE),Y(NBD,NGP),ES(NGP),ET(NBD,NGP)

      DO IB = 1,NBD
        DO IE = 1,NE
          ET(IB,IE)=ES(IE)
          Y(IB,IE)=YS(IS,IB,IE)
        ENDDO
      ENDDO

      RETURN
      END
