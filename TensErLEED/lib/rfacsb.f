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
      
C------------------------------------------------------------------------------
C  Subroutine SCRIPT creates a gnuplot script that plots the experimen-  110205
C  tal and the best-fit theoretical spectra of each beam provided by     110205
C  exp.column and theo.column to eps files. Furthermore a latex script   110205
C  is created that integrates these eps files and additional data into   110205
C  an postscript file                                                    110205

      SUBROUTINE SCRIPT(EMIN,EMAX,INTMAX,PLSIZE,XTICS,NBE,BENAME,RB,
     *                  NBED,TITLE,BRAV,IBP)

      REAL INTMAX,PLSIZE(2),RB(NBED)
      CHARACTER*4 BENAME(5,NBED)
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
