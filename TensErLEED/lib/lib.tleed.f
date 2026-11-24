C  Tensor LEED subroutines for reference & delta amplitude calculation
C  v1.72, LH,MR,AMI 13.10.2021
C  for use with ref-calc.f v1.72, delta.f v1.7+
C
C  as described in
C
C  V. Blum, K. Heinz, submitted to Comp. Phys. Comm. (2000).
C
C  Please cite this reference when using this program.
C
C  Welcome to the Erlangen Tensor LEED program package. Permission is granted
C  to use these programs freely, and to modify them as you may find appropriate.
C  However, we cannot take any responsibility for possible bugs or other
C  limitations you may encounter. If you find a bug or worthwhile improvement,
C  please notify us at
C
C                       tleed@fkp.physik.uni-erlangen.de
C
C  so we can make the updated code available to everyone. Likewise, in re-
C  distributing the code please refer back to us first to ensure the latest version
C  of the package is passed on.
C
C**********************************************************************************
C
C  Please read the comment in ref-calc.f, v1.72 .
C
C  Most subroutines are derived from those described in
C
C  Van Hove / Tong, "Surface Crystallography by LEED" (Springer, Berlin 1979),
C
C  with modifications for Tensor LEED where appropriate.
C  Initial modifications were made by P.J. Rous for Tensor LEED,
C  further ones by N. Bickel, W. Oed, R. Doell, U. Loeffler, M. Kottcke and others.
C
C  For this version, subroutines were modified resulting in some
C  loss of generality compared to original VHT subroutines. E.g., the MPERTI
C  subroutine is no longer compatible with the current version.
C
C  One word of advice: If you find a bug in these subroutines, do not ignore
C  it or just patch it somehow. Look into the code, find the bug and fix it
C  once and for all (and let others using the code know). It'll help a lot
C  of people avoid hitting the same trap over and over again.
C  Another one: Please write your comments in English!!
C
C  Subroutines are included in alphabetical order

!-------------------------------------------------------------------
!  SUBROUTINE GET_TOPLAY_PROPAGATORS computes propagators for
!  the complex amplitudes of each beam:
!    WK(i, 1)
!      beam i to/from 'vacuum' from/to topmost surface layer (i.e.,
!      where the solid starts).
!    WK(i, 2)/WK(i, 3)
!      beam i going towards (2) and coming from (3) the underlying stack
!
!    Parameters
!    ----------
!    WK : complex, OUTPUT
!        The propagators are stored here.
!    PQ : real
!        Lateral momentum for all the beams in the beam list in
!        atomic units (i.e., the reciprocal-space vectors). Only
!        the first N are to be considered.
!    N : integer
!        Total number of beams included in the calculation at the
!        current energy
!    INTL_VEC : real
!        Interlayer vector (z, x, y) between top layer and the
!        underlying stack
!
!    From common
!    -----------
!    E : real
!        Current energy in the solid (including real part of inner
!        potential)
!    AK2, AK3 : real
!        Lateral momentum of incident beam in vacuum
!    VPI : real, unused
!        Imaginary part of the inner potential
!    VPIO, VPIS : real                                                 << Would actually be nice to use VPI straight (everywhere)
!        Imaginary part of the inner potential in the top layer
!        and in the rest of the stack
!    ASE : real
!        Onset 'height' of inner potential, i.e., where the solid starts
!        with respect to the geometrical position of the topmost layer.
!
!    Author: Michele Riva
!    Added: 2021-08-26

      SUBROUTINE  GET_TOPLAY_PROPAGATORS (WK, PQ, N, INTL_VEC)

      INTEGER I, N
      COMPLEX WK, XX, YY, IU
      DIMENSION WK(N, 3)
      COMPLEX SQRT, EXP
      REAL AK, BK2, BK3, PQ, INTL_VEC, X
      DIMENSION INTL_VEC(3), PQ(2, N)

      COMMON  E, AK2, AK3    !, VPI
      COMMON  /ADS/ ASE, VPIS, VPIO, VV  ! VV necessary (although unused) as this is a named COMMON block from main

      IU = CMPLX(0.0,1.0)
      AK = 2.0 * E

      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)

!       YY is k_perp in topmost layer
        YY = CMPLX(AK - BK2 * BK2 - BK3 * BK3,
     +             -2.0 * VPIO + 0.000001)
        YY = SQRT(YY)

!       X is the phase shift due to in-plane propagation between top layer
!       and stack underneath
        X = BK2 * INTL_VEC(2) + BK3 * INTL_VEC(3)

        WK(I, 1) = EXP(IU * (YY * ASE))
        WK(I, 2) = EXP(IU * (YY * INTL_VEC(1) + X))
        WK(I, 3) = EXP(IU * (YY * INTL_VEC(1) - X))
      ENDDO
      RETURN
      END

!-------------------------------------------------------------------
!  SUBROUTINE GET_TOPLAY_MULTSCATT computes the LU factorization of
!  the complex multiple-scattering matrix between the top layer and
!  the rest of the stack.
!  The LU part of the matrix is returned into MULT_SCATT, while PIVOT
!  stores the permutation pivots of the LU decomposition.
!  The return values should be used in a subsequent call to CGETRS.
!  The call to CGETRS is then equivalent to computing
!  AMP1 = MULT_SCATT*AMP2,
!  with MULT_SCATT(i, j) corresponding to the multiple scattering of
!  beam j (exiting the top layer) into beam i, which originates from
!  an infinite number of reflections between the rest of the stack
!  (RREST_ABOVE) and the top layer itself (RTOP_BELOW) through all
!  beams. Beam i is then right below the top layer.
!
!  The matrix is computed without approximations relying on the fact
!  that (see Neumann Series)
!
!          (Id - A)^(-1) = SUM(A**k) .
!
!  Thus, SUM(A**k) can be calculated by inverting Id - A.
!
!  Author: Michele Riva
!  Added: 2021-08-27
!
!  Parameters
!  ----------
!  RTOP_BELOW : complex
!     Matrix of reflection for beams reaching the top layer from below
!  RREST_ABOVE : complex
!     Matrix of reflection for beams reaching the rest of the stack
!     from above
!  WK : complex
!     Propagators between top layer and rest, as calculated via
!     GET_TOPLAY_PROPAGATORS
!  MULT_SCATT : complex
!     The multiple scattering matrix output
!  N : integer
!     Total number of beams at the current energy

      SUBROUTINE  GET_TOPLAY_MULTSCATT(RTOP_BELOW, RREST_ABOVE, WK,
     +                                 MULT_SCATT, INV_PIVOT, N)

      INTEGER  N
      COMPLEX  RTOP_BELOW, RREST_ABOVE, WK, MULT_SCATT, RU, CZ
      DIMENSION RTOP_BELOW(N, N), RREST_ABOVE(N, N)
      DIMENSION WK(N, 3)
      DIMENSION MULT_SCATT(N, N)
      INTEGER   INV_PIVOT, INV_INFO   ! For matrix inversion
      DIMENSION INV_PIVOT(N)

      RU = CMPLX(1.0, 0.0)            ! Complex one
      CZ = CMPLX(0.0, 0.0)            ! Complex zero

!     Prepare the matrix (Id - SINGLE_SCATTERING) that will later be inverted

      DO J = 1, N        ! outer loop on second index
        DO I = 1, N
          MULT_SCATT(I, J) = CZ
          DO K = 1, N
            MULT_SCATT(I, J) =
     +          MULT_SCATT(I, J)
     +          - (WK(J, 2)                ! Propagate j toward rest
     +             * RREST_ABOVE(K, J)     ! Reflect j at rest into k
     +             * WK(K, 3)              ! Propagate k back toward top
     +             * RTOP_BELOW(I, K))     ! And reflect k at top into i
          ENDDO
        ENDDO
        MULT_SCATT(J, J) = MULT_SCATT(J, J) + RU
      ENDDO

!     Gaussian-elimination step of matrix inversion (LU decomposition)
      CALL CGETRF(N, N, MULT_SCATT, N, INV_PIVOT, INV_INFO)

      RETURN
      END

!-------------------------------------------------------------------
!  SUBROUTINE ADREF2T computes the complex amplitudes for beams
!  as a result of stacking the top layer above the rest of the slab.
!  The rest of the slab already includes the substrate bulk and all
!  surface layers except the topmost one.
!  Amplitudes are calculated for ONE INCIDENT BEAM from (1) reflection
!  and transmission matrices of the top layer for beams coming from above
!  and below, (2) reflection matrix of the rest, (3) multiple-scattering
!  matrix between top and rest, and (4) propagators.
!  Only the amplitudes necessary to output Tensors (AMP0, AMP1, AMP2) are
!  calculated for a generic beam. In addition, the outgoing wavefield XI
!  is calculated when beam (0, 0) is the current beam.
!
!  2021-08-27 Michele Riva: Reimplement subroutine and CHANGE ITS
!  SIGNATURE. The previous signature was:
!     (RA, TA, RAM, TAM, RB, ST, WK, U, INT, XI, PQ, N, EMACH, S,
!      AMP0, AMP1, AMP2, PSQ1, PSQ2)
! and was called in ref-calc with
!     (RAMP, TAPP, RAPM, TAMM, REABOV, S1, WK, XS, JNT, XI, PQ, NT,
!      EMACH, ASV, APLUS, AMINUS, ADOWN, PSQ1, PSQ2)
!
!  Parameters
!  -----------
!  RA, TA : complex
!      Reflection and transmission matrices for top layer and beams
!      coming from above. RA(i, j) is reflection of incident beam j
!      (from above) into beam i (moving up). TA(i, j) is transmission
!      of incident beam j (from above) into beam i (going down).
!  RAM, TAM : complex
!      Reflection and transmission matrices for top layer and beams
!      coming from below. RAM(i, j) is reflection of incident beam j
!      (from below) into beam i (going down). TAM(i, j) is transmission
!      of incident beam j (from below) into beam i (going up). RAM is
!      not used anymore (its effect is already included in the multiple
!      scattering matrix ST), but is kept in the signature for symmetry.
!  RB : complex
!      Reflection matrix for the rest of the slab for beams impinging
!      from above (i.e., from the region between top layer and rest)
!  ST : complex
!      Multiple-scattering matrix between top layer and the rest.
!      ST(i, j) is scattering of beam j -- exiting the top layer toward
!      the rest of the slab, right below the top layer -- into beam i,
!      which is also exiting the top layer propagating toward the rest
!      of the slab. ST should e calculated (once per each 'unknown' beam)
!      with GET_TOPLAY_MULTSCATT
!  WK : complex
!      Propagators for beams. WK(i, 1) propagates beam i from/to the
!      vacuum to/from the geometrical 'height' of the topmost layer.
!      WK(i, 2)/WK(i, 3) propagate beam i down/up from/to the top
!      layer to/from the rest of the slab. WK should be calculated
!      (once per each 'unknown' beam) with GET_TOPLAY_PROPAGATORS
!  XI : complex, OUTPUT
!      The wavefield in vacuum when the (0, 0) beam impinges on the
!      whole stack from vacuum will be stored in here. While the amplitudes
!      AMP0, AMP1, and AMP2 are always computed for all the requested
!      outgoing beams (needed for Tensors), XI is computed only once, i.e.,
!      when calculating incidence of the (0, 0) beam (IIN == 1).
!  N : integer
!      Total number of beams in the current beam set (i.e., at the
!      current energy)
!  AMP0, AMP1, AMP2 : complex, OUTPUT
!      Wavefield amplitudes as per following scheme:
!      ---- VACUUM ---------------------------------------------------
!                         |                           ^ XI
!                         v AMP0                      |
!      ---- TOP LAYER ------------------------------------------------
!                                  |           ^ AMP1
!                                  v AMP2      |
!      ---------------------------------------------------------------
!      //// REST OF SLAB /////////////////////////////////////////////
!  IIN : integer
!      Index of current beam in the beam set

      SUBROUTINE  ADREF2T (RA, TA, RAM, TAM, RB, ST, ST_PIVOT, WK, XI,
     +                     N, AMP0, AMP1, AMP2, IIN)

!   Alex: moved ST_PIVOT from complex to int -> should be?
      INTEGER    I, J, N, IIN, ST_PIVOT
      COMPLEX    RA, TA, RAM, TAM, RB, ST, WK, XI, CZ
      DIMENSION  RA(N, N), TA(N, N), RAM(N, N), TAM(N,N), RB(N, N)
      DIMENSION  ST(N, N), ST_PIVOT(N), XI(N)
      DIMENSION  WK(N, 3)
      COMPLEX    AMP0(N), AMP1(N), AMP2(N)
      COMPLEX    temp

      CZ = CMPLX(0.0, 0.0)

!     2021-06-27 Michele Riva: clean up and reimplement.
!     Now uses the propagators (WK) and the multiple-scattering matrix
!     computed (only once for each 'unknown' beam) via GET_TOPLAY_PROPAGATORS
!     and GET_TOPLAY_MULTSCATT.
!     AMP2(I) is the complex amplitude of beam I, propagating towards the
!     'rest of the slab' and right above it, when beam IIN comes in with unit
!     amplitude from the vacuum.

!     Construct the amplitudes of the incoming beams, right below the top
!     layer: beam IIN produces contributions in all other beams. In the
!     meanwhile, also zero the other amplitudes.
      temp = WK(IIN, 1)             ! Current beam from vacuum
      DO I = 1, N
        AMP2(I) = temp              ! Current beam from vacuum
     +            * TA(I, IIN)      ! Transmitted into beam I at top layer
        AMP1(I) = CZ
        AMP0(I) = CZ
      ENDDO

!     Now each beam is multiple-scattered between the top and rest. AMP2(I)
!     contains then the amplitude of beam I right below the top layer.
      call CGETRS('N', N, 1, ST, N, ST_PIVOT, AMP2, N, INFO)

!     The quantity we're interested in is the amplitude right above the rest,
!     though --> Propagate it down.
      DO I = 1, N
        AMP2(I) = AMP2(I) * WK(I, 2)
      ENDDO

!     AMP1(I) is the complex amplitude of beam I, propagating towards the
!     top layer from below and right below it, after reflection with the
!     rest of the slab of AMP2.
      DO J = 1, N                     ! Outer loop on second index
        temp = AMP2(J)                ! Beam J right above the 'rest'
        DO I = 1, N
          AMP1(I) = AMP1(I)
     +              + (temp           ! Beam J right above the 'rest'
     +                 * RB(I, J)     ! Reflected at the rest into beam I
     +                 * WK(I, 3))    ! Propagated right below the top layer
        ENDDO
      ENDDO

!     AMP0 contains the actual complex amplitudes of all beams right above
!     the top layer and coming from the vacuum side, when beam IIN comes
!     from the vacuum. All amplitudes are zero (zeroed above), except for
!     the one of IIN that is just propagated from vacuum to solid.

      AMP0(IIN) = WK(IIN, 1)

!     Now, only for IIN == 1, i.e., when calculating the (0, 0) beam, we also
!     compute the total wave-field exiting the solid for all beams, as this
!     is what is used to eventually calculate the IV curves (with a call to
!     RINT). This is stored into XI. We do the calculation only once, as we
!     only need the actual intensities when the incident beam is the (0, 0)
!     one. However, the calculation of the amplitudes done above is necessary
!     for computing the Tensors, and must be done for each calculated beam.
      IF (IIN.EQ.1) THEN
        DO I = 1, N
          ! First get the straight reflection of the (0, 0) beam at
          ! the top layer, right above it
          XI(I) = WK(IIN, 1)     ! (0, 0) beam from vacuum to top layer
     +            * RA(I, IIN)   ! Reflected at top layer into beam I

          ! Then prepare the contributions of the scattering from below
          temp = CZ
          DO J = 1, N
            temp = temp
     +             + (AMP1(J)      ! Beam J right below the top layer
     +                * TAM(I, J)) ! Transmitted into beam I while going
!                                    through the top layer from below
          ENDDO
          XI(I) = (temp + XI(I)) * WK(I, 1)  ! Finally sum up, and propagate to vacuum
        ENDDO
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
C  Subroutine Bessel generates the spherical Bessel functions by direct
C  evaluation of Pendry's expression (A.27). Since the argument is complex,
C  downward recursion is said to be more stable. For the range of interest
C  here, I have tested downward recursion vs. the current subroutine, and
C  find no deviations whatsoever. This one is more straightforward though
C  less elegant. VB.

      subroutine BESSEL(BJ, Z, N1)

C  note N1 is highest angular momentum for which BJ are calculated - here
C  2 * LMAX + 1. For low LMAX, the sums required when calculating temperature-
C  dependent phase shifts do not converge fast enough. For high LMAX, 2*LMAX+1
C  is a bit of an overkill, I guess ...

      complex Z
      complex*16 BJ(N1)

      complex*16 PRE,TERM,SUM,ZSQ
      integer K,L

C  begin computation

      ZSQ = Z*Z/2.

C  L is indeed the angular momentum in question. PRE is the L-dependent
C  pre-factor.

      PRE = 1.

      DO L = 0,N1-1,1

C  K is index of series expansion.

        TERM = 1.
        SUM = 1.
        K = 1

C  repeat the loop defined by 100 ...

 100    CONTINUE

          TERM = - TERM * (ZSQ/FLOAT(K*(2*L+2*K+1)))
          SUM = SUM + TERM

          K = K+1

C  ... until TERM really small!
        IF (CABS(CMPLX(TERM)).gt.1.E-16) GO TO 100

C  evaluate j_L(Z), i.e. BJ(L+1)

        BJ(L+1) = PRE * SUM

        PRE = PRE * Z / (2.*FLOAT(L+1)+1.0)
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C  Copypos copies global array SUBPOS (containing atomic positions in
C  composite layer) to local array POS (used in RTINV), ensuring that
C  the RTINV array dimensions are correct.

      SUBROUTINE  COPYPOS(POS,SUBPOS,NLTYPE,NNSUB,NLAY,ILTYPE)

      INTEGER NLTYPE,MNSUB,NLAY
      REAL POS,SUBPOS
      DIMENSION POS(NLAY,3),SUBPOS(NLTYPE,NNSUB,3)

      DO IPOS = 1,3
        DO ILAY = 1,NLAY
          POS(ILAY,IPOS) = SUBPOS(ILTYPE,ILAY,IPOS)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
COMMENT CPPP TABULATES THE FUNCTION PPP(I1,I2,I3), EACH ELEMENT        ! TODO: some code-style edits can be incorporated from MRiva_0.56_rtinv
C       CONTAINING THE INTEGRAL OF THE PRODUCT OF THREE LEGENDRE
C       FUNCTIONS P(I1),P(I2),P(I3). THE INTEGRALS ARE CALCULATED
C       FOLLOWING GAUNT'S SUMMATION SCHEME SET OUT BY SLATER
C       ATOMIC STRUCTURE, VOL1, 309,310
C  PPP IS USED BY SUBROUTINE PSTEMP IN COMPUTING TEMPERATURE-DEPENDENT
C  PHASE SHIFTS
C  AUTHOR  PENDRY
C  DIMENSION 42 REQUIRES N1+N2+N3-1 = 4*LMAX+2 .LE. 42
      SUBROUTINE  CPPP (PPP, N1, N2, N3, NFAC,F)
      ! DOUBLE PRECISION F
      REAL*8 F, SUM  ! MRiva 2021-10-14: Fixes divergence to +-Infinity, leading in turn to NaNs in spectra
      DIMENSION  PPP(N1,N2,N3), F(NFAC)
      F(1) = 1.0
      DO 370 I = 1, NFAC-1
  370 F(I + 1) = F(I) * DBLE(I)
      DO 460 I1 = 1, N1
      DO 460 I2 = 1, N2
      DO 460 I3 = 1, N3
      IM1 = I1
      IM2 = I2
      IM3 = I3
      IF (I1-I2)  380, 390, 390
  380 IM1 = I2
      IM2 = I1
  390 IF (IM2-I3)  400, 420, 420
  400 IM3 = IM2
      IM2 = I3
      IF (IM1-IM2)  410, 420, 420
  410 J = IM1
      IM1 = IM2
      IM2 = J
  420 A = 0.0
      IS = I1 + I2 + I3 - 3
      IF (MOD(IS,2)-1)  430, 460, 430
  430 IF (IABS(IM2-IM1)-IM3+1)  440, 440, 460
  440 SUM = 0.0
      IS = IS/2
      SIGN = 1.0
      DO 450 IT = 1, IM3
      SIGN =  - SIGN
      IA1 = IM1 + IT - 1
      IA2 = IM1 - IT + 1
      IA3 = IM3 - IT + 1
      IA4 = IM2 + IM3 - IT
      IA5 = IM2 - IM3 + IT
  450 SUM = SUM - SIGN * F(IA1) * F(IA4)/(F(IA2) * F(IA3) * F(IA5) * F(
     1IT))
      IA1 = 2 + IS - IM1
      IA2 = 2 + IS - IM2
      IA3 = 2 + IS - IM3
      IA4 = 3 + 2 * (IS - IM3)
      A =  - ( - 1.0)**(IS - IM2) * F(IA4) * F(IS + 1) * F(IM3) * SUM/(
     1F(IA1) * F(IA2) * F(IA3) * F(2 * IS + 2))
  460 PPP(I1,I2,I3) = A
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE DBGT_MOD COMPUTES THE REFLECTION MATRIX FOR A PAIR OF LAYERS   ! TODO: optimize?
!  FOR INCIDENCE FROM ONE SIDE ONLY. LAYER DOUBLING IS USED. NO
!  OVERWRITING OF INPUT TAKES PLACE.
!  NOTE  SUBROUTINE DBLG PRODUCES ADDITIONNALLY THE REFLECTION FROM
!  THE OTHER SIDE OF THE LAYERS AND THE TRANSMISSIONS IN BOTH DIRECTIONS
!   RA1,TA1,RA2,TA2= REFLECTION AND TRANSMISSION MATRICES FOR FIRST
!    LAYER (ON THE SIDE OF THE INCIDENT BEAMS). 1 AND 2 REFER TO
!    INCIDENCE ON FIRST LAYER FROM SIDE OPPOSITE SECOND LAYER AND SIDE
!    TOWARDS SECOND LAYER, RESP.
!   RB1= REFLECTION MATRIX OF SECOND LAYER (E.G. SUBSTRATE).
!   RAB= OUTPUT REFLECTION MATRIX OF COMBINED LAYER.
!   ASD= INTERLAYER VECTOR FROM FIRST LAYER TO SECOND LAYER.
!   N=NO. OF BEAMS USED.
!   NA= OFFSET IN LIST PQ OF PARTICULAR SUBSET OF BEAMS USED.
!   PQ= LIST OF BEAMS G.
!   NT= TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
!   S1,S2,XS,PP,IPL= WORKING SPACES.
!   NP=DIMENSION FOR WORKING SPACES, NP.GE.N.
!   EMACH= MACHINE ACCURACY.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY.
!   VPI= IMAGINARY PART OF ENERGY.
!   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!
!  ** SPECIAL VERSION FOR TENSOR-LEED
!  ** THIS SUBROUTINE PRODUCES ALL MATRICES NEEDED FOR CALCULATION OF
!  ** AMPLITUDES AT THE CURRENT LAYER IN -X DIRECTION  AND THOSE AT
!  ** THE LAYER BELOW IN +X DIRECTION BY MULTIPLYING THEM WITH THE
!  ** AMPLITUDES OF CURRENT LAYER IN +X DIRECTION (E.G. PRODUCED
!  ** BY SUBROUTINE ADREF1T). USE MULTAMP FOR THIS PURPOSE.
!  ** DIAGRAM OF ADDITIONAL INFORMATION:
!  **
!  ** LAYER   INTERL VEC   DIFFR MATRICES     AMPLITUDES  MATRICES
!  ** --------------------------------------------------------------
!  **   >> N-1 TOP LAYERS
!  **          ASE         + RA1,TA1            + AMP0
!  **   N ----------------------------------------------------------
!  **          ASD         - RA2,TA2            - AMP0  *  S1
!  **                      + RAB                + AMP0  *  S2
!  **   N+1 --------------------------------------------------------
!  **   >> PREVIOUS CALCULATION BY LAYER-DOUBLING
!  ** --------------------------------------------------------------
!  **
!  ** TO CALCULATE AMPLITUDES OF DEEPER LAYERS USE DBGT INSTEAD OF DBG
!  ** WHICH WILL PRODUCE THE MATRICES TO BE MULTIPLIED WITH AMP'S
!  ** CALCULATED BY ADREF1T
!  **
      SUBROUTINE DBGT_MOD(RA1, TA1, RA2, TA2, RB1, RAB, ASD, N,                     **T
     &                    PQ, S1, S2, XS, PP, IPL)
      COMPLEX RAB, RA1, TA1, RB1
      COMPLEX  S1, S2, PP, XS, XX, CZ, IU
      COMPLEX TEMP
      COMPLEX RA2(N, N), TA2(N, N)
      DIMENSION RA1(N, N), TA1(N, N), RB1(N, N), RAB(N, N),
     &          S1(N, N), S2(N, N)
      DIMENSION PP(N, 2), XS(N), PQ(2, N), IPL(N)
      DIMENSION  ASD(3)
      COMMON  E, AK2, AK3, VPI

      CZ = CMPLX(0.0,0.0)
      IU = CMPLX(0.0,1.0)
      AK = 2.0 * E

!     COMPUTE INTERLAYER PROPAGATORS PP
      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3, - 2.0 * VPI + 0.000001)
        XX = SQRT(XX) * ASD(1)           ! Phase shift due to perpendicular propagation
        X = BK2 * ASD(2) + BK3 * ASD(3)  ! Phase shift due to parallel propagation
        PP(I, 1) = EXP(IU * (XX + X))
        PP(I, 2) = EXP(IU * (XX - X))
      ENDDO

!     PRODUCE multiple-scattering MATRIX TO BE INVERTED
      DO K = 1, N
        TEMP = PP(K, 1)
        DO J = 1, N
          XX = CZ
          DO L = 1, N
            XX = XX - (RA2(J, L)
     +                 * PP(L, 2)
     +                 * RB1(L, K)
     +                 * TEMP)
          ENDDO
          S1(J, K) = XX
        ENDDO
        S1(K, K) = S1(K, K) + 1.0
      ENDDO

!     GAUSSIAN ELIMINATION STEP OF MATRIX INVERSION
      CALL CGETRF(N, N, S1, N, IPL, INFO)

!     SUBSTITUTION STEP OF MATRIX INVERSION
      CALL CGETRS('N', N, N, S1, N, IPL, TA1, N, INFO)

      DO K = 1, N
        DO L = 1, N
          XS(L) = PP(L, 1) * TA1(L, K)
        ENDDO
        DO J = 1, N  ! This loop can be replaced with a call to CGEMV??
          S1(J, K) = CZ
          DO L = 1, N
            S1(J, K) = S1(J, K) + RB1(J, L) * XS(L)
          ENDDO
        ENDDO

        DO L = 1, N
          S2(L, K) = XS(L)
          S1(L, K) = PP(L, 2) * S1(L, K)
        ENDDO
        DO J = 1, N  ! TODO: This loop, with the outer K, can be replaced with a call to CGEMM('N', 'N', N, N, N, (1.0, 0.0), TA2, N, S1, N, (1.0, 0.0), RAB, N), after assigning RAB = RA1
          RAB(J, K) = RA1(J, K)
          DO L = 1, N
            RAB(J, K) = RAB(J, K) + TA2(J, L) * S1(L, K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

!---------------------------------------------------------------------------
!  SUBROUTINE DBLGAS_MOD computes the combined reflection/transmission
!  matrices of two, possibly different, bulk layers from their
!  respective reflection/transmission matrices. Both reflection and
!  transmission may be different for incidence toward the solid and
!  away from it.
!
!  Edits: MRiva 2021-09-15 -- remove GOTOs, remove useless arguments
!         NT, NP (both == N) and EAMCH (unused). Rework code to use
!         effectively LAPACK (BLAS) routine CGEMM for matrix-matrix
!         multiplication. Get rid of S3, S4, and XS, unused as a
!         result of the reworking of the code.
!
!   RAMP,TAPP= REFLECTION AND TRANSMISSION MATRICES OF FIRST LAYER FOR
!    INCIDENCE FROM SIDE OF FIRST LAYER (OVERWRITTEN BY REFLECTION AND
!    TRANSMISSION MATRICES FOR INCIDENCE FROM SIDE OF FIRST LAYER),
!    i.e., beams propagating toward the solid
!   RAPM,TAMM= REFLECTION AND TRANSMISSION MATRICES OF FIRST LAYER FOR
!    INCIDENCE FROM SIDE OF SECOND LAYER, i.e., beams propagating away
!    from the solid.
!   RBPM,TBMM= REFLECTION AND TRANSMISSION MATRICES OF SECOND LAYER FOR
!    INCIDENCE FROM SIDE OF SECOND LAYER (OVERWRITTEN BY REFLECTION AND
!    TRANSMISSION MATRICES FOR INCIDENCE FROM SIDE OF SECOND LAYER),
!    i.e., beams propagating toward the solid
!   RBMP,TBPP= REFLECTION AND TRANSMISSION MATRICES OF SECOND LAYER FOR
!    INCIDENCE FROM SIDE OF FIRST LAYER, i.e., beams propagating away
!    from the solid.
!   ASD= INTERLAYER VECTOR FROM FIRST LAYER TO SECOND LAYER.
!   N=TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
!   PQ= LIST OF BEAMS G.
!   S1,S2,IPL= WORKING SPACES.
!   PP= interlayer propagators. PP(...,1) is propagation toward the
!       solid, PP(...,2) away from it. Calculated and used internally.
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY.
!   VPI= IMAGINARY PART OF ENERGY.
!   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
      SUBROUTINE DBLGAS_MOD(RAMP, TAPP, RAPM, TAMM,
     &                      RBMP, TBPP, RBPM, TBMM,
     &                      ASD, N, PQ, S1, S2, PP,
     &                      IPL)
      COMPLEX RAMP, TAPP, RAPM, TAMM, RBPM, TBMM, RBMP, TBPP
      COMPLEX  S1, S2, PP, XX, CZ, IU, C_ONE,
     &         C_M_ONE, P_UP, P_DOWN
      DIMENSION  RAMP(N, N), TAPP(N, N), RAPM(N, N), TAMM(N, N)
      DIMENSION  RBPM(N, N), TBMM(N, N), RBMP(N, N), TBPP(N, N)
      DIMENSION  S1(N, N), S2(N, N)
      DIMENSION  PP(N, 2), PQ(2, N), IPL(N)
      DIMENSION  ASD(3)
      COMMON  E, AK2, AK3, VPI

      CZ = CMPLX(0.0, 0.0)
      IU = CMPLX(0.0, 1.0)
      C_ONE = CMPLX(1.0, 0.0)
      C_M_ONE = CMPLX(-1.0, 0.0)
      AK = 2.0 * E

!     Compute interlayer propagators PP
      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3, - 2.0 * VPI + 0.000001)
        XX = SQRT(XX) * ASD(1)            ! out-of-plane, propagation phase shift
        X = BK2 * ASD(2) + BK3 * ASD(3)   ! in-plane, propagation phase shift
        PP(I, 1) = EXP(IU * (XX + X))     ! Beam i propagating toward solid
        PP(I, 2) = EXP(IU * (XX - X))     ! Beam i propagating away from solid
      ENDDO

!     Modify some of the reflection/transmission matrices to also
!     include propagation: RAPM(i, j), RBMP(i, j), TAMM(i, j)
!     and TBPP(i, j) are now:
!
!                                                       ^ i
!                                                  TAMM |
!     --<<TOP>>------------------------------------------------------
!                  j ^ RAPM |                         j ^
!                    |      v i                         |
!                                     |      ^ i               |
!                                   j v RBMP |               j v
!     --<<BOTTOM>>---------------------------------------------------
!                                                              | TBPP
!                                                              v i
!
!     Include for all the propagation of the incoming beam (j) from
!     the other layer. This means that the matrices become then:
!
!                                                       ^ i
!                                                  TAMM |
!     --<<TOP>>------------------------------------------------------
!                    ^ RAPM |       j |                 ^    j |
!                    |      v i       |                 |      |
!                    |                |      ^ i        |      |
!                  j |                v RBMP |        j |      v
!     --<<BOTTOM>>---------------------------------------------------
!                                                              | TBPP
!                                                              v i
      DO J = 1, N
        P_UP = PP(J, 2)
        P_DOWN = PP(J, 1)
        DO I = 1, N
          RAPM(I, J) = RAPM(I, J) * P_UP
          RBMP(I, J) = RBMP(I, J) * P_DOWN
          TAMM(I, J) = TAMM(I, J) * P_UP
          TBPP(I, J) = TBPP(I, J) * P_DOWN
        ENDDO
      ENDDO

!     Prepare first multiple-scattering matrix S1, used for calculating
!     reflection/transmission for incidence from above. Multiple
!     scattering is computed exactly from S1 = Id - X, since
!     (Id - X)**(-1) == sum_n X**n. Each single-pass scattering X(i, j)
!     is beam j coming from right below top layer, scattered at bottom,
!     propagated up, and scattered at top again. This means
!     X(i, j) = sum_k [RBMP(k, j)     j from above, incl. propagation, scattered into k at bottom
!                     * RAPM(i, k)]   k then propagated to top and scattered there into i, moving down
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S1 = Id - X
     &           RAPM, N, RBMP, N,
     &           CZ, S1, N)
      DO I = 1, N
        S1(I, I) = C_ONE + S1(I, I)  ! Id - X, diagonal
      ENDDO

!     LU decomposition for matrix inversion
      CALL CGETRF(N, N, S1, N, IPL, INFO)

!     Make TAPP = inv(S1)*TAPP, i.e., TAPP(i, j) is then beam j from
!     above, transmitted through top layer, and multiple-scattered
!     into beam i, which is right below the top layer, propagating down
      CALL CGETRS('N', N, N, S1, N, IPL, TAPP, N, INFO)

!     Prepare to compute the total reflection from above (stored in RAMP).
!     First, propagate-reflect the transmitted beam at the bottom layer
!     by taking the product S1(i, j) = sum_k RBMP(i, k) TAPP(k, j),
!     corresponding to incoming beam j from above top layer, and outgoing
!     beam i right above bottom layer (propagating up).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, RBMP, N, TAPP, N, CZ, S1, N)

!     Then propagate-transmit this at the top layer by multiplying
!     TAMM(i, k) with S1(k, j). Also, add to this 'reflection due to
!     multiple scattering in between' the 'specular' reflection RAMP(i, j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TAMM, N, S1, N, C_ONE, RAMP, N)

!     Finally, also compute the total transmission from above top layer
!     to below bottom one: Propagate outgoing beam k of TAPP(k, j) to
!     bottom layer, and transmit it there. Both done by multiplying
!     with TBPP(i, k), which already includes the propagation.
!     Notice that the calculation is done in S1, and not directly in
!     TAPP, as the CZ==0 term would zero it as a buffer before it can
!     be used in the calculation. TAPP = S1 is done later at the end.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TBPP, N, TAPP, N, CZ, S1, N)


!     Prepare second multiple-scattering matrix, used for calculating
!     reflection/transmission for incidence from below. Multiple
!     scattering is computed exactly from S2 = Id - X, since
!     (Id - X)**(-1) == sum_n X**n. Each single-pass scattering is
!     beam j coming from right below top layer, scattered at bottom,
!     propagated up, and scattered at top again. This means
!     X(i,j) = sum_k [RAPM(k, j)      j from below, incl. propagation, scattered into k at top
!                     * RBMP(i, k)]   k then propagated to bottom and scattered there into i, propagating up
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S2 = Id - X
     &           RBMP, N, RAPM, N,
     &           CZ, S2, N)
      DO I = 1, N
        S2(I, I) = C_ONE + S2(I, I)  ! Id - X, diagonal
      ENDDO

!     LU decomposition for matrix inversion
      CALL CGETRF(N, N, S2, N, IPL, INFO)

!     Make TBMM = inv(S2)*TBMM, i.e., TBMM(i,j) is then beam j from
!     below, transmitted through bottom layer, and multiple-scattered
!     into beam i, which is above the bottom layer, propagating up
      CALL CGETRS('N', N, N, S2, N, IPL, TBMM, N, INFO)

!     Prepare to compute the total reflection from below (stored in RBPM).
!     First, propagate-reflect the transmitted beam at the top layer
!     by taking the product S2(i, j) = sum_k RAPM(i, k) TBMM(k, j),
!     corresponding to incoming beam j from below bottom layer, and outgoing
!     beam i right below top layer (propagating down).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, RAPM, N, TBMM, N, CZ, S2, N)

!     Then propagate-transmit it at the bottom layer by multiplying
!     TBPP(i, k) with S2(k, j). Also, add to this 'reflection due to
!     multiple scattering in between' the 'specular' reflection RBPM(i, j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TBPP, N, S2, N, C_ONE, RBPM, N)

!     Finally, also compute the total transmission from below bottom
!     layer to above top one: Propagate outgoing beam k of TBMM(k, j) to
!     top layer, and transmit it there. Both done by multiplying
!     with TAMM(i, k), which already includes the propagation.
!     Notice that the calculation is done in S2, and not directly in
!     TBMM, as the CZ==0 term would zero it as a buffer before it can
!     be used in the calculation. TBMM = S2 is done later at the end.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, TAMM, N, TBMM, N, CZ, S2, N)

      DO K = 1, N
        DO J = 1, N
          TAPP(J, K) = S1(J, K)
          TBMM(J, K) = S2(J, K)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE DEBWAL COMPUTES DEBYE-WALLER FACTORS FOR
C  DIFFRACTION FROM BEAM G(PRIME) INTO BEAM G (AND ITS SYMMETRICAL
C  COUNTERPARTS) FOR, IN GENERAL, ANISOTROPIC ATOMIC VIBRATION,
C  INCLUDING ZERO-TEMPERATURE VIBRATION.
C   NG= NO. OF BEAMS CORRESPONDING TO G.
C   G= SET OF SYMMETRICAL BEAMS G.
C   GP= INCIDENT BEAM G(PRIME).
C   E= ENERGY.
C   VPI= OPTICAL POTENTIAL
C   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
C   T= ACTUAL TEMPERATURE.
C   T0= REFERENCE TEMPERATURE FOR VIBRATION AMPLITUDES.
C   DRX= RMS VIBRATION AMPLITUDE PERPENDICULAR TO SURFACE.
C   DRY= RMS VIBRATION AMPLITUDE PARALLEL TO SURFACE.
C   D04= FOURTH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE
C    (ISOTROPIC).
C   EDW= OUTPUT DEBYE-WALLER FACTOR (EDW(1,I) FOR REFLECTION, EDW(2,I)
C    FOR TRANSMISSION).
      SUBROUTINE  DEBWAL (NG, G, GP, E, VPI, AK2, AK3, T, T0, DRX, DRY,
     1D04, EDW)
      COMPLEX SQRT
      COMPLEX EDW
      DIMENSION  GP(2), G(2,12), EDW(2,NG)

      write(6,*) "Welcome to subroutine Debwal. If you really wish to"
      write(6,*) "compute anisotropic vibrations, you should consider"
      write(6,*) "using Tensor LEED. VB, 12.06.98."

      A1 = GP(1) + AK2
      A2 = GP(2) + AK3
      CC = REAL(SQRT(CMPLX(2.0 * E - A1 * A1 - A2 * A2, - 2.0 * VPI)))
      A1 = G(1,1) + AK2
      A2 = G(2,1) + AK3
      DD = REAL(SQRT(CMPLX(2.0 * E - A1 * A1 - A2 * A2, - 2.0 * VPI)))
C  C IS PERPENDICULAR COMPONENT OF SCATTERING VECTOR FOR REFLECTION,
C  D IS SAME FOR TRANSMISSION
      C = CC + DD
      D = CC - DD
      DO 330 I = 1, NG
      A1 = GP(1) - G(1,I)
      A2 = GP(2) - G(2,I)
C  D2 IS MEAN-SQUARE VIBRATION AMPLITUDE PARALLEL TO SURFACE AT ACTUAL
C  TEMPERATURE
      D2 = DRY * DRY * T/T0
C  ZERO-TEMPERATURE VIBRATION IS NOW MIXED IN
      D2 = SQRT(D2 * D2 + D04)
      A1 = (A1 * A1 + A2 * A2) * D2
C  D2 IS SAME AS ABOVE, BUT FOR PERPENDICULAR COMPONENTS
      D2 = DRX * DRX * T/T0
      D2 = SQRT(D2 * D2 + D04)
      EDW(1,I) = CMPLX(EXP( - 0.166667 * (A1 + C * C * D2)),0.0)
  330 EDW(2,I) = CMPLX(EXP( - 0.166667 * (A1 + D * D * D2)),0.0)
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE DECIDE determines which parts of the computation are
C  necessary for the current beam NEXIT and sets calculational
C  quantities appropriately (replaces old cryptographic masterpiece
C  wavenew). Note ICRSYM is not retained.

C  criteria: NEXIT = 0 -> do normal calculation
C            NEXIT > 0
C              evanescent beam -> skip complete calculation
C              non-evanescent beam
C                does (incident!) beam belong to those (exit!) beams for
C                which layer diffraction matrices were last calculated?
C                yes -> skip fd calculation up to ADREF2T
C                no (e.g. for first TLEED beam in case of off-normal incidence)
C                    -> perform entire calculation using this beam as inc.

      SUBROUTINE DECIDE(NEXIT,KNOWN,EMERGE,PSQ1,PSQ2,AK2,AK3,AK21,AK31,
     +                  KSQ,AK2M,AK3M,NT0,RBR1,RBR2)

C  NEXIT is current beam number
C  EMERGE is 1 if current beam NEXIT propagates in vacuum, 0 otherwise
C  KNOWN is 1 if scattering matrices for current beam have already been
C        calculated, 0 otherwise
C  PSQ1,PSQ2 is lateral momentum of current beam relative to last beam for
C        which layer diffraction matrices were computed (used in ADREF2T)
C  AK2,AK3 is lateral momentum of (true) incident beam
C  AK21,AK31 is lateral momentum of incident beam for which scattering
C        matrices were last computed (set to current value if computation required)
C  KSQ   is square of total momentum of wave
C  AK2M,AK3M is array containing negative lateral momentum of all NT0 Tensor LEED beams
C  NT0   is number of TLEED beams
C  RBR1,RBR2 are superlattice reciprocal vectors

      INTEGER NEXIT
      INTEGER EMERGE,KNOWN
      REAL PSQ1,PSQ2,AK2,AK3,AK21,AK31,KSQ
      REAL AK2M,AK3M
      DIMENSION AK2M(NT0),AK3M(NT0)
      INTEGER NT0
      REAL RBR1,RBR2
      DIMENSION RBR1(2),RBR2(2)

C  local variables

      REAL DET1,DET2

C  end declarations

      IF (NEXIT.eq.0) THEN

c  incident beam emerges, and scattering is certainly not yet known.

        EMERGE = 1
        KNOWN  = 0

        PSQ1=0.
        PSQ2=0.

        AK21 = AK2
        AK31 = AK3

        write(6,*) "Computation for incident beam started."

      ELSE

        IF ( (KSQ-AK2M(NEXIT)**2.-AK3M(NEXIT)**2.) .gt. 0. ) THEN

          EMERGE = 1

          write(6,*) "Beam NEXIT = ", NEXIT," can propagate in vacuum",
     +               " - computation required."

C  AK21,AK31 are lateral components of last incident beam for which entire
C  scattering was computed. Normally, that would be incident beam. In case
C  of off-normal incidence, recomputation is necessary for the first outgoing
C  beam as well. Only for beams not present in the reference structure at all,
C  another computation is necessary - this feature is currently not implemented.

          DET1 = (AK21-AK2M(NEXIT))*RBR2(2)
     +         - (AK31-AK3M(NEXIT))*RBR2(1)
          DET2 = (AK31-AK3M(NEXIT))*RBR1(1)
     +         - (AK21-AK2M(NEXIT))*RBR1(2)

          DET1 = DET1/(RBR1(1)*RBR2(2)-RBR1(2)*RBR2(1))
          DET2 = DET2/(RBR1(1)*RBR2(2)-RBR1(2)*RBR2(1))

          DET1 = ABS(ABS(DET1)-FLOAT(IFIX(ABS(DET1)+1.E-3)))
          DET2 = ABS(ABS(DET2)-FLOAT(IFIX(ABS(DET2)+1.E-3)))

C  now, if DET1 or DET2 .ne. 0, the current beam can not be related to the
C  one for which the last computation was performed by a superlattice vector,
C  i.e. a new computation is required
C  s.w. non-perpendicular incidence can be handled now

          IF ((DET1.gt.1.E-3).or.(DET2.gt.1.E-3).or.
     +        ((NEXIT.eq.1).and.((AK2.ne.0.).or.(AK3.ne.0.)))) THEN
            KNOWN = 0
            AK21 = AK2M(NEXIT)
            AK31 = AK3M(NEXIT)
            PSQ1 = 0.
            PSQ2 = 0.
            write(6,*) "Scattering matrices not yet known",
     +                 " - recompute entire scattering."
          ELSE
            KNOWN = 1
            PSQ1 = AK21 - AK2M(NEXIT)
            PSQ2 = AK31 - AK3M(NEXIT)
C            write(6,*) "AK21=",AK21,"AK2M=", AK2M(NEXIT)
C            write(6,*) "AK31=",AK31,"AK3M=", AK3M(NEXIT)
            write(6,*) "Scattering matrices known from previous beams",
     +                 " - only compute incoming amplitudes."
          END IF

        ELSE

          EMERGE = 0
          KNOWN = 0

          write(6,*) "Beam NEXIT = ", NEXIT," can't propagate at ",
     +               "current energy - skip computation."

        END IF

      END IF

      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE DELIMIT WRITES -1 TO OUTPUTFILE FOR AMPLITUDES
C
      SUBROUTINE DELIMIT(IFILE,IFORM)
C
      IF (IFORM.EQ.0) THEN
           WRITE (IFILE) -1
      ELSE
           WRITE (IFILE,10) -1
   10 FORMAT (I5) ! changed to I5 in v1.71
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
C      ROUTINE FINAL1 EVALUATES THE AMPLITUDE OF SPHERICAL
C      WAVES INCIDENT UPON THE ATOM AT THE ORIGIN OF THE
C      TOP LAYER GIVEN THE AMPLITUDE OF PLANE WAVES INCIDENT
C      ON EITHER SIDE OF THE TOP LAYER.
C
C  A0LM(LMMAX): INITIALLY INCIDENT SPH WAVE AMPLITUDES.
C  ALM(LMMAX) : FINALLY   INCIDENT SPH WAVE AMPLITUDES.
C
C  APLUS(NT)  : PLANE WAVES INCIDENT FROM -X (I.E. WITH
C               POSITIVE KZ- KZPLUS)
C  AMINUS(NT) : PLANE WAVES INCIDENT FROM +X (I.E. WITH
C               NEGATIVE KZ- KZMINUS)
C
C  CYLM(NT,LMMAX): SPHERICAL HARMONICS FOR EACH BEAM WITH
C                  POSITIVE KZ.FOR NEGATIVE KZ MULTIPLY BY
C                  (-1.)**(L+M)
C
C  XODST(LOD,LOD): TRANSPOSE OF ODD BLOCK OF MATRIX (1-X)
C  XEVST(LEV,LEV): TRANSPOSE OF EVEN BLOCK OF MATRIX (1-X)
C
C  IPLO(LOD),IPLE(LEV) : ELIMINATION INFORMATION FROM ZGE.
C
      SUBROUTINE FINAL1(ALM,A0LM,APLUS,AMINUS,
     1CYLM,XODST,XEVST,IPLO,IPLE,LMAX,LMMAX,
     1LOD,LEV,NT,LX,AEV,AOD,EMACH)
C
C
      COMPLEX AOD(LOD),AEV(LEV)
      COMPLEX ALM(LMMAX),A0LM(LMMAX)
      COMPLEX APLUS(NT),AMINUS(NT),CYLM(NT,LMMAX)
      COMPLEX XODST(LOD,LOD),XEVST(LEV,LEV)
      COMPLEX PRE,CZ,CI,ATEMP
      INTEGER IPLO(LOD),IPLE(LEV),LX(LMMAX)
C
C
      PI=4.0*ATAN(1.0)
      CI=CMPLX(0.0,1.0)
      CZ=CMPLX(0.0,0.0)
C
C
C   EVALUATE THE INTIALLY INCIDENT SPHERICAL WAVE AMPLITUDES
C   INCIDENT UPON THE ATOM AT THE ORIGIN OF THE TOP LAYER
C   PRODUCED BY THE INCIDENT PLANE WAVES APLUS & AMINUS
C
      DO 1 L=0,LMAX
      DO 1 M=-L,L
C
      LM=(L+1)**2-L+M
      LMM=LM-2*M
C
      PRE=4.0*PI*(CI**L)
      AL=((-1.)**L)
      AM=((-1.)**M)
      A0LM(LM)=CZ
C
      DO 1 IG=1,NT
C
      ATEMP=PRE*(AM*APLUS(IG)+AL*AMINUS(IG))*CYLM(IG,LMM)
C
   1  A0LM(LM)=A0LM(LM)+ATEMP
C
  32  FORMAT(I4,2E14.4)
C
C
C
C  EVALUATE FINAL SPHERICAL WAVE AMPLITUDES ALM FROM A0LM
C  BY MULTIPLYING ON THE LEFT BY (1-X) TRANSPOSED. I.E.
C  ALM=A0LM*(1-X) = (1-X)**T * A0LM
C
C  FIRST SORT A0LM INTO ODD & EVEN PARTS.
C
C
      DO 2 N=1,LEV
      AEV(N)=A0LM(LX(N))
    2  CONTINUE
C  2  WRITE(6,93)N,AEV(N)
C
      DO 3 N=(LEV+1),LMMAX
      AOD(N-LEV)=A0LM(LX(N))
   3  CONTINUE
C  3  WRITE(6,93)N,AOD(N)
  93  FORMAT(I4,2E14.4)
C
C
C  GAUSSIAN ELLIMINATE TRANSPOSED (1-X) MATRIX       ! Michele: calls can be replaced with CGETRF /TODO
C
      CALL ZGE(XODST,IPLO,LOD,LOD,EMACH)
      CALL ZGE(XEVST,IPLE,LEV,LEV,EMACH)
C
C  INVERT (1-X) & MULTIPLY BY A0LM ON THE RIGHT      ! Michele: calls can be replaced with CGETRS? Check what ZTU does differently than ZSU. Looks like it solves the system A**T X = B for X. /TODO
C
      CALL ZTU(XODST,IPLO,AOD,LOD,LOD,EMACH)
      CALL ZTU(XEVST,IPLE,AEV,LEV,LEV,EMACH)
C
C  STORE AOD & AEV
C
C
      DO N=1,LEV
          ALM(LX(N))=AEV(N)
      ENDDO
      DO N=(LEV+1),LMMAX
          ALM(LX(N))=AOD(N-LEV)
      ENDDO
C
C
      RETURN
      END

C---------------------------------------------------------------------
C  AUTHOR P. ROUS, MODIFIED W. OED 121190
      SUBROUTINE FINALOV(TSTORE,ALM,APLUS,AMINUS,NT,NLAYER,
     1LMN,LMMAX,CAF,LMAX,LXM,E,VPI)
C
      INTEGER LXM(LMMAX)
      COMPLEX ALM(LMMAX),CSUM,RSUM,AK,CI
      COMPLEX TSTORE(2,LMN,NT),CAF(LMAX+1)
      COMPLEX APLUS(NT), AMINUS(NT)
C
      CI=CMPLX(0.0,1.0)
      AK=-0.5/SQRT(CMPLX(2.0*E,-2.0*VPI+0.000001))
C
      IN=(NLAYER-1)*LMMAX
      K=0
      DO 1 L=0,LMAX
      DO 1 M=-L,L
      K=K+1
      KLM=LXM(K)
      CSUM=CMPLX(0.0,0.0)
      DO 2 JGP=1,NT
2     CSUM = CSUM + TSTORE(1,IN+KLM,JGP)*APLUS(JGP)+
     &TSTORE(2,IN+KLM,JGP)*AMINUS(JGP)
C
      IF( CABS(CAF(L+1)).GE.1.0E-06) THEN
      CSUM=(CI**L)*CSUM/(CAF(L+1)*AK)
      ELSE
      CSUM=CMPLX(0.0,0.0)
      ENDIF
C
      ALM(K)=CSUM

 1    CONTINUE
C
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE HEAD writes header for fd. output spectrum
C  enjoy the comments :-)

      SUBROUTINE HEAD(IFILE,TITLE,NPUN,NPU,KNT,SPQF,KSYM)

      CHARACTER*80 TITLE
      INTEGER NPUN,KNT
      INTEGER NPU,KSYM
      DIMENSION NPU(NPUN),KSYM(2,KNT)
      REAL SPQF
      DIMENSION SPQF(2,KNT)

  283 FORMAT(I5,2F10.5,I3)
C LH  Adjusted to allow for beams numbers until 99999
  284 FORMAT(I5,2F10.5,I3)

      WRITE(IFILE,'(A80)') TITLE
C  PUNCH A CARD WITH NPUN
      WRITE(IFILE,283)NPUN
      DO K=1,NPUN
      N=NPU(K)
C  PUNCH A CARD FOR EACH BEAM FOR WHICH PUNCH OUTPUT IS DESIRED (THIS
C  PASSES THE BEAM IDENTITY (0,0),(1,0),(0,1),ETC. TO THE DECK OF CARDS)
      WRITE(IFILE,284)N,(SPQF(I,N),I=1,2),KSYM(1,N)
      ENDDO
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE MFOLT_SIMPLE computes the probability of reflecting and
!  transmitting a beam into another for a composite layer.
!  This simplified version removes support for beam symmetry codes
!  and registry shifts. Uses input from RTINV_SIMPLE.
!
!  This stripped down version was created on 2021-09-02 by Michele Riva
!
!  Parameters
!  ----------
!  JG : real
!      Index (within this beam set) of scattered beam for which
!      reflection/transmission are to be computed
!  NA : real
!      Offset of the current beam set within PQ
!  AK2, AK3: real
!      In-plane components of wave vector of primary beam
!  CYLM : complex array, shape (NT, LMMAX)
!      Spherical harmonics for all beams (first index) at all
!      (L, M) combinations.
!  N : INTEGER
!      Total number of beams in the current beam set (same for both
!      incident and scattered)
!  LMMAX : INTEGER
!      Maximum number of (L, M) combinations
!  PQ : real array, shape (2, NT)
!      List of all beams, i.e., reciprocal lattice vectors
!  TS : complex vector, length LMN == NLAY*LMMAX
!      Multiple-scattered amplitudes from each of the (L, M)
!      combinations for each of the sublayers in the composite layer.
!      indices 1...LMMAX contain the (L, M) combinations for layer 1.
!      The other layers follow. This is calculated inside RTINV_SIMPLE.
!  LMN : INTEGER
!      == NLAY*LMMAX. Total number of (L, M) combinations for
!      all sublayers.
!  RG : complex array, shape (2, NLAY, NT)
!      Plane-wave propagators between subplanes. RG(1, ...)
!      corresponds to propagation into the solid; RG(2, ...)
!      towards vacuum;
!      RG(..., ILAY, IBEAM) corresponds to propagation of beam IBEAM
!      of the current beam set up to sublayer ILAY.
!  RG_R, RG_T : complex vectors, length NLAY
!      Contributions to reflection (_R) and transmission (_T) of
!      pure propagation to layer ILAY, relative to incidence/exit
!      at the extremes of the composite layer. These are calculated
!      and used internally in MFOLT_SIMPLE.
!  NLAY : INTEGER
!      Total number of sublayers in the current composite layer.
!  JGP : INTEGER
!      Index of incident beam for which reflection/transmission
!      are computed
!  LXM : integer vector, length LMMAX
!      Permutations of (L, M) sequence
!  RT_OUT : complex, OUTPUT
!      Reflection/transmission for the composite layer. Reflection
!      is RT_OUT(1), transmission is RT_OUT(2)
!  INC : {-1, +1}
!      Direction of the incident beam for which reflection and
!      transmission should be computed. +1 corresponds to
!      an incident beam moving toward the solid, -1 away from it.
!  POSS : real array, shape (NLAY, 3)
!      Atomic positions of each of the Bravais sublayers, sorted
!      from the surface towards the inside of the solid. Order
!      of coordinates is: POSS(..., 1) = z (perp, to surface),
!      POSS(..., 2/3) = x/y (in-plane).
      SUBROUTINE  MFOLT_SIMPLE (JG, NA, AK2, AK3, CYLM, N, LMMAX, PQ,
     +                          NT, TS, LMN, RG, RG_R, RG_T, NLAY, JGP,
     +                          LXM, RT_OUT, INC, POSS)
      INTEGER JG, NA, N, LMMAX, NT, LMN, NLAY, JGP, INC,  ! Input
     +        JGA                                         ! Internal use. TODO: update
      COMPLEX  CYLM, CZ, RU, CI, RT_OUT, ST, SM,
     +         CY, CTR, CTT, R, T, CR, CT, CS
      COMPLEX  EDW_R, EDW_T
      COMPLEX  EXP
      COMPLEX  TS(LMN), RG(2, NLAY, N)
      COMPLEX  RG_R(NLAY), RG_T(NLAY)
      REAL     K_IN, K_OUT          ! K_IN/K_OUT are total wave vectors of incident and scattered beams
      DIMENSION  GP(2), PQ(2, NT)   ! GP is the reciprocal-space part of the wave vector of the incident beam
      DIMENSION  CYLM(NT, LMMAX)
      DIMENSION  LXM(LMMAX), POSS(NLAY, 3), RT_OUT(2)
      DIMENSION  K_IN(2), K_OUT(2)

      REAL BR, BT

      COMMON /MFB/ GP, LM

CDIR$ NOVECTOR
      CZ = (0.0, 0.0)
      RU = (1.0, 0.0)
      CI = (0.0, 1.0)

      JGA = JG + NA
      K_OUT(1) = PQ(1, JGA) + AK2
      K_OUT(2) = PQ(2, JGA) + AK3
      K_IN(1) = GP(1) + AK2
      K_IN(2) = GP(2) + AK3

      IF (INC.EQ.1) THEN
         I1 = 1
         I2 = NLAY
      ELSE
         I1 = NLAY
         I2 = 1
      ENDIF

!     BR is the total in-plane propagation phase shift due to reflection
!     at the first layer encountered during propagation (I1): both K_IN
!     and K_OUT at first.
      BR =   (K_IN(1) - K_OUT(1)) * POSS(I1, 2)
     +     + (K_IN(2) - K_OUT(2)) * POSS(I1, 3)

!     BT is the total in-plane propagation phase shift due to transmission
!     between the first (I1) and last (I2) layers encountered: K_IN
!     at first, K_OUT at last
      BT =    K_IN(1)  * POSS(I1, 2) + K_IN(2)  * POSS(I1, 3)
     +     - (K_OUT(1) * POSS(I2, 2) + K_OUT(2) * POSS(I2, 3))

!     EDW_R/T are the in-plane propagators for reflection and
!     transmission
      EDW_R = EXP(CI * BR)
      EDW_T = EXP(CI * BT)

!     Prepare coefficients taking into account the effect of
!     propagation on the scattered amplitudes. These will then
!     be later combined with the relevant L,M contributions.
!        RG_R is effect on REFLECTION, RG_T effect on TRANSMISSION
!     For RG: RG(1, ...) is propagation towards solid, RG(2, ...) away
      DO J = 1, NLAY
        IF (INC.EQ.1) THEN
          ST = RG(1, J, JGP) / RG(1, 1, JGP)          ! Fraction of incident beam (JGP) reaching plane J
!                                                       (relative to the one reaching first layer) due to
!                                                       propagation. Both propagating down.
          RG_R(J) = ST
     +              * RG(2, J, JG) / RG(2, 1, JG)     ! Fraction of scattered beam (JG) leaving layer J
!                                                       (relative to the one reaching first layer) due to
!                                                       propagation. Both propagating up.
          RG_T(J) = ST
     +              * RG(1, NLAY, JG) / RG(1, J, JG)  ! Fraction of scattered beam (JG) leaving last layer
!                                                       (relative to the one reaching layer J) due to
!                                                       propagation. Both propagating down.
        ELSE
          ST = RG(2, NLAY, JGP) / RG(2, J, JGP)       ! Fraction of incident beam (JGP) reaching last
!                                                       layer relative to the one at plane J, due to
!                                                       propagation. Both propagating up.
          RG_R(J) = ST
     +              * RG(1, NLAY, JG) / RG(1, J, JG)  ! Fraction of scattered beam (JG) leaving last layer
!                                                       relative to the one reaching plane J, due to
!                                                       propagation. Both propagating up.
          RG_T(J) = ST
     +              * RG(2, J, JG) / RG(2, 1, JG)     ! Fraction of scattered beam (JG) leaving plane J
!                                                       relative to the one leaving first layer, due to
!                                                       propagation. Both propagating down.
        ENDIF
      ENDDO

!     Combine the pure effect of propagation with the (L,M)-dependent
!     scattering probabilities.
!       Some optimization from older code: it used to include products
!       of type SM*SL, with
!           SL = -1*(-1)**L,
!       and
!           SM = -SL * (-1)**M = (-1)**(L+M)
!       Thus, SM*SL == -1 * (-1)**(2L+M) == -1*(-1)**M. This quantity
!       is now called SM.
      R = CZ
      T = CZ
      JLM = 1
      L1 = 1
      DO L = 1, LM
        SM = RU
        LL = L + L - 1
        L1 = L1 + LL
        L2 = L1 - 2*L
        DO M = 1, LL
          KLM = LXM(JLM)
!         Sum over all sublayers
          CR = CZ
          CT = CZ
          DO J = 1, NLAY
            CS = TS((J - 1) * LMMAX + KLM)  ! Scattering at layer J with angular momentum numbers (L, M)
            CR = CR + RG_R(J) * CS          ! Including the propagation part on reflection
            CT = CT + RG_T(J) * CS          ! And the propagation part on transmission
          ENDDO

!         Finally include spherical harmonics and some more prefactors
          CY = CYLM(JG, L2 + M)
          IF (INC.EQ.1) THEN
            CTR = SM
            CTT = RU
          ELSE
            CTR = RU
            CTT = SM
          ENDIF
          R = R + CTR * CY * CR
          T = T + CTT * CY * CT
          JLM = JLM + 1
          SM =  - SM
        ENDDO
      ENDDO

      RT_OUT(1) = R * EDW_R
      RT_OUT(2) = T * EDW_T

      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE MULTAMP_OPT multiplies one vector into a matrix.
!  Optimized for Fortran memory access - by MR v1.72
!
      SUBROUTINE MULTAMP_OPT(INAMP,OUTAMP,MATRIX,N)
      COMPLEX INAMP(N),OUTAMP(N),MATRIX(N,N)
      COMPLEX CZ, TEMP

      CZ = (0.,0.)
      DO I=1,N
        OUTAMP(I) = CZ
      ENDDO

      DO J=1,N
        TEMP = INAMP(J)
        DO I=1,N
          OUTAMP(I) = OUTAMP(I) + MATRIX(I,J) * TEMP
        ENDDO
      ENDDO
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE OPENOUT OPENS THE OUTPUTFILE FOR THE AMPLITUDES
C  AND SUPPLIES ALL NEEDED GLOBAL INFORMATION
C
      SUBROUTINE OPENOUT(IFILE,FILENAM,IFORM)
      CHARACTER*(*) FILENAM
C
      IF (IFORM.EQ.0) THEN
      OPEN (IFILE,FILE=FILENAM,FORM='UNFORMATTED',STATUS='NEW')
      ELSE
      OPEN (IFILE,FILE=FILENAM,FORM='FORMATTED')
      ENDIF
      RETURN
      END

C---------------------------------------------------------------------
C  SUBROUTINE OUTAMP WRITES THE OUTPUT AMPLITUDES TO FILE# IFILE
C
      SUBROUTINE OUTAMP(IFILE,IFORM,NEXIT,PQ1,PQ2,ALM,LMMAX,AK21,AK31)
      COMPLEX ALM(LMMAX)
C
      IF (IFORM.EQ.0) THEN
C  UNFORMATTED OUTPUT
      WRITE (IFILE) NEXIT
      WRITE (IFILE) PQ1,PQ2,AK21,AK31,ALM
      ELSE
C  FORMATTED OUTPUT
      WRITE (IFILE,10) NEXIT
      WRITE (IFILE,20) PQ1,PQ2,AK21,AK31
      WRITE (IFILE,20) (ALM(I),I=1,LMMAX)
   10 FORMAT (I5)
   20 FORMAT(5E12.6)
      ENDIF
      RETURN
      END

C---------------------------------------------------------------------
C  OUTXIST STORES THE AMPLITUDES OF THE INTEGER ORDER
C  BEAMS AFTER THEY HAVE BEEN CALCULATED FOR THE INCIDENT BEAM
C
C  XI(NT): AMPLITUDES OF BEAMS FOR CURRENT ENERGY
C  XIST(NT0): ORIGINAL AMPLITUDES OF EACH EXIT BEAM
C                FOR EACH ENERGY: =0 IF IG=FRACTIONAL
C                ORDER BEAM.
C
      SUBROUTINE OUTXIST(IFILE,IFORM,E,PQF,SPQF,NPU,NT0,NT,XI,XIST,
     1                    L1,CAF)
C
      DIMENSION PQF(2,NT),SPQF(2,*),NPU(NT0)
      COMPLEX XI(NT), XIST(NT0), CAF(L1)
      COMMON /ADS/ ASE,VPIS,VPIO,VV ! note: V0 removed!

C   START LOOP OVER EXIT BEAMS.

       DO IG=1,NT0
         XIST(IG)=CMPLX(0.0,0.0)
         DO JG=1,NT
           P1=ABS(PQF(1,JG)-SPQF(1,NPU(IG)))
           P2=ABS(PQF(2,JG)-SPQF(2,NPU(IG)))
           IF((P1+P2).LT.0.000001) THEN
             XIST(IG)=XI(JG)
           ENDIF
         ENDDO
       ENDDO 

C  NOW WRITE ALL AMPLITUDES TO FILE # IFILE

      IF (IFORM.EQ.0) THEN
C  UNFORMATTED OUTPUT
        WRITE (IFILE) E,VPIS,VPIO,VV,L1,CAF,XIST
      ELSE
C  FORMATTED OUTPUT
        WRITE(IFILE,20) E,VPIS,VPIO,VV
        WRITE(IFILE,10)L1
        WRITE(IFILE,20)(CAF(I),I=1,L1)
        WRITE(IFILE,20)(XIST(I),I=1,NT0)
   10   FORMAT (I5)
   20   FORMAT (5E16.10)
      ENDIF

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE PSTEMP INCORPORATES THE THERMAL VIBRATION EFFECTS IN THE ! TODO: some code-style edits can be incorporated from MRiva_0.56_rtinv
C  PHASE SHIFTS, THROUGH A DEBYE-WALLER FACTOR. ISOTROPIC VIBRATION
C  AMPLITUDES ARE ASSUMED.
C   PPP= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C   N3= NO. OF INPUT PHASE SHIFTS.
C   N2= DESIRED NO. OF OUTPUT TEMPERATURE-DEPENDENT PHASE SHIFTS.
C   N1= N2+N3-1.
C   DR0= 4TH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDES.
C   DR= ISOTROPIC RMS VIBRATION AMPLITUDE AT REFERENCE TEMPERATURE T0.
C   T0= ARBITRARY REFERENCE TEMPERATURE FOR DR.
C   TEMP= ACTUAL TEMPERATURE.
C   E= CURRENT ENERGY (REAL NUMBER).
C   PHS= INPUT PHASE SHIFTS.
C   DEL= OUTPUT (COMPLEX) PHASE SHIFTS.
C   11.01.95 UL: uses BESSEL from TLEED package to calculate spherical
C   bessel-functions. Yield better convergence for higher vibration amplitudes
C   than original calculation scheme                                     110195
      SUBROUTINE  PSTEMP (PPP, N1, N2, N3, DR0, DR, T0, TEMP, E, PHS,
     1DEL,CTAB,SUM,BJ)
      COMPLEX  DEL, SUM, CTAB
      COMPLEX*16 BJ(N1)
      COMPLEX  Z, CI, CS, CL
      COMPLEX EXP,LOG
      DIMENSION  PPP(N1,N2,N3), PHS(N3), DEL(N2), SUM(N2)
      DIMENSION  CTAB(N3)

      PI = 3.14159265
      CI = CMPLX(0.0,1.0)
      DO 170 J = 1, N2
  170 DEL(J) = (0.0,0.0)
      ALFA = DR * DR * TEMP/T0
      ALFA = 0.166667 * SQRT(ALFA * ALFA + DR0)
      FALFE =  - 4.0 * ALFA * E
      IF (ABS(FALFE)-1.0E-3)  180, 200, 200
  180 DO 190 J = 1, N3
  190 DEL(J) = CMPLX(PHS(J),0.0)
      GO TO 360
C
COMMENT BJ(N1) IS LOADED WITH SPHERICAL BESSEL FUNCTIONS OF
C       THE FIRST KIND; ARGUMENT Z
  200 Z = FALFE * CI
c
c     use subroutine BESSEL
      call bessel(BJ,Z,N1)
c
  270 CS = CMPLX(1.0,0.0)
      FL = 1.0
      DO 280 I = 1, N1

      BJ(I) = EXP(FALFE) * FL * CS * BJ(I)

      FL = FL + 2.0
  280 CS = CS * CI
C
      FL = 1.0
      DO 290 I = 1, N3
      CTAB(I) = (EXP(2.0 * PHS(I) * CI) - (1.0,0.0)) * FL
  290 FL = FL + 2.0
C
      ITEST = 0
      LLLMAX = N2
      FL = 1.0
      DO 350 LLL = 1, N2
      SUM(LLL) = CMPLX(0.0,0.0)
      DO 300 L = 1, N3
      LLMIN = IABS(L - LLL) + 1
      LLMAX = L + LLL - 1
      DO 300 LL = LLMIN, LLMAX
  300 SUM(LLL) = SUM(LLL) + PPP(LL,LLL,L) * CTAB(L) * BJ(LL)

cvb  now, sum(LLL) is already the temperature-dependent t-matrix we were
C    looking for. It is next converted to a temp-dependent phase shift, only
C    to be converted back right after the pstemp call in tscatf. Kept for
C    the sake of compatibility with van Hove / Tong book only.

      DEL(LLL) =  - CI * LOG(SUM(LLL) + (1.0,0.0))/(2.0,0.0)
      ABSDEL = CABS(DEL(LLL))
      IL = LLL - 1
      IF (ABSDEL-1.0E-2)  320, 310, 310
  310 ITEST = 0
      GO TO 350
  320 IF (ITEST-1)  340, 330, 340
  330 LLLMAX = LLL
      GO TO 360
  340 ITEST = 1
  350 FL = FL + 2.0
  360 RETURN
      END
!-----------------------------------------------------------------------
!  SUBROUTINE RINT_SIMPLE COMPUTES THE REFLECTED BEAM INTENSITIES FROM THE
!  (COMPLEX) REFLECTED AMPLITUDES. INCLUDED ARE ANGULAR PREFACTORS AND A
!  SYMMETRY-RELATED PREFACTOR (SO THAT THE RESULTING INTENSITY IS
!  DIRECTLY COMPARABLE WITH EXPERIMENT). RINT PRINTS OUT NON-ZERO
!  INTENSITIES AND PUNCHES OUT (IF REQUESTED) THE INTENSITIES OF THE
!  BEAMS SPECIFIED IN NPUC.
!
!  This is a simplified version that removes support for symmetry codes
!  FK 2021-11-17: Add output of complex amplitudes in second file IFILE2
!
!   N= NO. OF BEAMS AT CURRENT ENERGY.
!   WV= INPUT REFLECTED AMPLITUDES.
!   AT= OUTPUT REFLECTED INTENSITIES.
!   PQ= LIST OF RECIPROCAL-LATTICE VECTORS G (BEAMS).
!   PQF= SAME AS PQ, BUT IN UNITS OF THE RECIPROCAL-LATTICE CONSTANTS.
!   SYM= LIST OF SYMMETRY CODES FOR ALL BEAMS.
!   VV= VACUUM LEVEL ABOVE SUBSTRATE MUFFIN-TIN CONSTANT.
!   THETA,FI= POLAR AND AZIMUTHAL ANGLES OF INCIDENT BEAM.
!   MPU= NO. OF BEAM INTENSITIES TO BE PUNCHED.
!   NPUC  INDICATES WHICH BEAM INTENSITIES ARE TO BE PUNCHED.
!   EEV= CURRENT ENERGY IN EV ABOVE VACUUM LEVEL.
!   A= STRUCTURAL PARAMETER OR OTHER IDENTIFIER TO BE PUNCHED ON CARDS.
!   NPNCH=0  NO PUNCH DESIRED.
!   NPNCH.NE.0  PUNCH DESIRED.
!   XIST= TEMPORARY STORAGE OF OUTPUT AMPLITUDES
!  IN COMMON BLOCKS
!   E= CURRENT ENERGY IN HARTREES ABOVE SUBSTRATE MUFFIN-TIN CONSTANT.
!   VPI= IMAGINARY PART OF CURRENT ENERGY.
!   CK2,CK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!   ANGSCALE= angular scaling factor per beam
      SUBROUTINE RINT_SIMPLE(N,WV,AT,ATP,PQ,PQF,VV,THETA,FI,
     &                       MPU,NPUC,EEV,AP,NPNCH,XIST,
     &                       IFILE,IFILE2)
      DIMENSION  WV(N), PQ(2,N), PQF(2,N), AT(N)
      DIMENSION  NPUC(MPU),ATP(MPU),XIST(MPU)
      COMPLEX WV, XIST
      REAL ANGSCALE
      DIMENSION ANGSCALE(N)
      COMMON  E, CK2, CK3, VPI

   35 FORMAT(1F7.2,1F7.4,4E14.5,/,100(5E14.5,/))                         140779
   36 FORMAT(1F7.2,1F7.4,8E18.10,/,100(10E18.10,/))
      AK = SQRT(AMAX1(2.0 * E - 2.0 * VV, 0.0))
      BK2 = AK * SIN(THETA) * COS(FI)
      BK3 = AK * SIN(THETA) * SIN(FI)

      C = AK * COS(THETA)  ! == k_perp in vacuum for incident beam
      DO J = 1, N
        AK2 = PQ(1, J) + BK2
        AK3 = PQ(2, J) + BK3
        A = 2.0 * E - 2.0 * VV - AK2 * AK2 - AK3 * AK3    ! == ||k_perp||^2
        IF (A.GT.0.0) THEN  ! beam is non-evanescent --> store its intensity for printing later
          A = SQRT(A)
          WR = REAL(WV(J))
          WI = AIMAG(WV(J))
!         AT IS REFLECTED INTENSITY (FOR UNIT INCIDENT CURRENT)
          AT(J) = (WR * WR + WI * WI) * A/C
          ANGSCALE(J) = SQRT(A/C)    ! for amplitude output below
        ELSE
          AT(J) = 0.0
          ANGSCALE(J) = 0.0   ! this will set printed complex amplitudes to zero for evanescent beams
        ENDIF
      ENDDO

      IF (NPNCH.EQ.0) RETURN  ! nothing to print to file

      DO K = 1,MPU
        NPC = NPUC(K)
        IF (NPC.EQ.0) THEN  ! This beam was not requested
          ATP(K) = 0.0
          XIST(K) = (0.0, 0.0)
        ELSE
          ATP(K) = AT(NPC)
          XIST(K) = WV(NPC) * ANGSCALE(NPC)
        ENDIF
      ENDDO

      WRITE(IFILE, 35) EEV, AP, (ATP(K), K=1,MPU)
      WRITE(IFILE2, 36) EEV, AP, (XIST(K), K=1,MPU)
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE RTINV_SIMPLE computes reflection and transmission matrices for
!  an atomic layer consisting of NLAY subplanes. A Beeby-type matrix
!  inversion is used (a perturbation treatment is done in subroutine
!  MPERTI -- which does not exist any longer).
!  This is a simplified version of the original RTINV that removes
!  support for (1) Beam symmetry codes (which were never used anyway),
!  (2) registry shifts SO(I) and SS(I) (which were also never fully
!  implemented), and (3) NOPT==1.
!
!  Parameters
!  ----------
!  R_TOP, T_TOP, R_BOT, T_BOT
!      Reflection and transmission matrices for beams incident from
!      above (_TOP) or below (_BOT) the composite layer. R/T(i, j) is
!      reflection/transmission of incident beam j into beam i. The
!      complex amplitudes of both beams are taken right next to the
!      composite layer, i.e., R_TOP(i, j): beams right above;
!      T_TOP(i, j): j right above, i right below; R_BOT(i, j): both
!      beams below; T_BOT(i, j): j right below, i right above.
!
!
!  SPECIFICATIONS NOT FOUND BELOW CAN BE FOUND IN THE SPECIFICATIONS
!  OF SUBROUTINE MPERTI (van Hove/Tong 1979), WHICH HAS THE SAME STRUCTURE
!  AS MTINV.
!  NOTE  REAL AND IMAGINARY PARTS OF MUFFIN-TIN CONSTANT ARE ASSUMED
!  CONSTANT THROUGHOUT LAYER (I.E. BETWEEN OUTERMOST NUCLEAR PLANES).
!   TAUG,TAUGM HAVE DIMENSION LMT=NTAU*LMMAX HERE.
!   TH(LMNI,LMNI2) IS SQUARE HERE WITH LMNI2=LMNI=NLAY*LMMAX.
!   IPL IS ADDITIONAL WORKING SPACE.
!
!  Author of edits: Michele Riva, 2021-09-02

!  REMOVED from arguments:
!    * extra reflection/transmission for other registry shifts
!    * SYM
!    * ID
!    * LX (unused)
!    * NM, NP --> both always == NT
!  (others)
!    * /SL/ common block
      SUBROUTINE RTINV_SIMPLE(
     &                 R_TOP,T_TOP,R_BOT,T_BOT, N,
     &                 AMULT, CYLM, PQ, NT,
     &                 LT, LXM, LMMAX, LEV,
     &                 TAU, LMT, TAUG, TAUGM, POSS, MGH,
     &                 NLAY, GH, LMG, RG, RG_PERP, TS, LMN, LM2N, VT,
     &                 TH, LMNI, IPL, LPSS, NORD,
     &                 NNSUB, RG1, RG2, TSTORE)
      COMPLEX R_TOP, T_TOP, R_BOT, T_BOT, CYLM, AMULT, XA, ST,
     +        CI, RU, CZ, AM, MFOLT_RT
      COMPLEX SQRT,EXP
      COMPLEX TH(LMNI, LMNI), TSTORE(2, LMN, NT)
      COMPLEX TAU(LMT, LEV), GH(LMG, LMMAX), RG(2, NLAY, N)            ! May be worth to have 2 as last?
      COMPLEX RG_PERP(N)
      COMPLEX TS(LMN), VT(LM2N), TAUG(LMT), TAUGM(LMT)
      DIMENSION  R_TOP(NT, NT), T_TOP(NT, NT), R_BOT(NT, NT),
     +           T_BOT(NT, NT), AMULT(N)
      DIMENSION  CYLM(NT, LMMAX), PQ(2, NT)                            ! Swap CYLM indices?
      DIMENSION  GP(2)
      DIMENSION  LT(LMMAX), LXM(LMMAX)
      DIMENSION  POSS(NLAY, 3)
      DIMENSION  MGH(NLAY, NLAY)
      DIMENSION  IPL(LMNI), MFOLT_RT(2)
CVB
      INTEGER NNSUB
      INTEGER LPSS, NORD
      DIMENSION LPSS(NNSUB), NORD(NNSUB)
      COMPLEX RG1(NLAY), RG2(NLAY), TEMP

      COMMON  E, AK2, AK3, VPI
      COMMON /MFB/ GP, L1
      COMMON /MS/  LMAX, EPS1, LITER
      COMMON /MPT/ NA, NS, ID, LAY, LM, NTAU, TST, TV, DCUT, NOPT, NEW ! NOPT always == 1 for Tensor LEED; ID == 1 (registry shifts), not used anymore

      PI = 3.14159265
      L1 = LM
      LOD = LMMAX - LEV
      LEE = (LMAX/2 + 1)**2
      LOE = ((LMAX - 1)/2 + 1)**2
      CI = CMPLX(0.0, 1.0)
      CZ = CMPLX(0.0, 0.0)
      RU = CMPLX(1.0, 0.0)

!     The subplanes are sorted (POSS, LPSS, NORD), the interlayer
!     propagators in (LM)-space (GH, MGH) are computed, and TAU is
!     computed for each chemical element, by the caller: SORT_SUBPLANES,
!     INTERPLANAR_PROPAGATORS, SUBPLANE_SCATTERING.

!     Create matrix TH to be inverted in Beeby\s scheme
      CALL THMAT(TH, LMNI, GH, LMG, LMMAX, MGH, NLAY, TAU, LMT,
     &           LEV, LPSS)

!     Prepare TH for inversion (Gaussian elimination)
      CALL CGETRF(LMNI, LMNI, TH, LMNI, IPL, INFO)

!     If NEW=1 (i.e. new energy) compute new spherical harmonics for
!     the beam directions. NB: NEW is always == 1! AMULT is overwritten
!     below.
      IF (NEW.EQ.1) THEN
        CALL BEAMS_YLM(CYLM, AMULT, PQ, NT, N, NA, E, AK2, AK3, VPI,
     &                 LMAX, LMMAX, VT)
      ENDIF
      DO IG = 1, N
!       Generate plane-wave propagators for use as prefactors for GH and TH
        JG = IG + NA
        BK2 = PQ(1, JG) + AK2
        BK3 = PQ(2, JG) + AK3
        C = BK2 * BK2 + BK3 * BK3
        XA = CMPLX(2.0 * E - C, - 2.0 * VPI + 0.000001)
        XA = SQRT(XA)
        RG_PERP(IG) = EXP(CI * XA * (POSS(NLAY, 1) - POSS(1, 1)))
        DO I = 1, NLAY
          X = BK2 * POSS(I, 2) + BK3 * POSS(I, 3)
          RG(1, I, IG) = EXP(CI * (XA * POSS(I, 1) + X))
          RG(2, I, IG) = EXP(CI * (XA * POSS(I, 1) - X))
        ENDDO
        IF (NEW.EQ.1) THEN
!         Generate prefactor of reflection and transmission matrix elements
          AMULT(IG) =  - 16.0 * PI * PI * CI / (TV * XA)
        ENDIF
      ENDDO

!     Start loop over incident beams
      DO JGP = 1, N
!       Generate quantities TAUG, TAUGM into which inverse of TH
!       will be multiplied
        CALL TAUY_SIMPLE(TAUG, TAUGM, LMT, TAU, LMT, LEV, CYLM, NT,
     &                   LMMAX, LT, NTAU, LOD, LEE, LOE, JGP)

        JGPS = JGP + NS
        GP(1) = PQ(1, JGP + NA)
        GP(2) = PQ(2, JGP + NA)
!       First consider incidence towards +x, i.e., beams propagating
!       towards the solid
        INC = 1
  285   CONTINUE  ! Will then continue from here with INC=-1

!       Include appropriate plane-wave propagating factors into TAUG/TAUGM
        DO I = 1, NLAY
          IN = (I - 1) * LMMAX
          LP = (LPSS(I) - 1) * LMMAX
          IF (INC.EQ.1) THEN
            ST = RG(1, I, JGP)
            DO K = 1, LMMAX
              TS(IN + K) = TAUG(LP + K) * ST
            ENDDO
          ELSE
            ST = RU / RG(2, I, JGP)
            DO K = 1, LMMAX
              TS(IN + K) = TAUGM(LP + K) * ST
            ENDDO
          ENDIF
        ENDDO

!       Do left multiplication with inverse of TH
        CALL CGETRS('N', LMNI, 1, TH, LMNI, IPL, TS, LMNI, INFO)

!       Include further plane-wave propagating factors, and store
!       the intra-layer scattering for computing Tensors.
        DO I = 1, NLAY
          IN = (I - 1) * LMMAX
          INSORT = (NORD(I) - 1) * LMMAX
          IF (INC.EQ.1) THEN
            ST = RU / RG(1, I, JGP)
          ELSE IF (INC.eq.-1) THEN
            ST = RG(2, I, JGP)
          END IF
!         121190 (ORIGINAL AUTHOR: P. ROUS)
!         take back sublayer reordering by SRTLAY for output via NORD here
          DO K = 1, LMMAX
            IF (INC.EQ.1) THEN
              TSTORE(1, INSORT + K, JGP) = 4.0 * PI * TS(IN + K)
            ELSE
              TSTORE(2, INSORT + K, JGP) = 4.0 * PI * TS(IN + K)
     +                                    * RG_PERP(JGP)
            ENDIF
            TS(IN + K) = TS(IN + K) * ST
          ENDDO
        ENDDO

!       Start loop over scattered beams
        DO JG = 1, N
!         Complete computation in MFOLT_SIMPLE
          CALL MFOLT_SIMPLE(JG, NA, AK2, AK3, CYLM, N, LMMAX, PQ,
     &                      NT, TS, LMN, RG, RG1, RG2, NLAY, JGP,
     &                      LXM, MFOLT_RT, INC, POSS)
          JGS = JG + NS
          AM = AMULT(JG)
!         Put final matrix elements in proper location
          IF (INC.EQ.-1) THEN
            R_BOT(JGS, JGPS) = MFOLT_RT(1) * AM
            T_BOT(JGS, JGPS) = MFOLT_RT(2) * AM
          ELSE
            R_TOP(JGS, JGPS) = MFOLT_RT(1) * AM
            T_TOP(JGS, JGPS) = MFOLT_RT(2) * AM
          ENDIF
        ENDDO

        IF (INC.EQ.1) THEN
!         Add unscattered plane wave to transmission from above
          T_TOP(JGPS, JGPS) = T_TOP(JGPS, JGPS) + RG_PERP(JGP)
!         And restart loop for beams moving from solid to vacuum
          INC = -1
          GO TO 285
        ELSE
!         Add unscattered plane wave to transmission from below
          T_BOT(JGPS, JGPS) = T_BOT(JGPS, JGPS) + RG_PERP(JGP)
        ENDIF
      ENDDO  ! END DO over incident beams (JGP)
      RETURN
      END


!-----------------------------------------------------------------------
!  SUBROUTINE BEAMS_YLM computes the spherical harmonics of the
!  propagation directions inside the crystal of the N beams
!  N_OFFSET+1..N_OFFSET+N in PQ and stores them in rows 1..N of CYLM.
!  Used by RSMF_SIMPLE and RTINV_SIMPLE for one beam set, and by the
!  main program of ref-calc for all NT beams before FINAL1.
!
!   CYLM= OUTPUT: CYLM(IG,K) IS THE K-TH SPHERICAL HARMONIC OF BEAM
!    N_OFFSET+IG; ROWS N+1..NT ARE LEFT UNCHANGED.
!   XA= OUTPUT: COMPLEX PERPENDICULAR WAVE-VECTOR COMPONENT OF BEAM
!    N_OFFSET+IG INSIDE THE CRYSTAL.
!   PQ= LIST OF BEAMS G.
!   NT= TOTAL NO. OF BEAMS IN MAIN PROGRAM AT CURRENT ENERGY.
!   N= NO. OF BEAMS TO COMPUTE.
!   N_OFFSET= OFFSET OF THE FIRST OF THEM IN LIST PQ.
!   E,VPI= CURRENT COMPLEX ENERGY.
!   AK2,AK3= PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!   LMAX= LARGEST VALUE OF L.
!   LMMAX= (LMAX+1)**2.
!   YLM= WORKING SPACE.
      SUBROUTINE BEAMS_YLM(CYLM, XA, PQ, NT, N, N_OFFSET, E, AK2, AK3,
     &                     VPI, LMAX, LMMAX, YLM)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NT, N, N_OFFSET, LMAX, LMMAX
      REAL, INTENT(IN) :: PQ(2, NT), E, AK2, AK3, VPI
      COMPLEX, INTENT(INOUT) :: CYLM(NT, LMMAX)
      COMPLEX, INTENT(OUT) :: XA(N), YLM(LMMAX)

      INTEGER IG, JG, K
      REAL BK2, BK3, C, B
      COMPLEX YA, CT, ST, CF

      YA = CMPLX(2.0 * E, - 2.0 * VPI + 0.000001)
      YA = SQRT(YA)
      DO IG = 1, N
        JG = IG + N_OFFSET
        BK2 = PQ(1, JG) + AK2
        BK3 = PQ(2, JG) + AK3
        C = BK2 * BK2 + BK3 * BK3
        XA(IG) = CMPLX(2.0 * E - C, - 2.0 * VPI + 0.000001)
        XA(IG) = SQRT(XA(IG))
        B = 0.0
        CF = CMPLX(1.0, 0.0)
        IF (C.GT.1.0E-7) THEN ! off-normal beam
          B = SQRT(C)
          CF = CMPLX(BK2/B, BK3/B)
        ENDIF
        CT = XA(IG)/YA
        ST = B/YA
        CALL SPHRM_MOD(LMAX, YLM, LMMAX, CT, ST, CF)
        DO K = 1, LMMAX
          CYLM(IG, K) = YLM(K)
        ENDDO
      ENDDO

      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE SPHRM_MOD COMPUTES SPHERICAL HARMONICS.
!  AUTHOR  PENDRY.
!
!  Edited: 2021-09-22 MRiva. Get rid of GOTOs and add ENDDOs
!
!
!   LMAX= LARGEST VALUE OF L.
!   YLM= OUTPUT COMPLEX SPHERICAL HARMONICS.
!   LMMAX= (LMAX+1)**2.
!   CT= COS(THETA) (COMPLEX).
!   ST= SIN(THETA) (COMPLEX).
!   CF= EXP(I*FI).
      SUBROUTINE  SPHRM_MOD(LMAX, YLM, LMMAX, CT, ST, CF)
      COMPLEX  YLM
      COMPLEX  CT, ST, CF, SF, SA
!  (21=LMAX+1; 441=(LMAX+1)**2   FOR LMAX=20 - if you need more (what for???),
!   please introduce variable dimensions, carrying the respective quantities through
!   all calls)

      DIMENSION  FAC1(21), FAC3(21), FAC2(441), YLM(LMMAX)

      IF (LMAX.gt.20) THEN
        write(6,*) "Dimensions for auxiliary arrays FAC1,FAC2,FAC3 in"
        write(6,*) "subroutine sphrm too low for LMAX.gt.20 . Remedy:"
        write(6,*) "increase dim. to"
        write(6,*) "FAC1(LMAX+1),FAC2((LMAX+1)**2),FAC3(LMAX+1)."
        STOP
      END IF

      PI = 3.14159265
      LM = 0
      CL = 0.0
      A = 1.0
      B = 1.0
      ASG = 1.0
      LL = LMAX + 1
      DO L = 1, LL
        FAC1(L) = ASG * SQRT((2.0 * CL + 1.0) * A/(4.0 * PI * B * B))
        FAC3(L) = SQRT(2.0 * CL)
        CM =  - CL
        LN = L + L - 1
CDIR$ SHORTLOOP
        DO M = 1, LN
          LO = LM + M
          FAC2(LO) = SQRT((CL + 1.0 + CM) * (CL + 1.0 - CM)
     +                    /((2.0 * CL + 3.0) * (2.0 * CL + 1.0)))
          CM = CM + 1.0
        ENDDO
        CL = CL + 1.0
        A = A * 2.0 * CL * (2.0 * CL - 1.0)/4.0
        B = B * CL
        ASG =  - ASG
        LM = LM + LN
      ENDDO

      LM = 1
      CL = 1.0
      ASG =  - 1.0
      SF = CF
      SA = CMPLX(1.0, 0.0)
      YLM(1) = CMPLX(FAC1(1), 0.0)
CDIR$ SHORTLOOP
      DO L = 1, LMAX
        LN = LM + L + L + 1
        YLM(LN) = FAC1(L + 1) * SA * SF * ST
        YLM(LM + 1) = ASG * FAC1(L + 1) * SA * ST/SF
        YLM(LN - 1) =  - FAC3(L + 1) * FAC1(L + 1) * SA * SF * CT/CF
        YLM(LM + 2) = ASG * FAC3(L + 1) * FAC1(L + 1) * SA * CT * CF/SF
        SA = ST * SA
        SF = SF * CF
        CL = CL + 1.0
        ASG =  - ASG
        LM = LN
      ENDDO

      LM = 1
      LL = LMAX - 1
      DO L = 1, LL
        LN = L + L - 1
        LM2 = LM + LN + 4
        LM3 = LM - LN
CDIR$ SHORTLOOP
        DO M = 1, LN
          LO = LM2 + M
          LP = LM3 + M
          LQ = LM + M + 1
          YLM(LO) =  - (FAC2(LP) * YLM(LP) - CT * YLM(LQ))/FAC2(LQ)
        ENDDO
        LM = LM + L + L + 1
      ENDDO
      RETURN
      END

!-----------------------------------------------------------------------------
!  SUBROUTINE SUBRAS_BLAS uses layer doubling to generate reflection and
!  transmission matrices for a stack of layers of sufficient thickness
!  to achieve convergence to the semi-infinite crystal limit.
!  The crystal is build from two 'layers' with identical reflection/
!  transmission matrices. The layers may be asymmetric, i.e., the
!  reflection/transmission matrices for beams propagating toward and
!  away from the solid may be different.
!
!  Modified: MRiva 2021-09-10. Get rid of M, get rid of most GOTOs.
!            Use optimized BLAS routine CGEMM for matrix multiplication
!            in TLRTA_BLAS.
!            TODO: The only useful output of this is RA, as it is used
!            only to stack the bulk. Figure out if there is a way to
!            do this exactly with a matrix inversion!
!
!  RA,TA,RB,TB
!    Input reflection and transmission matrices for reflection and
!    transmission by 'input layers' for incidence towards the solid
!    (RA, TA) and for incidence away from it (RB, TB).
!    These matrices will be overwritten and on output represent
!    reflection and transmission by a stack of layers for incidence
!    towards the solid (RA,TA) and for incidence away from it (RB,TB).
!
!   N= total NO. OF BEAMS at current energy.
!   S1,S2,S3,S4,S5,S6,PP,XS,INV_PIVOT= WORKING SPACES.
!   PQ= LIST OF RECIPROCAL LATTICE VECTORS G (BEAMS).
!   AS= INTERLAYER VECTOR (AS(1).GT.0.0).
!  IN COMMON BLOCKS
!   EPS= CONVERGENCE CRITERION FOR LAYER DOUBLING
!   LITER= LIMIT ON NUMBER OF DOUBLINGS.
!   E,VPI  COMPLEX CURRENT ENERGY.
!   AK2,AK3  PARALLEL COMPONENTS OF PRIMARY INCIDENT K-VECTOR.
!  NOTE  DO NOT IN GENERAL ASSIGN THE SAME STORAGE AREA FOR RA AND RB,
!  AND FOR TA AND TB, IN THE CALLING PROGRAM, EVEN THOUGH RA=RB, TA=TB
!  MAY HOLD FOR THE INDIVIDUAL LAYERS  THESE EQUALITIES WILL IN GENERAL
!  NO LONGER HOLD AFTER DOUBLING.
      SUBROUTINE SUBRAS_BLAS (RA, TA, RB, TB, N,
     &                        S1, S2, S3, S4, S5, S6,
     &                        PP, INV_PIVOT, PQ, AS)
      INTEGER  INV_PIVOT
      COMPLEX  RA, TA, RB, TB, PP, IU, XX
      COMPLEX  S1, S2, S3, S4, S5, S6
      COMPLEX  SQRT, EXP
      DIMENSION  RA(N,N), TA(N,N), RB(N,N), TB(N,N), PP(N,2)
      DIMENSION  S1(N,N), S2(N,N), S3(N,N), S4(N,N), S5(N,N), S6(N,N),
     &           PQ(2,N)
      DIMENSION  AS(3), INV_PIVOT(N)

      COMPLEX TEMP

! AMI note: BE CAREFULL, there are convergence criteria EPS and EPS1. A mixup between them caused trouble in v1.72 !
      COMMON  /MS/  LMAX, EPS, LITER
      COMMON  E, AK2, AK3, VPI

   10 FORMAT(4H X =,E15.4)
   20 FORMAT(24H NO CONV IN SUBREF AFTER,I3,2X,9HITER, X =,E15.4)
   30 FORMAT(I4,16H  ITER IN SUBREF)

      IU = CMPLX(0.0, 1.0)
      AK = 2.0 * E
!     Generate plane-wave propagators between layers
      DO I = 1, N
        BK2 = AK2 + PQ(1, I)
        BK3 = AK3 + PQ(2, I)
        XX = CMPLX(AK - BK2 * BK2 - BK3 * BK3, - 2.0 * VPI + 0.000001)
        XX = SQRT(XX) * AS(1)
        X = BK2 * AS(2) + BK3 * AS(3)
        PP(I, 1) = EXP(IU * (XX + X))
        PP(I, 2) = EXP(IU * (XX - X))
      ENDDO

      L = 0
      X = 0.0                                                           280193

!     START ITERATION OVER NUMBER OF DOUBLINGS. Layers may be asymmetric,
!     but should be identical. This means that they may have differences
!     in the transmission matrices from above and from below, but both
!     have the same reflection/transmission from above, and the same
!     from below.

  800 X0 = X
      IF (N.NE.1) THEN
        CALL TLRTA_BLAS (RA, TA, RB, TB, N, S1, S2, S3, S4, S5, S6, PP,
     &                   INV_PIVOT)
      ELSE
!       If only one beam, use direct formula rather than the
!       matrix inversion of TLRTA
!       XX now is the multiple scattering effect, from 1 / (1-x) = sum x**k:
!       propagate beam down [PP(1,1)], reflect there from above [RA(1,1)],
!       propagate up again [PP(1,2)], then reflect from below [RB(1,1)]
        XX = 1.0
     +       /(CMPLX(1.0, 0.0) - PP(1,1) * RA(1,1) * PP(1,2) * RB(1,1))
        TEMP = 1.0 + TA(1,1) * PP(1,1) * TB(1,1) * PP(1,2) * XX
        RA(1,1) = RA(1,1) * TEMP
        RB(1,1) = RB(1,1) * TEMP

        TA(1,1) = TA(1,1) * PP(1,1) * TA(1,1) * XX
        TB(1,1) = TB(1,1) * PP(1,2) * TB(1,1) * XX
      ENDIF

!     PERFORM ONE DOUBLING STEP
      X = 0.0
      DO J = 1, N
        TEMP = PP(J, 1)
        DO I = 1, N
          XX = RA(I, J) * PP(I, 2) * TEMP
          X = X + ABS(REAL(XX)) + ABS(AIMAG(XX))
        ENDDO
      ENDDO
      WRITE (6, 10)  X
      L = L + 1

!     Check number of iterations
      IF (L.GE.LITER) THEN
        WRITE (6, 20)  LITER, X
        write(6, *) "Make sure that interlayer spacing is ",
     +              "large enough for plane-wave space, else ",
     +              "use composite layer instead. Stopping program."
        STOP

!     Check convergence
      ELSEIF (ABS(X-X0)/X .LT. EPS) THEN
        WRITE (6, 30)  L
      ELSE
         GO TO 800
      ENDIF

      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE TAUY_SIMPLE PRODUCES THE QUANTITIES TAU*YLM(G) FOR ALL CHEMICAL
!  ELEMENTS AND, IF BEEBY-TYPE INVERSION IS CALLED FOR IN THE RSP
!  CALCULATION, TH*YLM(G) FOR ALL RELEVANT LAYERS. THESE ARE USED
!  EITHER FOR THE BEEBY-TYPE INVERSION IN MTINV OR AS INITIAL VALUES FOR
!  THE RSP ITERATION IN MPERTI.
!  FOR QUANTITIES NOT DESCRIBED BELOW SEE CALLING ROUTINES (MPERTI OR
!  MTINV).
!
!  Edits: MRiva 2021-09-06. Remove obsolete code, rewrite to remove
!         all GOTOs
!
!   TAUG,TAUGM= OUTPUT RESULTING VECTORS FOR K(G+) AND K(G-),RESP.
!   LTAUG= NTAU*LMMAX
!   LT=PERMUTATION OF (LM) ORDER, FROM LXGENT.
!   JGP= CURRENT INCIDENT BEAM.
      SUBROUTINE TAUY_SIMPLE(TAUG, TAUGM, LTAUG, TAU, LMT, LEV, CYLM,
     &                       NT, LMMAX, LT, NTAU, LOD, LEE, LOE, JGP)
      COMPLEX CZ, ST, SU, CF, CF1
      COMPLEX TAU(LMT, LEV), CYLM(NT, LMMAX)
      COMPLEX TAUG(LTAUG), TAUGM(LTAUG)
      DIMENSION LT(LMMAX)
      CZ = CMPLX(0.0, 0.0)

!
!  PERFORM MATRIX PRODUCT TAU*YLM(G+-) FOR EACH CHEMICAL ELEMENT
      DO I = 1, NTAU
        IS = (I - 1) * LMMAX
        DO JLM = 1, LEV
          ST = CZ
          JTAU = IS + JLM
          DO ILM = 1, LEV
            KLP = LT(ILM)
            IF (ILM.GT.LEE) THEN
              ST = ST - TAU(JTAU, ILM) * CYLM(JGP, KLP)
            ELSE
              ST = ST + TAU(JTAU, ILM) * CYLM(JGP, KLP)
            ENDIF
          ENDDO
          TAUG(JTAU) = ST
          TAUGM(JTAU) = ST
        ENDDO

        IS = IS + LEV
        DO JLM = 1, LOD
          ST = CZ
          SU = CZ
          JTAU = IS + JLM
          DO ILM = 1, LOD
            KLP = LT(ILM + LEV)
            CF = TAU(JTAU, ILM) * CYLM(JGP, KLP)
            IF (ILM.GT.LOE) THEN
              ST = ST - CF
              SU = SU + CF
            ELSE
              ST = ST + CF
              SU = SU - CF
            ENDIF
          ENDDO
          TAUG(JTAU) = ST
          TAUGM(JTAU) = SU
        ENDDO
      ENDDO
      END

!-----------------------------------------------------------------------
!  SUBROUTINE THMAT PRODUCES THE MATRIX THAT IS INVERTED IN THE BEEBY-
!  TYPE METHOD IN SUBROUTINE RTINV. FOR QUANTITIES NOT EXPLAINED
!  BELOW SEE RTINV.
!
!  Edited: 2021-10-07 MRiva add ENDDOs and indents. Remove LMNI2 (== LMNI)
!
!   TH= OUTPUT MATRIX, STILL TO BE INVERTED.
!   LMNI= NLAY*LMMAX.
      SUBROUTINE THMAT(TH, LMNI, GH, LMG, LMMAX, MGH, NLAY, TAU,
     1                 LMT, LEV, LPS)
      COMPLEX TH(LMNI, LMNI), GH(LMG, LMMAX), TAU(LMT, LEV)
      COMPLEX RU, CZ, ST, SU
      DIMENSION LPS(NLAY), MGH(NLAY, NLAY)

      LOD = LMMAX - LEV
      RU = (1.0, 0.0)
      CZ = (0.0, 0.0)

!     Start from identity matrix
      DO J = 1, LMNI
        DO I = 1, LMNI
          TH(I, J) = CZ
        ENDDO
        TH(J, J) = RU
      ENDDO

      IC = 0
      NLAY1 = NLAY - 1
      DO IN = 1, NLAY1
        IC = IC + 1
        ICL = (IC - 1) * LMMAX
        ID = IC
        IN1 = IN + 1
        DO INP = IN1, NLAY
          ID = ID + 1
          IDL = (ID - 1) * LMMAX
          M = (MGH(IN, INP) - 1) * LMMAX
          M1 = (MGH(INP, IN) - 1) * LMMAX
          LP = (LPS(IN) - 1) * LMMAX
          LP1 = (LPS(INP) - 1) * LMMAX

!         INSERT -TAU*GH FOR EVEN L+M (IN FIRST INDEX)
          DO I = 1, LEV
            DO J = 1, LMMAX
              ST = CZ
              SU = CZ
              DO K = 1,LEV
                ST = ST + TAU(LP + I, K) * GH(M + K, J)
                SU = SU + TAU(LP1 + I, K) * GH(M1 + K, J)
              ENDDO
              TH(ICL + I, IDL + J) = -ST
              TH(IDL + I, ICL + J) = -SU
            ENDDO
          ENDDO    ! I loop, over even L+M

          LP = LP + LEV
          LP1 = LP1 + LEV
          M = M + LEV
          M1 = M1 + LEV
          ICV = ICL + LEV
          IDV = IDL + LEV

!         INSERT -TAU*GH FOR ODD L+M (IN FIRST INDEX)
          DO I = 1, LOD
            DO J = 1, LMMAX
              ST = CZ
              SU = CZ
              DO K = 1, LOD
                ST = ST + TAU(LP + I, K) * GH(M + K, J)
                SU = SU + TAU(LP1 + I, K) * GH(M1 + K, J)
              ENDDO
              TH(ICV + I, IDL + J) = -ST
              TH(IDV + I, ICL + J) = -SU
            ENDDO
          ENDDO  ! I loop, over odd L+M
        ENDDO  ! INP loop
      ENDDO  ! IN loop
      RETURN
      END

!-----------------------------------------------------------------------
!  SUBROUTINE TLRTA_BLAS performs each layer doubling step for
!  subroutine SUBRAS_BLAS by matrix inversion. The two layers are
!  assumed to be identical, but can be asymmetric in the beam
!  propagation direction.
!
!  This modified version is optimized to use the LAPACK/BLAS routine
!  CGEMM for matrix multiplication, rather than doing explicit looping.
!
!   R1,T1,R2,T2
!    Input reflection and transmission matrices of the two layers.
!    Numbers 1 and 2 refer to incidence towards the solid and
!    away from it, respectively. They are overwritten and will
!    contain the corresponding matrices for the combination.
!
!   N= NO. OF BEAMS USED.
!   S1, S2, R1_P, R2_P, T1_P, T2_P, INV_PIVOT= WORKING SPACES.
!   PP= PLANE-WAVE PROPAGATORS, FROM SUBROUTINE SUBRAS_BLAS.

      SUBROUTINE  TLRTA_BLAS (R1, T1, R2, T2, N,
     &                        S1, S2, R1_P, R2_P, T1_P, T2_P,  ! storage only
     &                        PP, INV_PIVOT)
      COMPLEX  R1, T1, R2, T2, XX, YY, C_ZERO, PP
      COMPLEX  S1, S2, R1_P, R2_P, T1_P, T2_P
      COMPLEX  P_UP, P_DOWN, C_M_ONE, C_ONE
      DIMENSION  R1(N, N), T1(N, N), R2(N, N), T2(N, N)
      DIMENSION  R1_P(N, N), R2_P(N, N), T1_P(N, N), T2_P(N, N)
      DIMENSION  S1(N, N), S2(N, N), PP(N, 2)
      DIMENSION  INV_PIVOT(N)

      C_ZERO = CMPLX(0.0, 0.0)
      C_ONE = CMPLX(1.0, 0.0)
      C_M_ONE = CMPLX(-1.0, 0.0)

!     Store in R1_P, R2_P, T1_P, and T2_P modified versions of the
!     reflection/transmission matrices that include also propagation
!     between the layers for the incoming beam j.
!                                                        i ^
!                                          R2(i,j)         | T2(i,j)
!     --<<TOP>>--------------------------------------------------------
!                                          j ^  |        j ^
!                                            |  v i        |
!                  |  ^ i      |
!                j v  |      j v
!     --<<BOT>>--------------------------------------------------------
!                R1(i,j)       |  T1(i,j)
!                            i v
!                                                        i ^
!                                          R2_P(i,j)       |  T2_P(i,j)
!     --<<TOP>>--------------------------------------------------------
!                  |           |           j ^  |        j ^
!                  |           |             |  v i        |
!                  |  ^ i      |             |             |
!                j v  |      j v             |             |
!     --<<BOT>>--------------------------------------------------------
!                R1_P(i,j)     |  T1_P(i,j)
!                            i v
      DO J = 1, N
        P_UP = PP(J, 2)
        P_DOWN = PP(J, 1)
        DO I = 1, N
          R1_P(I, J) = R1(I, J) * P_DOWN
          R2_P(I, J) = R2(I, J) * P_UP
          T1_P(I, J) = T1(I, J) * P_DOWN
          T2_P(I, J) = T2(I, J) * P_UP
        ENDDO
      ENDDO

!     Prepare multiple-scattering matrix S1 for initial and final beams
!     propagating 'down' (i.e., toward the solid), and right below the
!     top layer. Use (Id-X)**-1 = sum_k X**k, with single-pass scattering
!     X(i,j) = sum_k R2_P(i,k) R1_P(k,j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S1 = Id - X
     &           R2_P, N, R1_P, N,
     &           C_ZERO, S1, N)
      DO I = 1, N
        S1(I, I) = S1(I, I) + C_ONE
      ENDDO

!     LU decomposition of matrix S1 (multiple scattering, beams moving down)
      CALL CGETRF(N, N, S1, N, INV_PIVOT, INFO)

!     Then multiply with transmission from above, such that T1(i, j) is:
!     beam j coming from above, transmitted through top and multiple-scattered
!     into i, which is right below the top layer and propagating down
      CALL CGETRS('N', N, N, S1, N, INV_PIVOT, T1, N, INFO)

!     Now calculate the total transmission from above into S1(i, j)
!     by taking S1(i, j) = sum_k [T1(k, j)     ! transmission through top of j into k, with multiple scattering; k right below top.
!                                 T1_P(i, k)]  ! transmission through bottom of k into i, including propagation down.
!     S1 will later be copied into T1.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T1_P, N, T1, N,
     &           C_ZERO, S1, N)

!     Go for the total reflection from above: (1) propagate-reflect
!     at the bottom layer [==multiply with R1_P(i,k)] the transmitted
!     and multiple-scattered beam k [==T1(k,j)]. Store this in S2(i,j).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, R1_P, N, T1, N,
     &           C_ZERO, S2, N)

!     Then (2) propagate-transmit at top [==multiply with T2_P(i,k)]
!     this beam k [==S2(k,j)]. Also, add to this the specular reflection
!     at the top layer R1(i,j), and store the whole thing in R1.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T2_P, N, S2, N,
     &           C_ONE, R1, N)

!     The original T1 is not needed anymore. Replace it with the total
!     transmission (temporarily in S1)
      DO J = 1, N
        DO I = 1, N
          T1(I, J) = S1(I, J)
        ENDDO
      ENDDO

!     Now repeat the whole thing for incidence from below.

!     Prepare multiple-scattering matrix S2 for initial and final beams
!     propagating 'up' (i.e., away from the solid), and right above the
!     bottom layer. Use (Id-X)**-1 = sum_k X**k, with single-pass scattering
!     X(i,j) = sum_k R1_P(i,k) R2_P(k,j)
      CALL CGEMM('N', 'N', N, N, N,
     &           C_M_ONE,  ! notice the C_M_ONE==-1.0 coefficient because of S2 = Id - X
     &           R1_P, N, R2_P, N,
     &           C_ZERO, S2, N)
      DO I = 1, N
        S2(I, I) = S2(I, I) + C_ONE
      ENDDO

!     LU decomposition of matrix S2 (multiple scattering, beams moving up)
      CALL CGETRF(N, N, S2, N, INV_PIVOT, INFO)

!     Then multiply with transmission from below, such that T2(i, j) is:
!     beam j coming from below, transmitted through bottom and multiple-
!     scattered into i, which is right above the bottom layer and
!     propagating up
      CALL CGETRS('N', N, N, S2, N, INV_PIVOT, T2, N, INFO)

!     Now calculate the total transmission from below into S2(i, j)
!     by taking S2(i, j) = sum_k [T2(k, j)     ! transmission through bottom of j into k, with multiple scattering; k right above bottom.
!                                 T2_P(i, k)]  ! transmission through top of k into i, including propagation up.
!     S2 will later be copied into T2.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T2_P, N, T2, N,
     &           C_ZERO, S2, N)

!     Go for the total reflection from below: (1) propagate-reflect
!     at the top layer [==multiply with R2_P(i,k)] the transmitted
!     and multiple-scattered beam k [==T2(k,j)]. Store this in S1(i,j).
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, R2_P, N, T2, N,
     &           C_ZERO, S1, N)

!     Then (2) propagate-transmit at bottom [==multiply with T1_P(i,k)]
!     this beam k [==S1(k,j)]. Also, add to this the specular reflection
!     at the top layer R2(i,j), and store the whole thing in R2.
      CALL CGEMM('N', 'N', N, N, N,
     &           C_ONE, T1_P, N, S1, N,
     &           C_ONE, R2, N)

!     Finally, store the total transmission in T2.
      DO J = 1, N
        DO I = 1, N
          T2(I, J) = S2(I, J)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE TSCATF INTERPOLATES TABULATED PHASE SHIFTS AND PRODUCES
C  THE ATOMIC T-MATRIX ELEMENTS (OUTPUT IN AF). THESE ARE ALSO
C  CORRECTED FOR THERMAL VIBRATIONS (OUTPUT IN CAF). AF AND CAF
C  ARE MEANT TO BE STORED IN ARRAY TMAT FOR LATER USE IN RSMF, RTINV.
C
C   IEL= CHEMICAL ELEMENT TO BE TREATED NOW, IDENTIFIED BY THE INPUT
C    SEQUENCE ORDER OF THE PHASE SHIFTS (IEL=1,2 OR 3).
C   L1= LMAX+1.
C   ES= LIST OF ENERGIES AT WHICH PHASE SHIFTS ARE TABULATED.
C   PHSS= TABULATED PHASE SHIFTS.
C   NPSI= NO. OF ENERGIES AT WHICH PHASE SHIFTS ARE GIVEN.
C   EB-V= CURRENT ENERGY (V CAN BE USED TO DESCRIBE LOCAL VARIATIONS
C    OF THE MUFFIN-TIN CONSTANT).
C   PPP= CLEBSCH-GORDON COEFFICIENTS FROM SUBROUTINE CPPP.
C   NN1= NN2+NN3-1.
C   NN2= NO. OF OUTPUT TEMPERATURE-CORRECTED PHASE SHIFTS DESIRED.
C   NN3= NO. OF INPUT PHASE SHIFTS.
C   DR0= FOURTH POWER OF RMS ZERO-TEMPERATURE VIBRATION AMPLITUDE.
C   DRPER= RMS VIBRATION AMPLITUDE PERPENDICULAR TO SURFACE.
C   DRPAR= RMS VIBRATION AMPLITUDE PARALLEL TO SURFACE.
C   T0= TEMPERATURE AT WHICH DRPER AND DRPAR HAVE BEEN COMPUTED.
C   T= CURRENT TEMPERATURE.
C   TSF0,TSF,AF,CAF  SEE ABOVE.
      SUBROUTINE TSCATF(IEL,L1,ES,PHSS,NPSI,EB,V,PPP,NN1,NN2,NN3,
     +DR0,DRPER,DRPAR,T0,T,AF,CAF,NEL,LMAX1,PHS,DEL,CTAB,SUM,BJ)

      DIMENSION PHSS(NPSI,NEL,LMAX1),PHS(LMAX1),ES(NPSI)
      DIMENSION PPP(NN1,NN2,NN3)
      COMPLEX CTAB,SUM
      DIMENSION CTAB(NN3),SUM(NN2)
      COMPLEX*16 BJ(NN1)
      COMPLEX CI,DEL(L1),CA,AF(L1),CAF(L1)

  700 FORMAT(42H TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS)
      CI=(0.0,1.0)
      E=EB-V
      IF (E.LT.ES(1)) GO TO 850
C  FIND SET OF PHASE SHIFTS APPROPRIATE TO DESIRED CHEMICAL ELEMENT
C  AND INTERPOLATE LINEARLY TO CURRENT ENERGY (OR EXTRAPOLATE TO ENERGIES
C  ABOVE THE RANGE GIVEN FOR THE PHASE SHIFTS)

  710 I = 1
  720 IF ((E-ES(I))*(E-ES(I+1)))  750, 750, 730
  730 I = I + 1
      IF (I-NPSI)  720, 740, 740
  740 I = I - 1
  750 FAC = (E - ES(I))/(ES(I + 1) - ES(I))
      DO 760 L = 1, L1
  760 PHS(L) = PHSS(I,IEL,L) + FAC * (PHSS(I + 1,IEL,L)
     +       - PHSS(I,IEL,L))
C  COMPUTE TEMPERATURE-INDEPENDENT T-MATRIX ELEMENTS
      DO 790 L = 1, L1
      A = PHS(L)
      AF(L) = A * CI
      AF(L) = EXP(AF(L))
      A = SIN(A)
      AF(L) = A * AF(L)
  790 CONTINUE
C  AVERAGE ANY ANISOTROPY OF RMS VIBRATION AMPLITUDES
      DR=SQRT((DRPER*DRPER+2.0*DRPAR*DRPAR)/3.0)
C  COMPUTE TEMPERATURE-DEPENDENT PHASE SHIFTS (DEL)
      CALL  PSTEMP (PPP,NN1,NN2,NN3,DR0,DR,T0,T,E,PHS,DEL,CTAB,SUM,BJ)
C  PRODUCE TEMPERATURE-DEPENDENT T-MATRIX ELEMENTS
      DO 840 L = 1, L1
      CA = DEL(L)
      CAF(L) = CA * CI
      CAF(L) = EXP(CAF(L))
      CA = SIN(CA)
      CAF(L) = CA * CAF(L)
  840 CONTINUE
      RETURN
  850 WRITE(6,700)
      STOP
      END

C-----------------------------------------------------------------------
C  SUBROUTINE ZGE PERFORMS GAUSSIAN ELIMINATION AS THE FIRST STEP IN THE
C  SOLUTION OF A SYSTEM OF LINEAR EQUATIONS. THIS IS USED TO MULTIPLY
C  THE INVERSE OF A MATRIX INTO A VECTOR, THE MULTIPLICATION BEING DONE
C  LATER BY SUBROUTINE ZSU.
C   A= INPUT (COMPLEX) MATRIX, WILL BE OVERWRITTEN BY A TRIANGULARIZED
C    MATRIX, TO BE TRANSMITTED TO SUBROUTINE ZSU.
C   INT= STORAGE FOR PERMUTATION OF MATRIX COLUMNS (AND ROWS), TO BE
C    TRANSMITTED TO SUBROUTINE ZSU.
C   NR= FIRST DIMENSION OF A (.GE.NC).
C   NC= ORDER OF A.
C   EMACH= MACHINE ACCURACY.
      SUBROUTINE  ZGE (A, INT, NR, NC, EMACH)
      COMPLEX  A, YR, DUM
      DIMENSION  A(NR,NC), INT(NC)
      N = NC
      DO 680 II = 2, N
      I = II - 1
      YR = A(I,I)
      IN = I
      DO 600 J = II, N
      IF (CABS(YR)-CABS(A(J,I)))  590, 600, 600
  590 YR = A(J,I)
      IN = J
  600 CONTINUE
      INT(I) = IN
      IF (IN-I)  610, 630, 610
  610 DO 620 J = I, N
      DUM = A(I,J)
      A(I,J) = A(IN,J)
  620 A(IN,J) = DUM
  630 IF (CABS(YR)-EMACH)  680, 680, 640
  640 DO 670 J = II, N
      IF (CABS(A(J,I))-EMACH)  670, 670, 650
  650 A(J,I) = A(J,I)/YR
CDIR$ IVDEP
      DO 660 K = II, N
  660 A(J,K) = A(J,K) - A(I,K) * A(J,I)
  670 CONTINUE
  680 CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
C  SUBROUTINE ZSU TERMINATES THE SOLUTION OF A SYSTEM OF LINEAR
C  EQUATIONS, INITIATED BY SUBROUTINE ZGE, BY BACK-SUBSTITUTING THE
C  CONSTANT VECTOR.
C   A= INPUT MATRIX, PREPARED BY SUBROUTINE ZGE.
C   INT= INPUT PERMUTATION FROM SUBROUTINE ZGE.
C   X= INPUT CONSTANT VECTOR AND OUTPUT RESULTING VECTOR.
C   NR= FIRST DIMENSION OF A (.GE.NC).
C   NC= ORDER OF A.
C   EMACH= MACHINE ACCURACY.
      SUBROUTINE  ZSU (A, INT, X, NR, NC, EMACH)
      COMPLEX  A, X, DUM
      DIMENSION  A(NR,NC), X(NC), INT(NC)
      N = NC
      DO 730 II = 2, N
      I = II - 1
      IF (INT(I)-I)  690, 700, 690
  690 IN = INT(I)
      DUM = X(IN)
      X(IN) = X(I)
      X(I) = DUM
CDIR$ IVDEP
  700 DO 720 J = II, N
      IF (CABS(A(J,I))-EMACH)  720, 720, 710
  710 X(J) = X(J) - A(J,I) * X(I)
  720 CONTINUE
  730 CONTINUE
      DO 780 II = 1, N
      I = N - II + 1
      IJ = I + 1
      IF (I-N)  740, 760, 740
CDIR$ IVDEP
  740 DO 750 J = IJ, N
  750 X(I) = X(I) - A(I,J) * X(J)
  760 IF (CABS(A(I,I))-EMACH*1.0E-5)  770, 780, 780
  770 A(I,I) = EMACH * 1.0E - 5 * (1.0,1.0)
  780 X(I) = X(I)/A(I,I)
      RETURN
      END

C---------------------------------------------------------------------
      SUBROUTINE ZTU(A,INT,X,ND,N,EMACH)
C     ==============
      COMPLEX A,X,Z1,Z2
      DIMENSION A(ND,ND),INT(ND),X(ND)
C
COMMENT ZTU IS A MODIFIED BACK SUBSTITUTION SUBROUTINE
C       USING THE OUTPUT OF ZGE TO CALCULATE X TIMES
C       A INVERSE, RETURNED IN X
      X(1)=X(1)/A(1,1)
      DO 32 I=2,N
      I1=I-1
      Z1=X(I)
      DO 31 J=1,I1
   31 Z1=Z1-X(J)*A(J,I)
   32 X(I)=Z1/A(I,I)
      DO 34 NI=2,N
      I=N-NI+1
      I1=I+1
      Z1=X(I)
      DO 33 J=I1,N
   33 Z1=Z1-X(J)*A(J,I)
      IN=INT(I)
      IF(I.NE.IN) X(I)=X(IN)
   34 X(IN)=Z1
      RETURN
      END
