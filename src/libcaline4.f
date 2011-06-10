      FUNCTION INTGR(ARG)
      REAL*4 ARG,INTGR
      INTGR=INT(ARG+SIGN(0.5,ARG))
      RETURN
      END

      SUBROUTINE CALINE4(
     +                   NR,      ! number of receptors 
     +                   XR,      ! x-coordinates of receptors
     +                   YR,      ! y-coordinates of receptors
     +                   ZR,      ! z-coordinates (heights) of receptors
     +                   NL,      ! number of links 
     +                   XL1,     ! x-coordinates of link vertices
     +                   YL1,     ! y-coordinates of link vertices
     +                   XL2,     ! x-coordinates of link vertices
     +                   YL2,     ! y-coordinates of link vertices
     +                   WL,      ! link widths
     +                   HL,      ! link heights
     +                   TYP,     ! link classifications *as integers*
     +                   VPHL,    ! traffic volume, per link
     +                   EFL,     ! emission factor, per link
     +                   U,       ! wind speed
     +                   BRG,     ! wind bearing
     +                   CLAS,    ! atmospheric stability class
     +                   MIXH,    ! mixing height
     +                   SIGTH,   ! sigma theta
     +                   TEMP,    ! temperature
     +                   Z0,      ! surface roughness
     +                   PTYP,    ! POLLUTANT ID
     +                   MOWT,    ! MOLECULAR WEIGHT
     +                   VS,      ! settling velocity
     +                   VD,      ! deposition velocity 
     +                   C)       ! resulting concentrations (per receptor)

C      IMPLICIT NONE
      INTEGER NR,NL
      REAL*4 U,BRG,SIGTH,TEMP,Z0,VS,VD

C     THESE APPEARED TO BE MISSING ("NO IMPLICIT TYPE")
C      REAL*4 AA,AD,ALT,ARG,ARG2,ASED,BASE,BB,BOT,BRG1,BSTO,
C     *  DPHI,EXP1,EXP2,ULIM2,
      
      INTEGER BLINE,CASE,CC,CLAS,CLINE,CLOSE,CNTR,COWNT,
     *  EFLCOD,METCOD,NCYC(20),NDLA(20),NUMLN(6),OH(6),OHS(6,8),
     *  PGCT,PTYP,RC,RR,RTYP,RUNOUT,SKIP,SPLIT,STOC(24),
     *  STOV(24,20),TYP(NL),VPHCOD,X1TOX2
      REAL*4 IDTI(20),IDT2(20),INC,KF,KR,KZ,LACC(20),LB,LBRG(20),
     *  LDCL(20),LIM1,LIM3,LLIM2,LQU(20),MIXH,MIXWL(20),MIXWLO(20),
     *  MIXWR(20),MIXWRO(20),MOWT,NE,
     *  NO,NOA,NO2,NO2A,NO2T
      DOUBLE PRECISION A,APRI,ARM1,ARM2,AXD1,AXD2,AYD1,AYD2,B,BPRI,
     *  D,DPRI,D1,D2,FAC2,HYP,INTG(7),L,LL(20),LPRI,
     *  PD,SIDE,XD,XPRI,YD,YPRI
      LOGICAL CANYON,COMPT,CYNBLF,CYNDUN,LEFTRF,NCTRIB,NOX,RTSIDE
C
C *****  DIMENSION ARRAYS  *****
C
      REAL*4 ACCR(20),ACCT(20),AZ(7),BRHOLD(20),ABRG1(20),
     *   ABRG2(20),C(NR),CADD(20),CHOLD(20),CLNK(20,20),
     *   CTOT(20),DCLR(20),DCLT(20),EFA(20),EFC(20),EFD(20),EFI(20),
     *   EFIO(20),EFL(NL),HL(NL),HLO(20),MODLNK(20),PRCLK(20),SCALU(2),
     *   SKR(24),SNO(24),SNO2(24),SO3(24),SPD(20),STOA(24),STOB(24),
     *   STOE(24,20),STOM(24),STOS(24),STOT(24),STOU(24),STPL(20),
     *   STPLO(20),VPHL(NL),WL(NL),WT(6),
     *   XL1(NL),XL2(NL),XR(NR),Y(8),YL1(NL),YL2(NL),YR(NR),ZR(NR)      

      CHARACTER COD(20),CRL(4),CRM(20,20),DASH,INLINE(80),LINE(80),
     *          SPACE,SPACEI(70),STB(7),TITLE1(9),TITLE2(5),TLINK(4)
      CHARACTER*2 TEYP(6),TP(20)
      CHARACTER*4 UNIT,UNTS(2)
      
C      COMMON PGCT,RR,RTYP,PTYP
C--      COMMON /DP/ JOB,NAME,RUN
C      COMMON /ENV/ INLINE,Z0,MOWT,VS,VD,NR,NL,SCAL,LC,RC,UALT
C      COMMON /M/ STPL,LQU,LACC,LDCL,NCYC,NDLA,VSP
C      COMMON /MT/ DCLT,ACCT,IDT1,IDT2,ACCR,DCLR,SPD
C      COMMON /ME/ EFA,EFC,EFD,EFI,VPHI,VPHO
      
      DATA AZ/1112.,566.,353.,219.,124.,56.,22./
C
C          STABILITY CLASS MODIFICATION MATRIX
C
      DATA COD/'A','B','C','D','E','F','G','H','I','J',
     *         'K','L','M','N','O','P','Q','R','S','T'/
      DATA CRL/'.','X','O',' '/
      DATA OHS /9,7,8,7,8,0,
     *          9,7,8,7,8,8,
     *          9,7,8,8,0,0,
     *          6,7,7,7,0,7,
     *          6,7,7,7,7,7,
     *          9,7,8,7,8,0,
     *          9,7,8,7,8,8,
     *          9,7,8,8,0,0/
      DATA SCALU /1.0000,0.3048/
      DATA STB/'A','B','C','D','E','F','G'/
      DATA TEYP /'AG','DP','FL','BG','PK','IN' /
      DATA UNTS /'(M) ','(FT)' /

C
C ***** INITIALIZATION OF CONSTANTS AND COUNTERS *****
C
      PGCT=1
      MLNUM=0
      PI=3.1415926
      RAD=PI/180.
      DEG=180./PI
      DREF=ALOG(10000.)
C     ...
      NOX=.FALSE.
      ICF1=0
      ICF2=0
      
      
      IF (PTYP.GT.1) GOTO 5
      MOWT=28.
      VS=0.
      GOTO 20
      
 5    IF (PTYP.GT.2) GOTO 10
      MOWT=46.
      VS=0.
      NOX=.TRUE.
      GOTO 20
      
 10   IF (PTYP.GT.3) GOTO 15
      VS=0.
      GOTO 20
      
 15   IF (PTYP.EQ.4) GOTO 20
      STOP                    ! PTYP OUT OF RANGE

 20   V1=VD-VS/2.
 
      SCAL=1.
      ALT=SCAL*UALT
      
C     TYP = HIGHWAY TYPE
C           1=AG: AT-GRADE
C           2=DP: DEPRESSED (CUT)
C           3=FL: FILL
C           4=BR: BRIDGE
C           5=PK: PARKING LOT
C           6:IN: INTERSECTION (MODAL LINK)

  995 RR=1
 1000 DO 1005 J=1,20
      CADD(J)=0.
      CTOT(J)=0
 1005 CONTINUE
      COWNT=1   
      RUNOUT=0
      NOPRINT=0
C 
 1010 CONTINUE     
C-- 1010 READ (5,1920,END=9990) RTYP,VPHCOD,EFLCOD,INTCOD,METCOD,RUN(RR)
C
C          RTYP = RUN TYPE
C          OPTIONS:  1 = STANDARD RUN
C                    2 = MULTI-RUN
C                    3 = WORST-CASE WIND ANGLE
C                    4 = MULTI-RUN/WORST-CASE
C                    9 = LAST RUN OF MULTI-RUN
C        VPHCOD = TRAFFIC VOL. CODE          0=USE PREVIOUS RUN DATA
C        EFLCOD = EMISSION FACTOR CODE                 "
C        INTCOD = INTERSECTION PARAMETER CODE
C        METCOD = METEOROLOGY CODE
C
C           RUN = RUN TITLE (12 CHAR. MAX.)
      IF (RTYP.LT.9) GOTO 1020
      RUNOUT=RTYP
      RTYP=SRTYP
      GOTO 1030
 1020 SRTYP=RTYP
 1030 IF (VPHCOD.EQ.0) GOTO 1040
      READ (5,*) (VPHL(ID),ID=1,NL)
 1040 IF (EFLCOD.EQ.0) GOTO 1050
      READ (5,*) (EFL(ID),ID=1,NL)
 1050 IF(INTCOD.EQ.0) GOTO 1070
C
      STOP                    ! UNIMPLEMENTED
C
 1070 DO 1080 J=1,NL
      STOV(COWNT,J)=VPHL(J)
      STOE(COWNT,J)=EFL(J)
 1080 CONTINUE
      IF (METCOD.EQ.0) GOTO 1105
      IF (PTYP.NE.2) GOTO 1090            ! IF NOT NOX THEN SKIP AHEAD
      READ (5,*) BRG,U,CLAS,MIXH,SIGTH,TEMP,O3,NOA,NO2A,KR
      KF=51.7*(EXP(-1450./(TEMP+273.)))   ! KF IN (1/(PPM*SEC))
      AMB=NO2A
      GOTO 1100
 1090 READ (5,*) BRG,U,CLAS,MIXH,SIGTH,AMB,TEMP
C
C        U = WIND SPEED (M/S)
C      BRG = WIND DIRECTION (DEGREES)
C    SIGTH = STANDARD DEV. OF WIND ANGLE (DEGREES)
C     CLAS = STABILITY CLASS (A-G)
C     MIXH = MIXING HEIGHT (M)
C      AMB = AMBIENT CONCENTRATION (PPM)
C     VPHL = TRAFFIC VOLUME (VEH/HR)
C     TEMP = TEMPERATURE (C)
C
C    *** NOX INPUTS ***
C
C       O3 = AMBIENT OZONE (PPM)
C      NOA = AMBIENT NO (PPM)
C     NO2A = AMBIENT NO2 (PPM)
C       KR = REVERSE REACTION RATE (1/SEC)
 1100 FPPM=(0.02241*((273.+TEMP)/273./
     *      EXP(-.03417*ALT/(TEMP+273.))))/MOWT
      IF (PTYP.EQ.4) FPPM=1.                     ! PARTICULATE OUTPUT IN MG/M3
      BSTO=BRG
      GOTO 1495
 1105 BRG=BSTO
C
C      ! USE 1ST ELEMENT TO HOLD INFO IF NOT WORST CASE
C
 1495 IR=1
      IF (RTYP.NE.3 .AND. RTYP.NE.4) GOTO 1700
C
C **** SELECTION OF RANGE FOR WORST CASE WIND SEARCH ****
C
      STOP                    ! UNIMPLEMENTED
C
 1700 SIGT=RAD*(SIGTH)
      BRG1=BRG                ! WIND ANGLE FOR OUTPUT
      BRG=BRG+180.
      IF (BRG.GE.360.) BRG=BRG-360.
      
      IF (RTYP.NE.3 .AND. RTYP.NE.4) GOTO 1710
C--      DO 1705 I=1,NL
C--      C(IR)=0.
C-- 1705 CONTINUE
      GOTO 1720
C
C-- 1710 DO 1715 I=1,NL
C--          DO 1715 J=1,NR
C--      C(I,J)=0.
C-- 1715 CONTINUE
C
C *****  LINK LOOP *****
C
 1710 CONTINUE
 1720 IM=1
      DO 6000 IL=1,NL
      VPH=VPHL(IL)
      EF=EFL(IL)
      IF (TYP(IL).NE.6) GOTO 1725
      ML=MODLNK(IM)
      IM=IM+1
 1725 IF (TYP(IL).EQ.2 .OR. TYP(IL).EQ.3) GOTO 1730
      H=HL(IL)
      GOTO 1735
 1730 H=0
 1735 W=WL(IL)
      CYNBLF=.FALSE.
      CANYON=.FALSE.
      IF (MIXWR(IL).NE.0 .OR. MIXWL(IL).NE.0) CYNBLF=.TRUE.
      IF (CYNBLF) ICF1=1
      IF (MIXWR(IL).NE.0 .AND. MIXWL(IL).NE.0) CANYON=.TRUE.
      IF (.NOT.CYNBLF .OR. (RTYP.NE.3 .AND. RTYP.NE.4)) GOTO 1800
C
C *****  LINK ROUTINE *****
C *************************
C
 1800 W2=W/2
      Q1=0.1726*VPH*EF
      XD=XL2(IL)-XL1(IL)
      YD=YL2(IL)-YL1(IL)
      IF (LL(IL).LT.DABS(XD)) LL(IL)=DABS(XD)
      LB=DEG*(DACOS(DABS(XD))/LL(IL))               ! LINK BEARING
      IF (XD.GT.0. .AND. YD.GE.0.) LB=90.-LB
      IF (XD.GT.0. .AND. YD.LT.0.) LB=90.+LB
      IF (XD.LT.0. .AND. YD.LE.0.) LB=270.-LB
      IF (XD.LT.0. .AND. YD.GT.0.) LB=270.+LB
      PHI=ABS(BRG-LB)                           ! WIND ANGLE W.R.T. LINK
      IF (CYNBLF) PHI=0.
      IF (PHI.LE.90) GOTO 1820
      IF (PHI.GE.270) GOTO 1810
      PHI=ABS(PHI-180.)
      GOTO 1820
 1810 PHI=ABS(PHI-360.)
 1820 BASE=1.1+PHI**3/2.5E5                     ! SET ELEMENT GROWTH BASE
      PHI=RAD*(PHI)                             ! CONVERSION OF PHI TO RADIANS
      IF (PHI.GT.1.5706) PHI=1.5706
      IF (PHI.LT.0.00017) PHI=0.00017
      DPHI=DEG*PHI
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      TPHI=TAN(PHI)
      IF (DPHI.GT.45) WMIX=W2/SPHI
      IF (DPHI.LE.45) WMIX=W2/SIN(RAD*45)
      DMIX=W2/SPHI
      IF (DMIX.GT.LL(IL)) DMIX=LL(IL)           ! CHECK FOR MAXIMUM DMIX
      SYCRIT=W2/0.6744
      TI=300.
      SYFAC=0.9*SYCRIT/SQRT(U*TI)
      DMIXM=((SYFAC+SQRT(SYFAC**2.+4.*SIGT*SYCRIT))/(2*SIGT))**2.
      IF ((DMIXM/U).LE.550) GOTO 1830
      CSA=-SYCRIT/SIGT
      CSB=-28.*SYCRIT*SQRT(U)/SIGT
      CSAB=CSB**2/4.+CSA**3/27.
      IF (CSAB.LT.0) GOTO 1830
      CSAB=SQRT(CSAB)
      CLA=(CSAB-CSB/2.)**(1./3.)
      CLB=-CSAB-CSB/2.
      IF (CLB.GE.0) CLB=CLB**(1./3.)
      IF (CLB.LT.0) CLB=-(-CLB)**(1./3.)
      DMIXM=(CLA+CLB)**2.
 1830 IF (DMIX.GT.DMIXM) DMIX=DMIXM
      IF (DMIX.LT.WMIX) DMIX=WMIX
      
      IF (COMPT) GOTO 1850                      ! WHAT IS COMPT?
      IF (RTYP.NE.2 .OR. RUNOUT.GE.6) GOTO 1000
      NOPRNT=1
      GOTO 1010
C
C *****  DEPRESSED SECTION *****
C
 1850 IF (HL(IL).LT.-1.5) GOTO 1855
      DSTR=1.
      HDS=1.
      GOTO 1860
 1855 HDS=HL(IL)
      DSTR=0.72*ABS(HDS)**0.83      ! RESIDENCE TIME FACTOR
 1860 TR=DSTR*WMIX/U                ! RESIDENCE TIME
      SGZI=1.7+0.1*TR 
      IF(TYP(IL).EQ.5) SGZI=1.      
      RFAC=(Z0/10.)**0.07           ! ROUGHNESS FACTOR
      IF (.NOT.CANYON) GOTO 1865
C     CALL MSGZ(...)      
      SGZF=SGZF*RFAC                ! SIGMA Z AT 10 KM, CANYON
      GOTO 1870
 1865 SGZF=AZ(CLAS)*RFAC            ! SIGMA Z AT 10 KM, PASSIVE RELEASE
 1870 SGZB=SGZB*RFAC                ! BASELINE SIGMA Z AT 10 KM    
      PZ2=(ALOG(SGZB)-ALOG(SGZI))/(DREF-ALOG(WMIX))
      PZ1=EXP((ALOG(SGZB)+ALOG(SGZI)-PZ2*(DREF+ALOG(WMIX)))/2.)
      PZ3=0.
      IF (DMIX.GE.10000. .OR. CLAS.EQ.1) GOTO 2995
      PZ3=ALOG((SGZF/SGZB)**(1./ALOG(10000./DMIX)))/ALOG(10000./DMIX)
      IF (-PZ2/(2.*PZ3).GT.ALOG(10000./DMIX)) GOTO 1880
      PZ3=-PZ2/(2.*ALOG(10000./DMIX))
C
C  **** INITIAL NOX CALCUALTIONS ****
C
 1880 IF (PTYP.NE.2) GOTO 2995
      NO=NOA+(0.925*Q1/(3.5*U))*FPPM*(46./30.)
      NO2=NOA+(0.075*Q1/(3.5*U))*FPPM
      AA=KF*03*NO-KR*NO2
      BB=-(KF*O3+KF*NO+KR)
      CCC=KF
      PP=SQRT(BB**2-4*AA*CCC)
C
C *****  END OF LINK ROUTINE  *****
C
C *****  RECEPTOR LOOP  *****
C
 2995 IF (RTYP.EQ.3 .OR.RTYP.EQ.4) GOTO 3010
      IR=0
 3000 IF (RTYP.EQ.3 .OR.RTYP.EQ.4) GOTO 3010
      IR=IR+1
      IF (IR.GT.NR) GOTO 5000
 3010 RDX1=XR(IR)-XL1(IL)
      RDY1=YR(IR)-YL1(IL)
      A=RDX1**2+RDY1**2
      B=(XR(IR)-XL2(IL))**2
      B=B+(YR(IR)-YL2(IL))**2
      L=(B-A-LL(IL)**2)/(2.*LL(IL))       ! OFFSET LENGTH
      IF (A.GT.L**2) D=DSQRT(A-L**2)      ! RECEPTOR DISTANCE
      IF (A.LE.L**2) D=0.
      IF (D.LT.W2) DFAC=W2
      IF (D.GT.W2) DFAC=D
      UWL=LL(IL)+L                        ! UPWIND LENGTH
      DWL=L                               ! DOWNWIND LENGTH
C--      CRM(IL,R)=CRL(4)
      IF (.NOT.CYNBLF) GOTO 3045
      TXVEC=XVEC
      TYVEC=YVEC
      IRLOC=0
      RLR=(XD*RDY1)-(YD*RDX1)
C      IF RLR < 0 RECEIVER ON RT. RLC > 0 RECEIVER ON LT.
      IF(RLR.GE.0) GOTO 3012
C      RECEIVER ON RIGHT. IN OR OUT OF CANYON?
      IF(MIXWR(IL).GT.0.AND.MIXWR(IL).LT.D) IRLOC=1
      GOTO 3014
C      RECEIVER ON LEFT. IN OR OUT OF CANYON?
 3012 IF(MIXWL(IL).LT.0.AND.MIXWL(IL)*(-1).LT.D) IRLOC=1
C      IF L<0 DOWN LINK. L+LL<0 UP LINK.
 3014 IF(L.GT.0.) IRLOC=IRLOC-2
      IF(L.LT.0.AND.LL(IL)+L.LT.0.) IRLOC=IRLOC+2
C      SET CAN/BLUFF MATRIX, PRINTFLAG, AND BYPASS CALC IF OUTSIDE
C--      CRM(IL,IR)=CRL(3)
      IF(IRLOC.NE.0) ICF2=1
      IF(IRLOC.EQ.1.OR.IRLOC.EQ.(-1).OR.IRLOC.EQ.3) GOTO 4000
C--      IF(IRLOC.EQ.2.OR.IRLOC.EQ.(-2)) CRM(IL,IR)=CRL(2)
C--      IF(IRLOC.EQ.0) CRM(IL,IR)=CRL(1)
C
C   CHECK USERS BRG FOR CANYON LINK
C
      DO 3040 IB=1,2
         UBRG=BRG1
         CB=LB-180.+180.*IB
         IF (CB.GT.360.) CB=CB-360.
         IF (CB.GT.1.) GOTO 3020 
         CB1=CB+361.
         CB2=CB+359.
         IF (UBRG.LT.359) UBRG=UBRG+360.
         GOTO 3030
 3020    IF (CB.LT.359.) GOTO 3025
         IF (UBRG.LT.358.) UBRG=UBRG+360
 3025    CB1=CB+1.
         CB2=CB-1.
 3030    IF (UBRG.GT.CB1 .OR. UBRG.LT.CB2) GOTO 3035
         BRG1=CB
         BRG=CB+180.
         IF (BRG.GE.360.) BRG=BRG-360.
         GOTO 3045
 3035    IF (IB.EQ.1) GOTO 3040
         STOP
 3040 CONTINUE
 3045 XVEC=COS(RAD*(450-BRG+0.01))        ! VIRTUAL DISPLACEMENT, X
      YVEC=SIN(RAD*(450.-BRG+0.01))       ! VIRTUAL DISPLACEMENT, Y
      IF (D.EQ.0) XPRI=XR(IR)+XVEC
      IF (D.NE.0) GOTO 3050
      YPRI=YR(IR)+YVEC
      GOTO 3055
 3050 XPRI=XR(IR)+D*XVEC
      YPRI=YR(IR)+D*YVEC
 3055 APRI=(XPRI-XL1(IL))**2+(YPRI-YL1(IL))**2
      BPRI=(XPRI-XL2(IL))**2+(YPRI-YL2(IL))**2
      LPRI=(BPRI-APRI-LL(IL)**2)/(2.*LL(IL))
      IF (APRI.GT.LPRI**2) DPRI=DSQRT(APRI-LPRI**2)
      IF (APRI.LE.LPRI**2) DPRI=0.
      IF (DPRI.LT.D) D=-D
      X1TOX2=0
      IF(LPRI-L) 3060,3065,3065
 3060 TMP=UWL
      UWL=-DWL
      DWL=-TMP
 3065 IF (TYP(IL).NE.2 .AND. 
     *    TYP(IL).NE.3) GOTO 3070
      D1=W2+2.*ABS(HL(IL))
      D2=W2
      IF(DABS(D).GE.D1) GOTO 3070
      IF(DABS(D).LE.D2) Z=ZR(IR)-HL(IL)
      IF(DABS(D).GT.D2) 
     *   Z=XR(IR)-HL(IL)*(1.-(DABS(D)-W2)/(2.*ABS(HL(IL))))
      GOTO 3080
 3070 Z=ZR(IR)
C
C *****  CALINE4 ROUTINE  *****
C *****************************
C
 3080 SGN=1.            ! SGN := +1 FOR UPWIND, -1 FOR DOWNWIND
      MFLG=1
      IF (TYP(IL).EQ.6) MFLG=0
      NE=0
      FLAGU=0.
      FLAGD=0.
      IF(MFLG.NE.0. .AND. UWL.LE.0. .AND. DWL.LT.0.) SGN=-1
C
C *****  ELEMENT LOOP *****
C
 3100 IF (DPHI.GT.45) ED1=DFAC*(1./TPHI)-W2
      IF (DPHI.LE.45) ED1=DFAC-W2
      SW=1.
      IF (X1TOX2.EQ.1) SW=-1
      IF (TYP(IL).EQ.6) ED1=(STPL(IL)+L)*SW
      IF (ED.GE.UWL) SGN=-1.
      IF(SQN.EQ.-1.AND.ED1.LE.DWL) GOTO 4000
      ED2=ED1+SGN*W
C      ! INITIALIZATION OF ELEMENT LIMITS
 3110 IF (SGN.EQ.-1.) GOTO 3160
      IF (ED1.LE.DWL .AND. ED2.LE.DWL) GOTO 3770
      IF (ED1.GT.DWL .AND. ED2.LT.UWL) GOTO 3250
      IF (ED1.LE.DWL) ED1=DWL
      IF (ED2.LT.UWL) GOTO 3250
      ED2=UWL
      FLAGU=1.
      GOTO 3250
 3160 IF (ED1.GE.DWL .AND. ED2.GE.DWL) GOTO 3770
      IF (ED1.LT.DWL .AND. ED2.GT.UWL) GOTO 3250
      IF (ED1.GE.UWL) ED1=UWL
      IF (ED2.GT.DWL) GOTO 3250
      ED2=DWL
      FLAGD=1.
 3250 EL2=ABS(ED2-ED1)/2.                       ! ELEMENT HALF-DISTANCE
      IF (EL2.EQ.0) GOTO 3762
      ECLD=(ED1+ED2)/2.                         ! ELEMENT CENTERLINE DISTANCE
      ELL2=W2/CPHI+(EL2-W2*TPHI)*SPHI           ! EQUIVALENT LINE HALF-LENGTH
      IF (PHI.GE.ATAN(W2/EL2)) CSL2=W2/SPHI     ! CENTRAL ZONE HALF-LENGTH
      IF (PHI.LT.ATAN(W2/EL2)) CSL2=EL2/CPHI    
      EM2=ABS((EL2-W2/TPHI)*SPHI)
      EN2=ELL2-EM2
C
C *****  RECEPTOR DISTANCE LOOP *****
C
      IF (TYP(IL).NE.6) GOTO 3335
      STOP  ! NOT IMPLEMENTED
C      IF (X1TOX2.EQ.1) GOTO 3305
C      ZD1=-L+ED1
C      ZD2=-L+ED2
C      GOTO 3315
C 3305 ZD1=-L-ED1
C      ZD2=-L-ED2
C 3315 CALL MODAL(ZD1,ZD2,IL,Q1)                 ! TODO: MODAL(...) (???)
 3335 QE=Q1*CSL2/W2
      FET=(ECLD+D*TPHI)*CPHI                    ! ELEMENT FETCH
      IF (FET.GT.10000 .AND. SGN.EQ.1) GOTO 3820  ! ELEMENT DOESNT CONTRIBUTE
      HYP=ECLD**2+D**2
      SIDE=FET**2
      IF (SIDE.GT.HYP) YE=0.                    ! DIST FROM CENTERLINE TO RCP
      IF (SIDE.LE.HYP) YE=DSQRT(HYP-SIDE)
      IF (ECLD.EQ.0) GOTO 3345
      IF (ECLD.GT.0 .AND. PHI.GT.DATAN(D/ECLD)) YE=-YE
      IF (ECLD.LT.0 .AND. PHI.LT.DATAN(D/ECLD)) YE=-YE
C
C *****  DETERMINE SIGMA Y AND SIGMA Z *****
C
 3345 IF (FET.GT.-CSL2) GOTO 3350               ! ELEMENT DOES NOT CONTRIBUTE
      IF (SGN.EQ.-1) GOTO 4000
      IF (SGN.NE.-1) GOTO 3768
 3350 IF (FET.GE.CSL2) GOTO 3355                ! RCP NOT COMPLETELY DOWNWIND
      FET=(CSL2+FET)/2.
 3355 SGZ=SGZI                                  ! CONSTANT OVER MIXING ZONE
      IF (FET.GT.WMIX) SGZ=PZ1*FET**PZ2         ! BASELINE CURVE
      IF (FET.GT.DMIX) SGZ=SGZ*((FET/DMIX)**PZ3)**(ALOG(FET/DMIX))
      KZ=SGZ**2.*U/(2.*FET)                     ! VERTICAL DIFFUSIVITY ESTIMATE
      TT=FET/U
      TI=300.
      IF (TT.GT.550.) TI=.001*(TT**2.)
      F1=1./(1.+.9*(TT/TI)**.5)
      SGY=FET*SIGT*F1                           ! SIGMA Y
      FAC1=0.399/(SGZ*U)                        ! SOURCE STRENGTH - WIND SPEED FACTOR
      IF (PTYP.NE.2) GOTO 3390
C
C  **** NOX CALCULATIONS ****
C
      TTPP=TT*PP
      IF (TTPP.GT.88) GOTO 3360
      TOP=2*AA*(EXP(TTPP)-1)
      BOT=BB*(1-EXP(TTPP))+PP*(1+EXP(TTPP))
      XCON=TOP/BOT
      GOTO 3365
 3360 XCON=(-BB-PP)/(2*CCC)
 3365 NO2T=NO2+XCON
      IF ((NO2T-NO2A).LE.0) QE=0
      IF ((NO2T-NO2A).GT.0) QE=((NO2T-NO2A)/(FPPM/(3.5*U)))*(CSL2/W2)

 3390 FAC2=0.
      TYE=YE
      RTSIDE=.TRUE.
      LEFTRF=.FALSE.
      CYNDUN=.FALSE.
      CNTR=1
 
 3400 YMIN=AMAX1((-YE-ELL2)/SGY,-3.)
      Y(1)=AMIN1((ELL2-YE)/SGY,3.)
      NCTRIB=(Y(1).LT.-3 .OR. YMIN.GT.3)
      IF (MFLG.EQ.0 .AND. NCTRIB) GOTO 3830
      IF (MFLG.EQ.0 .AND. .NOT.CYNBLF) GOTO 3405
      MFLG=1
      IF (NCTRIB .AND. CNTR.EQ.1) GOTO 3830     ! ELEMENT DOESNT CONTRIBUTE
      LEFTRF= NCTRIB .AND. CANYON .AND. .NOT.CYNDUN
      IF (LEFTRF) GOTO 3535                     ! REFLECTIONS ON LEFT SIDE
      IF (NCTRIB) GOTO 3550                     ! NO MORE CONTRIBUTION 
 3405 DO 3410 I=1,6
      WT(I)=0.
 3410 CONTINUE
      IF (FET.GE.CSL2) GOTO 3430                ! RCP NOT COMPLETELY DOWNWIND
      IF (ECLD.LT.-EL2) YV1=YE
      IF (ECLD.GE.-EL2) YV1=(ECLD+EL2)/SPHI+YE
      YV2=-(D+W2)/CPHI+YE
      IF (ECLD.LT.-EL2) YH=-(D+W2)*CPHI+YE
      IF (ECLD.GE.-EL2) YH=YV1-(YV1-YV2)*CPHI**2
      YMIN=AMAX1(YMIN,(YV2-YE)/SGY)
      Y(1)=AMIN1(Y(1),(YV1-YE)/SGY)
      NCTRIB=(Y(1).LT.-3 .OR. YMIN.GT.3.)
      IF (MFLG.EQ.0 .AND. NCTRIB) GOTO 3830
      IF (MFLG.EQ.0 .AND. .NOT.CYNBLF) GOTO 3430
      MFLG=1
      IF (NCTRIB .AND. CNTR.EQ.1) GOTO 3830     ! ELEMENT DOESNT CONTRIBUTE
      LEFTRF=NCTRIB .AND. CANYON .AND. .NOT.CYNDUN
      IF (LEFTRF) GOTO 3535                     ! REFLECTIONS ON LEFT SIDE OF CANYON
      IF (NCTRIB) GOTO 3550                     ! NO MORE CONTRIBUTION FROM BLUFF OR CANYON
 3430 INTG(1)=PDENS(Y(1))
      IF (Y(1)-INT(Y(1)).EQ.0) GOTO 3440
      IF (Y(1).LT.0) Y(2)=INT(Y(1))-1
      IF (Y(1).GE.0) Y(2)=INT(Y(1))
      GOTO 3445
 3440 Y(2)=Y(1)-1
 3445 DO 3500 I=2,7
      IF (Y(I).LT.YMIN) Y(I)=YMIN
      INTG(I)=PDENS(Y(I))                       ! ADJUSTMENT FOR ELEMENT END EFFECT
      SED1=Y(I-1)*SGY+YE
      SED2=Y(I)*SGY+YE
      IF (SED1.EQ.SED2) GOTO 3470
      LIM1=AMAX1(SED2,EM2)
      ULIM2=AMIN1(SED1,EM2)
      LLIM2=AMAX1(SED2,-EM2)
      LIM3=AMIN1(SED1,-EM2)
      ZON1=0.
      ZON2=0.
      ZON3=0.
      IF (FET.GE.CSL2) GOTO 3460                ! RCP NOT COMPLETELY DOWNWIND
      IF (SED1.LE.EM2) GOTO 3450
      ASED=(SED1+LIM1)/2.
      IF (ASED.GT.YH) PSSF=(YV1-ASED)*TPHI/(2.*CSL2)
      IF (ASED.LE.YH) PSSF=(ASED-YV2)/(TPHI*2.*CSL2)
      SSF=AMIN1(PSSF,(ELL2-ASED)/EN2)
      ZON1=(SED1-LIM1)*SSF
 3450 IF (SED1.LE.-EM2 .OR. SED2.GE.EM2) GOTO 3455
      ASED=(ULIM2-LLIM2)/2.
      IF (ASED.GT.YH) PSSF=(YV1-ASED)*TPHI/(2.*CSL2)
      IF (ASED.LE.YH) PSSF=(ASED-YV2)/(TPHI*2.*CSL2)
      ZON2=(ULIM2-LLIM2)*AMIN1(PSSF,1.)
 3455 IF(SED2.GE.-EM2) GOTO 3465
      ASED=(SED2+LIM3)/2.
      IF (ASED.GT.YH) PSSF=(YV1-ASED)*TPHI/(2.*CSL2)
      IF (ASED.LE.YH) PSSF=(ASED-YV2)/(TPHI*2.*CSL2)
      SSF=AMIN1(PSSF,(ELL2+ASED)/EN2)
      ZON3=(LIM3-SED2)*SSF
      GOTO 3465
 3460 IF (SED1.GT.EM2) ZON1=(SED1-LIM1)*((ELL2-(SED1+LIM1)/2.)/EN2)
      IF (SED1.GT.-EM2 .AND. SED2.LT.EM2) ZON2=ULIM2-LLIM2
      IF (SED2.LT.-EM2) ZON3=(LIM3-SED2)*((ELL2+(SED2+LIM3)/2.)/EN2)
 3465 WT(I-1)=(ZON1+ZON2+ZON3)/(SED1-SED2)
      IF (WT(I-1).LT.0.0) WT(I-1)=0.0
      IF (Y(I).EQ.YMIN) GOTO 3510
 3470 IF (Y(I).NE.YMIN) Y(I+1)=Y(I)-1
 3500 CONTINUE
 3510 DO 3560 I=1,6
      IF (WT(I).EQ.0) GOTO 3535
      IF ((SIGN(1.,Y(I))).EQ.(SIGN(1.,Y(I+1))))
     *    PD=DABS(INTG(I+1)-INTG(I)) 
      IF ((SIGN(1.,Y(I))).NE.(SIGN(1.,Y(I+1)))
     * .OR. (Y(I+1).EQ.0.))                           ! DH: UNREACHABLE?
     *    PD=1.-INTG(I)-INTG(I+1)
      FAC2=FAC2+PD*QE*WT(I)
 3530 CONTINUE
C
C *****  CANYON AND BLUFF CONTRIBUTIONS *****
C
 3535 IF (.NOT.CYNBLF) GOTO 3550
      IF (MIXWR(IL).EQ.0) GOTO 3546
      IF (MIXWL(IL).EQ.0) GOTO 3542
      IF (LEFTRF .OR. .NOT.RTSIDE) GOTO 3538
      YE=2*MIXWR(IL)-YE
      RTSIDE=.FALSE.
      CNTR=CNTR+1
      GOTO 3400
 3538 IF (.NOT.LEFTRF) GOTO 3540
      YE=TYE
      CYNDUN=.TRUE.
 3540 YE=2*MIXWL(IL)-YE
      RTSIDE=.TRUE.
      CNTR=CNTR+1
      GOTO 3400
 3542 IF (CNTR.GE.2) GOTO 3550
      YE=2*MIXWR(IL)-YE
      CNTR=2
      GOTO 3400
 3546 IF (CNTR.GE.2) GOTO 3550
      YE=2*MIXWL(IL)-YE
      CNTR=2
      GOTO 3400
      
 3550 FACT=FAC1*FAC2
C
C *****  DEPRESSED SECTION  *****
C
      IF (HDS.LT.-1.5 .AND.                                             
     *    DABS(D).LT.(W2-3.*HDS)) GO TO 3560                            
      GO TO 3580                                                        
 3560 IF (DABS(D).LE.W2) FACT=FACT*DSTR                                 
      IF (DABS(D).GT.W2) FACT=FACT*(DSTR-(DSTR-1.)*(DABS(D)-W2)/        
     *     (-3.*HDS))                                                   
C
C *****  DEPOSITION CORRECTION  *****
C
 3580 FAC3=0.                                                           
      IF (V1.EQ.0.) GO TO 3670                                          
      ARG=V1*SGZ/(KZ*SQRT(2.))+(Z+H)/(SGZ*SQRT(2.))                     
      IF (ARG.GT.9.) GO TO 3762                                         
      T=1./(1.+0.47047*ARG)                                             
      EFRC=(.3480242*T-.0958798*T**2+.7478556*T**3)*EXP(-1.*ARG**2)     
      FAC3=(SQRT(2.*PI)*V1*SGZ*EXP(V1*(Z+H)/KZ+.5*(V1*SGZ/KZ)**2)       
     *    *EFRC)/KZ                                                     
      IF (FAC3.GT.2.) FAC3=2.                                           
C
C *****  SETTLING CORRECTION  *****
C
 3670 IF (VS.EQ.0.) GO TO 3710                                          
      FAC4=EXP(-VS*(Z-H)/(2.*KZ)-(VS*SGZ/KZ)**2/8.)                     
      FACT=FACT*FAC4
C
C  *****  INCREMENTAL CONCENTRATION  *****
C
 3710 FAC5=0.                                                           
      CNT=0.                                                            
 3720 EXLS=0.                                                           
 3730 ARG1=-0.5*((Z+H+2.*CNT*MIXH)/SGZ)**2                              
      IF (ARG1.LT.-87.) EXP1=0.                                         
      IF (ARG1.GE.-87.) EXP1=EXP(ARG1)                                  
      ARG2=-0.5*((Z-H+2.*CNT*MIXH)/SGZ)**2                              
      IF (ARG2.LT.-87.) EXP2=0.                                         
      IF (ARG2.GE.-87.) EXP2=EXP(ARG2)                                  
      FAC5=FAC5+EXP1+EXP2                                               
      IF (MIXH.GE.1000. .OR. MIXH.EQ.0.) GOTO 3760                                     
      IF ((EXP1+EXP2+EXLS).EQ.0. .AND. CNT.LE.0.) GO TO 3760
      IF (CNT.GT.0) GOTO 3750
      CNT=ABS(CNT)+1.
      GOTO 3720
 3750 CNT=-1.*CNT
      EXLS=EXP1+EXP2
      GOTO 3730
 3760 INC=FACT*(FAC5-FAC3)                ! INCREMENTAL CONCENTRATION FROM ELEMENT
      C(IR)=C(IR)+INC               ! SUMMATION OF CONCENTRATIONS
      MFLG=1
 3762 IF (FLAGD.EQ.1.) GOTO 4000
 3768 IF (FLAGU.EQ.1.) GOTO 3820
 3770 NE=NE+1.
      STP=1.
      IF (TYP(IL).NE.6) STP=BASE**NE      ! STEP FACTOR
      ED1=ED2
      ED2=ED2+SGN*STP*W                   ! INCREMENT TO NEXT ELEMENT
      GOTO 3110
 3820 NE=0
      SGN=-1.
      FLAGU=0.
      GOTO 3100
 3830 IF (SGN.EQ.1. .AND. YMIN.GT.3. .AND. D.GT.-W2) GOTO 3820
      IF (SGN.EQ.-1. .AND. Y(1).LT.-3.) GOTO 4000
      GOTO 3762
C
C ***** END OF CALINE4 ROUTINE  *****
C
C
C *****  END LOOPS  *****
C
 4000 IF (RTYP.EQ.3 .OR. RTYP.EQ.4) GOTO 5000
      GOTO 3000
 5000 IF (.NOT.CYNBLF) GOTO 6000
      XVEC=TXVEC
      YVEC=TYVEC
 6000 CONTINUE

      DO 7010 I=1,NR
      C(IR)=C(IR)*FPPM
 7010 CONTINUE      

      END SUBROUTINE

C
C  ***** SGZ MODIFICATION SUBROUTINE *****
C
      SUBROUTINE MSGZ (VPH,U,CLAS,WIDTH,COMPT,SGZ)
      REAL AZ(7),SMOD(7,6)
      INTEGER CLAS
      LOGICAL COMPT
      DATA AZ/1112.,566.,353.,219.,124.,56.,22./
      DATA SMOD/0.2,0.27,0.52,99.,99.,99.,99.,
     *        0.38,0.45,0.58,0.7,99.,99.,99.,
     *        0.78,0.86,1.01,1.26,1.36,99.,99.,
     *        2.21,2.38,2.84,3.45,5.9,13.5,27.35,
     *        5.61,6.25,7.59,10.21,16.76,25.51,100.,
     *        9.34,10.82,14.12,21.11,100.,99.,99./
      HFF=6.82    ! HEAT FLUX FACTOR (MW*HR/CM/VEH)
      IF (CLAS.GT.1) GOTO 200
      MCLAS=1
      GOTO 250
 200  HF=HFF/VPH/(100*WIDTH)  ! ROADWAY HEAT FLUX  (MW/CM**2)
      WS=1.5
      DO 210 I=1,6
      IF (U.LT.WS) GOTO 220
      WS=WS+1.
 210  CONTINUE
 220  K=8-CLAS
      KNT=K
      IF (SMOD(I,K).EQ.99.) GOTO 260
      DO 230 J=K,6
      IF (SMOD(I,J).GE.99.) GOTO 240
      HF=HF-SMOD(I,J)
      IF (HF.LT.0) GOTO 240
      KNT=KNT+1
 230  CONTINUE
 240  MCLAS=8-KNT
 250  IF (MCLAS.EQ.1 .AND. U.LT.4.) SGZ=AZ(MCLAS)
      IF (MCLAS.EQ.1 .AND. U.GE.4.) GOTO 260
      IF (MCLAS.NE.1) SGZ=AZ(MCLAS)+(AZ(MCLAS-1)-AZ(MCLAS))*
     *       (HF+SMOD(I,KNT))/SMOD(I,KNT)
      COMPT=.TRUE.
      GOTO 270
 260  COMPT=.FALSE.
      SGZ=99999.
 270  RETURN
      END

C
C *****   PDENS FUNCTION   *****
C
      FUNCTION PDENS(ARG)
      REAL LIM
      LIM=ABS(ARG)
      T=1./(1.+0.23164*LIM)
      XX=LIM**2/(-2.)
      PDENS=0.3989*EXP(XX)*(0.3194*T-0.3566*T**2+1.7815*T**3
     *  -1.8213*T**4+1.3303*T**5)
      RETURN
      END FUNCTION

