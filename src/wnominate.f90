!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!*                 MULTI-DIMENSIONAL W-NOMINATE                        **
!*                                                                     **
!*           MULTI-DIMENSIONAL WEIGHTED NOMINAL THREE-STEP ESTIMATION  **
!*                             *        ******  *          *           **
!*           DEVELOPED BY KEITH T. POOLE AND HOWARD ROSENTHAL          **
!*           GRADUATE SCHOOL OF INDUSTRIAL ADMINISTRATION              **
!*           CARNEGIE-MELLON UNIVERSITY, PITTSBURGH, PA. 15213         **
!*                                                                     **
!*           INITIAL PROGRAM DEVELOPED 1982-1984                       **
!*           MULTIDIMENSIONAL PROGRAM DEVELOPED 1986-1987 IN           **
!*              CDC VECTOR FORTRAN                                     **
!*           OS2/SCALAR FORTRAN VERSION WRITTEN BY NOLAN MCCARTY       **
!*              AND KEITH T. POOLE 1991                                **
!*           THE MODIFIED NOMINATE ALGORITHM, W-NOMINATE, DEVELOPED BY **
!*              KEITH T. POOLE SEPTEMBER-NOVEMBER 1992.                **
!*                                                                     **
!*           PARAMETRIC BOOTSTRAP ALGORITHM DEVELOPED BY JEFF LEWIS    **
!*              AND KEITH T. POOLE 2003-04                             ** 
!*                                                                     **
!*           R IMPLEMENTATION DEVELOPED BY JEFF LEWIS AND KEITH T.     **
!*              POOLE JULY-AUGUST 2006                                 **
!*                                                                     **
!  FLOW DIAGRAM OF W-NOMINATE PROGRAM-- * INDICATES MAIN DO - LOOP
!
!     MAIN
!      |
!      |---CLEAN
!      |
!      |---KPASCORE
!      |
!      |
!      |            |--STAT 
!      |---WHOOPE---|--FOCUSW
!      |            |--FOCUS
!      |            |--STAT  
!      |
!      |---RECODE--JANICE
!      |
!      |            |--NORMZ
!     *|            |          |--FUNNEL
!     *|---MAXLNL---|--BHHH  --|--STEPR --FUNNEL--LOGLIK--ITHOBS
!     *|            |            
!     *|            |--RPRINT--|--GRID2
!     *|            |          |--CORR
!     *|            |--CROSS
!     *|
!     *|---OUTWRT--CORR
!      |
!      |---CROSS
!      |
!
!
!
!    PASSED BACK TO R:  $uweights
!                       $classify
!                       $fits
!                       $gmp
!                       $idealpoints(# dimensions * # members)
!                       $midpoints(# dimensions * # roll calls)
!                       $spreads(# dimensions * # roll calls)
!                       $eigenvalues(min(numvotes,nummembers)
!                       $exitstatus
!
!
!
! #$#$#$#$#$#$#$#$#$#$#
!
!   ARRAY SIZES -- FROM WNOMJLEWIS_REALLY_BIG.FOR -- INCLUDES OC
!
!              1719 = # ROWS OF INPUT DATA
!              3438 = 2 * # ROWS OF INPUT DATA
!              5157 = 3 * # ROWS OF INPUT DATA
!             15599 = # COLUMNS OF INPUT DATA
!             15999 = DUAL USE VECTORS =   MAX{# ROWS, # COLUMNS}
!     WS(.)   33398 > 2*MAX(P,Q)
!
! !!! FOR DEBUGGING PURPOSES THE MATRIX POOLE(NP,3*NRCALL)
! !!! IT IS DIMENSIONED POOLE(NP,NRCALL) IN THE CODE BELOW
! !!! TO SAVE SPACE
!
!
!  DO SEARCH-REPLACE ON THESE ARRAY SIZES BELOW
!
! #$#$#$#$#$#$#$#$#$#$#  
!
!    FROM OLD WNOM9707.FOR FORTRAN PROGRAM
!
!                       3709 = NP        
!                       3111 = NRCALL
!                       3711 = DUAL USE = MAX(NP,NRCALL)+11
!
!  ********%$%$%$  IN THE CODE BELOW:  DUAL USE = NP + NRCALL + 111
!   
        SUBROUTINE WNOM(LDATA,NUMMEMBERS,NUMVOTES,DIMS,TRIALS,          &
     &     POLARITY,REALDIMS,UBETA,UWEIGHTS,CLASSIFY,FITS,GMP,          &
     &     IDEALPOINTS,COVARIANCES,MIDPOINTS,SPREADS,EIGENVALUES,       &
     &     EXITSTATUS,LDATA2,GMPGMP,XFITS,ZMAT2,WVEC2,DSTAR,            &
     &     XDATA,XXX,XDATA3,ZMID,DYN,XSAVE,ZSAVE,CSAVE,KAV,KAY,         &
     &     KAN,XD,ISENS,LERIC,PSI,XMEANX,STDDEV,COVX,COVX2,KPJP,        &
     &     YBIGL,YYBIGL,STDDVZ,STDDVX,LMO,POOLE,TRUEX,TRUEX2,           &
     &     TRUEZMID,TRUEDYN,PROBMAT,XTARGET,XMAT0,XSAVE2,XSAVE3,        &
     &     XMAT)
!
!
        INTEGER LDATA(NUMMEMBERS,NUMVOTES),DIMS,TRIALS,POLARITY(DIMS),  &
     &          EXITSTATUS,LDATA2(NUMMEMBERS,NUMVOTES),REALDIMS
!       DOUBLE PRECISION UBETA,UWEIGHTS(DIMS),
        REAL UBETA,UWEIGHTS(DIMS),                                      &
     &          CLASSIFY((NUMMEMBERS+NUMVOTES)*4),                      &
     &          FITS(3*DIMS),                                           &
     &          GMP(NUMMEMBERS+NUMVOTES),                               &
     &          IDEALPOINTS(NUMMEMBERS*DIMS),                           &
     &          COVARIANCES(NUMMEMBERS*((DIMS*(DIMS+1))/2)),            &
     &          MIDPOINTS(NUMVOTES*DIMS),SPREADS(NUMVOTES*DIMS),        &
     &          EIGENVALUES(NUMMEMBERS)
        DIMENSION GMPGMP(NUMMEMBERS+NUMVOTES),                          &
     &            XFITS(3*DIMS),                                        &
     &            ZMAT2(NUMMEMBERS,NUMMEMBERS),WVEC2(NUMMEMBERS),       &
     &            DSTAR(NUMMEMBERS,NUMMEMBERS),                         &
     &            XDATA(NUMMEMBERS,DIMS),XXX(NUMMEMBERS),               &
     &            XDATA3(NUMMEMBERS,DIMS),SSS(100),                     &
     &            ZMID(NUMVOTES,DIMS),DYN(NUMVOTES,DIMS),               &
     &            XSAVE(NUMMEMBERS,2,2),ZSAVE(NUMVOTES,2,2),            &
     &            CSAVE(NUMVOTES,2),KAV(NUMVOTES),KAY(NUMVOTES),        &
     &            KAN(NUMVOTES),XD(NUMMEMBERS),BBB(50),BBBB(25),        &
     &            ISENS(NUMMEMBERS),LERIC(NUMVOTES),                    &
     &            PSI(NUMMEMBERS,NUMVOTES,2),                           &
     &            XMEANX(2*NUMMEMBERS,DIMS),                            &
     &            STDDEV(3*NUMMEMBERS,DIMS),                            &
     &        COVX(2*NUMMEMBERS,DIMS,DIMS),                             &
     &            COVX2(2*NUMMEMBERS,DIMS,DIMS),ALP(2),BTA(2),          &
     &            RR(2),KPJP(NUMMEMBERS+NUMVOTES+111,4),                &
     &            YBIGL(NUMMEMBERS+NUMVOTES+111),                       &
     &            YYBIGL(NUMMEMBERS+NUMVOTES+111),                      &
     &            STDDVZ(NUMVOTES,2,DIMS),STDDVX(NUMMEMBERS,DIMS),      &
     &            LMO(NUMMEMBERS+NUMVOTES+111),                         &
     &            LL(2,2),POOLE(NUMMEMBERS,NUMVOTES),                   &
     &            TRUEX(NUMMEMBERS,DIMS),TRUEX2(NUMMEMBERS,DIMS),       &
     &            TRUEZMID(NUMVOTES,DIMS),TRUEDYN(NUMVOTES,DIMS),       &
     &            PROBMAT(NUMMEMBERS,NUMVOTES),                         &
     &        XTARGET(NUMMEMBERS,DIMS),                                 &
     &            XMAT0(NUMMEMBERS,DIMS),XSAVE2(NUMMEMBERS,DIMS),       &
     &            XSAVE3(NUMMEMBERS,DIMS),FV1(50),FV2(50),              &
     &            XMATNEW(50,50),XMAT(NUMMEMBERS,DIMS),                 &
     &            Y16MIDP(50,50),YHAT(50),UUU(50,50),VVV(50,50),        &
     &            VCOV(50,50)
!
  100 FORMAT(6I6)  
  150 FORMAT(I5,75F10.4)
  200 FORMAT(I5,2X,1000I1)
  660 FORMAT(2I5,1X,36A1,2000F8.4)
 1015 FORMAT(' SINGULAR VALUE DECOMPOSITION',I5)
 1016 FORMAT(I5,50F7.3)
 1017 FORMAT(' ORTHOGONAL PROCRUSTES ROTATION MATRIX')
 1029 FORMAT(7X,'CONSTRAINED',7X,'UNCONSTRAINED'/1X,'OFF OPPOSITE',     &
     &3X,'OFF SAME'/2X,I5,9X,I5,9X,I5)
 1038 FORMAT(' NUMBER OF DIMENSIONS RESET TO=',I4)
 1039 FORMAT(' SAG IN RECOVERED LEGISLATOR COORDINATES')
 1040 FORMAT(I5,4I6,F10.4)
 1090 FORMAT(' R-SQUARES TRUE VS. REPRODUCED',2I4,2F7.3)
 2001 FORMAT(3I5,9F7.3)
!
!    SET VARIABLES AND ARRAYS PASSED IN FROM R TO THEIR NAMES
!        USED IN THE ORIGINAL NOMINATE AND W-NOMINATE FORTRAN
!        CODE
!
!  SET PRINTER SWITCH FOR DEBUGGING PURPOSES -- 1=PRINT, 0=NO PRINT
!
      IPRINT=0
!
!  INITIALIZE RANDOM NUMBER GENERATOR
!
      CALL RNDSTART()
      NTRIAL=TRIALS
      IF(NTRIAL.LT.5)NTRIAL=1
      NS=REALDIMS
      NP=NUMMEMBERS
      NRCALL=NUMVOTES
      NDUAL=NP+NRCALL+111
!
      JANLST=3
      MAXIT=30
      MAXSQZ=10
!
!      UBETA=15.0
!      UWEIGHTS(1)=0.5
!      UWEIGHTS(2)=0.5
!
      XVMIN=0.025
      KVMIN=10
      KSMIN=1
      KSMAX=1
!
!  TRANSPOSED CODE
!
!      DO 1 I=1,NP
!      DO 11 J=1,NRCALL
!      LDATA(I,J)=VOTES(J+(I-1)*NRCALL)
!  11  CONTINUE
!  1   CONTINUE
!
!  UN-TRANSPOSED CODE
!
!      DO 11 J=1,NRCALL
!      DO 191 I=1,NP
!      LDATA(I,J)=VOTES(I+(J-1)*NP)
!  191 CONTINUE
!  11  CONTINUE
!
      IF(IPRINT.EQ.1)THEN
         WRITE(11,100)NS,NP,NRCALL,NDUAL
         DO 291 I=1,NP
         WRITE(11,200)I,(LDATA(I,J),J=1,NRCALL)
  291    CONTINUE
      ENDIF
!
!  SUBROUTINE CLEAN--THROWS OUT ALL VOTES WITH LESS THAN XVMIN MAJ. AND
!    ALL LEGISLATORS VOTING LESS THAN KVMIN TIMES
!
      CALL CLEAN(NP,NRCALL,XVMIN,KVMIN,IPRINT,KPTSUM,LDATA,             &
     &           LDATA2,KAV,KAY,KAN)
!
!  RESET NUMBER OF DIMENSIONS IF IT EXCEEDS NP-1
!
      IF(NS.GT.NP-1)THEN
         NS=NP-1
         IF(IPRINT.EQ.1)WRITE(21,1038)NS
      ENDIF
!
!  ***************************
!
!  BOOTSTRAP LOOP BEGINS HERE
!
!  ***************************
!
!
      DO 77 K=1,NS
      DO 77 I=1,NP
      XMEANX(I,K)=0.0
      STDDEV(I,K)=0.0
      XMEANX(I+NP,K)=0.0
      STDDEV(I+NP,K)=0.0
      STDDEV(I+2*NP,K)=0.0
      DO 77 L=1,NS
      COVX(I,K,L)=0.0
      COVX(I+NP,K,L)=0.0
      COVX2(I,K,L)=0.0
      COVX2(I+NP,K,L)=0.0
  77  CONTINUE
!
      DO 9988 I999=1,NTRIAL
!
! **&*&*&*&*&*&*&*&*&*&*&*&*&
!  WARNING !!!!!!!!!!!!!
!  IN RECODE SUBROUTINE, THE IIII FLAG IS > 1 ON THE
!  2ND PASS, THAT IS, I999 > 1.  THIS FORTUNATELY DOES
!  NOT AFFECT THE CODE BECAUSE LDATA(.,.) USES "2" 
!  RATHER THAN "6" ON THE SECOND PASS!!!!!!!!!!!!!
!  CODE WAS CHANGED IN RECODE TO USE LDATA2(.,.) --
!  IT ALWAYS USES "1" AND "6"
! **&*&*&*&*&*&*&*&*&*&*&*&*&
!
!
!    INITIALIZE PSI MATRICES OF DIMENSIONAL WEIGHTING FACTORS
!
      BBB(1)=UBETA
      BBB(2)=UWEIGHTS(1)
      DO 69 K=1,2
      DO 69 I=1,NP
      ISENS(I)=1
      DO 69 J=1,NRCALL
      DO 68 L=1,NS
      BBBB(L)=BBB(2)
      ZMID(J,L)=0.0
      DYN(J,L)=0.0
  68  CONTINUE
      LERIC(J)=1
 69   PSI(I,J,K) = 1.0
!
!
!  CALCULATE AGREEMENT SCORE MATRIX, DOUBLE-CENTER IT, AND
!     EXTRACT EIGENVECTORS TO OBTAIN LEGISLATOR STARTS
!
      CALL KPASCORE(NP,NRCALL,NS,NDUAL,11,IPRINT,                       &
     &                  ZMAT2,WVEC2,DSTAR,LDATA)
!
      DO 44 I=1,NP
      EIGENVALUES(I)=WVEC2(NP+1-I)
      DO 58 K=1,NS
      XDATA(I,K)=ZMAT2(I,NP+1-K)*SQRT(ABS(WVEC2(NP+1-K)))
!      IDEALPOINTS(I+(K-1)*NP)=XDATA(I,K)
  58  CONTINUE
      IF(IPRINT.EQ.1)WRITE(11,150)I,EIGENVALUES(I),                     &
     &               (XDATA(I,K),K=1,NS)
  44  CONTINUE
!
!  NORMALIZE ESTIMATES OF LEGISLATORS TO BE WITHIN THE
!   UNIT HYPERSPHERE WITH CENTROID = 0
!
!
!  PERFORM METRIC SCALING TO INCREASE PRECISION OF STARTING COORDINATES
!
      CALL WHOOPE(NP,NS,DSTAR,XXX,XDATA,SSE1,SSE2,KTP,IPRINT)
      DO 71 K=1,NS
      SUM=0.0
      BB=-99.0
      DO 70 I=1,NP
      IF(NS.EQ.1)XDATA(I,1)=XXX(I)
  70  SUM=SUM+XDATA(I,K)
      DO 72 I=1,NP
      XDATA(I,K)=XDATA(I,K)-SUM/FLOAT(NP)
  72  BB=AMAX1(BB,XDATA(I,K))
      DO 73 I=1,NP
  73  XDATA(I,K)=XDATA(I,K)/BB 
  71  CONTINUE
!
!  NORMALIZE ESTIMATES OF LEGISLATORS TO BE WITHIN THE
!   UNIT HYPERSPHERE WITH CENTROID = 0
!
      DO 51 K=1,NS
      SUM=0.0
      DO 50 I=1,NP
  50  SUM=SUM+XDATA(I,K)
      DO 52 I=1,NP
      XDATA(I,K)=XDATA(I,K)-SUM/FLOAT(NP)
  52  CONTINUE
  51  CONTINUE
      BB=-99.0
      DO 56 I=1,NP
      SUM=0.0
      DO 54 K=1,NS
      SUM=SUM+XDATA(I,K)**2
  54  CONTINUE
      BB=AMAX1(BB,SUM)
  56  CONTINUE
      DO 55 I=1,NP
      DO 57 K=1,NS
      XDATA(I,K)=XDATA(I,K)*(1.0/SQRT(BB))
      XDATA3(I,K)=XDATA(I,K)
  57  CONTINUE
  55  CONTINUE
!
!  WRITE OUT METRIC SCALING COORDINATES USED AS STARTS
!
      DO 713 I=1,NP
!      WRITE(35,FMT3)I,(KNAME(I,J),J=1,JNE),
!     C                                 (XDATA(I,J),J=1,NS)
      DO 712 K=1,NS
      XMAT0(I,K)=XDATA(I,K)
  712 CONTINUE
  713 CONTINUE
!
!  SUBROUTINE RECODE--GENERATES STARTING VALUES FOR ROLL CALL MIDPOINTS
!    ON THE FIRST DIMENSION AND DETERMINES POLARITY OF ROLL CALL
!
      IIII=1
      KLSEN=POLARITY(1)
      KLSEN2=POLARITY(2)
      CALL RECODE(NS,NP,NRCALL,KLSEN,KLSEN2,11,KSMIN,KSMAX,KPTSUM,SSS,  &
     &            XDATA,ZMID,DYN,LDATA,LDATA2,XSAVE,ZSAVE,CSAVE,        &
     &            KAV,KAY,KAN,IIII,IPRINT)
!
!  OPTIONAL WRITE TO DISK
!
      DO 715 I=1,NP
      IF(IPRINT.EQ.1)WRITE(11,150)I,(XDATA(I,K),K=1,NS)
  715 CONTINUE
!
      DO 60 K=1,NS
      DO 61 J=1,NRCALL
!      MIDPOINTS(J+(K-1)*NRCALL)=ZMID(J,K)
!      SPREADS(J+(K-1)*NRCALL)=DYN(J,K)
!      IF(IPRINT.EQ.1)WRITE(11,150)J,MIDPOINTS(J+(K-1)*NRCALL),
!     C                            SPREADS(J+(K-1)*NRCALL)
  61  CONTINUE
  60  CONTINUE
!
!   DIMENSIONAL LOOP--ESTIMATES EACH DIMENSION SEPARATELY
!
      DO 2222 NNDIM=1,NS
      NDIM=NNDIM
!
!   FIND MAXIMUM AND MINIMUM LEGISLATORS FOR CURRENT DIMENSION STARTS AND
!       CALCULATE SHRINKAGE PARAMETER FOR CURRENT DIMENSION STARTS TO
!       CONFORM WITH UNIT HYPERSPHERE CONSTRAINT -- ALPHA IS THE SHRINKAGE
!       PARAMETER
!
      AA=99.0
      BB=-99.0
      ASUM=0.0
      BSUM=0.0
      DO 1 I=1,NP
      IF(NDIM.GE.2)THEN
          SUM=0.0
          DO 2 K=1,NDIM-1
            SUM=SUM+XDATA(I,K)**2
  2       CONTINUE
          XEPS=SQRT(ABS(1.0-SUM))
          ASUM=ASUM+ABS(XDATA(I,NDIM)*XEPS)
          BSUM=BSUM+XDATA(I,NDIM)**2
      ENDIF
      AA=AMIN1(AA,XDATA(I,NDIM))
      BB=AMAX1(BB,XDATA(I,NDIM))
      IF(ABS(AA-XDATA(I,NDIM)).LE..0001)KSMIN=I
      IF(ABS(BB-XDATA(I,NDIM)).LE..0001)KSMAX=I
  1   CONTINUE
!
!  CORRECT FOR SHRINKAGE
!
      ALPHA=1.0
      IF(NDIM.GE.2)THEN
           ALPHA=AMAX1(XD(1),ABS(XD(NP)))
           ALPHA=ALPHA/XDATA(KSMAX,NDIM)
!          ALPHA=ASUM/BSUM
          DO 22 I=1,NP
  22      XDATA(I,NDIM)=XDATA(I,NDIM)*ALPHA
      ENDIF
      SSS(1)=XDATA(KSMIN,NDIM)
      SSS(2)=XDATA(KSMAX,NDIM)
!
      IF(IPRINT.EQ.1)WRITE(21,2001)NDIM,KSMIN,KSMAX,ALPHA,SSS(1),SSS(2)
!
!  **UNIVERSAL LOOP**
!     3-STEP ESTIMATION; BETA, ROLL CALLS, LEGISLATORS
!     JANLST IS # OF UNIVERSAL ITERATIONS.
!
      KBAD1=0
      KBAD2=0
      JAN=0
      NITR=0
      DO 9999 IIII=1,JANLST
      JAN=JAN+1
      DO 6666 NSTEP=-1,2
!
!     UTILITY FUNCTION PHASE--BETA ESTIMATED
!
      IF(NSTEP.EQ.0.AND.IIII.GT.1.AND.NNDIM.EQ.1)THEN
           IF(I999.EQ.1)THEN
              CALL ECHOEVENT(4)
              CALL FLUSHCON()
              CALL PROCEVENT()
           ENDIF
           NOPAR=1
           NITR=NITR+1
           CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,  &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
      ENDIF
!
!     ROLL CALL PHASE--MIDPOINT AND DISTANCE BETWEEN TWO OUTCOMES ESTIMATED
!
      IF(NSTEP.EQ.1)THEN
           IF(I999.EQ.1)THEN
              CALL ECHOEVENT(2)
              CALL FLUSHCON()
              CALL PROCEVENT()
           ENDIF
           NOPAR=2
           NITR=NITR+1
           CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,  &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
!
           CALL OUTWRT(NRCALL,NSTEP,JAN,999,KPJP,YBIGL,NITR,            &
     &                  NDIM,NS,NP,NRCALL,NDUAL,IPRINT,XDATA,STDDVX,    &
     &                  KAV,DYN,ZMID,STDDVZ,XSAVE,ZSAVE,CSAVE,          &
     &                  KBAD1,KBAD2,SPRM,SPRSD,GMPGMP,XFITS)
!
      ENDIF
!
!     LEGISLATOR PHASE--COORDINATE ESTIMATED
!
      IF(NSTEP.EQ.2)THEN
           IF(I999.EQ.1)THEN
              CALL ECHOEVENT(1)
              CALL FLUSHCON()
              CALL PROCEVENT()
           ENDIF
           NOPAR=1
           NITR=NITR+1
           CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,  &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
!
           CALL OUTWRT(NP,NSTEP,JAN,999,KPJP,YYBIGL,NITR,               &
     &                  NDIM,NS,NP,NRCALL,NDUAL,IPRINT,XDATA,STDDVX,    &
     &                  KAV,DYN,ZMID,STDDVZ,XSAVE,ZSAVE,CSAVE,          &
     &                  KBAD1,KBAD2,SPRM,SPRSD,GMPGMP,XFITS)
      ENDIF
!
!     UTILITY FUNCTION PHASE--WEIGHT ESTIMATED (S => 2)
!
      IF(NSTEP.EQ.-1.AND.IIII.GT.1.AND.NNDIM.GT.1)THEN
           IF(I999.EQ.1)THEN
              CALL ECHOEVENT(5)
              CALL FLUSHCON()
              CALL PROCEVENT()
           ENDIF
           NOPAR=1
           NITR=NITR+1
           CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,  &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
      ENDIF
!
 6666 CONTINUE
!
!  CALCULATE NUMBER OF CONSTRAINED ROLL CALLS AND CORRELATION BETWEEN
!  CURRENT PARAMETER ESTIMATES AND THOSE FROM PREVIOUS ITERATIONS
!
      KRSCLL=NRCALL-KBAD1-KBAD2
      CALL CORR(NP,NRCALL,XSAVE,ZSAVE,CSAVE,ALP,BTA,RR,0)
      IF(IPRINT.EQ.1)WRITE(23,2001)KBAD2,KBAD1,KRSCLL,SPRM,SPRSD,       &
     &                  ALP(1),BTA(1),RR(1),                            &
     &                  SSS(1),SSS(3),SSS(4),SSS(2)
 9999 CONTINUE
!
!  END OF UNIVERSAL LOOP
!    SETUP PARAMETERS FOR FINAL ESTIMATE OF BETA/WEIGHT AND ROLL CALLS USING
!    FINAL ADJUSTED LEGISLATOR COORDINATES
!
!  FINAL ESTIMATION OF BETA (1ST DIMENSION ONLY)
!
      IF(NDIM.EQ.1)THEN
         JAN=JAN+1
         NSTEP=0
         NOPAR=1
         NITR=NITR+1
         CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,    &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
      ENDIF
!
!  FINAL ESTIMATION OF WEIGHT (2ND DIMENSION AND HIGHER)
!
      IF(NDIM.GT.1)THEN
         JAN=JAN+1
         NSTEP=-1
         NOPAR=1
         NITR=NITR+1
         CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,    &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
      ENDIF
!
!  FINAL ESTIMATION OF ROLL CALL PARAMETERS
!
      NSTEP=1
      NOPAR=2
      NITR=NITR+1
      CALL MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST,       &
     &                 NITR,KSMIN,KSMAX,KPTSUM,                         &
     &                 BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,           &
     &                 YBIGL,YYBIGL,XFITS,                              &
     &                 XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                 &
     &                 KPJP,LERIC,ISENS,IPRINT)
      CALL OUTWRT(NRCALL,NSTEP,JAN,9999,KPJP,YBIGL,NITR,                &
     &                  NDIM,NS,NP,NRCALL,NDUAL,IPRINT,XDATA,STDDVX,    &
     &                  KAV,DYN,ZMID,STDDVZ,XSAVE,ZSAVE,CSAVE,          &
     &                  KBAD1,KBAD2,SPRM,SPRSD,GMPGMP,XFITS)
!      CALL OUTWRT(NRCALL,NSTEP,JIO,KIO,JAN,9999,KPJP,YBIGL,
!     C     FMT2,NDIM,NS)
!
!
!
!  SETUP PARAMETERS TO CALCULATE CLASSIFICATIONS FOR LEGISLATORS
!
      NSTEP=2
      DO 3 II=1,NP
      LMO(II)=II
      JJ=II
!
!     SUBROUTINE CROSS COMPUTES 2X2 PREDICTIONS FOR ROLL CALLS AND LEGS.
!
      CALL CROSS(JJ,LL,                                                 &
     &                 NP,NRCALL,NDIM,NSTEP,JAN,JANLST,                 &
     &                 ZMID,XDATA,DYN,LDATA,PSI,BBB,BBBB)
!      CALL CROSS(JJ,LL)
      JJ=1
      DO 4 I=1,2
      DO 4 J=1,2
      KPJP(II,JJ)=LL(I,J)
      IF(I999.EQ.1)THEN
         CLASSIFY((II-1)*4 + JJ)=LL(I,J)
      ENDIF
  4   JJ=JJ+1
  3   CONTINUE
!
!  WRITE OUT FINAL ESTIMATE OF LEGISLATORS WITH CROSS CLASSIFICATIONS
!
      CALL OUTWRT(NP,NSTEP,999,999,KPJP,YYBIGL,NITR,                    &
     &                  NDIM,NS,NP,NRCALL,NDUAL,IPRINT,XDATA,STDDVX,    &
     &                  KAV,DYN,ZMID,STDDVZ,XSAVE,ZSAVE,CSAVE,          &
     &                  KBAD1,KBAD2,SPRM,SPRSD,GMPGMP,XFITS)
!      CALL OUTWRT(NP,NSTEP,KIO,KIO,999,999,KPJP,YYBIGL,FMT2,NDIM,NS)
!
!  SORT LEGISLATORS AND WRITE OUT THE 6 MOST EXTREME LEGISLATORS
!
      DO 1041 KK=1,NP
      XD(KK)=XDATA(KK,NDIM)
      LMO(I)=I
 1041 CONTINUE
      CALL RSORT(XD,NP,LMO)
      IF(IPRINT.EQ.1)WRITE(21,1039)
      DO 102 I=1,NP
!
!  WRITE OUT PARTIAL LIST OF LEGISLATOR COORDINATES--NOTE THAT
!    THIS FORMAT IS UNIQUE TO THE ROLL CALL DATA -- FOR DIFFERENT
!    MATRICES, THIS CODE WILL HAVE TO BE COMMENTED OUT -- SEE
!    EXAMPLE BELOW
!
      IF(I.GT.3.AND.I.LT.NP-2)GO TO 102
      IF(IPRINT.EQ.1)THEN
      WRITE(21,1040)ISENS(LMO(I)),                                      &
     &     (KPJP(LMO(I),J),J=1,4),XD(I)
      ENDIF
  102 CONTINUE
!
!  SETUP PARAMETERS TO CALCULATE CLASSIFICATIONS FOR ROLL CALLS
!
      NSTEP=1
      DO 5 II=1,NRCALL
      LMO(II)=II
      JJ=II
      CALL CROSS(JJ,LL,                                                 &
     &                 NP,NRCALL,NDIM,NSTEP,JAN,JANLST,                 &
     &                 ZMID,XDATA,DYN,LDATA,PSI,BBB,BBBB)
!      CALL CROSS(JJ,LL)
      JJ=1
      DO 6 I=1,2
      DO 6 J=1,2
      KPJP(II,JJ)=LL(I,J)
      IF(I999.EQ.1)THEN
         CLASSIFY(NP*4 + (II-1)*4 + JJ)=LL(I,J)
      ENDIF
  6   JJ=JJ+1
  5   CONTINUE
!
!  WRITE OUT FINAL ESTIMATES OF ROLL CALLS WITH CROSS CLASSIFICATIONS
!
      CALL OUTWRT(NRCALL,NSTEP,999,999,KPJP,YBIGL,NITR,                 &
     &                  NDIM,NS,NP,NRCALL,NDUAL,IPRINT,XDATA,STDDVX,    &
     &                  KAV,DYN,ZMID,STDDVZ,XSAVE,ZSAVE,CSAVE,          &
     &                  KBAD1,KBAD2,SPRM,SPRSD,GMPGMP,XFITS)
!      CALL OUTWRT(NRCALL,NSTEP,KIO,KIO,999,999,KPJP,YBIGL,FMT2,NDIM,NS)
!
!  WRITE OUT FINAL NUMBERS OF CONSTRAINED ROLL CALLS
!
      KRSCLL=NRCALL-KBAD1-KBAD2
      IF(IPRINT.EQ.1)WRITE(21,1029)KBAD2,KBAD1,KRSCLL
!
!    UPDATE THE PSI MATRICES FOLLOWING THE ESTIMATION OF COORDINATES
!       OF DIMENSION NDIM
!
      DO 1492 I=1,NP
      ISENS(I)=1
      DO 1492 J=1,NRCALL
      DY=0.0
      DN=0.0
!      DYY=-(BBB(2)*(XDATA(I,NDIM)-ZMID(J,NDIM)+DYN(J,NDIM)))**2
!      DNN=-(BBB(2)*(XDATA(I,NDIM)-ZMID(J,NDIM)-DYN(J,NDIM)))**2
      DO 1493 K=1,NDIM
      DY=DY+(BBBB(K)*(XDATA(I,K)-ZMID(J,K)+DYN(J,K)))**2
      DN=DN+(BBBB(K)*(XDATA(I,K)-ZMID(J,K)-DYN(J,K)))**2
 1493 CONTINUE
      DY=EXP(-DY/2.0)
      DN=EXP(-DN/2.0)
!      DYY=EXP(DYY/2.0)
!      DNN=EXP(DNN/2.0)
!      IF(DY.NE.DYY.OR.DN.NE.DNN)THEN
!         WRITE(*,1911)I,J,DY,DYY,DN,DNN
! 1911 FORMAT(2I5,4F10.4)
!         STOP
!      ENDIF
      PSI(I,J,1)=DY
      PSI(I,J,2)=DN
 1492 CONTINUE
      BBB(2)=.5
 2222 CONTINUE
!
!
!  COMPUTE PROBABILITY MATRIX FOR CHOICES
!
!
      IF(I999.EQ.1)THEN
         BETA=BBB(1)
!
!  PASS BACK THE BETA AND WEIGHT PARAMETERS TO R
!
         UBETA=BETA
         DO 9696 K=1,NS
         UWEIGHTS(K)=BBBB(K)
 9696    CONTINUE
!
!
!  TRANSFER ARRAYS THAT ARE PASSED BACK TO R
!
!
         DO 782 K=1,NS
         DO 781 I=1,NP
         GMP(I)=GMPGMP(I)
         IDEALPOINTS(I+(K-1)*NP)=XDATA(I,K)
  781    CONTINUE
  782    CONTINUE
         DO 786 K=1,NS
         DO 783 I=1,NRCALL
         GMP(I+NP)=GMPGMP(I+NP)
         MIDPOINTS(I+(K-1)*NRCALL)=ZMID(I,K)
         SPREADS(I+(K-1)*NRCALL)=DYN(I,K)
  783    CONTINUE
  786    CONTINUE
         DO 787 I=1,3*NS
         FITS(I)=XFITS(I)
  787    CONTINUE
!
           DO 109 J=1,NRCALL
           DO 109 I=1,NP
           ICH=LDATA(I,J)
           POOLE(I,J)=-1.0
           IF(ICH.LE.0)GO TO 109
!
           SUMY=0.0
           SUMN=0.0
           DO 96 K=1,NS
           SUMY=SUMY+(BBBB(K)*(XDATA(I,K)-ZMID(J,K)+DYN(J,K)))**2
           SUMN=SUMN+(BBBB(K)*(XDATA(I,K)-ZMID(J,K)-DYN(J,K)))**2
!
!  STORE ESTIMATED ("TRUE") COORDINATES
!
           TRUEX(I,K)=XDATA(I,K)
           TRUEX2(I,K)=XMAT0(I,K)
           TRUEZMID(J,K)=ZMID(J,K)
           TRUEDYN(J,K)=DYN(J,K)
!
   96      CONTINUE
           PART1=EXP(-SUMY/2.0)
           PART2=EXP(-SUMN/2.0)
           SUMY=BETA*EXP(-SUMY/2.0)
           SUMN=BETA*EXP(-SUMN/2.0)
           PIY=EXP(SUMY)/(EXP(SUMY)+EXP(SUMN))
           PROBMAT(I,J)=PIY
           IF(ICH.EQ.1)THEN
              POOLE(I,J)=PIY
!              POOLE(I,J+NRCALL)=PART1
!              POOLE(I,J+2*NRCALL)=PART2
           ENDIF
           IF(ICH.EQ.2)THEN
              POOLE(I,J)=1.0-PIY
!              POOLE(I,J+NRCALL)=PART2
!              POOLE(I,J+2*NRCALL)=PART1
           ENDIF
  109      CONTINUE
           DO 659 I=1,NP
           DO 659 J=1,NRCALL
!           WRITE(66,660)I,J,
!     C               POOLE(I,J),POOLE(I,J+NRCALL),POOLE(I,J+2*NRCALL)
  659      CONTINUE
!
           DO 1427 I=1,NP
           DO 1428 J=1,NRCALL
           ICH=LDATA(I,J)
           IF(ICH.LE.0)GO TO 1428
!
!           COINFLIP=RNUNF()
           COINFLIP=RANDOM()
!
           IF(COINFLIP.LE.PROBMAT(I,J))THEN
              LDATA(I,J)=1
              LDATA2(I,J)=1
           ENDIF
           IF(COINFLIP.GT.PROBMAT(I,J))THEN
              LDATA(I,J)=2
              LDATA2(I,J)=6
           ENDIF
 1428      CONTINUE
!           WRITE(66,200)I,
!     C                         (LDATA2(I,JJ),JJ=1,NRCALL)
 1427      CONTINUE           
      ENDIF
!
!  DRAW NEW ROLL CALL MATRIX
!
      IF(I999.GT.1)THEN
           IF(I999.EQ.2)THEN
              CALL ECHOEVENT(3)
              CALL FLUSHCON()
              CALL PROCEVENT()
           ENDIF
           IF(I999.GT.2)THEN
              CALL ECHOEVENT(0)
              CALL FLUSHCON()
              CALL PROCEVENT()
           ENDIF
!
           DO 108 J=1,NRCALL
           DO 108 I=1,NP
           ICH=LDATA(I,J)
           IF(ICH.LE.0)GO TO 108
!
!           COINFLIP=RNUNF()
           COINFLIP=RANDOM()
!
           IF(COINFLIP.LE.PROBMAT(I,J))THEN
              LDATA(I,J)=1
              LDATA2(I,J)=1
           ENDIF
           IF(COINFLIP.GT.PROBMAT(I,J))THEN
              LDATA(I,J)=2
              LDATA2(I,J)=6
           ENDIF
  108      CONTINUE           
           DO 107 I=1,NP
           DO 107 K=1,NS
           XTARGET(I,K)=XDATA(I,K)
  107      CONTINUE
           CALL ROTATE(NP,NS,XTARGET,TRUEX,IPRINT)
           DO 654 I=1,NP
           DO 654 K=1,NS
           STDDEV(I,K)=STDDEV(I,K)+(TRUEX(I,K)-XTARGET(I,K))**2
           XMEANX(I,K)=XMEANX(I,K)+XTARGET(I,K)
           DO 6541 L=1,NS
           COVX(I,K,L)=COVX(I,K,L)+XTARGET(I,K)*XTARGET(I,L)
 6541      CONTINUE
  654      CONTINUE
           DO 1070 I=1,NP
           DO 1070 K=1,NS
           XTARGET(I,K)=XMAT0(I,K)
 1070      CONTINUE
           CALL ROTATE(NP,NS,XTARGET,TRUEX2,IPRINT)
           DO 6540 I=1,NP
           DO 6540 K=1,NS
           STDDEV(I+NP,K)=STDDEV(I+NP,K)+                               &
     &                          (TRUEX2(I,K)-XTARGET(I,K))**2
           STDDEV(I+2*NP,K)=STDDEV(I+2*NP,K)+                           &
     &                          XTARGET(I,K)**2
           XMEANX(I+NP,K)=XMEANX(I+NP,K)+XTARGET(I,K)
           DO 6542 L=1,NS
           COVX(I+NP,K,L)=COVX(I+NP,K,L)+XTARGET(I,K)*XTARGET(I,L)
 6542      CONTINUE
 6540      CONTINUE
      ENDIF
!
!
!  ROTATE HECKMAN-SNYDER COORDINATES AND METRIC SCALING COORDINATES TO
!     BEST FIT THE W-NOMINATE COORDINATES AND CALCULATE FITS
!
      DO 2221 JJJJ=2,2
      DO 785 I=1,NP
      DO 784 K=1,NS
      IF(JJJJ.EQ.2)XMAT(I,K)=XDATA3(I,K)
  784 CONTINUE
  785 CONTINUE
!
!  ROTATE TRUE COORDINATES TO BEST FIT ESTIMATED COORDINATES USING
!    SCHONEMANN'S ORTHOGONAL PROCRUSTES SOLUTION
!
      DO 376 L=1,NS
      DO 377 K=1,NS
      SUM=0.0
      DO 378 I=1,NP
      SUM=SUM+XMAT(I,K)*XDATA(I,L)
  378 CONTINUE
      Y16MIDP(K,L)=SUM
  377 CONTINUE
  376 CONTINUE
!
!  CALL SINGULAR VALUE DECOMPOSITION ROUTINE
!
      XTOL=.001
      CALL SVDSVD(NS,NS,Y16MIDP,YHAT,UUU,VVV,IRANK,IPRINT)
!
      IF(IPRINT.EQ.1)WRITE(23,1015)IRANK
!
!  COMPUTE ROTATION MATRIX
!
      IF(IPRINT.EQ.1)WRITE(23,1017)
      DO 379 K=1,NS
      DO 380 M=1,NS
      SUM=0.0
      DO 381 L=1,NS
      SUM=SUM+UUU(K,L)*VVV(M,L)
  381 CONTINUE
      VCOV(K,M)=SUM
  380 CONTINUE
      IF(IPRINT.EQ.1)WRITE(23,1016)K,YHAT(K),(VCOV(K,L),L=1,NS) 
  379 CONTINUE
!
!  ROTATE COORDINATES
!
      DO 382 I=1,NP
      DO 383 K=1,NS
      SUM=0.0
      DO 384 L=1,NS
      SUM=SUM+XMAT(I,L)*VCOV(L,K)
  384 CONTINUE
      IF(JJJJ.EQ.1)XSAVE2(I,K)=SUM
      IF(JJJJ.EQ.2)XSAVE3(I,K)=SUM
  383 CONTINUE
  382 CONTINUE
      DO 385 I=1,NP
      DO 386 K=1,NS
      IF(JJJJ.EQ.1)XMAT(I,K)=XSAVE2(I,K)
      IF(JJJJ.EQ.2)XMAT(I,K)=XSAVE3(I,K)
  386 CONTINUE
  385 CONTINUE
!
!  CALCULATE R-SQUARE BETWEEN ESTIMATED LEGISLATOR COORDINATES AND
!    TRUE COORDINATES
!
      DO 275 JJ=1,50
      FV1(JJ)=0.0
      FV2(JJ)=0.0
      DO 276 JK=1,25
      XMATNEW(JJ,JK)=0.0
  276 CONTINUE
  275 CONTINUE
      DO 270 I=1,NP
      DO 271 K=1,NS
      FV1(K)=FV1(K)+XDATA(I,K)
      FV2(K)=FV2(K)+XMAT(I,K)
      FV1(K+NS)=FV1(K+NS)+XDATA(I,K)**2
      FV2(K+NS)=FV2(K+NS)+XMAT(I,K)**2
      DO 272 L=1,NS
      XMATNEW(K,L)=XMATNEW(K,L)+XDATA(I,K)*XMAT(I,L)
  272 CONTINUE
  271 CONTINUE
  270 CONTINUE
      DO 273 K=1,NS
      DO 274 L=1,NS
      AA=FLOAT(NP)*XMATNEW(K,L)-FV1(K)*FV2(L)
      AB=FLOAT(NP)*FV1(K+NS)-FV1(K)*FV1(K)
      AC=FLOAT(NP)*FV2(L+NS)-FV2(L)*FV2(L)
      IF(AB*AC.GT.0.0)RSQR1=(AA*AA)/(AB*AC)
      IF(IPRINT.EQ.1)WRITE(23,1090)K,L,RSQR1
  274 CONTINUE
  273 CONTINUE
 2221 CONTINUE
      DO 2220 I=1,NP
      IF(IPRINT.EQ.1)THEN
         WRITE(36,150)I,(XDATA(I,J),J=1,NS),                            &
     &                 (XSAVE3(I,J),J=1,NS)
      ENDIF
 2220 CONTINUE
!
!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!    END OF PARAMETRIC BOOTSTRAP LOOP
!
!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
 9988 CONTINUE
!
      IF(NTRIAL.EQ.1)THEN
              CALL ECHOEVENT(7)
              CALL FLUSHCON()
              CALL PROCEVENT()
      ENDIF
      IF(NTRIAL.GT.1)THEN
              CALL ECHOEVENT(6)
              CALL FLUSHCON()
              CALL PROCEVENT()
!
!  WRITE OUT COORDINATES WITH BOOTSTRAP STANDARD ERRORS IF NUMBER OF TRIALS >= 5
!
      DO 6567 I=1,NP
      DO 6568 K=1,NS
      DO 6569 L=1,NS
      AA=FLOAT(NTRIAL-1)*COVX(I,K,L)-XMEANX(I,K)*XMEANX(I,L)
      BB=FLOAT(NTRIAL-1)*COVX(I,K,K)-XMEANX(I,K)*XMEANX(I,K)
      CC=FLOAT(NTRIAL-1)*COVX(I,L,L)-XMEANX(I,L)*XMEANX(I,L)
      AAA1=FLOAT(NTRIAL-1)*COVX(I+NP,K,L)-                              &
     &                     XMEANX(I+NP,K)*XMEANX(I+NP,L)
      BBB1=FLOAT(NTRIAL-1)*COVX(I+NP,K,K)-                              &
     &                     XMEANX(I+NP,K)*XMEANX(I+NP,K)
      CCC1=FLOAT(NTRIAL-1)*COVX(I+NP,L,L)-                              &
     &                     XMEANX(I+NP,L)*XMEANX(I+NP,L)
!      WRITE(99,6549)I,AA,BB,CC,AAA1,BBB1,CCC1
 6549 FORMAT(I6,16F10.4)
      COVX2(I,K,L)=AA/SQRT(BB*CC)
      COVX2(I+NP,K,L)=AAA1/SQRT(BBB1*CCC1)
 6569 CONTINUE
 6568 CONTINUE
 6567 CONTINUE
!
      DO 656 K=1,NS
      DO 656 I=1,NP
      STDDEV(I,K)=SQRT(STDDEV(I,K)/FLOAT(NTRIAL-2))
      XMEANX(I,K)=XMEANX(I,K)/FLOAT(NTRIAL-1)
      STDDEV(I+NP,K)=SQRT(STDDEV(I+NP,K)/FLOAT(NTRIAL-2))
      XMEANX(I+NP,K)=XMEANX(I+NP,K)/FLOAT(NTRIAL-1)
  656 CONTINUE
      DO 657 I=1,NP
      IF(IPRINT.EQ.1)WRITE(26,150)I,                                    &
     &                (TRUEX(I,K),K=1,NS),                              &
     &                (XMEANX(I,K),K=1,NS),                             &
     &                (STDDEV(I,K),K=1,NS),                             &
     &                ((COVX2(I,K,L),L=1,NS),K=1,NS)
      KK=0
      DO 661 K=1,NS
      COVARIANCES(I+(K-1)*NP)=STDDEV(I,K)
      IF(K+1.LE.NS)THEN
         DO 658 L=K+1,NS
         COVARIANCES(I+NS*NP+KK*NP)=COVX2(I,K,L)
         KK=KK+1
  658    CONTINUE
      ENDIF
  661 CONTINUE
  657 CONTINUE
      ENDIF
!
      CLOSE(11)
      CLOSE(21)
      CLOSE(23)
      CLOSE(26)
      CLOSE(31)
      CLOSE(33)
      CLOSE(36)
!
!  SHUT DOWN RANDOM NUMBER GENERATOR
!
      CALL RNDEND()
      EXITSTATUS=1
        RETURN
        END
!
! ****************************************************************************
!  SUBROUTINE CLEAN----CALLED BY MAIN.  THROWS OUT ALL ROLLCALLS BELOW A
!     CUTOFF POINT (XVMIN)--USUALLY .05 OR .025.  THROWS OUT ALL LEGISLATORS
!     WHO HAVE VOTED KVMIN TIMES OR LESS (KVMIN=10 USUALLY).  KNAME(,)--
!     CHARACTER ARRAY HOLDING INFORMATION ON LEGISLATORS IS REARRANGED TO
!     HOLD ONLY INCLUDED LEGISLATORS.  LDATA( , ) IS THEN FILLED WITH THE
!     INCLUDED ROLL CALL VOTES.
! ****************************************************************************
!
      SUBROUTINE CLEAN(NP,NRCALL,XVMIN,KVMIN,IPRINT,KPTSUM,LDATA,       &
     &                 LDATA2,KAV,KAY,KAN)
      DIMENSION KMARG(50),KAV(NRCALL),KAY(NRCALL),                      &
     &          KAN(NRCALL),LDATA(NP,NRCALL),LDATA2(NP,NRCALL)
      CHARACTER*10 LMARG(10)
      INTEGER, ALLOCATABLE :: KKSUM(:)
      INTEGER, ALLOCATABLE :: LLSUM(:)
      INTEGER, ALLOCATABLE :: MMSUM(:)
      INTEGER, ALLOCATABLE :: NNSUM(:)
      REAL, ALLOCATABLE :: KD(:)     
      ALLOCATE ( KKSUM(NRCALL) )
      ALLOCATE ( LLSUM(NRCALL) )
      ALLOCATE ( MMSUM(NP) )
      ALLOCATE ( NNSUM(NP) )
      ALLOCATE ( KD(NRCALL) )      
!
      DO 34 I=1,NRCALL
      KKSUM(I)=0
      LLSUM(I)=0
  34  CONTINUE
      DO 35 I=1,50
      KMARG(I)=0
  35  CONTINUE
!
      LMARG(1)= '50 - 55   '
      LMARG(2)= '56 - 60   '
      LMARG(3)= '61 - 65   '
      LMARG(4)= '66 - 70   '
      LMARG(5)= '71 - 75   '
      LMARG(6)= '76 - 80   '
      LMARG(7)= '81 - 85   '
      LMARG(8)= '86 - 90   '
      LMARG(9)= '91 - 95   '
      LMARG(10)='96 - 99.5 '
      KPTSUM=0
      NAS=0
      NRS=0
      NOB=0
      NMOB=0
      DO 1 I=1,NP
      ISUM=0
      JSUM=0
      DO 2 J=1,NRCALL
      KD(J)=LDATA(I,J)
      IF(KD(J).EQ.1.OR.KD(J).EQ.2.OR.KD(J).EQ.3)KKSUM(J)=KKSUM(J)+1
      IF(KD(J).EQ.1.OR.KD(J).EQ.2.OR.KD(J).EQ.3)ISUM=ISUM+1
      IF(KD(J).EQ.4.OR.KD(J).EQ.5.OR.KD(J).EQ.6)LLSUM(J)=LLSUM(J)+1
  2   IF(KD(J).EQ.4.OR.KD(J).EQ.5.OR.KD(J).EQ.6)JSUM=JSUM+1
      MMSUM(I)=ISUM
  1   NNSUM(I)=JSUM
      DO 3 I=1,NP
      LL=MMSUM(I)+NNSUM(I)
      IF(LL.GT.KVMIN)NOB=NOB+1
      IF(LL.LE.KVMIN)NMOB=NMOB+1
      IF(LL.LE.KVMIN)GO TO 3
      NAS=0
      NRS=0
      DO 33 J=1,NRCALL
      KD(J)=LDATA(I,J)
      KK=KKSUM(J)+LLSUM(J)
      LL=MIN0(KKSUM(J),LLSUM(J))
      AA=0.0
      IF(KK.GT.0)AA=FLOAT(LL)/FLOAT(KK)
      IF(AA.GT.XVMIN)NAS=NAS+1
      IF(AA.LE.XVMIN)NRS=NRS+1
      IF(AA.GT.XVMIN)KAV(NAS)=J
      IF(AA.GT.XVMIN)KAY(NAS)=KKSUM(J)
      IF(AA.GT.XVMIN)KAN(NAS)=LLSUM(J)
      IF(AA.LE.XVMIN)GO TO 33
      IF(KD(J).EQ.2.OR.KD(J).EQ.3)KD(J)=1
      IF(KD(J).EQ.4.OR.KD(J).EQ.5)KD(J)=6
      IF(KD(J).EQ.7.OR.KD(J).EQ.8.OR.KD(J).EQ.9)KD(J)=0
      LDATA(NOB,NAS)=KD(J)
      LDATA2(NOB,NAS)=KD(J)
      IF(KD(J).NE.0)KPTSUM=KPTSUM+1
 33   CONTINUE
  3   CONTINUE
      IF(IPRINT.EQ.1)WRITE(23,1000)NRCALL,NRS,NAS,XVMIN
 1000 FORMAT(' ROLL-CALLS READ=',I4,2X,'NUMBER REJECTED=',I4,2X,        &
     &'NUMBER ACCEPTED=',I4,2X,'CUTOFF=',F6.3)
      IF(IPRINT.EQ.1)WRITE(23,1001)NP,NMOB,NOB,KVMIN
 1001 FORMAT(' LEGISLATORS READ=',I4,2X,'NUMBER REJECTED=',I4,2X,       &
     &'NUMBER ACCEPTED=',I4,2X,'CUTOFF=',I4)
      NRCALL=NAS
      NP=NOB
!
!   CALCULATE DISTRIBUTION OF ROLL CALL MARGINS
!
      DO 40 J=1,NRCALL
      KAY(J)=0
      KAN(J)=0
      DO 41 I=1,NP
      IF(LDATA(I,J).EQ.1)KAY(J)=KAY(J)+1
      IF(LDATA(I,J).EQ.6)KAN(J)=KAN(J)+1
  41  CONTINUE
  40  CONTINUE
      LLTOT=0
      LLALL=0
      DO 4 J=1,NRCALL
      LL=MAX0(KAY(J),KAN(J))
      LLTOT=LLTOT+LL
      AA=FLOAT(LL)/FLOAT(KAY(J)+KAN(J))
      LLALL=LLALL+KAY(J)+KAN(J)
      IF(AA.GE..50.AND.AA.LE..55)KMARG(1)=KMARG(1)+1
      IF(AA.GT..55.AND.AA.LE..60)KMARG(2)=KMARG(2)+1
      IF(AA.GT..60.AND.AA.LE..65)KMARG(3)=KMARG(3)+1
      IF(AA.GT..65.AND.AA.LE..70)KMARG(4)=KMARG(4)+1
      IF(AA.GT..70.AND.AA.LE..75)KMARG(5)=KMARG(5)+1
      IF(AA.GT..75.AND.AA.LE..80)KMARG(6)=KMARG(6)+1
      IF(AA.GT..80.AND.AA.LE..85)KMARG(7)=KMARG(7)+1
      IF(AA.GT..85.AND.AA.LE..90)KMARG(8)=KMARG(8)+1
      IF(AA.GT..90.AND.AA.LE..95)KMARG(9)=KMARG(9)+1
      IF(AA.GT..95.AND.AA.LE..995)KMARG(10)=KMARG(10)+1
  4   CONTINUE
      IF(IPRINT.EQ.1)WRITE(23,1003)
 1003 FORMAT('  DISTRIBUTION OF SCALABLE ROLL CALLS')
      DO 5 I=1,10
      AA=FLOAT(KMARG(I))/FLOAT(NRCALL)
      IF(IPRINT.EQ.1)WRITE(23,1002)I,LMARG(I),KMARG(I),AA
  5   CONTINUE
 1002 FORMAT(I4,1X,A10,I5,F7.3)
      XSUM=FLOAT(LLTOT)/FLOAT(LLALL)
      IF(IPRINT.EQ.1)WRITE(23,1004)LLTOT,LLALL,XSUM
 1004 FORMAT(' AVERAGE MAJORITY MARGIN= ',2I8,F9.5)
      DEALLOCATE (KKSUM)
      DEALLOCATE (LLSUM)
      DEALLOCATE (MMSUM)
      DEALLOCATE (NNSUM)
      DEALLOCATE (KD)   
      RETURN
      END
!
! **************************************************************************
!
!  SUBROUTINE KPASCORE -- PERFORMS EIGENVALUE-EIGENVECTOR DECOMPOSITION OF
!                       THE MATRIX OF LEGISLATOR BY LEGISLATOR AGREEMENT
!                       SCORES
!
! **************************************************************************
!
      SUBROUTINE KPASCORE(NP,NRCALL,NS,NDUAL,KIO,IPRINT,                &
     &                  ZMAT2,WVEC2,DSTAR,LDATA)
!
      DIMENSION ZMAT2(NP,NP),WVEC2(NP),                                 &
     &          DSTAR(NP,NP),LDATA(NP,NRCALL)
      REAL, ALLOCATABLE :: XCOL(:)
      REAL, ALLOCATABLE :: XROW(:)
      INTEGER, ALLOCATABLE :: KROW(:)
      REAL, ALLOCATABLE :: FV1(:)
      REAL, ALLOCATABLE :: FV2(:)
      REAL, ALLOCATABLE :: XAGREE(:,:)
      REAL, ALLOCATABLE :: ROWMEAN(:)
      REAL, ALLOCATABLE :: YCENTER(:,:)
      ALLOCATE ( XCOL(NRCALL) )
      ALLOCATE ( XROW(NP) )
      ALLOCATE ( KROW(NP) )
      ALLOCATE ( FV1(NP) )
      ALLOCATE ( FV2(NP) )
      ALLOCATE ( XAGREE(NP,NP) )
      ALLOCATE ( ROWMEAN(NP) )
      ALLOCATE ( YCENTER(NP,NP) )
!
!
  101 FORMAT(' PERFORMANCE INDEX EIGENVALUE/VECTOR ROUTINE='            &
     &       ,3I5,I6)
  102 FORMAT(I4,9F10.4)
 1000 FORMAT(3F10.4)
      DO 3 I=1,NRCALL
      XCOL(I)=0.0
  3   CONTINUE
      DO 30 I=1,NP
      XROW(I)=0.0
      KROW(I)=0
      ROWMEAN(I)=0.0
  30  CONTINUE
!
!  PERFORM THE PSEUDO-DOUBLE CENTER
!
!  CALCULATE THE ROLL CALL MEANS (#YEAS/#VOTING)
!
      XMAT=0.0
      KMAT=0
      XCHK1=0.0
      DO 1 J=1,NRCALL
      SUM=0.0
      KK=0
      DO 2 I=1,NP
      IF(LDATA(I,J).NE.0)THEN
         KK=KK+1
         KROW(I)=KROW(I)+1
      ENDIF
      IF(LDATA(I,J).EQ.1)THEN
         SUM=SUM+1.0
         XROW(I)=XROW(I)+1.0
      ENDIF
  2   CONTINUE
      XMAT=XMAT+SUM
      KMAT=KMAT+KK
      XCOL(J)=SUM/FLOAT(KK)
      XCHK1=XCHK1+XCOL(J)
  1   CONTINUE
      XCHK2=0.0
      DO 4 I=1,NP
      XROW(I)=XROW(I)/FLOAT(KROW(I))
      XCHK2=XCHK2+XROW(I)
  4   CONTINUE
!
!  MATRIX MEAN
!
      XCHK1=XCHK1/FLOAT(NRCALL)
      XCHK2=XCHK2/FLOAT(NP)
      XMAT=XMAT/FLOAT(KMAT)
      IF(IPRINT.EQ.1)WRITE(23,1000)XCHK1,XCHK2,XMAT
!
!  COMPUTE HECKMAN-SNYDER COVARIANCE MATRIX
!
      ALLMEAN=0.0
      KZERO=0
      DO 7 I=1,NP
      RSUM=0.0
      DO 8 J=1,NP
      SUM=0.0
      KK=0
      KKK=0
      DO 9 K=1,NRCALL
      IF(LDATA(I,K).EQ.0)GO TO 9
      IF(LDATA(J,K).EQ.0)GO TO 9
      KK=KK+1
!
!  SET UP SYMMETRIC MATRIX OF AGREEMENT SCORES
!
      IF(LDATA(I,K).EQ.LDATA(J,K))THEN
         KKK=KKK+1
      ENDIF
  9   CONTINUE
      IF(KK.EQ.0)THEN
         XAGREE(I,J)=0.25
         DSTAR(I,J)=1.0
         GO TO 88
      ENDIF
      XAGREE(I,J)=(1.0 - (FLOAT(KKK)/FLOAT(KK)))**2
      DSTAR(I,J)=(100.0 - (FLOAT(KKK)/FLOAT(KK))*100.0)/50.0
  88  CONTINUE
      RSUM=RSUM+XAGREE(I,J)
  8   CONTINUE
      ROWMEAN(I)=RSUM/FLOAT(NP)
      ALLMEAN=ALLMEAN+ROWMEAN(I)
  7   CONTINUE
      ALLMEAN=ALLMEAN/FLOAT(NP)
!
!  SETUP DOUBLE-CENTERED AGREEMENT SCORE MATRIX
!
      DO 33 I=1,NP
      DO 34 J=1,NP
      YCENTER(I,J)=(XAGREE(I,J)-ROWMEAN(I)-ROWMEAN(J)+ALLMEAN)/(-2.0)
  34  CONTINUE
  33  CONTINUE
!
!
!  EIGENVECTOR-EIGENVALUE DECOMPOSITION DOUBLE-CENTERED AGREEMENT
!     SCORE MATRIX
!
      CALL KPRS(NP,NP,YCENTER,WVEC2,1,ZMAT2,FV1,FV2,IER)
      IF(IPRINT.EQ.1)WRITE(23,101)NS,NP,IER,KZERO
      SUM2=0.0
      DO 20 I=1,NP
      SUM2=SUM2+ABS(WVEC2(I))
  20  CONTINUE
      YPER2=0.0
      IZULU=MIN(NP,20)
      DO 10 I=1,IZULU
      XPER2=(ABS(WVEC2(NP+1-I))/SUM2)*100.0
      YPER2=YPER2+XPER2
      IF(IPRINT.EQ.1)WRITE(23,102)I,WVEC2(NP+1-I),XPER2,YPER2
  10  CONTINUE
      DEALLOCATE (XCOL)
      DEALLOCATE (XROW)
      DEALLOCATE (KROW)
      DEALLOCATE (FV1)
      DEALLOCATE (FV2)
      DEALLOCATE (XAGREE)
      DEALLOCATE (ROWMEAN)
      DEALLOCATE (YCENTER)
      RETURN
      END
!
!
!  *************************************************************************
!   SUBROUTINE WHOOPE---IMPLEMENTS THE CONDITIONAL GLOBAL MINIMUM ALGORITHM
!   DEVELOPED BY POOLE, "LEAST SQUARES METRIC, UNIDIMENSIONAL UNFOLDING,"
!   PSYCHOMETRIKA, 1984.
!  *************************************************************************
!
!
      SUBROUTINE WHOOPE(NP,NS,DSTAR,ZZZ,XX,SSE1,SSE2,KTP,IPRINT)
      DIMENSION DSTAR(NP,NP),ZZZ(NP),DAT(20),XX(NP,25)
      REAL, ALLOCATABLE :: SAVEZ(:)
      REAL, ALLOCATABLE :: SAVED(:)
      REAL, ALLOCATABLE :: XXXX(:,:)
      ALLOCATE ( SAVEZ(NP) )
      ALLOCATE ( SAVED(NP) )
      ALLOCATE ( XXXX(NP,25) )
!
  100 FORMAT(I4,3F12.5,I8)
      KTP=1
      NPQ=NP-1
      CALL STATKP(NP,NS,DSTAR,ZZZ,XX,SSE1,RRSQ,KK)
      DAT(1)=SSE1
      II=0
      AKKK=0.0
      IF(IPRINT.EQ.1)WRITE(23,100)II,SSE1,RRSQ,AKKK,KK
      IF(SSE1.LE.0.001)SSE2=0.0
      IF(SSE1.LE.0.001)RETURN
      DO 99 II=1,10
      KTP=II
      DO 98 J=1,NP
      NPJ=J
      KK=0
      DO 918 JJ=1,NP
      IF(JJ.EQ.NPJ)GO TO 918
      KK=KK+1
      DO 919 K=1,NS
  919 XXXX(KK,K)=XX(JJ,K)
      SAVEZ(KK)=ZZZ(JJ)
      SAVED(KK)=DSTAR(NPJ,JJ)
  918 CONTINUE
      IF(NS.EQ.1)CALL FOCUSW(NP,NPQ,NPJ,SAVED,SAVEZ,ZZZ)
      IF(NS.GT.1)CALL FOCUS(NP,NPQ,NS,NPJ,SAVED,XX,XXXX)
  98  CONTINUE
      CALL STATKP(NP,NS,DSTAR,ZZZ,XX,SSE2,RRSQ,KK)
      DAT(II+1)=SSE2
      IF(SSE2.EQ.0.0)GO TO 9998
      AKKK=(DAT(II)-DAT(II+1))/DAT(II)
      IF(IPRINT.EQ.1)WRITE(23,100)II,SSE2,RRSQ,AKKK,KK
      IF(AKKK.LE..001)GO TO 9998
  99  CONTINUE
 9998 CONTINUE
      SUM=0.0
      DO 1 I=1,NP
  1   SUM=SUM+ZZZ(I)
      DO 2 I=1,NP
  2   ZZZ(I)=ZZZ(I)-(SUM/FLOAT(NP))
      DEALLOCATE (SAVEZ)
      DEALLOCATE (SAVED)
      DEALLOCATE (XXXX)
      RETURN
      END
!
!
!  *********************************************************************
!    SUBROUTINE FOCUSW---PERFORMS LEAST SQUARES METRIC SIMILARITIES
!    ANALYSIS USING THE CONDITIONAL GLOBAL MINIMUM ALGORITHM
!  *********************************************************************
!
!
      SUBROUTINE FOCUSW(NPT,NP,II,D,X,Z)
      DIMENSION D(NPT),X(NPT),Z(NPT)
      INTEGER, ALLOCATABLE :: LL(:)
      REAL, ALLOCATABLE :: Q(:)
      REAL, ALLOCATABLE :: XX(:,:)      
      ALLOCATE ( LL(NPT) )
      ALLOCATE ( Q(NPT) )
      ALLOCATE ( XX(NPT,2) )
!
      DO 303 I=1,NP
      LL(I)=I
  303 Q(I)=X(I)
      CALL RSORT(Q,NP,LL)
      ASUM=0.0
      BSUM=0.0
      WWSUM=0.0
      DO 66 I=1,NP
      IF(D(LL(I)).EQ.99.0)GO TO 66
      XX(I,1)=Q(I)-D(LL(I))
      XX(I,2)=Q(I)+D(LL(I))
      WWSUM=WWSUM+1.0
      ASUM=ASUM+XX(I,1)
      BSUM=BSUM+XX(I,1)**2
  66  CONTINUE
      AA=WWSUM*BSUM-ASUM*ASUM
      KK=1
      DO 77 I=1,NP
      IF(D(LL(I)).EQ.99.0)GO TO 77
      ASUM=ASUM-XX(I,1)+XX(I,2)
      BSUM=BSUM-XX(I,1)**2+XX(I,2)**2
      BB=WWSUM*BSUM-ASUM*ASUM
      CC=AMIN1(AA,BB)
      IF(ABS(CC-AA).LE..00001.AND.KK.GT.1)GO TO 88
      IF(ABS(CC-AA).LE..00001.AND.KK.EQ.1)Z(II)=                        &
     &                     (ASUM+XX(I,1)-XX(I,2))/WWSUM
      IF(ABS(CC-BB).LE..00001)Z(II)=ASUM/WWSUM
  88  AA=CC
      KK=KK+1
  77  CONTINUE
      DEALLOCATE (LL)
      DEALLOCATE (Q)
      DEALLOCATE (XX)
      RETURN
      END
!
!
!
!  **********************************************************************
!     SUBROUTINE STATKP--COMPUTES THE SUM OF SQUARED ERROR BETWEEN THE
!         THE INPUT DISTANCE MATRIX AND THE CURRENT DISTANCE MATRIX.
!  **********************************************************************
!
!
      SUBROUTINE STATKP(NP,NS,DSTAR,ZZZ,XX,SSE,RRSQ,KK)
      DIMENSION DSTAR(NP,NP),ZZZ(NP),XX(NP,25)
!
  100 FORMAT(2I5)
  125 FORMAT(I4,10F7.3)
      SSE=0.0
      ASUM=0.0
      BSUM=0.0
      CSUM=0.0
      DSUM=0.0
      ESUM=0.0
      KK=0
      DO 1 I=1,NP
!      WRITE(24,125)I,(XX(I,K),K=1,NS)
      DO 2 J=1,I
      IF(I.EQ.J)GO TO 2
      IF(DSTAR(I,J).EQ.99.0)GO TO 2
      KK=KK+1
      IF(NS.EQ.1)AA=ABS(ZZZ(I)-ZZZ(J))
      IF(NS.EQ.1)GO TO 10
      SSUMS=0.0
      DO 11 K=1,NS
  11  SSUMS=SSUMS+(XX(I,K)-XX(J,K))**2
      SSUMS=SQRT(SSUMS)
      AA=SSUMS
!      WRITE(24,125)KK,AA,DSTAR(I,J)
  10  CONTINUE
      ASUM=ASUM+AA
      BSUM=BSUM+DSTAR(I,J)
      CSUM=CSUM+AA*AA
      DSUM=DSUM+DSTAR(I,J)**2
      ESUM=ESUM+AA*DSTAR(I,J)
      SSE=SSE+(DSTAR(I,J)-AA)**2
  2   CONTINUE
  1   CONTINUE
      AA=FLOAT(KK)*ESUM-ASUM*BSUM
      BB=FLOAT(KK)*CSUM-ASUM*ASUM
      CC=FLOAT(KK)*DSUM-BSUM*BSUM
      RRSQ=(AA*AA)/(BB*CC)
      RETURN
      END
!
!
!  *********************************************************************
!    SUBROUTINE FOCUS IS A QUASI-GRADIENT ALGORITHM.  IT COMPUTES THE
!    COORDINATES
!  *********************************************************************
!
!
      SUBROUTINE FOCUS(NP,NPQ,NS,II,D,XX,XXXX)
      DIMENSION D(NP),XX(NP,25),XXXX(NP,25),ZZZ(100)
      DO 2 J=1,NS
      ZZZ(J)=0.0
  2   CONTINUE
      KK=0
      DO 4 J=1,NPQ
      IF(D(J).EQ.99.0)GO TO 4
      KK=KK+1
      SUM=0.0
      DO 7 K=1,NS
  7   SUM=SUM+(XXXX(J,K)-XX(II,K))**2
      IF(SUM.EQ.0.0E0) THEN 
      XC=1.0
      ELSE
      XC=D(J)/SQRT(SUM)
      ENDIF
  52  CONTINUE
      DO 8 K=1,NS
  8   ZZZ(K)=ZZZ(K)+XXXX(J,K)+XC*(XX(II,K)-XXXX(J,K))
  4   CONTINUE
      IF(KK.EQ.0)WRITE(23,310)II
      IF(KK.EQ.0)STOP
  310 FORMAT(' THIS IS YOUR PROBLEM STUPID!!!',I6)
      DO 1 K=1,NS
  1   XX(II,K)=ZZZ(K)/FLOAT(KK)
      RETURN
      END
!
!  **************************************************************************
!     SUBROUTINE RECODE--CALLED BY MAIN.  GENERATES STARTING VALUES FOR THE
!       ROLL CALL PARAMETERS.
!  **************************************************************************
!
      SUBROUTINE RECODE(NS,NP,NRCALL,KLSEN,KLSEN2,KIO,KSMIN,KSMAX,      &
     &                  KPTSUM,SSS,XDATA,ZMID,DYN,LDATA,LDATA2,XSAVE,   &
     &                  ZSAVE,CSAVE,KAV,KAY,KAN,IIII,IPRINT)
      DIMENSION SSS(100),ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),   &
     &          LDATA(NP,NRCALL),LDATA2(NP,NRCALL),                     &
     &          XSAVE(NP,2,2),ZSAVE(NRCALL,2,2),CSAVE(NRCALL,2),        &
     &          KAV(NRCALL),KAY(NRCALL),KAN(NRCALL)
      REAL, ALLOCATABLE :: X(:)
      INTEGER, ALLOCATABLE :: KV(:)
      INTEGER, ALLOCATABLE :: LL(:)
      ALLOCATE ( X(NP) )
      ALLOCATE ( KV(NP) )
      ALLOCATE ( LL(NP) )
!
 2000 FORMAT(' NUMBER LEGISLATORS=',I4,1X,'FURTHEST LEFT=',F10.4,1X,    &
     &'FURTHEST RIGHT=',F10.4/' POSITIONS: LEFT=',I4,2X,'RIGHT=',I4)
      AA=99.0
      BB=-99.0
!
!  SET FURTHEST RIGHT LEGISLATORS ON ALL DIMENSIONS TO POSITIVE VALUES 
!     IF NECESSARY AND FIND MAXIMUM AND MINIMUM LEGISLATORS
!
      XXX=XDATA(KLSEN,1)
      XXXX=XDATA(KLSEN2,2)
      IF(XXX.GT..0)GO TO 40
      DO 41 I=1,NP
  41  XDATA(I,1)=-XDATA(I,1)
  40  CONTINUE
      IF(XXXX.GT..0)GO TO 42
      DO 43 I=1,NP
  43  XDATA(I,2)=-XDATA(I,2)
  42  CONTINUE
      DO 1 I=1,NP
      AA=AMIN1(AA,XDATA(I,1))
      BB=AMAX1(BB,XDATA(I,1))
      IF(ABS(AA-XDATA(I,1)).LE..00001)KSMIN=I
      IF(ABS(BB-XDATA(I,1)).LE..00001)KSMAX=I
  1   CONTINUE
!
!  WRITE OUT MAXIMUM AND MINIMUM LEGISLATORS
!
      IF(IPRINT.EQ.1)WRITE(23,2000)NP,AA,BB,KSMIN,KSMAX
!
!  TRANSFORM LEGISLATOR CONFIGURATION SO THAT IT RUNS FROM -1 TO +1
!
      AA=(XDATA(KSMIN,1)+XDATA(KSMAX,1))/2.0
      BB=XDATA(KSMAX,1)-AA
      DO 10 I=1,NP
      LL(I)=I
      XDATA(I,1)=(XDATA(I,1)-AA)/BB
      XSAVE(I,1,1)=XDATA(I,1)
  10  X(I)=XDATA(I,1)
      CALL RSORT(X,NP,LL)
!
!  STORE 4 MOST EXTREME LEGISLATORS IN SSS(.)
      SSS(1)=X(1)
      SSS(2)=X(NP)
      SSS(3)=X(2)
      SSS(4)=X(NP-1)
!
!  CALL JANICE TO FIND OPTIMAL CUTTING LINE FOR EACH ROLL CALL AND
!    SET ZMID(.) TO THE OPTIMAL CUTTING LINE.  DETERMINE INITIAL
!    POLARITY OF ROLL CALL AND STORE -1.0 OR +1.0 IN DYN(.) DEPENDING
!    UPON POLARITY.
!
      KPTSUM=0
      DO 3 J=1,NRCALL
      SUMY=0.0
      SUMN=0.0
      DO 4 I=1,NP
      KV(I)=LDATA2(LL(I),J)
      IF(LDATA2(I,J).EQ.1)SUMY=SUMY+XDATA(I,1)
  4   IF(LDATA2(I,J).EQ.6)SUMN=SUMN+XDATA(I,1)
      SUMY=SUMY/FLOAT(KAY(J))
      SUMN=SUMN/FLOAT(KAN(J))
      SUMM=(SUMN+SUMY)/2.0
      KYES=KAY(J)
      KNO=KAN(J)
      JJJJ=J
!
!  JANICE FINDS OPTIMAL CUTTING POINT AND POLARITY
!
      CALL JANICE(NP,NRCALL,IIII,JJJJ,KAV,NP,KV,X,KYES,KNO,KIO,         &
     &            KFLIP,ZSAB,YPRE,SUMY,SUMN,SUMM,IPRINT)
      DO 7 I=1,NP
      IF(LDATA(I,J).NE.0)KPTSUM=KPTSUM+1
  7   IF(LDATA(I,J).EQ.6)LDATA(I,J)=2
      ZMID(J,1)=ZSAB
      DT=0.0
      IF(ZSAB.GT.0.0)DT=(ZSAB+1.0)/2.0
      IF(ZSAB.LE.0.0)DT=(1.0-ZSAB)/2.0
      IF(KFLIP.EQ.1)DYN(J,1)=DT
      IF(KFLIP.EQ.2)DYN(J,1)=-DT
      ZSAVE(J,1,1)=DYN(J,1)
      CSAVE(J,1)=ZMID(J,1)
  33  CONTINUE
  3   CONTINUE
      DEALLOCATE (X)
      DEALLOCATE (KV)
      DEALLOCATE (LL)
      RETURN
      END
!
!  ***************************************************************************
!    SUBROUTINE JANICE---CALLED BY RECODE.  SEARCHES EVERY PAIR OF ADJACENT
!      LEGISLATORS IN ASCENDING ORDER AND CALCULATES CLASSIFICATION ERROR
!      ASSOCIATED WITH USING THE MIDPOINT BETWEEN THE PAIR OF ADJACENT
!      LEGISLATORS AS THE ROLL CALL MIDPOINT.  THIS IS DONE FOR EACH POLARITY
!      OF THE ROLL CALL (LEFT = YES/RIGHT = YES).  PROGRAM RETURNS THE BEST
!      CUTPOINT AND POLARITY OF ROLL CALL.
!  ***************************************************************************
!
      SUBROUTINE JANICE(NP,NRCALL,IIII,JJJJ,KAV,NS,KV,X,KYES,KNO,JP,    &
     &                  KFLIP,ZSAB,YPRE,SUMY,SUMN,SUMM,IPRINT)
      DIMENSION KV(NP),X(NP),WK(10),KAV(NRCALL)
      REAL, ALLOCATABLE :: Y(:)
      REAL, ALLOCATABLE :: Z(:)
      INTEGER, ALLOCATABLE :: LV(:)
      INTEGER, ALLOCATABLE :: LVB(:)
      INTEGER, ALLOCATABLE :: LE(:)
      INTEGER, ALLOCATABLE :: LEB(:)
      ALLOCATE ( Y(NP))
      ALLOCATE ( Z(NP) )
      ALLOCATE ( LV(NP) )
      ALLOCATE ( LVB(NP) )
      ALLOCATE ( LE(NP) )
      ALLOCATE ( LEB(NP) )
  200 FORMAT('         KYES KNO   ZM   #ZM  CH  CL  EH  EL    %      PRE&
     &    YES    NO    MID')
  250 FORMAT(1X,4I4,F7.3,5I4,5F7.3)
  260 FORMAT(17X,F7.3,5I4,2F7.3)
!
      NSS=NS-1
      KCUT=1
      LCUT=6
      DO 999 III=1,2
      IF(III.EQ.2)KCUT=6
      IF(III.EQ.2)LCUT=1
!      IF(IIII.GT.1)KCUT=2
!      IF(IIII.GT.1)LCUT=1
      KMARK=1
      KSE=0
      KSV=0
      LSV=0
      LSE=0
      I=0
  10  I=I+1
      IF(I.EQ.NS)GO TO 12
      IF(I.LT.NS)Y(I)=(X(I)+X(I+1))/2.0
      IF(KV(I).EQ.KCUT)KSV=KSV+1
      IF(KV(I).EQ.LCUT)KSE=KSE+1
      IF(KMARK.EQ.1)GO TO 30
      IF(KV(I).EQ.KCUT)LSE=LSE-1
      IF(KV(I).EQ.LCUT)LSV=LSV-1
      GO TO 31
  30  KPDIS=I+1
      DO 3 J=KPDIS,NS
      IF(LCUT.EQ.KV(J)) THEN
      LSV=LSV+1
      ELSEIF(KCUT.EQ.KV(J)) THEN
      LSE=LSE+1
      ENDIF
  3   CONTINUE
      KMARK=0
  31  CONTINUE
      LV(I)=KSV
      LVB(I)=LSV
      LE(I)=KSE
      LEB(I)=LSE
      KT=LV(I)+LE(I)+LVB(I)+LEB(I)
      XKT=KT
      KTT=LE(I)+LEB(I)
      XKE=KTT
      Z(I)=(XKE/XKT)*100.0
      GO TO 10
  12  CONTINUE
      DO 64 J=1,5
      K=J
      AA=Z(J)
      AB=Y(J)
      LA=LV(J)
      LB=LE(J)
      LC=LVB(J)
       LD=LEB(J)
      DO 63 I=J,NSS
      IF(Z(I).GE.AA)GO TO 63
      K=I
      AA=Z(I)
      AB=Y(I)
      LA=LV(I)
      LB=LE(I)
      LC=LVB(I)
      LD=LEB(I)
  63  CONTINUE
      Z(K)=Z(J)
      Z(J)=AA
      Y(K)=Y(J)
      Y(J)=AB
      LV(K)=LV(J)
      LV(J)=LA
      LE(K)=LE(J)
      LE(J)=LB
      LVB(K)=LVB(J)
      LVB(J)=LC
      LEB(K)=LEB(J)
      LEB(J)=LD
  64  CONTINUE
      J=1
      DO 123 I=2,5
      IF(Z(1).EQ.Z(I)) J=J+1
  123 CONTINUE
      XMID=0.0
      DO 125 I=1,J
      XMID=XMID+Y(I)
      XFIT=FLOAT(LV(I))*FLOAT(LVB(I))
      YFIT=FLOAT(LE(I))*FLOAT(LEB(I))
      AFIT=FLOAT(LV(I))+FLOAT(LEB(I))
      BFIT=FLOAT(LVB(I))+FLOAT(LE(I))
      IF(AFIT.EQ.0.0.OR.BFIT.EQ.0.0)WK(I)=0.0
      IF(AFIT.EQ.0.0.OR.BFIT.EQ.0.0)GO TO 125
      WK(I)=(XFIT-YFIT)/(AFIT*BFIT)
  125 CONTINUE
      XMID=XMID/FLOAT(J)
      AA=Z(1)
      AB=Y(1)
      AC=WK(1)
      LA=LV(1)
      LB=LE(1)
      LC=LVB(1)
      LD=LEB(1)
      IF(J.EQ.1)GO TO 127
      DO 126 I=2,J
      IF(WK(I).LE.AC)GO TO 126
      AA=Z(I)
      AB=Y(I)
      AC=WK(I)
      LA=LV(I)
      LB=LE(I)
      LC=LVB(I)
      LD=LEB(I)
  126 CONTINUE
  127 JSAVE=J
      KPP=AMIN0(KYES,KNO)
      JKK=LB+LD
      XPRE=FLOAT(KPP-JKK)/FLOAT(KPP)
      NOWRIT=1
      IF(NOWRIT.EQ.1)GO TO 888
!      IF(IPRINT.EQ.0)GO TO 888
      IF(JJJJ.EQ.1.AND.III.EQ.1)WRITE(JP,200)
      IF(III.EQ.1)WRITE(JP,250)JJJJ,KAV(JJJJ),KYES,KNO,AB,JSAVE,LA,LC,  &
     &LB,LD,AA,XPRE,SUMY,SUMN,SUMM
      IF(III.EQ.2)WRITE(JP,260)AB,JSAVE,LA,LC,LB,LD,AA,XPRE
  888 IF(III.EQ.2)GO TO 999
      SAB=AB
      SXPRE=XPRE
  999 CONTINUE
      YPRE=AMAX1(SXPRE,XPRE)
      IF(ABS(YPRE-XPRE).LE..00001)THEN
         KFLIP=2
         ZSAB=AB
      ENDIF
      IF(ABS(YPRE-SXPRE).LE..00001)THEN
         KFLIP=1
         ZSAB=SAB
      ENDIF
      DEALLOCATE (Y)
      DEALLOCATE (Z)
      DEALLOCATE (LV)
      DEALLOCATE (LVB)
      DEALLOCATE (LE)
      DEALLOCATE (LEB)
      RETURN
      END
!
!
!
!  ************************************************************************
!    SORT SUBROUTINE--SORTS A VECTOR 'A' OF REAL ELEMENTS INTO ASCENDING
!    ORDER.  'LA' IS THE NUMBER OF ELEMENTS TO BE SORTED AND 'IR' IS A
!    VECTOR OF INTEGERS THAT RECORDS THE PERMUTATIONS--USUALLY SET TO
!    1,2,3,4,...
!  ************************************************************************
!
!
      SUBROUTINE RSORT(A,LA,IR)
      DIMENSION A(LA),IU(21),IL(21),IR(LA)
      IF (LA.LE.0) RETURN
      M = 1
      I = 1
      J = LA
      R = .375
    5 IF (I.EQ.J) GO TO 45
      IF (R.GT..5898437) GO TO 10
      R = R+3.90625E-2
      GO TO 15
   10 R = R-.21875
   15 K = I
!
! SELECT A CENTRAL ELEMENT OF THE  
! ARRAY AND SAVE IT IN LOCATION T  
!
      IJ = I+(J-I)*R 
      T = A(IJ)   
      IT = IR(IJ) 
!
! FIRST ELEMENT OF ARRAY IS GREATER
! THAN T, INTERCHANGE WITH T       
!
      IF (A(I).LE.T) GO TO 20  
      A(IJ) = A(I) 
      A(I) = T     
      T = A(IJ) 
      IR(IJ) = IR(I)  
      IR(I) = IT
      IT = IR(IJ) 
   20 L = J
!
! IF LAST ELEMENT OF ARRAY IS LESS THAN
! T, INTERCHANGE WITH T
!
      IF (A(J).GE.T) GO TO 30
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      IR(IJ) = IR(J)
      IR(J) = IT
      IT = IR(IJ)
!
! IF FIRST ELEMENT OF ARRAY IS GREATER
! THAN T, INTERCHANGE WITH T
!
      IF (A(I).LE.T) GO TO 30
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
      GO TO 30
   25 IF (A(L).EQ.A(K)) GO TO 30
      TT = A(L)
      A(L) = A(K)
      A(K) = TT
      ITT = IR(L)
      IR(L) = IR(K)
      IR(K) = ITT
!
! FIND AN ELEMENT IN THE SECOND HALF OF
! THE ARRAY WHICH IS SMALLER THAN T
!
   30 L = L-1
      IF (A(L).GT.T) GO TO 30
!
! FIND AN ELEMENT IN THE FIRST HALF OF
! THE ARRAY WHICH IS GREATER THAN T
!
   35 K = K+1
      IF (A(K).LT.T) GO TO 35
!
! INTERCHANGE THESE ELEMENTS
!
      IF (K.LE.L) GO TO 25
!
! SAVE UPPER AND LOWER SUBSCRIPTS OF
! THE ARRAY YET TO BE SORTED
!
      IF (L-I.LE.J-K) GO TO 40
      IL(M) = I
      IU(M) = L
      I = K
      M = M+1
      GO TO 50
   40 IL(M) = K
      IU(M) = J
      J = L
      M = M+1
      GO TO 50
!
! BEGIN AGAIN ON ANOTHER PORTION OF
! THE UNSORTED ARRAY
!
   45 M = M-1
      IF (M.EQ.0) RETURN
      I = IL(M)
      J = IU(M)
   50 IF (J-I.GE.11) GO TO 15
      IF (I.EQ.1) GO TO 5
      I = I-1
   55 I = I+1
      IF (I.EQ.J) GO TO 45
      T = A(I+1)
      IT = IR(I+1)
      IF (A(I).LE.T) GO TO 55
      K = I
   60 A(K+1) = A(K)
      IR(K+1) = IR(K)
      K = K-1
      IF (T.LT.A(K)) GO TO 60
      A(K+1) = T
      IR(K+1) = IT
      GO TO 55
      END
!
!
!  **************************************************************************
!    EIGENVECTOR/EIGENVALUE DECOMPOSITION SUBROUTINES FOR A SYMMETRIC MATRIX
!    SUBROUTINES ARE FROM EISPACK
!  **************************************************************************
!
      SUBROUTINE KPRS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
!
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
!
!     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
!     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
!     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
!     OF A REAL SYMMETRIC MATRIX.
!
!     ON INPUT
!
!        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
!        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!        DIMENSION STATEMENT.
!
!        N  IS THE ORDER OF THE MATRIX  A.
!
!        A  CONTAINS THE REAL SYMMETRIC MATRIX.
!
!        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
!        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
!        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
!
!     ON OUTPUT
!
!        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
!
!        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
!
!        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
!           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR KPTQLRAT
!           AND KPTQL2.  THE NORMAL COMPLETION CODE IS ZERO.
!
!        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
!
   10 IF (MATZ .NE. 0) GO TO 20
!     .......... FIND EIGENVALUES ONLY ..........
      CALL  KPTRED1(NM,N,A,W,FV1,FV2)
      CALL  KPTQLRAT(N,W,FV2,IERR)
      GO TO 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  KPTRED2(NM,N,A,W,FV1,Z)
      CALL  KPTQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
      SUBROUTINE KPTRED1(NM,N,A,D,E,E2)
!
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),E2(N)
      REAL F,G,H,SCALE
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
!     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
!          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(D(K))
!
         IF (SCALE .NE. 0.0E0) GO TO 140
!
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0E0
  125    CONTINUE
!
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 300
!
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
!
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0E0
!
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
!
  220       E(J) = G
  240    CONTINUE
!     .......... FORM P ..........
         F = 0.0E0
!
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
!
         H = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
!
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
!
  280    CONTINUE
!
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
!
  300 CONTINUE
!
      RETURN
      END
      SUBROUTINE KPTRED2(NM,N,A,D,E,Z)
!
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),Z(NM,N)
      REAL F,G,H,HH,SCALE
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
!     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!
!        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
!          PRODUCED IN THE REDUCTION.
!
!        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
!
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
!
         D(I) = A(N,I)
  100 CONTINUE
!
      IF (N .EQ. 1) GO TO 510
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 2) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(D(K))
!
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = D(L)
!
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0E0
            Z(J,I) = 0.0E0
  135    CONTINUE
!
         GO TO 290
!
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
!
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0E0
!
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
!
  220       E(J) = G
  240    CONTINUE
!     .......... FORM P ..........
         F = 0.0E0
!
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
!
         HH = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
!
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
!
            D(J) = Z(L,J)
            Z(I,J) = 0.0E0
  280    CONTINUE
!
  290    D(I) = H
  300 CONTINUE
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0E0
         H = D(I)
         IF (H .EQ. 0.0E0) GO TO 380
!
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
!
         DO 360 J = 1, L
            G = 0.0E0
!
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
!
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
!
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0E0
!
  500 CONTINUE
!
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0E0
  520 CONTINUE
!
      Z(N,N) = 1.0E0
      E(1) = 0.0E0
      RETURN
      END
      SUBROUTINE KPTQL2(NM,N,D,E,Z,IERR)
!
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,KPPYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  KPTRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  KPTRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1.
!
!        E HAS BEEN DESTROYED.
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS KPPYTHAG FOR  SQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      F = 0.0E0
      TST1 = 0.0E0
      E(N) = 0.0E0
!
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (TST1 .LT. H) TST1 = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = KPPYTHAG(P,1.0E0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
!
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = KPPYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + ABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE KPTQLRAT(N,D,E2,IERR)
!
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL D(N),E2(N)
      REAL B,C,F,G,H,P,R,S,T,KPEPSLON,KPPYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
!     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
!
!     ON INPUT
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
!          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!          THE SMALLEST EIGENVALUES.
!
!        E2 HAS BEEN DESTROYED.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS KPPYTHAG FOR  SQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
!
      F = 0.0E0
      T = 0.0E0
      E2(N) = 0.0E0
!
      DO 290 L = 1, N
         J = 0
         H = ABS(D(L)) + SQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = KPEPSLON(T)
         C = B * B
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * S)
         R = KPPYTHAG(P,1.0E0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
!
         DO 140 I = L1, N
  140    D(I) = D(I) - H
!
         F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0E0) G = B
         H = G
         S = 0.0E0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0E0) G = B
            H = G * P / R
  200    CONTINUE
!
         E2(L) = S * G
         D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0E0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0E0) GO TO 130
  210    P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
!
  250    I = 1
  270    D(I) = P
  290 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      REAL FUNCTION KPPYTHAG(A,B)
      REAL A,B
!
!     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
!
      REAL P,R,S,T,U
      P = AMAX1(ABS(A),ABS(B))
      IF (P .EQ. 0.0E0) GO TO 20
      R = (AMIN1(ABS(A),ABS(B))/P)**2
   10 CONTINUE
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         U = 1.0E0 + 2.0E0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 KPPYTHAG = P
      RETURN
      END
!
      REAL FUNCTION KPEPSLON (X)
      REAL X
      REAL A,B,C,EPS
      A = 4.0E0/3.0E0
  10  B = A - 1.0E0
      C = B + B + B
      EPS = ABS(C-1.0E0)
      IF(EPS .EQ. 0.0E0)GO TO 10
      KPEPSLON = EPS*ABS(X)
      RETURN
      END
!
!  **************************************************************************
!    SUBROUTINE MAXLNL---CALLED BY MAIN.  INITIALIZES VARIABLES FOR
!          ESTIMATION OF PARAMETERS
!  **************************************************************************
!
      SUBROUTINE MAXLNL(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,JANLST, &
     &                  NITR,KSMIN,KSMAX,KPTSUM,                        &
     &                  BBB,BBBB,SSS,ZMID,XDATA,DYN,LDATA,PSI,          &
     &                  YBIGL,YYBIGL,XFITS,                             &
     &                  XSAVE,ZSAVE,CSAVE,STDDVX,STDDVZ,                &
     &                  KPJP,LERIC,ISENS,IPRINT)                        &
     &                  
      DIMENSION LL(2,2),MM(2,2),BBB(50),BBBB(25),SSS(100),              &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),                      &
     &          KPJP(NDUAL,4),LERIC(NRCALL),ISENS(NP),                  &
     &          B(50),TESTB(50),DELTAB(50),GRAD(50),V(25,25),           &
     &          YBIGL(NDUAL),YYBIGL(NDUAL),XFITS(3*NS),                 &
     &          XSAVE(NP,2,2),ZSAVE(NRCALL,2,2),CSAVE(NRCALL,2),        &
     &          STDDVX(NP,25),STDDVZ(NRCALL,2,25)
       EXTERNAL FUNNEL
!
  23  FORMAT(' TOTAL VOTES PREDICTED FOR',I5,1X,'ROLL-CALLS'/13X,       &
     &'ACTUAL'/7X,'  YEA     NAY  ')
  24  FORMAT(' TOTAL VOTES PREDICTED FOR',I5,1X,'LEGISLATORS'/13X,      &
     &'ACTUAL'/7X,'  YEA     NAY  ')
  25  FORMAT('   YEA  ',I6,2X,I6/'   NAY  ',I6,2X,I6)
 2626 FORMAT(' PERCENT CLASSIFIED=',F7.3,'    AGGREGATE PRE=',F7.3)
!
!  ***** UTILITY FUNCTION PHASE *****
!    ESTIMATION OF BETA ON FIRST DIMENSION ONLY
!
      IF(NSTEP.EQ.0)THEN
           B(1)=BBB(1)
!
!  CALL SUBROUTINE NORMZ TO NORMALIZE COORDINATES AND CHECK FOR SAG
!    IN LEGISLATOR COORDINATES
!
           CALL NORMZ(NP,NRCALL,NS,NDIM,                                &
     &                  KSMIN,KSMAX,                                    &
     &                  BBB,BBBB,SSS,ZMID,XDATA,DYN,                    &
     &                  ISENS,IPRINT)
!
           NFEVAL=0
           CALL BHHH(NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,       &
     &                  JANLST,SVLNL,FNLNL,ICNVG,ITT,NVAL,              &
     &                  XLNL,B,TESTB,DELTAB,GRAD,V,                     &
     &                  ZMID,XDATA,DYN,LDATA,PSI,                       &
     &                  YBIGL,YYBIGL,ISENS,BBB,BBBB,IPRINT)
           CALL RPRINT(KPTSUM,NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,  &
     &                  JAN,JANLST,NITR,KSMIN,KSMAX,                    &
     &                  XLNL,B,V,SVLNL,FNLNL,ICNVG,ITT,NVAL,BBB,BBBB,   &
     &                  SSS,ZMID,XDATA,DYN,LDATA,PSI,XSAVE,ZSAVE,CSAVE, &
     &                  STDDVX,STDDVZ,YBIGL,YYBIGL,IPRINT)
      ENDIF
!
!  ***** ROLL CALL PHASE *****
!
      IF(NSTEP.EQ.1)THEN
           DO 22 I=1,2
           DO 22 J=1,2
  22       MM(I,J)=0
           NEQ=0
           JKP=0
           KTP=0
  1        CONTINUE
           NEQ=NEQ+1
           NFEVAL=0
           B(1)=DYN(NEQ,NDIM)
           B(2)=ZMID(NEQ,NDIM)
           IF(NDIM.GT.1)THEN
              B(1)=0.0
              B(2)=0.0
           ENDIF
           IF(LERIC(NEQ).NE.0)THEN
             CALL BHHH(NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,     &
     &                  JANLST,SVLNL,FNLNL,ICNVG,ITT,NVAL,              &
     &                  XLNL,B,TESTB,DELTAB,GRAD,V,                     &
     &                  ZMID,XDATA,DYN,LDATA,PSI,                       &
     &                  YBIGL,YYBIGL,ISENS,BBB,BBBB,IPRINT)
           ENDIF
             CALL RPRINT(KPTSUM,NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,&
     &                  JAN,JANLST,NITR,KSMIN,KSMAX,                    &
     &                  XLNL,B,V,SVLNL,FNLNL,ICNVG,ITT,NVAL,BBB,BBBB,   &
     &                  SSS,ZMID,XDATA,DYN,LDATA,PSI,XSAVE,ZSAVE,CSAVE, &
     &                  STDDVX,STDDVZ,YBIGL,YYBIGL,IPRINT)
             CALL CROSS(NEQ,LL,                                         &
     &                 NP,NRCALL,NDIM,NSTEP,JAN,JANLST,                 &
     &                 ZMID,XDATA,DYN,LDATA,PSI,BBB,BBBB)
           JJ=1
           DO 21 I=1,2
           DO 21 J=1,2
           KPJP(NEQ,JJ)=LL(I,J)
           JJ=JJ+1
  21       MM(I,J)=MM(I,J)+LL(I,J)
           LERIC(NEQ)=1
!           IF((LL(1,2)+LL(2,1)).EQ.0)LERIC(NEQ)=0
           LYES=LL(1,1)+LL(2,1)
           LNO =LL(1,2)+LL(2,2)
           JKP=JKP+AMIN0(LYES,LNO)
           KTP=KTP+LL(1,2)+LL(2,1)
           IF(NEQ.LT.NRCALL)GO TO 1
!
!  CALCULATE PERCENT CORRECTLY CLASSIFIED AND OVERALL PRE 
!
           XPER=FLOAT(MM(1,1)+MM(2,2))/                                 &
     &              FLOAT(MM(1,1)+MM(1,2)+MM(2,1)+MM(2,2))
           XPER=XPER*100.0
           XPRE=FLOAT(JKP-KTP)/FLOAT(JKP)
           XFITS(NDIM)=XPER
           XFITS(NS+NDIM)=XPRE
           IF(IPRINT.EQ.1)THEN
              WRITE(21,2626)XPER,XPRE
              WRITE(23,23)NRCALL
              WRITE(23,25)((MM(I,J),J=1,2),I=1,2)
           ENDIF
      ENDIF
!
!  ***** LEGISLATOR PHASE *****
!
      IF(NSTEP.EQ.2)THEN
           DO 26 I=1,2
           DO 26 J=1,2
  26       MM(I,J)=0
           NEQ=0
  10       CONTINUE
           NEQ=NEQ+1
           NFEVAL=0
           B(1)=XDATA(NEQ,NDIM)
           CALL BHHH(NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,       &
     &                  JANLST,SVLNL,FNLNL,ICNVG,ITT,NVAL,              &
     &                  XLNL,B,TESTB,DELTAB,GRAD,V,                     &
     &                  ZMID,XDATA,DYN,LDATA,PSI,                       &
     &                  YBIGL,YYBIGL,ISENS,BBB,BBBB,IPRINT)
           CALL RPRINT(KPTSUM,NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,  &
     &                  JAN,JANLST,NITR,KSMIN,KSMAX,                    &
     &                  XLNL,B,V,SVLNL,FNLNL,ICNVG,ITT,NVAL,BBB,BBBB,   &
     &                  SSS,ZMID,XDATA,DYN,LDATA,PSI,XSAVE,ZSAVE,CSAVE, &
     &                  STDDVX,STDDVZ,YBIGL,YYBIGL,IPRINT)
           CALL CROSS(NEQ,LL,                                           &
     &                 NP,NRCALL,NDIM,NSTEP,JAN,JANLST,                 &
     &                 ZMID,XDATA,DYN,LDATA,PSI,BBB,BBBB)
           JJ=1
           DO 27 I=1,2
           DO 27 J=1,2
           KPJP(NEQ,JJ)=LL(I,J)
           JJ=JJ+1
  27       MM(I,J)=MM(I,J)+LL(I,J)
           IF(NEQ.LT.NP)GO TO 10
           IF(IPRINT.EQ.1)THEN
              WRITE(23,24)NP
              WRITE(23,25)((MM(I,J),J=1,2),I=1,2)
           ENDIF
      ENDIF
!
!  ***** UTILITY FUNCTION PHASE *****
!     ESTIMATION OF WEIGHT FOR DIMENSIONS 2 AND ABOVE
!
      IF(NSTEP.EQ.-1)THEN
           B(1)=BBB(2)
           NFEVAL=0
           CALL BHHH(NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,       &
     &                  JANLST,SVLNL,FNLNL,ICNVG,ITT,NVAL,              &
     &                  XLNL,B,TESTB,DELTAB,GRAD,V,                     &
     &                  ZMID,XDATA,DYN,LDATA,PSI,                       &
     &                  YBIGL,YYBIGL,ISENS,BBB,BBBB,IPRINT)
           CALL RPRINT(KPTSUM,NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,  &
     &                  JAN,JANLST,NITR,KSMIN,KSMAX,                    &
     &                  XLNL,B,V,SVLNL,FNLNL,ICNVG,ITT,NVAL,BBB,BBBB,   &
     &                  SSS,ZMID,XDATA,DYN,LDATA,PSI,XSAVE,ZSAVE,CSAVE, &
     &                  STDDVX,STDDVZ,YBIGL,YYBIGL,IPRINT)
      ENDIF
      RETURN
      END
!
!  **************************************************************************
!    SUBROUTINE BHHH--CALLED BY MAXLIK.  IMPLEMENTS THE BHHH ALGORITHM FOR
!       MAXIMIZING LOG-LIKELIHOOD
!  **************************************************************************
!
      SUBROUTINE BHHH(NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR,JAN,      &
     &                JANLST,SVLNL,FNLNL,ICNVG,ITT,NVAL,                &
     &                XLNL,B,TESTB,DELTAB,GRAD,V,                       &
     &                  ZMID,XDATA,DYN,LDATA,PSI,                       &
     &                  YBIGL,YYBIGL,ISENS,BBB,BBBB,IPRINT)
      DIMENSION AV(25,25),KCONST(20),ZMAT(25,25),WVEC(25),FV1(25),      &
     &          FV2(25),B(50),TESTB(50),DELTAB(50),GRAD(50),V(25,25),   &
     &          BBB(50),BBBB(25),ZMID(NRCALL,25),                       &
     &          XDATA(NP,25),DYN(NRCALL,25),                            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),                      &
     &          YBIGL(NDUAL),YYBIGL(NDUAL),ISENS(NP)
!
      EXTERNAL FUNNEL
  101 FORMAT(' PERFORMANCE INDEX EIGENVECTOR/VALUE ROUTINE=',I4)
!
      NFEVAL=0
      TOLB=.4
      MAXIT=30
      MAXSQZ=10
      ICNVG=0
      DO 50 I=1,NOPAR
      KCONST(I)=1
      DELTAB(I)=0.0
   50 TESTB(I)=B(I)
!
!          BERNDT, HALL, HALL, AND HAUSMAN METHOD
!
      DO 1000 IT=1,MAXIT
      KFDRV=1
!
!
      XLNL=FUNNEL(0.0,                                                  &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!      XLNL=FUNNEL(0.0)
!
      IF(IT.EQ.1)SVLNL=-XLNL
      KFDRV=0
      IF(NOPAR.EQ.1)GO TO 20
      IF(NSTEP.EQ.0) KCONST(1)=1
      IF(NSTEP.EQ.1) KCONST(1)=1
      IF(NSTEP.EQ.1) KCONST(2)=1
      IF(NSTEP.EQ.1) KCONST(3)=0
      MSUM=KCONST(1)+KCONST(2)+KCONST(3)
      IF(MSUM.EQ.0)GO TO 1200
      IF(MSUM.EQ.1)GO TO 21
      IF(MSUM.EQ.2)GO TO 23
!
!         V11 V12 V13
!     V=  V21 V22 V23
!         V31 V32 V33
!
!     AV = INVERSE-V
!
!  CALCULATE MATRIX INVERSE
!
      CALL KPRS(25,3,V,WVEC,1,ZMAT,FV1,FV2,IER)
      DO 60 I=1,3
      DO 60 KX=1,3
      SUM=0.0
      DO 61 J=1,3
  61  SUM=SUM+ZMAT(KX,J)*(1.0/WVEC(J))*ZMAT(I,J)
  60  AV(I,KX)=SUM
      DET=V(1,1)*V(2,2)*V(3,3)-V(1,1)*V(2,3)*V(3,2)+V(1,2)*V(2,3)*V(3,1)
      DET=DET-V(1,2)*V(2,1)*V(3,3)+V(1,3)*V(2,1)*V(3,2)-V(1,3)*V(2,2)*  &
     &V(3,1)
      IF(ABS(DET).LE..00000001)ICNVG=9
      IF(ABS(DET).LE..00000001)RETURN
      DELTAB(1)=(AV(1,1)*GRAD(1)+AV(2,1)*GRAD(2)+AV(3,1)*GRAD(3))
      DELTAB(2)=(AV(1,2)*GRAD(1)+AV(2,2)*GRAD(2)+AV(3,2)*GRAD(3))
      DELTAB(3)=(AV(1,3)*GRAD(1)+AV(2,3)*GRAD(2)+AV(3,3)*GRAD(3))
      DO 77 I=1,3
      DO 77 J=1,3
      V(I,J)=AV(J,I)
77    CONTINUE
      GO TO 22
!
!     WHEN ONE OF THE THREE KCONST'S IS 1.
!
21    IF(KCONST(1).EQ.1.AND.KCONST(2).EQ.0.AND.KCONST(3).EQ.0)GO TO 30
      IF(KCONST(1).EQ.0.AND.KCONST(2).EQ.1.AND.KCONST(3).EQ.0)GO TO 31
      IF(KCONST(1).EQ.0.AND.KCONST(2).EQ.0.AND.KCONST(3).EQ.1)GO TO 32
30    IC=1
      JC=2
      KC=3
      GO TO 35
31    IC=2
      JC=1
      KC=3
      GO TO 35
32    IC=3
      JC=1
      KC=2
35    IF(ABS(V(IC,IC)).LE..00000001)ICNVG=9
      IF(ABS(V(IC,IC)).LE..00000001)RETURN
      DELTAB(IC)=GRAD(IC)/V(IC,IC)
      V(IC,IC)=1.0/V(IC,IC)
      DELTAB(JC)=0.0
      DELTAB(KC)=0.0
      GO TO 22
!
!     WHEN TWO OF THE THREE KCONST'S ARE 1'S.
!
23    IF(KCONST(1).EQ.1.AND.KCONST(2).EQ.1.AND.KCONST(3).EQ.0)GO TO 40
      IF(KCONST(1).EQ.1.AND.KCONST(2).EQ.0.AND.KCONST(3).EQ.1)GO TO 41
      IF(KCONST(1).EQ.0.AND.KCONST(2).EQ.1.AND.KCONST(3).EQ.1)GO TO 42
40    IC=1
      JC=2
      KC=3
      GO TO 45
41    IC=1
      JC=3
      KC=2
      GO TO 45
42    IC=2
      JC=3
      KC=1
45    DET=V(IC,IC)*V(JC,JC)-V(IC,JC)*V(JC,IC)
      IF(ABS(DET).LE..00000001)ICNVG=9
      IF(ABS(DET).LE..00000001)RETURN
      DELTAB(IC)=(V(JC,JC)*GRAD(IC)-V(IC,JC)*GRAD(JC))/DET
      DELTAB(JC)=(V(IC,IC)*GRAD(JC)-V(JC,IC)*GRAD(IC))/DET
      DELTAB(KC)=0.0
      S11=V(IC,IC)
      V(IC,IC)=V(JC,JC)/DET
      V(JC,JC)=S11/DET
      V(IC,JC)=-V(IC,JC)/DET
      V(JC,IC)=-V(JC,IC)/DET
      GO TO 22
  20  IF(ABS(V(1,1)).LE..00000001)ICNVG=9
      IF(ABS(V(1,1)).LE..00000001)RETURN
      DELTAB(1)=GRAD(1)/V(1,1)
      V(1,1)=1.0/V(1,1)
      DELTAB(2)=0.0
      DELTAB(3)=0.0
      GO TO 22
  22  CONTINUE
      COSDG=DELTAB(1)*GRAD(1)+DELTAB(2)*GRAD(2)+DELTAB(3)*              &
     &      GRAD(3)
!
!  AFTER 3 ITERATIONS STEPSIZE DEFAULTS TO PREVIOUS STEPSIZE
!
      IF(IT.LE.3)ALAMB=1.0
      IFBETR=1
      CALL STEPR(XLNL,XLNLNW,ALAMB,IFBETR,MAXSQZ,ISQZ,                  &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!
!          PRINT ITERATION RESULTS
!
      IF(NSTEP.EQ.0.AND.IPRINT.EQ.1)WRITE(23,1918)IT,                   &
     &                              (TESTB(J),J=1,NOPAR),               &
     &                              ALAMB,COSDG,XLNL,XLNLNW
 1918 FORMAT(I5,12F13.4/5X,12F13.4)
!
!          CHECK FOR CONVERGENCE OF LIKELIHOOD FUNCTION.
!
!
      IF(COSDG.LT.(TOLB**2))GO TO 1200
      IF(IFBETR.EQ.0)GO TO 1500
      DO 710 I=1,NOPAR
  710 B(I)=TESTB(I)
 1000 CONTINUE
      GO TO 1100
!
!          FAILURE TO CONVERGE
!
 1100 ICNVG=1
      ITT=IT
      GO TO 1206
!
!          CONVERGENCE
!
 1200 ICNVG=0
      ITT=IT
 1206 CONTINUE
 1207 DO 705 I=1,NOPAR
  705 B(I)=TESTB(I)
      XLNL=-FUNNEL(0.0,                                                 &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!      XLNL=-FUNNEL(0.0)
      FNLNL=XLNL
      ITT=IT
      NVAL=NFEVAL
      RETURN
!
!          FAILURE TO IMPROVE LNL
!
 1500 ICNVG=2
      DO 650 I=1,NOPAR
  650 TESTB(I)=B(I)
      GO TO 1207
      END
!
!  ***************************************************************************
!    SUBROUTINE STEPR---CALLED BY BHHH.  CALCULATES STEPSIZE FOR GRADIENT
!  ***************************************************************************
!
      SUBROUTINE STEPR(F,F0,R,IFOK,MAXSQZ,NSQZ,                         &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
      DIMENSION B(50),TESTB(50),DELTAB(50),GRAD(50),V(25,25),ISENS(NP), &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),BBB(50),              &
     &          YBIGL(NDUAL),YYBIGL(NDUAL)
!
!
      EXTERNAL FUNNEL
!
!
!        IF STEPPING FORWARD REDUCES -LNL, DOUBLE THE STEP.
!        IF NOT, HALVE THE STEP.  ONCE STEPPING TAKES PLACE
!        IN ONE DIRECTION, IT CONTINUES IN THAT DIRECTION
!        UNTIL NO FURTHER IMPROVEMENT OCCURS.
!
!
!
!  THIS SECTION INTIALIZES LOG-LIKELIHOOD--IF IT IS BETTER THAN OLD
!    VALUE THE STEP SIZE IS DOUBLED IN CODE SEGMENT BELOW.  IF IT IS
!    WORSE, THE STEP SIZE IN HALVED IN CODE SEGMENT IMMEDIATELY BELOW
!
      IFB = IFOK
      R0 = R
      NSQZ = 0
      IFOK=1
  2   CONTINUE
      F0 = FUNNEL(R0,                                                   &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!    2 F0 = FUNNEL(R0)
      IF(F0.LT.F.AND.IFB.EQ.1)GO TO 10
  99  R0=R0/2.
      NSQZ = NSQZ +1
      IF(NSQZ .LE. MAXSQZ) GO TO 2
      IFOK=0
      R = R0
      RETURN
   10 IF(NSQZ .GT. 0) GO TO 13
!
!  IF IMPROVEMENT, STEP SIZE IS DOUBLED
!
   14 R1 = R0 * 2.
      F1 = FUNNEL(R1,                                                   &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!      F1 = FUNNEL(R1)
      NSQZ = NSQZ +1
      IF(F1 .GT. F0) GO TO 11
      IF(F1.LT.0.0)GO TO 11
      R0 = R1
      F0 = F1
      IF(NSQZ .LT. MAXSQZ) GO TO 14
      R = R1
      F0 = F1
      RETURN
!
  11  CONTINUE
      F0 = FUNNEL(R0,                                                   &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!   11 F0 = FUNNEL(R0)
   13 R = R0
      RETURN
      END
!
!  **************************************************************************
!    FUNCTION FUNNEL---CALLED BY STEPR AND BY BHHH WHEN INITIALIZING.  IT
!      CALCULATES THE UPDATED PARAMETER VECTOR AS NOTED BELOW AND CALLS
!      SUBROUTINE LOGLIK WHICH CALCULATES THE CORRESPONDING LOG-LIKELIHOOD
!  **************************************************************************
!
      REAL FUNCTION FUNNEL(XLBA,                                        &
     &                     NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS, &
     &                  NOPAR,KFDRV,NFEVAL,DELTAB,B,                    &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
      DIMENSION B(50),TESTB(50),DELTAB(50),GRAD(50),V(25,25),ISENS(NP), &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),BBB(50),    &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),                      &
     &          YBIGL(NDUAL),YYBIGL(NDUAL)
!
!          UPDATE TRIAL PARAMETER VECTOR.
!
!        TESTB = B + XLBA * DELTAB
!
!           WHERE DELTAB = QMAT(INVERSE) * G
!           AND XLBA = STEPSIZE
!
      DO 10 I=1,NOPAR
  10  TESTB(I)=B(I)+XLBA*DELTAB(I)
      CALL LOGLIK(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS,          &
     &                  NOPAR,KFDRV,NFEVAL,                             &
     &                  TESTB,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,    &
     &                  YBIGL,YYBIGL)
!      CALL LOGLIK(TESTB,FLIKE,GRAD,V)
      FUNNEL=FLIKE
      RETURN
      END
!
!  ***************************************************************************
!    SUBROUTINE LOGLIK---CALLED BY FUNNEL.  CALCULATES LOG-LIKELIHOOD.  FOR
!       EACH APPROPRIATE OBSERVATION, IT CALLS SUBROUTINE ITHOBS (ITH
!       OBSERVATION) WHICH RETURNS THE LOG-LIKELIHOOD.  THE SUM IS CUMULATED
!       BELOW.
!  ***************************************************************************
!
      SUBROUTINE LOGLIK(NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NEQ,BBB,ISENS,    &
     &                  NOPAR,KFDRV,NFEVAL,                             &
     &                  B,FLIKE,GRAD,V,ZMID,XDATA,DYN,LDATA,PSI,        &
     &                  YBIGL,YYBIGL)
      DIMENSION B(50),GRAD(50),V(25,25),G(50),XBIGL(NRCALL),ISENS(NP),  &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),BBB(50),              &
     &          YBIGL(NDUAL),YYBIGL(NDUAL)
!
      IF(NSTEP.EQ.-1)NNBS=NP
      IF(NSTEP.EQ.0) NNBS=NP
      IF(NSTEP.EQ.1) NNBS=NP
      IF(NSTEP.EQ.2) NNBS=NRCALL
      NFEVAL=NFEVAL+1
      FLIKE=0.0
      IF(NSTEP.EQ.2)YYBIGL(NEQ)=0.0
      IF(NSTEP.NE.1)GO TO 3
      DO 2 J=NEQ,NRCALL
  2   YBIGL(J)=0.0
  3   CONTINUE
      IF(KFDRV.EQ.0)GO TO 202
      DO 201 J=1,NOPAR
      GRAD(J)=0.0
  201 CONTINUE
      DO 101 J=1,NOPAR
      DO 101 I=1,NOPAR
  101 V(I,J)=0.0
  202 CONTINUE
!      IF(NSTEP.LE.0)GO TO 301
!
!          LOOP OVER NUMBER OF OBSERVATIONS NNBS AND SUM LOG LIKELIHOOD
!     AND ITS GRADIENT WITH RESPECT TO B.
!
      I=0
  200 I=I+1
      IF(I.GT.NNBS)GO TO 211
      CALL ITHOBS(NP,NRCALL,NS,NDIM,NSTEP,NEQ,BBB,ISENS,                &
     &                  I,B,XLNL,G,XBIGL,ZMID,XDATA,DYN,LDATA,PSI)
!      CALL ITHOBS(I,B,XLNL,G,XBIGL)
      IF(NSTEP.EQ.2)YYBIGL(NEQ)=YYBIGL(NEQ)+XLNL
      IF(NSTEP.EQ.1)YBIGL(NEQ)=YBIGL(NEQ)+XLNL
      XLNL=-XLNL
      FLIKE=FLIKE+XLNL
      IF(KFDRV.EQ.0)GO TO 200
      DO 205 J=1,NOPAR
  205 GRAD(J)=GRAD(J)+G(J)
!
!          EVALUATE QMAT, THE APPROXIMATION TO THE HESSIAN, AS THE
!     OUTER PRODUCT OF THE GRADIENT VECTOR ONLY FOR THE BHHH METHOD.
!
      DO 206 L=1,NOPAR
      DO 206 J=1,NOPAR
  206 V(J,L)=V(J,L)+G(J)*G(L)
      GO TO 200
  211 RETURN
      END
!
!  ***************************************************************************
!    SUBROUTINE ITHOBS--CALLED BY LOGLIK.  CALCULATES THE LOG-LIKELIHOOD FOR
!      EACH OBSERVATION (I.E. THE LOG-LIKELIHOOD FOR ONE COLUMN/ROW IN THE
!      INPUT CHOICE MATRIX.
!  ***************************************************************************
!
      SUBROUTINE ITHOBS(NP,NRCALL,NS,NDIM,NSTEP,NEQ,BBB,ISENS,          &
     &                  JT,B,XL,G,XBIGL,ZMID,XDATA,DYN,LDATA,PSI)
      DIMENSION DSQ(4,NRCALL),XB(4,NRCALL),                             &
     &          EXB(4,NRCALL),DLDX(4),B(50),G(50),                      &
     &          PHI(NRCALL),XBIGL(NRCALL),BBB(50),ISENS(NP),            &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2)
!
!   SUBROUTINE TO ESTIMATE MULTINOMIAL LOGIT FOR SPATIAL MODELS
!
!   NP=NUMBER OF INDIVIDUALS=NP
!   NS=NUMBER OF DIMENSIONS
!   NQ=NUMBER OF STIMULI
!
      NQ=2
!
!
!   B (PARAMETER VECTOR IN MAXLIK) SHOULD HAVE THE FOLLOWING
!   FORMAT:
!
!    IF NSTEP EQUALS -1
!     B(1) IS THE KTH (K=2,...,S) DIMENSIONAL WEIGHT
!
!
!    IF NSTEP EQUALS  0
!     B(1) IS THE SIGNAL TO NOISE RATIO, BETA
!
!
!    IF NSTEP EQUALS  1
!     B(1) IS THE DISTANCE/2 BETWEEN THE TWO OUTCOMES
!     B(2) IS THE MIDPOINT OF THE OUTCOMES
!
!    IF NSTEP EQUALS  2
!     B(1) IS THE ITH LEGISLATOR'S COORDINATE
!
!
!    LDATA(JT, ... )IS THE CHOICE VARIABLE
!        (E.G. YES=1, NO=2)
!    XDATA(JT) ARE THE INDIVIDUAL COORDINATES
!
!
!
      IF(NSTEP.EQ.-1)THEN
         KRCALL=1
         LRCALL=NRCALL
         BETA=BBB(1)
         WEIGHT=B(1)
      ENDIF
      IF(NSTEP.EQ.0)THEN
         KRCALL=1
         LRCALL=NRCALL
         BETA=B(1)
         WEIGHT=BBB(2)
      ENDIF
      IF(NSTEP.GT.0)THEN
         BETA=BBB(1)
         WEIGHT=BBB(2)
      ENDIF
  11  CONTINUE
      DO 27 JJ=1,5
  27  G(JJ)=0.0
      XL=0.0
!
!  ***** LEGISLATOR PHASE *****
!    COMPUTE DERIVATIVE FOR LEGISLATORS
!
!
      IF(NSTEP.EQ.2)THEN
         WTSQ=WEIGHT**2
         ICH=LDATA(NEQ,JT)
         IF(ICH.LE.0)RETURN
         D2Y=(ABS(B(1)-ZMID(JT,NDIM)+DYN(JT,NDIM)))**2
         D2N=(ABS(B(1)-ZMID(JT,NDIM)-DYN(JT,NDIM)))**2
         XBY=-D2Y*WTSQ/2.0
         XBY=BETA*EXP(XBY)*PSI(NEQ,JT,1)
         EXBY=EXP(XBY)
         XBN=-D2N*WTSQ/2.0
         XBN=BETA*EXP(XBN)*PSI(NEQ,JT,2)
         EXBN=EXP(XBN)
         PHIXI=EXBY+EXBN
         AA=ALOG(PHIXI)
         XL=-AA
         IF(ICH.EQ.1)XL=XL+XBY
         IF(ICH.EQ.2)XL=XL+XBN
         DLDX(1)=B(1)-ZMID(JT,NDIM)+DYN(JT,NDIM)
         DLDX(2)=B(1)-ZMID(JT,NDIM)-DYN(JT,NDIM)
         DLDX(1)=DLDX(1)*WTSQ
         DLDX(2)=DLDX(2)*WTSQ
         IF(ICH.EQ.1)G(1)=-DLDX(1)*XBY
         IF(ICH.EQ.2)G(1)=-DLDX(2)*XBN
         G(1)=G(1)+(DLDX(1)*XBY*EXBY+DLDX(2)*XBN*EXBN)/PHIXI
!
!  IF LEGISLATOR IS CONSTRAINED DUE TO SAG PROBLEM, SET DERIVATIVE TO ZERO
!
         G(1)=G(1)*ISENS(NEQ)
      ENDIF
!
!  ***** ROLL CALL PHASE *****
!  ROLL CALL DERIVATIVES
!
!
      IF(NSTEP.EQ.1)THEN
         XBIGL(NEQ)=0.0
         ICH=LDATA(JT,NEQ)
         IF(ICH.LE.0)RETURN
         DSQ(1,NEQ)=(WEIGHT*(XDATA(JT,NDIM)-B(2)+B(1)))**2
         DSQ(2,NEQ)=(WEIGHT*(XDATA(JT,NDIM)-B(2)-B(1)))**2
!
!   COMPUTE PHI
!
         PHI(NEQ)=0.0
         DO 6 I=1,NQ
         XB(I,NEQ)=-DSQ(I,NEQ)/(2.0)
!        IF(XB(I,NEQ).LT.-10.0)GO TO 10
         XB(I,NEQ)=BETA*EXP(XB(I,NEQ))*PSI(JT,NEQ,I)
!  10    IF(XB(I,NEQ).LT.-10.0)XB(I,NEQ)=0.0
         EXB(I,NEQ)=EXP(XB(I,NEQ))
  6      PHI(NEQ)=PHI(NEQ)+EXB(I,NEQ)
!
!   COMPUTE LIKELIHOOD
!
         CC=ALOG(PHI(NEQ))
         XL=-CC
         XBIGL(NEQ)=-CC
         IF(ICH.GT.0)XL=XB(ICH,NEQ)+XL
         IF(ICH.GT.0)XBIGL(NEQ)=XBIGL(NEQ)+XB(ICH,NEQ)
!
!    COMPUTE DERIVATIVES FOR ROLL CALL OUTCOMES
!
         WTSQ=WEIGHT**2
         DLDX(1)=XDATA(JT,NDIM)-B(2)+B(1)
         DLDX(2)=XDATA(JT,NDIM)-B(2)-B(1)
         IF(ICH.EQ.1)G(1)=-WTSQ*DLDX(1)*XB(1,NEQ)
         IF(ICH.EQ.2)G(1)= WTSQ*DLDX(2)*XB(2,NEQ)
!         G(1)=WTSQ*(DLDX(2)*XB(2,NEQ)-DLDX(1)*XB(1,NEQ))
         G(1)=G(1)-WTSQ*(DLDX(2)*XB(2,NEQ)*EXB(2,NEQ) -                 &
     &                   DLDX(1)*XB(1,NEQ)*EXB(1,NEQ))/PHI(NEQ)
         G(2)=WTSQ*DLDX(ICH)*XB(ICH,NEQ)
         G(2)=G(2)-WTSQ*(DLDX(1)*XB(1,NEQ)*EXB(1,NEQ)+                  &
     &                   DLDX(2)*XB(2,NEQ)*EXB(2,NEQ))/PHI(NEQ)
      ENDIF
!
!  ***** UTILITY FUNCTION PHASE *****
!  DERIVATIVES FOR BETA
!
!
      IF(NSTEP.EQ.0)THEN
         DO 26 II=KRCALL,LRCALL
         XBIGL(II)=0.0
         ICH=LDATA(JT,II)
         IF(ICH.LE.0)GO TO 26
         DSQ(1,II)=(WEIGHT*(XDATA(JT,NDIM)-ZMID(II,NDIM)+               &
     &             DYN(II,NDIM)))**2
         DSQ(2,II)=(WEIGHT*(XDATA(JT,NDIM)-ZMID(II,NDIM)-               &
     &             DYN(II,NDIM)))**2
!
!   COMPUTE PHI
!
         PHI(II)=0.0
         DO 66 I=1,NQ
         XB(I,II)=-DSQ(I,II)/(2.0)
!        IF(XB(I,II).LT.-10.0)GO TO 67
         XB(I,II)=BETA*EXP(XB(I,II))*PSI(JT,II,I)
!  67    IF(XB(I,II).LT.-10.0)XB(I,II)=0.0
         EXB(I,II)=EXP(XB(I,II))
   66    PHI(II)=PHI(II)+EXB(I,II)
!
!   COMPUTE LIKELIHOOD
!
         CC=ALOG(PHI(II))
         XL=XL-CC
         XBIGL(II)=-CC
         IF(ICH.GT.0)XL=XB(ICH,II)+XL
         IF(ICH.GT.0)XBIGL(II)=XBIGL(II)+XB(ICH,II)
  26     CONTINUE
!
!   COMPUTE FIRST PARTIALS OF BETA
!
         DO 77 II=KRCALL,LRCALL
         ICH=LDATA(JT,II)
         IF(ICH.LE.0) GO TO 77
         DO 7 K=1,NQ
         T=-EXB(K,II)/PHI(II)
         IF(K.EQ.ICH)T=T+1.0
         G(1)=G(1)+T*(1.0/BETA)*XB(K,II)
  7      CONTINUE
  77     CONTINUE
      ENDIF
!
!
! ***UTILITY FUNCTION PHASE
!  DERIVATIVES FOR WEIGHT
!
!
      IF(NSTEP.EQ.-1)THEN
         DO 25 II=KRCALL,LRCALL
         XBIGL(II)=0.0
         ICH=LDATA(JT,II)
         IF(ICH.LE.0)GO TO 25
         DSQ(1,II)=(WEIGHT*(XDATA(JT,NDIM)-ZMID(II,NDIM)+               &
     &             DYN(II,NDIM)))**2
         DSQ(2,II)=(WEIGHT*(XDATA(JT,NDIM)-ZMID(II,NDIM)-               &
     &             DYN(II,NDIM)))**2
!
!   COMPUTE PHI
!
         PHI(II)=0.0
         DO 68 I=1,NQ
         XB(I,II)=-DSQ(I,II)/(2.0)
!        IF(XB(I,II).LT.-10.0)GO TO 69
         XB(I,II)=BETA*EXP(XB(I,II))*PSI(JT,II,I)
!  69    IF(XB(I,II).LT.-10.0)XB(I,II)=0.0
         EXB(I,II)=EXP(XB(I,II))
  68     PHI(II)=PHI(II)+EXB(I,II)
!
!   COMPUTE LIKELIHOOD
!
         CC=ALOG(PHI(II))
         XL=XL-CC
         XBIGL(II)=-CC
         IF(ICH.GT.0)XL=XB(ICH,II)+XL
         IF(ICH.GT.0)XBIGL(II)=XBIGL(II)+XB(ICH,II)
  25     CONTINUE
!
!   COMPUTE FIRST PARTIALS OF WEIGHT (STH DIMENSION)
!
         DO 78 II=KRCALL,LRCALL
         ICH=LDATA(JT,II)
         IF(ICH.LE.0) GO TO 78
         DO 8 K=1,NQ
         T=-EXB(K,II)/PHI(II)
         IF(K.EQ.ICH)T=T+1.0
         G(1)=G(1)-(1.0/WEIGHT)*T*XB(K,II)*DSQ(K,II)
  8      CONTINUE
  78     CONTINUE
      ENDIF
      RETURN
      END
!
!  **************************************************************************
!    SUBROUTINE RPRINT---CALLED BY MAXLIK.  STORES NEW PARAMETER ESTIMATES
!       AND (OPTIONALLY) WRITES PARAMETER ESTIMATES TO DISK.  IN ROLL CALL
!       PHASE, IF CONVERGED ESTIMATES ARE "UNREALISTIC", GRID2
!       SUBROUTINE IS CALLED TO PRODUCED "APPROPRIATE" CONSTRAINED
!       ESTIMATES
!  **************************************************************************
!
      SUBROUTINE RPRINT(KPTSUM,NEQ,NP,NRCALL,NDUAL,NS,NDIM,NSTEP,NOPAR, &
     &                  JAN,JANLST,NITR,KSMIN,KSMAX,                    &
     &                  XLNL,B,V,SVLNL,FNLNL,ICNVG,ITT,NVAL,BBB,BBBB,   &
     &                  SSS,ZMID,XDATA,DYN,LDATA,PSI,XSAVE,ZSAVE,CSAVE, &
     &                  STDDVX,STDDVZ,YBIGL,YYBIGL,IPRINT)
!
!
      DIMENSION SSSE(25),TTTE(25),B(50),V(25,25),ALP(2),BTA(2),RR(2),   &
     &          XSP(50),XNEW(50),BBB(50),BBBB(25),SSS(100),             &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),                      &
     &          XSAVE(NP,2,2),ZSAVE(NRCALL,2,2),CSAVE(NRCALL,2),        &
     &          STDDVX(NP,25),STDDVZ(NRCALL,2,25),                      &
     &          YBIGL(NDUAL),YYBIGL(NDUAL)
!
! 1001 FORMAT(I5,2F13.4,3I4,4F9.4/34X,5F9.4/25X,6F9.4)
! 1002 FORMAT(1X,3I4,2I7,F14.4,4F9.4)
! 1003 FORMAT(1X,I4,1X,19A1,1X,2F13.4,3I4,4F9.4)
! 1010 FORMAT(1X,I4,2F13.4,3I4,4F9.4/30X,F13.4,4F9.4)
! 1011 FORMAT(I5,2F13.4,3I4,4F9.4/34X,5F9.4)
! 1012 FORMAT(I5,5X,7F10.4)
! 1013 FORMAT(5X,I5,7F10.4)
! 1014 FORMAT(10X,7F10.4)
      DO 10 I=1,NOPAR
      SE=SQRT(ABS(V(I,I)))
      T=0.0
      IF(SE.GT.0) T=B(I)/SE
      SSSE(I)=SE
      TTTE(I)=T
   10 CONTINUE
!
!  ***** UTILITY FUNCTION PHASE *****
!    WEIGHT FOR DIMENSIONS 2 AND ABOVE
!
      IF(NSTEP.EQ.-1)THEN
!         IF(IPRINT.EQ.1)THEN
!           WRITE(23,1011)NEQ,SVLNL,FNLNL,ITT,
!     C        NVAL,ICNVG,BBB(2),B(1),SSSE(1),TTTE(1)
!         ENDIF
!
!     KPTSUM=TOTAL # OF CHOICES.
!
           AA=XLNL/FLOAT(KPTSUM)
!
!     AA=EXP(AA) IS GEOMETRIC MEAN PROBABILITY OF UTILITY FUNCTION.
!
           AA=EXP(AA)
           BBB(2)=B(1)
           BBBB(NDIM)=B(1)
           CALL CORR(NP,NRCALL,XSAVE,ZSAVE,CSAVE,ALP,BTA,RR,NSTEP)
!           IF(IPRINT.EQ.1)THEN
!           WRITE(21,1002)NDIM,NITR,NSTEP,IITIME,JITIME,XLNL,AA,
!     C       RR(1),BBB(1),BBB(2)
!           ENDIF
      ENDIF
!
!  ***** UTILITY FUNCTION PHASE *****
!     BETA FOR FIRST DIMENSION
!
      IF(NSTEP.EQ.0)THEN
!           IF(IPRINT.EQ.1)THEN
!           WRITE(23,1011)NEQ,SVLNL,FNLNL,ITT,
!     C     NVAL,ICNVG,BBB(1),B(1),SSSE(1),TTTE(1)
!           ENDIF
!
!     KPTSUM=TOTAL # OF CHOICES.
!
           AA=XLNL/FLOAT(KPTSUM)
!
!     AA=EXP(AA) IS GEOMETRIC MEAN PROBABILITY OF UTILITY FUNCTION.
!
           AA=EXP(AA)
           BBB(1)=B(1)
           CALL CORR(NP,NRCALL,XSAVE,ZSAVE,CSAVE,ALP,BTA,RR,NSTEP)
!           IF(IPRINT.EQ.1)THEN
!           WRITE(21,1002)NDIM,NITR,NSTEP,IITIME,JITIME,XLNL,AA,
!     C       RR(1),BBB(1),BBB(2)
!           ENDIF
      ENDIF
!
!  ***** ROLL CALL PHASE *****
!
      IF(NSTEP.EQ.1)THEN
           DYN(NEQ,NDIM)=B(1)
           ZMID(NEQ,NDIM)=B(2)
           IF(NDIM.EQ.1)THEN
!
!     IF THE MIDPOINT < FURTHEST LEFT LEGISLATOR, CALL GRID ROUTINE TO
!        CALCULATE CONSTRAINED ESTIMATES OF MIDPOINT AND YES OUTCOME
!        FOR ROLL CALL
!
!           IF(B(2).LT.SSS(1))CALL GRID2(NEQ,B,3,NDIM,
!     C                 NP,NRCALL,NDUAL,
!     C                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,
!     C                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
!
!     IF THE MIDPOINT > FURTHEST RIGHT LEGISLATOR, CALL GRID ROUTINE TO
!        CALCULATE CONSTRAINED ESTIMATES OF MIDPOINT AND YES OUTCOME
!        FOR ROLL CALL
!
!           IF(B(2).GT.SSS(2))CALL GRID2(NEQ,B,2,NDIM,
!     C                 NP,NRCALL,NDUAL,
!     C                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,
!     C                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
!
!  DO BOTH OF THE ABOVE AS IN THE SUPERCOMPUTER VERSION -- CONSTRAIN
!   AT THE ENDPOINT WITH THE BEST LOG-LIKELIHOOD
!
           IF(B(2).LT.SSS(1).OR.B(2).GT.SSS(2))THEN
              CALL GRID2(NEQ,B,2,NDIM,                                  &
     &                 NP,NRCALL,NDUAL,                                 &
     &                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,            &
     &                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
              SAVEB1=B(1)
              SAVEB2=B(2)
              SAVEFNL=FNLNL
              CALL GRID2(NEQ,B,3,NDIM,                                  &
     &                 NP,NRCALL,NDUAL,                                 &
     &                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,            &
     &                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
              IF(FNLNL.LT.SAVEFNL)THEN
                 B(1)=SAVEB1
                 B(2)=SAVEB2
                 FNLNL=SAVEFNL
              ENDIF
           ENDIF
!
!     IF THE YES/NO OUTCOME IS < FURTHEST LEFT LEGISLATOR AND YES/NO
!         OUTCOME IS > FURTHEST RIGHT LEGISLATOR THEN CALL GRID ROUTINE
!         TO CONSTRAIN THE YES/NO OUTCOMES BUT RETAIN MIDPOINT
!                  IF DYN(.) > 0
!
           IF((B(2)-B(1)).LT.SSS(1).AND.(B(2)+B(1)).GT.SSS(2))          &
     &                         CALL GRID2(NEQ,B,1,NDIM,                 &
     &                 NP,NRCALL,NDUAL,                                 &
     &                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,            &
     &                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
!
!                  IF DYN(.) < 0
!
           IF((B(2)+B(1)).LT.SSS(1).AND.(B(2)-B(1)).GT.SSS(2))          &
     &                         CALL GRID2(NEQ,B,1,NDIM,                 &
     &                 NP,NRCALL,NDUAL,                                 &
     &                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,            &
     &                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
!
!  CHECK IF PREVIOUS DIMENSION SET EQUAL TO +1/-1 AND DO GRID SEARCH
!     ALONG TRACKS THROUGH HYPERSHERE
!
           DYN(NEQ,NDIM)=B(1)
           ZMID(NEQ,NDIM)=B(2)
           ENDIF
!
           SUM=0.0
           DO 9 K=1,NDIM
  9        SUM=SUM+ZMID(NEQ,K)**2
           IF(NDIM.GT.1.AND.SUM.GT.1.0)                                 &
     &                          CALL GRID2(NEQ,B,4,NDIM,                &
     &                 NP,NRCALL,NDUAL,                                 &
     &                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,            &
     &                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
           DYN(NEQ,NDIM)=B(1)
           ZMID(NEQ,NDIM)=B(2)
           STDDVZ(NEQ,1,NDIM)=SSSE(1)
           STDDVZ(NEQ,2,NDIM)=SSSE(2)
           ZSAVE(NEQ,2,1)=ZSAVE(NEQ,1,1)
           ZSAVE(NEQ,1,1)=B(1)
           CSAVE(NEQ,2)=CSAVE(NEQ,1)
           CSAVE(NEQ,1)=B(2)
           ZSAVE(NEQ,2,2)=FNLNL
           ZSAVE(NEQ,1,2)=SVLNL
!
      ENDIF
!
!  ***** LEGISLATOR PHASE *****
!
      IF(NSTEP.EQ.2)THEN
           SUM=0.0
!
!  IMPOSE CONSTRAINT THAT LEGISLATOR COORDINATES HAVE TO LIE WITHIN A
!     UNIT HYPERSPHERE
!
           IF(NDIM.GT.1)THEN
              BMAX=-99999999999999.0
              XDATA(NEQ,NDIM)=B(1)
              DO 11 K=1,NDIM
              SUM=SUM+XDATA(NEQ,K)**2
  11          CONTINUE
              IF(SUM.GT.1.0)THEN
!                 WRITE(37,1012)NEQ,(XDATA(NEQ,K),K=1,NDIM),FNLNL
!
!  NORMALIZE LEGISLATOR COORDINATES TO BE ON SURFACE OF HYPERSPHERE
!
                 DO 13 K=1,NDIM
  13             XSP(K)=XDATA(NEQ,K)*(1.0/SQRT(SUM))
!                 WRITE(37,1014)(XSP(K),K=1,NDIM)
                 XDATA(NEQ,NDIM)=0.0
!
!  CREATE POINTS ON HYPERSPHERE SURFACE FROM NEW TO OLD POINT
!
                 WTSQ=BBB(2)**2
                 DO 14 JJ=1,21
                 XLL=0.0
                 XJJ=FLOAT(JJ)
                 ASUM=0.0
                 DO 15 K=1,NDIM
                 XNEW(K)=XSP(K)+((XJJ-1.0)/20.0)*(XDATA(NEQ,K)-XSP(K))
  15             ASUM=ASUM+XNEW(K)**2
                 DO 18 K=1,NDIM
  18             XNEW(K)=XNEW(K)*(1.0/SQRT(ASUM))
                 DO 12 JT=1,NRCALL
                 ICH=LDATA(NEQ,JT)
                 IF(ICH.LE.0)GO TO 12
                 D2Y=0.0
                 D2N=0.0
                 DO 16 K=1,NDIM
                 D2Y=D2Y+(BBBB(K)*(XNEW(K)-ZMID(JT,K)+DYN(JT,K)))**2
                 D2N=D2N+(BBBB(K)*(XNEW(K)-ZMID(JT,K)-DYN(JT,K)))**2
  16             CONTINUE
                 XBY=-D2Y/2.0
                 XBY=BBB(1)*EXP(XBY)
                 EXBY=EXP(XBY)
                 XBN=-D2N/2.0
                 XBN=BBB(1)*EXP(XBN)
                 EXBN=EXP(XBN)
                 PHIXI=EXBY+EXBN
                 AA=ALOG(PHIXI)
                 XL=-AA
                 IF(ICH.EQ.1)XL=XL+XBY
                 IF(ICH.EQ.2)XL=XL+XBN
                 XLL=XLL+XL
  12             CONTINUE
!                 WRITE(37,1013)JJ,(XNEW(K),K=1,NDIM),XLL
                 AMX=AMAX1(XLL,BMAX)
                 IF(ABS(AMX-XLL).LE..00001)KJJ=JJ
                 BMAX=AMX
  14             CONTINUE
                 XJJ=FLOAT(KJJ)
                 ASUM=0.0
                 DO 17 K=1,NDIM
                 XDATA(NEQ,K)=XSP(K)+                                   &
     &                    ((XJJ-1.0)/20.0)*(XDATA(NEQ,K)-XSP(K))
  17             ASUM=ASUM+XDATA(NEQ,K)**2
                 DO 19 K=1,NDIM
  19             XDATA(NEQ,K)=XDATA(NEQ,K)*(1.0/SQRT(ASUM))
                 B(1)=XDATA(NEQ,NDIM)
!                 WRITE(37,1013)KJJ,(XDATA(NEQ,K),K=1,NDIM),AMX
                 FNLNL=-AMX
              ENDIF
           ENDIF
           XDATA(NEQ,NDIM)=B(1)
           XSAVE(NEQ,2,1)=XSAVE(NEQ,1,1)
           XSAVE(NEQ,1,1)=B(1)
           XSAVE(NEQ,2,2)=FNLNL
           XSAVE(NEQ,1,2)=SVLNL
           STDDVX(NEQ,NDIM)=SSSE(1)
      ENDIF
      RETURN
      END
!
!  ***************************************************************************
!     SUBROUTINE CROSS---CALLED BY MAXLIK AND MAIN.  CALCULATES CROSS
!        CLASSIFICATIONS USING PREDICTED YES AND NO PROBABILITIES VERSUS
!        ACTUAL CHOICES.
!  ***************************************************************************
!
      SUBROUTINE CROSS(II,LL,                                           &
     &                 NP,NRCALL,NDIM,NSTEP,JAN,JANLST,                 &
     &                 ZMID,XDATA,DYN,LDATA,PSI,BBB,BBBB)
      DIMENSION L(2),LL(2,2),                                           &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),                      &
     &          BBB(50),BBBB(25)
      BETA=BBB(1)
      WEIGHT=BBB(2)
      WTSQ=WEIGHT*WEIGHT/2.0
      DO 5 I=1,2
      L(I)=0
      DO 5 J=1,2
  5   LL(I,J)=0
!
!  ***** UTILITY FUNCTION PHASE *****
!
      IF(NSTEP.EQ.0)THEN
           DO 100 I=1,NP
           ICH=LDATA(I,II)
           IF(ICH.LE.0)GO TO 100
           L(ICH)=L(ICH)+1
           SUMY=0.0
           SUMN=0.0
           DO 96 K=1,NDIM
!           SUMY=SUMY+(XDATA(I,K)-ZMID(II,K)+DYN(II,K))**2
!           SUMN=SUMN+(XDATA(I,K)-ZMID(II,K)-DYN(II,K))**2
           SUMY=SUMY+(BBBB(K)*(XDATA(I,K)-ZMID(II,K)+DYN(II,K)))**2
           SUMN=SUMN+(BBBB(K)*(XDATA(I,K)-ZMID(II,K)-DYN(II,K)))**2
   96      CONTINUE
           SUMY=BETA*EXP(-SUMY/2.0)
           SUMN=BETA*EXP(-SUMN/2.0)
           PIY=EXP(SUMY)/(EXP(SUMY)+EXP(SUMN))
           PIN=1.0-PIY
           IF(PIY.GE.PIN)LL(1,ICH)=LL(1,ICH)+1
           IF(PIY.LT.PIN)LL(2,ICH)=LL(2,ICH)+1
  100      CONTINUE
      ENDIF
!
!  ***** ROLL CALL PHASE *****
!
      IF(NSTEP.EQ.1)THEN
           DO 101 I=1,NP
           ICH=LDATA(I,II)
           IF(ICH.LE.0)GO TO 101
           L(ICH)=L(ICH)+1
           SUMY=0.0
           SUMN=0.0
           SUMYY=0.0
           SUMNN=0.0
           DO 97 K=1,NDIM
           SUMYY=SUMYY+((XDATA(I,K)-ZMID(II,K)+DYN(II,K)))**2
           SUMNN=SUMNN+((XDATA(I,K)-ZMID(II,K)-DYN(II,K)))**2
           SUMY=SUMY+(BBBB(K)*(XDATA(I,K)-ZMID(II,K)+DYN(II,K)))**2
           SUMN=SUMN+(BBBB(K)*(XDATA(I,K)-ZMID(II,K)-DYN(II,K)))**2
  97       CONTINUE
           IF(SUMY.LE.SUMN)LL(1,ICH)=LL(1,ICH)+1
           IF(SUMY.GT.SUMN)LL(2,ICH)=LL(2,ICH)+1
           SUMY=BETA*EXP(-SUMY/2.0)
           SUMN=BETA*EXP(-SUMN/2.0)
           PIY=EXP(SUMY)/(EXP(SUMY)+EXP(SUMN))
           PIN=1.0-PIY
!           IF(PIY.GE.PIN)LL(1,ICH)=LL(1,ICH)+1
!           IF(PIY.LT.PIN)LL(2,ICH)=LL(2,ICH)+1
  101      CONTINUE
      ENDIF
!
!  ***** LEGISLATOR PHASE *****
!
      IF(NSTEP.EQ.2)THEN
           DO 200 I=1,NRCALL
           ICH=LDATA(II,I)
           IF(ICH.LE.0)GO TO 200
!
!    USE THIS STATEMENT IF, ON THE LAST ITERATION, YOU DESIRE TO CALCULATE
!      THE CLASSIFICATIONS OF UNCONSTRAINED ROLL CALLS ONLY
!           IF(ABS(ZMID(I,K)).EQ.1.0.AND.JANLST.EQ.JAN)GO TO 200
!
           L(ICH)=L(ICH)+1
           SUMY=0.0
           SUMN=0.0
           DO 98 K=1,NDIM
           SUMY=SUMY+(BBBB(K)*(XDATA(II,K)-ZMID(I,K)+DYN(I,K)))**2
           SUMN=SUMN+(BBBB(K)*(XDATA(II,K)-ZMID(I,K)-DYN(I,K)))**2
  98       CONTINUE
           SUMY=BETA*EXP(-SUMY/2.0)
           SUMN=BETA*EXP(-SUMN/2.0)
           PIY=EXP(SUMY)/(EXP(SUMY)+EXP(SUMN))
           PIN=1.0-PIY
           IF(PIY.GE.PIN.AND.ICH.EQ.1)LL(1,1)=LL(1,1)+1
           IF(PIY.GE.PIN.AND.ICH.EQ.2)LL(1,2)=LL(1,2)+1
           IF(PIY.LT.PIN.AND.ICH.EQ.1)LL(2,1)=LL(2,1)+1
           IF(PIY.LT.PIN.AND.ICH.EQ.2)LL(2,2)=LL(2,2)+1
  200      CONTINUE
      ENDIF
      RETURN
      END
!
!  *************************************************************************
!    SUBROUTINE GRID2---CALLED BY RPRINT.  PERFORMS GRID SEARCH WHEN
!      ESTIMATES OF ROLL CALL PARAMETERS ARE UNREALISTIC.
!             IF KKS = 1 THEN THE YES AND NO OUTCOMES ARE OFF THE
!                        OPPOSITE ENDS OF THE DIMENSION
!             IF KKS = 2 THEN ZMID(.) > FURTHEST RIGHT LEGISLATOR
!             IF KKS = 3 THEN ZMID(.) < FURTHEST LEFT LEGISLATOR
!  *************************************************************************
!
      SUBROUTINE GRID2(NEQ,B,KKS,NDIM,                                  &
     &                 NP,NRCALL,NDUAL,                                 &
     &                 KSMIN,KSMAX,ZMID,XDATA,DYN,LDATA,PSI,            &
     &                 BBB,BBBB,SSS,YBIGL,YYBIGL,SVLNL,FNLNL)
      DIMENSION ZY(200),B(50),ZN(200),ZM(50,50),PTA(50),                &
     &          PTB(50),PTC(50),CMINUSZ(50),BMINUSA(50),                &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          LDATA(NP,NRCALL),PSI(NP,NRCALL,2),                      &
     &          BBB(50),BBBB(25),SSS(100),                              &
     &          YBIGL(NDUAL),YYBIGL(NDUAL)
 1000 FORMAT(' OFF OPPOSITE ENDS OUTCOMES'/1X,3I5,4F10.4)
 1001 FORMAT(' NEW ESTIMATES OF DYN AND ZM FROM GRID SEARCH'            &
     &        /6X,2I5,4F10.4)
 1002 FORMAT(' UNREALISTIC ESTIMATES OF DYN AND ZM'/1X,3I5,4F10.4)
 1003 FORMAT(' UNREALISTIC ESTIMATES OF DYN AND ZM'/1X,3I5/6F10.4/      &
     & 6F10.4)
 1004 FORMAT(' NEW ESTIMATES OF MIDPOINTS NDIM-1 AND NDIM'/             &
     &          1X,3I5,4F10.4/16X,6F10.4/16X,6F10.4/16X,6F10.4/         &
     &          16X,6F10.4)
 1005 FORMAT(' BIZARRE RESULTS'/1X,3I5,3F10.4/16X,6F10.4/               &
     &          16X,6F10.4/16X,6F10.4/16X,6F10.4)
 1006 FORMAT(' BIG ERROR',3I5,3F10.4/25X,4F10.4)
 1007 FORMAT(' BIZARRE RESULTS2'/1X,3I5,5F10.4)
      NCUT=100
      NCUTCUT=101
      BMAX=-999999999.0
      BMAXY=-999999999.0
      BMAXN=-999999999.0
      WEIGHT=BBB(2)
      KJ=1
!
!  TAKE INTO ACCOUNT POLARITY OF ROLL CALL
!
      XSIGN=+1.0
      IF(B(1).LE.0.0)XSIGN=-1.0
!
!  ***** OUTCOMES OFF OPPOSITE ENDS OF DIMENSION *****
!
      IF(KKS.EQ.1)THEN
!           WRITE(11,1000)KKS,NDIM,NEQ,B(1),B(2),YBIGL(NEQ),SVLNL
           DYNDYN=ABS(B(1))
           XINC=DYNDYN/FLOAT(NCUT)
!
!  SETUP ANCHOR POINT (ZF) TO START GRID SEARCH SUCH THAT AT LEAST ONE OF THE
!     OUTCOMES IS IN THE INTERIOR OF THE DIMENSION
!
!          ZF/KSMIN--------------0--ZMID-------KSMAX
!
           IF(B(2).GT.0.0)ZF=XDATA(KSMIN,NDIM)
!
!       ZF----KSMIN----ZMID------0-------------KSMAX
!
           IF(B(2).LE.0.0)ZF=B(2)-ABS(XDATA(KSMAX,NDIM)-B(2))
!
!  SETUP NCUT SEARCH POINTS FROM FAR LEFT ANCHOR POINT (ZF) INTO INTERIOR
!
           DO 1 J=1,NCUT
           XJ=FLOAT(J)
  1        ZY(J)=ZF+XINC*(XJ-1.0)
!
!  CALCULATE LOG LIKELIHOOD FOR EACH OF THE NCUT POINTS
!
           DO 2 J=1,NCUT
           XLX=0.0
           DO 3 K=1,NP
           ICH=LDATA(K,NEQ)
           IF(ICH.LE.0)GO TO 3
           DY=(WEIGHT*(XDATA(K,NDIM)-ZY(J)))**2
           DN=(WEIGHT*(XDATA(K,NDIM)-2.0*B(2)+ZY(J)))**2
           AA=BBB(1)*EXP(-DY/2.0)*PSI(K,NEQ,1)
           BB=BBB(1)*EXP(-DN/2.0)*PSI(K,NEQ,2)
           PHI=EXP(AA)+EXP(BB)
           CC=ALOG(PHI)
!  TO CATCH 1D NEAR-PERFECT VOTING ERRORS
           IF(CC.GT.100000000)THEN
                CALL ECHOEVENT(8)
                CALL FLUSHCON()
                CALL PROCEVENT()
           ENDIF
           IF(B(1).GE.0.0)THEN
                IF(ICH.EQ.1)XLX=XLX+AA-CC
                IF(ICH.EQ.2)XLX=XLX+BB-CC
           ENDIF
           IF(B(1).LT.0.0)THEN
                IF(ICH.EQ.2)XLX=XLX+AA-CC
                IF(ICH.EQ.1)XLX=XLX+BB-CC
           ENDIF
  3        CONTINUE
           AMX=AMAX1(XLX,BMAX)
           IF(ABS(AMX-XLX).LE..0001)KJ=J
           BMAX=AMX
  2        CONTINUE
           B(1)=(B(2)-ZY(KJ))*XSIGN
      ENDIF
!
!
! ***** ZMID(.) > FURTHEST RIGHT LEGISLATOR *****
!                    OR
! ***** ZMID(.) < FURTHEST LEFT LEGISLATOR *****
!
!
      IF(KKS.EQ.2.OR.KKS.EQ.3)THEN
!        IF(NDIM.EQ.2)WRITE(11,1002)KKS,NDIM,NEQ,B(1),B(2),YBIGL(NEQ),
!     C                             SVLNL
        IF(NDIM.EQ.1)THEN
!
!  SET MIDPOINT TO FURTHEST RIGHT LEGISLATOR
!
           IF(KKS.EQ.2)B(2)=XDATA(KSMAX,NDIM)
           IF(KKS.EQ.3)B(2)=XDATA(KSMIN,NDIM)
           DYNDYN=1.0
           XINC=DYNDYN/FLOAT(NCUT)
!
!  SET ANCHOR POINT (ZF) FOR GRID SEARCH TO 1 UNIT LEFT OF FURTHEST RIGHT
!     LEGISLATOR/FUTHEST LEFT LEGISLATOR
!
           ZF=B(2)-DYNDYN
           DO 11 J=1,NCUTCUT
           XJ=FLOAT(J)
  11       ZY(J)=ZF+XINC*(XJ-1.0)
           DO 12 J=1,NCUTCUT
           XLX=0.0
           XLX1=0.0
           XLX2=0.0
           DO 13 K=1,NP
           ICH=LDATA(K,NEQ)
           IF(ICH.LE.0)GO TO 13
           DY=(WEIGHT*(XDATA(K,NDIM)-ZY(J)))**2
           DN=(WEIGHT*(XDATA(K,NDIM)-2.0*B(2)+ZY(J)))**2
           AA=BBB(1)*EXP(-DY/2.0)*PSI(K,NEQ,1)
           BB=BBB(1)*EXP(-DN/2.0)*PSI(K,NEQ,2)
           PHI=EXP(AA)+EXP(BB)
           CC=ALOG(PHI)
!           IF(B(1).GE.0.0)THEN
                IF(ICH.EQ.1)XLX1=XLX1+AA-CC
                IF(ICH.EQ.2)XLX1=XLX1+BB-CC
!           ENDIF
!           IF(B(1).LT.0.0)THEN
                IF(ICH.EQ.2)XLX2=XLX2+AA-CC
                IF(ICH.EQ.1)XLX2=XLX2+BB-CC
!           ENDIF
  13       CONTINUE
           XLX=AMAX1(XLX1,XLX2)
           IF(J.EQ.1)AMX=XLX
           IF(J.GT.1)AMX=AMAX1(XLX,BMAX)
           IF(ABS(AMX-XLX).LE..0001)THEN
              KJ=J
              IF(ABS(XLX-XLX1).LE..0001)XSIGN=+1.0
              IF(ABS(XLX-XLX2).LE..0001)XSIGN=-1.0
           ENDIF
           BMAX=AMX
  12       CONTINUE
           B(1)=(B(2)-ZY(KJ))*XSIGN
        ENDIF
!
!
        IF(NDIM.GT.1)THEN
!
!  SET MAXIMUM VALUE OF MIDPOINT TO EDGE OF UNIT HYPERSPHERE
!
           SUM=0.0
           NDIMM1=NDIM-1
           DO 33 KK=1,NDIMM1
           SUM=SUM+ZMID(NEQ,KK)**2
  33       CONTINUE
           B2MAX=0.0
           IF(SUM.LT.1.0)THEN
               B2MAX=SQRT(ABS(1.0-SUM))
           ENDIF
           IF(B2MAX.EQ.0.0)THEN
               KKS=4
               GO TO 69
           ENDIF
           KKCUT=10
           KKKCUT=2*KKCUT + 1
           XINC=B2MAX/FLOAT(KKCUT)
!
!  CREATE YEA AND NAY TRACKS THROUGH THE HYPERSPHERE
!
           DO 34 J=1,KKKCUT
           XJ=FLOAT(J)
           ZY(J)=B2MAX-XINC*(XJ-1.0)
  34       ZN(J)=B2MAX-XINC*(XJ-1.0)
!
!  SEARCH FOR THE OPTIMAL PAIR OF OUTCOME COORDINATES
!
!
           XLDEF=0.0
           DO 35 J=1,KKKCUT
           DO 35 L=1,KKKCUT
           XLXY=0.0
           XLXN=0.0
           DO 36 K=1,NP
           ICH=LDATA(K,NEQ)
           IF(ICH.LE.0)GO TO 36
!
!  CALCULATE LOG-LIKELIHOOD DUE TO NDIM-1 PREVIOUS DIMENSIONS
!
           IF(J.EQ.KKKCUT.AND.L.EQ.KKKCUT)THEN
               AA=BBB(1)*PSI(K,NEQ,1)
               BB=BBB(1)*PSI(K,NEQ,2)
               PHI=EXP(AA)+EXP(BB)
               CC=ALOG(PHI)
               IF(ICH.EQ.1)XLDEF=XLDEF+AA-CC
               IF(ICH.EQ.2)XLDEF=XLDEF+BB-CC
           ENDIF
           DY=(WEIGHT*(XDATA(K,NDIM)-ZY(J)))**2
           DN=(WEIGHT*(XDATA(K,NDIM)-ZN(L)))**2
           AA=BBB(1)*EXP(-DY/2.0)*PSI(K,NEQ,1)
           BB=BBB(1)*EXP(-DN/2.0)*PSI(K,NEQ,2)
           PHI=EXP(AA)+EXP(BB)
           CC=ALOG(PHI)
!
!  ASSUME YEA ALTERNATIVE IS TO THE LEFT (POSITIVE SPREAD)
!
           IF(ICH.EQ.1)XLXY=XLXY+AA-CC
           IF(ICH.EQ.2)XLXY=XLXY+BB-CC
!
!  ASSUME YEA ALTERNATIVE IS TO THE RIGHT (NEGATIVE SPREAD)
!
           IF(ICH.EQ.2)XLXN=XLXN+AA-CC
           IF(ICH.EQ.1)XLXN=XLXN+BB-CC
  36       CONTINUE
           AMXY=AMAX1(XLXY,BMAXY)
           IF(ABS(AMXY-XLXY).LE..00001)THEN
              KJ=J
              KL=L
           ENDIF
           BMAXY=AMXY
           AMXN=AMAX1(XLXN,BMAXN)
           IF(ABS(AMXN-XLXN).LE..00001)THEN
              LJ=J
              LL=L
           ENDIF
           BMAXN=AMXN
  35       CONTINUE
           IF(AMXY.GT.AMXN)THEN
              B(1)=ABS(ZY(KJ)-ZN(KL))/2.0
              B(2)=(ZY(KJ)+ZN(KL))/2.0
              AMX=AMXY
           ENDIF
           IF(AMXY.LE.AMXN)THEN
              B(1)=ABS(ZY(LJ)-ZN(LL))/2.0
              B(2)=(ZY(LJ)+ZN(LL))/2.0
              AMX=AMXN
           ENDIF
           IF(XLDEF.GE.AMX)THEN
              B(1)=0.0
              B(2)=0.0
              AMX=SVLNL
           ENDIF
           IF(ABS(B(1)).LE.0.01)THEN
              B(1)=0.0
              B(2)=0.0
              AMX=SVLNL
           ENDIF
        ENDIF
      ENDIF
!
!
!  IF NDIM-1 IS AT EDGE OF HYPERSPHERE, CREATE TRACKS PARALLEL TO
!    CUTTING LINE THROUGH HYPERSPHERE AND COMPUTE NDIM-1 & NDIM
!    COORDINATES
!
!
  69  CONTINUE
      IF(KKS.EQ.4)THEN
        BMAX=-9999999999.0
!        IF(NDIM.EQ.2)WRITE(21,1003)KKS,NDIM,NEQ,
!     C                 (DYN(NEQ,JJ),ZMID(NEQ,JJ),JJ=1,NDIM),
!     C                 B(1),B(2),YBIGL(NEQ),SVLNL
!
!  COMPUTE TRACKS THROUGH CENTER OF HYPERSPHERE
!
           SUM=0.0
           DO 54 K=1,NDIM
  54       SUM=SUM+DYN(NEQ,K)**2
           T=1.0/SQRT(SUM)
           XPSI=0.0
           XGAM=0.0
!
!  POINTS A AND B ARE ON THE SURFACE OF THE UNIT HYPERSPHERE ONE UNIT APART
!
           DO 53 K=1,NDIM
           PTA(K)=T*(-DYN(NEQ,K))
           PTB(K)=T*( DYN(NEQ,K))
           XPSI=XPSI+(ZMID(NEQ,K)-PTA(K))**2
  53       XGAM=XGAM+(ZMID(NEQ,K)-PTB(K))**2
           COSB=(XPSI-XGAM-4.0)/(-4.0*SQRT(XGAM))
           COSA=(XGAM-XPSI-4.0)/(-4.0*SQRT(XPSI))
!
!  IF CUTTING LINE DOES NOT INTERSECT HYPERSPHERE SET CURRENT COORDINATES
!    TO ZERO AND RETURN
!
           IF(COSB.LT.0.0.OR.COSA.LT.0.0)THEN
!              WRITE(11,1007)KKS,NDIM,NEQ,T,XPSI,XGAM,COSB,COSA
              B(1)=0.0
              B(2)=0.0
              YBIGL(NEQ)=SVLNL
              RETURN
           ENDIF
!
!  CHECK RESULTS
!
           DALPHA=SQRT(XGAM)*COSB
           DBETA =SQRT(XPSI)*COSA
           IF(ABS(2.0-DALPHA-DBETA).GT..001)THEN
!              WRITE(11,1006)KKS,NDIM,NEQ,T,XPSI,XGAM,COSB,COSA,
!     C                       DALPHA,DBETA
              STOP
           ENDIF
!
!  FIND POINT C ALONG A TO B LINE THROUGH CENTER OF HYPERSPHERE
!
           DO 52 K=1,NDIM
  52       PTC(K)=PTA(K)+(DBETA/2.0)*(PTB(K)-PTA(K))
!
!  FIND POINTS D AND E -- INTERSECTION OF CUTTING LINE THROUGH
!       THE SURFACE OF THE HYPERSPHERE
!
           AA=0.0
           BB=0.0
           DO 51 K=1,NDIM
           CMINUSZ(K)=PTC(K)-ZMID(NEQ,K)
           BB=BB+2.0*CMINUSZ(K)*ZMID(NEQ,K)
  51       AA=AA+CMINUSZ(K)**2
           CC=0.0
           DO 58 K=1,NDIM-1
  58       CC=CC+ZMID(NEQ,K)**2
           CC=CC+B(2)**2
           CC=CC-1.0
           IF(BB*BB.LT.(4.0*AA*CC))THEN
!              WRITE(11,1005)KKS,NDIM,NEQ,AA,BB,CC,
!     C               (DYN(NEQ,JJ),ZMID(NEQ,JJ),JJ=1,NDIM),
!     C               (PTA(JJ),JJ=1,NDIM),(PTB(JJ),JJ=1,NDIM),
!     C               (PTC(JJ),JJ=1,NDIM),
!     C               DALPHA,DBETA,XPSI,XGAM
              B(1)=0.0
              B(2)=0.0
              YBIGL(NEQ)=SVLNL
              RETURN
           ENDIF
           IF(BB*BB.GT.(4.0*AA*CC))THEN
              T1=(-BB + SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
              T2=(-BB - SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
              DO 50 K=1,NDIM
              PTA(K)=ZMID(NEQ,K)+T1*CMINUSZ(K)
              PTB(K)=ZMID(NEQ,K)+T2*CMINUSZ(K)
              BMINUSA(K)=PTB(K)-PTA(K)
  50          CONTINUE
           ENDIF
           DO 55 J=1,11
           XJ=FLOAT(J)
           DO 55 K=1,NDIM
           ZM(J,K)=PTA(K)+(1.0-.1*(XJ-1.0))*BMINUSA(K)
  55       CONTINUE
!
!  CALCULATE LOG-LIKELIHOOD
!
           DO 56 L=1,11
           XLX=0.0
           DO 57 K=1,NP
           ICH=LDATA(K,NEQ)
           IF(ICH.LE.0)GO TO 57
           DY=0.0
           DN=0.0
           DO 48 KK=1,NDIM
           IF(KK.LE.NDIM-1)THEN
              DY=DY+(BBBB(KK)*(XDATA(K,KK)-ZM(L,KK)+DYN(NEQ,KK)))**2
              DN=DN+(BBBB(KK)*(XDATA(K,KK)-ZM(L,KK)-DYN(NEQ,KK)))**2
           ENDIF
           IF(KK.EQ.NDIM)THEN
              DY=DY+(WEIGHT*(XDATA(K,KK)-ZM(L,KK)+DYN(NEQ,KK)))**2
              DN=DN+(WEIGHT*(XDATA(K,KK)-ZM(L,KK)-DYN(NEQ,KK)))**2
           ENDIF
  48       CONTINUE
           AA=BBB(1)*EXP(-DY/2.0)
           BB=BBB(1)*EXP(-DN/2.0)
           PHI=EXP(AA)+EXP(BB)
           CC=ALOG(PHI)
           IF(ICH.EQ.1)XLX=XLX+AA-CC
           IF(ICH.EQ.2)XLX=XLX+BB-CC
  57       CONTINUE
           AMX=AMAX1(XLX,BMAX)
           IF(ABS(AMX-XLX).LE..00001)THEN
              KL=L
           ENDIF
           BMAX=AMX
  56       CONTINUE
           YBIGL(NEQ)=AMX
           FNLNL=AMX
!           WRITE(11,1004)KKS,NDIM,NEQ,(ZM(KL,JJ),JJ=1,NDIM),AMX,
!     C                    (PTC(JJ),JJ=1,NDIM),(CMINUSZ(JJ),JJ=1,NDIM),
!     C                    (PTA(JJ),JJ=1,NDIM),(PTB(JJ),JJ=1,NDIM),
!     C                    T1,T2,XPSI,XGAM,COSA,COSB,DALPHA,DBETA
           DO 47 K=1,NDIM
           ZMID(NEQ,K)=ZM(KL,K)
  47       CONTINUE
           B(2)=ZM(KL,NDIM)
!           IF(NDIM.EQ.2)WRITE(21,1001)NDIM,NEQ,B(1),B(2),AMX
           RETURN
      ENDIF
!
      YBIGL(NEQ)=AMX
      FNLNL=AMX
!      IF(NDIM.EQ.2)WRITE(21,1001)NDIM,NEQ,B(1),B(2),AMX
      RETURN
      END
!
!  ***************************************************************************
!    SUBROUTINE CORR--CALLED BY MAIN, RPRINT, AND OUTWRT.
!       CALCULATES CORRELATION BETWEEN CURRENT
!       SET OF PARAMETERS AND PARAMETERS FROM PREVIOUS ITERATION.
!       IF KHR = 1 CALCULATE CORRELATIONS FOR ROLL CALL PHASE.
!       IF KHR = 2 CALCULATE CORRELATIONS FOR LEGISLATOR PHASE.
!       IF KHR = 0 CALCULATE CORRELATIONS FOR ALL PARAMETERS
!  ***************************************************************************
!
      SUBROUTINE CORR(NP,NRCALL,XSAVE,ZSAVE,CSAVE,ALP,BTA,RR,KHR)
      DIMENSION ALP(2),BTA(2),RR(2),                                    &
     &          XSAVE(NP,2,2),ZSAVE(NRCALL,2,2),CSAVE(NRCALL,2)
!
!  CALCULATE CORRELATIONS FOR ROLL CALL PHASE
!
      IF(KHR.EQ.1)THEN
           ASUM=0.0
           BSUM=0.0
           CSUM=0.0
           DSUM=0.0
           ESUM=0.0
           AASUM=0.0
           BBSUM=0.0
           CCSUM=0.0
           DDSUM=0.0
           EESUM=0.0
           DO 3 I=1,NRCALL
           ASUM=ASUM+ZSAVE(I,1,1)
           AASUM=AASUM+CSAVE(I,1)
           BSUM=BSUM+ZSAVE(I,1,1)**2
           BBSUM=BBSUM+CSAVE(I,1)**2
           CSUM=CSUM+ZSAVE(I,2,1)
           CCSUM=CCSUM+CSAVE(I,2)
           DSUM=DSUM+ZSAVE(I,2,1)**2
           DDSUM=DDSUM+CSAVE(I,2)**2
           EESUM=EESUM+CSAVE(I,1)*CSAVE(I,2)
           ESUM=ESUM+ZSAVE(I,1,1)*ZSAVE(I,2,1)
  3        CONTINUE
           XCALL=FLOAT(NRCALL)
           AA=XCALL*ESUM-ASUM*CSUM
           BB=XCALL*BSUM-ASUM*ASUM
           CC=XCALL*DSUM-CSUM*CSUM
           DD=DSUM*ASUM-CSUM*ESUM
!      IF(NDIM.GE.8)WRITE(21,1994)XCALL,AA,BB,CC,DD,
!     C                          ASUM,BSUM,CSUM,DSUM,ESUM
! 1994 FORMAT(10F10.4)
           ALP(1)=0.0
           IF(CC.GT.0.0)ALP(1)=DD/CC
           BTA(1)=0.0
           IF(CC.GT.0.0)BTA(1)=AA/CC
           RR(1)=0.0
           IF((ABS(BB*CC)).GT.0.0)RR(1)=AA/SQRT(ABS(BB*CC))
           RR(1)=RR(1)*RR(1)
           AA=XCALL*EESUM-AASUM*CCSUM
           BB=XCALL*BBSUM-AASUM*AASUM
           CC=XCALL*DDSUM-CCSUM*CCSUM
           DD=DDSUM*AASUM-CCSUM*EESUM
!      IF(NDIM.GE.8)WRITE(21,1994)XCALL,AA,BB,CC,DD,
!     C                          AASUM,BBSUM,CCSUM,DDSUM,EESUM
           ALP(2)=0.0
           IF(CC.GT.0.0)ALP(2)=DD/CC
           BTA(2)=0.0
           IF(CC.GT.0.0)BTA(2)=AA/CC
           RR(2)=0.0
           IF((ABS(BB*CC)).GT.0.0)RR(2)=AA/SQRT(ABS(BB*CC))
           RR(2)=RR(2)*RR(2)
      ENDIF
!
!  CALCULATE CORRELATIONS FOR LEGISLATOR PHASE
!
      IF(KHR.EQ.2)THEN
           ASUM=0.0
           BSUM=0.0
           CSUM=0.0
           DSUM=0.0
           ESUM=0.0
           DO 11 I=1,NP
           ASUM=ASUM+XSAVE(I,1,1)
           BSUM=BSUM+XSAVE(I,1,1)**2
           CSUM=CSUM+XSAVE(I,2,1)
           DSUM=DSUM+XSAVE(I,2,1)**2
  11       ESUM=ESUM+XSAVE(I,1,1)*XSAVE(I,2,1)
           XNP=FLOAT(NP)
           AA=XNP*ESUM-ASUM*CSUM
           BB=XNP*BSUM-ASUM*ASUM
           CC=XNP*DSUM-CSUM*CSUM
           DD=DSUM*ASUM-CSUM*ESUM
           ALP(1)=DD/CC
           BTA(1)=AA/CC
           RR(1)=AA/SQRT(BB*CC)
           RR(1)=RR(1)*RR(1)
      ENDIF
!
!  CALCULATE CORRELATION OVER ALL PARAMETERS
!
      IF(KHR.LE.0)THEN
           AASUM=0.0
           BBSUM=0.0
           CCSUM=0.0
           DDSUM=0.0
           EESUM=0.0
           ASUM=0.0
           BSUM=0.0
           CSUM=0.0
           DSUM=0.0
           ESUM=0.0
           DO 21 I=1,NP
           ASUM=ASUM+XSAVE(I,1,1)
           BSUM=BSUM+XSAVE(I,1,1)**2
           CSUM=CSUM+XSAVE(I,2,1)
           DSUM=DSUM+XSAVE(I,2,1)**2
  21       ESUM=ESUM+XSAVE(I,1,1)*XSAVE(I,2,1)
           DO 23 I=1,NRCALL
           ASUM=ASUM+ZSAVE(I,1,1)
           AASUM=AASUM+CSAVE(I,1)
           BSUM=BSUM+ZSAVE(I,1,1)**2
           BBSUM=BBSUM+CSAVE(I,1)**2
           CSUM=CSUM+ZSAVE(I,2,1)
           CCSUM=CCSUM+CSAVE(I,2)
           DSUM=DSUM+ZSAVE(I,2,1)**2
           DDSUM=DDSUM+CSAVE(I,2)**2
           EESUM=EESUM+CSAVE(I,1)*CSAVE(I,2)
           ESUM=ESUM+ZSAVE(I,1,1)*ZSAVE(I,2,1)
  23       CONTINUE
           ASUM=ASUM+AASUM
           BSUM=BSUM+BBSUM
           CSUM=CSUM+CCSUM
           DSUM=DSUM+DDSUM
           ESUM=ESUM+EESUM
           XCALL=FLOAT(NP+2*NRCALL)
           AA=XCALL*ESUM-ASUM*CSUM
           BB=XCALL*BSUM-ASUM*ASUM
           CC=XCALL*DSUM-CSUM*CSUM
           DD=DSUM*ASUM-CSUM*ESUM
           ALP(1)=DD/CC
           BTA(1)=AA/CC
           RR(1)=AA/SQRT(ABS(BB*CC))
           RR(1)=RR(1)*RR(1)
      ENDIF
      RETURN
      END
!
!  ***************************************************************************
!    SUBROUTINE OUTWRT---CALLED BY MAIN DURING THE ROLL CALL AND LEGISLATOR
!      PHASES ONLY.  CALCULATES NUMBER OF CONSTRAINED ROLL CALLS AND GMP
!      FOR EACH ROLL CALL DURING ROLL CALL PHASE AND, AT THE LAST ITERATION,
!      WRITES ROLL CALL ESTIMATES TO DISK.  FOR LEGISLATOR PHASE, THE GMP
!      FOR EACH LEGISLATOR IS CALCULATED AND, AT THE LAST ITERATION, WRITES
!      LEGISLATOR ESTIMATES TO DISK.  FOR BOTH ROLL CALL AND LEGISLATOR
!      PHASES, THE CPU TIME OF THE PHASE IS CALCULATED, WRITTEN OUT, AND
!      STORED.
!  ***************************************************************************
!
      SUBROUTINE OUTWRT(NRX,NSTEP,JAN,JANLST,KPJP,YBIGL,NITR,           &
     &                  NDIM,NS,NP,NRCALL,NDUAL,IPRINT,XDATA,STDDVX,    &
     &                  KAV,DYN,ZMID,STDDVZ,XSAVE,ZSAVE,CSAVE,          &
     &                  KBAD1,KBAD2,SPRM,SPRSD,GMPGMP,XFITS)
!
      DIMENSION YBIGL(NDUAL),KPJP(NDUAL,4),ALP(2),BTA(2),RR(2),         &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          STDDVX(NP,25),STDDVZ(NRCALL,2,25),XFITS(3*NS),          &
     &          XSAVE(NP,2,2),ZSAVE(NRCALL,2,2),CSAVE(NRCALL,2),        &
     &          KAV(NRCALL),GMPGMP(NP+NRCALL)
!
 1002 FORMAT(1X,2I4,I2,4X,4I3,1X,101F7.3)
 1003 FORMAT(1X,3I4,2I7,F14.4,3F9.4)
 1004 FORMAT(20A4)
 1005 FORMAT(1X,I4,1X,19A1,4I4,4F10.4)
 1014 FORMAT(1X,2I7,2G20.13)
 1015 FORMAT(1X,I5,4I6,51F9.4)
!
!  ***** ROLL CALL PHASE *****
!
      IF(NSTEP.EQ.1)THEN
           KBAD1=0
           KBAD2=0
           SPRM=0.0
           SPRSD=0.0
           SUM=0.0
           KKSUM=0
           DO 11 I=1,NRX
           KSUM=0
           KVYBAD=0
!
!  KBAD1 IS THE NUMBER OF ROLL CALLS WITH MIDPOINT VECTORS OF LENGTH 1
!
           SUMA=0.0
           SUMY=0.0
           SUMN=0.0
           DO 1 K=1,NDIM
           SUMY=SUMY+(ZMID(I,K)-DYN(I,K))**2
           SUMN=SUMN+(ZMID(I,K)+DYN(I,K))**2
           SUMA=SUMA+ZMID(I,K)**2
  1        CONTINUE
           IF(ABS(1.0-SUMA).LT..001)KBAD1=KBAD1+1
           IF(ABS(1.0-SUMA).LT..001)KVYBAD=1
!
!  KBAD2 IS THE NUMBER OF ROLL CALLS WITH YES AND NO OUTCOMES EXTERIOR
!    TO THE UNIT HYPERSPHERE
!
           IF(SUMY.GT.1.0.AND.SUMN.GT.1.0)KBAD2=KBAD2+1
           SPRM=SPRM+ABS(2.0*DYN(I,NDIM))
           SPRSD=SPRSD+(2.0*DYN(I,NDIM))**2
           DO 12 JJ=1,4
  12       KSUM=KSUM+KPJP(I,JJ)
           KKSUM=KKSUM+KSUM
           SUM=SUM+YBIGL(I)
           IF(KSUM.EQ.0)AA=99.0
           IF(KSUM.EQ.0)GO TO 16
           AA=YBIGL(I)/FLOAT(KSUM)
           AA=EXP(AA)
  16       CONTINUE
           GMPGMP(NP+I)=AA
           IF(JAN.LT.JANLST.OR.NDIM.LT.NS)GO TO 11
           LYES=KPJP(I,1)+KPJP(I,3)
           LNO =KPJP(I,2)+KPJP(I,4)
           JKP =MIN0(LYES,LNO)
           KTP =KPJP(I,2)+KPJP(I,3)
           XPRE=0.0
           IF(JKP.GT.0)THEN
              XPRE=FLOAT(JKP-KTP)/FLOAT(JKP)
           ENDIF
           IF(IPRINT.EQ.1)THEN
             WRITE(33,1002)I,KAV(I),KVYBAD,(KPJP(I,J),J=1,4),AA,XPRE,   &
     &       (DYN(I,L),ZMID(I,L),L=1,NS),((STDDVZ(I,J,M),J=1,2),M=1,NS)
           ENDIF
!
!
  11       CONTINUE
!
!   SPRM IS THE AVERAGE DISTANCE BETWEEN THE YES AND NO OUTCOMES AND
!      SPRSD IS THE STANDARD DEVIATION OF THE DISTANCE
!
           SPRSD=SQRT(FLOAT(NRX)*SPRSD-SPRM*SPRM)/                      &
     &     FLOAT(NRX)
           SPRM=SPRM/FLOAT(NRX)
      ENDIF
!
! ***** LEGISLATOR PHASE *****
!
      IF(NSTEP.EQ.2)THEN
           SUM=0.0
           KKSUM=0
           DO 13 I=1,NRX
           KSUM=0
           DO 14 JJ=1,4
  14       KSUM=KSUM+KPJP(I,JJ)
           KKSUM=KKSUM+KSUM
           SUM=SUM+YBIGL(I)
           IF(KSUM.EQ.0)AA=99.0
           IF(KSUM.EQ.0)GO TO 166
           AA=YBIGL(I)/FLOAT(KSUM)
           AA=EXP(AA)
 166       CONTINUE
           GMPGMP(I)=AA
           IF(JAN.LT.JANLST.OR.NDIM.LT.NS)GO TO 13
!
!  WRITE OUT LEGISLATOR COORDINATES--NOTE THAT THIS FORMAT IS UNIQUE
!    TO THE ROLL CALL DATA -- FOR DIFFERENT MATRICES, THIS CODE
!    WILL HAVE TO BE COMMENTED OUT -- SEE EXAMPLE BELOW
!
!
           IF(IPRINT.EQ.1)THEN
           WRITE(31,1015)I,(KPJP(I,J),J=1,4),                           &
     &        AA,(XDATA(I,L),L=1,NS),                                   &
     &        (STDDVX(I,M),M=1,NS)
           ENDIF
  13       CONTINUE
      ENDIF
!
!  CALCULATE CPU TIME OF CURRENT PHASE
!
      AA=SUM/FLOAT(KKSUM)
      AA=EXP(AA)
!
!
      CALL CORR(NP,NRCALL,XSAVE,ZSAVE,CSAVE,ALP,BTA,RR,NSTEP)
      IF(IPRINT.EQ.1)WRITE(23,1014)NRX,KKSUM,SUM,AA
      IF(NSTEP.EQ.2.AND.IPRINT.EQ.1)THEN
           WRITE(21,1003)NDIM,NITR,NSTEP,IITIME,JITIME,SUM,AA,RR(1)
      ENDIF
      IF(NSTEP.EQ.1)THEN
         XFITS(2*NS+NDIM)=AA
         IF(IPRINT.EQ.1)THEN
           WRITE(21,1003)NDIM,NITR,NSTEP,IITIME,JITIME,SUM,AA,          &
     &          RR(1),RR(2)
         ENDIF
      ENDIF
      RETURN
      END
!
!  *************************************************************************
!    SUBROUTINE NORMZ---CALLED BY MAXLIK DURING UTILITY FUNCTION PHASE.
!       CHECKS FOR SAG IN LEGISLATOR COORDINATES AND CONSTRAINS "SAGGED"
!       LEGISLATORS FROM BEING FURTHER ESTIMATED.  READJUSTS LEGISLATOR
!       COORDINATES TO LIE BETWEEN -1 AND +1.  READJUSTS ROLL CALL OUTCOME
!       COORDINATES SO THAT SUM OF DYN**2 = SUM OF ABS(XI-ZMID)*DYN.
!  *************************************************************************
!
      SUBROUTINE NORMZ(NP,NRCALL,NS,NDIM,                               &
     &                  KSMIN,KSMAX,                                    &
     &                  BBB,BBBB,SSS,ZMID,XDATA,DYN,                    &
     &                  ISENS,IPRINT)
      DIMENSION XD(NP),LMO(NP),BBB(50),BBBB(25),SSS(100),               &
     &          ZMID(NRCALL,25),XDATA(NP,25),DYN(NRCALL,25),            &
     &          ISENS(NP)
!
 1000 FORMAT(' FURTHEST LEFT LEGISLATOR=',F10.4,2X,                     &
     &'FURTHEST RIGHT LEGISLATOR=',F10.4/14X,'OLD LEFT=',F10.4,         &
     &15X,'OLD RIGHT=',F10.4/' POSITIONS: LEFT=',I5,3X,'RIGHT=',I5)
 1001 FORMAT('    ADJUSTMENT TO OUTCOME COORDINATES=',2F10.4)
 1002 FORMAT(I6,4F10.4,4I6)
 1111 FORMAT(4F10.4)
!
!  ERROR CATCH FOR NP < 20
!
      K7MOA10=10
!
      CC=XDATA(KSMIN,NDIM)
      DD=XDATA(KSMAX,NDIM)
      DO 3 I=1,NP
      XD(I)=XDATA(I,NDIM)
  3   LMO(I)=I
      CALL RSORT(XD,NP,LMO)
!
!  EXAMINE LEGISLATOR CONFIGURATION FOR SAG.  IF SAG EXISTS, CONSTRAIN
!    ESTIMATION OF LEGISLATOR BY SETTING VALUE OF ISENS(.)=0.
!
!  IF 2ND DIMENSION OR HIGHER, DO NOT INVOKE CONSTRAINTS
!
      IF(NDIM.GT.1.OR.NP.LT.20)GO TO 444
!
      DO 99 JJ=1,K7MOA10
      IF(ABS(XD(JJ)-XD(JJ+1)).GT..1)ISENS(LMO(JJ))=0
      IF(ABS(XD(NP+1-JJ)-XD(NP-JJ)).GT..1)ISENS(LMO(NP+1-JJ))=0
      IF(IPRINT.EQ.1)THEN
      WRITE(23,1002)JJ,XD(JJ),XD(JJ+1),XD(NP+1-JJ),XD(NP-JJ),           &
     &ISENS(LMO(JJ)),ISENS(LMO(JJ+1)),ISENS(LMO(NP+1-JJ)),              &
     &ISENS(LMO(NP-JJ))
      ENDIF
      IF(ISENS(LMO(JJ)).EQ.1.AND.ISENS(LMO(NP+1-JJ)).EQ.1)GO TO 999
  99  CONTINUE
  999 CONTINUE
!
  444 CONTINUE
!
!
!  SET MIN AND MAX LEGISLATORS TO FURTHEST LEFT/RIGHT UNCONSTRAINED
!    LEGISLATORS
!
      KSMIN=LMO(1)
      KSMAX=LMO(NP)
      DO 77 J=1,K7MOA10
      IF(ISENS(LMO(J)).EQ.0)THEN
         KSMIN=LMO(J+1)
      ENDIF
      IF(ISENS(LMO(NP+1-J)).EQ.0)THEN
         KSMAX=LMO(NP+1-J-1)
      ENDIF
  77  CONTINUE
!
!  IF 2ND DIMENSION OR HIGHER, DO NOT INVOKE CONSTRAINTS
!
      IF(NDIM.GT.1)GO TO 445
!
      AA=(XDATA(KSMIN,NDIM)+XDATA(KSMAX,NDIM))/2.0
      BB=XDATA(KSMAX,NDIM)-AA
      IF(IPRINT.EQ.1)THEN
      WRITE(23,1111)AA,BB,XDATA(KSMIN,NDIM),XDATA(KSMAX,NDIM)
      ENDIF
!
!  SET MIN LEGISLATOR VALUE TO -1 AND MAX LEGISLATOR VALUE TO +1
!    AND NORMALIZE COORDINATES TO RUN FROM -1 TO +1
!
      DO 1 I=1,NP
      IF(ISENS(I).EQ.0.AND.XDATA(I,NDIM).LT.0.0)XDATA(I,NDIM)=-1.0
      IF(ISENS(I).EQ.0.AND.XDATA(I,NDIM).GT.0.0)XDATA(I,NDIM)= 1.0
      IF(ISENS(I).EQ.0)GO TO 1
      XDATA(I,NDIM)=(XDATA(I,NDIM)-AA)/BB
      IF(XDATA(I,NDIM).GT.1.0)XDATA(I,NDIM)=1.0
      IF(XDATA(I,NDIM).LT.-1.0)XDATA(I,NDIM)=-1.0
  1   CONTINUE
!
!  APPLY LINEAR TRANSFORMATION TO ROLL CALL COORDINATES
!
      DO 2 J=1,NRCALL
      ZMID(J,NDIM)=(ZMID(J,NDIM)-AA)/BB
  2   CONTINUE
!
!  MULTIPLY WEIGHT BY SCALING FACTOR SO AS TO ACHIEVE CANCELLATION
!    IN EXPONENT OF UTILITY FUNCTION
!
      BBB(2)=BBB(2)*BB
      BBBB(1)=BBB(2)
!
!  SET NEW MIN AND MAX VALUES
!
!
  445 CONTINUE
!
      SSS(1)=XDATA(KSMIN,NDIM)
      SSS(2)=XDATA(KSMAX,NDIM)
      SSS(3)=XDATA(LMO(2),NDIM)
      SSS(4)=XDATA(LMO(NP-1),NDIM)
!
      IF(IPRINT.EQ.1)THEN
      WRITE(23,1000)XDATA(KSMIN,NDIM),XDATA(KSMAX,NDIM),CC,DD,          &
     &    KSMIN,KSMAX
      ENDIF
      RETURN
      END
!
! **************************************************************************
!
!  SUBROUTINE SVDSVD -- PERFORMS SINGULAR VALUE DECOMPOSITION
!                       
!                       
!
! **************************************************************************
!
      SUBROUTINE SVDSVD(NROW,NCOL,Y16MIDP,YHAT,UUU,VVV,IRANK,IPRINT)
      DIMENSION Y16MIDP(50,50),UUU(50,50),VVV(50,50),YHAT(50),          &
     &          XCOV(50,50),WVEC(50),ZMAT(50,50),FV1(50),FV2(50)
!
  101 FORMAT(' PERFORMANCE INDEX EIGENVALUE/VECTOR ROUTINE='            &
     &       ,3I5,I6)
  102 FORMAT(' ACCURACY OF SVD',2F10.4)
  103 FORMAT(10F7.3)
 1001 FORMAT(10F10.4)
!   
      DO 3 I=1,NROW
      DO 2 L=1,NROW
      SUM=0.0
      DO 1 J=1,NCOL
      SUM=SUM+Y16MIDP(I,J)*Y16MIDP(L,J)
  1   CONTINUE
      XCOV(I,L)=SUM
  2   CONTINUE
  3   CONTINUE
!
!  GET U MATRIX
!
      CALL KPRS(50,NROW,XCOV,WVEC,1,ZMAT,FV1,FV2,IER)
      IF(IPRINT.EQ.1)WRITE(23,101)IER
      DO 4 I=1,NROW
      DO 5 J=1,NROW
      UUU(I,J)=ZMAT(I,NROW+1-J)
  5   CONTINUE
  4   CONTINUE
!
      DO 6 I=1,NCOL
      DO 7 L=1,NCOL
      SUM=0.0
      DO 8 J=1,NROW
      SUM=SUM+Y16MIDP(J,I)*Y16MIDP(J,L)
  8   CONTINUE
      XCOV(I,L)=SUM
  7   CONTINUE
  6   CONTINUE
!
!  GET V MATRIX
!
      CALL KPRS(50,NCOL,XCOV,WVEC,1,ZMAT,FV1,FV2,IER)
      IF(IPRINT.EQ.1)WRITE(23,101)IER
      IRANK=0
      DO 9 I=1,NCOL
      YHAT(I)=SQRT(ABS(WVEC(NCOL+1-I)))
      IF(YHAT(I).GE..001)IRANK=IRANK+1
      DO 10 J=1,NCOL
      VVV(I,J)=ZMAT(I,NCOL+1-J)
  10  CONTINUE
  9   CONTINUE
!
!  FORM U'XV TO GET SIGNS
!
      DO 16 L=1,NCOL
      DO 18 J=1,NCOL
      SUM=0.0
      DO 17 I=1,NROW
      SUM=SUM+UUU(I,J)*Y16MIDP(I,L)
  17  CONTINUE
      ZMAT(J,L)=SUM
  18  CONTINUE
  16  CONTINUE
      DO 19 K=1,NCOL
      DO 20 J=1,NCOL
      SUM=0.0
      DO 21 L=1,NCOL
      SUM=SUM+ZMAT(J,L)*VVV(L,K)
  21  CONTINUE
      XCOV(J,K)=SUM
  20  CONTINUE
  19  CONTINUE
      DO 22 I=1,NCOL
      WVEC(I)=1.0
      IF(XCOV(I,I).LT.0.0)WVEC(I)=-1.0
      IF(IPRINT.EQ.1)WRITE(23,103)(XCOV(I,J),J=1,NCOL)
  22  CONTINUE
!
!  DO SIGN CHANGE
!
      DO 23 I=1,NROW
      DO 24 J=1,NCOL
      UUU(I,J)=WVEC(J)*UUU(I,J)
  24  CONTINUE
  23  CONTINUE
!
!  CHECK ACCURACY OF DECOMPOSITION
!
      ESUM=0.0
      DO 11 L=1,NCOL
      FSUM=0.0
      GSUM=0.0
      DO 13 I=1,NROW
      SUM=0.0
      DO 12 J=1,NCOL
      SUM=SUM+UUU(I,J)*YHAT(J)*VVV(L,J)
  12  CONTINUE
      ESUM=ESUM+(Y16MIDP(I,L)-SUM)**2
      XCOV(I,L)=SUM
  13  CONTINUE
  11  CONTINUE
!      DO 15 I=1,NROW
!      WRITE(23,103)(Y16MIDP(I,J),J=1,NCOL),(XCOV(I,J),J=1,NCOL),
!     C             (UUU(I,J),J=1,NCOL),(VVV(I,J),J=1,NCOL),YHAT(I)
!  15  CONTINUE
!
      IF(IPRINT.EQ.1)WRITE(23,102)ESUM,ESUM/FLOAT(NROW*NCOL)
!
      RETURN
      END
!
! *********************************************************************
!   SUBROUTINE ROTATE -- IMPLEMENTS SCHONEMANN'S ORTHOGONAL PROCRUSTES
!                             ROTATION SOLUTION
!
!                          XMAT(,) IS ROTATED TO BEST MATCH XTRUE(,)
!
! *********************************************************************
!
      SUBROUTINE ROTATE(NP,NS,XMAT,XTRUE,IPRINT)
      DIMENSION XMAT(NP,25),XTRUE(NP,25),Y16MIDP(50,50),                &
     &          YHAT(50),UUU(50,50),VVV(50,50),VCOV(50,50),             &
     &          FV1(50),FV2(50),XMATNEW(50,50)
      REAL, ALLOCATABLE :: XSAVE(:,:)
      ALLOCATE ( XSAVE(NP,25))
 1015 FORMAT(' SINGULAR VALUE DECOMPOSITION',I5)
 1016 FORMAT(I5,99F7.3)
 1017 FORMAT(' ORTHOGONAL PROCRUSTES ROTATION MATRIX')
 1090 FORMAT(' R-SQUARES TRUE VS. REPRODUCED',2I4,2F7.3)
!
!  ROTATE ESTIMATED COORDINATES TO BEST FIT TRUE COORDINATES USING
!    SCHONEMANN'S ORTHOGONAL PROCRUSTES SOLUTION
!
      DO 376 L=1,NS
      DO 377 K=1,NS
      SUM=0.0
      DO 378 I=1,NP
      SUM=SUM+XMAT(I,K)*XTRUE(I,L)
  378 CONTINUE
      Y16MIDP(K,L)=SUM
  377 CONTINUE
  376 CONTINUE
!
!    CALL SINGULAR VALUE DECOMPOSITION ROUTINE
!
      XTOL=.001
!      CALL LSVRR(NS,NS,Y16MIDP,461,21,XTOL,IRANK,YHAT,UUU,
!     C           461,VVV,25)
      CALL SVDSVD(NS,NS,Y16MIDP,YHAT,UUU,VVV,IRANK,IPRINT)
!
      IF(IPRINT.EQ.1)WRITE(23,1015)IRANK
!
!  COMPUTE ROTATION MATRIX
!
      IF(IPRINT.EQ.1)WRITE(23,1017)
      DO 379 K=1,NS
      DO 380 M=1,NS
      SUM=0.0
      DO 381 L=1,NS
      SUM=SUM+UUU(K,L)*VVV(M,L)
  381 CONTINUE
      VCOV(K,M)=SUM
  380 CONTINUE
      IF(IPRINT.EQ.1)WRITE(23,1016)K,YHAT(K),(VCOV(K,L),L=1,NS) 
  379 CONTINUE
!
!  ROTATE OBSERVED COORDINATES
!
      DO 382 I=1,NP
      DO 383 K=1,NS
      SUM=0.0
      DO 384 L=1,NS
      SUM=SUM+XMAT(I,L)*VCOV(L,K)
  384 CONTINUE
      XSAVE(I,K)=SUM
  383 CONTINUE
  382 CONTINUE
      DO 385 I=1,NP
      DO 386 K=1,NS
      XMAT(I,K)=XSAVE(I,K)
  386 CONTINUE
!
!      IF(IPRINT.EQ.1)WRITE(23,1016)I,(XMAT(I,K),K=1,NS),
!     C                               (XTRUE(I,K),K=1,NS)
  385 CONTINUE
!
!  CALCULATE R-SQUARES FOR ROTATED MATRICES
!
      DO 275 JJ=1,50
      FV1(JJ)=0.0
      FV2(JJ)=0.0
      DO 276 JK=1,25
      XMATNEW(JJ,JK)=0.0
  276 CONTINUE
  275 CONTINUE
      DO 270 I=1,NP
      DO 271 K=1,NS
      FV1(K)=FV1(K)+XTRUE(I,K)
      FV2(K)=FV2(K)+XMAT(I,K)
      FV1(K+NS)=FV1(K+NS)+XTRUE(I,K)**2
      FV2(K+NS)=FV2(K+NS)+XMAT(I,K)**2
      DO 272 L=1,NS
      XMATNEW(K,L)=XMATNEW(K,L)+XTRUE(I,K)*XMAT(I,L)
  272 CONTINUE
  271 CONTINUE
  270 CONTINUE
      DO 273 K=1,NS
      DO 274 L=1,NS
      AA=FLOAT(NP)*XMATNEW(K,L)-FV1(K)*FV2(L)
      AB=FLOAT(NP)*FV1(K+NS)-FV1(K)*FV1(K)
      AC=FLOAT(NP)*FV2(L+NS)-FV2(L)*FV2(L)
      IF(AB*AC.GT.0.0)RSQR1=(AA*AA)/(AB*AC)
      IF(IPRINT.EQ.1)WRITE(23,1090)K,L,RSQR1
  274 CONTINUE
  273 CONTINUE
      DEALLOCATE (XSAVE)
      RETURN
      END
