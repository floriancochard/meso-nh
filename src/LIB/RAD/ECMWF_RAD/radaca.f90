!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
!OPTIONS XOPT(HSFUN)
SUBROUTINE RADACA ( KIDIA , KFDIA , KLON , KTDIA , KLEV &
                  &, PAPRS , PGELAM, PSIN  , PCLON, PSLON , PTH &
                  &, PAER  , POZON                         )

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *RADACA  - COMPUTES DISTRIBUTION OF AEROSOLS AND OZONE
!                >>> FOR ASSIMILATION PURPOSES <<

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *RADACA* FROM *RADINA* 

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------

!     J.-F. MAHFOUF    E.C.M.W.F.    97/07/04

!     Aadapted from

!     J.-J. MORCRETTE  E.C.M.W.F.    91/03/15

!     MODIFICATIONS.
!     --------------
!     J.-J. MORCRETTE  E.C.M.W.F.    93/03/15   OPERATIONAL CLIMATOLOGY
!     JJMorcrette  99-05-25     Revised aerosols
!     JJMorcrette 98-12-21 GISS volcanic aerosol climatology
!     JJMorcrette 99-09    monthly climatology of tropospheric aerosols
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOMCST   , ONLY : R        ,RPI
USE YOEAERD  , ONLY : CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED   ,&
            &RCAEOPS  ,RCAEOPL  ,RCAEOPU  ,RCAEOPD  ,RCTRBGA  ,&
            &RCVOBGA  ,RCSTBGA  ,RCTRPT   ,RAESC    ,RAESS    ,&
            &RAELC    ,RAELS    ,RAEUC    ,RAEUS    ,RAEDC    ,&
            &RAEDS
USE YOEOZOC  , ONLY : COZQC    ,COZQS    ,COZHC    ,COZHS
USE OYOERAD   , ONLY : LHVOLCA  ,LNEWAER
USE YOEAERC  , ONLY : RSINCT   ,RSINCV   ,REPAER   ,&
            &RTAEBC  ,RTAEOR   ,RTAESD   ,RTAESS   ,RTAESU   , &
            &RTAEVO 


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KTDIA




!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

REAL_B :: PAPRS(KLON,KLEV+1), PGELAM(KLON), PSIN(KLON) &
  &, PCLON(KLON), PSLON(KLON), PTH  (KLON,KLEV+1)

REAL_B :: PAER(KLON,6,KLEV),POZON(KLON,KLEV)
!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------

INTEGER_M :: IINLA1(KLON), IINLA2(KLON)
INTEGER_M :: IINLO1(KLON), IINLO2(KLON)

REAL_B :: ZAED  (KLON), ZAEL  (KLON), ZAES  (KLON), ZAEU  (KLON)
REAL_B :: ZAEQDN(KLON), ZAEQDO(KLON), ZAEQLN(KLON), ZAEQLO(KLON)
REAL_B :: ZAEQSN(KLON), ZAEQSO(KLON), ZAEQUN(KLON), ZAEQUO(KLON)

REAL_B :: ZAERBC(KLON), ZAEROR(KLON), ZAERSD(KLON)
REAL_B :: ZAERSS(KLON), ZAERSU(KLON), ZAERVO(KLON)

REAL_B :: ZAETRN(KLON),ZAETRO(KLON)

REAL_B :: ZALP(66)
REAL_B :: ZDPN(KLON)  , ZDPO(KLON)
REAL_B :: ZFAED(21)    , ZFAEL(21)    , ZFAES(21)    , ZFAEU(21)
REAL_B :: ZFOZQ(11)    , ZFOZH(11)
REAL_B :: ZGRTH(KLON)
REAL_B :: ZLON(KLON)  , ZLONR(72)    , ZNLO1(KLON) , ZNLO2(KLON)

REAL_B :: ZOZH (KLON) , ZOZQ (KLON)
REAL_B :: ZQOFN(KLON) , ZQOFO(KLON)
REAL_B :: ZSILAT(KLON), ZSINR(46)


!REAL_B :: ZALP(66)
!REAL_B :: ZFOZQ(11),ZFOZH(11),ZFAES(21),ZFAEL(21),  ZFAEU(21),ZFAED(21)
!
!REAL_B :: ZDPN (KLON) ,ZDPO (KLON) ,ZQOFN(KLON) ,ZQOFO(KLON)
!REAL_B :: ZOZH (KLON) ,ZOZQ (KLON)
!
!REAL_B :: ZAED  (KLON), ZAEL  (KLON), ZAES  (KLON), ZAEU  (KLON)
!REAL_B :: ZAEQSN(KLON),ZAEQSO(KLON),ZAEQLN(KLON),ZAEQLO(KLON)
!REAL_B :: ZAEQUN(KLON),ZAEQUO(KLON),ZAEQDN(KLON),ZAEQDO(KLON)
!
!REAL_B :: ZAERBC(KLON), ZAEROR(KLON), ZAERSD(KLON)
!REAL_B :: ZAERSS(KLON), ZAERSU(KLON), ZAERVO(KLON)
!
!REAL_B :: ZAETRN(KLON),ZAETRO(KLON)
!
!REAL_B :: ZSILAT(KLON), ZSINR(46)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IL, IMM, IMNC, IMNS, INLA, INLA1, INLA2, INLO1, INLO2, &
          &ITOTPT, JK, JL, JMM, JNN, NLATR, NLONR, JAER, JEND, &
          &JIL, JJL, IPRINT, ITOT, KRINT, KSHIFT, JLR, KCF

!     LOCAL REAL SCALARS
REAL_B :: ZAETR, ZCOS1, ZCOS10, ZCOS2, ZCOS3, ZCOS4,&
          &ZCOS5, ZCOS6, ZCOS7, ZCOS8, ZCOS9, ZCPHN3, &
          &ZCPHO3, ZDPNMO, ZGRIDR, ZLATR, ZSDPN3, ZSDPO3, &
          &ZSIN, ZSIN1, ZSIN10, ZSIN2, ZSIN3, ZSIN4, &
          &ZSIN5, ZSIN6, ZSIN7, ZSIN8, ZSIN9
REAL_B :: ZAERBC1, ZAERBC2, ZAEROR1, ZAEROR2, ZAERSD1, ZAERSD2, &
          &ZAERSS1, ZAERSS2, ZAERSU1, ZAERSU2

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

!*         1.     "NEW AEROSOL DISTRIBUTION" PARAMETERS COMPUTATIONS
!                 --------------------------------------------------


!*         1.1    VOLCANIC AEROSOL DISTRIBUTION PARAMETERS
!                 ----------------------------------------

!                 GISS CLIMATOLOGY 
!                 ----------------

KSHIFT=0
KRINT=1
KCF=0

IF (LHVOLCA) THEN
  NLATR=46
  ZGRIDR=180._JPRB/(NLATR-1)
  DO JLR=1,NLATR
    ZLATR=90._JPRB-(JLR-1)*ZGRIDR
    ZSINR(JLR)=SIN(ZLATR*RPI/180._JPRB)
  ENDDO

  IL=KSHIFT
  DO JL=KIDIA,KFDIA,KRINT
    IL=IL+1
    INLA=0
    ZSILAT(IL)=-9999._JPRB
    ZSIN=PSIN(JL)
    DO JLR=1,NLATR-1
      IF (ZSIN <= ZSINR(JLR) .AND. ZSIN > ZSINR(JLR+1)) THEN
        INLA=JLR
        ZSILAT(IL)=(ZSIN-ZSINR(INLA))/(ZSINR(INLA+1)-ZSINR(INLA))
        ZAERVO(IL)=RTAEVO(INLA)+ZSILAT(IL)*(RTAEVO(INLA+1)-RTAEVO(INLA))
      ENDIF
    ENDDO
    IF (ZSIN <= ZSINR(NLATR-1) .AND. ZSIN >= ZSINR(NLATR))THEN
      INLA=NLATR
      ZSILAT(IL)=(ZSIN-ZSINR(INLA-1))/(ZSINR(INLA)-ZSINR(INLA-1))
      ZAERVO(IL)=RTAEVO(INLA-1)&
       &+ZSILAT(IL)*(RTAEVO(INLA)-RTAEVO(INLA-1))
    ENDIF
    IF (INLA == 0) THEN
!      CALL ABOR1(' Problem with lat. interpolation in RADACA!')
      STOP ' Problem with lat. interpolation in RADACA!'
    ENDIF
  ENDDO

!                 TANRE ET AL. CLIMATOLOGY 
!                 ------------------------

ELSE
  IL = KSHIFT
  DO JL=KIDIA,KFDIA,KRINT
    IL = IL+1
    ZAERVO(IL)=RCVOBGA
  ENDDO
ENDIF
ITOTPT=IL



!*         1.2    TROPOSPHERIC AEROSOL DISTRIBUTION PARAMETERS
!                 --------------------------------------------

IF (LNEWAER) THEN
!  print *,'LNEWAER= ',LNEWAER
!  print *,'              SINLAT             LONGITUDE'
!  DO JL=KIDIA,KFDIA,KRINT
!    print 9001,JL,PSIN(JL),PGELAM(JL)
!9001    format(1x,'RADACA ',1x,I5,1x,2E15.8)
!  END DO         

!-- latitude       
  NLATR=46
  ZGRIDR=180._JPRB/(NLATR-1)
  DO JLR=1,NLATR
    ZLATR=90._JPRB-(JLR-1)*ZGRIDR
    ZSINR(JLR)=SIN(ZLATR*RPI/180._JPRB)
  END DO
  NLONR=72
  DO JLR=1,NLONR
    ZLONR(JLR)=(JLR-1)*2._JPRB*RPI/NLONR
  END DO  
        
!  print 9121,(ZSINR(JLR),JLR=1,NLATR)
!9121     format(1x,'ZSINR ',8E15.7)
!  print 9122,(ZLONR(JLR),JLR=1,NLONR)
!9122     format(1x,'ZLONR ',8E15.7)
           
           
  IL=KSHIFT
  DO JL=KIDIA,KFDIA,KRINT         
    IL=IL+1
    IINLA1(IL)=0
    IINLA2(IL)=0
    ZSILAT(IL)=-9999._JPRB
    ZSIN=PSIN(JL)
    DO JLR=1,NLATR-1
      IF (ZSIN <= ZSINR(JLR) .AND. ZSIN > ZSINR(JLR+1)) THEN
        INLA=JLR
        IINLA1(IL)=JLR
        IINLA2(IL)=JLR+1
        ZSILAT(IL)=(ZSIN-ZSINR(INLA))/(ZSINR(INLA+1)-ZSINR(INLA))
      ENDIF
    ENDDO
    IF (ZSIN <= ZSINR(NLATR-1) .AND. ZSIN >= ZSINR(NLATR))THEN
      INLA=NLATR
      IINLA1(IL)=NLATR-1
      IINLA2(IL)=NLATR
      ZSILAT(IL)=(ZSIN-ZSINR(INLA-1))/(ZSINR(INLA)-ZSINR(INLA-1))
    END IF    
    IF (INLA.EQ.0) THEN 
!      CALL ABOR1(' Problem with lat. interpolation in RADACA!')
      STOP ' Problem with lat. interpolation in RADACA!'
    ENDIF
!    print 9123,JL,IL,PSIN(JL),INLA,ZSINR(INLA),ZSILAT(IL)
!9123     format(1x,'Interp.Latit.',2I4,F10.7,I4,2F10.7)          
  END DO  

!-- longitude
  IL=KSHIFT
  DO JL=KIDIA,KFDIA,KRINT
    IL=IL+1
    IINLO1(IL)=0
    IINLO2(IL)=0
    ZLON(IL)=-9999.
    DO JLR=1,71
      IF (PGELAM(JL) < ZLONR(JLR+1) .AND. PGELAM(JL) >= ZLONR(JLR)) &
     &  THEN
        IINLO1(IL)=JLR
        IINLO2(IL)=JLR+1
        ZNLO1(IL)=ZLONR(JLR)
        ZNLO2(IL)=ZLONR(JLR+1)
      END IF
    END DO      
    IF (PGELAM(JL) >= ZLONR(72)) THEN
      IINLO1(IL)=72
      IINLO2(IL)= 1
      ZNLO1(IL)=ZLONR(72)
      ZNLO2(IL)=ZLONR(72)+2.*RPI
    ENDIF
!    print 9124,JL,IL,PGELAM(JL),IINLO1(IL),IINLO2(IL) &
!     & ,ZNLO1(IL),ZNLO2(IL)
!9124     format(1x,'Interp.Longit0.',2I4,F10.7,2I5,2F10.7)          
  END DO  
 
  IL=KSHIFT
  DO JL=KIDIA,KFDIA,KRINT        
    IL=IL+1
    IF (IINLO1(IL).EQ.0 .OR. IINLO2(IL).EQ.0) THEN 
!      CALL ABOR1(' Problem with long. interpolation in RADACA!')
      STOP ' Problem with long. interpolation in RADACA!'
    ENDIF
    ZLON(IL)=(PGELAM(JL)-ZNLO1(IL))/(ZNLO2(IL)-ZNLO1(IL))
    INLO1=IINLO1(IL)
    INLO2=IINLO2(IL)
    INLA1=IINLA1(IL)
    INLA2=IINLA2(IL)
    
    ZAERBC1=RTAEBC(INLA1,INLO1) &
     &      +ZSILAT(IL)*(RTAEBC(INLA2,INLO1)-RTAEBC(INLA1,INLO1))  
    ZAERBC2=RTAEBC(INLA1,INLO2) &
     &      +ZSILAT(IL)*(RTAEBC(INLA2,INLO2)-RTAEBC(INLA1,INLO2))  
    ZAERBC(IL)=ZAERBC1+ZLON(IL)*(ZAERBC2-ZAERBC1)
        
    ZAEROR1=RTAEOR(INLA1,INLO1) &
     &      +ZSILAT(IL)*(RTAEOR(INLA2,INLO1)-RTAEOR(INLA1,INLO1))  
    ZAEROR2=RTAEOR(INLA1,INLO2) &
     &      +ZSILAT(IL)*(RTAEOR(INLA2,INLO2)-RTAEOR(INLA1,INLO2))  
    ZAEROR(IL)=ZAEROR1+ZLON(IL)*(ZAEROR2-ZAEROR1)
        
    ZAERSD1=RTAESD(INLA1,INLO1) &
     &      +ZSILAT(IL)*(RTAESD(INLA2,INLO1)-RTAESD(INLA1,INLO1))  
    ZAERSD2=RTAESD(INLA1,INLO2) &
     &      +ZSILAT(IL)*(RTAESD(INLA2,INLO2)-RTAESD(INLA1,INLO2))  
    ZAERSD(IL)=ZAERSD1+ZLON(IL)*(ZAERSD2-ZAERSD1)
        
    ZAERSS1=RTAESS(INLA1,INLO1) &
     &      +ZSILAT(IL)*(RTAESS(INLA2,INLO1)-RTAESS(INLA1,INLO1))  
    ZAERSS2=RTAESS(INLA1,INLO2) &
     &      +ZSILAT(IL)*(RTAESS(INLA2,INLO2)-RTAESS(INLA1,INLO2))  
    ZAERSS(IL)=ZAERSS1+ZLON(IL)*(ZAERSS2-ZAERSS1)
        
    ZAERSU1=RTAESU(INLA1,INLO1) &
     &      +ZSILAT(IL)*(RTAESU(INLA2,INLO1)-RTAESU(INLA1,INLO1))  
    ZAERSU2=RTAESU(INLA1,INLO2) &
     &      +ZSILAT(IL)*(RTAESU(INLA2,INLO2)-RTAESU(INLA1,INLO2))  
    ZAERSU(IL)=ZAERSU1+ZLON(IL)*(ZAERSU2-ZAERSU1)
          
!    print 9125,JL,IL,PSIN(JL),PGELAM(JL),ZSILAT(IL) &
!     &      ,RTAESU(INLA2,INLO1),RTAESU(INLA1,INLO1),ZAERSU1 &
!     &      ,RTAESU(INLA2,INLO2),RTAESU(INLA1,INLO2),ZAERSU2 &
!     &                    ,INLA1,INLA2,INLO1,INLO2
!9125     format(1x,'Interp.Longit1.',2I4,9F10.7,4I5)          
!    print 9126,JL,IL,PSIN(JL),PGELAM(JL),ZSILAT(IL),ZLON(IL) &
!     &      ,ZNLO1(IL),ZNLO2(IL),INLA1,INLA2,INLO1,INLO2
!9126     format(1x,'Interp.Longit2.',2I4,6F10.7,4I5)          
!    print 9127,JL,IL,ZAERBC(IL),ZAEROR(IL),ZAERSD(IL),ZAERSS(IL) &
!     &      ,ZAERSU(IL)
!9127     format(1x,'Interp.Longit3.',2I4,5F10.7)          
  END DO       
END IF

!     ------------------------------------------------------------------   

!*       2.     OZONE
!               -----

ZSIN=PSIN(KIDIA)

!*       2.1     CALL TO LEGTRI.
!                ---------------
!***
CALL LEGTRI (ZSIN,6,66,ZALP)
!***

!*       2.2     LEGENDRE TRANSFORM FOR OZONE.
!                -----------------------------

DO JMM=1,11
  ZFOZQ(JMM)=_ZERO_
  ZFOZH(JMM)=_ZERO_
ENDDO
IMM=0
IMNC=0
IMNS=0
DO JMM=1,6
  IMM=IMM+1
  DO JNN=JMM,6
    IMNC=IMNC+1
    ZFOZQ(IMM)=ZFOZQ(IMM)+ZALP(IMNC)*COZQC(IMNC)
    ZFOZH(IMM)=ZFOZH(IMM)+ZALP(IMNC)*COZHC(IMNC)
  ENDDO
  IF(JMM /= 1) THEN
    IMM=IMM+1
    DO JNN=JMM,6
      IMNS=IMNS+1
      ZFOZQ(IMM)=ZFOZQ(IMM)+ZALP(IMNS+6)*COZQS(IMNS)
      ZFOZH(IMM)=ZFOZH(IMM)+ZALP(IMNS+6)*COZHS(IMNS)
    ENDDO
  ENDIF
ENDDO


!*       2.3     FOURIER TRANSFORM FOR OZONE.
!                ----------------------------

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  ZCOS1=PCLON(JL)
  ZSIN1=PSLON(JL)
  ZCOS2=ZCOS1*ZCOS1-ZSIN1*ZSIN1
  ZSIN2=ZSIN1*ZCOS1+ZCOS1*ZSIN1
  ZCOS3=ZCOS2*ZCOS1-ZSIN2*ZSIN1
  ZSIN3=ZSIN2*ZCOS1+ZCOS2*ZSIN1
  ZCOS4=ZCOS3*ZCOS1-ZSIN3*ZSIN1
  ZSIN4=ZSIN3*ZCOS1+ZCOS3*ZSIN1
  ZCOS5=ZCOS4*ZCOS1-ZSIN4*ZSIN1
  ZSIN5=ZSIN4*ZCOS1+ZCOS4*ZSIN1
  ZOZQ(IL)=&
   &ZFOZQ(1)+_TWO_*(ZFOZQ(2)*ZCOS1+ZFOZQ(3)*ZSIN1+ZFOZQ(4)*ZCOS2 &
   &+ZFOZQ(5)*ZSIN2+ZFOZQ(6)*ZCOS3+ZFOZQ(7)*ZSIN3+ZFOZQ(8)&
   &*ZCOS4+ZFOZQ(9)*ZSIN4+ZFOZQ(10)*ZCOS5+ZFOZQ(11)*ZSIN5)
  ZOZH(IL)=&
   &ZFOZH(1)+_TWO_*(ZFOZH(2)*ZCOS1+ZFOZH(3)*ZSIN1+ZFOZH(4)*ZCOS2 &
   &+ZFOZH(5)*ZSIN2+ZFOZH(6)*ZCOS3+ZFOZH(7)*ZSIN3+ZFOZH(8)&
   &*ZCOS4+ZFOZH(9)*ZSIN4+ZFOZH(10)*ZCOS5+ZFOZH(11)*ZSIN5)
  ZOZH(IL)=SQRT(ZOZH(IL))**3
ENDDO

!     ------------------------------------------------------------------

!       3.     AEROSOLS
!              --------
!***
!       3.1     CALL TO LEGTRI

!***
CALL LEGTRI (ZSIN,11,66,ZALP)
!***


!       3.2     LEGENDRE TRANSFORM FOR AEROSOLS
!               -------------------------------

DO JMM=1,21
  ZFAES(JMM) = _ZERO_
  ZFAEL(JMM) = _ZERO_
  ZFAEU(JMM) = _ZERO_
  ZFAED(JMM) = _ZERO_
ENDDO
IMM  = 0
IMNC = 0
IMNS = 0
DO JMM=1,11
  IMM  = IMM+1
  DO JNN=JMM,11
    IMNC = IMNC+1
    ZFAES(IMM) = ZFAES(IMM)+ZALP(IMNC)*RAESC(IMNC)
    ZFAEL(IMM) = ZFAEL(IMM)+ZALP(IMNC)*RAELC(IMNC)
    ZFAEU(IMM) = ZFAEU(IMM)+ZALP(IMNC)*RAEUC(IMNC)
    ZFAED(IMM) = ZFAED(IMM)+ZALP(IMNC)*RAEDC(IMNC)
  ENDDO
  IF(JMM /= 1) THEN
    IMM  = IMM+1
    DO JNN=JMM,11
      IMNS = IMNS+1
      ZFAES(IMM) = ZFAES(IMM)+ZALP(IMNS+11)*RAESS(IMNS)
      ZFAEL(IMM) = ZFAEL(IMM)+ZALP(IMNS+11)*RAELS(IMNS)
      ZFAEU(IMM) = ZFAEU(IMM)+ZALP(IMNS+11)*RAEUS(IMNS)
      ZFAED(IMM) = ZFAED(IMM)+ZALP(IMNS+11)*RAEDS(IMNS)
    ENDDO
  ENDIF
ENDDO

!       3.3     FOURIER TRANSFORM FOR AEROSOLS
!               ------------------------------

IL = KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL = IL+1
  ZCOS1    = PCLON(JL)
  ZSIN1    = PSLON(JL)
  ZCOS2    = ZCOS1*ZCOS1-ZSIN1*ZSIN1
  ZSIN2    = ZSIN1*ZCOS1+ZCOS1*ZSIN1
  ZCOS3    = ZCOS2*ZCOS1-ZSIN2*ZSIN1
  ZSIN3    = ZSIN2*ZCOS1+ZCOS2*ZSIN1
  ZCOS4    = ZCOS3*ZCOS1-ZSIN3*ZSIN1
  ZSIN4    = ZSIN3*ZCOS1+ZCOS3*ZSIN1
  ZCOS5    = ZCOS4*ZCOS1-ZSIN4*ZSIN1
  ZSIN5    = ZSIN4*ZCOS1+ZCOS4*ZSIN1
  ZCOS6    = ZCOS5*ZCOS1-ZSIN5*ZSIN1
  ZSIN6    = ZSIN5*ZCOS1+ZCOS5*ZSIN1
  ZCOS7    = ZCOS6*ZCOS1-ZSIN6*ZSIN1
  ZSIN7    = ZSIN6*ZCOS1+ZCOS6*ZSIN1
  ZCOS8    = ZCOS7*ZCOS1-ZSIN7*ZSIN1
  ZSIN8    = ZSIN7*ZCOS1+ZCOS7*ZSIN1
  ZCOS9    = ZCOS8*ZCOS1-ZSIN8*ZSIN1
  ZSIN9    = ZSIN8*ZCOS1+ZCOS8*ZSIN1
  ZCOS10   = ZCOS9*ZCOS1-ZSIN9*ZSIN1
  ZSIN10   = ZSIN9*ZCOS1+ZCOS9*ZSIN1
  ZAES(IL) = ZFAES(1) + _TWO_*&
   &( ZFAES(2)*ZCOS1  + ZFAES(3)*ZSIN1  + ZFAES(4)*ZCOS2 &
   &+ ZFAES(5)*ZSIN2  + ZFAES(6)*ZCOS3  + ZFAES(7)*ZSIN3 &
   &+ ZFAES(8)*ZCOS4  + ZFAES(9)*ZSIN4  + ZFAES(10)*ZCOS5 &
   &+ ZFAES(11)*ZSIN5 + ZFAES(12)*ZCOS6 + ZFAES(13)*ZSIN6 &
   &+ ZFAES(14)*ZCOS7 + ZFAES(15)*ZSIN7 + ZFAES(16)*ZCOS8 &
   &+ ZFAES(17)*ZSIN8 + ZFAES(18)*ZCOS9 + ZFAES(19)*ZSIN9 &
   &+ ZFAES(20)*ZCOS10+ ZFAES(21)*ZSIN10                 )
  ZAEL(IL) = ZFAEL(1) + _TWO_*&
   &( ZFAEL(2)*ZCOS1  + ZFAEL(3)*ZSIN1  + ZFAEL(4)*ZCOS2 &
   &+ ZFAEL(5)*ZSIN2  + ZFAEL(6)*ZCOS3  + ZFAEL(7)*ZSIN3 &
   &+ ZFAEL(8)*ZCOS4  + ZFAEL(9)*ZSIN4  + ZFAEL(10)*ZCOS5 &
   &+ ZFAEL(11)*ZSIN5 + ZFAEL(12)*ZCOS6 + ZFAEL(13)*ZSIN6 &
   &+ ZFAEL(14)*ZCOS7 + ZFAEL(15)*ZSIN7 + ZFAEL(16)*ZCOS8 &
   &+ ZFAEL(17)*ZSIN8 + ZFAEL(18)*ZCOS9 + ZFAEL(19)*ZSIN9 &
   &+ ZFAEL(20)*ZCOS10+ ZFAEL(21)*ZSIN10                 )
  ZAEU(IL) = ZFAEU(1) + _TWO_*&
   &( ZFAEU(2)*ZCOS1  + ZFAEU(3)*ZSIN1  + ZFAEU(4)*ZCOS2 &
   &+ ZFAEU(5)*ZSIN2  + ZFAEU(6)*ZCOS3  + ZFAEU(7)*ZSIN3 &
   &+ ZFAEU(8)*ZCOS4  + ZFAEU(9)*ZSIN4  + ZFAEU(10)*ZCOS5 &
   &+ ZFAEU(11)*ZSIN5 + ZFAEU(12)*ZCOS6 + ZFAEU(13)*ZSIN6 &
   &+ ZFAEU(14)*ZCOS7 + ZFAEU(15)*ZSIN7 + ZFAEU(16)*ZCOS8 &
   &+ ZFAEU(17)*ZSIN8 + ZFAEU(18)*ZCOS9 + ZFAEU(19)*ZSIN9 &
   &+ ZFAEU(20)*ZCOS10+ ZFAEU(21)*ZSIN10                 )
  ZAED(IL) = ZFAED(1) + _TWO_*&
   &( ZFAED(2)*ZCOS1  + ZFAED(3)*ZSIN1  + ZFAED(4)*ZCOS2 &
   &+ ZFAED(5)*ZSIN2  + ZFAED(6)*ZCOS3  + ZFAED(7)*ZSIN3 &
   &+ ZFAED(8)*ZCOS4  + ZFAED(9)*ZSIN4  + ZFAED(10)*ZCOS5 &
   &+ ZFAED(11)*ZSIN5 + ZFAED(12)*ZCOS6 + ZFAED(13)*ZSIN6 &
   &+ ZFAED(14)*ZCOS7 + ZFAED(15)*ZSIN7 + ZFAED(16)*ZCOS8 &
   &+ ZFAED(17)*ZSIN8 + ZFAED(18)*ZCOS9 + ZFAED(19)*ZSIN9 &
   &+ ZFAED(20)*ZCOS10+ ZFAED(21)*ZSIN10                 )
ENDDO


!     ------------------------------------------------------------------

!*       4.      VERTICAL DISTRIBUTION
!*               ---------------------


IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  ZDPO(IL)=PAPRS (JL,1)
  ZCPHO3=PAPRS (JL,1)**3
  ZSDPO3=SQRT  (ZCPHO3)
  IF (LNEWAER) THEN
    ZAEQSO(IL)= ZAERSS(IL)*CVDAES(1)
    ZAEQLO(IL)=(ZAEROR(IL)+ZAERSU(IL))*CVDAEL(1)
    ZAEQUO(IL)= ZAERBC(IL)*CVDAEU(1)
    ZAEQDO(IL)= ZAERSD(IL)*CVDAED(1)
  ELSE
    ZAEQSO(IL)=RCAEOPS*ZAES(IL)*CVDAES(1)
    ZAEQLO(IL)=RCAEOPL*ZAEL(IL)*CVDAEL(1)
    ZAEQUO(IL)=RCAEOPU*ZAEU(IL)*CVDAEU(1)
    ZAEQDO(IL)=RCAEOPD*ZAED(IL)*CVDAED(1)
  END IF  
  ZAETRO(IL)=_ONE_
  ZQOFO(IL)=ZOZQ(IL)*ZSDPO3 / (ZSDPO3 + ZOZH(IL))
ENDDO

DO JK=1,KLEV
  IL=KSHIFT
  IF (KCF == 0) THEN
    DO JL=KIDIA,KFDIA,KRINT
      IL=IL+1
      ZGRTH(IL)= PTH(JL,JK)/PTH(JL,JK+1)
    ENDDO
  ELSEIF (KCF == 1) THEN
    DO JL=KIDIA,KFDIA,KRINT
      IL=IL+1
      ZGRTH(IL)= PTH(IL,JK)/PTH(IL,JK+1)
    ENDDO
  ENDIF


  IL=KSHIFT
  DO JL=KIDIA,KFDIA,KRINT
    IL=IL+1
    ZDPN(IL)=PAPRS (JL,JK+1)
    ZCPHN3=PAPRS (JL,JK+1)**3
    ZSDPN3=SQRT  (ZCPHN3)
    IF (LNEWAER) THEN
      ZAEQSN(IL)= ZAERSS(IL)*CVDAES(JK+1)
      ZAEQLN(IL)=(ZAEROR(IL)+ZAERSU(IL))*CVDAEL(JK+1)
      ZAEQUN(IL)= ZAERBC(IL)*CVDAEU(JK+1)
      ZAEQDN(IL)= ZAERSD(IL)*CVDAED(JK+1)
    ELSE
      ZAEQSN(IL)=RCAEOPS*ZAES(IL)*CVDAES(JK+1)
      ZAEQLN(IL)=RCAEOPL*ZAEL(IL)*CVDAEL(JK+1)
      ZAEQUN(IL)=RCAEOPU*ZAEU(IL)*CVDAEU(JK+1)
      ZAEQDN(IL)=RCAEOPD*ZAED(IL)*CVDAED(JK+1)
    END IF  

    IF (_HALF_*(PAPRS(JL,JK)+PAPRS(JL,JK+1)) < 999._JPRB) THEN
! for models with top above 10hPa
      ZAETRN(IL)=_ONE_
      ZAETRO(IL)=_ONE_
    ELSE
      ZAETRN(IL)=ZAETRO(IL)*(MIN(_ONE_,     ZGRTH(IL)         ))**RCTRPT
    ENDIF

    ZAETR=SQRT  (ZAETRN(IL)*ZAETRO(IL))
    ZQOFN(IL)=ZOZQ(IL)*ZSDPN3/(ZSDPN3+ZOZH(IL))
    ZDPNMO    =ZDPN(IL)-ZDPO(IL)
    
    PAER(IL,1,JK)=(_ONE_-ZAETR)*(RCTRBGA*ZDPNMO+ ZAEQLN(IL)-ZAEQLO(IL))
    PAER(IL,2,JK)=(_ONE_-ZAETR)*(ZAEQSN(IL)-ZAEQSO(IL))
    PAER(IL,3,JK)=(_ONE_-ZAETR)*(ZAEQDN(IL)-ZAEQDO(IL))
    PAER(IL,4,JK)=(_ONE_-ZAETR)*(ZAEQUN(IL)-ZAEQUO(IL))
!old volc  PAER(IL,5,JK)= ZAETR * RCVOBGA*ZDPNMO
    PAER(IL,5,JK)=   ZAETR  * ZAERVO(IL) * ZDPNMO
    PAER(IL,6,JK)=   ZAETR  * RCSTBGA*ZDPNMO
    
!    print 9011,ZAETR,ZDPNMO,RCSTBGA,ZAERVO(IL),(PAER(IL,JAER,JK),JAER=1,6)
9011 format(1x,'ZAETR= ',10(1x,E12.5))

!old RH dependence          
!         AADS(IL,JK)=MAX(RCAEADM, (RCAEADK(1)*PAER(IL,1,JK)
!           + RCAEADK(2)*PAER(IL,2,JK)+RCAEADK(3)*PAER(IL,3,JK))/ZDPNMO)
    POZON(IL,JK)=ZQOFN(IL)-ZQOFO(IL)
!**** **************************************************
!**** **************************************************
  ENDDO
  IL=KSHIFT
  DO JL=KIDIA,KFDIA,KRINT
    IL=IL+1
    ZDPO(IL)=ZDPN(IL)
    ZQOFO(IL)=ZQOFN(IL)

    ZAEQSO(IL)=ZAEQSN(IL)
    ZAEQLO(IL)=ZAEQLN(IL)
    ZAEQUO(IL)=ZAEQUN(IL)
    ZAEQDO(IL)=ZAEQDN(IL)
    ZAETRO(IL)=ZAETRN(IL)
  ENDDO
  
  
!-- diagnostics in case of problem  
  DO JAER=1,6
    IL=KSHIFT
    DO JL=KIDIA,KFDIA,KRINT
      IL=IL+1
      PAER(IL,JAER,JK)=MAX(PAER(IL,JAER,JK),REPAER)
    END DO
    itot=il
  END DO   
!--
   
ENDDO


!IF (LNEWAER) THEN
!  print 9100,kshift+1,itot,ZSIN,REPAER
!9100   format(1x,'RADACA kshift/itot ',2I4,F10.6,E12.4)      
!  do jk=1,klev
!    do jl=kshift+1,itot,3
!      jend=min(jl+2,itot)
!      iprint=0
!      do jjl=jl,jend
!        do JAER=1,6
!          if (paer(jjl,jk,JAER).gt.1000.) then
!            iprint=1
!          end if
!        end do
!      end do    
!      if (iprint.eq.1) then
!        print 9102,JK,JL,((PAER(JIL,JAER,JK),JAER=1,6),JIL=JL,JEND)
!9102    format(1x,'RADACA ',I2,I4,18E9.2)
!      end if
!    end do
!  end do
!END IF  
!IF (IPRINT.EQ.1) THEN
!!  CALL ABOR1(' Problem with aerosols in RADACA')
!  STOP ' Problem with aerosols in RADACA'
!END IF  

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

RETURN
END SUBROUTINE RADACA
