SUBROUTINE SWR &
 &( KIDIA , KFDIA , KLON , KLEV  , KNU &
 &, PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU &
 &, PCGAZ , PPIZAZ, PRAY1, PRAY2 , PREFZ, PRJ  , PRK , PRMUE &
 &, PTAUAZ, PTRA1 , PTRA2, PTRCLD &
 &)

!**** *SWR* - CONTINUUM SCATTERING COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CONTINUUM SCATTERING

!**   INTERFACE.
!     ----------

!          *SWR* IS CALLED EITHER FROM *SW1S*
!                              OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
!     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)

!     EXTERNALS.
!     ----------

!          *SWDE*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        Ph. DANDIN Meteo-France 05-96 : Effect of cloud layer
!        JJMorcrette 990128 : sunshine duration
!        JJMorcrette 001218 : 6 spectral intervals
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE OYOERAD   , ONLY : NOVLP    ,NSW
USE YOECLD   , ONLY : REPSEC
USE YOEOVLP  , ONLY : RA1OVLP
!
USE MODI_SWDE
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KNU



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PALBD(KLON,NSW)      , PCG(KLON,NSW,KLEV)&
  &,  PCLD(KLON,KLEV)&
  &,  POMEGA(KLON,NSW,KLEV)&
  &,  PSEC(KLON)           , PTAU(KLON,NSW,KLEV)

REAL_B :: PRAY1(KLON,KLEV+1)   , PRAY2(KLON,KLEV+1)&
  &,  PREFZ(KLON,2,KLEV+1) , PRJ(KLON,6,KLEV+1)&
  &,  PRK(KLON,6,KLEV+1)   , PRMUE(KLON,KLEV+1)&
  &,  PCGAZ(KLON,KLEV)     , PPIZAZ(KLON,KLEV)&
  &,  PTAUAZ(KLON,KLEV)&
  &,  PTRA1(KLON,KLEV+1)   , PTRA2(KLON,KLEV+1)&
  &,  PTRCLD(KLON)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZC1I(KLON,KLEV+1)    , ZCLEQ(KLON,KLEV)&
  &,  ZCLEAR(KLON)         , ZCLOUD(KLON) &
  &,  ZGG(KLON)            , ZREF(KLON)&
  &,  ZRE1(KLON)           , ZRE2(KLON)&
  &,  ZRMUZ(KLON)          , ZRNEB(KLON)&
  &,  ZR21(KLON)           , ZR22(KLON)&
  &,  ZR23(KLON)           , ZSS1(KLON)&
  &,  ZTO1(KLON)           , ZTR(KLON,2,KLEV+1)&
  &,  ZTR1(KLON)           , ZTR2(KLON)&
  &,  ZW(KLON)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IKL, IKLP1, JA, JAJ, JK, JKM1, JL, INU1

!     LOCAL REAL SCALARS
REAL_B :: ZBMU0, ZBMU1, ZCORAE, ZCORCD, ZDEN, ZDEN1,&
          &ZFACOA, ZFACOC, ZGAP, ZMU1, ZMUE, ZRE11, &
          &ZTO, ZWW, ZALPHA1




!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------


DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = _ZERO_
      PRK(JL,JA,JK) = _ZERO_
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------


DO JL = KIDIA,KFDIA
  ZR23(JL) = _ZERO_
  ZC1I(JL,KLEV+1) = _ZERO_
  ZCLEAR(JL) = _ONE_
  ZCLOUD(JL) = _ZERO_
ENDDO

JK = 1
IKL = KLEV+1 - JK
IKLP1 = IKL + 1
DO JL = KIDIA,KFDIA
!++MODIF_MESONH
  IF (NOVLP.GE.5) THEN !MESONH VERSION
   ZFACOA =PTAUAZ(JL,IKL) 
   ZFACOC = _ONE_ - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
   ZCORAE = ZFACOA * PSEC(JL)
   ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
  ELSE !ECMWF VERSION
   ZFACOA = _ONE_ - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)
   ZFACOC = _ONE_ - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
   ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
   ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
  ENDIF
!--MODIF_MESONH
  ZR21(JL) = EXP(-ZCORAE   )
  ZR22(JL) = EXP(-ZCORCD   )
  ZSS1(JL) = PCLD(JL,IKL)*(_ONE_-ZR21(JL)*ZR22(JL))&
   &+ (_ONE_-PCLD(JL,IKL))*(_ONE_-ZR21(JL))
  ZCLEQ(JL,IKL) = ZSS1(JL)

!++MODIF_MESONH
  IF ((NOVLP == 1).OR.(NOVLP == 8)) THEN
!--MODIF_MESONH
!* maximum-random      
    ZCLEAR(JL) = ZCLEAR(JL)&
     &*(_ONE_-MAX(ZSS1(JL),ZCLOUD(JL)))&
     &/(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPSEC))
    ZC1I(JL,IKL) = _ONE_ - ZCLEAR(JL)
    ZCLOUD(JL) = ZSS1(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
    ZC1I(JL,IKL) = ZCLOUD(JL)
!++MODIF_MESONH
  ELSEIF ((NOVLP == 3).OR.((NOVLP .GE. 5).AND.(NOVLP.NE.8))) THEN
!--MODIF_MESONH
!* random
    ZCLEAR(JL) = ZCLEAR(JL)*(_ONE_ - ZSS1(JL))
    ZCLOUD(JL) = _ONE_ - ZCLEAR(JL)
    ZC1I(JL,IKL) = ZCLOUD(JL)
  ELSEIF (NOVLP == 4) THEN
!* Hogan & Illingworth, 2001  
    ZALPHA1=RA1OVLP(KLEV+1-JK)
    ZCLEAR(JL)=ZCLEAR(JL)*( &
      & ZALPHA1*(_ONE_-MAX(ZSS1(JL),ZCLOUD(JL))) &
      &        /(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPSEC)) &
      & +(_ONE_-ZALPHA1)*(_ONE_-ZSS1(JL)) )
    ZC1I(JL,IKL) = _ONE_ - ZCLEAR(JL) 
    ZCLOUD(JL) = ZSS1(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  IKL = KLEV+1 - JK
  IKLP1 = IKL + 1
  DO JL = KIDIA,KFDIA
!++MODIF_MESONH
    IF (NOVLP.GE.5) THEN !MESONH VERSION
     ZFACOA =PTAUAZ(JL,IKL) 
     ZFACOC = _ONE_ - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
     ZCORAE = ZFACOA * PSEC(JL)
     ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
    ELSE !ECMWF VERSION
     ZFACOA = _ONE_ - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)
     ZFACOC = _ONE_ - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
     ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
     ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
    ENDIF
!--MODIF_MESONH
    ZR21(JL) = EXP(-ZCORAE   )
    ZR22(JL) = EXP(-ZCORCD   )
    ZSS1(JL) = PCLD(JL,IKL)*(_ONE_-ZR21(JL)*ZR22(JL))&
     &+ (_ONE_-PCLD(JL,IKL))*(_ONE_-ZR21(JL))
    ZCLEQ(JL,IKL) = ZSS1(JL)

!++MODIF_MESONH
    IF ((NOVLP == 1).OR.(NOVLP == 8)) THEN
!--MODIF_MESONH
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       &*(_ONE_-MAX(ZSS1(JL),ZCLOUD(JL)))&
       &/(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPSEC))
      ZC1I(JL,IKL) = _ONE_ - ZCLEAR(JL)
      ZCLOUD(JL) = ZSS1(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
      ZC1I(JL,IKL) = ZCLOUD(JL)
!++MODIF_MESONH
  ELSEIF ((NOVLP == 3).OR.((NOVLP .GE. 5).AND.(NOVLP.NE.8))) THEN
!--MODIF_MESONH
!* random
      ZCLEAR(JL) = ZCLEAR(JL)*(_ONE_ - ZSS1(JL))
      ZCLOUD(JL) = _ONE_ - ZCLEAR(JL)
      ZC1I(JL,IKL) = ZCLOUD(JL)
    ELSEIF (NOVLP == 4) THEN
!* Hogan & Illingworth, 2001  
      ZALPHA1=RA1OVLP(KLEV+1-JK)
      ZCLEAR(JL)=ZCLEAR(JL)*( &
        & ZALPHA1*(_ONE_-MAX(ZSS1(JL),ZCLOUD(JL))) &
        &        /(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPSEC)) &
        & +(_ONE_-ZALPHA1)*(_ONE_-ZSS1(JL)) )
      ZC1I(JL,IKL) = _ONE_ - ZCLEAR(JL) 
      ZCLOUD(JL) = ZSS1(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------


DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = _ZERO_
  PRAY2(JL,KLEV+1) = _ZERO_
  PREFZ(JL,2,1) = PALBD(JL,KNU)
  PREFZ(JL,1,1) = PALBD(JL,KNU)
  PTRA1(JL,KLEV+1) = _ONE_
  PTRA2(JL,KLEV+1) = _ONE_
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA
    ZRNEB(JL)= PCLD(JL,JKM1)
    ZRE1(JL)=_ZERO_
    ZTR1(JL)=_ZERO_
    ZRE2(JL)=_ZERO_
    ZTR2(JL)=_ZERO_


!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------


    ZMUE = (_ONE_-ZC1I(JL,JK)) * PSEC(JL)+ ZC1I(JL,JK) * 1.66_JPRB
!-- just to test Box-type computations
!    ZMUE = PSEC(JL)
    PRMUE(JL,JK) = _ONE_/ZMUE


!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------


    ZGAP = PCGAZ(JL,JKM1)
    ZBMU0 = _HALF_ - 0.75_JPRB * ZGAP / ZMUE
    ZWW = PPIZAZ(JL,JKM1)
    ZTO = PTAUAZ(JL,JKM1)
    ZDEN = _ONE_ + (_ONE_ - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE &
     &+ (1-ZWW) * (_ONE_ - ZWW +_TWO_*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
    PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
    PTRA1(JL,JKM1) = _ONE_ / ZDEN

    ZMU1 = _HALF_
    ZBMU1 = _HALF_ - 0.75_JPRB * ZGAP * ZMU1
    ZDEN1= _ONE_ + (_ONE_ - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1 &
     &+ (1-ZWW) * (_ONE_ - ZWW +_TWO_*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
    PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
    PTRA2(JL,JKM1) = _ONE_ / ZDEN1


!     ------------------------------------------------------------------

!*         3.3  EFFECT OF CLOUD LAYER
!               ---------------------


!++MODIF_MESONH
    IF (NOVLP.GE.5)THEN !MESONH VERSION
     ZW(JL) =PCG(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)
     ZTO1(JL) = PTAU(JL,KNU,JKM1)*(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
     ZW(JL) =POMEGA(JL,KNU,JKM1)*(1-ZW(JL))/(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
     ZGG(JL) = PCG(JL,KNU,JKM1)/(1+PCG(JL,KNU,JKM1))
     ZGG(JL)=ZTO1(JL)*ZW(JL)*ZGG(JL)+PTAUAZ(JL,JKM1)*PPIZAZ(JL,JKM1)*PCGAZ(JL,JKM1)
     ZW(JL) =ZTO1(JL)*ZW(JL)+PTAUAZ(JL,JKM1)*PPIZAZ(JL,JKM1)
     ZTO1(JL) = ZTO1(JL) +  PTAUAZ(JL,JKM1)
     ZGG(JL)=ZGG(JL)/ZW(JL)
     ZW(JL) =ZW(JL)/ZTO1(JL)
    ELSE !ECMWF VERSION
     ZW(JL) = POMEGA(JL,KNU,JKM1)
     ZTO1(JL) = PTAU(JL,KNU,JKM1)/ZW(JL)+ PTAUAZ(JL,JKM1)/PPIZAZ(JL,JKM1)
     ZR21(JL) = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
     ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
     ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
      &+ (_ONE_ - ZR22(JL)) * PCGAZ(JL,JKM1)
     IF (ZW(JL) == _ONE_ .AND. PPIZAZ(JL,JKM1) == _ONE_) THEN
      ZW(JL)=_ONE_
     ELSE
      ZW(JL) = ZR21(JL) / ZTO1(JL)
     ENDIF
    ENDIF
!--MODIF_MESONH
     
    ZREF(JL) = PREFZ(JL,1,JKM1)
    ZRMUZ(JL) = PRMUE(JL,JK)
  ENDDO

  CALL SWDE ( KIDIA, KFDIA , KLON &
   &, ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW &
   &, ZRE1 , ZRE2  , ZTR1  , ZTR2      )

  DO JL = KIDIA,KFDIA

    PREFZ(JL,1,JK) = (_ONE_-ZRNEB(JL)) * (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1)&
     &/ (_ONE_-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))&
     &+ ZRNEB(JL) * ZRE2(JL)

    ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKM1)&
     &/ (_ONE_-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))&
     &* (_ONE_-ZRNEB(JL))

    PREFZ(JL,2,JK) = (_ONE_-ZRNEB(JL)) * (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1) )&
     &+ ZRNEB(JL) * ZRE1(JL)

    ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)+ PTRA1(JL,JKM1) * (_ONE_-ZRNEB(JL))

  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (_ONE_-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66_JPRB
!-- just to test Box-type computations
!  ZMUE = PSEC(JL)
  PRMUE(JL,1)=_ONE_/ZMUE
  PTRCLD(JL)=_ONE_-ZC1I(JL,1)
ENDDO


!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------


IF (NSW <= 4) THEN
  INU1=1
ELSE IF (NSW == 6) THEN
  INU1=3
END IF    

IF (KNU <= INU1) THEN
  JAJ = 2
  DO JL = KIDIA,KFDIA
    PRJ(JL,JAJ,KLEV+1) = _ONE_
    PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
  ENDDO

  DO JK = 1 , KLEV
    IKL = KLEV+1 - JK
    IKLP1 = IKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,  1,IKL)
      PRJ(JL,JAJ,IKL) = ZRE11
      PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,  1,IKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = _ONE_
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      IKL = KLEV+1 - JK
      IKLP1 = IKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,JAJ,IKL)
        PRJ(JL,JAJ,IKL) = ZRE11
        PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,JAJ,IKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWR

