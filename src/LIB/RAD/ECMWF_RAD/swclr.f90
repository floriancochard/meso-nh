SUBROUTINE SWCLR &
  &( KIDIA , KFDIA , KLON  , KLEV  , KAER  , KNU &
  &, PAER  , PALBP , PDSIG , PRAYL , PSEC &
  &, PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  &
  &, PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2 , PTRCLR &
!++MODIF_MESONH
  &, ODUST,PPIZA_DST, PCGA_DST, PTAUREL_DST )
!--MODIF_MESONH

!**** *SWCLR* - CLEAR-SKY COLUMN COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CLEAR-SKY COLUMN

!**   INTERFACE.
!     ----------

!          *SWCLR* IS CALLED EITHER FROM *SW1S*
!                                OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-11-15
!        Modified : 96-03-19 JJM-PhD (loop 107 in absence of aerosols)
!        JJMorcrette 990128 : sunshine duration
!        JJMorcrette 990128 : sunshine duration
!        99-05-25   JJMorcrette    Revised aerosols
!        JJMorcrette 001218 : 6 spectral intervals
!        D.St Martin 11/2015 : bug on ZFACOA for NOVLP>= 5
   
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE OYOESW    , ONLY : RTAUA    ,RPIZA    ,RCGA
USE OYOERAD   , ONLY : NOVLP    ,NSW
USE OYOERDI   , ONLY : REPCLC
USE YOERDU   , ONLY : REPSCT


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KNU



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAER(KLON,6,KLEV), PALBP(KLON,NSW)&
  &,  PDSIG(KLON,KLEV)&
  &,  PRAYL(KLON)&
  &,  PSEC(KLON)

REAL_B ::&
     &PCGAZ(KLON,KLEV)     &
  &,  PPIZAZ(KLON,KLEV)&
  &,  PRAY1(KLON,KLEV+1)  , PRAY2(KLON,KLEV+1)&
  &,  PREFZ(KLON,2,KLEV+1), PRJ(KLON,6,KLEV+1)&
  &,  PRK(KLON,6,KLEV+1)  , PRMU0(KLON,KLEV+1)&
  &,  PTAUAZ(KLON,KLEV)&
  &,  PTRA1(KLON,KLEV+1)  , PTRA2(KLON,KLEV+1)&
  &,  PTRCLR(KLON)

!
!++MODIF_MESONH
LOGICAL :: ODUST                   ! flag for DUST
REAL_B  :: PPIZA_DST(KLON,KLEV)    !ssa dust for current wavelength
REAL_B  :: PCGA_DST(KLON,KLEV)     !assym.fact dust for current wavelength
REAL_B  :: PTAUREL_DST(KLON,KLEV)  !tau/tau_{550} for current wavelength
!--MODIF_MESONH


!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZC0I(KLON,KLEV+1)&
  &,  ZCLE0(KLON,KLEV), ZCLEAR(KLON) &
  &,  ZR21(KLON)&
  &,  ZR23(KLON) , ZSS0(KLON) , ZSCAT(KLON)&
  &,  ZTR(KLON,2,KLEV+1)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IKL, JA, JAE, JAJ, JK, JKL, JKLP1, JKM1, JL, INU1

!     LOCAL REAL SCALARS
REAL_B :: ZBMU0, ZBMU1, ZCORAE, ZDEN, ZDEN1, ZFACOA,&
          &ZFF, ZGAP, ZGAR, ZMU1, ZMUE, ZRATIO, ZRE11, &
          &ZTO, ZTRAY, ZWW
!++MODIF_MESONH
REAL_B ::ZFACOA_NEW(KLON,KLEV)
!--MODIF_MESONH


!     ------------------------------------------------------------------

!*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
!                --------------------------------------------
ZFACOA_NEW(:,:)=0.0
DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = _ZERO_
      PRK(JL,JA,JK) = _ZERO_
    ENDDO
  ENDDO
ENDDO

! ------   NB: 'PAER' AEROSOLS ARE ENTERED FROM TOP TO BOTTOM
DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    PCGAZ(JL,JK) = _ZERO_
    PPIZAZ(JL,JK) =  _ZERO_
    PTAUAZ(JL,JK) = _ZERO_
  ENDDO

!++MODIF_MESONH  
  IF(NOVLP.LT.5)THEN !ECMWF VERSION
     DO JAE=1,6
        DO JL = KIDIA,KFDIA
           PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL, JAE, IKL)*RTAUA(KNU,JAE)
           PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL, JAE, IKL)&
                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)
           PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL, JAE, IKL)&
                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
        ENDDO
     ENDDO
  ELSE ! MESONH VERSION
     DO JAE=1,6
        DO JL = KIDIA,KFDIA
           !Special optical properties for dust
           IF (ODUST.AND.(JAE==3)) THEN
           !Ponderation of aerosol optical properties:first step 
           !ti
             PTAUAZ(JL,JK)=PTAUAZ(JL,JK) + PAER(JL,JAE,IKL) * PTAUREL_DST(JL,IKL)
           !wi*ti
             PPIZAZ(JL,JK)=PPIZAZ(JL,JK) + PAER(JL,JAE,IKL)  &
                    *PTAUREL_DST(JL,IKL)*PPIZA_DST(JL,IKL)
           !wi*ti*gi
             PCGAZ(JL,JK) = PCGAZ(JL,JK) + PAER(JL,JAE,IKL) &
                  *PTAUREL_DST(JL,IKL)*PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)
           !wi*ti*(gi**2)
             ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+PAER(JL, JAE, IKL)&
                &*PTAUREL_DST(JL,IKL) *PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)*PCGA_DST(JL,IKL)
           ELSE
           !Ponderation of aerosol optical properties:first step 
           !ti
             PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL, JAE, IKL)*RTAUA(KNU,JAE)
           !wi*ti
             PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL, JAE, IKL)&
                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)
           !wi*ti*gi
             PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL, JAE, IKL)&
                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
           !wi*ti*(gi**2)
             ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+PAER(JL, JAE, IKL)&
                &* RTAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)*RCGA(KNU,JAE)
           ENDIF
        ENDDO
     ENDDO
  ENDIF
!--MODIF_MESONH  

!++MODIF_MESONH  
  IF (NOVLP.LT.5) then !ECMWF VERSION
   DO JL = KIDIA,KFDIA
    IF (KAER /= 0) THEN
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
      ZGAR = PCGAZ(JL,JK)
      ZFF = ZGAR * ZGAR
      PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(_ONE_-PPIZAZ(JL,JK)*ZFF)
      PCGAZ(JL,JK) = ZGAR * (_ONE_ - ZRATIO) / (_ONE_ + ZGAR)
      PPIZAZ(JL,JK) =ZRATIO+(_ONE_-ZRATIO)*PPIZAZ(JL,JK)*(_ONE_-ZFF)&
       &/ (_ONE_ - PPIZAZ(JL,JK) * ZFF)
    ELSE
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
            PTAUAZ(JL,JK) = ZTRAY
      PCGAZ(JL,JK) = _ZERO_
      PPIZAZ(JL,JK) = _ONE_-REPSCT
    ENDIF
   ENDDO
  ELSE !MESONH VERSION
   DO JL = KIDIA,KFDIA
    IF (KAER /= 0) THEN
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      ZRATIO =PPIZAZ(JL,JK)+ZTRAY  ! wi*ti+1*tr
      !Ponderation G**2
      ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)/ZRATIO
      !Ponderation w
      PPIZAZ(JL,JK)=ZRATIO/(PTAUAZ(JL,JK)+ZTRAY)
      !Ponderation g
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/ZRATIO
      !Ponderation+delta-modified parameters tau - applies delta-Eddington to the complete phase function 
      PTAUAZ(JL,JK)=(ZTRAY+PTAUAZ(JL,JK))*&
          (_ONE_-PPIZAZ(JL,JK)*ZFACOA_NEW(JL,JK))
      !delta-modified parameters w
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)*(1-ZFACOA_NEW(JL,JK))/&
           (1-ZFACOA_NEW(JL,JK)*PPIZAZ(JL,JK))     
      !delta-modified parameters g
      PCGAZ(JL,JK)=PCGAZ(JL,JK)/(1+PCGAZ(JL,JK))
      
    ELSE
      ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
      ZFACOA_NEW(JL,JK)=_ZERO_
      PTAUAZ(JL,JK) = ZTRAY
      PCGAZ(JL,JK) = _ZERO_
      PPIZAZ(JL,JK) = _ONE_-REPSCT
    ENDIF
   ENDDO    
  ENDIF
!--MODIF_MESONH  


ENDDO
!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------


DO JL = KIDIA,KFDIA
  ZR23(JL) = _ZERO_
  ZC0I(JL,KLEV+1) = _ZERO_
  ZCLEAR(JL) = _ONE_
  ZSCAT(JL) = _ZERO_
ENDDO

JK = 1
JKL = KLEV+1 - JK
JKLP1 = JKL + 1
DO JL = KIDIA,KFDIA
!++MODIF_MESONH
  IF (NOVLP.GE.5) THEN
   ZFACOA = PTAUAZ(JL,JKL)
   ZCORAE = ZFACOA *  PSEC(JL)
  ELSE
   ZFACOA = _ONE_ - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
   ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
  ENDIF
!--MODIF_MESONH
  ZR21(JL) = EXP(-ZCORAE   )
  ZSS0(JL) = _ONE_-ZR21(JL)
  ZCLE0(JL,JKL) = ZSS0(JL)

!++MODIF_MESONH
  IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!--MODIF_MESONH
!* maximum-random      
    ZCLEAR(JL) = ZCLEAR(JL)&
     &*(_ONE_-MAX(ZSS0(JL),ZSCAT(JL)))&
     &/(_ONE_-MIN(ZSCAT(JL),_ONE_-REPCLC))
    ZC0I(JL,JKL) = _ONE_ - ZCLEAR(JL)
    ZSCAT(JL) = ZSS0(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
    ZC0I(JL,JKL) = ZSCAT(JL)
!++MODIF_MESONH
  ELSEIF ((NOVLP == 3).OR.(NOVLP .ge. 5)) THEN
!--MODIF_MESONH
!* random
    ZCLEAR(JL)=ZCLEAR(JL)*(_ONE_-ZSS0(JL))
    ZSCAT(JL) = _ONE_ - ZCLEAR(JL)
    ZC0I(JL,JKL) = ZSCAT(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  JKL = KLEV+1 - JK
  JKLP1 = JKL + 1
  DO JL = KIDIA,KFDIA
!++MODIF_MESONH
    IF (NOVLP.GE.5) THEN
     ZFACOA = PTAUAZ(JL,JKL)
     ZCORAE = ZFACOA *  PSEC(JL)
    ELSE
     ZFACOA = _ONE_ - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
     ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
    ENDIF
!--MODIF_MESONH
    ZR21(JL) = EXP(-ZCORAE   )
    ZSS0(JL) = _ONE_-ZR21(JL)
    ZCLE0(JL,JKL) = ZSS0(JL)

    IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       &*(_ONE_-MAX(ZSS0(JL),ZSCAT(JL)))&
       &/(_ONE_-MIN(ZSCAT(JL),_ONE_-REPCLC))
      ZC0I(JL,JKL) = _ONE_ - ZCLEAR(JL)
      ZSCAT(JL) = ZSS0(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
      ZC0I(JL,JKL) = ZSCAT(JL)
!++MODIF_MESONH
    ELSEIF ((NOVLP == 3).OR.(NOVLP.ge.5)) THEN
!--MODIF_MESONH
!* random
      ZCLEAR(JL)=ZCLEAR(JL)*(_ONE_-ZSS0(JL))
      ZSCAT(JL) = _ONE_ - ZCLEAR(JL)
      ZC0I(JL,JKL) = ZSCAT(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------


DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = _ZERO_
  PRAY2(JL,KLEV+1) = _ZERO_
  PREFZ(JL,2,1) = PALBP(JL,KNU)
  PREFZ(JL,1,1) = PALBP(JL,KNU)
  PTRA1(JL,KLEV+1) = _ONE_
  PTRA2(JL,KLEV+1) = _ONE_
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA


!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------


    ZMUE = (_ONE_-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_JPRB
    PRMU0(JL,JK) = _ONE_/ZMUE


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



    PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1)&
     &/ (_ONE_-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

    ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1)&
     &/ (_ONE_-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

    PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)&
     &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     &* PTRA2(JL,JKM1) )

    ZTR(JL,2,JKM1) = PTRA1(JL,JKM1)

  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (_ONE_-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66_JPRB
  PRMU0(JL,1)=_ONE_/ZMUE
  PTRCLR(JL)=_ONE_-ZC0I(JL,1)
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
    JKL = KLEV+1 - JK
    JKLP1 = JKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = _ONE_
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
        PRJ(JL,JAJ,JKL) = ZRE11
        PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------
RETURN
END SUBROUTINE SWCLR



