SUBROUTINE LWC &
 &( KIDIA , KFDIA, KLON  , KLEV &
 &, PBINT , PBSUI, PCLDLD, PCLDLU &
 &, PCNTRB, PEMIT, PFLUC &
 &, PFLUX                              &
 &)

!**** *LWC*   - LONGWAVE RADIATION, CLOUD EFFECTS

!     PURPOSE.
!     --------
!           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
!           RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PBINT  : (KLON,KLEV+1)       ; HALF LEVEL PLANCK FUNCTION
! PBSUI  : (KLON)              ; SURFACE PLANCK FUNCTION
! PCLDLD : (KLON,KLEV)         ; DOWNWARD EFFECTIVE CLOUD FRACTION
! PCLDLU : (KLON,KLEV)         ; UPWARD EFFECTIVE CLOUD FRACTION
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PEMIT  : (KLON)              ; SURFACE TOTAL LW EMISSIVITY
! PFLUC  : (KLON,2,KLEV+1)     ; CLEAR-SKY LW RADIATIVE FLUXES
!     ==== OUTPUTS ===
! PFLUX  : (KLON,2,KLEV+1)     ; TOTAL SKY LW RADIATIVE FLUXES :
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
!          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
!          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
!     CLOUDS

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        JJ Morcrette 97-04-18   Cleaning

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOERAD   , ONLY : NOVLP
USE OYOERDI   , ONLY : REPCLC
USE YOEOVLP  , ONLY : RA1OVLP


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B ::&
     &PBINT(KLON,KLEV+1), PBSUI(KLON)    &
  &,  PCLDLD(KLON,KLEV) , PCLDLU(KLON,KLEV)&
  &,  PCNTRB(KLON,KLEV+1,KLEV+1)&
  &,  PEMIT(KLON)&
  &,  PFLUC(KLON,2,KLEV+1) 

REAL_B :: PFLUX(KLON,2,KLEV+1)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZCLEAR(KLON)            , ZCLOUD(KLON)&
  &,  ZCLM(KLON,KLEV+1,KLEV+1), ZDNF(KLON,KLEV+1,KLEV+1)&
  &,  ZFD(KLON)               , ZFU(KLON)&
  &,  ZUPF(KLON,KLEV+1,KLEV+1)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IKCP1, IKM1, IKP1, IMAXC, IMXM1, IMXP1, JCLOUD,&
             &JK, JK1, JK2, JKJ, JL

!     LOCAL REAL SCALARS
REAL_B :: ZALPHA1, ZCFRAC


!     ------------------------------------------------------------------

!*         1.     INITIALIZATION
!                 --------------

!100  CONTINUE

!print *,' Enter LWC '
DO JL = KIDIA,KFDIA
  ZCLOUD(JL) = _ZERO_
ENDDO

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFLUX(JL,1,JK) = PFLUC(JL,1,JK)
    PFLUX(JL,2,JK) = PFLUC(JL,2,JK)
  ENDDO
ENDDO

!GM*******
IMAXC=KLEV
!GM*******


!     ------------------------------------------------------------------

!*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
!                  ---------------------------------------


IMXP1 = IMAXC + 1
IMXM1 = IMAXC - 1

!*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
!                  ------------------------------

!200  CONTINUE

DO JK1=1,KLEV+1
  DO JK2=1,KLEV+1
    DO JL = KIDIA,KFDIA
      ZUPF(JL,JK2,JK1)=PFLUC(JL,1,JK1)
      ZDNF(JL,JK2,JK1)=PFLUC(JL,2,JK1)
    ENDDO
  ENDDO
ENDDO
!print *,' LWC after Initialisation to clear-sky fluxes'

!*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
!                  ----------------------------------------------

!210  CONTINUE

DO JCLOUD = 1 , IMAXC
  IKCP1=JCLOUD+1

!*         2.1.1   ABOVE THE CLOUD
!                  ---------------

!2110 CONTINUE

  DO JK=IKCP1,KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZFU(JL)=_ZERO_
    ENDDO

    IF (JK  >  IKCP1) THEN
      DO JKJ=IKCP1,IKM1
        DO JL = KIDIA,KFDIA
          ZFU(JL) = ZFU(JL) + PCNTRB(JL,JK,JKJ)
        ENDDO
      ENDDO
    ENDIF

    DO JL = KIDIA,KFDIA
      ZUPF(JL,IKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
    ENDDO
  ENDDO

!*         2.1.2   BELOW THE CLOUD
!                  ---------------

!2120 CONTINUE

  DO JK=1,JCLOUD
    IKP1=JK+1
    DO JL = KIDIA,KFDIA
      ZFD(JL)=_ZERO_
    ENDDO

    IF (JK  <  JCLOUD) THEN
      DO JKJ=IKP1,JCLOUD
        DO JL = KIDIA,KFDIA
          ZFD(JL) = ZFD(JL) + PCNTRB(JL,JK,JKJ)
        ENDDO
      ENDDO
    ENDIF

    DO JL = KIDIA,KFDIA
      ZDNF(JL,IKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
    ENDDO
  ENDDO

ENDDO
!print *,' LWC after 213: Fluxes for unity emissivity'


!*         2.2     CLOUD COVER MATRIX
!                  ------------------

!*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
!     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1

!220  CONTINUE

DO JK1 = 1 , KLEV+1
  DO JK2 = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZCLM(JL,JK1,JK2) = _ZERO_
    ENDDO
  ENDDO
ENDDO
!print *,' LWC after Initialisation CC matrix'



!*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
!                  ------------------------------------------

!240  CONTINUE

DO JK1 = 2 , KLEV+1
  DO JL = KIDIA,KFDIA
    ZCLEAR(JL)=_ONE_
    ZCLOUD(JL)=_ZERO_
  ENDDO

  DO JK = JK1 - 1 , 1 , -1
    DO JL = KIDIA,KFDIA
!++MODIF_MESONH
      IF ((NOVLP==1).OR.(NOVLP==6).OR.(NOVLP==8)) THEN
!--MODIF_MESONH
!* maximum-random       
        ZCLEAR(JL)=ZCLEAR(JL)*(_ONE_-MAX(PCLDLU(JL,JK),ZCLOUD(JL)))&
         &/(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPCLC))
        ZCLM(JL,JK1,JK) = _ONE_ - ZCLEAR(JL)
        ZCLOUD(JL) = PCLDLU(JL,JK)
!++MODIF_MESONH
      ELSEIF ((NOVLP==2).OR.(NOVLP==7)) THEN
!--MODIF_MESONH
!* maximum      
        ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLU(JL,JK))
        ZCLM(JL,JK1,JK) = ZCLOUD(JL)
!++MODIF_MESONH
      ELSEIF ((NOVLP == 3).OR.(NOVLP==5)) THEN
!--MODIF_MESONH
!* random      
        ZCLEAR(JL) = ZCLEAR(JL)*(_ONE_ - PCLDLU(JL,JK))
        ZCLOUD(JL) = _ONE_ - ZCLEAR(JL)
        ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSEIF (NOVLP == 4) THEN
!** Hogan & Illingworth (2001)      
        ZALPHA1=RA1OVLP(KLEV+1-JK)
        ZCLEAR(JL)=ZCLEAR(JL)*( &
         & ZALPHA1*(_ONE_-MAX(PCLDLU(JL,JK),ZCLOUD(JL))) &
         &        /(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPCLC)) &
         & +(_ONE_-ZALPHA1)*(_ONE_-PCLDLU(JL,JK)) )
        ZCLM(JL,JK1,JK) = _ONE_ - ZCLEAR(JL) 
        ZCLOUD(JL) = PCLDLU(JL,JK)
      ENDIF
    ENDDO
  ENDDO

ENDDO
!print *,' LWC after 244: CC below level of calculation'


!*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
!                  ------------------------------------------

!250  CONTINUE

DO JK1 = 1 , KLEV
  DO JL = KIDIA,KFDIA
    ZCLEAR(JL)=_ONE_
    ZCLOUD(JL)=_ZERO_
  ENDDO

  DO JK = JK1 , KLEV
    DO JL = KIDIA,KFDIA
!++MODIF_MESONH
      IF ((NOVLP == 1).OR.(NOVLP==6).OR.(NOVLP==8)) THEN
!--MODIF_MESONH
!* maximum-random       
        ZCLEAR(JL)=ZCLEAR(JL)*(_ONE_-MAX(PCLDLD(JL,JK),ZCLOUD(JL)))&
         &/(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPCLC))
        ZCLM(JL,JK1,JK) = _ONE_ - ZCLEAR(JL)
        ZCLOUD(JL) = PCLDLD(JL,JK)
!++MODIF_MESONH
      ELSEIF ((NOVLP == 2).OR.(NOVLP==7)) THEN
!--MODIF_MESONH
!* maximum      
        ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLD(JL,JK))
        ZCLM(JL,JK1,JK) = ZCLOUD(JL)
!++MODIF_MESONH
      ELSEIF ((NOVLP == 3).OR.(NOVLP==5)) THEN
!--MODIF_MESONH
!* random      
        ZCLEAR(JL) = ZCLEAR(JL)*(_ONE_ - PCLDLD(JL,JK))
        ZCLOUD(JL) = _ONE_ - ZCLEAR(JL)
        ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSEIF (NOVLP == 4) THEN
!** Hogan & Illingworth (2001)      
        ZALPHA1=RA1OVLP(KLEV+1-JK)
        ZCLEAR(JL)=ZCLEAR(JL)*( &
         & ZALPHA1*(_ONE_-MAX(PCLDLD(JL,JK),ZCLOUD(JL))) &
         &        /(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPCLC)) &
         & +(_ONE_-ZALPHA1)*(_ONE_ - PCLDLD(JL,JK)) )
        ZCLM(JL,JK1,JK) = _ONE_ - ZCLEAR(JL)
        ZCLOUD(JL) = PCLDLD(JL,JK)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!print *,' LWC after 254: CC above level of calculation'



!*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
!                  ----------------------------------------------

!300  CONTINUE

!*         3.1     DOWNWARD FLUXES
!                  ---------------

!310  CONTINUE

DO JL = KIDIA,KFDIA
  PFLUX(JL,2,KLEV+1) = _ZERO_
ENDDO

DO JK1 = KLEV , 1 , -1

!*                 CONTRIBUTION FROM CLEAR-SKY FRACTION

  DO JL = KIDIA,KFDIA
    ZFD (JL) = (_ONE_ - ZCLM(JL,JK1,KLEV)) * ZDNF(JL,1,JK1)

!*                 CONTRIBUTION FROM ADJACENT CLOUD

    ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
  ENDDO

!*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

  DO JK = KLEV-1 , JK1 , -1
    DO JL = KIDIA,KFDIA
      ZCFRAC = ZCLM(JL,JK1,JK+1) - ZCLM(JL,JK1,JK)
      ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
    ENDDO
  ENDDO

  DO JL = KIDIA,KFDIA
    PFLUX(JL,2,JK1) = ZFD (JL)
  ENDDO

ENDDO
!      print *,' LWC after 317: Downward fluxes'




!*         3.2     UPWARD FLUX AT THE SURFACE
!                  --------------------------

!320  CONTINUE

DO JL = KIDIA,KFDIA
  PFLUX(JL,1,1) = PEMIT(JL)*PBSUI(JL)-(_ONE_-PEMIT(JL))*PFLUX(JL,2,1)
ENDDO



!*         3.3     UPWARD FLUXES
!                  -------------

!330  CONTINUE

DO JK1 = 2 , KLEV+1

!*                 CONTRIBUTION FROM CLEAR-SKY FRACTION

  DO JL = KIDIA,KFDIA
    ZFU (JL) = (_ONE_ - ZCLM(JL,JK1,1)) * ZUPF(JL,1,JK1)

!*                 CONTRIBUTION FROM ADJACENT CLOUD

    ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)
  ENDDO

!*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS

  DO JK = 2 , JK1-1
    DO JL = KIDIA,KFDIA
      ZCFRAC = ZCLM(JL,JK1,JK-1) - ZCLM(JL,JK1,JK)
      ZFU(JL) =  ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
    ENDDO
  ENDDO

  DO JL = KIDIA,KFDIA
    PFLUX(JL,1,JK1) = ZFU (JL)
  ENDDO

ENDDO
!      print *,' LWC after 337: Upward fluxes'

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE LWC

