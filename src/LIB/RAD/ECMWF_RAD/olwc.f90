!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWC ( KIDIA,KFDIA,KLON,KLEV &
     &  , PBINT,PBSUIN,PCLDLD,PCLDLU,PCNTRB,PEMIS,PFDN,PFUP &
     &  , PFLUX                                                 )
!
!**** *LWC*   - LONGWAVE RADIATION, CLOUD EFFECTS
!
!     PURPOSE.
!     --------
!           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
!           RADIANCES
!
!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PBINT  : (KLON,KLEV+1)     ; HALF LEVEL PLANCK FUNCTION
! PBSUIN : (KLON)            ; SURFACE PLANCK FUNCTION
! PCLDLD : (KLON,KLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
! PCLDLU : (KLON,KLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE
! PEMIS  : (KLON)            ; SURFACE EMISSIVITY
! PFDN   : (KLON,KLEV+1)     ; CLEAR-SKY DOWNWARD FLUX
! PFUP   : (KLON,KLEV+1)     ; CLEAR-SKY UPWARD FLUX
!     ==== OUTPUTS ===
! PFLUX(KLON,2,KLEV)         ; RADIATIVE FLUXES :
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
!          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
!          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
!     CLOUDS
!
!     EXTERNALS.
!     ----------
!
!          NONE
!
!     REFERENCE.
!     ----------
!
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
!
!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOERAD   , ONLY : NOVLP
USE OYOERDI   , ONLY : REPCLC
USE YOEDBUG  , ONLY : LDEBUG


IMPLICIT NONE



!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON

!-----------------------------------------------------------------------
!
!*       0.1   ARGUMENTS
!              ---------
!
REAL_B :: PBINT(KLON,KLEV+1),PBSUIN(KLON),PCLDLD(KLON,KLEV) &
     &  ,  PCLDLU(KLON,KLEV) &
     &  ,  PCNTRB(KLON,KLEV+1,KLEV+1) &
     &  ,  PFDN(KLON,KLEV+1),PFUP(KLON,KLEV+1) &
     &  ,  PEMIS(KLON)
!
REAL_B :: PFLUX(KLON,2,KLEV+1)
!
!-----------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
!
INTEGER_M :: IMX(KLON), IMXP(KLON)
!
REAL_B :: ZCLEAR(KLON),ZCLOUD(KLON),ZDNF(KLON,KLEV+1,KLEV+1) &
     &  , ZFD(KLON), ZFN10(KLON), ZFU(KLON) &
     &  , ZUPF(KLON,KLEV+1,KLEV+1) &
     &  , ZCLM(KLON,KLEV+1,KLEV+1)
!
!     LOCAL INTEGER SCALARS
INTEGER_M :: IKCP1, IKM1, IKP1, IMAXC, IMXM1, IMXP1, JCLOUD,&
             &JK, JK1, JK2, JKJ, JL, IMX1, IMX2, JKC, JKCP1, JKM1, JKP1

!     LOCAL REAL SCALARS
REAL_B :: ZCFRAC

!     ------------------------------------------------------------------
!
!*         1.     INITIALIZATION
!                 --------------
!
IMAXC = 0
!
DO JL = KIDIA,KFDIA
  IMX(JL)=0
  IMXP(JL)=0
  ZCLOUD(JL) = 0.
END DO
!
!*         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
!                 -------------------------------------------
!
DO JK = 1 , KLEV
  DO JL = KIDIA,KFDIA
    IMX1=IMX(JL)
    IMX2=JK
    IF (PCLDLU(JL,JK).GT.REPCLC) THEN
      IMXP(JL)=IMX2
    ELSE
      IMXP(JL)=IMX1
    END IF
    IMAXC=MAX(IMXP(JL),IMAXC)
    IMX(JL)=IMXP(JL)
  END DO
END DO
!CGM*******
IMAXC=KLEV
!CGM*******
!
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFLUX(JL,1,JK) = PFUP(JL,JK)
    PFLUX(JL,2,JK) = PFDN(JL,JK)
  END DO
END DO
!
!     ------------------------------------------------------------------
!
!*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
!                  ---------------------------------------
!
IF (IMAXC.GT.0) THEN
!
  IMXP1 = IMAXC + 1
  IMXM1 = IMAXC - 1
!
!*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
!                  ------------------------------
!
  DO JK1=1,KLEV+1
    DO JK2=1,KLEV+1
      DO JL = KIDIA,KFDIA
        ZUPF(JL,JK2,JK1)=PFUP(JL,JK1)
        ZDNF(JL,JK2,JK1)=PFDN(JL,JK1)
      END DO
    END DO
  END DO
!
!*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
!                  ----------------------------------------------
!
  DO JKC = 1 , IMAXC
    JCLOUD=JKC
    JKCP1=JCLOUD+1
!
!*         2.1.1   ABOVE THE CLOUD
!                  ---------------
!
    DO JK=JKCP1,KLEV+1
      JKM1=JK-1
      DO JL = KIDIA,KFDIA
        ZFU(JL)=0.
      END DO
      IF (JK .GT. JKCP1) THEN
        DO JKJ=JKCP1,JKM1
          DO JL = KIDIA,KFDIA
            ZFU(JL) = ZFU(JL) + PCNTRB(JL,JK,JKJ)
          END DO
        END DO
      END IF
!
      DO JL = KIDIA,KFDIA
        ZUPF(JL,JKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
      END DO
    END DO
!
!*         2.1.2   BELOW THE CLOUD
!                  ---------------
!
    DO JK=1,JCLOUD
      JKP1=JK+1
      DO JL = KIDIA,KFDIA
        ZFD(JL)=0.
      END DO
!
      IF (JK .LT. JCLOUD) THEN
        DO JKJ=JKP1,JCLOUD
          DO JL = KIDIA,KFDIA
            ZFD(JL) = ZFD(JL) + PCNTRB(JL,JK,JKJ)
          END DO
        END DO
      END IF
      DO JL = KIDIA,KFDIA
        ZDNF(JL,JKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
      END DO
    END DO
!
  END DO
!
!
!*         2.2     CLOUD COVER MATRIX
!                  ------------------
!
!*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
!     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1
!
  DO JK1 = 1 , KLEV+1
    DO JK2 = 1 , KLEV+1
      DO JL = KIDIA,KFDIA
        ZCLM(JL,JK1,JK2) = 0.
      END DO
    END DO
  END DO
!
!
!
!*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
!                  ------------------------------------------
!
  DO JK1 = 2 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
    END DO
    DO JK = JK1 - 1 , 1 , -1
      DO JL = KIDIA,KFDIA
        IF (NOVLP.EQ.1) THEN
!* maximum-random       
          ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLU(JL,JK),ZCLOUD(JL))) &
          &                       /(1.0-MIN(ZCLOUD(JL),1.-REPCLC))
          ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
          ZCLOUD(JL) = PCLDLU(JL,JK)
        ELSE IF (NOVLP.EQ.2) THEN 
!* maximum      
          ZCLOUD(JL) = AMAX1(ZCLOUD(JL) , PCLDLU(JL,JK))
          ZCLM(JL,JK1,JK) = ZCLOUD(JL)
        ELSE IF (NOVLP.EQ.3) THEN
!* random      
          ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLU(JL,JK))
          ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
          ZCLM(JL,JK1,JK) = ZCLOUD(JL)
        END IF
      END DO
    END DO
  END DO
!
!
!*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
!                  ------------------------------------------
!
  DO JK1 = 1 , KLEV
    DO JL = KIDIA,KFDIA
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
    END DO
    DO JK = JK1 , KLEV
      DO JL = KIDIA,KFDIA
        IF (NOVLP.EQ.1) THEN
!* maximum-random       
           ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLD(JL,JK),ZCLOUD(JL))) &
           &                    /(1.0-MIN(ZCLOUD(JL),1.-REPCLC))
           ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
           ZCLOUD(JL) = PCLDLD(JL,JK)
         ELSE IF (NOVLP.EQ.2) THEN 
!* maximum      
           ZCLOUD(JL) = AMAX1(ZCLOUD(JL) , PCLDLD(JL,JK))
           ZCLM(JL,JK1,JK) = ZCLOUD(JL)
         ELSE IF (NOVLP.EQ.3) THEN
!* random      
           ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLD(JL,JK))
           ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
           ZCLM(JL,JK1,JK) = ZCLOUD(JL)
         END IF
       END DO
     END DO
   END DO
!
!
!
!*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
!                  ----------------------------------------------
!
!*         3.1     DOWNWARD FLUXES
!                  ---------------
!
  DO JL = KIDIA,KFDIA
    PFLUX(JL,2,KLEV+1) = 0.
  END DO
!
  DO JK1 = KLEV , 1 , -1
!
!*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
!
    DO JL = KIDIA,KFDIA
      ZFD (JL) = (1. - ZCLM(JL,JK1,KLEV)) * ZDNF(JL,1,JK1)
!
!*                 CONTRIBUTION FROM ADJACENT CLOUD
!
      ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
    END DO
!
!*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
!
    DO JK = KLEV-1 , JK1 , -1
      DO JL = KIDIA,KFDIA
        ZCFRAC = ZCLM(JL,JK1,JK+1) - ZCLM(JL,JK1,JK)
        ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
      END DO
    END DO
!
    DO JL = KIDIA,KFDIA
      PFLUX(JL,2,JK1) = ZFD (JL)
    END DO
!
  END DO
!
!
!
!
!*         3.2     UPWARD FLUX AT THE SURFACE
!                  --------------------------
!
  DO JL = KIDIA,KFDIA
    PFLUX(JL,1,1) = PEMIS(JL)*PBSUIN(JL)-(1.-PEMIS(JL))*PFLUX(JL,2,1)
  END DO
!
!
!
!*         3.3     UPWARD FLUXES
!                  -------------
!
  DO JK1 = 2 , KLEV+1
!
!*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
!
    DO JL = KIDIA,KFDIA
      ZFU (JL) = (1. - ZCLM(JL,JK1,1)) * ZUPF(JL,1,JK1)
!
!*                 CONTRIBUTION FROM ADJACENT CLOUD
!
      ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)
    END DO
!
!*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
!
    DO JK = 2 , JK1-1
      DO JL = KIDIA,KFDIA
        ZCFRAC = ZCLM(JL,JK1,JK-1) - ZCLM(JL,JK1,JK)
        ZFU(JL) =  ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
      END DO
    END DO
!
    DO JL = KIDIA,KFDIA
      PFLUX(JL,1,JK1) = ZFU (JL)
    END DO
!
  END DO
!
!
END IF
!
!
!*         2.3     END OF CLOUD EFFECT COMPUTATIONS
!
!----------------------------------------------------------------------
RETURN
END SUBROUTINE OLWC
