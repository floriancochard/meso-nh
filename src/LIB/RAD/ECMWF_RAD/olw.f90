!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLW &
     & ( KIDIA, KFDIA , KLON  , KLEV &
     & , PCCO2, PCLDLD, PCLDLU &
     & , PDP  , PDT0  , PEMIS  &
     & , PAPH , PQOF  , PTH &
     & , PAER , PT    , PVIEW , PWV &
     & , PCOLR, PCOLC , PFLUX, PFLUC &
     & )
!
!**** *LW*   - ORGANIZES THE LONGWAVE CALCULATIONS
!
!     PURPOSE.
!     --------
!           DEPENDING ON KMODE, COMPUTES LONGWAVE FLUXES AND/OR
!           RADIANCES
!
!**   INTERFACE.
!     ----------
!
!        *LW* IS CALLED FROM *RADLSW*
!
!        EXPLICIT ARGUMENTS :
!        --------------------
! PAER   : (KLON,KLEV,6)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PCCO2  :                   ; CONCENTRATION IN CO2 (PA/PA)
! PCLDLD : (KLON,KLEV)       ; DOWNWARD EFFECTIVE FRACTIONAL COVER
! PCLDLU : (KLON,KLEV)       ; UPWARD EFFECTIVE FRACTIONAL COVER
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PDT0   : (KLON)            ; SURFACE TEMPERATURE DISCONTINUITY  
! PEMIS  : (KLON)            ; SURFACE EMISSIVITY
! PAPH   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
! PT     : (KLON,KLEV)       ; TEMPERATURE
! PTH    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
! PVIEW  : (KLON)            ; COSECANT OF VIEWING ANGLE
! PWV    : (KLON,KLEV)       ; SPECIFIC HUMIDITY  (PA/PA)
!     ==== OUTPUTS ===
!  IF KMODE = 0, 1, 2
! PFLUX(KLON,2,KLEV)         ; RADIATIVE FLUXES :
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL
! PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR SKY:
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL
! PCOLR(KLON,KLEV)           ; LONG-WAVE TENDENCY
! PCOLC(KLON,KLEV)           ; LONG-WAVE TENDENCY CLEAR SKY
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
!     ABSORBERS.
!          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
!     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
!          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
!     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
!     BOUNDARIES.
!          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
!          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.
!
!     EXTERNALS.
!     ----------
!
!          *LWU*, *LWBV*, *LWC*
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

USE YOEOLW   , ONLY : NUA
USE YOERDU   , ONLY : RCDAY
USE YOEDBUG  , ONLY : LDEBUG


IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON

!     DUMMY REAL SCALARS
REAL_B :: PCCO2, ZDFNET

!-----------------------------------------------------------------------
!
!*       0.1   ARGUMENTS
!              ---------
!
REAL_B :: PCLDLD(KLON,KLEV)      , PCLDLU(KLON,KLEV) &
     &  ,  PDP(KLON,KLEV)   , PDT0(KLON) &
     &  ,  PEMIS(KLON)      , PAPH(KLON,KLEV+1) &
     &  ,  PQOF(KLON,KLEV)  , PTH(KLON,KLEV+1) &
     &  ,  PAER(KLON,KLEV,6), PT(KLON,KLEV) &
     &  ,  PVIEW(KLON)      , PWV(KLON,KLEV)
!
REAL_B :: PCOLR(KLON,KLEV)       , PCOLC(KLON,KLEV) &
     &  ,  PFLUX(KLON,2,KLEV+1), PFLUC(KLON,2,KLEV+1)
!
!-------------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
REAL_B :: ZABCU(KLON,NUA,3*KLEV+1) &
     &  ,  ZBINT(KLON,KLEV+1) &
     &  ,  ZBSUI(KLON) &
     &  ,  ZCNTRB(KLON,KLEV+1,KLEV+1) &
     &  ,  ZFDN(KLON,KLEV+1) &
     &  ,  ZFUP(KLON,KLEV+1)
     
!     LOCAL INTEGER SCALARS
INTEGER_M :: JK, JKL, JL, ILIM, IUA, KLEVT, JK1, JK2

     
!
!     ------------------------------------------------------------------
!
!*         1.    INITIALIZATION
!                --------------
!
!     ------------------------------------------------------------------
!
!*         1.1   COMPUTES ABSORBER AMOUNTS
!                -------------------------
!

CALL OLWU ( KIDIA, KFDIA, KLON, KLEV &
     &  ,  PAER, PCCO2, PDP, PAPH, PQOF, PT, PVIEW, PWV &
     &  ,  ZABCU   )

!
!     ------------------------------------------------------------------
!
!*         2.    COMPUTES PLANCK FUNCTIONS
!                -------------------------
!                PERFORMS THE VERTICAL INTEGRATION
!                ---------------------------------

DO JK1=1,KLEV+1
  DO JK2=1,KLEV+1
    DO JL=KIDIA,KFDIA
      ZCNTRB(JL,JK1,JK2)=0.
    END DO
  END DO
END DO      
!
CALL OLWBV( KIDIA,KFDIA,KLON,KLEV &
     &          , PDP,PDT0,PEMIS,PTH &
     &          , PT &
     &          , PCOLC,PFLUC &
     &          , ZABCU,ZBINT,ZBSUI,ZCNTRB,ZFDN,ZFUP)
     
!     ------------------------------------------------------------------
!
!*         4.    INTRODUCES THE EFFECTS OF CLOUDS
!                --------------------------------
!
CALL OLWC ( KIDIA,KFDIA,KLON,KLEV &
     &  , ZBINT,ZBSUI,PCLDLD,PCLDLU,ZCNTRB,PEMIS,ZFDN,ZFUP &
     &  , PFLUX                                                )
!
DO JKL = 1 , KLEV
  JK = KLEV+1 - JKL

  DO JL = KIDIA,KFDIA
    ZDFNET = PFLUX(JL,1,JK+1) + PFLUX(JL,2,JK+1) &
     &        -PFLUX(JL,1,JK  ) - PFLUX(JL,2,JK  )
    PCOLR(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
    
    ZDFNET = PFLUC(JL,1,JK+1) + PFLUC(JL,2,JK+1) &
    &              -PFLUC(JL,1,JK  ) - PFLUC(JL,2,JK  )
    PCOLC(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
  END DO
  
END DO
!
!     ------------------------------------------------------------------
!
RETURN
END SUBROUTINE OLW
