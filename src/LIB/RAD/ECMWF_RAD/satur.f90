!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:38
!-----------------------------------------------------------------
SUBROUTINE SATUR ( KIDIA , KFDIA , KLON  , KTDIA , KLEV &
                 &, PAPRSF, PT    , PQSAT , KFLAG)

!***

! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION


!       J.F. MAHFOUF       E.C.M.W.F.     15/05/96


!       PURPOSE.
!       --------

!       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
!       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
!       AND LIQUID WATER CONTENT  

!       INTERFACE
!       ---------

!       THIS ROUTINE IS CALLED FROM *CALLPAR*.


!       PARAMETER     DESCRIPTION                                 UNITS
!       ---------     -----------                                 -----
!       INPUT PARAMETERS (INTEGER):

!      *KIDIA*        START POINT
!      *KFDIA*        END POINT
!      *KLON*         NUMBER OF GRID POINTS PER PACKET
!      *KTDIA*        START OF THE VERTICAL LOOP
!      *KLEV*         NUMBER OF LEVELS


!       INPUT PARAMETERS (REAL):

!      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
!      *PT*            TEMPERATURE AT T-DT                          K

!       INPUT PARAMETERS (INTEGER):

!      *KFLAG*         FLAG TO DETECT CALL FROM

!                      CONVECTION  KFLAG=1
!                      OTHER       KFLAG=2

!       OUTPUT PARAMETER (REAL):

!      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG

!-------------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
            &R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
            &RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU
USE YOEPHLI  , ONLY : LPHYLIN  ,RLPTRC   ,RLPAL1   ,RLPAL2


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KFLAG
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KTDIA




REAL_B :: PAPRSF(KLON,KLEV), PT(KLON,KLEV) , PQSAT(KLON,KLEV)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JK, JL

!     LOCAL REAL SCALARS
REAL_B :: Z3ES, Z4ES, ZCOR, ZEW, ZFOEEW, ZQMAX, ZQS, ZTARG


!!!!!!!   !DIR$ VFUNCTION EXPHF

#include "fcttre.h"

!----------------------------------------------------------------------

!*    1.           DEFINE CONSTANTS
!                  ----------------

ZQMAX=_HALF_

!     *
!----------------------------------------------------------------------

!     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
!                       --------------------------------------

DO JK=KTDIA,KLEV
  DO JL=KIDIA, KFDIA

    IF (LPHYLIN) THEN
      ZTARG = PT(JL,JK)
      IF (ZTARG > RTT) THEN
        Z3ES=R3LES
        Z4ES=R4LES
      ELSE
        Z3ES=R3IES
        Z4ES=R4IES
      ENDIF
      ZFOEEW = R2ES*EXP(Z3ES*(ZTARG-RTT)/(ZTARG-Z4ES))
      ZQS    = ZFOEEW/PAPRSF(JL,JK)
      IF (ZQS > ZQMAX) THEN
        ZQS=ZQMAX
      ENDIF
    ELSE
      IF(KFLAG == 1) THEN
        ZEW  = FOEEWMCU(PT(JL,JK))
      ELSE
        ZEW  = FOEEWM(PT(JL,JK))
      ENDIF
      ZQS  = ZEW/PAPRSF(JL,JK)
      ZQS  = MIN(ZQMAX,ZQS)
    ENDIF

    ZCOR = _ONE_/(_ONE_-RETV*ZQS)
    PQSAT(JL,JK)=ZQS*ZCOR
  ENDDO
ENDDO

RETURN
END SUBROUTINE SATUR
