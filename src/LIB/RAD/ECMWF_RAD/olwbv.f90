!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWBV &
     &  ( KIDIA, KFDIA, KLON , KLEV &
     &  , PDP  , PDT0 , PEMIS, PTH &
     &  , PT   &
     &  , PCOLC, PFLUC &
     &  , PABCU, PBINT, PBSUI, PCNTRB, PFDN, PFUP )
!
!**** *LWBV*   - COMPUTE PLANCK FUNC., PERF. VERT. INTEGRATION
!
!     PURPOSE.
!     --------
!           TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
!           VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
!           SAVING
!
!**   INTERFACE.
!     ----------
!
!        *LWVB* IS CALLED FROM *LW*
!
!        EXPLICIT ARGUMENTS :
!        --------------------
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PDT0   : (KLON)            ; SURFACE TEMPERATURE DISCONTINUITY
! PEMIS  : (KLON)            ; SURFACE EMISSIVITY
! PT     : (KLON,KLEV)       ; TEMPERATURE
! PTH    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
!     ==== OUTPUTS ===
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
!     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
!          2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
!     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
!     BOUNDARIES.
!          3. COMPUTES THE CLEAR-SKY COOLING RATES.
!
!     EXTERNALS.
!     ----------
!
!          *LWB*, *LWV*
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
!        MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
!                                          MEMORY)
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOEOLW   , ONLY : NISP     ,NIPD     ,NUA
USE YOERDU   , ONLY : NUAER    ,NTRAER   ,RCDAY
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
REAL_B :: PDP(KLON,KLEV) &
     &  ,  PDT0(KLON)        ,PEMIS(KLON) &
     &  ,  PTH(KLON,KLEV+1) &
     &  ,  PT(KLON,KLEV)
!
REAL_B :: PCOLC(KLON,KLEV), PFLUC(KLON,2,KLEV+1)
!     
REAL_B :: PABCU(KLON,NUA,3*KLEV+1) &
     &  ,  PBINT(KLON,KLEV+1) &
     &  ,  PBSUI(KLON) &
     &  ,  PCNTRB(KLON,KLEV+1,KLEV+1) &
     &  ,  PFDN(KLON,KLEV+1) &
     &  ,  PFUP(KLON,KLEV+1)
!
!-------------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
REAL_B :: ZB(KLON,NISP,KLEV+1) &
     &  ,  ZBSUR(KLON,NISP), ZBTOP(KLON,NISP) &
     &  ,  ZDBSL(KLON,NISP,KLEV*2) &
     &  ,  ZGA(KLON,8,2,KLEV), ZGB(KLON,8,2,KLEV) &
     &  ,  ZGASUR(KLON,8,2)   , ZGBSUR(KLON,8,2) &
     &  ,  ZGATOP(KLON,8,2)   , ZGBTOP(KLON,8,2)
     
REAL_B :: ZDFNET     
!
!     LOCAL INTEGER SCALARS
INTEGER_M :: JK, JKL, JL, JLW

!     ------------------------------------------------------------------
!
!*         2.    COMPUTES PLANCK FUNCTIONS
!                -------------------------
!
if (LDEBUG) print *, 'CALL OLWB'
CALL OLWB ( KIDIA, KFDIA, KLON  , KLEV  &
     &    , PDT0 , PT   , PTH &
     &    , ZB   , PBINT, PBSUI , ZBSUR , ZBTOP , ZDBSL &
     &    , ZGA  , ZGB  , ZGASUR, ZGBSUR, ZGATOP, ZGBTOP     )
!
!     ------------------------------------------------------------------
!
!*         3.    PERFORMS THE VERTICAL INTEGRATION
!                ---------------------------------
!
if (LDEBUG) print *, 'CALL OLWV'
CALL OLWV ( KIDIA, KFDIA, KLON , KLEV , NUAER, NTRAER &
     &    , PABCU, ZB   , PBINT, PBSUI, ZBSUR, ZBTOP, ZDBSL &
     &    , PEMIS &
     &    , ZGA  , ZGB  , ZGASUR,ZGBSUR,ZGATOP,ZGBTOP &
     &    , PCNTRB,PCOLC, PFLUC )
!
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFDN(JL,JK) = PFLUC(JL,2,JK)
    PFUP(JL,JK) = PFLUC(JL,1,JK)
  END DO
END DO
!
!
DO JKL = 1 , KLEV
  JK = KLEV+1 - JKL
  DO JL = KIDIA,KFDIA
    ZDFNET = PFLUC(JL,1,JK+1) + PFLUC(JL,2,JK+1) &
    &        -PFLUC(JL,1,JK  ) - PFLUC(JL,2,JK  )
    PCOLC(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
  END DO
END DO
!
!     ------------------------------------------------------------------
!
RETURN
END SUBROUTINE OLWBV
