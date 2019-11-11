!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWB &
     &  ( KIDIA, KFDIA, KLON  , KLEV  &
     &  , PDT0 , PT   , PTH &
     &  , PB   , PBINT, PBSUIN, PBSUR , PBTOP , PDBSL &
     &  , PGA  , PGB  , PGASUR, PGBSUR, PGATOP, PGBTOP )
!
!**** *LWB*   - COMPUTES BLACK-BODY FUNCTIONS FOR LONGWAVE CALCULATIONS
!
!     PURPOSE.
!     --------
!           COMPUTES PLANCK FUNCTIONS
!
!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PDT0   : (KLON)             ; SURFACE TEMPERATURE DISCONTINUITY
! PT     : (KLON,KLEV)        ; TEMPERATURE
! PTH    : (KLON,KLEV+1)      ; HALF LEVEL TEMPERATURE
!     ==== OUTPUTS ===
! PB     : (KLON,NISP,KLEV+1) ; SPECTRAL HALF LEVEL PLANCK FUNCTION
! PBINT  : (KLON,KLEV+1)      ; HALF LEVEL PLANCK FUNCTION
! PBSUIN : (KLON)             ; SURFACE PLANCK FUNCTION
! PBSUR  : (KLON,NISP)        ; SURFACE SPECTRAL PLANCK FUNCTION
! PBTOP  : (KLON,NISP)        ; TOP SPECTRAL PLANCK FUNCTION
! PDBSL  : (KLON,NISP,KLEV*2); SUB-LAYER PLANCK FUNCTION GRADIENT
! PGA    : (KLON,8,2,KLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
! PGB    : (KLON,8,2,KLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
! PGASUR, PGBSUR (KLON,8,2)   ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP (KLON,8,2)   ; T.O.A. PADE APPROXIMANTS
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
!     FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION
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
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS           "
!
!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOEOLW   , ONLY : MXIXT    ,NISP     ,NIPD     ,  GA     ,&
            & GB     ,TINTP    ,TSTAND   ,TSTP     ,XP

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
REAL_B :: PDT0(KLON), PT(KLON,KLEV), PTH(KLON,KLEV+1)
!
REAL_B :: PB(KLON,NISP,KLEV+1), PBINT(KLON,KLEV+1) &
     &  ,  PBSUIN(KLON)       , PBSUR(KLON,NISP) &
     &  ,  PBTOP(KLON,NISP)   , PDBSL(KLON,NISP,KLEV*2) &
     &  ,  PGA(KLON,8,2,KLEV) , PGB(KLON,8,2,KLEV) &
     &  ,  PGASUR(KLON,8,2)   , PGBSUR(KLON,8,2) &
     &  ,  PGATOP(KLON,8,2)   , PGBTOP(KLON,8,2)
!
!-------------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
INTEGER_M :: INDB(KLON),INDS(KLON)

REAL_B :: ZBLAY(KLON,KLEV),ZBLEV(KLON,KLEV+1) &
     &  , ZRES(KLON),ZRES2(KLON),ZTI(KLON),ZTI2(KLON)
     
REAL_B :: ZDST1, ZDSTO1, ZDSTX, ZDSTOX
     
INTEGER_M :: ILEV2, INDSU, INDT, INDTP, INDTO, INUS, INUE, IXTOX, IXTX &
     &  , JF, JG, JNU, JK, JK1, JK2, JL, IKL
!
!     ------------------------------------------------------------------
!
!
!*         1.0     PLANCK FUNCTIONS AND GRADIENTS
!                  ------------------------------
!
ILEV2=2*KLEV
INUS=1
INUE=NISP
       
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PBINT(JL,JK) = 0.
  END DO
END DO
DO JNU=1,NISP
  DO JL=KIDIA,KFDIA
    PBSUR(JL,JNU)=0.
    PBTOP(JL,JNU)=0.
  END DO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PB(JL,JNU,JK)=0.
    END DO
  END DO
  DO JK=1,ILEV2
    DO JL=KIDIA,KFDIA
      PDBSL(JL,JNU,JK)=0.
    END DO  
  END DO
END DO
DO JL = KIDIA,KFDIA
  PBSUIN(JL) = 0.
END DO
!
DO JNU=INUS,INUE
!
!
!*         1.1   LEVELS FROM SURFACE TO KLEV
!                ----------------------------
!                TEMPERATURE ENTERED FROM TOP TO BOTTOM
!
  DO JK = 1 , KLEV
    IKL=KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZTI(JL)=(PTH(JL,IKL+1)-TSTAND)/TSTAND
      ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU) &
      &       +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU) &
      &       )))))
      PBINT(JL,JK)=PBINT(JL,JK)+ZRES(JL)
      PB(JL,JNU,JK)= ZRES(JL)
      ZBLEV(JL,JK) = ZRES(JL)
      
      ZTI2(JL)=(PT(JL,IKL)-TSTAND)/TSTAND
      ZRES2(JL)=XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU) &
      &     +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU) &
      &       )))))
      ZBLAY(JL,JK) = ZRES2(JL)
    END DO
  END DO
!
!
!*         1.2   TOP OF THE ATMOSPHERE AND SURFACE
!                ---------------------------------
!                TEMPERATURE ENTERED FROM TOP TO BOTTOM
!
  DO JL = KIDIA,KFDIA
    ZTI(JL)=(PTH(JL,1)-TSTAND)/TSTAND
    ZTI2(JL) = (PTH(JL,KLEV+1) + PDT0(JL) - TSTAND) / TSTAND
    ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU) &
    &   +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU) &
    &      )))))
    ZRES2(JL) = XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU) &
    &    +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU) &
    &       )))))
    PBINT(JL,KLEV+1) = PBINT(JL,KLEV+1)+ZRES(JL)
    PB(JL,JNU,KLEV+1)= ZRES(JL)
    ZBLEV(JL,KLEV+1) = ZRES(JL)
    PBTOP(JL,JNU) = ZRES(JL)
    PBSUR(JL,JNU) = ZRES2(JL)
    PBSUIN(JL) = PBSUIN(JL) + ZRES2(JL)
  END DO
!
!
!*         1.3   GRADIENTS IN SUB-LAYERS
!                -----------------------
!
  DO JK = 1 , KLEV
    JK2 = 2 * JK
    JK1 = JK2 - 1
    DO JL = KIDIA,KFDIA
      PDBSL(JL,JNU,JK1) = ZBLAY(JL,JK  ) - ZBLEV(JL,JK)
      PDBSL(JL,JNU,JK2) = ZBLEV(JL,JK+1) - ZBLAY(JL,JK)
    END DO
  END DO
!
END DO
!
!*         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
!                ---------------------------------------------
!
DO JL=KIDIA,KFDIA
  ZDSTO1 = (PTH(JL,1)-TINTP(1)) / TSTP
  IXTOX = MAX( 1, MIN( MXIXT, INT( ZDSTO1 + 1. ) ) )
  ZDSTOX = (PTH(JL,1)-TINTP(IXTOX))/TSTP
  IF (ZDSTOX.LT.0.5) THEN
    INDTO=IXTOX
  ELSE
    INDTO=IXTOX+1
  END IF
  INDB(JL)=INDTO
  
  ZDST1 = (PTH(JL,KLEV+1)-TINTP(1)) / TSTP
  IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
  ZDSTX = (PTH(JL,KLEV+1)-TINTP(IXTX))/TSTP
  IF (ZDSTX.LT.0.5) THEN
    INDT=IXTX
  ELSE
    INDT=IXTX+1
  END IF
  INDS(JL)=INDT
END DO
!
DO JF=1,2
  DO JG=1, 8
    DO JL=KIDIA,KFDIA
      INDSU=INDS(JL)
      PGASUR(JL,JG,JF)=GA(INDSU,2*JG-1,JF)
      PGBSUR(JL,JG,JF)=GB(INDSU,2*JG-1,JF)
      INDTP=INDB(JL)
      PGATOP(JL,JG,JF)=GA(INDTP,2*JG-1,JF)
      PGBTOP(JL,JG,JF)=GB(INDTP,2*JG-1,JF)
    END DO
  END DO
END DO
!
DO JK=1,KLEV
  IKL=KLEV+1-JK
  DO JL=KIDIA,KFDIA
    ZDST1 = (PT(JL,IKL)-TINTP(1)) / TSTP
    IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
    ZDSTX = (PT(JL,IKL)-TINTP(IXTX))/TSTP
    IF (ZDSTX.LT.0.5) THEN
      INDT=IXTX
    ELSE
      INDT=IXTX+1
    END IF
    INDB(JL)=INDT
  END DO
!
  DO JF=1,2
    DO JG=1, 8
      DO JL=KIDIA,KFDIA
        INDT=INDB(JL)
        PGA(JL,JG,JF,JK)=GA(INDT,2*JG,JF)
        PGB(JL,JG,JF,JK)=GB(INDT,2*JG,JF)
      END DO
    END DO
  END DO
END DO
!
!     ------------------------------------------------------------------
!
RETURN
END SUBROUTINE OLWB
