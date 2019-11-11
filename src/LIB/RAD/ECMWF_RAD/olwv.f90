!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWV ( KIDIA,KFDIA,KLON,KLEV &
     &  , KUAER,KTRAER &
     &  , PABCU,PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL,PEMIS &
     &  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP &
     &  , PCNTRB,PCOLC,PFLUC                 )
!
!**** *LWV*   - LONGWAVE RADIATION, VERTICAL INTEGRATION
!
!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
!           FLUXES OR RADIANCES
!
!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PB     : (KLON,NISP,KLEV+1); SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
! PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
! PBSUIN : (KLON)            ; SURFACE PLANCK FUNCTION
! PBSUR  : (KLON,NISP)       ; SURFACE SPECTRAL PLANCK FUNCTION
! PBTOP  : (KLON,NISP)       ; T.O.A. SPECTRAL PLANCK FUNCTION
! PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PEMIS  : (KLON)            ; SURFACE EMISSIVITY
! PGA, PGB                   ; PADE APPROXIMANTS
! PGASUR, PGBSUR             ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP             ; T.O.A. PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR-SKY:
!                     1  ==>  UPWARD   FLUX TOTAL
! PCOLC(KLON,KLEV)           ; LONG-WAVE TENDENCY CLEAR SKY
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
!     CONTRIBUTIONS BY -  THE NEARBY LAYERS
!                      -  THE DISTANT LAYERS
!                      -  THE BOUNDARY TERMS
!          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
!
!     EXTERNALS.
!     ----------
!
!          *LWVN*, *LWVD*, *LWVB*
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

USE YOEOLW   , ONLY : NISP     ,NIPD     ,NUA


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KUAER, KTRAER

!-----------------------------------------------------------------------
!
!*       0.1   ARGUMENTS
!              ---------
!
!
REAL_B :: PABCU(KLON,NUA,3*KLEV+1), PB(KLON,NISP,KLEV+1) &
     &  ,  PBINT(KLON,KLEV+1) &
     &  ,  PBSUR(KLON,NISP), PBSUIN(KLON), PBTOP(KLON,NISP) &
     &  ,  PDBSL(KLON,NISP,KLEV*2), PEMIS(KLON) &
     &  ,  PGA(KLON,8,2,KLEV), PGB(KLON,8,2,KLEV) &
     &  ,  PGASUR(KLON,8,2)  , PGBSUR(KLON,8,2) &
     &  ,  PGATOP(KLON,8,2)  , PGBTOP(KLON,8,2)
!
REAL_B :: PCNTRB(KLON,KLEV+1,KLEV+1), PCOLC(KLON,KLEV) &
     &  ,  PFLUC(KLON,2,KLEV+1)
!
!-----------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
!
INTEGER_M :: ITX(KLON)
!
REAL_B :: ZADJD(KLON,KLEV+1), ZADJU(KLON,KLEV+1) &
     &  ,  ZDBDT(KLON,NISP,KLEV) &
     &  ,  ZDISD(KLON,KLEV+1), ZDISU(KLON,KLEV+1) &
     &  ,  ZFD(KLON), ZFDN(KLON,KLEV+1), ZFU(KLON) &
     &  ,  ZFUP(KLON,KLEV+1),ZGLAYD(KLON),ZGLAYU(KLON)
     
REAL_B :: ZDFNET     
!
!     LOCAL INTEGER SCALARS
INTEGER_M :: JA, JK, JL, JKL

!-----------------------------------------------------------------------
!
!*         1.    INITIALIZATION
!                --------------
!
!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------
!
DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    ZADJD(JL,JK)=0.
    ZADJU(JL,JK)=0.
    ZDISD(JL,JK)=0.
    ZDISU(JL,JK)=0.
  END DO
END DO
!
!     ------------------------------------------------------------------
!
!*         2.      VERTICAL INTEGRATION
!                  --------------------
!
!     ------------------------------------------------------------------
!
!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------
!
CALL OLWVN ( KIDIA,KFDIA,KLON,KLEV,KUAER &
     &  , PABCU,PDBSL,PGA,PGB &
     &  , ZADJD,ZADJU,PCNTRB,ZDBDT )
!
!     ------------------------------------------------------------------
!
!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------
!
CALL OLWVD ( KIDIA,KFDIA,KLON,KLEV,KTRAER &
     &  , PABCU,ZDBDT,PGA,PGB &
     &  , PCNTRB,ZDISD,ZDISU )
!
!     ------------------------------------------------------------------
!
!*         2.3     EXCHANGE WITH THE BOUNDARIES
!                  ----------------------------
!
CALL OLWVB ( KIDIA,KFDIA,KLON,KLEV,KUAER &
     &  , PABCU,ZADJD,ZADJU,PB,PBINT,PBSUIN,PBSUR,PBTOP &
     &  , ZDISD,ZDISU,PEMIS &
     &  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP &
     &  , PFLUC                     )
!
!
!     ------------------------------------------------------------------
RETURN
END SUBROUTINE OLWV
