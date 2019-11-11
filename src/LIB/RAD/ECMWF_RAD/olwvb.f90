!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWVB ( KIDIA,KFDIA,KLON,KLEV,KUAER &
     &  , PABCU,PADJD,PADJU,PB,PBINT,PBSUI,PBSUR,PBTOP &
     &  , PDISD,PDISU,PEMIS &
     &  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP &
     &  , PFLUC )
!
!**** *LWVB*   - L.W., VERTICAL INTEGRATION, EXCHANGE WITH BOUNDARIES
!
!     PURPOSE.
!     --------
!           INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
!           INTEGRATION
!
!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PADJ.. : (KLON,KLEV+1)     ; CONTRIBUTION BY ADJACENT LAYERS
! PB     : (KLON,NISP,KLEV+1); SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
! PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
! PBSUI  : (KLON)            ; SURFACE PLANCK FUNCTION
! PBSUR  : (KLON,NISP)       ; SPECTRAL SURFACE PLANCK FUNCTION
! PBTOP  : (KLON,NISP)       ; SPECTRAL T.O.A. PLANCK FUNCTION
! PDIS.. : (KLON,KLEV+1)     ; CONTRIBUTION BY DISTANT LAYERS
! PEMIS  : (KLON)            ; SURFACE EMISSIVITY
! PGA, PGB                   ; PADE APPROXIMANTS
! PGASUR, PGBSUR             ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP             ; T.O.A. PADE APPROXIMANTS 
!     ==== OUTPUTS ===
! PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR-SKY:
!                     1  ==>  UPWARD   FLUX TOTAL
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
!     ATMOSPHERE
!          2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
!     TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
!          3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES
!
!     EXTERNALS.
!     ----------
!
!          *LWTT*
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

USE YOEOLW   , ONLY : NISP     ,NIPD     ,NTRA     ,NUA      ,NG1P1


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KUAER


!-----------------------------------------------------------------------
!
!*       0.1   ARGUMENTS
!              ---------
!
!
REAL_B :: PABCU(KLON,NUA,3*KLEV+1) &
     &  ,  PADJD(KLON,KLEV+1), PADJU(KLON,KLEV+1) &
     &  ,  PB(KLON,NISP,KLEV+1), PBINT(KLON,KLEV+1) &
     &  ,  PBSUR(KLON,NISP), PBSUI(KLON), PBTOP(KLON,NISP) &
     &  ,  PDISD(KLON,KLEV+1), PDISU(KLON,KLEV+1) &
     &  ,  PEMIS(KLON) &
     &  ,  PGA(KLON,8,2,KLEV), PGB(KLON,8,2,KLEV) &
     &  ,  PGASUR(KLON,8,2), PGBSUR(KLON,8,2) &
     &  ,  PGATOP(KLON,8,2), PGBTOP(KLON,8,2)
!
REAL_B :: PFLUC(KLON,2,KLEV+1) 
!
!-----------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
!
INTEGER_M :: ITX(KLON)
!
REAL_B :: ZBGND(KLON), ZFD(KLON), ZFDN(KLON,KLEV+1) &
     &  ,  ZFN10(KLON), ZFU(KLON), ZFUP(KLON,KLEV+1) &
     &  ,  ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA) &
     &  ,  ZUU(KLON,NUA) , ZCNSOL(KLON), ZCNTOP(KLON)
!

!     LOCAL INTEGER SCALARS
INTEGER_M :: IN, JA, JK, JL, IND1, IND2, IND3, IND4, JLIM

!     LOCAL REAL SCALARS
REAL_B :: ZCNTOP1, ZCNTOP2, ZCNTOP3, ZCNTOP4, ZCNTOP5, ZCNTOP6

!-----------------------------------------------------------------------
!
!*         1.    INITIALIZATION
!                --------------
!
!
!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------
!
DO JA=1,NTRA
  DO JL=KIDIA,KFDIA
    ZTT (JL,JA)=1.0
    ZTT1(JL,JA)=1.0
    ZTT2(JL,JA)=1.0
  END DO
END DO
!
DO JA=1,NUA
  DO JL=KIDIA,KFDIA
    ZUU(JL,JA)=1.0
  END DO
END DO
!
!     ------------------------------------------------------------------
!
!*         2.      VERTICAL INTEGRATION
!                  --------------------
!
IND1=0
IND3=0
IND4=1
IND2=1
!
!
!*         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
!                  -----------------------------------
!
DO JK = 1 , KLEV
  IN=(JK-1)*NG1P1+1
!
  DO JA=1,KUAER
    DO JL=KIDIA,KFDIA
      ZUU(JL,JA)=PABCU(JL,JA,IN)
    END DO
  END DO
!
!
  CALL LWTT ( KIDIA,KFDIA,KLON &
     &          , PGATOP(1,1,1), PGBTOP(1,1,1), ZUU, ZTT )
!
  DO JL = KIDIA,KFDIA
    ZCNTOP(JL)=PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10) &
    &      +PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11) &
    &      +PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12) &
    &      +PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13) &
    &      +PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14) &
    &      +PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
    ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
    ZFDN(JL,JK)=ZFD(JL)
    PFLUC(JL,2,JK)=ZFD(JL)
  END DO
!
END DO
!
JK = KLEV+1
IN=(JK-1)*NG1P1+1
!
DO JL = KIDIA,KFDIA
  ZCNTOP(JL)= PBTOP(JL,1) &
  &   + PBTOP(JL,2) &
  &   + PBTOP(JL,3) &
  &   + PBTOP(JL,4) &
  &   + PBTOP(JL,5) &
  &   + PBTOP(JL,6)
  ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
  ZFDN(JL,JK)=ZFD(JL)
  PFLUC(JL,2,JK)=ZFD(JL)
END DO
!
!*         2.4     COOLING-TO-SPACE OF LAYERS ABOVE 10 HPA
!                  ---------------------------------------
!
!
!*         2.4.1   INITIALIZATION
!                  --------------
!
!
!*         2.5     EXCHANGE WITH LOWER LIMIT
!                  -------------------------
!
DO JL = KIDIA,KFDIA
  ZBGND(JL)=PBSUI(JL)*PEMIS(JL)-(1.-PEMIS(JL)) &
  &               *PFLUC(JL,2,1)-PBINT(JL,1)
END DO
!
JK = 1
IN=(JK-1)*NG1P1+1
!
DO JL = KIDIA,KFDIA
  ZCNSOL(JL)=PBSUR(JL,1) &
  & +PBSUR(JL,2) &
  & +PBSUR(JL,3) &
  & +PBSUR(JL,4) &
  & +PBSUR(JL,5) &
  & +PBSUR(JL,6)
  ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
  ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
  ZFUP(JL,JK)=ZFU(JL)
  PFLUC(JL,1,JK)=ZFU(JL)
END DO
!
DO JK = 2 , KLEV+1
  IN=(JK-1)*NG1P1+1
!
!
  DO JA=1,KUAER
    DO JL=KIDIA,KFDIA
      ZUU(JL,JA)=PABCU(JL,JA,1)-PABCU(JL,JA,IN)
    END DO
  END DO
!
!
  CALL LWTT ( KIDIA,KFDIA,KLON &
     &          , PGASUR(1,1,1), PGBSUR(1,1,1), ZUU, ZTT )
!
  DO JL = KIDIA,KFDIA
    ZCNSOL(JL)=PBSUR(JL,1)*ZTT(JL,1)          *ZTT(JL,10) &
    &      +PBSUR(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11) &
    &      +PBSUR(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12) &
    &      +PBSUR(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13) &
    &      +PBSUR(JL,5)*ZTT(JL,3)          *ZTT(JL,14) &
    &      +PBSUR(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
    ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
    ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
    ZFUP(JL,JK)=ZFU(JL)
    PFLUC(JL,1,JK)=ZFU(JL)
  END DO
!
!
END DO
!
!
!
!*         2.7     CLEAR-SKY FLUXES
!                  ----------------
!
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFLUC(JL,1,JK) = ZFUP(JL,JK)
    PFLUC(JL,2,JK) = ZFDN(JL,JK)
  END DO
END DO
!
!     ------------------------------------------------------------------
!
RETURN
END SUBROUTINE OLWVB
