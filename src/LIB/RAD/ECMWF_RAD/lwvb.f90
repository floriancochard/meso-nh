!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE LWVB &
  &( KIDIA , KFDIA , KLON  , KLEV  , KUAER &
  &, PABCU , PADJD , PADJU &
  &, PB    , PBINT , PBSUR , PBTOP &
  &, PDISD , PDISU , PEMIS , PEMIW &
  &, PGASUR, PGBSUR, PGATOP, PGBTOP &
  &, PDWFSU,PFLUC                  &
  &)

!**** *LWVB*   - L.W., VERTICAL INTEGRATION, EXCHANGE WITH BOUNDARIES

!     PURPOSE.
!     --------
!           INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
!           INTEGRATION

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PADJ.. : (KLON,KLEV+1)     ; CONTRIBUTION BY ADJACENT LAYERS
! PB     : (KLON,NSIL,KLEV+1); SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
! PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
! PBSUR  : (KLON,NSIL)       ; SPECTRAL SURFACE PLANCK FUNCTION
! PBTOP  : (KLON,NSIL)       ; SPECTRAL T.O.A. PLANCK FUNCTION
! PDIS.. : (KLON,KLEV+1)     ; CONTRIBUTION BY DISTANT LAYERS
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PGASUR, PGBSUR             ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP             ; T.O.A. PADE APPROXIMANTS 
!     ==== OUTPUTS ===
! PDWFSU : (KLON,NSIL)       ; SPECTRAL DOWNWARD FLUX AT SURFACE
! PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR-SKY:
!                     1  ==>  UPWARD   FLUX TOTAL

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
!     ATMOSPHERE AND ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY 
!     FLUXES

!     EXTERNALS.
!     ----------

!          *LWTT*

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
!        JJ Morcrette 96-06-07  Surface LW Window Emissivity

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : NSIL     ,NIPD     ,NTRA     ,NUA      ,NG1P1
!
USE MODI_LWTT
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KUAER



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------


REAL_B :: PABCU(KLON,NUA,3*KLEV+1)&
  &,  PADJD(KLON,KLEV+1) , PADJU(KLON,KLEV+1)&
  &,  PB(KLON,NSIL,KLEV+1),PBINT(KLON,KLEV+1)&
  &,  PBSUR(KLON,NSIL)   , PBTOP(KLON,NSIL)&
  &,  PDISD(KLON,KLEV+1) , PDISU(KLON,KLEV+1)&
  &,  PEMIS(KLON)        , PEMIW(KLON)&
  &,  PGASUR(KLON,NIPD,2), PGBSUR(KLON,NIPD,2)&
  &,  PGATOP(KLON,NIPD,2), PGBTOP(KLON,NIPD,2)

REAL_B :: PDWFSU(KLON,NSIL)  , PFLUC(KLON,2,KLEV+1)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZBSUR(KLON,NSIL)&
  &,  ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA)&
  &,  ZUU(KLON,NUA) , ZCNSOL(KLON)   , ZCNTOP(KLON)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IN, JA, JK, JL

!     LOCAL REAL SCALARS
REAL_B :: ZCNTOP1, ZCNTOP2, ZCNTOP3, ZCNTOP4, ZCNTOP5, ZCNTOP6


!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------


!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------

DO JA=1,NTRA
  DO JL=KIDIA,KFDIA
    ZTT (JL,JA)=_ONE_
    ZTT1(JL,JA)=_ONE_
    ZTT2(JL,JA)=_ONE_
  ENDDO
ENDDO

DO JA=1,NUA
  DO JL=KIDIA,KFDIA
    ZUU(JL,JA)=_ONE_
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------


!*         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
!                  -----------------------------------

DO JK = 1 , KLEV
  IN=(JK-1)*NG1P1+1

  DO JA=1,KUAER
    DO JL=KIDIA,KFDIA
      ZUU(JL,JA)=PABCU(JL,JA,IN)
    ENDDO
  ENDDO


  CALL LWTT &
   &( KIDIA        , KFDIA        , KLON &
   &, PGATOP(1,1,1), PGBTOP(1,1,1)&
   &, ZUU          , ZTT &
   &)

  DO JL = KIDIA,KFDIA
    ZCNTOP1=PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
    ZCNTOP2=PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
    ZCNTOP3=PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
    ZCNTOP4=PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
    ZCNTOP5=PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
    ZCNTOP6=PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
    ZCNTOP(JL)=ZCNTOP1+ZCNTOP2+ZCNTOP3+ZCNTOP4+ZCNTOP5+ZCNTOP6
    PFLUC(JL,2,JK)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
    IF (JK == 1) THEN
      PDWFSU(JL,1)=ZCNTOP1-PB(JL,1,JK)-PDWFSU(JL,1)
      PDWFSU(JL,2)=ZCNTOP2-PB(JL,2,JK)-PDWFSU(JL,2)
      PDWFSU(JL,3)=ZCNTOP3-PB(JL,3,JK)-PDWFSU(JL,3)
      PDWFSU(JL,4)=ZCNTOP4-PB(JL,4,JK)-PDWFSU(JL,4)
      PDWFSU(JL,5)=ZCNTOP5-PB(JL,5,JK)-PDWFSU(JL,5)
      PDWFSU(JL,6)=ZCNTOP6-PB(JL,6,JK)-PDWFSU(JL,6)
    ENDIF
  ENDDO

ENDDO

JK = KLEV+1
IN=(JK-1)*NG1P1+1

DO JL = KIDIA,KFDIA
  ZCNTOP(JL)= PBTOP(JL,1)&
   &+ PBTOP(JL,2)&
   &+ PBTOP(JL,3)&
   &+ PBTOP(JL,4)&
   &+ PBTOP(JL,5)&
   &+ PBTOP(JL,6)
  PFLUC(JL,2,JK)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
ENDDO


!*         2.5     EXCHANGE WITH LOWER LIMIT
!                  -------------------------

JK = 1
IN=(JK-1)*NG1P1+1

DO JL = KIDIA,KFDIA
  ZBSUR(JL,1)=PBSUR(JL,1)*PEMIS(JL) -(_ONE_-PEMIS(JL))*PDWFSU(JL,1)
  ZBSUR(JL,2)=PBSUR(JL,2)*PEMIS(JL) -(_ONE_-PEMIS(JL))*PDWFSU(JL,2)
  ZBSUR(JL,3)=PBSUR(JL,3)*PEMIW(JL) -(_ONE_-PEMIW(JL))*PDWFSU(JL,3)
  ZBSUR(JL,4)=PBSUR(JL,4)*PEMIW(JL) -(_ONE_-PEMIW(JL))*PDWFSU(JL,4)
  ZBSUR(JL,5)=PBSUR(JL,5)*PEMIS(JL) -(_ONE_-PEMIS(JL))*PDWFSU(JL,5)
  ZBSUR(JL,6)=PBSUR(JL,6)*PEMIS(JL) -(_ONE_-PEMIS(JL))*PDWFSU(JL,6)

  PFLUC(JL,1,JK) = ZBSUR(JL,1)&
   &+ ZBSUR(JL,2)&
   &+ ZBSUR(JL,3)&
   &+ ZBSUR(JL,4)&
   &+ ZBSUR(JL,5)&
   &+ ZBSUR(JL,6)

  ZBSUR(JL,1)=ZBSUR(JL,1)-PB(JL,1,1)
  ZBSUR(JL,2)=ZBSUR(JL,2)-PB(JL,2,1)
  ZBSUR(JL,3)=ZBSUR(JL,3)-PB(JL,3,1)
  ZBSUR(JL,4)=ZBSUR(JL,4)-PB(JL,4,1)
  ZBSUR(JL,5)=ZBSUR(JL,5)-PB(JL,5,1)
  ZBSUR(JL,6)=ZBSUR(JL,6)-PB(JL,6,1)
ENDDO

DO JK = 2 , KLEV+1
  IN=(JK-1)*NG1P1+1

  DO JA=1,KUAER
    DO JL=KIDIA,KFDIA
      ZUU(JL,JA)=PABCU(JL,JA,1)-PABCU(JL,JA,IN)
    ENDDO
  ENDDO


  CALL LWTT &
   &( KIDIA        , KFDIA        , KLON &
   &, PGASUR(1,1,1), PGBSUR(1,1,1)&
   &, ZUU, ZTT &
   &)

  DO JL = KIDIA,KFDIA
    ZCNSOL(JL)=ZBSUR(JL,1)*ZTT(JL,1)          *ZTT(JL,10)&
     &+ZBSUR(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)&
     &+ZBSUR(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)&
     &+ZBSUR(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)&
     &+ZBSUR(JL,5)*ZTT(JL,3)          *ZTT(JL,14)&
     &+ZBSUR(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
    PFLUC(JL,1,JK)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
  ENDDO

ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE LWVB
