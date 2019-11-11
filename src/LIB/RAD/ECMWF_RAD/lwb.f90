!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:34
!-----------------------------------------------------------------
SUBROUTINE LWB &
 &( KIDIA, KFDIA, KLON  , KLEV  , KMODE &
 &, PDT0 , PTAVE, PTL &
 &, PB   , PBINT, PBSUR , PBTOP , PDBSL &
 &, PGA  , PGB  , PGASUR, PGBSUR, PGATOP, PGBTOP    &
 &)

!**** *LWB*   - COMPUTES BLACK-BODY FUNCTIONS FOR LONGWAVE CALCULATIONS

!     PURPOSE.
!     --------
!           COMPUTES PLANCK FUNCTIONS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PDT0   : (KLON)            ; SURFACE TEMPERATURE DISCONTINUITY
! PTAVE  : (KLON,KLEV)       ; TEMPERATURE
! PTL    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
!     ==== OUTPUTS ===
! PB     : (KLON,NSIL,KLEV+1); SPECTRAL HALF LEVEL PLANCK FUNCTION
! PBINT  : (KLON,KLEV+1)     ; HALF LEVEL PLANCK FUNCTION
! PBSUR  : (KLON,NSIL)       ; SURFACE SPECTRAL PLANCK FUNCTION
! PBTOP  : (KLON,NSIL)       ; TOP SPECTRAL PLANCK FUNCTION
! PDBSL  : (KLON,NSIL,KLEV*2); SUB-LAYER PLANCK FUNCTION GRADIENT
! PGA    : (KLON,8,2,KLEV)   ; dB/dT-weighted LAYER PADE APPROXIMANTS
! PGB    : (KLON,8,2,KLEV)   ; dB/dT-weighted LAYER PADE APPROXIMANTS
! PGASUR, PGBSUR (KLON,8,2)  ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP (KLON,8,2)  ; T.O.A. PADE APPROXIMANTS

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
!     FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS           "

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        MODIFIED : 99-06-14  D.SALMOND  Optimisation

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : MXIXT    ,NSIL     ,NIPD     ,PDGA     ,&
            &PDGB     ,TINTP    ,TSTAND   ,TSTP     ,XP


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KMODE



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PDT0(KLON),PTAVE(KLON,KLEV),PTL(KLON,KLEV+1)

REAL_B :: PB(KLON,NSIL,KLEV+1)   , PBINT(KLON,KLEV+1)&
  &,  PBSUR(KLON,NSIL)       ,  PBTOP(KLON,NSIL)   &
  &,  PDBSL(KLON,NSIL,KLEV*2)&
  &,  PGA(KLON,NIPD,2,KLEV)  , PGB(KLON,NIPD,2,KLEV)&
  &,  PGASUR(KLON,NIPD,2)    , PGBSUR(KLON,NIPD,2)&
  &,  PGATOP(KLON,NIPD,2)    , PGBTOP(KLON,NIPD,2)

!-------------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------
INTEGER_M :: INDB(KLON)   , INDS(KLON)
REAL_B :: ZBLAY(KLON,KLEV), ZBLEV(KLON,KLEV+1)&
  &,  ZRES(KLON)      , ZRES2(KLON)&
  &,  ZTI(KLON)       , ZTI2(KLON)

!     LOCAL INTEGER SCALARS
INTEGER_M :: ILEV2, INDSU, INDT, INDTO, INDTP, INUE, INUS,&
             &IXTOX, IXTX, JF, JG, JK, JK1, JK2, JL, JNU

!     LOCAL REAL SCALARS
REAL_B :: ZDST1, ZDSTO1, ZDSTOX, ZDSTX


!     ------------------------------------------------------------------


!*         1.0     PLANCK FUNCTIONS AND GRADIENTS
!                  ------------------------------

ILEV2=2*KLEV
INUS=1
INUE=NSIL
IF (KMODE == 2) THEN
  INUS=3
  INUE=4
ENDIF

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PBINT(JL,JK) = _ZERO_
  ENDDO
ENDDO

DO JNU=1,NSIL
  DO JL=KIDIA,KFDIA
    PBSUR(JL,JNU)=_ZERO_
    PBTOP(JL,JNU)=_ZERO_
  ENDDO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PB(JL,JNU,JK)=_ZERO_
    ENDDO
  ENDDO
  DO JK=1,ILEV2
    DO JL=KIDIA,KFDIA
      PDBSL(JL,JNU,JK)=_ZERO_
    ENDDO
  ENDDO
ENDDO

DO JNU=INUS,INUE


!*         1.1   LEVELS FROM SURFACE TO KLEV
!                ----------------------------

  DO JK = 1 , KLEV
    DO JL = KIDIA,KFDIA
      ZTI(JL)=(PTL(JL,JK)-TSTAND)/TSTAND
      ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU)&
       &+ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU)&
       &)))))
      PBINT(JL,JK)=PBINT(JL,JK)+ZRES(JL)
      PB(JL,JNU,JK)= ZRES(JL)
      ZBLEV(JL,JK) = ZRES(JL)

      ZTI2(JL)=(PTAVE(JL,JK)-TSTAND)/TSTAND
      ZRES2(JL)=XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU)&
       &+ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,&
       &JNU)&
       &)))))
      ZBLAY(JL,JK) = ZRES2(JL)
    ENDDO
  ENDDO


!*         1.2   TOP OF THE ATMOSPHERE AND SURFACE
!                ---------------------------------

  DO JL = KIDIA,KFDIA
    ZTI(JL)=(PTL(JL,KLEV+1)-TSTAND)/TSTAND
    ZTI2(JL) = (PTL(JL,1) + PDT0(JL) - TSTAND) / TSTAND
    ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU)&
     &+ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU)&
     &)))))
    ZRES2(JL) = XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU)&
     &+ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU)&
     &)))))
    PBINT(JL,KLEV+1) = PBINT(JL,KLEV+1)+ZRES(JL)
    PB(JL,JNU,KLEV+1)= ZRES(JL)
    ZBLEV(JL,KLEV+1) = ZRES(JL)
    PBTOP(JL,JNU) = ZRES(JL)
    PBSUR(JL,JNU) = ZRES2(JL)
  ENDDO


!*         1.3   GRADIENTS IN SUB-LAYERS
!                -----------------------

  DO JK = 1 , KLEV
    JK2 = 2 * JK
    JK1 = JK2 - 1
    DO JL = KIDIA,KFDIA
      PDBSL(JL,JNU,JK1) = ZBLAY(JL,JK  ) - ZBLEV(JL,JK)
      PDBSL(JL,JNU,JK2) = ZBLEV(JL,JK+1) - ZBLAY(JL,JK)
    ENDDO
  ENDDO

ENDDO

!*         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
!                ---------------------------------------------

DO JL=KIDIA,KFDIA
  ZDSTO1 = (PTL(JL,KLEV+1)-TINTP(1)) / TSTP
  IXTOX = MAX( 1, MIN( MXIXT, INT( ZDSTO1 + _ONE_ ) ) )
  ZDSTOX = (PTL(JL,KLEV+1)-TINTP(IXTOX))/TSTP
  IF (ZDSTOX < _HALF_) THEN
    INDTO=IXTOX
  ELSE
    INDTO=IXTOX+1
  ENDIF
  INDB(JL)=INDTO
  ZDST1 = (PTL(JL,1)-TINTP(1)) / TSTP
  IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + _ONE_ ) ) )
  ZDSTX = (PTL(JL,1)-TINTP(IXTX))/TSTP
  IF (ZDSTX < _HALF_) THEN
    INDT=IXTX
  ELSE
    INDT=IXTX+1
  ENDIF
  INDS(JL)=INDT
ENDDO

DO JF=1,2
  DO JG=1,NIPD
    DO JL=KIDIA,KFDIA
      INDSU=INDS(JL)
      PGASUR(JL,JG,JF)=PDGA(INDSU,2*JG-1,JF)
      PGBSUR(JL,JG,JF)=PDGB(INDSU,2*JG-1,JF)
      INDTP=INDB(JL)
      PGATOP(JL,JG,JF)=PDGA(INDTP,2*JG-1,JF)
      PGBTOP(JL,JG,JF)=PDGB(INDTP,2*JG-1,JF)
    ENDDO
  ENDDO
ENDDO


DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDST1 = (PTAVE(JL,JK)-TINTP(1)) / TSTP
    IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + _ONE_ ) ) )
    ZDSTX = (PTAVE(JL,JK)-TINTP(IXTX))/TSTP
    IF (ZDSTX < _HALF_) THEN
      INDT=IXTX
    ELSE
      INDT=IXTX+1
    ENDIF
    INDB(JL)=INDT
  ENDDO

  DO JF=1,2
    DO JL=KIDIA,KFDIA
      INDT=INDB(JL)
      DO JG=1,NIPD
        PGA(JL,JG,JF,JK)=PDGA(INDT,2*JG,JF)
        PGB(JL,JG,JF,JK)=PDGB(INDT,2*JG,JF)
      ENDDO
    ENDDO
  ENDDO


ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE LWB
