!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:34
!-----------------------------------------------------------------
SUBROUTINE LW &
 &( KIDIA, KFDIA , KLON  , KLEV  , KMODE &
 &, PCCO2, PCLDLD, PCLDLU &
 &, PDP  , PDT0  , PEMIS , PEMIW &
 &, PPMB , PQOF  , PTL &
 &, PAER , PTAVE , PVIEW , PWV &
 &, PCOLR, PCOLC , PEMIT , PFLUX , PFLUC &
 &)

!**** *LW*   - ORGANIZES THE LONGWAVE CALCULATIONS

!     PURPOSE.
!     --------
!           COMPUTES LONGWAVE FLUXES 

!**   INTERFACE.
!     ----------

!        *LW* IS CALLED FROM *RADLSW*

!        EXPLICIT ARGUMENTS :
!        --------------------
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PCCO2  :                   ; CONCENTRATION IN CO2 (PA/PA)
! PCLDLD : (KLON,KLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
! PCLDLU : (KLON,KLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PDT0   : (KLON)            ; SURFACE TEMPERATURE DISCONTINUITY  
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PPMB   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
! PTAVE  : (KLON,KLEV)       ; TEMPERATURE
! PTL    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
! PVIEW  : (KLON)            ; COSECANT OF VIEWING ANGLE
! PWV    : (KLON,KLEV)       ; SPECIFIC HUMIDITY  (PA/PA)
!     ==== OUTPUTS ===
! PCOLR(KLON,KLEV)           ; LONG-WAVE TENDENCY
! PCOLC(KLON,KLEV)           ; LONG-WAVE TENDENCY CLEAR SKY
! PEMIT(KLON)                ; SURFACE TOTAL LW EMISSIVITY
! PFLUX(KLON,2,KLEV)         ; RADIATIVE FLUXES :
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL
! PFLUC(KLON,2,KLEV)         ; RADIATIVE FLUXES CLEAR SKY:
!                     1  ==>  UPWARD   FLUX TOTAL
!                     2  ==>  DOWNWARD FLUX TOTAL

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
!     ABSORBERS.
!          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
!     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
!          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
!     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
!     BOUNDARIES.
!          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
!          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.

!     EXTERNALS.
!     ----------

!          *LWU*, *LWBV*, *LWC*

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
!        99-05-25   JJMorcrette    Revised aerosols

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : NUA
USE YOERDU   , ONLY : RCDAY
!
USE MODI_LWU
USE MODI_LWBV
USE MODI_LWC

!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KMODE

!     DUMMY REAL SCALARS
REAL_B :: PCCO2



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAER(KLON,6,KLEV)&
  &,  PCLDLD(KLON,KLEV)   , PCLDLU(KLON,KLEV)&
  &,  PDP(KLON,KLEV)      , PDT0(KLON)       &
  &,  PEMIS(KLON)         , PEMIW(KLON)&
  &,  PPMB(KLON,KLEV+1)&
  &,  PQOF(KLON,KLEV)  &
  &,  PTL(KLON,KLEV+1)    , PTAVE(KLON,KLEV) &
  &,  PVIEW(KLON)         , PWV(KLON,KLEV)

REAL_B :: PCOLR(KLON,KLEV)    , PCOLC(KLON,KLEV)&
  &,  PEMIT(KLON) &
  &,  PFLUX(KLON,2,KLEV+1), PFLUC(KLON,2,KLEV+1)

!-------------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------
REAL_B :: ZABCU(KLON,NUA,3*KLEV+1)&
  &,  ZBINT(KLON,KLEV+1)        , ZBSUI(KLON)&
  &,  ZCNTRB(KLON,KLEV+1,KLEV+1)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JK, JKL, JL

!     LOCAL REAL SCALARS
REAL_B :: ZDCNET, ZDFNET


!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!100  CONTINUE

!     ------------------------------------------------------------------

!*         1.1   COMPUTES ABSORBER AMOUNTS
!                -------------------------

CALL LWU &
  &(  KIDIA, KFDIA, KLON, KLEV &
  &,  PAER , PCCO2, PDP , PPMB, PQOF , PTAVE, PVIEW, PWV &
  &,  ZABCU &
  &)

!     ------------------------------------------------------------------

!*         2.    COMPUTES PLANCK FUNCTIONS
!                -------------------------
!                PERFORMS THE VERTICAL INTEGRATION
!                ---------------------------------

CALL LWBV &
   &( KIDIA, KFDIA, KLON , KLEV  , KMODE &
   &, PDT0 , PEMIS, PEMIW, PTL   , PTAVE &
   &, PEMIT, PFLUC &
   &, ZABCU, ZBINT, ZBSUI, ZCNTRB &
   &)

!     ------------------------------------------------------------------

!*         4.    INTRODUCES THE EFFECTS OF CLOUDS
!                --------------------------------

!print *,'Just before LWC'
CALL LWC &
  &( KIDIA , KFDIA, KLON  , KLEV &
  &, ZBINT , ZBSUI, PCLDLD, PCLDLU &
  &, ZCNTRB, PEMIT, PFLUC &
  &, PFLUX    &
  &)

DO JKL = 1 , KLEV
  JK = KLEV+1 - JKL
  DO JL = KIDIA,KFDIA
    ZDCNET = PFLUC(JL,1,JK+1) + PFLUC(JL,2,JK+1)&
     &-PFLUC(JL,1,JK  ) - PFLUC(JL,2,JK  )
    PCOLC(JL,JK) = RCDAY * ZDCNET / PDP(JL,JKL)
    ZDFNET = PFLUX(JL,1,JK+1) + PFLUX(JL,2,JK+1)&
     &-PFLUX(JL,1,JK  ) - PFLUX(JL,2,JK  )
    PCOLR(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE LW
