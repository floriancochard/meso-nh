!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE LWV &
 &( KIDIA, KFDIA, KLON , KLEV , KUAER , KTRAER &
 &, PABCU, PB   , PBINT, PBSUR, PBTOP , PDBSL &
 &, PEMIS, PEMIW &
 &, PGA  , PGB  , PGASUR,PGBSUR,PGATOP, PGBTOP &
 &, PCNTRB,PFLUC &
 &)

!**** *LWV*   - LONGWAVE RADIATION, VERTICAL INTEGRATION

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
!           FLUXES OR RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PB     : (KLON,NSIL,KLEV+1); SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
! PBINT  : (KLON,KLEV+1)     ; HALF-LEVEL PLANCK FUNCTIONS
! PBSUR  : (KLON,NSIL)       ; SURFACE SPECTRAL PLANCK FUNCTION
! PBTOP  : (KLON,NSIL)       ; T.O.A. SPECTRAL PLANCK FUNCTION
! PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PGA, PGB                   ; PADE APPROXIMANTS
! PGASUR, PGBSUR             ; SURFACE PADE APPROXIMANTS
! PGATOP, PGBTOP             ; T.O.A. PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PFLUC(KLON,2,KLEV)           ; RADIATIVE FLUXES CLEAR-SKY

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
!     CONTRIBUTIONS BY -  THE NEARBY LAYERS
!                      -  THE DISTANT LAYERS
!                      -  THE BOUNDARY TERMS
!          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.

!     EXTERNALS.
!     ----------

!          *LWVN*, *LWVD*, *LWVB*

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
!        JJ Morcrette 96-06-07 Surface LW window emissivity
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : NSIL     ,NIPD     ,NUA
!
USE MODI_LWVN
USE MODI_LWVD
USE MODI_LWVB
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KTRAER
INTEGER_M :: KUAER



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------


REAL_B :: PABCU(KLON,NUA,3*KLEV+1)&
  &,  PB(KLON,NSIL,KLEV+1)   , PBINT(KLON,KLEV+1)&
  &,  PBSUR(KLON,NSIL)       , PBTOP(KLON,NSIL)&
  &,  PDBSL(KLON,NSIL,KLEV*2)&
  &,  PEMIS(KLON)            , PEMIW(KLON)&
  &,  PGA(KLON,NIPD,2,KLEV)  , PGB(KLON,NIPD,2,KLEV)&
  &,  PGASUR(KLON,NIPD,2)    , PGBSUR(KLON,NIPD,2)&
  &,  PGATOP(KLON,NIPD,2)    , PGBTOP(KLON,NIPD,2)

REAL_B :: PCNTRB(KLON,KLEV+1,KLEV+1),  PFLUC(KLON,2,KLEV+1)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZADJD(KLON,KLEV+1)  , ZADJU(KLON,KLEV+1)&
  &,  ZDBDT(KLON,NSIL,KLEV)&
  &,  ZDISD(KLON,KLEV+1)  , ZDISU(KLON,KLEV+1)&
  &,  ZDWFSU(KLON,NSIL)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JA, JK, JL


!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    ZADJD(JL,JK)=_ZERO_
    ZADJU(JL,JK)=_ZERO_
    ZDISD(JL,JK)=_ZERO_
    ZDISU(JL,JK)=_ZERO_
  ENDDO
ENDDO
DO JA=1,NSIL
  DO JL=KIDIA,KFDIA
    ZDWFSU(JL,JA)=_ZERO_
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------

!     ------------------------------------------------------------------

!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------

CALL LWVN &
  &( KIDIA, KFDIA, KLON  , KLEV , KUAER &
  &, PABCU, PDBSL, PGA   , PGB &
  &, ZADJD, ZADJU, PCNTRB, ZDBDT, ZDWFSU  &
  &)

!     ------------------------------------------------------------------

!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------

CALL LWVD &
  &( KIDIA , KFDIA, KLON , KLEV  , KTRAER &
  &, PABCU , ZDBDT, PGA  , PGB &
  &, PCNTRB, ZDISD, ZDISU, ZDWFSU &
  &)

!     ------------------------------------------------------------------

!*         2.3     EXCHANGE WITH THE BOUNDARIES
!                  ----------------------------

CALL LWVB &
  &( KIDIA , KFDIA , KLON  , KLEV  , KUAER &
  &, PABCU , ZADJD , ZADJU &
  &, PB    , PBINT , PBSUR , PBTOP &
  &, ZDISD , ZDISU , PEMIS , PEMIW &
  &, PGASUR, PGBSUR, PGATOP, PGBTOP &
  &, ZDWFSU,PFLUC  &
  &)

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE LWV
