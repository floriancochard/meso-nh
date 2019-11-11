!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOETHF


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!     ------------------------------------------------------------------

REAL_B :: R2ES
REAL_B :: R3LES
REAL_B :: R3IES
REAL_B :: R4LES
REAL_B :: R4IES
REAL_B :: R5LES
REAL_B :: R5IES
REAL_B :: RVTMP2
REAL_B :: RHOH2O
REAL_B :: R5ALVCP
REAL_B :: R5ALSCP
REAL_B :: RALVDCP
REAL_B :: RALSDCP
REAL_B :: RALFDCP
REAL_B :: RTWAT
REAL_B :: RTBER
REAL_B :: RTBERCU
REAL_B :: RTICE
REAL_B :: RTICECU

!     J.-J. MORCRETTE                   91/07/14  ADAPTED TO I.F.S.

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *R__ES*   REAL      *CONSTANTS USED FOR COMPUTATION OF SATURATION
!                         MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
!                         ICE(*R_IES*).
!     *RVTMP2*  REAL      *RVTMP2=RCPV/RCPD-1.
!     *RHOH2O*  REAL      *DENSITY OF LIQUID WATER.   (RATM/100.)
!     *R5ALVCP* REAL      *R5LES*RLVTT/RCPD
!     *R5ALSCP* REAL      *R5IES*RLSTT/RCPD
!     *RALVDCP* REAL      *RLVTT/RCPD
!     *RALSDCP* REAL      *RLSTT/RCPD
!     *RALFDCP* REAL      *RLMLT/RCPD
!     *RTWAT*   REAL      *RTWAT=RTT
!     *RTBER*   REAL      *RTBER=RTT-0.05
!     *RTBERCU  REAL      *RTBERCU=RTT-5.0
!     *RTICE*   REAL      *RTICE=RTT-0.1
!     *RTICECU* REAL      *RTICECU=RTT-23.0

!       ----------------------------------------------------------------
END MODULE YOETHF
