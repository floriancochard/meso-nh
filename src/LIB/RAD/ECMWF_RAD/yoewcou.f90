!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOEWCOU


#include "tsmbkind.h"

IMPLICIT NONE

SAVE


!*    ** *YOEWCOU* - VARIABLES FOR COUPLING WITH THE WAVE MODEL

!     P. VITERBO     E.C.M.W.F.       07/10/88
!     P. VITERBO     E.C.M.W.F.       03/02/92
!     J. DOYLE       E.C.M.W.F.       21/11/96 
!     J. BIDLOT      E.C.M.W.F.       13/06/97 
!     J. BIDLOT      E.C.M.W.F.       11/08/98 

!      NAME      TYPE       PURPOSE
!      ----      ----       -------

!     *NLONW*    INTEGER    NUMBER OF POINTS IN A LATITUDE LINE IN
!                           THE WAVE MODEL
!     *NLATW*    INTEGER    NUMBER OF LATITUDES IN THE WAVE MODEL.
!     *NLON1W*   INTEGER    *NLONW*
!     *NLAT1W*   INTEGER    TOTAL NUMBER OF LATITUDES WITH THE WAVE
!                           MODEL RESOLUTION.
!     *NNORXW*   INTEGER    NUMBER OF EXTRA POINTS NORTHWARDS OF THE
!                           NORTHERN BOUNDARY OF THE WAVE MODEL.
!     *CBEGDAT*  CHARACTER  INITIAL DATE OF FORECAST (YYYYMMDDHHmm)
!     *NDURAT*   INTEGER    DURATION (IN MINUTES) OF TOTAL WAVE INTEGRATION
!                              (OR FORECAST TIME IN MINUTES).
!     *NSTPW*    INTEGER    FREQUENCY OF CALL TO THE WAVE MODEL.
!     *NRESUM*   INTEGER    TIME STEP OF A RESTART EVENT
!     *LWCOU*    LOGICAL    TRUE IF THE WAVE MODEL IS TO BE RUN.
!     *LWCOU2W*  LOGICAL    TRUE IF TWO-WAY INTERACTION WITH THE WAVE MODEL.
!                           FALSE IF ONE-WAY INTERACTION WITH THE WAVE MODEL
!     *RSOUTW*   REAL       SOUTH BOUNDARY OF THE WAVE MODEL.
!     *RNORTW*   REAL       NORTH BOUNDARY OF THE WAVE MODEL.
!     *RDEGREW*  REAL       RESOLUTION OF THE WAVE MODEL (DEGREES).


INTEGER_M :: NLONW
INTEGER_M :: NLATW
INTEGER_M :: NLON1W
INTEGER_M :: NLAT1W
INTEGER_M :: NNORXW
INTEGER_M :: NDURAT
INTEGER_M :: NSTPW
INTEGER_M :: NRESUM


REAL_B :: RSOUTW
REAL_B :: RNORTW
REAL_B :: RDEGREW

LOGICAL LWCOU
LOGICAL LWCOU2W
CHARACTER CBEGDAT*12

!     ------------------------------------------------------------------
END MODULE YOEWCOU
