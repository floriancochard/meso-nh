!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE YOEPHLI


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEPHLI* CONTAINS CONSTANTS FOR THE LINEARIZED PHYSICS
!     ------------------------------------------------------------------

LOGICAL LPHYLIN

REAL_B :: RLPTRC
REAL_B :: RLPAL1
REAL_B :: RLPAL2
REAL_B :: RLPBB
REAL_B :: RLPCC
REAL_B :: RLPDD
REAL_B :: RLPMIXL
REAL_B :: RLPBETA
REAL_B :: RLPDRAG
REAL_B :: RLPEVAP
REAL_B :: RLPP00

!*     *YOEPHLI* CONTAINS CONSTANTS NEEDED BY 
!     THE LINEARIZED PHYSICS


!     J.F. MAHFOUF        E.C.M.W.F.    23/06/96


!     NAME        TYPE     DESCRIPTION
!     ----        ----     -----------

!     *RLPTRC*    REAL     CRITICAL TEMPERATURE FOR MIXED PHASE PROPERTIES
!                          OF WATER 
!     *RLPAL1*    REAL     SMOOTHING COEFFICIENT
!     *RLPAL2*    REAL     SMOOTHING COEFFICIENT
!     *RLPBB*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
!     *RLPCC*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
!     *RLPDD*     REAL     CONSTANT FROM THE LOUIS ET AL. FORMULATION
!     *RLPMIXL*   REAL     PSEUDO DEPTH OF THE PLANETARY BOUNDARY LAYER
!     *RLPBETA*   REAL     REDUCTION FACTOR OF THE ASYMPTOTIC MIXING LENGTH
!     *RLPDRAG*   REAL     COEFFICIENT FOR THE ESTIMATION OF SURFACE DRAG
!     *RLPEVAP*   REAL     FRACTION OF POSSIBLE RAINFALL EVAPORATION
!     *RLPP00*    REAL     PRESSURE ABOVE WHICH RADIATION IS NOT APPLIED
!     *LPHYLIN*   LOGICAL  TRUE WHEN LINEARIZED PHYSICS IS ACTIVATED 

!     ------------------------------------------------------------------
END MODULE YOEPHLI
