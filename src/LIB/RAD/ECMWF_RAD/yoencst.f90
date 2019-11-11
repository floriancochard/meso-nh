!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOENCST


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*     OF PHYSICAL CONSTANTS
!     YOU WILL FIND THE MEANINGS IN THE ANNEX 1 OF THE DOCUMENTATION

! A1.0 FUNDAMENTAL CONSTANTS
REAL_B :: RPI
REAL_B :: RCLUM
REAL_B :: RHPLA
REAL_B :: RKBOL
REAL_B :: RNAVO
! A1.1 ASTRONOMICAL CONSTANTS
REAL_B :: RDAY
REAL_B :: REA
REAL_B :: REPSM
REAL_B :: RSIYEA
REAL_B :: RSIDAY
REAL_B :: ROMEGA
! A1.2 GEOIDE
REAL_B :: RA
REAL_B :: RG
REAL_B :: R1SA
! A1.3 RADIATION
REAL_B :: RSIGMA
REAL_B :: RI0
! A1.4 THERMODYNAMIC GAS PHASE
REAL_B :: R
REAL_B :: RMD
REAL_B :: RMV
REAL_B :: RD
REAL_B :: RV
REAL_B :: RCPD
REAL_B :: RCPV
REAL_B :: RCVD
REAL_B :: RCVV
REAL_B :: RKAPPA
! A1.5,6 THERMODYNAMIC LIQUID,SOLID PHASES
REAL_B :: RCW
REAL_B :: RCS
! A1.7 THERMODYNAMIC TRANSITION OF PHASE
REAL_B :: RLVTT
REAL_B :: RLSTT
REAL_B :: RTT
REAL_B :: RATM
! A1.8 CURVE OF SATURATION
REAL_B :: RESTT
REAL_B :: RALPW
REAL_B :: RBETW
REAL_B :: RGAMW
REAL_B :: RALPS
REAL_B :: RBETS
REAL_B :: RGAMS

!    ------------------------------------------------------------------
END MODULE YOENCST
