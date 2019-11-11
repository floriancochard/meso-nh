!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARRTM1D

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     Parameters for 1-D radiation only computations from operational
!      library routines
!     991007  JJMorcrette
!     ------------------------------------------------------------------

INTEGER_M, PARAMETER :: JP_LON  = 1
INTEGER_M, PARAMETER :: JP_IDIA = 1
INTEGER_M, PARAMETER :: JP_FDIA = JP_LON-JP_IDIA+1
INTEGER_M, PARAMETER :: JP_TDIA = 1

!-- standard tropical INTEGER_M, PARAMETER :: JP_LEV  = 64
!-- ATEX              INTEGER_M, PARAMETER :: JP_LEV  = 83
!-- BOMEX             INTEGER_M, PARAMETER :: JP_LEV  = 84
!-- OPEN_CELLS        INTEGER_M, PARAMETER :: JP_LEV  = 63
!-- GATE_A,B,C        INTEGER_M, PARAMETER :: JP_LEV  = 46

INTEGER_M, PARAMETER :: JP_LEV = 100

INTEGER_M, PARAMETER :: JP_LW   = 6
INTEGER_M, PARAMETER :: JP_SW   = 6
INTEGER_M, PARAMETER :: JP_NUA  = 24
INTEGER_M, PARAMETER :: JP_MODE = 1
INTEGER_M, PARAMETER :: JP_AER  = 6
INTEGER_M, PARAMETER :: JP_LEVP1= JP_LEV+1

!     ------------------------------------------------------------------
END MODULE PARRTM1D
