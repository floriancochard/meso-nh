!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARRINT


#include "tsmbkind.h"

IMPLICIT NONE

SAVE


!     PARAMETERS REQUIRED FOR RADIATION INTERPOLATION

INTEGER_M, PARAMETER :: JPRADCW=2
INTEGER_M, PARAMETER :: JPRADCE=2
INTEGER_M, PARAMETER :: JPRADCT=JPRADCW+JPRADCE
INTEGER_M, PARAMETER :: JPRADFW=1
INTEGER_M, PARAMETER :: JPRADFE=1
END MODULE PARRINT
