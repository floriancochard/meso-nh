!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARDIM


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------


!     JPMXLE : MAXIMUM NUMBER OF LEVELS
!     JPMXGL : MAXIMUM NUMBER OF GAUSSIAN LATITUDES 

!     JPCSS : MAXIMUM NUMBER OF LEVELS IN SOIL         

!     JPFPPYX : Maximum number of fields in catalogue cilipdy

!     JPNULNAM : Fortran unit number used for namelist reading

INTEGER_M, PARAMETER :: JPMXLE=200
INTEGER_M, PARAMETER :: JPMXGL=1001
INTEGER_M, PARAMETER :: JPCSS=5
INTEGER_M, PARAMETER :: JPFPPYX=17
INTEGER_M, PARAMETER :: JPNULNAM=4

!     ------------------------------------------------------------------
END MODULE PARDIM
