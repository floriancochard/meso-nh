!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARCLI


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!------------------------------------------------------------------------------
!     JPNAV : Number of latitude of the NAVY grid.
!     JPPOCP: Number of latitude of the NAVY grid in polar caps.
!     JPAUT : Number of latitude of the N108 grids.
!     JPNMOI: Number of months to be processed.
!     JPNFIN: Number of fixed fields in the *NAVY* data
!     JPNFIS: Number of fixed fields in the *N108* data
!     JPNMOS: Number of monthly fields in the *N108* data
!     JPNFIF: Number of fixed fields in the final data
!     JPNMOF: Number of monthly fields in the final data
!     JPNFCS: Number of fixed fields in the *N108* data which are indifferent
!             to the land/sea distribution
!     JPNFCF: Number of fixed fields in the final data which are indifferent
!             to the land/sea distribution
!     PPTSHT: Shift (in radians) of the deep yearly temperature wave.
!     PPWSHT: Shift (in radians) of the deep yearly water content wave.
!     JPLDAT: Dimension of the array *IDATEF*.

INTEGER_M, PARAMETER :: JPNAV=1080
INTEGER_M, PARAMETER :: JPPOCP=120
INTEGER_M, PARAMETER :: JPAUT=216
INTEGER_M, PARAMETER :: JPINT=2*JPAUT+4
INTEGER_M, PARAMETER :: JPJNT=JPAUT+4
INTEGER_M, PARAMETER :: JPNFIN=7
INTEGER_M, PARAMETER :: JPNFIS=13
INTEGER_M, PARAMETER :: JPNMOS=5
INTEGER_M, PARAMETER :: JPNFIF=11
INTEGER_M, PARAMETER :: JPNMOF=7
INTEGER_M, PARAMETER :: JPNMOI=12
INTEGER_M, PARAMETER :: JPNFCS=7
INTEGER_M, PARAMETER :: JPNFCF=5
REAL_B   , PARAMETER :: PPTSHT=0.0786_JPRB
REAL_B   , PARAMETER :: PPWSHT=0.1969_JPRB
INTEGER_M, PARAMETER :: JPLDAT=11

END MODULE PARCLI
