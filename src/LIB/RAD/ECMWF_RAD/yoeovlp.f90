!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOEOVLP

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOEOVLP* VERTICAL DISTRIBUTION OF CLOUD OVERLAP PARAMETER
!     ------------------------------------------------------------------

!REAL_B,ALLOCATABLE:: RA1OVLP(:)
REAL_B :: RA1OVLP(61)

!     J.-J. MORCRETTE    E.C.M.W.F.     01/02/16

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *RA1OVLP* REAL      Alpha1 (Hogan, Illingworth, 2001)

!     ------------------------------------------------------------------
END MODULE YOEOVLP
