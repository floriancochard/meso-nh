!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOEOZOC


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOEOZOC* SPECTRAL DISTRIBUTION OF OZONE
!     ------------------------------------------------------------------


REAL_B :: COZQC(21)
REAL_B :: COZQS(15)
REAL_B :: COZHC(21)
REAL_B :: COZHS(15)

REAL_B :: RSINC(19)
REAL_B :: ROZT(19,0:35)
REAL_B :: RPROC(0:35)

!*     *YOEOZOC* SPECTRAL DISTRIBUTION OF OZONE
!                     (TRIANGULAR *T5* TRUNCATION FOR OZONE).

!     R.G AND M.J        E.C.M.W.F.     29/11/82.

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *COZ__*   REAL      *REFERS TO *OZONE.
!     *C___C*   REAL      *REFERS TO *COS COMPONENT.
!     *C___S*   REAL      *REFERS TO *SIN COMPONENT.
!       ----------------------------------------------------------------
END MODULE YOEOZOC
