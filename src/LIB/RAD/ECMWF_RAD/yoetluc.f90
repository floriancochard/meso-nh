!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOETLUC


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    *YOETLUC* - Look-up tables for thermodynamic functions
!     ------------------------------------------------------------------
!     PARAMETER( JPTLUCU1=100000 , JPTLUCU2=370000 )

!     REAL TLUCUA(JPTLUCU1:JPTLUCU2),TLUCUB(JPTLUCU1:JPTLUCU2)
!    +  ,  TLUCUC(JPTLUCU1:JPTLUCU2)

!     /POETLUC/ MTLUCUA,MTLUCUB,MTLUCUC

!     POINTER (MTLUCUA,TLUCUA),(MTLUCUB,TLUCUB),(MTLUCUC,TLUCUC)


!          D.SALMOND  CRAY (UK)   12/8/91
!          J.J. MORCRETTE ECMWF   92-09-18  Update to Cy44

!     NAME      TYPE     DEFINITION
!     ----      ----     ----------

!     TLUCUA    REAL     TABLE 1
!     TLUCUB    REAL     TABLE 1
!     TLUCUC    REAL     TABLE 1

!     ------------------------------------------------------------------
END MODULE YOETLUC
