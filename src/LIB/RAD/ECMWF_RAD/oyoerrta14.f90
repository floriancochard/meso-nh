!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA14


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA14* - RRTM COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG14 = 2

REAL_B , DIMENSION(NG14) :: FRACREFA
REAL_B , DIMENSION(NG14) :: FRACREFB

REAL_B :: KA(5,13,NG14)   ,ABSA(65,NG14)
REAL_B :: KB(5,13:59,NG14),ABSB(235,NG14)
REAL_B :: SELFREF(10,NG14)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA14


