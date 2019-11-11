!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA10


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA14* - RRTM COEFFICIENTS FOR INTERVAL 10
!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG10 = 6

REAL_B , DIMENSION(NG10) :: FRACREFA
REAL_B , DIMENSION(NG10) :: FRACREFB

REAL_B :: KA(5,13,NG10)   , ABSA(65,NG10)
REAL_B :: KB(5,13:59,NG10), ABSB(235,NG10)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA10
