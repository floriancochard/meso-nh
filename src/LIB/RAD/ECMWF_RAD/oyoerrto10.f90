!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO10


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 10
!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO10 = 16

REAL_B , DIMENSION(NO10) :: FRACREFAO
REAL_B , DIMENSION(NO10) :: FRACREFBO

REAL_B :: KAO(5,13,NO10)
REAL_B :: KBO(5,13:59,NO10)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO10
