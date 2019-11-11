!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO5


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO5* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO5  = 16

REAL_B :: FRACREFAO(NO5,9) ,FRACREFBO(NO5,5)

REAL_B , DIMENSION(NO5) :: CCL4O

REAL_B :: KAO(9,5,13,NO5)
REAL_B :: KBO(5,5,13:59,NO5)
REAL_B :: SELFREFO(10,NO5)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO5
