!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA5


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA5* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG5  = 16

REAL_B :: FRACREFA(NG5,9) ,FRACREFB(NG5,5)

REAL_B , DIMENSION(NG5) :: CCL4

REAL_B :: KA(9,5,13,NG5)   ,ABSA(585,NG5)
REAL_B :: KB(5,5,13:59,NG5),ABSB(1175,NG5)
REAL_B :: SELFREF(10,NG5)
REAL_B :: STRRAT1
REAL_B :: STRRAT2

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! CCL4    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
! STRRAT1 : REAL
! STRRAT2 : REAL    
!     -----------------------------------------------------------------
END MODULE OYOERRTA5
