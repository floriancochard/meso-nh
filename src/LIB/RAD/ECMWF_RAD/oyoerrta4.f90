!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA4


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA4* - RRTM COEFFICIENTS FOR INTERVAL 5
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG4  = 14

REAL_B :: FRACREFA(NG4,9)  ,FRACREFB(NG4,6)
REAL_B :: KA(9,5,13,NG4)   ,ABSA(585,NG4)
REAL_B :: KB(6,5,13:59,NG4),ABSB(1410,NG4)
REAL_B :: SELFREF(10,NG4)
REAL_B :: STRRAT1
REAL_B :: STRRAT2

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA4
