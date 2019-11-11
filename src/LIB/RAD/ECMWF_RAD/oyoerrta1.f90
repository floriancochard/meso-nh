!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA1


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA1* - RRTM COEFFICIENTS FOR INTERVAL 1
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG1  = 8

REAL_B :: FRACREFA(NG1)  , FRACREFB(NG1)
REAL_B :: KA(5,13,NG1)   , ABSA(65,NG1)
REAL_B :: KB(5,13:59,NG1), ABSB(235,NG1)
REAL_B :: SELFREF(10,NG1), FORREF(NG1)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! FORREF  : REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA1
