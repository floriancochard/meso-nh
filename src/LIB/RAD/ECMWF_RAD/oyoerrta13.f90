!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA13


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA13* - RRTM COEFFICIENTS FOR INTERVAL 13
!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG13 = 4

REAL_B :: FRACREFA(NG13,9)

REAL_B :: KA(9,5,13,NG13) ,ABSA(585,NG13)
REAL_B :: SELFREF(10,NG13)
REAL_B :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
! STRRAT1 : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA13
