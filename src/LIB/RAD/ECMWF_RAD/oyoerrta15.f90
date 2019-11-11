!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA15


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA15* - RRTM COEFFICIENTS FOR INTERVAL 15
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG15 = 2

REAL_B :: FRACREFA(NG15,9)

REAL_B :: KA(9,5,13,NG15) ,ABSA(585,NG15)
REAL_B :: SELFREF(10,NG15)
REAL_B :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL 
! STRRAT  : REAL    
!     -----------------------------------------------------------------
END MODULE OYOERRTA15
