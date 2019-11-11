!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA7


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA7* - RRTM COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG7  = 12

REAL_B :: FRACREFA(NG7,9)

REAL_B , DIMENSION(NG7) :: FRACREFB
REAL_B , DIMENSION(NG7) :: ABSCO2
REAL_B :: KA(9,5,13,NG7) ,ABSA(585,NG7)
REAL_B :: KB(5,13:59,NG7),ABSB(235,NG7)
REAL_B :: SELFREF(10,NG7)

REAL_B :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL 
! ABSCO2  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL  
! STRRAT  : REAL   
!     -----------------------------------------------------------------
END MODULE OYOERRTA7
