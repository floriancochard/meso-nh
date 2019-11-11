!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA6


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA6* - RRTM COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG6  = 8

REAL_B , DIMENSION(NG6) :: FRACREFA

REAL_B , DIMENSION(NG6) :: CFC11ADJ
REAL_B , DIMENSION(NG6) :: CFC12
REAL_B , DIMENSION(NG6) :: ABSCO2

REAL_B :: KA(5,13,NG6),ABSA(65,NG6)
REAL_B :: SELFREF(10,NG6)

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2  : REAL 
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! CFC11ADJ: REAL
! CFC12   : REAL
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA6
