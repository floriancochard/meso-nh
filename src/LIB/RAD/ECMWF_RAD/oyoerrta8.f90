!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA8


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA8* - RRTM COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG8  = 8

REAL_B , DIMENSION(NG8) :: FRACREFA
REAL_B , DIMENSION(NG8) :: FRACREFB
REAL_B , DIMENSION(NG8) :: CFC12
REAL_B , DIMENSION(NG8) :: CFC22ADJ
REAL_B , DIMENSION(NG8) :: ABSCO2A
REAL_B , DIMENSION(NG8) :: ABSCO2B
REAL_B , DIMENSION(NG8) :: ABSN2OA
REAL_B , DIMENSION(NG8) :: ABSN2OB
REAL_B , DIMENSION(59)  :: H2OREF
REAL_B , DIMENSION(59)  :: N2OREF
REAL_B , DIMENSION(59)  :: O3REF

REAL_B :: KA(5,7,NG8)    ,ABSA(35,NG8)
REAL_B :: KB(5,7:59,NG8) ,ABSB(265,NG8)
REAL_B :: SELFREF(10,NG8)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,7,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! H2OREF  : REAL    
! KA      : REAL     
! KB      : REAL     
! N2OREF  : REAL    
! O3REF   : REAL    
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA8
