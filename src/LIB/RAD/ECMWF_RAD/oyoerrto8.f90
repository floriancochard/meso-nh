!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO8


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO8* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 8
!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO8  = 16

REAL_B , DIMENSION(NO8) :: FRACREFAO
REAL_B , DIMENSION(NO8) :: FRACREFBO
REAL_B , DIMENSION(NO8) :: CFC12O
REAL_B , DIMENSION(NO8) :: CFC22ADJO
REAL_B , DIMENSION(NO8) :: ABSCO2AO
REAL_B , DIMENSION(NO8) :: ABSCO2BO
REAL_B , DIMENSION(NO8) :: ABSN2OAO
REAL_B , DIMENSION(NO8) :: ABSN2OBO

REAL_B :: KAO(5,7,NO8)
REAL_B :: KBO(5,7:59,NO8)
REAL_B :: SELFREFO(10,NO8)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2A : REAL     
! ABSCO2B : REAL     
! ABSN2OA : REAL     
! ABSN2OB : REAL 
! CFC12   : REAL     
! CFC22ADJ: REAL     
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO8
