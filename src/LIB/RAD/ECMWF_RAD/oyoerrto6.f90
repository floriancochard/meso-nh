!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO6


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO6* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO6  = 16

REAL_B , DIMENSION(NO6) :: FRACREFAO

REAL_B , DIMENSION(NO6) :: CFC11ADJO
REAL_B , DIMENSION(NO6) :: CFC12O
REAL_B , DIMENSION(NO6) :: ABSCO2O

REAL_B :: KAO(5,13,NO6)
REAL_B :: SELFREFO(10,NO6)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO6
