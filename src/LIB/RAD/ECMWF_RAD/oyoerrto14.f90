!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO14


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO14* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 14
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO14 = 16

REAL_B , DIMENSION(NO14) :: FRACREFAO
REAL_B , DIMENSION(NO14) :: FRACREFBO

REAL_B :: KAO(5,13,NO14)
REAL_B :: KBO(5,13:59,NO14)
REAL_B :: SELFREFO(10,NO14)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO14
