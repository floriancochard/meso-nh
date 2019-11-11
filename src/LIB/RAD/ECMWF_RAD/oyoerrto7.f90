!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO7


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO7* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 7
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO7  = 16

REAL_B :: FRACREFAO(NO7,9)

REAL_B , DIMENSION(NO7) :: FRACREFBO
REAL_B , DIMENSION(NO7) :: ABSCO2O
REAL_B :: KAO(9,5,13,NO7)
REAL_B :: KBO(5,13:59,NO7)
REAL_B :: SELFREFO(10,NO7)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL    
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL  
!     -----------------------------------------------------------------
END MODULE OYOERRTO7
