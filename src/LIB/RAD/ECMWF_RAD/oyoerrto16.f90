!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO16


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO16* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO16 = 16

REAL_B :: FRACREFAO(NO16,9)

REAL_B :: KAO(9,5,13,NO16)
REAL_B :: SELFREFO(10,NO16)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO16
