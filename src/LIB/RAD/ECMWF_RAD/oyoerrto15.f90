!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO15


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO15* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 15
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO15 = 16

REAL_B :: FRACREFAO(NO15,9)

REAL_B :: KAO(9,5,13,NO15)
REAL_B :: SELFREFO(10,NO15)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL 
!     -----------------------------------------------------------------
END MODULE OYOERRTO15
