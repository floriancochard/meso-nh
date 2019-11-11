!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO12


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO12* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 12
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO12 = 16

REAL_B :: FRACREFAO(NO12,9)
REAL_B :: KAO(9,5,13,NO12)
REAL_B :: SELFREFO(10,NO12)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
!     -----------------------------------------------------------------
END MODULE OYOERRTO12
