!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO1


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO1* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 1
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO1  = 16

REAL_B :: FRACREFAO(NO1)  , FRACREFBO(NO1)
REAL_B :: KAO(5,13,NO1)
REAL_B :: KBO(5,13:59,NO1)
REAL_B :: SELFREFO(10,NO1), FORREFO(NO1)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
!FRACREFAO: REAL    
!FRACREFBO: REAL
! FORREFO : REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO1
