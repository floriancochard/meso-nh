!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO2


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO2* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 2
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO2  = 16

!     The ith set of reference fractions are from the ith reference
!     pressure level.
REAL_B :: FRACREFAO(NO2,13), FRACREFBO(NO2)
REAL_B :: KAO(5,13,NO2)
REAL_B :: KBO(5,13:59,NO2)
REAL_B :: SELFREFO(10,NO2) , FORREFO(NO2)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
!FRACREFAO: REAL    
!FRACREFBO: REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL
! FORREFO : REAL  
!     -----------------------------------------------------------------
END MODULE OYOERRTO2
