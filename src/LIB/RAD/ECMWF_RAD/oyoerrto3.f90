!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO3


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO3* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 3
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO3  = 16

REAL_B :: FRACREFAO(NO3,10) ,FRACREFBO(NO3,5)

REAL_B , DIMENSION(NO3) :: FORREFO
REAL_B , DIMENSION(NO3) :: ABSN2OAO
REAL_B , DIMENSION(NO3) :: ABSN2OBO

REAL_B :: KAO(10,5,13,NO3)
REAL_B :: KBO(5,5,13:59,NO3)
REAL_B :: SELFREFO(10,NO3)

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSN2OAO: REAL
! ABSN2OBO: REAL
!FRACREFAO: REAL    
!FRACREFBO: REAL
! KAO     : REAL     
! KBO     : REAL     
! SELFREFO: REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO3
