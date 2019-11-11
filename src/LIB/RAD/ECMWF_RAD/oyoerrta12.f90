!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA12


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA12* - RRTM COEFFICIENTS FOR INTERVAL 12
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG12 = 8

REAL_B :: FRACREFA(NG12,9)
REAL_B :: KA(9,5,13,NG12) ,ABSA(585,NG12)
REAL_B :: SELFREF(10,NG12)

REAL_B :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL
! STRRAT1 : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA12
