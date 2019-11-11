!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA2


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA2* - RRTM COEFFICIENTS FOR INTERVAL 2
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG2  = 14

!     The ith set of reference fractions are from the ith reference
!     pressure level.
REAL_B :: FRACREFA(NG2,13), FRACREFB(NG2), REFPARAM(13)
REAL_B :: KA(5,13,NG2)   , ABSA(65,NG2)
REAL_B :: KB(5,13:59,NG2), ABSB(235,NG2)
REAL_B :: SELFREF(10,NG2), FORREF(NG2)

EQUIVALENCE (KA(1,1,1),ABSA(1,1)),(KB(1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! REFPARAM: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL
! FORREF  : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTA2


