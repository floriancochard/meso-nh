!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTO9


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO9* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 9
!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO9  = 16

REAL_B :: FRACREFAO(NO9,9)

REAL_B , DIMENSION(NO9) :: FRACREFBO
! 48 = 3*NO9      
REAL_B , DIMENSION(48) :: ABSN2OO

REAL_B :: KAO(11,5,13,NO9)
REAL_B :: KBO(5,13:59,NO9)
REAL_B :: SELFREFO(10,NO9)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSN2O  : REAL    
! CH4REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! N2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE OYOERRTO9
