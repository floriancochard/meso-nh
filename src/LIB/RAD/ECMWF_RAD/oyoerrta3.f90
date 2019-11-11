!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE OYOERRTA3


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA3* - RRTM COEFFICIENTS FOR INTERVAL 3
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG3  = 16

REAL_B :: FRACREFA(NG3,10) ,FRACREFB(NG3,5)

REAL_B , DIMENSION(16) :: FORREF
REAL_B , DIMENSION(16) :: ABSN2OA
REAL_B , DIMENSION(16) :: ABSN2OB
REAL_B , DIMENSION(10) :: ETAREF
REAL_B , DIMENSION(59) :: H2OREF
REAL_B , DIMENSION(59) :: N2OREF
REAL_B , DIMENSION(59) :: CO2REF

REAL_B :: KA(10,5,13,NG3)  ,ABSA(650,NG3)
REAL_B :: KB(5,5,13:59,NG3),ABSB(1175,NG3)
REAL_B :: SELFREF(10,NG3)
REAL_B :: STRRAT

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)),(KB(1,1,13,1),ABSB(1,1))

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL
! ABSB    : REAL
! ABSN2OA : REAL
! ABSN2OB : REAL
! CO2REF  : REAL
! ETAREF  : REAL
! FRACREFA: REAL    
! FRACREFB: REAL
! H2OREF  : REAL
! KA      : REAL     
! KB      : REAL     
! N2OREF  : REAL
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
END MODULE OYOERRTA3
