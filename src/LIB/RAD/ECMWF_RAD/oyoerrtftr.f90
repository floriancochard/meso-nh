!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:49
!-----------------------------------------------------------------
MODULE OYOERRTFTR

#include "tsmbkind.h"

USE OPARRRTM


IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

INTEGER_M :: NGC(JPBAND)
INTEGER_M :: NGS(JPBAND)
INTEGER_M :: NGN(JPGPT)
INTEGER_M :: NGB(JPGPT)

INTEGER_M :: NGM(JPG*JPBAND)
REAL_B ::    WT(JPG)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!  NGC   : INTEGER :
!  NGS   : INTEGER :
!  NGN   : INTEGER :
!  NGB   : INTEGER :
!  NGM   : INTEGER :
!  WT    : REAL    :
!    -------------------------------------------------------------------
END MODULE OYOERRTFTR


