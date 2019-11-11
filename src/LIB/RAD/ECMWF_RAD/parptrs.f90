!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARPTRS


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

INTEGER_M,parameter:: max_satellites=6
INTEGER_M,parameter:: max_ddrs=2
INTEGER_M,parameter:: max_statid_len=8
INTEGER_M,parameter:: NSYNOP=1
INTEGER_M,parameter:: NAIREP=2
INTEGER_M,parameter:: NSATOB=3
INTEGER_M,parameter:: NDRIBU=4
INTEGER_M,parameter:: NTEMP=5
INTEGER_M,parameter:: NPILOT=6
INTEGER_M,parameter:: NSATEM=7
INTEGER_M,parameter:: NPAOB=8
INTEGER_M,parameter:: NSCATT=9
INTEGER_M,parameter:: NTOVS=10
END MODULE PARPTRS
