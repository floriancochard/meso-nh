!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOEAERC


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERC* - TEGEN/GISS AEROSOL CLIMATOLOGY
!     ------------------------------------------------------------------

REAL_B :: RAERBC(72,46,12), RAEROR(72,46,12), RAERSD(72,46,12)
REAL_B :: RAERSS(72,46,12), RAERSU(72,46,12)
REAL_B :: REPAER
REAL_B :: RSINCT(46)      , RSINCV(24)
REAL_B :: RTAEBC(72,46)   , RTAEOR(72,46)   , RTAESD(72,46)
REAL_B :: RTAESS(72,46)   , RTAESU(72,46)   , RTAEVO(46)


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/12/21

!  NAME     TYPE    PURPOSE
!  ----  :  ----  : ---------------------------------------------------
! RAERBC : REAL   : OPTICAL THICKNESS BLACK CARBON TYPE AEROSOL
! RAEROR : REAL   : OPTICAL THICKNESS ORGANIC TYPE AEROSOL
! RAERSD : REAL   : OPTICAL THICKNESS SOIL-DUST TYPE AEROSOL
! RAERSS : REAL   : OPTICAL THICKNESS SEA-SALT TYPE AEROSOL
! RAERSU : REAL   : OPTICAL THICKNESS SULFATE-TYPE AEROSOL
! REPAER : REAL   : SECURITY PARAMETER FOR AEROSOLS
! RSINCT : REAL   : SINE OF LATITUDE IN TEGEN CLIMATOLOGY
! RSINCV : REAL   : SINE OF LATITUDE IN GISS VOLCANIC CLIMATOLOGY
! RTAEBC : REAL   : TIME-INTERPOLATED BLACK CARBON OPTICAL THICKNESS
! RTAEOR : REAL   : TIME-INTERPOLATED ORGANIC OPTICAL THICKNESS
! RTAESD : REAL   : TIME-INTERPOLATED SOIL-DUST OPTICAL THICKNESS
! RTAESS : REAL   : TIME-INTERPOLATED SEA-SALT OPTICAL THICKNESS
! RTAESU : REAL   : TIME-INTERPOLATED SULFATE CARBON OPTICAL THICKNESS
! RTAEVO : REAL   : TIME-INTERPOLATED STRATOS.VOLCANIC OPTICAL THICKNESS
!     ------------------------------------------------------------------
END MODULE YOEAERC


