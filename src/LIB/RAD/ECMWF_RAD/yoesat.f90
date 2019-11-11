!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOESAT


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESAT* - CONTROL OF SIMULATION OF SATELLITE RADIANCES
!     -----------------------------------------------------------------

INTEGER_M :: NGEO
INTEGER_M :: NPOLO
REAL_B :: RGALT(5)
REAL_B :: RGNAD(5)
REAL_B :: RGNOR(5)
REAL_B :: RGSOU(5)
REAL_B :: RGWST(5)
REAL_B :: RGEAS(5)
LOGICAL LGEOSE
LOGICAL LGEOSW
LOGICAL LGMS
LOGICAL LINDSA
LOGICAL LMTO
LOGICAL LNOAA
LOGICAL LNOAB
LOGICAL LNOAC
LOGICAL LNOAD


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME    TYPE     PURPOSE
!  ---- :  ----   : ---------------------------------------------------
! NGEO  : INTEGER : NUMBER OF GEOSTATIONARY SATELLITES
! NPOLO : INTEGER : NUMBER OF POLAR-ORBITER SATELLITES
! LGEOSE: LOGICAL : SWITCH FOR GEOS-EAST
! LGEOSW: LOGICAL : SWITCH FOR GEOS-WEST
! LGMS  : LOGICAL : SWITCH FOR GMS
! LINDSA: LOGICAL : SWITCH FOR INDSAT
! LMTO  : LOGICAL : SWITCH FOR METEOSAT
! LNOAA : LOGICAL : SWITCH FOR NOAA-A
! LNOAB : LOGICAL : SWITCH FOR NOAA-B
! LNOAC : LOGICAL : SWITCH FOR NOAA-C
! LNOAD : LOGICAL : SWITCH FOR NOAA-D
! RGALT : REAL    : ALTITUDE ABOVE EARTH'S SURFACE OF GEOSTAT. SATELLITE
! RGNAD : REAL    : LONGITUDE OF NADIR OF GEOSTATIONARY SATELLITE
! RGNOR : REAL    : LATITUDE OF NORTH LIMIT OF FIELD OF VIEW
! RGSOU : REAL    : LATITUDE OF SOUTH LIMIT OF FIELD OF VIEW
! RGWST : REAL    : LATITUDE OF WEST  LIMIT OF FIELD OF VIEW
! RGEAS : REAL    : LATITUDE OF EAST  LIMIT OF FIELD OF VIEW
!     -----------------------------------------------------------------
END MODULE YOESAT
