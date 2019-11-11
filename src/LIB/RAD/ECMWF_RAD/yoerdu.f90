!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE YOERDU


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOERDU* - CONTROL, PARAMETERS AND SECURITY IN RADIATION
!     ------------------------------------------------------------------

INTEGER_M :: NUAER
INTEGER_M :: NTRAER
INTEGER_M :: NIMP
INTEGER_M :: NOUT
REAL_B :: RCDAY
REAL_B :: R10E
REAL_B :: REPLOG
REAL_B :: REPSC
REAL_B :: REPSCO
REAL_B :: REPSCQ
REAL_B :: REPSCT
REAL_B :: REPSCW
REAL_B :: DIFF


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! NUAER  : INTEGER   NUMBER OF ABSORBER AMOUNTS        W OR W/O AEROSOLS
! NTRAER : INTEGER   NUMBER OF TRANSMISSION FUNCTIONS  W OR W/O AEROSOLS
! NIMP   : INTEGER   INDEX FOR EXTRA PRINTS WITHIN RADIATION CODE
! NOUT   : INTEGER   UNIT NUMBER FOR THE EXTRA PRINTS
! RCDAY  : REAL
! CCO2   : REAL      CONVERSION COEFFICIENT FOR CO2 IN S.W. CODE
! CH2O   : REAL      CONVERSION COEFFICIENT FOR H2O IN S.W. CODE
! R10E   : REAL      DECIMAL/NATURAL LOG.FACTOR
! DIFF   : REAL      DIFFUSIVITY FACTOR
!-SECURITY THRESHOLDS
! REPLOG : REAL      SEC. EPSILON FOR ABS.AMOUNT IN LAPLACE TRANSFORM
! REPSC  : REAL      SEC. EPSILON FOR CLOUD COVER
! REPSCO : REAL      SEC. EPSILON FOR OZONE AMOUNT
! REPSCQ : REAL      SEC. EPSILON FOR WATER VAPOR
! REPSCT : REAL      SEC. EPSILON FOR SHORTWAVE OPTICAL THICKNESS
! REPSCW : REAL      SEC. EPSILON FOR CLOUD LIQUID WATER PATH

!     -----------------------------------------------------------------
END MODULE YOERDU
