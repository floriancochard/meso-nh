!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOEMETH


#include "tsmbkind.h"

IMPLICIT NONE

SAVE


!     -----------------------------------------------------------------
!     ** YOEMETH - CONTROL PARAMETERS FOR METHANE OXIDATION
!     -----------------------------------------------------------------

!     * E.C.M.W.F. PHYSICS PACKAGE *

!     C. JAKOB        E.C.M.W.F.          98/04/07


!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *RALPHA1* REAL      CONSTANT IN TIME SCALE 1 CALCULATIONS
!     *RALPHA2* REAL      CONSTANT IN TIME SCALE 2 CALCULATIONS
!     *RQLIM*   REAL      UPPERLIMIT OD SPECIFIC HUMIDITY FOR 
!                         METHANE OXIDATION TO BE ACTIVE
!     *RPBOTOX* REAL      PRESSURE BELOW WHICH METHANE OXIDATION
!                         IS CONSIDERED ACTIVE
!     *RPBOTPH* REAL      PRESSURE BELOW WHICH H2O PHOTOLYSIS
!                         IS CONSIDERED ACTIVE
!     *RPTOPOX* REAL      PRESSURE BELOW WHICH SHORTEST TIME SCALE
!                         IS USED IN METHANE OXIDATION
!     *RPTOPPH* REAL      PRESSURE BELOW WHICH SHORTEST TIME SCALE
!                         IS USED IN H2O PHOTOLYISIS
!     *RALPHA3* REAL      CONSTANT IN TIME SCALE 2 CALCULATIONS
!     *RLOGPPH* REAL      CONSTANT IN TIME SCALE 2 CALCULATIONS

REAL_B :: RALPHA1
REAL_B :: RALPHA2
REAL_B :: RQLIM
REAL_B :: RPBOTOX
REAL_B :: RPBOTPH
REAL_B :: RPTOPOX
REAL_B :: RPTOPPH
REAL_B :: RALPHA3
REAL_B :: RLOGPPH

END MODULE YOEMETH
