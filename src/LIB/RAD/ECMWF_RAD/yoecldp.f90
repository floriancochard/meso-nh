!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOECLDP


#include "tsmbkind.h"

IMPLICIT NONE

SAVE


!     -----------------------------------------------------------------
!     ** YOECLDP - CONTROL PARAMETERS FOR PROGNOSTIC CLOUD SCHEME
!     -----------------------------------------------------------------

!     * E.C.M.W.F. PHYSICS PACKAGE *

!     C. JAKOB        E.C.M.W.F.          94/02/07


!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *RAMID*   REAL      BASE VALUE FOR CALCULATION OF RELATIVE 
!                         HUMIDITY THRESHOLD FOR ONSET OF STRATIFORM
!                         CONDENSATION (TIEDTKE, 1993, EQUATION 24)
!     *RCLDIFF* REAL      DIFFUSION-COEFFICIENT FOR EVAPORATION BY
!                         TURBULENT MIXING (IBID., EQU. 30)
!     *RCLCRIT* REAL      BASE VALUE OF CRITICAL CLOUD WATER CONTENT 
!                         FOR CONVERSION TO RAIN (SUNDQUIST, 1988)
!     *RKCONV*  REAL      BASE VALUE FOR CONVERSION COEFFICIENT (IBID.)
!     *RPRC1*   REAL      COALESCENCE CONSTANT (IBID.)
!     *RPRC2*   REAL      BERGERON-FINDEISEN CONSTANT (IBID.)
!     *RCLDMAX* REAL      MAXIMUM CLOUD WATER CONTENT
!     *RPECONS* REAL      EVAPORATION CONSTANT AFTER KESSLER 
!                         (TIEDTKE, 1993, EQU.35)
!     *RTAUMEL* REAL      RELAXATION TIME FOR MELTING OF SNOW
!     *RENTRTU* REAL      ENTRAINMENT PARAMETER (DEARDOFF, 1976)
!     *RENTRRA* REAL      ENTRAINMENT PARAMETER FOR RADIATION EFFECT
!     *RAMIN*   REAL      LIMIT FOR A
!     *RLMIN*   REAL      LIMIT FOR L
!     *RASMICE*   REAL    COEFFICIENT FOR "SMALL ICE" CALCULATION 
!                         (MACFARQUHAR AND HEYMSFIELD, JAS, 1997)
!     *RBSMICE*   REAL    COEFFICIENT FOR "SMALL ICE" CALCULATION 
!                         (MACFARQUHAR AND HEYMSFIELD, JAS, 1997)
REAL_B :: RAMID
REAL_B :: RCLDIFF
REAL_B :: RCLCRIT
REAL_B :: RKCONV
REAL_B :: RPRC1
REAL_B :: RPRC2
REAL_B :: RCLDMAX
REAL_B :: RPECONS
REAL_B :: RTAUMEL
REAL_B :: RENTRTU
REAL_B :: RENTRRA
REAL_B :: RAMIN
REAL_B :: RLMIN
REAL_B :: RASMICE
REAL_B :: RBSMICE
END MODULE YOECLDP
