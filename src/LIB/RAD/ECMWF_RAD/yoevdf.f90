!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOEVDF


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEVDF* CONTAINS CONSTANTS NEEDED BY *VDF....*
!     ------------------------------------------------------------------

REAL_B :: RLAM
REAL_B :: RKAP
REAL_B :: RCHAR
REAL_B :: RVDIFTS
REAL_B :: RTFREEZ
REAL_B :: RZ0ICE
REAL_B :: REPDU2
REAL_B :: REPUST
REAL_B :: RSEZ0H
REAL_B :: RSEZ0Q
REAL_B :: RNUM
REAL_B :: RNUH
REAL_B :: RNUQ
REAL_B :: RENTR
REAL_B :: RPAR
REAL_B :: RPAR1
REAL_B :: RPARSRF
REAL_B :: RPARZI
REAL_B :: RLAMSK
LOGICAL LELWDD

!*     *YOEVDF* CONTAINS CONSTANTS NEEDED BY *VDF....*
!     FOR THE COMPUTATION OF VERTICAL DIFFUSION


!     A.C.M. BELJAARS      E.C.M.W.F.    14/12/89

!     OBUKHOV-L UPDATE     ACMB          26/03/90.   

!     NAME        TYPE     DESCRIPTION
!     ----        ----     -----------

!     *RLAM*      REAL     *ASYMPTOTIC MIXING LENGTH FOR MOMENTUM
!     *RKAP*      REAL     *VONKARMAN CONSTANT
!     *RCHAR*     REAL     *CHARNOCK CONSTANT
!     *RVDIFTS*   REAL     *FACTOR FOR TIME STEP WEIGHTING IN *VDF....*
!     *RTFREEZ*   REAL     *TEMPERATURE AT WHICH SEA STARTS
!                           MELTING/FREEZING
!     *RZ0ICE*    REAL     *ROUGHNESS OVER SEA ICE
!     *REPDU2*    REAL     *MINIMUM VELOCITY DIFFERENCE IN RI-NUMBER
!     *REPDU2*    REAL     *MINIMUM VELOCITY DIFFERENCE IN RI-NUMBER  
!     *REPUST*    REAL     *MINIMUM FRICTION VELOCITY (SECURITY PARAMETER)
!     *RSEZ0H*    REAL     *ROUGHNESS LENGTH FOR HEAT OVER ROUGH SEA   
!     *RSEZ0Q*    REAL     *ROUGHNESS LENGTH FOR MOISTURE OVER ROUGH SEA 
!     *RNUM*      REAL     *SMOOTH SURFACE CONSTANT IN Z0M=RNUM/U*    
!     *RNUH*      REAL     *SMOOTH SURFACE CONSTANT IN Z0H=RNUH/U*  
!     *RNUQ*      REAL     *SMOOTH SURFACE CONSTANT IN Z0Q=RNUQ/U*    
!     *RENTR*     REAL     *ENTRAINMENT CONSTANT          
!     *RPAR*      REAL     *PARAMETER FOR TEMPERATURE EXCESS IN THERMAL 
!                           AT BOUNDARY LAYER TOP      
!     *RPAR1*     REAL     *COEFFICIENT OF (W*)**3 IN WS         
!     *RPARSRF*   REAL     *DEPTH OF SURFACE LAYER AS FRACTION OF PBL-H 
!     *RPARZI*    REAL     *ANSATZ FOR PBL-H IN W* COMPUTATION  
!     *RLAMSK*    REAL     *SKIN LAYER THERMAL BULK CONDUCTIVITY (W/M2/K)
!     *LELWDD*    LOGICAL  *TRUE WHEN LONGWAVE DOWNWARD DERIVATIVE IS
!                           FOR SKIN TEMPERATURE  
!     ------------------------------------------------------------------
END MODULE YOEVDF
