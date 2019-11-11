!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOEVEG

#include "tsmbkind.h"

USE PARDIM


IMPLICIT NONE

SAVE

!       ----------------------------------------------------------------
!**    ** *YOEVEG* CONTAINS CONSTANTS VEGETATION PARAMETERS
!       ----------------------------------------------------------------

REAL_B :: RCVA
REAL_B :: RCVB
REAL_B :: RCVC
REAL_B :: RCVBC
REAL_B :: RVLT
REAL_B :: RVXPKLT
REAL_B :: RVXMKLT
REAL_B :: RVKC
REAL_B :: RVABC
REAL_B :: RVRAD
REAL_B :: RVINTER
REAL_B :: RVZ0SN
REAL_B :: RCEPSW
REAL_B :: REPSR
REAL_B :: REPEVAP

REAL_B :: RVROOTSA(JPCSS)
!*     *YOEVEG* CONTAINS VEGETATION PARAMETERS,
!     USED IN *VDF...* AND *SRF...*.

!     A.C.M. BELJAARS      E.C.M.W.F.    14/12/89

!     NAME        TYPE     DESCRIPTION
!     ----        ----     -----------

!     *RCVA*      REAL     *CONSTANT TO DEFINE THE STOMATAL RESISTANCE
!     *RCVB*      REAL     *CONSTANT TO DEFINE THE STOMATAL RESISTANCE
!     *RCVC*      REAL     *MINIMUM STOMATAL RESITANCE
!     *RCVBC*     REAL     * RCVB X RCVC
!     *RVLT*      REAL     *LEAF AREA INDEX
!     *RVXPKLT*   REAL     *EXP(RCVKLT)
!     *RVXMKLT*   REAL     *EXP(-RCVKLT)
!     *RVKC*      REAL     * CVK X RCVC
!     *RVABC*     REAL     * (RCVA+RCVBC)/RCVC
!     *RVROOTSA*  REAL     *PERCENTAGE OF ROOTS IN EACH SOIL LAYER
!     *RVRAD*     REAL     *FRACTION OF THE NET S.W. RADIATION
!                           CONTRIBUTING TO P.A.R.
!     *RVINTER*   REAL     *EFFICIENCY OF INTERCEPTION OF PRECIPITATION
!     *RVZ0SN*    REAL     *FRACTION OF VEGETATION COVERED BY SNOW
!     *RCEPSW*     REAL     *MINIMUM RELATIVE HUMIDITY
!     *REPSR*     REAL     *MINIMUM VALUE FOR SHORT WAVE RADIATION IN
!                          THE CANOPY RESISTANCE COMPUTATION
!     *REPEVAP*   REAL     *MINIMUM ATMOSPHERIC DEMAND
!     ------------------------------------------------------------------
END MODULE YOEVEG
