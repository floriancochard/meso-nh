!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ###############
      MODULE MODD_LES_BUDGET
!     ###############
!
!!****  *MODD_LES_BUDGET* - declaration of LES variables for budget computations
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Sept. 29 2002
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
REAL :: XTIME_LES_BU
! time spent in subgrid LES computations in this time-step for budgets
!
REAL :: XTIME_LES_BU_PROCESS
! time spent in subgrid LES computations in this time-step for budgets
! for this process only (advection, microphysics, etc...)
!
REAL :: XCURRENT_TSTEP
!-------------------------------------------------------------------------------
!
!* variables of current model needed in budget computations
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_L_O_EXN_CP
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RHODJ
!
LOGICAL :: LCURRENT_USERV
LOGICAL :: LCURRENT_USERC
LOGICAL :: LCURRENT_USERR
LOGICAL :: LCURRENT_USERI
LOGICAL :: LCURRENT_USERS
LOGICAL :: LCURRENT_USERG
LOGICAL :: LCURRENT_USERH
!
INTEGER :: NCURRENT_RR
!
!-------------------------------------------------------------------------------
!
!* source of U, V, W, Th, R* from advection term
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RUS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RVS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RWS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RTHS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RTKES
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XCURRENT_RRS
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XCURRENT_RSVS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RTHLS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XCURRENT_RRTS
!
!-------------------------------------------------------------------------------
!
!* anomaly fields of all variables on MESONH grid
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XU_ANOM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XV_ANOM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XW_ANOM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XTHL_ANOM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XRT_ANOM
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XSV_ANOM
!
!-------------------------------------------------------------------------------
!
!* budget terms fot the current time-step
!
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_RES_KE    ! total production terms
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_RES_WThl  ! (including resolved
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_RES_Thl2  !  turbulent fluxes)
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_RES_WRt   ! for the resolved
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_RES_Rt2   ! fluxes, variances
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_RES_ThlRt ! and covariances
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_BU_RES_Sv2 !
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X_LES_BU_RES_WSv !
REAL, DIMENSION(:,:), ALLOCATABLE :: X_LES_BU_SBG_Tke
!
!-------------------------------------------------------------------------------
!
!* index for each processus taken into account in the budgets
!
INTEGER :: NLES_TOTADV
INTEGER :: NLES_RELA
INTEGER :: NLES_RAD
INTEGER :: NLES_GRAV
INTEGER :: NLES_COR
INTEGER :: NLES_MICR
INTEGER :: NLES_HTURB
INTEGER :: NLES_VTURB
INTEGER :: NLES_FORC
INTEGER :: NLES_PRES
INTEGER :: NLES_DIFF
INTEGER :: NLES_CURV
INTEGER :: NLES_PREF
INTEGER :: NLES_DP
INTEGER :: NLES_TP
INTEGER :: NLES_TR
INTEGER :: NLES_DISS
INTEGER :: NLES_TEND
INTEGER :: NLES_MISC
INTEGER :: NLES_ADVR
INTEGER :: NLES_ADVM
INTEGER :: NLES_NEST
!
INTEGER :: NLES_TOT
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_LES_BUDGET
