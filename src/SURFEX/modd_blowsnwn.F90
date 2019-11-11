MODULE MODD_BLOWSNW_n

!Purpose: 
!Declare variables and constants necessary to do the blowing snow calculations
!Here are only the variables which depend on the grid!

!Author: Vincent Vionnet
! based on modd_dstn.F90
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE BLOWSNW_t
  
  REAL, DIMENSION(:,:),POINTER    :: XSFSNW        ! Blowing snow fluxes
  ! 1: Number deposition flux 2: Mass deposition flux  3: Streamwise saltation flux
  REAL, DIMENSION(:,:),POINTER    :: XSNW_SUBL       !Sublimation rate
  ! 1: Instantaneous number 2: Instantaneous mass 3: Accumulated Mass
  REAL, DIMENSION(:,:),POINTER    :: XSNW_FTURB       ! Snow surface turbulent flux
  ! 1: Instantaneous number 2: Instantaneous mass 3: Accumulated Mass
  REAL, DIMENSION(:,:),POINTER    :: XSNW_FSED       ! Snow surface sedimentation flux
  ! 1: Instantaneous number 2: Instantaneous mass 3: Accumulated Mass
  REAL, DIMENSION(:,:),POINTER    :: XSNW_FNET       ! Total surface net flux
  !                                                        (salt+ susp)
  ! 1: Instantaneous number 2: Instantaneous mass 3: Accumulated Mass
  REAL, DIMENSION(:,:),POINTER    :: XSNW_FSALT      ! Saltation deposition/erosion flux
  ! 1: streamwise flux 2: Instantaneous 3: Accumulated
  REAL, DIMENSION(:,:),POINTER    :: XSNW_CANO_RGA      ! !Blowing snow radius at canopy levels (m)
  REAL, DIMENSION(:,:,:),POINTER    :: XSNW_CANO_VAR      ! !Blowing snow variables at canopy levels (_/m3)


END TYPE BLOWSNW_t

CONTAINS

SUBROUTINE BLOWSNW_INIT(YBLOWSNW)
TYPE(BLOWSNW_t), INTENT(INOUT) :: YBLOWSNW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BLOWSNW_N:BLOWSNW_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YBLOWSNW%XSFSNW)
  NULLIFY(YBLOWSNW%XSNW_SUBL)
  NULLIFY(YBLOWSNW%XSNW_FTURB)
  NULLIFY(YBLOWSNW%XSNW_FSED)
  NULLIFY(YBLOWSNW%XSNW_FNET)
  NULLIFY(YBLOWSNW%XSNW_FSALT)
IF (LHOOK) CALL DR_HOOK("MODD_BLOWSNW_N:BLOWSNW_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE BLOWSNW_INIT

END MODULE MODD_BLOWSNW_n
