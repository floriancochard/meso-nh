!!   #######################################
SUBROUTINE BLOWSNW_DIFFK(PK, PTKE, PVGK,KSNW,PKSNW)

!!   Compute particle eddy diffusivity for number and momentum
!!     ZZETA is the ratio between eddy diffusivity for momentum
!!          znd eddy diffusivity for number and mass.
!!    The value 0.25 is based on a sensititivy analysis using
!!    vertical profile of blowing snow fluxes and concentration
!!    measured using SPC and mechanical snow trap at Col du Lac
!!    Blanc experimental site. More in details in the PhD thesis
!!    of V. Vionnet : Etudes du transport de la neige par le vent en
!!    conditions alpines : observations et modelisation a l'aide d'un
!!    modele couple atmosphere/manteau neigeux

USE MODD_BLOWSNW_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE
  !
  !*       0.1   Declarations of dummy arguments :
  !
REAL, DIMENSION(:,:),   INTENT(IN)    :: PK           ! flow eddy diffusivity
REAL, DIMENSION(:,:),   INTENT(IN)    :: PTKE         ! turbulent kineti energy 
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PVGK         ! particle fall speed
INTEGER,                INTENT(IN)    :: KSNW         ! number of blowing snow variables
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PKSNW      ! particle eddy diffusivity
                                                    ! including inertial effects
  !
  !
  !*       0.2   Declarations of local variables :
  !   
REAL                              :: ZBETAC
REAL                              :: ZZETA
INTEGER                           :: JSV ! Loop counter
LOGICAL                           ::  ONEW_K
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
  !
  !
  !*       1.0   Initialize constants
  !
IF (LHOOK) CALL DR_HOOK('BLOWSNW_DIFFK',0,ZHOOK_HANDLE)


ZBETAC = 1.
ONEW_K=.FALSE.

  !
  !
  !*       2.0 Compute eddy diffusivity for number and mass
  !
IF(ONEW_K) THEN      
   DO JSV =1,KSNW
            PKSNW(:,:,JSV) = PK(:,:)*(1.+(ZBETAC*PVGK(:,:,JSV))**2./(1./3.*PTKE(:,:)))**(-1.)
   ENDDO
ELSE
!
   ZZETA = 0.25
   DO JSV =1,KSNW
       PKSNW(:,:,JSV)=PK(:,:)/XRSNOW_SBL
   ENDDO   
ENDIF

IF (LHOOK) CALL DR_HOOK('BLOWSNW_DIFFK',1,ZHOOK_HANDLE)

END SUBROUTINE BLOWSNW_DIFFK
