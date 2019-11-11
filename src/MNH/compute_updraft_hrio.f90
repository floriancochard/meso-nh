!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_COMPUTE_UPDRAFT_HRIO
!    ###########################
!
INTERFACE
!
!     #################################################################
      SUBROUTINE COMPUTE_UPDRAFT_HRIO(KKA,KKB,KKE,KKU,KKL, HFRAC_ICE,   &
                                 OENTR_DETR,OMIXUV,                &
                                 ONOMIXLG,KSV_LGBEG,KSV_LGEND,     &
                                 PZZ,PDZZ,                         &
                                 PSFTH,PSFRV,                      &
                                 PPABSM,PRHODREF,PUM,PVM,PTKEM,PWM,&
                                 PTHM,PRVM,PTHLM,PRTM,             &
                                 PSVM,PTHL_UP,PRT_UP,              &
                                 PRV_UP,PRC_UP,PRI_UP,PTHV_UP,     &
                                 PW_UP,PU_UP, PV_UP, PSV_UP,       &
                                 PFRAC_UP,PFRAC_ICE_UP,PRSAT_UP,   &
                                 PTHL_DO, PTHV_DO, PRT_DO,         &
                                 PU_DO, PV_DO, PSV_DO,             &
                                 PEMF,PDETR,PENTR,                 &
                                 PBUO_INTEG,KKLCL,KKETL,KKCTL,     &
                                 PDEPTH)
!     #################################################################
!
!*                    1.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER*1,            INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN) :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL, DIMENSION(:,:), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(:),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography
!
REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PWM        ! w mean wind
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRV_UP,PRC_UP, & ! updraft rv, rc
                                         PRI_UP,PTHV_UP,& ! updraft ri, THv
                                         PW_UP,PFRAC_UP,& ! updraft w, fraction
                                         PFRAC_ICE_UP,&   ! liquid/solid fraction in updraft
                                         PRSAT_UP         ! Rsat
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PTHL_DO,PTHV_DO,PRT_DO,PU_DO,PV_DO

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_DO           ! downdraft scalar var. 
                                         
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! entrainment, detrainment
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(:),  INTENT(INOUT)::  KKLCL,KKETL,KKCTL! LCL, ETL, CTL                                           
REAL, DIMENSION(:),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud


END SUBROUTINE COMPUTE_UPDRAFT_HRIO

END INTERFACE
!
END MODULE MODI_COMPUTE_UPDRAFT_HRIO!     ######spl
      SUBROUTINE COMPUTE_UPDRAFT_HRIO(KKA,KKB,KKE,KKU,KKL,HFRAC_ICE, &
                                 OENTR_DETR,OMIXUV,                  &
                                 ONOMIXLG,KSV_LGBEG,KSV_LGEND,       &
                                 PZZ,PDZZ,                           &
                                 PSFTH,PSFRV,                        &
                                 PPABSM,PRHODREF,PUM,PVM,PTKEM,PWM,  &                                      
                                 PTHM,PRVM,PTHLM,PRTM,               &
                                 PSVM,PTHL_UP,PRT_UP,                &
                                 PRV_UP,PRC_UP,PRI_UP,PTHV_UP,       &
                                 PW_UP,PU_UP, PV_UP, PSV_UP,         &
                                 PFRAC_UP,PFRAC_ICE_UP,PRSAT_UP,     &
                                 PTHL_DO, PTHV_DO, PRT_DO, &
                                 PU_DO, PV_DO, PSV_DO,         &
                                 PEMF,PDETR,PENTR,                   &
                                 PBUO_INTEG,KKLCL,KKETL,KKCTL,       &
                                 PDEPTH )
!     #################################################################
!!
!!****  *COMPUTE_UPDRAFT_HRIO* - calculates caracteristics of the updraft 
!!                         
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to build the updraft model 
!!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      !!     REFERENCE
!!     ---------
!!       Book 1 of Meso-NH documentation (chapter Turbulence)
!!       Soares et al. 2004 QJ
!!
!!     AUTHOR
!!     ------
!!     J.Pergaud
!!     V.Masson : Optimization 07/2010
!!     S. Riette : 07/2010 : modification for reproducibility  
!!     S. Riette may 2011: ice added, interface modified
!!     S. Riette Jan 2012: support for both order of vertical levels
!!     V.Masson, C.Lac : 02/2011 : SV_UP initialized by a non-zero value
!!     Q.Rodier  01/2019 : support RM17 mixing length 
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_PARAM_MFSHALL_n, ONLY : XPRES_UV,XALP_PERT,XCMF,XFRAC_UP_MAX,XA1,XB,&
                          XC,XBETA1
!USE MODD_CMFSHALL, ONLY : XPRES_UV,XALP_PERT,XCMF,XFRAC_UP_MAX,XA1,XALPHA_SEUIL,LNORM_RESOL,&
!                          XCOEF1,XCOEF2,XCOEF3,NBUOY,XB,&
!                          XC,XBETA1
USE MODD_GRID_n, ONLY : XDXHAT, XDYHAT
USE MODD_BLANK
USE MODD_TURB_n, ONLY :CTURBLEN

!USE MODI_COMPUTE_ENTR_DETR
USE MODI_TH_R_FROM_THL_RT_1D
USE MODI_SHUMAN_MF

USE MODI_COMPUTE_BL89_ML

IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER*1,            INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN) :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
! FAUX dans AROME, mais pas grave dans la paramétrisation
REAL, DIMENSION(:,:), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(:),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography
!
REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PWM        ! w mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRV_UP,PRC_UP, & ! updraft rv, rc
                                         PRI_UP,PTHV_UP,& ! updraft ri, THv
                                         PW_UP,PFRAC_UP,& ! updraft w, fraction
                                         PFRAC_ICE_UP,&   ! liquid/solid fraction in updraft
                                         PRSAT_UP         ! Rsat
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PTHL_DO,PTHV_DO,PRT_DO,PU_DO,PV_DO ! environment var.

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_DO           ! environment scalar var. 
                                         
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(:),  INTENT(INOUT) :: KKLCL,KKETL,KKCTL! LCL, ETL, CTL
REAL, DIMENSION(:),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
!                       1.2  Declaration of local variables
!
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::    &
                        ZTHM_F,ZRVM_F                 ! Theta,rv of
                                                      ! updraft environnement
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::    &
                        ZRTM_F, ZTHLM_F, ZTKEM_F,&    ! rt, thetal,TKE,pressure,
                        ZUM_F,ZVM_F,ZRHO_F,      &    ! density,momentum
                        ZPRES_F,ZTHVM_F,ZTHVM,   &    ! interpolated at the flux point
                        ZG_O_THVREF,             &    ! g*ThetaV ref
                        ZW_UP2                        ! w**2  of the updraft
!==================================================================================                        
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZW_UP   ! w  of the updraft
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZWM_M   ! w at mass levels

REAL, DIMENSION(SIZE(PSVM,1),SIZE(PTHM,2),SIZE(PSVM,3)) :: &
                        ZSVM_F ! scalar variables 

                        
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::  &
                        ZTH_UP,                  &    ! updraft THETA 
                        ZRC_MIX, ZRI_MIX              ! guess of Rc and Ri for KF mixture
!==================================================================================                        
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: PTHVREF                       ! THv de référence 
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZBUO                          ! Buoyancy 
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZZDZ                          ! Dz 
!==================================================================================                        

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(SIZE(PSFTH,1) )            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(SIZE(PTHM,1)) :: ZMIX1,ZMIX2

REAL, DIMENSION(SIZE(PTHM,1)) :: ZLUP         ! Upward Mixing length from the ground

INTEGER  :: ISV                ! Number of scalar variables                               
INTEGER  :: JK,JI,JSV          ! loop counters

LOGICAL, DIMENSION(SIZE(PTHM,1)) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(SIZE(PTHM,1))              :: GWORK1
LOGICAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: GWORK2

INTEGER  :: ITEST

REAL, DIMENSION(SIZE(PTHM,1)) :: ZRC_UP, ZRI_UP, ZRV_UP, ZRSATW, ZRSATI
!==================================================================================                        
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZWUP_MEAN    ! 
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZDETR_BUO, ZDETR_RT
!==================================================================================                        

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process
REAL  :: XFRAC_LIM ! surface maximale du thermique

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear
! pour le calcul de la resolution normalisée
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZA1, ZRESOL_NORM 
REAL :: ZRESOL_GRID   
REAL, DIMENSION(SIZE(PTHM,1)) :: MODIF
!

INTEGER :: KTCOUNT_MF ! current model time-step
REAL, DIMENSION (:), ALLOCATABLE :: ZWORK
!
!*    0.3    Declaration of namelists
!            ------------------------
!
!----------------------------------------------------------------------------

! Thresholds for the  perturbation of
! theta_l and r_t at the first level of the updraft

ZTMAX=2.0
ZRMAX=1.E-3
ZEPS=1.E-15
XFRAC_LIM=0.5
!
!------------------------------------------------------------------------

!                     INITIALISATION

! Initialisation of the constants   
ZRDORV   = XRD / XRV   !=0.622
ZRVORD   = (XRV / XRD) 

ZDEPTH_MAX1=3000. ! clouds with depth inferior to this value are keeped untouched
ZDEPTH_MAX2=4000. ! clouds with depth superior to this value are suppressed

!                 Local variables, internal domain
! Internal Domain

!number of scalar variables
ISV=SIZE(PSVM,3)

IF (OENTR_DETR) THEN
   ! si on prend en compte la résolution
   ! dans le calcul de l'entrainement et du détrainement
   !IF(LNORM_RESOL) THEN
        !grid resolution in the AROME CASE
        ZRESOL_GRID=sqrt(XDXHAT(1)*XDYHAT(1))
   !ENDIF

! Initialisation of intersesting Level :LCL,ETL,CTL
  KKLCL(:)=KKE
  KKETL(:)=KKE
  KKCTL(:)=KKE

  !
  ! Initialisation
  !* udraft governing variables
  PEMF(:,:)=0.
  PDETR(:,:)=0.
  PENTR(:,:)=0.

  ! Initialisation
  !* updraft core variables
  PRV_UP(:,:)=0.
  PRC_UP(:,:)=0.
  PRI_UP(:,:)=0.
  PW_UP(:,:)=0.
  ZTH_UP(:,:)=0.
  PFRAC_UP(:,:)=0.
  PTHV_UP(:,:)=0.

  PBUO_INTEG=0.
  ZBUO      =0.

  PFRAC_ICE_UP(:,:)=0.
  PRSAT_UP(:,:)=PRVM(:,:) ! should be initialised correctly but is (normaly) not used

  !cloud/dry air mixture cloud content
  ZRC_MIX = 0.
  ZRI_MIX = 0.

END IF

! Initialisation of environment variables at t-dt
! variables at flux level
ZTHLM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTHLM(:,:))
ZRTM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRTM(:,:))
ZUM_F  (:,:) = MZM_MF(KKA,KKU,KKL,PUM(:,:))
ZVM_F  (:,:) = MZM_MF(KKA,KKU,KKL,PVM(:,:))
ZTKEM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTKEM(:,:))

DO JSV=1,ISV
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  ZSVM_F(:,:,JSV) = MZM_MF(KKA,KKU,KKL,PSVM(:,:,JSV))
END DO
!                     
!          Initialisation of updraft characteristics 
PTHL_UP(:,:)=ZTHLM_F(:,:)
PRT_UP(:,:)=ZRTM_F(:,:)
PU_UP(:,:)=ZUM_F(:,:)
PV_UP(:,:)=ZVM_F(:,:)
PSV_UP(:,:,:)=ZSVM_F(:,:,:)
PSV_DO(:,:,:)=ZSVM_F(:,:,:)
PTHL_DO(:,:)=ZTHLM_F(:,:)
PRT_DO(:,:)=ZRTM_F(:,:)
PU_DO(:,:)=ZUM_F(:,:)
PV_DO(:,:)=ZVM_F(:,:)
PSV_DO(:,:,:)=0.

   ! initiation de l'équation de la dynamique
        !vertical velocity at mass level
        ZWM_M(:,:)=MZF_MF(KKA,KKU,KKL,PWM(:,:))
! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w2,Buoyancy term and mass flux (PEMF)

PTHL_UP(:,KKB)= ZTHLM_F(:,KKB)+MAX(0.,MIN(ZTMAX,(PSFTH(:)/SQRT(ZTKEM_F(:,KKB)))*XALP_PERT))
PRT_UP(:,KKB) = ZRTM_F(:,KKB)+MAX(0.,MIN(ZRMAX,(PSFRV(:)/SQRT(ZTKEM_F(:,KKB)))*XALP_PERT)) 
!------------------------
print*,OENTR_DETR
!------------------------
IF (OENTR_DETR) THEN
  ZTHM_F (:,:) = MZM_MF(KKA,KKU,KKL,PTHM (:,:))
  ZPRES_F(:,:) = MZM_MF(KKA,KKU,KKL,PPABSM(:,:))
  ZRHO_F (:,:) = MZM_MF(KKA,KKU,KKL,PRHODREF(:,:))
  ZRVM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRVM(:,:))

  ! thetav at mass and flux levels
  ZTHVM_F(:,:)=ZTHM_F(:,:)*((1.+ZRVORD*ZRVM_F(:,:))/(1.+ZRTM_F(:,:)))
  ZTHVM(:,:)=PTHM(:,:)*((1.+ZRVORD*PRVM(:,:))/(1.+PRTM(:,:)))

  PTHV_UP(:,:)=ZTHVM_F(:,:)

  ZW_UP2(:,:)=ZEPS
  ZW_UP2(:,KKB) = MAX(0.0001,(2./3.)*ZTKEM_F(:,KKB))

   ! initiation de l'équation de la dynamique
        ! initialisation du vent de l'updraft pour la zone grise
        ZW_UP(:,:)=SQRT(ZW_UP2)

  ! Computation of non conservative variable for the KKB level of the updraft
  ! (all or nothing ajustement)
  PRC_UP(:,KKB)=0.
  PRI_UP(:,KKB)=0.
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,KKB),ZPRES_F(:,KKB), &
             PTHL_UP(:,KKB),PRT_UP(:,KKB),ZTH_UP(:,KKB), &
             PRV_UP(:,KKB),PRC_UP(:,KKB),PRI_UP(:,KKB),ZRSATW(:),ZRSATI(:))

  ! compute updraft thevav and buoyancy term at KKB level
  PTHV_UP(:,KKB) = ZTH_UP(:,KKB)*((1+ZRVORD*PRV_UP(:,KKB))/(1+PRT_UP(:,KKB)))
  ! compute mean rsat in updraft
  PRSAT_UP(:,KKB) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,KKB)) + ZRSATI(:)*PFRAC_ICE_UP(:,KKB)
                                                            
  ! Closure assumption for mass flux at KKB level
  !
  ! calcul différent de la flottabilité
  PTHVREF=300.
  ! c'est la meilleure flottabilité ! 
  !NBUOY=0.
   !IF(NBUOY == 0.) THEN
      ZG_O_THVREF=XG/ZTHVM_F
   !ELSE
   !   ! in AROME XTHVREF does not exist ici 300
   !   ZG_O_THVREF=XG/PTHVREF(1,1) ! on revient à l'état de référence et pas ZTHVM_F
   !ENDIF
   ! Calcul de la fermeture de Julien Pergaut comme limite max de PHY

   DO JK=KKB,KKE-KKL,KKL   !  Vertical loop
    ZZDZ(:,JK)   = MAX(ZEPS,PZZ(:,JK+KKL)-PZZ(:,JK))  ! <== Delta Z between two flux level
  ENDDO

  ! compute L_up
  GLMIX=.TRUE.
  ZTKEM_F(:,KKB)=0.
  IF(CTURBLEN=='RM17') THEN
    ZDUDZ = MZF_MF(KKA,KKU,KKL,GZ_M_W_MF(KKA,KKU,KKL,PUM,PDZZ))
    ZDVDZ = MZF_MF(KKA,KKU,KKL,GZ_M_W_MF(KKA,KKU,KKL,PVM,PDZZ))
    ZSHEAR = SQRT(ZDUDZ*ZDUDZ + ZDVDZ*ZDVDZ)
  ELSE
    ZSHEAR = 0. !no shear in bl89 mixing length
  END IF  
  CALL COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ,ZTKEM_F(:,KKB),ZG_O_THVREF(:,KKB), &
                       ZTHVM_F,KKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
  ZLUP(:)=MAX(ZLUP(:),1.E-10)

  ! Compute Buoyancy flux at the ground
  ZWTHVSURF(:) = (ZTHVM_F(:,KKB)/ZTHM_F(:,KKB))*PSFTH(:)+     &
                (0.61*ZTHM_F(:,KKB))*PSFRV(:)

  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
    !elefant>
       MODIF(:)=tanh(1.83*sqrt(XDXHAT(1)*XDYHAT(1))/ZLUP)
      WHERE (ZWTHVSURF(:)>0.)
        PEMF(:,KKB) = XCMF * MODIF(:) * ZRHO_F(:,KKB) *&
        ((ZG_O_THVREF(:,KKB))*ZWTHVSURF*ZLUP)**(1./3.)
        PFRAC_UP(:,KKB)=MIN(PEMF(:,KKB)/(SQRT(ZW_UP2(:,KKB))*ZRHO_F(:,KKB)),XFRAC_UP_MAX)
        ZW_UP2(:,KKB)=(PEMF(:,KKB)/(PFRAC_UP(:,KKB)*ZRHO_F(:,KKB)))**2
        GTEST(:)=.TRUE.
      ELSEWHERE
        PEMF(:,KKB) =0.
        GTEST(:)=.FALSE.
      ENDWHERE
    !elefant<
ELSE
  GTEST(:)=PEMF(:,KKB+KKL)>0.
END IF

!--------------------------------------------------------------------------

!                        3. Vertical ascending loop
!                           -----------------------
!
! If GTEST = T the updraft starts from the KKB level and stops when GTEST becomes F
!
!
GTESTLCL(:)=.FALSE.
GTESTETL(:)=.FALSE.

   ! CALCUL DE LA RESOLUTION NORMALISEE
   
   !IF(NLES_DTCOUNT .EQ. 0.)THEN        
   !  print*,"PROBLEME NLES_DTCOUNT"
   !  print*,"NLES_DTCOUNT doit être 1."   
   !STOP
   !ENDIF

       WHERE (GTEST .AND. ZLUP(:)>0. )
        ! hauteur des thermiques du pas de temps precedent
        ! dans cette version ZLUP depend du flux de masse dans chaque
        ! colonne, il serait peut-être plus simple de faire une variable
        ! globale sur le domaine       
         ZRESOL_NORM(:)=ZRESOL_GRID/(ZLUP(:))
       ELSEWHERE
         ZRESOL_NORM(:)=0.5 ! ARBITRAIRE
       ENDWHERE
 
!       Loop on vertical level
DO JK=KKB,KKE-KKL,KKL

! IF the updraft top is reached for all column, stop the loop on levels
  ITEST=COUNT(GTEST)
  !IF (ITEST==0) print*,"cycle"
  IF (ITEST==0) CYCLE
!       Computation of entrainment and detrainment with KF90
!       parameterization in clouds and LR01 in subcloud layer


! to find the LCL (check if JK is LCL or not)

  WHERE ((PRC_UP(:,JK)+PRI_UP(:,JK)>0.).AND.(.NOT.(GTESTLCL)))
      KKLCL(:) = JK           
      GTESTLCL(:)=.TRUE.
  ENDWHERE
    
!  Buoyancy is computed on "flux" levels where updraft variables are known
   !================================================================================
   ! CALCUL DE LA FLOTTABILITE
   !================================================================================  

  ! Compute theta_v of updraft at flux level JK    
    
    ZRC_UP(:)   =PRC_UP(:,JK) ! guess
    ZRI_UP(:)   =PRI_UP(:,JK) ! guess 
    ZRV_UP(:)   =PRV_UP(:,JK)
    
    CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
               PPABSM(:,JK),PTHL_UP(:,JK),PRT_UP(:,JK),&
               ZTH_UP(:,JK),ZRV_UP,ZRC_UP,ZRI_UP,ZRSATW(:),ZRSATI(:))            
    
    WHERE (GTEST)
      PTHV_UP   (:,JK) = ZTH_UP(:,JK)*(1.+ZRVORD*ZRV_UP(:))/(1.+PRT_UP(:,JK))
    ENDWHERE ! fin temporaire de GTEST 

   ! test sur le calcul de la flottabilité de 2 manières 
   ! le calcul de la flottabilité est super important
   ! c'est ce qui decide de Wup
   !IF(NBUOY==1) THEN
   ! !WHERE (GTEST) ! comme dans EDKF, mais c'est toujours positif
   !               ! du coup le thermique ne s'arrète pas
   !   ZBUO      (:,JK) = ZG_O_THVREF(:,JK)*(PTHV_UP(:,JK) - PTHVREF(:,JK))    
   ! !ENDWHERE ! fin temporaire de GTEST 
   !  IF(LDUMMY1) THEN
   ! print*,"ZBUO_1(:,",JK,")=",minval(ZBUO(:,JK)),maxval(ZBUO(:,JK))
   !END IF
   !ELSEIF(NBUOY==2) THEN
   ! !WHERE (GTEST)! comme ca devrait être mais trop restrictif, le thermique
   ! !s'arrete d'emblée
   !   ZBUO      (:,JK) = ZG_O_THVREF(:,JK)*(PTHV_UP(:,JK) -2.*ZTHVM_F(:,JK) + PTHVREF(:,JK))
   ! !ENDWHERE ! fin temporaire de GTEST 
   ! IF(LDUMMY1) THEN 
   ! print*,"ZBUO_2(:,",JK,")=",minval(ZBUO(:,JK)),maxval(ZBUO(:,JK))
   !END IF
  !ELSEIF(NBUOY==3) THEN
   ! !WHERE (GTEST) ! un peu plus rigoureux que celle du dessus, mais trop restrictif, le thermique
   ! !s'arrete d'emblée 
   !   ZBUO      (:,JK) = XG/ZTHVM_F(:,JK)*(PTHV_UP(:,JK)-ZTHVM_F(:,JK))-&
   !                      XG/PTHVREF(:,JK)*(ZTHVM_F(:,JK)-PTHVREF(:,JK))    
   ! !ENDWHERE ! fin temporaire de GTEST 
   !IF(LDUMMY1) THEN
   ! print*,"ZBUO_3(:,",JK,")=",minval(ZBUO(:,JK)),maxval(ZBUO(:,JK))
   !END IF
   !ELSEIF(NBUOY==4) THEN
   ! devrait être plus proche de EDKF que NBUOY==0, mais ce n'est pas le cas. 
   ! IF(JK/=KKB) THEN
   !   ZRC_MIX(:,JK) = ZRC_MIX(:,JK-KKL) ! guess of Rc of mixture
   !   ZRI_MIX(:,JK) = ZRI_MIX(:,JK-KKL) ! guess of Ri of mixture
   ! ENDIF
   ! CALL COMPUTE_ENTR_DETR(JK,KKB,KKE,KKL,GTEST,GTESTLCL,HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
   !                        PPABSM(:,:),PZZ(:,:),PDZZ(:,:),ZTHVM(:,:),  &
   !                        PTHLM(:,JK),PRTM(:,JK),ZW_UP2(:,:),         &
   !                        PTHL_UP(:,JK),PRT_UP(:,JK),ZLUP(:),         &
   !                        PRC_UP(:,JK),PRI_UP(:,JK),ZRC_MIX(:,JK),ZRI_MIX(:,JK),                 &
   !                        PENTR(:,JK),PDETR(:,JK),ZBUO(:,JK)    )
   !

   !IF(LDUMMY1) THEN
   ! print*,"ZBUO_4(:,",JK,")=",minval(PBUO_INTEG(:,JK)),maxval(PBUO_INTEG(:,JK))
   !END IF

   !ELSE
    !WHERE (GTEST) ! une dernière possibilité
      ZBUO      (:,JK) = ZG_O_THVREF(:,JK)*(PTHV_UP(:,JK) - ZTHVM_F(:,JK))    
    !ENDWHERE ! fin temporaire de GTEST 
   !IF(LDUMMY1) THEN
   ! print*,"ZBUO_0(:,",JK,")=",minval(ZBUO(:,JK)),maxval(ZBUO(:,JK))
   !END IF
   !ENDIF

    !WHERE (GTEST)      
      ! anomalie de flottabilité * DZ
       PBUO_INTEG(:,JK) = ZBUO(:,JK)*(PZZ(:,JK+KKL)-PZZ(:,JK))      
      ! uniquement pour le calcul de l'ancienne vitesse verticale 
    !ENDWHERE ! fin temporaire de GTEST
 
   !================================================================================
   ! CALCUL DE ZA1
   !================================================================================  

   ! si on a besoin d'un calcul spécifique de ZA1 vs XA1 (constante)
   ! - LNORM_RESOL = calcul de ZA1 en prenant en compte la résolution normalisée

      ZA1(:)=XA1    
      !IF (LNORM_RESOL) THEN
      ! WHERE (GTEST)
      !   ZA1(:)=XCOEF1+XCOEF2*ZRESOL_NORM(:)**(-1*XCOEF3)
      ! ENDWHERE
      !ELSEIF (.NOT. XDUMMY1 == 0. ) THEN 
      ! WHERE (GTEST)
      !   ZA1(:)=XDUMMY1
      ! ENDWHERE
      !END IF
       !
      
   !================================================================================
   ! EQUATION DE LA DYNAMIQUE  
   !================================================================================  
      ! on calcule le vent vertical du thermique dans la zone grise
      ! ZA1 depend de la méthode utilisée
      ALLOCATE(ZWORK(SIZE(ZW_UP,1)))
       ZWORK(:)=0.
       !WHERE (GTEST .AND. PFRAC_UP(:,JK)>0 .AND. PFRAC_UP(:,JK)<=XFRAC_LIM)
       WHERE (GTEST .AND. PFRAC_UP(:,JK)>0 )
         !hypothse => alpha(k)=alpha(k+1)
         ZWCOE(:)  = (1.-PFRAC_UP(:,JK))/(1.-PFRAC_UP(:,JK)+XBETA1) 
         ZBUCOE(:) = 2.*ZWCOE(:)*ZA1
         ZWORK(:)=(ZW_UP(:,JK)-(1.-ZWCOE(:))*PWM(:,JK)-ZWCOE(:)*ZWM_M(:,JK))*&
                  (ZW_UP(:,JK)-(1.-ZWCOE(:))*PWM(:,JK)-ZWCOE(:)*ZWM_M(:,JK))+&
                   ZBUCOE(:)*PBUO_INTEG(:,JK)
       ELSEWHERE
        ZWORK(:)=0
       ENDWHERE        
       !
       WHERE (GTEST)    
       WHERE(ZWORK>=0.)
         ! il s'agit bien de Wu et non pas de Wu-Wm
         ZW_UP(:,JK+KKL)=(1.-ZWCOE(:))*PWM(:,JK+KKL)+ZWCOE(:)*ZWM_M(:,JK)+SQRT(ZWORK(:))
       ELSEWHERE
         ZW_UP(:,JK+KKL)=0.
       ENDWHERE   
        ZW_UP2(:,JK+KKL) = MAX(0.,ZW_UP(:,JK+KKL)*ZW_UP(:,JK+KKL))
        ZWUP_MEAN(:)     = MAX(ZEPS,0.5*(ZW_UP(:,JK+KKL)+ZW_UP(:,JK))-ZWM_M(:,JK))
       ENDWHERE
      DEALLOCATE(ZWORK)
   !
   !================================================================================
   ! CALCUL DE L ENTRAINEMENT ET DU DETRAINEMENT  
   !================================================================================  
   ! Hypothèse alpha jk= alpha point de masse
       !WHERE (GTEST .AND. PFRAC_UP(:,JK)>0. .AND. PFRAC_UP(:,JK)<XFRAC_LIM )
       WHERE (GTEST .AND. PFRAC_UP(:,JK)>0. .AND. ZWUP_MEAN(:)>0 )
   ! ZWUP_MEAN est Wu-Wm au point de masse
        PENTR(:,JK)  = MAX(0., ((1.-PFRAC_UP(:,JK))*XBETA1)/(1.-PFRAC_UP(:,JK)+XBETA1)*(ZA1*ZBUO(:,JK)/&
                               (ZWUP_MEAN(:)*ZWUP_MEAN(:))&
                              - XB &
                              - 1./(ZWUP_MEAN(:))*(PWM(:,JK+KKL)-PWM(:,JK))/(PZZ(:,JK+KKL)-PZZ(:,JK)))&
                          )
        ZDETR_BUO(:) = MAX(0., -((1.-PFRAC_UP(:,JK))*XBETA1)/(1.-PFRAC_UP(:,JK)+XBETA1)*(ZA1*ZBUO(:,JK)/&
                                (ZWUP_MEAN(:)*ZWUP_MEAN(:))))
        ZDETR_RT(:)  = XC*SQRT(MAX(0.,(PRT_UP(:,JK) - ZRTM_F(:,JK)))/MAX(ZEPS,ZRTM_F(:,JK))/ ZWUP_MEAN(:))
        PDETR(:,JK)  = ZDETR_RT(:)+ZDETR_BUO(:)
       ENDWHERE
!
   !================================================================================
   ! VARIABLES THERMODYNAMIQUES   
   !================================================================================  
   ! computation of the updraft characteritics at jk+1
   ! WHERE (GTEST .AND. PFRAC_UP(:,JK)>0. .AND. PFRAC_UP(:,JK)<XFRAC_LIM)
    WHERE (GTEST .AND. PFRAC_UP(:,JK)>0.)
     ZMIX2(:) = (PZZ(:,JK+KKL)-PZZ(:,JK))*PENTR(:,JK)/(1.-PFRAC_UP(:,JK)) !&
    ELSEWHERE
     ZMIX2(:) = 0. !&
    ENDWHERE  ! GTEST


! If the updraft did not stop, compute cons updraft characteritics at jk+KKL
   WHERE (GTEST)  
    PTHL_UP(:,JK+KKL)=(PTHL_UP(:,JK)*(1.-0.5*ZMIX2(:)) + PTHLM(:,JK)*ZMIX2(:)) &
                          /(1.+0.5*ZMIX2(:))   
    PRT_UP(:,JK+KKL) =(PRT_UP (:,JK)*(1.-0.5*ZMIX2(:)) + PRTM(:,JK)*ZMIX2(:))  &
                          /(1.+0.5*ZMIX2(:))
  ENDWHERE
  

  IF(OMIXUV) THEN
    IF(JK/=KKB) THEN
      WHERE(GTEST)
        PU_UP(:,JK+KKL) = (PU_UP (:,JK)*(1-0.5*ZMIX2(:)) + PUM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+KKL)-PUM(:,JK))/PDZZ(:,JK+KKL)+&
                           (PUM(:,JK)-PUM(:,JK-KKL))/PDZZ(:,JK))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+KKL) = (PV_UP (:,JK)*(1-0.5*ZMIX2(:)) + PVM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+KKL)-PVM(:,JK))/PDZZ(:,JK+KKL)+&
                           (PVM(:,JK)-PVM(:,JK-KKL))/PDZZ(:,JK))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE
    ELSE
      WHERE(GTEST)
        PU_UP(:,JK+KKL) = (PU_UP (:,JK)*(1-0.5*ZMIX2(:)) + PUM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+KKL)-PUM(:,JK))/PDZZ(:,JK+KKL))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+KKL) = (PV_UP (:,JK)*(1-0.5*ZMIX2(:)) + PVM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+KKL)-PVM(:,JK))/PDZZ(:,JK+KKL))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE

    ENDIF
  ENDIF
  DO JSV=1,ISV 
     IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
      WHERE(GTEST) 
           PSV_UP(:,JK+KKL,JSV) = (PSV_UP (:,JK,JSV)*(1-0.5*ZMIX2(:)) + &
                        PSVM(:,JK,JSV)*ZMIX2(:))  /(1+0.5*ZMIX2(:))
      ENDWHERE                        
  END DO  
  
! Compute non cons. var. at level JK+KKL
  ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
  ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,JK+KKL),ZPRES_F(:,JK+KKL), &
          PTHL_UP(:,JK+KKL),PRT_UP(:,JK+KKL),ZTH_UP(:,JK+KKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:))
  WHERE(GTEST)
    PRC_UP(:,JK+KKL)=ZRC_UP(:)
    PRV_UP(:,JK+KKL)=ZRV_UP(:)
    PRI_UP(:,JK+KKL)=ZRI_UP(:)
    PRSAT_UP(:,JK+KKL) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,JK+KKL)) + ZRSATI(:)*PFRAC_ICE_UP(:,JK+KKL)
  ENDWHERE
  

! Compute the updraft theta_v, buoyancy and w**2 for level JK+KKL
  WHERE(GTEST)
      PTHV_UP(:,JK+KKL) = ZTH_UP(:,JK+KKL)*((1+ZRVORD*PRV_UP(:,JK+KKL))/(1+PRT_UP(:,JK+KKL)))
 ENDWHERE
   !================================================================================
   ! CALCUL DU FLUX DE MASSE   
   !================================================================================
 ZMIX1(:)=0.
 WHERE(GTEST)
   ! identique dans Rio, EDKF ou dans la zone grise      
    ZMIX1(:)=ZZDZ(:,JK)*(PENTR(:,JK)-PDETR(:,JK))
    PEMF(:,JK+KKL)=PEMF(:,JK)*EXP(ZMIX1(:))
 ENDWHERE  
   !================================================================================
   ! FRACTION DE THERMIQUE   
   !================================================================================
        ! on cherche à savoir s'il y a des vitesses verticales non définies
        ! je n'utilise que ZW_UP2 pour pouvoir avoir une valeur si ZW_UP
        ! n'est pas défini 
IF (maxval(ZW_UP2(:,JK+KKL)) .NE. maxval(ZW_UP2(:,JK+KKL))) STOP 'probleme ici'
! si on est dans la zone grise la définition du flux de masse change
! donc celle de alpha aussi
WHERE(GTEST)
       !attention je ne suis pas au point de masse
WHERE ((SQRT(ZW_UP2(:,JK+KKL))-PWM(:,JK+KKL))>ZEPS)
  PFRAC_UP(:,JK+KKL)=PEMF(:,JK+KKL)/((SQRT(ZW_UP2(:,JK+KKL))-PWM(:,JK+KKL))*ZRHO_F(:,JK+KKL))
ELSEWHERE
  PFRAC_UP(:,JK+KKL)=0.
ENDWHERE
  PFRAC_UP(:,JK+KKL)=MIN(PFRAC_UP(:,JK+KKL),XFRAC_LIM)
ENDWHERE


! calcul des termes environmentaux au point de flux
  WHERE(GTEST) 
     !WHERE(PFRAC_UP(:,JK+KKL)>0 .AND. PFRAC_UP(:,JK+KKL)< XFRAC_LIM) 
     WHERE( PFRAC_UP(:,JK+KKL)>0 ) 
     PTHL_DO(:,JK+KKL)=((PTHLM(:,JK)+PTHLM(:,JK+KKL))*0.5-PFRAC_UP(:,JK+KKL)*PTHL_UP(:,JK+KKL)) &
                         /(1.-PFRAC_UP(:,JK+KKL))
     PRT_DO(:,JK+KKL) =((PRTM(:,JK)+PRTM(:,JK+KKL))*0.5-PFRAC_UP(:,JK+KKL)*PRT_UP(:,JK+KKL)) &
                         /(1.-PFRAC_UP(:,JK+KKL))
     PU_DO(:,JK+KKL)  = ((PUM(:,JK)+PUM(:,JK+KKL))*0.5-PFRAC_UP(:,JK+KKL)*PU_UP(:,JK+KKL)) &
                         /(1.-PFRAC_UP(:,JK+KKL))
     PV_DO(:,JK+KKL)  = ((PVM(:,JK)+PVM(:,JK+KKL))*0.5-PFRAC_UP(:,JK+KKL)*PV_UP(:,JK+KKL)) &
                         /(1.-PFRAC_UP(:,JK+KKL))
     PTHV_DO(:,JK+KKL)=(ZTHVM_F(:,JK+KKL)-PFRAC_UP(:,JK+KKL)*PTHV_UP(:,JK+KKL))&
                          /(1.-PFRAC_UP(:,JK+KKL))
    ENDWHERE  
  ENDWHERE  
  DO JSV=1,ISV 
     IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
      WHERE(GTEST) 
      !WHERE(PFRAC_UP(:,JK+KKL)>0 .AND. PFRAC_UP(:,JK+KKL)<= XFRAC_LIM) 
      WHERE(PFRAC_UP(:,JK+KKL)>0 ) 
           PSV_DO(:,JK+KKL,JSV) = ((PSVM(:,JK,JSV)+PSVM(:,JK+1,JSV))*0.5-PFRAC_UP(:,JK+KKL)*PSV_UP(:,JK+KKL,JSV))&
                                  /(1.-PFRAC_UP(:,JK+KKL))
      ENDWHERE  
      ENDWHERE                        
  ENDDO  

  WHERE(GTEST) 
  PFRAC_UP(:,JK+KKL)=MIN(XFRAC_UP_MAX,PFRAC_UP(:,JK+KKL))
  ENDWHERE  
! Test if the updraft has reach the ETL
  GTESTETL(:)=.FALSE.
  WHERE (GTEST.AND.(PBUO_INTEG(:,JK)<=0.))
      KKETL(:) = JK+KKL
      GTESTETL(:)=.TRUE.
  ENDWHERE

! Test is we have reached the top of the updraft

!  WHERE (GTEST.AND.((ZW_UP2(:,JK+KKL)<=ZEPS).OR.(PEMF(:,JK+KKL)<=ZEPS) .OR. PFRAC_UP(:,JK+KKL)<= XALPHA_SEUIL))
  WHERE ( GTEST .AND. ( (ZW_UP2(:,JK+KKL)<=10E-5) .OR. (PEMF(:,JK+KKL)<=10E-5)) )
      ZW_UP2   (:,JK+KKL)=0.
      PEMF     (:,JK+KKL)=0.
      GTEST    (:)       =.FALSE.
      PTHL_UP  (:,JK+KKL)=ZTHLM_F(:,JK+KKL)
      PRT_UP   (:,JK+KKL)=ZRTM_F(:,JK+KKL)
      PRC_UP   (:,JK+KKL)=0.
      PRI_UP   (:,JK+KKL)=0.
      PRV_UP   (:,JK+KKL)=ZRVM_F (:,JK+KKL)
      PTHV_UP  (:,JK+KKL)=ZTHVM_F(:,JK+KKL)
      
      PFRAC_UP (:,JK+KKL)=0.
      KKCTL    (:)       =JK+KKL
      PTHL_DO  (:,JK+KKL)=ZTHLM_F(:,JK+KKL)
      PRT_DO   (:,JK+KKL)=ZRTM_F(:,JK+KKL)
      PTHV_DO  (:,JK+KKL)=ZTHVM_F(:,JK+KKL)
  ENDWHERE

ENDDO ! boucle JK
!STOP
  PW_UP(:,:)=SQRT(ZW_UP2(:,:))

  PEMF(:,KKB) =0.

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of ZDEPTH_MAX1 m of depth) to 0 (for clouds of ZDEPTH_MAX2 m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
  DO JI=1,SIZE(PTHM,1) 
     PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
  END DO

GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
GWORK2(:,:) = SPREAD( GWORK1(:), DIM=2, NCOPIES=MAX(KKU,KKA) )
ZCOEF(:,:) = SPREAD( (1.-(PDEPTH(:)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1)), DIM=2, NCOPIES=SIZE(ZCOEF,2))
print*,"je sors de compute_updraft"

END SUBROUTINE COMPUTE_UPDRAFT_HRIO
