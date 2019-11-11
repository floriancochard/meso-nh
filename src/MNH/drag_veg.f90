!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #######################
       MODULE MODI_DRAG_VEG
!     #######################
!
INTERFACE

SUBROUTINE DRAG_VEG(PTSTEP,PUT,PVT,PTKET,ODEPOTREE, PVDEPOTREE, &
                    HCLOUD,PPABST,PTHT,PRT,PSVT,         &
                    PRHODJ,PZZ,PRUS, PRVS, PRTKES,       &
                    PTHS,PRRS,PSVS)
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT   ! variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTKET           !   at t
LOGICAL,                  INTENT(IN)    :: ODEPOTREE ! Droplet deposition on tree
REAL,                     INTENT(IN)    :: PVDEPOTREE! Velocity deposition on tree
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD       ! Kind of microphysical scheme
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST          !   at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT          !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT            !   at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ       ! Height (z)
!

!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS       ! Sources of Momentum
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTKES           ! Sources of Tke
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS         
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS       
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PTHS          
!
!

END SUBROUTINE DRAG_VEG

END INTERFACE

END MODULE MODI_DRAG_VEG
!
!     ###################################################################
SUBROUTINE DRAG_VEG(PTSTEP,PUT,PVT,PTKET,ODEPOTREE, PVDEPOTREE, &
                    HCLOUD,PPABST,PTHT,PRT,PSVT,         &
                    PRHODJ,PZZ,PRUS, PRVS, PRTKES,       &
                    PTHS,PRRS,PSVS)
!     ###################################################################
!
!!****  *DRAG_VEG_n * -
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     P. Aumond 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2009
!!       C.Lac      07/2011 : Add budgets
!!       S. Donier  06/2015 : bug surface aerosols
!!       C.Lac      07/2016 : Add droplet deposition
!!       C.Lac      10/2017 : Correction on deposition
!!---------------------------------------------------------------
!
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_DYN
USE MODD_DYN_n
USE MODD_VEG_n
USE MODD_BUDGET
USE MODD_PARAM_C2R2
USE MODD_NSV

!
USE MODI_SHUMAN
USE MODD_PGDFIELDS
USE MODD_GROUND_PAR
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_BUDGET

!  
IMPLICIT NONE
!  
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT   ! variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTKET           !   at t
LOGICAL,                  INTENT(IN)    :: ODEPOTREE ! Droplet deposition on tree
REAL,                     INTENT(IN)    :: PVDEPOTREE! Velocity deposition on tree
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD       ! Kind of microphysical scheme
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST          !   at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT          !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT            !   at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ       ! Height (z)
!

!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS       ! Sources of Momentum
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTKES           ! Sources of Tke
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS         
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS       
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS          
!
!
!*       0.2   Declarations of local variables :
!
INTEGER    ::  IIU,IJU,IKU,IKV         ! array size along the k direction 
INTEGER  :: JI, JJ, JK             ! loop index
!
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZH_TREE_PGD ! surface cover types
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLAI_PGD ! surface cover types
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) ::           &
                              ZWORK1, ZWORK2, ZWORK3, ZUT, ZVT,   &
                              ZUS, ZVS, ZTKES, ZTKET
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) ::           &
                              ZCDRAG, ZDENSITY
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2)) ::           &
                              ZVH,ZLAI           !  LAI, Vegetation height
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZT,ZEXN,ZLV,ZCPH                              
LOGICAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) &
            :: GDEP
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZWDEPR,ZWDEPS

!
!
IIU = SIZE(PUT,1)
IJU = SIZE(PUT,2)
IKU = SIZE(PUT,3)
!
!*       0.3     Initialisation de kelkes variables
!	       
ZVH(:,:)=0.
ZLAI(:,:)=0.
ZCDRAG(:,:,:)=0.
ZDENSITY(:,:,:)=0.
!
ALLOCATE(ZH_TREE_PGD(IIU,IJU))
ALLOCATE(ZLAI_PGD(IIU,IJU))
!
CALL MNHGET_SURF_PARAM_n(PH_TREE=ZH_TREE_PGD,PLAI_TREE=ZLAI_PGD)
!
ZVH(:,:)=ZH_TREE_PGD(:,:)
ZLAI(:,:)=ZLAI_PGD(:,:)
!
DEALLOCATE(ZH_TREE_PGD)
DEALLOCATE(ZLAI_PGD)
!
!-------------------------------------------------------------------------------
!
!
!*       1.     COMPUTES THE TRUE VELOCITY COMPONENTS
!	        -------------------------------------
!
ZUT(:,:,:) = PUT(:,:,:) 
ZVT(:,:,:) = PVT(:,:,:) 
ZTKET(:,:,:) = PTKET(:,:,:) 
!-------------------------------------------------------------------------------
!
!*      1.     Computations of wind tendency due to canopy drag
!              ------------------------------------------------
!
!
!
! Ext = - Cdrag  * u- * u- * Sv       tree canopy drag
!       - u'w'(ground)     * Sh       horizontal surfaces (ground)
!
!*      1.1    Drag coefficient by vegetation (Patton et al 2001)
!              ------------------------------
!
GDEP(:,:,:) = .FALSE.
DO JJ=2,(IJU-1)
 DO JI=2,(IIU-1)
   IF (ZVH(JI,JJ) /= 0) THEN
     DO JK=2,(IKU-1) 
         IF ((ZVH(JI,JJ)+PZZ(JI,JJ,2))<PZZ(JI,JJ,JK)) EXIT
         IF ((HCLOUD=='C2R2') .OR.  (HCLOUD=='KHKO')) THEN
           IF ((PRRS(JI,JJ,JK,2) >0.) .AND. (PSVS(JI,JJ,JK,NSV_C2R2BEG+1) >0.)) &
                   GDEP(JI,JJ,JK) = .TRUE.
         ELSE IF (HCLOUD /= 'NONE' .AND. HCLOUD /= 'REVE') THEN
           IF (PRRS(JI,JJ,JK,2) >0.) GDEP(JI,JJ,JK) = .TRUE.
         END IF
         ZCDRAG(JI,JJ,JK)  = 0.2 !0.075
         ZDENSITY(JI,JJ,JK) = MAX((4 * (ZLAI(JI,JJ) *&
                              (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                              (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                              (ZVH(JI,JJ)-(PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)))/&
                              ZVH(JI,JJ)**3)-&
                              (0.30*((ZLAI(JI,JJ) *&
                              (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                              (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                              (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) /&
                              (ZVH(JI,JJ)**3))-ZLAI(JI,JJ))))/&
                              ZVH(JI,JJ), 0.)

                                            
     END DO
   END IF
 END DO
END DO
! To exclude the first vertical level already dealt in rain_ice or rain_c2r2_khko
GDEP(:,:,2) = .FALSE.
!
!*      1.2    Drag force by wall surfaces
!              ---------------------------
!
!* drag force by vertical surfaces
!
ZUS(:,:,:)=  ZUT(:,:,:)/(1 + ZCDRAG(:,:,:)* ZDENSITY(:,:,:)*PTSTEP &
            *SQRT(ZUT(:,:,:)**2+ZVT(:,:,:)**2))
!
ZVS(:,:,:)=  ZVT(:,:,:)/(1 + ZCDRAG(:,:,:)* ZDENSITY(:,:,:)*PTSTEP &
            *SQRT(ZUT(:,:,:)**2+ZVT(:,:,:)**2))
!
PRUS(:,:,:)=PRUS(:,:,:)+((ZUS(:,:,:)-ZUT(:,:,:))*PRHODJ(:,:,:))/PTSTEP
!
PRVS(:,:,:)=PRVS(:,:,:)+((ZVS(:,:,:)-ZVT(:,:,:))*PRHODJ(:,:,:))/PTSTEP
!
IF (ODEPOTREE) THEN
  ZEXN(:,:,:)= (PPABST(:,:,:)/XP00)**(XRD/XCPD)
  ZT(:,:,:)= PTHT(:,:,:)*ZEXN(:,:,:)
  ZLV(:,:,:)=XLVTT +(XCPV-XCL) *(ZT(:,:,:)-XTT)
  ZCPH(:,:,:)=XCPD +XCPV*PRT(:,:,:,1)
  ZWDEPR(:,:,:)= 0.
  ZWDEPS(:,:,:)= 0.
  WHERE (GDEP)
   ZWDEPR(:,:,:)= PVDEPOTREE * PRT(:,:,:,2) * PRHODJ(:,:,:)
  END WHERE
  IF ((HCLOUD=='C2R2') .OR.  (HCLOUD=='KHKO')  .OR.  (HCLOUD=='LIMA')) THEN
   WHERE (GDEP)
   ZWDEPS(:,:,:)= PVDEPOTREE * PSVT(:,:,:,NSV_C2R2BEG+1) * PRHODJ(:,:,:)
   END WHERE
  END IF
  DO JJ=2,(IJU-1)
   DO JI=2,(IIU-1)
     DO JK=2,(IKU-2) 
       IF (GDEP(JI,JJ,JK)) THEN
          PRRS(JI,JJ,JK,2) = PRRS(JI,JJ,JK,2) + (ZWDEPR(JI,JJ,JK+1)-ZWDEPR(JI,JJ,JK))/ &
                             (PZZ(JI,JJ,JK+1)-PZZ(JI,JJ,JK))
         IF ((HCLOUD=='C2R2') .OR.  (HCLOUD=='KHKO').OR.  (HCLOUD=='LIMA')) THEN
          PSVS(JI,JJ,JK,NSV_C2R2BEG+1) =  PSVS(JI,JJ,JK,NSV_C2R2BEG+1) + &
                   (ZWDEPS(JI,JJ,JK+1)-ZWDEPS(JI,JJ,JK))/(PZZ(JI,JJ,JK+1)-PZZ(JI,JJ,JK))
         END IF
       END IF
     END DO
    END DO
   END DO
!
!
END IF
!
IF (LBUDGET_U) CALL BUDGET (PRUS,1,'DRAG_BU_RU')
IF (LBUDGET_V) CALL BUDGET (PRVS,2,'DRAG_BU_RV')
IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7,'DEPOTR_BU_RRC')
IF (LBUDGET_SV) CALL BUDGET (PSVS(:,:,:,NSV_C2R2BEG+1),14+(NSV_C2R2BEG-1),'DEPOTR_BU_RSV')
!
!
!*      3.     Computations of TKE  tendency due to canopy drag
!              ------------------------------------------------

!*      3.1    Creation of TKE by wake
!              -----------------------
!
! from Kanda and Hino (1994)
!
! Ext = + Cd * u+^3  * Sv/Vair        vertical surfaces or trees             
! Ext = - Cd * e * u  * Sv        trees Destruction of TKE due to 
!   small-scale motions forced by leaves from Kanda and Hino (1994)
!
! with Vair = Vair/Vtot * Vtot = (Vair/Vtot) * Stot * Dz
! and  Sv/Vair = (Sv/Stot) * Stot/Vair = (Sv/Stot) / (Vair/Vtot) / Dz
!
!ZTKES(:,:,:)=  (ZTKET(:,:,:) + (ZCDRAG(:,:,:)* ZDENSITY(:,:,:) &
!            *(SQRT(ZUT(:,:,:)**2+ZVT(:,:,:)**2))**3))   /&
!            (1.+(2.*ZCDRAG(:,:,:)* ZDENSITY(:,:,:)*SQRT(ZUT(:,:,:)**2+ZVT(:,:,:)**2)))
ZTKES(:,:,:)=  (ZTKET(:,:,:) + (ZCDRAG(:,:,:)* ZDENSITY(:,:,:) &
            *(SQRT(ZUT(:,:,:)**2+ZVT(:,:,:)**2))**3))*PTSTEP   /&
            (1.+PTSTEP*ZCDRAG(:,:,:)* ZDENSITY(:,:,:)*SQRT(ZUT(:,:,:)**2+ZVT(:,:,:)**2))
!
PRTKES(:,:,:)=PRTKES(:,:,:)+((ZTKES(:,:,:)-ZTKET(:,:,:))*PRHODJ(:,:,:)/PTSTEP)
!
IF (LBUDGET_TKE) CALL BUDGET (PRTKES(:,:,:),5,'DRAG_BU_RTKE')
!
END SUBROUTINE DRAG_VEG
