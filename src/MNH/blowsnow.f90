!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ########################
     MODULE MODI_BLOWSNOW
!    ########################
!
!
INTERFACE
!
      SUBROUTINE BLOWSNOW(HLBCX,HLBCY,PTSTEP,KRR,PPABST,PTHT,PRT,PZZ,PRHODREF,  &
                           PRHODJ,PEXNREF,PRS,PTHS,PSVT,PSVS,PSNWSUBL3D)
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
REAL,                     INTENT(IN)   :: PTSTEP ! Time step :XTSTEP in namelist
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     !  Moist  variable at time t

REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function

REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT) :: PRS   !  Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT  ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS  ! Scalar variable sources

REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSNWSUBL3D  ! Blowing snow sublimation profile 
!
END SUBROUTINE BLOWSNOW
!
END INTERFACE
END MODULE MODI_BLOWSNOW
!
!     ######################################################################
      SUBROUTINE BLOWSNOW(HLBCX,HLBCY,PTSTEP,KRR,PPABST,PTHT,PRT,PZZ,PRHODREF, &
                           PRHODJ,PEXNREF,PRS,PTHS,PSVT,PSVS,PSNWSUBL3D)
!     ######################################################################
!     ##########################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the evolution of blowing snow
!:      particles in Meso-NH and the related variables in Canopy 
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      Subroutine BLOWSNOW_VELGRAV: Computes settling velocity of blown snow
!                                     particles
!!      Subroutine SUBL_BLOWSNOW   : Computes sublimation of blown snow
!!                                    particles
!!      Subroutine SEDIM_BLOWSNOW  : Computes sedimentation of blown snow
!!                                    particles
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    Vionnet, V., Martin, E., Masson, V., Guyomarc’h, G., Naaim-Bouvet, F., Prokop, A., 
!!    Durand, Y. and Lac, C. : 
!!    Simulation of wind-induced snow transport and sublimation in alpine terrain 
!!    using a fully coupled snowpack/atmosphere model, The Cryosphere, 8, 395-415, 2014  
!!
!!    AUTHOR
!!    ------
!!      V. Vionnet      * CNRM/GAME*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/10/11
!!      Implementation in MNH 53 07/2017
!
!*       0.    DECLARATIONS
!  
USE MODE_ll
!
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_BLOWSNOW_n
USE MODD_BLOWSNOW
!
USE MODI_SUBL_BLOWSNOW
USE MODI_SEDIM_BLOWSNOW
USE MODI_BLOWSNOW_VELGRAV
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
REAL,                     INTENT(IN)   :: PTSTEP ! Time step :XTSTEP in namelist
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     !  Moist  variable at time t

REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function

REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT) :: PRS   !  Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT  ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS  ! Scalar variable sources

REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSNWSUBL3D  ! Blowing snow sublimation profile  
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR,JSV       ! Loop index for the moist and scalar variables
INTEGER :: IIB           !  Define the physical domain
INTEGER :: IIE           !
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !

REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSVT   ! scalar variable for microphysics only
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSVS   ! scalar tendency for microphysics only
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZVGK   ! settling velocity for blowing snow variables
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT


ALLOCATE(ZSVT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NSV_SNWEND - NSV_SNWBEG + 1))
ALLOCATE(ZSVS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NSV_SNWEND - NSV_SNWBEG + 1))
ALLOCATE(ZVGK(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NSV_SNWEND - NSV_SNWBEG + 1))

ZSVT = PSVT(:,:,:,NSV_SNWBEG:NSV_SNWEND) 
ZSVS = PSVS(:,:,:,NSV_SNWBEG:NSV_SNWEND)
ZVGK = 0.
!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
PTHS(:,:,:) = PTHS(:,:,:) / PRHODJ(:,:,:)
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) / PRHODJ(:,:,:)
END DO

DO  JSV = 1,SIZE(ZSVS,4)
    ZSVS(:,:,:,JSV) = ZSVS(:,:,:,JSV) / PRHODJ(:,:,:)
END DO

!
!  complete the vertical boundaries
!
PTHS(:,:,IKB-1) = PTHS(:,:,IKB)
PTHS(:,:,IKE+1) = PTHS(:,:,IKE)
!
PRS(:,:,IKB-1,1) = PRS(:,:,IKB,1)
PRS(:,:,IKE+1,1) = PRS(:,:,IKE,1)
PRS(:,:,IKB-1,2:) = 0.0
PRS(:,:,IKE+1,2:) = 0.0
!
PRT(:,:,IKB-1,1) = PRT(:,:,IKB,1)
PRT(:,:,IKE+1,1) = PRT(:,:,IKE,1)
PRT(:,:,IKB-1,2:) = 0.0
PRT(:,:,IKE+1,2:) = 0.0
!
ZSVS(:,:,IKB-1,:) = 0.0
ZSVS(:,:,IKE+1,:) = 0.0
ZSVT(:,:,IKB-1,:) = 0.0
ZSVT(:,:,IKE+1,:) = 0.0
!
!------------------------------------------------------------------------------
!
!*       3.     Settling velocity
!               ------------------------
!
! Compute number-averaged and mass-averaged settling velocity. Used later for: 
!      - sublimation as ventilation velocity 
!      - sedimentation
!
CALL BLOWSNOW_VELGRAV(ZSVT(:,:,1:IKE+1,:),PTHT(:,:,1:IKE+1),    &
                      PPABST(:,:,1:IKE+1),                       &
                      PRHODREF(:,:,1:IKE+1),ZVGK(:,:,1:IKE+1,:))

!------------------------------------------------------------------------------
!
!*       4.     Sublimation (optional)
!               ------------------------
!
IF(LSNOWSUBL) THEN 
! Initialize blowing snow sublimation profile         
     PSNWSUBL3D(:,:,:) = 0.        
!    Compute sublimation for MNH levels        
     CALL SUBL_BLOWSNOW(PZZ, PRHODJ,PRHODREF, PEXNREF, PPABST,          &
                         PTHT, PRT(:,:,:,1), PRT(:,:,:,2),PRT(:,:,:,3),  &
                         PRT(:,:,:,4), PRT(:,:,:,5),PRT(:,:,:,6),        & 
                         ZSVT,PTHS,PRS(:,:,:,1),ZSVS,PSNWSUBL3D,ZVGK(:,:,:,2)  )  
END IF
!------------------------------------------------------------------------------
!
!*       5.     Sedimentation 
!               ------------------------
!
CALL SEDIM_BLOWSNOW(PTHT(IIB:IIE,IJB:IJE,IKB:IKE), PTSTEP,&
                  PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),         &
                  PPABST(IIB:IIE,IJB:IJE,IKB:IKE),           &
                  PZZ(IIB:IIE,IJB:IJE,IKB:IKE+1),            &
                  ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:),           &
                  ZSVS(IIB:IIE,IJB:IJE,IKB:IKE,:),ZVGK(IIB:IIE,IJB:IJE,IKB:IKE,:))
!                  
!-------------------------------------------------------------------------------
!
!
!*      6.     SWITCH BACK TO THE PROGNOSTIC VARIABLES
!               ---------------------------------------
!
PTHS(:,:,:) = PTHS(:,:,:) * PRHODJ(:,:,:)
!
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) * PRHODJ(:,:,:)
END DO

DO  JSV = 1,SIZE(ZSVS,4)
    PSVS(:,:,:,JSV+NSV_SNWBEG-1) = ZSVS(:,:,:,JSV) * PRHODJ(:,:,:)
END DO
DO JSV = NSV_SNWBEG, NSV_SNWEND
     PSVT(:,:,:,JSV) = PSVS(:,:,:,JSV) * PTSTEP / PRHODJ(:,:,:)
END DO


DEALLOCATE(ZSVS)
DEALLOCATE(ZSVT)

END SUBROUTINE BLOWSNOW
