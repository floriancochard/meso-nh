!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      MODULE MODI_GRAVITY_IMPL
!     #####################
!
INTERFACE
      SUBROUTINE GRAVITY_IMPL (HLBCX, HLBCY,KRR,KRRL,KRRI, PTSTEP,             &
                            PTHT, PRT, PTHVREF, PRHODJ,                        &
                            PRWS, PRTHS, PRRS, PRTHS_CLD, PRRS_CLD             )
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL,                     INTENT(IN)    :: PTSTEP
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRWS
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRTHS
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRRS 
                                          ! Sources terms 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRTHS_CLD
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRRS_CLD
                                          ! tendencies from previous time-step
                                          ! cloud processes
!
END SUBROUTINE GRAVITY_IMPL
!
END INTERFACE
!
END MODULE MODI_GRAVITY_IMPL 
!     ##########################################################################
      SUBROUTINE GRAVITY_IMPL (HLBCX, HLBCY,KRR,KRRL,KRRI, PTSTEP,             &
                            PTHT, PRT, PTHVREF, PRHODJ,                        &
                            PRWS, PRTHS, PRRS, PRTHS_CLD, PRRS_CLD             )
!     ##########################################################################
!
!!****  *GRAVITY_IMPL * - routine to estimate gravity term using future buoyancy
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson        * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2011 
!!      Q.Rodier 06/15 correction on budget
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GRAVITY  
USE MODI_ADV_BOUNDARIES
USE MODD_BUDGET
USE MODI_BUDGET
!
!-------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL,                     INTENT(IN)    :: PTSTEP
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRWS
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRTHS
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRRS 
                                          ! Sources terms 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRTHS_CLD
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRRS_CLD
                                          ! tendencies from previous time-step
                                          ! cloud processes
!
!
!*       0.2   declarations of local variables
!
!  
! Tendencies of W due to gravity
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZRWS_GRAV
! Guess of future theta
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZTH
! Guess of future mixing ratios
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),SIZE(PRT,4)) :: ZR
!
INTEGER :: JR
!
!-------------------------------------------------------------------------------
!
ZRWS_GRAV = 0.
ZR = 0.
!
! guess of Theta at future time-step
ZTH(:,:,:) = (PRTHS(:,:,:) + PRTHS_CLD(:,:,:)) / PRHODJ * PTSTEP
DO JR = 1, KRR
 ZR(:,:,:,JR) = (PRRS(:,:,:,JR) + PRRS_CLD(:,:,:,JR)) / PRHODJ * PTSTEP
END DO
!
CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZTH, PTHT )    
DO JR = 1, KRR
  CALL ADV_BOUNDARIES (HLBCX, HLBCY, ZR(:,:,:,JR), PRT(:,:,:,JR))
END DO
!
! ====> A vérifier si c'est nécessaire d'échanger les champs
!       A priori, je dirai non.
!
! gravity effect on vertical speed
CALL GRAVITY ( KRR,KRRL, KRRI, ZTH, ZR, PRHODJ, PTHVREF, ZRWS_GRAV(:,:,:) )
!
PRWS(:,:,:) = PRWS(:,:,:) + ZRWS_GRAV(:,:,:)
!
IF (LBUDGET_W) CALL BUDGET (PRWS,3,'GRAV_BU_RW')
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GRAVITY_IMPL
