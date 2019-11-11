!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ######################################
     MODULE MODI_THRVRCRI_FROM_THLRTRCRI
!    ######################################
!
INTERFACE
      
!     #################################################################
      SUBROUTINE THRVRCRI_FROM_THLRTRCRI( KRR,                         &
                                       PTHL, PR, PLVOCPEXN,PLSOCPEXN,  &
                                       PTH, PRV                        )
!     #################################################################

!*               1.1  Declaration of Arguments
!

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.

REAL, DIMENSION(:,:,:), INTENT(IN)   :: PTHL     ! theta
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLVOCPEXN,PLSOCPEXN     ! L/(cp*Pi)

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PTH, PRV
!
END SUBROUTINE THRVRCRI_FROM_THLRTRCRI

END INTERFACE
!
END MODULE MODI_THRVRCRI_FROM_THLRTRCRI

                                       
!     #################################################################
      SUBROUTINE THRVRCRI_FROM_THLRTRCRI( KRR,                             &
                                       PTHL, PR, PLVOCPEXN,PLSOCPEXN,  &
                                       PTH, PRV                        )
!     #################################################################
!
!!
!!****  *THRVRCRI_FROM_THLRTRCRI* - 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!     
!!     V. Masson     *CNRM*
!!
!!    MODIFICATIONS
!!    -------------
!!     Original   06/2011
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.

REAL, DIMENSION(:,:,:), INTENT(IN)   :: PTHL     ! theta
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLVOCPEXN,PLSOCPEXN     ! L/(cp*Pi)

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PTH, PRV
!
!-------------------------------------------------------------------------------
!
!
IF ( KRR == 0 ) THEN
  PTH(:,:,:) = PTHL(:,:,:)
  PRV(:,:,:) = 0.
ELSE IF (KRR==1) THEN
  PTH(:,:,:) = PTHL(:,:,:)
  PRV(:,:,:) = PR  (:,:,:,1)
ELSE IF (KRR>=4) THEN
    ! Rnp at t
    PRV(:,:,:)  = PR(:,:,:,1)  - PR(:,:,:,2)  - PR(:,:,:,4)
    ! Theta_l at t
    PTH(:,:,:)  = PTHL(:,:,:)  + PLVOCPEXN(:,:,:) * PR(:,:,:,2) &
                               + PLSOCPEXN(:,:,:) * PR(:,:,:,4)
ELSE IF (KRR>=2) THEN
    ! Rnp at t
    PRV(:,:,:)  = PR(:,:,:,1)  - PR(:,:,:,2)
    ! Theta_l at t
    PTH(:,:,:)  = PTHL(:,:,:)  + PLVOCPEXN(:,:,:) * PR(:,:,:,2)
END IF
! 
!-------------------------------------------------------------------------------
!
END SUBROUTINE THRVRCRI_FROM_THLRTRCRI
