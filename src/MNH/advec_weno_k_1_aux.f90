!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##############################
      MODULE MODI_ADVEC_WENO_K_1_AUX
!     ##############################
!
INTERFACE
!
FUNCTION UP_UX(PSRC, PRUCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR ! output src term
END FUNCTION UP_UX
!
FUNCTION UP_MX(PSRC, PRUCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR ! output src term
END FUNCTION UP_MX
!
FUNCTION UP_VY(PSRC, PRVCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR ! output src term
END FUNCTION UP_VY
!
FUNCTION UP_MY(PSRC, PRVCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR ! output src term
END FUNCTION UP_MY
!
FUNCTION UP_WZ(PSRC, PRWCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR ! output src term
END FUNCTION UP_WZ
!
FUNCTION UP_MZ(PSRC, PRWCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR ! output src term
END FUNCTION UP_MZ
!
END INTERFACE
!
END MODULE MODI_ADVEC_WENO_K_1_AUX
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION UP_UX(PSRC, PRUCT) RESULT(PR)
!     ########################################################################
!!
!!****  UP_UX - upstream fluxes of U in X direction
!!              input variable PSRC is on U grid, and output PR is on mass grid
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
! upstream flux on mass points
!
PR(IIB:IIE,:,:) = PSRC(IIB:IIE,:,:)     * (0.5+SIGN(0.5,PRUCT(IIB:IIE,:,:))) +&
                  PSRC(IIB+1:IIE+1,:,:) * (0.5-SIGN(0.5,PRUCT(IIB:IIE,:,:)))
!
PR(IIB-1,:,:) = PR(IIE,:,:)
PR(IIE+1,:,:) = PR(IIB,:,:)
!
PR = PR * PRUCT
!
END FUNCTION UP_UX
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION UP_MX(PSRC, PRUCT) RESULT(PR)
!     ########################################################################
!!
!!****  UP_MX - upstream fluxes of variable in X direction
!!              input variable PSRC is on MASS grid, and output PR is on U grid
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS GRID at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on U GRID
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
! upstream flux on mass points
!
PR(IIB:IIE,:,:) = PSRC(IIB-1:IIE-1,:,:) * (0.5 + SIGN(0.5,PRUCT(IIB:IIE,:,:))) &
                  + PSRC(IIB:IIE,:,:)   * (0.5 - SIGN(0.5,PRUCT(IIB:IIE,:,:)))
!
PR(IIB-1,:,:) = PR(IIE,:,:)
PR(IIE+1,:,:) = PR(IIB,:,:)
!
PR = PR * PRUCT
!
END FUNCTION UP_MX
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION UP_VY(PSRC, PRVCT) RESULT(PR)
!     ########################################################################
!!
!!****  UP_VY - upstream fluxes of V in Y direction
!!              input variable PSRC is on V grid, and output PR is on MASS grid
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on V grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
! upstream flux on mass points
!
PR(:,IJB:IJE,:) = PSRC(:,IJB:IJE,:)     * (0.5+SIGN(0.5,PRVCT(:,IJB:IJE,:))) +&
                  PSRC(:,IJB+1:IJE+1,:) * (0.5-SIGN(0.5,PRVCT(:,IJB:IJE,:)))
!
PR(:,IJB-1,:) = PR(:,IJE,:)
PR(:,IJE+1,:) = PR(:,IJB,:)
!
PR = PR * PRVCT
!
END FUNCTION UP_VY
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION UP_MY(PSRC, PRVCT) RESULT(PR)
!     ########################################################################
!!
!!****  UP_MY - upstream fluxes of variable in Y direction
!!              input variable PSRC is on MASS grid, and output PR is on V grid
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on V GRID
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
! upstream flux on mass points
!
PR(:,IJB:IJE,:) = PSRC(:,IJB-1:IJE-1,:) * (0.5+SIGN(0.5,PRVCT(:,IJB:IJE,:))) +&
                  PSRC(:,IJB:IJE,:)     * (0.5-SIGN(0.5,PRVCT(:,IJB:IJE,:)))
!
PR(:,IJB-1,:) = PR(:,IJE,:)
PR(:,IJE+1,:) = PR(:,IJB,:)
!
PR = PR * PRVCT
!
END FUNCTION UP_MY
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION UP_WZ(PSRC, PRWCT) RESULT(PR)
!     ########################################################################
!!
!!****  UP_WZ - upstream fluxes of W in Z direction
!!              input variable PSRC is on W grid, and output PR is on MASS grid
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_PARAMETERS,ONLY: JPVEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB    ! Begining useful area in x,y,z directions
INTEGER :: IKE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
!
! upstream flux on mass points
!
PR(:,:,IKB:IKE) = PSRC(:,:,IKB:IKE)     * (0.5+SIGN(0.5,PRWCT(:,:,IKB:IKE))) +&
                  PSRC(:,:,IKB+1:IKE+1) * (0.5-SIGN(0.5,PRWCT(:,:,IKB:IKE)))
!
PR(:,:,IKB-1) = PR(:,:,IKB)
PR(:,:,IKE+1) = PR(:,:,IKE)
!
PR = PR * PRWCT
!
END FUNCTION UP_WZ
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION UP_MZ(PSRC, PRWCT) RESULT(PR)
!     ########################################################################
!!
!!****  UP_MZ - upstream fluxes of variable in Z direction
!!              input variable PSRC is on MASS grid, and output PR is on W grid
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_PARAMETERS,ONLY: JPVEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB    ! Begining useful area in x,y,z directions
INTEGER :: IKE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
!
! upstream flux on mass points
!
PR(:,:,IKB:IKE) = PSRC(:,:,IKB-1:IKE-1) * (0.5+SIGN(0.5,PRWCT(:,:,IKB:IKE))) +&
                  PSRC(:,:,IKB:IKE)     * (0.5-SIGN(0.5,PRWCT(:,:,IKB:IKE)))
!
PR(:,:,IKB-1) = PR(:,:,IKB)
PR(:,:,IKE+1) = PR(:,:,IKE)
!
PR = PR * PRWCT
!
END FUNCTION UP_MZ
