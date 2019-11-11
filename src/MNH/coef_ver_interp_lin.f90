!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:17
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_COEF_VER_INTERP_LIN
!     ###############################
INTERFACE COEF_VER_INTERP_LIN
!     ##############################################################
      SUBROUTINE COEF_VER_INTERP_LIN3D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ##############################################################
!
! third dimension of the arrays is vertical
!
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                             ! initial grid 
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                             ! target grid 
INTEGER, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
END SUBROUTINE COEF_VER_INTERP_LIN3D
!     ###############################################################
      SUBROUTINE COEF_VER_INTERP_LIN2D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ###############################################################
!
! second dimension of the arrays is vertical
!
REAL,   DIMENSION(:,:),   INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                             ! initial grid 
REAL,   DIMENSION(:,:),   INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                             ! target grid 
INTEGER, DIMENSION(:,:),              OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:,:),              OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
END SUBROUTINE COEF_VER_INTERP_LIN2D
!     ###############################################################
      SUBROUTINE COEF_VER_INTERP_LIN1D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ###############################################################
!
! first dimension of the arrays is vertical
!
REAL,   DIMENSION(:), INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                         ! initial grid 
REAL,   DIMENSION(:), INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                         ! target grid 
INTEGER, DIMENSION(:),              OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:),              OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
END SUBROUTINE COEF_VER_INTERP_LIN1D
!
!
END INTERFACE
END MODULE MODI_COEF_VER_INTERP_LIN
!     #################################
      MODULE MODI_COEF_VER_INTERP_LIN3D
!     #################################
INTERFACE
!     ###############################################################
      SUBROUTINE COEF_VER_INTERP_LIN3D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ###############################################################
!
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                             ! initial grid 
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                             ! target grid 
INTEGER, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
END SUBROUTINE COEF_VER_INTERP_LIN3D
END INTERFACE
END MODULE MODI_COEF_VER_INTERP_LIN3D
!     ###############################################################
      SUBROUTINE COEF_VER_INTERP_LIN3D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ###############################################################
!
!!****  *VER_INTERP_LIN* - vertical linear interpolation
!!
!!    PURPOSE
!!    -------
!     This function computes the interpolation coefficient XCOEFLIN
!     of the level XKLIN of grid PZ1 which is just under the points of
!     grid PZ2 (respectively called hereafter 'initial' and 'target'), 
!     in order to perform linear interpolations between these 2 grids.
!
!     CAUTION:
!     * The interpolation occurs on the WHOLE grid. Therefore, one must
!     only give as argument to this function the inner points of the domain,
!     particularly for the vertical grid, where there is no physical information
!     under the ground or at and over H.
!     * The level numbers must increase from bottom to top.
!!
!!**  METHOD
!!    ------
!!    two extrapolations are possible: with the two or four nearest points.
!!
!!   Interpolation with 2 points:
!!
!!    If there is less than two points on one side, the interpolation is linear.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	
!     V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/07/97
!!                  20/01/98  use explicit arguments
!!      P Jabouille 20/12/02  no extrapolation under the ground
!!      S. Malardel 11/2003   bug of no extrapolation under the ground
!!      V. Masson   10/2003   no extrapolation above top
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_VER_INTERP_LIN
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                             ! initial grid 
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                             ! target grid 
INTEGER, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
LOGICAL                                               :: GLEUG
LOGICAL,DIMENSION(SIZE(PZ1,1),SIZE(PZ1,2),SIZE(PZ1,3)):: GLEVEL
INTEGER                                               :: JK2,JI,JJ
INTEGER,DIMENSION(SIZE(PZ1,1),SIZE(PZ1,2))            :: ILEVEL
INTEGER,DIMENSION(SIZE(PZ1,1),SIZE(PZ1,2))            :: IUNDER
REAL                                                  :: ZEPS ! a small number
!-------------------------------------------------------------------------------
!
ZEPS=1.E-12
!
!
GLEUG = .FALSE.
IF (PRESENT(OLEUG)) THEN
  GLEUG = OLEUG
END IF
!
!*       1.    ALLOCATIONS
!              -----------
!
IF ( .NOT. PRESENT(PCOEFLIN) ) THEN
  IF (ALLOCATED(NKLIN   )) DEALLOCATE(NKLIN   )
  IF (ALLOCATED(XCOEFLIN)) DEALLOCATE(XCOEFLIN)
  !
  ALLOCATE(NKLIN   (SIZE(PZ2,1),SIZE(PZ2,2),SIZE(PZ2,3)))
  ALLOCATE(XCOEFLIN(SIZE(PZ2,1),SIZE(PZ2,2),SIZE(PZ2,3)))
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    LOOP ON THE TARGET VERTICAL GRID
!              --------------------------------
!
DO JK2=1,SIZE(PZ2,3)
!
!-------------------------------------------------------------------------------
!
!*       3.    Determination of the initial level under the target level JK2
!              -------------------------------------------------------------
!
  GLEVEL(:,:,:)=PZ1(:,:,:)<=SPREAD(PZ2(:,:,JK2),3,SIZE(PZ1,3)) *(1.-ZEPS)
  ILEVEL(:,:)  =COUNT(GLEVEL(:,:,:),3)
!
!* linear extrapolation under the ground
  IUNDER=ILEVEL
  ILEVEL(:,:)=MAX(ILEVEL(:,:),1)
!
!* linear extrapolation above the uppest level
  ILEVEL(:,:)=MIN(ILEVEL(:,:),SIZE(PZ1,3)-1)
!
IF ( .NOT. PRESENT(PCOEFLIN) ) THEN
!* save level in module
  NKLIN(:,:,JK2)=ILEVEL(:,:)
  NKLIN(:,:,JK2) = MAX(NKLIN(:,:,JK2), 1+JPVEXT)
ELSE
!* save level in argument
  KKLIN(:,:,JK2)=ILEVEL(:,:)
  KKLIN(:,:,JK2) = MAX(KKLIN(:,:,JK2), 1+JPVEXT)
END IF

!-------------------------------------------------------------------------------
!
!*       4.    Linear interpolation coefficients
!              ---------------------------------
!
IF ( .NOT. PRESENT(PCOEFLIN) ) THEN
  DO JI=1,SIZE(PZ1,1)
    DO JJ=1,SIZE(PZ1,2)
      XCOEFLIN(JI,JJ,JK2)=(PZ2(JI,JJ,JK2)-PZ1(JI,JJ,ILEVEL(JI,JJ)+1))              &
                         /(PZ1(JI,JJ,ILEVEL(JI,JJ))-PZ1(JI,JJ,ILEVEL(JI,JJ)+1)) 
      IF (.NOT. GLEUG) THEN ! no linear extrapolation under the ground
        IF (IUNDER(JI,JJ) < 1 .OR. ILEVEL(JI,JJ) <= JPVEXT) XCOEFLIN(JI,JJ,JK2)=1.
      END IF
      !
      IF (ILEVEL(JI,JJ)==SIZE(PZ1,3)-1) &  ! no extrapolation above
      XCOEFLIN(JI,JJ,JK2)=MAX(XCOEFLIN(JI,JJ,JK2),0.)
    ENDDO
  ENDDO
ELSE
  DO JI=1,SIZE(PZ1,1)
    DO JJ=1,SIZE(PZ1,2)
      PCOEFLIN(JI,JJ,JK2)=(PZ2(JI,JJ,JK2)-PZ1(JI,JJ,ILEVEL(JI,JJ)+1))              &
                         /(PZ1(JI,JJ,ILEVEL(JI,JJ))-PZ1(JI,JJ,ILEVEL(JI,JJ)+1)) 
      IF (.NOT. GLEUG) THEN ! no linear extrapolation under the ground
        IF (IUNDER(JI,JJ) < 1 .OR. ILEVEL(JI,JJ) <= JPVEXT) PCOEFLIN(JI,JJ,JK2)=1.
      END IF
      !
      IF (ILEVEL(JI,JJ)==SIZE(PZ1,3)-1) &  ! no extrapolation above
      PCOEFLIN(JI,JJ,JK2)=MAX(PCOEFLIN(JI,JJ,JK2),0.)
    ENDDO
  ENDDO
END IF
!
!-------------------------------------------------------------------------------
!
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COEF_VER_INTERP_LIN3D
!     ###############################################################
      SUBROUTINE COEF_VER_INTERP_LIN2D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ###############################################################
!
!!****  *VER_INTERP_LIN* - vertical linear interpolation
!!
!!    PURPOSE
!!    -------
!
!!
!!**  METHOD
!!    ------
!!
!!    This routine calls the 3D version ofCOEF_VER_INTERP_LIN after rewritting of
!!    the fields under 3D form.
!!
!!    EXTERNAL
!!    --------
!!
!!    SUBROUTINE COEF_VER_INTERP_LIN3D
!!    module   MODI_COEF_VER_INTERP_LIN3D
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    Book 2
!!
!!    AUTHOR
!!    ------
!!	
!     V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/07/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_COEF_VER_INTERP_LIN3D
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:), INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                           ! initial grid 
REAL,   DIMENSION(:,:), INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                           ! target grid 
INTEGER, DIMENSION(:,:),                OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:,:),                OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL,   DIMENSION(1,SIZE(PZ1,1),SIZE(PZ1,2))     :: ZZ1   ! altitudes of the points of the
!                                                         ! initial grid 
REAL,   DIMENSION(1,SIZE(PZ2,1),SIZE(PZ2,2))     :: ZZ2   ! altitudes of the points of the
!                                                         ! target grid 
INTEGER, ALLOCATABLE ,DIMENSION(:,:,:):: IKLIN 
!
REAL, ALLOCATABLE,  DIMENSION(:,:,:) :: ZCOEFLIN
!-------------------------------------------------------------------------------
!
ZZ1(1,:,:)  =PZ1(:,:)
ZZ2(1,:,:)  =PZ2(:,:)
!
IF ( PRESENT(PCOEFLIN) ) THEN
  ALLOCATE ( IKLIN ( 1,SIZE(KKLIN,1),SIZE(KKLIN,2)) )
  ALLOCATE ( ZCOEFLIN (1,SIZE(PCOEFLIN,1),SIZE(PCOEFLIN,2)) )
  !
  IF ( PRESENT(OLEUG) ) THEN
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:),IKLIN(:,:,:),ZCOEFLIN(:,:,:),OLEUG)
  ELSE
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:),IKLIN(:,:,:),ZCOEFLIN(:,:,:))
  ENDIF
  PCOEFLIN(:,:)=ZCOEFLIN(1,:,:)
  KKLIN(:,:)=IKLIN(1,:,:)
  DEALLOCATE(IKLIN,ZCOEFLIN)
ELSE
  IF ( PRESENT(OLEUG) ) THEN
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:),OLEUG=OLEUG)
  ELSE
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:))
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COEF_VER_INTERP_LIN2D
!     ###############################################################
      SUBROUTINE COEF_VER_INTERP_LIN1D(PZ1,PZ2,KKLIN,PCOEFLIN,OLEUG)
!     ###############################################################
!
!!****  *VER_INTERP_LIN* - vertical linear interpolation
!!
!!    PURPOSE
!!    -------
!
!!
!!**  METHOD
!!    ------
!!
!!    This routine calls the 3D version ofCOEF_VER_INTERP_LIN after rewritting of
!!    the fields under 3D form.
!!
!!    EXTERNAL
!!    --------
!!
!!    SUBROUTINE COEF_VER_INTERP_LIN3D
!!    module   MODI_COEF_VER_INTERP_LIN3D
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    Book 2
!!
!!    AUTHOR
!!    ------
!!	
!     V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/07/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_COEF_VER_INTERP_LIN3D
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:), INTENT(IN) :: PZ1   ! altitudes of the points of the
!                                         ! initial grid 
REAL,   DIMENSION(:), INTENT(IN) :: PZ2   ! altitudes of the points of the
!                                         ! target grid 
INTEGER, DIMENSION(:),              OPTIONAL :: KKLIN ! number of the level
                                              ! of the data to be interpolated
!
REAL,    DIMENSION(:),              OPTIONAL :: PCOEFLIN ! interpolation
                                              ! coefficient
LOGICAL, OPTIONAL, INTENT(IN) ::   OLEUG      ! for linear extrapolation under the ground
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL,   DIMENSION(1,1,SIZE(PZ1))   :: ZZ1   ! altitudes of the points of the
!                                           ! initial grid 
REAL,   DIMENSION(1,1,SIZE(PZ2))   :: ZZ2   ! altitudes of the points of the
!                                           ! target grid
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)   :: IKLIN 
!
REAL, ALLOCATABLE,    DIMENSION(:,:,:)   :: ZCOEFLIN
!
!-------------------------------------------------------------------------------
!
ZZ1(1,1,:)  =PZ1(:)
ZZ2(1,1,:)  =PZ2(:)
IF ( PRESENT(PCOEFLIN) ) THEN
ALLOCATE( IKLIN (1,1,SIZE(KKLIN))   )
ALLOCATE( ZCOEFLIN (1,1,SIZE(PCOEFLIN)) )
  !
  IF ( PRESENT(OLEUG) ) THEN
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:),IKLIN(:,:,:),ZCOEFLIN(:,:,:),OLEUG)
  ELSE
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:),IKLIN(:,:,:),ZCOEFLIN(:,:,:))
  ENDIF
  PCOEFLIN(:)=ZCOEFLIN(1,1,:)
  KKLIN(:)   =IKLIN(1,1,:)
  DEALLOCATE(IKLIN,ZCOEFLIN)
ELSE
  !
  IF ( PRESENT(OLEUG) ) THEN
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:),OLEUG=OLEUG)
  ELSE
    CALL COEF_VER_INTERP_LIN3D(ZZ1(:,:,:),ZZ2(:,:,:))
  ENDIF
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COEF_VER_INTERP_LIN1D
