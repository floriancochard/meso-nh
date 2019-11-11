!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_VER_INTERP_LIN
!     ######################
INTERFACE VER_INTERP_LIN
!     ##############################################
      FUNCTION VER_INTERP_LIN3D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
!
! third dimension of the arrays is vertical
!
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PVAR1 ! variable values on the initial
!                                             ! grid
INTEGER,DIMENSION(:,:,:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(SIZE(KKLIN,1),SIZE(KKLIN,2),SIZE(KKLIN,3)) &
                                     :: PVAR2 ! variable values on target
!                                             ! grid 
END FUNCTION VER_INTERP_LIN3D
!     ##############################################
      FUNCTION VER_INTERP_LIN2D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
!
! second dimension of the arrays is vertical
!
REAL,   DIMENSION(:,:),   INTENT(IN) :: PVAR1 ! variable values on the initial
!                                             ! grid
INTEGER,DIMENSION(:,:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:,:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(SIZE(KKLIN,1),SIZE(KKLIN,2))                               &
                                     :: PVAR2 ! variable values on target
!                                             ! grid 
END FUNCTION VER_INTERP_LIN2D
!     ##############################################
      FUNCTION VER_INTERP_LIN1D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
!
! first dimension of the arrays is vertical
!
REAL,   DIMENSION(:), INTENT(IN) :: PVAR1 ! variable values on the initial
!                                         ! grid
INTEGER,DIMENSION(:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(SIZE(KKLIN)) :: PVAR2 ! variable values on target
!                                         ! grid 
END FUNCTION VER_INTERP_LIN1D
!
!
END INTERFACE
END MODULE MODI_VER_INTERP_LIN
!     ############################
      MODULE MODI_VER_INTERP_LIN3D
!     ############################
INTERFACE
!     ##############################################
      FUNCTION VER_INTERP_LIN3D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
!
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PVAR1 ! variable values on the initial
!                                             ! grid
INTEGER,DIMENSION(:,:,:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(SIZE(KKLIN,1),SIZE(KKLIN,2),SIZE(KKLIN,3))                   &
                                     :: PVAR2 ! variable values on target
!                                             ! grid 
END FUNCTION VER_INTERP_LIN3D
END INTERFACE
END MODULE MODI_VER_INTERP_LIN3D
!     ##############################################
      FUNCTION VER_INTERP_LIN3D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
!
!!****  *VER_INTERP_LIN* - vertical linear interpolation
!!
!!    PURPOSE
!!    -------
!     This function interpolates the 3D fields from one grid
!     to another using linear interpolation cofficients stored in module
!     MODD_VER_INTERP_LIN.
!
!!
!!**  METHOD
!!    ------
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
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PVAR1 ! variable values on the initial
!                                             ! grid
INTEGER,DIMENSION(:,:,:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:,:,:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(SIZE(KKLIN,1),SIZE(KKLIN,2),SIZE(KKLIN,3))                   &
                                     :: PVAR2 ! variable values on target
!                                             ! grid 
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER                                               :: JI,JJ,JK2
!-------------------------------------------------------------------------------
!
DO JK2=1,SIZE(KKLIN,3)
  DO JJ=1,SIZE(KKLIN,2)
    DO JI=1,SIZE(KKLIN,1)
      PVAR2(JI,JJ,JK2)=    PCOEFLIN(JI,JJ,JK2) *PVAR1(JI,JJ,KKLIN(JI,JJ,JK2)  )&
                      +(1.-PCOEFLIN(JI,JJ,JK2))*PVAR1(JI,JJ,KKLIN(JI,JJ,JK2)+1)
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END FUNCTION VER_INTERP_LIN3D
!     ##############################################
      FUNCTION VER_INTERP_LIN2D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
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
!!    This routine calls the 3D version of VER_INTERP_LIN after rewritting of
!!    the fields under 3D form.
!!
!!    EXTERNAL
!!    --------
!!
!!    function VER_INTERP_LIN3D
!!    module   MODI_VER_INTERP_LIN3D
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
!!      Original    17/07/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_VER_INTERP_LIN3D
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:), INTENT(IN) :: PVAR1 ! variable values on the initial
!                                           ! grid
INTEGER,DIMENSION(:,:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:,:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(SIZE(KKLIN,1),SIZE(KKLIN,2)) :: PVAR2 ! variable values on
!                                                       ! target grid 
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL,   DIMENSION(1,SIZE(PVAR1,1),SIZE(PVAR1,2)) :: ZVAR1 ! variable values on the initial
!                                                         ! grid
REAL,   DIMENSION(1,SIZE(KKLIN,1),SIZE(KKLIN,2)) :: ZVAR2 ! variable values on target
!
INTEGER,DIMENSION(1,SIZE(KKLIN,1),SIZE(KKLIN,2)) :: IKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(1,SIZE(PCOEFLIN,1),SIZE(PCOEFLIN,2)):: ZCOEFLIN ! coefficient for level KKLIN
!
!-------------------------------------------------------------------------------
!
ZVAR1(1,:,:)=PVAR1(:,:)
IKLIN(1,:,:)=KKLIN(:,:)
ZCOEFLIN(1,:,:)=PCOEFLIN(:,:)
!
ZVAR2(:,:,:)=VER_INTERP_LIN3D(ZVAR1(:,:,:),IKLIN(:,:,:),ZCOEFLIN(:,:,:))
!
PVAR2(:,:)  =ZVAR2(1,:,:)
!
!-------------------------------------------------------------------------------
!
END FUNCTION VER_INTERP_LIN2D
!     ##############################################
      FUNCTION VER_INTERP_LIN1D(PVAR1,KKLIN,PCOEFLIN) RESULT(PVAR2)
!     ##############################################
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
!!    This routine calls the 3D version of VER_INTERP_LIN after rewritting of
!!    the fields under 3D form.
!!
!!    EXTERNAL
!!    --------
!!
!!    function VER_INTERP_LIN3D
!!    module   MODI_VER_INTERP_LIN3D
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
!!      Original    17/07/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_VER_INTERP_LIN3D
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:), INTENT(IN) :: PVAR1 ! variable values on the initial
!                                         ! grid
INTEGER,DIMENSION(:), INTENT(IN) :: KKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(:), INTENT(IN) :: PCOEFLIN ! coefficient for level KKLIN

REAL,   DIMENSION(SIZE(KKLIN)) :: PVAR2 ! variable values on target
!                                         ! grid 
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER,DIMENSION(1,1,SIZE(KKLIN)) :: IKLIN ! lower interpolating level of
!                                             ! grid 1 for each level of grid 2 
REAL,   DIMENSION(1,1,SIZE(PCOEFLIN)) :: ZCOEFLIN ! coefficient for level KKLIN
!
REAL,   DIMENSION(1,1,SIZE(PVAR1)) :: ZVAR1 ! variable values on the initial
!                                           ! grid
REAL,   DIMENSION(1,1,SIZE(KKLIN)) :: ZVAR2 ! variable values on target
!
!-------------------------------------------------------------------------------
!
ZVAR1(1,1,:)=PVAR1(:)
IKLIN(1,1,:)=KKLIN(:)
ZCOEFLIN(1,1,:)=PCOEFLIN(:)
!
ZVAR2(:,:,:)=VER_INTERP_LIN3D(ZVAR1(:,:,:),IKLIN(:,:,:),ZCOEFLIN(:,:,:))
!
PVAR2(:)    =ZVAR2(1,1,:)
!
!-------------------------------------------------------------------------------
!
END FUNCTION VER_INTERP_LIN1D
