!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 diag 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##########################################
      MODULE MODI_COMPUTE_MEAN_PRECIP
!     ##########################################
INTERFACE
      SUBROUTINE COMPUTE_MEAN_PRECIP(PVAR1,PMEAN_PR,PVAR2,KGRID)
REAL, DIMENSION(:,:), INTENT(IN) :: PVAR1    ! input array to be averaged
REAL, DIMENSION(2), INTENT(IN)   :: PMEAN_PR ! number of grid points of the
                                             ! small-scale model inside a
                                             ! large-scale mesh-grid
REAL, DIMENSION(:,:), INTENT(OUT):: PVAR2    ! output array 
INTEGER,              INTENT(OUT):: KGRID    ! localization of grid point
!
END SUBROUTINE COMPUTE_MEAN_PRECIP
END INTERFACE
END MODULE MODI_COMPUTE_MEAN_PRECIP
!
!     ##########################################
      SUBROUTINE COMPUTE_MEAN_PRECIP(PVAR1,PMEAN_PR,PVAR2,KGRID)
!     ##########################################
!
!!****  *COMPUTE_MEAN_PRECIP* - routine to compute the mean of a 2D field in the
!!                          mesh-grid of a large-scale model
!!
!!    PURPOSE
!!    -------
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
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!     P.Hereil N Asencio   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    3/02/98
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN) :: PVAR1    ! input array to be averaged
REAL, DIMENSION(2), INTENT(IN)   :: PMEAN_PR ! number of grid points of the
                                             ! small-scale model inside a
                                             ! large-scale mesh-grid
REAL, DIMENSION(:,:), INTENT(OUT):: PVAR2    ! output array 
INTEGER,              INTENT(OUT):: KGRID    ! localization of grid point
!
!*       0.2   Declarations of local variables
!
INTEGER           :: IIMESH_RATIO,IJMESH_RATIO
INTEGER           :: JI,JJ,JI2,JJ2,IJU,IIU
REAL              :: ZSOM
!-------------------------------------------------------------------------------
!
!
IIU=SIZE(PVAR1,1)
IJU=SIZE(PVAR1,2)
!
PVAR2(:,:)=XUNDEF
IIMESH_RATIO=FLOOR(PMEAN_PR(1)) ! number of grid points of the 
                                ! small-scale model inside a mesh-grid
                                ! of the large-scale model along x
IF (IIMESH_RATIO<1) THEN
  PRINT * ,'YOU WANT TO COMPUTE THE MEAN ACCUMULATED PRECIPITATION BUT'
  PRINT * ,'THE ASPECT RATIO ALONG X IS SMALLER THEN 1'
  IIMESH_RATIO = 1
END IF
!
IJMESH_RATIO=FLOOR(PMEAN_PR(2)) ! number of grid points of the 
                                ! small-scale model inside a mesh-grid
                                ! of the large-scale model along y
IF (IJMESH_RATIO<1) THEN
  PRINT * ,'YOU WANT TO COMPUTE THE MEAN ACCUMULATED PRECIPITATION' 
  PRINT * ,' BUT THE ASPECT RATIO ALONG Y IS SMALLER THEN 1'
  IJMESH_RATIO = 1
END IF
!
! localization of the center point  of the  large-scale mesh-grid
IF ( MOD( IIMESH_RATIO,2)  .EQ. 0) THEN
  IF ( MOD(IJMESH_RATIO,2)  .EQ. 0) THEN
    KGRID=7
  ELSE
    KGRID=2
  END IF
ELSE
  IF ( MOD(IJMESH_RATIO,2)  .EQ. 0) THEN
    KGRID=3
  ELSE
    KGRID=1
  END IF
END IF
!  ratio must be odd
IF ( MOD( IIMESH_RATIO,2)  .EQ. 0) THEN
  IIMESH_RATIO = IIMESH_RATIO +1
END IF
IF ( MOD( IJMESH_RATIO,2)  .EQ. 0) THEN
  IJMESH_RATIO = IJMESH_RATIO +1
END IF
!
DO JJ=2,IJU-IJMESH_RATIO
  DO JI=2,IIU-IIMESH_RATIO
    ZSOM=0.
    DO JJ2=JJ,JJ+IJMESH_RATIO-1
      DO JI2=JI,JI+IIMESH_RATIO-1
        ZSOM=PVAR1(JI2,JJ2)+ZSOM
      END DO
    END DO 
    ZSOM=ZSOM/(IIMESH_RATIO*IJMESH_RATIO) 
    PVAR2(JI+IIMESH_RATIO/2 , JJ+IIMESH_RATIO/2) = ZSOM
  END DO 
END DO
!
END SUBROUTINE COMPUTE_MEAN_PRECIP
