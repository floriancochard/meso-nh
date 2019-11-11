!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 param 2006/05/18 13:07:25
!-----------------------------------------------------------------
!    ##########################
     MODULE MODI_SURF_SOLAR_SUM 
!    ##########################
!
INTERFACE 
!
      SUBROUTINE SURF_SOLAR_SUM ( PXHAT, PYHAT, PMAP, PDIRSWD, PENERGY )
!
!
REAL, DIMENSION(:),   INTENT(IN) :: PXHAT    ! X coordinate
REAL, DIMENSION(:),   INTENT(IN) :: PYHAT    ! Y coordinate
REAL, DIMENSION(:,:), INTENT(IN) :: PMAP     ! map factor
REAL, DIMENSION(:,:), INTENT(IN) :: PDIRSWD  ! direct SW flux on hor. surf.
REAL,                 INTENT(OUT):: PENERGY  ! energy received by the surface
!
!
END SUBROUTINE SURF_SOLAR_SUM  
!
END INTERFACE
!
END MODULE MODI_SURF_SOLAR_SUM
!
!     ##################################################################
      SUBROUTINE SURF_SOLAR_SUM ( PXHAT, PYHAT, PMAP, PDIRSWD, PENERGY )
!     ##################################################################
!
!!****  * SURF_SOLAR_SUM * - computes the sum of energy received by
!!                           the surface from direct solar radiation
!!
!!    PURPOSE
!!    -------
!!
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
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/01/02
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!JUAN
USE MODE_REPRO_SUM
!JUAN
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
!
REAL, DIMENSION(:),   INTENT(IN) :: PXHAT    ! X coordinate
REAL, DIMENSION(:),   INTENT(IN) :: PYHAT    ! Y coordinate
REAL, DIMENSION(:,:), INTENT(IN) :: PMAP     ! map factor
REAL, DIMENSION(:,:), INTENT(IN) :: PDIRSWD  ! direct SW flux on hor. surf.
REAL,                 INTENT(OUT):: PENERGY  ! energy received by the surface
!
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: JI, JJ
INTEGER :: INFO_ll
!JUAN16
!REAL                               :: ZENERGY  
REAL, ALLOCATABLE, DIMENSION (:,:) :: ZENERGY_2D
!JUAN16
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!-------------------------------------------------------------------------------
!
!*       1.    SUM OF ENERGY FOR THIS PROCESSOR
!              --------------------------------
!
!JUAN16
ALLOCATE(ZENERGY_2D(IIB:IIE,IJB:IJE))
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    ZENERGY_2D(JI,JJ) = PDIRSWD(JI,JJ)*(PXHAT(JI+1)-PXHAT(JI)) &
                                      *(PYHAT(JJ+1)-PYHAT(JJ)) &
                                      /PMAP(JI,JJ)**2
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.    SUM WITH OTHER PROCESSORS
!              -------------------------
!
!CALL REDUCESUM_ll(ZENERGY,INFO_ll)
PENERGY = SUM_DD_R2_ll(ZENERGY_2D)
!JUAN16
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURF_SOLAR_SUM
