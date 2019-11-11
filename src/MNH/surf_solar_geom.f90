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
!    ###########################
     MODULE MODI_SURF_SOLAR_GEOM 
!    ###########################
!
INTERFACE 
!
      SUBROUTINE SURF_SOLAR_GEOM ( PZS, PZS_XY )
!
!
REAL, DIMENSION(:,:), INTENT(IN) :: PZS       ! model orography
REAL, DIMENSION(:,:), INTENT(OUT):: PZS_XY    ! model orography at vort. points
!
!
END SUBROUTINE SURF_SOLAR_GEOM  
!
END INTERFACE
!
END MODULE MODI_SURF_SOLAR_GEOM
!     #######################################################
      SUBROUTINE SURF_SOLAR_GEOM ( PZS, PZS_XY)
!     #######################################################
!
!!****  * SURF_SOLAR_GEOM * - computes continuous shape of the orography
!!                           for the modifications to the downwards
!!                           direct solar flux at the surface, due to
!!                           orientation and shape of this surface.
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
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CONF
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
!
REAL, DIMENSION(:,:), INTENT(IN) :: PZS       ! model orography
REAL, DIMENSION(:,:), INTENT(OUT):: PZS_XY    ! model orography at vort. points
!
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: JI, JJ
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF SURFACE SHAPE (using triangles)
!              ---------------------------------------------
!
!* orography of corners of the grid meshes
!
PZS_XY=XUNDEF
!
DO JJ=IJB,IJE+1
  DO JI=IIB,IIE+1
    PZS_XY(JI,JJ) = 0.25 * (PZS(JI-1,JJ-1)  &
                           +PZS(JI  ,JJ-1)  &
                           +PZS(JI-1,JJ  )  &
                           +PZS(JI  ,JJ  ) )
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.    COASTAL EFFECTS
!              ---------------
!
DO JJ=IJB-1,IJE+1
  DO JI=IIB-1,IIE+1
    IF (PZS(JI,JJ)==0.) THEN
      PZS_XY(JI  ,JJ  ) = 0.
      IF (JI/=IIE+1) PZS_XY(JI+1,JJ  ) = 0.
      IF (JJ/=IJE+1) PZS_XY(JI  ,JJ+1) = 0.
      IF (JI/=IIE+1 .AND. JJ/=IJE+1) PZS_XY(JI+1,JJ+1) = 0.
    END IF
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURF_SOLAR_GEOM
