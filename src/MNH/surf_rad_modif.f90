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
     MODULE MODI_SURF_RAD_MODIF 
!    ##########################
!
INTERFACE 
!
      SUBROUTINE SURF_RAD_MODIF ( PMAP, PXHAT, PYHAT,            &
                  PCOSZEN, PSINZEN, PAZIMSOL,PZS,PZS_XY,         &
                  PDIRFLASWD, PDIRSRFSWD                         )
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PMAP       ! map factor
REAL, DIMENSION(:),       INTENT(IN) :: PXHAT      ! X coordinate
REAL, DIMENSION(:),       INTENT(IN) :: PYHAT      ! Y coordinate
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN    ! COS(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PSINZEN    ! SIN(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PAZIMSOL   ! azimuthal solar angle
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS        ! (resolved) model orography
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS_XY     ! orography at vort. points
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDIRFLASWD ! Downward DIR SW Flux on flat surf
REAL, DIMENSION(:,:,:),   INTENT(OUT):: PDIRSRFSWD ! Downward SuRF. DIRect    SW Flux
!
END SUBROUTINE SURF_RAD_MODIF  
!
END INTERFACE
!
END MODULE MODI_SURF_RAD_MODIF
!
!     ###################################################################
      SUBROUTINE SURF_RAD_MODIF ( PMAP, PXHAT, PYHAT,            &
                  PCOSZEN, PSINZEN, PAZIMSOL,PZS,PZS_XY,         &
                  PDIRFLASWD, PDIRSRFSWD                         )
!     ###################################################################
!
!!****  * SURF_RAD_MODIF * - computes the modifications to the downwards
!!                           radiative fluxes at the surface, due to
!!                           orientation and shape of this surface.
!!
!!    PURPOSE
!!    -------
!!
!!    1) defines a continuous shape of the orography using triangles
!!       (SURF_SOLAR_GEOM)
!!
!!    2) modification of direct SW downwards flux due to the
!!       slope and orientation of the surface (SURF_SOLAR_SLOPES).
!!       The surface characteristics are compared to the azimuthal 
!!       and zenithal solar angles.
!!
!!    3) modification of direct SW by shadowing from other grid points orography.
!!
!!    4) A procedure is added to insure energy conservation after these modifications.
!!
!!       Only the RESOLVED orography is taken into account for these (4) effects.
!!       Therefore, these modifications will have an impact only for fine
!!       resolutions (large resolved slopes).
!!
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
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/02/95 
!!      V. Masson   28/02/00 extract the surface modifications of the
!!                           RADIATIONS routine, and add the subgrid solar
!!                           computations and the resolved shadows.
!!      V. Masson   18/02/02 rewrites the routine to add shadows from
!!                           one grid point to another
!!      V. Masson   03/03/03 moves local computatitions to surface schemes
!!                           and add multiple SW wavelengths
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_SURF_SOLAR_SUM
USE MODI_SURF_SOLAR_SLOPES
USE MODI_SURF_SOLAR_SHADOWS
!
USE MODD_CONF, ONLY : LCARTESIAN
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PMAP       ! map factor
REAL, DIMENSION(:),       INTENT(IN) :: PXHAT      ! X coordinate
REAL, DIMENSION(:),       INTENT(IN) :: PYHAT      ! Y coordinate
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN    ! COS(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PSINZEN    ! SIN(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PAZIMSOL   ! azimuthal solar angle
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS        ! (resolved) model orography
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS_XY     ! orography at vort. points
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDIRFLASWD ! Downward DIR SW Flux on flat surf
REAL, DIMENSION(:,:,:),   INTENT(OUT):: PDIRSRFSWD ! Downward SuRF. DIRect    SW Flux
!
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
REAL, DIMENSION(SIZE(PZS,1),SIZE(PZS,2)) :: ZMAP  ! map factor
REAL, DIMENSION(SIZE(PZS,1),SIZE(PZS,2),SIZE(PDIRFLASWD,3)) :: ZDIRSWD 
                                                      ! down SW on grid mesh
REAL, DIMENSION(SIZE(PZS,1),SIZE(PZS,2),4,SIZE(PDIRFLASWD,3)) :: ZDIRSWDT
!                                                     ! down SW on triangles
!                                                     ! (4 per grid mesh)
!
REAL, DIMENSION(SIZE(PDIRFLASWD,3)) :: ZENERGY1 
! energy received by the surface by direct solar radiation
REAL, DIMENSION(SIZE(PDIRFLASWD,3)) :: ZENERGY2
! before and after modification of radiation by terrain slopes
REAL, DIMENSION(SIZE(PDIRFLASWD,3)) :: ZENERGYP
! idem except taking into account only positive variations of energy
!
INTEGER :: ISWB  ! number of SW spectral bands
INTEGER :: JSWB  ! loop on SW spectral bands
!-------------------------------------------------------------------------------
!
!* initializations
!
IF (LCARTESIAN) THEN
  ZMAP=1.
ELSE
  ZMAP=PMAP
END IF
!
ISWB = SIZE(PDIRFLASWD,3)
!
!-------------------------------------------------------------------------------
!
DO JSWB = 1, ISWB
  CALL SURF_SOLAR_SUM     (PXHAT, PYHAT, ZMAP, PDIRFLASWD(:,:,JSWB), ZENERGY1(JSWB) )
END DO
!
!
!*       2.    Slope direction direct SW effects
!              ---------------------------------
!
CALL SURF_SOLAR_SLOPES  (ZMAP, PXHAT, PYHAT,                 &
                         PCOSZEN, PSINZEN, PAZIMSOL,         &
                         PZS, PZS_XY, PDIRFLASWD, ZDIRSWDT   )

!
!*       3.    RESOLVED shadows for direct solar radiation
!              -------------------------------------------
!
CALL SURF_SOLAR_SHADOWS (ZMAP, PXHAT, PYHAT,           &
                         PCOSZEN, PSINZEN, PAZIMSOL,   &
                         PZS, PZS_XY, ZDIRSWDT, ZDIRSWD)
!
!
!*       4.    Energy conservation
!              -------------------
!
DO JSWB = 1, ISWB
  CALL SURF_SOLAR_SUM(PXHAT, PYHAT, ZMAP,  &
                      ZDIRSWD(:,:,JSWB),   &
                      ZENERGY2(JSWB)       )
  !
  CALL SURF_SOLAR_SUM(PXHAT, PYHAT, ZMAP,                             &
                      MAX(ZDIRSWD(:,:,JSWB)-PDIRFLASWD(:,:,JSWB),0.), &
                      ZENERGYP(JSWB)                                  )
  !
  IF (ZENERGYP(JSWB)>0.) THEN
    PDIRSRFSWD(:,:,JSWB) = ZDIRSWD(:,:,JSWB)                              &
                         + (ZENERGY1(JSWB)-ZENERGY2(JSWB))/ZENERGYP(JSWB) &
                          * MAX(ZDIRSWD(:,:,JSWB)-PDIRFLASWD(:,:,JSWB),0.)
  ELSE
    PDIRSRFSWD(:,:,JSWB) = PDIRFLASWD(:,:,JSWB)
  END IF
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURF_RAD_MODIF
