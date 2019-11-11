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
!    #############################
     MODULE MODI_SURF_SOLAR_SLOPES 
!    #############################
!
INTERFACE 
!
      SUBROUTINE SURF_SOLAR_SLOPES ( PMAP, PXHAT, PYHAT,                      &
                  PCOSZEN, PSINZEN, PAZIMSOL,                                 &
                  PZS, PZS_XY, PDIRSRFSWD, PDIRSWDT                           )
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PMAP         ! map factor
REAL, DIMENSION(:),       INTENT(IN) :: PXHAT        ! X coordinate
REAL, DIMENSION(:),       INTENT(IN) :: PYHAT        ! Y coordinate
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN ! COS(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PSINZEN ! SIN(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PAZIMSOL! azimuthal solar angle
!                                               ! (radian, from North, clockwise)
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS    ! (resolved) model orography
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS_XY ! model orography at vort. points
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDIRSRFSWD!Downward SuRF. DIRect SW Flux
REAL, DIMENSION(:,:,:,:), INTENT(OUT):: PDIRSWDT ! shortwave flux received by 
!                                                ! each subgrid triangle
!
END SUBROUTINE SURF_SOLAR_SLOPES  
!
END INTERFACE
!
END MODULE MODI_SURF_SOLAR_SLOPES
!     #########################################################################
      SUBROUTINE SURF_SOLAR_SLOPES ( PMAP, PXHAT, PYHAT,                      &
                  PCOSZEN, PSINZEN, PAZIMSOL,                                 &
                  PZS, PZS_XY, PDIRSRFSWD, PDIRSWDT                           )
!     #########################################################################
!
!!****  * SURF_SOLAR_SLOPES * - computes the modifications to the downwards
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
!!      V. Masson   01/03/03 add multiple wavelengths
!!      V. Masson   04/01/11 standard definition of azimuthal angle
!!      J.P. Chaboureau & Juan 21/08/2017 correction for tiny solar zenithal angle in R*4
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : XPI, XMNH_TINY
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PMAP         ! map factor
REAL, DIMENSION(:),       INTENT(IN) :: PXHAT        ! X coordinate
REAL, DIMENSION(:),       INTENT(IN) :: PYHAT        ! Y coordinate
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN ! COS(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PSINZEN ! SIN(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PAZIMSOL! azimuthal solar angle
!                                               ! (radian, from North, clockwise)
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS     ! (resolved) model orography
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS_XY ! model orography at vort. points
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDIRSRFSWD!Downward SuRF. DIRect SW Flux
REAL, DIMENSION(:,:,:,:), INTENT(OUT):: PDIRSWDT ! shortwave flux received by 
!                                                ! each subgrid triangle
!
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: JI, JJ
INTEGER :: JT
!
REAL                   :: ZDZSDX   ! slope in X and Y direction
REAL                   :: ZDZSDY   ! of a triangle surface
REAL                   :: ZSLOPAZI ! azimuthal slope angle
REAL                   :: ZSLOPANG ! vertical slope angle
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
PDIRSWDT(:,:,:,:)=0.
!
!-------------------------------------------------------------------------------
!
!*       1.    LOOP ON GRID MESHES
!              -------------------
!
!* discretization of the grid mesh in four triangles
!
DO JT=1,4
!
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
!
!* If zenithal angle greater than Pi/2, sun is down.
!
      IF (PCOSZEN(JI,JJ)<0.) CYCLE
!
!-------------------------------------------------------------------------------
!
!*       2.    MODIFICATION OF RADIATION DUE TO LOCAL SLOPE
!              --------------------------------------------
!
!* slopes in x and y
!
      SELECT CASE (JT)
        CASE (1)
          ZDZSDX=(    2.* PZS   (JI,JJ)                   &
                   - (PZS_XY(JI,JJ)+PZS_XY(JI,JJ+1)) )    &
                 / (PXHAT(JI+1)-PXHAT(JI))  * PMAP(JI,JJ)
          ZDZSDY=(  PZS_XY(JI,JJ+1) - PZS_XY(JI,JJ) )     &
                 / (PYHAT(JJ+1)-PYHAT(JJ))  * PMAP(JI,JJ)
        CASE (2)
           ZDZSDX=(  PZS_XY(JI+1,JJ+1) -PZS_XY(JI,JJ+1))  &
                 / (PXHAT(JI+1)-PXHAT(JI))  * PMAP(JI,JJ)
           ZDZSDY=(  (PZS_XY(JI+1,JJ+1)+PZS_XY(JI,JJ+1))  &
                     - 2.* PZS (JI,JJ) )                  &
                 / (PYHAT(JJ+1)-PYHAT(JJ))  * PMAP(JI,JJ)
        CASE (3)
          ZDZSDX=(  (PZS_XY(JI+1,JJ)+PZS_XY(JI+1,JJ+1))   &
                   - 2.* PZS(JI,JJ)                    )  &
                 / (PXHAT(JI+1)-PXHAT(JI))  * PMAP(JI,JJ)
          ZDZSDY=(  PZS_XY(JI+1,JJ+1) - PZS_XY(JI+1,JJ) ) &
                 / (PYHAT(JJ+1)-PYHAT(JJ))  * PMAP(JI,JJ)
        CASE (4)
           ZDZSDX=(  PZS_XY(JI+1,JJ) - PZS_XY(JI,JJ) )    &
                 / (PXHAT(JI+1)-PXHAT(JI))  * PMAP(JI,JJ)
           ZDZSDY=(  2.* PZS(JI,JJ)                       &
                   - (PZS_XY(JI+1,JJ)+PZS_XY(JI,JJ)) )    &
                 / (PYHAT(JJ+1)-PYHAT(JJ))  * PMAP(JI,JJ)
      END SELECT
!
!* slope angles
!
      ZSLOPANG = ATAN(SQRT(ZDZSDX**2+ZDZSDY**2))
      ZSLOPAZI = - 0.5*XPI - ATAN2( ZDZSDY, ZDZSDX + SIGN(XMNH_TINY,ZDZSDX) )
!
!* modification of radiation received by 1 square meter of surface 
! (of the triangle) because of its orientation relative to the sun
!
      PDIRSWDT(JI,JJ,JT,:) = MAX( 0.0 , PDIRSRFSWD(JI,JJ,:) * ( &
         COS(ZSLOPANG)                                          &
       + SIN(ZSLOPANG) * PSINZEN(JI,JJ) / MAX(PCOSZEN(JI,JJ), XMNH_TINY) &
         *  COS(PAZIMSOL(JI,JJ)-ZSLOPAZI)                     ) &
                                )
!
!* normalizes received radiation by the surface of the triangle to obtain
!  radiation representative of an horizontal surface.
!
        PDIRSWDT(JI,JJ,JT,:) = PDIRSWDT(JI,JJ,JT,:) &
                             * SQRT(1. + ZDZSDX**2 + ZDZSDY**2)
!
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURF_SOLAR_SLOPES

