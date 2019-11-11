!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!
!     #####################
      MODULE MODI_PPM_RHODJ  
!     #####################
!
INTERFACE
!
      SUBROUTINE PPM_RHODJ (HLBCX,HLBCY,                          &
                          PCRU, PCRV, PCRW, PTSTEP, PRHODJ,       &
                          PRHOX1, PRHOX2, PRHOY1, PRHOY2,         &
                          PRHOZ1, PRHOZ2                          )
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRU  ! Contravariants compon.
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRV  ! 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRW  ! 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! density
!
REAL,                     INTENT(IN)    :: PTSTEP ! Single Time step 
! Temporary advected rhodj
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRHOX1,PRHOX2
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRHOY1,PRHOY2
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRHOZ1,PRHOZ2
!
END SUBROUTINE PPM_RHODJ   
!
END INTERFACE
!
END MODULE MODI_PPM_RHODJ
!
!     ######################################################################
      SUBROUTINE PPM_RHODJ (HLBCX,HLBCY,                          &
                          PCRU, PCRV, PCRW, PTSTEP, PRHODJ,       &
                          PRHOX1, PRHOX2, PRHOY1, PRHOY2,         &
                          PRHOZ1, PRHOZ2                          )
!     ######################################################################
!
!!****  *PPM_RHODJ * 
!!
!!    PURPOSE
!!    -------
!!    Calculate the advection of the density RHODJ to pass to the algorithm PPM
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
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original 11.05.2006. T.Maric
!!      C.Lac 04.2011 Splitted from ppm_met.f90 and ppm_scalar.f90
!!                    to limit duplication in the time splitting
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_PPM             
!
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRU  ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRV  !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRW  ! of momentum
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! density
! Temporary advected rhodj
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRHOX1,PRHOX2
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRHOY1,PRHOY2
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRHOZ1,PRHOZ2
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step 
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IGRID ! localisation on the model grid
!
REAL, DIMENSION(SIZE(PCRU,1),SIZE(PCRU,2),SIZE(PCRU,3)) :: ZUNIT
!
!-------------------------------------------------------------------------------
!
!
IGRID = 1
!
ZUNIT = 1.0
PRHOX1 = PPM_S0_X(HLBCX, IGRID, ZUNIT, PCRU, PRHODJ, PTSTEP)
PRHOY1 = PPM_S0_Y(HLBCY, IGRID, ZUNIT, PCRV, PRHOX1, PTSTEP)
PRHOZ1 = PPM_S0_Z(IGRID, ZUNIT, PCRW, PRHOY1, PTSTEP)
PRHOZ2 = PPM_S0_Z(IGRID, ZUNIT, PCRW, PRHODJ, PTSTEP)
PRHOY2 = PPM_S0_Y(HLBCY, IGRID, ZUNIT, PCRV, PRHOZ2, PTSTEP)
PRHOX2 = PPM_S0_X(HLBCX, IGRID, ZUNIT, PCRU, PRHOY2, PTSTEP)
!
!
END SUBROUTINE PPM_RHODJ
