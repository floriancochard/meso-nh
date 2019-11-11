!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!
!     #####################
      MODULE MODI_PPM_MET  
!     #####################
!
INTERFACE
!
      SUBROUTINE PPM_MET (HLBCX,HLBCY, KRR, TPDTCUR,              &
                          PCRU, PCRV, PCRW, PTSTEP,PTSTEP_PPM, PRHODJ, &
                          PRHOX1, PRHOX2, PRHOY1, PRHOY2,         &
                          PRHOZ1, PRHOZ2, PTHT, PTKET, PRT,       &
                          PRTHS, PRTKES, PRRS, HMET_ADV_SCHEME    )
!
USE MODD_TYPE_DATE, ONLY: DATE_TIME
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
CHARACTER (LEN=6),               INTENT(IN) :: HMET_ADV_SCHEME
!
INTEGER,                  INTENT(IN)    :: KRR    ! Number of moist variables
TYPE (DATE_TIME),         INTENT(IN)    :: TPDTCUR ! current date and time
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRU  ! Courant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRV  ! numbers
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRW  ! 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! density
! Temporary advected rhodj
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHOX1,PRHOX2
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHOY1,PRHOY2
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHOZ1,PRHOZ2
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step model  
REAL,                     INTENT(IN)    :: PTSTEP_PPM ! Time Step PPM
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET        ! Vars at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PRTHS, PRTKES! Source terms
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PRRS 
!
END SUBROUTINE PPM_MET   
!
END INTERFACE
!
END MODULE MODI_PPM_MET
!
!     ######################################################################
      SUBROUTINE PPM_MET (HLBCX,HLBCY, KRR, TPDTCUR,              &
                          PCRU, PCRV, PCRW, PTSTEP,PTSTEP_PPM, PRHODJ, &
                          PRHOX1, PRHOX2, PRHOY1, PRHOY2,         &
                          PRHOZ1, PRHOZ2, PTHT, PTKET, PRT,       &
                          PRTHS, PRTKES, PRRS, HMET_ADV_SCHEME    )
!     ######################################################################
!
!!****  *PPM_MET * 
!!
!!    PURPOSE
!!    -------
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
!!    MODULE MODD_ARGSLIST
!!         HALO2LIST_ll : type for a list of "HALO2_lls"
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
!!      Modification : 11.2011 C.Lac, V.Masson : Advection of (theta_l,r_t) 
!!                  10/2016  (C.Lac) Correction on the flag for Strang splitting
!!                                  to insure reproducibility between START and RESTA 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODE_ll
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_TYPE_DATE, ONLY: DATE_TIME
!
USE MODI_SHUMAN
USE MODI_PPM
USE MODI_ADVEC_PPM_ALGO
!
! incorporate ADVEC_4TH_ORDER_ALG, MZF4 and MZM4
!USE MODI_ADVEC_4TH_ORDER_AUX
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
CHARACTER (LEN=6),               INTENT(IN) :: HMET_ADV_SCHEME
!
INTEGER,                  INTENT(IN)    :: KRR    ! Number of moist variables
TYPE (DATE_TIME),         INTENT(IN)    :: TPDTCUR ! current date and time
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRU  ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRV  !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRW  ! of momentum
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! density
! Temporary advected rhodj
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHOX1,PRHOX2
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHOY1,PRHOY2
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHOZ1,PRHOZ2
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step model  
REAL,                     INTENT(IN)    :: PTSTEP_PPM ! Time Step PPM
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET ! Vars at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PRTHS, PRTKES! Source terms
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PRRS 
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR           ! Loop index for  moist variables
!
LOGICAL :: GTKEALLOC     ! true if TKE arrays are not zero-sized
!
INTEGER :: IGRID ! localisation on the model grid
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTES THE DOMAIN DIMENSIONS
!               ------------------------------
!
GTKEALLOC = SIZE(PTKET,1) /= 0
!
!-------------------------------------------------------------------------------
!
!*       2.     CALL THE ADVEC_PPM_ALGO ROUTINE FOR EACH FIELD
!               -----------------------------------------------
!
IGRID = 1
!
!
! Potential temperature
!
CALL ADVEC_PPM_ALGO(HMET_ADV_SCHEME, HLBCX, HLBCY, IGRID, PTHT, PRHODJ, PTSTEP, &
                    PTSTEP_PPM,PRHOX1, PRHOX2, PRHOY1, PRHOY2, PRHOZ1, PRHOZ2, &
                    PRTHS, TPDTCUR, PCRU, PCRV, PCRW)
!
!
! Turbulence variables
!
IF (GTKEALLOC) THEN
   CALL ADVEC_PPM_ALGO(HMET_ADV_SCHEME, HLBCX, HLBCY, IGRID, PTKET,PRHODJ,PTSTEP, &
                       PTSTEP_PPM,PRHOX1, PRHOX2, PRHOY1, PRHOY2, PRHOZ1, PRHOZ2, &
                       PRTKES, TPDTCUR, PCRU, PCRV, PCRW)
!
!
END IF
!
!
!
! Case with KRR moist variables
!
DO JRR=1,KRR
   CALL ADVEC_PPM_ALGO(HMET_ADV_SCHEME, HLBCX, HLBCY, IGRID,           &
                       PRT(:,:,:,JRR), PRHODJ, PTSTEP, PTSTEP_PPM,     &
                       PRHOX1, PRHOX2, PRHOY1, PRHOY2, PRHOZ1, PRHOZ2, &
                       PRRS(:,:,:,JRR), TPDTCUR, PCRU, PCRV, PCRW                 )
END DO
!
!
END SUBROUTINE PPM_MET
