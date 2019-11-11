!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_INIT_GROUND_PARAM_n
!     #######################
INTERFACE
      SUBROUTINE INIT_GROUND_PARAM_n(HINIT,KSV,HSV,PCO2,                       &
                               PZENITH,PAZIM,PSW_BANDS,PLW_BANDS,PDIR_ALB,PSCA_ALB,  &
                               PEMIS,PTSRAD                                )
!
CHARACTER(LEN=3),                  INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                           INTENT(IN)  :: KSV       ! number of scalars
CHARACTER(LEN=6), DIMENSION(KSV),  INTENT(INOUT)::HSV       ! name of all scalar variables
REAL,             DIMENSION(:,:),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/kg)
REAL,             DIMENSION(:,:),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(:,:),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clockwise)
REAL,             DIMENSION(:),    INTENT(IN)  :: PSW_BANDS ! middle wavelength of each SW band
REAL,             DIMENSION(:),    INTENT(IN)  :: PLW_BANDS ! middle wavelength of each LW band
REAL,             DIMENSION(:,:,:),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(:,:,:),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(:,:,:),INTENT(OUT) :: PEMIS     ! spectral emissivity
REAL,             DIMENSION(:,:),  INTENT(OUT) :: PTSRAD    ! radiative temperature
!
END SUBROUTINE INIT_GROUND_PARAM_n
!
END INTERFACE
END MODULE MODI_INIT_GROUND_PARAM_n
!
!     #############################################################
      SUBROUTINE INIT_GROUND_PARAM_n(HINIT,KSV,HSV,PCO2,                       &
                               PZENITH,PAZIM,PSW_BANDS,PLW_BANDS,PDIR_ALB,PSCA_ALB,  &
                               PEMIS,PTSRAD                                )
!     #############################################################
!
!!****  *INIT_GROUND_PARAM_n* - call initialization of externalized surface
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!      01/12/03    (D.Gazen) change emissions handling for surf. externalization
!!      Nov.  2010  (J.Escobar) PGI BUG , add SIZE(CSV) to interface
!!  06/2016     (G.Delautier) phasage surfex 8
!!  01/2018      (G.Delautier) SURFEX 8.1
!!                   02/2018 Q.Libois ECRAD
!!  03/2018     (P.Wautelet)   replace ADD_FORECAST_TO_DATE_SURF by DATETIME_CORRECTDATE
!!   2017  V.Vionnet Blow snow
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_DATETIME
USE MODE_IO_ll
USE MODE_FIELD
USE MODE_ll
!
USE MODD_DYN_n,      ONLY : NSTOP, XTSTEP
USE MODD_REF_n,      ONLY : XRHODREF
USE MODD_CH_M9_n,    ONLY : CNAMES
USE MODD_NSV
USE MODD_DUST,       ONLY : CDUSTNAMES
USE MODD_SALT,       ONLY : CSALTNAMES
USE MODD_CH_AEROSOL, ONLY : CAERONAMES
USE MODD_BLOWSNOW
USE MODD_TYPE_DATE,  ONLY : DATE_TIME
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_SURF=>DATE
!
USE MODD_PARAMETERS, ONLY : XUNDEF, JPVEXT
!
USE MODI_INIT_SURF_ATM_N
!
USE MODD_MNH_SURFEX_n 
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=3),                  INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                           INTENT(IN)  :: KSV       ! number of scalars
CHARACTER(LEN=6), DIMENSION(KSV),  INTENT(INOUT)::HSV       ! name of all scalar variables
REAL,             DIMENSION(:,:),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/kg)
REAL,             DIMENSION(:,:),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(:,:),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clockwise)
REAL,             DIMENSION(:),    INTENT(IN)  :: PSW_BANDS ! middle wavelength of each SW band
REAL,             DIMENSION(:),    INTENT(IN)  :: PLW_BANDS ! middle wavelength of each LW band
REAL,             DIMENSION(:,:,:),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(:,:,:),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(:,:,:),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(:,:),  INTENT(OUT) :: PTSRAD    ! radiative temperature
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:),    ALLOCATABLE :: ZCO2      ! CO2 concentration
REAL, DIMENSION(:),    ALLOCATABLE :: ZRHODREF  ! air density at the surface
REAL, DIMENSION(:),    ALLOCATABLE :: ZZENITH   ! solar zenithal angle
REAL, DIMENSION(:),    ALLOCATABLE :: ZAZIM     ! solar azimuthal angle
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZDIR_ALB  ! direct albedo
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSCA_ALB  ! diffuse albedo
REAL, DIMENSION(:),  ALLOCATABLE :: ZEMIS     ! spectral emissivity
REAL, DIMENSION(:),    ALLOCATABLE :: ZTSRAD    ! radiative temperature
REAL, DIMENSION(:),    ALLOCATABLE :: ZTSURF
!
TYPE(DATE_SURF) :: TDATE_END
!
REAL :: ZDURATION
INTEGER :: ISWB  ! number of SW bands
INTEGER :: ILWB  ! number of LW bands
INTEGER :: IIU   ! 1st array size
INTEGER :: IJU   ! 2nd array size
INTEGER :: IIB   ! X array physical boundary
INTEGER :: IJB   ! Y array physical boundary
INTEGER :: IIE   ! X array physical boundary
INTEGER :: IJE   ! Y array physical boundary
INTEGER :: ILU   ! total physical size
INTEGER :: JLAYER! loop index
INTEGER :: ISV
INTEGER :: IID,IRESP
TYPE (DATE_TIME), POINTER :: TZTCUR=>NULL()
TYPE (DATE_TIME)          :: TZDATE
!
CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE :: YSV_SURF ! name of the scalar variables
                                                        ! sent to SURFEX
!-------------------------------------------------------------------------------
!
!
ISWB = SIZE(PSW_BANDS)
ILWB = SIZE(PLW_BANDS)
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
ILU = (IIE-IIB+1) * (IJE-IJB+1)
!
ALLOCATE(ZCO2    (ILU))
ALLOCATE(ZRHODREF(ILU))
ALLOCATE(ZZENITH (ILU))
ALLOCATE(ZAZIM   (ILU))
ALLOCATE(ZDIR_ALB(ILU,ISWB))
ALLOCATE(ZSCA_ALB(ILU,ISWB))
ALLOCATE(ZEMIS   (ILU))
ALLOCATE(ZTSRAD  (ILU))
ALLOCATE(ZTSURF  (ILU))
!-------------------------------------------------------------------------------
!
ZRHODREF= RESHAPE(XRHODREF(IIB:IIE,IJB:IJE,JPVEXT+1), (/ ILU /) )
ZCO2    = RESHAPE(PCO2    (IIB:IIE,IJB:IJE), (/ ILU /) )
ZZENITH = RESHAPE(PZENITH (IIB:IIE,IJB:IJE), (/ ILU /) )
ZAZIM   = RESHAPE(PAZIM   (IIB:IIE,IJB:IJE), (/ ILU /) )
!
DO JLAYER=NSV_DSTBEG,NSV_DSTEND
  HSV(JLAYER) = TRIM(CDUSTNAMES(JLAYER-NSV_DSTBEG+1))
END DO
DO JLAYER=NSV_SLTBEG,NSV_SLTEND
  HSV(JLAYER) = TRIM(CSALTNAMES(JLAYER-NSV_SLTBEG+1))
END DO
DO JLAYER=NSV_CHGSBEG,NSV_CHGSEND
  HSV(JLAYER) = '#'//TRIM(CNAMES(JLAYER-NSV_CHGSBEG+1))
END DO
DO JLAYER=NSV_AERBEG,NSV_AEREND
  HSV(JLAYER) = '@'//TRIM(CAERONAMES(JLAYER-NSV_AERBEG+1))
END DO
!
CALL FIND_FIELD_ID_FROM_MNHNAME('DTCUR',IID,IRESP)
TZTCUR=>TFIELDLIST(IID)%TFIELD_T0D(1)%DATA
!
TZDATE = TZTCUR
TZDATE%TIME = TZDATE%TIME + NSTOP * XTSTEP
CALL DATETIME_CORRECTDATE(TZDATE)
!Done field by field because TYPE(DATE) different in MesoNH and SURFEX
TDATE_END%YEAR  = TZDATE%TDATE%YEAR
TDATE_END%MONTH = TZDATE%TDATE%MONTH
TDATE_END%DAY   = TZDATE%TDATE%DAY
!
DO JLAYER=NSV_SNWBEG,NSV_SNWEND
  HSV(JLAYER) = TRIM(CSNOWNAMES(JLAYER-NSV_SNWBEG+1))
END DO
!
ISV = SIZE(HSV)
IF(LBLOWSNOW) THEN
    ISV = ISV+NBLOWSNOW_2D ! When blowing snow scheme is used
                  ! NBLOWSN0W_2D variables are sent to SURFEX through ZP_SV.
                  ! They refer to the 2D fields advected by MNH including:
                  !             - total number concentration in Canopy
                  !             - total mass concentration in Canopy
                  !             - equivalent concentration in the saltation layer

    ALLOCATE(YSV_SURF(ISV))
    YSV_SURF(1:KSV)          = HSV(:)
    YSV_SURF(NSV+1:ISV) = YPBLOWSNOW_2D(:)                  
ELSE
    ALLOCATE(YSV_SURF(ISV))
    YSV_SURF(:)     = HSV(:)
ENDIF
CALL INIT_SURF_ATM_n(YSURF_CUR,'MESONH',HINIT,.FALSE.,                  &
                     ILU,ISV,SIZE(PSW_BANDS),                           &
                     YSV_SURF,ZCO2,ZRHODREF,                                 &
                     ZZENITH,ZAZIM,PSW_BANDS,ZDIR_ALB,ZSCA_ALB,         &
                     ZEMIS,ZTSRAD,ZTSURF,                               &
                     TZTCUR%TDATE%YEAR, TZTCUR%TDATE%MONTH,             &
                     TZTCUR%TDATE%DAY, TZTCUR%TIME,                     &
                     TDATE_END,'                            ','      ', &
                     'OK'                                       )
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
!
PDIR_ALB(IIB:IIE,IJB:IJE,:) = RESHAPE(ZDIR_ALB, (/ IIE-IIB+1, IJE-IJB+1, ISWB /) )
PSCA_ALB(IIB:IIE,IJB:IJE,:) = RESHAPE(ZSCA_ALB, (/ IIE-IIB+1, IJE-IJB+1, ISWB /) )
DO JLAYER=1,SIZE(PEMIS,3)
  PEMIS   (IIB:IIE,IJB:IJE,JLAYER) = RESHAPE(ZEMIS,    (/ IIE-IIB+1, IJE-IJB+1 /)       )
END DO
PTSRAD  (IIB:IIE,IJB:IJE)   = RESHAPE(ZTSRAD,   (/ IIE-IIB+1, IJE-IJB+1 /)       )
!-------------------------------------------------------------------------------
DEALLOCATE(ZCO2    )
DEALLOCATE(ZRHODREF)
DEALLOCATE(ZZENITH )
DEALLOCATE(ZAZIM   )
DEALLOCATE(ZDIR_ALB)
DEALLOCATE(ZSCA_ALB)
DEALLOCATE(ZEMIS   )
DEALLOCATE(ZTSRAD  )
DEALLOCATE(YSV_SURF  )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_GROUND_PARAM_n
