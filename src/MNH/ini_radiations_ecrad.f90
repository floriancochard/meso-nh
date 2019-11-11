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
!     ##########################
      MODULE MODI_INI_RADIATIONS_ECRAD
!     ##########################
!
INTERFACE
!
    SUBROUTINE INI_RADIATIONS_ECRAD(HINIFILE,HLUOUT,                      &
         PZHAT,PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,              &
         HLW,KDLON,KFLEV,KFLUX,KRAD,KSWB_OLD,HAER,KAER,KSTATM,                &
         PSTATM,PSEA,PTOWN,PBARE,POZON, PAER,PDST_WL, OSUBG_COND                    )
!
USE MODD_TYPE_DATE

CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
CHARACTER (LEN=*),      INTENT(IN) :: HAER       ! aerosol optical thickness climatology
CHARACTER (LEN=4),      INTENT(IN) :: HLW        ! LW scheme used
!
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! height level without orography
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPABST! pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT  !Temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PTSRAD ! surface radiative temperature
REAL, DIMENSION (:,:),  INTENT(IN) :: PSEA   ! sea fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PTOWN  ! town fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PBARE  ! bare soil fraction
REAL, DIMENSION(:,:),   INTENT(IN) :: PLAT, PLON ! arrays of latitude-longitude
!
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
                                                 !  which must be communicated in INIT
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date&time. Ajout PP
!
INTEGER, INTENT(OUT) :: KDLON   ! number of columns where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLEV   ! number of vertical levels where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLUX   ! number of top and ground fluxes in the output
INTEGER, INTENT(OUT) :: KRAD    ! number of satellite radiances to synthesize
INTEGER, INTENT(OUT) :: KAER    ! number of AERosol classes
INTEGER, INTENT(OUT) :: KSWB_OLD    ! number of SW band
INTEGER, INTENT(OUT) :: KSTATM  ! index of the STAndard ATMosphere level just
                                ! above the model top
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PSTATM   ! working standard atmosphere
!
REAL, DIMENSION(:,:,:),   POINTER :: POZON      ! ozone mixing ratio ( from climato.)
REAL, DIMENSION(:,:,:,:), POINTER :: PAER       ! aerosols optical thickness (from climato)
REAL, DIMENSION(:,:,:,:), POINTER :: PDST_WL    ! aerosols optical thickness (from climato)
LOGICAL, INTENT(IN)  ::   OSUBG_COND ! Switch for sub-grid condensation
!
END SUBROUTINE INI_RADIATIONS_ECRAD
!
END INTERFACE
!
END MODULE MODI_INI_RADIATIONS_ECRAD
!
!
!   #######################################################################
    SUBROUTINE INI_RADIATIONS_ECRAD(HINIFILE,HLUOUT,                      &
         PZHAT,PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,              &
         HLW,KDLON,KFLEV,KFLUX,KRAD,KSWB_OLD,HAER,KAER,KSTATM,                &
         PSTATM,PSEA,PTOWN,PBARE,POZON, PAER,PDST_WL, OSUBG_COND            )
!   #######################################################################
!
! INI_RADIATIONS_ECRAD - Initialization of ECRAD code
!
! PURPOSE
! -------
!   Declarations and intialization of local variables
!   Inspired from IFS routine SUECRAD
!   Initializes all variables of MODD_PARAM_ECRAD_n not initialized in radiation_scheme or through namelist
!
! INTERFACE
! ---------
!    INI_RADIATIONS_ECRAD is called from INI_MODELN. 
!
! AUTHOR
! ------
!   Quentin Libois, CNRM
!   Original: 2018-02-16
!
! MODIFICATIONS
! -------------
!
! TO DO

!
!*       0.    DECLARATIONS
!              ------------
!ECMWF and ECRAD radiation scheme specific modules 
#ifdef MNH_ECRAD
USE RADIATION_SETUP, ONLY : SETUP_RADIATION_SCHEME
USE MODI_SURDI
USE MODI_SURRTAB
USE MODI_SURRTPK
USE MODI_SURRTRF
!USE MODI_RRTM_INIT_140GP
!USE MODI_SRTM_INIT
USE MODD_PARAM_ECRAD_n
USE YOERDI
USE MODD_LUNIT_n , ONLY : TLUOUT
USE YOMLUN    ,ONLY : NULOUT
#endif
USE MODD_PARAM_RAD_n, ONLY : XDTRAD
USE MODD_DYN_n, ONLY : XTSTEP
USE PARKIND1 , ONLY : JPRB

USE MODI_INI_RADIATIONS_ECMWF
USE MODD_TYPE_DATE

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models                                               
CHARACTER (LEN=*),      INTENT(IN) :: HAER       ! aerosol optical thickness climatology
CHARACTER (LEN=4),      INTENT(IN) :: HLW        ! LW scheme used
!
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! height level without orography
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPABST! pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT  !Temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PTSRAD ! surface radiative temperature
REAL, DIMENSION (:,:),  INTENT(IN) :: PSEA   ! sea fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PTOWN  ! town fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PBARE  ! bare soil fraction
REAL, DIMENSION(:,:),   INTENT(IN) :: PLAT, PLON ! arrays of latitude-longitude
!
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date&time. Ajout PP
!
INTEGER, INTENT(OUT) :: KDLON   ! number of columns where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLEV   ! number of vertical levels where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLUX   ! number of top and ground fluxes in the output
INTEGER, INTENT(OUT) :: KRAD    ! number of satellite radiances to synthesize
INTEGER, INTENT(OUT) :: KAER    ! number of AERosol classes
INTEGER, INTENT(OUT) :: KSWB_OLD    ! number of SW band
INTEGER, INTENT(OUT) :: KSTATM  ! index of the STAndard ATMosphere level just
                                ! above the model top
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PSTATM   ! working standard atmosphere
!
REAL, DIMENSION(:,:,:),   POINTER :: POZON      ! ozone mixing ratio ( from climato.)
REAL, DIMENSION(:,:,:,:), POINTER :: PAER       ! aerosols optical thickness (from climato)
REAL, DIMENSION(:,:,:,:), POINTER :: PDST_WL    ! aerosols optical thickness (from climato)
LOGICAL, INTENT(IN)  ::   OSUBG_COND ! Switch for sub-grid condensation

#ifdef MNH_ECRAD

! Redirect OUTPUT to MNH OUTPUT_LISTING
NULOUT = TLUOUT%NLU

! Initialization of ECMWF still neede because many things intialized through this routine
CALL INI_RADIATIONS_ECMWF (HINIFILE,HLUOUT,                                           &
                               PZHAT,PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,    &
                               HLW,KDLON,KFLEV,KFLUX,KRAD,KSWB_OLD,HAER,KAER,KSTATM,       &
                               PSTATM,PSEA,PTOWN,PBARE,POZON, PAER,PDST_WL, OSUBG_COND )

! ECRAD specific variables

LCCNL = .FALSE.       ! True if CCN over land diagnosed
LCCNO = .FALSE.       ! True if CCN over sea is diagnosed

! Constant cloud condensation nuclei over land and sea
! In ECMWF original code, those values were 900 and 150
XCCNLND = 900_JPRB     ! constant CCN over land in m-3 (needed for Martin et al., 1994 parameterization)
XCCNSEA = 50_JPRB      ! constant CCN over sea in m-3

! NAERMACC is in the namelist
! NAERMACC = 0  -> Use of Tegen aerosol climatology
! NAERMACC = 1  -> Use of MACC aerosol classification
                
IF (NAERMACC==0) THEN ! number of aerosol species in MACC climatology
    NMCVAR = 6
ELSE    
    NMCVAR = 12   
END IF 

! MACC Aerosol climatology parameters
NMCLAT = 61
NMCLON = 120
NMCLEV = 1 ! because optical depth only

! Characteristic height for aerosols vertical distribution
XEXPBC = 1000 ! Black Carbon
XEXPSF = 4000 ! Sulfates
XEXPSS = 1000 ! Sea Salts
XEXPDS = 3000 ! Dust
XEXPOM = 2000 ! Organic Matter

NMINICE = 1             ! latitude dependence of RMINICE (0 otherwise)
XRE2DE = 0.64952_JPRB   ! Conversion factor between effective radius and particle size for ice
XRMINICE = 60._JPRB     ! Minimum ice cloud radius

! Decorrelation lengths for cloud fraction and liquid water content
NDECOLAT = 2            ! "Improved" decorrelation length of cloud fraction
XDECORR_CF = 2.0_JPRB
XDECORR_CW = 1.0_JPRB

! Constant mixing ratios for trace gases
XCCH4   = 1.72E-06_JPRB
XCN2O   = 310.E-09_JPRB
XCNO2   = 500.E-13_JPRB
XCCFC11 = 280.E-12_JPRB
XCCFC12 = 484.E-12_JPRB
XCCFC22 =   0.E-12_JPRB
XCCCL4  =   0.E-12_JPRB
! initialisation YOERDI
RCCH4   = 1.72E-06_JPRB
RCN2O   = 310.E-09_JPRB
RCNO2   = 500.E-13_JPRB
RCCFC11 = 280.E-12_JPRB
RCCFC12 = 484.E-12_JPRB
RCCFC22 =   0.E-12_JPRB
RCCCL4  =   0.E-12_JPRB

! Radiation computed every NRADFR timesteps
NRADFR = INT(XDTRAD/XTSTEP)

LCentredTimeSZA = .TRUE. ! compute SZA as centered between radiation time steps
LAverageSZA = .TRUE. ! compute SZA as average between  radiation time steps.
                     ! If True, LCentredTimeSZA must be True as well

!- RRTM as LW and SW scheme
LRRTM = .TRUE.
LSRTM = .TRUE.   

NREDGLW = 1 ! 140 g-points configuration of RRTM LW
NREDGSW = 1 ! 112 g-points configuration of RRTM SW

LAPPROXSWUPDATE  = .TRUE.
LAPPROXLWUPDATE = .TRUE.  

! Directory where to read ecRad databases
CDATADIR = "."

! Load modules from ifsrrtm
CALL SURDI
CALL SURRTAB
CALL SURRTPK
CALL SURRTRF
IF (LRRTM) THEN
    CALL RRTM_INIT_140GP(trim(CDATADIR)) ! LW
ENDIF    
IF (LSRTM) THEN 
    CALL SRTM_INIT(trim(CDATADIR))       ! SW
ENDIF

! Use the above parameters to prepare ecRad setup
CALL SETUP_RADIATION_SCHEME()

#else
PRINT *, "ECRAD LIBRARY NOT AVAILABLE = ###INI_RADIATIONS_ECRAD###"
#endif

END SUBROUTINE INI_RADIATIONS_ECRAD
