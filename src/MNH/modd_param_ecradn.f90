!MNH_LIC Copyright 2017-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################## 
      MODULE MODD_PARAM_ECRAD_n
!     ########################
!
!!****  *MODD_PARAM_ECRAD_n* - declaration of the control parameters for
!!                           calling the ECRAD radiation scheme
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define the set of space
!!    and time control parameters for the ECRAD radiation computations.
!!      Equivalent to YOERAD type YRERAD of IFS
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!     Q. Libois   *CNRM-GMME*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      29/05/2017 add ECRAD parameters as namelist
!!      Q. Libois
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------

USE MODD_PARAMETERS, ONLY: JPMODELMAX
USE PARKIND1 , ONLY : JPIM,JPRB
#ifdef MNH_ECRAD
USE radiation_config, ONLY : config_type
#endif
IMPLICIT NONE

TYPE PARAM_ECRAD_t
  ! Parameters entered in namelist NAM_PARAM_ECRAD
  INTEGER(KIND=JPIM) :: NSWSOLVER        ! SW solver ; 0: 'McICA ; 1: 'SPARTACUS' 2: 'SPARTACUS' with 3D effects    
  INTEGER(KIND=JPIM) :: NLWSOLVER        ! LW solver ; 0: 'McICA ; 1: 'SPARTACUS' 2: 'SPARTACUS' with 3D effects   
  INTEGER(KIND=JPIM) :: NLIQOPT          ! optical properties of liquid particles ; 2 = Slingo ; 3 = SOCRATES
  INTEGER(KIND=JPIM) :: NICEOPT          ! optical properties of liquid particles ; 3 = Fu ; 4 = Baran
  INTEGER(KIND=JPIM) :: NGAS             ! gas model 1 = 'Monochromatic', 2 = 'RRTMG-IFS', 3 = 'RRTMG-PSRAD'
  INTEGER(KIND=JPIM) :: NOVLP            ! overlap assumption ; 0= 'Max-Ran' ; 1= 'Exp-Ran'; 2 = 'Exp-Exp'
  INTEGER(KIND=JPIM) :: NREG             ! Number of regions for Triple Clouds
  INTEGER(KIND=JPIM) :: NRADLP           ! 0: ERA-15, 1: Zhang and Rossow, 2: Martin (1994), 3: Martin (1994) and Woods (2000)
  INTEGER(KIND=JPIM) :: NRADIP           ! 0: 40 mum, 1: Liou and Ou (1994), 2: Liou and Ou (1994) improved,
                                         ! 3: Sun and Rikus (1999)
  REAL(KIND=JPRB)    :: XCLOUD_FRAC_STD  ! Cloud water content horizontal fractional standard deviation in a gridbox  
  INTEGER(KIND=JPIM) :: NLWSCATTERING    ! 0: No longwave scattering
                                         ! 1: Longwave scattering by clouds only
                                         ! 2: Longwave scattering by clouds and aerosols
                                
  ! Parameters initialized in ini_radiations_ecrad - present in yoerad of IFS
  LOGICAL            :: LCCNO            ! If True CCN over ocean is diagnosed, otherwise default value
  LOGICAL            :: LCCNL            ! If True CCN over land is diagnosed, otherwise default value
  REAL(KIND=JPRB)    :: XCCNSEA          ! default CCN over sea
  REAL(KIND=JPRB)    :: XCCNLND          ! default CCN over land
  
  ! Aerosols
  INTEGER(KIND=JPIM) :: NAERMACC         ! Use of Tegen (0) or MACC (1) aerosol classification
                                         ! (that is over 6 or 12 aerosols species)
  INTEGER(KIND=JPIM) :: NMCLAT           ! number of latitudes in MACC climatology
  INTEGER(KIND=JPIM) :: NMCLON           ! number of longitudes in MACC climatology
  INTEGER(KIND=JPIM) :: NMCLEV           ! number of vertical levels in MACC climatology              
  INTEGER(KIND=JPIM) :: NMCVAR           ! number of aerosol species in MACC
  REAL(KIND=JPRB)    :: XEXPBC            ! exponential height of BC (m)  
  REAL(KIND=JPRB)    :: XEXPSF            ! exponential height of sulfates (m) 
  REAL(KIND=JPRB)    :: XEXPSS            ! exponential height of sea salts (m) 
  REAL(KIND=JPRB)    :: XEXPDS            ! exponential height of dust (m)  
  REAL(KIND=JPRB)    :: XEXPOM            ! exponential height of organic matter (m)  
  
  INTEGER(KIND=JPIM) :: NMINICE          ! variation of RMINICE with latitude (1), or not (0)
  REAL(KIND=JPRB)    :: XRE2DE           ! CONVERSION FACTOR BETWWEN EFFECTIVE RADIUS AND PARTICLE SIZE
  REAL(KIND=JPRB)    :: XRMINICE         ! MINIMUM SIZE FOR ICE PARTICLES (um) FOR ICE
  INTEGER(KIND=JPIM) :: NDECOLAT         ! Latitude dependence of cloud correlation length
                                         ! 0: SPECIFIED INDEPENDENT OF LATITUDE, 1: SHONK-HOGAN, 2: IMPROVED
  REAL(KIND=JPRB)    :: XDECORR_CF       ! Decorrelation length (km) for cloud fraction
  REAL(KIND=JPRB)    :: XDECORR_CW       ! Decorrelation length (km) for cloud water content
  INTEGER(KIND=JPIM) :: NRADFR           ! radiation every "NRADFR" time steps
  LOGICAL            :: LCentredTimeSZA  ! Compute solar zenith angle in radiation scheme half way between
                                         ! calls to radiation scheme (rather than previous behaviour, 
                                         ! which is half way between calls plus half a model timestep)      
  LOGICAL            :: LAverageSZA      ! Compute an averaged solar zenith angle across the time interval
                                         ! required (either a model timestep or a radiationtimestep)
                                         ! Should be used with LCentredTimeSZA=TRUE.
  ! Trace gases constant mixing ratios
  REAL(KIND=JPRB)    :: XCCH4
  REAL(KIND=JPRB)    :: XCN2O
  REAL(KIND=JPRB)    :: XCNO2
  REAL(KIND=JPRB)    :: XCCFC11
  REAL(KIND=JPRB)    :: XCCFC12
  REAL(KIND=JPRB)    :: XCCFC22
  REAL(KIND=JPRB)    :: XCCCL4   
  
  ! Use of RRTM for SW and LW
  LOGICAL :: LRRTM   
  LOGICAL :: LSRTM  
  
  INTEGER(KIND=JPIM) :: NREDGSW          ! Number of gpoints of SW RRTM 
                                         ! 0 full resolution for RRTM_SW (224)
                                         ! 1 ECMWF High resolution model configuration (_SW: 112)
                                         ! 2 ECMWF EPS configuration (_SW: 56)
  INTEGER(KIND=JPIM) :: NREDGLW          ! Number of gpoints of LW RRTM
                                         ! 0 full resolution for RRTM_LW (256)
                                         ! 1 ECMWF High resolution model configuration (_LW: 140)
                                         ! 2 ECMWF EPS configuration (_LW: 70)
                                         
  LOGICAL :: LAPPROXSWUPDATE !Update the shortwave upwelling flux
                             !every gridpoint to account for the local
                             !value of surface albedo
  LOGICAL :: LAPPROXLWUPDATE ! Update the longwave upwelling flux every
                             !timestep/gridpoint using the stored rate
                             !of change of the fluxes with respect to
                             !the surface upwelling longwave flux                                       

  ! Parameters initialized in ini_radiations_ecrad - not present in yoerad of IFS
  CHARACTER (LEN=255):: CDATADIR! RRTM data directory

  
!   INTEGER(KIND=JPIM) :: NSW    ! number of SW spectral intervals for albedo
!   INTEGER(KIND=JPIM) :: NSW_EC ! number of SW spectral intervals for SRTM  
!  LOGICAL :: LPSRAD
!  LOGICAL :: LSW_ML_E  
!  LOGICAL :: LLW_ML_E
!  LOGICAL :: LEFF3D              ! accounting for 3D effects with Spartacus?
!  LOGICAL :: LSIDEM              ! Side effects
  INTEGER(KIND=JPIM) :: NACTAERO
  INTEGER(KIND=JPIM) :: NAERCLD ! INDEX TO CONTROL SWITCHES FOR AEROSOL-MICROPHYSICS INTERACTION
#ifdef MNH_ECRAD
  type(config_type) :: rad_config
#endif
!-------------------------------------------------------------------------------
!
END TYPE PARAM_ECRAD_t

TYPE(PARAM_ECRAD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PARAM_ECRAD_MODEL

! PRINT*,"MODD_PARAM_ECRAD_n",NSWSOLVER
INTEGER(KIND=JPIM), POINTER :: NSWSOLVER=>NULL()
INTEGER(KIND=JPIM), POINTER :: NLWSOLVER=>NULL()
INTEGER(KIND=JPIM), POINTER :: NLIQOPT=>NULL()
INTEGER(KIND=JPIM), POINTER :: NICEOPT=>NULL()
INTEGER(KIND=JPIM), POINTER :: NGAS=>NULL()
INTEGER(KIND=JPIM), POINTER :: NOVLP=>NULL()
INTEGER(KIND=JPIM), POINTER :: NREG=>NULL()
!LOGICAL, POINTER :: LEFF3D=>NULL()
!LOGICAL, POINTER :: LSIDEM=>NULL()
!LOGICAL, POINTER :: LLWCSCA=>NULL()
!LOGICAL, POINTER :: LLWASCA=>NULL()
LOGICAL, POINTER :: LPSRAD=>NULL()
!LOGICAL, POINTER :: LSW_ML_E=>NULL()
!LOGICAL, POINTER :: LLW_ML_E=>NULL()
INTEGER(KIND=JPIM), POINTER :: NAERMACC=>NULL()
REAL(KIND=JPRB), POINTER :: XDECORR_CF=>NULL() 
REAL(KIND=JPRB), POINTER :: XDECORR_CW=>NULL() 
INTEGER(KIND=JPIM), POINTER :: NRADFR=>NULL()
LOGICAL, POINTER :: LCentredTimeSZA=>NULL()
LOGICAL, POINTER :: LAverageSZA=>NULL()
INTEGER(KIND=JPIM), POINTER :: NLWSCATTERING=>NULL()
INTEGER(KIND=JPIM), POINTER :: NRADLP=>NULL()
INTEGER(KIND=JPIM), POINTER :: NRADIP=>NULL()
INTEGER(KIND=JPIM), POINTER :: NACTAERO=>NULL()
INTEGER(KIND=JPIM), POINTER :: NMCVAR=>NULL()
INTEGER(KIND=JPIM), POINTER :: NMCLAT=>NULL()
INTEGER(KIND=JPIM), POINTER :: NMCLON=>NULL()
INTEGER(KIND=JPIM), POINTER :: NMCLEV=>NULL()
REAL(KIND=JPRB), POINTER    :: XEXPBC=>NULL()
REAL(KIND=JPRB), POINTER    :: XEXPSF=>NULL()
REAL(KIND=JPRB), POINTER    :: XEXPDS=>NULL()
REAL(KIND=JPRB), POINTER    :: XEXPSS=>NULL()
REAL(KIND=JPRB), POINTER    :: XEXPOM=>NULL()
INTEGER(KIND=JPIM), POINTER :: NAERCLD=>NULL() 
LOGICAL, POINTER :: LCCNO=>NULL() 
LOGICAL, POINTER :: LCCNL=>NULL()
REAL(KIND=JPRB), POINTER :: XCCNSEA=>NULL()  
REAL(KIND=JPRB), POINTER :: XCCNLND=>NULL()  
REAL(KIND=JPRB), POINTER :: XRE2DE=>NULL()   
REAL(KIND=JPRB), POINTER :: XRMINICE=>NULL()
INTEGER(KIND=JPIM), POINTER :: NMINICE=>NULL() 
INTEGER(KIND=JPIM), POINTER :: NDECOLAT=>NULL()
REAL(KIND=JPRB), POINTER :: XCLOUD_FRAC_STD=>NULL()
!INTEGER, POINTER :: NSW=>NULL() 
!INTEGER, POINTER :: NSW_EC=>NULL() 
REAL(KIND=JPRB), POINTER :: XCCH4=>NULL()
REAL(KIND=JPRB), POINTER :: XCN2O=>NULL()
REAL(KIND=JPRB), POINTER :: XCNO2=>NULL()
REAL(KIND=JPRB), POINTER :: XCCFC11=>NULL()
REAL(KIND=JPRB), POINTER :: XCCFC12=>NULL()
REAL(KIND=JPRB), POINTER :: XCCFC22=>NULL()
REAL(KIND=JPRB), POINTER :: XCCCL4=>NULL()
LOGICAL, POINTER :: LRRTM=>NULL()
LOGICAL, POINTER :: LSRTM=>NULL()
INTEGER(KIND=JPIM), POINTER :: NREDGLW=>NULL()
INTEGER(KIND=JPIM), POINTER :: NREDGSW=>NULL()
LOGICAL, POINTER :: LAPPROXSWUPDATE=>NULL()
LOGICAL, POINTER :: LAPPROXLWUPDATE=>NULL()
CHARACTER (LEN=255), POINTER :: CDATADIR=>NULL()
#ifdef MNH_ECRAD
type(config_type), pointer :: rad_config => NULL()
#endif
CONTAINS

SUBROUTINE PARAM_ECRAD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
! Current model is set to model KTO
    NSWSOLVER=>PARAM_ECRAD_MODEL(KTO)%NSWSOLVER
    NLWSOLVER=>PARAM_ECRAD_MODEL(KTO)%NLWSOLVER
    NLIQOPT=>PARAM_ECRAD_MODEL(KTO)%NLIQOPT
    NICEOPT=>PARAM_ECRAD_MODEL(KTO)%NICEOPT
    NGAS=>PARAM_ECRAD_MODEL(KTO)%NGAS
    NOVLP=>PARAM_ECRAD_MODEL(KTO)%NOVLP
    NREG=>PARAM_ECRAD_MODEL(KTO)%NREG
!    LEFF3D=>PARAM_ECRAD_MODEL(KTO)%LEFF3D
!    LSIDEM=>PARAM_ECRAD_MODEL(KTO)%LSIDEM
!    LLWCSCA=>PARAM_ECRAD_MODEL(KTO)%LLWCSCA
!    LLWASCA=>PARAM_ECRAD_MODEL(KTO)%LLWASCA
!    LPSRAD=>PARAM_ECRAD_MODEL(KTO)%LPSRAD
!    LSW_ML_E=>PARAM_ECRAD_MODEL(KTO)%LSW_ML_E
!    LLW_ML_E=>PARAM_ECRAD_MODEL(KTO)%LLW_ML_E
    NAERMACC=>PARAM_ECRAD_MODEL(KTO)%NAERMACC
    XDECORR_CF=>PARAM_ECRAD_MODEL(KTO)%XDECORR_CF
    XDECORR_CW=>PARAM_ECRAD_MODEL(KTO)%XDECORR_CW
    NRADFR=>PARAM_ECRAD_MODEL(KTO)%NRADFR
    LCentredTimeSZA=>PARAM_ECRAD_MODEL(KTO)%LCentredTimeSZA
    LAverageSZA=>PARAM_ECRAD_MODEL(KTO)%LAverageSZA
    NLWSCATTERING=>PARAM_ECRAD_MODEL(KTO)%NLWSCATTERING
    NRADLP=>PARAM_ECRAD_MODEL(KTO)%NRADLP
    NRADIP=>PARAM_ECRAD_MODEL(KTO)%NRADIP
    NACTAERO=>PARAM_ECRAD_MODEL(KTO)%NACTAERO    
    NMCVAR=>PARAM_ECRAD_MODEL(KTO)%NMCVAR
    NMCLAT=>PARAM_ECRAD_MODEL(KTO)%NMCLAT
    NMCLON=>PARAM_ECRAD_MODEL(KTO)%NMCLON
    NMCLEV=>PARAM_ECRAD_MODEL(KTO)%NMCLEV
    XEXPBC=>PARAM_ECRAD_MODEL(KTO)%XEXPBC
    XEXPSF=>PARAM_ECRAD_MODEL(KTO)%XEXPSF
    XEXPSS=>PARAM_ECRAD_MODEL(KTO)%XEXPSS
    XEXPDS=>PARAM_ECRAD_MODEL(KTO)%XEXPDS
    XEXPOM=>PARAM_ECRAD_MODEL(KTO)%XEXPOM
    NAERCLD=>PARAM_ECRAD_MODEL(KTO)%NAERCLD
    LCCNO=>PARAM_ECRAD_MODEL(KTO)%LCCNO
    LCCNL=>PARAM_ECRAD_MODEL(KTO)%LCCNL
    XCCNSEA=>PARAM_ECRAD_MODEL(KTO)%XCCNSEA
    XCCNLND=>PARAM_ECRAD_MODEL(KTO)%XCCNLND
    XRE2DE=>PARAM_ECRAD_MODEL(KTO)%XRE2DE
    XRMINICE=>PARAM_ECRAD_MODEL(KTO)%XRMINICE
    NMINICE=>PARAM_ECRAD_MODEL(KTO)%NMINICE
    NDECOLAT=>PARAM_ECRAD_MODEL(KTO)%NDECOLAT
    XCLOUD_FRAC_STD=>PARAM_ECRAD_MODEL(KTO)%XCLOUD_FRAC_STD
!     NSW=>PARAM_ECRAD_MODEL(KTO)%NSW
!     NSW_EC=>PARAM_ECRAD_MODEL(KTO)%NSW_EC
    XCCH4=>PARAM_ECRAD_MODEL(KTO)%XCCH4   
    XCN2O=>PARAM_ECRAD_MODEL(KTO)%XCN2O
    XCNO2=>PARAM_ECRAD_MODEL(KTO)%XCNO2   
    XCCFC11=>PARAM_ECRAD_MODEL(KTO)%XCCFC11
    XCCFC12=>PARAM_ECRAD_MODEL(KTO)%XCCFC12  
    XCCFC22=>PARAM_ECRAD_MODEL(KTO)%XCCFC22
    XCCCL4=>PARAM_ECRAD_MODEL(KTO)%XCCCL4
    LRRTM=>PARAM_ECRAD_MODEL(KTO)%LRRTM
    LSRTM=>PARAM_ECRAD_MODEL(KTO)%LSRTM  
    NREDGLW=>PARAM_ECRAD_MODEL(KTO)%NREDGLW
    NREDGSW=>PARAM_ECRAD_MODEL(KTO)%NREDGSW
    LAPPROXSWUPDATE=>PARAM_ECRAD_MODEL(KTO)%LAPPROXSWUPDATE
    LAPPROXLWUPDATE=>PARAM_ECRAD_MODEL(KTO)%LAPPROXLWUPDATE
    CDATADIR=>PARAM_ECRAD_MODEL(KTO)%CDATADIR
#ifdef MNH_ECRAD
    rad_config=>PARAM_ECRAD_MODEL(KTO)%rad_config
#endif
END SUBROUTINE PARAM_ECRAD_GOTO_MODEL

END MODULE MODD_PARAM_ECRAD_n




