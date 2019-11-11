!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #######################
     MODULE MODI_CALL_RTTOV8
!    #######################
INTERFACE
!
     SUBROUTINE CALL_RTTOV8(KDLON, KFLEV, KSTATM, PEMIS, PTSRAD, PSTATM,   &
                PTHT, PRT, PPABST, PZZ, PMFCONV, PCLDFR, PULVLKB, PVLVLKB,  &
                OUSERI, KRTTOVINFO, TPFILE    )
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
INTEGER, INTENT(IN)   :: KDLON !number of columns where the
                               !radiation calculations are performed
INTEGER, INTENT(IN)   :: KFLEV !number of vertical levels where the
                               !radiation calculations are performed
INTEGER, INTENT(IN)   :: KSTATM  !index of the standard atmosphere level
                                 !just above the model top
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PEMIS  !Surface IR EMISsivity
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD !RADiative Surface Temperature
REAL, DIMENSION(:,:),     INTENT(IN) :: PSTATM !selected standard atmosphere
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT   !THeta at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT    !moist variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST !pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ    !Model level heights
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! cloud fraction
REAL, DIMENSION(:,:),     INTENT(IN) :: PULVLKB ! U-wind at KB level
REAL, DIMENSION(:,:),     INTENT(IN) :: PVLVLKB ! V-wind at KB level
!
LOGICAL, INTENT(IN)                  :: OUSERI ! logical switch to compute both
                                               ! liquid and solid condensate (OUSERI=.TRUE.)
                                               ! or only liquid condensate (OUSERI=.FALSE.)
!
INTEGER, DIMENSION(:,:), INTENT(IN) :: KRTTOVINFO ! platform, satelit, sensor,
                                                  ! and selection calculations
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! File characteristics
!
END SUBROUTINE CALL_RTTOV8
END INTERFACE
END MODULE MODI_CALL_RTTOV8
!    #####################################################################
SUBROUTINE CALL_RTTOV8(KDLON, KFLEV, KSTATM, PEMIS, PTSRAD, PSTATM,     &
           PTHT, PRT, PPABST, PZZ, PMFCONV, PCLDFR, PULVLKB, PVLVLKB,  &
           OUSERI, KRTTOVINFO, TPFILE    )
!    #####################################################################
!!
!!****  *CALL_RTTOV* - 
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
!!    REFERENCE
!!    ---------
!!    See Chaboureau and Pinty, 2006
!!    Validation of a cirrus parameterization with Meteosat Second Generation
!!    observations. Geophys. Res. Let., doi:10.1029/2005GL024725
!!
!!    AUTHOR
!!    ------
!!      J.-P. Chaboureau       *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/12/03
!!      JP Chaboureau 27/03/2008 Vectorization
!!      JP Chaboureau 02/11/2009 move GANGL deallocation outside the sensor loop
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!----------------------------------------------------------------------------
!!
!!*       0.    DECLARATIONS
!!              ------------
!!
USE MODD_CST
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_GRID_n
USE MODD_DEEP_CONVECTION_n
USE MODD_REF_n
USE MODD_RADIATIONS_n,  ONLY : XSEA
!
USE MODN_CONF
!                                
USE MODD_RAD_TRANSF
!                               
USE MODI_DETER_ANGLE
USE MODI_PINTER
!
USE MODE_FIELD
USE MODE_FMWRIT
USE MODE_FMREAD
USE MODE_ll
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
USE MODE_POS
!
#ifdef MNH_RTTOV_8
USE rttov_const, ONLY :  &
     &  gas_id_watervapour  ,&
     & errorstatus_success,&
     & errorstatus_warning,&
     & errorstatus_fatal,&
     & sensor_id_mw,        &
     & npolar_return,       &
     & npolar_compute


USE rttov_types, ONLY : &
     & geometry_type  ,&
     & rttov_coef     ,&
     & rttov_scatt_coef     ,&
     & profile_type   ,&
     & profile_cloud_type   ,&
     & transmission_type    ,&
     & radiance_cloud_type

USE MOD_CPARAM, ONLY : jppf   ! Max no. profiles

USE parkind1, ONLY : jpim     ,jprb
!
IMPLICIT NONE
!
! -----------------------------------------------------------------------------
INTERFACE

!!! #include "rttov_setupindex.interface
SUBROUTINE rttov_setupindex (&
     & mchan,           & ! in
     & nprofiles,       & ! in
     & nfrequencies,    & ! in
     & nchannels,       & ! in
     & nbtout,          & ! in
     & coef,            & ! in
     & surfem,          & ! in
     & lprofiles,       & ! out
     & channels,        & ! out
     & polarisations,   & ! out
     & emissivity)        ! out
! Imported Type Definitions:
USE rttov_types, ONLY : &
     rttov_coef
USE rttov_const, ONLY :   &
     sensor_id_mw,     &
     npolar_return,    &
     npolar_compute

USE parkind1, ONLY : jpim     ,jprb
IMPLICIT NONE
! Subroutine arguments
INTEGER(Kind=jpim),  INTENT(in)    :: nprofiles          ! Number of profiles
INTEGER(Kind=jpim),  INTENT(in)    :: mchan(nprofiles)   ! nfrequencies/nprofiles
INTEGER(Kind=jpim),  INTENT(in)    :: nchannels          ! Number of radiances computed
INTEGER(Kind=jpim),  INTENT(in)    :: nfrequencies       ! Number of frequencies
!  (= channels used * profiles)
INTEGER(Kind=jpim),  INTENT(in)    :: nbtout             ! Number of BTs returned
INTEGER(Kind=jpim),  INTENT(out)   :: channels(nfrequencies)      ! Channel indices
INTEGER(Kind=jpim),  INTENT(out)   :: polarisations(nchannels,3)  ! Channel indices
INTEGER(Kind=jpim),  INTENT(out)   :: lprofiles(nfrequencies)     ! Profiles indices
REAL(Kind=jprb),     INTENT(in)    :: surfem(nchannels)           ! Input surface emissivity
REAL(Kind=jprb),     INTENT(out)   :: emissivity(nchannels)       ! Surface emissivity array for RTTOV
TYPE( rttov_coef ),  INTENT (in)   :: coef               ! coefficients
END SUBROUTINE rttov_setupindex
!
!!! #include "rttov_setupchan.interface"
SUBROUTINE rttov_setupchan (&
     & nprofiles,       & ! in
     & nchan,           & ! in
     & coef,            & ! in
     & nfrequencies,    & ! out
     & nchannels,       & ! out
     & nbtout)            ! out
  ! Imported Type Definitions:
USE rttov_types, ONLY : &
     rttov_coef
USE rttov_const, ONLY :   &
     sensor_id_mw,     &
     npolar_return,    &
     npolar_compute
USE parkind1, ONLY : jpim
IMPLICIT NONE
! Subroutine arguments
INTEGER(Kind=jpim),  INTENT(in)    :: nprofiles        ! Number of profiles
INTEGER(Kind=jpim),  INTENT(in)    :: nchan(nprofiles) ! Number of channels requested
TYPE( rttov_coef ),  INTENT (in)   :: coef             ! coefficients
INTEGER(Kind=jpim),  INTENT(out)   :: nchannels        ! Number of radiances computed
INTEGER(Kind=jpim),  INTENT(out)   :: nfrequencies     ! Number of frequencies
!  (= channels used * profiles)
INTEGER(Kind=jpim),  INTENT(out)   :: nbtout           ! Number of BTs returned
END SUBROUTINE rttov_setupchan
!
!!! #include "rttov_scatt_setupindex.interface"
SUBROUTINE rttov_scatt_setupindex (nprofiles, n_chan, coef, nchannels, &
     & lsprofiles,lsprofiles2, frequencies, nbtout)
USE parkind1     , ONLY: jpim, jprb
USE rttov_const, ONLY : npolar_return, npolar_compute, &
     & inst_id_ssmi
USE rttov_types, ONLY : rttov_coef     
IMPLICIT NONE
INTEGER (kind=jpim), INTENT ( in) :: nprofiles
INTEGER (kind=jpim), INTENT ( in) :: nchannels
INTEGER (kind=jpim), INTENT ( in) :: nbtout
INTEGER (kind=jpim), INTENT ( in) :: n_chan (nprofiles)    
TYPE   (rttov_coef), INTENT ( in) :: coef   
INTEGER (kind=jpim), INTENT (out), DIMENSION (nchannels)    :: lsprofiles   
INTEGER (kind=jpim), INTENT (out), DIMENSION (nbtout)       :: lsprofiles2  
INTEGER (kind=jpim), INTENT (out), DIMENSION (nchannels)    :: frequencies
END SUBROUTINE rttov_scatt_setupindex
!
!!! #include "rttov_cld.interface"
SUBROUTINE rttov_cld( &
     errorstatus,     & ! out
     nfrequencies,    & ! in
     nchannels,       & ! in
     nbtout,          & ! in
     nprofiles,       & ! in
     channels,        & ! in
     polarisations,   & ! in
     lprofiles,       & ! in
     profiles,        & ! inout  (to invalid clw absorption)
     cld_profiles,    & ! in
     coef,            & ! in
     calcemis,        & ! in
     emissivity,      & ! inout
     cld_radiance )     ! inout
USE rttov_const, ONLY :   &
     errorstatus_success ,&
     errorstatus_fatal   ,&
     overlap_scheme
USE rttov_types, ONLY :    &
     rttov_coef           ,&
     geometry_Type        ,&
     profile_Type         ,&
     profile_cloud_Type   ,&
     transmission_Type    ,&
     radiance_Type        ,&
     radiance_cloud_Type
USE parkind1, ONLY : jpim     ,jprb
IMPLICIT NONE
INTEGER(Kind=jpim),        INTENT(in)    :: nbtout  ! Number of output radiances
INTEGER(Kind=jpim),        INTENT(in)    :: nfrequencies  ! Number of output radiances
INTEGER(Kind=jpim),        INTENT(in)    :: nchannels
INTEGER(Kind=jpim),        INTENT(in)    :: nprofiles
INTEGER(Kind=jpim),        INTENT(in)    :: channels(nfrequencies)
INTEGER(Kind=jpim),        INTENT(in)    :: polarisations(nchannels,3)      ! Channel indices
INTEGER(Kind=jpim),        INTENT(in)    :: lprofiles(nfrequencies)
TYPE(profile_Type),        INTENT(inout) :: profiles(nprofiles) ! Profiles on RTTOV levels
TYPE(profile_cloud_Type),  INTENT(in)    :: cld_profiles(nprofiles) ! Cloud profiles on NWP levels
TYPE(rttov_coef),          INTENT(in)    :: coef  ! Coefficients
LOGICAL,                   INTENT(in)    :: calcemis(nchannels)  ! switch for emmissivity calc.
REAL(Kind=jprb),           INTENT(inout) :: emissivity(nchannels) ! surface emmissivity
TYPE(radiance_cloud_Type), INTENT(inout) :: cld_radiance    ! radiances (mw/cm-1/ster/sq.m)
INTEGER(Kind=jpim),        INTENT(out)   :: errorstatus(nprofiles)  ! return flag
END SUBROUTINE rttov_cld

!!! #include "rttov_cld_k.interface"
SUBROUTINE Rttov_cld_k  ( &
     errorstatus,    & ! out
     nfrequencies,   & ! in
     nchannels,      & ! in
     nbtout,         & ! in
     nprofiles,      & ! in
     channels,       & ! in
     polarisations,  & ! in
     lprofiles,      & ! in
     profiles,       & ! in
     cld_profiles,   & ! in
     coef,           & ! in
     switchrad,      & ! in
     calcemis,       & ! in
     emissivity,     & ! inout
     profiles_k ,    & ! inout
     cld_profiles_k ,& ! inout
     emissivity_k ,  & ! inout
     cld_radiance)     ! inout
USE rttov_const, ONLY :   &
     errorstatus_success ,&
     errorstatus_fatal   ,&
     overlap_scheme
USE rttov_types, ONLY :    &
     rttov_coef           ,&
     geometry_Type        ,&
     profile_Type         ,&
     profile_cloud_Type   ,&
     radiance_cloud_Type
USE parkind1, ONLY : jpim     ,jprb
IMPLICIT NONE
INTEGER(Kind=jpim),        INTENT(in)    :: nfrequencies
INTEGER(Kind=jpim),        INTENT(in)    :: nchannels
INTEGER(Kind=jpim),        INTENT(in)    :: nbtout
INTEGER(Kind=jpim),        INTENT(in)    :: nprofiles
INTEGER(Kind=jpim),        INTENT(in)    :: channels(nfrequencies)
INTEGER(Kind=jpim),        INTENT(in)    :: polarisations(nchannels,3)
INTEGER(Kind=jpim),        INTENT(in)    :: lprofiles(nfrequencies)
LOGICAL,                   INTENT(in)    :: switchrad  ! true if input is BT
TYPE(profile_Type),        INTENT(inout) :: profiles(nprofiles)
TYPE(profile_cloud_Type),  INTENT(in)    :: cld_profiles(nprofiles)
TYPE(rttov_coef),          INTENT(in)    :: coef
LOGICAL,                   INTENT(in)    :: calcemis(nchannels)
REAL(Kind=jprb),           INTENT(inout) :: emissivity(nchannels)
TYPE(radiance_cloud_type), INTENT(inout) :: cld_radiance! in because of meme allocation
TYPE(profile_Type),        INTENT(inout) :: profiles_k(nchannels)
TYPE(profile_cloud_Type),  INTENT(inout) :: cld_profiles_k(nchannels)
REAL(Kind=jprb),           INTENT(inout) :: emissivity_k(nchannels)
INTEGER(Kind=jpim),        INTENT(out)   :: errorstatus(nprofiles)
END SUBROUTINE Rttov_cld_k


!!! #include "rttov_scatt.interface"
SUBROUTINE rttov_scatt(&
     & errorstatus,&
     & nwp_levels,&
     & nrt_levels,&
     & nfrequencies,&
     & nchannels,&
     & nbtout,&
     & nprofiles,&
     & polarisations,&
     & channels,&
     & frequencies,&
     & lprofiles,&
     & lsprofiles,&
     & profiles,&
     & cld_profiles,&
     & coef_rttov,&
     & coef_scatt,&
     & calcemiss,&
     & emissivity_in,&
     & cld_radiance ) 
USE rttov_types, ONLY :&
     & rttov_coef ,&
     & rttov_scatt_coef ,&
     & geometry_Type ,&
     & profile_Type ,&
     & profile_cloud_Type ,&
     & profile_scatt_aux ,&
     & transmission_Type ,&
     & radiance_Type ,&
     & radiance_cloud_Type 
USE parkind1, ONLY : jpim ,jprb
INTEGER (Kind=jpim), INTENT (in) :: nwp_levels
INTEGER (Kind=jpim), INTENT (in) :: nrt_levels
INTEGER (Kind=jpim), INTENT (in) :: nprofiles
INTEGER (Kind=jpim), INTENT (in) :: nfrequencies
INTEGER (Kind=jpim), INTENT (in) :: nchannels
INTEGER (Kind=jpim), INTENT (in) :: nbtout
INTEGER (Kind=jpim), INTENT (in) :: channels (nfrequencies)
INTEGER (Kind=jpim), INTENT (in) :: frequencies (nchannels)
INTEGER (Kind=jpim), INTENT (in) :: polarisations (nchannels,3)
INTEGER (Kind=jpim), INTENT (in) :: lprofiles (nfrequencies)
INTEGER (Kind=jpim), INTENT (in) :: lsprofiles (nchannels)
INTEGER (Kind=jpim), INTENT (out) :: errorstatus (nprofiles)
LOGICAL, INTENT (in) :: calcemiss (nchannels)
REAL (Kind=jprb), INTENT (in) :: emissivity_in (nchannels)
TYPE (profile_Type), INTENT (inout) :: profiles (nprofiles)
TYPE (rttov_coef), INTENT (in) :: coef_rttov
TYPE (rttov_scatt_coef), INTENT (in) :: coef_scatt
TYPE (profile_cloud_Type), INTENT (in) :: cld_profiles (nprofiles)
TYPE (radiance_cloud_Type), INTENT (inout) :: cld_radiance
END SUBROUTINE rttov_scatt

!!! #include "rttov_readcoeffs.interface"
SUBROUTINE rttov_readcoeffs  (&
     & errorstatus,  & ! out
     & coef,         & ! out
     & instrument,   & ! in Optional
     & kmyproc,      & ! in Optional
     & kioproc,      & ! in Optional
     & file_id,      & ! in Optional
     & channels      ) ! in Optional
USE rttov_const, ONLY :   &
     version             ,&
     release             ,&
     minor_version       ,&
     rttov_magic_string  ,&
     sensor_id_mw        ,&
     sensor_id_ir        ,&
     errorstatus_info    ,&
     errorstatus_success ,&
     errorstatus_fatal   ,&
     gas_id_mixed        ,&
     gas_id_watervapour  ,&
     gas_id_ozone        ,&
     gas_id_wvcont       ,&
     gas_id_co2          ,&
     gas_id_n2o          ,&
     gas_id_co           ,&
     gas_id_ch4          ,&
     gas_unit_specconc   ,&
     gas_unit_ppmv       ,&
     earthradius         ,&
     gas_name            ,&
     pressure_top
USE rttov_types, ONLY : &
     rttov_coef
USE parkind1, ONLY : jpim     ,jprb
IMPLICIT NONE
INTEGER(Kind=jpim), OPTIONAL, INTENT(in) :: kmyproc  ! logical processor id
INTEGER(Kind=jpim), OPTIONAL, INTENT(in) :: kioproc  ! processor dedicated for io
INTEGER(Kind=jpim), OPTIONAL, INTENT (in) :: instrument(3)  ! (platform, satellite identification, instrument) number
INTEGER(Kind=jpim), OPTIONAL, INTENT (in) :: file_id      ! file logical unit number
INTEGER(Kind=jpim), OPTIONAL, INTENT (in) :: channels(:)      ! list of channels to extract
INTEGER(Kind=jpim), INTENT (out) :: errorstatus       ! return code
TYPE( rttov_coef ), INTENT (out) :: coef   ! coefficients
END SUBROUTINE rttov_readcoeffs

!!! #include "rttov_initcoeffs.interface"
SUBROUTINE rttov_initcoeffs  (&
     & errorstatus,   &! out
     & coef,          &! out
     & knproc,        &! in Optional
     & kmyproc,       &! in Optional
     & kioproc        )! in Optional
USE rttov_const, ONLY :   &
     & sensor_id_mw        ,&
     & errorstatus_info    ,&
     & errorstatus_success ,&
     & errorstatus_fatal   ,&
     & gas_id_mixed        ,&
     & gas_id_watervapour  ,&
     & gas_id_ozone        ,&
     & gas_id_wvcont       ,&
     & gas_id_co2          ,&
     & gas_id_n2o          ,&
     & gas_id_co           ,&
     & gas_id_ch4          ,&
     & gas_unit_specconc   ,&
     & gas_unit_ppmv       ,&
     & earthradius         ,&
     & gas_name            ,&
     & pressure_top
! Imported Type Definitions:
USE rttov_types, ONLY : &
     & rttov_coef
USE parkind1, ONLY : jpim     ,jprb
IMPLICIT NONE
INTEGER(Kind=jpim), OPTIONAL, INTENT(in) :: knproc   ! number of procs
INTEGER(Kind=jpim), OPTIONAL, INTENT(in) :: kmyproc  ! logical processor id
INTEGER(Kind=jpim), OPTIONAL, INTENT(in) :: kioproc  ! procs dedicated for io
! scalar arguments with intent(out):
INTEGER(Kind=jpim), INTENT (out) :: errorstatus       ! return code
TYPE( rttov_coef ), INTENT (out) :: coef   ! coefficients
END SUBROUTINE rttov_initcoeffs

!!! #include "rttov_readscattcoeffs.interface"
SUBROUTINE rttov_readscattcoeffs  (&
     & errorstatus,   &! out
     & coef_rttov,    &! in
     & coef_scatt,    &! out
     & file_id       ) ! in Optional
! Imported Type Definitions:
USE rttov_types, ONLY : &
     & rttov_coef, &
     & rttov_scatt_coef
USE rttov_const, ONLY :   &
     & inst_name           ,&
     & platform_name       ,&
     & errorstatus_info    ,&
     & errorstatus_success ,&
     & errorstatus_fatal
USE parkind1, ONLY : jpim     ,jprb
IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
INTEGER(Kind=jpim), INTENT (out) :: errorstatus             ! return code
! scalar arguments with optional intent(in):
INTEGER(Kind=jpim), OPTIONAL, INTENT (in)  :: file_id       ! file logical unit number
! array arguments with intent(in):
TYPE( rttov_coef ), INTENT (in) :: coef_rttov               ! clear-sky coefficients
! array arguments with intent(out):
TYPE( rttov_scatt_coef ), INTENT (out) :: coef_scatt        ! coefficients
END SUBROUTINE rttov_readscattcoeffs

END INTERFACE
!!! #include "rttov_opencoeff.interface"
!!! #include "rttov_errorhandling.interface"
!!! #include "rttov_dealloc_coef.interface"
!!! #include "rttov_errorreport.interface"
#endif
!!!
!!!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!!!
INTEGER, INTENT(IN)   :: KDLON   !number of columns where the
! radiation calculations are performed
INTEGER, INTENT(IN)   :: KFLEV   !number of vertical levels where the
! radiation calculations are performed
INTEGER, INTENT(IN)   :: KSTATM  !index of the standard atmosphere level
                                !just above the model top
!!!
REAL, DIMENSION(:,:),     INTENT(IN) :: PEMIS  !Surface IR EMISsivity
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD !RADiative Surface Temperature
REAL, DIMENSION(:,:),     INTENT(IN) :: PSTATM !selected standard atmosphere
                                !
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT   !THeta at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT    !moist variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST !pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ    !Model level heights
!!!
!!!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! cloud fraction
REAL, DIMENSION(:,:),     INTENT(IN) :: PULVLKB ! U-wind at KB level
REAL, DIMENSION(:,:),     INTENT(IN) :: PVLVLKB ! V-wind at KB level
!!!
LOGICAL, INTENT(IN)                  :: OUSERI ! logical switch to compute both
! liquid and solid condensate (OUSERI=.TRUE.)
! or only liquid condensate (OUSERI=.FALSE.)
!!!
INTEGER, DIMENSION(:,:), INTENT(IN) :: KRTTOVINFO ! platform, satelit, sensor,
                                                  ! and selection calculations
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! File characteristics
!
#ifdef MNH_RTTOV_8
!!!
!!!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!!!
!!!
INTEGER, PARAMETER :: JPNSAT=3       ! No. of Satellite required
                                !
INTEGER :: JI,JJ,JK,JK1,JK2,JKRAD,JKF,JSAT,JC ! loop indexes
                                !
INTEGER :: IJSAT        ! number of columns/=NUNDEF which 
                        ! have to be treated in the table KRTTOVINFO(:,:)
INTEGER :: IIB,IIE        ! I index value of the first/last inner mass point
INTEGER :: IJB,IJE        ! J index value of the first/last inner mass point
INTEGER :: IKB,IKE        ! K index value of the first/last inner mass point
INTEGER :: IIU          ! array size for the first  index
INTEGER :: IJU          ! array size for the second index
INTEGER :: IKU          ! array size for the third  index
INTEGER :: IKR          ! real array size for the third  index
INTEGER (Kind=jpim) :: iwp_levels ! equal to IKR (call to rttov_scatt)
INTEGER :: IIJ          ! reformatted array index
INTEGER :: IKSTAE       ! level number of the STAndard atmosphere array
INTEGER :: IKUP         ! vertical level above which STAndard atmosphere data
INTEGER, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,3)) :: IKKOZ ! indice array used to
! vertically interpolate the ozone content on the model grid
                                !
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPAVE  ! mean-layer pressure
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTAVE  ! mean-layer temperature
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQVAVE ! mean-layer specific humidity
REAL, DIMENSION(:,:), ALLOCATABLE :: ZO3AVE ! mean-layer ozone content
REAL, DIMENSION(:),   ALLOCATABLE :: ZREMIS ! Reformatted PEMIS array
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEXNT ! Exner function
REAL, DIMENSION(SIZE(PSTATM,1)) :: ZSTAZZ,ZSTAOZ ! STAndard atmosphere height
!    and OZone content
REAL :: ZOZ ! variable used to interpolate the ozone profile

REAL, DIMENSION(:),   ALLOCATABLE   ::  ZULAT
REAL, DIMENSION(:),   ALLOCATABLE   ::  ZULON

REAL, DIMENSION(:,:,:), ALLOCATABLE   ::  ZTBTMP
REAL, DIMENSION(:,:), ALLOCATABLE :: ZANTMP, ZUTH
REAL :: ZZH, zdeg_to_rad, zrad_to_deg, zbeta, zalpha

! Other arrays for zenithal solar angle 
! REAL, DIMENSION(:,:),   ALLOCATABLE :: ZCOSZEN, ZSINZEN, ZAZIMSOL

! Other arrays for condensation
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZTEMP  ! Temperature
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZNCLD  ! grid scale cloud fraction
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRC    ! grid scale r_c (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRR    ! grid scale r_r (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRS    ! grid scale r_s (kg/kg)
! -----------------------------------------------------------------------------
INTEGER, PARAMETER :: JPLEV=43, JPNAV=3, JPNSAV=5, JPNSSV=6, JPNCVCLD=6  

REAL, DIMENSION(JPLEV) ::  ZPRES !Fixed level pressures used in RTTOV
REAL, DIMENSION(JPLEV) ::  ZPRES_INV
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZAV    !Profile array content
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZSAV   !Surface array content
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZSSV   !Surface Skin array content
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZAP    !Full-level Model Pressure (hPa)
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZAP_HL !Half-level Model Pressure (hPa)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCV    !Temperature and cloud variable
                                !on full-level model
REAL, DIMENSION(:), ALLOCATABLE :: ZANGL   !Satellite zenith angle (deg)
REAL, DIMENSION(:), ALLOCATABLE :: ZANGS   !Solar zenith angle (deg)
! -----------------------------------------------------------------------------
! Jacobian fields
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTEMPK, ZWVAPK
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTEMPKP, ZTEMPKPP, ZWVAPKP, ZWVAPKPP
! -----------------------------------------------------------------------------
! INDEXES AND TEMPORAL ARRAYS FOR VECTORIZATION
INTEGER :: JIS, IBEG, IEND, IDIM, ICPT
INTEGER, DIMENSION(:),   ALLOCATABLE :: IMSURFP
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZAVP, ZCVP
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZSAVP, ZSSVP, ZAPP, ZAP_HLP
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZZTMP, ZZTMPP
REAL, DIMENSION(:),   ALLOCATABLE :: ZANGLP, ZREMISP
LOGICAL, DIMENSION(:), ALLOCATABLE ::  GANGL
! -----------------------------------------------------------------------------
INTEGER :: INRAD = 2 ! INRAD=1 RADIANCE; INRAD=2 BRIGHTNESS TEMPERATURE
! -----------------------------------------------------------------------------
! Realistic maximum values for hydrometeor content in kg/kg
REAL :: ZRCMAX = 5.0E-03, ZRRMAX = 5.0E-03, ZRIMAX = 2.0E-03, ZRSMAX = 5.0E-03
! -----------------------------------------------------------------------------
INTEGER, DIMENSION(:), ALLOCATABLE :: IMSURF   !Surface type index
                                
INTEGER :: IKFBOT, IKFTOP, INDEX, ISUM, JLEV, JCH, IWATER, ICAN
REAL, DIMENSION(:), ALLOCATABLE :: ZTEXTR, ZQVEXTR !Array used in interpolation
REAL, DIMENSION(:), ALLOCATABLE :: ZQVSAT, ZVINT !Array used in interpolation
REAL, DIMENSION(:), ALLOCATABLE :: ZPSUM, ZTSUM, ZQVSUM, ZO3SUM !Array used in interpolation
REAL :: zconst, ZPS, ZTGRAD, ZQGRAD, ZOGRAD !variables used in interpolation
REAL, DIMENSION(:), ALLOCATABLE :: ZPIN, ZFIN, ZOUT
!  at the open of the file LFI routines 
CHARACTER(LEN=8)  :: YINST  
CHARACTER(LEN=4)  :: YBEG, YEND
CHARACTER(LEN=2)  :: YCHAN, YTWO   
CHARACTER(LEN=1)  :: YONE   
                               
INTEGER, PARAMETER :: JPPLAT=16

CHARACTER(LEN=3), DIMENSION(JPPLAT) :: YPLAT= (/ &
     'N  ','D  ','MET','GO ','GMS','FY2','TRM','ERS', &
     'EOS','MTP','ENV','MSG','FY1','ADS','MTS','CRL' /)
CHARACTER(LEN=2), DIMENSION(2) :: YLBL_MVIRI = (/ 'WV', 'IR'/)
CHARACTER(LEN=3), DIMENSION(7) :: YLBL_SSMI = (/ &
     '19V','19H','22V','37V','37H','85V','85H'/)
CHARACTER(LEN=3), DIMENSION(9) :: YLBL_TMI = (/ &
     '10V','10H','19V','19H','22V','37V','37H','85V','85H'/)
CHARACTER(LEN=3), DIMENSION(8) :: YLBL_SEVIRI = (/ &
     '039', '062','073','087','097','108','120','134'/)
CHARACTER(LEN=3), DIMENSION(4) :: YLBL_GOESI = (/ &
     '039', '067','107','120'/)

! -----------------------------------------------------------------------------
!*JPC*VECTORIZATION
!! One profile per run
!! INTEGER (Kind=jpim) :: nprofiles = 1
INTEGER (Kind=jpim) :: nprofiles, ntruepro
!*JPC*VECTORIZATION

! RTTOV_readcoeffs interface
! ====================
INTEGER(Kind=jpim) :: errorstatus
INTEGER(Kind=jpim) :: instrument(3)
TYPE( rttov_coef ) :: coef         ! coefficients
TYPE( rttov_scatt_coef ) :: coef_scatt

! RTTOV interface
! ====================
INTEGER(Kind=jpim), ALLOCATABLE :: rttov_errorstatus(:)  ! rttov error return code
INTEGER(Kind=jpim) :: nfrequencies
INTEGER(Kind=jpim) :: nchannels
INTEGER(Kind=jpim) :: nbtout
INTEGER(Kind=jpim), ALLOCATABLE :: channels   (:), n_chan(:)
INTEGER(Kind=jpim), ALLOCATABLE :: polarisations   (:,:)
INTEGER(Kind=jpim), ALLOCATABLE :: frequencies   (:)
INTEGER(Kind=jpim), ALLOCATABLE :: lprofiles  (:),lsprofiles(:),lsprofiles2(:)
TYPE(profile_Type), ALLOCATABLE  :: profiles(:)
TYPE(profile_cloud_type), ALLOCATABLE :: cld_profiles(:)
TYPE(transmission_type)               :: transmission
LOGICAL              :: addcloud = .FALSE.
LOGICAL, ALLOCATABLE         :: calcemis(:)
REAL(Kind=jprb), ALLOCATABLE :: emissivity (:)
TYPE(radiance_cloud_type)             :: radiance

REAL(Kind=jprb),    ALLOCATABLE :: input_emissivity (:)
CHARACTER (len=6)  :: NameOfRoutine = 'tstrad'
! RTTOV K/AD interface
! ====================
LOGICAL    :: switchrad  ! true if input is BT
TYPE(profile_Type), ALLOCATABLE :: profiles_k(:)
TYPE(profile_cloud_Type), ALLOCATABLE :: cld_profiles_k(:)
REAL(Kind=jprb),    ALLOCATABLE       :: emissivity_k (:)

! variables for input
! ====================
! Parameter for WV conversion used in all tstrad suite
REAL(Kind=jprb), PARAMETER :: q_mixratio_to_ppmv  = 1.60771704e+6_JPRB
REAL(Kind=jprb), PARAMETER :: o3_mixratio_to_ppmv = 6.03504e+5_JPRB
INTEGER(Kind=jpim) :: alloc_status(40)

TYPE(TFIELDDATA) :: TZFIELD

! - End of header --------------------------------------------------------
!!!----------------------------------------------------------------------------
!!!
!!!*       1.    INITIALIZATION OF CONSTANTS FOR TRANSFERT CODE
!!!              ----------------------------------------------
!!!

! JPC from refprof.dat
ZPRES=(/   0.100,    0.290,    0.690,    1.420,    2.611,    4.407, &
     6.950,   10.370,   14.810,   20.400,   27.260,   35.510, &
     45.290,   56.730,   69.970,   85.180,  102.050,  122.040, &
     143.840,  167.950,  194.360,  222.940,  253.710,  286.600, &
     321.500,  358.280,  396.810,  436.950,  478.540,  521.460, &
     565.540,  610.600,  656.430,  702.730,  749.120,  795.090, &
     839.950,  882.800,  922.460,  957.440,  985.880, 1005.430, &
     1013.250 /)

DO JK=1,JPLEV
  JKRAD=JPLEV-JK+1
  ZPRES_INV(JK)=ZPRES(JKRAD)*100. ! Conversion from hPa to Pa
END DO

errorstatus     = 0
alloc_status(:) = 0

PRINT *,'NB OF SAT SIZE(KRTTOVINFO,1)=',SIZE(KRTTOVINFO,1)
PRINT *,'NB OF SAT SIZE(KRTTOVINFO,2)=',SIZE(KRTTOVINFO,2)
DO JSAT=1,SIZE(KRTTOVINFO,2)
  IF (KRTTOVINFO(1,JSAT) /= NUNDEF) THEN
    IJSAT = JSAT
  END IF
END DO

JSAT=1
instrument(1)=KRTTOVINFO(1,JSAT)
instrument(2)=KRTTOVINFO(2,JSAT)
instrument(3)=KRTTOVINFO(3,JSAT)
PRINT *,'range(KRTTOVINFO(3,JSAT)) ',range(KRTTOVINFO(3,JSAT))
PRINT *,'range(instrument(3)) ',range(instrument(3))
CALL rttov_readcoeffs (errorstatus, coef, instrument)
CALL rttov_initcoeffs (errorstatus, coef)

switchrad = INRAD == 2
PRINT *,' RADIANCE OR TB CALCULATION: INRAD=',INRAD,' switchrad=',switchrad

!!!----------------------------------------------------------------------------
!!!
!!!*       2.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES
!!!              ----------------------------------------------
                               
IIU = SIZE(PTHT,1)
IJU = SIZE(PTHT,2)
IKU = SIZE(PTHT,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
IKR = IKE - IKB +1
                              
IKSTAE = SIZE(PSTATM,1)
IKUP   = IKE-JPVEXT+1

!*JPC*VECTORIZATION
! Determine the number of profiles per RTTOV run
nprofiles = JPPF
!*JPC*VECTORIZATION

  
!!!----------------------------------------------------------------------------
!!!
!!!*       3.    INITIALIZES THE MEAN-LAYER VARIABLES
!!!              ------------------------------------
                               
ALLOCATE(ZEXNT(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
ZEXNT(:,:,:)= ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
                              
! Pressure
ALLOCATE(ZPAVE(KDLON,KFLEV))
DO JK=IKB,IKE
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZPAVE(IIJ,JKRAD)  = PPABST(JI,JJ,JK)*0.01  !Pressure in hPa
    END DO
  END DO
END DO
                               
! Temperature
ALLOCATE(ZTEMP(IIU,IJU,IKU))
ZTEMP=PTHT*ZEXNT
ALLOCATE(ZTAVE(KDLON,KFLEV))
DO JK=IKB,IKE
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZTAVE(IIJ,JKRAD)  = ZTEMP(JI,JJ,JK)
    END DO
  END DO
END DO

! Water vapor
ALLOCATE(ZQVAVE(KDLON,KFLEV))
ZQVAVE(:,:) = 0.0
IF( SIZE(PRT(:,:,:,:),4) >= 1 ) THEN
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
        ZQVAVE(IIJ,JKRAD) = PRT(JI,JJ,JK,1)
      END DO
    END DO
  END DO
END IF

! Ozone
ALLOCATE(ZO3AVE(KDLON,KFLEV))
    
ZSTAOZ(:) = PSTATM(:,6)/PSTATM(:,4)
ZSTAZZ(:) = 1000.0*PSTATM(:,1)

DO JJ = IJB,IJE
  DO JK2 = IKB,IKE
    JKRAD = JK2-JPVEXT
    IKKOZ(:,JK2) = IKB-1
    DO JK1 = 1,IKSTAE
      DO JI = IIB,IIE
        IKKOZ(JI,JK2)=IKKOZ(JI,JK2) + NINT(0.5 + SIGN(0.5,    &
             -ZSTAZZ(JK1)+0.5*(PZZ(JI,JJ,JK2)+PZZ(JI,JJ,JK2+1)) ))
      END DO
    END DO
    DO JI = IIB,IIE
      ZOZ=(0.5*(PZZ(JI,JJ,JK2)+PZZ(JI,JJ,JK2+1))- ZSTAZZ(IKKOZ(JI,JK2))) &
           /( ZSTAZZ(IKKOZ(JI,JK2)+1)           - ZSTAZZ(IKKOZ(JI,JK2)))
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZO3AVE(IIJ,JKRAD) =( (1.- ZOZ) * ZSTAOZ(IKKOZ(JI,JK2))    &
           +  ZOZ     * ZSTAOZ(IKKOZ(JI,JK2)+1))
    END DO
  END DO
END DO
                                
!  Standard atmosphere extension                               
DO JK=IKUP,KFLEV
  JK1 = (KSTATM-1)+(JK-IKUP)
  JK2 = JK1+1
  ZPAVE(:,JK)  = 0.5*( PSTATM(JK1,2)+PSTATM(JK2,2) )
  ZTAVE(:,JK)  = 0.5*( PSTATM(JK1,3)+PSTATM(JK2,3) )
  ZQVAVE(:,JK) = 0.5*( PSTATM(JK1,5)/PSTATM(JK1,4)+PSTATM(JK2,5)/PSTATM(JK2,4))
  JK1 = (KSTATM)+(JK-IKUP)
  ZO3AVE(:,JK) = ZSTAOZ(JK1)
END DO
!!!
!!!----------------------------------------------------------------------------
!!!
!!!*       4.    INTERPOLATES THE ATMOSPHERIC VARIABLES ONTO THE RTTOV GRID
!              ----------------------------------------------------------
!!!WITH INVERSION OF VERTICAL LEVELS!
                               
ALLOCATE(ZAV(JPLEV,JPNAV,KDLON))

ALLOCATE(ZTEXTR(JPLEV))
ALLOCATE(ZQVEXTR(JPLEV))
ALLOCATE(ZVINT(JPLEV))
ISUM=JPLEV+KFLEV
ALLOCATE(ZPSUM(ISUM))
ALLOCATE(ZTSUM(ISUM))
ALLOCATE(ZQVSUM(ISUM))
ALLOCATE(ZO3SUM(ISUM))
ALLOCATE(ZQVSAT(ISUM))
ZPSUM(:)=0.
ZTSUM(:)=0.
ZQVSUM(:)=0.
ZO3SUM(:)=0.
ZQVSAT(:)=0.
zconst= 287./1005.
IWATER = coef % fmv_gas_pos( gas_id_watervapour )
DO JI=IIB,IIE
  DO JJ=IJB,IJE
    IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
    ZPS=XP00*0.01 * & !Surface Pressure in hPa
         (0.5*(ZEXNT(JI,JJ,IKB)+ZEXNT(JI,JJ,IKB-1)))**(XCPD/XRD)
    DO JK=1,KFLEV
      JKRAD = KFLEV-JK+1 !INVERSION OF VERTICAL LEVELS!
      ZPSUM(JKRAD)=ZPAVE(IIJ,JK)
      ZTSUM(JKRAD)=ZTAVE(IIJ,JK)
      ZQVSUM(JKRAD)=ZQVAVE(IIJ,JK)
      ZO3SUM(JKRAD)=ZO3AVE(IIJ,JK)
    END DO
    ZTSUM(KFLEV+1)=ZTAVE(IIJ,1)
    ZQVSUM(KFLEV+1)=ZQVAVE(IIJ,1)
    IKFBOT=0
    DO JKF=1,JPLEV
      IF (ZPRES(JKF) > ZPS) THEN
        IKFBOT=JKF
        EXIT
      END IF
    END DO
    INDEX = KFLEV
    IF (IKFBOT /= 0) THEN
!!!-----Extrapolates temperature below surface pressure-------------------
      INDEX=JPLEV-IKFBOT+1
      INDEX=INDEX+KFLEV+1
      ZTSUM((KFLEV+2):INDEX) = PTSRAD(JI,JJ)
      ZPSUM((KFLEV+1))=ZPS
      ZPSUM((KFLEV+2):INDEX)=ZPRES(IKFBOT:JPLEV)
      ZO3SUM((KFLEV+2):INDEX)=ZO3AVE(IIJ,1)
    END IF
!!!-----Extrapolates profile above highest declared level-----------------
!!!----- => linear extrapolation -----------------------------------------
    IKFTOP = 1
    DO JLEV=1,INDEX
      IF(ZPRES(JLEV) >= ZPAVE(IIJ,KFLEV) ) EXIT
      IKFTOP = IKFTOP + 1
    END DO
    IF (IKFTOP /= 1) THEN
      ZTGRAD = (ZTSUM(1) - ZTSUM(2)) / (ZPSUM(1)-ZPSUM(2))
      ZQGRAD = (ZQVSUM(1) - ZQVSUM(2)) / (ZPSUM(1)-ZPSUM(2))
      ZOGRAD = (ZO3SUM(1) - ZO3SUM(2)) / (ZPSUM(1)-ZPSUM(2))
      DO JLEV=INDEX, 1, -1
        ZTSUM(JLEV+IKFTOP-1) = ZTSUM(JLEV)
        ZQVSUM(JLEV+IKFTOP-1) = ZQVSUM(JLEV)
        ZO3SUM(JLEV+IKFTOP-1) = ZO3SUM(JLEV)
        ZPSUM(JLEV+IKFTOP-1) = ZPSUM(JLEV)
      END DO
      INDEX = INDEX + IKFTOP-1
      DO JLEV=1,IKFTOP-1
        ZPSUM(JLEV) = ZPRES(JLEV)
        ZTSUM(JLEV) = ZTSUM(IKFTOP)  &
             + ZTGRAD * (ZPSUM(JLEV) - ZPSUM(IKFTOP))
        ZQVSUM(JLEV) = ZQVSUM(IKFTOP)  &
             + ZQGRAD * (ZPSUM(JLEV) - ZPSUM(IKFTOP))
        ZO3SUM(JLEV) = ZO3SUM(IKFTOP)  &
             + ZOGRAD * (ZPSUM(JLEV) - ZPSUM(IKFTOP))
      END DO
    ENDIF
!!!-----Interpolates to given pressure grid-------------------------------
    ALLOCATE(ZPIN(INDEX))
    ALLOCATE(ZFIN(INDEX))
    ALLOCATE(ZOUT(JPLEV))
    DO JLEV=1,INDEX
      JKRAD=INDEX-JLEV+1
      ZPIN(JKRAD) = ZPSUM(JLEV)*100.
      ZFIN(JKRAD) = ZTSUM(JLEV)
    END DO
    CALL PINTER(ZFIN, ZPIN, ZFIN, ZFIN, ZOUT, ZPRES_INV, &
         1, 1, INDEX, 1, JPLEV, 'LOG', 'RHU.')
    DO JLEV=1,JPLEV
      JKRAD=JPLEV-JLEV+1
      ZVINT(JKRAD) = ZOUT(JLEV)
    END DO
    ZAV(:,1,IIJ)= ZVINT(:) ! temperature K
    DO JLEV=1,INDEX
      JKRAD=INDEX-JLEV+1
      ZFIN(JKRAD) = ZQVSUM(JLEV)
    END DO
    CALL PINTER(ZFIN, ZPIN, ZFIN, ZFIN, ZOUT, ZPRES_INV, &
         1, 1, INDEX, 1, JPLEV, 'LOG', 'RHU.')
    DO JLEV=1,JPLEV
      JKRAD=JPLEV-JLEV+1
      ZVINT(JKRAD) = ZOUT(JLEV)
    END DO
    ZAV(:,2,IIJ)= ZVINT(:)*q_mixratio_to_ppmv ! water vapor mr ppmv
    DO JLEV=1,INDEX
      JKRAD=INDEX-JLEV+1
      ZFIN(JKRAD) = ZO3SUM(JLEV)
    END DO
    CALL PINTER(ZFIN, ZPIN, ZFIN, ZFIN, ZOUT, ZPRES_INV, &
         1, 1, INDEX, 1, JPLEV, 'LOG', 'RHU.')
    DO JLEV=1,JPLEV
      JKRAD=JPLEV-JLEV+1
      ZVINT(JKRAD) = ZOUT(JLEV)
    END DO
    ZAV(:,3,IIJ)= ZVINT(:)*o3_mixratio_to_ppmv ! ozone mixing ratio  ppmv
     DO JLEV=1,JPLEV
       ZAV(JLEV,1,IIJ)= MAX(coef%lim_prfl_tmin(JLEV), &
            MIN(coef%lim_prfl_tmax(JLEV),ZAV(JLEV,1,IIJ)))
       ZAV(JLEV,2,IIJ)= MAX(coef%lim_prfl_gmin(JLEV,IWATER), &
            MIN(coef%lim_prfl_gmax(JLEV,IWATER),ZAV(JLEV,2,IIJ)))
     END DO
    DEALLOCATE(ZPIN,ZFIN,ZOUT)
  END DO
END DO
DEALLOCATE(ZVINT)
DEALLOCATE(ZPAVE,ZTAVE,ZQVAVE,ZO3AVE)
!
!--------------------------------------------------------------------------
!
!*       6.    CALLS THE RTTOV RADIATION CODE
!	       ------------------------------
!
!*       6.1   INITIALIZES 2D AND SURFACE FIELDS
!
!
ALLOCATE(ZANGS(KDLON))
ZANGS(:)=0. ! zenithal solar angle not used
!
ALLOCATE(IMSURF(KDLON))
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
    IMSURF(IIJ) = NINT(XSEA(JI,JJ)) ! Surface Mask 0=land, 1=sea, 2=sea-ice
  END DO
END DO
!
ALLOCATE(ZSAV(JPNSAV,KDLON)) ! Surface 2m array contents
! fields taken at first level rather than at 2m
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
    ZSAV(1,IIJ) = ZTEMP(JI,JJ,IKB) ! 2m temperature (K)
    ZSAV(2,IIJ) = PRT(JI,JJ,IKB,1)*q_mixratio_to_ppmv ! 2m water vapor (ppmv)
    ZSAV(3,IIJ) = XP00*0.01 * & !Surface Pressure in hPa
         (0.5*(ZEXNT(JI,JJ,IKB)+ZEXNT(JI,JJ,IKB-1)))**(XCPD/XRD)
    ZSAV(4,IIJ) = PULVLKB(JI,JJ) ! 2m wind speed u (m/s)
    ZSAV(5,IIJ) = PVLVLKB(JI,JJ) ! 2m wind speed v (m/s)
  END DO
END DO
!
ALLOCATE(ZSSV(JPNSSV,KDLON)) !Surface skin array contents
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
    ZSSV(1,IIJ) =  PTSRAD(JI,JJ)
    ZSSV(2,IIJ) = 2.3  ! FASTEM-2 land coef (Bare soil see Table 3 svr)
    ZSSV(3,IIJ) = 1.9  ! FASTEM-2 land coef
    ZSSV(4,IIJ) = 21.8 ! FASTEM-2 land coef
    ZSSV(5,IIJ) = 0.0  ! FASTEM-2 land coef
    ZSSV(6,IIJ) = 0.5  ! FASTEM-2 land coef
  END DO
END DO
!
!
ALLOCATE(ZAP(KDLON,IKR))
DO JK=IKB,IKE
  JKRAD = IKE-JK+1 !INVERSION OF VERTICAL LEVELS!
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZAP(IIJ,JKRAD)=PPABST(JI,JJ,JK)*0.01 !Pressure in hPa
    END DO
  END DO
END DO
!
!
ALLOCATE(ZAP_HL(KDLON,IKR+1))
DO JK=IKB,IKE+1
  JKRAD = IKE-JK+2 !INVERSION OF VERTICAL LEVELS!
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZAP_HL(IIJ,JKRAD)=XP00*0.01 * & !Pressure in hPa
           (0.5*(ZEXNT(JI,JJ,JK)+ZEXNT(JI,JJ,JK-1)))**(XCPD/XRD)
    END DO
  END DO
END DO
DEALLOCATE(ZEXNT)
!
!
ALLOCATE(ZNCLD(IIU,IJU,IKU))
ZNCLD=0.
ALLOCATE(ZRC(IIU,IJU,IKU))
ZRC=0.
ALLOCATE(ZRI(IIU,IJU,IKU))
ZRI=0.
ALLOCATE(ZRR(IIU,IJU,IKU))
ZRR=0.
ALLOCATE(ZRS(IIU,IJU,IKU))
ZRS=0.
IF( SIZE(PRT(:,:,:,:),4) >= 3 ) THEN
  ZRC=PRT(:,:,:,2)
  ZRR=PRT(:,:,:,3)
  IF( OUSERI ) THEN
! ice
    ZRI=PRT(:,:,:,4)
    ZRS=PRT(:,:,:,5)+PRT(:,:,:,6)
  END IF
  ZNCLD=PCLDFR
END IF

! temperature and cloud field on full-model levels
ALLOCATE(ZCV(KDLON,IKR,JPNCVCLD))
ZCV = 0.

DO JK=IKB,IKE
  JKRAD = IKE-JK+1 !INVERSION OF VERTICAL LEVELS!
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZCV(IIJ,JKRAD,1)=ZTEMP(JI,JJ,JK) !Temperature (K)
      ZCV(IIJ,JKRAD,2)=ZNCLD(JI,JJ,JK) !Cloud cover (fraction)
      ZCV(IIJ,JKRAD,3)=MIN(ZRCMAX,ZRC(JI,JJ,JK))   !Cloud liquid water (kg/kg)
      ZCV(IIJ,JKRAD,4)=MIN(ZRIMAX,ZRI(JI,JJ,JK))   !Cloud ice water (kg/kg)
! rttov_iniscatt modified
!      ZCV(IIJ,JKRAD,5)=ZRR(JI,JJ,JK)   !rain (kg/m2/s)
!      ZCV(IIJ,JKRAD,6)=ZRS(JI,JJ,JK)   !solid precipitation (kg/m2/s)
      ZCV(IIJ,JKRAD,5)=MIN(ZRRMAX,ZRR(JI,JJ,JK))  !rain (kg/kg)
      ZCV(IIJ,JKRAD,6)=MIN(ZRSMAX,ZRS(JI,JJ,JK))   !solid precipitation (kg/kg)
    END DO
  END DO
END DO
DEALLOCATE(ZTEMP,ZNCLD,ZRC,ZRI,ZRR,ZRS)
!
!
ALLOCATE(ZREMIS(KDLON))
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
    ZREMIS(IIJ)   = PEMIS(JI,JJ)
  END DO
END DO
!
ALLOCATE(ZULAT(KDLON))
ALLOCATE(ZULON(KDLON))
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
    ZULON(IIJ) = XLON(JI,JJ)
    ZULAT(IIJ) = XLAT(JI,JJ)
  END DO
END DO
!
!*       6.2   CALLS THE RTTOV ROUTINES
!
!
ALLOCATE( rttov_errorstatus(nprofiles))

! Profiles on RTTOV pressure levels
ALLOCATE( profiles(nprofiles))
DO JI = 1, nprofiles
! allocate model profiles atmospheric arrays with model levels dimension
  profiles(JI) % nlevels =  coef % nlevels
  ALLOCATE( profiles(JI) % p  ( coef % nlevels ) )
  ALLOCATE( profiles(JI) % t  ( coef % nlevels ) )
  ALLOCATE( profiles(JI) % q  ( coef % nlevels ) )
  ALLOCATE( profiles(JI) % o3 ( coef % nlevels ) )
  ALLOCATE( profiles(JI) % clw( coef % nlevels ) )
  profiles(JI) % p(:) = coef % ref_prfl_p(:)
END DO
! Cloud additional profiles
ALLOCATE( cld_profiles(nprofiles))
DO JI = 1, nprofiles
! allocate model profiles atmospheric arrays with model levels dimension
  cld_profiles(JI) % nlevels =  IKR
  ALLOCATE( cld_profiles(JI) % p  ( IKR ) )
  ALLOCATE( cld_profiles(JI) % ph ( IKR+1 ) )
  ALLOCATE( cld_profiles(JI) % t  ( IKR ) )
  ALLOCATE( cld_profiles(JI) % cc ( IKR ) )
  ALLOCATE( cld_profiles(JI) % clw( IKR ) )
  ALLOCATE( cld_profiles(JI) % ciw( IKR ) )
  ALLOCATE( cld_profiles(JI) % rain( IKR ) )
  ALLOCATE( cld_profiles(JI) % sp( IKR ) )
END DO

! -----------------------------------------------------------------------------
!              *** LOOP OVER SENSORS ***
! -----------------------------------------------------------------------------
DO JSAT=1,IJSAT ! loop over sensors
 
  instrument(1)=KRTTOVINFO(1,JSAT)
  instrument(2)=KRTTOVINFO(2,JSAT)
  instrument(3)=KRTTOVINFO(3,JSAT)
  PRINT *,' JSAT=',JSAT, instrument

! Read and initialise coefficients
! -----------------------------------------------------------------------------
  CALL rttov_readcoeffs (errorstatus, coef, instrument)
  IF(errorstatus /= 0) THEN
    WRITE(*,*) 'error rttov_readcoeffs :',errorstatus
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','CALL_RTTOV8','error rttov_readcoeffs')
  ENDIF
  CALL rttov_initcoeffs (errorstatus,coef)
  IF(errorstatus /= 0) THEN
    WRITE(*,*) 'error rttov_initcoeffs :',errorstatus
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','CALL_RTTOV8','error rttov_initcoeffs')
  ENDIF

  ! Read coef file for cloud/rain absorption/scattering
  IF( coef% id_sensor == sensor_id_mw) THEN
    CALL rttov_readscattcoeffs (errorstatus, coef, coef_scatt)
  ENDIF

  ALLOCATE(ZANGL(KDLON))
  ZANGL=XUNDEF
  IF (KRTTOVINFO(1,JSAT) == 1) THEN ! NOAA PLATFORM
    ZANGL=0.
  ELSEIF (KRTTOVINFO(1,JSAT) == 2) THEN ! DMSP PLATFORM
    ZANGL=53.1 ! see Saunders, 2002, RTTOV7 - science/validation rep, page 8
! METEOSAT PLATFORM
  ELSEIF (KRTTOVINFO(1,JSAT) == 3) THEN 
    CALL DETER_ANGLE(5, KDLON, ZULAT, ZULON, ZANGL)
! Conversion from cosecant to angle (deg)
    WHERE (ZANGL /= XUNDEF .AND. ZANGL /=0.) ZANGL=ACOS(1./ZANGL)*180./XPI
! MSG PLATFORM
  ELSEIF (KRTTOVINFO(1,JSAT) == 12) THEN
    CALL DETER_ANGLE(6, KDLON, ZULAT, ZULON, ZANGL)
! Conversion from cosecant to angle (deg)
    WHERE (ZANGL /= XUNDEF .AND. ZANGL /=0.) ZANGL=ACOS(1./ZANGL)*180./XPI
  ELSEIF (KRTTOVINFO(1,JSAT) == 4) THEN ! GOES-E PLATFORM
    CALL DETER_ANGLE(1, KDLON, ZULAT, ZULON, ZANGL)
! Conversion from cosecant to angle (deg)
    WHERE (ZANGL /= XUNDEF .AND. ZANGL /=0.) ZANGL=ACOS(1./ZANGL)*180./XPI
  ELSEIF (KRTTOVINFO(1,JSAT) == 7) THEN ! TRMM PLATFORM
    ZANGL=52.3 ! see Kummerow et al., J. Appl. Meteorol., Dec. 2000
  ENDIF
! Coefficients computed from transmittances for 6 viewing angles in the range 
! 0 to 63.6 deg (Saunders, 2002, RTTOV7 - science/validation rep., page 3)
  WHERE (ZANGL > 65.) ZANGL=65.

  ALLOCATE(n_chan(nprofiles))
  n_chan=coef%fmv_chn
  CALL rttov_setupchan(nprofiles,n_chan,coef,nfrequencies,nchannels,nbtout)

  ALLOCATE( channels   ( nfrequencies ) )
  ALLOCATE( lprofiles  ( nfrequencies ) )
  ALLOCATE( lsprofiles  ( nchannels ) )
  ALLOCATE( lsprofiles2  ( nbtout ) )
  ALLOCATE( emissivity ( nchannels ) )
  ALLOCATE( frequencies ( nchannels ) ) 
  ALLOCATE( polarisations ( nchannels ,3) )
  ALLOCATE( input_emissivity ( nchannels ) )
  ALLOCATE( calcemis  ( nchannels ) )

  ALLOCATE( transmission % tau_surf      ( nchannels ) )
  ALLOCATE( transmission % tau_layer     ( coef % nlevels, nchannels ) )
  ALLOCATE( transmission % od_singlelayer( coef % nlevels, nchannels ) )

  calcemis(1:nchannels)              = .TRUE.
  input_emissivity(1:nchannels)      = 0.5
  emissivity(1:nchannels)            = 0.

! allocate radiance results arrays with number of channels
  ALLOCATE( radiance % clear    ( nchannels ) )
  ALLOCATE( radiance % cloudy   ( nchannels ) )
  ALLOCATE( radiance % total    ( nchannels ) )
  ALLOCATE( radiance % bt       ( nchannels ) )
  ALLOCATE( radiance % bt_clear ( nchannels ) )
  ALLOCATE( radiance % upclear  ( nchannels ) )
  ALLOCATE( radiance % dnclear  ( nchannels ) )
  ALLOCATE( radiance % reflclear( nchannels ) )
  ALLOCATE( radiance % overcast ( IKR, nchannels ) )
  ALLOCATE( radiance % downcld  ( IKR, nchannels ) )
  ALLOCATE( radiance % cldemis  ( IKR, nchannels ) )
  ALLOCATE( radiance % wtoa     ( IKR, nchannels ) )
  ALLOCATE( radiance % wsurf    ( IKR, nchannels ) )
  ALLOCATE( radiance % cs_wtoa  ( nchannels ) )
  ALLOCATE( radiance % cs_wsurf ( nchannels ) )
  ALLOCATE( radiance % out  ( nbtout ) )
  ALLOCATE( radiance % out_clear( nbtout ) )
  ALLOCATE( radiance % total_out( nbtout ) )
  ALLOCATE( radiance % clear_out( nbtout ) )
  ALLOCATE( radiance % freq_used( nchannels) )

! Allocate new profiles for K code
  IF ( KRTTOVINFO(4,JSAT) == 1 .OR. KRTTOVINFO(4,JSAT) == 3) THEN
! Profiles on RTTOV pressure levels
    ALLOCATE( profiles_k(nchannels))
    DO JI = 1, nchannels
! allocate model profiles atmospheric arrays with model levels dimension
      profiles_k(JI) % nlevels =  coef % nlevels
      ALLOCATE( profiles_k(JI) % p  ( coef % nlevels ) )
      ALLOCATE( profiles_k(JI) % t  ( coef % nlevels ) )
      ALLOCATE( profiles_k(JI) % q  ( coef % nlevels ) )
      ALLOCATE( profiles_k(JI) % o3 ( coef % nlevels ) )
      ALLOCATE( profiles_k(JI) % clw( coef % nlevels ) )
      profiles_k(JI) % p(:) = coef % ref_prfl_p(:)
    END DO
! Cloud additional profiles
    ALLOCATE( cld_profiles_k(nchannels))
    DO JI = 1, nchannels
! allocate model profiles atmospheric arrays with model levels dimension
      cld_profiles_k(JI) % nlevels =  IKR
      ALLOCATE( cld_profiles_k(JI) % p  ( IKR ) )
      ALLOCATE( cld_profiles_k(JI) % ph ( IKR+1 ) )
      ALLOCATE( cld_profiles_k(JI) % t  ( IKR ) )
      ALLOCATE( cld_profiles_k(JI) % cc ( IKR ) )
      ALLOCATE( cld_profiles_k(JI) % clw( IKR ) )
      ALLOCATE( cld_profiles_k(JI) % ciw( IKR ) )
    END DO
    ALLOCATE( emissivity_k( nchannels ))
  END IF


! fixed values
  profiles(1:nprofiles) % ozone_data = .TRUE. 
  profiles(1:nprofiles) % co2_data   = .FALSE.
  profiles(1:nprofiles) % clw_data   = .FALSE.
  profiles(1:nprofiles) % s2m % o    = 0.
  profiles(1:nprofiles) % azangle    = 0. !!!!!! WARNING
  profiles(1:nprofiles) % ctp        = 500._JPRB  ! default value
  profiles(1:nprofiles) % cfraction  = 0._JPRB    ! default value
! See rttov_emiscld.F90
  cld_profiles(1:nprofiles) % kice   = 0          ! Hexagonal columns
!  cld_profiles(1:nprofiles) % kice   = 1          ! Aggregates
!  cld_profiles(1:nprofiles) % kradip = 0          ! Ou-Liou
!  cld_profiles(1:nprofiles) % kradip = 1          ! Wyser
!  cld_profiles(1:nprofiles) % kradip = 2          ! Boudala et al.
  cld_profiles(1:nprofiles) % kradip = 3          ! McFarquhar
    
  PRINT *,'cld_profiles % kice = ',cld_profiles(1) % kice
  PRINT *,'cld_profiles % kradip = ',cld_profiles(1) % kradip
   
  CALL rttov_setupindex (n_chan,nprofiles,nfrequencies,nchannels,nbtout,coef, &
       & input_emissivity,lprofiles,channels,polarisations,emissivity) 

!!! Set up remaining indices
  IF( coef% id_sensor == sensor_id_mw) &
       CALL rttov_scatt_setupindex (nprofiles,n_chan,coef,nchannels, &
       & lsprofiles, lsprofiles2, frequencies,nbtout)

!!! METEOSAT, GOES, OR MSG PLATFORM
  IF (KRTTOVINFO(1,JSAT) == 3 .OR. KRTTOVINFO(1,JSAT) == 4 &
       .OR. KRTTOVINFO(1,JSAT) == 12) &
       calcemis(1:nchannels) = .FALSE.


  ALLOCATE(GANGL(KDLON))
  GANGL(:) = .TRUE.
  WHERE( ZANGL(:) == XUNDEF) 
    GANGL(:) = .FALSE.
  END WHERE

  IDIM = COUNT( GANGL(:) )  ! number of columns with a defined sat angle

  ALLOCATE(ZANGLP(IDIM))
  ZANGLP  = PACK( ZANGL,MASK=GANGL )

  ALLOCATE(ZAVP(JPLEV,JPNAV,IDIM))
  DO JC=1,JPNAV
    DO JK=1,JPLEV
      ZAVP(JK,JC,:)  = PACK( ZAV(JK,JC,:),MASK=GANGL )
    END DO
  END DO

  ALLOCATE(ZSAVP(JPNSAV,IDIM)) 
  DO JK=1,JPNSAV
    ZSAVP(JK,:) = PACK( ZSAV(JK,:),MASK=GANGL )
  END DO

  ALLOCATE(IMSURFP(IDIM))
  IMSURFP  = PACK( IMSURF,MASK=GANGL )

  ALLOCATE(ZSSVP(JPNSSV,IDIM)) 
  DO JK=1,JPNSSV
    ZSSVP(JK,:)  = PACK( ZSSV(JK,:),MASK=GANGL )
  END DO

  ALLOCATE(ZCVP(IDIM,IKR,JPNCVCLD))
  DO JC=1,JPNCVCLD
    DO JK=1,IKR
      ZCVP(:,JK,JC)  = PACK( ZCV(:,JK,JC),MASK=GANGL )
    END DO
  END DO

  ALLOCATE(ZAPP(IDIM,IKR))
  DO JK=1,IKR
    ZAPP(:,JK)  = PACK( ZAP(:,JK),MASK=GANGL )
  END DO

  ALLOCATE(ZAP_HLP(IDIM,IKR+1))
  DO JK=1,IKR+1
    ZAP_HLP(:,JK)  = PACK( ZAP_HL(:,JK),MASK=GANGL )
  END DO

  ALLOCATE(ZREMISP(IDIM))
  ZREMISP  = PACK( ZREMIS,MASK=GANGL )

  ALLOCATE(ZZTMP(coef%fmv_chn,KDLON))
  ALLOCATE(ZZTMPP(coef%fmv_chn,IDIM))
  ZZTMP=XUNDEF
  ZZTMPP=XUNDEF

  IF ( KRTTOVINFO(4,JSAT) == 1 .OR. KRTTOVINFO(4,JSAT) == 3) THEN
    ALLOCATE(ZTEMPKP(coef%fmv_chn,KDLON,JPLEV))
    ALLOCATE(ZTEMPKPP(coef%fmv_chn,IDIM,JPLEV))
    ALLOCATE(ZWVAPKP(coef%fmv_chn,KDLON,JPLEV))
    ALLOCATE(ZWVAPKPP(coef%fmv_chn,IDIM,JPLEV))
    ZTEMPKP=XUNDEF
    ZTEMPKPP=XUNDEF
    ZWVAPKP=XUNDEF
    ZWVAPKPP=XUNDEF
  ENDIF
    
  DO JIS=1,IDIM,nprofiles
    IBEG = JIS
    IEND = MIN(JIS+nprofiles-1,IDIM)
    ntruepro=IEND-IBEG+1

    ICPT=IBEG
    DO JI=1,ntruepro
      profiles(JI) % t(:) = ZAVP(:,1,ICPT)
      profiles(JI) % q(:) = ZAVP(:,2,ICPT)
      profiles(JI) % o3(:) = ZAVP(:,3,ICPT)
! Surface
      profiles(JI) % s2m % p = ZSAVP(3,ICPT)
      profiles(JI) % s2m % q = ZSAVP(2,ICPT)
      profiles(JI) % s2m % t = ZSAVP(1,ICPT)
      profiles(JI) % s2m % u = ZSAVP(4,ICPT)
      profiles(JI) % s2m % v = ZSAVP(5,ICPT)
      profiles(JI) % skin % surftype = IMSURFP(ICPT)
      profiles(JI) % skin % t = ZSSVP(1,ICPT)
      profiles(JI) % skin % fastem(:) = &
! RTTOV 8.5 example
!        (/ 3.0_JPRB, 5.0_JPRB, 15.0_JPRB, 0.1_JPRB, 0.3_JPRB /)
! Bare soil see Table 3 svr rttov7)
         (/ 2.3_JPRB, 1.9_JPRB, 21.8_JPRB, 0.0_JPRB, 0.5_JPRB /)
! Angles
      profiles(JI) % zenangle   = ZANGLP(ICPT)
! Cloudy atmosphere on Meso-NH levels
      cld_profiles(JI) % p  (:) = ZAPP(ICPT,:)
      cld_profiles(JI) % ph (:) = ZAP_HLP(ICPT,:)
      cld_profiles(JI) % t  (:) = ZCVP(ICPT,:,1)
      cld_profiles(JI) % cc (:) = ZCVP(ICPT,:,2)
      cld_profiles(JI) % clw(:) = ZCVP(ICPT,:,3)
      cld_profiles(JI) % ciw(:) = ZCVP(ICPT,:,4)
      cld_profiles(JI) % rain(:) = ZCVP(ICPT,:,5)
      cld_profiles(JI) % sp(:) = ZCVP(ICPT,:,6)
      ICPT=ICPT+1
    END DO

    ICAN=0
    ICPT=IBEG
    DO JI=1,ntruepro
      DO JCH=1,coef%fmv_chn
        ICAN=ICAN+1
        IF (.NOT.calcemis(ICAN)) emissivity(ICAN) = ZREMISP(ICPT)
      END DO
      ICPT=ICPT+1
    END DO

    IF( coef% id_sensor /= sensor_id_mw) THEN
      CALL rttov_cld( &
               & rttov_errorstatus,  &! out
               & nfrequencies,   &! in
               & nchannels,      &! in
               & nbtout,         &! in
               & nprofiles,      &! in
               & channels,       &! in
               & polarisations,  &! in
               & lprofiles,      &! in
               & profiles,       &! inout  (to invalid clw absorption)
               & cld_profiles,   &! in
               & coef,           &! in
               & calcemis,       &! in
               & emissivity,     &! inout
               & radiance )       ! inout
    ELSE
      iwp_levels=IKR
      CALL rttov_scatt( &
               & rttov_errorstatus,  &! out
               & iwp_levels, & ! in
               & coef%nlevels, & ! in
               & nfrequencies,   &! in
               & nchannels,      &! in
               & nbtout,         &! in
               & nprofiles,      &! in
               & polarisations,  &! in
               & channels,           & ! in
               & frequencies,        & ! in
               & lprofiles,      &! in
               & lsprofiles,         & ! in
               & profiles,       &! inout  (to invalid clw absorption)
               & cld_profiles,   &! in
               & coef,           &! in
               & coef_scatt,     &! in
               & calcemis,       &! in
               & emissivity,     &! inout
               & radiance )       ! inout
    END IF

    IF (INRAD==1) THEN
! cloudy radiance for given cloud
      ICAN=0
      ICPT=IBEG
      DO JI=1,ntruepro
        DO JCH=1,coef%fmv_chn
          ICAN=ICAN+1
          ZZTMPP(JCH,ICPT) = radiance%total_out (ICAN) 
        END DO
        ICPT=ICPT+1
      END DO
    ELSE
! BT equivalent to total radiance
      ICAN=0
      ICPT=IBEG
      DO JI=1,ntruepro
        DO JCH=1,coef%fmv_chn
          ICAN=ICAN+1
          ZZTMPP(JCH,ICPT) = radiance%out (ICAN) 
        END DO
        ICPT=ICPT+1
      END DO
    ENDIF
!   PRINT *,'size',coef%fmv_chn,IDIM,KDLON,SIZE(ZZTMPP,1),SIZE(ZZTMPP,2)
!   PRINT *,'ZZTMP min/max ',MINVAL(ZZTMPP(:,:)),MAXVAL(ZZTMPP(:,:))


! Calling for K code
    IF ( KRTTOVINFO(4,JSAT) == 1 .OR. KRTTOVINFO(4,JSAT) == 3) THEN
!!!    IF (JIS==1) THEN
!!      IF( coef% id_sensor /= sensor_id_mw) THEN
      CALL rttov_cld_k  ( &
           & rttov_errorstatus,     &! out
           & nfrequencies,    &! in
           & nchannels,       &! in
           & nbtout,          &! in
           & nprofiles,       &! in
           & channels,        &! in
           & polarisations,   &! in
           & lprofiles,       &! in
           & profiles,        &! in
           & cld_profiles,    &! in
           & coef,            &! in
           & switchrad,       &! in
           & calcemis,        &! in
           & emissivity,      &! inout
           & profiles_k ,     &! inout
           & cld_profiles_k , &! inout
           & emissivity_k ,   &! inout
           & radiance)         ! inout
!!!      ENDIF

      ICAN=0
      ICPT=IBEG
      DO JI=1,ntruepro
        DO JCH=1,coef%fmv_chn
          ICAN=ICAN+1
          DO JK=1,JPLEV
            ZTEMPKPP(JCH,ICPT,JK) = profiles_k(ICAN) % t (JK) 
            ZWVAPKPP(JCH,ICPT,JK) = profiles_k(ICAN) % q (JK) 
          END DO
        END DO
        ICPT=ICPT+1
      END DO
!       DO JK=1,JPLEV
!         PRINT *,JK,' temp ',MINVAL(ZTEMPKPP(:,:,JK)),MAXVAL(ZTEMPKPP(:,:,JK))
!         PRINT *,JK,' vap ',MINVAL(ZWVAPKPP(:,:,JK)),MAXVAL(ZWVAPKPP(:,:,JK))
!       END DO
    END IF
  END DO
! Unpack the vector
  DO JCH=1,coef%fmv_chn
    ZZTMP(JCH,:)  = UNPACK( ZZTMPP(JCH,:), MASK=GANGL, FIELD=XUNDEF )
  END DO
  DEALLOCATE(ZZTMPP,ZANGLP)
  DEALLOCATE(ZAVP,ZSAVP,IMSURFP,ZSSVP,ZCVP,ZAPP,ZAP_HLP,ZREMISP)
! -----------------------------------------------------------------------------
! Generate angle and BT images
  ALLOCATE(ZANTMP(IIU,IJU))
  ZANTMP = XUNDEF
  ALLOCATE(ZTBTMP(IIU,IJU,coef%fmv_chn))
  ZTBTMP = XUNDEF
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
      ZANTMP(JI,JJ) = ZANGL(IIJ)
      ZTBTMP(JI,JJ,:) = ZZTMP(:,IIJ)
    END DO
  END DO
  DEALLOCATE(ZANGL,ZZTMP)
! -----------------------------------------------------------------------------
  IF ( KRTTOVINFO(4,JSAT) == 1 .OR. KRTTOVINFO(4,JSAT) == 3) THEN
    DO JCH=1,coef%fmv_chn
      DO JK=1,JPLEV
        ZTEMPKP(JCH,:,JK)=UNPACK(ZTEMPKPP(JCH,:,JK),MASK=GANGL,FIELD=XUNDEF )
        ZWVAPKP(JCH,:,JK)=UNPACK(ZWVAPKPP(JCH,:,JK),MASK=GANGL,FIELD=XUNDEF )
      END DO
    END DO
    DEALLOCATE(ZTEMPKPP,ZWVAPKPP)
  ENDIF
! -----------------------------------------------------------------------------
  IF (KRTTOVINFO(3,JSAT) == 20) THEN ! MVIRI
    YINST='MVIRI'
!    YINST=inst_name(KRTTOVINFO(3,JSAT))
!    DO JK1=1,LEN_TRIM(inst_name(KRTTOVINFO(3,JSAT)))
!      YINST(JK1:JK1)=CHAR(ICHAR(YINST(JK1:JK1))-32)
!    END DO
    TZFIELD%CMNHNAME   = TRIM(YINST)//'_ANGL'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'degree'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = TRIM(YINST)//' ANGLE'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    PRINT *,TZFIELD%CMNHNAME//TZFIELD%CCOMMENT
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZANTMP)
  END IF
  DEALLOCATE(ZANTMP)
! -----------------------------------------------------------------------------
  YBEG='    '
  IF (KRTTOVINFO(1,JSAT) <= 2 .OR. KRTTOVINFO(1,JSAT) == 4) THEN ! NOAA
    WRITE(YTWO,'(I2.2)') KRTTOVINFO(2,JSAT)
    YBEG=TRIM(YPLAT(KRTTOVINFO(1,JSAT)))//YTWO
  ELSEIF (KRTTOVINFO(1,JSAT) <= JPPLAT) THEN
    WRITE(YONE,'(I1.1)') KRTTOVINFO(2,JSAT)
    YBEG=TRIM(YPLAT(KRTTOVINFO(1,JSAT)))//YONE
  ELSE
    YBEG='XXXX'
  END IF
  WRITE(YTWO,'(I2.2)') KRTTOVINFO(3,JSAT)
!*JPC*VECTORIZATION
!  DO JCH=1,nbtout
  DO JCH=1,coef%fmv_chn
!*JPC*VECTORIZATION
    YEND='    '
    WRITE(YCHAN,'(I2.2)') JCH
    IF (KRTTOVINFO(3,JSAT) == 0) THEN ! HIRS
      YEND='H'//YCHAN
    ELSEIF (KRTTOVINFO(3,JSAT) == 3) THEN ! AMSU-A
      YEND='A'//YCHAN
    ELSEIF (KRTTOVINFO(3,JSAT) == 4) THEN ! AMSU-B
      YEND='B'//YCHAN
    ELSEIF (KRTTOVINFO(3,JSAT) == 6) THEN ! SSMI
      YEND=YLBL_SSMI(JCH)
    ELSEIF (KRTTOVINFO(3,JSAT) == 9) THEN ! TMI
      YEND=YLBL_TMI(JCH)
    ELSEIF (KRTTOVINFO(3,JSAT) == 20) THEN ! MVIRI
      YEND=YLBL_MVIRI(JCH)
    ELSEIF (KRTTOVINFO(3,JSAT) == 21) THEN ! SEVIRI
      YEND=YLBL_SEVIRI(JCH)
    ELSEIF (KRTTOVINFO(3,JSAT) == 22) THEN ! GOES-I
      YEND=YLBL_GOESI(JCH)
    ELSE
      YEND=YTWO//YCHAN
    END IF
    IF (INRAD==1) THEN
      TZFIELD%CMNHNAME   = TRIM(YBEG)//'_'//TRIM(YEND)//'rad'
      TZFIELD%CUNITS     = 'mw/cm-1/ster/sq.m'
      TZFIELD%CCOMMENT   = TRIM(YBEG)//'_'//TRIM(YEND)//' rad'
    ELSE
      TZFIELD%CMNHNAME   = TRIM(YBEG)//'_'//TRIM(YEND)//'BT'
      TZFIELD%CUNITS     = 'K'
      TZFIELD%CCOMMENT   = TRIM(YBEG)//'_'//TRIM(YEND)//' BT'
    ENDIF
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    PRINT *,TZFIELD%CMNHNAME//TZFIELD%CCOMMENT, &
         MINVAL(ZTBTMP(:,:,JCH),ZTBTMP(:,:,JCH)/=XUNDEF), &
         MAXVAL(ZTBTMP(:,:,JCH),ZTBTMP(:,:,JCH)/=XUNDEF)
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZTBTMP(:,:,JCH))
    IF (KRTTOVINFO(3,JSAT) == 4.AND. JCH==3 ) THEN ! AMSU-B
      TZFIELD%CMNHNAME   = TRIM(YBEG)//'_UTH'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'percent'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = TRIM(YBEG)//'_UTH'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 2
      TZFIELD%LTIMEDEP   = .TRUE.
! UTH computation from Buehler and John JGR 2005
      ZZH= 833000. ! (m) nominal altitude of the satellite
      zdeg_to_rad = XPI / 180.0
      zrad_to_deg = 180.0 / XPI
      zbeta = zdeg_to_rad*0.55 ! angle of incident radiation
! viewing angle alpha
      zalpha = zrad_to_deg*ASIN(XRADIUS/(XRADIUS+zzh)*SIN(zbeta))
      ALLOCATE(ZUTH(IIU,IJU))
      ZUTH = XUNDEF
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IF (ZTBTMP(JI,JJ,JCH)/=XUNDEF) THEN
            ZUTH(JI,JJ) = 100.*COS(zdeg_to_rad*zalpha)   &
                 *EXP(18.341-0.0764737*ZTBTMP(JI,JJ,JCH))
          END IF
        END DO
      END DO
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZUTH)
      DEALLOCATE(ZUTH)
    END IF
  END DO
! -----------------------------------------------------------------------------
! Jacobian fields
  IF ( KRTTOVINFO(4,JSAT) == 1 .OR. KRTTOVINFO(4,JSAT) == 3) THEN
    ALLOCATE(ZTEMPK(IIU,IJU,IKU))
    ALLOCATE(ZWVAPK(IIU,IJU,IKU))
    ALLOCATE(ZFIN(JPLEV))
    DO JCH=1,coef%fmv_chn
      YEND='    '
      WRITE(YCHAN,'(I2.2)') JCH
      IF (KRTTOVINFO(3,JSAT) == 0) THEN ! HIRS
        YEND='H'//YCHAN
      ELSEIF (KRTTOVINFO(3,JSAT) == 3) THEN ! AMSU-A
        YEND='A'//YCHAN
      ELSEIF (KRTTOVINFO(3,JSAT) == 4) THEN ! AMSU-B
        YEND='B'//YCHAN
      ELSEIF (KRTTOVINFO(3,JSAT) == 6) THEN ! SSMI
        YEND=YLBL_SSMI(JCH)
      ELSEIF (KRTTOVINFO(3,JSAT) == 9) THEN ! TMI
        YEND=YLBL_TMI(JCH)
      ELSEIF (KRTTOVINFO(3,JSAT) == 20) THEN ! MVIRI
        YEND=YLBL_MVIRI(JCH)
      ELSEIF (KRTTOVINFO(3,JSAT) == 21) THEN ! SEVIRI
        YEND=YLBL_SEVIRI(JCH)
      ELSEIF (KRTTOVINFO(3,JSAT) == 22) THEN ! GOES-I
        YEND=YLBL_GOESI(JCH)
      ELSE
        YEND=YTWO//YCHAN
      END IF
      ZTEMPK = XUNDEF
      ZWVAPK = XUNDEF
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = (JI-JPHEXT) + (IIE-IIB+1)*(JJ-IJB)
          DO JK=1,JPLEV
            JKRAD=JPLEV-JK+1
            ZFIN(JK)=ZTEMPKP(JCH,IIJ,JKRAD)
          END DO
          CALL PINTER(ZFIN, ZPRES_INV, ZFIN, ZFIN, &
               ZTEMPK(JI,JJ,IKB:IKE), PPABST(JI,JJ,IKB:IKE), &
               1, 1, JPLEV, 1, IKR, 'LOG', 'RHU.')
          DO JK=1,JPLEV
            JKRAD=JPLEV-JK+1
            ZFIN(JK)=ZWVAPKP(JCH,IIJ,JKRAD)
          END DO
          CALL PINTER(ZFIN, ZPRES_INV, ZFIN, ZFIN, &
               ZWVAPK(JI,JJ,IKB:IKE), PPABST(JI,JJ,IKB:IKE), &
               1, 1, JPLEV, 1, IKR, 'LOG', 'RHU.')
        END DO
      END DO
      !
      TZFIELD%CMNHNAME   = TRIM(YBEG)//'_'//TRIM(YEND)//'JAT'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'K K-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = TRIM(YBEG)//'_'//TRIM(YEND)//' JATEMP'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 2
      TZFIELD%LTIMEDEP   = .TRUE.
      PRINT *,TZFIELD%CMNHNAME//TZFIELD%CCOMMENT, &
           MINVAL(ZTEMPK(:,:,:),ZTEMPK(:,:,:)/=XUNDEF), &
           MAXVAL(ZTEMPK(:,:,:),ZTEMPK(:,:,:)/=XUNDEF)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZTEMPK(:,:,:))
      !
      TZFIELD%CMNHNAME   = TRIM(YBEG)//'_'//TRIM(YEND)//'JAV'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'K'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = TRIM(YBEG)//'_'//TRIM(YEND)//' JAWVAP'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 2
      TZFIELD%LTIMEDEP   = .TRUE.
      WHERE (ZWVAPK(:,:,:) /= XUNDEF) &
           ZWVAPK(:,:,:)=ZWVAPK(:,:,:)*(-0.1*PRT(:,:,:,1))
      PRINT *,TZFIELD%CMNHNAME//TZFIELD%CCOMMENT, &
           MINVAL(ZWVAPK(:,:,:),ZWVAPK(:,:,:)/=XUNDEF), &
           MAXVAL(ZWVAPK(:,:,:),ZWVAPK(:,:,:)/=XUNDEF)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWVAPK(:,:,:))
    END DO
    DEALLOCATE(ZTEMPKP,ZWVAPKP,ZFIN)
  ENDIF
! -----------------------------------------------------------------------------
  DEALLOCATE(GANGL,ZTBTMP)
  DEALLOCATE(channels,lprofiles,lsprofiles,lsprofiles2,emissivity,frequencies)
  DEALLOCATE(n_chan,polarisations,input_emissivity,calcemis)
  DEALLOCATE( transmission % tau_surf   )
  DEALLOCATE( transmission % tau_layer  )
  DEALLOCATE( transmission % od_singlelayer )
  DEALLOCATE(radiance % clear)
  DEALLOCATE( radiance % cloudy)
  DEALLOCATE( radiance % total )
  DEALLOCATE( radiance % bt )
  DEALLOCATE( radiance % bt_clear )
  DEALLOCATE( radiance % upclear )
  DEALLOCATE( radiance % dnclear )
  DEALLOCATE( radiance % reflclear )
  DEALLOCATE( radiance % overcast )
  DEALLOCATE( radiance % downcld )
  DEALLOCATE( radiance % cldemis  )
  DEALLOCATE( radiance % wtoa   )
  DEALLOCATE( radiance % wsurf   )
  DEALLOCATE( radiance % cs_wtoa  )
  DEALLOCATE( radiance % cs_wsurf )
  DEALLOCATE( radiance % out )
  DEALLOCATE( radiance % out_clear )
  DEALLOCATE( radiance % total_out )
  DEALLOCATE( radiance % clear_out )
  IF ( KRTTOVINFO(4,JSAT) == 1 .OR. KRTTOVINFO(4,JSAT) == 3) THEN
    DEALLOCATE(ZTEMPK,ZWVAPK)
    DEALLOCATE( profiles_k)
    DEALLOCATE( cld_profiles_k)
    DEALLOCATE( emissivity_k)
  ENDIF
END DO
DEALLOCATE(ZULAT,ZULON,ZANGS,IMSURF)
DEALLOCATE(ZAV,ZSAV,ZSSV,ZCV,ZAP,ZAP_HL)
#else
PRINT *, "RTTOV 8.7 LIBRARY NOT AVAILABLE = ###CALL_RTTOV8####"
#endif
!
END SUBROUTINE CALL_RTTOV8
