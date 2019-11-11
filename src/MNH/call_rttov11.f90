!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ########################
     MODULE MODI_CALL_RTTOV11
!    ########################
INTERFACE
!
     SUBROUTINE CALL_RTTOV11(KDLON, KFLEV, PEMIS, PTSRAD,   &
                PTHT, PRT, PPABST, PZZ, PMFCONV, PCLDFR, PULVLKB, PVLVLKB,  &
                OUSERI, KRTTOVINFO, TPFILE    )
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
INTEGER, INTENT(IN)   :: KDLON !number of columns where the
                               !radiation calculations are performed
INTEGER, INTENT(IN)   :: KFLEV !number of vertical levels where the
                               !radiation calculations are performed
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PEMIS  !Surface IR EMISsivity
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD !RADiative Surface Temperature
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
END SUBROUTINE CALL_RTTOV11
END INTERFACE
END MODULE MODI_CALL_RTTOV11
!    #####################################################################
SUBROUTINE CALL_RTTOV11(KDLON, KFLEV, PEMIS, PTSRAD,     &
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
!!      JP Chaboureau 30/05/2017 exclude the first layer when considering clouds
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!----------------------------------------------------------------------------
!!
!!*       0.    DECLARATIONS
!!              ------------
!!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_GRID_n
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LUNIT_n
USE MODD_DEEP_CONVECTION_n
USE MODD_REF_n
USE MODD_RADIATIONS_n,  ONLY : XSEA
!
USE MODN_CONF
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
#ifdef MNH_RTTOV_11
USE rttov_const, ONLY :  &
       & sensor_id_ir, sensor_id_hi, sensor_id_mw, &
       & q_mixratio_to_ppmv, tmin, tmax, qmin, qmax, pmin, pmax
USE rttov_types
USE parkind1, ONLY: jpim, jprb, jplm
!
IMPLICIT NONE
!
! -----------------------------------------------------------------------------
#include "rttov_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_scatt.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_init_rad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#endif
!!!
!!!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!!!
INTEGER, INTENT(IN)   :: KDLON   !number of columns where the
! radiation calculations are performed
INTEGER, INTENT(IN)   :: KFLEV   !number of vertical levels where the
! radiation calculations are performed
!!!
REAL, DIMENSION(:,:),     INTENT(IN) :: PEMIS  !Surface IR EMISsivity
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD !RADiative Surface Temperature
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
#ifdef MNH_RTTOV_11
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

REAL, DIMENSION(:,:,:), ALLOCATABLE   ::  ZBT
REAL, DIMENSION(:,:), ALLOCATABLE :: ZANTMP, ZUTH
REAL :: ZZH, zdeg_to_rad, zrad_to_deg, zbeta, zalpha

! Other arrays for zenithal solar angle 
! REAL, DIMENSION(:,:),   ALLOCATABLE :: ZCOSZEN, ZSINZEN, ZAZIMSOL

! -----------------------------------------------------------------------------
REAL, DIMENSION(:), ALLOCATABLE :: ZANGL   !Satellite zenith angle (deg)
REAL, DIMENSION(:), ALLOCATABLE :: ZANGS   !Solar zenith angle (deg)
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
LOGICAL (kind=jplm)       , ALLOCATABLE :: calcemis   (:)
INTEGER (kind=jpim)       , ALLOCATABLE :: frequencies (:)
TYPE (rttov_chanprof)     , ALLOCATABLE :: chanprof    (:) ! Channel and profile indices
TYPE (profile_type)       , ALLOCATABLE :: profiles    (:), profiles_k    (:)
TYPE (profile_cloud_type) , ALLOCATABLE :: cld_profiles(:), cld_profiles_k(:)
TYPE(rttov_emissivity),  ALLOCATABLE :: emissivity(:)  ! Input/output surface emissivity
LOGICAL(KIND=jplm),      ALLOCATABLE :: calcrefl(:)    ! Flag to indicate calculation of BRDF within RTTOV
TYPE(rttov_reflectance), ALLOCATABLE :: reflectance(:) ! Input/output surface BRDF
  TYPE(transmission_type)              :: transmission   ! Output transmittances
  INTEGER(KIND=jpim) :: asw

integer (kind=jpim)        :: errorstatus
type (radiance_type)       :: radiance, radiance_k  
type (rttov_options)       :: opts     ! Defaults to everything optional switched off
type (rttov_options_scatt) :: opts_scatt
type (rttov_coefs     )    :: coef_rttov
type (rttov_scatt_coef)    :: coef_scatt

integer (kind=jpim) :: instrument (3)
integer (kind=jpim) :: ilev, iprof, ichan, nprof, nchan, nlev, nchannels
real    (kind=jprb) :: zenangle
integer (kind=jpim), parameter :: fin = 10
character (len=256) :: outstring
! -----------------------------------------------------------------------------
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZTEMP
TYPE(TFIELDDATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
!*       0.     ARRAYS BOUNDS INITIALIZATION
!
IIU=SIZE(PTHT,1)
IJU=SIZE(PTHT,2)
IKU=SIZE(PTHT,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=IKU-JPVEXT

errorstatus     = 0
nlev=IKE-IKB+1
nprof=1
ZTEMP = PTHT * ( PPABST/XP00 ) ** (XRD/XCPD)                   
DO JSAT=1,SIZE(KRTTOVINFO,2)
  IF (KRTTOVINFO(1,JSAT) /= NUNDEF) THEN
    IJSAT = JSAT
  END IF
END DO

! -----------------------------------------------------------------------------
!              *** LOOP OVER SENSORS ***
! -----------------------------------------------------------------------------
DO JSAT=1,IJSAT ! loop over sensors
 
  instrument(1)=KRTTOVINFO(1,JSAT)
  instrument(2)=KRTTOVINFO(2,JSAT)
  instrument(3)=KRTTOVINFO(3,JSAT)
! PRINT *,' JSAT=',JSAT, instrument

!!! METEOSAT, GOES, OR MSG PLATFORM
  IF (KRTTOVINFO(1,JSAT) == 3 .OR. KRTTOVINFO(1,JSAT) == 4 &
       .OR. KRTTOVINFO(1,JSAT) == 12) THEN
    opts % rt_ir % addsolar = .FALSE.         ! Do not include solar radiation
    opts % interpolation % addinterp  = .TRUE.  ! Allow interpolation of input profile
    opts % interpolation % interp_mode = 1       ! Set interpolation method
    opts % rt_all % addrefrac          = .FALSE.  ! Do not include refraction in path calc
    opts % rt_ir % addclouds          = .TRUE.  ! Include cloud effects
    opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects
    opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
    opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
    opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
    opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
    opts % rt_ir % co_data             = .FALSE. !
!   opts % rt_mw % clw_data            = .FALSE. !
!   opts%rt_ir%user_cld_opt_param = .FALSE.
  ELSE
    opts % rt_ir % addclouds          = .FALSE.  ! Include cloud effects
  END IF
  opts % config % verbose            = .FALSE.  ! Enable printing of warnings
  opts % config % do_checkinput      = .FALSE.


! Read and initialise coefficients
! -----------------------------------------------------------------------------
  CALL rttov_read_coefs (errorstatus, coef_rttov, opts, instrument=instrument)
  IF(errorstatus /= 0) THEN
    WRITE(*,*) 'error rttov_readcoeffs :',errorstatus
    CALL PRINT_MSG(NVERB_FATAL,'GEN','CALL_RTTOV11','error rttov_readcoeffs')
  ENDIF
!  CALL rttov_initcoeffs (errorstatus,coef_rttov)
!  IF( errorstatus/= 0) THEN
!    WRITE(*,*) 'error rttov_initcoeffs :',errorstatus
!    CALL PRINT_MSG(NVERB_FATAL,'GEN','CALL_RTTOV11','error rttov_initcoeffs')
!  ENDIF
  
! Read coef file for cloud/rain absorption/scattering
  IF( coef_rttov%coef%id_sensor == sensor_id_mw) THEN
    CALL rttov_read_scattcoeffs (errorstatus, coef_rttov%coef, coef_scatt)
  END IF

  nchan     = coef_rttov%coef%fmv_chn   ! number of channels on instrument
  nchannels = nprof * nchan             ! total channels to simulate

  ALLOCATE(ZBT(IIU,IJU,nchannels))
  ZBT(:,:,:)=999.
!  PRINT *,'ncan=',nchan,' nchannels=',nchannels

  ALLOCATE (chanprof     (nchannels))
  ALLOCATE (frequencies  (nchannels))
  ALLOCATE (emissivity   (nchannels))
  ALLOCATE (calcemis    (nchannels))
  ALLOCATE (profiles     (nprof))
  ALLOCATE (cld_profiles (nprof))
! Request RTTOV / FASTEM to calculate surface emissivity
  calcemis  = .TRUE.
  emissivity % emis_in = 0.0_JPRB

!!! METEOSAT, GOES, OR MSG PLATFORM
  IF (KRTTOVINFO(1,JSAT) == 3 .OR. KRTTOVINFO(1,JSAT) == 4 &
       .OR. KRTTOVINFO(1,JSAT) == 12) calcemis = .FALSE.

!  IF( coef_rttov%coef% id_sensor /= sensor_id_mw) THEN
!  ! Allocate arrays for surface reflectance
!    ALLOCATE(calcrefl(nchannels))
!    ALLOCATE(reflectance(nchannels))
!  END IF

! Setup indices
  IF( coef_rttov%coef% id_sensor /= sensor_id_mw) THEN
    DO JCH=1,nchannels
      chanprof(JCH)%prof = 1
      chanprof(JCH)%chan = JCH
    END DO
  ELSE
    CALL rttov_scatt_setupindex ( &
         & nprof,           & ! in  
         & nchan,           & ! in 
         & coef_rttov%coef, & ! in
         & nchannels,       & ! in
         & chanprof,        & ! out
       & frequencies)       ! out
  END IF

  asw = 1_jpim ! Switch for allocation passed into RTTOV subroutines

! Allocate profiles (input) and radiance (output) structures
  CALL rttov_alloc_prof(errorstatus, nprof,     profiles, nlev, opts,asw,coef_rttov,init = .TRUE._jplm)
  IF( coef_rttov%coef% id_sensor == sensor_id_mw) THEN
! CALL rttov_alloc_opt_param( &
!     & errorstatus,          &
!     & cld_opt_param,        &
!     & nchanprof,            &
!     & nlevels-1_jpim,       &
!     & nphangle,             &
!     & asw)
! ELSE
    CALL rttov_alloc_scatt_prof(      nprof, cld_profiles, nlev, .FALSE._jplm, 1_jpim, init = .TRUE._jplm)
  END IF

  CALL rttov_alloc_rad       (errorstatus, nchannels, radiance, nlev-1_jpim,asw)
!    WRITE(*,*) 'error rttov_alloc_rad :',errorstatus
  ! Allocate transmittance structure
  CALL rttov_alloc_transmission( &
      & errorstatus,             &
      & transmission,            &
      & nlev-1_jpim,          &
      & nchannels,               &
      & asw,                     &
      & init=.TRUE.)

  profiles(1) % zenangle    = 0. ! zenith
  cld_profiles(1) % use_totalice = .FALSE.
  profiles(1) % skin % fastem(:) = &
! RTTOV 8.5 example
!        (/ 3.0_JPRB, 5.0_JPRB, 15.0_JPRB, 0.1_JPRB, 0.3_JPRB /)
! Bare soil see Table 3 svr rttov7)
       (/ 2.3_JPRB, 1.9_JPRB, 21.8_JPRB, 0.0_JPRB, 0.5_JPRB /)

  profiles(1) % nlevels =  nlev
  profiles(1) % nlayers =  nlev-1

 ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
  IF (errorstatus /= 0) THEN
    WRITE(*,*) 'error in rttov options'
    STOP
  ENDIF

!!  opts%interpolation%reg_limit_extrap = .TRUE.
!!  profiles(1)%gas_units = 1 ! kg/kg over moist air
!PRINT *,'nlev=',nlev,' tmax=',tmax,' tmin=',tmin,' qmax=',qmax,' qmin=',qmin
!PRINT *, coef_rttov%coef % nlevels
  DO JI=IIB,IIE
    DO JJ=IJB,IJE      
      DO JK=IKB,IKE ! nlevels
        JKRAD = nlev-JK+2 !INVERSION OF VERTICAL LEVELS!
!PRINT *,'jk=',jk,' jkrad=',jkrad
        profiles(1) % p(JKRAD) = PPABST(JI,JJ,JK)*0.01
        profiles(1) % t(JKRAD) = MIN(tmax,MAX(tmin,ZTEMP(JI,JJ,JK)))
!PRINT *,'jk=',JK,' ZTEMP=',ZTEMP(JI,JJ,JK),' t=',profiles(1) % t(JKRAD)
        profiles(1) % q(JKRAD) = MIN(qmax,MAX(qmin,PRT(JI,JJ,JK,1)*q_mixratio_to_ppmv))
!        PRINT *,JK,profiles(1) % p(JKRAD) ,profiles(1) % t(JKRAD) ,profiles(1) % q(JKRAD) 
      END DO
      profiles(1) % elevation = 0.5*( PZZ(JI,JJ,1)+PZZ(JI,JJ,IKB) )
      profiles(1) % skin % t = MIN(tmax,MAX(tmin,PTSRAD(JI,JJ)))
      profiles(1) % s2m % t = MIN(tmax,MAX(tmin,ZTEMP(JI,JJ,IKB)))
      profiles(1) % s2m % q = MIN(qmax,MAX(qmin,PRT(JI,JJ,1,IKB)*q_mixratio_to_ppmv))
      profiles(1) % s2m % u = PULVLKB(JI,JJ) ! 2m wind speed u (m/s)
      profiles(1) % s2m % v = PVLVLKB(JI,JJ) ! 2m wind speed v (m/s)
      profiles(1) % s2m % p = PPABST(JI,JJ,IKB)*0.01
      IF (NINT(XSEA(JI,JJ)).EQ.0.) THEN
        profiles(1) % skin % surftype = 0 ! Surface Mask 0=land, 1=sea, 2=sea-ice
      ELSE
        profiles(1) % skin % surftype = 1
        profiles(1) % skin % watertype = 1 ! Ocean water
      END IF
      profiles(1) % ctp = 500.0_JPRB     ! Not used but still required by RTTOV
      IF( coef_rttov%coef% id_sensor /= sensor_id_mw) THEN
        profiles(1)%ish = 2 ! Aggregates
        profiles(1)%idg = 4 ! McFarquar et al (2003)
        DO JK=IKB+1,IKE-1 ! nlayers
          JKRAD = nlev-JK+1 !INVERSION OF VERTICAL LEVELS!
          profiles(1) %cfrac(JKRAD) = PCLDFR(JI,JJ,JK)
          profiles(1) %cloud(1,JKRAD) = PRT(JI,JJ,JK,2)*XRHODREF(JI,JJ,JK)*1.0E03
          IF( OUSERI ) THEN
            profiles(1) %cloud(6,JKRAD) = (PRT(JI,JJ,JK,4)+PRT(JI,JJ,JK,5)) \
            *XRHODREF(JI,JJ,JK)*1.0E03
          END IF
        END DO
      ELSE
        DO JK=IKB,IKE
          JKRAD = nlev-JK+2 !INVERSION OF VERTICAL LEVELS!
          cld_profiles(1) % ph (JKRAD) = 0.5*( PPABST(JI,JJ,JK) + PPABST(JI,JJ,JK+1) )*0.01
          cld_profiles(1) %cc(JKRAD) = PCLDFR(JI,JJ,JK)
          cld_profiles(1) %clw(JKRAD) = MIN(ZRCMAX,PRT(JI,JJ,JK,2))
          cld_profiles(1) %rain(JKRAD) = MIN(ZRRMAX,PRT(JI,JJ,JK,3))
          IF( OUSERI ) THEN
            cld_profiles(1) %ciw(JKRAD) = MIN(ZRIMAX,PRT(JI,JJ,JK,4))
            cld_profiles(1) %sp(JKRAD) = MIN(ZRSMAX,PRT(JI,JJ,JK,5)+PRT(JI,JJ,JK,6))
          END IF
        END DO
        cld_profiles (1) % ph (nlev+1) =   profiles (1) % s2m % p
!          PRINT *,nlev+1,' cld_profiles(1) % ph (nlev+1) =',cld_profiles(1) % ph (nlev+1) 
      END IF

      DO JCH=1,nchannels
        IF (.NOT.calcemis(JCH)) emissivity(JCH)%emis_in = PEMIS(JI,JJ)
      END DO

!write(*,*) 'Calling forward model' 

! Forward model run
      IF ( coef_rttov%coef% id_sensor /= sensor_id_mw) THEN
        CALL rttov_direct(                &
             & errorstatus,              &! out   error flag
             & chanprof,                 &! in    channel and profile index structure
             & opts,                     &! in    options structure
             & profiles,                 &! in    profile array
             & coef_rttov,               &! in    coefficients strucutre
             & transmission,             &! inout compscauted transmittances
             & radiance,                 &! inout computed radiances
             & calcemis    = calcemis,   &! in    flag for internal emissivity calcs
             & emissivity  = emissivity) !, &! inout input/output emissivities per channel
!             & calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
!             & reflectance = reflectance) ! inout input/output BRDFs per channel
      ELSE
        CALL rttov_scatt ( &
             & errorstatus,         &! out
             & opts_scatt,          &! in
             & nlev,                &! in
             & chanprof,            &! in
             & frequencies,         &! in
             & profiles,            &! in  
             & cld_profiles,        &! in
             & coef_rttov,          &! in
             & coef_scatt,          &! in
             & calcemis,           &! in
             & emissivity,          &! in
             & radiance)             ! out
      END IF
! STOP
      DO JCH=1,nchannels
        ZBT(JI,JJ,JCH)= radiance % bt(JCH)
      END DO
    END DO
  END DO
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

  DO JCH=1,nchannels
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

!    IF (INRAD==1) THEN
!    TZFIELD%CMNHNAME   = TRIM(YBEG)//'_'//TRIM(YEND)//'rad'
!    TZFIELD%CUNITS     = 'mw/cm-1/ster/sq.m'
!    TZFIELD%CCOMMENT   = TRIM(YBEG)//'_'//TRIM(YEND)//' rad'
!    ELSE
    TZFIELD%CMNHNAME   = TRIM(YBEG)//'_'//TRIM(YEND)//'BT'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CCOMMENT   = TRIM(YBEG)//'_'//TRIM(YEND)//' BT'
!    ENDIF
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .FALSE.
!    PRINT *,'YRECFM='//TRIM(TZFIELD%CMNHNAME)
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZBT(:,:,JCH))
  END DO
  DEALLOCATE(chanprof,frequencies,emissivity,calcemis,profiles,cld_profiles)
  DEALLOCATE(ZBT)
!  IF( coef_rttov%coef% id_sensor /= sensor_id_mw) THEN
!    DEALLOCATE(calcrefl,reflectance)
!  END IF
END DO

#else
PRINT *, "RTTOV 11.1 LIBRARY NOT AVAILABLE = ###CALL_RTTOV11####"
#endif
!
END SUBROUTINE CALL_RTTOV11
