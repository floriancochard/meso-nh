!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      PROGRAM CH_MODEL0D
!!    ##################
!!
!!***  *CH_MODEL0D*
!!
!!    PURPOSE
!!    -------
!     Monitor of the zero-dimensional box model
!!
!!**  METHOD
!!    ------
!!    This monitor makes all necessary initial calls,
!!    controls the timestep variables of the solver and the
!!    diagnostic routines and calls the different utilities
!!    (e.g. to write results to disk at a given interval)
!!
!!    REFERENCES
!!    ----------
!!    none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!    01/08/96 (K. Suhre) fix namelist filename to CHCONTROL1.nam
!!    16/02/99 (K. Suhre) add on-line photolysis rates TUV
!!    16/02/99 (V. Crassier & K. Suhre) insert vectorization
!!    18/02/99 (K. Suhre) take out timing commands
!!    23/03/99 (K. Suhre) update reaction rates at every time step
!!                        (for liquid phase chemistry)
!!    21/09/04 (P. Tulet) update for MASDEV44 bug2
!!    21/03/06 (P. Tulet) update for MASDEV46 and add ORILAM aerosol scheme
!!    24/24/14 (M. Leriche) add ReLACS3
!!    M.Leriche 2015 : masse molaire Black carbon à 12 g/mol
!!    Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_INIT_MODEL0D
USE MODI_CH_INIT_CCS
USE MODI_CH_READ_CHEM
USE MODI_CH_SHOW_CHEM
USE MODI_CH_INIT_METEO
USE MODI_CH_UPDATE_METEO
USE MODI_CH_INIT_OUTPUT
USE MODI_CH_OUTPUT
USE MODI_CH_DIAGNOSTICS
USE MODI_CH_SET_RATES
USE MODI_CH_SET_PHOTO_RATES
USE MODI_CH_SOLVER_n
USE MODI_CH_WRITE_CHEM
USE MODI_CH_UPDATE_JVALUES
USE MODI_CH_INIT_SCHEME
USE MODI_CH_ORILAM
USE MODI_CH_INI_ORILAM
USE MODI_CH_AER_EQM_INIT0d
USE MODD_CH_AEROSOL
!USE MODD_CH_AERO_n
USE MODI_CH_AER_SURF
USE MODI_CH_INIT_JVALUES
USE MODI_CH_AER_INIT_SOA
USE MODI_CH_AER_MOD_INIT
USE MODD_CST, ONLY : XMNH_TINY
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_INIT_JVALUES
USE MODD_CH_MODEL0D
USE MODD_CH_M9_n,     ONLY: NEQ, NREAC,      & ! no. of species & reactions
                            NMETEOVARS,      & ! no. of meteo variables
                            CNAMES,          & ! names of chem. species
                            METEOTRANSTYPE     ! TYPE of meteo struct variable

USE MODD_CONF, ONLY : CPROGRAM 
USE MODD_CH_M9_SCHEME, ONLY :  CCSTYPE,TACCS

USE MODE_IO_ll
USE MODE_MODELN_HANDLER
!!
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
CHARACTER*256      :: YNAMELISTFILE = "CHCONTROL1.nam" !  namelist input file
!
! reaction rates and auxiliary variables  
TYPE(CCSTYPE), POINTER        :: TZK 
TYPE(METEOTRANSTYPE), DIMENSION(1) :: TZM 
!
! meteo variables to be transferred into CCS
!
REAL, DIMENSION(1,1,1)         :: ZRHODREF ! air density
REAL, DIMENSION(1)             :: ZDENAIR ! air density
REAL, DIMENSION(1)             :: ZTHT     ! potential temperature (K)
REAL, DIMENSION(1)             :: ZTEMP    ! air temperature (K)
REAL, DIMENSION(1)             :: ZPABST   ! pressure
REAL, DIMENSION(1,1,1)         :: ZZZ      ! height of each model level
REAL, DIMENSION(1,1,1,2)       :: ZRT      ! moist variables (1= vapor, 2= cloud water)
REAL, DIMENSION(1)             :: ZRC      ! vapor water
REAL, DIMENSION(1)             :: ZRV      ! cloud water
REAL, DIMENSION(1,1)           :: ZZS      ! height of surface 
REAL, DIMENSION(1,1)           :: ZALBUV   ! surface UV aldedo

! Variables for TUV (to compute solar zenith angle)
INTEGER, DIMENSION(0:11) :: IBIS, INOBIS ! Cumulative number of days per month
                                         ! for bissextile and regular years
REAL  :: ZTIME, ZTIMERAD, ZUT, ZTUT, ZDATE, ZAD, ZTSIDER, ZA1, ZA2 ! variables for solar angle
REAL  :: ZLAT, ZLON   ! latitude and longitude in radians
REAL  :: XLAT, XLON ! latitude and longitude in degrees
REAL  :: ZCOSZEN, ZSOLANG
REAL  :: ZSINDEL, ZCOSDEL, ZDECSOL
REAL, DIMENSION(1,1)           :: ZZENITH  ! solar zenith angle


REAL, DIMENSION(1,1,1,JPJVMAX) :: XJVALUES ! TUV coeff of photodissociation
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZCONC, ZAERO    
                                  ! prognostic chem./aerosol species at time XTSIMUL
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZNEWCONC , ZNEWAERO
                                  ! prognostic chem./aerosol species after one timestep
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZFRAC 
!
INTEGER               :: ITCOUNT = 0   ! loop counter for JVALUE update
INTEGER, DIMENSION(4,1) :: KMASK         ! mask to stop microphysic at
                                         ! low aerosol values
REAL, DIMENSION(1,JPIN)   :: ZM          ! aerosol moments (JPIN = 3*JPMODE)
REAL, DIMENSION(1,JPMODE) :: ZSIG0       ! log of standard deviation (LN(SIG)) per mode
REAL, DIMENSION(1,JPMODE) :: ZN0         ! particles number pr mode
REAL, DIMENSION(1,JPMODE) :: ZRG0        ! mean radius (in um) per mode
REAL, DIMENSION(1,JPMODE) :: ZRHOP0      ! aerosol density per mode
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZCTOTG     ! condensable gas (microgram/m3)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZCTOTA ! aerosol composition per mode  (microgram/m3)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZCCTOT ! fraction of aerosol composition per mode  
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZMI ! molecular mass of each aerosol compounds
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZSOLORG
REAL, DIMENSION(1,JPIN)   :: ZSEDA       ! aerosol moment tendency for
                                         ! sedimentation (not use in 0D equal to 0)
REAL, DIMENSION(1)        :: ZMU         ! Gas viscosity (kg/(ms))
REAL, DIMENSION(1)        :: ZLAMBDA     ! Mean free path of background
REAL, DIMENSION(1,JPMODE) :: ZOM         ! Omega factor
REAL, DIMENSION(1)        :: ZSO4RAT     ! sulfuric acid tendency
! 
! Constant values
REAL :: ZVALBC, ZVALOC,  ZCOEFAEROBC, ZCOEFAEROOC
REAL :: ZP00, ZRD, ZCPD, ZPI, ZAVOGADRO, ZBOLTZ, ZMD
INTEGER :: IMONTH, IYEAR, IDAY
REAL :: ZDEN2MOL
CHARACTER(LEN=10)    :: CCH_SCHEME

! Index
INTEGER :: JI, JJ, JN

!------------------------------------------------------------------------------
!
CALL GOTO_MODEL(1)

! 0)  Initialization of constants

CPROGRAM = "MOD0D"

ZPI = 2.*ASIN(1.) ! Pi
ZP00 = 1.E5 ! Reference pressure
ZAVOGADRO = 6.0221367E+23
ZBOLTZ    = 1.380658E-23
ZMD    = 28.9644E-3
ZRD    = ZAVOGADRO * ZBOLTZ / ZMD 
ZCPD   = 7.* ZRD /2.
ZDEN2MOL = 1E-6 * ZAVOGADRO / ZMD

XCH_TUV_ALBNEW = 0.02  ! surface UV albedo (to be modified)
XCH_TUV_DOBNEW = 320.  ! O3 dobson (to be modified)
!
!*       1.   INITIALISATION
!        -------------------
CALL INITIO_ll()
!
!*       1.1  read namelist and initialize time control variables
!
! read name of namelist file from standard input
PRINT *, "This is the chemistry box model Version 1.2"
PRINT *, "Run parameters will be read from namelistfile file ", YNAMELISTFILE
CALL CH_INIT_MODEL0D(YNAMELISTFILE)
! 
!*       1.2  initialize chemical core system (CCS)
!
CALL CH_INIT_SCHEME(6)
IF (.NOT. ASSOCIATED(TACCS(1)%NVERB)) THEN
   CALL CH_ALLOCATE_TACCS(1,1)
 END IF
TZK=>TACCS(1)
KMASK(:,:)=1
CALL CH_INIT_CCS(NEQ,6,NVERB) ! output channel of the box-model is stdout (6)

! Allocation of gaseous chemical variables 
ALLOCATE(ZCONC(1,NEQ))
ALLOCATE(ZNEWCONC(1,NEQ))
ALLOCATE(ZFRAC(1,NEQ))
!

IF (LORILAM) THEN
! Initialisation of SOA
  CALL CH_AER_INIT_SOA(6, NVERB)
! Initialisation of aerosols tables
  CALL CH_AER_MOD_INIT
  DO JN=1, SIZE(CNAMES)
  IF (TRIM(CNAMES(JN)) .EQ. "ALKA") CCH_SCHEME = "RELACS"
  IF (TRIM(CNAMES(JN)) .EQ. "HC3")  CCH_SCHEME = "RACM"
  IF (TRIM(CNAMES(JN)) .EQ. "URG1") CCH_SCHEME = "RELACS2"
  IF (TRIM(CNAMES(JN)) .EQ. "GLY")  CCH_SCHEME = "RELACS3"
  IF (TRIM(CNAMES(JN)) .EQ. "UR29") CCH_SCHEME = "CACM"
  ENDDO


! Allocation of aerosol chemical variables 
IF (LVARSIGI.AND.LVARSIGJ) ALLOCATE(ZAERO(1,(NSP+NCARB+NSOA+2)*JPMODE))
IF (LVARSIGI.AND..NOT.LVARSIGJ) ALLOCATE(ZAERO(1,(NSP+NCARB+NSOA+1)*JPMODE+1))
IF (.NOT.LVARSIGI.AND.LVARSIGJ) ALLOCATE(ZAERO(1,(NSP+NCARB+NSOA+1)*JPMODE+1))
IF (.NOT.LVARSIGI.AND..NOT.LVARSIGJ) ALLOCATE(ZAERO(1,(NSP+NCARB+NSOA+1)*JPMODE))
ALLOCATE(ZNEWAERO(1,SIZE(ZAERO,2)))
ALLOCATE(XSURF(1,JPMODE))       ! aerosol surface per mode
ALLOCATE(XDP(1,JPMODE))         ! aerosol diameter per mode
ALLOCATE(ZCTOTG(1,NSP+NCARB+NSOA))        ! condensable gas (microgram/m3)
ALLOCATE(ZCTOTA(1,NSP+NCARB+NSOA,JPMODE)) ! aerosol composition per mode  (microgram/m3)
ALLOCATE(ZCCTOT(1,NSP+NCARB+NSOA,JPMODE)) ! fraction of aerosol composition per mode  
ALLOCATE(ZMI(1,NSP+NCARB+NSOA)) ! molecular mass in g/mol-1
ALLOCATE(ZSOLORG(1,10)) ! solubility of soa (mass fraction)

! Constants initialization
ZMI(:,:) = 250.
ZMI(1,JP_AER_SO4)  = 98.
ZMI(1,JP_AER_NO3)  = 63.
ZMI(1,JP_AER_NH3)  = 17.
ZMI(1,JP_AER_H2O)  = 18.
ZMI(1,JP_AER_BC)   = 12.
IF (NSOA == 10) THEN
ZMI(1,JP_AER_SOA1) = 88. 
ZMI(1,JP_AER_SOA2) = 180.
ZMI(1,JP_AER_SOA3) = 1.5374857E+02
ZMI(1,JP_AER_SOA4) = 1.9586780E+02
ZMI(1,JP_AER_SOA5) = 1951.
ZMI(1,JP_AER_SOA6) = 195.
ZMI(1,JP_AER_SOA7) = 165.
ZMI(1,JP_AER_SOA8) = 195.
ZMI(1,JP_AER_SOA9) = 270.
ZMI(1,JP_AER_SOA10) = 210.
END IF
ELSE
ALLOCATE(ZAERO(1,0))
ALLOCATE(ZNEWAERO(1,SIZE(ZAERO,2)))
ALLOCATE(XSURF(1,0))       ! aerosol surface per mode
ALLOCATE(XDP(1,0))         ! aerosol diameter per mode
ALLOCATE(ZCTOTG(1,0))        ! condensable gas (microgram/m3)
ALLOCATE(ZCTOTA(1,0,0)) ! aerosol composition per mode  (microgram/m3)
ALLOCATE(ZCCTOT(1,0,0)) ! fraction of aerosol composition per mode  
ALLOCATE(ZMI(1,0))
END IF

!*       1.3  prepare reading of meteo-variables and read first set of data

CALL CH_INIT_METEO(TZM(1))
CALL CH_UPDATE_METEO(TZM(1),XTSIMUL)
TZM(:)=TZM(1)

ZALBUV= XCH_TUV_ALBNEW
ZZZ(1,1,1) = TZM(1)%XMETEOVAR(1)
ZZS(:,:) = 0.
ZRHODREF(1,1,1)= TZM(1)%XMETEOVAR(2)
ZDENAIR(1) = ZRHODREF(1,1,1)
ZPABST(1)  = TZM(1)%XMETEOVAR(2) * TZM(1)%XMETEOVAR(3) * 288.290947
ZTHT(1) = TZM(1)%XMETEOVAR(3) / ((ZPABST(1)/ZP00)**(ZRD/ZCPD))
ZRT(1,1,1,1)= TZM(1)%XMETEOVAR(4)
ZRT(1,1,1,2)= TZM(1)%XMETEOVAR(5)
!Temperature (K)
ZTEMP(1)   = ZTHT(1)*((ZPABST(1)/ZP00)**(ZRD/ZCPD))
!Water vapor (kg/kg)
ZRV(1)     = ZRT(1,1,1,1)
!Cloud vapor (kg/kg)
ZRC(1)     = ZRT(1,1,1,2)

XLAT = TZM(1)%XMETEOVAR(6)
XLON = TZM(1)%XMETEOVAR(7)
IYEAR = INT(TZM(1)%XMETEOVAR(8))
IMONTH = INT(TZM(1)%XMETEOVAR(9))
IDAY = INT(TZM(1)%XMETEOVAR(10))
! 
!*       1.4  initialize vector of chemical species
!
CALL CH_READ_CHEM(ZCONC(1,:), ZAERO(1,:), ZRHODREF(:,1,1), CINITFILE)
!
IF (NVERB >= 5) THEN
  PRINT *, NEQ, ' chemical concentrations initialized to:'
  CALL CH_SHOW_CHEM(ZCONC(1,:),CNAMES)
  IF (LORILAM) THEN
    PRINT *, SIZE(CAERONAMES), ' aerosol concentrations initialized to:'
    CALL CH_SHOW_CHEM(ZAERO(1,:),CAERONAMES)
  END IF
END IF
! 
ZCONC(1,:) = ZCONC(1,:) * ZDEN2MOL * ZRHODREF(1,1,1) ! convert ppp to molec.cm-3
ZNEWCONC(1,:) = ZCONC(1,:)
IF (LORILAM) THEN
ZAERO(1,:) = ZAERO(1,:) * ZDEN2MOL * ZRHODREF(1,1,1) ! convert ppp to molec.cm-3
ZNEWAERO(1,:) = ZAERO(1,:)  
END IF
!*       1.5  initialize data for jvalues

CALL CH_INIT_JVALUES(IDAY, IMONTH, IYEAR, 6, XCH_TUV_DOBNEW)

ZCONC(:,:)  = MAX(ZCONC(:,:), XMNH_TINY)
IF (LORILAM) THEN
ZSIG0(:,1) = LOG(XINISIGI)
ZSIG0(:,2) = LOG(XINISIGJ)

IF (CRGUNIT=="MASS") THEN
ZRG0(:,1) = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
ZRG0(:,2) = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
ZRG0(:,1) = XINIRADIUSI 
ZRG0(:,2) = XINIRADIUSJ
ENDIF

ZAERO(:,:)  = MAX(ZAERO(:,:), XMNH_TINY)
END IF


IF (LORILAM) CALL CH_AER_EQM_INIT0d(ZMI, ZAERO,ZM, ZRHOP0, ZSIG0, ZRG0, ZN0, ZCTOTA)

!*       1.6    COMPUTES THE DAY OF THE YEAR
!
INOBIS(:) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
IBIS(0) = INOBIS(0)
DO JI=1,11
  IBIS(JI) = INOBIS(JI)+1
END DO

ZLAT = XLAT*(ZPI/180.)
ZLON = XLON*(ZPI/180.)
!
ZTIMERAD     = XTSIMUL + 0.5*XCH_TUV_TUPDATE
IBIS(0) = INOBIS(0)
DO JI=1,11
  IBIS(JI) = INOBIS(JI)+1
END DO
IF( MOD(IYEAR,4).EQ.0 ) THEN
  ZDATE = FLOAT(IDAY +   IBIS(IMONTH-1)) - 1
  ZAD = 2.0*ZPI*ZDATE/366.0
ELSE
  ZDATE = FLOAT(IDAY + INOBIS(IMONTH-1)) - 1
  ZAD = 2.0*ZPI*ZDATE/365.0
END IF
ZDECSOL = 0.006918-0.399912*COS(ZAD)   +0.070257*SIN(ZAD)    &
         -0.006758*COS(2.*ZAD)+0.000907*SIN(2.*ZAD) &
         -0.002697*COS(3.*ZAD)+0.00148 *SIN(3.*ZAD)
ZSINDEL = SIN(ZDECSOL)
ZCOSDEL = COS(ZDECSOL)
ZA1 = (1.00554*ZDATE- 6.28306)*(ZPI/180.0)
ZA2 = (1.93946*ZDATE+23.35089)*(ZPI/180.0)
ZTSIDER = (7.67825*SIN(ZA1)+10.09176*SIN(ZA2)) / 60.0
         print*,'DMOD(XTSIMUL/3600.,24.0) =',DMOD(XTSIMUL/3600.,24.0)
ZTIME = XTSIMUL
ZUT   = DMOD( 24.0 + DMOD(ZTIME/3600.,24.0),24.0)
ZTUT = ZUT - ZTSIDER + ZLON*((180./ZPI)/15.0)
ZSOLANG = (ZTUT-12.0)*15.0*(ZPI/180.)          ! hour angle in radians
ZCOSZEN = SIN(ZLAT)*ZSINDEL +                 &! Cosine of the zenithal
              COS(ZLAT)*ZCOSDEL*COS(ZSOLANG)  !       solar angle
ZZENITH(:,:) = ACOS(ZCOSZEN)

!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTES THE SOLAR DECLINATION ANGLE
!	        ------------------------------------
!

!*       2.1  set reaction and photolysis rates
!   conversion from part/part to molec./cm3


IF (LORILAM) THEN
!   Compute aerosol surface (in m2 / cc)
   XSURF(:,1) =  ZPI * XN0IMIN * (ZRG0(:,1))**2 * EXP(2.*ZSIG0(:,1)**2)
   XSURF(:,2) =  ZPI * XN0JMIN * (ZRG0(:,2))**2 * EXP(2.*ZSIG0(:,2)**2)
      ! Mean diameter in meter 
   XDP(:,:) = 2.E-6 * ZRG0(:,:)
END IF

CALL CH_SET_RATES(XTSIMUL,ZCONC,TZM,1,6,NVERB,1,NEQ)

TZK%MODELLEVEL = 1
!

CALL CH_UPDATE_JVALUES(6, ZZENITH, ZRT,              &
       ZALBUV, ZZS, ZZZ, XLAT, XLON,                &
       SIZE(ZZZ,1), SIZE(ZZZ,2), SIZE(ZZZ,3), SIZE(ZRT,4), &
       IDAY, IMONTH, IYEAR, XTSIMUL,&
       LCH_TUV_ONLINE, CCH_TUV_CLOUDS, &
       XCH_TUV_ALBNEW, XCH_TUV_DOBNEW, ZRHODREF, XJVALUES,&
       1,1,1,1,3,3,10)


CALL CH_SET_PHOTO_RATES(XTSIMUL,ZCONC,1,TZM,1,6,NVERB,1,KMASK,NEQ, XJVALUES )

IF (LORILAM) THEN


          CALL CH_INI_ORILAM(ZM, ZSIG0, ZRG0, ZN0, ZCTOTG, ZCTOTA, ZCCTOT,&
                             ZSEDA, ZOM, ZRHOP0, ZAERO, ZCONC, ZRV, ZDENAIR, &
                             ZPABST, ZTEMP, ZRC, ZFRAC, ZMI, CCH_SCHEME )

ZCONC(:,JP_CH_H2SO4) =  ZCONC(:,JP_CH_H2SO4) + ZAERO(:,JP_CH_SO4i) + ZAERO(:,JP_CH_SO4j)

END IF

!*       2.2  prepare result output
!
CALL CH_INIT_OUTPUT(TZM(1))
! 
!*       2.3  prepare diagnostics output
!
CALL CH_INIT_DIAGNOSTICS
!
!*       3.   TIME LOOP FOR INTEGRATION IN TIME
!        --------------------------------------
!
!PRINT *, "CH_MODEL0D: simulation started, XTSIMUL = ", XTSIMUL
time_loop: DO
! 
!*       3.1  update meteorological data and reaction rates
!
CALL CH_UPDATE_METEO(TZM(1),XTSIMUL)
TZM(:)=TZM(1)
ZZZ(1,1,1) = TZM(1)%XMETEOVAR(1)
ZRHODREF(1,1,1)= TZM(1)%XMETEOVAR(2)
ZDENAIR(1) = ZRHODREF(1,1,1)
ZPABST(1)  = TZM(1)%XMETEOVAR(2) * TZM(1)%XMETEOVAR(3) * 288.290947 
ZTHT(1) = TZM(1)%XMETEOVAR(3) / ((ZPABST(1)/ZP00)**(ZRD/ZCPD))
ZRT(1,1,1,1)= TZM(1)%XMETEOVAR(4)
ZRT(1,1,1,2)= TZM(1)%XMETEOVAR(5)
XLAT = TZM(1)%XMETEOVAR(6)
XLON = TZM(1)%XMETEOVAR(7)
IYEAR  = INT(TZM(1)%XMETEOVAR(8))
IMONTH = INT(TZM(1)%XMETEOVAR(9))
IDAY   = INT(TZM(1)%XMETEOVAR(10))
!Temperature (K)
ZTEMP(1)   = ZTHT(1)*((ZPABST(1)/ZP00)**(ZRD/ZCPD))
!Water vapor (kg/kg)
ZRV(1)     = ZRT(1,1,1,1)
!Cloud vapor (kg/kg)
ZRC(1)     = ZRT(1,1,1,2)
!
IF (LORILAM) THEN
!   Compute aerosol surface (in m2 / cc)
   CALL CH_AER_SURF(ZM, ZRG0, ZSIG0, XSURF)
      ! Mean diameter in meter 
   XDP(:,:) = 2.E-6 * ZRG0(:,:)
END IF

  ! ZCONC(1,:) = ZCONC(1,:) * ZDEN2MOL * ZRHODREF(1,1,1) !convert ppp to  molec.cm-3

   CALL CH_SET_RATES(XTSIMUL,ZCONC,TZM,1,6,NVERB,1,NEQ)
   TZK%MODELLEVEL = 1
! 
!*       3.2  update photolysis rates
ZUT       = DMOD( 24.0+DMOD(XTSIMUL/3600.,24.0),24.0 )
ZTUT = ZUT - ZTSIDER + ZLON*((180./ZPI)/15.0)
ZSOLANG = (ZTUT-12.0)*15.0*(ZPI/180.)          ! hour angle in radians
ZCOSZEN = SIN(ZLAT)*ZSINDEL +                 &! Cosine of the zenithal
              COS(ZLAT)*ZCOSDEL*COS(ZSOLANG)  !       solar angle
ZZENITH(:,:) = ACOS(ZCOSZEN)

CALL CH_UPDATE_JVALUES(6, ZZENITH, ZRT,              &
       ZALBUV, ZZS, ZZZ, XLAT, XLON,                &
       SIZE(ZZZ,1), SIZE(ZZZ,2), SIZE(ZZZ,3), SIZE(ZRT,4), &
       IDAY, IMONTH, IYEAR, XTSIMUL,&
       LCH_TUV_ONLINE, CCH_TUV_CLOUDS, &
       XCH_TUV_ALBNEW, XCH_TUV_DOBNEW, ZRHODREF, XJVALUES,&
       1,1,1,1,3,3,10)
!
  ITCOUNT = ITCOUNT + 1
  IF ( MOD(ITCOUNT, MAX(1, INT(XCH_TUV_TUPDATE/MAX(XDTACT,1.))) ) .EQ. 0 ) THEN
!
    CALL CH_SET_PHOTO_RATES(XTSIMUL,ZCONC,1,TZM,1,6,NVERB,1,KMASK,NEQ, XJVALUES)
  ENDIF
! 

!*       3.3  write results to disk
!
  write_to_disk : IF (XTSIMUL >= XTNEXTOUT) THEN

    ZCONC(1,:) = ZCONC(1,:)*1.E9 / (ZDEN2MOL * ZRHODREF(1,1,1)) ! convert molec.cm-3 to ppb
    IF (LORILAM) ZAERO(1,:) = ZAERO(1,:) / (ZDEN2MOL * ZRHODREF(1,1,1)) ! convert molec.cm-3 to ppp
    CALL CH_OUTPUT(ZCONC,ZAERO, ZMI, TZM, 1, 1)
    ZCONC(1,:) = ZCONC(1,:)*1.E-9 * (ZDEN2MOL * ZRHODREF(1,1,1)) ! convert  ppb to molec.cm-3
    IF (LORILAM) ZAERO(1,:) = ZAERO(1,:) * (ZDEN2MOL * ZRHODREF(1,1,1)) ! convert  ppp to molec.cm-3

    XTNEXTOUT = XTSIMUL + XDTOUT
  ENDIF write_to_disk
! 
!*       3.4  calculate diagnostics and write them to disk
!
  diagnostics : IF (XTSIMUL >= XTNEXTDIAG) THEN
    CALL CH_DIAGNOSTICS(ZCONC, 1, 1)
    XTNEXTDIAG = XTSIMUL + XDTDIAG
  ENDIF diagnostics
! 
!*       3.5  exit timeloop when end of simulation has been reached
!
  IF (XTSIMUL >= XTEND) EXIT time_loop
! 
!*       3.6  solve differential eqn. for next timestep
!

   CALL CH_SOLVER_n(XTSIMUL, XDTACT, ZCONC, ZNEWCONC, NEQ, 1, 1)


IF (LORILAM)  ZSO4RAT(:) = (ZNEWCONC(:,JP_CH_H2SO4) - ZCONC(:,JP_CH_H2SO4) ) / XDTACT

!
!        3.7 integrate aerosol variables
IF (LORILAM) THEN
  CALL CH_ORILAM(ZAERO,ZNEWCONC, ZM, ZSIG0, ZRG0, ZN0,  ZCTOTG, ZCTOTA,      &
                 ZCCTOT,  XDTACT , ZSEDA,ZMU, ZLAMBDA, ZRHOP0, ZOM, ZSO4RAT, &
                 ZRV, ZDENAIR, ZPABST, ZTEMP, ZRC, ZFRAC, ZMI, XTSIMUL,      &
                 CCH_SCHEME, ZSOLORG)

  ZNEWAERO(1,:) = ZAERO(1,:)


END IF

! 
!*       3.8  terminate the timestep by updating the variables
!
  ZCONC(1,:) = ZNEWCONC(1,:)

  XTSIMUL  = XTSIMUL + XDTACT
  XDTACT   = MIN(XDTACT,XTEND-XTSIMUL)
!
!
ENDDO time_loop

!
!*       4.   TERMINATE
!        --------------
!   
!*       4.1  write final result to disk (restart file)
!
! convert molec.cm-3 to ppb
ZCONC(1,:) = ZCONC(1,:) / (ZDEN2MOL * ZRHODREF(1,1,1)) ! convert molec.cm-3 to ppp
IF (LORILAM) ZAERO(1,:) = ZAERO(1,:) / (ZDEN2MOL * ZRHODREF(1,1,1)) ! convert molec.cm-3 to ppp
CALL CH_WRITE_CHEM(ZCONC(1,:), ZAERO(1,:), ZRHODREF(:,1,1), COUTFILE)
!
!*       4.2  finish
!
PRINT *, "CH_MODEL0D: simulation finished correctly, XTSIMUL = ", XTSIMUL
!
END PROGRAM CH_MODEL0D
