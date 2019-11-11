!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/Attic/radar_simulator.f90,v $ $Revision: 1.1.2.3.2.1.12.2.2.2 $ $Date: 2015/09/16 14:31:20 $
!-----------------------------------------------------------------
!      ###########################
       MODULE MODI_RADAR_SIMULATOR 
!      ###########################
!
INTERFACE
    SUBROUTINE RADAR_SIMULATOR(PUM,PVM,PWM,PRT,PCIT,PRHODREF,PTEMP,PPABSM,PREFL_CART,PLATLON,PCRT)
! variables en entree
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET :: PUM,PVM,PWM ! wind components
REAL,DIMENSION(:,:,:,:),INTENT(IN), TARGET  :: PRT  ! microphysical  mix. ratios at t
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET  :: PCIT ! pristine ice concentration at t
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET  :: PRHODREF ! density of the ref. state
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET  :: PTEMP    ! air temperature
REAL,DIMENSION(:,:,:),  INTENT(IN)  :: PPABSM   ! Absolute pressure
! variables en sortie
REAL,DIMENSION(:,:,:,:,:),INTENT(OUT) :: PREFL_CART ! radar variables (including reflectivity in dBZ) on observation cartesian grid
REAL,DIMENSION(:,:,:),  INTENT(OUT) :: PLATLON! latlon of cartesian grid points
REAL,DIMENSION(:,:,:),OPTIONAL,INTENT(IN), TARGET  :: PCRT  ! rain concentration at t

!
END SUBROUTINE RADAR_SIMULATOR
!
END INTERFACE
!
END MODULE MODI_RADAR_SIMULATOR
!
!     #########################################################################
      SUBROUTINE RADAR_SIMULATOR(PUM,PVM,PWM,PRT,PCIT,PRHODREF,PTEMP,PPABSM, &
           PREFL_CART,PLATLON,PCRT)
!     #########################################################################
!
!!****  *RADAR_SIMULATOR * - computes some pertinent radar parameters on PPIs
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the equivalent reflectivity
!!     of a mixed phase cloud on a PPI and then project it on a cartesian grid.
!!
!!**  METHOD
!!    ------
!!      First the geometry of radar data (radar, elevation, azimuth, range bin)
!!    is defined. Then ray paths are determined. Necessary model variables are 
!!    interpolated on each bin.
!!    A call to RADAR_SCATTERING yields reflectivities for each bin. In the end,
!!    reflectivities are interpolated on a cartesian grid.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!      Module MODD_REF
!!      Module MODD_RAIN_ICE_DESCR
!!      Module MODD_RAIN_ICE_PARAM
!!      Module MODD_PARAMETERS
!!      Module MODD_LUNIT
!
!!      Module MODE_IO_ll
!!      Module MODE_FM
!
!!      Module MODD_GR_FIELD_n
!!      Module MODD_GRID_n
!!      Module MODD_CONF_n
!!      Module MODD_GRID 
!!      Module MODE_GRIDPROJ
!!      Module MODD_RADAR 
!!      Module MODE_INTERPOL_BEAM
!!      Module MODE_FGAU  
!!      Module MODI_RADAR_SCATTERING
!!      Module MODI_SET_MSK 
!!      Module MODD_BUDGET
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RADAR_SIMULATOR )
!!
!!    AUTHOR
!!    ------
!!      O. Caumont & V. Ducrocq      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    26/03/2004
!!    O. Caumont  14/09/2009 modifications to allow for polar outputs
!!    O. Caumont 11/02/2010 thresholding and conversion from linear to log values after interpolation instead of before.
!!    O. Caumont 01/2011 gate-to-gate path computations revised (formulation+efficiency); comments in outputs revised
!!    C. Augros  22/02/2012 add of comments
!!    Attention : dans cette version, la temperature est forcee a 15° et l'elevation
!!    a la valeur donnee dans DIAG (pas de prise en compte de la courbure du faisceau)
!!
!!    ------------- MY_MODIF 8 ------
!!    C. Augros 02/2013
!!    All thresholding is done in radar_scattering.
!!
!!    ------------- MY_MODIF 9 ------
!!    C. Augros 20/02/2013
!!    Calculation of RHV, PDP and DHV
!!
!!    ------------- MY_MODIF 10 ------
!!    C. Augros 28/02/2013
!!    add of Avv (specific vertical attenuation) and T° in output files
!!
!!    ------------- MY_MODIF 11 ------
!!    C. Augros 7/03/2013
!!    Add of NDIFF=7 (TmatInt) for snow and graupel
!!    RHR-RHG, ZDA-ZDG, KDR-KDG (RhoHV, ZDR and Kdp for rain, snow and graupel)
!!    in output files
!!
!!    C. Augros 13/03/2013
!!    Thresholding of radar variables for each specie (specific SNR value for rain, snow, ice, graupel)
!!    
!!    C. Augros 22/03/2013
!!    Correction of interpolation part: add of "IF LATT" to set variables (AER-AEG, ATR-ATG)
!!    to xvalground or xundef
!!
!!    C. Augros 27/03/2013
!!    for polar output:
!!    NBAZIM set in nameliste (720) 
!!    ZAZIM_BASE(JAZ)=(0.5+JAZ-1)*ZZSTEP
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!       
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST              , ONLY: XPI,XRD,XRV,XRADIUS
USE MODD_REF
USE MODD_PARAMETERS
USE MODD_LUNIT
!
USE MODE_IO_ll
USE MODE_MSG
!
USE MODD_GR_FIELD_n
USE MODD_GRID_n
USE MODD_CONF_n
USE MODD_GRID             , ONLY: XLATORI,XLONORI,XRPK,XLAT0,XLON0,XBETA
USE MODE_GRIDPROJ
USE MODD_RADAR            , ONLY: XLAT_RAD,XLON_RAD,XALT_RAD,XDT_RAD,XELEV,&
     XX_INI,XY_INI,XZ_INI,XSTEP_RAD,NBRAD,NBELEV,NBAZIM,NBSTEPMAX,&
     NCURV_INTERPOL,LATT,LCART_RAD,NPTS_H,NPTS_V,LQUAD,XGRID,XVALGROUND,&
     NMAX,LREFR,LDNDZ,XREFLMIN
!
USE MODE_INTERPOL_BEAM
USE MODE_FGAU             , ONLY: GAULEG,GAUHER
USE MODI_RADAR_SCATTERING
! convective/stratiform
USE MODI_SET_MSK 
! /convective/stratiform
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET :: PUM,PVM,PWM ! wind components
REAL,DIMENSION(:,:,:,:),INTENT(IN), TARGET :: PRT  ! microphysical  mix. ratios at t
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET :: PCIT ! pristine ice concentration at t
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET :: PRHODREF ! density of the ref. state
REAL,DIMENSION(:,:,:),  INTENT(IN), TARGET :: PTEMP    ! air temperature
REAL,DIMENSION(:,:,:),  INTENT(IN)  :: PPABSM     ! Absolute pressure
!
REAL,DIMENSION(:,:,:,:,:),INTENT(OUT) :: PREFL_CART ! radar reflectivity in dBZ and other parameters on observation cartesian or polar grid
REAL,DIMENSION(:,:,:),    INTENT(OUT) :: PLATLON! latlon of cartesian grid points
REAL,DIMENSION(:,:,:),OPTIONAL,INTENT(IN), TARGET  :: PCRT  ! rain concentration at t

!
!*       0.2   Declarations of local variables :
!
!
TYPE(PAMOD),DIMENSION(:),ALLOCATABLE :: TVARMOD ! array of pointers to grid-point model fields
TYPE(PARAD),DIMENSION(:),ALLOCATABLE :: TVARRAD ! array of pointers to model fields interpolated onto radar PPIs
REAL :: ZRDSRV      ! XRD/XRV
REAL :: ZRDSDG      ! PI/180
INTEGER :: IIB,IIE          ! Loop limits for coordinate X
INTEGER :: IJB,IJE          ! Loop limits for coordinate Y
INTEGER :: IKB,IKE          ! Loop limits for coordinate Z
INTEGER :: IIU,IJU,IKU      ! Loop variables of model 
INTEGER :: IIELV ! maximum number of elevations 
INTEGER :: ILUOUT0 ! Logical unit number for output-listing
INTEGER :: IRESP   ! Return code of FM-routines
INTEGER :: JI,JL,JEL,JAZ,JH,JV ! Loop variables of control
INTEGER :: IEL,IND,INDV
INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE :: IREFL_CART_NB,IREFL_CART_NBR,IREFL_CART_NBI,IREFL_CART_NBS,IREFL_CART_NBG!
INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE :: IVDOP_CART_NB!
REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE  :: ZZE ! gate equivalent reflectivity factor, ZDR, KDP, 
REAL,DIMENSION(:,:,:,:),ALLOCATABLE    :: ZELEV        ! elevation in rad.
REAL,DIMENSION(:),  ALLOCATABLE        :: ZAZIM_BASE   ! azimuth in rad. of the beam centre
REAL,DIMENSION(:,:,:),  ALLOCATABLE    :: ZAZIM        ! azimuth in rad. of discretized beam
!
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZX_RAY   ! x positions of the points along the ray-tracing
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZY_RAY   ! y positions of the points along the ray-tracing
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZZ_RAY   ! z positions of the points along the ray_tracing
REAL, DIMENSION(:,:),        ALLOCATABLE :: ZLAT     ! latitude of the points 
REAL, DIMENSION(:,:),        ALLOCATABLE :: ZLON     ! longitude of the points 
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZDX_NAT ! x increment on the natural referential
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZDY_NAT ! y increment on the natural referential
REAL                                     :: ZRPK,ZLAT0,ZBETA,ZLON0 ! projection characteristics
REAL                                     :: ZCLAT0,ZSLAT0          ! cos and sin of lat0
REAL                                     :: ZMAP     ! Map factor
REAL                                     :: ZGAMMA,ZCOSG,ZSING  ! angle of projection and its cos and sin values
!
REAL, DIMENSION(:),          ALLOCATABLE :: ZXHATM ! X values of the mass points
REAL, DIMENSION(:),          ALLOCATABLE :: ZYHATM ! Y values of the mass points
REAL, DIMENSION(:,:,:),      ALLOCATABLE :: ZZM    ! Z values of the mass points
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZT_RAY  ! temperature interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZR_RAY  ! rain            mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZI_RAY  ! pristine ice    mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZS_RAY  ! snow/aggregates mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZG_RAY  ! graupel         mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZH_RAY  ! hail            mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZCIT_RAY  ! pristine ice concentration interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZRHODREF_RAY
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZVDOP_RAY  
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZUM_RAY,ZVM_RAY,ZWM_RAY 
!
REAL :: ZKE=4./3. ! COEFF OF EFFECTIVE RADIUS IN AIR REFRACTIVE INDEX COMPUTATION
REAL :: ZN0,ZN1 ! REFRACTIVE INDEX OF AIR
REAL :: ZDNDZ1 ! vertical gradient of refractive index of air
!

REAL,DIMENSION(:),ALLOCATABLE :: ZX_H,ZW_H ! Gauss-Hermite horizontal points & weights
REAL,DIMENSION(:),ALLOCATABLE :: ZX_V,ZW_V ! Gauss-Hermite vertical points & weights
 
INTEGER :: IXGRID,IYGRID
INTEGER :: INVAR
REAL :: r,h,alph

! convective/stratiform
LOGICAL,DIMENSION(:,:,:),ALLOCATABLE             :: GBU_MSK
REAL, DIMENSION(:,:,:),ALLOCATABLE, TARGET       :: ZBU_MASK
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZBU_MASK_RAY
! refractivity
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZN_RAY,ZDNDZ_RAY ! refractivity and its vertical gradient in radar coordinates
REAL, DIMENSION(:,:,:),ALLOCATABLE       :: ZN,ZDNDZ ! index of ZN_RAY,ZDNDZ_RAY in ZZE

INTEGER :: IZER,IZEI,IZES,IZEG
INTEGER :: IVDOP
INTEGER :: IAER,IAEI,IAES,IAEG
INTEGER :: IAVR,IAVI,IAVS,IAVG
INTEGER :: IATR,IATI,IATS,IATG
INTEGER :: IRHV,IPDP,IDHV
INTEGER :: IRHR, IRHS, IRHG, IZDA, IZDS, IZDG, IKDR, IKDS, IKDG
INTEGER :: IHAS,IRFR,IDNZ   
REAL :: ZZSTEP
INTEGER :: IZEH, IRHH,IKDH,IZDH ! hail
INTEGER :: IAEH,IAVH,IATH
INTEGER :: IMR,IMI,IMS,IMG,IMH,ICIT,ITEM,ICRT
LOGICAL :: GHAIL
INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE ::IREFL_CART_NBH
!
!Modif pour LIMA
!
LOGICAL :: GLIMA
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE, TARGET :: ZCRT_RAY  ! rain concentration interpolated along the rays!
REAL, DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZWORK
!-------------------------------------------------------------------------------
!
!
!*       1.     INITIALIZATION 
!   	        --------------
!
!
!*       1.1 IO and dimensions initialization
!   
ILUOUT0 = TLUOUT0%NLU
!
IIU=SIZE(PTEMP,1)
IJU=SIZE(PTEMP,2)
IKU=SIZE(PTEMP,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = JPVEXT + 1
IKE = IKU - JPVEXT
! convective/stratiform
ALLOCATE(GBU_MSK(IIU,IJU,4))
CALL SET_MSK(PRT,PRHODREF,GBU_MSK) ! on récupère GBU_MSK
ALLOCATE(ZBU_MASK(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)))
ZBU_MASK(:,:,1)=0.
WHERE(GBU_MSK(:,:,2)) ! stratiform
   ZBU_MASK(:,:,1)=1.
END WHERE
WHERE(GBU_MSK(:,:,1)) ! convective
   ZBU_MASK(:,:,1)=2.
END WHERE
DEALLOCATE(GBU_MSK)
ZBU_MASK(:,:,:)=SPREAD(ZBU_MASK(:,:,1),DIM=3,NCOPIES=SIZE(PTEMP,3))
IIELV=MAXVAL(NBELEV(1:NBRAD))
! LIMA
IF (PRESENT(PCRT)) THEN
  GLIMA=.TRUE.
ELSE
  GLIMA=.FALSE.
ENDIF
IF (SIZE(PRT,4)== 7 ) THEN
  GHAIL=.TRUE.
ELSE
  GHAIL=.FALSE.
ENDIF

!
!*       1.2 Some constants and parameters
!   
ZRDSRV=XRD/XRV       ! XRD/XRV
ZRDSDG=XPI/180.       ! PI/180
!
!*       1.3 beam characteristics initialization
!  
! azimuths 0=N 90=E
WRITE(ILUOUT0,*) "NBAZIM",NBAZIM
ALLOCATE(ZAZIM_BASE(NBAZIM),ZAZIM(NBRAD,NBAZIM,NPTS_H)) 
!calculation of the azimut of the center of the beam (ZAZIM_BASE)
!so that each pixel of the square grid circling the PPI (which number of pixel in the radius is NMAX)
!NMAX=INT(NBSTEPMAX*XSTEP_RAD/XGRID) : number of range 
!contains one azimut only

IF (LCART_RAD) THEN  
  DO JAZ=1,NBAZIM ! NBAZIM defined in mode_interpol_beam
    IF(JAZ<=NMAX) THEN
      ZAZIM_BASE(JAZ)=ATAN((JAZ-.5)/NMAX)
    ELSE IF(JAZ<=3*NMAX) THEN
      ZAZIM_BASE(JAZ)=XPI/2+ATAN((-2*NMAX+JAZ-.5)/NMAX)
    ELSE IF(JAZ<=5*NMAX) THEN
      ZAZIM_BASE(JAZ)=XPI-ATAN((4*NMAX-JAZ+.5)/NMAX)
    ELSE IF(JAZ<=7*NMAX) THEN
      ZAZIM_BASE(JAZ)=3*XPI/2-ATAN((6*NMAX-JAZ+.5)/NMAX)
    ELSE 
      ZAZIM_BASE(JAZ)=2*XPI-ATAN((8*NMAX-JAZ+.5)/NMAX)
    END IF
  END DO
ELSE
  ZZSTEP=2*XPI/NBAZIM
  DO JAZ=1,NBAZIM
      ZAZIM_BASE(JAZ)=(0.5+JAZ-1)*ZZSTEP
  END DO
END IF

WRITE(ILUOUT0,*) "ZAZIM_BASE(1)",ZAZIM_BASE(1)
WRITE(ILUOUT0,*) "ZAZIM_BASE(NBAZIM/2.)",ZAZIM_BASE(NINT(NBAZIM/2.))
WRITE(ILUOUT0,*) "ZAZIM_BASE(NBAZIM)",ZAZIM_BASE(NBAZIM)

!copy in the 3D matrix (ZAZIM) containing the horizontal discretization of the beam for all azimut of all radars 
ZAZIM(:,1:NBAZIM,:)=SPREAD(SPREAD(ZAZIM_BASE(1:NBAZIM),DIM=1,NCOPIES=NBRAD),DIM=3,NCOPIES=NPTS_H)
!  
! elevations 
ALLOCATE(ZELEV(NBRAD,IIELV,NBSTEPMAX+1,NPTS_V))
! 4D matrix containing the vertical discretization of the beam for all elevations (with a value for each range step)
ZELEV(:,:,:,:)=SPREAD(SPREAD(XELEV(1:NBRAD,1:IIELV),DIM=3,NCOPIES=NBSTEPMAX+1),&
     DIM=4,NCOPIES=NPTS_V)
! 
! Discretization of the gate 
! Calculation of the position ZX_H and ZX_V and weights ZW_H and ZW_V for each discretization of the beam (horizontally and vertically)
ALLOCATE(ZX_H((NPTS_H+1)/2),ZW_H((NPTS_H+1)/2))
ALLOCATE(ZX_V((NPTS_V+1)/2),ZW_V((NPTS_V+1)/2))
IF(LQUAD) THEN
   CALL GAULEG(NPTS_H,ZX_H,ZW_H)
   CALL GAULEG(NPTS_V,ZX_V,ZW_V)
   XDT_RAD(:)=XDT_RAD(:)/2.
ELSE
   CALL GAUHER(NPTS_H,ZX_H,ZW_H)
   CALL GAUHER(NPTS_V,ZX_V,ZW_V)
   XDT_RAD(:)=XDT_RAD(:)/SQRT(8.*LOG(2.)) ! variable change
END IF
!
DO JI=1,NBRAD
  IEL=NBELEV(JI)
  DO JV=1,NPTS_V
   ! to change if NPTS_V even
    ZELEV(1:NBRAD,1:IEL,:,JV)=ZELEV(1:NBRAD,1:IEL,:,JV)&
         +SIGN(ZX_V(ABS((2*JV-NPTS_V-1)/2)+1),2.*JV-NPTS_V-1.)* &
         SPREAD(SPREAD(XDT_RAD(1:NBRAD),DIM=2,NCOPIES=IIELV),DIM=3,NCOPIES=NBSTEPMAX+1)
  END DO
END DO

DO JH=1,NPTS_H
   ! to change if NPTS_H even
  ZAZIM(1:NBRAD,:,JH)=ZAZIM(1:NBRAD,:,JH)+SIGN(ZX_H(ABS((2*JH-NPTS_H-1)/2)+1),2.*JH-NPTS_H-1.)* &
       SPREAD(XDT_RAD(1:NBRAD),DIM=2,NCOPIES=NBAZIM)*ZRDSDG
END DO
ZELEV(:,:,:,:)=ZELEV(:,:,:,:)*ZRDSDG ! in radian
! initialisation of refractivity indices
! 1 : ZHH
! 2 : ZDR
! 3 : KDP
! 4 : CSR
IZER=5 ! ZER
IZEI=IZER+1 ! ZEI
IZES=IZEI+1 ! ZES
IZEG=IZES+1 ! ZEG
IF (GHAIL) THEN
  IZEH=IZEG+1 !ZEH
  IVDOP=IZEH+1 !VRU
ELSE
  IVDOP=IZEG+1 !VRU
END IF
IF (LATT) THEN
  IF (GHAIL) THEN
    IAER=IVDOP+1
    IAEI=IAER+1
    IAES=IAEI+1
    IAEG=IAES+1
    IAEH=IAEG+1
    IAVR=IAEH+1
    IAVI=IAVR+1
    IAVS=IAVI+1
    IAVG=IAVS+1
    IAVH=IAVG+1
    IATR=IAVH+1
    IATI=IATR+1
    IATS=IATI+1
    IATG=IATS+1
    IATH=IATG+1
    IRHV=IATH+1
  ELSE
    IAER=IVDOP+1
    IAEI=IAER+1
    IAES=IAEI+1
    IAEG=IAES+1
    IAVR=IAEG+1
    IAVI=IAVR+1
    IAVS=IAVI+1
    IAVG=IAVS+1
    IATR=IAVG+1
    IATI=IATR+1
    IATS=IATI+1
    IATG=IATS+1
    IRHV=IATG+1
  ENDIF
ELSE
    IRHV=IVDOP+1        
ENDIF
IPDP=IRHV+1
IDHV=IPDP+1
IRHR=IDHV+1
IRHS=IRHR+1
IRHG=IRHS+1
IF (GHAIL) THEN
  IRHH=IRHG+1
  IZDA=IRHH+1
ELSE
  IZDA=IRHG+1
ENDIF
IZDS=IZDA+1
IZDG=IZDS+1
IF (GHAIL) THEN
  IZDH=IZDG+1
  IKDR=IZDH+1
ELSE
  IKDR=IZDG+1
ENDIF
IKDS=IKDR+1
IKDG=IKDS+1
IF (GHAIL) THEN
  IKDH=IKDG+1
  IHAS=IKDH+1
ELSE
  IHAS=IKDG+1
ENDIF
IMR=IHAS+1
IMI=IMR+1
IMS=IMI+1
IMG=IMS+1
IF (GHAIL) THEN
  IMH=IMG+1
  ICIT=IMH+1
ELSE
  ICIT=IMG+1
ENDIF
ITEM=ICIT+1
IF (GLIMA) THEN
  ICRT=ITEM+1
  IF(LREFR) THEN
    IRFR=ICRT+1
    IF(LDNDZ) IDNZ=IRFR+1
  ELSE
    IF(LDNDZ) IDNZ=ICRT+1
  ENDIF
ELSE
  IF(LREFR) THEN
    IRFR=ITEM+1
    IF(LDNDZ) IDNZ=IRFR+1
  ELSE
    IF(LDNDZ) IDNZ=ITEM+1
  ENDIF
ENDIF
!
!----------------------------------------------------------------------------------------
!*       2.    RAY TRACING DEFINITION
!              ----------------------
!
!*       2.1 some initializations for MESO-NH conformal projection
!
IF (XRPK<0.) THEN     ! projection from north pole 
  ZRPK=-XRPK
  ZLAT0=-XLAT0
  ZBETA=-XBETA
  ZLON0=XLON0+180.
  WRITE(ILUOUT0,*) 'projection from north pole'
  WRITE(ILUOUT0,*) 'ZRPK',ZRPK
  WRITE(ILUOUT0,*) 'ZLAT0',ZLAT0
  WRITE(ILUOUT0,*) 'ZBETA',ZBETA
  WRITE(ILUOUT0,*) 'ZLON0',ZLON0
ELSE                  ! projection from south pole
  ZRPK=XRPK
  ZLAT0=XLAT0
  ZBETA=XBETA
  ZLON0=XLON0
  WRITE(ILUOUT0,*) 'projection from south pole'
  WRITE(ILUOUT0,*) 'ZRPK:',ZRPK
  WRITE(ILUOUT0,*) 'ZLAT0:',ZLAT0
  WRITE(ILUOUT0,*) 'ZBETA:',ZBETA
  WRITE(ILUOUT0,*) 'ZLON0:',ZLON0
ENDIF
ZCLAT0 = COS(ZRDSDG*ZLAT0)
ZSLAT0 = SIN(ZRDSDG*ZLAT0)
!
!  Positions of the mass points in the MESO-NH conformal projection 
!
ALLOCATE(ZXHATM(IIU))
ALLOCATE(ZYHATM(IJU)) 
ALLOCATE(ZZM(IIU,IJU,IKU))
!
ZXHATM(1:IIU-1) = .5*(XXHAT(1:IIU-1)+XXHAT(2:IIU))
ZXHATM(IIU)     = 2.*XXHAT(IIU)-ZXHATM(IIU-1)
!
ZYHATM(1:IJU-1) = .5*(XYHAT(1:IJU-1)+XYHAT(2:IJU))
ZYHATM(IJU)     = 2.*XYHAT(IJU)-ZYHATM(IJU-1)
!
ZZM(:,:,1:IKU-1)= .5*(XZZ(:,:,1:IKU-1)+XZZ(:,:,2:IKU))
ZZM(:,:,IKU)= 2. * XZZ(:,:,IKU) - ZZM(:,:,IKU-1)
!
!*       2.2 initialization of the ray beginning
!
!
! position of the radar on the MESO-NH conformal projection grid
!
DO JI=1,NBRAD
   CALL SM_XYHAT(XLATORI,XLONORI,   & 
       XLAT_RAD(JI),XLON_RAD(JI),XX_INI(JI),XY_INI(JI))  
END DO
XZ_INI(:)=XALT_RAD(:)       ! z positions of the ground source signal 
!
ALLOCATE(ZX_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V),& !6D matrix : X position of the pixel in the model grid for each range
         ZY_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V),& !of each discretisation of each elevation and each azimut of each radar 
         ZZ_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V),&
         ZLAT(NPTS_H,NPTS_V),&
         ZLON(NPTS_H,NPTS_V),&
         ZDX_NAT(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V),&
         ZDY_NAT(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
!
ZX_RAY(:,:,:,:,:,:)=0.
ZY_RAY(:,:,:,:,:,:)=0.
ZZ_RAY(:,:,:,:,:,:)=0.
ZDX_NAT(:,:,:,:,:,:)=0.
ZDY_NAT(:,:,:,:,:,:)=0.
!
ZX_RAY(1:NBRAD,:,:,1,:,:)=SPREAD(SPREAD(SPREAD(SPREAD(XX_INI(1:NBRAD),DIM=2,NCOPIES=IIELV),& !initialization with the X position of the radar
     DIM=3,NCOPIES=NBAZIM),DIM=4,NCOPIES=NPTS_H),DIM=5,NCOPIES=NPTS_V)
ZY_RAY(1:NBRAD,:,:,1,:,:)=SPREAD(SPREAD(SPREAD(SPREAD(XY_INI(1:NBRAD),DIM=2,NCOPIES=IIELV), & !initialization with the Y position of the radar
     DIM=3,NCOPIES=NBAZIM),DIM=4,NCOPIES=NPTS_H),DIM=5,NCOPIES=NPTS_V)
ZZ_RAY(1:NBRAD,:,:,1,:,:)=SPREAD(SPREAD(SPREAD(SPREAD(XZ_INI(1:NBRAD),DIM=2,NCOPIES=IIELV), & !initialization with the Z position of the radar
     DIM=3,NCOPIES=NBAZIM),DIM=4,NCOPIES=NPTS_H),DIM=5,NCOPIES=NPTS_V)
! refractivity
IF(LREFR) ALLOCATE(ZN_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
IF(LDNDZ) ALLOCATE(ZDNDZ_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
IF(NCURV_INTERPOL == 1) THEN
  ALLOCATE(ZN(IIU,IJU,IKU),ZDNDZ(IIU,IJU,IKU))
  ZN(:,:,:)=1.+.776E-6*PPABSM(:,:,:)*( 1. + 4810.*PRT(:,:,:,1)/(ZRDSRV+PRT(:,:,:,1))/PTEMP(:,:,:) )/PTEMP(:,:,:)
  ! vertical gradient of n approximated using a centred finite difference scheme
  ZDNDZ(:,:,IKB+1:IKE-1)=((ZZM(:,:,IKB+1:IKE-1)-ZZM(:,:,IKB:IKE-2))*2*(ZN(:,:,IKB+2:IKE)-ZN(:,:,IKB+1:IKE-1)) &
    +(ZZM(:,:,IKB+2:IKE)-ZZM(:,:,IKB+1:IKE-1))**2*(ZN(:,:,IKB+1:IKE-1)-ZN(:,:,IKB:IKE-2))) &
    /(ZZM(:,:,IKB+1:IKE-1)-ZZM(:,:,IKB:IKE-2))/(ZZM(:,:,IKB+2:IKE)-ZZM(:,:,IKB+1:IKE-1))/(ZZM(:,:,IKB+2:IKE)-ZZM(:,:,IKB:IKE-2))
! + valeurs limites
  ZDNDZ(:,:,IKB)=(ZN(:,:,IKB+1)-ZN(:,:,IKB))/(ZZM(:,:,IKB+1)-ZZM(:,:,IKB))
  ZDNDZ(:,:,IKE)=(ZN(:,:,IKE)-ZN(:,:,IKE-1))/(ZZM(:,:,IKE)-ZZM(:,:,IKE-1))
ENDIF

!
!*       2.3  positions of the rays in the MESO-NH conformal projection (calculation of ZX_RAY, ZY_RAY and ZZ_RAY)
!
DO JI=1,NBRAD
  IEL=NBELEV(JI)
  WRITE(ILUOUT0,*) 'RADAR #',JI,'Number of ELEVATIONS: ',NBELEV(JI)
  WRITE(ILUOUT0,*) '  Elevations used:'
  DO JEL=1,IEL
    WRITE(ILUOUT0,*) "    ",ZELEV(JI,JEL,1,:)/ZRDSDG
    DO JAZ=1,NBAZIM
      label: DO JL=1,NBSTEPMAX
        ! SM_LATLON takes bidimensional arrays as arguments
        CALL SM_LATLON(XLATORI,XLONORI,   &   
               ZX_RAY(JI,JEL,JAZ,JL,:,:), ZY_RAY(JI,JEL,JAZ,JL,:,:),ZLAT(:,:),ZLON(:,:))
        DO JH=1,NPTS_H
          DO JV=1,NPTS_V
            !   
            !*        Compute positions of the gates
            !      
            ! Compute local Map factor and other projection factors
            IF(XRPK<0.)  ZLAT(JH,JV)=-ZLAT(JH,JV)     ! projection from north pole 
               
            IF(ABS(ZRPK-1.)>1.E-10 .AND. ABS(COS(ZRDSDG*ZLAT(JH,JV)))<1.E-10) THEN
              WRITE(ILUOUT0,*) 'Error in projection : '
              WRITE(ILUOUT0,*) 'pole in the domain, but not with stereopolar projection'
              !callabortstop
              CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SIMULATOR',&
                             'Error in projection: pole in the domain, but not with stereopolar projection')
            ENDIF
            !
            IF(ABS(ZCLAT0)<1.E-10 .AND. ABS(ZRPK-1.)<1.E-10) THEN
              ZMAP = (1.+ZSLAT0)/(1.+SIN(ZRDSDG*ZLAT(JH,JV)))
            ELSE IF(ABS(COS(ZRDSDG*ZLAT(JH,JV)))>1.E-10) THEN
              ZMAP = ((ZCLAT0/COS(ZRDSDG*ZLAT(JH,JV)))**(1.-ZRPK))      &
                   * ((1.+ZSLAT0)/(1.+SIN(ZRDSDG*ZLAT(JH,JV))))**ZRPK
            ELSE
              ZMAP = (1.+ZSLAT0)/(1.+SIN(ZRDSDG*ZLAT(JH,JV)))
            END IF
            ZGAMMA=(ZRPK*(ZLON(JH,JV)-ZLON0)-ZBETA)*ZRDSDG
            ZCOSG=COS(ZGAMMA)
            ZSING=SIN(ZGAMMA)
                  
            ! compute positions of radar gates (2 methods) : 
            !  First method : gate-to-gate computation using the model's index of refraction
            IF(NCURV_INTERPOL == 1) THEN 
              ! first compute vertical position (height)          
              ! compute the index of refraction at the radar gate boundaries
              CALL INTERPOL_BEAM(ZN(:,:,:),ZN1,ZX_RAY(JI,JEL,JAZ,JL,JH,JV),&
              ZY_RAY(JI,JEL,JAZ,JL,JH,JV),ZZ_RAY(JI,JEL,JAZ,JL,JH,JV),ZXHATM(:),ZYHATM(:),ZZM(:,:,:))
              IF(LREFR) ZN_RAY(JI,JEL,JAZ,JL,JH,JV)=(ZN1-1.)*1.E6 !LREFR: if true writes out refractivity (N ≡ (n − 1) × 106)
              IF(LDNDZ) THEN !LDNDZ: if true writes out vertical gradient of refractivity
                IF(JL==1) THEN
                  ZDNDZ_RAY(JI,JEL,JAZ,JL,JH,JV)=0. ! this is not true, this is set to XVALGROUND afterwards
                ELSE
                  ZDNDZ_RAY(JI,JEL,JAZ,JL,JH,JV)=(ZN1-ZN0)*1.E6/(ZZ_RAY(JI,JEL,JAZ,JL,1,1)-ZZ_RAY(JI,JEL,JAZ,JL-1,1,1))
                END IF
              END IF
              IF(ZN1==-XUNDEF) THEN ! we are underground
                ZZ_RAY(JI,JEL,JAZ,JL:NBSTEPMAX+1,:,:)=-XUNDEF ! rest of the ray is flagged undefined
                EXIT label
              ELSE
                IF(JL > 1) THEN 
                  ! next line to comment (std refraction)
                  !ZN1=ZN0-(ZZ_RAY(JI,JEL,JAZ,JL,1,1)-ZZ_RAY(JI,JEL,JAZ,JL-1,1,1))/(4.*XRADIUS)
                  IF(ZN0/ZN1*(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL-1,JH,JV))/(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL,JH,JV)) &
                       *COS(ZELEV(JI,JEL,JL-1,JV)) >= 1.) THEN ! it means the slope of the ray path is 0 relative to the Earth
                    ZELEV(JI,JEL,JL,JV)=-ACOS(2.-ZN0/ZN1*(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL-1,JH,JV)) &
                                /(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL,JH,JV))*COS(ZELEV(JI,JEL,JL-1,JV)))
                  ELSE ! usual formula
                    ZELEV(JI,JEL,JL,JV)=ZELEV(JI,JEL,JL-1,JV)/ABS(ZELEV(JI,JEL,JL-1,JV))* &
                                        ACOS(ZN0/ZN1*(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL-1,JH,JV))/            &
                                       (XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL,JH,JV))*COS(ZELEV(JI,JEL,JL-1,JV)))
                  END IF
                  ZDNDZ1=(ZN1-ZN0)/(ZZ_RAY(JI,JEL,JAZ,JL,1,1)-ZZ_RAY(JI,JEL,JAZ,JL-1,1,1))
                ELSE ! for first gate DNDZ1 is the local value at radar
                  CALL INTERPOL_BEAM(ZDNDZ(:,:,:),ZDNDZ1,ZX_RAY(JI,JEL,JAZ,JL,JH,JV),&
                            ZY_RAY(JI,JEL,JAZ,JL,JH,JV),ZZ_RAY(JI,JEL,JAZ,JL,JH,JV),ZXHATM(:),ZYHATM(:),ZZM(:,:,:))
                END IF
                IF(ZDNDZ1>-ZN1/XRADIUS/COS(ZELEV(JI,JEL,JL,JV))) THEN
                  ZKE=1./(1.+XRADIUS/ZN1*ZDNDZ1*COS(ZELEV(JI,JEL,JL,JV)))
                ELSE
                  ZKE=1./(1.-XRADIUS/ZN1*ZDNDZ1*COS(ZELEV(JI,JEL,JL,JV)))
                END IF
                ! elements finis
                !ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV)=SQRT(XSTEP_RAD**2+(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL,JH,JV))**2 &
                !+2.*XSTEP_RAD*(XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL,JH,JV))*SIN(ZELEV(JI,JEL,JL,JV)))-XRADIUS
                ! Doviak & Zrnic
                ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV)=ZZ_RAY(JI,JEL,JAZ,JL,JH,JV)+SQRT(XSTEP_RAD**2+(ZKE*XRADIUS)**2 &
                           +2.*XSTEP_RAD*ZKE*XRADIUS*SIN(ZELEV(JI,JEL,JL,JV)))-ZKE*XRADIUS
                ZN0=ZN1
                ! then compute horizontal position
                ZDX_NAT(JI,JEL,JAZ,JL,JH,JV)=XRADIUS*ASIN(XSTEP_RAD/ &
                        (XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV))*COS(ZELEV(JI,JEL,JL,JV)))*SIN(ZAZIM(JI,JAZ,JH))
                ZDY_NAT(JI,JEL,JAZ,JL,JH,JV)=XRADIUS*ASIN(XSTEP_RAD/ &
                        (XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV))*COS(ZELEV(JI,JEL,JL,JV)))*COS(ZAZIM(JI,JAZ,JH))
                ZX_RAY(JI,JEL,JAZ,JL+1,JH,JV)=ZX_RAY(JI,JEL,JAZ,JL,JH,JV) & !!!
                        +  (ZMAP* XRADIUS *((ZDX_NAT(JI,JEL,JAZ,JL,JH,JV) * ZCOSG) - &
                        (ZDY_NAT(JI,JEL,JAZ,JL,JH,JV)* ZSING)  ) &
                        /(ZZ_RAY(JI,JEL,JAZ,JL,JH,JV) + XRADIUS))
                ZY_RAY(JI,JEL,JAZ,JL+1,JH,JV)=ZY_RAY(JI,JEL,JAZ,JL,JH,JV) + &
                        (ZMAP* XRADIUS *((ZDX_NAT(JI,JEL,JAZ,JL,JH,JV) * ZSING) +  &
                        (ZDY_NAT(JI,JEL,JAZ,JL,JH,JV)* ZCOSG)  ) &
                        /(ZZ_RAY(JI,JEL,JAZ,JL,JH,JV) + XRADIUS))
              END IF
            ELSE 
              ! 2nd method : effective Earth radius model Doviak & Zrnic 1993 (2.28b) p. 21
              ! vertical position
              ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV)=SQRT((JL*XSTEP_RAD)**2+(ZKE*XRADIUS)**2+ &
                   2.*JL*XSTEP_RAD*ZKE*XRADIUS*SIN(ZELEV(JI,JEL,1,JV)))-ZKE*XRADIUS+ &
                   ZZ_RAY(JI,JEL,JAZ,1,JH,JV)
              ! This formula is given by Doviak & Zrnic 1993 (9.9 p. 307) 
              ZELEV(JI,JEL,JL+1,JV)=ZELEV(JI,JEL,1,JV)+ATAN(JL*XSTEP_RAD*COS(ZELEV(JI,JEL,1,JV))&
                  /(ZKE*XRADIUS+JL*XSTEP_RAD*SIN(ZELEV(JI,JEL,1,JV))))
              ! horizontal position (Doviak & Zrnic)
              ZDX_NAT(JI,JEL,JAZ,JL,JH,JV)=ZKE*XRADIUS*ASIN(JL*XSTEP_RAD*COS(ZELEV(JI,JEL,JL,JV)) &
                  /(ZKE*XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV)))*SIN(ZAZIM(JI,JAZ,JH))
              ZDY_NAT(JI,JEL,JAZ,JL,JH,JV)=ZKE*XRADIUS*ASIN(JL*XSTEP_RAD*COS(ZELEV(JI,JEL,JL,JV)) &
                  /(ZKE*XRADIUS+ZZ_RAY(JI,JEL,JAZ,JL+1,JH,JV)))*COS(ZAZIM(JI,JAZ,JH))
              ZX_RAY(JI,JEL,JAZ,JL+1,JH,JV)=ZX_RAY(JI,JEL,JAZ,1,JH,JV) & !!!
                  +  (ZMAP* XRADIUS *((ZDX_NAT(JI,JEL,JAZ,JL,JH,JV) * ZCOSG) - &
                  (ZDY_NAT(JI,JEL,JAZ,JL,JH,JV)* ZSING)  ) &
                  /(ZZ_RAY(JI,JEL,JAZ,JL,JH,JV) + XRADIUS))
              ZY_RAY(JI,JEL,JAZ,JL+1,JH,JV)=ZY_RAY(JI,JEL,JAZ,1,JH,JV) + &
                  (ZMAP* XRADIUS *((ZDX_NAT(JI,JEL,JAZ,JL,JH,JV) * ZSING) +  &
                  (ZDY_NAT(JI,JEL,JAZ,JL,JH,JV)* ZCOSG)  ) &
                  /(ZZ_RAY(JI,JEL,JAZ,JL,JH,JV) + XRADIUS))
            END IF
          END DO
        END DO
      END DO label
    END DO
  END DO
END DO
DEALLOCATE(ZLAT,ZLON)
DEALLOCATE(ZDX_NAT,ZDY_NAT)
IF(NCURV_INTERPOL == 1) DEALLOCATE(ZN,ZDNDZ)
! end of geometrical part ; I determined z[xyz]_ray
WRITE(ILUOUT0,*) 'BEAM DEFINITION DONE'
!
!-------------------------------------------------------------------------------
!*       3.    INTERPOLATION OF THE MODEL VARIABLES ON THE RAYS 
!              ------------------------------------------------
!  
!
!*       3.1  allocation of the arrays and initialization of the arrays of pointers 
!              (to avoid multiple calls to interpol_beam)
!
! 1: temperature; 2: rhodref, 3: rain mixing ratio; 4: r_i; 5: CIT; 6: r_s; 7: r_g; 8: convective/stratiform; 9: u; 10: v; 11: w; 12: rain concentration (LIMA only)
  ALLOCATE(TVARMOD(NRR+5)) ! pointer toward the matrix of model variables (Tempe, rhodref,mixing rations...) in the model projection (X,Y,Z)
  ALLOCATE(TVARRAD(NRR+5)) ! pointer toward the matrix of model variables (Tempe, rhodref,mixing rations...) interpolated in the
                         ! radar projection (iradar,ielev,iazim,irangestep,idiscretH,idiscretV)
       
TVARMOD(1)%P=>PTEMP(:,:,:)
TVARMOD(2)%P=>PRHODREF(:,:,:)

ALLOCATE(ZT_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V),& ! temperature and reference density interpolated along the ray 
     ZRHODREF_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
TVARRAD(1)%P=>ZT_RAY(:,:,:,:,:,:)
TVARRAD(2)%P=>ZRHODREF_RAY(:,:,:,:,:,:)
INVAR=2
! raindrops
IF(SIZE(PRT,4)>2) THEN ! PRT : 4D matrix containing the mixing ratios of different species in the model grid (X,Y,Z)
                       !SIZE(PRT,4) : number of hydrometeor species
  INVAR=INVAR+1
  TVARMOD(INVAR)%P=>PRT(:,:,:,3)
  ALLOCATE(ZR_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))  
  TVARRAD(INVAR)%P=>ZR_RAY(:,:,:,:,:,:)
END IF
! pristine ice
IF (SIZE(PRT,4)>3) THEN
  INVAR=INVAR+1
  TVARMOD(INVAR)%P=>PRT(:,:,:,4)
  ALLOCATE(ZI_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
  TVARRAD(INVAR)%P=>ZI_RAY(:,:,:,:,:,:)
  INVAR=INVAR+1
  TVARMOD(INVAR)%P=>PCIT(:,:,:)
  ALLOCATE(ZCIT_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
  TVARRAD(INVAR)%P=>ZCIT_RAY(:,:,:,:,:,:)
END IF
! snow
IF (SIZE(PRT,4)>4) THEN
  INVAR=INVAR+1
  TVARMOD(INVAR)%P=>PRT(:,:,:,5)
  ALLOCATE(ZS_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
  TVARRAD(INVAR)%P=>ZS_RAY(:,:,:,:,:,:)
END IF
! graupel
IF (SIZE(PRT,4)>5) THEN
  INVAR=INVAR+1
  TVARMOD(INVAR)%P=>PRT(:,:,:,6)
  ALLOCATE(ZG_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
  TVARRAD(INVAR)%P=>ZG_RAY(:,:,:,:,:,:)
END IF
! HAIL
IF (SIZE(PRT,4)>6) THEN
  INVAR=INVAR+1
  TVARMOD(INVAR)%P=>PRT(:,:,:,7)
  ALLOCATE(ZH_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
  TVARRAD(INVAR)%P=>ZH_RAY(:,:,:,:,:,:)
END IF
! convective/stratiform
TVARMOD(INVAR+1)%P=>ZBU_MASK(:,:,:)
ALLOCATE(ZBU_MASK_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
TVARRAD(INVAR+1)%P=>ZBU_MASK_RAY(:,:,:,:,:,:)
! wind components
TVARMOD(INVAR+2)%P=>PUM(:,:,:)
TVARMOD(INVAR+3)%P=>PVM(:,:,:)
TVARMOD(INVAR+4)%P=>PWM(:,:,:)
ALLOCATE(ZUM_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
ALLOCATE(ZVM_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
ALLOCATE(ZWM_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
TVARRAD(INVAR+2)%P=>ZUM_RAY(:,:,:,:,:,:)
TVARRAD(INVAR+3)%P=>ZVM_RAY(:,:,:,:,:,:)
TVARRAD(INVAR+4)%P=>ZWM_RAY(:,:,:,:,:,:)

IF (GLIMA) THEN
  ALLOCATE(ZCRT_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
  TVARMOD(INVAR+5)%P=>PCRT(:,:,:)
  TVARRAD(INVAR+5)%P=>ZCRT_RAY(:,:,:,:,:,:)
ENDIF
!*       3.2   interpolation of all model variables
!interpolation from TVARMOD to TVARRAD of the model variables in the radar projection, using the position (ZX_RAY, ZY_RAY, ZZ_RAY)
!of the beam in the model grid 
CALL INTERPOL_BEAM(TVARMOD,TVARRAD,ZX_RAY(:,:,:,:,:,:),&
     ZY_RAY(:,:,:,:,:,:),ZZ_RAY(:,:,:,:,:,:),ZXHATM(:),ZYHATM(:),ZZM(:,:,:)) 
!
DEALLOCATE(ZBU_MASK)
DEALLOCATE(ZXHATM,ZYHATM,ZZM)
DEALLOCATE(ZX_RAY,ZY_RAY)
DEALLOCATE(TVARMOD,TVARRAD)
!
!Doppler velocities (unfolded): wind contribution
!Calculation of the radial velocity in the radar projection from U,V and W model (ZUM_RAY, ZVM_RAY and ZWM_RAY model)
ALLOCATE(ZVDOP_RAY(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,NPTS_H,NPTS_V))
DO JI=1,NBRAD
  IEL=NBELEV(JI)
  DO JEL=1,IEL
    DO JAZ=1,NBAZIM
      DO JL=1,NBSTEPMAX+1
        DO JH=1,NPTS_H
          DO JV=1,NPTS_V
            IF(ZUM_RAY(JI,JEL,JAZ,JL,JH,JV)/=-XUNDEF) THEN
              ZVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)=(ZUM_RAY(JI,JEL,JAZ,JL,JH,JV)*SIN(ZAZIM(JI,JAZ,JH))&
                       +ZVM_RAY(JI,JEL,JAZ,JL,JH,JV)*COS(ZAZIM(JI,JAZ,JH)))*COS(ZELEV(JI,JEL,JL,JV))&
                       +ZWM_RAY(JI,JEL,JAZ,JL,JH,JV)*SIN(ZELEV(JI,JEL,JL,JV))
            ELSE
              ZVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)=-XUNDEF
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END DO
DEALLOCATE(ZAZIM,ZUM_RAY,ZVM_RAY,ZWM_RAY)
!
WRITE(ILUOUT0,*) 'INTERPOLATION OF MODEL VARIABLES DONE'
!
!-----------------------------------------------------------------------------------------
!*       4.    COMPUTING REFLECTIVITIES ALONG THE RAY BEAM (BACKSCATTERING  + ATTENUATION) 
!              ---------------------------------------------------------------------------    
ALLOCATE(ZZE(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,SIZE(PREFL_CART(:,:,:,:,:),5)))
!ZZE : 5D matrix (iradar, ielev, iaz, irangestep, ivar) containing the radar variables that will be calculated
!in polar or cartesian projection (same projection as the observation grid)
!PREFL_CART is the same matrix, but is the name of the variable which is given as output of radar_simulator
!whereas ZZE is a intermediate variable

!calculation of ZZE and ZBU_MASK_RAY
!-----------------------------------------------------------------------

IF (GLIMA) THEN
  IF (GHAIL) THEN
    CALL RADAR_SCATTERING(ZT_RAY,ZRHODREF_RAY,ZR_RAY,ZI_RAY,ZCIT_RAY,ZS_RAY,ZG_RAY,ZVDOP_RAY, &
     ZELEV,ZX_H,ZX_V,ZW_H,ZW_V,ZZE(:,:,:,:,1:IHAS-1),ZBU_MASK_RAY,PCR_RAY=ZCRT_RAY,PH_RAY=ZH_RAY)
  ELSE
    CALL RADAR_SCATTERING(ZT_RAY,ZRHODREF_RAY,ZR_RAY,ZI_RAY,ZCIT_RAY,ZS_RAY,ZG_RAY,ZVDOP_RAY, &
     ZELEV,ZX_H,ZX_V,ZW_H,ZW_V,ZZE(:,:,:,:,1:IHAS-1),ZBU_MASK_RAY,PCR_RAY=ZCRT_RAY)
  ENDIF
ELSE
  IF (GHAIL) THEN
    CALL RADAR_SCATTERING(ZT_RAY,ZRHODREF_RAY,ZR_RAY,ZI_RAY,ZCIT_RAY,ZS_RAY,ZG_RAY,ZVDOP_RAY, &
     ZELEV,ZX_H,ZX_V,ZW_H,ZW_V,ZZE(:,:,:,:,1:IHAS-1),ZBU_MASK_RAY,PH_RAY=ZH_RAY)
  ELSE
    CALL RADAR_SCATTERING(ZT_RAY,ZRHODREF_RAY,ZR_RAY,ZI_RAY,ZCIT_RAY,ZS_RAY,ZG_RAY,ZVDOP_RAY, &
     ZELEV,ZX_H,ZX_V,ZW_H,ZW_V,ZZE(:,:,:,:,1:IHAS-1),ZBU_MASK_RAY)
  ENDIF
ENDIF
DEALLOCATE(ZVDOP_RAY)
! convective/stratiform
DEALLOCATE(ZBU_MASK_RAY)
! /convective/stratiform
! conversion discretised gates -> single point gates for other output fields
! model variables (beam height, mixing ratios, CIT, refractivity, refractivity gradient) in the radar projection are also added to ZZE
ZZE(:,:,:,:,IHAS)=ZZ_RAY(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2) ! beam height
ALLOCATE(ZWORK(SIZE(ZRHODREF_RAY,1),SIZE(ZRHODREF_RAY,2),&
               SIZE(ZRHODREF_RAY,3),SIZE(ZRHODREF_RAY,4),&
               SIZE(ZRHODREF_RAY,5),SIZE(ZRHODREF_RAY,6)))
! M_r
WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZR_RAY/=-XUNDEF)
  ZWORK=ZRHODREF_RAY*ZR_RAY
ELSEWHERE
  ZWORK=-XUNDEF
END WHERE
WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZR_RAY==XVALGROUND)
  ZWORK=XVALGROUND
END WHERE
ZZE(:,:,:,:,IMR)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
! M_i
WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZI_RAY/=-XUNDEF)
  ZWORK=ZRHODREF_RAY*ZI_RAY
ELSEWHERE
  ZWORK=-XUNDEF
END WHERE
WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZI_RAY==XVALGROUND)
  ZWORK=XVALGROUND
END WHERE
ZZE(:,:,:,:,IMI)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
! M_s
WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZS_RAY/=-XUNDEF)
  ZWORK=ZRHODREF_RAY*ZS_RAY
ELSEWHERE
  ZWORK=-XUNDEF
END WHERE
WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZS_RAY==XVALGROUND)
  ZWORK=XVALGROUND
END WHERE
ZZE(:,:,:,:,IMS)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
! M_g
WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZG_RAY/=-XUNDEF)
  ZWORK=ZRHODREF_RAY*ZG_RAY
ELSEWHERE
  ZWORK=-XUNDEF
END WHERE
WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZG_RAY==XVALGROUND)
  ZWORK=XVALGROUND
END WHERE
ZZE(:,:,:,:,IMG)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
DEALLOCATE(ZWORK)
IF (GHAIL) THEN 
  ALLOCATE(ZWORK(SIZE(ZRHODREF_RAY,1),SIZE(ZRHODREF_RAY,2),&
               SIZE(ZRHODREF_RAY,3),SIZE(ZRHODREF_RAY,4),&
               SIZE(ZRHODREF_RAY,5),SIZE(ZRHODREF_RAY,6)))
  ! M_h
  WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZH_RAY/=-XUNDEF)
    ZWORK=ZRHODREF_RAY*ZH_RAY
  ELSEWHERE
    ZWORK=-XUNDEF
  END WHERE
  WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZH_RAY==XVALGROUND)
    ZWORK=XVALGROUND
  END WHERE
  ZZE(:,:,:,:,IMH)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
  DEALLOCATE(ZWORK)
ENDIF
! CIT
IF (GLIMA)THEN
  ALLOCATE(ZWORK(SIZE(ZRHODREF_RAY,1),SIZE(ZRHODREF_RAY,2),&
                 SIZE(ZRHODREF_RAY,3),SIZE(ZRHODREF_RAY,4),&
                 SIZE(ZRHODREF_RAY,5),SIZE(ZRHODREF_RAY,6)))
  WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZCIT_RAY/=-XUNDEF)
    ZWORK=ZRHODREF_RAY*ZCIT_RAY
  ELSEWHERE
    ZWORK=-XUNDEF
  END WHERE
  WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZCIT_RAY==XVALGROUND)
    ZWORK=XVALGROUND
  END WHERE

  ZZE(:,:,:,:,ICIT)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
  DEALLOCATE(ZWORK)
ELSE
  ZZE(:,:,:,:,ICIT)=ZCIT_RAY(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
ENDIF
! temperature
ZZE(:,:,:,:,ITEM)=ZT_RAY(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2) 
!CRT 
IF (GLIMA)THEN
  ALLOCATE(ZWORK(SIZE(ZRHODREF_RAY,1),SIZE(ZRHODREF_RAY,2),&
                 SIZE(ZRHODREF_RAY,3),SIZE(ZRHODREF_RAY,4),&
                 SIZE(ZRHODREF_RAY,5),SIZE(ZRHODREF_RAY,6)))
  WHERE(ZRHODREF_RAY/=-XUNDEF .AND. ZCRT_RAY/=-XUNDEF)
    ZWORK=ZRHODREF_RAY*ZCRT_RAY
  ELSEWHERE
    ZWORK=-XUNDEF
  END WHERE
  WHERE(ZRHODREF_RAY==XVALGROUND .OR. ZCRT_RAY==XVALGROUND)
    ZWORK=XVALGROUND
  END WHERE
  ZZE(:,:,:,:,ICRT)=ZWORK(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2)
  DEALLOCATE(ZWORK)
ENDIF
IF(LREFR) THEN
   ZZE(:,:,:,:,IRFR)=ZN_RAY(:,:,:,:,(NPTS_H+1)/2,(NPTS_V+1)/2) ! refractivity
   ZZE(:,:,:,NBSTEPMAX+1,IRFR)=XVALGROUND
   DEALLOCATE(ZN_RAY)
END IF
IF(LDNDZ) THEN
   ZZE(:,:,:,:,IDNZ)=ZDNDZ_RAY(:,:,:,:,(NPTS_H+1)/2,(NPTS_v+1)/2) ! refractivity vertical gradient
   DEALLOCATE(ZDNDZ_RAY)
   ZZE(:,:,:,1,IDNZ)=XVALGROUND ! we can do it now
END IF
!
DEALLOCATE(ZELEV)
DEALLOCATE(ZW_H,ZW_V)
DEALLOCATE(ZX_H,ZX_V)
DEALLOCATE(ZT_RAY,ZRHODREF_RAY,ZZ_RAY)
IF(ALLOCATED(ZR_RAY)) DEALLOCATE(ZR_RAY)
IF(ALLOCATED(ZI_RAY)) DEALLOCATE(ZI_RAY,ZCIT_RAY)
IF(ALLOCATED(ZS_RAY)) DEALLOCATE(ZS_RAY)
IF(ALLOCATED(ZG_RAY)) DEALLOCATE(ZG_RAY)
!
!----------------------------------------------------------------------------------------------
!*       5.    INTERPOLATION ON THE CARTESIAN GRID 
!              -----------------------------------
IF (LCART_RAD) THEN !if cartesian interpolation
  ALLOCATE(IREFL_CART_NB(NBRAD,IIELV,2*NMAX,2*NMAX),IVDOP_CART_NB(NBRAD,IIELV,2*NMAX,2*NMAX))
  ALLOCATE(IREFL_CART_NBR(NBRAD,IIELV,2*NMAX,2*NMAX))
  ALLOCATE(IREFL_CART_NBI(NBRAD,IIELV,2*NMAX,2*NMAX))
  ALLOCATE(IREFL_CART_NBS(NBRAD,IIELV,2*NMAX,2*NMAX))
  ALLOCATE(IREFL_CART_NBG(NBRAD,IIELV,2*NMAX,2*NMAX))
  IF (GHAIL) ALLOCATE(IREFL_CART_NBH(NBRAD,IIELV,2*NMAX,2*NMAX))
  PREFL_CART(:,:,:,:,:)=0.
  IREFL_CART_NB(:,:,:,:)=0
  IREFL_CART_NBR(:,:,:,:)=0
  IREFL_CART_NBS(:,:,:,:)=0
  IREFL_CART_NBI(:,:,:,:)=0
  IREFL_CART_NBG(:,:,:,:)=0
  IF (GHAIL)   IREFL_CART_NBH(:,:,:,:)=0
  IVDOP_CART_NB(:,:,:,:)=0
!
!*       5.1  reflectivity on a cartesian grid (this is the way DSO/CMR creates BUFRs)
!  
  DO JI=1,NBRAD
    IEL=NBELEV(JI)
    DO JEL=1,IEL
      DO JAZ=1,NBAZIM
        DO JL=1,NBSTEPMAX+1
          IF ((ZZE(JI,JEL,JAZ,JL,ITEM)/=XVALGROUND).AND.(ZZE(JI,JEL,JAZ,JL,ITEM)/=-XUNDEF)) THEN  !conversion en °C
            ZZE(JI,JEL,JAZ,JL,ITEM)=ZZE(JI,JEL,JAZ,JL,ITEM)-273.15
          ENDIF
          IXGRID=CEILING(NMAX+((JL-1)*XSTEP_RAD*SIN(ZAZIM_BASE(JAZ))/XGRID))
          IYGRID=CEILING(NMAX+((JL-1)*XSTEP_RAD*COS(ZAZIM_BASE(JAZ))/XGRID))
          ! assigning polar grid values to cartesian grid

          !*************************************************
          !****            RADIAL VELOCITY              ****
          !*************************************************
          !XVALGROUND for VRU
          IF(ZZE(JI,JEL,JAZ,JL,IVDOP)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)==XVALGROUND &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
               ) THEN ! if any XVALGROUND in the pixel for ZHH -> pixel set to XVALGROUND for all variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)=XVALGROUND
            IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=1

          !-xundef for VRU
          ELSE IF(ZZE(JI,JEL,JAZ,JL,IVDOP)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)==-XUNDEF &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
               ) THEN ! if any -XUNDEF in the pixel for ZHH-> pixel set to -XUNDEF for all general variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)=-XUNDEF
            IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=1
          
          !if no xvalground and no -xundef for VRU
          ELSE 
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP) &
                 +ZZE(JI,JEL,JAZ,JL,IVDOP)
            IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)+1
          END IF
          
          !*************************************************
          !****            GENERAL VARIABLES            ****
          !*************************************************
          !Keeping -XUNDEF and XVALGROUND values for ZHH and "general" variables
          IF(ZZE(JI,JEL,JAZ,JL,1)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,1)==XVALGROUND &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
               ) THEN ! if any XVALGROUND in the pixel for ZHH -> pixel set to XVALGROUND for all variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,1:4)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHV:IDHV)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IHAS:)=XVALGROUND
            IREFL_CART_NB(JI,JEL,IXGRID,IYGRID)=1
            !IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=1

          !-xundef for ZHH
          ELSE IF(ZZE(JI,JEL,JAZ,JL,1)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,1)==-XUNDEF &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
               ) THEN ! if any -XUNDEF in the pixel for ZHH-> pixel set to -XUNDEF for all general variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,1:4)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHV:IDHV)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IHAS:)=-XUNDEF
            IREFL_CART_NB(JI,JEL,IXGRID,IYGRID)=1
            !IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=1

          !if no xvalground and no -xundef for REFL: incrementation of all polar pixels inside the cartesian pixel
          ELSE 
            PREFL_CART(JI,JEL,IXGRID,IYGRID,1:4)=PREFL_CART(JI,JEL,IXGRID,IYGRID,1:4) &
                 +ZZE(JI,JEL,JAZ,JL,1:4)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHV:IDHV)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHV:IDHV) &
                 +ZZE(JI,JEL,JAZ,JL,IRHV:IDHV)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IHAS:)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IHAS:) &
                 +ZZE(JI,JEL,JAZ,JL,IHAS:)
            IREFL_CART_NB(JI,JEL,IXGRID,IYGRID)=IREFL_CART_NB(JI,JEL,IXGRID,IYGRID)+1
          
            !if no xundef for IVDOP  
            !IF(ZZE(JI,JEL,JAZ,JL,IVDOP)/=-XUNDEF.AND.PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)/=-XUNDEF) THEN
            !  !if no xvalground for IVDOP               
            !  IF (ZZE(JI,JEL,JAZ,JL,IVDOP)/=XVALGROUND.AND.PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)/=XVALGROUND)THEN
            !    PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)+ZZE(JI,JEL,JAZ,JL,IVDOP)
            !    IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)+1
            !  ELSE
            !    PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)=XVALGROUND
            !    IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=1
            !  END IF
            !ELSE
             !  PREFL_CART(JI,JEL,IXGRID,IYGRID,IVDOP)=-XUNDEF
            !   IVDOP_CART_NB(JI,JEL,IXGRID,IYGRID)=1
            !END IF !end if no xundef for IVDOP
          END IF !END IF ZZE(JI,JEL,JAZ,JL,1).OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,1)==XVALGROUND

          !**********************************************
          !****           RAIN VARIABLES             ****
          !**********************************************
          IF(ZZE(JI,JEL,JAZ,JL,IZER)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZER)==XVALGROUND &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
               ) THEN ! if any XVALGROUND in the pixel for ZHH -> pixel set to XVALGROUND for all variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZER)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHR)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDA)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDR)=XVALGROUND
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVR)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATR)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAER)=XVALGROUND
            END IF
            IREFL_CART_NBR(JI,JEL,IXGRID,IYGRID)=1
 
          !-xundef for ZER
          ELSE IF(ZZE(JI,JEL,JAZ,JL,IZER)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZER)==-XUNDEF &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
               ) THEN ! if any -XUNDEF in the pixel for ZER-> pixel set to -XUNDEF for all rain variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZER)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHR)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDA)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDR)=-XUNDEF
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAER)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVR)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATR)=-XUNDEF
            END IF
            IREFL_CART_NBR(JI,JEL,IXGRID,IYGRID)=1

          !if no xvalground and no -xundef for REFL: incrementation of all polar pixels inside the cartesian pixel
          ELSE 
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZER)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZER) &
                 +ZZE(JI,JEL,JAZ,JL,IZER)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHR)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHR) &
                 +ZZE(JI,JEL,JAZ,JL,IRHR)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDA)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDA) &
                 +ZZE(JI,JEL,JAZ,JL,IZDA)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDR)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDR) &
                 +ZZE(JI,JEL,JAZ,JL,IKDR)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAER)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAER) &
                 +ZZE(JI,JEL,JAZ,JL,IAER)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVR)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVR) &
                 +ZZE(JI,JEL,JAZ,JL,IAVR)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATR)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IATR) &
                 +ZZE(JI,JEL,JAZ,JL,IATR)
            END IF
            IREFL_CART_NBR(JI,JEL,IXGRID,IYGRID)=IREFL_CART_NBR(JI,JEL,IXGRID,IYGRID)+1
          END IF


          !**********************************************
          !****           ICE VARIABLES              ****
          !**********************************************
          IF(ZZE(JI,JEL,JAZ,JL,IZEI)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEI)==XVALGROUND &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
               ) THEN ! if any XVALGROUND in the pixel for ZHH -> pixel set to XVALGROUND for all variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEI)=XVALGROUND
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEI)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVI)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATI)=XVALGROUND
            ENDIF
            IREFL_CART_NBI(JI,JEL,IXGRID,IYGRID)=1
 
          !-xundef for ZEI
          ELSE IF(ZZE(JI,JEL,JAZ,JL,IZEI)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEI)==-XUNDEF &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
               ) THEN ! if any -XUNDEF in the pixel for ZER-> pixel set to -XUNDEF for all rain variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEI)=-XUNDEF
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEI)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVI)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATI)=-XUNDEF
            ENDIF
            IREFL_CART_NBI(JI,JEL,IXGRID,IYGRID)=1

          !if no xvalground and no -xundef for REFL: incrementation of all polar pixels inside the cartesian pixel
          ELSE 
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEI)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEI) &
                 +ZZE(JI,JEL,JAZ,JL,IZEI)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEI)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEI) &
                 +ZZE(JI,JEL,JAZ,JL,IAEI)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVI)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVI) &
                 +ZZE(JI,JEL,JAZ,JL,IAVI)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATI)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IATI) &
                 +ZZE(JI,JEL,JAZ,JL,IATI)
            ENDIF
            IREFL_CART_NBI(JI,JEL,IXGRID,IYGRID)=IREFL_CART_NBI(JI,JEL,IXGRID,IYGRID)+1
          END IF

          !**********************************************
          !****           SNOW VARIABLES             ****
          !********************************************** 
          IF(ZZE(JI,JEL,JAZ,JL,IZES)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZES)==XVALGROUND &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
               ) THEN ! if any XVALGROUND in the pixel for ZES -> pixel set to XVALGROUND for all snow variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZES)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHS)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDS)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDS)=XVALGROUND
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAES)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVS)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATS)=XVALGROUND
            END IF
            IREFL_CART_NBS(JI,JEL,IXGRID,IYGRID)=1
 
          !-xundef for ZES
          ELSE IF(ZZE(JI,JEL,JAZ,JL,IZES)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZES)==-XUNDEF &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
               ) THEN ! if any -XUNDEF in the pixel for ZHH-> pixel set to -XUNDEF for all general variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZES)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHS)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDS)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDS)=-XUNDEF
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAES)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVS)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATS)=-XUNDEF
            ENDIF
            IREFL_CART_NBS(JI,JEL,IXGRID,IYGRID)=1

          !if no xvalground and no -xundef for REFL: incrementation of all polar pixels inside the cartesian pixel
          ELSE 
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZES)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZES) &
                 +ZZE(JI,JEL,JAZ,JL,IZES)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHS)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHS) &
                 +ZZE(JI,JEL,JAZ,JL,IRHS)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDS)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDS) &
                 +ZZE(JI,JEL,JAZ,JL,IZDS)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDS)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDS) &
                 +ZZE(JI,JEL,JAZ,JL,IKDS)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAES)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAES) &
                 +ZZE(JI,JEL,JAZ,JL,IAES)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVS)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAES) &
                 +ZZE(JI,JEL,JAZ,JL,IAES)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATS)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IATS) &
                 +ZZE(JI,JEL,JAZ,JL,IATS)
            END IF
            IREFL_CART_NBS(JI,JEL,IXGRID,IYGRID)=IREFL_CART_NBS(JI,JEL,IXGRID,IYGRID)+1
          END IF            

          !**********************************************
          !****          GRAUPEL VARIABLES           ****
          !********************************************** 
          IF(ZZE(JI,JEL,JAZ,JL,IZEG)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEG)==XVALGROUND &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
               ) THEN ! if any XVALGROUND in the pixel for ZES -> pixel set to XVALGROUND for all snow variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEG)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHG)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDG)=XVALGROUND
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDG)=XVALGROUND
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEG)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVG)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATG)=XVALGROUND
            END IF

            IREFL_CART_NBG(JI,JEL,IXGRID,IYGRID)=1
 
          !-xundef for ZEG
          ELSE IF(ZZE(JI,JEL,JAZ,JL,IZEG)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEG)==-XUNDEF &
               .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
               .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
               ) THEN ! if any -XUNDEF in the pixel for ZHH-> pixel set to -XUNDEF for all general variables
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEG)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHG)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDG)=-XUNDEF
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDG)=-XUNDEF
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEG)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVG)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATG)=-XUNDEF
            END IF
            IREFL_CART_NBG(JI,JEL,IXGRID,IYGRID)=1

          !if no xvalground and no -xundef for REFL: incrementation of all polar pixels inside the cartesian pixel
          ELSE 
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEG) &
                 +ZZE(JI,JEL,JAZ,JL,IZEG)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHG) &
                 +ZZE(JI,JEL,JAZ,JL,IRHG)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDG) &
                 +ZZE(JI,JEL,JAZ,JL,IZDG)
            PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDG) &
                 +ZZE(JI,JEL,JAZ,JL,IKDG)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEG) &
                 +ZZE(JI,JEL,JAZ,JL,IAEG)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVG) &
                 +ZZE(JI,JEL,JAZ,JL,IAVG)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IATG)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IATG) &
                 +ZZE(JI,JEL,JAZ,JL,IATG)
            END IF
            IREFL_CART_NBG(JI,JEL,IXGRID,IYGRID)=IREFL_CART_NBG(JI,JEL,IXGRID,IYGRID)+1
          END IF           
          !**********************************************
          !****          HAIL VARIABLES           ****
          !********************************************** 
          IF (GHAIL) THEN
            IF(ZZE(JI,JEL,JAZ,JL,IZEH)==XVALGROUND.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEH)==XVALGROUND &
                 .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==XVALGROUND) &    ! case for refractivity at boundaries
                 .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==XVALGROUND) & ! case for refractivity gradient at origin
                 ) THEN ! if any XVALGROUND in the pixel for ZES -> pixel set to XVALGROUND for all snow variables
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEH)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHH)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDH)=XVALGROUND
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDH)=XVALGROUND
              IF (LATT) THEN
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEH)=XVALGROUND
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVH)=XVALGROUND
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IATH)=XVALGROUND
              END IF

              IREFL_CART_NBG(JI,JEL,IXGRID,IYGRID)=1
 
            !-xundef for ZEG
            ELSE IF(ZZE(JI,JEL,JAZ,JL,IZEH)==-XUNDEF.OR.PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEH)==-XUNDEF &
                 .OR.(LREFR.AND.ZZE(JI,JEL,JAZ,JL,IRFR)==-XUNDEF) &    ! case for refractivity at boundaries
                 .OR.(LDNDZ.AND.ZZE(JI,JEL,JAZ,JL,IDNZ)==-XUNDEF) & ! case for refractivity gradient at origin
                 ) THEN ! if any -XUNDEF in the pixel for ZHH-> pixel set to -XUNDEF for all general variables
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEH)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHH)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDH)=-XUNDEF
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDH)=-XUNDEF
              IF (LATT) THEN
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEH)=-XUNDEF
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVH)=-XUNDEF
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IATH)=-XUNDEF
              END IF
              IREFL_CART_NBG(JI,JEL,IXGRID,IYGRID)=1

            !if no xvalground and no -xundef for REFL: incrementation of all polar pixels inside the cartesian pixel
            ELSE 
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZEH) &
                   +ZZE(JI,JEL,JAZ,JL,IZEH)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IRHH) &
                   +ZZE(JI,JEL,JAZ,JL,IRHH)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IZDH) &
                   +ZZE(JI,JEL,JAZ,JL,IZDH)
              PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IKDH) &
                   +ZZE(JI,JEL,JAZ,JL,IKDH)
              IF (LATT) THEN
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAEH) &
                   +ZZE(JI,JEL,JAZ,JL,IAEH)
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IAVH) &
                   +ZZE(JI,JEL,JAZ,JL,IAVH)
                PREFL_CART(JI,JEL,IXGRID,IYGRID,IATH)=PREFL_CART(JI,JEL,IXGRID,IYGRID,IATH) &
                   +ZZE(JI,JEL,JAZ,JL,IATH)
              END IF
              IREFL_CART_NBH(JI,JEL,IXGRID,IYGRID)=IREFL_CART_NBH(JI,JEL,IXGRID,IYGRID)+1
            END IF
          END IF 
        END DO !JL
      END DO !JAZ
    END DO !JEL
  END DO !JI

DEALLOCATE(ZZE)
!*       5.2  writing cartesian grid output (averaging) and dB conversion
  DO JI=1,NBRAD
    IEL=NBELEV(JI)
    DO JEL=1,IEL
      DO JV=2*NMAX,1,-1
        DO JH=1,2*NMAX

          !****   RADIAL VELOCITY
          IF(IVDOP_CART_NB(JI,JEL,JH,JV) == 0) THEN
            IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
              ! out of range
              PREFL_CART(JI,JEL,JH,JV,IVDOP)=XVALGROUND
            ELSE
              PREFL_CART(JI,JEL,JH,JV,IVDOP)=0
              WRITE(ILUOUT0,*) "Warning: some pixels have no radial velocity; increase XGRID or decrease XSTEP_RAD"
            END IF
          ELSE            
            PREFL_CART(JI,JEL,JH,JV,IVDOP)=PREFL_CART(JI,JEL,JH,JV,IVDOP)/IVDOP_CART_NB(JI,JEL,JH,JV)
          END IF

          !****   GENERAL VARIABLES
          IF(IREFL_CART_NB(JI,JEL,JH,JV) == 0) THEN ! no values inside the cartesian pixel
            IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
               ! out of range
              PREFL_CART(JI,JEL,JH,JV,1:4)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IRHV:IDHV)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IHAS:)=XVALGROUND
            ELSE
              PREFL_CART(JI,JEL,JH,JV,1:4)=0
              PREFL_CART(JI,JEL,JH,JV,IRHV:IDHV)=0
              PREFL_CART(JI,JEL,JH,JV,IHAS:)=0
              WRITE(ILUOUT0,*) "Warning: some pixels have no reflectivity; increase XGRID or decrease XSTEP_RAD"
            END IF
          ELSE
            !PREFL_CART(JI,JEL,JH,JV,:)=PREFL_CART(JI,JEL,JH,JV,:)/IREFL_CART_NB(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,1:4)=PREFL_CART(JI,JEL,JH,JV,1:4)/IREFL_CART_NB(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IRHV:IDHV)=PREFL_CART(JI,JEL,JH,JV,IRHV:IDHV)/IREFL_CART_NB(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IHAS:)=PREFL_CART(JI,JEL,JH,JV,IHAS:)/IREFL_CART_NB(JI,JEL,JH,JV)
                     
            ! --------- Converting to dB -----------
            IF(PREFL_CART(JI,JEL,JH,JV,1)> 0) THEN
              PREFL_CART(JI,JEL,JH,JV,1)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,1)) ! Z_equiv in dBZ
              IF(PREFL_CART(JI,JEL,JH,JV,2) > 0.) THEN
                PREFL_CART(JI,JEL,JH,JV,2)=PREFL_CART(JI,JEL,JH,JV,1) &
                                          -10.*LOG10(PREFL_CART(JI,JEL,JH,JV,2)) ! Zdr=Z_HH-Z_VV  
              !ELSE
              !  PREFL_CART(JI,JEL,JH,JV,2)=-XUNDEF Zvv<0
              ENDIF                      
            ELSE IF (PREFL_CART(JI,JEL,JH,JV,1)== 0) THEN
              PREFL_CART(JI,JEL,JH,JV,1)=-XUNDEF 
              !PREFL_CART(JI,JEL,JH,JV,IZER:IZEG)=-XUNDEF 
            END IF    
          END IF !***** END GENERAL VARIABLES

          !****   RAIN VARIABLES
          IF(IREFL_CART_NBR(JI,JEL,JH,JV) == 0) THEN ! no values inside the cartesian pixel
            IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
              ! out of range
              PREFL_CART(JI,JEL,JH,JV,IZER)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IRHR)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IZDA)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IKDR)=XVALGROUND
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAER)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IAVR)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IATR)=XVALGROUND
              END IF
            ELSE
              PREFL_CART(JI,JEL,JH,JV,IZER)=0
              PREFL_CART(JI,JEL,JH,JV,IRHR)=0
              PREFL_CART(JI,JEL,JH,JV,IZDA)=0
              PREFL_CART(JI,JEL,JH,JV,IKDR)=0
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAER)=0
                PREFL_CART(JI,JEL,JH,JV,IAVR)=0
                PREFL_CART(JI,JEL,JH,JV,IATR)=0
              END IF
              WRITE(ILUOUT0,*) "Warning: some pixels have no reflectivity; increase XGRID or decrease XSTEP_RAD"
            END IF
          ELSE
            PREFL_CART(JI,JEL,JH,JV,IZER)=PREFL_CART(JI,JEL,JH,JV,IZER)/IREFL_CART_NBR(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IRHR)=PREFL_CART(JI,JEL,JH,JV,IRHR)/IREFL_CART_NBR(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IZDA)=PREFL_CART(JI,JEL,JH,JV,IZDA)/IREFL_CART_NBR(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IKDR)=PREFL_CART(JI,JEL,JH,JV,IKDR)/IREFL_CART_NBR(JI,JEL,JH,JV)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,JH,JV,IAER)=PREFL_CART(JI,JEL,JH,JV,IAER)/IREFL_CART_NBR(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IAVR)=PREFL_CART(JI,JEL,JH,JV,IAVR)/IREFL_CART_NBR(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IATR)=PREFL_CART(JI,JEL,JH,JV,IATR)/IREFL_CART_NBR(JI,JEL,JH,JV)
            END IF
            
            ! --------- Converting to dB -----------
            IF(PREFL_CART(JI,JEL,JH,JV,IZER)> 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZER)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZER)) ! Z_equiv in dBZ
              IF(PREFL_CART(JI,JEL,JH,JV,IZDA) > 0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IZDA)=PREFL_CART(JI,JEL,JH,JV,IZER) &
                                              -10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZDA)) ! Zdr=Z_HH-Z_VV  
              ENDIF                      
            ELSE IF (PREFL_CART(JI,JEL,JH,JV,IZER)== 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZER)=-XUNDEF 
            END IF    

            IF(LATT) THEN
              IF(PREFL_CART(JI,JEL,JH,JV,IATR)>0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IATR)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IATR))
              ENDIF
            ENDIF
          END IF !****** END RAIN

          !****   ICE VARIABLES
          IF(IREFL_CART_NBI(JI,JEL,JH,JV) == 0) THEN ! no values inside the cartesian pixel
            IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
              ! out of range
              PREFL_CART(JI,JEL,JH,JV,IZEI)=XVALGROUND
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAEI)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IAVI)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IATI)=XVALGROUND
              END IF
            ELSE
              PREFL_CART(JI,JEL,JH,JV,IZEI)=0
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAEI)=0
                PREFL_CART(JI,JEL,JH,JV,IAVI)=0
                PREFL_CART(JI,JEL,JH,JV,IATI)=0
              END IF
              WRITE(ILUOUT0,*) "Warning: some pixels have no reflectivity; increase XGRID or decrease XSTEP_RAD"
            END IF
          ELSE
            PREFL_CART(JI,JEL,JH,JV,IZEI)=PREFL_CART(JI,JEL,JH,JV,IZEI)/IREFL_CART_NBI(JI,JEL,JH,JV)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,JH,JV,IAEI)=PREFL_CART(JI,JEL,JH,JV,IAEI)/IREFL_CART_NBI(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IAVI)=PREFL_CART(JI,JEL,JH,JV,IAVI)/IREFL_CART_NBI(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IATI)=PREFL_CART(JI,JEL,JH,JV,IATI)/IREFL_CART_NBI(JI,JEL,JH,JV)
            END IF
            ! --------- Converting to dB -----------
            IF(PREFL_CART(JI,JEL,JH,JV,IZEI)> 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZEI)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZEI)) ! Z_equiv in dB                    
            ELSE IF (PREFL_CART(JI,JEL,JH,JV,IZEI)== 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZEI)=-XUNDEF 
            END IF    

            IF(LATT) THEN
              IF(PREFL_CART(JI,JEL,JH,JV,IATI)>0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IATI)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IATI))
              ENDIF
            ENDIF
          END IF !****** END ICE

          !****   SNOW VARIABLES
          IF(IREFL_CART_NBS(JI,JEL,JH,JV) == 0) THEN ! no values inside the cartesian pixel
            IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
               ! out of range
              PREFL_CART(JI,JEL,JH,JV,IZES)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IRHS)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IZDS)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IKDS)=XVALGROUND
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAES)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IAVS)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IATS)=XVALGROUND
              ENDIF
            ELSE
              PREFL_CART(JI,JEL,JH,JV,IZES)=0
              PREFL_CART(JI,JEL,JH,JV,IRHS)=0
              PREFL_CART(JI,JEL,JH,JV,IZDS)=0
              PREFL_CART(JI,JEL,JH,JV,IKDS)=0
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAES)=0
                PREFL_CART(JI,JEL,JH,JV,IAVS)=0
                PREFL_CART(JI,JEL,JH,JV,IATS)=0
              END IF
              WRITE(ILUOUT0,*) "Warning: some pixels have no reflectivity; increase XGRID or decrease XSTEP_RAD"
            END IF
          ELSE
            PREFL_CART(JI,JEL,JH,JV,IZES)=PREFL_CART(JI,JEL,JH,JV,IZES)/IREFL_CART_NBS(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IRHS)=PREFL_CART(JI,JEL,JH,JV,IRHS)/IREFL_CART_NBS(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IZDS)=PREFL_CART(JI,JEL,JH,JV,IZDS)/IREFL_CART_NBS(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IKDS)=PREFL_CART(JI,JEL,JH,JV,IKDS)/IREFL_CART_NBS(JI,JEL,JH,JV)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,JH,JV,IAES)=PREFL_CART(JI,JEL,JH,JV,IAES)/IREFL_CART_NBS(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IAVS)=PREFL_CART(JI,JEL,JH,JV,IAVS)/IREFL_CART_NBS(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IATS)=PREFL_CART(JI,JEL,JH,JV,IATS)/IREFL_CART_NBS(JI,JEL,JH,JV)
            END IF
                      
            ! --------- Converting to dB -----------
            IF(PREFL_CART(JI,JEL,JH,JV,IZES)> 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZES)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZES)) ! Z_equiv in dBZ
              IF(PREFL_CART(JI,JEL,JH,JV,IZDS) > 0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IZDS)=PREFL_CART(JI,JEL,JH,JV,IZES) &
                                             -10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZDS)) ! Zdr=Z_HH-Z_VV  
              ENDIF                      
            ELSE IF (PREFL_CART(JI,JEL,JH,JV,IZES)== 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZES)=-XUNDEF 
            END IF    

            IF(LATT) THEN
              IF(PREFL_CART(JI,JEL,JH,JV,IATS)>0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IATS)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IATS))
              ENDIF
            ENDIF
          END IF !******** END SNOW

          !****   GRAUPEL VARIABLES
          IF(IREFL_CART_NBG(JI,JEL,JH,JV) == 0) THEN ! no values inside the cartesian pixel
            IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
               ! out of range
              PREFL_CART(JI,JEL,JH,JV,IZEG)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IRHG)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IZDG)=XVALGROUND
              PREFL_CART(JI,JEL,JH,JV,IKDG)=XVALGROUND
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAEG)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IAVG)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IATG)=XVALGROUND
              END IF
            ELSE
              PREFL_CART(JI,JEL,JH,JV,IZEG)=0
              PREFL_CART(JI,JEL,JH,JV,IRHG)=0
              PREFL_CART(JI,JEL,JH,JV,IZDG)=0
              PREFL_CART(JI,JEL,JH,JV,IKDG)=0
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAEG)=0
                PREFL_CART(JI,JEL,JH,JV,IAVG)=0
                PREFL_CART(JI,JEL,JH,JV,IATG)=0
              ENDIF
              WRITE(ILUOUT0,*) "Warning: some pixels have no reflectivity; increase XGRID or decrease XSTEP_RAD"
            END IF
          ELSE
            PREFL_CART(JI,JEL,JH,JV,IZEG)=PREFL_CART(JI,JEL,JH,JV,IZEG)/IREFL_CART_NBG(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IRHG)=PREFL_CART(JI,JEL,JH,JV,IRHG)/IREFL_CART_NBG(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IZDG)=PREFL_CART(JI,JEL,JH,JV,IZDG)/IREFL_CART_NBG(JI,JEL,JH,JV)
            PREFL_CART(JI,JEL,JH,JV,IKDG)=PREFL_CART(JI,JEL,JH,JV,IKDG)/IREFL_CART_NBG(JI,JEL,JH,JV)
            IF (LATT) THEN
              PREFL_CART(JI,JEL,JH,JV,IAEG)=PREFL_CART(JI,JEL,JH,JV,IAEG)/IREFL_CART_NBG(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IAVG)=PREFL_CART(JI,JEL,JH,JV,IAVG)/IREFL_CART_NBG(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IATG)=PREFL_CART(JI,JEL,JH,JV,IATG)/IREFL_CART_NBG(JI,JEL,JH,JV)
            END IF
            
            ! --------- Converting to dB -----------
            IF(PREFL_CART(JI,JEL,JH,JV,IZEG)> 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZEG)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZEG)) ! Z_equiv in dBZ
              IF(PREFL_CART(JI,JEL,JH,JV,IZDG) > 0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IZDG)=PREFL_CART(JI,JEL,JH,JV,IZEG) &
                                             -10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZDG)) ! Zdr=Z_HH-Z_VV  
              ENDIF                      
            ELSE IF (PREFL_CART(JI,JEL,JH,JV,IZEG)== 0) THEN
              PREFL_CART(JI,JEL,JH,JV,IZES)=-XUNDEF 
            END IF    

            IF(LATT) THEN
              IF(PREFL_CART(JI,JEL,JH,JV,IATG)>0.) THEN
                PREFL_CART(JI,JEL,JH,JV,IATG)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IATG))
              ENDIF
            ENDIF
          END IF !******** END GRAUPEL
           !****   HAIL VARIABLES
          IF (GHAIL) THEN
            IF(IREFL_CART_NBH(JI,JEL,JH,JV) == 0) THEN ! no values inside the cartesian pixel
              IF((JH+SIGN(.5,JH-.5-NMAX)-.5-NMAX)**2+(JV+SIGN(.5,JV-.5-NMAX)-.5-NMAX)**2>=NMAX**2) THEN 
                ! out of range
                PREFL_CART(JI,JEL,JH,JV,IZEH)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IRHH)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IZDH)=XVALGROUND
                PREFL_CART(JI,JEL,JH,JV,IKDH)=XVALGROUND
                IF (LATT) THEN
                  PREFL_CART(JI,JEL,JH,JV,IAEH)=XVALGROUND
                  PREFL_CART(JI,JEL,JH,JV,IAVH)=XVALGROUND
                  PREFL_CART(JI,JEL,JH,JV,IATH)=XVALGROUND
                END IF
              ELSE
                PREFL_CART(JI,JEL,JH,JV,IZEH)=0
                PREFL_CART(JI,JEL,JH,JV,IRHH)=0
                PREFL_CART(JI,JEL,JH,JV,IZDH)=0
                PREFL_CART(JI,JEL,JH,JV,IKDH)=0
                IF (LATT) THEN
                  PREFL_CART(JI,JEL,JH,JV,IAEH)=0
                  PREFL_CART(JI,JEL,JH,JV,IAVH)=0
                  PREFL_CART(JI,JEL,JH,JV,IATH)=0
                ENDIF
                WRITE(ILUOUT0,*) "Warning: some pixels have no reflectivity; increase XGRID or decrease XSTEP_RAD"
              END IF
            ELSE
              PREFL_CART(JI,JEL,JH,JV,IZEH)=PREFL_CART(JI,JEL,JH,JV,IZEH)/IREFL_CART_NBH(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IRHH)=PREFL_CART(JI,JEL,JH,JV,IRHH)/IREFL_CART_NBH(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IZDH)=PREFL_CART(JI,JEL,JH,JV,IZDH)/IREFL_CART_NBH(JI,JEL,JH,JV)
              PREFL_CART(JI,JEL,JH,JV,IKDH)=PREFL_CART(JI,JEL,JH,JV,IKDH)/IREFL_CART_NBH(JI,JEL,JH,JV)
              IF (LATT) THEN
                PREFL_CART(JI,JEL,JH,JV,IAEH)=PREFL_CART(JI,JEL,JH,JV,IAEH)/IREFL_CART_NBH(JI,JEL,JH,JV)
                PREFL_CART(JI,JEL,JH,JV,IAVH)=PREFL_CART(JI,JEL,JH,JV,IAVH)/IREFL_CART_NBH(JI,JEL,JH,JV)
                PREFL_CART(JI,JEL,JH,JV,IATH)=PREFL_CART(JI,JEL,JH,JV,IATH)/IREFL_CART_NBH(JI,JEL,JH,JV)
              END IF
            
              ! --------- Converting to dB -----------
              IF(PREFL_CART(JI,JEL,JH,JV,IZEH)> 0) THEN
                PREFL_CART(JI,JEL,JH,JV,IZEH)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZEH)) ! Z_equiv in dBZ
                IF(PREFL_CART(JI,JEL,JH,JV,IZDH) > 0.) THEN
                  PREFL_CART(JI,JEL,JH,JV,IZDH)=PREFL_CART(JI,JEL,JH,JV,IZEH) &
                                             -10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IZDH)) ! Zdr=Z_HH-Z_VV  
                ENDIF                      
              ELSE IF (PREFL_CART(JI,JEL,JH,JV,IZEH)== 0) THEN
                PREFL_CART(JI,JEL,JH,JV,IZES)=-XUNDEF 
              END IF    

              IF(LATT) THEN
                IF(PREFL_CART(JI,JEL,JH,JV,IATH)>0.) THEN
                  PREFL_CART(JI,JEL,JH,JV,IATH)=10.*LOG10(PREFL_CART(JI,JEL,JH,JV,IATH))
                ENDIF
              ENDIF
            END IF !******** END HAIL
          END IF
        END DO
      END DO
    END DO
  END DO
  WRITE(ILUOUT0,*) 'CARTESIAN GRID INTERPOLATION DONE'
  DEALLOCATE(IREFL_CART_NB,IVDOP_CART_NB)
!
!*       5.3  positions of the cartesian grid (as in R2, provided 'as is')
!
  DO JI=1,NBRAD
     DO JEL=1,IEL
        DO JV=1,2*NMAX
           DO JH=1,2*NMAX
              r=SQRT((JH-.5-NMAX)**2+(JV-.5-NMAX)**2)*XGRID
              h=XRADIUS+r*XELEV(JI,JEL)*ZRDSDG+r*r/(2.*ZKE*XRADIUS)
              alph=ACOS((XRADIUS*XRADIUS+h*h-r*r)/(2.*XRADIUS*h))
              
              PLATLON(JI,2*JH-1,JV)=ASIN(SIN(XLAT_RAD(JI)*ZRDSDG)*COS(alph)+COS(XLAT_RAD(JI)*ZRDSDG)*SIN(alph)* & ! LAT
                   (JV-.5-NMAX)/SQRT((JH-.5-NMAX)**2+(JV-.5-NMAX)**2))/ZRDSDG
              PLATLON(JI,2*JH,JV)=XLON_RAD(JI)+ASIN(SIN(alph)* & ! lon
                   (JH-.5-NMAX)/SQRT((JH-.5-NMAX)**2+(JV-.5-NMAX)**2)/ &
                   COS(PLATLON(JI,2*JH-1,JV)*ZRDSDG))/ZRDSDG
           END DO
        END DO
     END DO
  END DO


ELSE ! if polar output
  PREFL_CART(:,:,:,:,:)=ZZE(:,:,:,:,:)

  DEALLOCATE(ZZE)

  DO JI=1,NBRAD
     DO JAZ=1,NBAZIM
        PLATLON(1,JAZ,1)=ZAZIM_BASE(JAZ) ! pourquoi PLATLON(1,JAZ,1) et pas PLATLON(JI,JAZ,1)? 
     END DO ! en coordonnees polaires on indique la position centrale (en radian) de chaque azimut
  END DO


  !--------- Converting to dB for polar output -----------
  DO JI=1,NBRAD
    IEL=NBELEV(JI)
    DO JEL=1,IEL
      DO JAZ=1,NBAZIM
        DO JL=1,NBSTEPMAX+1
          !conversion en deg celcius
          IF (PREFL_CART(JI,JEL,JAZ,JL,ITEM)/=-XUNDEF .AND. PREFL_CART(JI,JEL,JAZ,JL,ITEM)/=XVALGROUND) &
          PREFL_CART(JI,JEL,JAZ,JL,ITEM)=PREFL_CART(JI,JEL,JAZ,JL,ITEM)-273.15
          
          !------ ZHH and ZDR  
          IF(PREFL_CART(JI,JEL,JAZ,JL,1)> 0) THEN
            PREFL_CART(JI,JEL,JAZ,JL,1)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,1)) ! Z_equiv in dBZ
            IF(PREFL_CART(JI,JEL,JAZ,JL,2) > 0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,2)=PREFL_CART(JI,JEL,JAZ,JL,1) &
                                          -10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,2)) ! Zdr=Z_HH-Z_VV  
            ENDIF                      
          ELSE IF (PREFL_CART(JI,JEL,JAZ,JL,1)== 0) THEN
              PREFL_CART(JI,JEL,JAZ,JL,1)=-XUNDEF 
          END IF
          
          !------ RAIN : ZER, ZDA, ATR
          IF(PREFL_CART(JI,JEL,JAZ,JL,IZER)> 0) THEN
            PREFL_CART(JI,JEL,JAZ,JL,IZER)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZER)) ! Z_equiv in dBZ
            IF(PREFL_CART(JI,JEL,JAZ,JL,IZDA) > 0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IZDA)=PREFL_CART(JI,JEL,JAZ,JL,IZER) &
                                            -10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZDA)) ! Zdr=Z_HH-Z_VV  
            ENDIF                      
          ELSE IF (PREFL_CART(JI,JEL,JAZ,JL,IZER)== 0) THEN
             PREFL_CART(JI,JEL,JAZ,JL,IZER)=-XUNDEF 
          END IF    

          IF(LATT) THEN
            IF(PREFL_CART(JI,JEL,JAZ,JL,IATR)>0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IATR)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IATR))
            ENDIF 
          ENDIF !------- END RAIN

! ------- ICE  -----------
          IF(PREFL_CART(JI,JEL,JAZ,JL,IZEI)> 0) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IZEI)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZEI)) ! Z_equiv in dB                    
          ELSE IF (PREFL_CART(JI,JEL,JAZ,JL,IZEI)== 0) THEN
             PREFL_CART(JI,JEL,JAZ,JL,IZEI)=-XUNDEF 
          END IF    

          IF(LATT) THEN
            IF(PREFL_CART(JI,JEL,JAZ,JL,IATI)>0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IATI)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IATI))
            ENDIF
          END IF !------ END ICE

! --------- SNOW   -----------
          IF(PREFL_CART(JI,JEL,JAZ,JL,IZES)> 0) THEN
            PREFL_CART(JI,JEL,JAZ,JL,IZES)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZES)) ! Z_equiv in dBZ
            IF(PREFL_CART(JI,JEL,JAZ,JL,IZDS) > 0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IZDS)=PREFL_CART(JI,JEL,JAZ,JL,IZES) &
                                             -10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZDS)) ! Zdr=Z_HH-Z_VV  
            ENDIF                      
          ELSE IF (PREFL_CART(JI,JEL,JAZ,JL,IZES)== 0) THEN
            PREFL_CART(JI,JEL,JAZ,JL,IZES)=-XUNDEF 
          END IF    

          IF(LATT) THEN
            IF(PREFL_CART(JI,JEL,JAZ,JL,IATS)>0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IATS)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IATS))
            ENDIF
          END IF !------ END SNOW   

! --------- GRAUPEL  -----------
          IF(PREFL_CART(JI,JEL,JAZ,JL,IZEG)> 0) THEN
            PREFL_CART(JI,JEL,JAZ,JL,IZEG)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZEG)) ! Z_equiv in dBZ
            IF(PREFL_CART(JI,JEL,JAZ,JL,IZDG) > 0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IZDG)=PREFL_CART(JI,JEL,JAZ,JL,IZEG) &
                                            -10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZDG)) ! Zdr=Z_HH-Z_VV  
            ENDIF                      
          ELSE IF (PREFL_CART(JI,JEL,JAZ,JL,IZEG)== 0) THEN
            PREFL_CART(JI,JEL,JAZ,JL,IZES)=-XUNDEF 
          END IF    

          IF(LATT) THEN
            IF(PREFL_CART(JI,JEL,JAZ,JL,IATG)>0.) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IATG)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IATG))
            END IF
          END IF !------ END GRAUPEL
! --------- HAIL  -----------
          IF (GHAIL) THEN
            IF(PREFL_CART(JI,JEL,JAZ,JL,IZEH)> 0) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IZEH)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZEH)) ! Z_equiv in dBZ
              IF(PREFL_CART(JI,JEL,JAZ,JL,IZDH) > 0.) THEN
                PREFL_CART(JI,JEL,JAZ,JL,IZDH)=PREFL_CART(JI,JEL,JAZ,JL,IZEH) &
                                            -10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IZDH)) ! Zdr=Z_HH-Z_VV  
              ENDIF                      
            ELSE IF (PREFL_CART(JI,JEL,JAZ,JL,IZEH)== 0) THEN
              PREFL_CART(JI,JEL,JAZ,JL,IZES)=-XUNDEF 
            END IF    

            IF(LATT) THEN
              IF(PREFL_CART(JI,JEL,JAZ,JL,IATH)>0.) THEN
                PREFL_CART(JI,JEL,JAZ,JL,IATH)=10.*LOG10(PREFL_CART(JI,JEL,JAZ,JL,IATH))
              END IF
            END IF !------ END HAIL
          END IF
        END DO
      END DO
    END DO
  END DO


END IF ! end condition on cartesian or polar outout
DEALLOCATE(ZAZIM_BASE)
WRITE(ILUOUT0,*) 'ROUTINE RADAR_SIMULATOR COMPLETED'
END SUBROUTINE RADAR_SIMULATOR

