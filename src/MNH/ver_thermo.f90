!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_VER_THERMO
!     ######################
INTERFACE
      SUBROUTINE VER_THERMO(TPFILE,OSHIFT,                                                 &
                            PTHV_MX,PR_MX,PZS_LS,PZSMT_LS,PZMASS_MX,PZFLUX_MX,PPMHP_MX,PJ, &
                            PDXX,PDYY,PEXNTOP2D,PPSURF,PDIAG,                              &
                            PLSTH_MX,PLSRV_MX                                              )
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),          INTENT(IN)  :: TPFILE     ! File characteristics
LOGICAL,                  INTENT(IN)  :: OSHIFT     ! T: vertical shift of BL (used for GRIB file data)
!                                                   ! F: no vertical shift (used for MESONH data)
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PTHV_MX   ! thetav on mixed grid
REAL,   DIMENSION(:,:,:,:), INTENT(IN)   :: PR_MX     ! r on mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! mass point altitudes on
!                                                  ! mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX_MX ! flux point altitudes on
!                                                  ! mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PPMHP_MX  ! pressure minus hyd. pressure
REAL,   DIMENSION(:,:)  , INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:)  , INTENT(IN)  :: PZSMT_LS  ! large scale smooth orography
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! Jacobian
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! metric coefficient dyy
REAL,   DIMENSION(:,:)  , INTENT(IN)     :: PEXNTOP2D ! top Exner function
REAL,   DIMENSION(:,:)  , INTENT(OUT)     :: PPSURF    ! Surface pressure
REAL,                     INTENT(OUT) :: PDIAG     ! diagnostics computing time
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSTH_MX ! large scale potential temperature
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSRV_MX ! large scale vapor mixing ratios
!
END SUBROUTINE VER_THERMO
END INTERFACE
END MODULE MODI_VER_THERMO
!     ######################################################################
      SUBROUTINE VER_THERMO(TPFILE,OSHIFT,                                                 &
                            PTHV_MX,PR_MX,PZS_LS,PZSMT_LS,PZMASS_MX,PZFLUX_MX,PPMHP_MX,PJ, &
                            PDXX,PDYY,PEXNTOP2D,PPSURF,PDIAG,                              &
                            PLSTH_MX,PLSRV_MX                                              )
!     ######################################################################
!
!!****  *VER_THERMO* - initializes the thermodynamic and reference state
!!                     fields in MESO-NH for a real case from Aladin fields.
!!
!!    PURPOSE
!!    -------
!!    This routine initializes the potential temperature and the vapor mixing
!!    ratio on the MESO-NH grid from the fields on the mixed grid
!!    defined by the altitudes of its mass points.
!!    The other mixing ratio, if any, are initialized to zero.
!!    The reference state variables and the total dry mass are also computed.
!!
!!**  METHOD
!!    ------
!!
!!  1 The initialization of thetav and rv is performed in VER_INT_THERMO
!!
!!  2 theta is deduced:
!!                   1 + rv
!!   theta=thetav *------------------
!!                   1 + Rv/Rd*rv
!!
!!  3 The reference anelastic state variables are computed in SET_REFZ
!!
!!  4 The total dry mass is computed in DRY_MASS
!!    ____
!!  5 rhod J* (XRHODJ) is computed in SET_REF
!!
!!  6 theta and rv are stored in XTHT and XRT
!!
!!    EXTERNAL
!!    --------
!!    subroutine VER_INT_THERMO : to initialize thetav and rv
!!    subroutine SET_REFZ       : to initialize the reference state 1D variables
!!    subroutine TOTAL_DMASS    : to compute the total dry mass
!!    subroutine SET_REF        : to compute  rhoJ
!!
!!    Module MODI_VER_INT_THERMO: interface for subroutine VER_INT_THERMO
!!    Module MODI_SET_REFZ      : interface for subroutine SET_REFZ
!!    Module MODI_TOTAL_DMASS   : interface for subroutine TOTAL_DMASS
!!    Module MODI_SET_REF       : interface for subroutine SET_REF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB   : verbosity level for output-listing
!!      Module MODD_CONF1     : contains configuration variables for model 1.
!!         NRR     : number of moist variables
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_CST       : contains physical constants
!!         XRD : gas constant for dry air
!!         XRV : gas constant for vapor
!!         XG  : gravity constant
!!      Module MODD_GRID1     : contains grid variables
!!         XZZ :
!!         XZHAT:
!!      Module MODD_REF1      : contains 3D reference state variables for model1
!!         XRHODJ      :
!!         XTHVREF     :
!!         XEXNREF     :
!!         XREFMASS    : Mass of the ref. atmosphere contained in the simulation
!                        domain
!!         XMASS_O_PHI0: normalization constant used in the PHI0 computation
!!         XLINMASS    : Lineic mass through open boundaries
!!      Module MODD_FIELD1    : contains prognostics  variables
!!         XTHM : theta (:,:,:)         at t-dt
!!         XRM  : r (:,:,:,:)           at t-dt
!!      Module MODD_LBC1 : contains declaration of lateral boundary conditions
!!         CLBCX   : X-direction LBC type at left(1) and right(2) boundaries
!!         CLBCY   : Y-direction LBC type at left(1) and right(2) boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/12/94
!!                  Sept. 21, 1995  (J.Stein and V.Masson) surface pressure
!!                  Jan.  09, 1996  (V. Masson) pressure function deduced from
!!                                  hydrostatic pressure
!!                  Aug.  20, 1996  (V. Masson) change call to DRY_MASS,
!!                                  correction in XTHM computation
!!                  Oct.  25, 1996  (V. Masson) add deallocations
!!                  Jan   15, 1997  (Stein,Lafore) Durran anelastic equation
!!                  Jun   10, 1997  (V. Masson) add NH pressure
!!                  Jan.  25, 1998  (Stein) add the LS fields' treatment
!!                  Jun.  06, 2006  (Mallet) replace DRY_MASS by TOTAL_DMASS
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                            2014  (M.Faivre)
!!                         08/2015  (M.Moge)    removing UPDATE_HALO_ll on 
!!                                              XRHODREF, XTHVREF, XRVREF, XEXNREF, XRHODJ
!!                                              because we now do it in SET_REF
!!                     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ARGSLIST_ll, ONLY: LIST_ll
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_DYN_n
USE MODD_FIELD_n,     ONLY: XTHT,XRT,XPABST,XDRYMASST
USE MODD_GRID_n
USE MODD_IO_ll,       ONLY: TFILEDATA,TFILE_DUMMY
USE MODD_LBC_n
USE MODD_LSFIELD_n
USE MODD_LUNIT,       ONLY: CLUOUT0,TLUOUT0
USE MODD_LUNIT_n,     ONLY: CLUOUT
USE MODD_PARAMETERS
USE MODD_REF_n
!
USE MODD_DIM_n
USE MODE_EXTRAPOL
USE MODE_FIELD,       ONLY: TFIELDDATA,TYPEREAL
USE MODE_FMWRIT
USE MODE_ll
USE MODE_MPPDB
!
USE MODI_COMPUTE_EXNER_FROM_TOP
USE MODI_SET_REF
USE MODI_SET_REFZ
USE MODI_TOTAL_DMASS
USE MODI_VER_INT_THERMO
USE MODI_WATER_SUM
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
TYPE(TFILEDATA),          INTENT(IN)  :: TPFILE     ! File characteristics
LOGICAL,                  INTENT(IN)  :: OSHIFT     ! T: vertical shift of BL (used for GRIB file data)
!                                                   ! F: no vertical shift (used for MESONH data)
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PTHV_MX   ! thetav on mixed grid
REAL,   DIMENSION(:,:,:,:), INTENT(IN)   :: PR_MX     ! r on mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! mass point altitudes on
!                                                  ! mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX_MX ! flux point altitudes on
!                                                  ! mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PPMHP_MX  ! pressure minus hyd. pressure
REAL,   DIMENSION(:,:)  , INTENT(IN)  :: PZS_LS    ! mixed grid orography
REAL,   DIMENSION(:,:)  , INTENT(IN)  :: PZSMT_LS  ! large scale smooth orography
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! Jacobian
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! metric coefficient dyy
REAL,   DIMENSION(:,:)  , INTENT(IN)     :: PEXNTOP2D ! top Exner function
REAL,   DIMENSION(:,:)  , INTENT(OUT)     :: PPSURF    ! Surface pressure
REAL,                     INTENT(OUT) :: PDIAG     ! diagnostics computing time
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSTH_MX ! large scale potential temperature
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSRV_MX ! large scale vapor mixing ratios
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER                                          :: ILBX,ILBY
INTEGER                                          :: IIB,IIE,IIU
INTEGER                                          :: IJB,IJE,IJU
INTEGER                                          :: IKB,IKE,IKU
INTEGER                                          :: JRR
REAL, DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)):: ZTHV
!                                                  ! virtual potential temperature
!                                                  ! on MESONH grid
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))   :: ZHEXNFLUX,ZHEXNMASS,ZPMHP
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))   :: ZRHOD,ZSUMRT
!
INTEGER      :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll => NULL()  ! list of fields to exchange
!
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
TYPE(TFIELDDATA) :: TZFIELD

!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIU=SIZE(PJ,1)
IJU=SIZE(PJ,2)
IKB=JPVEXT+1
IKE=SIZE(PJ,3)-JPVEXT
IKU=SIZE(PJ,3)
!
!-------------------------------------------------------------------------------
!
!*       1.    SHIFT AND INTERPOLATION TO MESONH GRID
!              --------------------------------------
!
ALLOCATE(XTHT(IIU,IJU,IKU))
ALLOCATE(XRT(IIU,IJU,IKU,NRR))

CALL MPPDB_CHECK3D(PTHV_MX,"ver_thermo:PTHV_MX",PRECISION)
CALL MPPDB_CHECK3D(PR_MX(:,:,:,1),"ver_thermo:PR_MX",PRECISION)
CALL MPPDB_CHECK2D(PZS_LS,"ver_thermo:PZS_LS",PRECISION)
CALL MPPDB_CHECK2D(PZSMT_LS,"ver_thermo:PZSMT_LS",PRECISION)
CALL MPPDB_CHECK3D(PZMASS_MX,"ver_thermo:PZMASS_MX",PRECISION)
CALL MPPDB_CHECK3D(PZFLUX_MX,"ver_thermo:PZFLUX_MX",PRECISION)
CALL MPPDB_CHECK3D(PPMHP_MX,"ver_thermo:PPMHP_MX",PRECISION)
CALL MPPDB_CHECK2D(PEXNTOP2D,"ver_thermo:PEXNTOP2D",PRECISION)


IF ( PRESENT(PLSTH_MX)) THEN
  ALLOCATE(XLSTHM(IIU,IJU,IKU))
  ALLOCATE(XLSRVM(IIU,IJU,IKU))
  CALL MPPDB_CHECK3D(PLSTH_MX,"PLSTH_MX",PRECISION)
  CALL MPPDB_CHECK3D(PLSRV_MX,"PLSRV_MX",PRECISION)
  !
  CALL VER_INT_THERMO(TPFILE,OSHIFT,PTHV_MX,PR_MX,PZS_LS,PZSMT_LS,PZMASS_MX,PZFLUX_MX,PPMHP_MX,PEXNTOP2D, &
                      ZTHV,XRT,ZPMHP,PDIAG,PLSTH_MX,PLSRV_MX,XLSTHM,XLSRVM)
ELSE
  CALL VER_INT_THERMO(TPFILE,OSHIFT,PTHV_MX,PR_MX,PZS_LS,PZSMT_LS,PZMASS_MX,PZFLUX_MX,PPMHP_MX,PEXNTOP2D, &
                      ZTHV,XRT,ZPMHP,PDIAG)
END IF
!
XTHT(:,:,:)=ZTHV(:,:,:)*(1.+WATER_SUM(XRT(:,:,:,:)))/(1.+XRV/XRD*XRT(:,:,:,1))
!
!20131113 add update_halo here
CALL ADD3DFIELD_ll(TZFIELDS_ll,XTHT )
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
CALL MPPDB_CHECK3D(XTHT,"PGDFILTER9:XTHT",PRECISION)
!
ZTHV(:,:,1)=ZTHV(:,:,2)
XTHT(:,:,1)=XTHT(:,:,2)
XRT(:,:,1,:)=XRT(:,:,2,:)
!
IF (NRR>=3) THEN
  WHERE  (XRT(:,:,:,3)<1.E-20)
    XRT(:,:,:,3)=0.
  END WHERE
END IF
CALL EXTRAPOL('E',XTHT)

CALL MPPDB_CHECK3D(XTHT,"VERTHERMO::XTHT",PRECISION)

DO JRR=1,SIZE(XRT,4)
  CALL EXTRAPOL('E',XRT(:,:,:,JRR))
END DO
!
IF (NVERB>=10) THEN
  TZFIELD%CMNHNAME   = 'THV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'THV'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_THV'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZTHV)
END IF
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION OF 1D REFERENCE STATE VARIABLES
!              -------------------------------------------
!
CALL MPPDB_CHECK3D(ZTHV,"VERTHERMO bef set_refz::ZTHV",PRECISION)
CALL MPPDB_CHECK3D(XRT(:,:,:,1),"VERTHERMO bef set_refz::XRT",PRECISION)
!
CALL SET_REFZ(ZTHV,XRT(:,:,:,1))

CALL MPPDB_CHECK3D(ZTHV,"VERTHERMO aft set_refz::ZTHV",PRECISION)
CALL MPPDB_CHECK3D(XRT(:,:,:,1),"VERTHERMO: aft set_refz:XRT",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF 3D REFERENCE STATE VARIABLES
!              -------------------------------------------
!
ALLOCATE(XRHODREF(IIU,IJU,IKU))
ALLOCATE(XTHVREF(IIU,IJU,IKU))
ALLOCATE(XRVREF(IIU,IJU,IKU))
ALLOCATE(XEXNREF(IIU,IJU,IKU))
ALLOCATE(XRHODJ(IIU,IJU,IKU))
XRVREF(:,:,:) = 0.
CALL SET_REF(0,TFILE_DUMMY,XZZ,XZHAT,PJ,PDXX,PDYY,CLBCX,CLBCY,       &
             XREFMASS,XMASS_O_PHI0,XLINMASS,XRHODREF,XTHVREF,XRVREF, &
             XEXNREF,XRHODJ)

CALL MPPDB_CHECK3D(XRHODREF,"VERTHERMO::XRHODREF",PRECISION)
CALL MPPDB_CHECK3D(XTHVREF,"VERTHERMO::XTHVREF",PRECISION)
CALL MPPDB_CHECK3D(XRVREF,"VERTHERMO::XRVREF",PRECISION)
CALL MPPDB_CHECK3D(XEXNREF,"VERTHERMO::XEXNREF",PRECISION)
CALL MPPDB_CHECK3D(XRHODJ,"VERTHERMO::XRHODJ",PRECISION)

!
!-------------------------------------------------------------------------------
!
!*       4.    PRESSURE
!              --------
!
CALL COMPUTE_EXNER_FROM_TOP(ZTHV,XZZ,PEXNTOP2D,ZHEXNFLUX,ZHEXNMASS)
!
PPSURF(:,:) = 1.5*ZPMHP(:,:,JPVEXT+1) - 0.5*ZPMHP(:,:,JPVEXT+2) &
             + XP00*ZHEXNFLUX(:,:,JPVEXT+1) ** (XCPD/XRD)
!
ALLOCATE(XPABST(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
XPABST(:,:,:)=ZPMHP(:,:,:) + XP00*ZHEXNMASS(:,:,:) ** (XCPD/XRD)

CALL EXTRAPOL('E',XPABST)
!
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTATION OF TOTAL DRY MASS
!              -----------------------------
!
ZSUMRT(:,:,:) = 0.
DO JRR=1,SIZE(XRT,4)
  ZSUMRT(:,:,:) = ZSUMRT(:,:,:) + XRT(:,:,:,JRR)
END DO
!
ZRHOD(:,:,:)=XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**(XRD/XCPD) &
            /(XRD*ZTHV(:,:,:)*(1.+ZSUMRT(:,:,:)))
!
CALL TOTAL_DMASS(PJ,ZRHOD,XDRYMASST)
!
!-------------------------------------------------------------------------------
!
!*       7.    LARGE SCALE FIELDS INITIALIZATIONS
!              ----------------------------------
!
                                 ! 3D case
!
  CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,               &
                    IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
                    IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  CALL GET_SIZEY_LB(NIMAX_ll,NJMAX_ll,NRIMY,               &
                    IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV, &
                    IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)

IF ( .NOT. PRESENT(PLSTH_MX) ) THEN
  ALLOCATE(XLSTHM(IIU,IJU,IKU))
  ALLOCATE(XLSRVM(IIU,IJU,IKU))
  XLSTHM=XTHT
  XLSRVM=XRT(:,:,:,1)
END IF
! copy at the external levels
XLSTHM(:,:,IKB-1)=XLSTHM(:,:,IKB)
XLSTHM(:,:,IKE+1)=XLSTHM(:,:,IKE)
XLSRVM(:,:,IKB-1)=XLSRVM(:,:,IKB)
XLSRVM(:,:,IKE+1)=XLSRVM(:,:,IKE)
!
CALL EXTRAPOL('E',XLSTHM,XLSRVM)
!
IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2*JPHEXT
    NSIZELBXU_ll=2*NRIMX+2*JPHEXT
    NSIZELBY_ll=2*NRIMY+2*JPHEXT
    NSIZELBYV_ll=2*NRIMY+2*JPHEXT
   ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,IKU))
   ALLOCATE(XLBYTHM(IISIZEYF,IJSIZEYF,IKU))
ELSE
    NSIZELBX_ll=2*JPHEXT     ! 2
    NSIZELBXU_ll=2*(JPHEXT+1) ! 4
    NSIZELBY_ll=2*JPHEXT     ! 2
    NSIZELBYV_ll=2*(JPHEXT+1) ! 4
   ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,IKU))
   ALLOCATE(XLBYTHM(IISIZEY2,IJSIZEY2,IKU))
END IF
!
!ILBX=SIZE(XLBXTHM,1)/2-1
!XLBXTHM(1:ILBX+1,:,:)         = XTHT(IIB-1:IIB-1+ILBX,:,:)
!XLBXTHM(ILBX+2:2*ILBX+2,:,:)  = XTHT(IIE+1-ILBX:IIE+1,:,:)
!ILBY=SIZE(XLBYTHM,2)/2-1
!XLBYTHM(:,1:ILBY+1,:)        = XTHT(:,IJB-1:IJB-1+ILBY,:)
!XLBYTHM(:,ILBY+2:2*ILBY+2,:) = XTHT(:,IJE+1-ILBY:IJE+1,:)
!
IF ( NRR > 0 ) THEN
  IF (       LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
        .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH                     &
     ) THEN
!   ALLOCATE(XLBXRM(2*NRIMX+2*JPHEXT,IJU,IKU,NRR))
!   ALLOCATE(XLBYRM(IIU,2*NRIMY+2*JPHEXT,IKU,NRR))
! ELSE
!   ALLOCATE(XLBXRM(2,IJU,IKU,NRR))
!   ALLOCATE(XLBYRM(IIU,2,IKU,NRR))
      NSIZELBXR_ll=2*NRIMX+2*JPHEXT
      NSIZELBYR_ll=2*NRIMY+2*JPHEXT
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,IKU,NRR))
      ALLOCATE(XLBYRM(IISIZEYF,IJSIZEYF,IKU,NRR))
    ELSE
      NSIZELBXR_ll=2*JPHEXT     !2
      NSIZELBYR_ll=2*JPHEXT     !2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,IKU,NRR))
      ALLOCATE(XLBYRM(IISIZEY2,IJSIZEY2,IKU,NRR))
  ENDIF
  !
  IF (SIZE(XLBXRM) .NE. 0 ) THEN
     ILBX=SIZE(XLBXRM,1)/2-JPHEXT
     XLBXRM(1:ILBX+JPHEXT,:,:,:)                  = XRT(1:ILBX+JPHEXT,:,:,:)
     XLBXRM(ILBX+JPHEXT+1:2*ILBX+2*JPHEXT,:,:,:)  = XRT(IIE+1-ILBX:IIE+JPHEXT,:,:,:)
  ENDIF
  IF (SIZE(XLBYRM) .NE. 0 ) THEN
     ILBY=SIZE(XLBYRM,2)/2-JPHEXT
     XLBYRM(:,1:ILBY+JPHEXT,:,:)        = XRT(:,1:ILBY+JPHEXT,:,:)
     XLBYRM(:,ILBY+JPHEXT+1:2*ILBY+2*JPHEXT,:,:) = XRT(:,IJE+1-ILBY:IJE+JPHEXT,:,:)
  ENDIF
ELSE
   NSIZELBXR_ll=0
   NSIZELBYR_ll=0
   ALLOCATE(XLBXRM(0,0,0,0))
   ALLOCATE(XLBYRM(0,0,0,0))
END IF
!
!NSIZELBXR_ll=SIZE(XLBXRM,1)
!NSIZELBYR_ll=SIZE(XLBYRM,2)   !! coding for one processor

ILBX=SIZE(XLBXTHM,1)
ILBY=SIZE(XLBYTHM,2)
IF(LWEST_ll() .AND. .NOT. L1D) THEN
  XLBXTHM(1:NRIMX+JPHEXT,        :,:)   = XTHT(1:NRIMX+JPHEXT,        :,:)
  XLBXRM(1:NRIMX+JPHEXT,        :,:,:)   = XRT(1:NRIMX+JPHEXT,        :,:,:)
ENDIF
IF(LEAST_ll() .AND. .NOT. L1D) THEN
  XLBXTHM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)   = XTHT(IIU-NRIMX-JPHEXT+1:IIU,    :,:)
  XLBXRM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:,:)   = XRT(IIU-NRIMX-JPHEXT+1:IIU,    :,:,:)
ENDIF
IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYTHM(:,1:NRIMY+JPHEXT,        :)    = XTHT(:,1:NRIMY+JPHEXT,      :)
  XLBYRM(:,1:NRIMY+JPHEXT,        :,:)   = XRT(:,1:NRIMY+JPHEXT,      :,:)
ENDIF
IF(LNORTH_ll().AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYTHM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)    = XTHT(:,IJU-NRIMY-JPHEXT+1:IJU,  :)
  XLBYRM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:,:)   = XRT(:,IJU-NRIMY-JPHEXT+1:IJU,  :,:)
ENDIF
!
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT0%NLU,*) 'Routine VER_THERMO completed'
!
END SUBROUTINE VER_THERMO
