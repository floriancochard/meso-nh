!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_WRITE_LFIFM_n
!     #########################
!
INTERFACE
!
SUBROUTINE WRITE_LFIFM_n(TPFILE,HDADFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
IMPLICIT NONE
!
TYPE(TFILEDATA), INTENT(IN) :: TPFILE   ! File characteristics
CHARACTER(LEN=*),INTENT(IN) :: HDADFILE ! Corresponding FM-file name of its DAD model
END SUBROUTINE WRITE_LFIFM_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LFIFM_n
!
!
!     ##########################################
      SUBROUTINE WRITE_LFIFM_n(TPFILE,HDADFILE)
!     ##########################################
!
!!****  *WRITE_LFIFM_n* - routine to write a LFIFM file for model $n
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write an initial LFIFM File 
!     of name YFMFILE//'.lfi' with the FM routines.  
!
!!**  METHOD
!!    ------
!!      The data are written in the LFIFM file :
!!        - dimensions
!!        - grid variables
!!        - configuration variables
!!        - prognostic variables at time t and t-dt
!!        - 1D anelastic reference state
!!
!!      The localization on the model grid is also indicated :
!!
!!        IGRID = 1 for mass grid point
!!        IGRID = 2 for U grid point
!!        IGRID = 3 for V grid point
!!        IGRID = 4 for w grid point
!!        IGRID = 0 for meaningless case
!!          
!!
!!    EXTERNAL
!!    --------
!!      WRITE_BALLOON_n : routine to write balloon records
!!      WRITE_LB_n : routine to write LB fields
!!      FMWRIT     : FM-routine to write a record
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DIM_n   : contains dimensions
!!      Module MODD_TIME    : contains time variables for all models
!!      Module MODD_TIME_n   : contains time variables 
!!      Module MODD_GRID    : contains spatial grid variables for all models
!!      Module MODD_GRID_n : contains spatial grid variables
!!      Module MODD_REF     : contains reference state variables
!!      Module MODD_LUNIT_n: contains logical unit variables.
!!      Module MODD_CONF    : contains configuration variables for all models
!!      Module MODD_CONF_n  : contains configuration variables
!!      Module MODD_FIELD_n  : contains prognostic variables
!!      Module MODD_GR_FIELD_n : contains surface prognostic variables
!!      Module MODD_LSFIELD_n  : contains Larger Scale variables
!!      Module MODD_PARAM_n    : contains parameterization options
!!      Module MODD_TURB_n    : contains turbulence options
!!      Module MODD_FRC    : contains forcing variables
!!      Module MODD_DEEP_CONVECTION_n : contains deep convection tendencies
!!      Module MODD_PARAM_KAFR_n : contains configuration
!!      Module MODD_AIRCRAFT_BALLOON : contains balloon and aircraft variables
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq   *Meteo France* 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/05/94 
!!       V. Ducrocq    27/06/94                  
!!       J.Stein       20/10/94 (name of the FMFILE)
!!       J.Stein       06/12/94 add the LS fields   
!!       J.P. Lafore   09/01/95 add the DRYMASST
!!       J.Stein       20/01/95 add TKE and change the ycomment for the water 
!!                              variables       
!!       J.Stein       23/01/95 add a TKE switch and MODD_PARAM_n
!!       J.Stein       16/03/95 remove R from the historical variables
!!       J.Stein       20/03/95 add the EPS var.  
!!       J.Stein       30/06/95 add the variables related to the subgrid condens
!!       S. Belair     01/09/95 add surface variables and ground parameters
!!       J.-P. Pinty   15/09/95 add the radiation parameters
!!       J.Stein       23/01/96 add the TSZ0 option for the surface scheme
!!       M.Georgelin   13/12/95 add the forcing variables 
!!       J.-P. Pinty   15/02/96 add external control for the forcing
!!       J.Stein P.Bougeault  15/03/96 add the cloud fraction and change the
!!                                     surface parameters for TSZ0 option
!!       J.Stein P.Jabouille  30/04/96 add the storage type
!!       J.Stein P.Jabouille  20/05/96 switch for XSIGS and XSRC 
!!       J.Stein              10/10/96 change Xsrc into XSRCM and XRCT
!!       J.P. Lafore          30/07/96 add YFMFILE and HDADFILE writing
!!                                     corresponding to MY_NAME and DAD_NAME (for nesting)
!!       V.Masson      08/10/96 add LTHINSHELL
!!       J.-P. Pinty   15/12/96 add the microphysics (ice)
!!       J.-P. Pinty   11/01/97 add the deep convection
!!       J.-P. Pinty   27/01/97 split the recording of the SV array
!!       J.-P. Pinty   29/01/97 set recording of PRCONV and PACCONV in mm/h and
!!                                                         mm respectively
!!       J. Viviand    04/02/97 convert precipitation rates in mm/h
!!       J.P. Lafore   25/11/96 resolution ratio and position for nesting
!!       J.P. Lafore   26/02/97 adding of "surfacic" LS fields
!!       J.Stein       22/06/97 use the absolute pressure
!!       V.Masson      09/07/97 add directional z0 and Subgrid-Scale Orography
!!       V.Masson      18/08/97 call to fmwrit directly with dates and strings
!!       J.Stein       22/10/97 add the LB fields for U,V,W, THETA, RV....
!!       P.Bechtold    24/01/98 add convective tracer tendencies
!!       P.Jabouille   15/10/98 //
!!       P.Jabouille   25/05/99 replace 'DTRAD_CLONLY' by 'DTRAD_CLLY' (size too long)
!!       J. Stein      20/05/98 remove NXEND and NYEND
!!       V. Masson     04/01/00 remove TSZ0 option
!!       P. Jabouille  03/04/00 write XCIT only for MESONH program
!!       K. Suhre      03/12/99 add chemical variable names                         
!        F.solmon /V.Masson   06/00 adapt for patch surface variables
!!       D.Gazen       22/01/01 use MODD_NSV and add names to scalar variables
!!       G.Jaubert     06/06/01 add Balloon current positions
!!       P.Jabouille   10/04/02 extra radiative surface flux
!!       J.-P. Pinty   29/11/02 add C3R5, ICE2, ICE4, CELEC
!!       V. Masson     01/2004  removes surface (externalization)
!!                     05/2006  Remove KEPS
!!       J. escobar    02/09/2009 missing YDIR for CLDFR variable
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!       P. Aumond     12/2009 Mean_UM,...
!!       M. Leriche    16/07/10 add ice phase chemical species
!!       C. Barthe     Jan. 2011  add diagnostics for elec
!!       J. Escobar    Feb. 2012  replace MINVAL/MAXVAL by MIN_ll/MAX_ll in OUTPUT_LISTING
!!       P.Peyrille    06/12 2D west african monsoon: ADV forcing and fluxes writing
!!                     AEROSOLS and ozone vertical distribution are also written
!!       M.Tomasini    06/12 2D west african monsoon: nesting for ADV forcing writing
!!       Pialat/Tulet  15/02/2012 add ForeFire variables
!!       J. Escobar    Mars 2014 , missing YDIR="XY" in 1.6 for tendencies fields 
!!       J.escobar & M.Leriche 23/06/2014 Pb with JSA increment versus ini_nsv order initialization 
!!       P. Tulet      Nov 2014 accumulated moles of aqueous species that fall at the surface
!!       M.Faivre      2014
!!       C.Lac         Dec.2014 writing past wind fields for centred advection
!!       J.-P. Pinty   Jan 2015 add LNOx and flash map diagnostics
!!       J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!       P. Tulet & M. Leriche    Nov 2015 add mean pH value in the rain at the surface
!!       J.escobar     04/08/2015 suit Pb with writ_lfin JSA increment , modif in ini_nsv to have good order initialization
!!       Modification    01/2016  (JP Pinty) Add LIMA
!!       M.Mazoyer     04/16 : Add supersaturation fields
!!       P.Wautelet    11/07/2016 removed MNH_NCWRIT define
!!       JP Chaboureau 27/11/2017 add wind tendency forcing
!!                   02/2018 Q.Libois move Diagnostic related to the radiations in radiations.f90
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!       V. Vionnet    07/2017, add blowing snow variables
!!       P.Wautelet    11/01/2019: bug correction in write XBL_DEPTH->XSBL_DEPTH
!!       C.Lac         18/02/2019: add rain fraction as an output field              
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_TIME
USE MODD_TIME_n
USE MODD_FIELD_n
USE MODD_MEAN_FIELD_n
USE MODD_DUMMY_GR_FIELD_n
USE MODD_LSFIELD_n
USE MODD_DYN_n
USE MODD_PARAM_n
USE MODD_REF
USE MODD_LUNIT_n
USE MODD_TURB_n
USE MODD_RADIATIONS_n,   ONLY : XDTHRAD, NCLEARCOL_TM1, XFLALWD, &
                                XZENITH, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, &
                                XDIRSRFSWD, XSCAFLASWD, XDIRFLASWD, XAZIM
USE MODD_REF_n,  ONLY : XRHODREF
USE MODD_FRC
USE MODD_PRECIP_n
USE MODD_ELEC_n
USE MODD_CST
USE MODD_CLOUDPAR
USE MODD_DEEP_CONVECTION_n
USE MODD_PARAM_KAFR_n
USE MODD_NESTING
USE MODD_PARAMETERS
USE MODD_GR_FIELD_n
USE MODD_CH_MNHC_n,       ONLY: LUSECHEM,LCH_CONV_LINOX, &
                                LUSECHAQ,LUSECHIC,LCH_PH, XCH_PHINIT
USE MODD_CH_PH_n
USE MODD_CH_M9_n
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_ICE_C1R3_DESCR,  ONLY: C1R3NAMES
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES, LLNOX_EXPLICIT
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_NSV
USE MODD_AIRCRAFT_BALLOON
USE MODD_HURR_CONF, ONLY: LFILTERING,CFILTERING,NDIAG_FILT
USE MODD_HURR_FIELD_n
USE MODD_PREP_REAL, ONLY: CDUMMY_2D, XDUMMY_2D
USE MODD_DUST
USE MODD_SALT
USE MODD_PASPOL
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE
#endif
USE MODD_CONDSAMP
USE MODD_CH_AEROSOL
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_PAST_FIELD_n
USE MODD_ADV_n, ONLY: CUVW_ADV_SCHEME,XRTKEMS,CTEMP_SCHEME,LSPLIT_CFL
USE MODD_ELEC_FLASH
!
USE MODD_PARAM_LIMA     , ONLY: NMOD_CCN, LSCAV, LAERO_MASS,                &
                                NMOD_IFN, NMOD_IMM, NINDICE_CCN_IMM, LHHONI
USE MODD_PARAM_LIMA_WARM, ONLY: CLIMA_WARM_NAMES, CAERO_MASS
USE MODD_PARAM_LIMA_COLD, ONLY: CLIMA_COLD_NAMES
USE MODD_LIMA_PRECIP_SCAVENGING_n
!
USE MODE_FM,    ONLY: IO_FILE_CLOSE_ll
USE MODE_FMWRIT
USE MODE_ll
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_IO_ll, ONLY: UPCASE,CLOSE_ll
USE MODE_FIELD
USE MODE_GATHER_ll
USE MODE_GRIDPROJ
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
USE MODI_WRITE_LB_n
USE MODI_WRITE_BALLOON_n
USE MODI_DUSTLFI_n
USE MODI_SALTLFI_n
USE MODI_CH_AER_REALLFI_n
!
!20131128
USE MODE_MPPDB
USE MODE_EXTRAPOL
! Modif Eddy fluxes
USE MODD_DEF_EDDY_FLUX_n       ! Ajout PP
USE MODD_DEF_EDDYUV_FLUX_n     ! Ajout PP
USE MODD_LATZ_EDFLX            ! Ajout PP
!
USE MODD_2D_FRC                  ! Ajout PP
USE MODD_ADVFRC_n              ! Modif PP ADV FRC
USE MODD_RELFRC_n
!
USE MODD_PARAM_C2R2
! 
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA), INTENT(IN) :: TPFILE   ! File characteristics
CHARACTER(LEN=*),INTENT(IN) :: HDADFILE ! Corresponding FM-file name of its DAD model
!
!*       0.2   Declarations of local variables
!
INTEGER           :: ILUOUT         ! logical unit
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears 
                                    !in LFI subroutines at the open of the file              
!
INTEGER           :: JSV            ! loop index for scalar variables
INTEGER           :: JSA            ! beginning of chemical-aerosol variables

! 
CHARACTER(LEN=3)  :: YFRC           ! to mark the time of the forcing
INTEGER           :: JT             ! loop index
!
INTEGER           :: JMOM, IMOMENTS, JMODE, ISV_NAME_IDX  ! dust modes
! 
REAL,DIMENSION(:,:), ALLOCATABLE  :: ZWORK2D     ! Working array
REAL,DIMENSION(:,:,:), ALLOCATABLE  :: ZWORK3D     ! Working array
!
REAL                              :: ZLATOR, ZLONOR ! geographical coordinates of 1st mass point
REAL                              :: ZXHATM, ZYHATM ! conformal    coordinates of 1st mass point
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
INTEGER :: IMI ! Current model index
!
INTEGER           :: ICH_NBR        ! to write number and names of scalar 
INTEGER,DIMENSION(:),ALLOCATABLE :: ICH_NAMES !(chem+aero+dust) variables
CHARACTER(LEN=NMNHNAMELGTMAX),DIMENSION(:),ALLOCATABLE :: YDSTNAMES,YCHNAMES, YSLTNAMES
INTEGER           :: ILREC,ILENG    !in NSV.DIM and NSV.TITRE
INTEGER           :: INFO_ll
INTEGER :: IKRAD
INTEGER           :: JI,JJ,JK   ! loop index
INTEGER           :: IIU,IJU,IKU,IIB,IJB,IKB,IIE,IJE,IKE ! Arrays bounds
!
CHARACTER(LEN=2)  :: INDICE
INTEGER           :: I, IID
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
!
!*	0. Initialization
!
IMI = GET_CURRENT_MODEL_INDEX()
!
ILUOUT=TLUOUT%NLU
!
ALLOCATE(ZWORK2D(SIZE(XTHT,1),SIZE(XTHT,2)))
ALLOCATE(ZWORK3D(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)))
!
!*       0.2     ARRAYS BOUNDS INITIALIZATION
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
!*       1.     WRITES IN THE LFI FILE
! 
!
!*       1.0    File and HDADFILE writing :
!
CALL IO_WRITE_FIELD(TPFILE,'FILETYPE',TPFILE%CTYPE)
!
IF (LEN_TRIM(HDADFILE)>0) THEN
  CALL IO_WRITE_FIELD(TPFILE,'DXRATIO',NDXRATIO_ALL(IMI))
  CALL IO_WRITE_FIELD(TPFILE,'DYRATIO',NDYRATIO_ALL(IMI))
  CALL IO_WRITE_FIELD(TPFILE,'XOR',    NXOR_ALL(IMI))
  CALL IO_WRITE_FIELD(TPFILE,'YOR',    NYOR_ALL(IMI))
END IF
!
!*       1.1    Type and Dimensions :
!
CALL IO_WRITE_FIELD(TPFILE,'IMAX',NIMAX_ll)
CALL IO_WRITE_FIELD(TPFILE,'JMAX',NJMAX_ll)
CALL IO_WRITE_FIELD(TPFILE,'KMAX',NKMAX)
!
CALL IO_WRITE_FIELD(TPFILE,'JPHEXT',JPHEXT)
!
!*       1.2    Grid variables :
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_WRITE_FIELD(TPFILE,'RPK',   XRPK)
  CALL IO_WRITE_FIELD(TPFILE,'LONORI',XLONORI)
  CALL IO_WRITE_FIELD(TPFILE,'LATORI',XLATORI)
! 
!* diagnostic of 1st mass point
!
  ALLOCATE(ZXHAT_ll(NIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(NJMAX_ll+2 * JPHEXT))
  CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP) !//
  CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP) !//
  ZXHATM = 0.5 * (ZXHAT_ll(1)+ZXHAT_ll(2))
  ZYHATM = 0.5 * (ZYHAT_ll(1)+ZYHAT_ll(2))
  CALL SM_LATLON(XLATORI,XLONORI,ZXHATM,ZYHATM,ZLATOR,ZLONOR)
  DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
!
  CALL IO_WRITE_FIELD(TPFILE,'LONOR',ZLONOR)
  CALL IO_WRITE_FIELD(TPFILE,'LATOR',ZLATOR)
END IF 
!
CALL IO_WRITE_FIELD(TPFILE,'THINSHELL',LTHINSHELL)
CALL IO_WRITE_FIELD(TPFILE,'LAT0',XLAT0)
CALL IO_WRITE_FIELD(TPFILE,'LON0',XLON0)
CALL IO_WRITE_FIELD(TPFILE,'BETA',XBETA)
!
CALL IO_WRITE_FIELD(TPFILE,'XHAT',XXHAT)
CALL IO_WRITE_FIELD(TPFILE,'YHAT',XYHAT)
CALL IO_WRITE_FIELD(TPFILE,'ZHAT',XZHAT)
CALL IO_WRITE_FIELD(TPFILE,'ZTOP',XZTOP)
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_WRITE_FIELD(TPFILE,'LAT',XLAT)
  CALL IO_WRITE_FIELD(TPFILE,'LON',XLON)
END IF
!
CALL IO_WRITE_FIELD(TPFILE,'ZS',   XZS)
IF(ASSOCIATED(XZWS)) THEN
  CALL IO_WRITE_FIELD(TPFILE,'ZWS',  XZWS)
END IF
CALL IO_WRITE_FIELD(TPFILE,'ZSMT', XZSMT)
CALL IO_WRITE_FIELD(TPFILE,'SLEVE',LSLEVE)
!
IF (LSLEVE) THEN
  CALL IO_WRITE_FIELD(TPFILE,'LEN1',XLEN1)
  CALL IO_WRITE_FIELD(TPFILE,'LEN2',XLEN2)
END IF
!
!
CALL IO_WRITE_FIELD(TPFILE,'DTMOD',TDTMOD)
CALL IO_WRITE_FIELD(TPFILE,'DTCUR',TDTCUR)
CALL IO_WRITE_FIELD(TPFILE,'DTEXP',TDTEXP)
CALL IO_WRITE_FIELD(TPFILE,'DTSEG',TDTSEG)
!
!*       1.3    Configuration  variables :
!
CALL IO_WRITE_FIELD(TPFILE,'L1D',      L1D)
CALL IO_WRITE_FIELD(TPFILE,'L2D',      L2D)
CALL IO_WRITE_FIELD(TPFILE,'PACK',     LPACK)
CALL IO_WRITE_FIELD(TPFILE,'CARTESIAN',LCARTESIAN)
CALL IO_WRITE_FIELD(TPFILE,'LBOUSS',   LBOUSS)
!
CALL IO_WRITE_FIELD(TPFILE,'SURF',     CSURF)
CALL IO_WRITE_FIELD(TPFILE,'CPL_AROME',LCPL_AROME)
CALL IO_WRITE_FIELD(TPFILE,'COUPLING', LCOUPLING)
!
!*       1.4    Prognostic variables :
!
!
!*       1.4.1  Time t:
!
!20131128 check XUT-> X_Y_W_U wind component for PRC
!  CALL EXTRAPOL('W',XUT)
!  CALL EXTRAPOL('E',XUT)
!  CALL EXTRAPOL('N',XUT)
!  CALL EXTRAPOL('S',XUT)
CALL MPPDB_CHECK3D(XUT,"write_lfifmn before IO_WRITE_FIELD::XUT",PRECISION)
CALL IO_WRITE_FIELD(TPFILE,'UT',XUT)
CALL MPPDB_CHECK3D(XUT,"write_lfifmn after IO_WRITE_FIELD::XUT",PRECISION)
!
!20131128 check XVT-> X_Y_W_V wind component for PRC
CALL MPPDB_CHECK3D(XVT,"write_lfifmn::XVT",PRECISION)
!
CALL IO_WRITE_FIELD(TPFILE,'VT',XVT)
CALL IO_WRITE_FIELD(TPFILE,'WT',XWT)
!
CALL IO_WRITE_FIELD(TPFILE,'THT',XTHT)
!
!*       1.4.2  Time t-dt:
!
IF ( (CUVW_ADV_SCHEME == 'CEN4TH') .AND. (CTEMP_SCHEME == 'LEFR') ) THEN
  CALL IO_WRITE_FIELD(TPFILE,'UM', XUM)
  CALL IO_WRITE_FIELD(TPFILE,'VM', XVM)
  CALL IO_WRITE_FIELD(TPFILE,'WM', XWM)
  CALL IO_WRITE_FIELD(TPFILE,'DUM',XDUM)
  CALL IO_WRITE_FIELD(TPFILE,'DVM',XDVM)
  CALL IO_WRITE_FIELD(TPFILE,'DWM',XDWM)
END IF
!
IF (MEAN_COUNT /= 0) THEN
!
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
!
  TZFIELD%NGRID      = 2
!
  TZFIELD%CMNHNAME   = 'UMME'
  TZFIELD%CLONGNAME  = 'UMME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U component of mean wind'
  ZWORK3D = XUM_MEAN/MEAN_COUNT
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'U2ME'
  TZFIELD%CLONGNAME  = 'U2ME'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U component of mean wind variance'
  ZWORK3D = XU2_MEAN/MEAN_COUNT-XUM_MEAN**2/MEAN_COUNT**2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'UMMA'
  TZFIELD%CLONGNAME  = 'UMMA'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U component of max wind'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XUM_MAX)
!
  TZFIELD%NGRID      = 3
!
  TZFIELD%CMNHNAME   = 'VMME'
  TZFIELD%CLONGNAME  = 'VMME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V component of mean wind'
  ZWORK3D = XVM_MEAN/MEAN_COUNT
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'V2ME'
  TZFIELD%CLONGNAME  = 'V2ME'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V component of mean wind variance'
  ZWORK3D = XV2_MEAN/MEAN_COUNT-XVM_MEAN**2/MEAN_COUNT**2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'VMMA'
  TZFIELD%CLONGNAME  = 'VMMA'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V component of max wind'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XVM_MAX)
!
  TZFIELD%NGRID      = 4
!
  TZFIELD%CMNHNAME   = 'WMME'
  TZFIELD%CLONGNAME  = 'WMME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical mean wind'
  ZWORK3D = XWM_MEAN/MEAN_COUNT
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'W2ME'
  TZFIELD%CLONGNAME  = 'W2ME'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical mean wind variance'
  ZWORK3D = XW2_MEAN/MEAN_COUNT-XWM_MEAN**2/MEAN_COUNT**2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'WMMA'
  TZFIELD%CLONGNAME  = 'WMMA'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical max wind'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XWM_MAX)
!
  TZFIELD%NGRID      = 1
!
  TZFIELD%CMNHNAME   = 'THMME'
  TZFIELD%CLONGNAME  = 'THMME'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean potential temperature'
  ZWORK3D = XTHM_MEAN/MEAN_COUNT
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'TH2ME'
  TZFIELD%CLONGNAME  = 'TH2ME'
  TZFIELD%CUNITS     = 'K2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean potential temperature variance'
  ZWORK3D = XTH2_MEAN/MEAN_COUNT-XTHM_MEAN**2/MEAN_COUNT**2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'THMMA'
  TZFIELD%CLONGNAME  = 'THMMA'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_max potential temperature'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTHM_MAX)
!
  TZFIELD%CMNHNAME   = 'TEMPMME'
  TZFIELD%CLONGNAME  = 'TEMPMME'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean temperature'
  ZWORK3D = XTEMPM_MEAN/MEAN_COUNT
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'TEMP2ME'
  TZFIELD%CLONGNAME  = 'TEMP2ME'
  TZFIELD%CUNITS     = 'K2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean temperature variance'
  ZWORK3D = XTEMP2_MEAN/MEAN_COUNT-XTEMPM_MEAN**2/MEAN_COUNT**2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'TEMPMMA'
  TZFIELD%CLONGNAME  = 'TEMPMMA'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_max temperature'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTEMPM_MAX)
!
  TZFIELD%CMNHNAME   = 'PABSMME'
  TZFIELD%CLONGNAME  = 'PABSMME'
  TZFIELD%CUNITS     = 'Pa'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean ABSolute Pressure'
  ZWORK3D = XPABSM_MEAN/MEAN_COUNT
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'PABS2ME'
  TZFIELD%CLONGNAME  = 'PABS2ME'
  TZFIELD%CUNITS     = 'Pa2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean ABSolute Pressure variance'
  ZWORK3D = XPABS2_MEAN/MEAN_COUNT-XPABSM_MEAN**2/MEAN_COUNT**2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'PABSMMA'
  TZFIELD%CLONGNAME  = 'PABSMMA'
  TZFIELD%CUNITS     = 'Pa'
  TZFIELD%CCOMMENT   = 'X_Y_Z_max ABSolute Pressure'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XPABSM_MAX)
!
  IF (CTURB /= 'NONE') THEN
    TZFIELD%CMNHNAME   = 'TKEMME'
    TZFIELD%CLONGNAME  = 'TKEMME'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_mean kinetic energy'
    ZWORK3D= XTKEM_MEAN/MEAN_COUNT
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
!
    TZFIELD%CMNHNAME   = 'TKEMMA'
    TZFIELD%CLONGNAME  = 'TKEMMA'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_max kinetic energy'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTKEM_MAX)
  END IF
!
END IF
!
!
IF (CTURB /= 'NONE') THEN
  CALL IO_WRITE_FIELD(TPFILE,'TKET',XTKET)
  IF (CPROGRAM == 'MESONH' .AND. LSPLIT_CFL) CALL IO_WRITE_FIELD(TPFILE,'TKEMS',XRTKEMS)
END IF
!
!
CALL IO_WRITE_FIELD(TPFILE,'PABST',XPABST)
!
IF (NRR >=1) THEN
  IF (LUSERV) CALL IO_WRITE_FIELD(TPFILE,'RVT',XRT(:,:,:,IDX_RVT))
  IF (LUSERC) THEN
    CALL IO_WRITE_FIELD(TPFILE,'RCT',XRT(:,:,:,IDX_RCT))
    WRITE (ILUOUT,*) IDX_RCT,' RC min-max ',MIN_ll(XRT(:,:,:,IDX_RCT),INFO_ll),MAX_ll(XRT(:,:,:,IDX_RCT),INFO_ll)
  END IF
  IF (LUSERR) THEN
    CALL IO_WRITE_FIELD(TPFILE,'RRT',XRT(:,:,:,IDX_RRT))
    WRITE (ILUOUT,*) IDX_RRT,' RR min-max ',MIN_ll(XRT(:,:,:,IDX_RRT),INFO_ll),MAX_ll(XRT(:,:,:,IDX_RRT),INFO_ll)
  END IF 
  IF (LUSERI) THEN
    CALL IO_WRITE_FIELD(TPFILE,'RIT',XRT(:,:,:,IDX_RIT))
    WRITE (ILUOUT,*) IDX_RIT,' RI min-max ',MIN_ll(XRT(:,:,:,IDX_RIT),INFO_ll),MAX_ll(XRT(:,:,:,IDX_RIT),INFO_ll)
    IF ( CPROGRAM == 'MESONH' .AND. CCLOUD(1:3) == 'ICE') THEN
      CALL IO_WRITE_FIELD(TPFILE,'CIT',XCIT(:,:,:))
    END IF
  END IF 
  IF (LUSERS) THEN
    CALL IO_WRITE_FIELD(TPFILE,'RST',XRT(:,:,:,IDX_RST))
    WRITE (ILUOUT,*) IDX_RST,' RS min-max ',MINVAL(XRT(:,:,:,IDX_RST)),MAXVAL(XRT(:,:,:,IDX_RST))
  END IF
  IF (LUSERG) THEN
    CALL IO_WRITE_FIELD(TPFILE,'RGT',XRT(:,:,:,IDX_RGT))
    WRITE (ILUOUT,*) IDX_RGT,' RG min-max ',MINVAL(XRT(:,:,:,IDX_RGT)),MAXVAL(XRT(:,:,:,IDX_RGT))
  END IF 
  IF (LUSERH) CALL IO_WRITE_FIELD(TPFILE,'RHT',XRT(:,:,:,IDX_RHT))
END IF
!
IF (NSV >=1) THEN
  JSA=0
  ! User scalar variables
  IF (NSV_USER>0) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = 1,NSV_USER
      WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'SVT',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  END IF
  ! microphysical C2R2 scheme scalar variables
  IF (NSV_C2R2END>=NSV_C2R2BEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'm-3'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_C2R2BEG,NSV_C2R2END
      TZFIELD%CMNHNAME   = TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  END IF
  ! microphysical C3R5 scheme additional scalar variables
  IF (NSV_C1R3END>=NSV_C1R3BEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'm-3'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_C1R3BEG,NSV_C1R3END
      TZFIELD%CMNHNAME   = TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  END IF
!
! microphysical LIMA variables
!
  IF (NSV_LIMA_END>=NSV_LIMA_BEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
  END IF
  !
  DO JSV = NSV_LIMA_BEG,NSV_LIMA_END
    !
    TZFIELD%CUNITS     = 'kg-1'
    WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
    !
! Nc
    IF (JSV .EQ. NSV_LIMA_NC) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_NAMES(1))//'T'
    END IF
! Nr
    IF (JSV .EQ. NSV_LIMA_NR) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_NAMES(2))//'T'
    END IF
! N CCN free
    IF (JSV .GE. NSV_LIMA_CCN_FREE .AND. JSV .LT. NSV_LIMA_CCN_ACTI) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_FREE + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_NAMES(3))//INDICE//'T'
    END IF
! N CCN acti
    IF (JSV .GE. NSV_LIMA_CCN_ACTI .AND. JSV .LT. NSV_LIMA_CCN_ACTI + NMOD_CCN) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_ACTI + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_NAMES(4))//INDICE//'T'
    END IF
! Scavenging
    IF (JSV .EQ. NSV_LIMA_SCAVMASS) THEN
      TZFIELD%CMNHNAME   = TRIM(CAERO_MASS(1))//'T'
      TZFIELD%CUNITS     = 'kg kg-1'
    END IF
! Ni
    IF (JSV .EQ. NSV_LIMA_NI) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_NAMES(1))//'T'
    END IF
! N IFN free
    IF (JSV .GE. NSV_LIMA_IFN_FREE .AND. JSV .LT. NSV_LIMA_IFN_NUCL) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_FREE + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_NAMES(2))//INDICE//'T'
    END IF
! N IFN nucl
    IF (JSV .GE. NSV_LIMA_IFN_NUCL .AND. JSV .LT. NSV_LIMA_IFN_NUCL + NMOD_IFN) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_NUCL + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_NAMES(3))//INDICE//'T'
    END IF
! N IMM nucl
    I = 0
    IF (JSV .GE. NSV_LIMA_IMM_NUCL .AND. JSV .LT. NSV_LIMA_IMM_NUCL + NMOD_IMM) THEN
      I = I + 1
      WRITE(INDICE,'(I2.2)')(NINDICE_CCN_IMM(I))
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_NAMES(4))//INDICE//'T'
    END IF
! Hom. freez. of CCN
    IF (JSV .EQ. NSV_LIMA_HOM_HAZE) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_NAMES(5))//'T'
    END IF
    !
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
!
    JSA=JSA+1
  END DO
!
  IF (LSCAV .AND. LAERO_MASS) THEN
  IF (ASSOCIATED(XINPAP)) THEN
  IF (SIZE(XINPAP) /= 0 ) THEN
    CALL IO_WRITE_FIELD(TPFILE,'INPAP',XINPAP)
    !
    ZWORK2D(:,:)  = XRHOLW*XINPRR(:,:)*XSVT(:,:,2,NSV_LIMA_SCAVMASS)/ &
                                        max( 1.e-20,XRT(:,:,2,3) ) !~2=at ground level
    TZFIELD%CMNHNAME   = 'INPBP'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'INPBP'
    TZFIELD%CUNITS     = 'kg m-2 s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_INstantaneous Precipitating Aerosol Rate'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2D)
    !
    CALL IO_WRITE_FIELD(TPFILE,'ACPAP',XACPAP)
  END IF
  END IF
  END IF
!
!
  ! electrical scalar variables
  IF (NSV_ELECEND>=NSV_ELECBEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_ELECBEG,NSV_ELECEND
      TZFIELD%CMNHNAME   = TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      IF (JSV .GT. NSV_ELECBEG .AND. JSV .LT. NSV_ELECEND) THEN 
        TZFIELD%CUNITS     = 'C m-3'
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      ELSE
        TZFIELD%CUNITS     = 'm-3'
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3,A8)')'X_Y_Z_','SVT',JSV,' (nb ions/m3)'
      END IF
      ZWORK3D(:,:,:) = XSVT(:,:,:,JSV) * XRHODREF(:,:,:) ! C/kg --> C/m3
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
      JSA=JSA+1
    END DO
  END IF
  !
  IF (CELEC /= 'NONE') THEN
    CALL IO_WRITE_FIELD(TPFILE,'EFIELDU',XEFIELDU)
    CALL IO_WRITE_FIELD(TPFILE,'EFIELDV',XEFIELDV)
    CALL IO_WRITE_FIELD(TPFILE,'EFIELDW',XEFIELDW)
 !
    TZFIELD%CMNHNAME   = 'EMODULE'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'V m-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    ZWORK3D(:,:,:) = (XEFIELDU**2 + XEFIELDV**2 + XEFIELDW**2)**0.5
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK3D)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('NI_IAGGS',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XNI_IAGGS*1.E12)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('NI_IDRYG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XNI_IDRYG*1.E12)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('NI_SDRYG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XNI_SDRYG*1.E12)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('INDUC_CG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XIND_RATE*1.E12)
 !
    CALL IO_WRITE_FIELD(TPFILE,'TRIG_IC',   NMAP_TRIG_IC)
    CALL IO_WRITE_FIELD(TPFILE,'IMPACT_CG', NMAP_IMPACT_CG)
    CALL IO_WRITE_FIELD(TPFILE,'AREA_CG',   NMAP_2DAREA_CG)
    CALL IO_WRITE_FIELD(TPFILE,'AREA_IC',   NMAP_2DAREA_IC)
    CALL IO_WRITE_FIELD(TPFILE,'FLASH_3DCG',NMAP_3DCG)
    CALL IO_WRITE_FIELD(TPFILE,'FLASH_3DIC',NMAP_3DIC)
 !
    IF (LLNOX_EXPLICIT) THEN
      TZFIELD%CMNHNAME   = 'LINOX'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'mol mol-1'
      TZFIELD%CDIR       = 'XY'
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,NSV_LNOXEND))
      JSA=JSA+1
    END IF
  END IF
  ! lagrangian variables
  IF (NSV_LGEND>=NSV_LGBEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_LGBEG,NSV_LGEND
      TZFIELD%CMNHNAME   = TRIM(CLGNAMES(JSV-NSV_LGBEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  END IF
  ! Passive scalar variables        
  IF (LPASPOL) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_PPBEG,NSV_PPEND
      WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'SVT',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  END IF
!
  IF ( ((CCLOUD == 'KHKO') .OR.(CCLOUD == 'C2R2')) .AND. (.NOT. LSUPSAT)) THEN
    CALL IO_WRITE_FIELD(TPFILE,'SUPSATMAX',XSUPSAT(:,:,:))
    CALL IO_WRITE_FIELD(TPFILE,'NACT',     XNACT(:,:,:))
  END IF
  IF ( ((CCLOUD == 'KHKO') .OR.(CCLOUD == 'C2R2')) .AND. LSUPSAT) THEN
    CALL IO_WRITE_FIELD(TPFILE,'SSPRO',XSSPRO(:,:,:))
    CALL IO_WRITE_FIELD(TPFILE,'NPRO', XNPRO(:,:,:))
  END IF
!
#ifdef MNH_FOREFIRE
  ! ForeFire scalar variables
  IF ( LFOREFIRE ) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_FFBEG,NSV_FFEND
      WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'SVT',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  END IF
#endif
! Blowing snow variables
  IF (LBLOWSNOW) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    DO JSV = NSV_SNWBEG,NSV_SNWEND
      TZFIELD%CMNHNAME=TRIM(CSNOWNAMES(JSV-NSV_SNWBEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    DO JSV = 1,(NSV_SNW)
      WRITE(TZFIELD%CMNHNAME,'(A10,I3.3)')'SNOWCANO_M',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A8,I3.3)')'X_Y_Z_','SNOWCANO',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSNWCANO(:,:,JSV))
      JSA=JSA+1
    END DO
  ENDIF
  ! Conditional sampling variables  
  IF (LCONDSAMP) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_CSBEG,NSV_CSEND
      WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'SVT',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      JSA=JSA+1
    END DO
  !
  END IF
  ! number of chemical variables (chem+aero+dust)
  ICH_NBR = 0
  IF (LUSECHEM) ICH_NBR = ICH_NBR +NSV_CHEMEND-NSV_CHEMBEG+1 
  IF (LUSECHIC) ICH_NBR = ICH_NBR +NSV_CHICEND-NSV_CHICBEG+1
  IF (.NOT.LUSECHEM.AND.LCH_CONV_LINOX) ICH_NBR = ICH_NBR + &
                                                  NSV_LNOXEND-NSV_LNOXBEG+1 
  IF (LORILAM)  ICH_NBR = ICH_NBR +NSV_AEREND -NSV_AERBEG+1 
  IF (LDUST)    ICH_NBR = ICH_NBR +NSV_DSTEND -NSV_DSTBEG+1
  IF (LDEPOS_DST(IMI))  ICH_NBR = ICH_NBR +NSV_DSTDEPEND -NSV_DSTDEPBEG+1 
  IF (LDEPOS_SLT(IMI))  ICH_NBR = ICH_NBR +NSV_SLTDEPEND -NSV_SLTDEPBEG+1 
  IF (LDEPOS_AER(IMI))  ICH_NBR = ICH_NBR +NSV_AERDEPEND -NSV_AERDEPBEG+1 
  IF (LSALT)    ICH_NBR = ICH_NBR +NSV_SLTEND -NSV_SLTBEG+1
  IF (ICH_NBR /=0) ALLOCATE(YCHNAMES(ICH_NBR))
  ! chemical scalar variables
  IF (LUSECHEM) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_CHEMBEG,NSV_CHEMEND
      TZFIELD%CMNHNAME   = TRIM(UPCASE(CNAMES(JSV-NSV_CHEMBEG+1)))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ppp'
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      !
      YCHNAMES(JSV-JSA)=TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1) ! without T
    END DO
    !
    IF (LUSECHIC) THEN
      DO JSV = NSV_CHICBEG,NSV_CHICEND
        TZFIELD%CMNHNAME   = TRIM(UPCASE(CICNAMES(JSV-NSV_CHICBEG+1)))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CUNITS     = 'ppp'
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        !
        YCHNAMES(JSV-JSA)=TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1) ! without M
      END DO
    ENDIF
    IF (LUSECHAQ.AND.NRR>=3) THEN ! accumulated moles of aqueous species that fall at the surface (mol i/m2) 
      TZFIELD%NDIMS = 2
      DO JSV = NSV_CHACBEG+NSV_CHAC/2,NSV_CHACEND
        TZFIELD%CMNHNAME   = 'ACPR_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHEMBEG+1)))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CUNITS     = 'mol i m-2'
        TZFIELD%CCOMMENT   = 'X_Y_Accumulated moles of aqueous species at the surface'
        ZWORK2D(:,:)  = XACPRAQ(:,:,JSV-NSV_CHACBEG-NSV_CHAC/2+1)
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2D)
      END DO
      TZFIELD%NDIMS = 3
    END IF
    IF (LUSECHAQ.AND.LCH_PH) THEN  ! pH values in cloud
      CALL IO_WRITE_FIELD(TPFILE,'PHC',XPHC)
      IF (NRR>=3) THEN
        CALL IO_WRITE_FIELD(TPFILE,'PHR',XPHR)
        ! compute mean pH in accumulated surface water
        !ZWORK2D(:,:) = 10**(-XCH_PHINIT)
        WHERE (XACPRR > 0.)
          ZWORK2D(:,:) =  XACPHR(:,:) *1E3 / XACPRR(:,:) ! moles of H+ / l of water 
        ELSE WHERE
          ZWORK2D(:,:) = XUNDEF
        END WHERE
        WHERE ((ZWORK2D(:,:) < 1E-1).AND.(ZWORK2D(:,:) > 1E-14))
          ZWORK2D(:,:) = -LOG10(ZWORK2D(:,:))           ! mean pH of surface water
        END WHERE
        TZFIELD%CMNHNAME   = 'MEANPHR'
        TZFIELD%CSTDNAME   = ''
        TZFIELD%CLONGNAME  = 'MEANPHR'
        TZFIELD%CUNITS     = '1'
        TZFIELD%CDIR       = 'XY'
        TZFIELD%CCOMMENT   = 'X_Y_MEAN_PH'
        TZFIELD%NGRID      = 1
        TZFIELD%NTYPE      = TYPEREAL
        TZFIELD%NDIMS      = 2
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2D)
      ENDIF
    ENDIF
  ELSE IF (LCH_CONV_LINOX) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'ppp'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_LNOXBEG,NSV_LNOXEND
      TZFIELD%CMNHNAME   = 'LINOXT'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)') 'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
      YCHNAMES(JSV-JSA)=TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
    END DO
  ENDIF  
  ! aerosol scalar variables
  IF (LORILAM) THEN
    IF ((CPROGRAM == 'REAL  ').AND.(NSV_AER > 1).AND.(IMI==1).AND.(LAERINIT))  &
      CALL CH_AER_REALLFI_n(XSVT(:,:,:,NSV_AERBEG:NSV_AEREND),XSVT(:,:,:,NSV_CHEMBEG-1+JP_CH_CO), XRHODREF)
    IF ((CPROGRAM == 'IDEAL ').AND.(NSV_AER > 1).AND.(IMI==1))  &
      CALL CH_AER_REALLFI_n(XSVT(:,:,:,NSV_AERBEG:NSV_AEREND),XSVT(:,:,:,NSV_CHEMBEG-1+JP_CH_CO),  XRHODREF)
    IF (NSV_AEREND>=NSV_AERBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_AERBEG,NSV_AEREND
        TZFIELD%CMNHNAME   = TRIM(UPCASE(CAERONAMES(JSV-NSV_AERBEG+1)))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        IF (JSV==NSV_AERBEG) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_AERBEG ',JSV
        IF (JSV==NSV_AEREND) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_AEREND ',JSV
        YCHNAMES(JSV-JSA)=  TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
      END DO
    END IF
    IF (LDEPOS_AER(IMI)) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_AERDEPBEG,NSV_AERDEPEND
        TZFIELD%CMNHNAME   = TRIM(CDEAERNAMES(JSV-NSV_AERDEPBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        IF (JSV==NSV_AERDEPBEG) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_AERDEPBEG ',JSV
        IF (JSV==NSV_AERDEPEND) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_AERDEPEND ',JSV
        YCHNAMES(JSV-JSA) = TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
      END DO   ! Loop on aq dust scalar variables      
    ENDIF
  END IF
  ! dust scalar variables
  IF (LDUST) THEN
    IF ((CPROGRAM == 'REAL  ').AND.(NSV_DST > 1).AND.(IMI==1).AND.(LDSTINIT)) &
      CALL DUSTLFI_n(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), XRHODREF)
    IF ((CPROGRAM == 'IDEAL ').AND.(NSV_DST > 1).AND.(IMI==1)) &
      CALL DUSTLFI_n(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), XRHODREF)
    !At this point, we have the tracer array in order of importance, i.e.
    !if mode 2 is most important it will occupy place 1-3 of XSVT  
    IF ((CPROGRAM == 'REAL  ').AND.((LDSTINIT).OR.(LDSTPRES)).OR.&
       (CPROGRAM == 'IDEAL ')               ) THEN
      ! In this case CDUSTNAMES is not allocated. We will use YPDUST_INI,
      !but remember that this variable does not follow JPDUSTORDER
      IMOMENTS = INT(NSV_DSTEND - NSV_DSTBEG+1)/NMODE_DST  
      !Should equal 3 at this point
      IF (IMOMENTS > 3) THEN
        WRITE(ILUOUT,*) 'Error in write_lfin: number of moments must equal or inferior to 3'
        WRITE(ILUOUT,*) NSV_DSTBEG, NSV_DSTEND,NMODE_DST,IMOMENTS
 !callabortstop
        CALL IO_FILE_CLOSE_ll(TLUOUT)
        CALL ABORT
        STOP
      END IF ! Test IMOMENTS
      ALLOCATE(YDSTNAMES(NSV_DSTEND - NSV_DSTBEG+1))
      !
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      IF (IMOMENTS == 1) THEN
        DO JMODE=1, NMODE_DST
          ISV_NAME_IDX = (JPDUSTORDER(JMODE) - 1)*3 + 2
          JSV = (JMODE-1)*IMOMENTS  & !Number of moments previously counted
                  +  1              & !Number of moments in this mode
                  + (NSV_DSTBEG -1)      !Previous list of tracers
          TZFIELD%CMNHNAME   = TRIM(YPDUST_INI(ISV_NAME_IDX))//'T'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
          CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
          YDSTNAMES((JMODE-1)*IMOMENTS+1)=TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
        END DO ! Loop on mode
      ELSE
        DO JMODE=1,NMODE_DST
          DO JMOM=1,IMOMENTS
            ISV_NAME_IDX = (JPDUSTORDER(JMODE) - 1)*3 + JMOM
            JSV = (JMODE-1)*IMOMENTS  & !Number of moments previously counted
                 + JMOM               & !Number of moments in this mode
                 + (NSV_DSTBEG -1)
            TZFIELD%CMNHNAME   = TRIM(YPDUST_INI(ISV_NAME_IDX))//'T'  !The refererence which will be written to file
            TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
            WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
            CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
            YDSTNAMES((JMODE-1)*IMOMENTS+JMOM)=TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
          END DO ! Loop on moment
        END DO ! loop on mode
      END IF ! Valeur IMOMENTS
!
      DO JSV = NSV_DSTBEG,NSV_DSTEND
        YCHNAMES(JSV-JSA) = YDSTNAMES(JSV-NSV_DSTBEG+1)
      END DO   
      DEALLOCATE(YDSTNAMES)
    ELSE 
      ! We are in the subprogram MESONH, CDUSTNAMES are allocated and are 
      !in the same order as the variables in XSVT (i.e. following JPDUSTORDER)
      !
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_DSTBEG,NSV_DSTEND
        TZFIELD%CMNHNAME   = TRIM(CDUSTNAMES(JSV-NSV_DSTBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        IF (JSV==NSV_DSTBEG) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_DSTBEG ',JSV
        IF (JSV==NSV_DSTEND) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_DSTEND ',JSV
        YCHNAMES(JSV-JSA) = TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
      END DO   ! Loop on dust scalar variables
    END IF 
    IF (LDEPOS_DST(IMI)) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_DSTDEPBEG,NSV_DSTDEPEND
        TZFIELD%CMNHNAME   = TRIM(CDEDSTNAMES(JSV-NSV_DSTDEPBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        IF (JSV==NSV_DSTDEPBEG) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_DSTDEPBEG ',JSV
        IF (JSV==NSV_DSTDEPEND) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_DSTDEPEND ',JSV
        YCHNAMES(JSV-JSA) = TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
      END DO   ! Loop on aq dust scalar variables      
    ENDIF
  ENDIF  
  ! sea salt scalar variables
  IF (LSALT) THEN
    IF ((CPROGRAM == 'REAL  ').AND.(NSV_SLT > 1).AND.(IMI==1).AND.(LSLTINIT)) &
      CALL SALTLFI_n(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND), XRHODREF, XZZ)
    IF ((CPROGRAM == 'IDEAL ').AND.(NSV_SLT > 1).AND.(IMI==1)) &
      CALL SALTLFI_n(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND), XRHODREF, XZZ)
    !At this point, we have the tracer array in order of importance, i.e.
    !if mode 2 is most important it will occupy place 1-3 of XSVT  
    IF (((CPROGRAM == 'REAL  ').AND.(LSLTINIT)).OR.&
        (CPROGRAM == 'IDEAL ')           ) THEN
      ! In this case CSALTNAMES is not allocated. We will use YPSALT_INI,
      !but remember that this variable does not follow JPSALTORDER
      IMOMENTS = INT(NSV_SLTEND - NSV_SLTBEG+1)/NMODE_SLT  
      !Should equal 3 at this point
      IF (IMOMENTS .NE. 3) THEN
        WRITE(ILUOUT,*) 'Error in write_lfin: number of moments must be 3'
        WRITE(ILUOUT,*) NSV_SLTBEG, NSV_SLTEND,NMODE_SLT,IMOMENTS
 !callabortstop
        CALL IO_FILE_CLOSE_ll(TLUOUT)
        CALL ABORT
        STOP
      END IF
      ALLOCATE(YSLTNAMES(NSV_SLTEND - NSV_SLTBEG+1))
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      DO JMODE=1, NMODE_SLT
        DO JMOM = 1, IMOMENTS
          !Index from which names are picked
          ISV_NAME_IDX = (JPSALTORDER(JMODE)-1)*IMOMENTS + JMOM 
          !Index which counts in the XSVT
          JSV = (JMODE-1)*IMOMENTS      & !Number of moments previously counted
               + JMOM                   & !Number of moments in this mode
               + (NSV_SLTBEG -1)          !Previous list of tracers 

          TZFIELD%CMNHNAME   = TRIM(YPSALT_INI(ISV_NAME_IDX))//'T'  !The refererence which will be written to file
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
          CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
          YSLTNAMES((JMODE-1)*IMOMENTS+JMOM)=TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
        END DO ! Loop on moments
      END DO   ! Loop on modes
      !
      DO JSV = NSV_SLTBEG,NSV_SLTEND
        YCHNAMES(JSV-JSA) = YSLTNAMES(JSV-NSV_SLTBEG+1)
      END DO   
      DEALLOCATE(YSLTNAMES)
    ELSE 
      ! We are in the subprogram MESONH, CSALTNAMES are allocated and are 
      !in the same order as the variables in XSVT (i.e. following JPSALTORDER)
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_SLTBEG,NSV_SLTEND
        TZFIELD%CMNHNAME   = TRIM(CSALTNAMES(JSV-NSV_SLTBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        IF (JSV==NSV_SLTBEG) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_SLTBEG ',JSV
        IF (JSV==NSV_SLTEND) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_SLTEND ',JSV
        YCHNAMES(JSV-JSA) = TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
      END DO   ! Loop on sea salt scalar variables
    END IF 
    IF (LDEPOS_SLT(IMI)) THEN        
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_SLTDEPBEG,NSV_SLTDEPEND
        TZFIELD%CMNHNAME   = TRIM(CDESLTNAMES(JSV-NSV_SLTDEPBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
        IF (JSV==NSV_SLTDEPBEG) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_SLTDEPBEG ',JSV
        IF (JSV==NSV_SLTDEPEND) WRITE(ILUOUT,*)'MNHC: write_lfin:NSV_SLTDEPEND ',JSV
        YCHNAMES(JSV-JSA) = TZFIELD%CMNHNAME(1:LEN_TRIM(TZFIELD%CMNHNAME)-1)
      END DO   ! Loop on aq dust scalar variables      
    ENDIF
  ENDIF  
  !
  DO JSV=1,ICH_NBR
    WRITE(ILUOUT,*)JSV,TRIM(YCHNAMES(JSV))
  END DO
  TZFIELD%CMNHNAME   = 'NSV.DIM'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'NSV.DIM'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = 'Number of chemical variables'
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ICH_NBR)
  !
  IF (ICH_NBR/=0) THEN
    TZFIELD%CMNHNAME   = 'NSV.TITRE'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'NSV.TITRE'
    TZFIELD%CUNITS     = ''
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 0
    TZFIELD%NTYPE      = TYPEINT
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    ILREC=LEN(YCHNAMES(1))
    ILENG=ILREC*ICH_NBR
    ALLOCATE(ICH_NAMES(ILENG))
    DO JSV = 1,ICH_NBR
      DO JT = 1,ILREC
        ICH_NAMES(ILREC*(JSV-1)+JT) = ICHAR(YCHNAMES(JSV)(JT:JT))
      ENDDO
    ENDDO
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ICH_NAMES)
    DEALLOCATE(YCHNAMES,ICH_NAMES)
  END IF 
  !
  ! lagrangian variables
  IF (NSV_LGEND>=NSV_LGBEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_LGBEG,NSV_LGEND
      TZFIELD%CMNHNAME   = TRIM(CLGNAMES(JSV-NSV_LGBEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
    END DO
  END IF
END IF
!
!
CALL IO_WRITE_FIELD(TPFILE,'LSUM', XLSUM)
CALL IO_WRITE_FIELD(TPFILE,'LSVM', XLSVM)
CALL IO_WRITE_FIELD(TPFILE,'LSWM', XLSWM)
CALL IO_WRITE_FIELD(TPFILE,'LSTHM',XLSTHM)
IF (LUSERV) CALL IO_WRITE_FIELD(TPFILE,'LSRVM',XLSRVM)
!
CALL WRITE_LB_n(TPFILE)
!
!
CALL IO_WRITE_FIELD(TPFILE,'DRYMASST',XDRYMASST)
!
IF( CTURB /= 'NONE' .AND. CTOM=='TM06') THEN
  CALL IO_WRITE_FIELD(TPFILE,'BL_DEPTH',XBL_DEPTH)
END IF
!
IF( CTURB /= 'NONE' .AND. LRMC01) THEN
  CALL IO_WRITE_FIELD(TPFILE,'SBL_DEPTH',XSBL_DEPTH)
END IF
!
IF( CTURB /= 'NONE' .AND. CSCONV == 'EDKF' .AND.(CPROGRAM == 'MESONH' .OR. CPROGRAM == 'DIAG')) THEN
  CALL IO_WRITE_FIELD(TPFILE,'WTHVMF',XWTHVMF)
END IF
!
IF( NRR > 1 .AND. CTURB /= 'NONE' ) THEN
  CALL IO_WRITE_FIELD(TPFILE,'SRCT',XSRCT)
  CALL IO_WRITE_FIELD(TPFILE,'SIGS',XSIGS)
END IF
!
!*       1.5    Reference state variables :
!
CALL IO_WRITE_FIELD(TPFILE,'RHOREFZ',XRHODREFZ)
CALL IO_WRITE_FIELD(TPFILE,'THVREFZ',XTHVREFZ)
CALL IO_WRITE_FIELD(TPFILE,'EXNTOP', XEXNTOP)
!
!
!*       1.6  Tendencies                                         
!
IF (CPROGRAM == 'MESONH') THEN
  IF (CTEMP_SCHEME/='LEFR') THEN
    CALL IO_WRITE_FIELD(TPFILE,'US_PRES',XRUS_PRES)
    CALL IO_WRITE_FIELD(TPFILE,'VS_PRES',XRVS_PRES)
    CALL IO_WRITE_FIELD(TPFILE,'WS_PRES',XRWS_PRES)
  END IF
  IF (LSPLIT_CFL) THEN
    CALL IO_WRITE_FIELD(TPFILE,'THS_CLD',XRTHS_CLD)
!
    IF (NRR >=1) THEN
      IF (LUSERV) CALL IO_WRITE_FIELD(TPFILE,'RVS_CLD',XRRS_CLD(:,:,:,IDX_RVT))
      IF (LUSERC) CALL IO_WRITE_FIELD(TPFILE,'RCS_CLD',XRRS_CLD(:,:,:,IDX_RCT))
      IF (LUSERR) CALL IO_WRITE_FIELD(TPFILE,'RRS_CLD',XRRS_CLD(:,:,:,IDX_RRT))
      IF (LUSERI) CALL IO_WRITE_FIELD(TPFILE,'RIS_CLD',XRRS_CLD(:,:,:,IDX_RIT))
      IF (LUSERS) CALL IO_WRITE_FIELD(TPFILE,'RSS_CLD',XRRS_CLD(:,:,:,IDX_RST))
      IF (LUSERG) CALL IO_WRITE_FIELD(TPFILE,'RGS_CLD',XRRS_CLD(:,:,:,IDX_RGT))
      IF (LUSERH) CALL IO_WRITE_FIELD(TPFILE,'RHS_CLD',XRRS_CLD(:,:,:,IDX_RHT))
    END IF 
  END IF
END IF 
!
!IF (LSPLIT_CFL) THEN
! IF (NSV >=1) THEN
!    DO JSV = NSV_C2R2BEG,NSV_C2R2END
!     IF (JSV == NSV_C2R2BEG ) THEN
!       TZFIELD%CMNHNAME   = 'RSVS_CLD1'
!       TZFIELD%CSTDNAME   = ''
!       TZFIELD%CLONGNAME  = 'RSVS_CLD1'
!       TZFIELD%CUNITS     = '1'
!       TZFIELD%CDIR       = 'XY'
!       TZFIELD%CCOMMENT   = 'X_Y_Z_RHS_CLD'
!       TZFIELD%NGRID      = 1
!       TZFIELD%NTYPE      = TYPEREAL
!       TZFIELD%NDIMS      = 3
!       TZFIELD%LTIMEDEP   = .TRUE.
!       CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XRRS_CLD(:,:,:,IRR))
!     END IF
!     IF (JSV == NSV_C2R2END ) THEN
!       TZFIELD%CMNHNAME   = 'RSVS_CLD2'
!       TZFIELD%CSTDNAME   = ''
!       TZFIELD%CLONGNAME  = 'RSVS_CLD2'
!       TZFIELD%CUNITS     = '1'
!       TZFIELD%CDIR       = 'XY'
!       TZFIELD%CCOMMENT   = 'X_Y_Z_RHS_CLD'
!       TZFIELD%NGRID      = 1
!       TZFIELD%NTYPE      = TYPEREAL
!       TZFIELD%NDIMS      = 3
!       TZFIELD%LTIMEDEP   = .TRUE.
!       CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XRRS_CLD(:,:,:,IRR))
!     END IF
!    END DO
! END IF
!ENDIF
!
!*       1.8    Diagnostic variables related to the radiations
!
!
IF (CRAD /= 'NONE') THEN
  CALL IO_WRITE_FIELD(TPFILE,'DTRAD_FULL',TDTRAD_FULL)
  CALL IO_WRITE_FIELD(TPFILE,'DTRAD_CLLY',TDTRAD_CLONLY)
!
  CALL IO_WRITE_FIELD(TPFILE,'DTHRAD',      XDTHRAD)
  CALL IO_WRITE_FIELD(TPFILE,'FLALWD',      XFLALWD)
  CALL IO_WRITE_FIELD(TPFILE,'DIRFLASWD',   XDIRFLASWD)
  CALL IO_WRITE_FIELD(TPFILE,'SCAFLASWD',   XSCAFLASWD)
  CALL IO_WRITE_FIELD(TPFILE,'DIRSRFSWD',   XDIRSRFSWD)
  CALL IO_WRITE_FIELD(TPFILE,'CLEARCOL_TM1',NCLEARCOL_TM1)
  CALL IO_WRITE_FIELD(TPFILE,'ZENITH',      XZENITH)
  CALL IO_WRITE_FIELD(TPFILE,'AZIM',        XAZIM)
  CALL IO_WRITE_FIELD(TPFILE,'DIR_ALB',     XDIR_ALB)
  CALL IO_WRITE_FIELD(TPFILE,'SCA_ALB',     XSCA_ALB)
  !
  CALL PRINT_MSG(NVERB_INFO,'IO','WRITE_LFIFM_n','EMIS: writing only first band')
  CALL FIND_FIELD_ID_FROM_MNHNAME('EMIS',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%NDIMS = 2
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XEMIS(:,:,1))
  !
  CALL IO_WRITE_FIELD(TPFILE,'TSRAD',       XTSRAD)
ENDIF
!
IF (NRR > 1 .AND. CPROGRAM == 'MESONH') THEN
  CALL IO_WRITE_FIELD(TPFILE,'CLDFR',XCLDFR)
  CALL IO_WRITE_FIELD(TPFILE,'RAINFR',XRAINFR)
END IF
!
!
!*       1.9     Diagnostic variables related to deep convection
!
!
IF (CDCONV /= 'NONE' .OR. CSCONV == 'KAFR') THEN
!
! 
!
  CALL IO_WRITE_FIELD(TPFILE,'DTDCONV',  TDTDCONV)
  CALL IO_WRITE_FIELD(TPFILE,'COUNTCONV',NCOUNTCONV)
  CALL IO_WRITE_FIELD(TPFILE,'DTHCONV',  XDTHCONV)
  CALL IO_WRITE_FIELD(TPFILE,'DRVCONV',  XDRVCONV)
  CALL IO_WRITE_FIELD(TPFILE,'DRCCONV',  XDRCCONV)
  CALL IO_WRITE_FIELD(TPFILE,'DRICONV',  XDRICONV)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PRCONV',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XPRCONV*3.6E6)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PACCONV',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XPACCONV*1.0E3)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PRSCONV',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XPRSCONV*3.6E6)
!
  IF ( LCH_CONV_LINOX ) THEN
    CALL IO_WRITE_FIELD(TPFILE,'IC_RATE',    XIC_RATE)
    CALL IO_WRITE_FIELD(TPFILE,'CG_RATE',    XCG_RATE)
    CALL IO_WRITE_FIELD(TPFILE,'IC_TOTAL_NB',XIC_TOTAL_NUMBER)
    CALL IO_WRITE_FIELD(TPFILE,'CG_TOTAL_NB',XCG_TOTAL_NUMBER)
  END IF
!
  IF ( LCHTRANS .AND. NSV > 0 ) THEN
   ! scalar variables are recorded
   ! individually in the file
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = 1, NSV_USER
      WRITE(TZFIELD%CMNHNAME,'(A7,I3.3)')'DSVCONV',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_C2R2BEG, NSV_C2R2END
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_C1R3BEG, NSV_C1R3END
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_ELECBEG, NSV_ELECEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_PPBEG, NSV_PPEND
      WRITE(TZFIELD%CMNHNAME,'(A7,I3.3)')'DSVCONV',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
#ifdef MNH_FOREFIRE
    IF (LFOREFIRE) THEN
      DO JSV = NSV_FFBEG, NSV_FFEND
        WRITE(TZFIELD%CMNHNAME,'(A7,I3.3)')'DSVCONV',JSV
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
#endif
    IF (LUSECHEM) THEN
      DO JSV = NSV_CHEMBEG, NSV_CHEMEND
        TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHEMBEG+1)))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
      IF (LORILAM) THEN
        DO JSV = NSV_AERBEG, NSV_AEREND
          TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(UPCASE(CAERONAMES(JSV-NSV_AERBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
          CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
        END DO
      END IF
! linox scalar variables
    ELSE IF (LCH_CONV_LINOX) THEN
      DO JSV = NSV_LNOXBEG,NSV_LNOXEND
        TZFIELD%CMNHNAME   = 'DSVCONV_LINOX'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
        CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
    DO JSV = NSV_LGBEG, NSV_LGEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CLGNAMES(JSV-NSV_LGBEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_DSTBEG, NSV_DSTEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CDUSTNAMES(JSV-NSV_DSTBEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_SLTBEG, NSV_SLTEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CSALTNAMES(JSV-NSV_SLTBEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
    END DO
  END IF
!
END IF
!
!
!*       1.10   Diagnostic variables related to the precipitations
!
IF (CPROGRAM /= 'IDEAL') THEN
  IF (ASSOCIATED(XINPRC)) THEN
  IF (SIZE(XINPRC) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRC',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XINPRC*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRC',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XACPRC*1.0E3)
!
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINDEP)) THEN
  IF (SIZE(XINDEP) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INDEP',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XINDEP*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACDEP',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XACDEP*1.0E3)
!
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINPRR)) THEN
  IF (SIZE(XINPRR) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRR',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XINPRR*3.6E6)
!
    CALL IO_WRITE_FIELD(TPFILE,'INPRR3D',XINPRR3D)
    CALL IO_WRITE_FIELD(TPFILE,'EVAP3D', XEVAP3D)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRR',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XACPRR*1.0E3)
!
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINPRS)) THEN
  IF (SIZE(XINPRS) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRS',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XINPRS*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRS',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XACPRS*1.0E3)
  END IF
  END IF
!
  IF (ASSOCIATED(XINPRG)) THEN
  IF (SIZE(XINPRG) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XINPRG*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XACPRG*1.0E3)
  END IF
  END IF
!
  IF (ASSOCIATED(XINPRH)) THEN
  IF (SIZE(XINPRH) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRH',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XINPRH*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRH',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XACPRH*1.0E3)
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINPRS)) THEN
  IF (SIZE(XINPRS) /= 0 ) THEN
    ZWORK2D = XINPRR + XINPRS
    IF (SIZE(XINPRG) /= 0 ) ZWORK2D = ZWORK2D + XINPRG
    IF (SIZE(XINPRH) /= 0 ) ZWORK2D = ZWORK2D + XINPRH
    IF (SIZE(XINPRC) /= 0 ) ZWORK2D = ZWORK2D + XINPRC
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRT',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2D*3.6E6)
!
    ZWORK2D = XACPRR + XACPRS
    IF (SIZE(XINPRG) /= 0 ) ZWORK2D = ZWORK2D + XACPRG
    IF (SIZE(XINPRH) /= 0 ) ZWORK2D = ZWORK2D + XACPRH
    IF (SIZE(XINPRC) /= 0 ) ZWORK2D = ZWORK2D + XACPRC
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRT',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2D*1.0E3)
  END IF
  END IF
!
END IF
!
IF(LBLOWSNOW) THEN
  IF (ASSOCIATED(XSNWSUBL3D)) THEN
    IF (SIZE(XSNWSUBL3D) /= 0 ) THEN
      TZFIELD%CMNHNAME   = 'SNWSUBL3D'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'kg m-3 s-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_INstantaneous 3D Drifting snow sublimation flux'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XSNWSUBL3D(:,:,:))
      ZWORK2D(:,:) = 0.
      DO JK = IKB,IKE
        ZWORK2D(:,:) = ZWORK2D(:,:)+XSNWSUBL3D(:,:,JK) * &
                    (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW*3600*24
      END DO
      ZWORK2D(:,:) = ZWORK2D(:,:)*1000. ! vapor water in mm unit
      !
      TZFIELD%CMNHNAME   = 'COL_SNWSUBL'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'mm day-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Column Sublimation Rate (mmSWE/day)'
      TZFIELD%NGRID      = 4
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 2
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2D(:,:))
    END IF
  END IF
ENDIF
!
!*       1.11   Forcing variables
!
!
IF (LFORCING) THEN
!
  CALL IO_WRITE_FIELD(TPFILE,'FRC',NFRC)
!
  DO JT=1,NFRC
!
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD%CMNHNAME   = 'DTFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Date of forcing profile '//YFRC
    TZFIELD%NGRID      = 0
    TZFIELD%NTYPE      = TYPEDATE
    TZFIELD%NDIMS      = 0
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,TDTFRC(JT))
!
    TZFIELD%CMNHNAME   = 'UFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Zonal component of horizontal forcing wind'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XUFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'VFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Meridian component of horizontal forcing wind'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XVFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'WFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Vertical forcing wind'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XWFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'THFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Forcing potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'RVFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Forcing vapor mixing ratio'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XRVFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'TENDTHFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale potential temperature tendency for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTENDTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'TENDRVFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'kg kg-1 s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale vapor mixing ratio tendency for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTENDRVFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'GXTHFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K m-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale potential temperature gradient for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XGXTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'GYTHFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K m-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale potential temperature gradient for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XGYTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'PGROUNDFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'Pa'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Forcing ground pressure'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 0
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XPGROUNDFRC(JT))
!
    TZFIELD%CMNHNAME   = 'TENDUFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale U tendency for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTENDUFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'TENDVFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale V tendency for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTENDVFRC(:,JT))
!
  END DO
!
!
END IF
!
! -------------------------------------------------------------------------
IF ( L2D_ADV_FRC ) THEN
!
  TZFIELD%CMNHNAME   = 'NADVFRC1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'NADVFRC1'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = 'Number of forcing profiles'
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,NADVFRC)
!
  DO JT=1,NADVFRC
!
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD%CMNHNAME   = 'DTADV'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Date and time of the advecting forcing '//YFRC
    TZFIELD%NGRID      = 0
    TZFIELD%NTYPE      = TYPEDATE
    TZFIELD%NDIMS      = 0
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,TDTADVFRC(JT))
!                                                                
    TZFIELD%CMNHNAME   = 'TH_ADV'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDTHFRC(:,:,:,JT))
!    
    TZFIELD%CMNHNAME   = 'Q_ADV'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'kg kg-1 s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDRVFRC(:,:,:,JT))
!
  ENDDO
ENDIF
!
IF ( L2D_REL_FRC ) THEN
!
  TZFIELD%CMNHNAME   = 'NRELFRC1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'NRELFRC1'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = 'Number of forcing profiles'
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,NRELFRC)
!
  DO JT=1,NRELFRC
!
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD%CMNHNAME   = 'DTREL'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Date and time of the relaxation forcing '//YFRC
    TZFIELD%NGRID      = 0
    TZFIELD%NTYPE      = TYPEDATE
    TZFIELD%NDIMS      = 0
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,TDTRELFRC(JT))
!                                                                
    TZFIELD%CMNHNAME   = 'TH_REL'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XTHREL(:,:,:,JT))
!    
    TZFIELD%CMNHNAME   = 'Q_REL'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XRVREL(:,:,:,JT))
!
  ENDDO
ENDIF
!
!*       1.11bis   Eddy Fluxes variables    ! Modif PP
!
IF ( LTH_FLX ) THEN
   CALL IO_WRITE_FIELD(TPFILE,'VT_FLX',XVTH_FLUX_M)
   CALL IO_WRITE_FIELD(TPFILE,'WT_FLX',XWTH_FLUX_M)
END IF
!
IF ( LUV_FLX) CALL IO_WRITE_FIELD(TPFILE,'VU_FLX',XVU_FLUX_M)
!
!*       1.12   Balloon variables
!
!
IF (LFLYER) CALL WRITE_BALLOON_n(TPFILE)
!
!
!*       1.13    Filtered variables for hurricane initialization
!
!
IF ( CPROGRAM=='REAL  ' ) THEN
  IF (LFILTERING) THEN
  !
    IF (NDIAG_FILT >=0) THEN
!
!             i) Total fields (TOT=BASIC+TOTDIS)
!
      CALL IO_WRITE_FIELD(TPFILE,'UT15',   XUTOT)
      CALL IO_WRITE_FIELD(TPFILE,'VT15',   XVTOT)
      CALL IO_WRITE_FIELD(TPFILE,'TEMPTOT',XTTOT)
      IF (INDEX(CFILTERING,'P')/=0) CALL IO_WRITE_FIELD(TPFILE,'PRESTOT',XPTOT)
      IF (INDEX(CFILTERING,'Q')/=0) CALL IO_WRITE_FIELD(TPFILE,'HUMTOT', XQTOT)
!
!             ii) Environmental fields (ENV=TOT-VORDIS)
!
      CALL IO_WRITE_FIELD(TPFILE,'UT16',   XUENV)
      CALL IO_WRITE_FIELD(TPFILE,'VT16',   XVENV)
      CALL IO_WRITE_FIELD(TPFILE,'TEMPENV',XTENV)
      IF (INDEX(CFILTERING,'P')/=0) CALL IO_WRITE_FIELD(TPFILE,'PRESENV',XPENV)
      IF (INDEX(CFILTERING,'Q')/=0) CALL IO_WRITE_FIELD(TPFILE,'HUMENV', XQENV)
!
    END IF
    IF (NDIAG_FILT >=1) THEN
!
!             iii) Basic (filtered) fields
!
      CALL IO_WRITE_FIELD(TPFILE,'UT17',   XUBASIC)
      CALL IO_WRITE_FIELD(TPFILE,'VT17',   XVBASIC)
      CALL IO_WRITE_FIELD(TPFILE,'TEMPBAS',XTBASIC)
      IF (INDEX(CFILTERING,'P')/=0) CALL IO_WRITE_FIELD(TPFILE,'PRESBAS',XPBASIC)
      IF (INDEX(CFILTERING,'Q')/=0) CALL IO_WRITE_FIELD(TPFILE,'HUMBAS', XQBASIC)
    END IF
    IF (NDIAG_FILT >=2) THEN
!
!             iv) Total disturbance tangential wind
!
      CALL IO_WRITE_FIELD(TPFILE,'VTDIS',XVTDIS)
!
    END IF
!
  END IF
!
!*       1.14    Dummy variables in PREP_REAL_CASE
!
  IF (ALLOCATED(CDUMMY_2D)) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSA=1,SIZE(XDUMMY_2D,3)
      TZFIELD%CMNHNAME   = ADJUSTL(CDUMMY_2D(JSA))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
      CALL IO_WRITE_FIELD(TPFILE,TZFIELD,XDUMMY_2D(:,:,JSA))
    END DO
  END IF
!
END IF
!
!
DEALLOCATE(ZWORK2D,ZWORK3D)
!
!-------------------------------------------------------------------------------!
!
END SUBROUTINE WRITE_LFIFM_n  
