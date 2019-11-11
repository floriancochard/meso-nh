!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! BUG2 masdev4_7 2007/09/19 12:46:39
!-----------------------------------------------------------------
!     #################
      MODULE MODI_INI_MODEL_n
!     #################
!
INTERFACE
!
       SUBROUTINE INI_MODEL_n(KMI,PTSTEP_OLD,HLUOUT,HINIFILE)
!
       INTEGER, INTENT(IN)              :: KMI      ! Model index 
       REAL,              INTENT(IN)    :: PTSTEP_OLD ! OLD Time STEP 
       CHARACTER (LEN=*), INTENT(IN)    :: HLUOUT   ! name for output-listing
       !  of nested models
       CHARACTER (LEN=28), INTENT(IN)   :: HINIFILE ! name of
!
END SUBROUTINE INI_MODEL_n
!
END INTERFACE
!
END MODULE MODI_INI_MODEL_n
!     ######################################################
      SUBROUTINE INI_MODEL_n(KMI,PTSTEP_OLD,HLUOUT,HINIFILE)
!     ######################################################
!
!!****  *INI_MODEL_n* - routine to initialize the nested model _n
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the variables
!     of the nested model _n.
!
!!**  METHOD
!!    ------
!!      The initialization of the model _n is performed as follows :
!!       - Memory for arrays are then allocated :
!!            * If  turbulence kinetic energy variable is not needed
!!    (CTURB='NONE'),  XTKET, XTKEM and XTKES are zero-size arrays.
!!            * If  dissipation of TKE variable is not needed
!!    (CTURBLEN /='KEPS'),  XEPST, XEPSM and XREPSS are zero-size arrays.
!!            * Memory for mixing ratio arrays is allocated according to the
!!     value of logicals LUSERn (the number NRR of moist variables is deduced).
!!            * The latitude (XLAT), longitude (XLON) and map factor (XMAP)
!!     arrays are zero-size arrays if Cartesian geometry (LCARTESIAN=.TRUE.)
!!            * Memory for reference state without orography ( XRHODREFZ and
!!     XTHVREFZ) is only allocated in INI_MODEL1
!!            * The horizontal Coriolis parameters (XCORIOX and XCORIOY) arrays
!!     are  zero-size arrays if thinshell approximation (LTHINSHELL=.TRUE.)
!!            * The Curvature coefficients (XCURVX and XCURVY) arrays
!!     are zero-size arrays if Cartesian geometry (LCARTESIAN=.TRUE.)
!!            * Memory for the Jacobian (ZJ) local array is allocated
!!     (This variable is computed in SET_GRID and used in SET_REF).
!!       - The spatial and temporal grid variables are initialized by SET_GRID.
!!       - The metric coefficients are computed by METRICS (they are using in
!!     the SET-REF call).
!!       - The prognostic variables and are read in initial
!!    LFIFM file (in READ_FIELD)
!!       - The reference state variables  are initialized by SET_REF.
!!       - The temporal indexes of the outputs are computed by SET_OUTPUT_TIMES
!!       - The large scale sources are computed in case of coupling case by
!!    INI_CPL.
!!       - The initialization of the parameters needed for the dynamics
!!         of the model n is realized in INI_DYNAMICS.
!!       - Then the initial file (DESFM+LFIFM files) is closed by FMCLOS.
!!       - The initialization of the parameters needed for the ECMWF radiation
!!         code is realized in INI_RADIATIONS.
!!       - The contents of the scalar variables are overwritten by
!!         the chemistry initialization subroutine CH_INIT_FIELDn when
!!         the flags LUSECHEM and LCH_INIT_FIELD are set to TRUE.
!!         This allows easy initialization of the chemical fields at a
!!         restart of the model.
!!
!!    EXTERNAL
!!    --------
!!      FMLOOK      : to retrieve a logical unit number associated with a file
!!      FMREAD      : to read a LFIFM file
!!      FMFREE      : to release a logical unit number
!!      SET_DIM     : to initialize dimensions
!!      SET_GRID    : to initialize grid
!!      METRICS     : to compute metric coefficients
!!      READ_FIELD  : to initialize field
!!      FMCLOS      : to close a FM-file
!!      SET_REF     : to initialize reference state for anelastic approximation
!!      INI_DYNAMICS: to initialize parameters for the dynamics
!!      INI_TKE_EPS : to initialize the TKE 
!!      SET_DIRCOS  : to compute the director cosinus of the orography
!!      INI_RADIATIONS : to initialize radiation computations
!!      CH_INIT_CCS: to initialize the chemical core system
!!      CH_INIT_FIELDn: to (re)initialize the scalar variables
!!      INI_DEEP_CONVECTION : to initialize the deep convection scheme
!!      CLEANLIST_ll : deaalocate a list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS : contains declaration of parameter variables
!!         JPHEXT : Horizontal external points number
!!         JPVEXT : Vertical external points number
!!
!!      Module MODD_MODD_DYN   : contains declaration of parameters
!!                               for the dynamics
!!      Module MODD_CONF       : contains declaration of configuration variables
!!                               for all models
!!         NMODEL     : Number of nested models
!!         NVERB      : Level of informations on output-listing
!!                          0 for minimum  prints
!!                          5 for intermediate level of prints
!!                         10 for maximum  prints
!!
!!      Module MODD_REF        : contains declaration of reference state
!!                               variables for all models
!!      Module MODD_FIELD_n    : contains declaration of prognostic fields
!!      Module MODD_LSFIELD_n  : contains declaration of Larger Scale fields
!!      Module MODD_GRID_n     : contains declaration of spatial grid variables
!!      Module MODD_TIME_n     : contains declaration of temporal grid variables
!!      Module MODD_REF_n      : contains declaration of reference state
!!                               variables
!!      Module MODD_CURVCOR_n  : contains declaration of curvature and Coriolis
!!                               variables
!!      Module MODD_BUDGET     : contains declarations of the budget parameters
!!      Module MODD_RADIATIONS_n:contains declaration of the variables of the
!!                               radiation interface scheme
!!      Module MODD_STAND_ATM  : contains declaration of the 5 standard
!!                               atmospheres used for the ECMWF-radiation code
!!      Module MODD_FRC        : contains declaration of the control variables
!!                               and of the forcing fields
!!      Module MODD_CH_MNHC_n   : contains the control parameters for chemistry
!!      Module MODD_DEEP_CONVECTION_n: contains declaration of the variables of
!!                                     the deep convection scheme
!!        
!!        
!!        
!!        
!!      Module MODN_CONF_n     : contains declaration of namelist NAM_CONFn and 
!!                             uses module MODD_CONF_n (configuration variables)
!!      Module MODN_LUNIT_n    : contains declaration of namelist NAM_LUNITn and
!!                             uses module MODD_LUNIT_n (Logical units)
!!      Module MODN_DYN_n      : contains declaration of namelist NAM_DYNn and
!!                             uses module MODD_DYN_n (control of dynamics)
!!      Module MODN_PARAM_n    : contains declaration of namelist NAM_PARAMn and
!!                             uses module MODD_PARAM_n (control of physical
!!                             parameterization)
!!      Module MODN_LBC_n      : contains declaration of namelist NAM_LBCn and 
!!                             uses module MODD_LBC_n (lateral boundaries)
!!      Module MODN_TURB_n     : contains declaration of namelist NAM_TURBn and 
!!                             uses module MODD_TURB_n (turbulence scheme)
!!      Module MODN_PARAM_RAD_n: contains declaration of namelist NAM_PARAM_RADn
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_MODEL_n)
!!      
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/06/94
!!      Modification 17/10/94  (Stein)  For LCORIO
!!      Modification 20/10/94  (Stein)  For SET_GRID and NAMOUTN
!!      Modification 26/10/94  (Stein)  Modifications of the namelist names
!!      Modification 10/11/94  (Lafore) allocatation of tke fields
!!      Modification 22/11/94  (Stein)  change the READ_FIELDS call ( add
!!                                      pressure function
!!      Modification 06/12/94  (Stein)  add the LS fields
!!                   12/12/94  (Stein)  rename END_INI in INI_DYNAMICS
!!      Modification 09/01/95  (Stein)  add the turbulence scheme
!!      Modification Jan 19, 1995 (J. Cuxart) add the TKE initialization
!!                   Jan 23, 1995 (J. Stein ) remove the condition
!!                             LTHINSHELL=T LCARTESIAN=T => stop
!!      Modification Feb 16, 1995 (I.Mallet) add the METRICS call and
!!                                      change the SET_REF call (add
!!                                      the lineic mass)
!!      Modification Mar 10, 1995 (I. Mallet) add the COUPLING initialization
!!                   June 29,1995 (Ph. Hereil, J. Stein) add the budget init.
!!      Modification Sept. 1, 1995 (S. Belair) Reading of the surface variables
!!                                 and parameters for ISBA (i.e., add a
!!                                 CALL READ_GR_FIELD)
!!      Modification 18/08/95     (J.P.Lafore)   time step change case
!!                   25/09/95     (J. Cuxart and J.Stein)   add LES variables
!!                                and the diachronic file initialization
!!      Modification Sept 20,1995 (Lafore) coupling for the dry mass Md
!!      Modification Sept. 12, 1995 (J.-P. Pinty) add the initialization of
!!                                      the ECMWF radiation code
!!      Modification Sept. 13, 1995 (J.-P. Pinty) control the allocation of the
!!                                      arrays of MODD_GR_FIELD_n
!!      Modification Nove. 17, 1995 (J.Stein) control of the control !! 
!!                   March 01, 1996 (J. Stein) add the cloud fraction    
!!                   April 03, 1996 (J. Stein) unify the ISBA and TSZ0 cases
!!      Modification 13/12/95 (M. Georgelin) add the forcing variables in
!!                                           the call read_field, and their
!!                                           allocation.
!!                   Mai   23, 1996 (J. Stein) allocate XSEA in the TSZ0 case
!!                   June  11, 1996 (V. Masson) add XSILT and XLAKE of
!!                                              MODD_GR_FIELD_n
!!                   August 7, 1996 (K. Suhre)  add (re)initialization of
!!                                              chemistry
!!                   Octo. 11, 1996 (J. Stein ) add XSRCT and XSRCM
!!                   October 8, 1996 (J. Cuxart, E. Sanchez) Moist LES diagnostics
!!                                     and control on TKE initialization.
!!      Modification 19/12/96 (J.-P. Pinty) add the ice parameterization and
!!                                          the precipitation fields
!!      Modification 11/01/97 (J.-P. Pinty) add the deep convection
!!                   Nov.   1, 1996 (V. Masson) Read the vertical grid kind
!!                   Nov.  20, 1996 (V. Masson) control of convection calling time
!!                   July  16, 1996 (J.P.Lafore) update of EXSEG file reading
!!                   Oct.  08, 1996 (J.P.Lafore, V.Masson)
!!                                       MY_NAME and DAD_NAME reading and check
!!                   Oct.  30, 1996 (J.P.Lafore) resolution ratio reading for nesting
!!                                       and Bikhardt interpolation coef. initialization
!!                   Nov.  22, 1996 (J.P.Lafore) allocation of LS sources for nesting
!!                   Feb.  26, 1997 (J.P.Lafore) allocation of "surfacic" LS fields
!!                   March 10, 1997 (J.P.Lafore) forcing only for model 1
!!                   June  22, 1997 (J. Stein)   add the absolute pressure
!!                   July  09, 1997 (V. Masson)  add directional z0 and SSO
!!                   Aug.  18, 1997 (V. Masson)  consistency between storage
!!                                               type and CCONF
!!                   Dec.  22, 1997 (J. Stein)   add the LS field spawning
!!                   Jan.  24, 1998 (P.Bechtold) change MODD_FRC and MODD_DEEP_CONVECTION
!!                   Dec.  24, 1997 (V.Masson)   directional z0 parameters
!!                   Aug.  13, 1998 (V. Ducrocq P Jabouille)   //
!!                   Mai.  26, 1998 (J. Stein)   remove NXEND,NYEND
!!                   Feb.   1, 1999 (J. Stein)   compute the Bikhardt
!!                                       interpolation coeff. before the call to set_grid
!!                   April  5, 1999 (V. Ducrocq) change the DXRATIO_ALL init.
!!                   April  12, 1999 (J. Stein)  cleaning + INI_SPAWN_LS
!!                   Apr.   7, 1999 (P Jabouille) store the metric coefficients
!!                                                in modd_metrics_n
!!                   Jui.   15,1999 (P Jabouille) split the routines in two parts
!!                   Jan.   04,2000 (V. Masson)   removes the TSZ0 case
!!                   Apr.   15,2000 (P Jabouille) parallelization of grid nesting
!!                   Aug.   20,2000 (J Stein    ) tranpose XBFY
!!                   Jui    01,2000 (F.solmon )   adapatation for patch approach
!!                   Jun.   15,2000 (J.-P. Pinty) add C2R2 initialization
!!                   Nov.  15,2000 (V.Masson) use of ini_modeln in prep_real_case
!!                   Nov.  15,2000 (V.Masson) call of LES routines
!!                   Nov.  15,2000 (V.Masson) aircraft and balloon initialization routines
!!                   Jan.  22,2001 (D.Gazen) update_nsv set NSV_* var. for current model
!!                   Mar.  04,2002 (V.Ducrocq) initialization to temporal series
!!                   Mar.  15,2002 (F.Solmon) modification of ini_radiation interface
!!                   Nov.  29,2002 (JP Pinty) add C3R5, ICE2, ICE4, ELEC
!!                   Jan.  2004    (V.Masson) externalization of surface
!!                   May   2006    Remove KEPS
!---------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODE_IO_ll
USE MODE_FM
USE MODE_FMREAD
USE MODE_TYPE_ZDIFFU
!
USE MODD_NSV, ONLY : NSV,CSV
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CONF
USE MODD_DUST
USE MODD_DYN
USE MODD_DYNZD
USE MODD_FRC
USE MODD_REF
USE MODD_SERIES, ONLY: LSERIES
USE MODD_TIME
USE MODD_TURB_CLOUD, ONLY: NMODEL_CLOUD, CTURBLEN_CLOUD,XCEI
USE MODD_NESTING
USE MODD_METRICS_n
USE MODD_DYN_n
USE MODD_DYNZD_n
USE MODD_FIELD_n
USE MODD_ADV_n
USE MODD_LSFIELD_n
USE MODD_GRID_n
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_TIME_n
USE MODD_REF_n
USE MODD_CURVCOR_n
USE MODD_DIM_n
USE MODD_BUDGET
USE MODD_RADIATIONS_n,  ONLY : XSLOPANG, XSLOPAZI, XDTHRAD, NCLEARCOL_TM1, XSTATM,    &
                               NDLON,NFLEV,NFLUX,NRAD,NSWB,NAER,NSTATM,               &
                               XCCO2,XSINDEL,XCOSDEL,XTSIDER,XCORSOL,                 &
                               XOZON, XFLALWD, XAER, &
                               XDIRFLASWD, XSCAFLASWD, XDIRSRFSWD,                     &
                               NSWB_MNH, XSW_BANDS, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, &
                               XZENITH, NSWB_MNH, XALBUV, XAZIM, XSEA, XZS_XY
USE MODD_PARAM_RAD_n,   ONLY : CLW, CAER
!
USE MODD_STAND_ATM,     ONLY : XSTROATM, XSMLSATM, XSMLWATM, XSPOSATM, XSPOWATM
USE MODD_CH_MNHC_n, ONLY: LUSECHEM, LCH_INIT_FIELD, CCHEM_INPUT_FILE, LCH_CONV_LINOX
USE MODD_CH_AEROSOL, ONLY : LORILAM
USE MODD_PARAM_CONVECT_n
USE MODD_DEEP_CONVECTION_n
USE MODD_OUT_n
USE MODD_BIKHARDT_n
USE MODD_NUDGING_n, ONLY : LNUDGING
USE MODD_DIAG_FLAG, ONLY : LCHEMDIAG
!
!USE MODN_TURB_n
!USE MODN_CONF_n
!USE MODN_ADV_n
!USE MODN_LUNIT_n
!USE MODN_DYN_n
!USE MODN_PARAM_n
!USE MODN_LBC_n
!USE MODN_PARAM_RAD_n
!USE MODN_PARAM_CONVECT_n
!USE MODN_CH_MNHC_n

USE MODD_LUNIT_n
USE MODD_CONF_n
USE MODD_GET_n
USE MODD_TURB_n
USE MODD_LBC_n
!
!
USE MODI_INI_BUDGET
USE MODI_INI_SW_SETUP
USE MODI_SET_GRID
USE MODI_METRICS
USE MODI_UPDATE_METRICS
USE MODI_READ_FIELD
USE MODI_SET_REF
USE MODI_INI_DYNAMICS
USE MODI_INI_TKE_EPS
USE MODI_SET_DIRCOS
USE MODI_INI_CPL
USE MODI_INI_RADIATIONS
USE MODI_INI_RADIATIONS_ECMWF
USE MODI_CH_INIT_CCS
USE MODI_CH_INIT_FIELD_n
USE MODI_INI_DEEP_CONVECTION
USE MODI_INI_BIKHARDT_n
USE MODI_INI_ONE_WAY_n
USE MODI_GET_SIZEX_LB
USE MODI_GET_SIZEY_LB
USE MODI_INI_SPAWN_LS_n
USE MODI_INI_AIRCRAFT_BALLOON
USE MODI_UPDATE_NSV
USE MODI_INI_ELEC
USE MODI_INI_LG
USE MODI_SURF_SOLAR_GEOM
USE MODI_SUNPOS_n
USE MODI_INI_SURF_RAD
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_MNHREAD_ZS_DUMMY_n
USE MODI_INIT_GROUND_PARAM_n
USE MODI_INI_AIRCRAFT_BALLOON
USE MODI_INI_SURFSTATION_n
USE MODI_INI_POSPROFILER_n
USE MODI_CH_INIT_JVALUES
USE MODI_CH_AER_MOD_INIT
!
USE MODD_PARAM_n
USE MODE_MODELN_HANDLER

USE MODI_TEMPORAL_DIST
!JUAN P1/P2 SPLITTING
USE MODE_SPLITTINGZ_ll , ONLY : GET_DIM_EXTZ_ll
!JUAN
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER, INTENT(IN)              :: KMI      ! Model Index 
REAL,    INTENT(IN)              :: PTSTEP_OLD ! OLD Time STEP 

CHARACTER (LEN=*), INTENT(IN)    :: HLUOUT   ! name for output-listing
                                             !  of nested models
CHARACTER (LEN=28),  INTENT(IN)   :: HINIFILE ! name of
                                             ! the initial file
!
!*       0.2   declarations of local variables
!
INTEGER             :: IRESP   ! Return code of FM routines
INTEGER             :: IMASDEV ! version of MESOHN in the input file
INTEGER             :: ILUOUT  ! Logical unit number of output-listing
CHARACTER(LEN=2)    :: YDIR   ! Type  of the data field in LFIFM file
INTEGER             :: IGRID   ! C-grid indicator in LFIFM file
INTEGER             :: ILENCH  ! Length of comment string in LFIFM file
CHARACTER (LEN=100) :: YCOMMENT!comment string in LFIFM file
CHARACTER (LEN=16)  :: YRECFM  ! Name of the desired field in LFIFM file
INTEGER             :: IIU     ! Upper dimension in x direction (local)
INTEGER             :: IJU     ! Upper dimension in y direction (local)
INTEGER             :: IIU_ll  ! Upper dimension in x direction (global)
INTEGER             :: IJU_ll  ! Upper dimension in y direction (global)
INTEGER             :: IKU     ! Upper dimension in z direction
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZJ ! Jacobian
LOGICAL             :: GINIDCONV ! logical switch for the deep convection
                               ! initialization
LOGICAL             :: GINIRAD ! logical switch for the radiation
                               ! initialization
!
!
TYPE(LIST_ll), POINTER :: TZINITHALO2D_ll ! pointer for the list of 2D fields
                                      !  which must be communicated in INIT
TYPE(LIST_ll), POINTER :: TZINITHALO3D_ll ! pointer for the list of 3D fields
                                      !  which must be communicated in INIT
!
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
INTEGER :: IINFO_ll  ! Return code of //routines
INTEGER :: IIY,IJY
!JUAN Z_SPLITING
INTEGER :: IIU_B,IJU_B
INTEGER :: IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll
!JUAN Z_SPLITING
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZCO2   ! CO2 concentration near the surface
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCOVER ! surface cover types
INTEGER                             :: ICOVER ! number of cover types
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDIR_ALB ! direct albedo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSCA_ALB ! diffuse albedo
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZEMIS    ! emissivity
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTSRAD   ! surface temperature
!------------------------------------------
! Dummy pointers needed to correct an ifort Bug
REAL, DIMENSION(:), POINTER :: DPTR_XZHAT
REAL, DIMENSION(:), POINTER :: DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4
REAL, DIMENSION(:), POINTER :: DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4
REAL, DIMENSION(:), POINTER :: DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4
REAL, DIMENSION(:), POINTER :: DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4
CHARACTER(LEN=4), DIMENSION(:), POINTER :: DPTR_CLBCX,DPTR_CLBCY
INTEGER, DIMENSION(:,:,:), POINTER :: DPTR_NKLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_NKLIN_LBXV,DPTR_NKLIN_LBYV
INTEGER, DIMENSION(:,:,:), POINTER :: DPTR_NKLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_NKLIN_LBXM,DPTR_NKLIN_LBYM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXU,DPTR_XCOEFLIN_LBYU
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXV,DPTR_XCOEFLIN_LBYV
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXW,DPTR_XCOEFLIN_LBYW
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXM,DPTR_XCOEFLIN_LBYM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXWM,DPTR_XLBYWM,DPTR_XLBXTHM,DPTR_XLBYTHM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXTKEM,DPTR_XLBYTKEM
REAL, DIMENSION(:,:,:,:),   POINTER :: DPTR_XLBXSVM,DPTR_XLBYSVM                
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XLBXRM,DPTR_XLBYRM
REAL, DIMENSION(:,:,:),   POINTER ::  DPTR_XZZ
REAL, DIMENSION(:,:,:), POINTER ::   DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSTHM,DPTR_XLSRVM
REAL, DIMENSION(:,:,:), POINTER ::   DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS,DPTR_XLSTHS,DPTR_XLSRVS

!-------------------------------------------------------------------------------
!
!*       0.    PROLOGUE
!              --------
!
NULLIFY(TZINITHALO2D_ll)
NULLIFY(TZINITHALO3D_ll)
!
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
CALL FMLOOK_ll(HLUOUT,HLUOUT,ILUOUT,IRESP)
CLUOUT = HLUOUT
CINIFILE=HINIFILE
!
!-------------------------------------------------------------------------------
!
!*       2.   END OF READING
!             --------------
!*       2.1  Read number of forcing fields
!
IF (LFORCING) THEN ! Retrieve the number of time-dependent forcings.
  YRECFM='FRC'
  YDIR='--'
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,NFRC,IGRID,ILENCH,YCOMMENT,IRESP)
  IF ( (IRESP /= 0) .OR. (NFRC <=0) ) THEN
    WRITE(ILUOUT,'(A/A)') &
     "INI_MODEL_n ERROR: you want to read forcing variables from FMfile", &
     "                   but no fields have been found by FMREAD"
    STOP 1
  END IF
END IF
!
!*       2.2  Checks the position of vertical absorbing layer
!
IKU=NKMAX+2*JPVEXT
!
YRECFM = 'ZHAT'
ALLOCATE(XZHAT(IKU))
 YDIR='--'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,XZHAT,IGRID,ILENCH,YCOMMENT,IRESP)
IF (XALZBOT>=XZHAT(IKU) .AND. LVE_RELAX) THEN
  WRITE(ILUOUT,FMT=*) "INI_MODEL_n ERROR: you want to use vertical relaxation"
  WRITE(ILUOUT,FMT=*) "                  but bottom of layer XALZBOT(",XALZBOT,")"
  WRITE(ILUOUT,FMT=*) "                  is upper than model top    (",XZHAT(IKU),")"
  STOP
END IF
IF (XALZBOT>=XZHAT(IKU-4) .AND. LVE_RELAX) THEN
  WRITE(ILUOUT,FMT=*) "INI_MODEL_n WARNING: you want to use vertical relaxation"
  WRITE(ILUOUT,FMT=*) "                    but the layer defined by XALZBOT(",XALZBOT,")"
  WRITE(ILUOUT,FMT=*) "                    contains less than 5 model levels"
END IF
DEALLOCATE(XZHAT)
!
!*       2.3  Compute sizes of arrays of the extended sub-domain
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IIU_ll=NIMAX_ll + 2 * JPHEXT
IJU_ll=NJMAX_ll + 2 * JPHEXT
! initialize NIMAX and NJMAX for not updated versions regarding the parallelism
! spawning,...
CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
!
NRR=0
NRRL=0
NRRI=0
IF ((CGETRVM /= 'SKIP' ).OR.(CGETRVT /= 'SKIP' )) THEN
  NRR = NRR+1
END IF
IF ((CGETRCM /= 'SKIP' ).OR.(CGETRCT /= 'SKIP' )) THEN
  NRR = NRR+1
  NRRL = NRRL+1
END IF
IF ((CGETRRM /= 'SKIP').OR.(CGETRRT /= 'SKIP' )) THEN
  NRR = NRR+1
  NRRL = NRRL+1
END IF
IF ((CGETRIM /= 'SKIP').OR.(CGETRIT /= 'SKIP' )) THEN
  NRR = NRR+1
  NRRI = NRRI+1
END IF
IF ((CGETRSM /= 'SKIP').OR.(CGETRST /= 'SKIP' )) THEN
  NRR = NRR+1
  NRRI = NRRI+1
END IF
IF ((CGETRGM /= 'SKIP').OR.(CGETRGT /= 'SKIP' )) THEN
  NRR = NRR+1
  NRRI = NRRI+1
END IF
IF ((CGETRHM /= 'SKIP').OR.(CGETRHT /= 'SKIP' )) THEN
  NRR = NRR+1
  NRRI = NRRI+1
END IF
IF (NVERB >= 5) THEN
  WRITE (UNIT=ILUOUT,FMT='("THERE ARE ",I2," WATER VARIABLES")') NRR
  WRITE (UNIT=ILUOUT,FMT='("THERE ARE ",I2," LIQUID VARIABLES")') NRRL
  WRITE (UNIT=ILUOUT,FMT='("THERE ARE ",I2," SOLID VARIABLES")') NRRI
END IF
!
!*       2.3  Update NSV and floating indices for the current model
!
! 
CALL UPDATE_NSV(KMI) 
!
!-------------------------------------------------------------------------------
!
!*       3.    ALLOCATE  MEMORY
!              -----------------
!
!*       3.1   Module MODD_FIELD_n
!
ALLOCATE(XUM(IIU,IJU,IKU))      ; XUM  = 0.0
ALLOCATE(XVM(IIU,IJU,IKU))      ; XVM  = 0.0
ALLOCATE(XWM(IIU,IJU,IKU))      ; XWM  = 0.0
ALLOCATE(XTHM(IIU,IJU,IKU))     ; XTHM = 0.0
ALLOCATE(XUT(IIU,IJU,IKU))      ; XUT  = 0.0
ALLOCATE(XVT(IIU,IJU,IKU))      ; XVT  = 0.0
ALLOCATE(XWT(IIU,IJU,IKU))      ; XWT  = 0.0
ALLOCATE(XTHT(IIU,IJU,IKU))     ; XTHT = 0.0
ALLOCATE(XRUS(IIU,IJU,IKU))     ; XRUS = 0.0
ALLOCATE(XRVS(IIU,IJU,IKU))     ; XRVS = 0.0
ALLOCATE(XRWS(IIU,IJU,IKU))     ; XRWS = 0.0
ALLOCATE(XRTHS(IIU,IJU,IKU))
IF (CTURB /= 'NONE') THEN
  ALLOCATE(XTKEM(IIU,IJU,IKU))
  ALLOCATE(XTKET(IIU,IJU,IKU))
  ALLOCATE(XRTKES(IIU,IJU,IKU))
ELSE
  ALLOCATE(XTKEM(0,0,0))
  ALLOCATE(XTKET(0,0,0))
  ALLOCATE(XRTKES(0,0,0))
END IF
IF (CTOM == 'TM06') THEN
  ALLOCATE(XBL_DEPTH(IIU,IJU))
ELSE
  ALLOCATE(XBL_DEPTH(0,0))
END IF
IF (LRMC01) THEN
  ALLOCATE(XSBL_DEPTH(IIU,IJU))
ELSE
  ALLOCATE(XSBL_DEPTH(0,0))
END IF
!
ALLOCATE(XPABSM(IIU,IJU,IKU)) ; XPABSM = 0.0
ALLOCATE(XPABST(IIU,IJU,IKU)) ; XPABST = 0.0
!
ALLOCATE(XRM(IIU,IJU,IKU,NRR))
ALLOCATE(XRT(IIU,IJU,IKU,NRR))
ALLOCATE(XRRS(IIU,IJU,IKU,NRR))
!
IF (CTURB /= 'NONE' .AND. NRR>1) THEN
  ALLOCATE(XSRCM(IIU,IJU,IKU))
  ALLOCATE(XSRCT(IIU,IJU,IKU))
  ALLOCATE(XSIGS(IIU,IJU,IKU))
ELSE
  ALLOCATE(XSRCM(0,0,0))
  ALLOCATE(XSRCT(0,0,0))
  ALLOCATE(XSIGS(0,0,0))
END IF
!
IF (NRR>1) THEN
  ALLOCATE(XCLDFR(IIU,IJU,IKU))
ELSE
  ALLOCATE(XCLDFR(0,0,0))
END IF
!
ALLOCATE(XSVM(IIU,IJU,IKU,NSV))
ALLOCATE(XSVT(IIU,IJU,IKU,NSV))
ALLOCATE(XRSVS(IIU,IJU,IKU,NSV))
!
!
!*       3.2   Module MODD_ADV_n      
!
IF (CMET_ADV_SCHEME(1:3) == 'PPM' .OR. CSV_ADV_SCHEME(1:3) == 'PPM') THEN
  ALLOCATE(XRTHMS(IIU,IJU,IKU)) ; XRTHMS = 0.
  IF (CTURB /= 'NONE') THEN
    ALLOCATE(XRTKEMS(IIU,IJU,IKU))
    XRTKEMS = 0.
  ENDIF
  ALLOCATE(XRRMS(IIU,IJU,IKU,NRR)); XRRMS = 0.
  ALLOCATE(XRSVMS(IIU,IJU,IKU,NSV)) ; XRSVMS = 0.
END IF
   
!
!*       3.3   Module MODD_GRID_n and MODD_METRICS_n
!
IF (LCARTESIAN) THEN
  ALLOCATE(XLON(0,0))
  ALLOCATE(XLAT(0,0))
  ALLOCATE(XMAP(0,0))
ELSE
  ALLOCATE(XLON(IIU,IJU))
  ALLOCATE(XLAT(IIU,IJU))
  ALLOCATE(XMAP(IIU,IJU))
END IF
ALLOCATE(XXHAT(IIU))
ALLOCATE(XDXHAT(IIU))
ALLOCATE(XYHAT(IJU))
ALLOCATE(XDYHAT(IJU))
ALLOCATE(XZS(IIU,IJU))
ALLOCATE(XZSMT(IIU,IJU))
ALLOCATE(XZZ(IIU,IJU,IKU))
ALLOCATE(XZHAT(IKU))
ALLOCATE(XDIRCOSZW(IIU,IJU))
ALLOCATE(XDIRCOSXW(IIU,IJU))
ALLOCATE(XDIRCOSYW(IIU,IJU))
ALLOCATE(XCOSSLOPE(IIU,IJU))
ALLOCATE(XSINSLOPE(IIU,IJU))
!
ALLOCATE(XDXX(IIU,IJU,IKU))
ALLOCATE(XDYY(IIU,IJU,IKU))
ALLOCATE(XDZX(IIU,IJU,IKU))
ALLOCATE(XDZY(IIU,IJU,IKU))
ALLOCATE(XDZZ(IIU,IJU,IKU))
!
!*       3.4   Modules MODD_REF and  MODD_REF_n
!
IF (KMI == 1) THEN
  ALLOCATE(XRHODREFZ(IKU),XTHVREFZ(IKU))
END IF
ALLOCATE(XRHODREF(IIU,IJU,IKU))
ALLOCATE(XTHVREF(IIU,IJU,IKU))
ALLOCATE(XEXNREF(IIU,IJU,IKU))
ALLOCATE(XRHODJ(IIU,IJU,IKU))
IF (CEQNSYS=='DUR' .AND. LUSERV) THEN
  ALLOCATE(XRVREF(IIU,IJU,IKU))
ELSE
  ALLOCATE(XRVREF(0,0,0))
END IF
!
!*       3.5   Module MODD_CURVCOR_n
!
IF (LTHINSHELL) THEN
  ALLOCATE(XCORIOX(0,0))
  ALLOCATE(XCORIOY(0,0))
ELSE
  ALLOCATE(XCORIOX(IIU,IJU))
  ALLOCATE(XCORIOY(IIU,IJU))
END IF
  ALLOCATE(XCORIOZ(IIU,IJU))
IF (LCARTESIAN) THEN
  ALLOCATE(XCURVX(0,0))
  ALLOCATE(XCURVY(0,0))
ELSE
  ALLOCATE(XCURVX(IIU,IJU))
  ALLOCATE(XCURVY(IIU,IJU))
END IF
!
!*       3.6   Module MODD_DYN_n
!
CALL GET_DIM_EXT_ll('Y',IIY,IJY)
IF (L2D) THEN
  ALLOCATE(XBFY(IIY,IJY,IKU))
ELSE
  ALLOCATE(XBFY(IJY,IIY,IKU)) ! transposition needed by the optimisition of the
                              ! FFT solver
END IF
!JUAN Z_SPLITING
CALL GET_DIM_EXT_ll('B',IIU_B,IJU_B)
ALLOCATE(XBFB(IIU_B,IJU_B,IKU))
CALL GET_DIM_EXTZ_ll('SXP2_YP1_Z',IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll)
ALLOCATE(XBF_SXP2_YP1_Z(IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll))
!JUAN Z_SPLITING
ALLOCATE(XAF(IKU),XCF(IKU))
ALLOCATE(XTRIGSX(3*IIU_ll))
ALLOCATE(XTRIGSY(3*IJU_ll))
ALLOCATE(XRHOM(IKU))
ALLOCATE(XALK(IKU))
ALLOCATE(XALKW(IKU))
IF ( LHORELAX_UVWTH .OR. LHORELAX_RV .OR.                                  &
     LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI  .OR. LHORELAX_RS  .OR. &
     LHORELAX_RG .OR. LHORELAX_RH .OR. LHORELAX_TKE .OR.                   &
     ANY(LHORELAX_SV) ) THEN
  ALLOCATE(XKURELAX(IIU,IJU))
  ALLOCATE(XKVRELAX(IIU,IJU))
  ALLOCATE(XKWRELAX(IIU,IJU))
  ALLOCATE(LMASK_RELAX(IIU,IJU))
ELSE
  ALLOCATE(XKURELAX(0,0))
  ALLOCATE(XKVRELAX(0,0))
  ALLOCATE(XKWRELAX(0,0))
  ALLOCATE(LMASK_RELAX(0,0))
END IF
!
! Additional fields for truly horizontal diffusion (Module MODD_DYNZD$n)
IF (LZDIFFU) THEN
  CALL INIT_TYPE_ZDIFFU_HALO2(XZDIFFU_HALO2)
ELSE
  CALL INIT_TYPE_ZDIFFU_HALO2(XZDIFFU_HALO2,0)
ENDIF
!
!*       3.7   Larger Scale variables (Module MODD_LSFIELD$n)
!
!
! upper relaxation part
!
ALLOCATE(XLSUM(IIU,IJU,IKU))    ; XLSUM  = 0.0
ALLOCATE(XLSVM(IIU,IJU,IKU))    ; XLSVM  = 0.0
ALLOCATE(XLSWM(IIU,IJU,IKU))    ; XLSWM  = 0.0
ALLOCATE(XLSTHM(IIU,IJU,IKU))   ; XLSTHM = 0.0
IF ( NRR > 0 ) THEN
  ALLOCATE(XLSRVM(IIU,IJU,IKU)) ; XLSRVM = 0.0
ELSE
  ALLOCATE(XLSRVM(0,0,0))
END IF
!
!  lbc part
!
IF ( L1D) THEN                         ! 1D case
!
  NSIZELBX_ll=0
  NSIZELBXU_ll=0
  NSIZELBY_ll=0
  NSIZELBYV_ll=0
  NSIZELBXTKE_ll=0
  NSIZELBXR_ll=0
  NSIZELBXSV_ll=0
  NSIZELBYTKE_ll=0
  NSIZELBYR_ll=0
  NSIZELBYSV_ll=0
  ALLOCATE(XLBXUM(0,0,0))
  ALLOCATE(XLBYUM(0,0,0))
  ALLOCATE(XLBXVM(0,0,0))
  ALLOCATE(XLBYVM(0,0,0))
  ALLOCATE(XLBXWM(0,0,0))
  ALLOCATE(XLBYWM(0,0,0))
  ALLOCATE(XLBXTHM(0,0,0))
  ALLOCATE(XLBYTHM(0,0,0))
  ALLOCATE(XLBXTKEM(0,0,0))
  ALLOCATE(XLBYTKEM(0,0,0))
  ALLOCATE(XLBXRM(0,0,0,0))
  ALLOCATE(XLBYRM(0,0,0,0))
  ALLOCATE(XLBXSVM(0,0,0,0))
  ALLOCATE(XLBYSVM(0,0,0,0))
!
ELSEIF( L2D ) THEN                         ! 2D case
!
  NSIZELBY_ll=0
  NSIZELBYV_ll=0
  NSIZELBYTKE_ll=0
  NSIZELBYR_ll=0
  NSIZELBYSV_ll=0
  ALLOCATE(XLBYUM(0,0,0))
  ALLOCATE(XLBYVM(0,0,0))
  ALLOCATE(XLBYWM(0,0,0))
  ALLOCATE(XLBYTHM(0,0,0))
  ALLOCATE(XLBYTKEM(0,0,0))
  ALLOCATE(XLBYRM(0,0,0,0))
  ALLOCATE(XLBYSVM(0,0,0,0))
!
  CALL GET_SIZEX_LB(HLUOUT,NIMAX_ll,NJMAX_ll,NRIMX,   &
       IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU,         &
       IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
!
  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2
    NSIZELBXU_ll=2*NRIMX+2
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,IKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,IKU))
  ELSE
    NSIZELBX_ll=2
    NSIZELBXU_ll=4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,IKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,IKU))
  END IF
!
  IF (CTURB /= 'NONE') THEN
    IF ( LHORELAX_TKE) THEN
      NSIZELBXTKE_ll=2* NRIMX+2
      ALLOCATE(XLBXTKEM(IISIZEXF,IJSIZEXF,IKU))
    ELSE
      NSIZELBXTKE_ll=2
      ALLOCATE(XLBXTKEM(IISIZEX2,IJSIZEX2,IKU))
    END IF
  ELSE
    NSIZELBXTKE_ll=0
    ALLOCATE(XLBXTKEM(0,0,0))
  END IF
  !
  IF ( NRR > 0 ) THEN
    IF (LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
         .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH               &
       ) THEN
      NSIZELBXR_ll=2* NRIMX+2
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,IKU,NRR))
    ELSE
      NSIZELBXR_ll=2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,IKU,NRR))
    ENDIF
  ELSE
    NSIZELBXR_ll=0
    ALLOCATE(XLBXRM(0,0,0,0))
  END IF
  !
  IF ( NSV > 0 ) THEN
    IF ( ANY( LHORELAX_SV(:)) ) THEN
      NSIZELBXSV_ll=2* NRIMX+2
      ALLOCATE(XLBXSVM(IISIZEXF,IJSIZEXF,IKU,NSV))
    ELSE
      NSIZELBXSV_ll=2
      ALLOCATE(XLBXSVM(IISIZEX2,IJSIZEX2,IKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
  END IF
!
ELSE                                   ! 3D case
!
!
  CALL GET_SIZEX_LB(HLUOUT,NIMAX_ll,NJMAX_ll,NRIMX,   &
       IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU,          &
       IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  CALL GET_SIZEY_LB(HLUOUT,NIMAX_ll,NJMAX_ll,NRIMY,   &
       IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV,          &
       IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)
!
  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2
    NSIZELBXU_ll=2*NRIMX+2
    NSIZELBY_ll=2*NRIMY+2
    NSIZELBYV_ll=2*NRIMY+2
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,IKU))
    ALLOCATE(XLBYUM(IISIZEYF,IJSIZEYF,IKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYVM(IISIZEYFV,IJSIZEYFV,IKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYWM(IISIZEYF,IJSIZEYF,IKU))
    ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYTHM(IISIZEYF,IJSIZEYF,IKU))
  ELSE
    NSIZELBX_ll=2
    NSIZELBXU_ll=4
    NSIZELBY_ll=2
    NSIZELBYV_ll=4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,IKU))
    ALLOCATE(XLBYUM(IISIZEY2,IJSIZEY2,IKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYVM(IISIZEY4,IJSIZEY4,IKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYWM(IISIZEY2,IJSIZEY2,IKU))
    ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYTHM(IISIZEY2,IJSIZEY2,IKU))
  END IF
  !
  IF (CTURB /= 'NONE') THEN
    IF ( LHORELAX_TKE) THEN
      NSIZELBXTKE_ll=2*NRIMX+2
      NSIZELBYTKE_ll=2*NRIMY+2
      ALLOCATE(XLBXTKEM(IISIZEXF,IJSIZEXF,IKU))
      ALLOCATE(XLBYTKEM(IISIZEYF,IJSIZEYF,IKU))
    ELSE
      NSIZELBXTKE_ll=2
      NSIZELBYTKE_ll=2
      ALLOCATE(XLBXTKEM(IISIZEX2,IJSIZEX2,IKU))
      ALLOCATE(XLBYTKEM(IISIZEY2,IJSIZEY2,IKU))
    END IF
  ELSE
    NSIZELBXTKE_ll=0
    NSIZELBYTKE_ll=0
    ALLOCATE(XLBXTKEM(0,0,0))
    ALLOCATE(XLBYTKEM(0,0,0))
  END IF
  !
  IF ( NRR > 0 ) THEN
    IF (LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
          .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH              &
       ) THEN
      NSIZELBXR_ll=2*NRIMX+2
      NSIZELBYR_ll=2*NRIMY+2
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,IKU,NRR))
      ALLOCATE(XLBYRM(IISIZEYF,IJSIZEYF,IKU,NRR))
    ELSE
      NSIZELBXR_ll=2
      NSIZELBYR_ll=2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,IKU,NRR))
      ALLOCATE(XLBYRM(IISIZEY2,IJSIZEY2,IKU,NRR))
    ENDIF
  ELSE
    NSIZELBXR_ll=0
    NSIZELBYR_ll=0
    ALLOCATE(XLBXRM(0,0,0,0))
    ALLOCATE(XLBYRM(0,0,0,0))
  END IF
  !
  IF ( NSV > 0 ) THEN
    IF ( ANY( LHORELAX_SV(:)) ) THEN
      NSIZELBXSV_ll=2*NRIMX+2
      NSIZELBYSV_ll=2*NRIMY+2
      ALLOCATE(XLBXSVM(IISIZEXF,IJSIZEXF,IKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEYF,IJSIZEYF,IKU,NSV))
    ELSE
      NSIZELBXSV_ll=2
      NSIZELBYSV_ll=2
      ALLOCATE(XLBXSVM(IISIZEX2,IJSIZEX2,IKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEY2,IJSIZEY2,IKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    NSIZELBYSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
    ALLOCATE(XLBYSVM(0,0,0,0))
  END IF
END IF      ! END OF THE IF STRUCTURE ON THE MODEL DIMENSION
!
!
IF ( KMI > 1 ) THEN 
  ! it has been assumed that the THeta field used the largest rim area compared
  ! to the others prognostic variables, if it is not the case, you must change
  ! these lines
  ALLOCATE(XCOEFLIN_LBXM(SIZE(XLBXTHM,1),SIZE(XLBXTHM,2),SIZE(XLBXTHM,3)))
  ALLOCATE(   NKLIN_LBXM(SIZE(XLBXTHM,1),SIZE(XLBXTHM,2),SIZE(XLBXTHM,3)))
  ALLOCATE(XCOEFLIN_LBYM(SIZE(XLBYTHM,1),SIZE(XLBYTHM,2),SIZE(XLBYTHM,3)))
  ALLOCATE(   NKLIN_LBYM(SIZE(XLBYTHM,1),SIZE(XLBYTHM,2),SIZE(XLBYTHM,3)))
  ALLOCATE(XCOEFLIN_LBXU(SIZE(XLBXUM,1),SIZE(XLBXUM,2),SIZE(XLBXUM,3)))
  ALLOCATE(   NKLIN_LBXU(SIZE(XLBXUM,1),SIZE(XLBXUM,2),SIZE(XLBXUM,3)))
  ALLOCATE(XCOEFLIN_LBYU(SIZE(XLBYUM,1),SIZE(XLBYUM,2),SIZE(XLBYUM,3)))
  ALLOCATE(   NKLIN_LBYU(SIZE(XLBYUM,1),SIZE(XLBYUM,2),SIZE(XLBYUM,3)))
  ALLOCATE(XCOEFLIN_LBXV(SIZE(XLBXVM,1),SIZE(XLBXVM,2),SIZE(XLBXVM,3)))
  ALLOCATE(   NKLIN_LBXV(SIZE(XLBXVM,1),SIZE(XLBXVM,2),SIZE(XLBXVM,3)))
  ALLOCATE(XCOEFLIN_LBYV(SIZE(XLBYVM,1),SIZE(XLBYVM,2),SIZE(XLBYVM,3)))
  ALLOCATE(   NKLIN_LBYV(SIZE(XLBYVM,1),SIZE(XLBYVM,2),SIZE(XLBYVM,3)))
  ALLOCATE(XCOEFLIN_LBXW(SIZE(XLBXWM,1),SIZE(XLBXWM,2),SIZE(XLBXWM,3)))
  ALLOCATE(   NKLIN_LBXW(SIZE(XLBXWM,1),SIZE(XLBXWM,2),SIZE(XLBXWM,3)))
  ALLOCATE(XCOEFLIN_LBYW(SIZE(XLBYWM,1),SIZE(XLBYWM,2),SIZE(XLBYWM,3)))
  ALLOCATE(   NKLIN_LBYW(SIZE(XLBYWM,1),SIZE(XLBYWM,2),SIZE(XLBYWM,3)))
END IF
!
!  allocation of the LS fields for vertical relaxation and numerical diffusion
IF( .NOT. LSTEADYLS )  THEN
!
  ALLOCATE(XLSUS(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
  ALLOCATE(XLSVS(SIZE(XLSVM,1),SIZE(XLSVM,2),SIZE(XLSVM,3)))
  ALLOCATE(XLSWS(SIZE(XLSWM,1),SIZE(XLSWM,2),SIZE(XLSWM,3)))
  ALLOCATE(XLSTHS(SIZE(XLSTHM,1),SIZE(XLSTHM,2),SIZE(XLSTHM,3)))
  ALLOCATE(XLSRVS(SIZE(XLSRVM,1),SIZE(XLSRVM,2),SIZE(XLSRVM,3)))
!
ELSE
!
  ALLOCATE(XLSUS(0,0,0))
  ALLOCATE(XLSVS(0,0,0))
  ALLOCATE(XLSWS(0,0,0))
  ALLOCATE(XLSTHS(0,0,0))
  ALLOCATE(XLSRVS(0,0,0))
!
END IF
!  allocation of the LB fields for horizontal relaxation and Lateral Boundaries
IF( .NOT. ( LSTEADYLS .AND. KMI==1 ) )  THEN
!
  ALLOCATE(XLBXTKES(SIZE(XLBXTKEM,1),SIZE(XLBXTKEM,2),SIZE(XLBXTKEM,3)))
  ALLOCATE(XLBYTKES(SIZE(XLBYTKEM,1),SIZE(XLBYTKEM,2),SIZE(XLBYTKEM,3)))
  ALLOCATE(XLBXUS(SIZE(XLBXUM,1),SIZE(XLBXUM,2),SIZE(XLBXUM,3)))
  ALLOCATE(XLBYUS(SIZE(XLBYUM,1),SIZE(XLBYUM,2),SIZE(XLBYUM,3)))
  ALLOCATE(XLBXVS(SIZE(XLBXVM,1),SIZE(XLBXVM,2),SIZE(XLBXVM,3)))
  ALLOCATE(XLBYVS(SIZE(XLBYVM,1),SIZE(XLBYVM,2),SIZE(XLBYVM,3)))
  ALLOCATE(XLBXWS(SIZE(XLBXWM,1),SIZE(XLBXWM,2),SIZE(XLBXWM,3)))
  ALLOCATE(XLBYWS(SIZE(XLBYWM,1),SIZE(XLBYWM,2),SIZE(XLBYWM,3)))
  ALLOCATE(XLBXTHS(SIZE(XLBXTHM,1),SIZE(XLBXTHM,2),SIZE(XLBXTHM,3)))
  ALLOCATE(XLBYTHS(SIZE(XLBYTHM,1),SIZE(XLBYTHM,2),SIZE(XLBYTHM,3)))
  ALLOCATE(XLBXRS(SIZE(XLBXRM,1),SIZE(XLBXRM,2),SIZE(XLBXRM,3),SIZE(XLBXRM,4)))
  ALLOCATE(XLBYRS(SIZE(XLBYRM,1),SIZE(XLBYRM,2),SIZE(XLBYRM,3),SIZE(XLBYRM,4)))
  ALLOCATE(XLBXSVS(SIZE(XLBXSVM,1),SIZE(XLBXSVM,2),SIZE(XLBXSVM,3),SIZE(XLBXSVM,4)))
  ALLOCATE(XLBYSVS(SIZE(XLBYSVM,1),SIZE(XLBYSVM,2),SIZE(XLBYSVM,3),SIZE(XLBYSVM,4)))
!
ELSE
!
  ALLOCATE(XLBXTKES(0,0,0))
  ALLOCATE(XLBYTKES(0,0,0))
  ALLOCATE(XLBXUS(0,0,0))
  ALLOCATE(XLBYUS(0,0,0))
  ALLOCATE(XLBXVS(0,0,0))
  ALLOCATE(XLBYVS(0,0,0))
  ALLOCATE(XLBXWS(0,0,0))
  ALLOCATE(XLBYWS(0,0,0))
  ALLOCATE(XLBXTHS(0,0,0))
  ALLOCATE(XLBYTHS(0,0,0))
  ALLOCATE(XLBXRS(0,0,0,0))
  ALLOCATE(XLBYRS(0,0,0,0))
  ALLOCATE(XLBXSVS(0,0,0,0))
  ALLOCATE(XLBYSVS(0,0,0,0))
!
END IF
!
!
!*       4.8   Module MODD_RADIATIONS_n (except XOZON and XAER)
!
!

NSWB_MNH = 6
ALLOCATE(XSW_BANDS (NSWB_MNH))
ALLOCATE(XZENITH   (IIU,IJU))
ALLOCATE(XAZIM     (IIU,IJU))
ALLOCATE(XALBUV    (IIU,IJU))
ALLOCATE(XDIRSRFSWD(IIU,IJU,NSWB_MNH))
ALLOCATE(XSCAFLASWD(IIU,IJU,NSWB_MNH))
ALLOCATE(XFLALWD   (IIU,IJU))
!
IF (CRAD /= 'NONE') THEN
  ALLOCATE(XSLOPANG(IIU,IJU))
  ALLOCATE(XSLOPAZI(IIU,IJU))
  ALLOCATE(XDTHRAD(IIU,IJU,IKU))
  ALLOCATE(XDIRFLASWD(IIU,IJU,NSWB_MNH))
  ALLOCATE(XDIR_ALB(IIU,IJU,NSWB_MNH))
  ALLOCATE(XSCA_ALB(IIU,IJU,NSWB_MNH))
  ALLOCATE(XEMIS  (IIU,IJU))
  ALLOCATE(XTSRAD (IIU,IJU))    ; XTSRAD = 0.0
  ALLOCATE(XSEA (IIU,IJU))
  ALLOCATE(XZS_XY (IIU,IJU))
  ALLOCATE(NCLEARCOL_TM1(IIU,IJU))
ELSE
  ALLOCATE(XSLOPANG(0,0))
  ALLOCATE(XSLOPAZI(0,0))
  ALLOCATE(XDTHRAD(0,0,0))
  ALLOCATE(XDIRFLASWD(0,0,0))
  ALLOCATE(XDIR_ALB(0,0,0))
  ALLOCATE(XSCA_ALB(0,0,0))
  ALLOCATE(XEMIS  (0,0))
  ALLOCATE(XTSRAD (0,0))
  ALLOCATE(XSEA (0,0))
  ALLOCATE(XZS_XY (0,0))
  ALLOCATE(NCLEARCOL_TM1(0,0))
END IF

IF (CRAD == 'ECMW') THEN
  ALLOCATE(XSTROATM(31,6))
  ALLOCATE(XSMLSATM(31,6))
  ALLOCATE(XSMLWATM(31,6))
  ALLOCATE(XSPOSATM(31,6))
  ALLOCATE(XSPOWATM(31,6))
  ALLOCATE(XSTATM(31,6))
ELSE
  ALLOCATE(XSTROATM(0,0))
  ALLOCATE(XSMLSATM(0,0))
  ALLOCATE(XSMLWATM(0,0))
  ALLOCATE(XSPOSATM(0,0))
  ALLOCATE(XSPOWATM(0,0))
  ALLOCATE(XSTATM(0,0))
END IF
!
!*       4.9   Module MODD_DEEP_CONVECTION_n
!
IF (CDCONV /= 'NONE') THEN
  ALLOCATE(NCOUNTCONV(IIU,IJU))
  ALLOCATE(XDTHCONV(IIU,IJU,IKU))
  ALLOCATE(XDRVCONV(IIU,IJU,IKU))
  ALLOCATE(XDRCCONV(IIU,IJU,IKU))
  ALLOCATE(XDRICONV(IIU,IJU,IKU))
  ALLOCATE(XPRCONV(IIU,IJU))
  ALLOCATE(XPACCONV(IIU,IJU))
  ALLOCATE(XPRSCONV(IIU,IJU))
  ! diagnostics
  IF (LCH_CONV_LINOX) THEN
    ALLOCATE(XIC_RATE(IIU,IJU))
    ALLOCATE(XCG_RATE(IIU,IJU))
    ALLOCATE(XIC_TOTAL_NUMBER(IIU,IJU))
    ALLOCATE(XCG_TOTAL_NUMBER(IIU,IJU))
  ELSE
    ALLOCATE(XIC_RATE(0,0))
    ALLOCATE(XCG_RATE(0,0))
    ALLOCATE(XIC_TOTAL_NUMBER(0,0))
    ALLOCATE(XCG_TOTAL_NUMBER(0,0))
  END IF
  IF ( LDIAGCONV )  THEN
    ALLOCATE(XUMFCONV(IIU,IJU,IKU))
    ALLOCATE(XDMFCONV(IIU,IJU,IKU))
    ALLOCATE(XPRLFLXCONV(IIU,IJU,IKU))
    ALLOCATE(XPRSFLXCONV(IIU,IJU,IKU))
    ALLOCATE(XCAPE(IIU,IJU))
    ALLOCATE(NCLTOPCONV(IIU,IJU))
    ALLOCATE(NCLBASCONV(IIU,IJU))
  ELSE
    ALLOCATE(XUMFCONV(0,0,0))
    ALLOCATE(XDMFCONV(0,0,0))
    ALLOCATE(XPRLFLXCONV(0,0,0))
    ALLOCATE(XPRSFLXCONV(0,0,0))
    ALLOCATE(XCAPE(0,0))
    ALLOCATE(NCLTOPCONV(0,0))
    ALLOCATE(NCLBASCONV(0,0))
  END IF
ELSE
  ALLOCATE(XPRCONV(0,0))
  ALLOCATE(XPACCONV(0,0))
  ALLOCATE(XPRSCONV(0,0))
END IF
!
IF (CDCONV == 'KAFR' .AND. LSUBG_COND .AND. LSIG_CONV) THEN
  ALLOCATE(XMFCONV(IIU,IJU,IKU))
ELSE
  ALLOCATE(XMFCONV(0,0,0))
ENDIF
!
IF (CDCONV == 'KAFR' .AND. LCHTRANS .AND. NSV > 0 ) THEN
  ALLOCATE(XDSVCONV(IIU,IJU,IKU,NSV))
ELSE
  ALLOCATE(XDSVCONV(0,0,0,0))
END IF
!
!*       4.10  Local variables
!
ALLOCATE(ZJ(IIU,IJU,IKU))
!
!*      4.11 Forcing variables (Module MODD_FRC)
!
IF (KMI == 1) THEN
  IF ( LFORCING ) THEN
    ALLOCATE(TDTFRC(NFRC))
    ALLOCATE(XUFRC(IKU,NFRC))
    ALLOCATE(XVFRC(IKU,NFRC))
    ALLOCATE(XWFRC(IKU,NFRC))
    ALLOCATE(XTHFRC(IKU,NFRC))
    ALLOCATE(XRVFRC(IKU,NFRC))
    ALLOCATE(XTENDTHFRC(IKU,NFRC))
    ALLOCATE(XTENDRVFRC(IKU,NFRC))
    ALLOCATE(XGXTHFRC(IKU,NFRC))
    ALLOCATE(XGYTHFRC(IKU,NFRC))
    ALLOCATE(XPGROUNDFRC(NFRC))
  ELSE
    ALLOCATE(TDTFRC(0))
    ALLOCATE(XUFRC(0,0))
    ALLOCATE(XVFRC(0,0))
    ALLOCATE(XWFRC(0,0))
    ALLOCATE(XTHFRC(0,0))
    ALLOCATE(XRVFRC(0,0))
    ALLOCATE(XTENDTHFRC(0,0))
    ALLOCATE(XTENDRVFRC(0,0))
    ALLOCATE(XGXTHFRC(0,0))
    ALLOCATE(XGYTHFRC(0,0))
    ALLOCATE(XPGROUNDFRC(0))
  END IF
  IF ( LFORCING ) THEN
    ALLOCATE(XWTFRC(IIU,IJU,IKU))
  ELSE
    ALLOCATE(XWTFRC(0,0,0))
  END IF
END IF
!
!*      4.12   Module MODD_ICE_CONC_n
!
IF (     (CCLOUD == 'ICE3'.OR.CCLOUD == 'ICE4') .AND.   &
     (CPROGRAM == 'DIAG  '.OR.CPROGRAM == 'MESONH')) THEN
  ALLOCATE(XCIT(IIU,IJU,IKU))
ELSE
  ALLOCATE(XCIT(0,0,0))
END IF
!
!*      4.13   Module MODD_TURB_CLOUD
!
IF (KMI == NMODEL_CLOUD .AND. CTURBLEN_CLOUD/='NONE' ) THEN
  IF (ALLOCATED(XCEI)) THEN
    DEALLOCATE(XCEI)
    ALLOCATE(XCEI(IIU,IJU,IKU))
  END IF
ELSE
  IF (.NOT.(ALLOCATED(XCEI))) ALLOCATE(XCEI(0,0,0))
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       5.    INITIALIZE BUDGET VARIABLES
!              ---------------------------
!
IF ( CBUTYPE /= "NONE" .AND. NBUMOD == KMI ) THEN
  CALL INI_BUDGET(HLUOUT,XTSTEP,NSV,NRR,                       &
             LNUMDIFF,LHORELAX_UVWTH,LHORELAX_RV, LHORELAX_RC,LHORELAX_RR,    &
             LHORELAX_RI,LHORELAX_RS,LHORELAX_RG, LHORELAX_RH,LHORELAX_TKE,   &
             LHORELAX_SV,LVE_RELAX,LCHTRANS,LNUDGING,                         &
             CRAD,CDCONV,CTURB,CTURBDIM,CCLOUD,CMET_ADV_SCHEME,CSV_ADV_SCHEME )
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       6.    INITIALIZE INTERPOLATION COEFFICIENTS
!
CALL INI_BIKHARDT_n (NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),KMI)
!
!-------------------------------------------------------------------------------
!
!*       7.    INITIALIZE GRIDS AND METRIC COEFFICIENTS
!              ----------------------------------------
!
CALL SET_GRID(KMI,HINIFILE,HLUOUT,IIU,IJU,IKU,NIMAX_ll,NJMAX_ll,         &
              XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4,           &
              XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4,           &
              NXOR_ALL(KMI),NYOR_ALL(KMI),NXEND_ALL(KMI),NYEND_ALL(KMI), &
              NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),                       &
              CLBCX,CLBCY,                                               &
              XTSTEP,XSEGLEN,                                            &
              XLONORI,XLATORI,XLON,XLAT,                                 &
              XXHAT,XYHAT,XDXHAT,XDYHAT, XMAP,                           &
              XZS,XZZ,XZHAT,LSLEVE,XLEN1,XLEN2,XZSMT,                    &
              ZJ,                                                        &
              TDTMOD,TDTCUR,NSTOP,NOUT_TIMES,NOUT_NUMB)
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!* update halos of metric coefficients
!
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!
CALL SET_DIRCOS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,TZINITHALO2D_ll,   &
                XDIRCOSXW,XDIRCOSYW,XDIRCOSZW,XCOSSLOPE,XSINSLOPE  )
!
! grid nesting initializations
IF ( KMI == 1 ) THEN
  XTSTEP_MODEL1=XTSTEP
END IF
!
NDT_2_WAY(KMI)=4
!
!-------------------------------------------------------------------------------
!
!*      8.    INITIALIZE DATA FOR JVALUES AND AEROSOLS 
!
IF ( LUSECHEM .OR. LCHEMDIAG ) THEN
  IF ((KMI==1).AND.(CPROGRAM == "MESONH".OR.CPROGRAM == "DIAG  "))  &
    CALL CH_INIT_JVALUES(TDTCUR%TDATE%DAY, TDTCUR%TDATE%MONTH,      &
                                           TDTCUR%TDATE%YEAR, ILUOUT)
!
  IF (LORILAM) CALL CH_AER_MOD_INIT
END IF
!
!-------------------------------------------------------------------------------
!
!*       9.    INITIALIZE THE PROGNOSTIC FIELDS
!              --------------------------------
!
CALL READ_FIELD(HINIFILE,HLUOUT,IIU,IJU,IKU,PTSTEP_OLD,XTSTEP,                &
                CGETTKEM,CGETRVM,CGETRCM,CGETRRM,CGETRIM,CGETRSM,             &
                CGETRGM,CGETRHM,CGETSVM,CGETSRCM,                             &
                CGETTKET,CGETRVT,CGETRCT,CGETRRT,CGETRIT,CGETCIT,             &
                CGETRST,CGETRGT,CGETRHT,CGETSVT,CGETSRCT,CGETSIGS,CGETCLDFR,  &
                CGETBL_DEPTH,CGETSBL_DEPTH,                                   &
                NSIZELBX_ll,NSIZELBXU_ll,NSIZELBY_ll,NSIZELBYV_ll,            &
                NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                &
                NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,        &
                XUM,XVM,XWM,XTHM,XPABSM,XTKEM,XRM,XSVM,XSRCM,                 &
                XUT,XVT,XWT,XTHT,XPABST,XTKET,XRT,XSVT,XCIT,XDRYMASST,        &
                XSIGS,XSRCT,XCLDFR,XBL_DEPTH,XSBL_DEPTH,                      &
                XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,                              &
                XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,                        &
                XLBXRM,XLBXSVM,                                               &
                XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,                        &
                XLBYRM,XLBYSVM,                                               &
                NFRC,TDTFRC,XUFRC,XVFRC,XWFRC,XTHFRC,XRVFRC,                  &
                XTENDTHFRC,XTENDRVFRC,XGXTHFRC,XGYTHFRC,                      &
                XPGROUNDFRC                                                   )
!
!-------------------------------------------------------------------------------
!
!
!*       10.   INITIALIZE REFERENCE STATE
!              ---------------------------
!
!
CALL SET_REF(KMI,HINIFILE,HLUOUT,                                &
             XZZ,XZHAT,ZJ,XDXX,XDYY,CLBCX,CLBCY,                 &
             XREFMASS,XMASS_O_PHI0,XLINMASS,                      &
             XRHODREF,XTHVREF,XRVREF,XEXNREF,XRHODJ              )
!
!-------------------------------------------------------------------------------
!
!*       11.1    INITIALIZE THE TURBULENCE VARIABLES
!               -----------------------------------
!
IF ((CTURB == 'TKEL').AND.(CCONF=='START')) THEN
  CALL INI_TKE_EPS(CGETTKEM,CGETTKET,XTHVREF,XZZ, &
                   XUM,XVM,XTHM,                  &
                   XUT,XVT,XTHT,                  &
                   XTKEM,XTKET,TZINITHALO3D_ll    )
END IF
!
!
!*       11.2   INITIALIZE THE LES VARIABLES
!               ----------------------------
!
CALL INI_LES_n
!
!-------------------------------------------------------------------------------
!
!*       12.    INITIALIZE THE SOURCE OF TOTAL DRY MASS Md
!               ------------------------------------------
!
IF((KMI==1).AND.LSTEADYLS) THEN
  XDRYMASSS = 0.
END IF
!
!-------------------------------------------------------------------------------
!
!*       13.    INITIALIZE THE MICROPHYSICS                   
!               ----------------------------
!
CALL INI_MICRO_n
!
!-------------------------------------------------------------------------------
!
!*       14.   INITIALIZE THE LARGE SCALE SOURCES
!              ----------------------------------
!
IF ((KMI==1).AND.(.NOT. LSTEADYLS)) THEN
  CALL INI_CPL(HLUOUT,NSTOP,XTSTEP,LSTEADYLS,                                 &
               CGETTKEM,                                                      &
               CGETRVM,CGETRCM,CGETRRM,CGETRIM,                               &
               CGETRSM,CGETRGM,CGETRHM,CGETSVM,                               &
               NSV,NIMAX_ll,NJMAX_ll,                                         &
               NSIZELBX_ll,NSIZELBXU_ll,NSIZELBY_ll,NSIZELBYV_ll,             &
               NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                 &
               NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,         &
               XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XDRYMASST,                     &
               XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,          &
               XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,          &
               XLSUS,XLSVS,XLSWS,XLSTHS,XLSRVS,XDRYMASSS,                     &
               XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,          &
               XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS           )
END IF
!
IF ( KMI > 1) THEN
  ! Use dummy pointers to correct an ifort BUG
  DPTR_XBMX1=>XBMX1
  DPTR_XBMX2=>XBMX2
  DPTR_XBMX3=>XBMX3
  DPTR_XBMX4=>XBMX4
  DPTR_XBMY1=>XBMY1
  DPTR_XBMY2=>XBMY2
  DPTR_XBMY3=>XBMY3
  DPTR_XBMY4=>XBMY4
  DPTR_XBFX1=>XBFX1
  DPTR_XBFX2=>XBFX2
  DPTR_XBFX3=>XBFX3
  DPTR_XBFX4=>XBFX4
  DPTR_XBFY1=>XBFY1
  DPTR_XBFY2=>XBFY2
  DPTR_XBFY3=>XBFY3
  DPTR_XBFY4=>XBFY4
  DPTR_CLBCX=>CLBCX
  DPTR_CLBCY=>CLBCY
  !
  DPTR_XZZ=>XZZ
  DPTR_XZHAT=>XZHAT
  DPTR_XLSUM=>XLSUM
  DPTR_XLSVM=>XLSVM
  DPTR_XLSWM=>XLSWM
  DPTR_XLSTHM=>XLSTHM
  DPTR_XLSRVM=>XLSRVM
  DPTR_XLSUS=>XLSUS
  DPTR_XLSVS=>XLSVS
  DPTR_XLSWS=>XLSWS
  DPTR_XLSTHS=>XLSTHS
  DPTR_XLSRVS=>XLSRVS
  !
  DPTR_NKLIN_LBXU=>NKLIN_LBXU
  DPTR_XCOEFLIN_LBXU=>XCOEFLIN_LBXU
  DPTR_NKLIN_LBYU=>NKLIN_LBYU
  DPTR_XCOEFLIN_LBYU=>XCOEFLIN_LBYU
  DPTR_NKLIN_LBXV=>NKLIN_LBXV
  DPTR_XCOEFLIN_LBXV=>XCOEFLIN_LBXV
  DPTR_NKLIN_LBYV=>NKLIN_LBYV
  DPTR_XCOEFLIN_LBYV=>XCOEFLIN_LBYV
  DPTR_NKLIN_LBXW=>NKLIN_LBXW
  DPTR_XCOEFLIN_LBXW=>XCOEFLIN_LBXW
  DPTR_NKLIN_LBYW=>NKLIN_LBYW
  DPTR_XCOEFLIN_LBYW=>XCOEFLIN_LBYW
  DPTR_NKLIN_LBXM=>NKLIN_LBXM
  DPTR_XCOEFLIN_LBXM=>XCOEFLIN_LBXM
  DPTR_NKLIN_LBYM=>NKLIN_LBYM
  DPTR_XCOEFLIN_LBYM=>XCOEFLIN_LBYM
  !
  CALL INI_SPAWN_LS_n(NDAD(KMI),XTSTEP,KMI,                                 &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,      &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,      &
       NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),                  &
       DPTR_CLBCX,DPTR_CLBCY,DPTR_XZZ,DPTR_XZHAT,                                &
       LSLEVE,XLEN1,XLEN2,                                   &
       DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSTHM,DPTR_XLSRVM,                      &
       DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS,DPTR_XLSTHS,DPTR_XLSRVS,                      &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,    &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,    &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,    &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM     )
  !
  DPTR_XLBXUM=>XLBXUM
  DPTR_XLBYUM=>XLBYUM
  DPTR_XLBXVM=>XLBXVM
  DPTR_XLBYVM=>XLBYVM
  DPTR_XLBXWM=>XLBXWM
  DPTR_XLBYWM=>XLBYWM
  DPTR_XLBXTHM=>XLBXTHM
  DPTR_XLBYTHM=>XLBYTHM
  DPTR_XLBXTKEM=>XLBXTKEM
  DPTR_XLBYTKEM=>XLBYTKEM
  DPTR_XLBXRM=>XLBXRM
  DPTR_XLBYRM=>XLBYRM
  DPTR_XLBXSVM=>XLBXSVM
  DPTR_XLBYSVM=>XLBYSVM
  CALL INI_ONE_WAY_n(NDAD(KMI),CLUOUT,XTSTEP,KMI,1,                         &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
       NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),NDTRATIO(KMI),      &
       DPTR_CLBCX,DPTR_CLBCY,NRIMX,NRIMY,                                &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,      &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,      &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,      &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM,      &
       CCLOUD,                                                 &
       DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM,DPTR_XLBXWM,DPTR_XLBYWM,    &
       DPTR_XLBXTHM,DPTR_XLBYTHM,                                                  &
       DPTR_XLBXTKEM,DPTR_XLBYTKEM,                                                &
       DPTR_XLBXRM,DPTR_XLBYRM,DPTR_XLBXSVM,DPTR_XLBYSVM                           )
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       15.    INITIALIZE THE SCALAR VARIABLES
!               -------------------------------
!
IF (LLG .AND. LINIT_LG .AND. CPROGRAM=='MESONH') &
  CALL INI_LG(XXHAT,XYHAT,XZZ,XSVM,XSVT,XLBXSVM,XLBYSVM)

!
!*       16.    BUILT THE GENERIC OUTPUT NAME
!               ----------------------------
!
!!$COUTFILE=CEXP//'.$n.'//TRIM(ADJUSTL(CSEG))
!!$CFMDIAC =CEXP//'.$n.'//TRIM(ADJUSTL(CSEG))//'.000'
WRITE(COUTFILE,'(A,".",I1,".",A)') CEXP,KMI,TRIM(ADJUSTL(CSEG))
WRITE(CFMDIAC, '(A,".",I1,".",A)') CEXP,KMI,TRIM(ADJUSTL(CSEG))//'.000'
IF (CPROGRAM=='MESONH') THEN
  IF ( NDAD(KMI) == 1)  CDAD_NAME(KMI) = CEXP//'.1.'//CSEG
  IF ( NDAD(KMI) == 2)  CDAD_NAME(KMI) = CEXP//'.2.'//CSEG
  IF ( NDAD(KMI) == 3)  CDAD_NAME(KMI) = CEXP//'.3.'//CSEG
  IF ( NDAD(KMI) == 4)  CDAD_NAME(KMI) = CEXP//'.4.'//CSEG
  IF ( NDAD(KMI) == 5)  CDAD_NAME(KMI) = CEXP//'.5.'//CSEG
  IF ( NDAD(KMI) == 6)  CDAD_NAME(KMI) = CEXP//'.6.'//CSEG
  IF ( NDAD(KMI) == 7)  CDAD_NAME(KMI) = CEXP//'.7.'//CSEG
  IF ( NDAD(KMI) == 8)  CDAD_NAME(KMI) = CEXP//'.8.'//CSEG
END IF
!
!-------------------------------------------------------------------------------
!
!*       17.    INITIALIZE THE PARAMETERS FOR THE DYNAMICS
!               ------------------------------------------
!
CALL INI_DYNAMICS(HLUOUT,XLON,XLAT,XRHODJ,XTHVREF,XMAP,XZZ,XDXHAT,XDYHAT,     &
             XZHAT,CLBCX,CLBCY,XTSTEP,                                        &
             LVE_RELAX,LHORELAX_UVWTH,LHORELAX_RV,                            &
             LHORELAX_RC,LHORELAX_RR,LHORELAX_RI,LHORELAX_RS,LHORELAX_RG,     &
             LHORELAX_RH,LHORELAX_TKE,LHORELAX_SV,                            &
             LHORELAX_SVC2R2,LHORELAX_SVC1R3,LHORELAX_SVELEC,LHORELAX_SVLG,   &
             LHORELAX_SVCHEM,LHORELAX_SVAER,LHORELAX_SVDST,LHORELAX_SVSLT,    &
             XRIMKMAX,NRIMX,NRIMY,XALKTOP,XALZBOT,                            &
             XT4DIFF,                                                         &
             XCORIOX,XCORIOY,XCORIOZ,XCURVX,XCURVY,                           &
             XDXHATM,XDYHATM,XRHOM,XAF,XBFY,XCF,XTRIGSX,XTRIGSY,NIFAXX,NIFAXY,&
             XALK,XALKW,NALBOT,LMASK_RELAX,XKURELAX,XKVRELAX,XKWRELAX,        &
             XDK2,XDK4,LZDIFFU,XZDIFFU_HALO2,                                 &
             XBFB,&
             XBF_SXP2_YP1_Z ) !JUAN Z_SPLITING
!
!-------------------------------------------------------------------------------
!
!*      18.    SURFACE FIELDS
!              --------------
!
!*      18.1   Radiative setup
!              ---------------
!
IF (CRAD   /= 'NONE') THEN
  IF (CGETRAD =='INIT') THEN
    GINIRAD  =.TRUE.
  ELSE
    GINIRAD  =.FALSE.
  END IF
  CALL INI_RADIATIONS(HINIFILE,HLUOUT,GINIRAD,TDTCUR,XZZ, &
                      XDXX, XDYY,                         &
                      XSINDEL,XCOSDEL,XTSIDER,XCORSOL,    &
                      XSLOPANG,XSLOPAZI,                  &
                      XDTHRAD,XDIRFLASWD,XSCAFLASWD,      &
                      XFLALWD,XDIRSRFSWD,NCLEARCOL_TM1,   &
                      XZENITH,XAZIM,                      &
                      TDTRAD_FULL,TDTRAD_CLONLY,          &
                      TZINITHALO2D_ll                     )
  !
  IF (GINIRAD) CALL SUNPOS_n(XZENITH,PAZIMSOL=XAZIM)
  CALL SURF_SOLAR_GEOM    (XZS, XZS_XY)
ELSE
  XAZIM       = 0.
  XZENITH     = XPI/2.
  XDIRSRFSWD  = 0.
  XSCAFLASWD  = 0.
  XFLALWD     = 300.  ! W/m2
  XTSIDER     = 0.
END IF
!
!
CALL INI_SW_SETUP (CRAD,NSWB_MNH,XSW_BANDS)
!
!
!       18.1.1 Special initialisation for CO2 content
!              CO2 (molar mass=44) horizontally and vertically homogeneous at 360 ppm
!
XCCO2 = 360.0E-06 * 44.0E-03 / XMD
!
!
!*      18.2   Externalized surface fields
!              ---------------------------
!
ALLOCATE(ZCO2(IIU,IJU))
ZCO2(:,:) = XCCO2
!

ALLOCATE(ZDIR_ALB(IIU,IJU,NSWB_MNH))
ALLOCATE(ZSCA_ALB(IIU,IJU,NSWB_MNH))
ALLOCATE(ZEMIS  (IIU,IJU))
ALLOCATE(ZTSRAD (IIU,IJU))
!
CALL FMREAD(HINIFILE,'MASDEV',HLUOUT,'--',IMASDEV,IGRID,ILENCH,YCOMMENT,IRESP)
IF (IMASDEV>=46) THEN
  CALL FMREAD(HINIFILE,'SURF',HLUOUT,'--',CSURF,IGRID,ILENCH,YCOMMENT,IRESP)
ELSE
  CSURF = "EXTE"
END IF
!
!
IF (CSURF=='EXTE' .AND. (CPROGRAM=='MESONH' .OR. CPROGRAM=='DIAG  ')) THEN
  CALL GOTO_SURFEX(KMI)
  !* initialization of surface
  CALL INIT_GROUND_PARAM_n ('ALL',CSV,ZCO2,                             &
                            XZENITH,XAZIM,XSW_BANDS,ZDIR_ALB,ZSCA_ALB,  &
                            ZEMIS,ZTSRAD                                )
  !
  IF (SIZE(XEMIS)>0) THEN
    XDIR_ALB = ZDIR_ALB
    XSCA_ALB = ZSCA_ALB
    XEMIS    = ZEMIS
    XTSRAD   = ZTSRAD
    CALL MNHGET_SURF_PARAM_n (PSEA=XSEA)                          
  END IF
ELSE
  !* fields not physically necessary, but must be initialized
  IF (SIZE(XEMIS)>0) THEN
    XDIR_ALB = 0.
    XSCA_ALB = 0.
    XEMIS    = 1.
    XTSRAD   = XTT
    XSEA     = 1.
  END IF
END IF
!
  !* special case after spawning in prep_real_case
IF (CSURF=='EXRM' .AND. CPROGRAM=='REAL  ') CSURF = 'EXTE'
!
DEALLOCATE(ZDIR_ALB)
DEALLOCATE(ZSCA_ALB)
DEALLOCATE(ZEMIS   )
DEALLOCATE(ZTSRAD  )
!
DEALLOCATE(ZCO2)
!
!
!* in a RESTART case, reads surface radiative quantities in the MESONH file
!
IF (CRAD   == 'ECMW' .AND. CGETRAD=='READ') THEN
  CALL INI_SURF_RAD(HINIFILE, CLUOUT, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD)
END IF
!
!
!*      18.3   Mesonh fields
!              -------------
!
IF (CPROGRAM/='REAL  ') CALL MNHREAD_ZS_DUMMY_n(CINIFILE)
!
!-------------------------------------------------------------------------------
!
!*       19.    INITIALIZE THE PARAMETERS FOR THE PHYSICS
!               -----------------------------------------
!
IF (CRAD   == 'ECMW') THEN
!
!* get cover mask for aerosols
!
  IF (CPROGRAM=='MESONH' .OR. CPROGRAM=='DIAG  ') THEN
    IF (CSURF=='EXTE') THEN
      CALL GOTO_SURFEX(KMI)
      CALL MNHGET_SURF_PARAM_n(KCOVER=ICOVER)
      ALLOCATE(ZCOVER(IIU,IJU,ICOVER))
      CALL MNHGET_SURF_PARAM_n(PCOVER=ZCOVER)
    ELSE
      ALLOCATE(ZCOVER(IIU,IJU,255))
      ZCOVER(:,:,:) = 0.
      ZCOVER(:,:,1) = 1.
    END IF
!
    CALL INI_RADIATIONS_ECMWF (HINIFILE,HLUOUT,                                           &
                               XZHAT,XPABST,XTHT,XTSRAD,XLAT,XLON,TDTCUR,                 &
                               CLW,NDLON,NFLEV,NFLUX,NRAD,NSWB,CAER,NAER,NSTATM,          &
                               XSTATM,ZCOVER,XOZON, XAER, LSUBG_COND                      )
!
    DEALLOCATE(ZCOVER)
!
  END IF
ELSE
  ALLOCATE (XOZON(0,0,0))
  ALLOCATE (XAER(0,0,0,0))
END IF
!
!
!
IF (CDCONV /= 'NONE') THEN
  IF (CGETCONV=='INIT') THEN
    GINIDCONV=.TRUE.
  ELSE
    GINIDCONV=.FALSE.
  END IF
!
!  commensurability between convection calling time and time step
!
  XDTCONV=XTSTEP*REAL( INT( (MIN(XDTCONV,1800.)+1.E-10)/XTSTEP ) )
  XDTCONV=MAX( XDTCONV, XTSTEP )
  IF (NVERB>=10) THEN
    WRITE(ILUOUT,*) 'XDTCONV has been set to : ',XDTCONV
  END IF
  CALL INI_DEEP_CONVECTION (HINIFILE,HLUOUT,GINIDCONV,TDTCUR,        &
                           NCOUNTCONV,XDTHCONV,XDRVCONV,XDRCCONV,    &
                           XDRICONV,XPRCONV,XPRSCONV,XPACCONV,       &
                           XUMFCONV,XDMFCONV,XMFCONV,XPRLFLXCONV,XPRSFLXCONV,&
                           XCAPE,NCLTOPCONV,NCLBASCONV,              &
                           TDTDCONV,LCHTRANS,XDSVCONV,               &
                           LCH_CONV_LINOX, XIC_RATE, XCG_RATE,       &
                           XIC_TOTAL_NUMBER, XCG_TOTAL_NUMBER        )

END IF
!
!
!-------------------------------------------------------------------------------
!
!
!*      20.    ALLOCATION OF THE TEMPORAL SERIES
!              ---------------------------------
!
IF (LSERIES .AND. CPROGRAM/='DIAG  ') CALL INI_SERIES_n
!
!-------------------------------------------------------------------------------
!
!
!*      21.   (re)initialize scalar variables
!             -------------------------------
!
!
IF ( LUSECHEM .OR. LCHEMDIAG ) THEN
  IF (CPROGRAM=='MESONH'.OR. CPROGRAM=='DIAG  ' .OR. CPROGRAM=='IDEAL ') &
        CALL CH_INIT_FIELD_n(KMI, ILUOUT, NVERB)
END IF

!
!-------------------------------------------------------------------------------
!
!*      22.    UPDATE HALO
!              -----------
!
!
CALL UPDATE_HALO_ll(TZINITHALO3D_ll,IINFO_ll)
CALL UPDATE_HALO_ll(TZINITHALO2D_ll,IINFO_ll)
CALL CLEANLIST_ll(TZINITHALO3D_ll)
CALL CLEANLIST_ll(TZINITHALO2D_ll)
!
!
!-------------------------------------------------------------------------------
!
!*      23.    DEALLOCATION
!              -------------
!
DEALLOCATE(ZJ)
!
DEALLOCATE(XSTROATM)
DEALLOCATE(XSMLSATM)
DEALLOCATE(XSMLWATM)
DEALLOCATE(XSPOSATM)
DEALLOCATE(XSPOWATM)
!
!-------------------------------------------------------------------------------
!
!*      24.     BALLOON and AIRCRAFT initializations
!              ------------------------------------
!
CALL INI_AIRCRAFT_BALLOON(HINIFILE,CLUOUT,XTSTEP, TDTSEG, XSEGLEN, NRR, NSV,  &
                          CTURB=="TKEL" ,                                     &
                          XLATORI, XLONORI                                    )
!
!-------------------------------------------------------------------------------
!
!*      25.     STATION initializations
!              -----------------------
!
CALL INI_SURFSTATION_n(CLUOUT,XTSTEP, TDTSEG, XSEGLEN, NRR, NSV,  &
                       CTURB=="TKEL" ,                            &
                       XLATORI, XLONORI                           )
!
!-------------------------------------------------------------------------------
!
!*      26.     PROFILER initializations
!              ------------------------
!
CALL INI_POSPROFILER_n(CLUOUT,XTSTEP, TDTSEG, XSEGLEN, NRR, NSV,  &
                       CTURB=="TKEL",                             &
                       XLATORI, XLONORI                           )
!
!-------------------------------------------------------------------------------
!
! 
END SUBROUTINE INI_MODEL_n

