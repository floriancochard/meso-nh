!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG3 2007/11/20 12:37:45
!-----------------------------------------------------------------
!     #################
      MODULE MODI_MODEL_n
!     #################
!
INTERFACE
!
       SUBROUTINE MODEL_n(KTCOUNT,OEXIT)
!
INTEGER, INTENT(IN)   :: KTCOUNT  ! temporal loop index of model KMODEL
LOGICAL, INTENT(INOUT):: OEXIT    ! switch for the end of the temporal loop
!
END SUBROUTINE MODEL_n
!
END INTERFACE
!
END MODULE MODI_MODEL_n

!     ################################### 
      SUBROUTINE MODEL_n(KTCOUNT, OEXIT) 
!     ###################################
!
!!****  *MODEL_n * -monitor of the model version _n 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to build up a typical model version
!     by sequentially calling the specialized routines.
!
!!**  METHOD
!!    ------
!!      Some preliminary initializations are performed in the first section.
!!    Then, specialized routines are called to update the guess of the future
!!    instant XRxxS of the variable xx by adding the effects of all the
!!    different sources of evolution.
!!
!!              (guess of xx at t+dt) * Rhod_ref * Jacobian
!!      XRxxS = -------------------------------------------
!!                           2 dt
!!
!!      At this level, the informations are transferred with a USE association
!!    from the INIT step, where the modules have been previously filled. The
!!    transfer to the subroutines computing each source term is performed by
!!    argument in order to avoid repeated compilations of these subroutines.
!!      This monitor model_n, must therefore be duplicated for each model,
!!    model1 corresponds in this case to the outermost model, model2 is used
!!    for the first level of gridnesting,....  
!!      The effect of all parameterizations is computed in PHYS_PARAM_n, which
!!    is itself a monitor. This is due to a possible large number of
!!    parameterizations, which can be activated and therefore, will require a
!!    very large list of arguments. To circumvent this problem, we transfer by
!!    a USE association, the necessary informations in this monitor, which will
!!    dispatch the pertinent information to every parametrization.
!!      Some elaborated diagnostics, LES tools, budget storages are also called
!!    at this level because they require informations about the fields at every
!!    timestep.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine FMLOOK: to recover the logical unit number linked to a FMfile
!!      Subroutine FMOPEN: to open a FMfile
!!      Subroutine WRITE_DESFM: to write the descriptive part of a FMfile
!!      Subroutine WRITE_LFIFM: to write the binary part of a FMfile
!!      Subroutine SET_MASK   : to compute all the masks selected for budget
!!                         computations
!!      Subroutine BOUNDARIES   : set the fields at the marginal points in every
!!                         directions according the selected boundary conditions
!!      Subroutine INITIAL_GUESS: initializes the guess of the future instant
!!      Subroutine LES_FLX_SPECTRA: computes the resolved fluxes and the
!!                     spectra of some quantities when running in LES mode.
!!      Subroutine ADVECTION: computes the advection terms.
!!      Subroutine DYN_SOURCES: computes the curvature, Coriolis, gravity terms.
!!      Subroutine NUM_DIFF: applies the fourth order numerical diffusion.
!!      Subroutine RELAXATION: performs the relaxation to Larger Scale fields
!!                             in the upper levels and outermost vertical planes
!!      Subroutine PHYS_PARAM_n : computes the parameterized physical terms
!!      Subroutine RAD_BOUND: prepares the velocity normal components for the bc.
!!      Subroutine RESOLVED_CLOUD : computes the sources terms for water in any
!!                                  form
!!      Subroutine PRESSURE : computes the pressure gradient term and the
!!                            absolute pressure
!!      Subroutine EXCHANGE : updates the halo of each subdomains
!!      Subroutine ENDSTEP : advances in time the  fields.
!!      Subroutines UVW_LS_COUPLING and SCALAR_LS_COUPLING:
!!                                 compute the large scale fields, used to
!!                                 couple Model_n with outer informations.
!!      Subroutine ENDSTEP_BUDGET: writes the budget informations.
!!      Subroutine FMCLOS        : closes a FM file
!!      Subroutine ADD_FORECAST_TO_DATE : transform the current time in GMT
!!      Subroutine FORCING : computes forcing terms
!!      Subroutine ADD3DFIELD_ll : add a field to 3D-list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!          MODD_DYN
!!          MODD_CONF
!!          MODD_NESTING
!!          MODD_BUDGET
!!          MODD_PARAMETERS
!!          MODD_CONF_n
!!          MODD_CURVCOR_n
!!          MODD_DYN_n
!!          MODD_DIM_n
!!          MODD_ADV_n
!!          MODD_FIELD_n
!!          MODD_LSFIELD_n
!!          MODD_GRID_n
!!          MODD_METRICS_n
!!          MODD_LBC_n
!!          MODD_PARAM_n
!!          MODD_REF_n
!!          MODD_LUNIT_n
!!          MODD_OUT_n
!!          MODD_TIME_n
!!          MODD_TURB_n
!!          MODD_CLOUDPAR_n
!!          MODD_TIME
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty                  * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/09/94
!!      Modification 20/10/94  (J.Stein) for the outputs and abs_layers routines
!!      Modification 10/11/94  (J.Stein) change ABS_LAYER_FIELDS call
!!      Modification 16/11/94  (J.Stein) add call to the renormalization
!!      Modification 17/11/94  (J.-P. Lafore and J.-P. Pinty) call NUM_DIFF
!!      Modification 08/12/94  (J.Stein) cleaning + remove (RENORM + ABS_LAYER..
!!                             ..) + add RELAXATION + LS fiels in the arguments
!!      Modification 19/12/94  (J.Stein) switch for the num diff
!!      Modification 22/12/94  (J.Stein) update tdtcur + change dyn_source call
!!      Modification 05/01/95  (J.Stein) add the parameterization monitor
!!      Modification 09/01/95  (J.Stein) add the 1D switch
!!      Modification 10/01/95  (J.Stein) displace the TDTCUR computation
!!      Modification 03/01/95  (J.-P. Lafore) Absolute pressure diagnosis
!!      Modification Jan 19, 1995 (J. Cuxart) Shunt the DYN_SOURCES in 1D cases.
!!      Modification Jan 24, 1995 (J. Stein)  Interchange Boundaries and
!!                           Initial_guess to correct a bug in 2D configuration
!!      Modification Feb 02, 1995 (I.Mallet) update BOUNDARIES and RAD_BOUND
!!                                           calls
!!      Modification Mar 10, 1995 (I.Mallet) add call to SET_COUPLING
!!                   March,21, 1995 (J. Stein) remove R from the historical var.
!!                   March,26, 1995 (J. Stein) add the EPS variable
!!                   April 18, 1995 (J. Cuxart) add the LES call
!!                   Sept 20,1995 (Lafore) coupling for the dry mass Md
!!                   Nov   2,1995 (Stein) displace the temporal counter increase
!!                   Jan   2,1996 (Stein) rm the test on the temporal counter
!!      Modification Feb   5,1996 (J. Vila) implementation new advection
!!                                          schemes for scalars
!!      Modification Feb  20,1996 (J.Stein) doctor norm
!!                   Dec95 - Jul96 (Georgelin, Pinty, Mari, Suhre) FORCING
!!                   June 17,1996 (Vincent, Lafore, Jabouille)
!!                                        statistics of computing time
!!                   Aug 8, 1996 (K. Suhre) add chemistry
!!                   October 12, 1996 (J. Stein) save the PSRC value
!!                   Sept 05,1996 (V.Masson) print of loop index for debugging
!!                                           purposes
!!                   July 22,1996 (Lafore) improve write of computing time statistics
!!                   July 29,1996 (Lafore) nesting introduction
!!                   Aug.  1,1996 (Lafore) synchronization between models
!!                   Sept. 4,1996 (Lafore) modification of call to routine SET_COUPLING
!!                                         now splitted in 2 routines
!!                                         (UVW_LS_COUPLING and SCALAR_LS_COUPLING)
!!                   Sept  5,1996 (V.Masson) print of loop index for debugging
!!                                           purposes
!!                   Sept 25,1996 (V.Masson) test for coupling performed here
!!                   Oct. 29,1996 (Lafore)   one-way nesting implementation
!!                   Oct. 12,1996 (J. Stein) save the PSRC value
!!                   Dec. 12,1996 (Lafore)   change call to RAD_BOUND
!!                   Dec. 21,1996 (Lafore)   two-way nesting implementation
!!                   Mar. 12,1997 (Lafore)   introduction of "surfacic" LS fields
!!                   Nov 18, 1996 (J.-P. Pinty) FORCING revisited (translation)
!!                   Dec 04, 1996 (J.-P. Pinty) include mixed-phase clouds
!!                   Dec 20, 1996 (J.-P. Pinty) update the budgets
!!                   Dec 23, 1996 (J.-P. Pinty) add the diachronic file control
!!                   Jan 11, 1997 (J.-P. Pinty) add the deep convection control
!!                   Dec  20,1996 (V.Masson) call boundaries before the writing
!!                   Fev 25, 1997 (P.Jabouille) modify the LES tools
!!                   April 3,1997 (Lafore)      merging of the nesting
!!                                              developments on MASTER3
!!                   Jul.  8,1997 (Lafore)  print control for nesting (NVERB>=7)
!!                   Jul. 28,1997 (Masson)  supress LSTEADY_DMASS
!!                   Aug. 19,1997 (Lafore)  full Clark's formulation introduction
!!                   Sept 26,1997 (Lafore)  LS source calculation at restart
!!                                          (temporarily test to have LS at instant t)
!!                   Jan. 28,1998 (Bechtold) add SST forcing
!!                   fev. 10,1998 (Lafore)  RHODJ computation and storage for budget
!!                   Jul. 10,1998 (Stein )  sequentiel loop for nesting
!!                   Apr. 07,1999 (Stein )  cleaning of the nesting subroutines
!!                   oct. 20,1998 (Jabouille) //
!!                   oct. 20,2000 (J.-P. Pinty) add the C2R2 scheme
!!                   fev. 01,2001 (D.Gazen) add module MODD_NSV for NSV variables
!!                   mar,  4,2002 (V.Ducrocq) call to temporal series
!!                   mar, 8, 2001 (V. Masson) advection of perturbation of theta in neutral cases.
!!                   Nov, 6, 2002 (V. Masson) time counters for budgets & LES
!!                   mars 20,2001 (Pinty)   add ICE4 and C3R5 options
!!                   jan. 2004    (Masson)  surface externalization
!!                   sept 2004 (M. Tomasini) Cloud mixing length modification
!!                   june 2005 (P. Tulet)  add aerosols / dusts
!!                   Jul. 2005 (N. Asencio)  two_way and phys_param calls: 
!!                             Add the surface parameters : precipitating 
!!                             hydrometeors, Short and Long Wave , MASKkids array 
!!                   april 2006 (T.Maric) Add halo related to 4th order advection scheme
!!                   May 2006 Remove KEPS
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODE_ll
!
USE MODE_FM
!
USE MODD_TIME
USE MODD_DYN
USE MODD_DYNZD
USE MODD_CONF
USE MODD_NESTING
USE MODD_FMOUT
USE MODD_BUDGET
USE MODD_PARAMETERS
USE MODD_PARAM_ICE,        ONLY : LWARM,LSEDIC
USE MODD_FRC
USE MODD_AIRCRAFT_BALLOON
USE MODD_STATION_n
USE MODD_PROFILER_n
USE MODD_PARAM_C2R2,      ONLY : LSEDC, LRAIN, LACTIT
USE MODD_PARAM_C1R3,      ONLY : LSEDI, LHHONI
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_LUNIT
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_SERIES, ONLY: LSERIES
USE MODD_TURB_CLOUD, ONLY: NMODEL_CLOUD,CTURBLEN_CLOUD,XCEI
!
USE MODD_SUB_MODEL_n
USE MODD_GET_n
USE MODD_CONF_n
USE MODD_CURVCOR_n
USE MODD_DIM_n
USE MODD_DYN_n
USE MODD_DYNZD_n
USE MODD_ADV_n
USE MODD_FIELD_n
USE MODD_LSFIELD_n
USE MODD_GRID_n
USE MODD_METRICS_n
USE MODD_LBC_n
USE MODD_PARAM_n
USE MODD_REF_n
USE MODD_LUNIT_n
USE MODD_OUT_n
USE MODD_TIME_n 
USE MODD_TURB_n
USE MODD_CLOUDPAR_n
USE MODD_PRECIP_n
USE MODD_BIKHARDT_n
USE MODD_DEEP_CONVECTION_n
USE MODD_NSV
USE MODD_RADIATIONS_n, ONLY : XTSRAD,XSCAFLASWD,XDIRFLASWD,XDIRSRFSWD
USE MODD_SERIES_n, ONLY: NFREQSERIES
USE MODD_CH_MNHC_n,    ONLY: LUSECHEM,LCH_CONV_LINOX
USE MODD_NUDGING_n
!
USE MODI_INITIAL_GUESS
USE MODI_BOUNDARIES
USE MODI_ADVECTION
USE MODI_DYN_SOURCES
USE MODI_RELAXATION
USE MODI_NUM_DIFF
USE MODI_PHYS_PARAM_n
USE MODI_RAD_BOUND
USE MODI_PRESSUREZ
USE MODI_ENDSTEP
USE MODI_EXCHANGE
USE MODI_RESOLVED_CLOUD
USE MODI_LS_COUPLING
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_LFIFM_n
USE MODI_MNHWRITE_ZS_DUMMY_n
USE MODI_ENDSTEP_BUDGET
USE MODI_BUDGET_FLAGS
USE MODI_ADD_FORECAST_TO_DATE
USE MODI_FORCING
USE MODI_NUDGING
USE MODI_TEMPORAL_DIST
USE MODI_WRITE_LFIFMN_FORDIACHRO_n
USE MODI_MENU_DIACHRO
USE MODI_MASK_COMPRESS
USE MODI_CART_COMPRESS
USE MODI_SHUMAN
USE MODI_ONE_WAY_n
USE MODI_TWO_WAY
USE MODI_SPAWN_LS_n
USE MODI_LES_INI_TIMESTEP_n
USE MODI_WRITE_LES_n
USE MODI_AIRCRAFT_BALLOON
USE MODI_WRITE_AIRCRAFT_BALLOON
USE MODI_UPDATE_NSV
USE MODI_PROFILER_n
USE MODI_STATION_n
USE MODI_WRITE_SERIES_n
USE MODI_WRITE_PROFILER_n
USE MODI_WRITE_STATION_n
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_INI_DIAG_IN_RUN
USE MODI_END_DIAG_IN_RUN
USE MODI_TURB_CLOUD_INDEX
USE MODI_INI_LG
USE MODE_MODELN_HANDLER
!
!JUAN Z_SPLITTING
USE MODD_TIMEZ
USE MODE_MNH_TIMING
!JUAN Z_SPLITTING
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
INTEGER, INTENT(IN)   :: KTCOUNT
LOGICAL, INTENT(INOUT):: OEXIT
!
!*       0.2   declarations of local variables
!
INTEGER :: ILUOUT      ! Logical unit number for the output listing
INTEGER :: IIU,IJU,IKU ! array size in first, second and third dimensions
INTEGER :: JSV,JRR     ! Loop index for scalar and moist variables
CHARACTER (LEN=28) :: YFMFILE   ! name of the OUTPUT FM-file
CHARACTER (LEN=28) :: YDADFILE  ! name of the corresponding DAD model OUTPUT FM-file
CHARACTER (LEN=4)  :: YNUMBER   ! character string for the OUTPUT FM-file number
CHARACTER (LEN=4)  :: YDADNUMBER! character string for the DAD model OUTPUT FM-file number
CHARACTER (LEN=32) :: YDESFM    ! name of the desfm part of this FM-file
INTEGER  :: INBVAR              ! number of HALO2_lls to allocate
INTEGER  :: IRESP               ! return code in FM routines
INTEGER  :: IINFO_ll            ! return code of parallel routine
INTEGER :: INPRAR               ! number of articles predicted  in
                                !  the LFIFM file
INTEGER :: ININAR               ! number of articles  present in
                                !  the LFIFM file
INTEGER :: ITYPE                ! type of file (cpio or not)
INTEGER :: JOUT                 ! loop index on the output instant list
INTEGER :: IOUTDAD              ! numero of the OUTPUT FM-file of DAD model
INTEGER :: JOUTDAD              ! loop index on the output instant list for DAD model
LOGICAL :: GSTEADY_DMASS        ! conditional call to mass computation
!
                                ! for computing time analysis
!JUAN
REAL,DIMENSION(2)         :: ZTIME,ZTIME1,ZTIME2,ZEND,ZTOT,ZALL,ZTOT_PT
!
REAL,DIMENSION(2)         :: ZTIME_STEP,ZTIME_STEP_PTS
CHARACTER                 :: YMI
INTEGER                   :: IPOINTS
CHARACTER(len=12)         :: YTCOUNT,YPOINTS
!JUAN

REAL         :: ZSTAT_CSTORE,ZSTAT_CBOUND,ZSTAT_CGUESS,ZSTAT_CADV,ZSTAT_CSOURCES
REAL         :: ZSTAT_CDIFF,ZSTAT_CRELAX,ZSTAT_CPARAM
REAL         :: ZSTAT_CSPECTRA,ZSTAT_CRAD_BOUND,ZSTAT_CPRESS
REAL         :: ZSTAT_CCLOUD,ZSTAT_CSTEP_SWA,ZSTAT_CSTEP_MISC
REAL         :: ZSTAT_CCOUPL,ZSTAT_CSTEP_BUD
REAL         :: ZSTAT_CTURB,ZSTAT_C1WAY,ZSTAT_C2WAY
REAL         :: ZSTAT_CRAD,ZSTAT_CDCONV,ZSTAT_CGROUND,ZSTAT_CHALO
REAL         :: ZSTAT_CFORCING,ZSTAT_CNUDGING,ZSTAT_CCHEM
!
REAL         :: ZTSTEP       ! Double timestep except for cold start (single)
REAL         :: ZTSTEP_MET,ZTSTEP_SV ! Effective time step for advection
!
INTEGER :: ISYNCHRO          ! model synchronic index relative to its father
                             ! = 1  for the first time step in phase with DAD
                             ! = 0  for the last  time step (out of phase)
INTEGER      :: IMI ! Current model index
REAL, DIMENSION(:,:),ALLOCATABLE          :: ZSEA
REAL, DIMENSION(:,:),ALLOCATABLE          :: ZTOWN
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
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSTHM,DPTR_XLSRVM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS,DPTR_XLSTHS,DPTR_XLSRVS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLBXUS,DPTR_XLBYUS,DPTR_XLBXVS,DPTR_XLBYVS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLBXWS,DPTR_XLBYWS,DPTR_XLBXTHS,DPTR_XLBYTHS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XLBXTKES,DPTR_XLBYTKES
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XLBXRS,DPTR_XLBYRS,DPTR_XLBXSVS,DPTR_XLBYSVS
!
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XRHODJ,DPTR_XUM,DPTR_XVM,DPTR_XWM,DPTR_XTHM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XTKEM,DPTR_XRUS,DPTR_XRVS,DPTR_XRWS,DPTR_XRTHS
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XRTKES,DPTR_XDIRFLASWD,DPTR_XSCAFLASWD,DPTR_XDIRSRFSWD
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XRM,DPTR_XSVM,DPTR_XRRS,DPTR_XRSVS
REAL, DIMENSION(:,:), POINTER :: DPTR_XINPRC,DPTR_XINPRR,DPTR_XINPRS,DPTR_XINPRG
REAL, DIMENSION(:,:), POINTER :: DPTR_XINPRH,DPTR_XPRCONV,DPTR_XPRSCONV
LOGICAL, DIMENSION(:,:),POINTER :: DPTR_GMASKkids

!
!-------------------------------------------------------------------------------
!
!*        1    PRELIMINARY
!              ------------
ITYPE = 1
IMI = GET_CURRENT_MODEL_INDEX()
!
IF ((CCONF=='START').AND.(KTCOUNT==1)) THEN
  ZTSTEP = XTSTEP
ELSE
  ZTSTEP = 2.*XTSTEP
END IF
!
IF (CMET_ADV_SCHEME(1:3) == 'PPM') THEN
  ZTSTEP_MET = XTSTEP
ELSE
  ZTSTEP_MET = ZTSTEP
END IF
!
IF (CSV_ADV_SCHEME(1:3) == 'PPM') THEN
  ZTSTEP_SV = XTSTEP
ELSE
  ZTSTEP_SV = ZTSTEP
END IF
!
!*       1.0   update NSV_* variables for current model
!              ----------------------------------------
CALL UPDATE_NSV(IMI)
!
!*       1.1   RECOVER THE LOGICAL UNIT NUMBER FOR THE OUTPUT PRINTS
!
CALL FMLOOK_ll(CLUOUT,CLUOUT,ILUOUT,IRESP)
!
!*       1.2   SET ARRAY SIZE
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=NKMAX+2*JPVEXT
!
IF (IMI==1) THEN
  GSTEADY_DMASS=LSTEADYLS
ELSE
  GSTEADY_DMASS=.FALSE.
END IF
!
!*       1.3   OPEN THE DIACHRONIC FILE
!
IF (KTCOUNT == 1) THEN
!
  NULLIFY(TZFIELDS_ll,TZLSFIELD_ll,TZFIELDM_ll)
  NULLIFY(TZHALO2M_ll)
  NULLIFY(TZLSHALO2_ll)
  NULLIFY(TZFIELDT_ll)
  NULLIFY(TZHALO2T_ll)
  NULLIFY(TZFIELDMT_ll)
  NULLIFY(TZHALO2MT_ll)
  NULLIFY(TZFIELDSC_ll)
  NULLIFY(TZHALO2SC_ll)
!
  ALLOCATE(ZWT_ACT_NUC(SIZE(XWT,1),SIZE(XWT,2),SIZE(XWT,3)))
  ALLOCATE(GMASKkids(SIZE(XWT,1),SIZE(XWT,2)))
!
! initialization of the FM file output number
  IOUT=0
!
  INPRAR = 50
  CALL FMOPEN_ll(CFMDIAC,'WRITE',CLUOUT,INPRAR,ITYPE,NVERB,ININAR,IRESP)
  YDESFM=ADJUSTL(ADJUSTR(CFMDIAC)//'.des')
  CALL WRITE_DESFM_n(IMI,YDESFM,CLUOUT)
  CALL WRITE_LFIFMN_FORDIACHRO_n(CFMDIAC)
!
!*       1.4   Initialization of the list of fields for the halo updates
!
!                 a) Sources terms
!
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRUS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRVS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRWS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRTHS)
  IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDS_ll, XRTKES)
  DO JRR=1,NRR
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XRRS(:,:,:,JRR))
  ENDDO
  DO JSV=1,NSV
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XRSVS(:,:,:,JSV))
  ENDDO
  IF (SIZE(XSRCM,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDS_ll, XSRCM)
  !
  IF (LNUMDIFF .AND. NHALO==1 ) THEN
  !
  !                 b) LS fields
  !
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSUM)
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSVM)
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSWM)
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSTHM)
    IF (NRR >= 1) THEN
      CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSRVM)
    ENDIF
  !
  !                 c) Fields at t-dt
  !
    CALL ADD3DFIELD_ll(TZFIELDM_ll, XUM)
    CALL ADD3DFIELD_ll(TZFIELDM_ll, XVM)
    CALL ADD3DFIELD_ll(TZFIELDM_ll, XWM)
    CALL ADD3DFIELD_ll(TZFIELDM_ll, XTHM)
    IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDM_ll, XTKEM)
    DO JRR=1,NRR
      CALL ADD3DFIELD_ll(TZFIELDM_ll, XRM(:,:,:,JRR))
    ENDDO
    DO JSV=1,NSV
      CALL ADD3DFIELD_ll(TZFIELDM_ll, XSVM(:,:,:,JSV))
    ENDDO
  !
  !*       1.5   Initialize the list of fields for the halo updates (2nd layer)
  !
    INBVAR = 4+NRR+NSV
    IF (SIZE(XRTKES,1) /= 0) INBVAR=INBVAR+1
    IF( NHALO==1 ) CALL INIT_HALO2_ll(TZHALO2M_ll,INBVAR,IIU,IJU,IKU)
    IF( NHALO==1 ) CALL INIT_HALO2_ll(TZLSHALO2_ll,4+MIN(1,NRR),IIU,IJU,IKU)
  !
  !*       1.6   Initialise the 2nd layer of the halo of the LS fields
  !
    IF ( LSTEADYLS .AND. NHALO==1 ) CALL UPDATE_HALO2_ll(TZLSFIELD_ll, TZLSHALO2_ll, IINFO_ll)
  END IF
  !
  !        1.7  Initialize the list of fields (2nd layer) time t if 4th order
!             advection schemes are used
!
   IF( CUVW_ADV_SCHEME == "CEN4TH" .AND. NHALO==1 ) THEN
!
      CALL ADD3DFIELD_ll(TZFIELDMT_ll, XUT)
      CALL ADD3DFIELD_ll(TZFIELDMT_ll, XVT)
      CALL ADD3DFIELD_ll(TZFIELDMT_ll, XWT)
!
      INBVAR = 3
      IF( NHALO==1 ) CALL INIT_HALO2_ll(TZHALO2MT_ll,INBVAR,IIU,IJU,IKU)
!
   END IF
!
   IF( CMET_ADV_SCHEME == "CEN4TH" .AND. NHALO==1 ) THEN 
!
      CALL ADD3DFIELD_ll(TZFIELDT_ll, XTHT)
      IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDT_ll, XTKET)
      DO JRR=1,NRR
         CALL ADD3DFIELD_ll(TZFIELDT_ll, XRT(:,:,:,JRR))
      ENDDO
!
      INBVAR = 1+NRR
      IF (SIZE(XRTKES,1) /= 0) INBVAR=INBVAR+1
      IF( NHALO==1 ) CALL INIT_HALO2_ll(TZHALO2T_ll,INBVAR,IIU,IJU,IKU)
!
   END IF
!
   IF( CSV_ADV_SCHEME == "CEN4TH" .AND. NHALO==1 ) THEN 
!
      DO JSV=1,NSV
         CALL ADD3DFIELD_ll(TZFIELDSC_ll, XSVT(:,:,:,JSV))
      ENDDO

      INBVAR = NSV
      IF( NHALO==1 ) CALL INIT_HALO2_ll(TZHALO2SC_ll,INBVAR,IIU,IJU,IKU)
!
   END IF
!
  !
  XT_START = 0.0
  ! 
  XT_STORE     = 0.0
  XT_BOUND     = 0.0
  XT_GUESS     = 0.0
  XT_FORCING   = 0.0
  XT_NUDGING   = 0.0
  XT_ADV       = 0.0
  XT_SOURCES   = 0.0
  !
  XT_DIFF      = 0.0
  XT_RELAX     = 0.0
  XT_PARAM     = 0.0
  XT_SPECTRA   = 0.0
  XT_HALO      = 0.0
  XT_RAD_BOUND = 0.0
  XT_PRESS     = 0.0
  !
  XT_CLOUD     = 0.0
  XT_STEP_SWA  = 0.0
  XT_STEP_MISC = 0.0
  XT_COUPL     = 0.0
  XT_1WAY      = 0.0
  XT_STEP_BUD  = 0.0
  !
  XT_RAD       = 0.0
  XT_DCONV     = 0.0
  XT_GROUND    = 0.0
  XT_TURB      = 0.0
  XT_CHEM      = 0.0
  XT_2WAY      = 0.0
  !
END IF
!
!*       1.7   Allocation of arrays for observation diagnostics
!
CALL INI_DIAG_IN_RUN(IIU,IJU,IKU,LFLYER,LSTATION,LPROFILER)
!
!
CALL SECOND_MNH2(ZEND)
!
!-------------------------------------------------------------------------------
!
!*       2.    ONE-WAY NESTING AND LARGE SCALE FIELD REFRESH
!              ---------------------------------------------
!
!
CALL SECOND_MNH2(ZTIME1)
!
ISYNCHRO = MODULO (KTCOUNT, NDTRATIO(IMI) )      ! test of synchronisation
!


IF (IMI/=1 .AND. NDAD(IMI)/=IMI .AND. (ISYNCHRO==1 .OR. NDTRATIO(IMI) == 1) ) THEN     
!                                                                        
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
  DPTR_XCOEFLIN_LBXM=>XCOEFLIN_LBXM
  DPTR_XLSTHM=>XLSTHM
  DPTR_XLSRVM=>XLSRVM
  DPTR_XLSUM=>XLSUM
  DPTR_XLSVM=>XLSVM
  DPTR_XLSWM=>XLSWM
  DPTR_XLSTHS=>XLSTHS
  DPTR_XLSRVS=>XLSRVS
  DPTR_XLSUS=>XLSUS
  DPTR_XLSVS=>XLSVS
  DPTR_XLSWS=>XLSWS
  !
  IF ( LSTEADYLS                     ) THEN
    NCPL_CUR=0
  ELSE
    IF (NCPL_CUR/=1) THEN
      IF ( KTCOUNT+1 == NCPL_TIMES(NCPL_CUR-1,IMI)  ) THEN
        !
        !  LS sources are interpolated from the LS field 
        ! values of model DAD(IMI)
        CALL SPAWN_LS_n(NDAD(IMI),XTSTEP,IMI,                        &
             DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
             DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
             NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),                    &
             DPTR_CLBCX,DPTR_CLBCY,DPTR_XZZ,DPTR_XZHAT,LSLEVE,XLEN1,XLEN2,DPTR_XCOEFLIN_LBXM, &
             DPTR_XLSTHM,DPTR_XLSRVM,DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,                        &
             DPTR_XLSTHS,DPTR_XLSRVS,DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS                         )
      END IF
    END IF
    !
  END IF
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
  !
  DPTR_NKLIN_LBXM=>NKLIN_LBXM
  DPTR_XCOEFLIN_LBXM=>XCOEFLIN_LBXM
  DPTR_NKLIN_LBYM=>NKLIN_LBYM
  DPTR_XCOEFLIN_LBYM=>XCOEFLIN_LBYM
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
  !  
  DPTR_XLBXUS=>XLBXUS
  DPTR_XLBYUS=>XLBYUS
  DPTR_XLBXVS=>XLBXVS
  DPTR_XLBYVS=>XLBYVS
  DPTR_XLBXWS=>XLBXWS
  DPTR_XLBYWS=>XLBYWS
  DPTR_XLBXTHS=>XLBXTHS
  DPTR_XLBYTHS=>XLBYTHS
  DPTR_XLBXTKES=>XLBXTKES
  DPTR_XLBYTKES=>XLBYTKES
  DPTR_XLBXRS=>XLBXRS
  DPTR_XLBYRS=>XLBYRS
  DPTR_XLBXSVS=>XLBXSVS
  DPTR_XLBYSVS=>XLBYSVS
  !
  CALL ONE_WAY_n(NDAD(IMI),CLUOUT,XTSTEP,IMI,KTCOUNT,                   &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
       NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),NDTRATIO(IMI),         &
       DPTR_CLBCX,DPTR_CLBCY,NRIMX,NRIMY,                                &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,      &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,      &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,      &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM,      &
       GSTEADY_DMASS,CCLOUD,                                   &
       DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM,DPTR_XLBXWM,DPTR_XLBYWM,              &
       DPTR_XLBXTHM,DPTR_XLBYTHM,                                        &
       DPTR_XLBXTKEM,DPTR_XLBYTKEM,                                      &
       DPTR_XLBXRM,DPTR_XLBYRM,DPTR_XLBXSVM,DPTR_XLBYSVM,                          &
       XDRYMASST,XDRYMASSS,                                    &
       DPTR_XLBXUS,DPTR_XLBYUS,DPTR_XLBXVS,DPTR_XLBYVS,DPTR_XLBXWS,DPTR_XLBYWS,              &
       DPTR_XLBXTHS,DPTR_XLBYTHS,                                        &
       DPTR_XLBXTKES,DPTR_XLBYTKES,                                      &
       DPTR_XLBXRS,DPTR_XLBYRS,DPTR_XLBXSVS,DPTR_XLBYSVS                           )
  !
END IF
!
CALL SECOND_MNH2(ZTIME2)                                                  
XT_1WAY = XT_1WAY + ZTIME2 - ZTIME1                                
!
!-------------------------------------------------------------------------------
!
!*       3.    LATERAL BOUNDARY CONDITIONS EXCEPT FOR NORMAL VELOCITY
!              ------------------------------------------------------
!
ZTIME1=ZTIME2
!
!*       3.1   Set the lagragian variables values at the LB
!
IF( LLG .AND. IMI==1 ) CALL SETLB_LG
!
CALL BOUNDARIES (                                                   &
            XTSTEP,CLBCX,CLBCY,NRR,NSV,KTCOUNT,                     &
            XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
            XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
            XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,   &
            XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS,   &
            XRHODJ,                                                 &
            XUM, XVM, XWM, XTHM, XTKEM, XRM, XSVM,XSRCM,            &
            XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT                   )
!
CALL SECOND_MNH2(ZTIME2)
!
XT_BOUND = XT_BOUND + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!* initializes surface number
IF (CSURF=='EXTE') CALL GOTO_SURFEX(IMI)
!-------------------------------------------------------------------------------
!
!*       4.    STORAGE IN A SYNCHRONOUS FILE
!              -----------------------------
!
ZTIME1 = ZTIME2
!
DO JOUT = 1,NOUT_NUMB
  IF (KTCOUNT == NOUT_TIMES(JOUT)) THEN
    IOUT=IOUT+1
    GCLOSE_OUT=.TRUE.
    INPRAR = 22 +2*(4+NRR+NSV)
    WRITE (YNUMBER,FMT="('.',I3.3)") IOUT
    YFMFILE=ADJUSTL(ADJUSTR(COUTFILE)//YNUMBER)
!
!        search for the corresponding Output of the DAD model
!
    IF (NDAD(IMI) == IMI .OR.  IMI == 1) THEN
      YDADFILE=YFMFILE
    ELSE
      IOUTDAD=0
      DO JOUTDAD =1,JPOUTMAX
        IF ( XFMOUT(NDAD(IMI),JOUTDAD) /= XUNDEF .AND.                 &
             XFMOUT(NDAD(IMI),JOUTDAD) <= (XFMOUT(IMI,JOUT)+1.E-10) )   &
                     IOUTDAD=IOUTDAD+1
      END DO
      IF(IOUTDAD>0) THEN
        WRITE (YDADNUMBER,FMT="('.',I3.3)") IOUTDAD
        YDADFILE=ADJUSTL(ADJUSTR(CDAD_NAME(IMI))//YDADNUMBER)
      ELSE
        WRITE (YDADFILE,FMT="('NO_DAD_FILE')")
      END IF
    END IF
!
    CALL FMOPEN_ll(YFMFILE,'WRITE',CLUOUT,INPRAR,ITYPE,NVERB,ININAR,IRESP)
    YDESFM=ADJUSTL(ADJUSTR(YFMFILE)//'.des')
!  
    CALL WRITE_DESFM_n(IMI,YDESFM,CLUOUT)
    CALL WRITE_LFIFM_n(YFMFILE,YDADFILE)
    COUTFMFILE = YFMFILE
    CALL MNHWRITE_ZS_DUMMY_n(CPROGRAM)
    IF (CSURF=='EXTE') THEN
      CALL GOTO_SURFEX(IMI)
      CALL WRITE_SURF_ATM_n('MESONH','ALL')
      CALL DIAG_SURF_ATM_n('MESONH')
      CALL WRITE_DIAG_SURF_ATM_n('MESONH','ALL')
    END IF
    !
    ! Reinitialize Lagragian variables at every model output
    IF (LLG .AND. LINIT_LG .AND. CINIT_LG=='FMOUT') THEN
      CALL INI_LG(XXHAT,XYHAT,XZZ,XSVM,XSVT,XLBXSVM,XLBYSVM)
      IF (NVERB>=5) THEN
        WRITE(UNIT=ILUOUT,FMT=*) '************************************'
        WRITE(UNIT=ILUOUT,FMT=*) '*** Lagrangian variables refreshed after ',TRIM(YFMFILE),' output'
        WRITE(UNIT=ILUOUT,FMT=*) '************************************'
      END IF
    END IF
!
  END IF
!
END DO
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STORE = XT_STORE + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       5.    INITIALIZATION OF THE BUDGET VARIABLES
!              --------------------------------------
!
IF (NBUMOD==IMI) THEN
  LBU_ENABLE = CBUTYPE /='NONE'.AND. CBUTYPE /='SKIP' .AND. MODULO(KTCOUNT,2)==1
ELSE
  LBU_ENABLE = .FALSE.
END IF
!
IF (NBUMOD==IMI .AND. CBUTYPE=='MASK' .AND. MODULO(KTCOUNT,2)==1) THEN
  CALL SET_MASK
  IF (LBU_RU)   XBURHODJU(:,NBUTIME,:) = XBURHODJU(:,NBUTIME,:)    &
                            + 2.*MASK_COMPRESS(MXM(XRHODJ))
  IF (LBU_RV)   XBURHODJV(:,NBUTIME,:) = XBURHODJV(:,NBUTIME,:)    &
                            + 2.*MASK_COMPRESS(MYM(XRHODJ))
  IF (LBU_RW)   XBURHODJW(:,NBUTIME,:) = XBURHODJW(:,NBUTIME,:)    &
                            + 2.*MASK_COMPRESS(MZM(XRHODJ))
  IF (ALLOCATED(XBURHODJ))                                         &
                XBURHODJ (:,NBUTIME,:) = XBURHODJ (:,NBUTIME,:)    &
                              + 2.*MASK_COMPRESS(XRHODJ)
END IF
!
IF (NBUMOD==IMI .AND. CBUTYPE=='CART' .AND. MODULO(KTCOUNT,2)==1) THEN
  IF (LBU_RU)   XBURHODJU(:,:,:) = XBURHODJU(:,:,:)    &
                + 2.*CART_COMPRESS(MXM(XRHODJ))
  IF (LBU_RV)   XBURHODJV(:,:,:) = XBURHODJV(:,:,:)    &
                + 2.*CART_COMPRESS(MYM(XRHODJ))
  IF (LBU_RW)   XBURHODJW(:,:,:) = XBURHODJW(:,:,:)    &
                + 2.*CART_COMPRESS(MZM(XRHODJ))
  IF (ALLOCATED(XBURHODJ))                             &
                XBURHODJ (:,:,:) = XBURHODJ (:,:,:)    &
                + 2.*CART_COMPRESS(XRHODJ)
END IF
!
CALL BUDGET_FLAGS(LUSERV, LUSERC, LUSERR,         &
                  LUSERI, LUSERS, LUSERG, LUSERH  )
!
XTIME_BU   = 0.0
!
!-------------------------------------------------------------------------------
!
!*       6.    INITIALIZATION OF THE FIELD TENDENCIES
!              --------------------------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
!
CALL INITIAL_GUESS ( NRR, NSV, KTCOUNT, XRHODJ,IMI,                         &
                     XUM, XVM, XWM, XTHM, XRM, XTKEM, XSVM,                 &
                     ZTSTEP,                                                &
                     XRUS, XRVS, XRWS, XRTHS, XRRS, XRTKES, XRSVS           )

CALL SECOND_MNH2(ZTIME2)
!
XT_GUESS = XT_GUESS + ZTIME2 - ZTIME1 - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       7.    INITIALIZATION OF THE LES FOR CURRENT TIME-STEP
!              -----------------------------------------------
!
XTIME_LES_BU   = 0.0
XTIME_LES      = 0.0
IF (LLES) CALL LES_INI_TIMESTEP_n(KTCOUNT)
!
!-------------------------------------------------------------------------------
!
!*       8.    TWO-WAY INTERACTIVE GRID-NESTING
!              --------------------------------
!
!
CALL SECOND_MNH2(ZTIME1)
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
GMASKkids(:,:)=.FALSE.
!
IF (NMODEL>1) THEN
  ! correct an ifort bug
  DPTR_XRHODJ=>XRHODJ
  DPTR_XUM=>XUM
  DPTR_XVM=>XVM
  DPTR_XWM=>XWM
  DPTR_XTHM=>XTHM
  DPTR_XRM=>XRM
  DPTR_XTKEM=>XTKEM
  DPTR_XSVM=>XSVM
  DPTR_XRUS=>XRUS
  DPTR_XRVS=>XRVS
  DPTR_XRWS=>XRWS
  DPTR_XRTHS=>XRTHS
  DPTR_XRRS=>XRRS
  DPTR_XRTKES=>XRTKES
  DPTR_XRSVS=>XRSVS
  DPTR_XINPRC=>XINPRC
  DPTR_XINPRR=>XINPRR
  DPTR_XINPRS=>XINPRS
  DPTR_XINPRG=>XINPRG
  DPTR_XINPRH=>XINPRH
  DPTR_XPRCONV=>XPRCONV
  DPTR_XPRSCONV=>XPRSCONV
  DPTR_XDIRFLASWD=>XDIRFLASWD
  DPTR_XSCAFLASWD=>XSCAFLASWD
  DPTR_XDIRSRFSWD=>XDIRSRFSWD
  DPTR_GMASKkids=>GMASKkids
  !
  CALL TWO_WAY(     CLUOUT,NRR,NSV,KTCOUNT,DPTR_XRHODJ,IMI,XTSTEP,        &        
       DPTR_XUM ,DPTR_XVM ,DPTR_XWM , DPTR_XTHM, DPTR_XRM, DPTR_XTKEM, DPTR_XSVM,              &        
       DPTR_XRUS,DPTR_XRVS,DPTR_XRWS,DPTR_XRTHS,DPTR_XRRS,DPTR_XRTKES,DPTR_XRSVS,              &
       DPTR_XINPRC,DPTR_XINPRR,DPTR_XINPRS,DPTR_XINPRG,DPTR_XINPRH,DPTR_XPRCONV,DPTR_XPRSCONV, &
       DPTR_XDIRFLASWD,DPTR_XSCAFLASWD,DPTR_XDIRSRFSWD,DPTR_GMASKkids           )
END IF
!                                                       
CALL SECOND_MNH2(ZTIME2)
XT_2WAY = XT_2WAY + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       9.    ADVECTION
!              ---------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (LNEUTRAL) XTHT=XTHT-XTHVREF
!
IF (CUVW_ADV_SCHEME == "CEN4TH" .AND. NHALO == 1 ) THEN
   CALL UPDATE_HALO2_ll(TZFIELDMT_ll, TZHALO2MT_ll, IINFO_ll)
ENDIF
!
IF (CMET_ADV_SCHEME == "CEN4TH" .AND. NHALO == 1 ) THEN
   CALL UPDATE_HALO2_ll(TZFIELDT_ll,  TZHALO2T_ll,  IINFO_ll)
ENDIF
!
IF (CSV_ADV_SCHEME == "CEN4TH" .AND. NHALO == 1 ) THEN
   CALL UPDATE_HALO2_ll(TZFIELDSC_ll, TZHALO2SC_ll, IINFO_ll)
ENDIF
!
IF (CMET_ADV_SCHEME(1:3) == 'PPM' .OR. CSV_ADV_SCHEME(1:3) == 'PPM') THEN
   CALL ADVECTION ( CUVW_ADV_SCHEME, CMET_ADV_SCHEME, CSV_ADV_SCHEME,   &
                    NLITER, CLBCX, CLBCY,                               &
                    NRR, NSV, KTCOUNT, ZTSTEP_MET, ZTSTEP_SV,           &
                    XUM, XVM, XWM, XTHM, XRM, XTKEM, XSVM,              &
                    XUT, XVT, XWT, XTHT, XRT, XTKET, XSVT,              &
                    XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,               &
                    XRUS, XRVS, XRWS, XRTHS, XRRS, XRTKES,              &
                    XRSVS, TZHALO2MT_ll, TZHALO2T_ll, TZHALO2SC_ll,     &
                    XRTHMS, XRRMS, XRTKEMS, XRSVMS                      )
ELSE
   CALL ADVECTION ( CUVW_ADV_SCHEME, CMET_ADV_SCHEME, CSV_ADV_SCHEME,   &
                    NLITER, CLBCX, CLBCY,                               &
                    NRR, NSV, KTCOUNT, ZTSTEP_MET, ZTSTEP_SV,           &
                    XUM, XVM, XWM, XTHM, XRM, XTKEM, XSVM,              &
                    XUT, XVT, XWT, XTHT, XRT, XTKET, XSVT,              &
                    XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,               &
                    XRUS, XRVS, XRWS, XRTHS, XRRS, XRTKES,              &
                    XRSVS, TZHALO2MT_ll, TZHALO2T_ll, TZHALO2SC_ll      )
ENDIF
!
IF (LNEUTRAL) XTHT=XTHT+XTHVREF
!
IF (NMODEL_CLOUD==IMI .AND. CTURBLEN_CLOUD/='NONE') THEN
  CALL TURB_CLOUD_INDEX(ZTSTEP,YFMFILE,CLUOUT,                    &
                        LTURB_DIAG,GCLOSE_OUT,NRRI,               &
                        XRRS,XRM,XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY, &
                        XCEI )
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ADV = XT_ADV + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       10.    FORCING
!               -------
!
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF ( LFORCING ) THEN
  CALL FORCING(LUSERV,XRHODJ,XCORIOZ,XZHAT,XZZ,TDTCUR,&
               XUT,XVT,XWT,XTHT,XTKET,XRT,XSVT,       &
               XUM,XVM,XWM,XTHM,XTKEM,XRM,XSVM,       &
               XRUS,XRVS,XRWS,XRTHS,XRTKES,XRRS,XRSVS,IMI)
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_FORCING = XT_FORCING + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       11.    NUDGING
!               -------
!
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF ( LNUDGING ) THEN
  CALL NUDGING(LUSERV,XRHODJ,XTNUDGING,         &
               XUM,XVM,XWM,XTHM,XRM,            &
               XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM, &
               XRUS,XRVS,XRWS,XRTHS,XRRS)

END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_NUDGING = XT_NUDGING + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       12.    DYNAMICAL SOURCES
!               -----------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF( LTRANS ) THEN
  XUT(:,:,:) = XUT(:,:,:) + XUTRANS
  XVT(:,:,:) = XVT(:,:,:) + XVTRANS
END IF
!
CALL DYN_SOURCES( NRR,NRRL, NRRI,IMI,                           &
                  XUT, XVT, XWT, XTHT, XRT,                    &
                  XCORIOX, XCORIOY, XCORIOZ, XCURVX, XCURVY,   &
                  XRHODJ, XZZ, XTHVREF, XEXNREF,               &
                  XRUS, XRVS, XRWS, XRTHS                      )
!
IF( LTRANS ) THEN
  XUT(:,:,:) = XUT(:,:,:) - XUTRANS
  XVT(:,:,:) = XVT(:,:,:) - XVTRANS
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_SOURCES = XT_SOURCES + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       13.    NUMERICAL DIFFUSION
!               -------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF ( LNUMDIFF ) THEN
!
  IF( NHALO==1 ) THEN
    CALL UPDATE_HALO2_ll(TZFIELDM_ll, TZHALO2M_ll, IINFO_ll)
    IF ( .NOT. LSTEADYLS) CALL UPDATE_HALO2_ll(TZLSFIELD_ll, TZLSHALO2_ll, IINFO_ll)
  ENDIF
  CALL NUM_DIFF ( CLBCX, CLBCY, NRR, NSV, XDK2, XDK4,IMI,               &
                  XUM, XVM, XWM, XTHM, XTKEM, XRM, XSVM,                &
                  XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XRHODJ,               &
                  XRUS, XRVS, XRWS, XRTHS, XRTKES, XRRS, XRSVS,         &
                  LZDIFFU,TZHALO2M_ll, TZLSHALO2_ll,XZDIFFU_HALO2      )
END IF
!
DO JSV = NSV_CHEMBEG,NSV_CHEMEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_AERBEG,NSV_AEREND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_LNOXBEG,NSV_LNOXEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_DSTBEG,NSV_DSTEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_SLTBEG,NSV_SLTEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
!
CALL SECOND_MNH2(ZTIME2)
!
XT_DIFF = XT_DIFF + ZTIME2 - ZTIME1 &
          - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       14.    UPPER AND LATERAL RELAXATION
!               ----------------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF(LVE_RELAX .OR. LHORELAX_UVWTH .OR. LHORELAX_RV .OR.                   &
   LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI .OR. LHORELAX_RS .OR.   &
   LHORELAX_RG .OR. LHORELAX_RH .OR. LHORELAX_TKE .OR.                   &
   ANY(LHORELAX_SV)) THEN
  CALL RELAXATION (LVE_RELAX,LHORELAX_UVWTH,LHORELAX_RV,LHORELAX_RC,   &
                   LHORELAX_RR,LHORELAX_RI,LHORELAX_RS,LHORELAX_RG,    &
                   LHORELAX_RH,LHORELAX_TKE,LHORELAX_SV,               &
                   LHORELAX_SVC2R2,LHORELAX_SVC1R3,                    &
                   LHORELAX_SVELEC,LHORELAX_SVLG,                      &
                   LHORELAX_SVCHEM,LHORELAX_SVAER,                     &
                   LHORELAX_SVDST,LHORELAX_SVSLT,                      &
                   KTCOUNT,NRR,NSV,XTSTEP,XRHODJ,                      &
                   XUM, XVM, XWM, XTHM, XRM, XSVM, XTKEM,              &
                   XLSUM, XLSVM, XLSWM, XLSTHM,                        &
                   XLBXUM, XLBXVM, XLBXWM, XLBXTHM,                    &
                   XLBXRM, XLBXSVM, XLBXTKEM,                          &
                   XLBYUM, XLBYVM, XLBYWM, XLBYTHM,                    &
                   XLBYRM, XLBYSVM, XLBYTKEM,                          &
                   NALBOT, XALK, XALKW,                                &
                   LMASK_RELAX,XKURELAX, XKVRELAX, XKWRELAX,           &
                   NRIMX,NRIMY,                                        &
                   XRUS, XRVS, XRWS, XRTHS, XRRS, XRSVS, XRTKES        )
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_RELAX = XT_RELAX + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       15.    PARAMETRIZATIONS' MONITOR
!               -------------------------
!
ZTIME1 = ZTIME2
!
CALL PHYS_PARAM_n(KTCOUNT,ZTSTEP,YFMFILE,GCLOSE_OUT,              &
                  XT_RAD,XT_DCONV,XT_GROUND,XT_TURB,XT_CHEM,                &
                  ZTIME,GMASKkids)
!
IF (CDCONV/='NONE') THEN
  XPACCONV = XPACCONV + XPRCONV * XTSTEP
  IF (LCH_CONV_LINOX) THEN
    XIC_TOTAL_NUMBER = XIC_TOTAL_NUMBER + XIC_RATE * XTSTEP
    XCG_TOTAL_NUMBER = XCG_TOTAL_NUMBER + XCG_RATE * XTSTEP
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_PARAM = XT_PARAM + ZTIME2 - ZTIME1 - XTIME_LES - ZTIME
!
!-------------------------------------------------------------------------------
!
!*       16.    TEMPORAL SERIES
!               ---------------
!
ZTIME1 = ZTIME2
!
IF (LSERIES) THEN
  IF ( MOD (KTCOUNT-1,NFREQSERIES) == 0 ) CALL SERIES_n
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_MISC = XT_STEP_MISC + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       17.    LARGE SCALE FIELD REFRESH
!               -------------------------
!
ZTIME1 = ZTIME2
!
IF (.NOT. LSTEADYLS) THEN
  IF (  IMI==1                             .AND.      &
    NCPL_CUR < NCPL_NBR                              ) THEN
    IF (KTCOUNT+1 == NCPL_TIMES(NCPL_CUR,1)          ) THEN
                                  ! The next current time reachs a
      NCPL_CUR=NCPL_CUR+1         ! coupling one, LS sources are refreshed
      !
      CALL LS_COUPLING(CLUOUT,XTSTEP,GSTEADY_DMASS,                         &
             CGETTKEM,                                                      &
             CGETRVM,CGETRCM,CGETRRM,CGETRIM,                               &
             CGETRSM,CGETRGM,CGETRHM,CGETSVM,NSV,                           &
             NIMAX_ll,NJMAX_ll,                                             &
             NSIZELBX_ll,NSIZELBXU_ll,NSIZELBY_ll,NSIZELBYV_ll,             &
             NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                 &
             NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,         &
             XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XDRYMASST,                     &
             XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,          &
             XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,          &
             XLSUS,XLSVS,XLSWS,XLSTHS,XLSRVS,XDRYMASSS,                     &
             XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,          &
             XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS           )
      !
      DO JSV=NSV_CHEMBEG,NSV_CHEMEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_LNOXBEG,NSV_LNOXEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_AERBEG,NSV_AEREND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_DSTBEG,NSV_DSTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
     END IF
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_COUPL = XT_COUPL + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       18.    LATERAL BOUNDARY CONDITION FOR THE NORMAL VELOCITY
!               --------------------------------------------------
!
ZTIME1 = ZTIME2
!
CALL RAD_BOUND (CLBCX,CLBCY,XRIMKMAX,                    &
                ZTSTEP, XDXHAT, XDYHAT, XZZ ,            &
                XUM, XVM, XUT, XVT,                      &
                XLBXUM, XLBYVM, XLBXUS, XLBYVS,          &
                XCPHASE, XRHODJ,                         &
                XRUS, XRVS, XRWS                         )
!
CALL SECOND_MNH2(ZTIME2)
!
XT_RAD_BOUND = XT_RAD_BOUND + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       19.    PRESSURE COMPUTATION
!               --------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
!
IF(.NOT. L1D) THEN
  CALL PRESSUREZ ( GCLOSE_OUT, YFMFILE, CLUOUT,                           &
                  CLBCX,CLBCY,CPRESOPT,NITR,LITRADJ,KTCOUNT, XRELAX,IMI, &
                  XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY,XDXHATM,XDYHATM,XRHOM, &
                  XAF,XBFY,XCF,XTRIGSX,XTRIGSY,NIFAXX,NIFAXY,XPABSM,     &
                  NRR,NRRL,NRRI,XDRYMASST,XREFMASS,XMASS_O_PHI0,         &
                  XTHT,XRT,XRHODREF,XTHVREF,XRVREF,XEXNREF, XLINMASS,    &
                  XRUS, XRVS, XRWS, XPABST,                              &
                  XBFB,&
                  XBF_SXP2_YP1_Z) !JUAN Z_SPLITING
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_PRESS = XT_PRESS + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       20.    WATER MICROPHYSICS
!               ------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (CCLOUD /= 'NONE') THEN
  IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO' .OR. CCLOUD == 'C3R5' ) THEN
    IF ( LFORCING ) THEN
      ZWT_ACT_NUC(:,:,:) = XWT(:,:,:) + XWTFRC(:,:,:)
    ELSE
      ZWT_ACT_NUC(:,:,:) = XWT(:,:,:)
    END IF
    IF (CTURB /= 'NONE' ) THEN
      ZWT_ACT_NUC(:,:,:) = ZWT_ACT_NUC(:,:,:) + (2./3. * XTKET(:,:,:))**0.5
    ENDIF
  ELSE
    ZWT_ACT_NUC(:,:,:) = 0.
  END IF
!
  IF (CSURF=='EXTE') THEN
    ALLOCATE (ZSEA(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ALLOCATE (ZTOWN(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ZSEA(:,:) = 0.
    ZTOWN(:,:)= 0.
    CALL MNHGET_SURF_PARAM_n (PSEA=ZSEA(:,:),PTOWN=ZTOWN(:,:))
    CALL RESOLVED_CLOUD ( CCLOUD, NRR, NSPLITR, NSPLITG, IMI, KTCOUNT,         &
                          CLBCX,CLBCY,YFMFILE, CLUOUT, CRAD, CTURBDIM,         &
                          GCLOSE_OUT, LSUBG_COND,LSIGMAS,LSUBG_AUCV,ZTSTEP,    &
                          XZZ, XRHODJ, XRHODREF, XEXNREF,                      &
                          XPABSM, XPABST, XTHM, XTHT,XRM,XRT,XSIGS,XMFCONV,    &
                          ZWT_ACT_NUC, XRTHS, XRRS,                            &
                          XSVM, XSVT, XRSVS,                                   &
                          XSRCM, XCLDFR,XCIT,                                  &
                          LSEDIC,LACTIT, LSEDC, LSEDI, LRAIN, LWARM, LHHONI,   &
                          XINPRC,XINPRR, XINPRR3D, XEVAP3D,                    &
                          XINPRS, XINPRG, XINPRH, ZSEA, ZTOWN                  )
    DEALLOCATE(ZSEA)
    DEALLOCATE(ZTOWN)
  ELSE
    CALL RESOLVED_CLOUD ( CCLOUD, NRR, NSPLITR, NSPLITG, IMI, KTCOUNT,         &
                          CLBCX,CLBCY,YFMFILE, CLUOUT, CRAD, CTURBDIM,         &
                          GCLOSE_OUT, LSUBG_COND,LSIGMAS,LSUBG_AUCV,ZTSTEP,    &
                          XZZ, XRHODJ, XRHODREF, XEXNREF,                      &
                          XPABSM, XPABST, XTHM, XTHT,XRM,XRT,XSIGS,XMFCONV,    &
                          ZWT_ACT_NUC, XRTHS, XRRS,                            &
                          XSVM, XSVT, XRSVS,                                   &
                          XSRCM, XCLDFR,XCIT,                                  &
                          LSEDIC, LACTIT, LSEDC, LSEDI, LRAIN, LWARM, LHHONI,  &
                          XINPRC,XINPRR, XINPRR3D, XEVAP3D,                    &
                          XINPRS, XINPRG, XINPRH                               )
  END IF
!
  IF (CCLOUD /= 'REVE' ) THEN
    XACPRR = XACPRR + XINPRR * XTSTEP
    IF ((CCLOUD(1:3) == 'ICE' .AND. LSEDIC ) .OR.                       &
        ((CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO') &
                              .AND. LSEDC  )      )                     &
      XACPRC = XACPRC + XINPRC * XTSTEP
    IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5') THEN
      XACPRS = XACPRS + XINPRS * XTSTEP
      XACPRG = XACPRG + XINPRG * XTSTEP
      IF (CCLOUD == 'ICE4') XACPRH = XACPRH + XINPRH * XTSTEP          
    END IF
  END IF
!
! It is necessary that SV_C2R2 and SV_C1R3 are contiguous in the preceeding CALL
!
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_CLOUD = XT_CLOUD + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       21.    L.E.S. COMPUTATIONS
!               -------------------
!
ZTIME1 = ZTIME2
!
CALL LES_n
!
CALL SECOND_MNH2(ZTIME2)
!
XT_SPECTRA = XT_SPECTRA + ZTIME2 - ZTIME1 + XTIME_LES_BU + XTIME_LES
!
!-------------------------------------------------------------------------------
!
!*       22.    UPDATE HALO OF EACH SUBDOMAINS FOR TIME T+DT
!               --------------------------------------------
!
ZTIME1 = ZTIME2
!
CALL EXCHANGE (ZTSTEP,NRR,NSV,XRHODJ,TZFIELDS_ll,                      &
               XRUS, XRVS,XRWS,XRTHS,XRRS,XRTKES,XRSVS                 )
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_HALO = XT_HALO + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       23.    TEMPORAL SWAPPING
!               -----------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
!
CALL ENDSTEP  ( XTSTEP,NRR,NSV,KTCOUNT,IMI,XRHODJ,        &
                XRUS,XRVS,XRWS,XDRYMASSS,                 &
                XRTHS,XRRS,XRTKES,XRSVS,                  &
                XLSUS,XLSVS,XLSWS,                        &
                XLSTHS,XLSRVS,                            &
                XLBXUS,XLBXVS,XLBXWS,                     &
                XLBXTHS,XLBXRS,XLBXTKES,XLBXSVS,          &
                XLBYUS,XLBYVS,XLBYWS,                     &
                XLBYTHS,XLBYRS,XLBYTKES,XLBYSVS,          &
                XUM,XVM,XWM,XPABSM,                       &
                XTHM,XRM,XTKEM,XSVM,XSRCM,                &
                XUT,XVT,XWT,XPABST,XDRYMASST,             &
                XTHT, XRT, XTKET, XSVT,XSRCT,             &
                XLSUM,XLSVM,XLSWM,                        &
                XLSTHM,XLSRVM,                            &
                XLBXUM,XLBXVM,XLBXWM,                     &
                XLBXTHM,XLBXRM,XLBXTKEM,XLBXSVM,          &
                XLBYUM,XLBYVM,XLBYWM,                     &
                XLBYTHM,XLBYRM,XLBYTKEM,XLBYSVM,          &
                CMET_ADV_SCHEME, CSV_ADV_SCHEME           )
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_SWA = XT_STEP_SWA + ZTIME2 - ZTIME1 - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       24.1    BALLOON and AIRCRAFT
!               --------------------
!
ZTIME1 = ZTIME2
!
IF (LFLYER)                                                                   &
  CALL AIRCRAFT_BALLOON(CLUOUT, XTSTEP,                                       &
                      TDTEXP, TDTMOD, TDTSEG, TDTCUR,                         &
                      XXHAT, XYHAT, XZZ, XMAP, XLONORI, XLATORI,              &
                      XUT, XVT, XWT, XPABST, XTHT, XRT, XSVT, XTKET, XTSRAD   )


!-------------------------------------------------------------------------------
!
!*       24.2    STATION (observation diagnostic)
!               --------------------------------
!
IF (LSTATION)                                                            &
  CALL STATION_n(CLUOUT, XTSTEP,                                         &
                 TDTEXP, TDTMOD, TDTSEG, TDTCUR,                         &
                 XXHAT, XYHAT, XZZ,                                      &
                 XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET, XTSRAD, XPABST   )
!
!---------------------------------------------------------
!
!*       24.3    PROFILER (observation diagnostic)
!               ---------------------------------
!
IF (LPROFILER)                                                           &
  CALL PROFILER_n(CLUOUT, XTSTEP,                                        &
                  TDTEXP, TDTMOD, TDTSEG, TDTCUR,                        &
                  XXHAT, XYHAT, XZZ,                                     &
                  XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET, XTSRAD, XPABST  )
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_MISC = XT_STEP_MISC + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       24.4   deallocation of observation diagnostics
!               ---------------------------------------
!
CALL END_DIAG_IN_RUN
!
!-------------------------------------------------------------------------------
!
!
!*       25.    STORAGE OF BUDGET FIELDS
!               ------------------------
!
ZTIME1 = ZTIME2
!
IF (NBUMOD==IMI .AND. CBUTYPE/='NONE') THEN
  CALL ENDSTEP_BUDGET(CFMDIAC,CLUOUT,KTCOUNT,TDTCUR,TDTMOD,XTSTEP,NSV)
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STEP_BUD = XT_STEP_BUD + ZTIME2 - ZTIME1 + XTIME_BU
!
!-------------------------------------------------------------------------------
!
!*       26.    FM FILE CLOSURE
!               ---------------
!
IF (GCLOSE_OUT) THEN
  GCLOSE_OUT=.FALSE.
  CALL FMCLOS_ll(YFMFILE,'KEEP',CLUOUT,IRESP)
END IF
!
!-------------------------------------------------------------------------------
!
!*       27.    CURRENT TIME REFRESH
!               --------------------
!
TDTCUR%TIME=TDTCUR%TIME + XTSTEP
CALL ADD_FORECAST_TO_DATE(TDTCUR%TDATE%YEAR, &
                          TDTCUR%TDATE%MONTH,&
                          TDTCUR%TDATE%DAY,  &
                          TDTCUR%TIME        )
!
!-------------------------------------------------------------------------------
!
!*       28.    CPU ANALYSIS
!               ------------
!
CALL SECOND_MNH2(ZTIME2)
XT_START=XT_START+ZTIME2-ZEND
!
!
IF ( KTCOUNT == NSTOP .AND. IMI==1) THEN
  OEXIT=.TRUE.  
END IF
!
IF (OEXIT) THEN
!
  IF (LSERIES) CALL WRITE_SERIES_n(CFMDIAC,CLUOUT )
  CALL WRITE_AIRCRAFT_BALLOON(CFMDIAC)
  CALL WRITE_STATION_n(CFMDIAC)
  CALL WRITE_PROFILER_n(CFMDIAC)
  CALL WRITE_LES_n(' ')
  CALL WRITE_LES_n('A')
  CALL WRITE_LES_n('E')
  CALL WRITE_LES_n('H')
  CALL MENU_DIACHRO(CFMDIAC,CLUOUT,'END')
  CALL FMCLOS_ll(CFMDIAC,'KEEP',CLUOUT,IRESP)
  !
  CALL FMCLOS_ll(CINIFILE,'KEEP',CLUOUT,IRESP)
!
!*       28.1   print statistics!
!
  ! Set File Timing OUTPUT
  !
  CALL SET_ILUOUT_TIMING(ILUOUT)
  !
  ! Compute global time
  !
  CALL TIME_STAT_ll(XT_START,ZTOT)
  !
  CALL TIME_HEADER_ll(IMI)
  !
  CALL TIME_STAT_ll(XT_1WAY,ZTOT,       ' ONE WAY','=')
  CALL TIME_STAT_ll(XT_BOUND,ZTOT,      ' BOUNDARIES','=')
  CALL TIME_STAT_ll(XT_STORE,ZTOT,      ' STORE-FIELDS','=')
  CALL TIME_STAT_ll(XT_GUESS,ZTOT,      ' INITIAL_GUESS','=')
  CALL TIME_STAT_ll(XT_2WAY,ZTOT,       ' TWO WAY','=')
  CALL TIME_STAT_ll(XT_ADV,ZTOT,        ' ADVECTION','=')
  CALL TIME_STAT_ll(XT_FORCING,ZTOT,    ' FORCING','=')
  CALL TIME_STAT_ll(XT_NUDGING,ZTOT,    ' NUDGING','=')
  CALL TIME_STAT_ll(XT_SOURCES,ZTOT,    ' DYN_SOURCES','=')
  CALL TIME_STAT_ll(XT_DIFF,ZTOT,       ' NUM_DIFF','=')
  CALL TIME_STAT_ll(XT_RELAX,ZTOT,      ' RELAXATION','=')
  !
  CALL  TIMING_LEGEND() 
  !
  CALL TIME_STAT_ll(XT_PARAM,ZTOT,      ' PHYS_PARAM','=')
    CALL TIME_STAT_ll(XT_RAD,ZTOT,      '   RAD       = '//CRAD  ,'-')
    CALL TIME_STAT_ll(XT_DCONV,ZTOT,    '   DEEP CONV = '//CDCONV,'-')
    CALL TIME_STAT_ll(XT_GROUND,ZTOT,   '   GROUND'              ,'-')
    CALL TIME_STAT_ll(XT_TURB,ZTOT,     '   TURB      = '//CTURB ,'-')
    CALL TIME_STAT_ll(XT_CHEM,ZTOT,     '   CHIMIE'              ,'-')
  CALL  TIMING_LEGEND()
  CALL TIME_STAT_ll(XT_COUPL,ZTOT,      ' SET_COUPLING','=')
  CALL TIME_STAT_ll(XT_RAD_BOUND,ZTOT,  ' RAD_BOUND','=')
  !
  CALL  TIMING_LEGEND()
  ! 
  CALL TIME_STAT_ll(XT_PRESS,ZTOT,      ' PRESSURE ','=')
  !JUAN Z_SPLITTING
    CALL TIME_STAT_ll(T_MAP_B_SX_YP2_ZP1,ZTOT,          '   REMAP      B=>FFTX'  ,'-')
    CALL TIME_STAT_ll(T_MAP_SX_YP2_ZP1_SXP2_Y_ZP1,ZTOT, '   REMAP   FFTX=>FFTY'  ,'-')
    CALL TIME_STAT_ll(T_MAP_SXP2_Y_ZP1_B,ZTOT,          '   REMAP   FTTY=>B'     ,'-')
    CALL TIME_STAT_ll(T_MAP_SXP2_Y_ZP1_SXP2_YP1_Z,ZTOT, '   REMAP   FFTY=>SUBZ'  ,'-')
    CALL TIME_STAT_ll(T_MAP_B_SXP2_Y_ZP1,ZTOT,          '   REMAP      B=>FFTY-1','-')
    CALL TIME_STAT_ll(T_MAP_SXP2_YP1_Z_SXP2_Y_ZP1,ZTOT, '   REMAP   SUBZ=>FFTY-1','-')
    CALL TIME_STAT_ll(T_MAP_SXP2_Y_ZP1_SX_YP2_ZP1,ZTOT, '   REMAP FFTY-1=>FFTX-1','-')
    CALL TIME_STAT_ll(T_MAP_SX_YP2_ZP1_B,ZTOT,          '   REMAP FFTX-1=>B     ','-')
  ! JUAN P1/P2
  CALL TIME_STAT_ll(XT_CLOUD,ZTOT,      ' RESOLVED_CLOUD','=')
  CALL TIME_STAT_ll(XT_HALO,ZTOT,       ' EXCHANGE_HALO','=')
  CALL TIME_STAT_ll(XT_STEP_SWA,ZTOT,   ' ENDSTEP','=')
  CALL TIME_STAT_ll(XT_STEP_BUD,ZTOT,   ' BUDGETS','=')
  CALL TIME_STAT_ll(XT_SPECTRA,ZTOT,    ' LES','=')
  CALL TIME_STAT_ll(XT_STEP_MISC,ZTOT,  ' MISCELLANEOUS','=')
  !
  ! sum of call subroutine
  !
  ZALL   = XT_1WAY + XT_BOUND   + XT_STORE   + XT_GUESS    +  XT_2WAY   + &
           XT_ADV  + XT_FORCING + XT_NUDGING + XT_SOURCES  +  XT_DIFF   + &
           XT_RELAX+ XT_PARAM   + XT_COUPL   + XT_RAD_BOUND+XT_PRESS    + &
           XT_CLOUD+  XT_HALO   + XT_SPECTRA + XT_STEP_SWA +XT_STEP_MISC+ &
           XT_STEP_BUD
  CALL TIME_STAT_ll(ZALL,ZTOT,          ' SUM(CALL)','=')
  CALL  TIMING_SEPARATOR('=')
  !
  ! Gobale Stat
  !
  WRITE(ILUOUT,FMT=*)
  WRITE(ILUOUT,FMT=*)
  CALL  TIMING_LEGEND() 
  !
  ! MODELN all included
  !
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')  
  WRITE(YMI,FMT="(I0)") IMI
  CALL TIME_STAT_ll(XT_START,ZTOT,      ' MODEL'//YMI,'+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
  !
  ! Timing/ Steps
  !
  ZTIME_STEP     =  XT_START / FLOAT(KTCOUNT)
  WRITE(YTCOUNT,FMT="(I0)") KTCOUNT
  CALL TIME_STAT_ll(ZTIME_STEP,ZTOT,     ' SECOND/STEP='//YTCOUNT,'=')
  !
  ! Timing/Step/Points
  !
  IPOINTS = NIMAX_ll*NJMAX_ll*NKMAX
  WRITE(YPOINTS,FMT="(I0)") IPOINTS
  ZTIME_STEP_PTS =  ZTIME_STEP / FLOAT(IPOINTS) * 1e6
  CALL TIME_STAT_ll(ZTIME_STEP_PTS,ZTOT_PT)
  CALL TIME_STAT_ll(ZTIME_STEP_PTS,ZTOT_PT,  ' MICROSEC/STP/PT='//YPOINTS,'-')
  !
  CALL  TIMING_SEPARATOR('=')
  !
  !
  !
END IF



END SUBROUTINE MODEL_n
