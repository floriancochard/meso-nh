!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_MODEL_n
!     ###################
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
!!      Subroutine IO_FILE_OPEN_ll: to open a file
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
!!      Subroutine IO_FILE_CLOSE_ll: closes a file
!!      Subroutine DATETIME_CORRECTDATE: transform the current time in GMT
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
!!                   Fev. 2006 (M. Leriche) add aqueous phase chemistry
!!                   april 2006 (T.Maric) Add halo related to 4th order advection scheme
!!                   May 2006 Remove KEPS
!!                   Oct 2008 (C.Lac) FIT for variables advected with PPM
!!                   July 2009 : Displacement of surface diagnostics call to be
!!                               coherent with  surface diagnostics obtained with DIAG
!!                   10/11/2009 (P. Aumond) Add mean moments
!!                   Nov, 12, 2009 (C. Barthe) add cloud electrification and lightning flashes
!!                   July 2010 (M. Leriche) add ice phase chemical species
!!                   April 2011 (C.Lac) : Remove instant M 
!!                   April 2011 (C.Lac, V.Masson) : Time splitting for advection
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!       P. Tulet      Nov 2014 accumulated moles of aqueous species that fall at the surface   
!!                   Dec 2014 (C.Lac) : For reproducibility START/RESTA
!!      J.Escobar 20/04/2015: missing UPDATE_HALO before UPDATE_HALO2
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!       C.Lac    11/09/2015: correction of the budget due to FIT temporal scheme
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!                   Sep 2015 (S. Bielli) : Remove YDADFILE from argument call 
!                              of write_phys_param
!!      J.Escobar : 19/04/2016 : Pb IOZ/NETCDF , missing OPARALLELIO=.FALSE. for PGD files
!!      M.Mazoyer : 04/2016      DTHRAD used for radiative cooling when LACTIT
!!!      Modification    01/2016  (JP Pinty) Add LIMA
!!  06/2016     (G.Delautier) phasage surfex 8
!!      M.Leriche : 03/2016 Move computation of accumulated chem. in rain to ch_monitor
!!                  09/2016 Add filter on negative values on AERDEP SV before relaxation
!!                  10/2016  (C.Lac) _ Correction on the flag for Strang splitting
!!                                  to insure reproducibility between START and RESTA
!!                                  _  Add OSPLIT_WENO
!!                                  _ Add droplet deposition 
!!                   10/2016 (M.Mazoyer) New KHKO output fields
!!      P.Wautelet : 11/07/2016 : removed MNH_NCWRIT define
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!                   10/2017 (C.Lac) Necessity to have chemistry processes as
!!                            the las process modifying XRSVS
!!  01/2018      (G.Delautier) SURFEX 8.1
!!  03/2018     (P.Wautelet)   replace ADD_FORECAST_TO_DATE by DATETIME_CORRECTDATE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                  07/2017  (V. Vionnet) : Add blowing snow scheme 
!!      S. Riette : 11/2016 Add ZPABST to keep pressure constant during timestep
!!                   01/2018 (C.Lac) Add VISCOSITY
!!  Philippe Wautelet: 21/01/2019: add LIO_ALLOW_NO_BACKUP and LIO_NO_WRITE to modd_io_ll
!                                  to allow to disable writes (for bench purposes)
!!                   02/2019 C.Lac add rain fraction as an output field
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_2D_FRC
USE MODD_ADV_n
USE MODD_AIRCRAFT_BALLOON
USE MODD_BAKOUT
USE MODD_BIKHARDT_n
USE MODD_BLANK 
USE MODD_BUDGET
USE MODD_CH_AERO_n,      ONLY: XSOLORG, XMI
USE MODD_CH_MNHC_n,      ONLY: LUSECHEM,LCH_CONV_LINOX,LUSECHAQ,LUSECHIC, &
                               LCH_INIT_FIELD
USE MODD_CLOUD_MF_n  
USE MODD_VISCOSITY
USE MODD_DRAG_n
USE MODD_CLOUDPAR_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST, ONLY: XMD
USE MODD_CURVCOR_n
USE MODD_DEEP_CONVECTION_n
USE MODD_DIM_n
USE MODD_DUST,           ONLY: LDUST
USE MODD_DYN
USE MODD_DYN_n
USE MODD_DYNZD
USE MODD_DYNZD_n
USE MODD_ELEC_DESCR
USE MODD_FIELD_n
USE MODD_FRC
USE MODD_FRC_n
USE MODD_GET_n
USE MODD_GRID,           ONLY: XLONORI,XLATORI
USE MODD_GRID_n
USE MODD_ICE_C1R3_DESCR, ONLY: XRTMIN_C1R3=>XRTMIN
USE MODD_IO_ll,          ONLY: LIO_NO_WRITE, TFILEDATA,TFILE_SURFEX,TFILE_DUMMY
USE MODD_LBC_n
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_LSFIELD_n
USE MODD_LUNIT,          ONLY: TLUOUT0,TOUTDATAFILE
USE MODD_LUNIT_n,        ONLY: TDIAFILE,TINIFILE,TINIFILEPGD,TLUOUT
USE MODD_MEAN_FIELD
USE MODD_MEAN_FIELD_n
USE MODD_METRICS_n
USE MODD_MNH_SURFEX_n
USE MODD_NESTING
USE MODD_NSV
USE MODD_NUDGING_n
USE MODD_OUT_n
USE MODD_PARAM_C1R3,     ONLY: NSEDI => LSEDI, NHHONI => LHHONI
USE MODD_PARAM_C2R2,     ONLY: NSEDC => LSEDC, NRAIN => LRAIN, NACTIT => LACTIT,LACTTKE,LDEPOC
USE MODD_PARAMETERS
USE MODD_PARAM_ICE,      ONLY: LWARM,LSEDIC,LCONVHG,LDEPOSC
USE MODD_PARAM_LIMA,     ONLY: MSEDC => LSEDC, MWARM => LWARM, MRAIN => LRAIN, LACTI,  &
                               MACTIT => LACTIT, LSCAV, NMOD_CCN, LCOLD,               &
                               MSEDI => LSEDI, MHHONI => LHHONI, NMOD_IFN, LHAIL,      &
                               XRTMIN_LIMA=>XRTMIN, MACTTKE=>LACTTKE
USE MODD_BLOWSNOW_n
USE MODD_BLOWSNOW                       
USE MODD_PARAM_MFSHALL_n
USE MODD_PARAM_n
USE MODD_PAST_FIELD_n
USE MODD_PRECIP_n
USE MODD_PROFILER_n
USE MODD_RADIATIONS_n,   ONLY: XTSRAD,XSCAFLASWD,XDIRFLASWD,XDIRSRFSWD, XAER, XDTHRAD
USE MODD_RAIN_ICE_DESCR, ONLY: XRTMIN
USE MODD_REF_n
USE MODD_SALT,           ONLY: LSALT
USE MODD_SERIES,         ONLY: LSERIES
USE MODD_SERIES_n,       ONLY: NFREQSERIES
USE MODD_STATION_n
USE MODD_SUB_MODEL_n
USE MODD_TIME
USE MODD_TIME_n 
USE MODD_TIMEZ
USE MODD_TURB_CLOUD,     ONLY: NMODEL_CLOUD,CTURBLEN_CLOUD,XCEI
USE MODD_TURB_n
!
USE MODE_DATETIME
USE MODE_ELEC_ll
USE MODE_FM
USE MODE_GRIDCART         
USE MODE_GRIDPROJ
USE MODE_IO_ll
USE MODE_IO_WRITE_FIELD
USE MODE_ll
USE MODE_MNH_TIMING
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
!
USE MODI_ADVECTION_METSV
USE MODI_ADVECTION_UVW     
USE MODI_ADVECTION_UVW_CEN 
USE MODI_ADV_FORCING_n
USE MODI_AER_MONITOR_n
USE MODI_AIRCRAFT_BALLOON
USE MODI_BOUNDARIES
USE MODI_BUDGET_FLAGS
USE MODI_CART_COMPRESS
USE MODI_CH_MONITOR_n
USE MODI_DIAG_SURF_ATM_N
USE MODI_DYN_SOURCES
USE MODI_END_DIAG_IN_RUN
USE MODI_ENDSTEP
USE MODI_ENDSTEP_BUDGET
USE MODI_EXCHANGE
USE MODI_FORCING
USE MODI_FORC_SQUALL_LINE
USE MODI_FORC_WIND
USE MODI_GET_HALO
USE MODI_GRAVITY_IMPL         
USE MODI_INI_DIAG_IN_RUN
USE MODI_INI_LG
USE MODI_INI_MEAN_FIELD
USE MODI_INITIAL_GUESS
USE MODI_LES_INI_TIMESTEP_n
USE MODI_LES_N
USE MODI_VISCOSITY
USE MODI_LIMA_PRECIP_SCAVENGING
USE MODI_LS_COUPLING
USE MODI_MASK_COMPRESS
USE MODI_MEAN_FIELD
USE MODI_MENU_DIACHRO
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_MNHWRITE_ZS_DUMMY_n
USE MODI_NUDGING
USE MODI_NUM_DIFF
USE MODI_ONE_WAY_n
USE MODI_PHYS_PARAM_n
USE MODI_PRESSUREZ
USE MODI_PROFILER_n
USE MODI_RAD_BOUND
USE MODI_RELAX2FW_ION 
USE MODI_RELAXATION
USE MODI_REL_FORCING_n
USE MODI_RESOLVED_CLOUD
USE MODI_RESOLVED_ELEC_n
USE MODI_SERIES_N
USE MODI_SETLB_LG
USE MODI_SET_MASK
USE MODI_SHUMAN
USE MODI_SPAWN_LS_n
USE MODI_STATION_n
USE MODI_TURB_CLOUD_INDEX
USE MODI_TWO_WAY
USE MODI_UPDATE_NSV
USE MODI_WRITE_AIRCRAFT_BALLOON
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_DIAG_SURF_ATM_N
USE MODI_WRITE_LES_n
USE MODI_WRITE_LFIFM_n
USE MODI_WRITE_LFIFMN_FORDIACHRO_n
USE MODI_WRITE_PROFILER_n
USE MODI_WRITE_SERIES_n
USE MODI_WRITE_STATION_n
USE MODI_BLOWSNOW

USE MODI_WRITE_SURF_ATM_N
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
INTEGER :: IIB,IIE,IJB,IJE,IKB,IKE ! index values for the physical subdomain
INTEGER :: JSV,JRR     ! Loop index for scalar and moist variables
INTEGER  :: INBVAR              ! number of HALO2_lls to allocate
INTEGER  :: IRESP               ! return code in FM routines
INTEGER  :: IINFO_ll            ! return code of parallel routine
INTEGER :: IVERB                ! LFI verbosity level
LOGICAL :: GSTEADY_DMASS        ! conditional call to mass computation
!
                                ! for computing time analysis
REAL*8,DIMENSION(2)         :: ZTIME,ZTIME1,ZTIME2,ZEND,ZTOT,ZALL,ZTOT_PT
!
REAL*8,DIMENSION(2)         :: ZTIME_STEP,ZTIME_STEP_PTS
CHARACTER                 :: YMI
INTEGER                   :: IPOINTS
CHARACTER(len=16)         :: YTCOUNT,YPOINTS

REAL         :: ZSTAT_CSTORE,ZSTAT_CBOUND,ZSTAT_CGUESS,ZSTAT_CADV,ZSTAT_CSOURCES
REAL         :: ZSTAT_CDIFF,ZSTAT_CRELAX,ZSTAT_CPARAM
REAL         :: ZSTAT_CSPECTRA,ZSTAT_CRAD_BOUND,ZSTAT_CPRESS
REAL         :: ZSTAT_CCLOUD,ZSTAT_CSTEP_SWA,ZSTAT_CSTEP_MISC
REAL         :: ZSTAT_CCOUPL,ZSTAT_CSTEP_BUD,ZSTAT_CSTEP_CDRAG
REAL         :: ZSTAT_CSTEP_CTRACER,ZSTAT_CSTEP_CELEC
REAL         :: SCONV_CTURB,ZSTAT_C1WAY,ZSTAT_C2WAY,ZSTAT_CMAFL
REAL         :: ZSTAT_CRAD,ZSTAT_CDCONV,ZSTAT_CGROUND,ZSTAT_CHALO
REAL         :: ZSTAT_CFORCING,ZSTAT_CNUDGING,ZSTAT_CCHEM
!
REAL         :: ZPERCALL,ZPRICE
REAL         :: ZPERCSTORE,ZPERCBOUND,ZPERCGUESS,ZPERCADV,ZPERCSOURCES,ZPERCDRAG
REAL         :: ZPERCDIFF,ZPERCRELAX,ZPERCPARAM
REAL         :: ZPERCSPECTRA,ZPERCRAD_BOUND,ZPERCPRESS
REAL         :: ZPERCCLOUD,ZPERCSTEP_SWA,ZPERCSTEP_MISC
REAL         :: ZPERCELEC
REAL         :: ZPERCCOUPL,ZPERCSTEP_BUD
REAL         :: ZPERCTURB,ZPERC1WAY,ZPERC2WAY
REAL         :: ZPERCRAD,ZPERCSHADOWS,ZPERCKAFR,ZPERCGROUND,ZPERCHALO,ZPERCMAFL,ZPERTRACER
REAL         :: ZPERCFORCING,ZPERCNUDGING,ZPERCCHEM
REAL         :: ZTSTEP_UVW   ! Double timestep except for cold start (single)
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
REAL, DIMENSION(:,:),   POINTER :: DPTR_XLSZWSM,DPTR_XLSZWSS
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
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDC
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDR
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDS
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDG
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZSPEEDH
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRC3D
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRS3D
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRG3D
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))           :: ZINPRH3D
!
LOGICAL :: KWARM        
LOGICAL :: KRAIN        
LOGICAL :: KSEDC  
LOGICAL :: KACTIT
LOGICAL :: KSEDI
LOGICAL :: KHHONI
REAL :: TEMPS
INTEGER :: NSV_END
CHARACTER (LEN=100) :: YCOMMENT   ! Comment string in LFIFM file
CHARACTER (LEN=LEN_HREC)  :: YRECFM     ! Name of the desired field in LFIFM file
!
INTEGER             :: ILENG      ! Length of comment string in LFIFM file
INTEGER             :: IGRID      ! C-grid indicator in LFIFM file
INTEGER             :: ILENCH     ! Length of comment string in LFIFM file
!
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZRUS,ZRVS,ZRWS
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZPABST !To give pressure at t 
                                                     ! (and not t+1) to resolved_cloud
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZJ
!
! for various testing
INTEGER :: IK
REAL, DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZTMP
!
TYPE(LIST_ll), POINTER :: TZFIELDC_ll   ! list of fields to exchange
TYPE(HALO2LIST_ll), POINTER :: TZHALO2C_ll   ! list of fields to exchange
LOGICAL :: GCLD                     ! conditionnal call for dust wet deposition
LOGICAL :: GCLOUD_ONLY              ! conditionnal radiation computations for
                                !      the only cloudy columns
REAL, DIMENSION(SIZE(XRSVS,1), SIZE(XRSVS,2), SIZE(XRSVS,3), NSV_AER)  :: ZWETDEPAER


!
TYPE(TFILEDATA),POINTER :: TZBAKFILE, TZOUTFILE
! TYPE(TFILEDATA),SAVE    :: TZDIACFILE
!-------------------------------------------------------------------------------
!
TZBAKFILE=> NULL()
TZOUTFILE=> NULL()
!
!*       0.    MICROPHYSICAL SCHEME
!              ------------------- 
SELECT CASE(CCLOUD)
CASE('C2R2','KHKO','C3R5')
  KWARM  = .TRUE.        
  KRAIN  = NRAIN
  KSEDC  = NSEDC
  KACTIT = NACTIT
!
  KSEDI  = NSEDI
  KHHONI = NHHONI
CASE('LIMA')
  KWARM  = MWARM        
  KRAIN  = MRAIN
  KSEDC  = MSEDC
  KACTIT = MACTIT
!
  KSEDI  = MSEDI
  KHHONI = MHHONI
CASE('ICE3','ICE4') !default values
  KWARM  = LWARM        
  KRAIN  = .TRUE.
  KSEDC  = .TRUE.
  KACTIT = .FALSE.
!
  KSEDI  = .TRUE.
  KHHONI = .FALSE.
END SELECT
!
!
!*        1    PRELIMINARY
!              ------------
IMI = GET_CURRENT_MODEL_INDEX()
!
!*       1.0   update NSV_* variables for current model
!              ----------------------------------------
!
CALL UPDATE_NSV(IMI)
!
!*       1.1   RECOVER THE LOGICAL UNIT NUMBER FOR THE OUTPUT PRINTS
!
ILUOUT = TLUOUT%NLU
!
!*       1.2   SET ARRAY SIZE
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=NKMAX+2*JPVEXT
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
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
  NULLIFY(TZFIELDS_ll,TZLSFIELD_ll,TZFIELDT_ll)
  NULLIFY(TZHALO2T_ll)
  NULLIFY(TZLSHALO2_ll)
  NULLIFY(TZFIELDSC_ll)
!
  ALLOCATE(ZWT_ACT_NUC(SIZE(XWT,1),SIZE(XWT,2),SIZE(XWT,3)))
  ALLOCATE(GMASKkids(SIZE(XWT,1),SIZE(XWT,2)))
!
! initialization of the FM file backup/output number
  IBAK=0
  IOUT=0
!
  IF ( .NOT. LIO_NO_WRITE ) THEN
    CALL IO_FILE_OPEN_ll(TDIAFILE)
!
    CALL IO_WRITE_HEADER(TDIAFILE)
    CALL WRITE_DESFM_n(IMI,TDIAFILE)
    CALL WRITE_LFIFMN_FORDIACHRO_n(TDIAFILE)
  END IF
!
!*       1.4   Initialization of the list of fields for the halo updates
!
!                 a) Sources terms
!
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRUS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRVS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRWS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRTHS)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRUS_PRES)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRVS_PRES)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRWS_PRES)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XRTHS_CLD)
  IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDS_ll, XRTKES)
  DO JRR=1,NRR
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XRRS(:,:,:,JRR))
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XRRS_CLD(:,:,:,JRR))
  ENDDO
  DO JSV=1,NSV
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XRSVS(:,:,:,JSV))
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XRSVS_CLD(:,:,:,JSV))
  ENDDO
  IF (SIZE(XSRCT,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDS_ll, XSRCT)
  !
  IF ((LNUMDIFU .OR. LNUMDIFTH .OR. LNUMDIFSV) ) THEN
  !
  !                 b) LS fields
  !
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSUM)
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSVM)
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSWM)
    CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSTHM)
    CALL ADD2DFIELD_ll(TZLSFIELD_ll, XLSZWSM)
    IF (NRR >= 1) THEN
      CALL ADD3DFIELD_ll(TZLSFIELD_ll, XLSRVM)
    ENDIF
  !
  !                 c) Fields at t
  !
    CALL ADD3DFIELD_ll(TZFIELDT_ll, XUT)
    CALL ADD3DFIELD_ll(TZFIELDT_ll, XVT)
    CALL ADD3DFIELD_ll(TZFIELDT_ll, XWT)
    CALL ADD3DFIELD_ll(TZFIELDT_ll, XTHT)
    IF (SIZE(XRTKES,1) /= 0) CALL ADD3DFIELD_ll(TZFIELDT_ll, XTKET)
    DO JRR=1,NRR
      CALL ADD3DFIELD_ll(TZFIELDT_ll, XRT(:,:,:,JRR))
    ENDDO
    DO JSV=1,NSV
      CALL ADD3DFIELD_ll(TZFIELDT_ll, XSVT(:,:,:,JSV))
    ENDDO
  !
  !*       1.5   Initialize the list of fields for the halo updates (2nd layer)
  !
    INBVAR = 4+NRR+NSV
    IF (SIZE(XRTKES,1) /= 0) INBVAR=INBVAR+1
    CALL INIT_HALO2_ll(TZHALO2T_ll,INBVAR,IIU,IJU,IKU)
    CALL INIT_HALO2_ll(TZLSHALO2_ll,4+MIN(1,NRR),IIU,IJU,IKU)
  !
  !*       1.6   Initialise the 2nd layer of the halo of the LS fields
  !
    IF ( LSTEADYLS ) THEN
       CALL UPDATE_HALO_ll(TZLSFIELD_ll, IINFO_ll)
       CALL DEL2DFIELD_ll(TZLSFIELD_ll,XLSZWSM,IINFO_ll) 
       CALL UPDATE_HALO2_ll(TZLSFIELD_ll, TZLSHALO2_ll, IINFO_ll)
    END IF
  END IF
  !
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
  XT_ADVUVW    = 0.0
  XT_GRAV      = 0.0
  XT_SOURCES   = 0.0
  !
  XT_DIFF      = 0.0
  XT_RELAX     = 0.0
  XT_PARAM     = 0.0
  XT_SPECTRA   = 0.0
  XT_HALO      = 0.0
  XT_VISC      = 0.0  
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
  XT_MAFL      = 0.0
  XT_DRAG      = 0.0
  XT_TRACER    = 0.0
  XT_SHADOWS   = 0.0
  XT_ELEC      = 0.0
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
  DPTR_XLSZWSM=>XLSZWSM
  DPTR_XLSTHS=>XLSTHS
  DPTR_XLSRVS=>XLSRVS
  DPTR_XLSUS=>XLSUS
  DPTR_XLSVS=>XLSVS
  DPTR_XLSWS=>XLSWS
  DPTR_XLSZWSS=>XLSZWSS
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
             DPTR_XLSTHM,DPTR_XLSRVM,DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSZWSM,                        &
             DPTR_XLSTHS,DPTR_XLSRVS,DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS, DPTR_XLSZWSS                         )
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
  CALL ONE_WAY_n(NDAD(IMI),XTSTEP,IMI,KTCOUNT,                   &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
       NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),NDTRATIO(IMI),         &
       DPTR_CLBCX,DPTR_CLBCY,NRIMX,NRIMY,                                &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,      &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,      &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,      &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM,      &
       GSTEADY_DMASS,CCLOUD,LUSECHAQ,LUSECHIC,                           &
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
IF (CCONF == "START" .OR. (CCONF == "RESTA" .AND. KTCOUNT /= 1 )) THEN
CALL MPPDB_CHECK3DM("before BOUNDARIES:XUT, XVT, XWT, XTHT, XTKET",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET)
CALL BOUNDARIES (                                                   &
            XTSTEP,CLBCX,CLBCY,NRR,NSV,KTCOUNT,                     &
            XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
            XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
            XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,   &
            XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS,   &
            XRHODJ,                                                 &
            XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT, XSRCT            )
CALL MPPDB_CHECK3DM("after  BOUNDARIES:XUT, XVT, XWT, XTHT, XTKET",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET)
END IF
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
IF (IBAK < NBAK_NUMB ) THEN
  IF (KTCOUNT == TBACKUPN(IBAK+1)%NSTEP) THEN
    IBAK=IBAK+1
    GCLOSE_OUT=.TRUE.
    !
    TZBAKFILE => TBACKUPN(IBAK)%TFILE
    IVERB    = TZBAKFILE%NLFIVERB
    !
    CALL IO_FILE_OPEN_ll(TZBAKFILE)
    !
    CALL WRITE_DESFM_n(IMI,TZBAKFILE)
    CALL IO_WRITE_HEADER(TBACKUPN(IBAK)%TFILE)
    CALL WRITE_LFIFM_n(TBACKUPN(IBAK)%TFILE,TBACKUPN(IBAK)%TFILE%TDADFILE%CNAME)
    TOUTDATAFILE => TZBAKFILE
    CALL MNHWRITE_ZS_DUMMY_n(TZBAKFILE)
    IF (CSURF=='EXTE') THEN
      TFILE_SURFEX => TZBAKFILE
      CALL GOTO_SURFEX(IMI)
      CALL WRITE_SURF_ATM_n(YSURF_CUR,'MESONH','ALL',.FALSE.)
      NULLIFY(TFILE_SURFEX)
    END IF
    !
    ! Reinitialize Lagragian variables at every model backup
    IF (LLG .AND. LINIT_LG .AND. CINIT_LG=='FMOUT') THEN
      CALL INI_LG(XXHAT,XYHAT,XZZ,XSVT,XLBXSVM,XLBYSVM)
      IF (IVERB>=5) THEN
        WRITE(UNIT=ILUOUT,FMT=*) '************************************'
        WRITE(UNIT=ILUOUT,FMT=*) '*** Lagrangian variables refreshed after ',TRIM(TZBAKFILE%CNAME),' backup'
        WRITE(UNIT=ILUOUT,FMT=*) '************************************'
      END IF
    END IF
    ! Reinitialise mean variables
    IF (LMEAN_FIELD) THEN
       CALL INI_MEAN_FIELD
    END IF
!
  ELSE
    !Necessary to have a 'valid' CNAME when calling some subroutines
    TZBAKFILE => TFILE_DUMMY
  END IF
ELSE
  !Necessary to have a 'valid' CNAME when calling some subroutines
  TZBAKFILE => TFILE_DUMMY
END IF
!
IF (IOUT < NOUT_NUMB ) THEN
  IF (KTCOUNT == TOUTPUTN(IOUT+1)%NSTEP) THEN
    IOUT=IOUT+1
    !
    TZOUTFILE => TOUTPUTN(IOUT)%TFILE
    !
    CALL IO_FILE_OPEN_ll(TZOUTFILE)
    !
    CALL IO_WRITE_HEADER(TZOUTFILE)
    CALL IO_WRITE_FIELDLIST(TOUTPUTN(IOUT))
    CALL IO_WRITE_FIELD_USER(TOUTPUTN(IOUT))
    !
    CALL IO_FILE_CLOSE_ll(TZOUTFILE)
    !
  END IF
END IF
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
  LBU_ENABLE = CBUTYPE /='NONE'.AND. CBUTYPE /='SKIP' 
ELSE
  LBU_ENABLE = .FALSE.
END IF
!
IF (NBUMOD==IMI .AND. CBUTYPE=='MASK' ) THEN
  CALL SET_MASK
  IF (LBU_RU)   XBURHODJU(:,NBUTIME,:) = XBURHODJU(:,NBUTIME,:)    &
                            + MASK_COMPRESS(MXM(XRHODJ))
  IF (LBU_RV)   XBURHODJV(:,NBUTIME,:) = XBURHODJV(:,NBUTIME,:)    &
                            + MASK_COMPRESS(MYM(XRHODJ))
  IF (LBU_RW)   XBURHODJW(:,NBUTIME,:) = XBURHODJW(:,NBUTIME,:)    &
                            + MASK_COMPRESS(MZM(1,IKU,1,XRHODJ))
  IF (ALLOCATED(XBURHODJ))                                         &
                XBURHODJ (:,NBUTIME,:) = XBURHODJ (:,NBUTIME,:)    &
                              + MASK_COMPRESS(XRHODJ)
END IF
!
IF (NBUMOD==IMI .AND. CBUTYPE=='CART' ) THEN
  IF (LBU_RU)   XBURHODJU(:,:,:) = XBURHODJU(:,:,:)    &
                + CART_COMPRESS(MXM(XRHODJ))
  IF (LBU_RV)   XBURHODJV(:,:,:) = XBURHODJV(:,:,:)    &
                + CART_COMPRESS(MYM(XRHODJ))
  IF (LBU_RW)   XBURHODJW(:,:,:) = XBURHODJW(:,:,:)    &
                + CART_COMPRESS(MZM(1,IKU,1,XRHODJ))
  IF (ALLOCATED(XBURHODJ))                             &
                XBURHODJ (:,:,:) = XBURHODJ (:,:,:)    &
                + CART_COMPRESS(XRHODJ)
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
CALL INITIAL_GUESS ( NRR, NSV, KTCOUNT, XRHODJ,IMI, XTSTEP,                 &
                     XRUS, XRVS, XRWS, XRTHS, XRRS, XRTKES, XRSVS,          &
                     XUT, XVT, XWT, XTHT, XRT, XTKET, XSVT )
!
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
  DPTR_XUM=>XUT
  DPTR_XVM=>XVT
  DPTR_XWM=>XWT
  DPTR_XTHM=>XTHT
  DPTR_XRM=>XRT
  DPTR_XTKEM=>XTKET
  DPTR_XSVM=>XSVT
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
  CALL TWO_WAY(     NRR,NSV,KTCOUNT,DPTR_XRHODJ,IMI,XTSTEP,                                    &
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
!-------------------------------------------------------------------------------
!
!*       10.    FORCING
!               -------
!
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XDXHAT,XDYHAT,XZZ,ZJ)
  XMAP=1.
ELSE
  CALL SM_GRIDPROJ(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XLATORI,XLONORI, &
                   XMAP,XLAT,XLON,XDXHAT,XDYHAT,XZZ,ZJ)
END IF
!
IF ( LFORCING ) THEN
  CALL FORCING(XTSTEP,LUSERV,XRHODJ,XCORIOZ,XZHAT,XZZ,TDTCUR,&
               XUFRC_PAST, XVFRC_PAST,                &
               XUT,XVT,XWT,XTHT,XTKET,XRT,XSVT,       &
               XRUS,XRVS,XRWS,XRTHS,XRTKES,XRRS,XRSVS,IMI,ZJ)
END IF
!
IF ( L2D_ADV_FRC ) THEN 
  CALL ADV_FORCING_n(XRHODJ,TDTCUR,XTHT,XRT,XZZ,XRTHS,XRRS)
END IF
IF ( L2D_REL_FRC ) THEN 
  CALL REL_FORCING_n(XRHODJ,TDTCUR,XTHT,XRT,XZZ,XRTHS,XRRS)
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
               XUT,XVT,XWT,XTHT,XRT,            &
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
CALL DYN_SOURCES( NRR,NRRL, NRRI,                              &
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
IF ( LNUMDIFU .OR. LNUMDIFTH .OR. LNUMDIFSV ) THEN
!
  CALL UPDATE_HALO_ll(TZFIELDT_ll, IINFO_ll)
  CALL UPDATE_HALO2_ll(TZFIELDT_ll, TZHALO2T_ll, IINFO_ll)
  IF ( .NOT. LSTEADYLS ) THEN
     CALL UPDATE_HALO_ll(TZLSFIELD_ll, IINFO_ll)
     CALL UPDATE_HALO2_ll(TZLSFIELD_ll, TZLSHALO2_ll, IINFO_ll)
  END IF
  CALL NUM_DIFF ( CLBCX, CLBCY, NRR, NSV,                               &
                  XDK2U, XDK4U, XDK2TH, XDK4TH, XDK2SV, XDK4SV, IMI,    &
                  XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT,                &
                  XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XRHODJ,               &
                  XRUS, XRVS, XRWS, XRTHS, XRTKES, XRRS, XRSVS,         &
                  LZDIFFU,LNUMDIFU, LNUMDIFTH, LNUMDIFSV,               &
                  TZHALO2T_ll, TZLSHALO2_ll,XZDIFFU_HALO2      )
END IF
!
DO JSV = NSV_CHEMBEG,NSV_CHEMEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_CHICBEG,NSV_CHICEND
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
DO JSV = NSV_PPBEG,NSV_PPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
#ifdef MNH_FOREFIRE
DO JSV = NSV_FFBEG,NSV_FFEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
#endif
DO JSV = NSV_CSBEG,NSV_CSEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_DSTDEPBEG,NSV_DSTDEPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_SLTDEPBEG,NSV_SLTDEPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_AERDEPBEG,NSV_AERDEPEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
DO JSV = NSV_SNWBEG,NSV_SNWEND
  XRSVS(:,:,:,JSV) = MAX(XRSVS(:,:,:,JSV),0.)
END DO
IF (CELEC .NE. 'NONE') THEN
  XRSVS(:,:,:,NSV_ELECBEG) = MAX(XRSVS(:,:,:,NSV_ELECBEG),0.)
  XRSVS(:,:,:,NSV_ELECEND) = MAX(XRSVS(:,:,:,NSV_ELECEND),0.)
END IF
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
IF(LVE_RELAX .OR. LVE_RELAX_GRD .OR. LHORELAX_UVWTH .OR. LHORELAX_RV .OR.&
   LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI .OR. LHORELAX_RS .OR.   &
   LHORELAX_RG .OR. LHORELAX_RH .OR. LHORELAX_TKE .OR.                   &
   ANY(LHORELAX_SV)) THEN
  CALL RELAXATION (LVE_RELAX,LVE_RELAX_GRD,LHORELAX_UVWTH,LHORELAX_RV,LHORELAX_RC,   &
                   LHORELAX_RR,LHORELAX_RI,LHORELAX_RS,LHORELAX_RG,    &
                   LHORELAX_RH,LHORELAX_TKE,LHORELAX_SV,               &
                   LHORELAX_SVC2R2,LHORELAX_SVC1R3,                    &
                   LHORELAX_SVELEC,LHORELAX_SVLG,                      &
                   LHORELAX_SVCHEM,LHORELAX_SVCHIC,LHORELAX_SVAER,     &
                   LHORELAX_SVDST,LHORELAX_SVSLT,LHORELAX_SVPP,        &
                   LHORELAX_SVCS,LHORELAX_SVSNW,                       &
#ifdef MNH_FOREFIRE
                   LHORELAX_SVFF,                                      &
#endif
                   KTCOUNT,NRR,NSV,XTSTEP,XRHODJ,                      &
                   XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET,              &
                   XLSUM, XLSVM, XLSWM, XLSTHM,                        &
                   XLBXUM, XLBXVM, XLBXWM, XLBXTHM,                    &
                   XLBXRM, XLBXSVM, XLBXTKEM,                          &
                   XLBYUM, XLBYVM, XLBYWM, XLBYTHM,                    &
                   XLBYRM, XLBYSVM, XLBYTKEM,                          &
                   NALBOT, XALK, XALKW,                                &
                   NALBAS, XALKBAS, XALKWBAS,                          &
                   LMASK_RELAX,XKURELAX, XKVRELAX, XKWRELAX,           &
                   NRIMX,NRIMY,                                        &
                   XRUS, XRVS, XRWS, XRTHS, XRRS, XRSVS, XRTKES        )
END IF

IF (CELEC.NE.'NONE' .AND. LRELAX2FW_ION) THEN
   CALL RELAX2FW_ION (KTCOUNT, IMI, XTSTEP, XRHODJ, XSVT, NALBOT,      &
                      XALK, LMASK_RELAX, XKWRELAX, XRSVS )   
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
CALL PHYS_PARAM_n(KTCOUNT,TZBAKFILE, GCLOSE_OUT,                &
                  XT_RAD,XT_SHADOWS,XT_DCONV,XT_GROUND,XT_MAFL, &
                  XT_DRAG,XT_TURB,XT_TRACER,                    &
                  ZTIME,ZWETDEPAER,GMASKkids,GCLOUD_ONLY)
!
IF (CDCONV/='NONE') THEN
  XPACCONV = XPACCONV + XPRCONV * XTSTEP
  IF (LCH_CONV_LINOX) THEN
    XIC_TOTAL_NUMBER = XIC_TOTAL_NUMBER + XIC_RATE * XTSTEP
    XCG_TOTAL_NUMBER = XCG_TOTAL_NUMBER + XCG_RATE * XTSTEP
  END IF
END IF
!
IF (IBAK>0 .AND. IBAK <= NBAK_NUMB ) THEN
  IF (KTCOUNT == TBACKUPN(IBAK)%NSTEP) THEN
    IF (CSURF=='EXTE') THEN
      CALL GOTO_SURFEX(IMI)
      CALL DIAG_SURF_ATM_n(YSURF_CUR,'MESONH')
      TFILE_SURFEX => TZBAKFILE
      CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR,'MESONH','ALL')
      NULLIFY(TFILE_SURFEX)
    END IF
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
      CALL LS_COUPLING(XTSTEP,GSTEADY_DMASS,CCONF,                          &
             CGETTKET,                                                      &
             CGETRVT,CGETRCT,CGETRRT,CGETRIT,                               &
             CGETRST,CGETRGT,CGETRHT,CGETSVT,LCH_INIT_FIELD, NSV,           &
             NIMAX_ll,NJMAX_ll,                                             &
             NSIZELBX_ll,NSIZELBXU_ll,NSIZELBY_ll,NSIZELBYV_ll,             &
             NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                 &
             NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,         &
             XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XLSZWSM,XDRYMASST,             &
             XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,          &
             XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,          &
             XLSUS,XLSVS,XLSWS,XLSTHS,XLSRVS,XLSZWSS,XDRYMASSS,             &
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
      DO JSV=NSV_DSTDEPBEG,NSV_DSTDEPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTDEPBEG,NSV_SLTDEPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_PPBEG,NSV_PPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
#ifdef MNH_FOREFIRE
      DO JSV=NSV_FFBEG,NSV_FFEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
#endif
      DO JSV=NSV_CSBEG,NSV_CSEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SNWBEG,NSV_SNWEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
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
!
!
!*      8 Bis . Blowing snow scheme
!              ---------
!
IF(LBLOWSNOW) THEN
 CALL BLOWSNOW(CLBCX,CLBCY,XTSTEP,NRR,XPABST,XTHT,XRT,XZZ,XRHODREF,  &
                     XRHODJ,XEXNREF,XRRS,XRTHS,XSVT,XRSVS,XSNWSUBL3D )
ENDIF
!
!-----------------------------------------------------------------------
!
!*       8 Ter  VISCOSITY (no-slip condition inside)
!              ---------
!
!
IF ( LVISC ) THEN
!
ZTIME1 = ZTIME2
!
   CALL VISCOSITY(CLBCX, CLBCY, NRR, NSV, XMU_V,XPRANDTL,         &
                  LVISC_UVW,LVISC_TH,LVISC_SV,LVISC_R,            &
                  LDRAG,    &
                  XUT, XVT, XWT, XTHT, XRT, XSVT,                 &
                  XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,           &
                  XRUS, XRVS, XRWS, XRTHS, XRRS, XRSVS,XDRAG )
!
ENDIF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_VISC = XT_VISC + ZTIME2 - ZTIME1
!!
!-------------------------------------------------------------------------------
!
!*       9.    ADVECTION
!              ---------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
!
!
CALL MPPDB_CHECK3DM("before ADVEC_METSV:XU/V/W/TH/TKE/T,XRHODJ",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET,XRHODJ)
 CALL ADVECTION_METSV ( TZBAKFILE, GCLOSE_OUT,CUVW_ADV_SCHEME,         &
                 CMET_ADV_SCHEME, CSV_ADV_SCHEME, CCLOUD, NSPLIT,      &
                 LSPLIT_CFL, XSPLIT_CFL, LCFL_WRIT,                    &
                 CLBCX, CLBCY, NRR, NSV, TDTCUR, XTSTEP,               &
                 XUT, XVT, XWT, XTHT, XRT, XTKET, XSVT, XPABST,        &
                 XTHVREF, XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,        &
                 XRTHS, XRRS, XRTKES, XRSVS,                           &
                 XRTHS_CLD, XRRS_CLD, XRSVS_CLD, XRTKEMS               )
CALL MPPDB_CHECK3DM("after  ADVEC_METSV:XU/V/W/TH/TKE/T,XRHODJ ",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET,XRHODJ)
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ADV = XT_ADV + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
ZRWS = XRWS
!
CALL GRAVITY_IMPL ( CLBCX, CLBCY, NRR, NRRL, NRRI,XTSTEP,            &
                 XTHT, XRT, XTHVREF, XRHODJ, XRWS, XRTHS, XRRS,      &
                 XRTHS_CLD, XRRS_CLD                                 )   
!
! At the initial instant the difference with the ref state creates a 
! vertical velocity production that must not be advected as it is 
! compensated by the pressure gradient
!
IF (KTCOUNT == 1 .AND. CCONF=='START') XRWS_PRES = - (XRWS - ZRWS) 
!
CALL SECOND_MNH2(ZTIME2)
!
XT_GRAV = XT_GRAV + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
! 
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
!MPPDB_CHECK_LB=.TRUE.
CALL MPPDB_CHECK3DM("before ADVEC_UVW:XU/V/W/TH/TKE/T,XRHODJ,XRU/V/Ws",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET,XRHODJ,XRUS,XRVS,XRWS)
IF ((CUVW_ADV_SCHEME(1:3)=='CEN') .AND. (CTEMP_SCHEME == 'LEFR')) THEN
  IF (CUVW_ADV_SCHEME=='CEN4TH') THEN
    NULLIFY(TZFIELDC_ll)
    NULLIFY(TZHALO2C_ll)
      CALL ADD3DFIELD_ll(TZFIELDC_ll, XUT)
      CALL ADD3DFIELD_ll(TZFIELDC_ll, XVT)
      CALL ADD3DFIELD_ll(TZFIELDC_ll, XWT)
      CALL INIT_HALO2_ll(TZHALO2C_ll,3,IIU,IJU,IKU)
      CALL UPDATE_HALO_ll(TZFIELDC_ll,IINFO_ll)
      CALL UPDATE_HALO2_ll(TZFIELDC_ll, TZHALO2C_ll, IINFO_ll)
  END IF
 CALL ADVECTION_UVW_CEN(CUVW_ADV_SCHEME,                &
                           CLBCX, CLBCY,                           &
                           XTSTEP, KTCOUNT,                        &
                           XUM, XVM, XWM, XDUM, XDVM, XDWM,        &
                           XUT, XVT, XWT,                          &
                           XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,   &
                           XRUS,XRVS, XRWS,                        &
                           TZHALO2C_ll                             )
  IF (CUVW_ADV_SCHEME=='CEN4TH') THEN
    CALL CLEANLIST_ll(TZFIELDC_ll)
    NULLIFY(TZFIELDC_ll)
    CALL  DEL_HALO2_ll(TZHALO2C_ll)
    NULLIFY(TZHALO2C_ll)
  END IF
ELSE

  CALL ADVECTION_UVW(CUVW_ADV_SCHEME, CTEMP_SCHEME,                  &
                 NWENO_ORDER, LSPLIT_WENO,                           &
                 CLBCX, CLBCY, XTSTEP,                               &
                 XUT, XVT, XWT,                                      &
                 XRHODJ, XDXX, XDYY, XDZZ, XDZX, XDZY,               &
                 XRUS, XRVS, XRWS,                                   &
                 XRUS_PRES, XRVS_PRES, XRWS_PRES                     )
END IF
!
CALL MPPDB_CHECK3DM("after  ADVEC_UVW:XU/V/W/TH/TKE/T,XRHODJ,XRU/V/Ws",PRECISION,&
                   &  XUT, XVT, XWT, XTHT, XTKET,XRHODJ,XRUS,XRVS,XRWS)
!MPPDB_CHECK_LB=.FALSE.
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ADVUVW = XT_ADVUVW + ZTIME2 - ZTIME1 - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
IF (NMODEL_CLOUD==IMI .AND. CTURBLEN_CLOUD/='NONE') THEN
  CALL TURB_CLOUD_INDEX(XTSTEP,TZBAKFILE,                         &
                        LTURB_DIAG,GCLOSE_OUT,NRRI,               &
                        XRRS,XRT,XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY, &
                        XCEI )
END IF
!
!-------------------------------------------------------------------------------
!
!*       18.    LATERAL BOUNDARY CONDITION FOR THE NORMAL VELOCITY
!               --------------------------------------------------
!
ZTIME1 = ZTIME2
!
CALL MPPDB_CHECK3DM("before RAD_BOUND :XRU/V/WS",PRECISION,XRUS,XRVS,XRWS)
ZRUS=XRUS
ZRVS=XRVS
ZRWS=XRWS
!
  CALL RAD_BOUND (CLBCX,CLBCY,CTURB,XCARPKMAX,           &
                XTSTEP,                                  &
                XDXHAT, XDYHAT, XZHAT,                   &
                XUT, XVT,                                &
                XLBXUM, XLBYVM, XLBXUS, XLBYVS,          &
                XCPHASE, XCPHASE_PBL, XRHODJ,            &
                XTKET,XRUS, XRVS, XRWS                   )
ZRUS=XRUS-ZRUS
ZRVS=XRVS-ZRVS
ZRWS=XRWS-ZRWS
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
ZPABST = XPABST
!
IF(.NOT. L1D) THEN
!
CALL MPPDB_CHECK3DM("before pressurez:XRU/V/WS",PRECISION,XRUS,XRVS,XRWS)
  XRUS_PRES = XRUS
  XRVS_PRES = XRVS
  XRWS_PRES = XRWS
!
  CALL PRESSUREZ( CLBCX,CLBCY,CPRESOPT,NITR,LITRADJ,KTCOUNT, XRELAX,IMI, &
                  XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY,XDXHATM,XDYHATM,XRHOM, &
                  XAF,XBFY,XCF,XTRIGSX,XTRIGSY,NIFAXX,NIFAXY,            &
                  NRR,NRRL,NRRI,XDRYMASST,XREFMASS,XMASS_O_PHI0,         &
                  XTHT,XRT,XRHODREF,XTHVREF,XRVREF,XEXNREF, XLINMASS,    &
                  XRUS, XRVS, XRWS, XPABST,                              &
                  XBFB,&
                  XBF_SXP2_YP1_Z) !JUAN Z_SPLITING
!
  XRUS_PRES = XRUS - XRUS_PRES + ZRUS
  XRVS_PRES = XRVS - XRVS_PRES + ZRVS
  XRWS_PRES = XRWS - XRWS_PRES + ZRWS
!
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_PRESS = XT_PRESS + ZTIME2 - ZTIME1 &
           - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
!-------------------------------------------------------------------------------
!
!*       20.    CHEMISTRY/AEROSOLS
!               ------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (LUSECHEM) THEN
  CALL CH_MONITOR_n(ZWETDEPAER,KTCOUNT,XTSTEP, ILUOUT, NVERB)
END IF
!
! For inert aerosol (dust and sea salt) => aer_monitor_n
IF ((LDUST).OR.(LSALT)) THEN
!
! tests to see if any cloud exists
!   
    GCLD=.TRUE.
    IF (GCLD .AND. NRR.LE.3 ) THEN 
      IF( MAXVAL(XCLDFR(:,:,:)).LE. 1.E-10 .AND. GCLOUD_ONLY ) THEN
          GCLD = .FALSE.                ! only the cloudy verticals would be 
                                        ! refreshed but there is no clouds 
      END IF
    END IF
!
    IF (GCLD .AND. NRR.GE.4 ) THEN 
      IF( CCLOUD(1:3)=='ICE' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN(4) .AND. GCLOUD_ONLY ) THEN
            GCLD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
      IF( CCLOUD=='C3R5' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN_C1R3(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN_C1R3(4) .AND. GCLOUD_ONLY ) THEN
            GCLD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
      IF( CCLOUD=='LIMA' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN_LIMA(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN_LIMA(4) .AND. GCLOUD_ONLY ) THEN
            GCLD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
    END IF

!
        CALL AER_MONITOR_n(KTCOUNT,XTSTEP, ILUOUT, NVERB, GCLD)
END IF
!
!
CALL SECOND_MNH2(ZTIME2)
!
XT_CHEM = XT_CHEM + ZTIME2 - ZTIME1 &
      - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
ZTIME = ZTIME + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS

!-------------------------------------------------------------------------------
!
!*       20.    WATER MICROPHYSICS
!               ------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (CCLOUD /= 'NONE' .AND. CELEC == 'NONE') THEN
!
  IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO' .OR. CCLOUD == 'C3R5' &
                                             .OR. CCLOUD == "LIMA" ) THEN
    IF ( LFORCING ) THEN
      ZWT_ACT_NUC(:,:,:) = XWT(:,:,:) + XWTFRC(:,:,:)
    ELSE
      ZWT_ACT_NUC(:,:,:) = XWT(:,:,:)
    END IF
    IF (CTURB /= 'NONE' ) THEN
     IF ( ((CCLOUD=='C2R2'.OR.CCLOUD=='KHKO').AND.LACTTKE) .OR. (CCLOUD=='LIMA'.AND.MACTTKE) ) THEN 
       ZWT_ACT_NUC(:,:,:) = ZWT_ACT_NUC(:,:,:) +  (2./3. * XTKET(:,:,:))**0.5
     ELSE
       ZWT_ACT_NUC(:,:,:) = ZWT_ACT_NUC(:,:,:) 
     ENDIF
    ENDIF
  ELSE
    ZWT_ACT_NUC(:,:,:) = 0.
  END IF
!
  XRTHS_CLD  = XRTHS
  XRRS_CLD   = XRRS
  XRSVS_CLD  = XRSVS
  IF (CSURF=='EXTE') THEN
    ALLOCATE (ZSEA(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ALLOCATE (ZTOWN(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ZSEA(:,:) = 0.
    ZTOWN(:,:)= 0.
    CALL MNHGET_SURF_PARAM_n (PSEA=ZSEA(:,:),PTOWN=ZTOWN(:,:))
    CALL RESOLVED_CLOUD ( CCLOUD, CACTCCN, CSCONV, CMF_CLOUD, NRR, NSPLITR,    &
                          NSPLITG, IMI, KTCOUNT,                               &
                          CLBCX,CLBCY,TZBAKFILE, CRAD, CTURBDIM,               &
                          GCLOSE_OUT, LSUBG_COND,LSIGMAS,CSUBG_AUCV,XTSTEP,    &
                          XZZ, XRHODJ, XRHODREF, XEXNREF,                      &
                          ZPABST, XTHT,XRT,XSIGS,VSIGQSAT,XMFCONV,XTHM,XRCM,   &
                          XPABSM, ZWT_ACT_NUC,XDTHRAD, XRTHS, XRRS,            &
                          XSVT, XRSVS,                                         &
                          XSRCT, XCLDFR,XCIT,                                  &
                          LSEDIC,KACTIT, KSEDC, KSEDI, KRAIN, KWARM, KHHONI,   &
                          LCONVHG, XCF_MF,XRC_MF, XRI_MF,                      &
                          XINPRC,ZINPRC3D,XINPRR, XINPRR3D, XEVAP3D,             &
                          XINPRS,ZINPRS3D, XINPRG,ZINPRG3D, XINPRH,ZINPRH3D,   &
                          XSOLORG, XMI,ZSPEEDC, ZSPEEDR, ZSPEEDS, ZSPEEDG, ZSPEEDH, &
                          XINDEP, XSUPSAT, XNACT, XNPRO,XSSPRO, XRAINFR,       &
                          ZSEA, ZTOWN                                          )
    DEALLOCATE(ZTOWN)
  ELSE
    CALL RESOLVED_CLOUD ( CCLOUD, CACTCCN, CSCONV, CMF_CLOUD, NRR, NSPLITR,    &
                          NSPLITG, IMI, KTCOUNT,                               &
                          CLBCX,CLBCY,TZBAKFILE, CRAD, CTURBDIM,               &
                          GCLOSE_OUT, LSUBG_COND,LSIGMAS,CSUBG_AUCV,           &
                          XTSTEP,XZZ, XRHODJ, XRHODREF, XEXNREF,               &
                          ZPABST, XTHT,XRT,XSIGS,VSIGQSAT,XMFCONV,XTHM,XRCM,   &
                          XPABSM, ZWT_ACT_NUC,XDTHRAD, XRTHS, XRRS,            &
                          XSVT, XRSVS,                                         &
                          XSRCT, XCLDFR,XCIT,                                  &
                          LSEDIC,KACTIT, KSEDC, KSEDI, KRAIN, KWARM, KHHONI,   &
                          LCONVHG, XCF_MF,XRC_MF, XRI_MF,                      &
                          XINPRC,ZINPRC3D,XINPRR, XINPRR3D, XEVAP3D,             &
                          XINPRS,ZINPRS3D, XINPRG,ZINPRG3D, XINPRH,ZINPRH3D,   &
                          XSOLORG, XMI,ZSPEEDC, ZSPEEDR, ZSPEEDS, ZSPEEDG, ZSPEEDH, &
                          XINDEP, XSUPSAT, XNACT, XNPRO,XSSPRO, XRAINFR        )
  END IF
  XRTHS_CLD  = XRTHS - XRTHS_CLD
  XRRS_CLD   = XRRS  - XRRS_CLD
  XRSVS_CLD  = XRSVS - XRSVS_CLD
!
  IF (CCLOUD /= 'REVE' ) THEN
    XACPRR = XACPRR + XINPRR * XTSTEP
          IF ( (CCLOUD(1:3) == 'ICE' .AND. LSEDIC ) .OR.                     &
        ((CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO' &
                           .OR. CCLOUD == 'LIMA' ) .AND. KSEDC ) )    THEN
      XACPRC = XACPRC + XINPRC * XTSTEP
      IF (LDEPOSC .OR. LDEPOC) XACDEP = XACDEP + XINDEP * XTSTEP
    END IF
    IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. &
                                 (CCLOUD == 'LIMA' .AND. LCOLD ) ) THEN
      XACPRS = XACPRS + XINPRS * XTSTEP
      XACPRG = XACPRG + XINPRG * XTSTEP
      IF (CCLOUD == 'ICE4' .OR. (CCLOUD == 'LIMA' .AND. LHAIL)) XACPRH = XACPRH + XINPRH * XTSTEP          
    END IF
!
! Lessivage des CCN et IFN nucléables par Slinn
!
    IF (LSCAV .AND. (CCLOUD == 'LIMA')) THEN
      CALL LIMA_PRECIP_SCAVENGING(CCLOUD, ILUOUT, KTCOUNT,XTSTEP,XRT(:,:,:,3), &
                              XRHODREF, XRHODJ, XZZ, XPABST, XTHT,             &
                              XSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END),           &
                              XRSVS(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END), XINPAP   )
!
      XACPAP(:,:) = XACPAP(:,:) + XINPAP(:,:) * XTSTEP
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
!*       21.    CLOUD ELECTRIFICATION AND LIGHTNING FLASHES
!               -------------------------------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (CELEC /= 'NONE' .AND. (CCLOUD(1:3) == 'ICE')) THEN
  ZWT_ACT_NUC(:,:,:) = 0.
!
  XRTHS_CLD = XRTHS
  XRRS_CLD  = XRRS
  XRSVS_CLD = XRSVS
  IF (CSURF=='EXTE') THEN
    ALLOCATE (ZSEA(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ALLOCATE (ZTOWN(SIZE(XRHODJ,1),SIZE(XRHODJ,2)))
    ZSEA(:,:) = 0.
    ZTOWN(:,:)= 0.
    CALL MNHGET_SURF_PARAM_n (PSEA=ZSEA(:,:),PTOWN=ZTOWN(:,:))
    CALL RESOLVED_ELEC_n (CCLOUD, CSCONV, CMF_CLOUD,                     &
                          NRR, NSPLITR, IMI, KTCOUNT, OEXIT,             &
                          CLBCX, CLBCY, CRAD, CTURBDIM,                  &
                          LSUBG_COND, LSIGMAS,VSIGQSAT,CSUBG_AUCV,       &
                          XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF,        &
                          ZPABST, XTHT, XRTHS, XWT,  XRT, XRRS,          &
                          XSVT, XRSVS, XCIT,                             &
                          XSIGS, XSRCT, XCLDFR, XMFCONV, XCF_MF, XRC_MF, &
                          XRI_MF, LSEDIC, LWARM,                         &
                          XINPRC, XINPRR, XINPRR3D, XEVAP3D,             &
                          XINPRS, XINPRG, XINPRH,                        &
                          ZSEA, ZTOWN                                    )
    DEALLOCATE(ZTOWN)
  ELSE
    CALL RESOLVED_ELEC_n (CCLOUD, CSCONV, CMF_CLOUD,                     &
                          NRR, NSPLITR, IMI, KTCOUNT, OEXIT,             &
                          CLBCX, CLBCY, CRAD, CTURBDIM,                  &
                          LSUBG_COND, LSIGMAS,VSIGQSAT, CSUBG_AUCV,      &
                          XTSTEP, XZZ, XRHODJ, XRHODREF, XEXNREF,        &
                          ZPABST, XTHT, XRTHS, XWT,                      &
                          XRT, XRRS, XSVT, XRSVS, XCIT,                  &
                          XSIGS, XSRCT, XCLDFR, XMFCONV, XCF_MF, XRC_MF, &
                          XRI_MF, LSEDIC, LWARM,                         &
                          XINPRC, XINPRR, XINPRR3D, XEVAP3D,             &
                          XINPRS, XINPRG, XINPRH                         )
  END IF
  XRTHS_CLD = XRTHS - XRTHS_CLD
  XRRS_CLD  = XRRS  - XRRS_CLD
  XRSVS_CLD = XRSVS - XRSVS_CLD
!
  XACPRR = XACPRR + XINPRR * XTSTEP
  IF ((CCLOUD(1:3) == 'ICE' .AND. LSEDIC)) & 
       XACPRC = XACPRC + XINPRC * XTSTEP
  IF (CCLOUD(1:3) == 'ICE') THEN
    XACPRS = XACPRS + XINPRS * XTSTEP
    XACPRG = XACPRG + XINPRG * XTSTEP
    IF (CCLOUD == 'ICE4') XACPRH = XACPRH + XINPRH * XTSTEP          
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
XT_ELEC = XT_ELEC + ZTIME2 - ZTIME1 &
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
!*       21. bis    MEAN_UM
!               --------------------
!
IF (LMEAN_FIELD) THEN
   CALL MEAN_FIELD(XUT, XVT, XWT, XTHT, XTKET, XPABST)
END IF
!
!-------------------------------------------------------------------------------
!
!*       22.    UPDATE HALO OF EACH SUBDOMAINS FOR TIME T+DT
!               --------------------------------------------
!
ZTIME1 = ZTIME2
!
CALL EXCHANGE (XTSTEP,NRR,NSV,XRHODJ,TZFIELDS_ll,     &
               XRUS, XRVS,XRWS,XRTHS,XRRS,XRTKES,XRSVS)
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
CALL ENDSTEP  ( XTSTEP,NRR,NSV,KTCOUNT,IMI,               &
                CUVW_ADV_SCHEME,CTEMP_SCHEME,XRHODJ,      &
                XRUS,XRVS,XRWS,XDRYMASSS,                 &
                XRTHS,XRRS,XRTKES,XRSVS,                  &
                XLSUS,XLSVS,XLSWS,                        &
                XLSTHS,XLSRVS,XLSZWSS,                    &
                XLBXUS,XLBXVS,XLBXWS,                     &
                XLBXTHS,XLBXRS,XLBXTKES,XLBXSVS,          &
                XLBYUS,XLBYVS,XLBYWS,                     &
                XLBYTHS,XLBYRS,XLBYTKES,XLBYSVS,          &
                XUM,XVM,XWM,XZWS,                         &
                XUT,XVT,XWT,XPABST,XDRYMASST,             &
                XTHT, XRT, XTHM, XRCM, XPABSM,XTKET, XSVT,&
                XLSUM,XLSVM,XLSWM,                        &
                XLSTHM,XLSRVM,XLSZWSM,                    &
                XLBXUM,XLBXVM,XLBXWM,                     &
                XLBXTHM,XLBXRM,XLBXTKEM,XLBXSVM,          &
                XLBYUM,XLBYVM,XLBYWM,                     &
                XLBYTHM,XLBYRM,XLBYTKEM,XLBYSVM           )
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
  CALL AIRCRAFT_BALLOON(XTSTEP,                                               &
                      TDTEXP, TDTMOD, TDTSEG, TDTCUR,                         &
                      XXHAT, XYHAT, XZZ, XMAP, XLONORI, XLATORI,              &
                      XUT, XVT, XWT, XPABST, XTHT, XRT, XSVT, XTKET, XTSRAD,  &
                      XRHODREF,XCIT,PSEA=ZSEA(:,:))


!-------------------------------------------------------------------------------
!
!*       24.2    STATION (observation diagnostic)
!               --------------------------------
!
IF (LSTATION)                                                            &
  CALL STATION_n(XTSTEP,                                                 &
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
  CALL PROFILER_n(XTSTEP,                                                &
                  TDTEXP, TDTMOD, TDTSEG, TDTCUR,                        &
                  XXHAT, XYHAT, XZZ,XRHODREF,                            &
                  XUT, XVT, XWT, XTHT, XRT, XSVT, XTKET, XTSRAD, XPABST, &
                  XAER, XCLDFR, XCIT)
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
IF ( .NOT. LIO_NO_WRITE ) THEN
  IF (NBUMOD==IMI .AND. CBUTYPE/='NONE') THEN
    CALL ENDSTEP_BUDGET(TDIAFILE,KTCOUNT,TDTCUR,TDTMOD,XTSTEP,NSV)
  END IF
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
  CALL IO_FILE_CLOSE_ll(TZBAKFILE)
END IF
!
!-------------------------------------------------------------------------------
!
!*       27.    CURRENT TIME REFRESH
!               --------------------
!
TDTCUR%TIME=TDTCUR%TIME + XTSTEP
CALL DATETIME_CORRECTDATE(TDTCUR)
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
  IF ( .NOT. LIO_NO_WRITE ) THEN
    IF (LSERIES) CALL WRITE_SERIES_n(TDIAFILE)
    CALL WRITE_AIRCRAFT_BALLOON(TDIAFILE)
    CALL WRITE_STATION_n(TDIAFILE)
    CALL WRITE_PROFILER_n(TDIAFILE)
    CALL WRITE_LES_n(TDIAFILE,' ')
    CALL WRITE_LES_n(TDIAFILE,'A')
    CALL WRITE_LES_n(TDIAFILE,'E')
    CALL WRITE_LES_n(TDIAFILE,'H')
    CALL MENU_DIACHRO(TDIAFILE,'END')
    CALL IO_FILE_CLOSE_ll(TDIAFILE)
  END IF
  !
  CALL IO_FILE_CLOSE_ll(TINIFILE)
  IF (CSURF=="EXTE") CALL IO_FILE_CLOSE_ll(TINIFILEPGD,OPARALLELIO=.FALSE.)
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
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_SEND,ZTOT,    '   W3D_SEND ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_RECV,ZTOT,    '   W3D_RECV ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_WRIT,ZTOT,    '   W3D_WRIT ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_WAIT,ZTOT,    '   W3D_WAIT ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT3D_ALL ,ZTOT,    '   W3D_ALL ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT2D_GATH,ZTOT,    '   W2D_GATH ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT2D_WRIT,ZTOT,    '   W2D_WRIT ','-')
    CALL TIME_STAT_ll(TIMEZ%T_WRIT2D_ALL ,ZTOT,    '   W2D_ALL ','-')
  CALL TIME_STAT_ll(XT_GUESS,ZTOT,      ' INITIAL_GUESS','=')
  CALL TIME_STAT_ll(XT_2WAY,ZTOT,       ' TWO WAY','=')
  CALL TIME_STAT_ll(XT_ADV,ZTOT,        ' ADVECTION MET','=')
  CALL TIME_STAT_ll(XT_ADVUVW,ZTOT,     ' ADVECTION UVW','=')
  CALL TIME_STAT_ll(XT_GRAV,ZTOT,       ' GRAVITY','=')
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
    CALL TIME_STAT_ll(XT_SHADOWS,ZTOT,  '   SHADOWS'             ,'-')
    CALL TIME_STAT_ll(XT_DCONV,ZTOT,    '   DEEP CONV = '//CDCONV,'-')
    CALL TIME_STAT_ll(XT_GROUND,ZTOT,   '   GROUND'              ,'-')
    CALL TIME_STAT_ll(XT_TURB,ZTOT,     '   TURB      = '//CTURB ,'-')
    CALL TIME_STAT_ll(XT_MAFL,ZTOT,     '   MAFL      = '//CSCONV,'-')
    CALL TIME_STAT_ll(XT_CHEM,ZTOT,     '   CHIMIE'              ,'-')
  CALL  TIMING_LEGEND()
  CALL TIME_STAT_ll(XT_COUPL,ZTOT,      ' SET_COUPLING','=')
  CALL TIME_STAT_ll(XT_RAD_BOUND,ZTOT,  ' RAD_BOUND','=')
  !
  CALL  TIMING_LEGEND()
  ! 
  CALL TIME_STAT_ll(XT_PRESS,ZTOT,      ' PRESSURE ','=','F')
  !JUAN Z_SPLITTING
    CALL TIME_STAT_ll(TIMEZ%T_MAP_B_SX_YP2_ZP1,ZTOT,          '   REMAP       B=>FFTXZ'  ,'-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_SX_YP2_ZP1_SXP2_Y_ZP1,ZTOT, '   REMAP   FFTXZ=>FFTYZ'  ,'-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_SXP2_Y_ZP1_B,ZTOT,          '   REMAP   FTTYZ=>B'      ,'-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_SXP2_Y_ZP1_SXP2_YP1_Z,ZTOT, '   REMAP   FFTYZ=>SUBZ'   ,'-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_B_SXP2_Y_ZP1,ZTOT,          '   REMAP       B=>FFTYZ-1','-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_SXP2_YP1_Z_SXP2_Y_ZP1,ZTOT, '   REMAP    SUBZ=>FFTYZ-1','-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_SXP2_Y_ZP1_SX_YP2_ZP1,ZTOT, '   REMAP FFTYZ-1=>FFTXZ-1','-','F')
    CALL TIME_STAT_ll(TIMEZ%T_MAP_SX_YP2_ZP1_B,ZTOT,          '   REMAP FFTXZ-1=>B     ' ,'-','F')
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
           XT_ADVUVW  + XT_GRAV +                                         &
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
  CALL IO_FILE_CLOSE_ll(TLUOUT)
  IF (IMI==NMODEL) CALL IO_FILE_CLOSE_ll(TLUOUT0)
END IF
!
END SUBROUTINE MODEL_n
