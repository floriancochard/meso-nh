!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! BUG2 masdev4_7 2007/10/04 11:51:06
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_DEFAULT_DESFM_n
!     ###########################
!
INTERFACE
!
SUBROUTINE DEFAULT_DESFM_n(KMI)
INTEGER,         INTENT(IN)  :: KMI       ! Model index
END SUBROUTINE DEFAULT_DESFM_n
!
END INTERFACE
!
END MODULE MODI_DEFAULT_DESFM_n
!
!
!
!     ###############################
      SUBROUTINE DEFAULT_DESFM_n(KMI)
!     ###############################
!
!!****  *DEFAULT_DESFM_n * - set default values for descriptive variables of 
!!                         model KMI
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set default values for the variables
!     in descriptor files by filling the corresponding variables which
!     are stored in modules.
!       
!
!!**  METHOD
!!    ------
!!      Each variable in modules, which can be initialized by  reading its 
!!    value in the descriptor file is set to a default value.    
!!     When this routine is used during INIT, the modules of the first model 
!!   are used to temporarily store  the variables associated with a nested 
!!   model.  
!!     When this routine is used during  SPAWNING, the modules of a second 
!!   model must be initialized. 
!!     Default values for variables common to all models are set only
!!   at the first call of DEFAULT_DESFM_n (i.e. when KMI=1)
!!
!!
!!    EXTERNAL
!!    --------
!!     NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS : JPHEXT,JPVEXT 
!!
!!      Module MODD_CONF       : CCONF,L2D,L1D,LFLAT,NMODEL,NVERB 
!!
!!      Module MODD_DYN        : XSEGLEN,XASSELIN,LCORIO,LNUMDIFF
!!                               XALKTOP,XALZBOT
!!
!!      Module MODD_FMOUT      : XOUT
!!
!!      Module MODD_NESTING    : NDAD(m),NDTRATIO(m),XWAY(m) 
!!
!!      Module MODD_CONF_n    : LUSERV,LUSERC,LUSERR,LUSERI,LUSERS  
!!                              LUSERG,LUSERH,CSEG,CEXP 
!!
!!      Module MODD_LUNIT_n   : CINIFILE,CCPLFILE
!!
!!
!!      Module MODD_DYN_n     : XTSTEP,CPRESOPT,NITR,XRELAX,LHO_RELAX 
!!         LVE_RELAX,XRIMKMAX,NRIMX,NRIMY
!!
!!      Module MODD_ADV_n : CUVW_ADV_SCHEME,CMET_ADV_SCHEME,CSV_ADV_SCHEME,NLITER
!!
!!      Module MODD_PARAM_n : CTURB,CRAD,CDCONV
!!
!!      Module MODD_LBC_n : CLBCX, CLBCY,NLBLX,NLBLY,XCPHASE 
!!
!!      Module MODD_TURB_n : XIMPL,CTURBLEN,CTURBDIM,LTURB_FLX,LTURB_DIAG,LSUBG_COND
!!                           LTGT_FLX
!!
!!
!!      Module MODD_PARAM_RAD_n:
!!          XDTRAD,XDTRAD_CLONLY,LCLEAR_SKY,NRAD_COLNBR, NRAD_DIAG
!!
!!      Module MODD_BUDGET : CBUTYPE,NBUMOD,XBULEN,NBUKL, NBUKH,LBU_KCP,XBUWRI
!!         NBUIL, NBUIH,NBUJL, NBUJH,LBU_ICP,LBU_JCP,NBUMASK 
!!         LBU_RU,LBU_RV,LBU_RW,LBU_RTH,LBU_RTKE,LBU_RRV,LBU_RRC,LBU_RRR  
!!         LBU_RRI,LBU_RRS,LBU_RRG,LBU_RRH,LBU_RSVx
!!         NADVXU, NADVYU, NADVZU, NCURVU, NCORU, NDIFU, NRELU, NHTURBU,  
!!         NVTURBU, NPRESU
!!         NADVXV, NADVYV, NADVZV, NCURVV, NCORV, NDIFV, NRELV, NHTURBV,  
!!         NVTURBV, NPRESV      
!!         NADVXW, NADVYW, NADVZW, NCURVW, NCORW, NGRAVW, NDIFW, NRELW,  
!!         NHTURBW, NVTURBW, NPRESW
!!         NADVXTH, NADVYTH, NADVZTH, NPREFTH, NDIFTH, NRELTH, NHTURBTH, 
!!         NVTURBTH, NREVATH, NCONDTH        
!!         NADVXTKE, NADVYTKE, NADVZTKE, NDIFTKE, NDPTKE, NTPTKE, NDISSTKE,
!!         NTRTKE
!!         NADVXRV, NADVYRV, NADVZRV, NDIFRV, NRELRV, NHTURBRV, NVTURBRV,
!!         NREVARV, NCONDRV
!!         NADVXRC, NADVYRC, NADVZRC, NDIFRC, NHTURBRC, NVTURBRC, NACCRRC,
!!         NAUTORC, NCONDRC
!!         NADVXRR, NADVYRR, NADVZRR, NDIFRR, NACCRRR, NAUTORR, NREVARR, 
!!         NSEDIRR
!!         NADVXRI, NADVYRI, NADVZRI, NDIFRI
!!         NADVXRS, NADVYRS, NADVZRS, NDIFRS
!!         NADVXRG, NADVYRG, NADVZRG, NDIFRG
!!         NADVXRH, NADVYRH, NADVZRH, NDIFRH
!!         NADVXSVx, NADVYSVx,  NADVZSVx, NDIFSVx, NHTURBSVx, NVTURBSVx 
!!
!!      Module MODD_BLANK :
!!
!!          XDUMMYi, NDUMMYi, LDUMMYi, CDUMMYi 
!!
!!      Module MODD_FRC :
!!
!!          LGEOST_UV_FRC,LGEOST_TH_FRC,LTEND_THRV_FRC
!!          LVERT_MOTION_FRC,LRELAX_THRV_FRC,LRELAX_UV_FRC,XRELAX_TIME_FRC
!!          XRELAX_HEIGHT_FRC,CRELAX_HEIGHT_TYPE,LTRANS,XUTRANS,XVTRANS,
!!          LPGROUND_FRC
!!
!!      Module MODD_PARAM_ICE :
!!
!!          LWARM,CPRISTINE_ICE
!!
!!      Module MODD_PARAM_CONVECT_n :
!!
!!          XDTCONV,LDEEP,LSHAL,LREFRESH_ALL,LDOWN,NICE,LCHTRANS  
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine DEFAULT_DESFM_n)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      02/06/94
!!      Modifications 17/10/94  (Stein)  For LCORIO
!!      Modifications 06/12/94  (Stein)  remove LBOUSS+add LABSLAYER, LNUMDIFF
!!                                       ,LSTEADYLS
!!      Modifications 06/12/94  (Stein)  remove LABSLAYER, add LHO_RELAX, 
!!                                       LVE_RELAX, NRIMX, NRIMY, XRIMKMAX
!!      Modifications 09/01/95  (Lafore) add LSTEADY_DMASS 
!!      Modifications 09/01/95  (Stein)  add the turbulence scheme namelist
!!      Modifications 09/01/95  (Stein)  add the 1D switch
!!      Modifications 10/03/95  (Mallet) add the coupling files
!!                    29/06/95  ( Stein, Nicolau, Hereil) add the budgets 
!!      Modifications 25/09/95  ( Stein )add the LES tools
!!      Modifications 25/10/95  ( Stein )add the radiations
!!      Modifications 23/10/95  (Vila, lafore) new scalar advection scheme
!!      Modifications 24/02/96  (Stein)  change the default value for CCPLFILE
!!      Modifications 12/02/96  (Lafore) transformation to DEFAULT_DESFM_n for 
!!                                       spawning
!!      Modifications 25/04/96  (Suhre)  add the blank module
!!      Modifications 29/07/96  (Pinty&Suhre) add module MODD_FRC
!!      Modifications 11/04/96  (Pinty)  add the rain-ice scheme and modify
!!                                       the splitted arrays in MODD_PARAM_RAD_n
!!      Modifications 11/01/97  (Pinty)  add the deep convection scheme
!!      Modifications 24/11/96  (Masson)  add LREFRESH_ALL in deep convection
!!      Modifications 12/02/96  (Lafore) transformation to DEFAULT_DESFM_n for spawning
!!      Modifications 22/07/96  (Lafore) gridnesting implementation
!!      Modifications 29/07/96  (Lafore) add the module MODD_FMOUT
!!      Modifications 23/06/97  (Stein)  add the equation system name
!!      Modifications 10/07/97  (Masson) add MODD_PARAM_GROUNDn : CROUGH
!!      Modifications 28/07/97  (Masson) remove LREFRESH_ALL and LSTEADY_DMASS
!!      Modifications 08/10/97  (Stein)  switch (_n=1) to initialize the
!!                                       parameters common to all models
!!      Modifications 24/01/98 (Bechtold) add LREFRESH_ALL, LCHTRANS, 
!!                                         LTEND_THRV_FR and LSST_FRC
!!      Modifications 18/07/99  (Stein)  add LRAD_DIAG
!!      Modification  15/03/99 (Masson)  use of XUNDEF
!!      Modification  11/12/00 (Tomasini) Add CSEA_FLUX to MODD_PARAMn
!!      Modification  22/01/01 (Gazen) delete NSV and add LHORELAX_SVC2R2
!!                                     LHORELAX_SVCHEM,LHORELAX_SVLG
!!      Modification 15/03/02 (Solmon) radiation scheme: remove NSPOT and add
!!                                   default for aerosol and cloud rad. prop. control
!!      Modification 22/05/02 (Jabouille) put chimical default here
!!      Modification 01/2004  (Masson) removes surface (externalization)
!!                      09/04 (M. Tomasini) New namelist to modify the 
!!                                             Cloud mixing length 
!!                   07/05 (P.Tulet) New namelists for dust and aerosol
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_PARAMETERS
USE MODD_CONF             !        For INIT only DEFAULT_DESFM1
USE MODD_DYN
USE MODD_NESTING
USE MODD_FMOUT
USE MODD_SERIES
USE MODD_CONF_n           ! modules used to set the default values is only      
USE MODD_LUNIT_n          ! the one corresponding to model 1. These memory 
USE MODD_DIM_n            ! addresses will then be filled by the values read in 
USE MODD_DYN_n            ! the DESFM corresponding to model n which may have    
USE MODD_ADV_n            ! missing values. This is why we affect default values.
USE MODD_PARAM_n          !      For SPAWNING DEFAULT_DESFM2 is also used 
USE MODD_LBC_n           
USE MODD_OUT_n
USE MODD_TURB_n
USE MODD_BUDGET
USE MODD_LES
USE MODD_PARAM_RAD_n
USE MODD_BLANK
USE MODD_FRC
USE MODD_PARAM_ICE
USE MODD_PARAM_C2R2
USE MODD_TURB_CLOUD
USE MODD_PARAM_CONVECT_n
USE MODD_CH_MNHC_n
USE MODD_SERIES_n
USE MODD_NUDGING_n
USE MODD_CH_AEROSOL
USE MODD_DUST
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,         INTENT(IN)  :: KMI       ! Model index 
!
!*       0.2   declaration of local variables
!
INTEGER             :: JM      ! loop index
!
!-------------------------------------------------------------------------------
!
!*      1.    SET DEFAULT VALUES FOR MODD_LUNIT_n :
!             ----------------------------------
!
CINIFILE='INIFILE'
CCPLFILE(:)='    '
!
!-------------------------------------------------------------------------------
!
!*      2.    SET DEFAULT VALUES FOR MODD_CONF AND MODD_CONF_n :
!             ------------------------------------------------
!
IF (KMI == 1) THEN                                             
  CCONF      ='START'
  LTHINSHELL = .FALSE.
  L2D        = .FALSE.
  L1D        = .FALSE.
  LFLAT      = .FALSE.
  NMODEL     = 1
  CEQNSYS    = 'DUR'
  NVERB      = 5
  CEXP       = 'EXP01'
  CSEG       = 'SEG01'
  LFORCING   = .FALSE.
  LPACK      = .TRUE.
  NHALO      = 1
  CSPLIT     ='YSPLITTING'
  NZ_PROC    = 1               !JUAN Z_SPLITTING :: number of proc in Z splitting
  NZ_SPLITTING = 10            !JUAN Z_SPLITTING :: for debug NZ=1=flat_inv;  NZ=10=flat_invz; NZ=1+2 the two
  LLG        = .FALSE.
  LINIT_LG   = .TRUE.
  CINIT_LG   = 'FMOUT'
  LNOMIXLG   = .FALSE.
END IF                                                           
!
CCLOUD    = 'NONE'
LUSERV    = .TRUE.
LUSERC    = .FALSE.
LUSERR    = .FALSE.  
LUSERI    = .FALSE.   
LUSERS    = .FALSE.  
LUSERG    = .FALSE.   
LUSERH    = .FALSE.  
!NSV       = 0
!NSV_USER   = 0
LUSECI    = .FALSE.   
!
!-------------------------------------------------------------------------------
!
!*      3.    SET DEFAULT VALUES FOR MODD_DYN AND MODD_DYN_n :
!             -----------------------------------------------
!
IF (KMI == 1) THEN                                             
  XSEGLEN   = 43200.
  XASSELIN  = 0.2
  XASSELIN_SV = 0.02
  LCORIO    = .TRUE.
  LNUMDIFF  = .FALSE.
  XALZBOT   = 4000.
  XALKTOP   = 0.01
END IF                                                        
!
XTSTEP    = 60.
CPRESOPT  = 'CRESI'
NITR      = 4
LITRADJ   = .FALSE.
XRELAX    = 1.
LVE_RELAX = .FALSE.
XRIMKMAX  = 0.01 / XTSTEP
XT4DIFF   = 1800.          
!
IF (KMI == 1) THEN   ! for model 1 we have a Large scale information       
  NRIMX = JPRIMMAX   ! for U,V,W,TH,Rv used for the hor. relaxation
  NRIMY = JPRIMMAX 
ELSE
  NRIMX = 0          ! for inner models we use only surfacic fields to
  NRIMY = 0          ! give the lbc and no hor. relaxation is used
END IF
!
LHORELAX_UVWTH = .FALSE.
LHORELAX_RV = .FALSE.
LHORELAX_RC = .FALSE. ! for all these fields, no large scale is usally available
LHORELAX_RR = .FALSE. ! for model 1 and for inner models, we only use surfacic
LHORELAX_RS = .FALSE. ! fiels ( no hor. relax. )
LHORELAX_RI = .FALSE.
LHORELAX_RG = .FALSE.
LHORELAX_RH = .FALSE.
LHORELAX_TKE = .FALSE.
LHORELAX_SV(:) = .FALSE.
LHORELAX_SVC2R2 = .FALSE.
LHORELAX_SVC1R3 = .FALSE.
LHORELAX_SVELEC = .FALSE.
LHORELAX_SVLG   = .FALSE.
LHORELAX_SVCHEM = .FALSE.
LHORELAX_SVDST  = .FALSE.
LHORELAX_SVSLT  = .FALSE.
LHORELAX_SVAER  = .FALSE.
!
!-------------------------------------------------------------------------------
!
!*      4.    SET DEFAULT VALUES FOR MODD_NESTING :
!             -----------------------------------
!
IF (KMI == 1) THEN
  NDAD(1)=1
  DO JM=2,JPMODELMAX
    NDAD(JM)  = JM - 1
  END DO
  NDTRATIO(:) = 1
  XWAY(:)     = 2.      ! two-way interactive gridnesting
  XWAY(1)     = 0.      ! except for model 1
END IF
!
!-------------------------------------------------------------------------------
!
!*      5.    SET DEFAULT VALUES FOR MODD_ADV_n :
!             ----------------------------------
!
CUVW_ADV_SCHEME =  'CEN4TH'
CMET_ADV_SCHEME =  'FCT2ND'
CSV_ADV_SCHEME  =  'FCT2ND'
NLITER    = 2
!
!-------------------------------------------------------------------------------
!
!*      6.    SET DEFAULT VALUES FOR MODD_PARAM_n :
!             -----------------------------------
!
CTURB   = 'NONE'
CRAD    = 'NONE'
CDCONV  = 'NONE'
CELEC   = 'NONE'
!
!-------------------------------------------------------------------------------
!
!*      7.    SET DEFAULT VALUES FOR MODD_LBC_n :
!             ---------------------------------
!
CLBCX(1) ='CYCL'
CLBCX(2) ='CYCL'
CLBCY(1) ='CYCL'
CLBCY(2) ='CYCL'
NLBLX(:) = 1
NLBLY(:) = 1
XCPHASE = 20.
!
!-------------------------------------------------------------------------------
!
!*      8.    SET DEFAULT VALUES FOR MODD_NUDGING_n :
!             ---------------------------------
!
LNUDGING = .FALSE.
XTNUDGING = 21600.
!
!-------------------------------------------------------------------------------
!
!*      9.    SET DEFAULT VALUES FOR MODD_FMOUT and MODD_OUT_n :
!             ------------------------------------------------
!
IF (KMI == 1) XFMOUT (:,:) = XUNDEF
!
!
!-------------------------------------------------------------------------------
!
!*      10.   SET DEFAULT VALUES FOR MODD_TURB_n :
!             ----------------------------------
!
XIMPL     = 1.
CTURBLEN  = 'BL89'
CTURBDIM  = '1DIM'
LTURB_FLX =.FALSE.
LTURB_DIAG=.FALSE.
LSUBG_COND=.FALSE.
LSUBG_AUCV=.FALSE.
LSIGMAS   =.TRUE. 
LSIG_CONV=.FALSE.
LSIG_CONV=.TRUE.
LRMC01   =.FALSE.
CTOM     ='NONE'
!
!-------------------------------------------------------------------------------
!
!*      11.   SET DEFAULT VALUES FOR MODD_BUDGET :
!             ------------------------------------
!
!       11.1 General  budget variables
!
IF (KMI == 1) THEN 
  CBUTYPE = 'NONE'
  NBUMOD = 1
  XBULEN = XSEGLEN
  XBUWRI = XSEGLEN
  NBUKL  = 1
  NBUKH  = 0         
  LBU_KCP = .TRUE.
!
!       11.2 Variables for the cartesian box
!
  NBUIL = 1
  NBUIH = 0      
  NBUJL = 1
  NBUJH = 0          
  LBU_ICP = .TRUE.
  LBU_JCP = .TRUE.
!
!       11.3 Variables for the mask
!
  NBUMASK = 1
!
!       11.4 Variables for budget and processes choice
  LBU_RU = .FALSE.
  NASSEU  = 0
  NNESTU  = 0
  NADVXU  = 0
  NADVYU  = 0
  NADVZU  = 0
  NFRCU   = 0
  NNUDU   = 0
  NCURVU  = 0
  NCORU   = 0
  NDIFU   = 0
  NRELU   = 0
  NVTURBU = 0
  NHTURBU = 0
  NPRESU  = 0
!
!                    Budget of RV
  LBU_RV = .FALSE.
  NASSEV  = 0
  NNESTV  = 0
  NADVXV  = 0
  NADVYV  = 0
  NADVZV  = 0
  NFRCV   = 0
  NNUDV   = 0
  NCURVV  = 0
  NCORV   = 0
  NDIFV   = 0
  NRELV   = 0
  NVTURBV = 0
  NHTURBV = 0
  NPRESV  = 0
!
!                    Budget of RW
  LBU_RW = .FALSE.
  NASSEW  = 0
  NNESTW  = 0
  NADVXW  = 0
  NADVYW  = 0
  NADVZW  = 0
  NFRCW   = 0
  NNUDW   = 0
  NCURVW  = 0
  NCORW   = 0
  NGRAVW  = 0
  NDIFW   = 0
  NRELW   = 0
  NVTURBW = 0
  NHTURBW = 0
  NPRESW  = 0
!
!                    Budget of RTH
  LBU_RTH = .FALSE.
  NASSETH  = 0
  NNESTTH  = 0
  NADVTH   = 0
  NADVXTH  = 0
  NADVYTH  = 0
  NADVZTH  = 0
  NFRCTH   = 0
  NNUDTH   = 0
  NPREFTH  = 0
  NDIFTH   = 0
  NRELTH   = 0
  NRADTH   = 0
  NDCONVTH = 0
  NVTURBTH = 0
  NHTURBTH = 0
  NDISSHTH = 0
  NNEGATH  = 0
  NREVATH  = 0
  NCONDTH  = 0
  NHENUTH  = 0
  NHONTH   = 0
  NSFRTH   = 0
  NDEPSTH  = 0
  NDEPGTH  = 0
  NREVATH  = 0
  NRIMTH   = 0
  NACCTH   = 0
  NCFRZTH  = 0
  NWETGTH  = 0
  NDRYGTH  = 0
  NGMLTTH  = 0
  NIMLTTH  = 0
  NBERFITH = 0
  NCDEPITH = 0
  NWETHTH = 0
  NHMLTTH = 0
!
!                    Budget of RTKE
  LBU_RTKE = .FALSE.
  NASSETKE = 0
  NADVTKE  = 0
  NADVXTKE = 0
  NADVYTKE = 0
  NADVZTKE = 0
  NFRCTKE  = 0
  NDIFTKE  = 0
  NRELTKE  = 0
  NDPTKE   = 0
  NTPTKE   = 0
  NDISSTKE = 0
  NTRTKE   = 0
!
!                    Budget of RRV
  LBU_RRV = .FALSE.
  NASSERV  = 0
  NNESTRV  = 0
  NADVRV   = 0
  NADVXRV  = 0
  NADVYRV  = 0
  NADVZRV  = 0
  NFRCRV   = 0
  NNUDRV   = 0
  NDIFRV   = 0
  NRELRV   = 0
  NDCONVRV = 0
  NVTURBRV = 0
  NHTURBRV = 0
  NNEGARV  = 0
  NREVARV  = 0
  NCONDRV  = 0
  NHENURV  = 0
  NDEPSRV  = 0
  NDEPGRV  = 0
  NREVARV  = 0
  NCDEPIRV = 0
!
!                    Budget of RRC
  LBU_RRC = .FALSE.
  NASSERC  = 0
  NNESTRC  = 0
  NADVRC   = 0
  NADVXRC  = 0
  NADVYRC  = 0
  NADVZRC  = 0
  NFRCRC   = 0
  NDIFRC   = 0
  NRELRC   = 0
  NDCONVRC = 0
  NVTURBRC = 0
  NHTURBRC = 0
  NNEGARC  = 0
  NACCRRC  = 0
  NAUTORC  = 0
  NCONDRC  = 0
  NAUTORC  = 0
  NACCRRC  = 0
  NHONRC   = 0
  NRIMRC   = 0
  NWETGRC  = 0
  NDRYGRC  = 0
  NIMLTRC  = 0
  NBERFIRC = 0
  NCDEPIRC = 0
  NWETHRC = 0
!
!                    Budget of RRR
  LBU_RRR = .FALSE.
  NASSERR  = 0
  NNESTRR  = 0
  NADVRR   = 0
  NADVXRR  = 0
  NADVYRR  = 0
  NADVZRR  = 0
  NFRCRR   = 0
  NDIFRR   = 0
  NRELRR   = 0
  NNEGARR  = 0
  NACCRRR  = 0
  NAUTORR  = 0
  NREVARR  = 0
  NSEDIRR  = 0
  NSFRRR   = 0
  NACCRR   = 0
  NCFRZRR  = 0
  NWETGRR  = 0
  NDRYGRR  = 0
  NGMLTRR  = 0
  NWETHRR  = 0
  NHMLTRR  = 0
!
!                    Budget of RRI
  LBU_RRI = .FALSE.
  NASSERI  = 0
  NNESTRI  = 0
  NADVRI   = 0
  NADVXRI  = 0
  NADVYRI  = 0
  NADVZRI  = 0
  NFRCRI   = 0
  NDIFRI   = 0
  NRELRI   = 0
  NDCONVRI = 0
  NVTURBRI = 0
  NHTURBRI = 0
  NNEGARI  = 0
  NSEDIRI  = 0
  NHENURI  = 0
  NHONRI   = 0
  NAGGSRI  = 0
  NAUTSRI  = 0
  NCFRZRI  = 0
  NWETGRI  = 0
  NDRYGRI  = 0
  NIMLTRI  = 0
  NBERFIRI = 0
  NCDEPIRI = 0
  NWETHRI = 0
!
!                    Budget of RRS
  LBU_RRS = .FALSE.
  NASSERS  = 0
  NNESTRS  = 0
  NADVRS   = 0
  NADVXRS  = 0
  NADVYRS  = 0
  NADVZRS  = 0
  NFRCRS   = 0
  NDIFRS   = 0
  NRELRS   = 0
  NNEGARS  = 0
  NSEDIRS  = 0
  NDEPSRS  = 0
  NAGGSRS  = 0
  NAUTSRS  = 0
  NRIMRS   = 0
  NACCRS   = 0
  NCMELRS  = 0
  NWETGRS  = 0
  NDRYGRS  = 0
  NWETHRS  = 0
!
!                    Budget of RRG
  LBU_RRG = .FALSE.
  NASSERG  = 0
  NNESTRG  = 0
  NADVRG   = 0
  NADVXRG  = 0
  NADVYRG  = 0
  NADVZRG  = 0
  NFRCRG   = 0
  NDIFRG   = 0
  NRELRG   = 0
  NNEGARG  = 0
  NSEDIRG  = 0
  NSFRRG   = 0
  NDEPGRG  = 0
  NRIMRG   = 0
  NACCRG   = 0
  NCMELRG  = 0
  NCFRZRG  = 0
  NWETGRG  = 0
  NDRYGRG  = 0
  NGMLTRG  = 0
  NWETHRG  = 0
!
!                    Budget of RRH
  LBU_RRH = .FALSE.
  NASSERH  = 0
  NNESTRH  = 0
  NADVRH   = 0
  NADVXRH  = 0
  NADVYRH  = 0
  NADVZRH  = 0
  NFRCRH   = 0
  NDIFRH   = 0
  NRELRH   = 0
  NNEGARH  = 0
  NWETGRH  = 0
  NWETHRH  = 0
  NHMLTRH  = 0
!
!                    Budget of RSVx
  LBU_RSV = .FALSE.
  NASSESV  = 0
  NNESTSV  = 0
  NADVSV   = 0
  NADVXSV  = 0
  NADVYSV  = 0
  NADVZSV  = 0
  NFRCSV   = 0
  NDIFSV   = 0
  NRELSV   = 0
  NDCONVSV = 0
  NVTURBSV = 0
  NHTURBSV = 0
  NCHEMSV  = 0
!
!
END IF
!
!-------------------------------------------------------------------------------
!
!*      12.    SET DEFAULT VALUES FOR MODD_LES :
!             ---------------------------------
!
IF (KMI == 1) THEN 
  LLES_MEAN               = .FALSE.
  LLES_RESOLVED           = .FALSE.
  LLES_SUBGRID            = .FALSE.
  LLES_UPDRAFT            = .FALSE.
  LLES_DOWNDRAFT          = .FALSE.
  LLES_SPECTRA            = .FALSE.
!
  NLES_LEVELS             = NUNDEF
  XLES_ALTITUDES          = XUNDEF
  NSPECTRA_LEVELS         = NUNDEF
  XSPECTRA_ALTITUDES      = XUNDEF
  NLES_TEMP_SERIE_I       = NUNDEF
  NLES_TEMP_SERIE_J       = NUNDEF
  NLES_TEMP_SERIE_Z       = NUNDEF
  CLES_NORM_TYPE          = 'NONE'
  CBL_HEIGHT_DEF          = 'KE'
  XLES_TEMP_SAMPLING      = XUNDEF
  XLES_TEMP_MEAN_START    = XUNDEF
  XLES_TEMP_MEAN_END      = XUNDEF
  XLES_TEMP_MEAN_STEP     = 3600.
  LLES_CART_MASK          = .FALSE.
  NLES_IINF               = NUNDEF
  NLES_ISUP               = NUNDEF
  NLES_JINF               = NUNDEF
  NLES_JSUP               = NUNDEF
  LLES_NEB_MASK           = .FALSE.
  LLES_CORE_MASK          = .FALSE.
  LLES_MY_MASK            = .FALSE.
!
  LLES_PDF               = .FALSE.
  NPDF                   = 1
  XTH_PDF_MIN            = 270.
  XTH_PDF_MAX            = 350.
  XW_PDF_MIN             = -10.
  XW_PDF_MAX             = 10.
  XTHV_PDF_MIN           = 270.
  XTHV_PDF_MAX           = 350.
  XRV_PDF_MIN            = 0.
  XRV_PDF_MAX            = 20.
  XRC_PDF_MIN            = 0.
  XRC_PDF_MAX            = 1.
  XRR_PDF_MIN            = 0.
  XRR_PDF_MAX            = 1.
  XRI_PDF_MIN            = 0.
  XRI_PDF_MAX            = 1.
  XRS_PDF_MIN            = 0.
  XRS_PDF_MAX            = 1.
  XRG_PDF_MIN            = 0.
  XRG_PDF_MAX            = 1.
  XRT_PDF_MIN            = 0.
  XRT_PDF_MAX            = 20.
  XTHL_PDF_MIN           = 270.
  XTHL_PDF_MAX           = 350.
END IF
!
!-------------------------------------------------------------------------------
!
!*      13.   SET DEFAULT VALUES FOR MODD_PARAM_RAD_n :
!             ---------------------------------------
!
XDTRAD        = XTSTEP
XDTRAD_CLONLY = XTSTEP
LCLEAR_SKY    =.FALSE.
NRAD_COLNBR   = 1000
NRAD_DIAG     = 0
CLW ='RRTM'
CAER='SURF'
CEFRADL='MART'
CEFRADI='LIOU'
COPWSW = 'FOUQ'
COPISW = 'EBCU'
COPWLW = 'SMSH'
COPILW = 'EBCU'
XFUDG = 1.
!
!-------------------------------------------------------------------------------
!
!*      14.   SET DEFAULT VALUES FOR MODD_BLANK :
!             -----------------------------------
!
IF (KMI == 1) THEN 
  XDUMMY1       = 0.
  XDUMMY2       = 0.
  XDUMMY3       = 0.
  XDUMMY4       = 0.
  XDUMMY5       = 0.
  XDUMMY6       = 0.
  XDUMMY7       = 0.
  XDUMMY8       = 0.
  XDUMMY=0.
!
  NDUMMY1       = 0
  NDUMMY2       = 0
  NDUMMY3       = 0
  NDUMMY4       = 0
  NDUMMY5       = 0
  NDUMMY6       = 0
  NDUMMY7       = 0
  NDUMMY8       = 0
  NDUMMY=0
!
  LDUMMY1       = .TRUE.
  LDUMMY2       = .TRUE.
  LDUMMY3       = .TRUE.
  LDUMMY4       = .TRUE.
  LDUMMY5       = .TRUE.
  LDUMMY6       = .TRUE.
  LDUMMY7       = .TRUE.
  LDUMMY8       = .TRUE.
  LDUMMY=.TRUE.
!
  CDUMMY1       = ' '
  CDUMMY2       = ' '
  CDUMMY3       = ' '
  CDUMMY4       = ' '
  CDUMMY5       = ' '
  CDUMMY6       = ' '
  CDUMMY7       = ' '
  CDUMMY8       = ' '
  CDUMMY= ' '
END IF
!
!------------------------------------------------------------------------------
!
!*      15.   SET DEFAULT VALUES FOR MODD_FRC :
!             ---------------------------------
!
IF (KMI == 1) THEN 
  LGEOST_UV_FRC      = .FALSE.
  LGEOST_TH_FRC      = .FALSE.
  LTEND_THRV_FRC      = .FALSE.
  LVERT_MOTION_FRC   = .FALSE.
  LRELAX_THRV_FRC    = .FALSE.
  LRELAX_UV_FRC      = .FALSE.
  XRELAX_TIME_FRC    = 10800.
  XRELAX_HEIGHT_FRC  = 0.
  CRELAX_HEIGHT_TYPE = "FIXE"
  LTRANS             = .FALSE.
  XUTRANS            = 0.0
  XVTRANS            = 0.0
  LPGROUND_FRC       = .FALSE.
END IF
!
!-------------------------------------------------------------------------------
!
!
!*      16.   SET DEFAULT VALUES FOR MODD_PARAM_ICE :
!             ---------------------------------------
!
IF (KMI == 1) THEN 
  LWARM = .TRUE.
  CPRISTINE_ICE = 'PLAT'
  LSEDIC = .FALSE.
END IF
!
!-------------------------------------------------------------------------------
!
!
!*      17.   SET DEFAULT VALUES FOR MODD_PARAM_CONVECT_n :
!             --------------------------------------------
!
XDTCONV       = MAX( 300.0,XTSTEP )
NICE          = 1
LREFRESH_ALL  = .TRUE.
LCHTRANS      = .FALSE.
LDOWN         = .TRUE.
LDEEP         = .TRUE.
LSHAL         = .FALSE.
LSETTADJ      = .FALSE.
XTADJD        = 3600.
XTADJS        = 10800.
LDIAGCONV     = .FALSE.
NENSM         = 0
!
!-------------------------------------------------------------------------------
!
!*      18.   SET DEFAULT VALUES FOR MODD_PARAM_C2R2 :
!             ----------------------------------------
!
IF (KMI == 1) THEN
  XNUC    = 1.0
  XALPHAC = 3.0
  XNUR    = 2.0
  XALPHAR = 1.0
!
  LRAIN = .TRUE.
  LSEDC = .TRUE.
  LACTIT = .FALSE.
!
  HPARAM_CCN = 'XXX'
  HINI_CCN   = 'XXX'
  HTYPE_CCN  = 'X'
!
  XCHEN      = 0.0
  XKHEN      = 0.0
  XMUHEN     = 0.0
  XBETAHEN   = 0.0
!
  XCONC_CCN   = 0.0
  XAERDIFF    = 0.0
  XAERHEIGHT  = 2000
  XR_MEAN_CCN = 0.0
  XLOGSIG_CCN = 0.0
  XFSOLUB_CCN = 1.0
  XACTEMP_CCN = 280.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      19.   SET DEFAULT VALUES FOR MODD_CH_MNHC_n
!             -------------------------------------
!
LUSECHEM            = .FALSE.
LCH_INIT_FIELD      = .FALSE.  
LCH_SURFACE_FLUX    = .FALSE.
LCH_CONV_SCAV       = .FALSE. 
LCH_CONV_LINOX      = .FALSE.
CCH_EXPLICIT_SCAV   = 'IHLEKESS'
CCHEM_INPUT_FILE    = 'EXSEG1.nam'
CCH_TDISCRETIZATION = 'SPLIT'  
NCH_SUBSTEPS        = 1
LCH_TUV_ONLINE      = .FALSE.      
CCH_TUV_LOOKUP      = 'PHOTO.TUV39' 
CCH_TUV_CLOUDS      = 'NONE'        
XCH_TUV_ALBNEW      = -1.  
XCH_TUV_DOBNEW      = -1.  
XCH_TUV_TUPDATE     = 600. 
CCH_VEC_METHOD      = 'MAX' 
NCH_VEC_LENGTH      = 1000   
XCH_TS1D_TSTEP      = 600.  
CCH_TS1D_COMMENT    = 'no comment' 
CCH_TS1D_FILENAME   = 'IO1D' 
!
!-------------------------------------------------------------------------------
!
!*      20.   SET DEFAULT VALUES FOR MODD_SERIES AND MODD_SERIE_n
!             ---------------------------------------------------
!
IF (KMI == 1) THEN
  LSERIES      = .FALSE.
  LMASKLANDSEA = .FALSE.
  LWMINMAX     = .FALSE.
ENDIF
!
NIBOXL = 1 + JPHEXT
NIBOXH = 1 + 2*JPHEXT  
NJBOXL = 1 + JPHEXT
NJBOXH = 1 + 2*JPHEXT  
NKCLS  = 1 + JPVEXT
NKLOW  = 1 + JPVEXT
NKMID  = 1 + JPVEXT
NKUP   = 1 + JPVEXT
NKCLA  = 1 + JPVEXT
NBJSLICE = 1
NJSLICEL(:) = 1 + JPHEXT
NJSLICEH(:) = 1 + 2*JPHEXT  
NFREQSERIES  = INT(XSEGLEN /(100.*XTSTEP) )
NFREQSERIES  = MAX(NFREQSERIES,1)
!
!-------------------------------------------------------------------------------
!
!*      21.   SET DEFAULT VALUES FOR MODD_TURB_CLOUD
!             --------------------------------------
!
IF (KMI == 1) THEN
  NMODEL_CLOUD = NUNDEF
  CTURBLEN_CLOUD = 'DELT'
  XCOEF_AMPL_SAT = 5.
  XCEI_MIN = 0.001E-06
  XCEI_MAX = 0.01E-06
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      22.   SET DEFAULT VALUES FOR MODD_AEROSOL
!             -----------------------------------
IF (KMI == 1) THEN ! other values are defined in modd_ch_aerosol
!
! aerosol lognormal parameterization
CRGUNIT   = 'MASS'    ! type of log-normal geometric mean radius given
!                                     ! in nameliste (mass on number)

LVARSIGI  = .FALSE.   ! switch to active pronostic dispersion for I mode
LVARSIGJ  = .FALSE.   ! switch to active pronostic dispersion for J mode
LHETEROSO4 = .FALSE.  ! switch to active sulfates heteronegeous
                      ! production
LSEDIMAERO = .FALSE.  ! switch to active aerosol sedimentation
LTHERMIJ   = .FALSE.  ! switch to active thermodynamics
                      ! mineral equilibrium for each mode
CMINERAL      = "NONE"   ! mineral equilibrium scheme
CORGANIC      = "NONE"   ! mineral equilibrium scheme
CNUCLEATION   = "NONE" ! sulfates nucleation scheme

ENDIF

!*      23.   SET DEFAULT VALUES FOR MODD_DUST
!             --------------------------------
!
IF (KMI == 1) THEN ! other values initialized in modd_dust
  LDUST      = .FALSE.
  NMODE_DST  = 3
  CRGUNITD   = 'MASS'
  LVARSIG    = .FALSE.
  LSEDIMDUST = .FALSE.
ENDIF
!

!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DESFM_n
