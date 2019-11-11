!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!#######################
MODULE MODI_SPAWN_FIELD2
!#######################
!
INTERFACE
!
      SUBROUTINE SPAWN_FIELD2(KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,HTURB,   &
               PUT,PVT,PWT,PTHVT,PRT,PHUT,PTKET,PSVT,PZWS,PATC,                &
               PSRCT,PSIGS,                                                    &
               PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,                        &
               PDTHFRC,PDRVFRC,PTHREL,PRVREL,                                  &
               PVU_FLUX_M,PVTH_FLUX_M,PWTH_FLUX_M,                             &
               TPSONFILE,KIUSON,KJUSON,                                        &
               KIB2,KJB2,KIE2,KJE2,                                            &
               KIB1,KJB1,KIE1,KJE1                                             )
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
CHARACTER (LEN=4), INTENT(IN) :: HTURB !  Kind of turbulence parameterization
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PUT,PVT,PWT        !  model 2
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PTKET              ! variables
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PRT,PSVT,PATC      !   at t
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PTHVT,PHUT         !
REAL, DIMENSION(:,:),     INTENT(OUT) :: PZWS
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PSRCT,PSIGS  ! secondary
                                                            ! prognostic variables
           ! Larger Scale fields for relaxation and diffusion
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSUM, PLSVM, PLSWM 
REAL, DIMENSION(:,:),            INTENT(OUT) :: PLSZWSM
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSTHM,  PLSRVM     
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PDTHFRC,PDRVFRC
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PTHREL,PRVREL
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PVU_FLUX_M,PVTH_FLUX_M,PWTH_FLUX_M
!
           ! Arguments for spawning with 2 input files (father+son1)
TYPE(TFILEDATA),   OPTIONAL, INTENT(IN) :: TPSONFILE ! input FM-file SON
INTEGER,           OPTIONAL, INTENT(IN) :: KIUSON  ! upper dimensions of the
INTEGER,           OPTIONAL, INTENT(IN) :: KJUSON  !input FM-file SON
INTEGER,           OPTIONAL, INTENT(IN) :: KIB2,KJB2 ! indexes for common
INTEGER,           OPTIONAL, INTENT(IN) :: KIE2,KJE2 !domain in model2
INTEGER,           OPTIONAL, INTENT(IN) :: KIB1,KJB1 !and in
INTEGER,           OPTIONAL, INTENT(IN) :: KIE1,KJE1 !SON
END SUBROUTINE SPAWN_FIELD2
!
END INTERFACE
!
END MODULE MODI_SPAWN_FIELD2
!     ##########################################################################
      SUBROUTINE SPAWN_FIELD2(KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,HTURB,   &
               PUT,PVT,PWT,PTHVT,PRT,PHUT,PTKET,PSVT, PZWS,PATC,                &
               PSRCT,PSIGS,                                                    &
               PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,                        &
               PDTHFRC,PDRVFRC,PTHREL,PRVREL,                                  &
               PVU_FLUX_M,PVTH_FLUX_M,PWTH_FLUX_M,                             &
               TPSONFILE,KIUSON,KJUSON,                                        &
               KIB2,KJB2,KIE2,KJE2,                                            &
               KIB1,KJB1,KIE1,KJE1                                             )
!     ##########################################################################
!
!!****  *SPAWN_FIELD2 * - subroutine generating the model 2 prognostic and LS
!!                      fields, consistently with the spawning model 1.
!!
!!    PURPOSE
!!    -------
!!
!!      The prognostic and LS fields are interpolated from the model 1, to 
!!    initialize the model 2.
!!
!!**  METHOD
!!    ------
!!
!!      The model 2 variables are transmitted by argument (P or K prefixes),
!!    while the ones of model 1 are declared through calls to MODD_... 
!!    (X or N prefixes)
!!
!!      For the case where the resolution ratio between models is 1, 
!!    the horizontal interpolation becomes a simple equality.
!!      For the general case where resolution ratio is not egal to one,
!!    fields are interpolated using 2 types of interpolations:
!!                 1. Clark and Farley (JAS 1984) on 9 points 
!!                 2. Bikhardt on 16 points
!!
!!    EXTERNAL
!!    --------
!!      
!!      Routine BIKHARDT      : to perform horizontal interpolations
!!      Routine CLARK_FARLEY  : to perform horizontal interpolations
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS : contains parameters 
!!      Module MODD_CONF       : contains NVERB
!!      Module MODD_CONF1      : contains CONF_MODEL(1)%NRR (total Number of moist variables)
!!      Module MODD_FIELD1     : contains pronostic variables of model 1
!!      Module MODD_LSFIELD1   : contains LB and LS variables of model 1
!!      Module MODD_REF1       : contains RHODJ of model 1
!!      Module MODD_GRID1      : contains grid variables
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       SUBROUTINE SPAWN_FIELD2 (Book2 of the documentation)
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       J.P. Lafore     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    12/01/95
!!      Modification 20/03/95 (I.Mallet) change Large Scale fields initialization 
!!      Modification 27/04/95 (    "   ) remove R from the historical variables 
!!      Modification 17/04/96  (Lafore) Different resolution ratio case introduction
!!      Modification 10/06/96 (V.Masson) remove the loops in case of no resolution change
!!                                       and bug in initialization of ZBFY
!!      Modification 10/06/96 (V.Masson) interpolation computations performed in
!!                                       independant routines
!!                   10/10/96 (J. Stein) add SRCM and SRCT
!!      Modification 21/11/96 (Lafore)   move from BIKHARDT2 to BIKHARDT routine
!!      Modification 21/11/96 (Lafore)   "surfacic" LS fields
!!      Modification 10/07/97 (Masson)   remove pressure interpolations
!!      Modification 17/07/97 (Masson)   add EPS and tests on other variables
!!      Modification 14/09/97 (Masson)   interpolation of relative humidity
!!      Modification 14/09/97 (J. Stein) add the LB and LS fields
!!      Modification 27/07/98 (P. Jabouille) compute HU for all the cases
!!      Modification 01/02/01 (D.Gazen)  add module MODD_NSV for NSV variable
!!      Modification 07/07/05 (D.Barbary) spawn with 2 input files (father+son1)
!!      Modification 05/06                Remove EPS, Clark and Farley
!!      Modification 06/12  (M.Tomasini)  Interpolation of turbulent fluxes (EDDY_FLUX)
!!                                        for 2D west african monsoon
!!      Modification 07/13  (Bosseur & Filippi) Adds Forefire
!!      Modification 2014 (M.Faivre)
!!      Modification 01/15  (C. Barthe)   add LNOx
!!      Modification 25/02/2015 (M.Moge) correction of the parallelization attempted by M.Faivre
!!      Modification 15/04/2016 (P.Tulet) bug allocation ZSVT_C
!!                   29/04/2016 (J.Escobar) bug in use of ZSVT_C in SET_LSFIELD_1WAY_ll        
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Modification 05/03/2018 (J.Escobar) bypass gridnesting special case KD(X/Y)RATIO == 1 not parallelized
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 14/03/2019: correct ZWS when variable not present in file
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_2D_FRC
USE MODD_ADVFRC_n
USE MODD_BIKHARDT_n
USE MODD_CH_AEROSOL,      ONLY: CAERONAMES
USE MODD_CH_M9_n,         ONLY: CNAMES, CICNAMES
USE MODD_CONF
USE MODD_CST
USE MODD_CONF_n,          ONLY:  CONF_MODEL
USE MODD_DUST,            ONLY: CDUSTNAMES
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES
USE MODD_FIELD_n,         ONLY: FIELD_MODEL, XZWS_DEFAULT
USE MODD_IO_ll,           ONLY : TFILEDATA
USE MODD_LATZ_EDFLX
USE MODD_LBC_n,           ONLY:  LBC_MODEL
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_LUNIT_n,         ONLY:  LUNIT_MODEL,TLUOUT
USE MODD_NSV
USE MODD_REF_n,           ONLY:  REF_MODEL
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA     , ONLY : NMOD_CCN, NMOD_IFN, NMOD_IMM, NINDICE_CCN_IMM,&
                                 LSCAV, LAERO_MASS, LHHONI
USE MODD_PARAM_LIMA_COLD, ONLY : CLIMA_COLD_NAMES
USE MODD_PARAM_LIMA_WARM, ONLY : CLIMA_WARM_NAMES, CAERO_MASS
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_RELFRC_n 
USE MODD_SALT,            ONLY: CSALTNAMES
USE MODD_SPAWN
!
USE MODE_FIELD,           ONLY: TFIELDDATA,TYPEREAL
USE MODE_FMREAD
USE MODE_IO_ll,           ONLY: UPCASE
USE MODE_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_THERMO
!
USE MODI_BIKHARDT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
CHARACTER (LEN=4), INTENT(IN) :: HTURB !  Kind of turbulence parameterization
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PUT,PVT,PWT        !  model 2
REAL, DIMENSION(:,:),     INTENT(OUT) :: PZWS
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PTKET              ! variables
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PRT,PSVT,PATC      !   at t
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PTHVT,PHUT         !
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PSRCT,PSIGS  ! secondary
                                                            ! prognostic variables
           ! Larger Scale fields for relaxation and diffusion
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSUM, PLSVM, PLSWM 
REAL, DIMENSION(:,:),            INTENT(OUT) :: PLSZWSM
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSTHM,  PLSRVM 
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PDTHFRC,PDRVFRC
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PTHREL,PRVREL
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PVU_FLUX_M,PVTH_FLUX_M,PWTH_FLUX_M
           ! Arguments for spawning with 2 input files (father+son1)
TYPE(TFILEDATA),   OPTIONAL, INTENT(IN) :: TPSONFILE ! input FM-file SON
INTEGER,           OPTIONAL, INTENT(IN) :: KIUSON  ! upper dimensions of the
INTEGER,           OPTIONAL, INTENT(IN) :: KJUSON  !input FM-file SON
INTEGER,           OPTIONAL, INTENT(IN) :: KIB2,KJB2 ! indexes for common
INTEGER,           OPTIONAL, INTENT(IN) :: KIE2,KJE2 !domain in model2
INTEGER,           OPTIONAL, INTENT(IN) :: KIB1,KJB1 !and in
INTEGER,           OPTIONAL, INTENT(IN) :: KIE1,KJE1 !SON
!
!*       0.2    Declarations of local variables 
!
INTEGER             :: ILUOUT    ! Logical unit number for the output listing 
INTEGER             :: IRESP     ! Return codes in FM routines
INTEGER             :: JRR,JSV   ! Loop index for moist and scalar variables 
INTEGER             :: IRR       ! Number of moist variables 
!
REAL, DIMENSION(SIZE(XRT1,1),SIZE(XRT1,2),SIZE(XRT1,3)) :: ZHUT ! relative humidity
                                                             ! (model 1)
REAL, DIMENSION(SIZE(XTHT1,1),SIZE(XTHT1,2),SIZE(XTHT1,3)) :: ZTHVT! virtual pot. T
                                                                ! (model 1)          
!$20140708
REAL, DIMENSION(:,:),   ALLOCATABLE   :: ZZWS_C, ZLSZWSM_C
!$***** 3D
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZUT_C, ZLSUM_C 
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZVT_C, ZLSVM_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZWT_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZTHVT_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZLSWM_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZLSTHM_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZLSRVM_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZTKET_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZHUT_C, ZSRCM_C, ZSRCT_C, ZSIGS_C
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZVU_FLUX_M_C, ZVTH_FLUX_M_C, ZWTH_FLUX_M_C
!$***** 4D
REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: ZSVT_C
REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: ZRT_C, ZDTHFRC_C, ZDRVFRC_C
REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: ZTHREL_C, ZRVREL_C
!$                    
INTEGER  :: IMI, JI,KI
!$20140708
INTEGER  :: IDIMX_C, IDIMY_C
INTEGER  :: IINFO_ll
!$
! Arrays for reading fields of input SON 1 file
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZWORK3D
REAL, DIMENSION(:,:),   ALLOCATABLE   :: ZWORK2D
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZTHT1,ZTHVT1
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZPABST1,ZHUT1
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZRT1
LOGICAL :: GUSERV
!
CHARACTER(LEN=15) :: YVAL
CHARACTER(LEN=2)  :: INDICE
INTEGER           :: I
TYPE(TFIELDDATA)             :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE:
!              ---------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(2)
CALL GO_TOMODEL_ll(2, IINFO_ll)
!
!*       1.0  recovers logical unit number of output listing
!
ILUOUT = TLUOUT%NLU
!
!*       1.1   Secondary variables
!
CALL COMPUTE_THV_HU(CONF_MODEL(1)%LUSERV,XRT1,XTHT1,XPABST1,ZTHVT,ZHUT)
!
!*       1.2   Working arrays for reading in SON input file
!
IF (PRESENT(TPSONFILE)) THEN
  ALLOCATE(ZWORK3D(KIUSON,KJUSON,SIZE(PUT,3)))
  ALLOCATE(ZWORK2D(KIUSON,KJUSON))
  ALLOCATE(ZPABST1(KIE1-KIB1+1,KJE1-KJB1+1,SIZE(PUT,3)))
  ALLOCATE(ZTHT1(KIE1-KIB1+1,KJE1-KJB1+1,SIZE(PUT,3)))
  ALLOCATE(ZTHVT1(KIE1-KIB1+1,KJE1-KJB1+1,SIZE(PUT,3)))
  IF (CONF_MODEL(1)%NRR /= 0) THEN
    ALLOCATE(ZHUT1(KIE1-KIB1+1,KJE1-KJB1+1,SIZE(PUT,3)))
    ALLOCATE(ZRT1(KIE1-KIB1+1,KJE1-KJB1+1, SIZE(PUT,3),SIZE(PRT,4)))
  END IF
END IF
! 
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZATION OF PROGNOSTIC AND LS VARIABLES OF MODEL 2:
!              ---------------------------------------------------------
! 
!
!!$IF (KDXRATIO == 1 .AND. KDYRATIO == 1 ) THEN
!!$!
!!$!*       2.1   special case of spawning - no change of resolution :
!!$!
!!$!*       2.1.1  variables which always exist
!!$!
!!$  PUT  (:,:,:)   =  FIELD_MODEL(1)%XUT  (KXOR:KXEND,KYOR:KYEND,:)
!!$  PVT  (:,:,:)   =  FIELD_MODEL(1)%XVT  (KXOR:KXEND,KYOR:KYEND,:)
!!$  PWT  (:,:,:)   =  FIELD_MODEL(1)%XWT  (KXOR:KXEND,KYOR:KYEND,:)
!!$  PTHVT(:,:,:)   =  ZTHVT(KXOR:KXEND,KYOR:KYEND,:)
!!$!
!!$  PLSUM (:,:,:)  =  FIELD_MODEL(1)%XUT (KXOR:KXEND,KYOR:KYEND,:)
!!$  PLSVM (:,:,:)  =  FIELD_MODEL(1)%XVT (KXOR:KXEND,KYOR:KYEND,:)
!!$  PLSWM (:,:,:)  =  FIELD_MODEL(1)%XWT (KXOR:KXEND,KYOR:KYEND,:)
!!$  PLSTHM(:,:,:)  =  FIELD_MODEL(1)%XTHT(KXOR:KXEND,KYOR:KYEND,:)
!!$!
!!$  PLSRVM(:,:,:)  = 0.
!!$!
!!$!$20140707
!!$CALL MPPDB_CHECK3D(PUT,"SPAWN_FIELD2:PUT",PRECISION)
!!$CALL MPPDB_CHECK3D(PVT,"SPAWN_FIELD2:PVT",PRECISION)
!!$!$
!!$!*       2.1.2  TKE variable
!!$!
!!$  IF (HTURB /= 'NONE') THEN
!!$    PTKET(:,:,:)   =  FIELD_MODEL(1)%XTKET(KXOR:KXEND,KYOR:KYEND,:)
!!$  ENDIF
!!$!
!!$!*       2.1.3  moist variables
!!$!
!!$  IF (CONF_MODEL(1)%NRR /= 0) THEN
!!$    PRT  (:,:,:,:) =  FIELD_MODEL(1)%XRT  (KXOR:KXEND,KYOR:KYEND,:,:)
!!$    PLSRVM(:,:,:)  =  FIELD_MODEL(1)%XRT  (KXOR:KXEND,KYOR:KYEND,:,1)
!!$    PHUT (:,:,:)   =  ZHUT (KXOR:KXEND,KYOR:KYEND,:)
!!$  ENDIF
!!$!
!!$!*       2.1.4  scalar variables
!!$!
!!$  IF (NSV /= 0) THEN
!!$    PSVT (:,:,:,:) =  FIELD_MODEL(1)%XSVT (KXOR:KXEND,KYOR:KYEND,:,:)
!!$  ENDIF
!!$!
!!$!*       2.1.5  secondary prognostic variables
!!$!
!!$  IF (CONF_MODEL(1)%NRR > 1) THEN
!!$    PSRCT (:,:,:) =  FIELD_MODEL(1)%XSRCT (KXOR:KXEND,KYOR:KYEND,:)
!!$    PSIGS(:,:,:) =  FIELD_MODEL(1)%XSIGS(KXOR:KXEND,KYOR:KYEND,:)
!!$  ENDIF
!!$!
!!$!*       2.1.6  Large scale variables
!!$!
!!$  PLSUM  (:,:,:)   =  LSFIELD_MODEL(1)%XLSUM  (KXOR:KXEND,KYOR:KYEND,:)
!!$  PLSVM  (:,:,:)   =  LSFIELD_MODEL(1)%XLSVM  (KXOR:KXEND,KYOR:KYEND,:)
!!$  PLSWM  (:,:,:)   =  LSFIELD_MODEL(1)%XLSWM  (KXOR:KXEND,KYOR:KYEND,:)
!!$  PLSTHM(:,:,:)    =  LSFIELD_MODEL(1)%XLSTHM (KXOR:KXEND,KYOR:KYEND,:)
!!$  IF ( CONF_MODEL(1)%NRR > 0 ) THEN
!!$    PLSRVM  (:,:,:)   =  LSFIELD_MODEL(1)%XLSRVM  (KXOR:KXEND,KYOR:KYEND,:) 
!!$  END IF
!!$!
!!$!*       2.1.7  Advective forcing fields for 2D (Modif MT)
!!$!
!!$  IF (L2D_ADV_FRC) THEN
!!$    PDTHFRC(:,:,:,:)= ADVFRC_MODEL(1)%XDTHFRC (KXOR:KXEND,KYOR:KYEND,:,:)
!!$    PDRVFRC(:,:,:,:)= ADVFRC_MODEL(1)%XDRVFRC (KXOR:KXEND,KYOR:KYEND,:,:)
!!$  ENDIF
!!$  IF (L2D_REL_FRC) THEN
!!$    PTHREL(:,:,:,:)= RELFRC_MODEL(1)%XTHREL (KXOR:KXEND,KYOR:KYEND,:,:)
!!$    PRVREL(:,:,:,:)= RELFRC_MODEL(1)%XRVREL (KXOR:KXEND,KYOR:KYEND,:,:)
!!$  ENDIF
!!$!
!!$!*       2.1.8  Turbulent fluxes for 2D (Modif MT)                                    
!!$!
!!$  IF (LUV_FLX) THEN
!!$    PVU_FLUX_M(:,:,:)= EDDYUV_FLUX_MODEL(1)%XVU_FLUX_M (KXOR:KXEND,KYOR:KYEND,:)
!!$  END IF
!!$!
!!$  IF (LTH_FLX) THEN
!!$    PVTH_FLUX_M(:,:,:)= EDDY_FLUX_MODEL(1)%XVTH_FLUX_M (KXOR:KXEND,KYOR:KYEND,:)
!!$    PWTH_FLUX_M(:,:,:)= EDDY_FLUX_MODEL(1)%XWTH_FLUX_M (KXOR:KXEND,KYOR:KYEND,:)
!!$  END IF
!!$!
!!$!-------------------------------------------------------------------------------
!!$!
!!$ELSE
!
!-------------------------------------------------------------------------------
!
!*       2.2  general case - change of resolution :
!             -----------------------------------
!
!$20140708 get XDIM, YDIM = G2^G1@resol1
  CALL GOTO_MODEL(1)
  CALL GO_TOMODEL_ll(1, IINFO_ll)
  CALL GET_CHILD_DIM_ll(2, IDIMX_C, IDIMY_C, IINFO_ll)
!
!$20140708 use  ZTHVM_C in BIKAT top cal PTHVM_C
  ALLOCATE(ZZWS_C(IDIMX_C,IDIMY_C))
  ALLOCATE(ZLSZWSM_C(IDIMX_C,IDIMY_C))
  !$**** 3D
  ALLOCATE(ZUT_C(IDIMX_C,IDIMY_C,SIZE(PUT,3)))
  ALLOCATE(ZLSUM_C(IDIMX_C,IDIMY_C,SIZE(PUT,3)))
  ALLOCATE(ZVT_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZLSVM_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZWT_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZLSWM_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZLSTHM_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZLSRVM_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  !$20140709
  ALLOCATE(ZHUT_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZTKET_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZSRCT_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZSIGS_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZTHVT_C(IDIMX_C,IDIMY_C,SIZE(PUT,3)))
  ALLOCATE(ZVU_FLUX_M_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZVTH_FLUX_M_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  ALLOCATE(ZWTH_FLUX_M_C(IDIMX_C,IDIMY_C,SIZE(PVT,3)))
  !$***** 4D
  ALLOCATE(ZRT_C(IDIMX_C,IDIMY_C,SIZE(PUT,3),SIZE(PRT,4)))
  ALLOCATE(ZSVT_C(IDIMX_C,IDIMY_C,SIZE(PUT,3),NSV))
  ALLOCATE(ZDRVFRC_C(IDIMX_C,IDIMY_C,SIZE(PUT,3),SIZE(PRT,4)))
  ALLOCATE(ZDTHFRC_C(IDIMX_C,IDIMY_C,SIZE(PUT,3),SIZE(PRT,4)))
  ALLOCATE(ZRVREL_C(IDIMX_C,IDIMY_C,SIZE(PUT,3),SIZE(PRT,4)))
  ALLOCATE(ZTHREL_C(IDIMX_C,IDIMY_C,SIZE(PUT,3),SIZE(PRT,4)))
  !$initialize
  !$***** 3D
  ZUT_C   =0.
  ZLSUM_C =0.
  ZVT_C   =0.
  ZWT_C   =0.
  ZTHVT_C =0.
  ZZWS_C  =0.
  ZLSZWSM_C=0.
  ZHUT_C  =0.
  ZTKET_C =0.
  ZSRCT_C =0.
  ZSIGS_C =0.
  ZVU_FLUX_M_C=0.
  ZVTH_FLUX_M_C=0.
  ZWTH_FLUX_M_C=0.
  !$***** 4D
  ZRT_C   =0.
  ZSVT_C  =0.
  ZDRVFRC_C=0.
  ZDTHFRC_C=0.
  ZRVREL_C=0.
  ZTHREL_C=00
!
    CALL SET_LSFIELD_1WAY_ll(XZWS1(:,:),ZZWS_C(:,:),2)
    CALL SET_LSFIELD_1WAY_ll(XLSZWSM1(:,:),ZLSZWSM_C(:,:),2)
    !
    CALL LS_FORCING_ll(2, IINFO_ll, .TRUE.)
    CALL GO_TOMODEL_ll(2, IINFO_ll)
    CALL GOTO_MODEL(2)
    CALL UNSET_LSFIELD_1WAY_ll()
    !
  !$***** 3D VARS
  DO JI=1,SIZE(PUT,3)
    CALL GOTO_MODEL(1)
    CALL GO_TOMODEL_ll(1, IINFO_ll)
    !
    !$series of SET_LSFIELD_1WAY_ll
    !$***** 3D VARS
    CALL SET_LSFIELD_1WAY_ll(XUT1(:,:,JI),ZUT_C(:,:,JI),2)
    CALL SET_LSFIELD_1WAY_ll(XLSUM1(:,:,JI), ZLSUM_C(:,:,JI),2)
    !
    CALL SET_LSFIELD_1WAY_ll(XVT1(:,:,JI),ZVT_C(:,:,JI),2)
    CALL SET_LSFIELD_1WAY_ll(XLSVM1(:,:,JI),ZLSVM_C(:,:,JI),2)
    !
    CALL SET_LSFIELD_1WAY_ll(XWT1(:,:,JI),ZWT_C(:,:,JI),2)
    CALL SET_LSFIELD_1WAY_ll(XLSWM1(:,:,JI),ZLSWM_C(:,:,JI),2)
    !
    CALL SET_LSFIELD_1WAY_ll(ZTHVT(:,:,JI), ZTHVT_C(:,:,JI),2)
    CALL SET_LSFIELD_1WAY_ll(XLSTHM1(:,:,JI),ZLSTHM_C(:,:,JI),2)
    !$conditionnal VARS
    IF (HTURB /= 'NONE') THEN
      CALL SET_LSFIELD_1WAY_ll(XTKET1(:,:,JI), ZTKET_C(:,:,JI),2)
    ENDIF
    IF (CONF_MODEL(1)%NRR>=1) THEN
      CALL SET_LSFIELD_1WAY_ll(XLSRVM1(:,:,JI), ZLSRVM_C(:,:,JI),2)
      CALL SET_LSFIELD_1WAY_ll(ZHUT(:,:,JI),ZHUT_C(:,:,JI),2)
    ENDIF
    IF (CONF_MODEL(1)%NRR>1 .AND. HTURB /='NONE') THEN
      CALL SET_LSFIELD_1WAY_ll(XSRCT1(:,:,JI),ZSRCT_C(:,:,JI),2)
      CALL SET_LSFIELD_1WAY_ll(XSIGS1(:,:,JI),ZSIGS_C(:,:,JI),2)
    ENDIF
    IF (LUV_FLX)                                    &
      CALL SET_LSFIELD_1WAY_ll(XVU_FLUX_M1(:,:,JI),ZVU_FLUX_M_C(:,:,JI),2)
    IF (LTH_FLX) THEN
      CALL SET_LSFIELD_1WAY_ll(XVTH_FLUX_M1(:,:,JI),ZVTH_FLUX_M_C(:,:,JI),2)
      CALL SET_LSFIELD_1WAY_ll(XWTH_FLUX_M1(:,:,JI),ZWTH_FLUX_M_C(:,:,JI),2) 
    ENDIF
    !
    CALL LS_FORCING_ll(2, IINFO_ll, .TRUE.)
    CALL GO_TOMODEL_ll(2, IINFO_ll)
    CALL GOTO_MODEL(2)
    CALL UNSET_LSFIELD_1WAY_ll()
!
  ENDDO
!if the child grid is the whole father grid, we first need to extrapolate
!the data on a "pseudo halo" before doing BIKHARDT interpolation
! -------> done in LS_FORCING_ll
  !$***** 4D VARS
  DO JI=1,SIZE(PUT,3)
    DO KI=1,SIZE(PRT,4)
      CALL GOTO_MODEL(1)
      CALL GO_TOMODEL_ll(1, IINFO_ll)
      IF (CONF_MODEL(1)%NRR>=1) THEN
        CALL SET_LSFIELD_1WAY_ll(XRT1(:,:,JI,KI),ZRT_C(:,:,JI,KI),2)
      ENDIF
      IF ( L2D_ADV_FRC ) THEN
        CALL SET_LSFIELD_1WAY_ll(ADVFRC_MODEL(1)%XDTHFRC(:,:,JI,KI),ZDTHFRC_C(:,:,JI,KI),2)
        CALL SET_LSFIELD_1WAY_ll(ADVFRC_MODEL(1)%XDRVFRC(:,:,JI,KI),ZDRVFRC_C(:,:,JI,KI),2)
      ENDIF
      IF (L2D_REL_FRC) THEN
        CALL SET_LSFIELD_1WAY_ll(RELFRC_MODEL(1)%XTHREL(:,:,JI,KI),ZTHREL_C(:,:,JI,KI),2)
        CALL SET_LSFIELD_1WAY_ll(RELFRC_MODEL(1)%XRVREL(:,:,JI,KI),ZRVREL_C(:,:,JI,KI),2)
      ENDIF
      !
      CALL LS_FORCING_ll(2, IINFO_ll, .TRUE.)
      CALL GO_TOMODEL_ll(2, IINFO_ll)
      CALL GOTO_MODEL(2)
      CALL UNSET_LSFIELD_1WAY_ll()
!
    ENDDO
  ENDDO
  !$***** 4D NSV 
  IF (NSV>=1) THEN
     DO JI=1,SIZE(PUT,3)
        DO KI=1,NSV
           CALL GOTO_MODEL(1)
           CALL GO_TOMODEL_ll(1, IINFO_ll)           
           CALL SET_LSFIELD_1WAY_ll(FIELD_MODEL(1)%XSVT(:,:,JI,KI),ZSVT_C(:,:,JI,KI),2)           
           CALL LS_FORCING_ll(2, IINFO_ll, .TRUE.)
           CALL GO_TOMODEL_ll(2, IINFO_ll)
           CALL GOTO_MODEL(2)
           CALL UNSET_LSFIELD_1WAY_ll()
           !
        ENDDO
     ENDDO
  ENDIF

!if the child grid is the whole father grid, we first need to extrapolate
!the data on a "pseudo halo" before doing BIKHARDT interpolation
! -------> done in LS_FORCING_ll
!
!                        Interpolation of the U variable at t
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,2,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZUT_C,PUT)
    CALL MPPDB_CHECK3D(PUT,"SPAWN_FIELD2:PUT",PRECISION)
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,2,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZLSUM_C,PLSUM)
    CALL MPPDB_CHECK3D(PLSUM,"SPAWN_FIELD2:PLSUM",PRECISION)
!
!                        Interpolation of the V variable at t
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,3,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZVT_C,PVT)
    CALL MPPDB_CHECK3D(PVT,"SPAWN_FIELD2:PVT",PRECISION)
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,3,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZLSVM_C,PLSVM)
    CALL MPPDB_CHECK3D(PLSVM,"SPAWN_FIELD2:PLSVM",PRECISION)

!                        Interpolation of the ZWS variable at t
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,3,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZZWS_C,PZWS)
    CALL MPPDB_CHECK2D(PZWS,"SPAWN_FIELD2:PZWS",PRECISION)
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,3,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZLSZWSM_C,PLSZWSM)
    CALL MPPDB_CHECK2D(PLSZWSM,"SPAWN_FIELD2:PLSZWSM",PRECISION)
!
!
!                        Interpolation of variables at t
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,4,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZWT_C,PWT)
    CALL MPPDB_CHECK3D(PWT,"SPAWN_FIELD2:PWT",PRECISION)
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,4,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZLSWM_C,PLSWM)
    CALL MPPDB_CHECK3D(PLSWM,"SPAWN_FIELD2:PLSWM",PRECISION)
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZLSTHM_C,PLSTHM)
    CALL MPPDB_CHECK3D(PLSTHM,"SPAWN_FIELD2:PLSTHM",PRECISION)
!
!
    IF (CONF_MODEL(1)%NRR>=1) THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZLSRVM_C,PLSRVM)
      CALL MPPDB_CHECK3D(PLSRVM,"SPAWN_FIELD2:PLSRVM",PRECISION)               
    ENDIF
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZTHVT_C,PTHVT)
    CALL MPPDB_CHECK3D(PTHVT,"SPAWN_FIELD2:PTHVT",PRECISION)
!
    IF (HTURB /= 'NONE') THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZTKET_C,PTKET)
      CALL MPPDB_CHECK3D(PTKET,"SPAWN_FIELD2:PTKET",PRECISION)
    ENDIF
!
    IF (CONF_MODEL(1)%NRR>=1) THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZHUT_C,PHUT)
      CALL MPPDB_CHECK3D(PHUT,"SPAWN_FIELD2:PHUT",PRECISION)
    ENDIF
!
    IF (CONF_MODEL(1)%NRR>=1) THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZRT_C,PRT)
      CALL MPPDB_CHECK3D(PRT(:,:,:,1),"SPAWN_FIELD2:PRT",PRECISION)
    ENDIF
!
    IF (NSV>=1) THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZSVT_C,PSVT)
      CALL MPPDB_CHECK3D(PSVT(:,:,:,1),"SPAWN_FIELD2:PSVT",PRECISION)
    ENDIF
!
    IF (CONF_MODEL(1)%NRR>1 .AND. HTURB /='NONE') THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZSRCT_C,PSRCT)
      CALL MPPDB_CHECK3D(PSRCT,"SPAWN_FIELD2:PSRCT",PRECISION)
    ENDIF
!
    IF (CONF_MODEL(1)%NRR>1 .AND. HTURB /='NONE') THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZSIGS_C,PSIGS)
      CALL MPPDB_CHECK3D(PSIGS,"SPAWN_FIELD2:PSIGS",PRECISION)
    ENDIF
!
    IF ( L2D_ADV_FRC ) THEN      ! MT adding for ADVFRC
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZDTHFRC_C,PDTHFRC)
!
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZDRVFRC_C,PDRVFRC)
    ENDIF
    IF (L2D_REL_FRC) THEN      ! MT adding for REL FRC
       WRITE(ILUOUT,FMT=*) 'SPAWN_FIELD2: Appel a BIKHARDT pour RELFRC'
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,       &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZTHREL_C,PTHREL)
!
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,       &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZRVREL_C,PRVREL)
    ENDIF
!
    IF ( LUV_FLX) THEN      ! MT adding for EDDY_FLUX
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,       &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZVU_FLUX_M_C,PVU_FLUX_M)
      CALL MPPDB_CHECK3D(PVU_FLUX_M,"SPAWN_FIELD2:PVU_FLUX_M",PRECISION)
    ENDIF
!
    IF (LTH_FLX) THEN
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,       &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZVTH_FLUX_M_C,PVTH_FLUX_M)
      CALL MPPDB_CHECK3D(PVTH_FLUX_M,"SPAWN_FIELD2:PVTH_FLUX_M",PRECISION)
!
      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                     2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,       &
                     LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                     ZWTH_FLUX_M_C,PWTH_FLUX_M)
      CALL MPPDB_CHECK3D(PWTH_FLUX_M,"SPAWN_FIELD2:PWTH_FLUX_M",PRECISION)
    ENDIF
!
!!$END IF
!
IF (CONF_MODEL(1)%NRR>=3) THEN
  WHERE  (PRT(:,:,:,3)<1.E-20)
    PRT(:,:,:,3)=0.
  END WHERE
END IF
!
!
!*       2.2.3  Informations from model SON1
! (LS fields are not treated because they are identical in the father file)
!
IF (PRESENT(TPSONFILE)) THEN
  !
  !variables which always exist
  !
  CALL IO_READ_FIELD(TPSONFILE,'UT',ZWORK3D) ! U wind component at time t
  PUT(KIB2:KIE2,KJB2:KJE2,:) = ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  CALL IO_READ_FIELD(TPSONFILE,'VT',ZWORK3D) ! V wind component at time t
  PVT(KIB2:KIE2,KJB2:KJE2,:) = ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  CALL IO_READ_FIELD(TPSONFILE,'WT',ZWORK3D) ! W wind component at time t
  PWT(KIB2:KIE2,KJB2:KJE2,:) = ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  CALL IO_READ_FIELD(TPSONFILE,'ZWS',ZWORK2D,IRESP) !
  !If the field ZWS is not in the file, set its value to XZWS_DEFAULT
  !ZWS is present in files since MesoNH 5.4.2
  IF ( IRESP/=0 ) THEN
    WRITE (YVAL,'( E15.8 )') XZWS_DEFAULT
    CALL PRINT_MSG(NVERB_WARNING,'IO','SPAWN_FIELD2','ZWS not found in file: using default value: '//TRIM(YVAL)//' m')
    ZWORK2D(:,:) = XZWS_DEFAULT
  END IF
  PZWS(KIB2:KIE2,KJB2:KJE2) = ZWORK2D(KIB1:KIE1,KJB1:KJE1)
  !
  ! moist variables
  !
  IRR=1
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    GUSERV=.TRUE.
    CALL IO_READ_FIELD(TPSONFILE,'RVT',ZWORK3D,IRESP) ! Vapor at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    CALL IO_READ_FIELD(TPSONFILE,'RCT',ZWORK3D,IRESP) ! Cloud at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    CALL IO_READ_FIELD(TPSONFILE,'RRT',ZWORK3D,IRESP) ! Rain at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    CALL IO_READ_FIELD(TPSONFILE,'RIT',ZWORK3D,IRESP) ! Ice at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    CALL IO_READ_FIELD(TPSONFILE,'RST',ZWORK3D,IRESP) ! Snow at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    CALL IO_READ_FIELD(TPSONFILE,'RGT',ZWORK3D,IRESP) ! Graupel at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IF (IRR<=CONF_MODEL(1)%NRR) THEN
    CALL IO_READ_FIELD(TPSONFILE,'HVT',ZWORK3D,IRESP) ! Hail at time t
    IF(IRESP==0) ZRT1(:,:,:,IRR)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    IF(IRESP==0) IRR=IRR+1
  END IF
  IRR=IRR-1
  WRITE(ILUOUT,FMT=*) 'SPAWN_FIELD2: spawing with a SON input file'
  WRITE(ILUOUT,FMT=*) '    ',CONF_MODEL(1)%NRR,' moist variables in model1 and model2, ',    &
                             IRR,' moist variables in input SON'
  CALL IO_READ_FIELD(TPSONFILE,'THT',ZWORK3D) ! Theta at time t
  ZTHT1(:,:,:)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  CALL IO_READ_FIELD(TPSONFILE,'PABST',ZWORK3D) ! Pressure at time t
  ZPABST1(:,:,:)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  !
  CALL COMPUTE_THV_HU(GUSERV,ZRT1,ZTHT1,ZPABST1,ZTHVT1,ZHUT1)
  !
  PTHVT(KIB2:KIE2,KJB2:KJE2,:) = ZTHVT1(:,:,:)
  IF (CONF_MODEL(1)%NRR /= 0) THEN
    PHUT(KIB2:KIE2,KJB2:KJE2,:) = ZHUT1(:,:,:)  
    PRT(KIB2:KIE2,KJB2:KJE2,:,:) = ZRT1(:,:,:,:)  
  END IF
  !
  ! TKE variables
  !
  IF (HTURB/='NONE') THEN
    CALL IO_READ_FIELD(TPSONFILE,'TKET',ZWORK3D,IRESP) ! Turbulence Kinetic Energy at time t
    IF(IRESP==0) PTKET(KIB2:KIE2,KJB2:KJE2,:)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  END IF
  !
  ! Scalar variables
  !
  IF (NSV /= 0) THEN
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
      DO JSV = 1, NSV_USER      ! Users Scalar Variables
        WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'SVT',JSV
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! LIMA variables
    !
    DO JSV = NSV_LIMA_BEG,NSV_LIMA_END
      TZFIELD%CSTDNAME   = ''
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CUNITS     = 'kg-1'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
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
      ! time t
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
      IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    END DO
    !
    ! ELEC Scalar Variables
    !
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
          TZFIELD%CUNITS   = 'C m-3'
          WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        ELSE
          TZFIELD%CUNITS   = 'm-3'
          WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3,A)')'X_Y_Z_','SVT',JSV,' (nb ions/m3)'
        END IF
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Chemical Scalar Variables
    !
    IF (NSV_CHEMEND>=NSV_CHEMBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_CHEMBEG,NSV_CHEMEND
        TZFIELD%CMNHNAME   = TRIM(CNAMES(JSV-NSV_CHEMBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Ice phase chemical Scalar Variables
    !
    IF (NSV_CHICEND>=NSV_CHICBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_CHICBEG,NSV_CHICEND
        CICNAMES(JSV-NSV_CHICBEG+1) = UPCASE(CICNAMES(JSV-NSV_CHICBEG+1))
        TZFIELD%CMNHNAME   = TRIM(CICNAMES(JSV-NSV_CHICBEG+1))//'T'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Orilam Scalar Variables
    !
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Dust Scalar Variables
    !
    IF (NSV_DSTEND>=NSV_DSTBEG) THEN
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Sea Salt Scalar Variables
    !
    IF (NSV_SLTEND>=NSV_SLTBEG) THEN
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! LG Scalar Variables
    !
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! LNOx Scalar Variables
    !
!PW:TODO/bug1?: LINOX or LINOXT?
!PW:TODO/bug2?: Same name of variable in a loop!
    IF (NSV_LNOXEND>=NSV_LNOXBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'ppp' !PW: TODO: not sure (depends if LINOX or LINOXT)
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_LNOXBEG,NSV_LNOXEND
        TZFIELD%CMNHNAME   = 'LINOX'
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Passive scalar variables
    !
    IF (NSV_PPEND>=NSV_PPBEG) THEN
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
#ifdef MNH_FOREFIRE
    !
    ! ForeFire variables
    !
    IF (NSV_FFEND>=NSV_FFBEG) THEN
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
#endif
    !
    ! Passive scalar variables
    !
    IF (NSV_CSEND>=NSV_CSBEG) THEN
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
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PSVT(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
    !
    ! Passive scalar variables
    !
    IF (NSV_PP>=1) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'm-3'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = 1,NSV_PP
        WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'ATC',JSV+NSV_PPBEG-1
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','ATC',JSV+NSV_PPBEG-1
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PATC(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
#ifdef MNH_FOREFIRE
    !
    ! ForeFire variables
    !
    IF (NSV_FF>=1) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 'm-3'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = 1,NSV_FF
        WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'ATC',JSV+NSV_FFBEG-1
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','ATC',JSV+NSV_FFBEG-1
        CALL IO_READ_FIELD(TPSONFILE,TZFIELD,ZWORK3D,IRESP)
        IF(IRESP==0) PATC(KIB2:KIE2,KJB2:KJE2,:,JSV)=ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
      END DO
    END IF
#endif
  END IF
  !
  ! Secondary pronostic variables
  !
  IF (HTURB /= 'NONE' .AND. IRR>1) THEN
    CALL IO_READ_FIELD(TPSONFILE,'SRCT',ZWORK3D,IRESP) ! turbulent flux SRC at time t
    IF(IRESP == 0) PSRCT(KIB2:KIE2,KJB2:KJE2,:) =                    &
                                        ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
    CALL IO_READ_FIELD(TPSONFILE,'SIGS',ZWORK3D,IRESP) ! subgrid condensation
    IF(IRESP == 0) PSIGS(KIB2:KIE2,KJB2:KJE2,:) =                    &
                                        ZWORK3D(KIB1:KIE1,KJB1:KJE1,:)
  END IF
END IF
!
!*       2.2.4  secondary prognostic variables correction
!
IF (CONF_MODEL(1)%NRR > 1 .AND. HTURB /= 'NONE')  PSRCT(:,:,:) = MIN( 1.0, MAX( 0.0, PSRCT(:,:,:)) )
!
IF ( CONF_MODEL(1)%NRR == 0 ) THEN
  PHUT (:,:,:)= 0.
END IF
!-------------------------------------------------------------------------------
!
CALL GOTO_MODEL(IMI)
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CONTAINS 
!      
      SUBROUTINE COMPUTE_THV_HU(OUSERV,PR,PTH,PPABS,PTHV,PHU)
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, INTENT(IN)   :: OUSERV
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTH,PPABS
REAL, DIMENSION(:,:,:,:), INTENT(IN)   :: PR
REAL, DIMENSION(:,:,:),   INTENT(OUT)  :: PTHV,PHU
!
!*       0.2    Declarations of local variables 
!
REAL, DIMENSION(SIZE(PR,1),SIZE(PR,2),SIZE(PR,3)) :: ZSUMR ! sum of water ratios
!
IF (OUSERV) THEN
  ZSUMR(:,:,:) = 0.
  IRR=SIZE(PR,4)
  DO JRR=1,IRR
    ZSUMR(:,:,:) = ZSUMR(:,:,:) + PR(:,:,:,JRR)
  END DO
  PTHV(:,:,:)=PTH(:,:,:)*(1.+XRV/XRD*PR(:,:,:,1))/(1.+ZSUMR(:,:,:))
  PHU (:,:,:)=100.*PPABS(:,:,:)/(XRD/XRV/MAX(PR(:,:,:,1),1.E-16)+1.) &
               /SM_FOES(PTH(:,:,:)*(PPABS(:,:,:)/XP00)**(XRD/XCPD))
ELSE
  PTHV(:,:,:)=PTH(:,:,:)
  PHU (:,:,:)=0.
END IF
!
!
END SUBROUTINE COMPUTE_THV_HU
!
END SUBROUTINE SPAWN_FIELD2
