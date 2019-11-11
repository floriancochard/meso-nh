!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!!    ########################
      MODULE MODI_CH_MONITOR_n
!!    ########################
!!
!
INTERFACE
!!
SUBROUTINE CH_MONITOR_n(PWETDEPAER, KTCOUNT,PTSTEP, KLUOUT, KVERB)
IMPLICIT NONE
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PWETDEPAER ! tendency of aerosol wet depostion
INTEGER, INTENT(IN) :: KTCOUNT    ! iteration count
REAL,  INTENT(IN)   :: PTSTEP    ! Double timestep except 
                                  ! for the first time step (single one)
INTEGER, INTENT(IN) :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN) :: KVERB      ! verbosity level
END SUBROUTINE CH_MONITOR_n
!!
END INTERFACE
!!
END MODULE MODI_CH_MONITOR_n
!!
!!    ####################################################### 
      SUBROUTINE CH_MONITOR_n(PWETDEPAER, KTCOUNT,PTSTEP, KLUOUT, KVERB)
!!    #######################################################
!!
!!*** *CH_MONITOR_n*  monitor of the chemical module
!!
!!    PURPOSE
!!    -------
!!       The purpose of this subroutine is to control the chemical module
!!    i.e. to pass the meteorological parameters from MesoNH to its chemical
!!    part and to call the different subroutines (calculation of rate constants,
!!    photolysis rates, stiff solver,..)
!!
!!    METHOD
!!    ------
!!       The calculation  of the chemical terms is performed using a loop
!!    over all spatial dimensions. 
!!
!!       For each single grid point, all necessary meteorological parameters are
!!    passed into the chemical core system (variable TZM). This variable is
!!    then passed on to the subroutines that calculate the reaction and
!!    photolysis rates. Then the chemical solver is called. As the chemistry
!!    part works with different units than MesoNH (MesoNH uses mixing ratio,
!!    the chemisty part uses molec/cm3) some unit conversion is also performed.
!!
!!       Temporal integration is performed over a double timestep 2*XTSTEP
!!    (except in the case of a cold start). If the timestep of MesoNH
!!    is too large for the chemical solver, several smaller steps can
!!    be taken using the NCH_SUBSTEPS parameter.
!!    Three options of temporal discretization are implemented:
!!    "SPLIT"  : from XRSVS the scalar variable at t+dt is calculated and
!!               given as input to the solver; the result is rewritten 
!!               into XRSVS; this corresponds to applying first only dynamics
!!               and then only chemistry; this option assures positivity, but
!!               degrades the order of the temporal integration.
!!               In fact, an overhead of a factor two is produced here.
!!               A future solution will be to calculate the dynamics
!!               of the scalar variables not using leapfrog, but forward
!!               temporal integration.
!!    "CENTER" : the scalar variables at t (XSVT) are taken in order to 
!!               calculate the tendencies for chemistry, that are then applied
!!               together with all other terms in parallel; this option
!!               is consistent with the MesoNH leapfrog scheme, but
!!               unfortunately it tends to be unstable due to the stiffness
!!               of the chemical system; thus this option is not recommended.
!!    "LAGGED" : the scalar variables at t-dt (XSVM) are taken in order to
!!               calculate the tendencies for chemistry, that are then applied
!!               together with all other terms in parallel; this option
!!               does not garantee positivity, but seems to be stable.
!!    The options "CENTER" and "LAGGED" are implemented more for test than
!!    for production purposes.
!!
!!    REFERENCE
!!    ---------
!!    Book 1, 2, 3 of MesoNH-chemistry
!!
!!    AUTHOR
!!    ------
!!    K. Suhre    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/05/95
!!    26/10/95 KS: add conversion mixing ratio -> concentration
!!                 and use first guess variable as input (split)
!!    27/10/95 KS: change parameterlist
!!    04/08/96 (K. Suhre) restructered in order to run with grid-nesting
!!    09/03/99 (V. Crassier & K. Suhre) vectorization
!!    09/03/99 (K. Suhre) TUV online
!!    06/06/00 (C. Mari) add 1-D timeseries for chemistry
!!    21/03/01 (C. Mari & J. Escobar) Code optimization
!!    01/08/01 (C. Mari)  change CH_SOLVER to $n 
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    01/01/04 (P. Tulet)   Bugs in conversion for CENTER and LAGGED options
!!    01/01/06 (P. Tulet)   ORILAM aerosol scheme
!!    04/06/07 (M. Leriche) add pH
!!    30/07/07 (JP Pinty) add Rosenbrock solver
!!    26/03/08 (M Leriche) add microphysical transfert from collision/coalescence
!!    10/11/08 (M Leriche) add microphysical transfert from rain sedimentation
!!    24/04/14 (M Leriche) Bugs in orilam transfert zsvt in xrsvs
!!                         + supress line transfer H2SO4 from AP to gas phase
!!                         imply transfer H2SO4 AP in aqueous phase if aq.chem.
!!    04/2014 (C.Lac) Remove GCENTER with FIT temporal scheme
!!    06/11/14 (M Leriche) Bug in pH computing
!!    11/12/15 (M. Leriche & P. Tulet) add ch_init_ice initialise index for ice chem.
!!    18/01/16 (M Leriche) for sedimentation fusion C2R2 and khko
!!    15/02/16 (M Leriche) call ch_init_rosenbrock only one time
!!    20/01/17 (G.Delautier) bug if CPROGRAM/=DIAG
!!    01/10/17 (C.Lac) add correction of negativity
!!    Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 12/02/2019: bugfix: ZINPRR was not initialized all the time
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_METEO_TRANS_KESS
USE MODI_CH_METEO_TRANS_C2R2
USE MODI_CH_SET_RATES
USE MODI_CH_SET_PHOTO_RATES
USE MODI_CH_SOLVER_n
USE MODI_CH_UPDATE_JVALUES
USE MODI_BUDGET
USE MODI_CH_INIT_ICE
USE MODI_CH_AQUEOUS_TMICICE
USE MODI_CH_AQUEOUS_TMICKESS
USE MODI_CH_AQUEOUS_TMICC2R2
USE MODI_CH_AQUEOUS_TMICKHKO
USE MODI_CH_AQUEOUS_SEDIM1MOM
USE MODI_CH_AQUEOUS_SEDIM2MOM
USE MODI_CH_AQUEOUS_CHECK
USE MODI_SUM_ll
USE MODI_CH_AER_SEDIM_n
USE MODI_CH_AER_WETDEP_n
USE MODI_CH_ORILAM
USE MODI_CH_INI_ORILAM
USE MODI_CH_AER_EQM_CORMASS
USE MODI_CH_AER_SURF
USE MODI_CH_AER_DEPOS
!
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODI_WRITE_TS1D
USE MODD_CST, ONLY : XMNH_TINY
!
USE MODI_CH_PRODLOSS
!     IMPLICIT ARGUMENTS
!     ------------------
! 
USE MODD_BUDGET
USE MODD_LUNIT_n
USE MODD_NSV, ONLY : NSV_CHEMBEG,NSV_CHEMEND,NSV_CHEM,& ! index for chemical SV
                     NSV_CHACBEG,NSV_CHACEND,NSV_CHAC,& ! index for aqueous SV
                     NSV_CHGSBEG,NSV_CHGSEND,         & ! index for gas phase SV
                     NSV_CHICBEG,NSV_CHICEND,         & ! index for ice phase SV
                     NSV_C2R2BEG,                     & ! index for number concentration
                     NSV_AERBEG, NSV_AEREND, NSV_AER, & ! index for aerosols SV
                     XSVMIN
!
USE MODD_CH_M9_n,   ONLY: NEQ,            &! number of prognostic chem. species
                          NEQAQ,          &! number of aqueous chem. species
                          NMETEOVARS,     &! number of meteorological variables
                          CNAMES,         &! names of the chem. species
                          CICNAMES,       &! names of the ice chem. species
                          METEOTRANSTYPE, &! type for meteo . transfer
                          NREAC,          &
                          NNONZEROTERMS,  &
                          CREACS          
!
USE MODI_CH_TERMS
USE MODI_CH_NONZEROTERMS
USE MODI_CH_GET_RATES
!
USE MODD_CH_MNHC_n, ONLY: CCH_TDISCRETIZATION
                  ! temporal discretization:
                  ! "SPLIT"  : use time-splitting, input fields for solver are
                  !            scalar variables at t+dt (derived from XRSVS)
                  ! "CENTER" : input fields for solver are
                  !            scalar variables at t (XSVT)
                  ! "LAGGED" : input fields for solver are
                  !            scalar variables at t-dt (XSVM)
USE MODD_CH_MNHC_n, ONLY: NCH_SUBSTEPS
                  ! number of chemical timesteps to be taken during one
                  ! double timestep of MesoNH (MesoNH integrates with timesteps
                  ! of lenght 2*XTSTEP using leapfrog), the timestep of the
                  ! solver will be calculated as
                  ! ZDTSOLVER = 2*XTSTEP/NCH_SUBSTEPS
USE MODD_CH_MNHC_n, ONLY: LCH_TUV_ONLINE, CCH_TUV_LOOKUP, CCH_TUV_CLOUDS,  &
                          XCH_TUV_ALBNEW, XCH_TUV_DOBNEW, XCH_TUV_TUPDATE, &
                          CCH_VEC_METHOD, NCH_VEC_LENGTH
                          ! used for vectorization and photolysis rates
USE MODD_CH_MNHC_n, ONLY: LUSECHAQ, LUSECHIC, LCH_PH, LCH_RET_ICE, XRTMIN_AQ
                  ! aqueous chemistry and pH
USE MODD_CH_SOLVER_n
!
USE MODD_CH_PH_n                      ! pH value in 3D
!
USE MODD_FIELD_n,   ONLY: XSVT,      &! scalar variable at t
                          XRSVS,     &! source of scalar variable
                          XRT,       &! water mixing ratio at t
                          XCIT,      &! pristine conc. at t
                          XRRS,      &! source of water mixing ratio 
                          XPABST,    &! pressure
                          XTHT        ! potential temperature
!
USE MODD_REF_n,     ONLY: XRHODREF,  &! dry density for ref. state
                          XRHODJ      ! ( rhod J ) = dry density
!
USE MODD_TIME,      ONLY: TDTEXP 
!
USE MODD_TIME_n,    ONLY: TDTCUR      ! Current Time and Date
!
USE MODD_CONF,      ONLY: CPROGRAM, L1D
USE MODD_PARAM_n,   ONLY: CCLOUD
!
USE MODD_PARAMETERS,ONLY: JPHEXT,    &! number of horizontal External points
                          JPVEXT      ! number of vertical External points
!
USE MODD_CST,       ONLY: XAVOGADRO, &! Avogadro number
                          XMD,       &! Molar mass of dry air
                          XP00, XRD, XCPD
!
USE MODD_CH_PRODLOSSTOT_n             ! Total production/loss for chemical
                                      ! species
USE MODD_CH_BUDGET_n                  ! Extended production/loss terms for
                                      ! chemical species             
!
USE MODD_DIAG_FLAG,     ONLY: CSPEC_BU_DIAG,CSPEC_DIAG
! variables used by TUV
!
USE MODD_GRID_n,    ONLY: XZZ,&       ! height z
                          XZS,&       ! orography
                          XLAT, XLON
USE MODD_GRID,      ONLY: XLAT0,XLON0 ! Reference longitude and latitude
USE MODD_CONF_n,    ONLY: LUSERV,&    ! Logical to use wapor water
                          LUSERC,&    ! Logical to use cloud water
                          LUSERR,&    ! Logical to use rain water
                          NRR,   &    ! Total number of water variables
                          NRRL        ! Number of liquid water variables
USE MODD_SUB_CH_MONITOR_n
USE MODD_DYN_n,     ONLY: XTSTEP      ! time step of MesoNH
!
! variables used by ORILAM
!
USE MODD_PRECIP_n, ONLY: XEVAP3D
USE MODD_CLOUDPAR_n, ONLY: NSPLITR  ! Nb of required small time step integration
!
!variables used by microphysical mass transfer - sedimentation
!
USE MODD_CLOUDPAR_n, ONLY: NSPLITR
!
!variables used by rosenbrock solver
!
USE MODD_CH_ROSENBROCK_n, ONLY: NSPARSEDIM,   & ! Dim of NSPARSE_xxx vectors
                              NSPARSE_IROW, & ! row index
                              NSPARSE_ICOL, & ! col index
                              NSPARSE_CROW, & ! first row element index
                              NSPARSE_DIAG, & ! diag index
                              NEQ_NAQ,          & ! number of Non-AQueous species
                              NSPARSEDIM_NAQ,   & ! Dim of NSPARSE_xxx vectors
                              NSPARSE_IROW_NAQ, & ! row index
                              NSPARSE_ICOL_NAQ, & ! col index
                              NSPARSE_CROW_NAQ, & ! first row element index
                              NSPARSE_DIAG_NAQ    ! diag index
                                      ! of the gridpoint sparse JACobian matrix
!
USE MODD_RBK90_JacobianSP_n         ! vectorized form of the sparse indexes
USE MODD_RBK90_Parameters_n, ONLY: NVAR, LU_NONZERO
!
! parameters of the namelist to come
!
USE MODD_VAR_ll
USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n

USE MODD_CH_INIT_JVALUES, ONLY: JPJVMAX ! number of photolysis reactions in TUV
USE MODD_CH_JVALUES_n,    ONLY: XJVALUES
USE MODD_CH_MNHC_n,       ONLY: CCH_SCHEME , LCH_CONV_SCAV
USE MODD_RADIATIONS_n,    ONLY: XZENITH, XALBUV
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PWETDEPAER ! tendency of aerosol wet depostion
INTEGER, INTENT(IN) :: KTCOUNT    ! iteration count
REAL,  INTENT(IN)   :: PTSTEP    ! Double timestep except 
                                  ! for the first time step (single one)
INTEGER, INTENT(IN) :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN) :: KVERB      ! verbosity level
!
!*      0.2    declarations of local variables
!
INTEGER :: JI,JJ,JK,JL,JM,JN   ! loop counters
REAL    :: ZDTSOLVER        ! timestep for the solver
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCHEM, ZOLDCHEM,  ZNEWCHEM 
REAL, DIMENSION(:,:), ALLOCATABLE :: ZAERO, ZOLDAERO,  ZNEWAERO
        ! arrays for parameter passage to solver
!
REAL, DIMENSION(:), ALLOCATABLE    :: ZCONV            
        ! conversion factor mixing ratio * RhoDJ ! to molec./cm3
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPH
!
!Varibales for integrated prod/loss for given species
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPRODTOT  ! Production/loss tables
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLOSSTOT  ! for all species
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPROD     ! Production/loss tables
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLOSS     ! for selected species
!
!Variables for detailed production/destruction terms for given species
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTCHEMTOT    ! detailed production/loss terms
INTEGER, DIMENSION(:,:),ALLOCATABLE   :: IINDEX   ! indices of non-zero terms
INTEGER                               :: IREAC    ! indices of reaction
INTEGER, DIMENSION(:),ALLOCATABLE   :: IIND   
TYPE REAC
  INTEGER, DIMENSION(:), POINTER  :: IB_REAC
  REAL   , DIMENSION(:,:), POINTER  :: ZB_REAC
END TYPE
TYPE (REAC), ALLOCATABLE, DIMENSION(:)  :: ZTCHEM
!     
INTEGER                                 :: JO
INTEGER                                 :: JR
INTEGER                                 :: JS
!
!
REAL    :: ZDEN2MOL
        !  ZDEN2MOL = 6.0221367E+23 * 1E-6 / 28.9644E-3
        !  conversion factor density to mol/cm3
        !  n_molec (moelc./cm3):  M = 1E-6*RHO(kg/m3) * XAVOGADRO / XMD
!
TYPE(METEOTRANSTYPE), DIMENSION(:), ALLOCATABLE  :: TZM 
        ! meteo variables to be transferred into CCS
!
!
LOGICAL :: GSPLIT           ! use timesplitting as temporal discretization
!
INTEGER             :: IIU  ! Upper dimension in x direction
INTEGER             :: IJU  ! Upper dimension in y direction
INTEGER             :: IKU  ! Upper dimension in z direction
INTEGER             :: IIB  ! indice I Beginning in x direction
INTEGER             :: IJB  ! indice J Beginning in y direction
INTEGER             :: IKB  ! indice K Beginning in z direction
INTEGER             :: IIE  ! indice I End       in x direction
INTEGER             :: IJE  ! indice J End       in y direction
INTEGER             :: IKE  ! indice K End       in z direction
!
!---------------------------------------------------------------------------
!   variables for the vectorization
!
INTEGER :: ITOT,IMAX,IDUM
LOGICAL :: GEND,GENDTOT
!
INTEGER :: IDTI,IDTJ,IDTK
INTEGER :: IDT1,IDT2,IDT3
INTEGER :: IDUMI,IDUMJ,IDUMK
INTEGER :: ITOTI,ITOTJ,ITOTK
!
!-------------------------------------------------------------------------------
!   variables for TUV
!
REAL                   :: ZRATIO, ZMASSTOT, ZMASSPOS
INTEGER                :: IINFO_ll       ! return code of parallel routine
INTEGER                :: JSV            ! loop index for SV
INTEGER                :: IMI            ! model index
!
!-------------------------------------------------------------------------------
!   variables for the aerosol module
!
REAL                   :: ZTIME                ! current time 
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZM, ZSIG0, ZN0, ZRG0, &   ! work array
                                       ZCTOTG, ZSEDA, ZFRAC, ZMI ! for aerosols
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZCTOTA, ZCCTOT
                                         ! first dimension is vectorization,
                                         ! second dim. are the modes*moments
REAL, DIMENSION(:),   ALLOCATABLE :: ZRV, ZDENAIR, ZPRESSURE, ZTEMP, ZRC
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRHOP0, ZOM, ZSOLORG
REAL, DIMENSION(:),   ALLOCATABLE :: ZLAMBDA, ZMU, ZSO4RAT

REAL,DIMENSION(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),SIZE(XSVT,4)) :: ZSVT
REAL,DIMENSION(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NSV_AER) :: ZCWETAERO
!
!-------------------------------------------------------------------------------
!   variables for AQueous/NAQueous cases
!
INTEGER                :: JRR          ! Loop index for the moist variables
REAL,DIMENSION(SIZE(XRT,1),SIZE(XRT,2),SIZE(XRT,3),SIZE(XRT,4))     :: ZRT_VOL
                                       ! liquid content in vol/vol
REAL, DIMENSION(SIZE(XRT,1), SIZE(XRT,2))     :: ZINPRR! Rain instant precip
!
!-------------------------------------------------------------------------------
!
! get model index
  IMI = GET_CURRENT_MODEL_INDEX()
!
!*       1.    PREPARE MONITOR
!              ---------------
!
!*       1.1   compute dimensions of arrays
!
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKU = SIZE(XRSVS,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
! Correction of negativity
!
DO JSV = 1, SIZE(XSVT,4)
 XRSVS(:,:,:,JSV)   = MAX((XRSVS(:,:,:,JSV)),XSVMIN(JSV))
END DO
!
!
IF (KTCOUNT == 1) THEN
!
!        1.1.1  determine mask to use for vectorisation at first step
!
  IDTI=(IIE-IIB+1)
  IDTJ=(IJE-IJB+1)
  IDTK=(IKE-IKB+1)
  ITOT=IDTI*IDTJ*IDTK
!
! the mask option will become a namelist parameter
!
  SELECT CASE (CCH_VEC_METHOD)
!
!***************************************************
! No mask (local)
!***************************************************
    CASE('LOC')
!
      ISVECNPT=1
      ISVECNMASK=ITOT
      IDT1=1
      IDT2=1
      IDT3=1
!
!***************************************************
! Horizontal mask
!***************************************************
    CASE('HOR')
!
      ISVECNPT=IDTI*IDTJ
      ISVECNMASK=IDTK
      IDT1=IDTI
      IDT2=IDTJ
      IDT3=1
!
!***************************************************
! Vertical mask
!***************************************************
    CASE('VER')
!
      ISVECNPT=IDTJ*IDTK
      ISVECNMASK=IDTI
      IDT1=1
      IDT2=IDTJ
      IDT3=IDTK
!
!***************************************************
! 1 mask with all grid points
! (no parallelisation)
!***************************************************
    CASE('TOT')
!
      ISVECNPT=IDTI*IDTJ*IDTK
      ISVECNMASK=1
      IDT1=IDTI
      IDT2=IDTJ
      IDT3=IDTK
!
!****************************************************
! Choice of a maximum number of points
!****************************************************
    CASE('MAX')
!
      GEND=.FALSE.
      GENDTOT=.FALSE.
!
      IMAX=MIN(NCH_VEC_LENGTH,ITOT)
!
      DO WHILE (.NOT.(GENDTOT))
!
        IDUM=IMAX
        DO WHILE (.NOT.(GEND))
          IF ((ITOT-IDUM*(ITOT/IDUM)) == 0) THEN
            GEND=.TRUE.
            ISVECNMASK=ITOT/IDUM
            ISVECNPT=IDUM
          ELSE
            IDUM=IDUM-1
          END IF
        END DO
!
        GEND=.FALSE.
        ITOTI=ISVECNPT
        IDUMI=IDTI
        DO WHILE (.NOT.(GEND) .AND. IDUMI >= 1)
          IF (     (ITOTI-IDUMI*(ITOTI/IDUMI)) == 0 &
             .AND. (IDTI-IDUMI*(IDTI/IDUMI)) == 0) THEN
            IDT1=IDUMI
            ITOTJ=ITOTI/IDUMI
            IDUMJ=IDTJ
            DO WHILE (.NOT.(GEND) .AND. IDUMJ >= 1)
              IF (     (ITOTJ-IDUMJ*(ITOTJ/IDUMJ)) == 0 &
                 .AND. (IDTJ-IDUMJ*(IDTJ/IDUMJ)) == 0) THEN
                IDT2=IDUMJ
                ITOTK=ITOTJ/IDUMJ
                IDUMK=IDTK
                DO WHILE (.NOT.(GEND) .AND. IDUMK >= 1)
                  IF (     (ITOTK-IDUMK*(ITOTK/IDUMK)) == 0 &
                     .AND. (IDTK-IDUMK*(IDTK/IDUMK)) == 0) THEN
                    IDT3=IDUMK
                    GEND=.TRUE.
                  ELSE
                    IDUMK=IDUMK-1
                  END IF
                END DO
              ELSE
                IDUMJ=IDUMJ-1
              END IF
            END DO      
          ELSE
            IDUMI=IDUMI-1
          END IF
        END DO
!
        GENDTOT=GEND
!
      END DO
!
  END SELECT
!
  ALLOCATE (ISVECMASK(6,ISVECNMASK))
!
!**********************************
! Compute mask boundaries
!**********************************
!
  ISVECMASK(1,1)=IIB
  ISVECMASK(2,1)=IIB+IDT1-1
  ISVECMASK(3,1)=IJB
  ISVECMASK(4,1)=IJB+IDT2-1
  ISVECMASK(5,1)=IKB
  ISVECMASK(6,1)=IKB+IDT3-1
!
  IF (ISVECNMASK .GE. 2) THEN
    DO JI=2,ISVECNMASK
      ISVECMASK(1,JI)=ISVECMASK(1,JI-1)+IDT1-IIB
      ISVECMASK(3,JI)=ISVECMASK(3,JI-1)
      ISVECMASK(5,JI)=ISVECMASK(5,JI-1)
!
      ISVECMASK(3,JI)=ISVECMASK(3,JI)+(ISVECMASK(1,JI)/IDTI)*IDT2-IJB
      ISVECMASK(5,JI)=ISVECMASK(5,JI)+(ISVECMASK(3,JI)/IDTJ)*IDT3-IKB
!
      ISVECMASK(1,JI)=ISVECMASK(1,JI)-IDTI*(ISVECMASK(1,JI)/IDTI)+IIB
      ISVECMASK(3,JI)=ISVECMASK(3,JI)-IDTJ*(ISVECMASK(3,JI)/IDTJ)+IJB
      ISVECMASK(5,JI)=ISVECMASK(5,JI)-IDTK*(ISVECMASK(5,JI)/IDTK)+IKB
!
      ISVECMASK(2,JI)=ISVECMASK(1,JI)+IDT1-1
      ISVECMASK(4,JI)=ISVECMASK(3,JI)+IDT2-1
      ISVECMASK(6,JI)=ISVECMASK(5,JI)+IDT3-1
    END DO
  END IF
!
!        1.1.2  determine sparse indexes to describe the jacobian matrix 
!               with vectorisation in a Rosenbrock solver without aqueous
!               chemistry
!
  IF (CSOLVER(1:2)=="RO" .AND. NEQAQ==0) THEN ! only for gaseous chemistry rosenbrock solver
    CALL PREPARE_LU_ROSENBROCK
  END IF
!
  ALLOCATE(LU_DIM_SPECIES(ISVECNPT))
  LU_DIM_SPECIES(:) = NEQ
!
!        1.1.3 determine index for ice phase chemistry or degassing with ICE3/4
  IF ((LUSECHAQ).AND.((CCLOUD=='ICE3' .OR. CCLOUD=='ICE4'))) THEN 
     CALL CH_INIT_ICE(LUSECHIC,LCH_RET_ICE,CNAMES,CICNAMES,NEQ,NEQAQ) 
  ENDIF
!
ENDIF  ! first time step
!
!*       1.2   calculate timestep variables
!
ZDTSOLVER = PTSTEP / NCH_SUBSTEPS
!
!*       1.3   give minimum value and conserve mass for aerosols
!
!
IF (LORILAM) THEN

  DO JSV = 1, SIZE(XSVT,4)
    ZSVT(:,:,:,JSV) =  XRSVS(:,:,:,JSV) *PTSTEP / XRHODJ(:,:,:) 
  END DO
  ZSVT(:,:,:,NSV_CHEMBEG:NSV_CHEMEND) = MAX(ZSVT(:,:,:,NSV_CHEMBEG:NSV_CHEMEND), XMNH_TINY)
  ZSVT(:,:,:,NSV_AERBEG:NSV_AEREND)   = MAX(ZSVT(:,:,:,NSV_AERBEG:NSV_AEREND), XMNH_TINY)
!
END IF
!
!*       1.4 compute conversion factor ppp/m3 --> molec/cm3
!
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
!
!*       1.5   set logical variables for temporal discretization
!
SELECT CASE (CCH_TDISCRETIZATION)
  CASE ("SPLIT")
    GSPLIT = .TRUE. 
    IF (KVERB >= 10) WRITE(KLUOUT,*) "CH_MONITOR_n: using SPLIT option"
  CASE ("CENTER")
    GSPLIT = .FALSE. 
    IF (KVERB >= 10) WRITE(KLUOUT,*) "CH_MONITOR_n: using CENTER option"
  CASE ("LAGGED")
    GSPLIT = .FALSE. 
    IF (KVERB >= 10) WRITE(KLUOUT,*) "CH_MONITOR_n: using LAGGED option"
  CASE DEFAULT
    ! the following line should never be reached:
    ! callabortstop
    CALL ABORT
    STOP "CH_MONITOR_n: CCH_TDISCRETIZATION option not valid"
END SELECT
!
!
IF (CPROGRAM=='DIAG  ') THEN
  IF (LEN_TRIM(CSPEC_BU_DIAG)/=0.OR.LEN_TRIM(CSPEC_DIAG)/=0) GSPLIT=.FALSE.  ! Modif. for DIAG
END IF
!
!
!*       1.6   allocate tables
!
ALLOCATE(TZM(ISVECNPT))
ALLOCATE(ZCHEM(ISVECNPT,NEQ))    !dimension of the 2nd row NEQ is provisional
ALLOCATE(ZNEWCHEM(ISVECNPT,NEQ)) !dimension of the 2nd row NEQ is provisional
ALLOCATE(ZOLDCHEM(ISVECNPT,NEQ)) !dimension of the 2nd row NEQ is provisional
ALLOCATE(ZCONV(ISVECNPT))
IF (LUSECHAQ.AND.LCH_PH) ALLOCATE(ZPH(ISVECNPT,NRRL))
IF (NEQ_PLT>0) THEN
  ALLOCATE(ZPRODTOT(ISVECNPT,NEQ))
  ALLOCATE(ZLOSSTOT(ISVECNPT,NEQ))
  ALLOCATE(ZPROD(ISVECNPT,NEQ_PLT))
  ALLOCATE(ZLOSS(ISVECNPT,NEQ_PLT))
END IF
IF (NEQ_BUDGET>0) THEN
  ALLOCATE(ZTCHEMTOT(ISVECNPT,NEQ,NREAC))
  ALLOCATE(ZTCHEM(NEQ_BUDGET))
  ALLOCATE(IIND(NEQ_BUDGET))
  ALLOCATE(IINDEX(2,NNONZEROTERMS))
  CALL CH_NONZEROTERMS(IMI,IINDEX,NNONZEROTERMS)
  DO JM=1,NEQ_BUDGET
      IIND(JM)=COUNT((IINDEX(1,:))==NSPEC_BUDGET(JM))
      ALLOCATE(ZTCHEM(JM)%IB_REAC(IIND(JM)))
      ALLOCATE(ZTCHEM(JM)%ZB_REAC(ISVECNPT,IIND(JM)))
  END DO
END IF
IF (LORILAM) THEN
  ALLOCATE(ZAERO(ISVECNPT,NSV_AER))
  ALLOCATE(ZNEWAERO(ISVECNPT,NSV_AER))
  ALLOCATE(ZOLDAERO(ISVECNPT,NSV_AER))
  ALLOCATE(ZM(ISVECNPT,JPIN))
  ALLOCATE(ZSEDA(ISVECNPT,JPIN))
  ALLOCATE(ZRHOP0(ISVECNPT,JPMODE))
  ALLOCATE(ZSIG0(ISVECNPT,JPMODE))
  ALLOCATE(ZRG0(ISVECNPT,JPMODE))
  ALLOCATE(ZN0(ISVECNPT,JPMODE))
  ALLOCATE(ZCTOTA(ISVECNPT,NSP+NCARB+NSOA,JPMODE))
  ALLOCATE(ZCCTOT(ISVECNPT,NSP+NCARB+NSOA,JPMODE))
  ALLOCATE(ZCTOTG(ISVECNPT,NSP+NCARB+NSOA))
  ALLOCATE(ZMU(ISVECNPT))
  ALLOCATE(ZLAMBDA(ISVECNPT))
  ALLOCATE(ZOM(ISVECNPT,JPMODE))
  ALLOCATE(ZSO4RAT(ISVECNPT))
  ALLOCATE(ZRV(ISVECNPT))
  ALLOCATE(ZRC(ISVECNPT))
  ALLOCATE(ZPRESSURE(ISVECNPT))
  ALLOCATE(ZTEMP(ISVECNPT))
  ALLOCATE(ZDENAIR(ISVECNPT))
  ALLOCATE(ZFRAC(ISVECNPT,NEQ))
  ALLOCATE(ZMI(ISVECNPT,NSP+NCARB+NSOA))
  ALLOCATE(ZSOLORG(ISVECNPT,NSOA))
  ALLOCATE(XSURF(ISVECNPT,JPMODE))
  ALLOCATE(XDP(ISVECNPT,JPMODE))
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    UPDATE PHOTOLYSIS RATES
!              -----------------------
!
IF (KTCOUNT==1 .OR. &
    (MOD(ISTCOUNT, MAX(1, INT(XCH_TUV_TUPDATE/XTSTEP)) ) .EQ. 0)) THEN
!
  WRITE(KLUOUT,*)"TIME call update jvalue: ",TDTCUR%TIME
!
  IF (.NOT.ASSOCIATED(XJVALUES)) &
             ALLOCATE(XJVALUES(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPJVMAX))
  XJVALUES(:,:,:,:) = 0.
  CALL CH_UPDATE_JVALUES(KLUOUT,  XZENITH, XRT,                  &
       XALBUV, XZS, XZZ, XLAT0, XLON0,                           &
       SIZE(XZZ,1), SIZE(XZZ,2), SIZE(XZZ,3), NRR,               &
       TDTCUR%TDATE%DAY, TDTCUR%TDATE%MONTH, TDTCUR%TDATE%YEAR, TDTCUR%TIME,&
       LCH_TUV_ONLINE,  CCH_TUV_CLOUDS,                          &
       XCH_TUV_ALBNEW, XCH_TUV_DOBNEW, XRHODREF, XJVALUES,       &
       IIB,IIE,IJB,IJE,IIU,IJU, KVERB   )
ENDIF
!
ISTCOUNT = ISTCOUNT + 1
!
!-------------------------------------------------------------------------------
!
!*       3.    MICROPHYSICS TERM FOR AEROSOL AND AQUEOUS CHEMISTRY
!              ---------------------------------------------------
!
!*       3.1 sedimentation term and wet deposition for aerosols tendency (XSEDA)
!
IF (LORILAM) THEN
  ZTIME  = TDTCUR%TIME ! need for ch_orilam
  XSEDA(:,:,:,:) = 0.
  ZSEDA(:,:) = 0.
! dry sedimentation
  IF ((LSEDIMAERO).AND.(CPROGRAM/='DIAG  ')) THEN
    CALL CH_AER_SEDIM_n(PTSTEP,                                                &
            ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_AERBEG:NSV_AEREND),               &
            XTHT(IIB:IIE,IJB:IJE,IKB:IKE),  XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), &
          XPABST(IIB:IIE,IJB:IJE,IKB:IKE), XVDEPAERO(IIB:IIE,IJB:IJE,:),       &
             XZZ(IIB:IIE,IJB:IJE,IKB:IKE),     XSEDA(IIB:IIE,IJB:IJE,IKB:IKE,:))
  ENDIF
! implicit  wet deposition
  IF ((LCH_CONV_SCAV).AND.(CPROGRAM/='DIAG  ')) THEN
    DO JN=1,NSV_AER
    ZCWETAERO(:,:,:,JN) =  (XRSVS(:,:,:,JN+NSV_AERBEG-1)+PWETDEPAER(:,:,:,JN))*PTSTEP / XRHODJ(:,:,:) 
    END DO
    ZCWETAERO(:,:,:,:)= MAX(ZCWETAERO(:,:,:,:), XMNH_TINY)

    CALL CH_AER_WETDEP_n(PTSTEP, ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_AERBEG:NSV_AEREND),             &
                         ZCWETAERO(IIB:IIE,IJB:IJE,IKB:IKE,:), XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), &
                         XSEDA(IIB:IIE,IJB:IJE,IKB:IKE,:))
  ENDIF
! explicit wet deposition
  IF ((LDEPOS_AER(IMI)).AND.(CPROGRAM/='DIAG  ')) THEN
    CALL CH_AER_DEPOS(NSPLITR, PTSTEP,                    &
                        XZZ(IIB:IIE,IJB:IJE,IKB:IKE),       &
                        XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),  &
                        XRT(IIB:IIE,IJB:IJE,IKB:IKE,:),     &
                        XRRS(IIB:IIE,IJB:IJE,IKB:IKE,:),    &
                        XRHODJ(IIB:IIE,IJB:IJE,IKB:IKE),    &
                        ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:),    &
                        XMI(IIB:IIE,IJB:IJE,IKB:IKE,:),     &
                        XTHT(IIB:IIE,IJB:IJE,IKB:IKE),      &
                        XPABST(IIB:IIE,IJB:IJE,IKB:IKE),    &
                        XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),   &
                        XSEDA(IIB:IIE,IJB:IJE,IKB:IKE,:))

  ENDIF
! Update aerosol tendency before aerosol solver
  DO JSV = 1, SIZE(XSVT,4)
    XRSVS(:,:,:,JSV) = ZSVT(:,:,:,JSV) * XRHODJ(:,:,:) / PTSTEP
  END DO
ENDIF
!
!*       3.2 check where aqueous concentration>0 + micropĥysics term
!            sedimentation, autoconversion and accretion
!
IF (LUSECHAQ.AND.(NRRL>=2) ) THEN
  DO JRR = 2, 3
    ZRT_VOL(:,:,:,JRR) = XRT(:,:,:,JRR)*XRHODREF(:,:,:)/1.e3
  END DO
  CALL CH_AQUEOUS_CHECK (PTSTEP, XRHODREF, XRHODJ, XRRS, XRSVS, NRRL, &
                         NRR, NEQ, NEQAQ, CNAMES, XRTMIN_AQ, LUSECHIC )
  IF (MAXVAL(ZRT_VOL(:,:,:,2))>XRTMIN_AQ) THEN
    SELECT CASE ( CCLOUD )
      CASE ('KESS')
        CALL CH_AQUEOUS_TMICKESS(PTSTEP, XRHODREF, XRHODJ, XRTMIN_AQ,     &
                          XRT(:,:,:,2), XRT(:,:,:,3),                     &
                          XRRS(:,:,:,2), XRRS(:,:,:,3),                   &
                          XSVT(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),  &
                          XRSVS(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2), &
                          XSVT(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND),    &
                          XRSVS(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND)    )

      CASE ('ICE3','ICE4')
        CALL CH_AQUEOUS_TMICICE(PTSTEP, XRHODREF, XRHODJ, XTHT, XPABST,     &
                          XRTMIN_AQ, LUSECHIC, LCH_RET_ICE, CNAMES,         &
                          CICNAMES, NEQ, NEQAQ,                             & 
                          XRT(:,:,:,1), XRT(:,:,:,2), XRT(:,:,:,3),         &
                          XRT(:,:,:,4), XRT(:,:,:,5), XRT(:,:,:,6),         &
                          XCIT(:,:,:), XRRS(:,:,:,2), XRRS(:,:,:,3),        &
                          XRRS(:,:,:,4), XRRS(:,:,:,5),XRRS(:,:,:,6),       &
                          XSVT(:,:,:,NSV_CHGSBEG:NSV_CHGSEND),              &
                          XRSVS(:,:,:,NSV_CHGSBEG:NSV_CHGSEND),             &
                          XSVT(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),    &
                          XRSVS(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),   &
                          XSVT(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND),      &
                          XRSVS(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND),     &
                          XSVT(:,:,:,NSV_CHICBEG:NSV_CHICEND),              &
                          XRSVS(:,:,:,NSV_CHICBEG:NSV_CHICEND)              )

      CASE ('C2R2','C3R5')
        CALL CH_AQUEOUS_TMICC2R2(PTSTEP, XRTMIN_AQ, XRHODREF, XRHODJ,          &
                          XRT(:,:,:,2), XRT(:,:,:,3),                          &
                          XRRS(:,:,:,2), XRRS(:,:,:,3),                        &
                          XSVT(:,:,:,NSV_C2R2BEG+1), XSVT(:,:,:,NSV_C2R2BEG+2),&
                          XSVT(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),       &
                          XRSVS(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),      &
                          XSVT(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND),         &
                          XRSVS(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND)         )
      CASE ('KHKO')
        CALL CH_AQUEOUS_TMICKHKO(PTSTEP, XRTMIN_AQ, XRHODREF, XRHODJ,          &
                          XRT(:,:,:,2), XRT(:,:,:,3),                          &
                          XRRS(:,:,:,2), XRRS(:,:,:,3),                        &
                          XRSVS(  :,:,:,NSV_C2R2BEG+1),                        &
                          XSVT(:,:,:,NSV_C2R2BEG+1), XSVT(:,:,:,NSV_C2R2BEG+2),&
                          XSVT(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),       &
                          XRSVS(:,:,:,NSV_CHACBEG:NSV_CHACBEG-1+NEQAQ/2),      &
                          XSVT(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND),         &
                          XRSVS(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND)         )
    END SELECT
  ENDIF
  IF (MAXVAL(ZRT_VOL(:,:,:,3))>XRTMIN_AQ) THEN
    SELECT CASE ( CCLOUD )
      CASE ('KESS','ICE3','ICE4')
        CALL CH_AQUEOUS_SEDIM1MOM(NSPLITR, CCLOUD, LUSECHIC,                    &
                                  PTSTEP , XZZ, XRHODREF,                       &
                                  XRHODJ, XRRS(:,:,:,3), XRRS(:,:,:,5),         &
                                  XRRS(:,:,:,6),                                &
                                  XRSVS(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND), &
                                  XRSVS(:,:,:,NSV_CHICBEG:NSV_CHICEND),         &
                                  ZINPRR(:,:)                                   )

      CASE ('C2R2','C3R5','KHKO')
        CALL CH_AQUEOUS_SEDIM2MOM(NSPLITR, CCLOUD, PTSTEP, XRTMIN_AQ,         &
                                  XZZ, XRHODREF, XRHODJ,                       &
                                  XRT(:,:,:,3),XRRS(:,:,:,3),                  &
                                  XSVT(:,:,:,NSV_C2R2BEG+2),                   &
                                  XRSVS(:,:,:,NSV_C2R2BEG+2),                  &
                                  XSVT(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND), &
                                  XRSVS(:,:,:,NSV_CHACBEG+NEQAQ/2:NSV_CHACEND),&
                                  ZINPRR(:,:)                                  )
    END SELECT
  ELSE
    ZINPRR(:,:) = 0.
  END IF
ELSE IF (LUSECHAQ.AND.(NRRL==1) ) THEN
  CALL CH_AQUEOUS_CHECK (PTSTEP, XRHODREF, XRHODJ, XRRS, XRSVS, NRRL, &
                         NRR, NEQ, NEQAQ, CNAMES, XRTMIN_AQ, LUSECHIC )
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.    INTEGRATE OVER ALL GRID POINTS
!              -------------------------------
!
DO JL=1,ISVECNMASK
!
!*       4.1   transfer chemical species from 4D into 1D array for solver
!              and convert from part/part to molec./cm3
!
  IDTI=ISVECMASK(2,JL)-ISVECMASK(1,JL)+1
  IDTJ=ISVECMASK(4,JL)-ISVECMASK(3,JL)+1
  IDTK=ISVECMASK(6,JL)-ISVECMASK(5,JL)+1
  IF (CSOLVER(1:2)=="RO" .AND. NEQAQ>0) THEN ! aqueous chemistry case rosenbrock solver
    CALL PREPARE_LU_AQUEOUS_ROSENBROCK !size of the jacobian matrix depending on
                                       !the presence of cloud and/or rain
  END IF
!
  IF (LORILAM) THEN
    ZRV(:) = 0.
    ZRC(:) = 0.
!ocl novrec
!cdir nodep
  DO JM=0,ISVECNPT-1
    JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
    JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
    JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
    ZSEDA(JM+1,:) = XSEDA(JI,JJ,JK,:)
    !Pressure (Pa)
    ZPRESSURE(JM+1) = XPABST(JI,JJ,JK)
    !Air density (kg/m3)
    ZDENAIR(JM+1) = XRHODREF(JI, JJ, JK)
    !Temperature (K)
    ZTEMP(JM+1)   = XTHT(JI,JJ,JK)*((XPABST(JI,JJ,JK)/XP00)**(XRD/XCPD))
    !Water vapor (kg/kg)
    IF (SIZE(XRT,4) .GE. 1) ZRV(JM+1) = XRT(JI, JJ, JK, 1)
    !Cloud vapor (kg/kg)
    IF (SIZE(XRT,4) .GE. 2) ZRC(JM+1) = XRT(JI, JJ, JK, 2)
    !Molar mass (kg/kg)
    ZMI(JM+1,:)     = XMI(JI, JJ, JK, :)
    !Moments (ppp)
    ZM(JM+1,:)      = XM3D(JI,JJ,JK,:)
    ZSIG0(JM+1,:)   = LOG(XSIG3D(JI,JJ,JK,:))
    ZRG0(JM+1,:)    =  XRG3D(JI,JJ,JK,:)  
    ZN0(JM+1,:)     =  XN3D(JI,JJ,JK,:)  
    IF (NSOA > 0) ZSOLORG(JM+1,:) = XSOLORG(JI,JJ,JK,:)
  ENDDO
  DO JN = 1, NSV_AER
!Vectorization:
!ocl novrec
!cdir nodep
      DO JM=0,ISVECNPT-1
        JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
        JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
        JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
        ZCONV(JM+1) = (XRHODREF(JI,JJ,JK)/XRHODJ(JI,JJ,JK))*ZDEN2MOL
        IF (GSPLIT) THEN
          ZAERO(JM+1,JN) = XRSVS(JI,JJ,JK,NSV_AERBEG+JN-1)*PTSTEP*ZCONV(JM+1)
        ELSE 
          ZAERO(JM+1,JN) = XSVT(JI,JJ,JK,NSV_AERBEG+JN-1)*ZDEN2MOL*XRHODREF(JI,JJ,JK)
        END IF
      END DO
    END DO
!
!* initialize  aerosol surface and aerosol diameter
!
    CALL CH_AER_SURF(ZM, ZRG0, ZSIG0, XSURF) ! Compute aerosol surface (m2/cc)
    XDP(:,:) = 2.E-6 * ZRG0(:,:) ! Mean diameter in meter 
  END IF
!
!
  IF (GSPLIT) THEN
    DO JM = 0, ISVECNPT-1
!Vectorization:
!ocl novrec
!cdir nodep
      JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
      JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
      JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
      ZCONV(JM+1) = (XRHODREF(JI,JJ,JK)/XRHODJ(JI,JJ,JK))*ZDEN2MOL
      DO JN = 1, LU_DIM_SPECIES(JM+1)
        ZCHEM(JM+1,JN) = XRSVS(JI,JJ,JK,NSV_CHEMBEG+JN-1) * PTSTEP &
                                                          * ZCONV(JM+1)
      END DO
      DO JN = 1, NEQAQ/2 ! set aqueous concentrations to zero where LW<XRTMIN_AQ
        IF (((((XRRS(JI,JJ,JK,2)/XRHODJ(JI,JJ,JK))*PTSTEP)*XRHODREF(JI,JJ,JK))/1.e3) &
                                 < XRTMIN_AQ) THEN ! cloud
          ZCHEM(JM+1,NEQ-NEQAQ+JN) = 0.
        ENDIF
        IF (((((XRRS(JI,JJ,JK,3)/XRHODJ(JI,JJ,JK))*PTSTEP)*XRHODREF(JI,JJ,JK))/1.e3) &
                                 < XRTMIN_AQ) THEN ! rain
          ZCHEM(JM+1,NEQ-NEQAQ/2+JN) = 0.
        ENDIF
      END DO 
    END DO
  ELSE 
    DO JM = 0, ISVECNPT-1
!Vectorization:
!ocl novrec
!cdir nodep
      JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
      JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
      JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
      ZCONV(JM+1) = (XRHODREF(JI,JJ,JK)/XRHODJ(JI,JJ,JK))*ZDEN2MOL
      DO JN = 1, LU_DIM_SPECIES(JM+1)
        ZCHEM(JM+1,JN) = XSVT(JI,JJ,JK,NSV_CHEMBEG+JN-1) * ZDEN2MOL &
                                                         * XRHODREF(JI,JJ,JK)
      END DO
      DO JN = 1, NEQAQ/2 ! set aqueous concentrations to zero where LW<XRTMIN_AQ
        IF (((XRT(JI,JJ,JK,2)*XRHODREF(JI,JJ,JK))/1.e3) < XRTMIN_AQ) THEN ! cloud
          ZCHEM(JM+1,NEQ-NEQAQ+JN) = 0.
        ENDIF
        IF (((XRT(JI,JJ,JK,3)*XRHODREF(JI,JJ,JK))/1.e3) < XRTMIN_AQ) THEN ! rain
          ZCHEM(JM+1,NEQ-NEQAQ/2+JN) = 0.
        ENDIF
      END DO 
    END DO
  END IF
!
!*       4.2   transfer meteo data into chemical core system
!
  SELECT CASE ( CCLOUD )
  CASE ('NONE','KESS','ICE3','ICE4')
    IF (GSPLIT) THEN ! LWC and LWR computed from tendencies
      CALL CH_METEO_TRANS_KESS(JL, XRHODJ, XRHODREF, XRRS, XTHT, XPABST, &
                             ISVECNPT, ISVECMASK, TZM, TDTCUR%TDATE%DAY, &
                             TDTCUR%TDATE%MONTH, TDTCUR%TDATE%YEAR,      &
                             XLAT, XLON, XLAT0, XLON0, LUSERV, LUSERC,   &
                             LUSERR, KLUOUT, CCLOUD, PTSTEP              )
    ELSE 
      CALL CH_METEO_TRANS_KESS(JL, XRHODJ, XRHODREF, XRT, XTHT, XPABST,  &
                             ISVECNPT, ISVECMASK, TZM, TDTCUR%TDATE%DAY, &
                             TDTCUR%TDATE%MONTH, TDTCUR%TDATE%YEAR,      &
                             XLAT, XLON, XLAT0, XLON0, LUSERV, LUSERC,   &
                             LUSERR, KLUOUT, CCLOUD                      )
    ENDIF

  CASE ('C2R2','KHKO','C3R5') !add cloud and rain C. for mean radius
    IF (GSPLIT) THEN ! LWC and LWR computed from tendencies
      CALL CH_METEO_TRANS_C2R2(JL, XRHODJ, XRHODREF, XRRS, XRSVS(:,:,:,NSV_C2R2BEG+1), &
                             XRSVS(:,:,:,NSV_C2R2BEG+2), XTHT, XPABST, ISVECNPT,       &
                             ISVECMASK, TZM, TDTCUR%TDATE%DAY, TDTCUR%TDATE%MONTH,     &
                             TDTCUR%TDATE%YEAR, XLAT,XLON, XLAT0, XLON0, LUSERV,       &
                             LUSERC, LUSERR, KLUOUT, CCLOUD,  PTSTEP                   )
    ELSE 
      CALL CH_METEO_TRANS_C2R2(JL, XRHODJ, XRHODREF, XRT, XSVT(:,:,:,NSV_C2R2BEG+1), &
                             XSVT(:,:,:,NSV_C2R2BEG+2), XTHT, XPABST, ISVECNPT,      &
                             ISVECMASK, TZM, TDTCUR%TDATE%DAY, TDTCUR%TDATE%MONTH,   &
                             TDTCUR%TDATE%YEAR, XLAT,XLON, XLAT0, XLON0, LUSERV,     &
                             LUSERC, LUSERR, KLUOUT, CCLOUD                          )
    ENDIF
  END SELECT
!
!*       4.3    calculate reaction and photolysis rates and current pH value
!
  IF (LUSECHAQ.AND.LCH_PH) THEN
    SELECT CASE(NRRL)
      CASE(1)
        DO JM=0,ISVECNPT-1
          JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
          JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
          JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
          ZPH(JM+1,1) = XPHC(JI,JJ,JK)
        END DO
      CASE(2)
        DO JM=0,ISVECNPT-1
          JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
          JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
          JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
          ZPH(JM+1,1) = XPHC(JI,JJ,JK)
          ZPH(JM+1,2) = XPHR(JI,JJ,JK)
        END DO
    END SELECT
    CALL CH_SET_RATES &
      (TDTCUR%TIME, ZCHEM, TZM, IMI, KLUOUT, KVERB, ISVECNPT, NEQ, NRRL, ZPH)
  ELSE
    CALL CH_SET_RATES &
      (TDTCUR%TIME, ZCHEM, TZM, IMI, KLUOUT, KVERB, ISVECNPT, NEQ, NRRL)
  ENDIF
!
  CALL CH_SET_PHOTO_RATES( TDTCUR%TIME, ZCHEM, JL, TZM, IMI, KLUOUT, KVERB, &
                           ISVECNPT, ISVECMASK, NEQ, XJVALUES)
!
!*       4.4 initialize aerosol parameters and moments of 0th,
!            6th, aerosol surface and aerosol diameter order
!
  IF (LORILAM) THEN
    IF (KTCOUNT == 1) THEN
          CALL CH_INI_ORILAM(ZM, ZSIG0, ZRG0, ZN0, ZCTOTG, ZCTOTA, ZCCTOT,   &
                             ZSEDA, ZOM, ZRHOP0, ZAERO, ZCHEM, ZRV, ZDENAIR, &
                             ZPRESSURE, ZTEMP, ZRC, ZFRAC, ZMI,CCH_SCHEME)
    END IF
! transfer non-volatile species from aerosol to gas-phase variables
! this line seems to be useless and transfer all H2SO4 from AP to cloud
! droplets is LUSECHAQ and LORILAM set to true
!   ZCHEM(:,JP_CH_H2SO4) = ZAERO(:,JP_CH_SO4i) + ZAERO(:,JP_CH_SO4j)
  END IF
!
!*       4.5   solve chemical system for the timestep of the monitor
!
  ZOLDCHEM(:,:) = ZCHEM(:,:)
  DO JM = 1, NCH_SUBSTEPS
    CALL CH_SOLVER_n &
          (TDTCUR%TIME, ZDTSOLVER, ZCHEM, ZNEWCHEM, NEQ, ISVECNPT, IMI)
    ZCHEM(:,:) = MAX(0.0,ZNEWCHEM(:,:))
  END DO
 IF (CSOLVER(1:2)=="RO" .AND. NEQAQ>0) THEN ! aqueous chemistry case rosenbrock solver
    DEALLOCATE(LU_IROW)
    DEALLOCATE(LU_ICOL)
    DEALLOCATE(LU_CROW)
    DEALLOCATE(LU_DIAG)
  END IF
!
!*       4.6   solve aerosol system
!
  IF (LORILAM) THEN
    ZSO4RAT(:)    = (ZNEWCHEM(:,JP_CH_H2SO4)-ZOLDCHEM(:,JP_CH_H2SO4)) / PTSTEP
    ZOLDAERO(:,:) = ZAERO(:,:)
    CALL CH_ORILAM(ZAERO,ZNEWCHEM, ZM, ZSIG0, ZRG0, ZN0,  ZCTOTG,  &
                   ZCTOTA, ZCCTOT, PTSTEP, ZSEDA,                  &
                   ZMU, ZLAMBDA, ZRHOP0, ZOM, ZSO4RAT,             &
                   ZRV, ZDENAIR,ZPRESSURE, ZTEMP, ZRC, ZFRAC, ZMI, &
                   ZTIME,CCH_SCHEME,ZSOLORG )
    ZNEWAERO(:,:) = ZAERO(:,:)
!
!*       4.7   return results to MesoNH scalar variables - aerosols
!
!Vectorization:
!ocl novrec
!cdir nodep
    DO JM=0,ISVECNPT-1
      JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
      JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
      JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
      XSIG3D(JI,JJ,JK,:)     = EXP(ZSIG0(JM+1,:))
      XRG3D(JI,JJ,JK,:)      = ZRG0(JM+1,:)
      XN3D(JI,JJ,JK,:)       = ZN0(JM+1,:)
      XRHOP3D(JI,JJ,JK,:)    = ZRHOP0(JM+1,:)
      XCTOTA3D(JI,JJ,JK,:,:) = ZCTOTA(JM+1,:,:)
      XM3D(JI,JJ,JK,:)       = ZM(JM+1,:)
      XFRAC(JI,JJ,JK,:)      = ZFRAC(JM+1,:)
      XMI(JI,JJ,JK,:)        = ZMI(JM+1,:)
    END DO
    DO JN = 1, NSV_AER
!Vectorization:
!ocl novrec
!cdir nodep
      DO JM=0,ISVECNPT-1
        JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
        JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
        JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
        IF (GSPLIT) THEN
          XRSVS(JI,JJ,JK,NSV_AERBEG+JN-1)=ZNEWAERO(JM+1,JN)/(PTSTEP*ZCONV(JM+1))
        ELSE
          XRSVS(JI,JJ,JK,NSV_AERBEG+JN-1) = XRSVS(JI,JJ,JK,NSV_AERBEG+JN-1) &
                                  + (ZNEWAERO(JM+1,JN) - ZOLDAERO(JM+1,JN)) &
                                      / (PTSTEP * ZCONV(JM+1))
        END IF
      END DO
    END DO
  END IF
!
!
!*       4.8.1  read production/loss terms for chemical species and filter
!               selected species
!
  IF (NEQ_PLT>0) THEN
    CALL CH_PRODLOSS(TDTCUR%TIME,ZCHEM,ZPRODTOT,ZLOSSTOT,IMI,ISVECNPT,NEQ)
    DO JM=1, NEQ_PLT
      DO JN=1,ISVECNPT
        ZPROD(JN,JM)=ZPRODTOT(JN,NIND_SPEC(JM))
        ZLOSS(JN,JM)=ZLOSSTOT(JN,NIND_SPEC(JM))*ZCHEM(JN,NIND_SPEC(JM))
      END DO
    END DO
  END IF
!
!
!*       4.8.2  read extended production/loss terms for chemical species and
!               filter selected species
!
  IF (NEQ_BUDGET>0) THEN
        CALL CH_TERMS(TDTCUR%TIME,ZCHEM,ZTCHEMTOT,IMI,ISVECNPT,NEQ,NREAC)
        DO JM=1,NEQ_BUDGET
          DO JN=1,ISVECNPT    
            JS=1      
            DO JO=1,NNONZEROTERMS
                IF(NSPEC_BUDGET(JM).EQ.IINDEX(1,JO)) THEN
                  ZTCHEM(JM)%ZB_REAC(JN,JS)=ZTCHEMTOT(JN,IINDEX(1,JO),IINDEX(2,JO))
                  ZTCHEM(JM)%IB_REAC(JS)=IINDEX(2,JO)
                  JS=JS+1
                END IF
            END DO
          END DO
        END DO
  END IF
!
!*       4.9   return result to MesoNH scalar variables - chemical species
!
  IF (GSPLIT) THEN
    DO JM = 0, ISVECNPT-1
!Vectorization:
!ocl novrec
!cdir nodep
      JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
      JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
      JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
      DO JN = 1, LU_DIM_SPECIES(JM+1)
        XRSVS(JI,JJ,JK,NSV_CHEMBEG+JN-1) = ZNEWCHEM(JM+1,JN) &
                            / (PTSTEP * ZCONV(JM+1))
      END DO
    END DO
  ELSE
    DO JM = 0, ISVECNPT-1
!Vectorization:
!ocl novrec
!cdir nodep
      JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
      JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
      JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
!
      DO JN = 1, LU_DIM_SPECIES(JM+1)
        XRSVS(JI,JJ,JK,NSV_CHEMBEG+JN-1) = XRSVS(JI,JJ,JK,NSV_CHEMBEG+JN-1) &
                          + (ZNEWCHEM(JM+1,JN) - ZOLDCHEM(JM+1,JN)) &
                            / (PTSTEP * ZCONV(JM+1))
      END DO
    END DO
  END IF
  IF (CSOLVER(1:2)=="RO" .AND. NEQAQ>0) THEN ! aqueous chemistry case rosenbrock solver
    DEALLOCATE(LU_DIM_SPECIES) 
  END IF
!
!*       4.10   return result to MesoNH scalar variables - pH values
!
  IF (LUSECHAQ.AND.LCH_PH) THEN
    SELECT CASE(NRRL)
      CASE(1)
        DO JM=0,ISVECNPT-1
          JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
          JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
          JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
          XPHC(JI,JJ,JK) = ZPH(JM+1,1)
        END DO
      CASE(2)
        DO JM=0,ISVECNPT-1
          JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
          JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
          JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
          XPHC(JI,JJ,JK) = ZPH(JM+1,1)
          XPHR(JI,JJ,JK) = ZPH(JM+1,2)
        END DO
    END SELECT
  ENDIF
!
!
!*      4.11 return result to MesoNH scalar variables - prod/loss terms
!
  IF (NEQ_PLT>0) THEN
        DO JM=0,ISVECNPT-1
          JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
          JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
          JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
          DO JN=1,NEQ_PLT
            XPROD(JI,JJ,JK,JN) = ZPROD(JM+1,JN)/(ZDEN2MOL*XRHODREF(JI,JJ,JK))
            XLOSS(JI,JJ,JK,JN) = ZLOSS(JM+1,JN)/(ZDEN2MOL*XRHODREF(JI,JJ,JK))
          END DO
        END DO
  END IF
!
!
!*      4.12 return result to MesoNH scalar variables - extended prod/loss terms
!
  IF (NEQ_BUDGET>0) THEN
        DO JM=0,ISVECNPT-1
          JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
          JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
          JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
          DO JN=1,NEQ_BUDGET
            DO JS=1,IIND(JN)
              XTCHEM(JN)%XB_REAC(JI,JJ,JK,JS)=(ZTCHEM(JN)%ZB_REAC(JM+1,JS))/(ZDEN2MOL*XRHODREF(JI,JJ,JK))
              XTCHEM(JN)%NB_REAC(JS)=ZTCHEM(JN)%IB_REAC(JS)
            END DO
          END DO
        END DO
  END IF
!
!
END DO
!
!*        4.13  compute accumalated concentrations in rain at the surface
!
IF (CCLOUD /= 'REVE' ) THEN
  IF (LUSECHAQ) THEN
    DO JSV=1,NSV_CHAC/2
      WHERE((XRRS(:,:,IKB,3) .GT. 0.).AND.(XRSVS(:,:,IKB,JSV+NSV_CHACBEG+NSV_CHAC/2-1).GT.0.))
          XACPRAQ(:,:,JSV) = XACPRAQ(:,:,JSV) + &
              (XRSVS(:,:,IKB,JSV+NSV_CHACBEG+NSV_CHAC/2-1))/ (XMD*XRRS(:,:,IKB,3))*& ! moles i  / kg eau
               1E3*ZINPRR(:,:) * XTSTEP ! moles i / m2
      END WHERE
    ENDDO
    IF (LCH_PH) THEN
      WHERE ((ZINPRR(:,:)>0.).AND.(XPHR(:,:,IKB)>0.))
      ! moles of H+ / m2
        XACPHR(:,:) =  XACPHR(:,:) + 1E3*ZINPRR(:,:) * XTSTEP * &
                     10**(-XPHR(:,:,IKB))
      END WHERE
    END IF
  END IF
END IF

! Correction of negativity
!
DO JSV = 1, SIZE(XSVT,4)
 XRSVS(:,:,:,JSV)   = MAX((XRSVS(:,:,:,JSV)),XSVMIN(JSV))
END DO
!
IF (LBUDGET_SV) THEN
  DO JSV=NSV_CHEMBEG,NSV_CHEMEND
    CALL BUDGET(XRSVS(:,:,:,JSV),JSV+12,'CHEM_BU_RSV')
  ENDDO
ENDIF
!
!----------------------------------------------------------------------
!
IF ((CPROGRAM =='DIAG  ').OR.(L1D)) THEN
  CALL WRITE_TS1D
END IF
!
DEALLOCATE(TZM)
DEALLOCATE(ZCHEM)
DEALLOCATE(ZNEWCHEM)
DEALLOCATE(ZOLDCHEM)
DEALLOCATE(ZCONV)
IF (LUSECHAQ.AND.LCH_PH) DEALLOCATE(ZPH)
!
IF (NEQ_PLT>0) THEN
  DEALLOCATE(ZPRODTOT)
  DEALLOCATE(ZLOSSTOT)
  DEALLOCATE(ZPROD)
  DEALLOCATE(ZLOSS)
END IF
IF (NEQ_BUDGET>0) THEN
  DEALLOCATE(ZTCHEMTOT)
  DEALLOCATE(ZTCHEM)
  DEALLOCATE(IIND)
  DEALLOCATE(IINDEX)
END IF
IF (LORILAM) THEN
  DEALLOCATE(ZAERO)
  DEALLOCATE(ZNEWAERO)
  DEALLOCATE(ZOLDAERO)
  DEALLOCATE(ZM)
  DEALLOCATE(ZSEDA)
  DEALLOCATE(ZN0)
  DEALLOCATE(ZRG0)
  DEALLOCATE(ZSIG0)
  DEALLOCATE(ZRHOP0)
  DEALLOCATE(ZCTOTA)
  DEALLOCATE(ZCCTOT)
  DEALLOCATE(ZCTOTG)
  DEALLOCATE(ZMU)
  DEALLOCATE(ZLAMBDA)
  DEALLOCATE(ZOM)
  DEALLOCATE(ZSO4RAT)
  DEALLOCATE(ZRV)
  DEALLOCATE(ZRC)
  DEALLOCATE(ZPRESSURE)
  DEALLOCATE(ZTEMP)
  DEALLOCATE(ZDENAIR)
  DEALLOCATE(ZFRAC)
  DEALLOCATE(ZMI)
  DEALLOCATE(ZSOLORG)
  DEALLOCATE(XDP)
  DEALLOCATE(XSURF)
END IF
!-------------------------------------------------------------------------------
!
CONTAINS
!
  SUBROUTINE PREPARE_LU_ROSENBROCK
!
USE MODI_CH_INIT_ROSENBROCK
!
!   local variables
!
INTEGER                :: JISHIFT        ! shift index in a loop
INTEGER                :: JILOCAL        ! shift index in a loop
INTEGER                :: ILAST          ! last elemnt of NSPARSE_DIAG vector
!
!-------------------------------------------------------------------------------
!
    CALL CH_INIT_ROSENBROCK(IMI,KLUOUT)
!
! add vectorization of the LU_arrays created by CH_INIT_ROSENBROCK
!
    LU_NONZERO = NSPARSEDIM*ISVECNPT
    ALLOCATE(LU_IROW(LU_NONZERO))
    ALLOCATE(LU_ICOL(LU_NONZERO))
    DO JI = 1, ISVECNPT
      JISHIFT = NSPARSEDIM*(JI-1)
      JILOCAL = NEQ*(JI-1)
      LU_IROW(JISHIFT+1:JISHIFT+NSPARSEDIM) = NSPARSE_IROW(1:NSPARSEDIM)+JILOCAL
      LU_ICOL(JISHIFT+1:JISHIFT+NSPARSEDIM) = NSPARSE_ICOL(1:NSPARSEDIM)+JILOCAL
    END DO
!
    NVAR = NEQ*ISVECNPT
    ALLOCATE(LU_CROW(NVAR+1))
    ALLOCATE(LU_DIAG(NVAR+1))
    ILAST = NSPARSE_DIAG(NEQ)
    DO JI = 1, ISVECNPT
      JISHIFT = NEQ*(JI-1)
      JILOCAL = ILAST*(JI-1)
      LU_CROW(JISHIFT+1:JISHIFT+NEQ) = NSPARSE_CROW(1:NEQ)+JILOCAL
      LU_DIAG(JISHIFT+1:JISHIFT+NEQ) = NSPARSE_DIAG(1:NEQ)+JILOCAL
    END DO
    LU_CROW(NVAR+1) = LU_NONZERO+1
    LU_DIAG(NVAR+1) = LU_NONZERO+1
  RETURN
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREPARE_LU_ROSENBROCK
!
!
  SUBROUTINE PREPARE_LU_AQUEOUS_ROSENBROCK
!
USE MODI_CH_INIT_ROSENBROCK
USE MODD_CH_MNHC_n, ONLY : XRTMIN_AQ
!
!   local variables
!
INTEGER                :: JISHIFT        ! shift index in a loop
INTEGER                :: JILOCAL        ! shift index in a loop
INTEGER                :: ILAST          ! last elemnt of NSPARSE_DIAG vector
REAL,    DIMENSION(SIZE(XRRS,1),SIZE(XRRS,2),SIZE(XRRS,3),SIZE(XRRS,4))   &
                                :: ZRRS      ! work array
LOGICAL, DIMENSION(:), ALLOCATABLE   :: GWATER
INTEGER, DIMENSION(:), ALLOCATABLE   :: IMASKAQ
INTEGER                :: IWATER
INTEGER                :: ISPARSEDIM
INTEGER                :: IEQ
INTEGER                :: ILAST_NAQ
INTEGER                :: JRR           ! Loop index for the moist variables
REAL                   :: ZRTMIN_AQ
!
!-------------------------------------------------------------------------------
!
DO JRR = 2, NRRL+1
  ZRRS(:,:,:,JRR)  = XRRS(:,:,:,JRR) / XRHODJ(:,:,:)
END DO
!
!-------------------------------------------------------------------------------
!
! Same as in PREPARE_LU_ROSENBROCK but in the case of non-homogeneous
! chemical systems with are put together, here a mixture of NEQ and NEQ_NAQ 
! system dimensions.
!
    IF (KTCOUNT == 1) THEN
      IF (JL==1) CALL CH_INIT_ROSENBROCK(IMI,KLUOUT)
      IF( ASSOCIATED(LU_DIM_SPECIES) ) THEN
        DEALLOCATE(LU_DIM_SPECIES)
      END IF
    END IF
!
! Create the GWATER mask
!
    ALLOCATE(GWATER(ISVECNPT))
    GWATER(:) = .FALSE.
    IF (GSPLIT) THEN
      ZRTMIN_AQ = XRTMIN_AQ/PTSTEP
      SELECT CASE ( CCLOUD )
        CASE('REVE')
          DO JM=0,ISVECNPT-1
            JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
            JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
            JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
            GWATER(JM+1) = XRRS(JI,JJ,JK,2)>(ZRTMIN_AQ*1.e3/XRHODREF(JI,JJ,JK))
          END DO
        CASE('KESS','ICE3','ICE4','C2R2','C3R5','KHKO')
          DO JM=0,ISVECNPT-1
            JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
            JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
            JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
            GWATER(JM+1) = XRRS(JI,JJ,JK,2)>(ZRTMIN_AQ*1.e3/XRHODREF(JI,JJ,JK)) &
                      .OR. XRRS(JI,JJ,JK,3)>(ZRTMIN_AQ*1.e3/XRHODREF(JI,JJ,JK))
          END DO
      END SELECT
    ELSE 
      SELECT CASE ( CCLOUD )
        CASE('REVE')
          DO JM=0,ISVECNPT-1
            JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
            JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
            JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
            GWATER(JM+1) = XRT(JI,JJ,JK,2)>(XRTMIN_AQ*1.e3/XRHODREF(JI,JJ,JK))
          END DO
        CASE('KESS','ICE3','ICE4','C2R2','C3R5','KHKO')
          DO JM=0,ISVECNPT-1
            JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
            JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+ISVECMASK(3,JL)
            JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
            GWATER(JM+1) = XRT(JI,JJ,JK,2)>(XRTMIN_AQ*1.e3/XRHODREF(JI,JJ,JK)) &
                      .OR. XRT(JI,JJ,JK,3)>(XRTMIN_AQ*1.e3/XRHODREF(JI,JJ,JK))
          END DO
      END SELECT
    END IF
    IWATER = COUNT(GWATER(:))
    ALLOCATE(IMASKAQ(ISVECNPT)); IMASKAQ(:) = 0
    IF( IWATER>=1 ) THEN
      WHERE( GWATER(:) ) 
        IMASKAQ(:) = 1
      END WHERE
    END IF
    DEALLOCATE(GWATER)
!
! add vectorization of the LU_arrays created by CH_INIT_ROSENBROCK
! but now taking into account a full system to solve with NEQ species
! (gazeous+aqueous species) and a reduced system with NEQ_NAQ<NEQ 
! species (pure gazeous case)
!
    ALLOCATE(LU_DIM_SPECIES(ISVECNPT))
    LU_DIM_SPECIES(:) = NEQ*IMASKAQ(:) + NEQ_NAQ*(1-IMASKAQ(:))
!
    LU_NONZERO = NSPARSEDIM*IWATER + NSPARSEDIM_NAQ*(ISVECNPT-IWATER)
    ALLOCATE(LU_IROW(LU_NONZERO))
    ALLOCATE(LU_ICOL(LU_NONZERO))
    JISHIFT = 0
    JILOCAL = 0
    DO JI = 1, ISVECNPT
      ISPARSEDIM = NSPARSEDIM*IMASKAQ(JI) + NSPARSEDIM_NAQ*(1-IMASKAQ(JI))
      IF( ISPARSEDIM==NSPARSEDIM ) THEN
        LU_IROW(JISHIFT+1:JISHIFT+ISPARSEDIM)=NSPARSE_IROW(1:ISPARSEDIM)+JILOCAL
        LU_ICOL(JISHIFT+1:JISHIFT+ISPARSEDIM)=NSPARSE_ICOL(1:ISPARSEDIM)+JILOCAL
      ELSE
        LU_IROW(JISHIFT+1:JISHIFT+ISPARSEDIM)=NSPARSE_IROW_NAQ(1:ISPARSEDIM)+ &
                                                                         JILOCAL
        LU_ICOL(JISHIFT+1:JISHIFT+ISPARSEDIM)=NSPARSE_ICOL_NAQ(1:ISPARSEDIM)+ &
                                                                         JILOCAL
      END IF
      JISHIFT = JISHIFT + ISPARSEDIM
      JILOCAL = JILOCAL + NEQ*IMASKAQ(JI) + NEQ_NAQ*(1-IMASKAQ(JI))
    END DO
!
    NVAR = NEQ*IWATER + NEQ_NAQ*(ISVECNPT-IWATER)
    ALLOCATE(LU_CROW(NVAR+1))
    ALLOCATE(LU_DIAG(NVAR+1))
    JISHIFT = 0
    JILOCAL = 0
    ILAST     = NSPARSE_DIAG(NEQ)
    ILAST_NAQ = NSPARSE_DIAG_NAQ(NEQ_NAQ)
    DO JI = 1, ISVECNPT
      IEQ = LU_DIM_SPECIES(JI)
      IF( IEQ==NEQ ) THEN
        LU_CROW(JISHIFT+1:JISHIFT+IEQ) = NSPARSE_CROW(1:IEQ)+JILOCAL
        LU_DIAG(JISHIFT+1:JISHIFT+IEQ) = NSPARSE_DIAG(1:IEQ)+JILOCAL
      ELSE
        LU_CROW(JISHIFT+1:JISHIFT+IEQ) = NSPARSE_CROW_NAQ(1:IEQ)+JILOCAL
        LU_DIAG(JISHIFT+1:JISHIFT+IEQ) = NSPARSE_DIAG_NAQ(1:IEQ)+JILOCAL
      END IF
      JISHIFT = JISHIFT + IEQ
      JILOCAL = JILOCAL + ILAST*IMASKAQ(JI) + ILAST_NAQ*(1-IMASKAQ(JI))
    END DO
    LU_CROW(NVAR+1) = LU_NONZERO+1
    LU_DIAG(NVAR+1) = LU_NONZERO+1
!
    DEALLOCATE(IMASKAQ)
  RETURN
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREPARE_LU_AQUEOUS_ROSENBROCK
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_MONITOR_n

