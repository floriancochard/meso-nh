!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!     #######################
      PROGRAM PREP_IDEAL_CASE
!     #######################
!
!!****  *PREP_IDEAL_CASE* - program to write an initial FM-file 
!!
!!    PURPOSE
!!    -------
!       The purpose of this program is to prepare an initial meso-NH file
!     (LFIFM and DESFM files) filled with some idealized fields.    
!
!      ---- The present version can provide two types of fields:
!
!      1) CIDEAL = 'CSTN' : 3D fields derived  from a vertical profile with
!         ---------------   n levels of constant moist Brunt Vaisala frequency
!             The vertical profile is read in EXPRE file.                 
!             These fields can be used for model runs 
!
!      2) CIDEAL = 'RSOU' : 3D fields derived from a radiosounding.
!          --------------- 
!             The radiosounding is read in EXPRE file. 
!             The following kind of data  is permitted :
!                  YKIND = 'STANDARD'  :   Zsol, Psol, Tsol, TDsol
!                                         (Pressure, dd, ff) , 
!                                         (Pressure, T, Td)
!                  YKIND = 'PUVTHVMR'  : zsol, Psol, Thvsol, Rsol
!                                        (Pressure, U, V) , 
!                                        (Pressure, THv, R)
!                  YKIND = 'PUVTHVHU'  :  zsol, Psol, Thvsol, Husol
!                                         (Pressure, U, V) , 
!                                         (Pressure, THv, Hu)
!                  YKIND = 'ZUVTHVHU'  :  zsol, Psol, Thvsol, Husol
!                                         (height, U, V) , 
!                                         (height, THv, Hu)
!                  YKIND = 'ZUVTHVMR'  :  zsol, Psol, Thvsol, Rsol
!                                         (height, U, V) , 
!                                         (height, THv, R)
!                  YKIND = 'PUVTHDMR'  : zsol, Psol, Thdsol, Rsol
!                                         (Pressure, U, V) , 
!                                         (Pressure, THd, R)
!                  YKIND = 'PUVTHDHU'  : zsol, Psol, Thdsol, Husol
!                                         (Pressure, U, V) , 
!                                         (Pressure, THd, Hu)
!                  YKIND = 'ZUVTHDMR'  :  zsol, Psol, Thdsol, Rsol
!                                         (height, U, V) , 
!                                         (height, THd, R)
!                  YKIND = 'ZUVTHLMR'  :  zsol, Psol, Thdsol, Rsol
!                                         (height, U, V) , 
!                                         (height, THl, Rt)
!
!             These fields can be used for model runs 
!
!      Cases (1) and (2) can be balanced
!      (geostrophic, hydrostatic  and anelastic balances) if desired.
!
!      ---- The orography can be flat (YZS='FLAT'), but also 
!      sine-shaped (YZS='SINE') or  bell-shaped (YZS='BELL')
!
!      ---- The U(z)  profile given in the RSOU and CSTN cases can
!      be multiplied (CUFUN="Y*Z") by a function of y (function FUNUY)  
!      The V(z) profile  given in the RSOU and CSTN cases can
!      be multiplied (CVFUN="X*Z") by a function of x (function FUNVX). 
!      If it is not the case, i.e. U(y,z)=U(z) then CUFUN="ZZZ" and 
!      CVFUN="ZZZ" for V(y,z)=V(z). Instead of these separable forms,
!      non-separables functions FUNUYZ (CUFUN="Y,Z")  and FUNVXZ (CVFUN="X,Z") 
!      can be used to specify the wind components.
!
!!**  METHOD
!!    ------
!!      The directives and data to perform the preparation of the initial FM
!!    file are stored in EXPRE file. This file is composed  of two parts : 
!!          - a namelists-format  part which is present in all cases
!!          - a free-format  part which contains data in cases 
!!       of discretised orography (CZS='DATA')
!!       of radiosounding (CIDEAL='RSOU') or Nv=cste  profile (CIDEAL='CSTN')
!!       of forced version (LFORCING=.TRUE.)
!!    
!!
!!      The following  PREP_IDEAL_CASE program  :
!!
!!             - initializes physical constants by calling INI_CST 
!!
!!             - sets default values for global variables which will be 
!!     written  in DESFM file and for variables in EXPRE file (namelists part)
!!     which will be written in LFIFM file.    
!!
!!             - reads the namelists part of EXPRE file which gives 
!!     informations about the preinitialization to perform,
!!
!!             - allocates memory for arrays, 
!!
!!             - initializes fields depending on the 
!!              directives  (CIDEAL in namelist NAM_CONF_PRE) :
!!  
!!                * grid variables : 
!!                  The gridpoints are regularly spaced by XDELTAX, XDELTAY.
!!               The grid is stretched along the z direction, the mesh varies 
!!               from XDZGRD near the ground to XDZTOP near the top and the 
!!               weigthing function is a TANH function characterized by its 
!!               center and width above and under this center
!!                  The orography is initialized following the kind of orography
!!               (YZS in namelist NAM_CONF_PRE) and the degrees of freedom :
!!                     sine-shape ---> ZHMAX, IEXPX,IEXPY
!!                     bell-shape ---> ZHMAX, ZAX,ZAY,IIZS,IJZS
!!                  The horizontal grid variables are initialized following
!!                the kind of geometry (LCARTESIAN in namelist NAM_CONF_PRE) 
!!                and the grid parameters XLAT0,XLON0,XBETA in both geometries
!!                and XRPK,XLONORI,XLATORI  in conformal projection.
!!                  In the  case of initialization from a radiosounding, the
!!                date and time is read in free-part of the EXPRE file. In other
!!                cases year, month and day are set to NUNDEF and time to 0.
!!
!!               * prognostic fields : 
!!
!!                     U,V,W, Theta and r. are first determined. They are
!!                multiplied by rhoj after the anelastic reference state 
!!                computation.
!!                     For the CSTN and RSOU cases, the determination of 
!!                Theta and rv is performed  respectively by SET_RSOU
!!                and by SET_CSTN which call the common routine SET_MASS. 
!!                These three routines have  the following actions :
!!          ---   The input vertical profile   is converted in 
!!                variables (U,V,thetav,r) and  interpolated
!!                on the vertical grid of the model without orography. 
!!          ---   thetav  is computed on the model grid without orography 
!!                by integration of thermal wind balance. ( A variation of
!!                the u-wind component( x-model axis component)  is possible in
!!                y direction, a variation of the v-wind component (y-model
!!                axis component) is possible in x direction.).
!!          ---   The mass fields (theta and r ) and the wind components are 
!!                then interpolated on the model grid with orography.         
!!          ---   An  anelastic correction is  applied in PRESSURE_IN_PREP in
!!                the case of non-vanishing orography.    
!!            
!!               * anelastic reference state variables :
!!
!!                   1D reference state : 
!!                     RSOU and CSTN cases : rhorefz and thvrefz are computed 
!!                         by  SET_MASS (called by SET_RSOU or SET_CSTN).
!!                         They are deduced from thetav and r on the model grid
!!                         without orography.
!!                   The 3D reference state is  computed by SET_REF   
!!            
!!               * The total mass of dry air is computed by TOTAL_DMASS              
!!
!!             - writes the DESFM file, 
!!
!!             - writes the LFIFM file . 
!!
!!    EXTERNAL
!!    --------
!!      DEFAULT_DESFM : to set default values for variables which can be 
!!                      contained in DESFM file
!!      DEFAULT_EXPRE : to  set default values for other global variables 
!!                      which can be contained in namelist-part of EXPRE file
!!      Module MODE_GRIDPROJ : contains conformal projection routines
!!           SM_GRIDPROJ   : to compute some grid variables, in
!!                           case of conformal projection.
!!      Module MODE_GRIDCART : contains cartesian geometry routines
!!           SM_GRIDCART   : to compute some grid variables, in
!!                           case of cartesian geometry.
!!      SET_RSOU      : to initialize mass fields from a radiosounding
!!      SET_CSTN      : to initialize mass fields from a vertical profile of 
!!                      n layers of Nv=cste 
!!      SET_REF       : to compute  rhoJ 
!!      RESSURE_IN_PREP : to apply an anelastic correction in the case of
!!                        non-vanishing orography 
!!      FMOPEN        : to open a FM-file (DESFM + LFIFM)
!!      WRITE_DESFM   : to write the  DESFM file
!!      WRI_LFIFM     : to write the   LFIFM file  
!!      FMCLOS        : to close a FM-file (DESFM + LFIFM)
!!
!!      MXM,MYM,MZM   : Shuman operators
!!      WGUESS        : to compute W with the continuity equation from 
!!                      the U,V values 
!!
!!
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains parameters
!!      Module MODD_DIM1      : contains dimensions 
!!      Module MODD_CONF       : contains  configuration variables for 
!!                                all models
!!      Module MODD_CST        : contains physical constants
!!      Module MODD_GRID       : contains grid variables  for all models
!!      Module MODD_GRID1     : contains grid variables
!!      Module MODD_TIME      : contains time variables for all models  
!!      Module MODD_TIME1     : contains time variables  
!!      Module MODD_REF        : contains reference state variables for
!!                               all models
!!      Module MODD_REF1      : contains reference state variables 
!!      Module MODD_LUNIT      : contains variables which concern names
!!                            and logical unit numbers of files  for all models
!!      Module MODD_FIELD1    : contains prognostics  variables
!!      Module MODD_GR_FIELD1 : contains the surface prognostic variables 
!!      Module MODD_LSFIELD1    : contains Larger Scale fields
!!      Module MODD_DYN1        : contains dynamic control variables for model 1
!!      Module MODD_LBC1        : contains lbc control variables for model 1
!!
!!
!!      Module MODN_CONF1    : contains  configuration variables for model 1
!!                               and the NAMELIST list
!!      Module MODN_LUNIT1    : contains variables which concern names
!!                               and logical unit numbers of files and 
!!                               the NAMELIST list
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (program PREP_IDEAL_CASE)
!!    
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                            05/05/94
!!      updated                V. Ducrocq   27/06/94   
!!      updated                P.M.         27/07/94
!!      updated                V. Ducrocq   23/08/94 
!!      updated                V. Ducrocq   01/09/94 
!!      namelist changes       J. Stein     26/10/94 
!!      namelist changes       J. Stein     04/11/94 
!!      remove the second step of the geostrophic balance 14/11/94 (J.Stein)
!!      add grid stretching in the z direction + Larger scale fields +
!!      cleaning                                           6/12/94 (J.Stein) 
!!      periodize the orography and the grid sizes in the periodic case
!!                                                        19/12/94 (J.Stein) 
!!      correct a bug in the Larger Scale Fields initialization
!!                                                        19/12/94 (J.Stein) 
!!      add the vertical grid stretching                  02/01/95 (J. Stein)
!!      Total mass of dry air computation                 02/01/95 (J.P.Lafore) 
!!      add the 1D switch                                 13/01/95 (J. Stein)
!!      enforce a regular vertical grid if desired        18/01/95 (J. Stein)
!!      add the tdtcur initialization                     26/01/95 (J. Stein)
!!      bug in the test of the type of RS localization    25/02/95 (J. Stein)
!!      remove R from the historical variables            16/03/95 (J. Stein)
!!      error on the grid stretching                      30/06/95 (J. Stein)
!!      add the soil fields                               01/09/95 (S.Belair)
!!      change the streching function  and the wind guess
!!        (J. Stein and V.Masson)                         21/09/95 
!!      reset to FALSE LUSERC,..,LUSERH                   12/12/95 (J. Stein)
!!      enforce the RS localization in 1D and 2D config.
!!      + add the 'TSZ0' option for the soil variables    28/01/96 (J. Stein)
!!      initialization of domain from center point        31/01/96 (V. Masson)
!!      add the constant file reading                     05/02/96 (J. Stein)
!!      enter vertical model levels values                20/10/95 (T.Montmerle)
!!      add LFORCING option                               19/02/96 (K. Suhre)
!!      modify structure of NAM_CONF_PRE                  20/02/96 (J.-P. Pinty)
!!      default of the domain center when use of pgd file 12/03/96 (V. Masson)
!!      change the surface initialization                 20/03/96 ( Stein,
!!                                                    Bougeault, Kastendeutsch )
!!      change the DEFAULT_DESFMN CALL                    17/04/96 ( Lafore )
!!      set the STORAGE_TYPE to 'TT' (a single instant)   30/04/96 (Stein, 
!!                                                    Jabouille)
!!      new wguess to spread  the divergence              15/05/96 (Stein)
!!      set LTHINSHELL to TRUE + return to the old wguess 29/08/96 (Stein)
!!      MY_NAME and DAD_NAME writing for nesting          30/07/96 (Lafore)
!!      MY_NAME and DAD_NAME reading in pgd file          26/09/96 (Masson)
!!       and reading of pgd grid in a new routine
!!      XXHAT and XYHAT are set to 0. at origine point    02/10/96 (Masson)
!!      add LTHINSHELL in namelist NAM_CONF_PRE           08/10/96 (Masson)
!!      restores use of TS and T2                         26/11/96 (Masson)
!!      value  XUNDEF for soil and vegetation fields on sea 27/11/96 (Masson)
!!      use of HUG and HU2 in both ISBA and TSZ0 cases    04/12/96 (Masson)
!!      add initialization of chemical variables          06/08/96 (K. Suhre)
!!      add MANUAL option for the terrain elevation       12/12/96 (J.-P. Pinty)
!!      set DATA instead of MANUAL for the terrain
!!      elevation option
!!      add new anelastic equations' systems              29/06/97 (Stein)
!!      split mode_lfifm_pgd                              29/07/97 (Masson)
!!      add directional z0 and subgrid scale orography    31/07/97 (Masson)
!!      separates surface treatment in PREP_IDEAL_SURF    15/03/99 (Masson)
!!      new PGD fields allocations                        15/03/99 (Masson)
!!      iterative call to pressure solver                 15/03/99 (Masson)
!!      removes TSZ0 case                                 04/01/00 (Masson)
!!      parallelization                                   18/06/00 (Pinty)
!!      adaptation for patch approach                     02/07/00 (Solmon/Masson)
!!      bug in W LB field on Y direction                  05/03/01 (Stein)
!!      add module MODD_NSV for NSV variable              01/02/01 (D. Gazen) 
!!      allow namelists in different orders               15/10/01 (I. Mallet)
!!      allow LUSERC and LUSERI in 1D configuration       05/06/02 (P. Jabouille)
!!      add  ZUVTHLMR case (move in set_rsou latter)      05/12/02 Jabouille/Masson
!!      move LHORELAX_SV (after INI_NSV)                  30/04/04 (Pinty)
!!      Correction Parallel bug IBEG & IDEND  evalution   13/11/08 J.Escobar
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS       ! Declarative modules
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CST
USE MODD_GRID         
USE MODD_GRID_n
USE MODD_METRICS_n
USE MODD_PGDDIM
USE MODD_PGDGRID
USE MODD_TIME
USE MODD_TIME_n
USE MODD_REF
USE MODD_REF_n
USE MODD_LUNIT
USE MODD_FIELD_n
USE MODD_DYN_n
USE MODD_LBC_n
USE MODD_LSFIELD_n
USE MODD_PARAM_n
USE MODD_CH_MNHC_n, ONLY:  LUSECHEM, LCH_INIT_FIELD, CCHEM_INPUT_FILE 
USE MODD_CH_M9,     ONLY:  NEQ
USE MODD_CH_AEROSOL,ONLY:  LORILAM, CORGANIC, LVARSIGI, LVARSIGJ
USE MODD_DUST,      ONLY:  LDUST
USE MODD_SALT,      ONLY:  LSALT
USE MODD_VAR_ll,    ONLY:  NPROC
USE MODD_LUNIT_n
USE MODD_CONF_n
USE MODD_NSV,      ONLY : NSV,NSV_CHEM,           &
                          NSV_DSTEND, NSV_DSTBEG
!
USE MODN_BLANK
!
USE MODE_THERMO
USE MODE_POS
USE MODE_GRIDCART         ! Executive modules
USE MODE_GRIDPROJ
USE MODE_FM
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODE_MODELN_HANDLER
!
USE MODI_DEFAULT_DESFM_n    ! Interface modules
USE MODI_DEFAULT_EXPRE
USE MODI_READ_HGRID
USE MODI_SHUMAN
USE MODI_SET_RSOU
USE MODI_SET_CSTN
USE MODI_SET_FRC
USE MODI_PRESSURE_IN_PREP
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_LFIFM_n
USE MODI_METRICS
USE MODI_UPDATE_METRICS
USE MODI_SET_REF
USE MODI_SET_PERTURB
USE MODI_TOTAL_DMASS
USE MODI_WGUESS
USE MODI_CH_INIT_SCHEME
USE MODI_CH_INIT_CCS
USE MODI_CH_INIT_FIELD_n
USE MODI_GATHER_ll
USE MODI_INI_NSV
USE MODI_READ_PRE_IDEA_NAM_n
USE MODI_CH_AER_INIT_SOA
USE MODI_ZSMT_PIC
USE MODI_READ_VER_GRID
!
!JUAN
USE MODE_SPLITTINGZ_ll
USE MODD_SUB_MODEL_n
USE MODE_MNH_TIMING
!JUAN
IMPLICIT NONE
!
!*       0.1  Declarations of global variables not declared in the modules
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XJ ! Jacobian
REAL :: XLATCEN=XUNDEF, XLONCEN=XUNDEF ! latitude and longitude of the center of
                                     ! the domain for initialization. This 
                                     ! point is vertical vorticity point
                                     !          ------------------------
REAL :: XDELTAX=0.5E4, XDELTAY=0.5E4 ! horizontal mesh lengths  
                                     !  used to determine  XXHAT,XYHAT
!
INTEGER :: NLUPRE,NLUOUT           ! Logical unit numbers for EXPRE file
                                   ! and for output_listing file
INTEGER :: NRESP                   ! return code in FM routines
INTEGER :: NTYPE                   ! type of file (cpio or not)
INTEGER :: NNPRAR                  ! number of articles predicted  in
                                   !  the LFIFM file
INTEGER :: NNINAR                  ! number of articles  present in
                                   !  the LFIFM file
LOGICAL :: GFOUND                  ! Return code when searching namelist
!
INTEGER :: JLOOP,JILOOP,JJLOOP     ! Loop indexes
!
INTEGER :: NIB,NJB,NKB             ! Begining useful area  in x,y,z directions
INTEGER :: NIE,NJE                 ! Ending useful area  in x,y directions
INTEGER :: NIU,NJU,NKU             ! Upper bounds in x,y,z directions
CHARACTER (LEN=32) ::  CEXPRE            ! name of the EXPRE file
CHARACTER (LEN=32) :: CDESFM             ! Name of DESFM file 
CHARACTER(LEN=4)   :: CIDEAL ='CSTN'     ! kind of idealized fields
                                         ! 'CSTN' : Nv=cste case 
                                         ! 'RSOU' : radiosounding case
CHARACTER(LEN=4)   :: CZS    ='FLAT'     ! orography selector
                                         ! 'FLAT' : zero orography
                                         ! 'SINE' : sine-shaped orography 
                                         ! 'BELL' : bell-shaped orography 
REAL    :: XHMAX=XUNDEF            ! Maximum height for orography
REAL    :: NEXPX=3,NEXPY=1         ! Exponents for  orography in case of CZS='SINE'
REAL    :: XAX= 1.E4, XAY=1.E4     ! Widths for orography in case CZS='BELL'
                                   ! along x and y 
INTEGER :: NIZS = 5, NJZS = 5      ! Localization of the center in 
                                   ! case CZS ='BELL' 
!
!*       0.1.1 Declarations of local variables for N=cste and 
!              radiosounding cases :
!
INTEGER            :: NYEAR,NMONTH,NDAY ! year, month and day in EXPRE file
REAL               :: XTIME             ! time in EXPRE file
LOGICAL            :: LPERTURB =.FALSE. ! Logical to add a perturbation to 
                                        ! a basic state 
LOGICAL            :: LBOUSS   =.FALSE. ! Logical to obtain a Boussinesq   
                                        ! version ( Nref = 0 + rho = cte)
LOGICAL            :: LGEOSBAL =.TRUE.  ! Logical to satisfy the geostrophic
                                        ! balance
                                        ! .TRUE. for geostrophic balance
                                        ! .FALSE. to ignore this balance
LOGICAL            :: LPV_PERT =.FALSE. ! Logical to add a PV pertubation
LOGICAL            :: LRMV_BL  =.FALSE. ! Logical to remove the boundary layer
                                        ! before PV inversion
CHARACTER(LEN=3)   :: CFUNU ='ZZZ'      ! CHARACTER STRING for variation of
                                        ! U in y direction
                                        ! 'ZZZ'  : U = U(Z)
                                        ! 'Y*Z'  : U = F(Y) * U(Z)
                                        ! 'Y,Z'  : U = G(Y,Z)
CHARACTER(LEN=3)   :: CFUNV ='ZZZ'      ! CHARACTER STRING for variation of
                                        ! V in x direction
                                        ! 'ZZZ'  : V = V(Z)
                                        ! 'Y*Z'  : V = F(X) * V(Z)
                                        ! 'Y,Z'  : V = G(X,Z)
CHARACTER(LEN=6)   :: CTYPELOC='IJGRID' ! Type of informations  used to give the
                                        ! localization of vertical profile
                                        ! 'IJGRID'  for (i,j) point  on index space
                                        ! 'XYHATM' for (x,y) coordinates on
                                        !  conformal or cartesian plane
                                        ! 'LATLON' for (latitude,longitude) on
                                        !   spherical earth  
REAL               :: XLATLOC= 45., XLONLOC=0.
                                        ! Latitude and longitude of the vertical
                                        ! profile localization  (used in case 
                                        ! CTYPELOC='LATLON') 
REAL               :: XXHATLOC=2.E4, XYHATLOC=2.E4 
                                        ! (x,y) of the vertical profile
                                        ! localization  (used in cases 
                                        ! CTYPELOC='LATLON' and 'XYHATM') 
INTEGER, DIMENSION(1) :: NILOC=4, NJLOC=4 
                                        ! (i,j) of the vertical profile
                                        ! localization 
!
!
REAL,DIMENSION(:,:,:),ALLOCATABLE   :: XCORIOZ ! Coriolis parameter (this
                                                 ! is exceptionnaly a 3D array
                                                 ! for computing needs)
!
!
!*       0.1.2 Declarations of local variables used when a PhysioGraphic Data
!              file is used :
!
INTEGER             :: JSV                      ! loop index on scalar var.
CHARACTER(LEN=28)   :: CPGD_FILE=' '            ! Physio-Graphic Data file name
LOGICAL  :: LREAD_ZS = .TRUE.,                & ! switch to use orography 
                                                ! coming from the PGD file
            LREAD_GROUND_PARAM = .TRUE.         ! switch to use soil parameters
                                                ! useful for the soil scheme
                                                ! coming from the PGD file

INTEGER           :: NSLEVE   =12         ! number of iteration for smooth orography
REAL              :: XSMOOTH_ZS = XUNDEF  ! optional uniform smooth orography for SLEVE coordinate
CHARACTER(LEN=28) :: YPGD_NAME, YPGD_DAD_NAME   ! general information
CHARACTER(LEN=8)  :: YKIND                      ! Kind of radiosounding data
CHARACTER(LEN=2)  :: YPGD_TYPE
!
INTEGER           :: IINFO_ll                   ! return code of // routines
TYPE(LIST_ll), POINTER :: TZ_FIELDS_ll           ! list of metric coefficient fields
!
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
INTEGER                           :: IBEG,IEND,IXOR,IXDIM,IYOR,IYDIM,ILBX,ILBY
REAL, DIMENSION(:),   ALLOCATABLE :: ZXHAT_ll, ZYHAT_ll
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTHL, ZRW, ZQSAT, ZRSAT, ZT
                                 ! variables to transform THETAL in THETA
!
INTEGER             :: ILENCH, IGRID, IRESP
CHARACTER (LEN=100) :: YCOMMENT
!
!JUAN TIMING
REAL,DIMENSION(2)         :: ZTIME1,ZTIME2,ZEND,ZTOT
CHARACTER                 :: YMI
INTEGER                   :: IMI
!JUAN TIMING
!
!*       0.2  Namelist declarations
!
NAMELIST/NAM_CONF_PRE/ LTHINSHELL,LCARTESIAN,    &! Declarations in MODD_CONF
                       LPACK,                    &!
                       NVERB,CIDEAL,CZS,         &!+global variables initialized
                       LBOUSS,LPERTURB,LPV_PERT, &! at their declarations
                       LRMV_BL,LFORCING,CEQNSYS   ! at their declarations
NAMELIST/NAM_GRID_PRE/ XLON0,XLAT0,            & ! Declarations in MODD_GRID
                       XBETA,XRPK,             & 
                       XLONORI,XLATORI
NAMELIST/NAM_GRIDH_PRE/ XLATCEN,XLONCEN,       & ! local variables  initialized
                 XDELTAX,XDELTAY,              & ! at their declarations
                 XHMAX,NEXPX,NEXPY,            &
                 XAX,XAY,NIZS,NJZS
NAMELIST/NAM_VPROF_PRE/LGEOSBAL, CFUNU,CFUNV,   &! global variables initialized
                     CTYPELOC,XLATLOC,XLONLOC,  &!  at their declarations
                     XXHATLOC,XYHATLOC,NILOC,NJLOC
NAMELIST/NAM_REAL_PGD/CPGD_FILE,                 & ! Physio-Graphic Data file
                                                   !  name
                      LREAD_ZS,                  & ! switch to use orography 
                                                   ! coming from the PGD file
                      LREAD_GROUND_PARAM
NAMELIST/NAM_SLEVE/NSLEVE, XSMOOTH_ZS
!
!*       0.3  Auxillary Namelist declarations
!
NAMELIST/NAM_DUST_PRE/ LDUST, LSALT
!
!-------------------------------------------------------------------------------
!
!*       0.    PROLOGUE
!              --------
!
CALL GOTO_MODEL(1)
!
CALL INITIO_ll()
NULLIFY(TZ_FIELDS_ll)
CALL VERSION
CPROGRAM='IDEAL '
!
!JUAN TIMING
  XT_START     = 0.0
  XT_STORE     = 0.0
!
  CALL SECOND_MNH2(ZEND)
!
!JUAN TIMING
!
!*       1.    INITIALIZE PHYSICAL CONSTANTS :         
!              ------------------------------
!
NVERB = 5
CALL INI_CST
!
!-------------------------------------------------------------------------------
!
!
!*  	 2.    SET DEFAULT VALUES  :  
!              --------------------
!
!
!*       2.1  For variables in DESFM file
!
CALL DEFAULT_DESFM_n(1)
!
CSURF = "NONE"
!
!
!*       2.2  For other global variables in EXPRE file
!
CALL DEFAULT_EXPRE
!-------------------------------------------------------------------------------
!
!*  	 3.    READ THE EXPRE FILE :  
!              --------------------
!
!*       3.1   initialize logical unit numbers (EXPRE and output-listing files)
!              and open these files :
! 
! 
CLUOUT  = 'OUTPUT_LISTING1'
CLUOUT0 = CLUOUT
CEXPRE  = 'PRE_IDEA1.nam'
CALL OPEN_ll(UNIT=NLUOUT,FILE=CLUOUT,IOSTAT=NRESP,FORM='FORMATTED',ACTION='WRITE', &
     MODE=GLOBAL)
CALL OPEN_ll(UNIT=NLUPRE,FILE=CEXPRE,IOSTAT=NRESP,ACTION='READ', &
     DELIM='QUOTE',MODE=GLOBAL)    
!
!*       3.2   read in NLUPRE the namelist informations
!
WRITE(NLUOUT,FMT=*) 'attempt to read ',TRIM(CEXPRE),' file'
CALL POSNAM(NLUPRE,'NAM_REAL_PGD',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_REAL_PGD)
!
!
CALL POSNAM(NLUPRE,'NAM_CONF_PRE',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_CONF_PRE)
CALL POSNAM(NLUPRE,'NAM_GRID_PRE',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_GRID_PRE)
CALL POSNAM(NLUPRE,'NAM_GRIDH_PRE',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_GRIDH_PRE)
CALL POSNAM(NLUPRE,'NAM_VPROF_PRE',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_VPROF_PRE)
CALL POSNAM(NLUPRE,'NAM_BLANK',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_BLANK)
CALL READ_PRE_IDEA_NAM_n(NLUPRE,NLUOUT)
CALL POSNAM(NLUPRE,'NAM_DUST_PRE',GFOUND,NLUOUT)
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_DUST_PRE)
!
!
IF( LEN_TRIM(CPGD_FILE) /= 0 ) THEN 
  ! open the PGD_FILE
  CALL FMOPEN_ll(CPGD_FILE,'READ',CLUOUT,NNPRAR,2,NVERB,NNINAR,NRESP)
  ! read the grid in the PGD file
  CALL FMREAD(CPGD_FILE,'IMAX',CLUOUT,'--',NIMAX,IGRID,ILENCH,YCOMMENT,IRESP)
  CALL FMREAD(CPGD_FILE,'JMAX',CLUOUT,'--',NJMAX,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
NIMAX_ll=NIMAX   !! _ll variables are global variables
NJMAX_ll=NJMAX   !! but the old names are kept in PRE_IDEA1.nam file
!
!*       3.3   check some parameters:
!
L1D=.FALSE. ; L2D=.FALSE.
!
IF ((NIMAX == 1).OR.(NJMAX == 1)) THEN 
  L2D=.TRUE.
  NJMAX_ll=1
  NIMAX_ll=MAX(NIMAX,NJMAX)
  WRITE(NLUOUT,FMT=*) ' NJMAX HAS BEEN SET TO 1 SINCE 2D INITIAL FILE IS REQUIRED &
                   & (L2D=TRUE) )' 
END IF
!
IF ((NIMAX == 1).AND.(NJMAX == 1)) THEN 
  L1D=.TRUE.
  NIMAX_ll = 1
  NJMAX_ll = 1
  WRITE(NLUOUT,FMT=*) ' 1D INITIAL FILE IS REQUIRED (L1D=TRUE) ' 
END IF
!
IF(.NOT. L1D) THEN
  LHORELAX_UVWTH=.TRUE.
  LHORELAX_RV=.TRUE.
ENDIF
!
NRIMX= MIN(JPRIMMAX,NIMAX_ll/2)
!
IF (L2D) THEN
  NRIMY=0
ELSE
  NRIMY= MIN(JPRIMMAX,NJMAX_ll/2)
END IF
!
IF (L1D) THEN
  NRIMX=0
  NRIMY=0
END IF
!
IF (L1D .AND. ( LPERTURB .OR. LGEOSBAL .OR.                &
               (.NOT. LCARTESIAN ) .OR. (.NOT. LTHINSHELL) ))THEN 
  LGEOSBAL   = .FALSE.
  LPERTURB   = .FALSE.
  LCARTESIAN = .TRUE. 
  LTHINSHELL = .TRUE. 
  WRITE(NLUOUT,FMT=*) ' LGEOSBAL AND LPERTURB HAVE BEEN SET TO FALSE &
                      & AND LCARTESIAN AND LTHINSHELL TO TRUE        &
                      & SINCE 1D INITIAL FILE IS REQUIRED (L1D=TRUE)' 
END IF
!
!*       3.4   compute the number of moist variables :
!
IF (.NOT.LUSERV) THEN
  LUSERV = .TRUE.
  WRITE(NLUOUT,FMT=*) ' LUSERV HAS BEEN RESET TO TRUE, SINCE A MOIST VARIABLE &
                   & IS PRESENT IN EXPRE FILE (CIDEAL = RSOU OR CSTN)' 
END IF
!
IF((LUSERI .OR. LUSERC).AND. (CIDEAL /= 'RSOU')) THEN
  WRITE(NLUOUT,FMT=*) 'USE OF HYDROMETEORS IS ONLY ALLOWED IN RSOU CASE'
  WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
  STOP
ENDIF
IF (LUSERI) THEN
  LUSERC =.TRUE.
  LUSERR =.TRUE.
  LUSERI =.TRUE.
  LUSERS =.TRUE.
  LUSERG =.TRUE.
  LUSERH =.FALSE.
  CCLOUD='ICE3'
ELSEIF(LUSERC) THEN
  LUSERR =.FALSE.
  LUSERI =.FALSE.
  LUSERS =.FALSE.
  LUSERG =.FALSE.
  LUSERH =.FALSE.
  CCLOUD='REVE'
ELSE
  LUSERC =.FALSE.
  LUSERR =.FALSE.
  LUSERI =.FALSE.
  LUSERS =.FALSE.
  LUSERG =.FALSE.
  LUSERH =.FALSE.
  LHORELAX_RC=.FALSE.
  LHORELAX_RR=.FALSE.
  LHORELAX_RI=.FALSE.
  LHORELAX_RS=.FALSE.
  LHORELAX_RG=.FALSE.
  LHORELAX_RH=.FALSE.
  CCLOUD='NONE'
!
END IF
!
NRR=0
IF (LUSERV) NRR=NRR+1
IF (LUSERC) NRR=NRR+1
IF (LUSERR) NRR=NRR+1
IF (LUSERI) NRR=NRR+1
IF (LUSERS) NRR=NRR+1
IF (LUSERG) NRR=NRR+1
IF (LUSERH) NRR=NRR+1
!
!
!*       3.5   Chemistry
!
IF (LORILAM .OR. LCH_INIT_FIELD) THEN
  ! Always initialize chemical scheme variables before INI_NSV call !
  CALL CH_INIT_SCHEME(NLUOUT)
  CALL CH_INIT_CCS(NEQ,NLUOUT,NVERB)
  LUSECHEM = .TRUE.
  IF (LORILAM) THEN
    CORGANIC = "MPMPO"
    LVARSIGI = .TRUE.
    LVARSIGJ = .TRUE.
    CALL CH_AER_INIT_SOA(NLUOUT, NVERB)
  END IF
END IF
! initialise NSV_* variables
CALL INI_NSV(1)
LHORELAX_SV(:)=.FALSE.
IF(.NOT. L1D) LHORELAX_SV(1:NSV)=.TRUE.
!
!-------------------------------------------------------------------------------
!
!*       4.    ALLOCATE MEMORY FOR ARRAYS :  
!   	       ----------------------------
!
!*       4.1  Vertical Spatial grid 
!
CALL READ_VER_GRID(CEXPRE)
!
!*       4.2  Initialize parallel variables and compute array's dimensions
!
!
IF(LGEOSBAL) THEN
  CALL SET_SPLITTING_ll('XSPLITTING')  ! required for integration of thermal wind balance
ELSE
  CALL SET_SPLITTING_ll('BSPLITTING')
ENDIF
CALL SET_JP_ll(1,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
CALL SET_FMPACK_ll(L1D,L2D,LPACK)
CALL SET_LBX_ll(CLBCX(1), 1)
CALL SET_LBY_ll(CLBCY(1), 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(NIMAX_ll+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(NJMAX_ll+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
!JUAN
! CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!JUAN
!
! sizes of arrays of the extended sub-domain
!
CALL GET_DIM_EXT_ll('B',NIU,NJU)
CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
CALL GET_INDICE_ll(NIB,NJB,NIE,NJE)
CALL GET_OR_ll('B',IXOR,IYOR)
NKB=1+JPVEXT
NKU=NKMAX+2*JPVEXT
!
!*       4.3  Global variables absent from the modules :
!
ALLOCATE(XJ(NIU,NJU,NKU))
SELECT CASE(CIDEAL)
  CASE('RSOU','CSTN')
    IF (LGEOSBAL) ALLOCATE(XCORIOZ(NIU,NJU,NKU))  ! exceptionally a 3D array  
  CASE DEFAULT                      ! undefined preinitialization
    WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : CIDEAL IS NOT CORRECTLY DEFINED'
    WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
    STOP
END SELECT 
!
!*       4.4   Prognostic variables at M instant (module MODD_FIELD1):
!
ALLOCATE(XUM(NIU,NJU,NKU))
ALLOCATE(XVM(NIU,NJU,NKU))
ALLOCATE(XWM(NIU,NJU,NKU))
ALLOCATE(XTHM(NIU,NJU,NKU))
ALLOCATE(XPABSM(NIU,NJU,NKU))
ALLOCATE(XRM(NIU,NJU,NKU,NRR))
ALLOCATE(XSVM(NIU,NJU,NKU,NSV))
!
!*       4.5   Grid variables (module MODD_GRID1 and MODD_METRICS1):
!
ALLOCATE(XMAP(NIU,NJU))
ALLOCATE(XLAT(NIU,NJU))
ALLOCATE(XLON(NIU,NJU))
ALLOCATE(XDXHAT(NIU),XDYHAT(NJU))
IF (LEN_TRIM(CPGD_FILE)==0) ALLOCATE(XZS(NIU,NJU))
IF (LEN_TRIM(CPGD_FILE)==0) ALLOCATE(XZSMT(NIU,NJU))
ALLOCATE(XZZ(NIU,NJU,NKU))
!
ALLOCATE(XDXX(NIU,NJU,NKU))
ALLOCATE(XDYY(NIU,NJU,NKU))
ALLOCATE(XDZX(NIU,NJU,NKU))
ALLOCATE(XDZY(NIU,NJU,NKU))
ALLOCATE(XDZZ(NIU,NJU,NKU))
!
!*       4.6   Reference state variables (modules MODD_REF and MODD_REF1):
!
ALLOCATE(XRHODREFZ(NKU),XTHVREFZ(NKU))
XTHVREFZ(:)=0.0
IF(CEQNSYS == 'DUR') THEN
  ALLOCATE(XRVREF(NIU,NJU,NKU))
ELSE
  ALLOCATE(XRVREF(0,0,0))
END IF
ALLOCATE(XRHODREF(NIU,NJU,NKU),XTHVREF(NIU,NJU,NKU),XEXNREF(NIU,NJU,NKU))
ALLOCATE(XRHODJ(NIU,NJU,NKU))
!
!*       4.7   Larger Scale fields (modules MODD_LSFIELD1):
!
ALLOCATE(XLSUM(NIU,NJU,NKU))
ALLOCATE(XLSVM(NIU,NJU,NKU))
ALLOCATE(XLSWM(NIU,NJU,NKU))
ALLOCATE(XLSTHM(NIU,NJU,NKU))
ALLOCATE(XLSRVM(NIU,NJU,NKU))
!
!  allocate lateral boundary field used for coupling
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
ELSEIF( L2D ) THEN             ! 2D case (not yet parallelized)
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
  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2
    NSIZELBXU_ll=2*NRIMX+2
    ALLOCATE(XLBXUM(2*NRIMX+2,NJU,NKU))
    ALLOCATE(XLBXVM(2*NRIMX+2,NJU,NKU))
    ALLOCATE(XLBXWM(2*NRIMX+2,NJU,NKU))
    ALLOCATE(XLBXTHM(2*NRIMX+2,NJU,NKU))
  ELSE
    NSIZELBX_ll=2
    NSIZELBXU_ll=4 
    ALLOCATE(XLBXUM(4,NJU,NKU))
    ALLOCATE(XLBXVM(2,NJU,NKU))
    ALLOCATE(XLBXWM(2,NJU,NKU))
    ALLOCATE(XLBXTHM(2,NJU,NKU))
  END IF  
  !
  IF ( NRR > 0 ) THEN
    IF (       LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
          .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH                     &
       ) THEN 
      NSIZELBXR_ll=2* NRIMX+2
      ALLOCATE(XLBXRM(2*NRIMX+2,NJU,NKU,NRR))
    ELSE
      NSIZELBXR_ll=2
      ALLOCATE(XLBXRM(2,NJU,NKU,NRR))
    ENDIF
  ELSE
    NSIZELBXR_ll=0
    ALLOCATE(XLBXRM(0,0,0,0))
  END IF
  !
  IF ( NSV > 0 ) THEN 
    IF ( ANY( LHORELAX_SV(:)) ) THEN
      NSIZELBXSV_ll=2* NRIMX+2
      ALLOCATE(XLBXSVM(2*NRIMX+2,NJU,NKU,NSV))
    ELSE
      NSIZELBXSV_ll=2
      ALLOCATE(XLBXSVM(2,NJU,NKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
  END IF
!
ELSE                                   ! 3D case
!
  CALL GET_SIZEX_LB(CLUOUT,NIMAX_ll,NJMAX_ll,NRIMX,   &
       IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU,         &
       IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  CALL GET_SIZEY_LB(CLUOUT,NIMAX_ll,NJMAX_ll,NRIMY,   &
       IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV,         &
       IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)
!
  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2
    NSIZELBXU_ll=2*NRIMX+2
    NSIZELBY_ll=2*NRIMY+2
    NSIZELBYV_ll=2*NRIMY+2
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,NKU))
    ALLOCATE(XLBYUM(IISIZEYF,IJSIZEYF,NKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,NKU))
    ALLOCATE(XLBYVM(IISIZEYFV,IJSIZEYFV,NKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,NKU))
    ALLOCATE(XLBYWM(IISIZEYF,IJSIZEYF,NKU))
    ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,NKU))
    ALLOCATE(XLBYTHM(IISIZEYF,IJSIZEYF,NKU))
  ELSE
    NSIZELBX_ll=2
    NSIZELBXU_ll=4
    NSIZELBY_ll=2
    NSIZELBYV_ll=4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,NKU))
    ALLOCATE(XLBYUM(IISIZEY2,IJSIZEY2,NKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,NKU))
    ALLOCATE(XLBYVM(IISIZEY4,IJSIZEY4,NKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,NKU))
    ALLOCATE(XLBYWM(IISIZEY2,IJSIZEY2,NKU))
    ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,NKU))
    ALLOCATE(XLBYTHM(IISIZEY2,IJSIZEY2,NKU))
  END IF  
  !
  IF ( NRR > 0 ) THEN
    IF (       LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
          .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH                     &
       ) THEN 
      NSIZELBXR_ll=2*NRIMX+2
      NSIZELBYR_ll=2*NRIMY+2
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,NKU,NRR))
      ALLOCATE(XLBYRM(IISIZEYF,IJSIZEYF,NKU,NRR))
    ELSE
      NSIZELBXR_ll=2
      NSIZELBYR_ll=2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,NKU,NRR))
      ALLOCATE(XLBYRM(IISIZEY2,IJSIZEY2,NKU,NRR))
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
      ALLOCATE(XLBXSVM(IISIZEXF,IJSIZEXF,NKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEYF,IJSIZEYF,NKU,NSV))
    ELSE
      NSIZELBXSV_ll=2
      NSIZELBYSV_ll=2
      ALLOCATE(XLBXSVM(IISIZEX2,IJSIZEX2,NKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEY2,IJSIZEY2,NKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    NSIZELBYSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
    ALLOCATE(XLBYSVM(0,0,0,0))
  END IF
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       5.     INITIALIZE ALL THE MODEL VARIABLES
!   	        ----------------------------------
!
!
!*       5.1    Grid variables and RS localization:
!
!*       5.1.1  Horizontal Spatial grid :
!
IF( LEN_TRIM(CPGD_FILE) /= 0 ) THEN 
!--------------------------------------------------------
! the MESONH horizontal grid will be read in the PGD_FILE 
!--------------------------------------------------------
  CALL READ_HGRID(1,CPGD_FILE,YPGD_NAME,YPGD_DAD_NAME,YPGD_TYPE)
! control the cartesian option
  IF( LCARTESIAN ) THEN
     WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : IN GENERAL, THE USE OF A PGD_FILE &
                & IMPLIES THAT YOU MUST TAKE INTO ACCOUNT THE EARTH SPHERICITY'
     WRITE(NLUOUT,FMT=*) 'NEVERTHELESS, LCARTESIAN HAS BEEN KEPT TO TRUE'
  END IF   
!
!* use of the externalized surface
!
  CSURF = "EXTE"
!
! determine whether the model is flat or no
!
  IF( ABS( MAXVAL(XZS(NIB:NIU-JPHEXT,NJB:NJU-JPHEXT)) ) < 1.E-10 ) THEN
    LFLAT=.TRUE.
  ELSE
    LFLAT=.FALSE.
  END IF
!

ELSE
!------------------------------------------------------------------------
! the MESONH horizontal grid is built from the PRE_IDEA1.nam informations
!------------------------------------------------------------------------
!
  ALLOCATE(XXHAT(NIU),XYHAT(NJU))
!
! define the grid localization at the earth surface by the central point
! coordinates
!
  IF (XLONCEN/=XUNDEF .OR. XLATCEN/=XUNDEF) THEN
    IF (XLONCEN/=XUNDEF .AND. XLATCEN/=XUNDEF) THEN 
!
! it should be noted that XLATCEN and XLONCEN refer to a vertical
! vorticity point and (XLATORI, XLONORI) refer to the mass point of
! conformal coordinates (0,0). This is to allow the centering of the model in
! a non-cyclic  configuration regarding to XLATCEN or XLONCEN.
!
      ALLOCATE(ZXHAT_ll(NIMAX_ll+2*JPHEXT),ZYHAT_ll(NJMAX_ll+2*JPHEXT))
      ZXHAT_ll=0.
      ZYHAT_ll=0.
      CALL SM_LATLON(XLATCEN,XLONCEN,                     &
                       -XDELTAX*(NIMAX_ll/2-0.5+JPHEXT),  &
                       -XDELTAY*(NJMAX_ll/2-0.5+JPHEXT),  &
                       XLATORI,XLONORI)
        DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
!
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : XLATORI=' , XLATORI, &
                          ' XLONORI= ', XLONORI
    ELSE
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : LATITUDE AND LONGITUDE OF THE CENTER &
                           & POINT MUST BE INITIALIZED ALL TOGETHER OR NOT'
      WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
      STOP
    END IF
  END IF
!
  IF (NPROC > 1) THEN
    CALL GET_DIM_EXT_ll('B',IXDIM,IYDIM)
    IBEG = IXOR-JPHEXT+1
    IEND = IBEG+IXDIM-1
    XXHAT(:) = (/ (FLOAT(JLOOP)*XDELTAX, JLOOP=IBEG,IEND) /)
    IBEG = IYOR-JPHEXT+1
    IEND = IBEG+IYDIM-1
    XYHAT(:) = (/ (FLOAT(JLOOP)*XDELTAY, JLOOP=IBEG,IEND) /)
!
  ELSE
    XXHAT(:) = (/ (FLOAT(JLOOP-NIB)*XDELTAX, JLOOP=1,NIU) /)
    XYHAT(:) = (/ (FLOAT(JLOOP-NJB)*XDELTAY, JLOOP=1,NJU) /)
  END IF
END IF
!
!*       5.1.2  Orography and Gal-Chen Sommerville transformation :
!
IF (    LEN_TRIM(CPGD_FILE) == 0  .OR. .NOT. LREAD_ZS) THEN
  SELECT CASE(CZS)     ! 'FLAT' or 'SINE' or 'BELL'
  CASE('FLAT')
    LFLAT = .TRUE.
    IF (XHMAX==XUNDEF) THEN
      XZS(:,:) = 0.
    ELSE
      XZS(:,:) = XHMAX
    END IF
  CASE('SINE')       ! sinus-shaped orography 
    IF (XHMAX==XUNDEF) XHMAX=300.
    LFLAT    =.FALSE.
    XZS(:,:) = XHMAX          &      ! three-dimensional case   
    *SPREAD((/((SIN((XPI/(NIMAX_ll+2*JPHEXT-1))*JLOOP)**2)**NEXPX,JLOOP=IXOR-1,IXOR+NIU-2)/),2,NJU) &
    *SPREAD((/((SIN((XPI/(NJMAX_ll+2*JPHEXT-1))*JLOOP)**2)**NEXPY,JLOOP=IYOR-1,IYOR+NJU-2)/),1,NIU)
    IF(L1D) THEN                     ! one-dimensional case
      XZS(:,:) = XHMAX 
    END IF        
  CASE('BELL')       ! bell-shaped orography 
    IF (XHMAX==XUNDEF) XHMAX=300.
    LFLAT = .FALSE.
    IF(.NOT.L2D) THEN                ! three-dimensional case
      XZS(:,:) = XHMAX  / ( 1.                                           &
        + ( (SPREAD(XXHAT(1:NIU),2,NJU) - FLOAT(NIZS) * XDELTAX) /XAX ) **2  &
        + ( (SPREAD(XYHAT(1:NJU),1,NIU) - FLOAT(NJZS) * XDELTAY) /XAY ) **2  ) **1.5
    ELSE                             ! two-dimensional case
      XZS(:,:) = XHMAX  / ( 1.                                          &
        + ( (SPREAD(XXHAT(1:NIU),2,NJU) - FLOAT(NIZS) * XDELTAX) /XAX ) **2 )
    ENDIF
    IF(L1D) THEN                     ! one-dimensional case
      XZS(:,:) = XHMAX 
    END IF        
  CASE('DATA')       ! discretized orography
    LFLAT    =.FALSE.
    WRITE(NLUOUT,FMT=*) 'CZS="DATA",   ATTEMPT TO READ ARRAY     &
                    &XZS(NIB:NIU-JPHEXT:1,NJU-JPHEXT:NJB:-1) &
                    &starting from the first index'
    CALL POSKEY(NLUPRE,NLUOUT,'ZS')
    DO JJLOOP = NJU-JPHEXT,NJB,-1    ! input like a map prior the sounding
      READ(NLUPRE,FMT=*) XZS(NIB:NIU-JPHEXT,JJLOOP)
    END DO
!
  CASE DEFAULT   ! undefined  shape of orography
    WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: ERRONEOUS TERRAIN TYPE'
    WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
    STOP
  END SELECT
!
  CALL ADD2DFIELD_ll(TZ_FIELDS_ll, XZS)
  CALL UPDATE_HALO_ll(TZ_FIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZ_FIELDS_ll)
!
END IF
!
IF (LWEST_ll())  THEN
  DO JILOOP = 1,JPHEXT
    XZS(JILOOP,:) = XZS(NIB,:)
  END DO
END IF
IF (LEAST_ll()) THEN
  DO JILOOP = NIU-JPHEXT+1,NIU
    XZS(JILOOP,:)=XZS(NIU-JPHEXT,:)
  END DO
END IF
IF (LSOUTH_ll()) THEN
  DO JJLOOP = 1,JPHEXT
    XZS(:,JJLOOP)=XZS(:,NJB)
  END DO
END IF
IF (LNORTH_ll()) THEN
  DO JJLOOP =NJU-JPHEXT+1,NJU
    XZS(:,JJLOOP)=XZS(:,NJU-JPHEXT)
  END DO
END IF
!
IF ( LEN_TRIM(CPGD_FILE) == 0  .OR. .NOT. LREAD_ZS) THEN
  IF (LSLEVE) THEN
    CALL ZSMT_PIC(NSLEVE,XSMOOTH_ZS)
  ELSE
    XZSMT(:,:) = 0.
  END IF
END IF
!
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(CLUOUT,XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XDXHAT,XDYHAT,XZZ,XJ)
  XMAP=1.
ELSE
  CALL SM_GRIDPROJ(CLUOUT,XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XLATORI,XLONORI, &
                   XMAP,XLAT,XLON,XDXHAT,XDYHAT,XZZ,XJ)
END IF
!
!*       5.1.3  Compute the localization in index space of the vertical profile
!               in CSTN and RSOU cases  :
!
IF (CTYPELOC =='LATLON' ) THEN  
  IF (.NOT.LCARTESIAN) THEN                            ! compute (x,y) if 
    CALL SM_XYHAT(XLATORI,XLONORI,                 &   ! the localization 
                  XLATLOC,XLONLOC,XXHATLOC,XYHATLOC)   ! is given in latitude 
  ELSE                                                 ! and longitude
    WRITE(NLUOUT,FMT=*) 'CTYPELOC CANNOT BE LATLON IN CARTESIAN GEOMETRY'
    WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
    STOP
  END IF 
END IF  
!
ALLOCATE(ZXHAT_ll(NIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(NJMAX_ll+2 * JPHEXT))
CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,NRESP) !//
CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,NRESP) !//
IF (CTYPELOC /= 'IJGRID') THEN                                               
  NILOC = MINLOC(ABS(XXHATLOC-ZXHAT_ll(:)))
  NJLOC = MINLOC(ABS(XYHATLOC-ZYHAT_ll(:)))
END IF
!
IF ( L1D .AND. (NILOC(1) /= NIB .OR. NJLOC(1) /= NJB) ) THEN
  NILOC = NIB
  NJLOC = NJB
  WRITE(NLUOUT,FMT=*) 'FOR 1D CONFIGURATION, THE RS INFORMATIONS ARE TAKEN AT &
                      & I=NIB AND J=NJB (CENTRAL VERTICAL)'
END IF
!
IF ( L2D .AND. ( NJLOC(1) /= NJB) ) THEN
  NJLOC = NJB
  WRITE(NLUOUT,FMT=*) 'FOR 2D CONFIGURATION, THE RS INFORMATIONS ARE TAKEN AT &
                      & J=NJB (CENTRAL PLANE)'
END IF
!
!*       5.2    Prognostic variables (not multiplied by  rhoJ) : u,v,w,theta,r
!               and 1D anelastic reference state
!
IF(LPV_PERT .AND. .NOT.(LGEOSBAL)) THEN
  WRITE(NLUOUT,FMT=*) 'FOR PV INVERSION, LGEOSBAL HAS TO BE TRUE'
  STOP
ENDIF
!
IF(LPV_PERT .AND. NPROC>1) THEN
    WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : THE USE OF A PV INVERSION HAS TO BE &
                        & PERFORMED WITH MONOPROCESSOR MODE'
    WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
    STOP
ENDIF
!
!*       5.2.1  Use a Radiosounding : CIDEAL='RSOU''
!
IF (CIDEAL == 'RSOU') THEN
  WRITE(NLUOUT,FMT=*) 'CIDEAL="RSOU", attempt to read DATE'
  CALL POSKEY(NLUPRE,NLUOUT,'RSOU')
  READ(NLUPRE,FMT=*)  NYEAR,NMONTH,NDAY,XTIME
  TDTCUR = DATE_TIME(DATE(NYEAR,NMONTH,NDAY),XTIME)
  TDTEXP = TDTCUR
  TDTSEG = TDTCUR
  TDTMOD = TDTCUR
  READ(NLUPRE,*) YKIND
  BACKSPACE(NLUPRE)    ! because YKIND read again in set_rsou
  WRITE(NLUOUT,FMT=*) 'CIDEAL="RSOU", ATTEMPT TO PROCESS THE SOUNDING DATA'
  IF (LGEOSBAL) THEN
    CALL SET_RSOU(CEXPRE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS,LPV_PERT,LRMV_BL,XCORIOZ)
  ELSE
    CALL SET_RSOU(CEXPRE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS,LPV_PERT,LRMV_BL)
  END IF 
!
!*       5.2.2  N=cste  and U(z) : CIDEAL='CSTN'
!
ELSE IF (CIDEAL == 'CSTN') THEN
  WRITE(NLUOUT,FMT=*) 'CIDEAL="CSTN", attempt to read DATE'
  CALL POSKEY(NLUPRE,NLUOUT,'CSTN')
  READ(NLUPRE,FMT=*)  NYEAR,NMONTH,NDAY,XTIME
  TDTCUR = DATE_TIME(DATE(NYEAR,NMONTH,NDAY),XTIME)
  TDTEXP = TDTCUR
  TDTSEG = TDTCUR
  TDTMOD = TDTCUR
  WRITE(NLUOUT,FMT=*) 'CIDEAL="CSTN", ATTEMPT TO PROCESS THE SOUNDING DATA'
  IF (LGEOSBAL) THEN
    CALL SET_CSTN(CEXPRE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS,LPV_PERT,LRMV_BL,XCORIOZ)
  ELSE
    CALL SET_CSTN(CEXPRE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS,LPV_PERT,LRMV_BL)
  END IF
!
END IF 
!
!*       5.3    Forcing variables
!
IF (LFORCING) THEN
  WRITE(NLUOUT,FMT=*) 'FORCING IS ENABLED, ATTEMPT TO SET FORCING FIELDS'
  CALL POSKEY(NLUPRE,NLUOUT,'ZFRC','PFRC')
  CALL SET_FRC(CEXPRE)
END IF
!
!*       5.4    3D Reference state variables :
!
!
!*       5.4.1  metrics coefficients and update halos:
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!*       5.4.2  3D reference state :
!
CALL SET_REF(0,'NIL',CLUOUT,                                     &
             XZZ,XZHAT,XJ,XDXX,XDYY,CLBCX,CLBCY,                 &
             XREFMASS,XMASS_O_PHI0,XLINMASS,                     &
             XRHODREF,XTHVREF,XRVREF,XEXNREF,XRHODJ              )
!
!
!*       5.5.1  Absolute pressure :
!
XPABSM= XP00 * XEXNREF**(XCPD/XRD)
!
!*       5.5.2  Total mass of dry air Md computation :
!
CALL TOTAL_DMASS(CLUOUT,XJ,XRHODREF,XDRYMASST)
!
!
!*       5.6    Complete prognostic variables (multipliy by  rhoJ) at time t :
!
! U grid   : gridpoint 2
IF (LWEST_ll())  XUM(1,:,:)    = 2.*XUM(2,:,:) - XUM(3,:,:)
! V grid   : gridpoint 3
!JUAN IF (LEAST_ll())  XVM(:,1,:)    = 2.*XVM(:,2,:) - XVM(:,3,:)
IF (LSOUTH_ll())  XVM(:,1,:)    = 2.*XVM(:,2,:) - XVM(:,3,:)
! SV : gridpoint 1
XSVM(:,:,:,:) = 0.
!
!
!*       5.7   Larger scale fields initialization :
!
XLSUM(:,:,:) = XUM(:,:,:)        ! these fields do not satisfy the 
XLSVM(:,:,:) = XVM(:,:,:)        ! lower boundary condition but are 
XLSWM(:,:,:) = XWM(:,:,:)        ! in equilibrium
XLSTHM(:,:,:)= XTHM(:,:,:)
XLSRVM(:,:,:)= XRM(:,:,:,1)
!
! enforce the vertical homogeneity under the ground and above the top of
! the model for the LS fields
!
XLSUM(:,:,NKB-1)=XLSUM(:,:,NKB)
XLSUM(:,:,NKU)=XLSUM(:,:,NKU-1)
XLSVM(:,:,NKB-1)=XLSVM(:,:,NKB)
XLSVM(:,:,NKU)=XLSVM(:,:,NKU-1)
XLSWM(:,:,NKB-1)=XLSWM(:,:,NKB)
XLSWM(:,:,NKU)=XLSWM(:,:,NKU-1)
XLSTHM(:,:,NKB-1)=XLSTHM(:,:,NKB)
XLSTHM(:,:,NKU)=XLSTHM(:,:,NKU-1)
IF ( NRR > 0 ) THEN
  XLSRVM(:,:,NKB-1)=XLSRVM(:,:,NKB)
  XLSRVM(:,:,NKU)=XLSRVM(:,:,NKU-1)
END IF
!
ILBX=SIZE(XLBXUM,1)
ILBY=SIZE(XLBYUM,2)
IF(LWEST_ll() .AND. .NOT. L1D) THEN
  XLBXUM(1:NRIMX+1,        :,:)     = XUM(2:NRIMX+2,        :,:)
  XLBXVM(1:NRIMX+1,        :,:)     = XVM(1:NRIMX+1,        :,:)
  XLBXWM(1:NRIMX+1,        :,:)     = XWM(1:NRIMX+1,        :,:)
  XLBXTHM(1:NRIMX+1,        :,:)   = XTHM(1:NRIMX+1,        :,:)
  XLBXRM(1:NRIMX+1,        :,:,:)   = XRM(1:NRIMX+1,        :,:,:)
ENDIF
IF(LEAST_ll() .AND. .NOT. L1D) THEN
  XLBXUM(ILBX-NRIMX:ILBX,:,:)     = XUM(NIU-NRIMX:NIU,    :,:)
  XLBXVM(ILBX-NRIMX:ILBX,:,:)     = XVM(NIU-NRIMX:NIU,    :,:)
  XLBXWM(ILBX-NRIMX:ILBX,:,:)     = XWM(NIU-NRIMX:NIU,    :,:)
  XLBXTHM(ILBX-NRIMX:ILBX,:,:)   = XTHM(NIU-NRIMX:NIU,    :,:)
  XLBXRM(ILBX-NRIMX:ILBX,:,:,:)   = XRM(NIU-NRIMX:NIU,    :,:,:)
ENDIF
IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYUM(:,1:NRIMY+1,        :)     = XUM(:,1:NRIMY+1,      :)
  XLBYVM(:,1:NRIMY+1,        :)     = XVM(:,2:NRIMY+2,      :)
  XLBYWM(:,1:NRIMY+1,        :)     = XWM(:,1:NRIMY+1,  :)
  XLBYTHM(:,1:NRIMY+1,        :)    = XTHM(:,1:NRIMY+1,      :)
  XLBYRM(:,1:NRIMY+1,        :,:)   = XRM(:,1:NRIMY+1,      :,:)
ENDIF
IF(LNORTH_ll().AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYUM(:,ILBY-NRIMY:ILBY,:)     = XUM(:,NJU-NRIMY:NJU,  :)
  XLBYVM(:,ILBY-NRIMY:ILBY,:)     = XVM(:,NJU-NRIMY:NJU,  :)
  XLBYWM(:,ILBY-NRIMY:ILBY,:)     = XWM(:,NJU-NRIMY:NJU,  :)
  XLBYTHM(:,ILBY-NRIMY:ILBY,:)    = XTHM(:,NJU-NRIMY:NJU,  :)
  XLBYRM(:,ILBY-NRIMY:ILBY,:,:)   = XRM(:,NJU-NRIMY:NJU,  :,:)
ENDIF
DO JSV = 1, NSV
  IF(LWEST_ll() .AND. .NOT. L1D) &
  XLBXSVM(1:NRIMX+1,        :,:,JSV)   = XSVM(1:NRIMX+1,        :,:,JSV)
  IF(LEAST_ll() .AND. .NOT. L1D) &
  XLBXSVM(ILBX-NRIMX:ILBX,:,:,JSV)   = XSVM(NIU-NRIMX:NIU,    :,:,JSV)
  IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
  XLBYSVM(:,1:NRIMY+1,        :,JSV)   = XSVM(:,1:NRIMY+1,      :,JSV)
  IF(LNORTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
  XLBYSVM(:,ILBY-NRIMY:ILBY,:,JSV)   = XSVM(:,NJU-NRIMY:NJU,  :,JSV)
END DO
!
!
!*       5.8   Add a perturbation to a basic state :
!
IF(LPERTURB) CALL SET_PERTURB(CEXPRE)
!
!
!*       5.9   Anelastic correction and pressure:
!
IF ( .NOT. L1D ) CALL PRESSURE_IN_PREP(XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!
!*       5.10  Compute THETA, vapor and cloud mixing ratio
!
IF (CIDEAL == 'RSOU') THEN
  IF(YKIND == 'ZUVTHLMR') THEN
    ALLOCATE(ZTHL (NIU,NJU,NKU))
    ALLOCATE(ZRW  (NIU,NJU,NKU))
    ALLOCATE(ZQSAT(NIU,NJU,NKU))
    ALLOCATE(ZRSAT(NIU,NJU,NKU))
    ALLOCATE(ZT   (NIU,NJU,NKU))
    ZTHL=XTHM
    ZRW(:,:,:)=XRM(:,:,:,1)
!
    XRM(:,:,:,2)=0.
!
    DO JLOOP=1,20
      ZT=XTHM*(XPABSM/XP00)**(XRD/XCPD)
      ZQSAT=QSAT(ZT,XPABSM)
      ZRSAT=1./(1./ZQSAT(:,:,:)-1.)
      WHERE(ZRW(:,:,:)>ZRSAT(:,:,:))
        XRM(:,:,:,2) = ZRW(:,:,:) - ZRSAT(:,:,:)
        XRM(:,:,:,1) = ZRSAT(:,:,:)
      ELSEWHERE
        XRM(:,:,:,2) = 0.
        XRM(:,:,:,1) = ZRW(:,:,:)
      END WHERE
      XTHM(:,:,:) = 0.5 * XTHM(:,:,:) &
                   +0.5 *(ZTHL(:,:,:) + XLVTT*XTHM(:,:,:)/XCPD/ZT(:,:,:)*XRM(:,:,:,2))
    END DO
    DEALLOCATE(ZTHL )
    DEALLOCATE(ZRW  )
    DEALLOCATE(ZQSAT)
    DEALLOCATE(ZRSAT)
    DEALLOCATE(ZT   )
!
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*  	 6.    INITIALIZE SCALAR VARIABLES FOR CHEMISTRY
!   	       -----------------------------------------
!
!  before calling chemistry
CCONF = 'START'
CSTORAGE_TYPE='TT'                  ! instant t and t-dt are the same
CALL CLOSE_ll(CEXPRE,IOSTAT=NRESP)  ! Close the EXPRE file 
!
IF ( LCH_INIT_FIELD ) CALL CH_INIT_FIELD_n(1, NLUOUT, NVERB)
!
!-------------------------------------------------------------------------------
!
!*   	 7.    WRITE THE FMFILE 
!   	       ----------------
!
CALL SECOND_MNH2(ZTIME1)
!
NNPRAR = 22 + 2*(NRR+NSV)   &    ! 22 = number of grid variables + reference 
       + 8 + 17                  ! state variables + dimension variables
                                 ! 2*(8+NRR+NSV) + 1 = number of prognostic
                                 ! variables at time t and t-dt
NTYPE=1
CDESFM=ADJUSTL(ADJUSTR(CINIFILE)//'.des')
!
CALL FMOPEN_ll(CINIFILE,'WRITE',CLUOUT,NNPRAR,NTYPE,NVERB,NNINAR,NRESP)
!
CALL WRITE_DESFM_n(1,CDESFM,CLUOUT)
!
CALL WRITE_LFIFM_n(CINIFILE,'                            ')  ! There is no DAD model for PREP_IDEAL_CASE 
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STORE = XT_STORE + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*     8.     EXTERNALIZED SURFACE
!             --------------------
!
!
IF (CSURF =='EXTE') THEN
  ! Switch to model 1 surface variables
  CALL GOTO_SURFEX(1)
  !* definition of physiographic fields
  ! computed ...
  IF (LEN_TRIM(CPGD_FILE)==0 .OR. .NOT. LREAD_GROUND_PARAM) THEN
    CPGDFILE = CINIFILE
    CALL PGD_SURF_ATM('MESONH',CINIFILE,'MESONH',.TRUE.,.TRUE.)
  ELSE
  ! ... or read from file.
    CPGDFILE = CPGD_FILE
    CALL INIT_PGD_SURF_ATM('MESONH','PGD',                         &
                            '                            ','      ',&
                            TDTCUR%TDATE%YEAR, TDTCUR%TDATE%MONTH,  &
                            TDTCUR%TDATE%DAY, TDTCUR%TIME           )
!
  END IF
  !
  !* forces orography from atmospheric file
  IF (.NOT. LREAD_ZS) CALL MNHPUT_ZS_n
  !
  !* writing of physiographic fields in the file
  COUTFMFILE = CINIFILE
  CALL WRITE_SURF_ATM_n('MESONH','PGD')
  !
  !* deallocation of physiographic fields
  CALL DEALLOC_SURF_ATM_n
  !
  !* rereading of physiographic fields and definition of prognostic fields
  !* writing of all surface fields
  CPGDFILE   = CINIFILE
  COUTFMFILE = CINIFILE
  CALL PREP_SURF_MNH('                            ','      ')
ELSE
  CSURF = "NONE"
END IF
!
!-------------------------------------------------------------------------------
!
!*     9.     CLOSES THE FILE
!             ---------------
!
CALL FMCLOS_ll(CINIFILE,'KEEP',CLUOUT,NRESP)
IF( LEN_TRIM(CPGD_FILE) /= 0 ) THEN
  CALL FMCLOS_ll(CPGD_FILE,'KEEP',CLUOUT,NRESP)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*      10.    PRINTS ON OUTPUT-LISTING
!              ------------------------
!
IF (NVERB >= 5) THEN
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: LCARTESIAN,CIDEAL,CZS=', &
                                    LCARTESIAN,CIDEAL,CZS 
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: LUSERV=',LUSERV
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: XLON0,XLAT0,XBETA,XRPK,XLONORI,XLATORI=', &
                                    XLON0,XLAT0,XBETA,XRPK,XLONORI,XLATORI
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: XDELTAX,XDELTAY=',XDELTAX,XDELTAY
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: NVERB=',NVERB
  IF(LCARTESIAN) THEN
    WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: No map projection used.'
  ELSE
    IF (XRPK == 1.) THEN
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: Polar stereo used.'
    ELSE IF (XRPK == 0.) THEN
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: Mercator used.'
    ELSE
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: Lambert used, cone factor=',XRPK 
    END IF
  END IF
END IF
!
IF (NVERB >= 5) THEN
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: IIB, IJB, IKB=',NIB,NJB,NKB
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: IIU, IJU, IKU=',NIU,NJU,NKU
END IF
!
!
!*       28.1   print statistics!
!
  !
  CALL SECOND_MNH2(ZTIME2)
  XT_START=XT_START+ZTIME2-ZEND
  !
  ! Set File Timing OUTPUT
  !
  CALL SET_ILUOUT_TIMING(NLUOUT)
  !
  ! Compute global time
  !
  CALL TIME_STAT_ll(XT_START,ZTOT)
  !
  !
  IMI = 1
  CALL TIME_HEADER_ll(IMI)
  !
  CALL TIME_STAT_ll(XT_STORE,ZTOT,      ' STORE-FIELDS','=')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')  
  WRITE(YMI,FMT="(I0)") IMI
  CALL TIME_STAT_ll(XT_START,ZTOT,      ' MODEL'//YMI,'+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
WRITE(NLUOUT,FMT=*) ' '
WRITE(NLUOUT,FMT=*) '****************************************************'
WRITE(NLUOUT,FMT=*) '* PREP_IDEAL_CASE: PREP_IDEAL_CASE ENDS CORRECTLY. *'
WRITE(NLUOUT,FMT=*) '****************************************************'
!
CALL CLOSE_ll(CLUOUT,IOSTAT=NRESP)
CALL END_PARA_ll(IINFO_ll)
!
!
STOP
!
END PROGRAM PREP_IDEAL_CASE
