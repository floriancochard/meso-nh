!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##################
      MODULE MODD_BUDGET
!     ##################
!
!!****  *MODD_BUDGET* - declaration of budget variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the budget 
!     variables.     
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_PARAMETERS: JPBUMAX, JPBUPROCMAX
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_BUDGET)
!!          
!!    AUTHOR
!!    ------
!!	P. Hereil   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        23/02/95 
!!      J.-P. Lafore    10/02/98    adding of rhodj declaration for budget  
!!      V. Ducrocq      4/06/99     //
!!      J.-P. Pinty     25/09/00    additional budget terms for C2R2 scheme
!!      D. Gazen        22/01/01    add NCHEMSV
!!      V. Masson       06/11/02    new flags for budget calls and time counters
!!      V. Masson       27/11/02    add 2way nesting effect
!!      P. Jabouille    07/07/04    add budget terms for microphysics
!!      C. Barthe       19/11/09    add budget terms for electricity          
!!      C.Lac           04/2016  negative contribution to the budget splitted between advection, turbulence and microphysics for KHKO/C2R2
!!      C. Barthe            /16    add budget terms for LIMA
!!      C. LAc          10/2016 add droplets deposition
!!      S. Riette       11/2016  New budgets for ICE3/ICE4
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
USE MODD_PARAMETERS, ONLY :JPBUMAX, JPBUPROMAX, NMNHNAMELGTMAX
!
IMPLICIT NONE
!
!                       General variables
LOGICAL, SAVE :: LBU_ENABLE
!
INTEGER, SAVE, DIMENSION(JPBUMAX,JPBUPROMAX)  &  ! number of processes to be
                             :: NBUINC           ! avoided for every budget
                                                 ! between one active
                                                 ! source to the next one
INTEGER, SAVE, DIMENSION(JPBUMAX)             &  ! counter for all the processes
                             :: NBUCTR_ACTV      ! activated or not
!
CHARACTER (LEN=4), SAVE :: CBUTYPE         ! type of desired budget 'CART'
                                           ! (cartesian box) or 'MASK' (budget
                                           ! zone defined by a mask) or 'NONE'
                                           ! (no budget)
INTEGER, SAVE :: NBUMOD                    ! model in which budget is 
                                           ! calculated
INTEGER, SAVE, DIMENSION(:),             & ! number of processes for each 
                 ALLOCATABLE :: NBUPROCNBR ! budget 
!
INTEGER, SAVE, DIMENSION(:),             & ! process counter linked to each 
                 ALLOCATABLE :: NBUPROCCTR ! budget 
!
CHARACTER(LEN=2), SAVE, DIMENSION(:,:),  & ! resulting string character of the 
        ALLOCATABLE :: CBUACTION           ! transcription of the budget actions 
                                           ! (integer) read in  namelists or 
                                           ! set by default
CHARACTER (LEN=NMNHNAMELGTMAX), SAVE, DIMENSION(:,:),& ! names of records on the FM file 
                 ALLOCATABLE :: CBURECORD  ! for the budgets 
!
CHARACTER (LEN=99), SAVE, DIMENSION(:,:),& ! name of a process for a budget. It
                 ALLOCATABLE :: CBUCOMMENT ! will appear in the comment part of 
                                           ! the previous record
!
LOGICAL, SAVE :: LBU_BEG                   ! switch for budget beginning
!
REAL, SAVE    :: XBULEN                    ! length in seconds of the budget 
                                           ! temporal average
!
INTEGER, SAVE :: NBUSTEP                   ! number of model timesteps required 
                                           ! for the budget time average
REAL, SAVE    :: XBUWRI                       ! period in seconds of
                                           ! budget writing on FM-files
INTEGER, SAVE :: NBUWRNB                   ! number of budget periods when storage
                                           ! arrays are written on FM-files
INTEGER, SAVE :: NBUTSHIFT                 ! temporal shift for budgets writing
!
INTEGER, SAVE :: NBUKL, NBUKH              ! lowest and highest K indice values 
                                           ! of the budget box 
LOGICAL, SAVE :: LBU_KCP                   ! switch for compression in K
                                           ! direction
!
!                Variables used by the cartesian box case ('CART') only
!
INTEGER, SAVE :: NBUIL, NBUIH              ! lowest and highest I indice values 
                                           ! of the cartesian box 
INTEGER, SAVE :: NBUJL, NBUJH              ! lowest and highest J indice values 
                                           ! of the cartesian box 
LOGICAL, SAVE :: LBU_ICP                   ! switch for compression in I
                                           ! direction
LOGICAL, SAVE :: LBU_JCP                   ! switch for comppression in J
                                           ! direction
!
!                Variables used by the  mask case ('MASK') only
!
INTEGER, SAVE :: NBUMASK                   ! number of MASK zones for which 
                                           ! budgets are performed 
LOGICAL, SAVE, DIMENSION(:,:,:),         & ! define the zone where the MASK 
           ALLOCATABLE :: LBU_MASK         ! is True 
!                                          
REAL, SAVE, DIMENSION(:,:,:,:),          & ! surface for each mask at each   
           ALLOCATABLE :: XBUSURF          ! budget step   
!             
INTEGER, SAVE :: NBUTIME                   ! number of budget time periods
!
!                       Variables for budget storage 
!
!                       General variables
INTEGER, SAVE :: NBUSIL, NBUSIH      ! lowest and highest I indices of the intersection
                                     ! of the cartesian box with the sub-domain
INTEGER, SAVE :: NBUSJL, NBUSJH      ! lowest and highest J indices of the intersection
                                     ! of the global cartesian box
INTEGER, SAVE :: NBUIMAX_ll                ! second dimension of the budget
INTEGER, SAVE :: NBUJMAX_ll                ! second dimension of the budget
                                           ! array in the global domain (in CART case)
!
INTEGER, SAVE :: NBUIMAX                   ! first dimension of the budget
                                           ! tabular
INTEGER, SAVE :: NBUJMAX                   ! second dimension of the budget
                                           ! tabular
INTEGER, SAVE :: NBUKMAX                   ! dimension along K of the budget
                                           ! tabular
REAL, SAVE, DIMENSION(:,:,:,:),          & ! budget arrays for RU, RV and
        ALLOCATABLE :: XBURU, XBURV, XBURW ! RW (wind components) respectively
REAL, SAVE, DIMENSION(:,:,:,:),          & ! budget arrays for RTH (potential 
        ALLOCATABLE :: XBURTH, XBURTKE     ! temperature) and RTKE (kinetic
                                           ! energy)
REAL, SAVE, DIMENSION(:,:,:,:),          & ! budget arrays for RRV (water vapor)
        ALLOCATABLE :: XBURRV, XBURRC      ! and RRC (cloud water)
REAL, SAVE, DIMENSION(:,:,:,:),          & ! budget arrays for RRR (rain water)
        ALLOCATABLE :: XBURRR, XBURRI      ! and RRI (ice)
REAL, SAVE, DIMENSION(:,:,:,:),          & ! budget arrays for RRS (snow)
        ALLOCATABLE :: XBURRS, XBURRG      ! and RRG (graupel)
REAL, SAVE, DIMENSION(:,:,:,:),          & ! budget array for RRH (hail)
        ALLOCATABLE :: XBURRH              ! 
REAL, SAVE, DIMENSION(:,:,:,:,:), &
                 ALLOCATABLE :: XBURSV       ! Budget of the SVx
REAL, SAVE, DIMENSION(:,:,:),            & ! budget arrays for RHODJ at
               ALLOCATABLE :: XBURHODJ , & !   scalar localization
                              XBURHODJU, & !        U localization
                              XBURHODJV, & !        V localization
                              XBURHODJW    !    and W localization
!
!      Allowed processes for the budget of the x scalar variables
!        (transport part only)
!
! For each budget, the switches values (from 0 to JPBUPROMAX) for budgets 
! activation may be set by the user in a namelist. Their default value is 0.
! In the following declaration, the corresponding process names  are given 
! beside as comments.
!     
!      Allowed processes for the budget of RU (wind component along x)
!
! Courant namelist: NAM_BURU
!
LOGICAL, SAVE :: LBU_RU     ! True when the budget of RU is performed
!                         
INTEGER, SAVE :: NASSEU     ! time filter
INTEGER, SAVE :: NNESTU     ! Efffect of 2way nesting on U
INTEGER, SAVE :: NADVU      ! advection 
INTEGER, SAVE :: NFRCU      ! forcing
INTEGER, SAVE :: NNUDU      ! nudging
INTEGER, SAVE :: NCURVU     ! curvature
INTEGER, SAVE :: NCORU      ! Coriolis terms 
INTEGER, SAVE :: NDIFU      ! numerical diffusion
INTEGER, SAVE :: NRELU      ! relaxation
INTEGER, SAVE :: NHTURBU    ! horizontal TURBulence
INTEGER, SAVE :: NVTURBU    ! vertical turbulence 
INTEGER, SAVE :: NDRAGU     ! vegetation drag        
INTEGER, SAVE :: NMAFLU     ! mass flux            
INTEGER, SAVE :: NPRESU     ! pressure term
!     
!      Allowed processes for the budget of RV (wind component along y)
!                                                  
! Courant namelist: NAM_BURV
!
LOGICAL, SAVE :: LBU_RV     ! True when the budget of RV is performed
!
INTEGER, SAVE :: NASSEV     ! time filter
INTEGER, SAVE :: NNESTV     ! Efffect of 2way nesting on V
INTEGER, SAVE :: NADVV      ! advection 
INTEGER, SAVE :: NFRCV      ! forcing
INTEGER, SAVE :: NNUDV      ! nudging
INTEGER, SAVE :: NCURVV     ! curvature
INTEGER, SAVE :: NCORV      ! Coriolis terms 
INTEGER, SAVE :: NDIFV      ! numerical diffusion
INTEGER, SAVE :: NRELV      ! relaxation
INTEGER, SAVE :: NHTURBV    ! horizontal turbulence
INTEGER, SAVE :: NVTURBV    ! vertical turbulence 
INTEGER, SAVE :: NDRAGV     ! vegetation drag         
INTEGER, SAVE :: NMAFLV     ! mass flux            
INTEGER, SAVE :: NPRESV     ! pressure term
!
!      Allowed processes for the budget of RW (wind vertical component)
!                                                  
! Courant namelist: NAM_BURW
!
LOGICAL, SAVE :: LBU_RW     ! True when the budget of RW is performed 
!                                                  
INTEGER, SAVE :: NASSEW     ! time filter
INTEGER, SAVE :: NNESTW     ! Efffect of 2way nesting on W
INTEGER, SAVE :: NADVW      ! advection
INTEGER, SAVE :: NFRCW      ! forcing
INTEGER, SAVE :: NNUDW      ! nudging
INTEGER, SAVE :: NCURVW     ! curvature
INTEGER, SAVE :: NCORW      ! Coriolis terms 
INTEGER, SAVE :: NGRAVW     ! gravity term
INTEGER, SAVE :: NDIFW      ! numerical diffusion
INTEGER, SAVE :: NRELW      ! relaxation
INTEGER, SAVE :: NHTURBW    ! horizontal turbulence 
INTEGER, SAVE :: NVTURBW    ! vertical turbulence 
INTEGER, SAVE :: NPRESW     ! pressure term
!
!      Allowed processes for the budget of RTH (potential temperature)
!                                                  
! Courant namelist: NAM_BURTH
!
LOGICAL, SAVE :: LBU_RTH    ! True when the budget of RTH is performed
!
INTEGER, SAVE :: NASSETH    ! time filter
INTEGER, SAVE :: NNESTTH    ! Efffect of 2way nesting on Th
INTEGER, SAVE :: NADVTH     ! Total advection for PPM
INTEGER, SAVE :: NFRCTH     ! forcing
INTEGER, SAVE :: N2DADVTH   ! 2d advecting forcing
INTEGER, SAVE :: N2DRELTH   ! 2d relaxation forcing
INTEGER, SAVE :: NNUDTH     ! nudging
INTEGER, SAVE :: NPREFTH    ! theta source term due to the reference pressure
                            ! (Dyn. Sources) only present if KRR>0
INTEGER, SAVE :: NDIFTH     ! numerical diffusion
INTEGER, SAVE :: NRELTH     ! relaxation
INTEGER, SAVE :: NRADTH     ! RADiation
INTEGER, SAVE :: NDCONVTH   ! KAFR CONVection
INTEGER, SAVE :: NMAFLTH    ! Mass flux              
INTEGER, SAVE :: NHTURBTH   ! horizontal turbulence
INTEGER, SAVE :: NVTURBTH   ! vertical turbulence
INTEGER, SAVE :: NDISSHTH   ! dissipative heating
INTEGER, SAVE :: NNEGATH    ! negative correction induced by hydrometeors
INTEGER, SAVE :: NNETURTH    ! negative correction induced by hydrometeors
INTEGER, SAVE :: NNEADVTH    ! negative correction induced by hydrometeors
INTEGER, SAVE :: NNECONTH    ! negative correction induced by hydrometeors
INTEGER, SAVE :: NREVATH    ! rain evaporation
INTEGER, SAVE :: NCONDTH    ! evaporation/condensation
INTEGER, SAVE :: NHENUTH    ! HEterogenous NUcleation ICE3
INTEGER, SAVE :: NHONTH     ! HOmogeneous Nucleation  ICE3
INTEGER, SAVE :: NSFRTH     ! Spontaneous FReezing    ICE3
INTEGER, SAVE :: NDEPSTH    ! DEPosition on Snow      ICE3
INTEGER, SAVE :: NDEPGTH    ! DEPosition on Graupel   ICE3
INTEGER, SAVE :: NRIMTH     ! RIMing of cloudwater    ICE3
INTEGER, SAVE :: NACCTH     ! ACCretion of rainwater  ICE3
INTEGER, SAVE :: NCFRZTH    ! Conversion FReeZing     ICE3
INTEGER, SAVE :: NWETGTH    ! WET Growth of graupel   ICE3
INTEGER, SAVE :: NDRYGTH    ! DRY Growth of graupel   ICE3
INTEGER, SAVE :: NGMLTTH    ! Graupel MeLTing         ICE3
INTEGER, SAVE :: NIMLTTH    ! Ice MeLTing             ICE3
INTEGER, SAVE :: NBERFITH   ! BERgeron-FIndeisen gth. ICE3
INTEGER, SAVE :: NCDEPITH   ! Cond./DEPosition on ice ICE3
INTEGER, SAVE :: NWETHTH    ! wet growth of hail      ICE4
INTEGER, SAVE :: NDRYHTH    ! dry growth of hail      ICE4
INTEGER, SAVE :: NHMLTTH    ! melting of hail         ICE4
INTEGER, SAVE :: NADJUTH    ! adjustement before rain_ice ICE3
INTEGER, SAVE :: NCORRTH    ! tendencies correction after ICE3
INTEGER, SAVE :: NHINDTH    ! Heterogeneous Nucleation by Deposition LIMA
INTEGER, SAVE :: NHINCTH    ! Heterogeneous Nucleation by Contact    LIMA
INTEGER, SAVE :: NHONHTH    ! Haze Homogeneous Nucleation            LIMA
INTEGER, SAVE :: NHONCTH    ! droplet homogeneous nucleation         LIMA
INTEGER, SAVE :: NHONRTH    ! drop homogeneous nucleation            LIMA
INTEGER, SAVE :: NCEDSTH    ! adjustment
INTEGER, SAVE :: NSEDITH    ! Temperature transport by hydrometeors sedimentation
!
!      Allowed processes for the budget of RTKE (kinetic energy)
!                                                  
! Courant namelist: NAM_BURTKE
!
LOGICAL, SAVE :: LBU_RTKE   ! True when the budget of RTKE is performed
!
INTEGER, SAVE :: NASSETKE   ! time filter
INTEGER, SAVE :: NADVTKE    ! Total advection for PPM
INTEGER, SAVE :: NFRCTKE    ! forcing
INTEGER, SAVE :: NDIFTKE    ! numerical diffusion
INTEGER, SAVE :: NRELTKE    ! relaxation
INTEGER, SAVE :: NDPTKE     ! dynamic production of TKE
INTEGER, SAVE :: NTPTKE     ! thermal production of TKE
INTEGER, SAVE :: NDRAGTKE   ! vegetation drag              
INTEGER, SAVE :: NDISSTKE   ! dissipation of TKE
INTEGER, SAVE :: NTRTKE     ! turbulent transport of TKE
!
!
!      Allowed processes for the budget of moist variable RRV (water vapor)
!                                                  
! Courant namelist: NAM_BURRV
!
LOGICAL, SAVE :: LBU_RRV   ! true when the budget of RRV is performed
!
INTEGER, SAVE :: NASSERV   ! time filter
INTEGER, SAVE :: NNESTRV   ! Effect of 2way nesting on Rv
INTEGER, SAVE :: NADVRV    ! Total advection for PPM
INTEGER, SAVE :: NFRCRV    ! forcing
INTEGER, SAVE :: N2DADVRV  ! 2d advecting forcing
INTEGER, SAVE :: N2DRELRV  ! 2d relaxation forcing
INTEGER, SAVE :: NNUDRV    ! nudging
INTEGER, SAVE :: NDIFRV    ! numerical diffusion
INTEGER, SAVE :: NRELRV    ! relaxation
INTEGER, SAVE :: NDCONVRV  ! KAFR CONVection
INTEGER, SAVE :: NMAFLRV   ! Mass flux           
INTEGER, SAVE :: NHTURBRV  ! horizontal turbulence 
INTEGER, SAVE :: NVTURBRV  ! vertical turbulence
INTEGER, SAVE :: NNEGARV   ! negative correction                            
INTEGER, SAVE :: NNETURRV   ! negative correction                            
INTEGER, SAVE :: NNECONRV   ! negative correction                            
INTEGER, SAVE :: NNEADVRV   ! negative correction                            
INTEGER, SAVE :: NREVARV   ! rain evaporation
INTEGER, SAVE :: NCONDRV   ! evaporation/condensation
INTEGER, SAVE :: NHENURV   ! HEterogenous NUcleation ICE3
INTEGER, SAVE :: NDEPSRV   ! DEPosition on Snow      ICE3
INTEGER, SAVE :: NDEPGRV   ! DEPosition on Graupel   ICE3
INTEGER, SAVE :: NCDEPIRV  ! Cond./DEPosition on ice ICE3
INTEGER, SAVE :: NADJURV   ! adjustement before rain_ice ICE3
INTEGER, SAVE :: NCORRRV    ! tendencies correction after ICE3
INTEGER, SAVE :: NHINDRV   ! Heterogeneous Nucleation by Deposition LIMA
INTEGER, SAVE :: NHONHRV   ! Haze Homogeneous Nucleation            LIMA
INTEGER, SAVE :: NCEDSRV   ! adjustement 
!
!      Allowed processes for the budget of moist variable RRC (cloud water)
!                                                  
! Courant namelist: NAM_BURRC
!
LOGICAL, SAVE :: LBU_RRC    ! True when the budget of RRC is performed
!
INTEGER, SAVE :: NASSERC    ! time filter
INTEGER, SAVE :: NNESTRC    ! Efffect of 2way nesting on Rc
INTEGER, SAVE :: NADVRC     ! Total advection for PPM
INTEGER, SAVE :: NFRCRC     ! forcing
INTEGER, SAVE :: NDIFRC     ! numerical diffusion
INTEGER, SAVE :: NRELRC     ! relaxation
INTEGER, SAVE :: NDCONVRC   ! Deep CONVection
INTEGER, SAVE :: NHTURBRC   ! horizontal turbulence 
INTEGER, SAVE :: NVTURBRC   ! vertical turbulence
INTEGER, SAVE :: NNEGARC    ! negative correction                            
INTEGER, SAVE :: NNETURRC    ! negative correction                            
INTEGER, SAVE :: NNECONRC    ! negative correction                            
INTEGER, SAVE :: NNEADVRC    ! negative correction                            
INTEGER, SAVE :: NACCRRC    ! accretion
INTEGER, SAVE :: NAUTORC    ! autoconversion
INTEGER, SAVE :: NCONDRC    ! evaporation/condensation
INTEGER, SAVE :: NHONRC     ! HOmogeneous Nucleation  ICE3
INTEGER, SAVE :: NRIMRC     ! RIMing of cloudwater    ICE3
INTEGER, SAVE :: NCMELRC    ! collection by snow and conversion into rain with T>XTT ICE3
INTEGER, SAVE :: NWETGRC    ! WET Growth of graupel   ICE3
INTEGER, SAVE :: NDRYGRC    ! DRY Growth of graupel   ICE3
INTEGER, SAVE :: NIMLTRC    ! Ice MeLTing             ICE3
INTEGER, SAVE :: NBERFIRC   ! BERgeron-FIndeisen gth. ICE3
INTEGER, SAVE :: NCDEPIRC   ! Cond./DEPosition on ice ICE3
INTEGER, SAVE :: NHENURC    ! CCN Activation C2R2
INTEGER, SAVE :: NSEDIRC    ! sedimentation  C2R2
INTEGER, SAVE :: NDEPORC    ! ground deposition     
INTEGER, SAVE :: NDEPOTRRC  ! deposition on tree
INTEGER, SAVE :: NWETHRC    ! wet growth of hail
INTEGER, SAVE :: NDRYHRC    ! dry growth of hail      ICE4
INTEGER, SAVE :: NADJURC    ! adjustement before rain_ice ICE3
INTEGER, SAVE :: NHINCRC    ! Heterogeneous Nucleation by Contact LIMA
INTEGER, SAVE :: NHONCRC    ! droplet homogeneous nucleation      LIMA
INTEGER, SAVE :: NCEDSRC    ! adjustment                          LIMA
INTEGER, SAVE :: NREVARC    ! evaporation of rain drops
INTEGER, SAVE :: NCORRRC    ! rain <-> cloud transfer at the beginning of LIMA
INTEGER, SAVE :: NR2C1RC    ! rain -> cloud change after sedimentation in LIMA
INTEGER, SAVE :: NCVRCRC    ! rain -> cloud change after other microphysical processes in LIMA
!
!      Allowed processes for the budget of moist variable RRR (rain water)
!
! Courant namelist: NAM_BURRR
!
LOGICAL, SAVE :: LBU_RRR    ! True when the budget of RRR is performed
!
INTEGER, SAVE :: NASSERR    ! time filter
INTEGER, SAVE :: NNESTRR    ! Efffect of 2way nesting on Rr
INTEGER, SAVE :: NADVRR     ! Total advection for PPM
INTEGER, SAVE :: NFRCRR     ! forcing
INTEGER, SAVE :: NDIFRR     ! numerical diffusion
INTEGER, SAVE :: NRELRR     ! relaxation
INTEGER, SAVE :: NNEGARR    ! negative correction                            
INTEGER, SAVE :: NACCRRR    ! accretion
INTEGER, SAVE :: NAUTORR    ! autoconversion
INTEGER, SAVE :: NREVARR    ! rain evaporation
INTEGER, SAVE :: NSEDIRR    ! sedimentation
INTEGER, SAVE :: NSFRRR     ! Spontaneous FReezing    ICE3
INTEGER, SAVE :: NACCRR     ! ACCretion of rainwater  ICE3
INTEGER, SAVE :: NCMELRR    ! collection of droplets by snow and conversion into rain with T>XTT ICE3
INTEGER, SAVE :: NCFRZRR    ! Conversion FReeZing     ICE3
INTEGER, SAVE :: NWETGRR    ! WET Growth of graupel   ICE3
INTEGER, SAVE :: NDRYGRR    ! DRY Growth of graupel   ICE3
INTEGER, SAVE :: NGMLTRR    ! Graupel MeLTing         ICE3
INTEGER, SAVE :: NWETHRR    ! wet growth of hail      ICE4
INTEGER, SAVE :: NDRYHRR    ! dry growth of hail      ICE4
INTEGER, SAVE :: NHMLTRR    ! melting of hail         ICE4
INTEGER, SAVE :: NCORRRR    ! tendencies correction after ICE3
INTEGER, SAVE :: NHONRRR    ! drop homogeneous nucleation LIMA
INTEGER, SAVE :: NR2C1RR    ! rain -> cloud change after sedimentation in LIMA
INTEGER, SAVE :: NCVRCRR    ! rain -> cloud change after other microphysical processes in LIMA
!
!      Allowed processes for the budget of moist variable RRI (ice)
!
! Courant namelist: NAM_BURRI
!
LOGICAL, SAVE :: LBU_RRI    ! True when the budget of RRI is performed
!
INTEGER, SAVE :: NASSERI    ! time filter
INTEGER, SAVE :: NNESTRI    ! Efffect of 2way nesting on Ri
INTEGER, SAVE :: NADVRI     ! Total advection for PPM
INTEGER, SAVE :: NFRCRI     ! forcing
INTEGER, SAVE :: NDIFRI     ! numerical diffusion
INTEGER, SAVE :: NRELRI     ! relaxation
INTEGER, SAVE :: NDCONVRI   ! Deep CONVection
INTEGER, SAVE :: NHTURBRI   ! horizontal turbulence
INTEGER, SAVE :: NVTURBRI   ! vertical turbulence
INTEGER, SAVE :: NNEGARI    ! negative correction                            
INTEGER, SAVE :: NSEDIRI    ! SEDImentation           ICE3
INTEGER, SAVE :: NHENURI    ! HEterogenous NUcleation ICE3
INTEGER, SAVE :: NHONRI     ! HOmogeneous Nucleation  ICE3
INTEGER, SAVE :: NAGGSRI    ! AGGregation of snow     ICE3
INTEGER, SAVE :: NAUTSRI    ! AUToconversion of ice   ICE3
INTEGER, SAVE :: NCFRZRI    ! Conversion FReeZing     ICE3
INTEGER, SAVE :: NWETGRI    ! WET Growth of graupel   ICE3
INTEGER, SAVE :: NDRYGRI    ! DRY Growth of graupel   ICE3
INTEGER, SAVE :: NIMLTRI    ! Ice MeLTing             ICE3
INTEGER, SAVE :: NBERFIRI   ! BERgeron-FIndeisen gth. ICE3
INTEGER, SAVE :: NCDEPIRI   ! Cond./DEPosition on ice ICE3
INTEGER, SAVE :: NWETHRI    ! wet growth of hail      ICE4
INTEGER, SAVE :: NDRYHRI    ! dry growth of hail      ICE4
INTEGER, SAVE :: NADJURI    ! adjustement before rain_ice ICE3
INTEGER, SAVE :: NHINDRI ! heterogeneous nucleation by deposition LIMA
INTEGER, SAVE :: NHINCRI ! heterogeneous nucleation by contact    LIMA
INTEGER, SAVE :: NHONHRI ! haze homogeneous nucleation source     LIMA
INTEGER, SAVE :: NHONCRI ! droplet homogeneous nucleation         LIMA
INTEGER, SAVE :: NCNVIRI ! Conversion of snow to r_i              LIMA
INTEGER, SAVE :: NCNVSRI ! Conversion of pristine ice to r_s      LIMA
INTEGER, SAVE :: NHMSRI  ! Hallett-Mossop ice multiplication process due to snow riming LIMA
INTEGER, SAVE :: NHMGRI  ! Hallett-Mossop ice multiplication process due to graupel riming LIMA
INTEGER, SAVE :: NCEDSRI ! adjustement LIMA
INTEGER, SAVE :: NCORRRI    ! ice <-> snow transfer at the beginning of LIMA
!
!      Allowed processes for the budget of moist variable RRS (snow)
!
! Courant namelist: NAM_BURRS
!
LOGICAL, SAVE :: LBU_RRS    ! True when the budget of RRS is performed
!
INTEGER, SAVE :: NASSERS    ! time filter
INTEGER, SAVE :: NNESTRS    ! Efffect of 2way nesting on Rs
INTEGER, SAVE :: NADVRS     ! Total advection for PPM
INTEGER, SAVE :: NFRCRS     ! forcing
INTEGER, SAVE :: NDIFRS     ! numerical diffusion
INTEGER, SAVE :: NRELRS     ! relaxation
INTEGER, SAVE :: NNEGARS    ! negative correction                            
INTEGER, SAVE :: NSEDIRS    ! SEDImentation           ICE3
INTEGER, SAVE :: NDEPSRS    ! DEPosition on Snow      ICE3
INTEGER, SAVE :: NAGGSRS    ! AGGregation of snow     ICE3
INTEGER, SAVE :: NAUTSRS    ! AUToconversion of ice   ICE3
INTEGER, SAVE :: NRIMRS     ! RIMing of cloudwater    ICE3
INTEGER, SAVE :: NACCRS     ! ACCretion of rainwater  ICE3
INTEGER, SAVE :: NCMELRS    ! Conversion MeLTing      ICE3
INTEGER, SAVE :: NWETGRS    ! WET Growth of graupel   ICE3
INTEGER, SAVE :: NDRYGRS    ! DRY Growth of graupel   ICE3
INTEGER, SAVE :: NWETHRS    ! wet growth of hail      ICE4
INTEGER, SAVE :: NDRYHRS    ! dry growth of hail      ICE4
INTEGER, SAVE :: NCNVIRS   ! Conversion of snow to r_i         LIMA
INTEGER, SAVE :: NCNVSRS   ! Conversion of pristine ice to r_s LIMA
INTEGER, SAVE :: NHMSRS    ! Hallett-Mossop ice multiplication process due to snow riming LIMA
INTEGER, SAVE :: NCORRRS    ! ice <-> snow transfer at the beginning of LIMA
!
!      Allowed processes for the budget of moist variable RRG (graupel)
!
! Courant namelist: NAM_BURRG
!
LOGICAL, SAVE :: LBU_RRG    ! True when the budget of RRG is performed
!
INTEGER, SAVE :: NASSERG    ! time filter
INTEGER, SAVE :: NNESTRG    ! Efffect of 2way nesting on Rg
INTEGER, SAVE :: NADVRG     ! Total advection for PPM
INTEGER, SAVE :: NFRCRG     ! forcing
INTEGER, SAVE :: NDIFRG     ! numerical diffusion
INTEGER, SAVE :: NRELRG     ! relaxation
INTEGER, SAVE :: NNEGARG    ! negative correction                            
INTEGER, SAVE :: NSEDIRG    ! SEDImentation           ICE3
INTEGER, SAVE :: NSFRRG     ! Spontaneous FReezing    ICE3
INTEGER, SAVE :: NDEPGRG    ! DEPosition on Snow      ICE3
INTEGER, SAVE :: NRIMRG     ! RIMing of cloudwater    ICE3
INTEGER, SAVE :: NACCRG     ! ACCretion of rainwater  ICE3
INTEGER, SAVE :: NCMELRG    ! Conversion MeLTing      ICE3
INTEGER, SAVE :: NCFRZRG    ! Conversion FReeZing     ICE3
INTEGER, SAVE :: NWETGRG    ! WET Growth of graupel   ICE3
INTEGER, SAVE :: NDRYGRG    ! DRY Growth of graupel   ICE3
INTEGER, SAVE :: NGMLTRG    ! Graupel MeLTing         ICE3
INTEGER, SAVE :: NWETHRG    ! wet growth of hail      ICE4
INTEGER, SAVE :: NDRYHRG    ! dry growth of hail      ICE4
INTEGER, SAVE :: NCORRRG    ! tendencies correction after ICE3
INTEGER, SAVE :: NHGCVRG    ! Hail to Graupel ConVersion ICE4
INTEGER, SAVE :: NGHCVRG    ! Graupel to Hail ConVersion ICE4
INTEGER, SAVE :: NHONRRG    ! drop homogeneous nucleation LIMA
INTEGER, SAVE :: NHMGRG     ! Hallett-Mossop ice multiplication process due to graupel riming
INTEGER, SAVE :: NCOHGRG    ! conversion of hail to graupel
!
!      Allowed processes for the budget of moist variable RRH (hail)
!
! Courant namelist: NAM_BURRH
!
LOGICAL, SAVE :: LBU_RRH    ! True when the budget of RRH is performed
!
INTEGER, SAVE :: NASSERH    ! time filter
INTEGER, SAVE :: NNESTRH    ! Efffect of 2way nesting on Rh
INTEGER, SAVE :: NADVRH     ! Total advection for PPM
INTEGER, SAVE :: NFRCRH     ! forcing
INTEGER, SAVE :: NDIFRH     ! numerical diffusion
INTEGER, SAVE :: NRELRH     ! relaxation
INTEGER, SAVE :: NNEGARH    ! negative correction 
INTEGER, SAVE :: NSEDIRH    ! sedimentation
INTEGER, SAVE :: NWETGRH    ! wet growth of graupel
INTEGER, SAVE :: NWETHRH    ! wet growth of hail
INTEGER, SAVE :: NCOHGRH    ! reconversion from hail to graupel LIMA
INTEGER, SAVE :: NDRYHRH    ! dry growth of hail      ICE4
INTEGER, SAVE :: NHMLTRH    ! melting                           
INTEGER, SAVE :: NCORRRH    ! tendencies correction after ICE3
INTEGER, SAVE :: NHGCVRH    ! Hail to Graupel ConVersion ICE4
INTEGER, SAVE :: NGHCVRH    ! Graupel to Hail ConVersion ICE4
!
! Courant namelist: NAM_BURSV
!
LOGICAL, SAVE :: LBU_RSV    ! True when the budget of RSVx is performed
!
INTEGER, SAVE :: NASSESV    ! Asselin-Robert time filter
INTEGER, SAVE :: NNESTSV    ! Efffect of 2way nesting on Sv
INTEGER, SAVE :: NADVSV     ! Total advection for PPM
INTEGER, SAVE :: NFRCSV     ! forcing
INTEGER, SAVE :: NDIFSV     ! numerical diffusion
INTEGER, SAVE :: NRELSV     ! relaxation
INTEGER, SAVE :: NDCONVSV   !  Deep CONVection
INTEGER, SAVE :: NMAFLSV    ! mass flux            
INTEGER, SAVE :: NDEPOSV    ! deposition on the ground
INTEGER, SAVE :: NDEPOTRSV  ! deposition on tree    
INTEGER, SAVE :: NHTURBSV   ! horizontal turbulence
INTEGER, SAVE :: NVTURBSV   ! vertical turbulence
INTEGER, SAVE :: NCHEMSV    ! chemistry activity
!
INTEGER, SAVE :: NNEGASV
!
! Allowed processes for the budget of electric charge carried by water vapor
INTEGER, SAVE :: NDEPSQV
INTEGER, SAVE :: NDEPGQV
INTEGER, SAVE :: NREVAQV
INTEGER, SAVE :: NDEPIQV
INTEGER, SAVE :: NNEUTQV
!
! Allowed processes for the budget of electric charge carried by cloud droplets
INTEGER, SAVE :: NAUTOQC
INTEGER, SAVE :: NACCRQC
INTEGER, SAVE :: NRIMQC
INTEGER, SAVE :: NWETGQC
INTEGER, SAVE :: NDRYGQC
INTEGER, SAVE :: NIMLTQC
INTEGER, SAVE :: NBERFIQC
INTEGER, SAVE :: NDEPIQC
INTEGER, SAVE :: NINDQC  ! inductive process
INTEGER, SAVE :: NSEDIQC
INTEGER, SAVE :: NNEUTQC
!
! Allowed processes for the budget of electric charge carried by rain drops
INTEGER, SAVE :: NAUTOQR
INTEGER, SAVE :: NACCRQR
INTEGER, SAVE :: NREVAQR
INTEGER, SAVE :: NACCQR
INTEGER, SAVE :: NCFRZQR
INTEGER, SAVE :: NWETGQR
INTEGER, SAVE :: NDRYGQR
INTEGER, SAVE :: NGMLTQR
INTEGER, SAVE :: NSEDIQR
INTEGER, SAVE :: NNEUTQR
!
! Allowed processes for the budget of electric charge carried by ice crystals
INTEGER, SAVE :: NAGGSQI
INTEGER, SAVE :: NAUTSQI
INTEGER, SAVE :: NCFRZQI
INTEGER, SAVE :: NWETGQI
INTEGER, SAVE :: NDRYGQI
INTEGER, SAVE :: NIMLTQI
INTEGER, SAVE :: NBERFIQI
INTEGER, SAVE :: NDEPIQI
INTEGER, SAVE :: NNIISQI  ! non-inductive I-S
INTEGER, SAVE :: NSEDIQI
INTEGER, SAVE :: NNEUTQI
!
! Allowed processes for the budget of electric charge carried by snow
INTEGER, SAVE :: NDEPSQS
INTEGER, SAVE :: NAGGSQS
INTEGER, SAVE :: NAUTSQS
INTEGER, SAVE :: NRIMQS
INTEGER, SAVE :: NACCQS
INTEGER, SAVE :: NCMELQS
INTEGER, SAVE :: NWETGQS
INTEGER, SAVE :: NDRYGQS
INTEGER, SAVE :: NNIISQS  ! non-inductive I-S
INTEGER, SAVE :: NSEDIQS
INTEGER, SAVE :: NNEUTQS
!
! Allowed processes for the budget of electric charge carried by graupel
INTEGER, SAVE :: NDEPGQG
INTEGER, SAVE :: NRIMQG
INTEGER, SAVE :: NACCQG
INTEGER, SAVE :: NCMELQG
INTEGER, SAVE :: NCFRZQG
INTEGER, SAVE :: NWETGQG
INTEGER, SAVE :: NDRYGQG
INTEGER, SAVE :: NGMLTQG
INTEGER, SAVE :: NINDQG  ! inductive process
INTEGER, SAVE :: NSEDIQG
INTEGER, SAVE :: NNEUTQG
!
! must add processes for electric charge carried by hail
!
!
REAL :: XTIME_BU          ! budget time in this time-step
REAL :: XTIME_BU_PROCESS  ! budget time per process for this time-step
!
LOGICAL :: LBUDGET_U  ! flag to compute budget of RhoJu  and/or LES budgets with u
LOGICAL :: LBUDGET_V  ! flag to compute budget of RhoJv  and/or LES budgets with u
LOGICAL :: LBUDGET_W  ! flag to compute budget of RhoJw  and/or LES budgets with u
LOGICAL :: LBUDGET_TH ! flag to compute budget of RhoJTh and/or LES budgets with th
LOGICAL :: LBUDGET_TKE! flag to compute budget of RhoJTke and/or LES budgets with Tke
LOGICAL :: LBUDGET_RV ! flag to compute budget of RhoJrv and/or LES budgets with rv
LOGICAL :: LBUDGET_RC ! flag to compute budget of RhoJrc and/or LES budgets with rc
LOGICAL :: LBUDGET_RR ! flag to compute budget of RhoJrr and/or LES budgets with rr
LOGICAL :: LBUDGET_RI ! flag to compute budget of RhoJri and/or LES budgets with ri
LOGICAL :: LBUDGET_RS ! flag to compute budget of RhoJrs and/or LES budgets with rs
LOGICAL :: LBUDGET_RG ! flag to compute budget of RhoJrg and/or LES budgets with rg
LOGICAL :: LBUDGET_RH ! flag to compute budget of RhoJrh and/or LES budgets with rh
LOGICAL :: LBUDGET_SV ! flag to compute budget of RhoJsv and/or LES budgets with sv
!
END MODULE MODD_BUDGET
