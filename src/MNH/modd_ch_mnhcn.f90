!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/modd_ch_mnhcn.f90,v $ $Revision: 1.2.4.1.2.1.12.2 $ $Date: 2014/01/09 15:01:56 $
!-----------------------------------------------------------------
!!    #####################
      MODULE MODD_CH_MNHC_n
!!    #####################
!!
!!*** *MODD_CH_MNHC$n*
!!
!!    PURPOSE
!!    -------
!       This module contains all parameters that control the 
!     chemical part of MesoNH.
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre      *Laboratoire d'Aerollogie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/05/95
!!    27/07/96 (K. Suhre) restructured
!!    30/11/99 (K. Suhre) add new parameters
!!    16/11/00 (C. Mari)  add new parameters
!!    28/05/02 (C. Mari)  move default values to default_desfmn
!!    13/07/03 (J.-P. Pinty) add flag for lightning production of NOx
!!    12/04/07 (M. Leriche) add flag for aqueous chemistry
!!    30/05/07 (M. Leriche) add flag and real for pH calculation
!!    25/04/08 (M. Leriche) add threshold for aqueous phase chemistry
!!    16/09/10 (M. Leriche) add flag for ice phase chemistry
!!    13/01/11 (M. Leriche) add flag for retention in ice 
!!    24/03/16 (M. Leriche) remove surface option -> manage them in SURFEX
!!    01/10/16 (F. Brosse)  add production/destruction terms computation
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_MNHC_t
!
! switch which indicates whether chemistry is used or not
!
  LOGICAL :: LUSECHEM
!
! switch which indicates whether aqueous chemistry is used or not
!
  LOGICAL :: LUSECHAQ
!
! switch which indicates whether ice phase chemistry is used or not
!
  LOGICAL :: LUSECHIC
!
!* Initialization
!
  LOGICAL :: LCH_INIT_FIELD ! flag indicating whether initialization
                 ! of chemical fields shall be done during MesoNH run using
                 ! CH_INIT_FIELD (overwrites initial values from FM-files)
                 ! or not
!
!* Scavenging in parameterized convective clouds
!
  LOGICAL :: LCH_CONV_SCAV 
                 ! flag for calculation of scavenging 
                 ! by convective precipitations (active only if LCHTRANS=.TRUE.)
!
!* pH calculation
!
  LOGICAL :: LCH_PH
               ! flag for calculation of pH
  REAL :: XCH_PHINIT
              ! initial pH value if pH is calculated or
              ! pH value if pH is fixed
!
!* threshold for aqueous phase chemistry
!
  REAL    :: XRTMIN_AQ
               ! minimal water mixing ratio (vol/vol) needed for the
               ! chemistry of the aqueous phases (cloud and rain) 
!
!* Retention in ice when lusechic = f
!
  LOGICAL  :: LCH_RET_ICE ! flag for all soluble gases retained in ice
                      
!*
  CHARACTER(LEN=10) :: CCH_SCHEME
                 ! name of chemical scheme
!
  CHARACTER(LEN=80) :: CCHEM_INPUT_FILE 
                 ! name of general 
                 ! purpose ASCII input file (handeled by CH_OPEN_INPUT)
!
  CHARACTER(LEN=10) :: CCH_TDISCRETIZATION 
                 ! temporal discretization:
                 ! "SPLIT"  : use time-splitting, input fields for solver are
                 !            scalar variables at t+dt (derived from XRSVS)
                 ! "CENTER" : input fields for solver are 
                 !            scalar variables at t (XSVT)
                 ! "LAGGED" : input fields for solver are 
                 !            scalar variables at t-dt (XSVM)
!
  INTEGER           :: NCH_SUBSTEPS
                 ! number of chemical timesteps to be taken during one 
                 ! double timestep of MesoNH (MesoNH integrates with timesteps
                 ! of lenght 2*XTSTEP using leapfrog), the timestep of the 
                 ! solver will be calculated as 
                 ! ZDTSOLVER = 2*XTSTEP / NCH_SUBSTEPS
!* LiNOx
!
  LOGICAL :: LCH_CONV_LINOX
                 ! flag for calculation of NOx production by lightnings
!
!* photolysis rates (TUV)
!
  LOGICAL      :: LCH_TUV_ONLINE  ! switch online/lookup table
  CHARACTER(LEN=80) :: CCH_TUV_LOOKUP  ! name of lookup table file
  CHARACTER(LEN=4)  :: CCH_TUV_CLOUDS  ! method for calculating the
                                ! impact of clouds on radiation
                                ! "FOUQ" (model clouds, only 1-D)
  REAL :: XCH_TUV_ALBNEW  ! surface albedo (if negative the albedo
                        ! will be read from DATAX/albedo.dat)
  REAL :: XCH_TUV_DOBNEW  ! scaling factor for ozone column dobson
                        ! (if negative, no scaling will be performed,
                        ! note: the O3 profile will be read from
                        ! DATAX/O3.profile, if this file is empty, the
                        ! US standard O3 profile will be used)
  REAL :: XCH_TUV_TUPDATE ! update frequency for TUV (in seconds)
!
!* vectorization
!
  CHARACTER(LEN=3) :: CCH_VEC_METHOD          ! type of vectorization mask
                                       ! 'MAX' take NCH_VEC_LENGTH points
                                       ! 'TOT' take all grid points
                                       ! 'HOR' take horizontal layers
                                       ! 'VER' take vertical columns
  INTEGER          :: NCH_VEC_LENGTH          ! number of points for 'MAX' option
!
!* 1-D time series
!
  REAL              :: XCH_TS1D_TSTEP         ! time between two call to write_ts1d
  CHARACTER(LEN=80) :: CCH_TS1D_COMMENT       ! comment for write_ts1d
  CHARACTER(LEN=80) :: CCH_TS1D_FILENAME      ! filename for write_ts1d files
!
!* total production/loss for chemical species
! 
  CHARACTER(LEN=1024)    :: CSPEC_PRODLOSS
!
!* extended production/loss terms for chemical species
!
  CHARACTER(LEN=1024)    :: CSPEC_BUDGET
!
!
END TYPE CH_MNHC_t

TYPE(CH_MNHC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_MNHC_MODEL

LOGICAL, POINTER :: LUSECHEM=>NULL()
LOGICAL, POINTER :: LUSECHAQ=>NULL()
LOGICAL, POINTER :: LUSECHIC=>NULL()
LOGICAL, POINTER :: LCH_INIT_FIELD=>NULL()
LOGICAL, POINTER :: LCH_CONV_SCAV=>NULL()
LOGICAL, POINTER :: LCH_PH=>NULL()
LOGICAL, POINTER :: LCH_RET_ICE=>NULL()
REAL, POINTER :: XCH_PHINIT=>NULL()
REAL, POINTER :: XRTMIN_AQ=>NULL()
CHARACTER(LEN=10), POINTER :: CCH_SCHEME=>NULL()
CHARACTER(LEN=80), POINTER :: CCHEM_INPUT_FILE=>NULL()
CHARACTER(LEN=10), POINTER :: CCH_TDISCRETIZATION=>NULL()
INTEGER, POINTER :: NCH_SUBSTEPS=>NULL()
LOGICAL, POINTER :: LCH_CONV_LINOX=>NULL()
LOGICAL, POINTER :: LCH_TUV_ONLINE=>NULL()
CHARACTER(LEN=80), POINTER :: CCH_TUV_LOOKUP=>NULL()
CHARACTER(LEN=4), POINTER :: CCH_TUV_CLOUDS=>NULL()
REAL, POINTER :: XCH_TUV_ALBNEW=>NULL()
REAL, POINTER :: XCH_TUV_DOBNEW=>NULL()
REAL, POINTER :: XCH_TUV_TUPDATE=>NULL()
CHARACTER(LEN=3), POINTER :: CCH_VEC_METHOD=>NULL()
INTEGER, POINTER :: NCH_VEC_LENGTH=>NULL()
REAL, POINTER :: XCH_TS1D_TSTEP=>NULL()
CHARACTER(LEN=80), POINTER :: CCH_TS1D_COMMENT=>NULL()
CHARACTER(LEN=80), POINTER :: CCH_TS1D_FILENAME=>NULL()
CHARACTER(LEN=1024), POINTER :: CSPEC_PRODLOSS=>NULL()
CHARACTER(LEN=1024), POINTER :: CSPEC_BUDGET=>NULL()

CONTAINS

SUBROUTINE CH_MNHC_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
LUSECHEM=>CH_MNHC_MODEL(KTO)%LUSECHEM
LUSECHAQ=>CH_MNHC_MODEL(KTO)%LUSECHAQ
LUSECHIC=>CH_MNHC_MODEL(KTO)%LUSECHIC
LCH_INIT_FIELD=>CH_MNHC_MODEL(KTO)%LCH_INIT_FIELD
LCH_CONV_SCAV=>CH_MNHC_MODEL(KTO)%LCH_CONV_SCAV
LCH_PH=>CH_MNHC_MODEL(KTO)%LCH_PH
LCH_RET_ICE=>CH_MNHC_MODEL(KTO)%LCH_RET_ICE
XCH_PHINIT=>CH_MNHC_MODEL(KTO)%XCH_PHINIT
XRTMIN_AQ=>CH_MNHC_MODEL(KTO)%XRTMIN_AQ
CCH_SCHEME=>CH_MNHC_MODEL(KTO)%CCH_SCHEME
CCHEM_INPUT_FILE=>CH_MNHC_MODEL(KTO)%CCHEM_INPUT_FILE
CCH_TDISCRETIZATION=>CH_MNHC_MODEL(KTO)%CCH_TDISCRETIZATION
NCH_SUBSTEPS=>CH_MNHC_MODEL(KTO)%NCH_SUBSTEPS
LCH_CONV_LINOX=>CH_MNHC_MODEL(KTO)%LCH_CONV_LINOX
LCH_TUV_ONLINE=>CH_MNHC_MODEL(KTO)%LCH_TUV_ONLINE
CCH_TUV_LOOKUP=>CH_MNHC_MODEL(KTO)%CCH_TUV_LOOKUP
CCH_TUV_CLOUDS=>CH_MNHC_MODEL(KTO)%CCH_TUV_CLOUDS
XCH_TUV_ALBNEW=>CH_MNHC_MODEL(KTO)%XCH_TUV_ALBNEW
XCH_TUV_DOBNEW=>CH_MNHC_MODEL(KTO)%XCH_TUV_DOBNEW
XCH_TUV_TUPDATE=>CH_MNHC_MODEL(KTO)%XCH_TUV_TUPDATE
CCH_VEC_METHOD=>CH_MNHC_MODEL(KTO)%CCH_VEC_METHOD
NCH_VEC_LENGTH=>CH_MNHC_MODEL(KTO)%NCH_VEC_LENGTH
XCH_TS1D_TSTEP=>CH_MNHC_MODEL(KTO)%XCH_TS1D_TSTEP
CCH_TS1D_COMMENT=>CH_MNHC_MODEL(KTO)%CCH_TS1D_COMMENT
CCH_TS1D_FILENAME=>CH_MNHC_MODEL(KTO)%CCH_TS1D_FILENAME
CSPEC_PRODLOSS=>CH_MNHC_MODEL(KTO)%CSPEC_PRODLOSS  
CSPEC_BUDGET=>CH_MNHC_MODEL(KTO)%CSPEC_BUDGET

END SUBROUTINE CH_MNHC_GOTO_MODEL

END MODULE MODD_CH_MNHC_n
