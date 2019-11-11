!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!####################################
MODULE MODI_READ_ALL_DATA_MESONH_CASE
!####################################
INTERFACE
      SUBROUTINE READ_ALL_DATA_MESONH_CASE(TZPRE_REAL1,HFMFILE,TPPGDFILE, &
                              HDAD_NAME                                   )
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),POINTER, INTENT(INOUT) :: TZPRE_REAL1 !PRE_REAL1 file
CHARACTER(LEN=28), INTENT(IN)    :: HFMFILE    ! name of the Mesonh input file
TYPE(TFILEDATA),   INTENT(IN)    :: TPPGDFILE  ! physiographic data file
CHARACTER(LEN=*),  INTENT(INOUT) :: HDAD_NAME  ! true name of the Mesonh input file
!
END SUBROUTINE READ_ALL_DATA_MESONH_CASE
!
END INTERFACE
!
END MODULE MODI_READ_ALL_DATA_MESONH_CASE
!
!     #####################################################################
      SUBROUTINE READ_ALL_DATA_MESONH_CASE(TZPRE_REAL1,HFMFILE,TPPGDFILE, &
                              HDAD_NAME                                   )
!     #####################################################################
!
!!****  *READ_ALL_DATA_MESONH_CASE* - reads data for the initialization of real cases.
!! 
!!    PURPOSE
!!    -------
!!    This routine calls the routines reading the different input files.
!!    It closes the PGD file.
!!
!!**  METHOD
!!    ------
!!    Projection is read in READ_LFIFM_PGD (MODD_GRID).
!!    Grid and definition of large domain are read in PGD file and MESONH file
!!     (READ_LFIFM_PGD and MESONH_INFO in READ_GENERAL) and coherence is tested.
!!    The PGD files are also read in READ_LFIFM_PGD.
!!    The PGD file is closed.
!!    The MESO-NH domain is defined from PRE_REAL1.nam inputs in SET_SUBDOMAIN.
!!    Vertical grid is defined in READ_VER_GRID.
!!    Variables fields are read in Mesonh file and stored on MESO-NH domain
!!     (READ_LFI_PRC).
!!
!!    EXTERNAL
!!    --------
!!    subroutine READ_LFIFM_PGD: to read PGD file
!!    subroutine READ_GRID_TIME_MESONH_CASE  : to read the geographic and time data.
!!    subroutine SET_SUBDOMAIN : to define the horizontal MESO-NH domain.
!!    subroutine READ_VER_GRID : to read the vertical grid in namelist file.
!!    subroutine READ_PRC_FMFILE : to read the large scale fields on the MESO-NH
!!                               domain.
!!    
!!    Module MODI_READ_GRID_TIME_MESONH_CASE  : interface for subroutine
!!                                              READ_GRID_TIME_MESONH_CASE
!!    Module MODI_SET_SUBDOMAIN : interface for subroutine SET_SUBDOMAIN
!!    Module MODI_READ_VER_GRID : interface for subroutine READ_VER_GRID
!!    Module MODI_READ_PRC_FMFILE  : interface for subroutine READ_PRC_FMFILE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     : contains logical unit names for all models
!!         TLUOUT0 : output-listing file
!!      Module MODD_PGDDIM    : contains dimension of PGD fields
!!         NPGDIMAX: dimension along x (no external point)
!!         NPGDJMAX: dimension along y (no external point)
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    16/12/94
!!                  24/09/96 (V. Masson) add HDAD_NAME
!!                  25/10/96 (V. Masson) remove 3D pressure from routine outputs
!!                  12/12/96 (V. Masson) add vertical velocity
!!                  07/05/97 (V. Masson) add tke
!!                  10/06/97 (V. Masson) add pressure
!!                  10/07/97 (V. Masson) add epsilon
!!                  11/07/97 (V. Masson) add READ_DESFM
!!                  11/07/97 (V. Masson) add passive scalars
!!                  20/01/98 (J. Stein)  add the Large scale fields reading
!!                  30/04/98 (V. Masson) Large scale VEG and LAI
!!                  02/07/00 (F.Solmon/V.Masson) adaptation for patch approach
!!                  01/2004  (V. Masson) removes surface (externalization)
!!                  01/06/02 (O.Nuissier) bogussing of tropical cyclone
!!                  Aou   09, 2005 (D.Barbary) call to compare_dad
!!                  19/03/2008 (J.Escobar) rename INIT to INIT_MNH --> grib problem
!!                  2014 (M.Faivre)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_FM
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_READ_GRID_TIME_MESONH_CASE ! interface modules
USE MODI_READ_HGRID
USE MODI_READ_VER_GRID
USE MODI_READ_PRC_FMFILE  
USE MODI_CHECK_ZS
USE MODI_CHECK_ZHAT
USE MODI_SET_BOGUS_VORTEX
USE MODI_COMPARE_DAD
USE MODI_ZS_BOUNDARY
!
USE MODD_CONF           ! declaration modules
USE MODD_CONF_n
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_PARAM_n
USE MODD_LUNIT
USE MODD_LUNIT_n
USE MODD_PGDDIM
USE MODD_PARAMETERS
USE MODD_REF_n
USE MODD_CST
USE MODD_TIME_n
USE MODD_DIM_n
USE MODD_DYN_n
USE MODD_GRID_n
USE MODD_GR_FIELD_n
USE MODD_FIELD_n
USE MODD_LBC_n
USE MODD_HURR_CONF, ONLY: LBOGUSSING, CDADATMFILE,CDADBOGFILE
USE MODD_PREP_REAL
!
USE MODI_INIT_MNH
!
!20131113 add modules for update_halo and check
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
TYPE(TFILEDATA),POINTER, INTENT(INOUT) :: TZPRE_REAL1 !PRE_REAL1 file
CHARACTER(LEN=28), INTENT(IN)    :: HFMFILE    ! name of the Mesonh input file
TYPE(TFILEDATA),   INTENT(IN)    :: TPPGDFILE  ! physiographic data file
CHARACTER(LEN=*),  INTENT(INOUT) :: HDAD_NAME  ! true name of the Mesonh input file
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IIINF_LS   ! lower boundary of the MESO-NH domain
!                     ! in the large domain (direction I).
INTEGER :: IISUP_LS   ! upper boundary of the MESO-NH domain
!                     ! in the large domain (direction I).
INTEGER :: IJINF_LS   ! lower boundary of the MESO-NH domain
!                     ! in the large domain (direction J).
INTEGER :: IJSUP_LS   ! upper boundary of the MESO-NH domain
!                     ! in the large domain (direction J).
INTEGER :: IXOR_LS    ! I shift between PGD file and LS atmospheric file
INTEGER :: IYOR_LS    ! J shift between PGD file and LS atmospheric file
!
INTEGER :: IRESP      ! return-code if problems occured
INTEGER :: ILUOUT0    ! logical unit for file TLUOUT0
!
CHARACTER(LEN=28) :: YPGD_NAME, YPGD_DAD_NAME
CHARACTER(LEN=28) :: YOUTFILE
CHARACTER(LEN=2)  :: YPGD_TYPE
!
!* temporary namelist configuration variables
!
INTEGER                           :: IVERB      ! verbosity level
CHARACTER(LEN=3)                  :: YEQNSYS    ! EQuatioN SYStem
CHARACTER(LEN=5)                  :: YPRESOPT   ! PRESsure OPTion
LOGICAL                           :: GRES
REAL                              :: ZRES
INTEGER                           :: IPRE_REAL1
!
!20131113 add vars related to ADD3DFIELD and UPDATE_HALO
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll => NULL() ! list of fields to exchange
!
!------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!-------------------------------------------------------------------------------
!
!*            1.  reading of the atm. fields, on the entire atm. file domain
!                 ----------------------------------------------------------
!
!*       1.1   save of some configuration variables before call to INIT
!
YEQNSYS = CEQNSYS
YPRESOPT= CPRESOPT
IVERB   = NVERB
GRES    = LRES
ZRES    = XRES
!
!* note that some quantities such as grid and domain definition
! will be erased by the following routines
!
YOUTFILE=CINIFILE
CINIFILE=HFMFILE
!
CALL IO_FILE_CLOSE_ll(TZPRE_REAL1)
CALL INIT_MNH
CALL IO_FILE_OPEN_ll(TZPRE_REAL1)
!
CINIFILE=YOUTFILE
!
CSTORAGE_TYPE='TT'       ! after prep_real_case, file type is 'TT'
!
!*       1.2   restores values after call to INIT
!
IF (LEN_TRIM(YEQNSYS)>0)  CEQNSYS = YEQNSYS
IF (LEN_TRIM(YPRESOPT)>0) CPRESOPT= YPRESOPT
IF (IVERB/=NUNDEF)        NVERB   = IVERB
LRES = GRES
XRES = ZRES
!
!*            2. Reading of physiographic data domain
!                ------------------------------------
!
CALL READ_HGRID(0,TPPGDFILE,YPGD_NAME,YPGD_DAD_NAME,YPGD_TYPE)
!
!*            3. Reading of large-scale grid and time
!                ------------------------------------
!
CALL READ_GRID_TIME_MESONH_CASE(HFMFILE,IXOR_LS,IYOR_LS,HDAD_NAME)
!
!*            4. Definition of horizontal domain
!                -------------------------------
!
IIINF_LS = 1
IISUP_LS = NIMAX + 2* JPHEXT
IJINF_LS = 1
IJSUP_LS = NJMAX + 2* JPHEXT
!
!
!*            5. Reading of vertical grid
!                ------------------------
!
CALL READ_VER_GRID(TZPRE_REAL1,XZHAT_LS,LSLEVE_LS,XLEN1_LS,XLEN2_LS)
!
!*            6. Add a bogus vortex from observations (radar, satellite,...)
!                -----------------------------------------------------------
!
IF (LBOGUSSING) THEN
  IF (LEN_TRIM(CDADBOGFILE) >0 ) THEN
    IF (LEN_TRIM(CDADATMFILE) == 0 ) THEN
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_ALL_DATA_MESONH_CASE','CDADATMFILE not initialized in nml NAM_HURR_CONF')
    ELSE 
      IF (LEN_TRIM(HDAD_NAME) >0 ) THEN
        IF (ADJUSTL(ADJUSTR(HDAD_NAME)) .NE. ADJUSTL(ADJUSTR(CDADATMFILE))) THEN
          WRITE(ILUOUT0,*) &
           'READ_ALL_DATA_MESONH_CASE: CDADATMFILE is NOT the DAD of input file'
          WRITE(ILUOUT0,*) ' CDADATMFILE='//TRIM(CDADATMFILE)
          WRITE(ILUOUT0,*) ' DAD_NAME of model1='//TRIM(HDAD_NAME)
          WRITE(ILUOUT0,*) '-> JOB ABORTED'
 !callabortstop
          CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_ALL_DATA_MESONH_CASE','')
        ELSE    
          CALL COMPARE_DAD(CDADATMFILE,CDADBOGFILE,IRESP)
          IF (IRESP .NE. 0) THEN
 !callabortstop
            CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_ALL_DATA_MESONH_CASE',&
                           'unable to replace the DAD of output file with CDADBOGFILE')
          ENDIF
          HDAD_NAME=CDADBOGFILE
        ENDIF
      ENDIF ! no DAD_NAME (bogussing of model1)
    ENDIF
  ENDIF
  CALL SET_BOGUS_VORTEX(XUT,XVT,XTHT)
ENDIF
!
!*            7.  truncation of the fields on the new domain
!                 ------------------------------------------
!
IF (IIINF_LS/=1 .OR. IISUP_LS/=SIZE(XTHT,1)) CLBCX(:) ='OPEN'
IF (IJINF_LS/=1 .OR. IJSUP_LS/=SIZE(XTHT,2)) CLBCY(:) ='OPEN'
!
CALL READ_PRC_FMFILE(IIINF_LS,IISUP_LS,IJINF_LS,IJSUP_LS                     )
!
!*            8.  Orography
!                 ---------
!
ALLOCATE(XZS(IISUP_LS-IIINF_LS+1,IJSUP_LS-IJINF_LS+1))
CALL IO_READ_FIELD(TPPGDFILE,'ZS',XZS)
CALL ZS_BOUNDARY(XZS,XZS_LS)
!
ALLOCATE(XZSMT(IISUP_LS-IIINF_LS+1,IJSUP_LS-IJINF_LS+1))
IF (TPPGDFILE%NMNHVERSION(1)<4 .OR. (TPPGDFILE%NMNHVERSION(1)==4 .AND. TPPGDFILE%NMNHVERSION(2)<=6)) THEN
  XZSMT = XZS
ELSE
  CALL IO_READ_FIELD(TPPGDFILE,'ZSMT',XZSMT)
END IF 
CALL ZS_BOUNDARY(XZSMT,XZSMT_LS)
!
!
!*            9. Surface pressure computation
!                ----------------------------
!
ALLOCATE(XPS_LS    (IISUP_LS-IIINF_LS+1,IJSUP_LS-IJINF_LS+1))
!
XPS_LS(:,:) = XP00* (                                                     &
                      0.5 * ( (XPMASS_LS(:,:,JPVEXT  )/XP00)**(XRD/XCPD)  &
                             +(XPMASS_LS(:,:,JPVEXT+1)/XP00)**(XRD/XCPD)  &
                            )                                             &
                     )**(XCPD/XRD)
!
!20131113 add update_halo
CALL ADD2DFIELD_ll(TZFIELDS_ll,XPS_LS )
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
CALL MPPDB_CHECK2D(XPS_LS,"PGDFILTER9:XPS_LS",PRECISION)
!
!
!*           10. Check coherence between the 2 orographies
!                -----------------------------------------
!
!20131023 mise en commentaire du check_zs et zhat
!
!IF (LEN_TRIM(HDAD_NAME)>0) CALL CHECK_ZS(HFMFILE,HDAD_NAME,IIINF_LS,IJINF_LS)
!IF (LEN_TRIM(HDAD_NAME)>0) CALL CHECK_ZHAT(HFMFILE,HDAD_NAME)
!
!-------------------------------------------------------------------------------
!
WRITE (ILUOUT0,*) 'Monitor READ_ALL_DATA_MESONH_CASE completed'
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_ALL_DATA_MESONH_CASE
