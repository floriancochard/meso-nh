!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################
      PROGRAM SPAWNING
!     ################
!
!!****  *SPAWNING * -general monitor to spawn a model from an other one.
!!
!!    PURPOSE
!!    -------
!! 
!!      This program is the general monitor to spawn a model. Firstly, it calls
!!    the subroutine INIT, which performs the initialization of the model 1.
!!      The domain size and configuration of model 2 are initialized in
!!     INIT before the data structures of ComLib set up.
!!      Then, the program calls the subroutine SPAWN_MODEL2, which initializes
!!    by horizontal interpolation, the model 2 in a sub-domain of model 1, 
!!    and writes the resulting fields in a FM-file.
!!
!!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!       subroutine INIT      : performs the initialization of the model 1
!!
!!       subroutine SPAWN_MODEL2 : performs horizontal interpolation and writes
!!                                 FM-file   
!!
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    NONE
!!
!!    REFERENCE
!!    ---------
!!
!!       NONE
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
!!      Original    10/01/95
!!      Modification 29/01/96  (Lafore) Update for MASDEV2_2 version
!!      Modification 19/03/96  (Lafore) Spawning for surface fields
!!      Modification 24/10/96  (Masson) Initialization of outer points
!!      Modification 19/11/96  (Masson) Add deep convection
!!      Modification 03/11/97  (Lafore) update of call to BOUNDARIES
!!      Modification 03/06/98  (Stein ) update of call to BOUNDARIES
!!      Modification 15/03/99  (Masson) call to READ_EXSPA and program name
!!      Modification 15/07/99  (Jabouille) create MODD_SPAWN
!!      Modification 14/12/00  (Jabouille) add NAM_BLANK reading
!!      Modification 15/10/01  (Mallet) allow namelists in different orders
!!      Modification 07/07/05  (Barbary) spawn with 2 input files (father+son1)
!!                                      to keep finest fields of son1
!!      Modification 05/06     Remove EPS
!!      Modification 19/03/2008 (J.Escobar) rename INIT to INIT_MNH --> grib problem
!!      Modification 05/02/2015 (M.Moge) read namelist NAM_CONFZ, before INIT_MNH
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 19/04/2016 : Pb IOZ/NETCDF , missing OPARALLELIO=.FALSE. for PGD files
!!  06/2016     (G.Delautier) phasage surfex 8
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!  
!*       0.1    Declarative modules common to all the models
!
USE MODD_CONF
USE MODD_CST
USE MODD_DYN
USE MODD_GRID
USE MODD_LUNIT
USE MODD_PARAMETERS
USE MODD_REF
USE MODD_SPAWN
USE MODN_BLANK
USE MODD_NSV
USE MODN_CONFZ
!  
!*       0.2    Declarative modules of model 1
!
USE MODD_CONF_n
USE MODD_CURVCOR_n
USE MODD_DIM_n
USE MODD_DYN_n, LRES_n=>LRES, XRES_n=>XRES 
USE MODD_FIELD_n
USE MODD_IO_ll, ONLY: NIO_VERB,NVERB_DEBUG,TFILEDATA
USE MODD_LSFIELD_n
USE MODD_LBC_n
USE MODD_LUNIT_n
USE MODD_OUT_n
USE MODD_PARAM_n
USE MODD_REF_n
USE MODD_TIME_n
USE MODD_CH_MNHC_n
USE MODD_GRID_n
! 
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME,IO_FILE_PRINT_LIST
USE MODE_ll
USE MODE_POS
USE MODE_FM
USE MODE_MODELN_HANDLER
!
USE MODI_SPAWN_MODEL2    
USE MODI_BOUNDARIES
!
USE MODI_VERSION
USE MODI_INIT_MNH
USE MODD_MNH_SURFEX_n
USE MODE_MPPDB
!
!
USE MODN_CONF, ONLY : JPHEXT , NHALO
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.3    Local variables
!
!
CHARACTER (LEN=28) :: YSONFILE = ' '  ! possible name of SON input FM-file
CHARACTER (LEN=28) :: YSPAFILE = ' '  ! possible name of the output FM-file
CHARACTER (LEN= 2) :: YSPANBR = '00'  ! NumBeR associated to the SPAwned file
INTEGER            :: IINFO_ll        ! return code of // routines
INTEGER :: IRESP                      ! Return codes in FM routines
INTEGER :: ILUSPA                     ! Logical unit number for the EXSPA file
LOGICAL :: GFOUND                     ! Return code when searching namelist
!
LOGICAL :: LSPAWN_SURF = .TRUE.  ! .TRUE. : surface fields are spawned
LOGICAL                           :: LRES
REAL                              :: XRES
TYPE(TFILEDATA),POINTER :: TZINIFILE  => NULL()
TYPE(TFILEDATA),POINTER :: TZEXPAFILE => NULL()
!
NAMELIST/NAM_SPAWN_SURF/LSPAWN_SURF, LRES, XRES  
NAMELIST/NAM_CONF_SPAWN/JPHEXT, NHALO
!
!-------------------------------------------------------------------------------
!
CALL MPPDB_INIT()
!
! First Switch to model 1 variables
CALL GOTO_MODEL(1,ONOFIELDLIST=.TRUE.)
!
CALL VERSION
CPROGRAM='SPAWN '
CDOMAIN= ''
!
CALL INITIO_ll()
!-------------------------------------------------------------------------------
!
!*       1.    SPAWNING INITIALIZATION 
!              -----------------------
!
CALL READ_EXSPA(CINIFILE,CINIFILEPGD,&
                NXOR,NYOR,NXSIZE,NYSIZE,NDXRATIO,NDYRATIO, &
                LBAL_ONLY, &
                CDOMAIN,YSPAFILE,YSPANBR,CDADINIFILE,CDADSPAFILE,YSONFILE)
!
!
!*       2.    NAM_BLANK, NAM_SPAWN_SURF and NAM_CONFZ READING AND EXSPA file CLOSURE
!              ----------------------------------------
!
CALL IO_FILE_FIND_BYNAME('SPAWN1.nam',TZEXPAFILE,IRESP)
CALL IO_FILE_OPEN_ll(TZEXPAFILE)
ILUSPA = TZEXPAFILE%NLU
!
CALL INIT_NMLVAR
CALL POSNAM(ILUSPA,'NAM_SPAWN_SURF',GFOUND)
IF (GFOUND) READ(UNIT=ILUSPA,NML=NAM_SPAWN_SURF)
CALL UPDATE_MODD_FROM_NMLVAR
CALL POSNAM(ILUSPA,'NAM_BLANK',GFOUND)
IF (GFOUND) READ(UNIT=ILUSPA,NML=NAM_BLANK)
CALL POSNAM(ILUSPA,'NAM_CONFZ',GFOUND)
IF (GFOUND) READ(UNIT=ILUSPA,NML=NAM_CONFZ)
CALL POSNAM(ILUSPA,'NAM_CONF_SPAWN',GFOUND)
IF (GFOUND) READ(UNIT=ILUSPA,NML=NAM_CONF_SPAWN)
CALL IO_FILE_CLOSE_ll(TZEXPAFILE)
!
!-------------------------------------------------------------------------------
!
!*       3.    MODEL 1 INITIALIZATION
!              ----------------------
!
CALL INIT_MNH
!
CALL IO_FILE_FIND_BYNAME(TRIM(CINIFILE),TZINIFILE,IRESP)
CALL IO_FILE_CLOSE_ll(TZINIFILE)
CALL IO_FILE_CLOSE_ll(TINIFILEPGD,OPARALLELIO=.FALSE.)
!-------------------------------------------------------------------------------
!
!*       4.    INITIALIZATION OF OUTER POINTS OF MODEL 1
!              -----------------------------------------
!
CALL BOUNDARIES                                                     & 
           (XTSTEP,CLBCX,CLBCY,NRR,NSV,1,                           &
            XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
            XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
            XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,   &
            XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS,   &
            XRHODJ,                                                 &
            XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT, XSRCT            )
CALL MPPDB_CHECK3D(XUT,"SPAWNING-after boundaries::XUT",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       5.    SPAWNING OF MODEL 2 FROM MODEL 1
!              --------------------------------
!
CALL IO_FILE_OPEN_ll(TZEXPAFILE)
ILUSPA = TZEXPAFILE%NLU
!
CALL SET_POINTERS_TO_MODEL1()
CALL GOTO_MODEL(2)
CALL INIT_NMLVAR
CALL POSNAM(ILUSPA,'NAM_SPAWN_SURF',GFOUND)
IF (GFOUND) READ(UNIT=ILUSPA,NML=NAM_SPAWN_SURF)
CALL UPDATE_MODD_FROM_NMLVAR
CALL GOTO_MODEL(1)
CALL IO_FILE_CLOSE_ll(TZEXPAFILE)
!
CALL GO_TOMODEL_ll(2,IINFO_ll)
!
CALL SPAWN_MODEL2 (NRR,NSV_USER,CTURB,CSURF,CCLOUD,                     &
                   CCHEM_INPUT_FILE,YSPAFILE,YSPANBR,YSONFILE,          &
                   CINIFILE, CINIFILEPGD, LSPAWN_SURF                   )
!
CALL SURFEX_DEALLO_LIST
!
IF(NIO_VERB>=NVERB_DEBUG) CALL IO_FILE_PRINT_LIST()
!
CALL IO_FILE_CLOSE_ll(TLUOUT)
!
CALL END_PARA_ll(IINFO_ll)
!JUAN CALL ABORT
STOP

CONTAINS 

SUBROUTINE INIT_NMLVAR
LRES=LRES_n
XRES=XRES_n
END SUBROUTINE INIT_NMLVAR

SUBROUTINE UPDATE_MODD_FROM_NMLVAR
LRES_n=LRES
XRES_n=XRES
END SUBROUTINE UPDATE_MODD_FROM_NMLVAR

SUBROUTINE SET_POINTERS_TO_MODEL1()
!
USE MODD_DEF_EDDY_FLUX_n
USE MODD_DEF_EDDYUV_FLUX_n
USE MODD_FIELD_n
USE MODD_GRID_n
USE MODD_LATZ_EDFLX
USE MODD_PRECIP_n
!
XXHAT1 => XXHAT
XYHAT1 => XYHAT
XZHAT1 => XZHAT
XZTOP1 => XZTOP
XZS1 => XZS
XZSMT1 => XZSMT
XZZ1 => XZZ
LSLEVE1 => LSLEVE
XLEN11 => XLEN1
XLEN21 => XLEN2
XTHT1 => XTHT
XTKET1 => XTKET
XPABST1 => XPABST
XRT1 => XRT
XUT1 => XUT
XVT1 => XVT
XWT1 => XWT
XZWS1 => XZWS
XSRCT1 => XSRCT
XSIGS1 => XSIGS
TDTCUR1 => TDTCUR
XLSUM1 => XLSUM
XLSVM1 => XLSVM
XLSWM1 => XLSWM
XLSZWSM1 => XLSZWSM
XLSTHM1 => XLSTHM
XLSRVM1 => XLSRVM
!
XINPRC1 => XINPRC
XACPRC1 => XACPRC
XINDEP1 => XINDEP
XACDEP1 => XACDEP
XINPRR1 => XINPRR
XINPRR3D1 => XINPRR3D
XEVAP3D1 => XEVAP3D
XACPRR1 => XACPRR
XINPRS1 => XINPRS
XACPRS1 => XACPRS
XINPRG1 => XINPRG
XACPRG1 => XACPRG
XINPRH1 => XINPRH
XACPRH1 => XACPRH
!
IF (LUV_FLX) THEN
  XVU_FLUX_M1  => XVU_FLUX_M
END IF
IF (LTH_FLX) THEN
  XVTH_FLUX_M1 => XVTH_FLUX_M
  XWTH_FLUX_M1 => XWTH_FLUX_M1
END IF
!
END SUBROUTINE SET_POINTERS_TO_MODEL1
!
!-------------------------------------------------------------------------------
!
END PROGRAM SPAWNING  
