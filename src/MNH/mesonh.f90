!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##############
      PROGRAM MESONH
!     ##############
!
!!****  *MESONH * -general monitor of the model
!!
!!    PURPOSE
!!    -------
!! 
!!      This program is the general monitor of the model. Firstly, it calls the 
!!    subroutine INIT, which performs the sequential initialization of the  
!!    nested models. Then, the program calls the temporal loops of all the 
!!    models, by calling a recursive function which make the temporal nesting
!!    of the different nested models. 
!!    
!!**  METHOD
!!    ------
!!
!!      The initialization is a sequentially performed together with the
!!    temporal loop of all the nested models. The spatial nesting can be
!!    performed in an arbitrary way, the only constrainst is for the first model
!!    which must contain all the others. For the moment, only 8 models can be
!!    runned at the time and the imbriquation level can also go to this upper
!!    value.
!!    
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!       subroutine INIT  : performs the sequential initialization of the nested 
!!                         models
!!       subroutine MODEL: choose the right MODELn to be called
!! 
!!       subroutine KID_MODEL: recursive function which calls the kid models
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    MODD_CONF:     NMODEL, CPROGRAM
!!    MODD_CONF_n:   CSTORAGE_TYPE
!!                   
!!
!!    REFERENCE
!!    ---------
!!
!!    NONE
!!
!!    AUTHOR
!!    ------
!!
!!       J. STEIN  * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    20/10/94 
!!      J.Stein                     08/12/94  clean the code and add the modules
!!      J.Stein and P.Jabouille     30/04/96  add the storage_type
!!      J.P.Lafore                  11/07/96  multi-tasking introduction for nesting
!!      J.P.Lafore                  01/08/96  events implementation
!!      J.P.Lafore                  17/11/97  events modification for two-way nesting
!!      J.Stein                     08/07/98  sequential form for the nesting
!!      J.Stein                     08/04/99  general case of the sequential form
!!      V. Masson                   15/03/99  MASDEV number and PROGRAM name
!!      J.P. Chaboureau             15/03/04  loop limited to 100000 iterations
!!                                            remplaced by infinite loop
!!      J.Escobar                 19/03/2008  rename INIT to INIT_MNH --> grib problem
!!      J.Escobar                  6/11/2014  remove test on LCHECK otherwise never call MPPDB_INIT
!!      J.Escobar                 15/09/2015  WENO5 & JPHEXT <> 1
!!      G.Delautier                  06/2016  phasage surfex 8
!!      J. Pianezze               01/08/2016  add sfxoasis coupling functions
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
#ifdef CPLOASIS
  USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD, IP
  USE MODD_DYN_n, ONLY : XTSTEP
  USE MODD_SFX_OASIS, ONLY : LOASIS, LOASIS_GRID
#endif
!
USE MODD_CONF
USE MODD_NESTING
USE MODD_CONF_n
USE MODD_IO_ll, ONLY: NIO_VERB,NVERB_DEBUG
!
USE MODI_MODEL_n
USE MODI_KID_MODEL
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY : IO_FILE_PRINT_LIST
USE MODE_MODELN_HANDLER
!
USE MODI_VERSION
USE MODI_INIT_MNH
USE MODD_MNH_SURFEX_n
!
#ifdef CPLOASIS
  USE MODI_SFX_OASIS_INIT
  USE MODI_MNH_OASIS_GRID
  USE MODI_MNH_OASIS_DEFINE
  USE MODI_SFX_OASIS_END
#endif
!
USE MODE_MPPDB
!
IMPLICIT NONE 
!
!*       0.1    declarations of local variables
!
INTEGER       :: JMODEL                       ! loop index 
INTEGER       :: ITEMP_MODEL1                 ! loop increment 
LOGICAL       :: GEXIT                        ! flag for the end of the
                                              ! temporal loop       
INTEGER       :: IINFO_ll                     ! return code of // routines
!
#ifdef CPLOASIS
  CHARACTER(LEN=28)  :: CNAMELIST
  LOGICAL            :: L_MASTER
#endif
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION
!              --------------
! Switch to model 1 variables
#ifndef CPLOASIS
CALL MPPDB_INIT()
#endif
!
CALL GOTO_MODEL(1,ONOFIELDLIST=.TRUE.)
!
#ifdef CPLOASIS
  CNAMELIST='EXSEG1.nam'
  CALL SFX_OASIS_INIT(CNAMELIST, NMNH_COMM_WORLD)
#endif
!
CALL INITIO_ll()
!
CALL VERSION
CPROGRAM='MESONH'
!
CALL INIT_MNH
!
!
GEXIT=.FALSE.
!
!
!*         1.1   INITIALIZATION GRID OASIS
!                -------------------------
!
!
#ifdef CPLOASIS
IF(IP==1) THEN
  L_MASTER=.TRUE.
ELSE
  L_MASTER=.FALSE.
END IF
!
IF (LOASIS_GRID) THEN
  CALL MNH_OASIS_GRID(L_MASTER,NMNH_COMM_WORLD)
ENDIF
#endif
!
!
!*         1.2   INITIALIZATION PARTITION OASIS
!                ------------------------------
!
#ifdef CPLOASIS
IF (LOASIS) THEN
  CALL MNH_OASIS_DEFINE(CPROGRAM,IP)
END IF
#endif
!
!-------------------------------------------------------------------------------
!
!*       2.    TEMPORAL LOOP 
!              -------------
!
DO JMODEL=1,NMODEL
  CALL GO_TOMODEL_ll(JMODEL,IINFO_ll)
  CALL GOTO_MODEL(JMODEL)
  CSTORAGE_TYPE='TT'
  CALL MODEL_n(1,GEXIT)
END DO
!
IF(GEXIT) THEN
   !callabortstop
  CALL ABORT
  STOP
ENDIF
!
ITEMP_MODEL1=1
DO
  ITEMP_MODEL1=ITEMP_MODEL1+1
  !
  CALL GO_TOMODEL_ll(1,IINFO_ll)
  CALL GOTO_MODEL(1)
  CALL MODEL_n(ITEMP_MODEL1,GEXIT)
  !
  CALL KID_MODEL(1,ITEMP_MODEL1,GEXIT)
  !
  IF(GEXIT) EXIT
  !
END DO
!
IF(NIO_VERB>=NVERB_DEBUG) CALL IO_FILE_PRINT_LIST()
!
!-------------------------------------------------------------------------------
!
!*       3.    FINALIZE THE PARALLEL SESSION
!              -----------------------------
!
IF (LCHECK) THEN
  CALL MPPDB_BARRIER()
ELSE
  CALL END_PARA_ll(IINFO_ll)
#ifdef CPLOASIS
IF (LOASIS) THEN
  CALL SFX_OASIS_END
END IF
#endif
END IF
!
!
CALL SURFEX_DEALLO_LIST
!
!-------------------------------------------------------------------------------
!
!callabortstop
!CALL ABORT
STOP
!
END PROGRAM MESONH  
