!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODI_INI_SIZE_n
!     #################
!
INTERFACE
!
SUBROUTINE INI_SIZE_n(KMI,HLUOUT,TPINIFILE,HINIFILEPGD)
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
INTEGER,            INTENT(IN)    :: KMI          !Model Index
CHARACTER (LEN=*),  INTENT(IN)    :: HLUOUT       !Name for output-listing of nested models
TYPE(TFILEDATA),    INTENT(IN)    :: TPINIFILE    !Initial file
CHARACTER (LEN=*),  INTENT(IN)    :: HINIFILEPGD
!
       END SUBROUTINE INI_SIZE_n
!
END INTERFACE
!
END MODULE MODI_INI_SIZE_n
!-----------------------------------------------------------------
!     #######################################################
      SUBROUTINE INI_SIZE_n(KMI,HLUOUT,TPINIFILE,HINIFILEPGD)
!     #######################################################
!
!!
!!****  *INI_SIZE_n* - routine to initialize the sizes ratio positions of nested model _n
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the sizes ratio positions
!     of the nested model _n.
!
!!**  METHOD
!!    ------
!!      The first part of the initialization of the model _n is performed as 
!!    follows :
!!       - The logical unit number associated to output_listing file HLUOUT is 
!!    retrieved and module MODD_LUNIT_n is initialized.
!!       -  Then the description of the segment to perform for the model _n is
!!    retrieved : 
!!            * If there is more than one model, the variables in EXSEG file 
!!    which have been updated in INI_SEG (WRITE_DESFM)  are read in order
!!    to initialize properly the corresponding variables in modules MODD_XXXX_n.
!!    (If there is only one model, the variables in modules MODD_XXXX1   
!!    have been  already properly initialized by the routine READ_EXSEG)
!!            * The kind of geometry (cartesian or spherical geometries) is
!!    also initialized by reading LFIFM file. 
!!       - The dimensions of arrays in initial file are initialized by SET_DIM.
!!   
!!    EXTERNAL
!!    --------
!!      IO_READ_FIELD : to read a LFIFM file
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!        
!!      Module MODD_PARAMETERS : JPHEXT
!!      Module MODD_CONF       : NVERB,CCONF,LCARTESIAN,LTHINSHELL
!!      Module MODD_DYN        : LCORIO  
!!      Module MODD_NESTING    : MY_NAME, DAD_NAME, NDXRATIO_ALL, NDYRATIO_ALL
!!                               NXOR_ALL,NYOR_ALL,NXEND_ALL,NYEND_ALL
!!      Module MODN_LUNIT_n     : contain the namelist NAM_LUNIT_n
!!      Module MODN_CONF_n      : idem...
!!      Module MODN_DYN_n       : NIMAX_ll,NJMAX_ll
!!      Module MODN_ADV_n
!!      Module MODN_DYN_n
!!      Module MODN_PARAM_n
!!      Module MODN_PARAM_RAD_n
!!      Module MODN_PARAM_KAFR_n
!!      Module MODD_LBC_n       : CLBCX,CLBCY ...
!!      Module MODN_TURB_n
!!      Module MODN_CH_MNHC_n
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France * (former ini_modeln)
!!
!!    MODIFICATIONS
!!    -------------
!!             Mar. 31 1999 (Gicquel N.) Part of model initialization necessary 
!!                                       for routine INI_PARA
!!             Apr. 04 2000 (P Jabouille) Halo size and kind of splitting choice
!!             Oct. 10 2001 (I. Mallet)  allow namelists in different orders
!!             Jan. 2004   (V. Masson)  externalization of surface
!!             June 2006   (D. Gazen) _n: no more read of updated var. 
!!             J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,       ONLY: CCONF, LCARTESIAN, NVERB, LTHINSHELL, NHALO, CSPLIT, &
                           L1D, L2D, LPACK
USE MODD_CONFZ,      ONLY: NZ_PROC
USE MODD_DIM_n,      ONLY: NIMAX_ll, NJMAX_ll, NKMAX
USE MODD_DYN,        ONLY: LCORIO
USE MODD_IO_ll,      ONLY: GSMONOPROC, TFILEDATA
USE MODD_LBC_n,      ONLY: CLBCX, CLBCY
USE MODD_LUNIT_n,    ONLY: CLUOUT, CINIFILE, CINIFILEPGD, TLUOUT
USE MODD_NESTING,    ONLY: CMY_NAME, CDAD_NAME, NDAD, NDXRATIO_ALL, NDYRATIO_ALL, &
                           NXOR_ALL, NYOR_ALL, NXEND_ALL,NYEND_ALL
USE MODD_PARAMETERS, ONLY: JPMODELMAX, JPHEXT,JPVEXT 
!
USE MODE_FM,         ONLY: SET_FMPACK_ll
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_ll
USE MODE_MSG
USE MODE_POS
USE MODE_SPLITTINGZ_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,            INTENT(IN)    :: KMI          !Model Index
CHARACTER (LEN=*),  INTENT(IN)    :: HLUOUT       !Name for output-listing of nested models
TYPE(TFILEDATA),    INTENT(IN)    :: TPINIFILE    !Initial file
CHARACTER (LEN=*),  INTENT(IN)    :: HINIFILEPGD
!
!*       0.2   declarations of local variables
!
INTEGER             :: IRESP   ! Return code of FM routines 
INTEGER             :: ILUOUT  ! Logical unit number of output-listing
INTEGER             :: IJPHEXT
!
!-------------------------------------------------------------------------------
!
!*       1.    RETRIEVE LOGICAL UNIT NUMBER AND INITIALIZE MODD_LUNIT_n
!              --------------------------------------------------------
!
ILUOUT = TLUOUT%NLU
CLUOUT = HLUOUT
CINIFILEPGD=HINIFILEPGD
!
!-------------------------------------------------------------------------------
!
!*       2.    RETRIEVE SEGMENT DESCRIPTION
!              ----------------------------
!
!*       2.0   Retrieve DAD_NAME and MY_NAME to check the DAD model identity
!
CALL IO_READ_FIELD(TPINIFILE,'MY_NAME',CMY_NAME(KMI),IRESP)
IF (IRESP /= 0)  THEN
  WRITE(ILUOUT,FMT=9000) 'MY_NAME',IRESP
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SIZE_n','')
END IF
!
CALL IO_READ_FIELD(TPINIFILE,'DAD_NAME',CDAD_NAME(KMI),IRESP)
IF (IRESP /= 0)  THEN
  WRITE(ILUOUT,FMT=9000) 'DAD_NAME',IRESP
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SIZE_n','')
END IF
!
IF ( KMI > 1 ) THEN
  IF ( TRIM(CDAD_NAME(KMI)) /= TRIM(CMY_NAME(NDAD(KMI))) ) THEN
    WRITE(UNIT=ILUOUT,FMT=9005) NDAD(KMI)
    WRITE(ILUOUT,FMT=*) ' THE INITIAL FM-File IS NOT CONSISTANT WITH THE ONE OF THE DAD MODEL!'
    WRITE(UNIT=ILUOUT,FMT=*) 'KMI=',KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'CDAD_NAME(KMI)='//TRIM(CDAD_NAME(KMI))
    WRITE(UNIT=ILUOUT,FMT=*) 'CMY_NAME(NDAD(KMI))='//TRIM(CMY_NAME(NDAD(KMI)))
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SIZE_n','')
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.    INITIALIZE MODEL DIMENSIONS
!              ---------------------------
!
!
!*       3.1  Read dimensions in initial file and initialize  subdomain 
!             dimensions and parallel variables
!
CALL IO_READ_FIELD(TPINIFILE,'IMAX',  NIMAX_ll)
CALL IO_READ_FIELD(TPINIFILE,'JMAX',  NJMAX_ll)
CALL IO_READ_FIELD(TPINIFILE,'KMAX',  NKMAX)
CALL IO_READ_FIELD(TPINIFILE,'JPHEXT',IJPHEXT)
!
IF ( IJPHEXT .NE. JPHEXT ) THEN
   WRITE(ILUOUT,FMT=*) ' INI_SIZE_N : JPHEXT in namelist NAM_CONF ( or default or .des value )&
        & JPHEXT=',JPHEXT
   WRITE(ILUOUT,FMT=*)' different from LFI file=',TPINIFILE%CNAME ,' value JPHEXT=',IJPHEXT
   WRITE(ILUOUT,FMT=*) '-> JOB ABORTED'
   CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SIZE_n','')
   !WRITE(NLUOUT,FMT=*) ' JPHEXT HAS BEEN SET TO ', IJPHEXT
   !IJPHEXT = JPHEXT
END IF
!
IF (KMI == 1) THEN   ! special initialisation for the outer model
  NDXRATIO_ALL(KMI) = 1
  NDYRATIO_ALL(KMI) =  1
  NXOR_ALL(KMI) =  1
  NYOR_ALL(KMI) =  1
  NXEND_ALL(KMI) = NIMAX_ll + 2*JPHEXT
  NYEND_ALL(KMI) = NJMAX_ll + 2*JPHEXT
  NDAD(1) = 0
  CALL SET_SPLITTING_ll(CSPLIT)
!!$  CALL SET_NZ_PROC_ll(NZ_PROC)
  CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT, NHALO)
  CALL SET_DAD0_ll()
  CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
ENDIF
!
!   read the nested model location in its father's grid
!   and compute the coordinates of the corner points
IF (LEN_TRIM(CDAD_NAME(KMI))>0) THEN
  CALL IO_READ_FIELD(TPINIFILE,'DXRATIO',NDXRATIO_ALL(KMI))
  CALL IO_READ_FIELD(TPINIFILE,'DYRATIO',NDYRATIO_ALL(KMI))
  CALL IO_READ_FIELD(TPINIFILE,'XOR',NXOR_ALL(KMI))
  CALL IO_READ_FIELD(TPINIFILE,'YOR',NYOR_ALL(KMI))
  NXEND_ALL(KMI)=NXOR_ALL(KMI)-1 + NIMAX_ll/NDXRATIO_ALL(KMI) +2*JPHEXT
  NYEND_ALL(KMI)=NYOR_ALL(KMI)-1 + NJMAX_ll/NDYRATIO_ALL(KMI) +2*JPHEXT
ELSE
  NDXRATIO_ALL(KMI)=1
  NDYRATIO_ALL(KMI)=1
END IF
!
CALL SET_LBX_ll(CLBCX(1), KMI)
CALL SET_LBY_ll(CLBCY(1), KMI)
CALL SET_XRATIO_ll(NDXRATIO_ALL(KMI), KMI)
CALL SET_YRATIO_ll(NDYRATIO_ALL(KMI), KMI)
CALL SET_XOR_ll(NXOR_ALL(KMI), KMI)
CALL SET_XEND_ll(NXEND_ALL(KMI), KMI)
CALL SET_YOR_ll(NYOR_ALL(KMI), KMI)
CALL SET_YEND_ll(NYEND_ALL(KMI), KMI)
CALL SET_DAD_ll(NDAD(KMI), KMI)
!
IF (KMI == 1) NDAD(KMI)=1  ! return to mesonh meaning
!  
IF (NVERB >= 5) THEN
  WRITE(UNIT=ILUOUT,    &
  FMT="(' DIMENSIONS INITIALIZED BY INI_SIZE_n :',/,'NIMAX =',I5,' NJMAX =',I5,' NKMAX =',I5)")  & 
  NIMAX_ll,NJMAX_ll,NKMAX
END IF
!
!*       3.2  Set the configuration (MODD_CONF)
!
IF (KMI == 1) THEN   
!
  IF( (NIMAX_ll == 1).AND.(NJMAX_ll == 1) .AND. .NOT.L1D) THEN
    L1D=.TRUE.
    WRITE(UNIT=ILUOUT,FMT=9002) KMI
    WRITE(ILUOUT,FMT=*) 'THIS IS A 1D CONFIGURATION : L1D is set to T'
  ENDIF
  IF (L1D .AND. .NOT.GSMONOPROC) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SIZE_n','this is a 1D simulation: it has to be performed in monoprocess mode')
  ENDIF
!
  IF( (NIMAX_ll /= 1).AND.(NJMAX_ll == 1) .AND. .NOT.L2D) THEN
    L2D=.TRUE.
    WRITE(UNIT=ILUOUT,FMT=9002) KMI
    WRITE(ILUOUT,FMT=*) 'THIS IS A 2D CONFIGURATION : L2D is set to T'
  ENDIF
  IF (L2D .AND. .NOT.GSMONOPROC) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SIZE_n','this is a 2D simulation: it has to be performed in monoprocess mode')
  ENDIF
!
  CALL SET_FMPACK_ll(L1D,L2D,LPACK)
!
END IF
!-------------------------------------------------------------------------------
!
!*       4.   FORMATS
!              -------
!
9000  FORMAT(/,'FATAL ERROR IN INI_SIZE_n: pb to read ',A16,' IRESP=',I2)
9002  FORMAT(/,'WARNING IN READ_EXSEG FOR MODEL ', I2,' : ',/, &
             '----------------------------------' )
9003  FORMAT(/,'FATAL ERROR IN INI_SIZE_n FOR MODEL ', I2,' : ',/, &
             '--------------------------------------' )
9005  FORMAT(/,'FATAL ERROR IN INI_SIZE_n FOR MODEL_n  AND ITS DAD MODEL ', I2,' : ',/, &
               '--------------------------------------------------------' )
!
!-------------------------------------------------------------------------------
END SUBROUTINE INI_SIZE_n
