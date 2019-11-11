!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:27
!-----------------------------------------------------------------
!     #################
      MODULE MODI_INI_SIZE_n
!     #################
!
INTERFACE
!
       SUBROUTINE INI_SIZE_n(KMI,HLUOUT,HINIFILE)
!
       INTEGER, INTENT(IN)              :: KMI      ! Model Index 
       CHARACTER (LEN=*), INTENT(IN)    :: HLUOUT   ! name for output-listing
       !  of nested models
       CHARACTER (LEN=*),  INTENT(IN)   :: HINIFILE ! name of
                                             ! the initial file
!
       END SUBROUTINE INI_SIZE_n
!
END INTERFACE
!
END MODULE MODI_INI_SIZE_n
!-----------------------------------------------------------------
!     ##########################################
      SUBROUTINE INI_SIZE_n(KMI,HLUOUT,HINIFILE)
!     ##########################################
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
!!      FMLOOK_ll   : to retrieve a logical unit number associated with a file 
!!      FMREAD      : to read a LFIFM file
!!      FMCLOS      : to close a FM-file
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
!!      Module MODN_PARAM_CONVECT_n
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX, JPHEXT,JPVEXT 
USE MODD_CONF, ONLY: CCONF, LCARTESIAN, NVERB, LTHINSHELL, NHALO, CSPLIT, &
                     NZ_PROC, L1D, L2D, LPACK
USE MODD_DYN, ONLY: LCORIO
USE MODD_NESTING, ONLY: CMY_NAME, CDAD_NAME, NDAD, NDXRATIO_ALL, NDYRATIO_ALL, &
                        NXOR_ALL, NYOR_ALL, NXEND_ALL,NYEND_ALL
USE MODD_DIM_n, ONLY: NIMAX_ll, NJMAX_ll, NKMAX
USE MODD_LBC_n, ONLY: CLBCX, CLBCY
USE MODD_LUNIT_n, ONLY: CLUOUT, CINIFILE
USE MODD_IO_ll,   ONLY : GSMONOPROC
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_FMREAD
USE MODE_FM
USE MODE_POS
!
!JUAN
USE MODE_SPLITTINGZ_ll
!JUAN
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER, INTENT(IN)              :: KMI      ! Model Index 
CHARACTER (LEN=*), INTENT(IN)    :: HLUOUT   ! name for output-listing
                                             !  of nested models
CHARACTER (LEN=*),  INTENT(IN)   :: HINIFILE ! name of
                                             ! the initial file
!
!*       0.2   declarations of local variables
!
INTEGER             :: IRESP   ! Return code of FM routines 
INTEGER             :: ILUOUT  ! Logical unit number of output-listing
CHARACTER(LEN=2)    :: YDIR    ! Type  of the data field in LFIFM file
INTEGER             :: IGRID   ! C-grid indicator in LFIFM file 
INTEGER             :: ILENCH  ! Length of comment string in LFIFM file
CHARACTER (LEN=100) :: YCOMMENT! comment string in LFIFM file
CHARACTER (LEN=16)  :: YRECFM  ! Name of the desired field in LFIFM file
!
!-------------------------------------------------------------------------------
!
!*       1.    RETRIEVE LOGICAL UNIT NUMBER AND INITIALIZE MODD_LUNIT_n
!              --------------------------------------------------------
!
CALL FMLOOK_ll(HLUOUT,HLUOUT,ILUOUT,IRESP)
CLUOUT = HLUOUT
CINIFILE=HINIFILE
!
!-------------------------------------------------------------------------------
!
!*       2.    RETRIEVE SEGMENT DESCRIPTION
!              ----------------------------
!
!*       2.0   Retrieve DAD_NAME and MY_NAME to check the DAD model identity
!
YRECFM = 'MY_NAME'
YDIR='--'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,CMY_NAME(KMI),IGRID,ILENCH,YCOMMENT,IRESP)
IF (IRESP /= 0)  THEN
  WRITE(ILUOUT,FMT=9000) YRECFM,IRESP
  STOP
END IF
!
YRECFM = 'DAD_NAME'
YDIR='--'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,CDAD_NAME(KMI),IGRID,ILENCH,YCOMMENT,IRESP)
IF (IRESP /= 0)  THEN
  WRITE(ILUOUT,FMT=9000) YRECFM,IRESP
  STOP
END IF
!
IF ( KMI > 1 ) THEN
  IF ( TRIM(CDAD_NAME(KMI)) /= TRIM(CMY_NAME(NDAD(KMI))) ) THEN
    WRITE(UNIT=ILUOUT,FMT=9005) NDAD(KMI)
    WRITE(ILUOUT,FMT=*) ' THE INITIAL FM-File IS NOT CONSISTANT WITH THE ONE OF THE DAD MODEL!'
    WRITE(UNIT=ILUOUT,FMT=*) 'KMI=',KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'CDAD_NAME(KMI)='//TRIM(CDAD_NAME(KMI))
    WRITE(UNIT=ILUOUT,FMT=*) 'CMY_NAME(NDAD(KMI))='//TRIM(CMY_NAME(NDAD(KMI)))
    STOP
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
YRECFM='IMAX'
YDIR='--'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,NIMAX_ll,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='JMAX'
YDIR='--'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,NJMAX_ll,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='KMAX'
YDIR='--'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,YDIR,NKMAX,IGRID,ILENCH,YCOMMENT,IRESP)
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
  CALL SET_NZ_PROC_ll(NZ_PROC)
  CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT, NHALO)
  CALL SET_DAD0_ll()
  CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
ENDIF
!
!   read the nested model location in its father's grid
!   and compute the coordinates of the corner points
IF (LEN_TRIM(CDAD_NAME(KMI))>0) THEN
  YDIR='--'
  CALL FMREAD(HINIFILE,'DXRATIO',HLUOUT,YDIR,NDXRATIO_ALL(KMI),IGRID,ILENCH,YCOMMENT,IRESP)
  CALL FMREAD(HINIFILE,'DYRATIO',HLUOUT,YDIR,NDYRATIO_ALL(KMI),IGRID,ILENCH,YCOMMENT,IRESP)
  CALL FMREAD(HINIFILE,'XOR',HLUOUT,YDIR,NXOR_ALL(KMI),IGRID,ILENCH,YCOMMENT,IRESP)
  CALL FMREAD(HINIFILE,'YOR',HLUOUT,YDIR,NYOR_ALL(KMI),IGRID,ILENCH,YCOMMENT,IRESP)
  NXEND_ALL(KMI)=NXOR_ALL(KMI)+NIMAX_ll/NDXRATIO_ALL(KMI)+JPHEXT
  NYEND_ALL(KMI)=NYOR_ALL(KMI)+NJMAX_ll/NDYRATIO_ALL(KMI)+JPHEXT
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
    WRITE(ILUOUT,FMT=*) 'THIS IS A 1D SIMULATION : IT  HAS TO BE &
                      & PERFORMED WITH MONOPROCESSOR MODE'
    STOP
  ENDIF
!
  IF( (NIMAX_ll /= 1).AND.(NJMAX_ll == 1) .AND. .NOT.L2D) THEN
    L2D=.TRUE.
    WRITE(UNIT=ILUOUT,FMT=9002) KMI
    WRITE(ILUOUT,FMT=*) 'THIS IS A 2D CONFIGURATION : L2D is set to T'
  ENDIF
  IF (L2D .AND. .NOT.GSMONOPROC) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(ILUOUT,FMT=*) 'THIS IS A 2D SIMULATION : IT  HAS TO BE &
                      & PERFORMED WITH MONOPROCESSOR MODE'
    STOP
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
