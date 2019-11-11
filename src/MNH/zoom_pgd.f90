!MNH_LIC Copyright 2005-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################
      PROGRAM ZOOM_PGD
!     ################
!!
!!    PURPOSE
!!    -------
!!   This program zooms the physiographic data fields.
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     march 2005
!!    10/10/2011  J.Escobar call INI_PARAZ_ll
!!    30/03/2012  S.Bielli  Add NAM_NCOUT
!!  06/2016     (G.Delautier) phasage surfex 8
!!    08/07/2016  P.Wautelet Removed MNH_NCWRIT define
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_CONF,   ONLY : CPROGRAM, L1D, L2D, LPACK
USE MODD_IO_ll,  ONLY:  NIO_VERB,NVERB_DEBUG,TFILE_OUTPUTLISTING,TFILEDATA
USE MODD_LUNIT,  ONLY : CLUOUT0, TLUOUT0, TOUTDATAFILE
USE MODD_PARAMETERS, ONLY : XUNDEF, NUNDEF, JPVEXT, JPHEXT, JPMODELMAX
USE MODD_PARAM_n,     ONLY : CSURF
USE MODD_DIM_n,       ONLY : NIMAX, NJMAX
USE MODD_CONF_n,   ONLY : CSTORAGE_TYPE
!
USE MODE_POS
USE MODE_FM
USE MODE_FMWRIT
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY : IO_FILE_ADD2LIST, IO_FILE_PRINT_LIST
USE MODE_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
USE MODI_READ_HGRID
USE MODI_WRITE_HGRID
USE MODI_SET_SUBDOMAIN
!JUANZ
USE MODE_SPLITTINGZ_ll
!JUANZ
!
USE MODI_VERSION
USE MODI_READ_ALL_NAMELISTS
USE MODI_ZOOM_PGD_SURF_ATM
USE MODI_WRITE_PGD_SURF_ATM_N
USE MODD_MNH_SURFEX_n
!
USE MODN_CONFIO, ONLY : NAM_CONFIO
!
IMPLICIT NONE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IRESP    ! return code for I/O
INTEGER :: ILUOUT0
INTEGER :: ILUNAM
INTEGER :: IINFO_ll
CHARACTER(LEN=28)  :: CPGDFILE              ! name of the PGD file
CHARACTER(LEN=28)  :: YZOOMFILE             ! name of the output file
CHARACTER(LEN=2)   :: YZOOMNBR
CHARACTER(LEN=28)  :: YMY_NAME,YDAD_NAME
CHARACTER(LEN=28)  :: YPGDFILE
CHARACTER(LEN=2)   :: YSTORAGE_TYPE
LOGICAL :: GFOUND
INTEGER            :: IXOR_DAD,IYOR_DAD  ! compared to Dad file, if any
INTEGER            :: IXOR,IYOR          ! given or computed
INTEGER            :: IDXRATIO,IDYRATIO
TYPE(TFILEDATA),POINTER :: TZNMLFILE  => NULL()
TYPE(TFILEDATA),POINTER :: TZPGDFILE  => NULL()
TYPE(TFILEDATA),POINTER :: TZZOOMFILE => NULL()
!
REAL,  DIMENSION(:,:), ALLOCATABLE :: ZZS1,ZZSMT1,ZZS2,ZZSMT2
!
NAMELIST/NAM_PGDFILE/CPGDFILE,YZOOMFILE,YZOOMNBR
!------------------------------------------------------------------------------
!
CALL GOTO_MODEL(1)
CALL VERSION
CPROGRAM='ZOOMPG'
CSTORAGE_TYPE = 'PG'
!
CALL INI_CST
! 
!
!*    1.      Set default names and parallelized I/O
!             --------------------------------------
!
CALL INITIO_ll()
!
CLUOUT0='OUTPUT_LISTING0'                    ! name of the output-listing
CALL IO_FILE_ADD2LIST(TLUOUT0,'OUTPUT_LISTING0','OUTPUTLISTING','WRITE')
CALL IO_FILE_OPEN_ll(TLUOUT0)
TFILE_OUTPUTLISTING => TLUOUT0
ILUOUT0=TLUOUT0%NLU
!
CALL IO_FILE_ADD2LIST(TZNMLFILE,'PRE_ZOOM1.nam','NML','READ')
CALL IO_FILE_OPEN_ll(TZNMLFILE)
ILUNAM = TZNMLFILE%NLU
!
CPGDFILE  = 'PGDFILE'                         ! name of the input file
YZOOMFILE = ''
YZOOMNBR  = '00'
CALL POSNAM(ILUNAM,'NAM_PGDFILE',GFOUND,ILUOUT0)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_PGDFILE)
CALL POSNAM(ILUNAM,'NAM_CONFIO',GFOUND,ILUOUT0)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_CONFIO)
CALL SET_CONFIO_ll()
!
!------------------------------------------------------------------------------
!
!*    2.      ZOOM OF PGD DOMAIN
!             ------------------
!
!*    2.1     Open PGD file
!             -------------
!
CALL IO_FILE_ADD2LIST(TZPGDFILE,TRIM(CPGDFILE),'UNKNOWN','READ',KLFINPRAR=INT(1,KIND=LFI_INT),KLFITYPE=2,KLFIVERB=5)
CALL IO_FILE_OPEN_ll(TZPGDFILE)
!
!*    2.2     Reading of initial grid
!             -----------------------
!
CALL READ_HGRID(1,TZPGDFILE,YMY_NAME,YDAD_NAME,YSTORAGE_TYPE)
!
! NIMAX, NJMAX: size of input domain
ALLOCATE(ZZS1  (NIMAX+2*JPHEXT,NJMAX+2*JPHEXT))
ALLOCATE(ZZSMT1(NIMAX+2*JPHEXT,NJMAX+2*JPHEXT))
CALL IO_READ_FIELD(TZPGDFILE,'ZS',ZZS1)
CALL IO_READ_FIELD(TZPGDFILE,'ZSMT',ZZSMT1)
!
!*    2.3     Define subdomain
!             ----------------
!
CALL SET_SUBDOMAIN(TZNMLFILE,TZPGDFILE,IXOR_DAD,IYOR_DAD,IXOR,IYOR,IDXRATIO,IDYRATIO)
!
CALL IO_FILE_CLOSE_ll(TZNMLFILE)
!
! NIMAX, NJMAX: size of output domain
!
CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
CALL SET_DIM_ll(NIMAX, NJMAX, 1)
CALL SET_LBX_ll('OPEN',1)
CALL SET_LBY_ll('OPEN', 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(NIMAX+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(NJMAX+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
!JUANZ CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!
!
!*    2.4     Writing of final grid
!             ---------------------
!
IF ( (LEN_TRIM(YZOOMFILE) == 0) .OR. (ADJUSTL(YZOOMFILE) == ADJUSTL(CPGDFILE)) ) THEN
  YZOOMFILE=ADJUSTL(ADJUSTR(CPGDFILE)//'.z'//ADJUSTL(YZOOMNBR))
END IF
!
CALL IO_FILE_ADD2LIST(TZZOOMFILE,TRIM(YZOOMFILE),'ZOOMPGD','WRITE',KLFINPRAR=INT(1,KIND=LFI_INT),KLFITYPE=1,KLFIVERB=5)
!PW: TODO: points to dad file (if existing) ! TZZOOMFILE%TDADFILE =>
!
CALL IO_FILE_OPEN_ll(TZZOOMFILE)
CALL WRITE_HGRID(1,TZZOOMFILE)
!
!*    2.5     Preparation of surface physiographic fields
!             -------------------------------------------
!
CALL IO_READ_FIELD(TZPGDFILE,'SURF',CSURF)
!
!
IF (CSURF=='EXTE') THEN
  CALL SURFEX_ALLOC_LIST(1)
  YSURF_CUR => YSURF_LIST(1)
  CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)      
  YPGDFILE   = CPGDFILE
  CPGDFILE   = YZOOMFILE
  TOUTDATAFILE => TZZOOMFILE
  CALL GOTO_SURFEX(1)
  CALL ZOOM_PGD_SURF_ATM(YSURF_CUR,'MESONH',YPGDFILE,'MESONH',YZOOMFILE,'MESONH')
!
!*    2.6     Writes the physiographic fields
!             -------------------------------
!
  CALL WRITE_PGD_SURF_ATM_n(YSURF_CUR,'MESONH')
ELSE
  ALLOCATE(ZZS2(NIMAX+2*JPHEXT,NJMAX+2*JPHEXT))
  ZZS2(:,:)=ZZS1(IXOR:IXOR+NIMAX+2*JPHEXT-1,IYOR:IYOR+NJMAX+2*JPHEXT-1)
  CALL IO_WRITE_FIELD(TZZOOMFILE,'ZS',ZZS2)
END IF
!
ALLOCATE(ZZSMT2(NIMAX+2*JPHEXT,NJMAX+2*JPHEXT))
ZZSMT2(:,:)=ZZSMT1(IXOR:IXOR+NIMAX+2*JPHEXT-1,IYOR:IYOR+NJMAX+2*JPHEXT-1)
CALL IO_WRITE_FIELD(TZZOOMFILE,'ZSMT',ZZSMT2)
!
!*    2.7     Write configuration variables in the output file
!             ------------------------------------------------
!
CALL IO_WRITE_HEADER(TZZOOMFILE)
CALL IO_WRITE_FIELD(TZZOOMFILE,'DXRATIO',IDXRATIO)
CALL IO_WRITE_FIELD(TZZOOMFILE,'DYRATIO',IDYRATIO)
CALL IO_WRITE_FIELD(TZZOOMFILE,'XOR',    IXOR_DAD)
CALL IO_WRITE_FIELD(TZZOOMFILE,'YOR',    IYOR_DAD)
CALL IO_WRITE_FIELD(TZZOOMFILE,'L1D',    L1D)
CALL IO_WRITE_FIELD(TZZOOMFILE,'L2D',    L2D)
CALL IO_WRITE_FIELD(TZZOOMFILE,'PACK',   LPACK)
CALL IO_WRITE_FIELD(TZZOOMFILE,'SURF',   CSURF)
CALL IO_FILE_CLOSE_ll(TZZOOMFILE)
!
!*    2.8     Shift to new PGD file
!             ---------------------
!
CPGDFILE = YZOOMFILE
!                       
!------------------------------------------------------------------------------
!
!*    3.     CLOSE PARALLELIZED I/O
!            ----------------------
!
CALL IO_FILE_CLOSE_ll(TZPGDFILE)
!
IF(NIO_VERB>=NVERB_DEBUG) CALL IO_FILE_PRINT_LIST()
!
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*) '***************************'
WRITE(ILUOUT0,*) '* ZOOM_PGD ends correctly *'
WRITE(ILUOUT0,*) '***************************'
!
CALL IO_FILE_CLOSE_ll(TLUOUT0)
!
CALL END_PARA_ll(IINFO_ll)

IF (CSURF=='EXTE') CALL SURFEX_DEALLO_LIST
!
!-------------------------------------------------------------------------------
!
END PROGRAM ZOOM_PGD
