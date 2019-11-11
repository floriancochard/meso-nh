!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  D.Gazen   : avril 2016 change error message
!  P. Wautelet : may 2016: use NetCDF Fortran module
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  Philippe Wautelet: 29/10/2018: better detection of older MNH version numbers
!  Philippe Wautelet: 13/12/2018: moved some operations to new mode_io_*_nc4 modules
!  Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN + move management
!                                 of NNCID and NLFIFLU to the nc4 and lfi subroutines
!  Philippe Wautelet: 21/01/2019: add LIO_ALLOW_NO_BACKUP and LIO_NO_WRITE to modd_io_ll
!                                 to allow to disable writes (for bench purposes)
!-----------------------------------------------------------------

MODULE MODE_FM
USE MODD_MPIF

USE MODE_MSG

IMPLICIT NONE 

PRIVATE 

PUBLIC SET_FMPACK_ll
PUBLIC IO_FILE_OPEN_ll, IO_FILE_CLOSE_ll

CONTAINS 

SUBROUTINE SET_FMPACK_ll(O1D,O2D,OPACK)
USE MODD_IO_ll, ONLY : LPACK,L1D,L2D
!JUAN
USE MODD_VAR_ll, ONLY : IP
!JUAN

IMPLICIT NONE 

LOGICAL, INTENT(IN) :: O1D,O2D,OPACK

LPACK = OPACK
L1D   = O1D
L2D   = O2D

IF ( IP .EQ. 1 ) PRINT *,'INIT L1D,L2D,LPACK = ',L1D,L2D,LPACK

END SUBROUTINE SET_FMPACK_ll

SUBROUTINE IO_FILE_OPEN_ll(TPFILE,KRESP,OPARALLELIO,HPOSITION,HSTATUS,HPROGRAM_ORIG)
!
USE MODD_CONF,  ONLY: CPROGRAM
USE MODD_IO_ll, ONLY: LIO_NO_WRITE, TFILEDATA
USE MODE_FMREAD
USE MODE_IO_ll, ONLY : OPEN_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST,IO_FILE_FIND_BYNAME
!
TYPE(TFILEDATA),POINTER,INTENT(INOUT)         :: TPFILE ! File structure
INTEGER,                INTENT(OUT), OPTIONAL :: KRESP  ! Return code
LOGICAL,                INTENT(IN),  OPTIONAL :: OPARALLELIO
CHARACTER(LEN=*),       INTENT(IN),  OPTIONAL :: HPOSITION
CHARACTER(LEN=*),       INTENT(IN),  OPTIONAL :: HSTATUS
CHARACTER(LEN=*),       INTENT(IN),  OPTIONAL :: HPROGRAM_ORIG !To emulate a file coming from this program
!
INTEGER :: IRESP
TYPE(TFILEDATA),POINTER :: TZFILE_DES
TYPE(TFILEDATA),POINTER :: TZFILE_DUMMY
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_FILE_OPEN_ll','opening '//TRIM(TPFILE%CNAME)//' for '//TRIM(TPFILE%CMODE)// &
               ' (filetype='//TRIM(TPFILE%CTYPE)//')')
!
IF (.NOT.ASSOCIATED(TPFILE)) CALL PRINT_MSG(NVERB_FATAL,'IO','IO_FILE_OPEN_ll','TPFILE is not associated')
!
IF ( LIO_NO_WRITE .AND. TPFILE%CMODE == 'WRITE' .AND. TPFILE%CTYPE/='OUTPUTLISTING') THEN
  CALL PRINT_MSG(NVERB_WARNING,'IO','IO_FILE_OPEN_ll','opening file '//TRIM(TPFILE%CNAME)//' in write mode but LIO_NO_WRITE is set')
END IF
!
TZFILE_DES   => NULL()
TZFILE_DUMMY => NULL()
!
TPFILE%NOPEN         = TPFILE%NOPEN + 1
TPFILE%NOPEN_CURRENT = TPFILE%NOPEN_CURRENT + 1
!
IF (TPFILE%LOPENED) THEN
  CALL PRINT_MSG(NVERB_INFO,'IO','IO_FILE_OPEN_ll','file '//TRIM(TPFILE%CNAME)//' is already in open state')
  RETURN
END IF
!
TPFILE%LOPENED       = .TRUE.
!
!Check if file is in filelist
CALL IO_FILE_FIND_BYNAME(TRIM(TPFILE%CNAME),TZFILE_DUMMY,IRESP)
IF (IRESP/=0) CALL PRINT_MSG(NVERB_ERROR,'IO','IO_FILE_OPEN_ll','file '//TRIM(TPFILE%CNAME)//' not in filelist')
!
SELECT CASE(TPFILE%CTYPE)
  !Chemistry input files
  CASE('CHEMINPUT')
    CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='FORMATTED',POSITION='REWIND',STATUS='OLD',MODE='GLOBAL')


  !Chemistry tabulation files
  CASE('CHEMTAB')
    CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='FORMATTED',MODE='GLOBAL')


  !GPS files
  CASE('GPS')
    CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='FORMATTED',MODE='SPECIFIC')


  !Meteo files
  CASE('METEO')
   CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='UNFORMATTED',MODE='GLOBAL',RECL=100000000)


  !Namelist files
  CASE('NML')
    CALL OPEN_ll(TPFILE,IOSTAT=IRESP,DELIM='QUOTE',MODE='GLOBAL')


  !OUTPUTLISTING files
  CASE('OUTPUTLISTING')
    CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='FORMATTED',MODE='GLOBAL')


  !SURFACE_DATA files
  CASE('SURFACE_DATA')
    IF (TPFILE%CFORM=='FORMATTED') THEN
      CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM=TPFILE%CFORM,MODE='GLOBAL')
    ELSE IF (TPFILE%CACCESS=='DIRECT') THEN
      CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM=TPFILE%CFORM,ACCESS=TPFILE%CACCESS,RECL=TPFILE%NRECL,MODE='GLOBAL')
    ELSE
      CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM=TPFILE%CFORM,MODE='GLOBAL')
    END IF


  !Text files
  CASE('TXT')
    IF(TPFILE%NRECL>0) THEN
      CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='FORMATTED',POSITION=HPOSITION,STATUS=HSTATUS,RECL=TPFILE%NRECL,MODE='GLOBAL')
    ELSE
      CALL OPEN_ll(TPFILE,IOSTAT=IRESP,FORM='FORMATTED',POSITION=HPOSITION,STATUS=HSTATUS,MODE='GLOBAL')
    END IF


  CASE DEFAULT
    !Do not open '.des' file if OUTPUT
    IF(TPFILE%CTYPE/='OUTPUT' .AND. CPROGRAM/='LFICDF') THEN
      CALL IO_FILE_ADD2LIST(TZFILE_DES,TRIM(TPFILE%CNAME)//'.des','DES',TPFILE%CMODE,TPDATAFILE=TPFILE,OOLD=.TRUE.) !OOLD=T because the file may already be in the list
      CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_FILE_OPEN_ll','OPEN_ll for '//TRIM(TPFILE%CNAME)//'.des')
      CALL OPEN_ll(TZFILE_DES,FORM='FORMATTED',DELIM='QUOTE',IOSTAT=IRESP,RECL=1024*8,OPARALLELIO=OPARALLELIO)
      TZFILE_DES%LOPENED       = .TRUE.
      TZFILE_DES%NOPEN_CURRENT = TZFILE_DES%NOPEN_CURRENT + 1
      TZFILE_DES%NOPEN         = TZFILE_DES%NOPEN + 1
    ENDIF
    !
    CALL FMOPEN_ll(TPFILE,IRESP,OPARALLELIO=OPARALLELIO,HPROGRAM_ORIG=HPROGRAM_ORIG)

END SELECT
!
IF (PRESENT(KRESP)) KRESP = IRESP
!
END SUBROUTINE IO_FILE_OPEN_ll

SUBROUTINE FMOPEN_ll(TPFILE,KRESP,OPARALLELIO,HPROGRAM_ORIG)
USE MODD_IO_ll,               ONLY: TFILEDATA
USE MODE_IO_ll,               ONLY: OPEN_ll,GCONFIO
!JUANZ
USE MODD_CONFZ,ONLY  : NB_PROCIO_R,NB_PROCIO_W
!JUANZ
#if defined(MNH_IOCDF4)
USE MODD_NETCDF, ONLY:IDCDF_KIND
use mode_io_file_nc4, only: io_create_file_nc4, io_open_file_nc4
#endif
use mode_io_file_lfi, only: io_create_file_lfi, io_open_file_lfi

TYPE(TFILEDATA), INTENT(INOUT) :: TPFILE ! File structure
INTEGER,         INTENT(OUT)   :: KRESP  ! return-code
LOGICAL,         INTENT(IN),  OPTIONAL :: OPARALLELIO
CHARACTER(LEN=*),INTENT(IN),  OPTIONAL :: HPROGRAM_ORIG !To emulate a file coming from this program
!
!   Local variables
!
INTEGER               :: IROWF, IRESP
CHARACTER(LEN=7)      :: YACTION ! Action upon the file ('READ' or 'WRITE')
CHARACTER(LEN=8)      :: YRESP
INTEGER               :: IERR
INTEGER               :: INB_PROCIO
LOGICAL               :: GPARALLELIO
LOGICAL               :: GEXIST_LFI, GEXIST_NC4

YACTION = TPFILE%CMODE

CALL PRINT_MSG(NVERB_DEBUG,'IO','FMOPEN_ll','opening '//TRIM(TPFILE%CNAME)//' for '//TRIM(YACTION))

IF ( PRESENT(OPARALLELIO) ) THEN
  GPARALLELIO = OPARALLELIO
ELSE  !par defaut on active les IO paralleles en Z si possible
  GPARALLELIO = .TRUE.
ENDIF

IF (.NOT. GCONFIO) THEN
   PRINT *, 'FMOPEN_ll Aborting... Please, ensure to call SET_CONFIO_ll before &
        &the first FMOPEN_ll call.'
   STOP
END IF

IROWF  = 0
IRESP  = 0

IROWF=LEN_TRIM(TPFILE%CNAME)

IF (IROWF.EQ.0) THEN
  IRESP=-45
  GOTO 1000
ENDIF

 SELECT CASE (YACTION)
 CASE('READ')
    INB_PROCIO = NB_PROCIO_R
 CASE('WRITE')
    INB_PROCIO = NB_PROCIO_W
 END SELECT
CALL OPEN_ll(TPFILE,STATUS="UNKNOWN",MODE='IO_ZSPLIT',IOSTAT=IRESP,     &
             KNB_PROCIO=INB_PROCIO,OPARALLELIO=GPARALLELIO,HPROGRAM_ORIG=HPROGRAM_ORIG)

IF (IRESP /= 0) GOTO 1000

IF (TPFILE%LMASTER) THEN
  ! Proc I/O case
  INQUIRE(FILE=TRIM(TPFILE%CNAME)//'.lfi',EXIST=GEXIST_LFI)
  INQUIRE(FILE=TRIM(TPFILE%CNAME)//'.nc',EXIST=GEXIST_NC4)

  IF (YACTION == 'READ') THEN
    IF (.NOT.GEXIST_LFI .AND. .NOT.GEXIST_NC4) &
      CALL PRINT_MSG(NVERB_FATAL,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)//': no .nc or .lfi file')

    SELECT CASE (TRIM(TPFILE%CFORMAT))
      CASE ('NETCDF4')
        IF (.NOT.GEXIST_NC4 .AND. GEXIST_LFI) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)// &
                         ': .nc file does not exist but .lfi exists -> forced to LFI')
          TPFILE%CFORMAT='LFI'
        END IF
      CASE ('LFI')
        IF (.NOT.GEXIST_LFI .AND. GEXIST_NC4) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)// &
                         ': .lfi file does not exist but .nc exists -> forced to NETCDF4')
          TPFILE%CFORMAT='NETCDF4'
        END IF
      CASE ('LFICDF4')
        IF (GEXIST_NC4) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)// &
                         ': LFICDF4 format is not allowed in READ mode -> forced to NETCDF4')
          TPFILE%CFORMAT='NETCDF4'
        ELSE IF (GEXIST_LFI) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)// &
                         ': LFICDF4 format is not allowed in READ mode -> forced to LFI')
          TPFILE%CFORMAT='LFI'
        END IF
      CASE DEFAULT
        IF (GEXIST_NC4) THEN
          CALL PRINT_MSG(NVERB_ERROR,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)// &
                         ': invalid fileformat (-> forced to NETCDF4 if no abort)')
          TPFILE%CFORMAT='NETCDF4'
        ELSE IF (GEXIST_LFI) THEN
          CALL PRINT_MSG(NVERB_ERROR,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)// &
                         ': invalid fileformat (-> forced to LFI if no abort)')
          TPFILE%CFORMAT='LFI'
        END IF
    END SELECT
  END IF
END IF

#if defined(MNH_IOCDF4)
IF (TPFILE%CFORMAT=='NETCDF4' .OR. TPFILE%CFORMAT=='LFICDF4') THEN
  SELECT CASE (YACTION)
    CASE('READ')
      call io_open_file_nc4(tpfile)
    CASE('WRITE')
      call io_create_file_nc4(TPFILE, hprogram_orig=HPROGRAM_ORIG)
  END SELECT
END IF
#endif

IF (TPFILE%CFORMAT=='LFI' .OR. TPFILE%CFORMAT=='LFICDF4') THEN
  SELECT CASE (YACTION)
    CASE('READ')
      call io_open_file_lfi(tpfile,iresp)
    CASE('WRITE')
      call io_create_file_lfi(tpfile,iresp)
  END SELECT
END IF

! Broadcast ERROR
CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
IF (IRESP /= 0) GOTO 1000


1000 CONTINUE

IF (IRESP.NE.0)  THEN
  WRITE(YRESP,"( I0 )") IRESP
  CALL PRINT_MSG(NVERB_ERROR,'IO','FMOPEN_ll',TRIM(TPFILE%CNAME)//': exit with IRESP='//TRIM(YRESP))
END IF

KRESP=IRESP

END SUBROUTINE FMOPEN_ll

SUBROUTINE IO_FILE_CLOSE_ll(TPFILE,KRESP,OPARALLELIO,HPROGRAM_ORIG)
!
USE MODD_CONF,  ONLY: CPROGRAM
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_IO_ll, ONLY : CLOSE_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
!
TYPE(TFILEDATA),  INTENT(INOUT)         :: TPFILE ! File structure
INTEGER,          INTENT(OUT), OPTIONAL :: KRESP  ! Return code
LOGICAL,          INTENT(IN),  OPTIONAL :: OPARALLELIO
CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: HPROGRAM_ORIG !To emulate a file coming from this program
!
INTEGER                 :: IRESP, JI
TYPE(TFILEDATA),POINTER :: TZFILE_DES
TYPE(TFILEDATA),POINTER :: TZFILE_IOZ
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_FILE_CLOSE_ll','closing '//TRIM(TPFILE%CNAME))
!
IF (.NOT.TPFILE%LOPENED) THEN
  CALL PRINT_MSG(NVERB_ERROR,'IO','IO_FILE_CLOSE_ll','trying to close a file not opened: '//TRIM(TPFILE%CNAME))
  RETURN
ENDIF
!
IF (TPFILE%NOPEN_CURRENT>1) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_FILE_CLOSE_ll',TRIM(TPFILE%CNAME)// &
                 ': decrementing NOPEN_CURRENT (still opened after this call)')
  TPFILE%NOPEN_CURRENT = TPFILE%NOPEN_CURRENT - 1
  TPFILE%NCLOSE        = TPFILE%NCLOSE        + 1
  !
  DO JI = 1,TPFILE%NSUBFILES_IOZ
    TZFILE_IOZ => TPFILE%TFILES_IOZ(JI)%TFILE
    TZFILE_IOZ%NOPEN_CURRENT = TZFILE_IOZ%NOPEN_CURRENT - 1
    TZFILE_IOZ%NCLOSE        = TZFILE_IOZ%NCLOSE        + 1
  END DO
  !
  RETURN
END IF
!
SELECT CASE(TPFILE%CTYPE)
  !Chemistry input files
  CASE('CHEMINPUT')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  !Chemistry tabulation files
  CASE('CHEMTAB')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  !GPS files
  CASE('GPS')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  !Meteo files
  CASE('METEO')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  !Namelist files
  CASE('NML')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  !OUTPUTLISTING files
  CASE('OUTPUTLISTING')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP,OPARALLELIO=.FALSE.)
    !
    TPFILE%NLU = -1


  !SURFACE_DATA files
  CASE('SURFACE_DATA')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  !Text files
  CASE('TXT')
    CALL CLOSE_ll(TPFILE,IOSTAT=IRESP)
    !
    TPFILE%NLU = -1


  CASE DEFAULT
    !Do not close (non-existing) '.des' file if OUTPUT
    IF(TPFILE%CTYPE/='OUTPUT' .AND. CPROGRAM/='LFICDF') THEN
      CALL IO_FILE_FIND_BYNAME(TRIM(TPFILE%CNAME)//'.des',TZFILE_DES,IRESP)
      IF (IRESP/=0) CALL PRINT_MSG(NVERB_ERROR,'IO','IO_FILE_CLOSE_ll','file '//TRIM(TPFILE%CNAME)//'.des not in filelist')
      !
      TZFILE_DES%NOPEN_CURRENT = TZFILE_DES%NOPEN_CURRENT - 1
      TZFILE_DES%NCLOSE        = TZFILE_DES%NCLOSE + 1
      !
      IF (TZFILE_DES%NOPEN_CURRENT==0) THEN
        CALL CLOSE_ll(TZFILE_DES,IOSTAT=IRESP)
        TZFILE_DES%LOPENED = .FALSE.
        TZFILE_DES%NLU     = -1
      END IF
    ENDIF
    !
    CALL FMCLOS_ll(TPFILE,KRESP=IRESP,OPARALLELIO=OPARALLELIO,HPROGRAM_ORIG=HPROGRAM_ORIG)
    !
    DO JI = 1,TPFILE%NSUBFILES_IOZ
      TZFILE_IOZ => TPFILE%TFILES_IOZ(JI)%TFILE
      IF (.NOT.TZFILE_IOZ%LOPENED) &
        CALL PRINT_MSG(NVERB_ERROR,'IO','IO_FILE_CLOSE_ll','file '//TRIM(TZFILE_IOZ%CNAME)//' is not opened')
      IF (TZFILE_IOZ%NOPEN_CURRENT/=1) &
        CALL PRINT_MSG(NVERB_WARNING,'IO','IO_FILE_CLOSE_ll','file '//TRIM(TZFILE_IOZ%CNAME)//&
                       ' is currently opened 0 or several times (expected only 1)')
      TZFILE_IOZ%LOPENED       = .FALSE.
      TZFILE_IOZ%NOPEN_CURRENT = 0
      TZFILE_IOZ%NCLOSE        = TZFILE_IOZ%NCLOSE + 1
    END DO
END SELECT
!
TPFILE%LOPENED       = .FALSE.
TPFILE%NOPEN_CURRENT = 0
TPFILE%NCLOSE        = TPFILE%NCLOSE + 1
!
IF (PRESENT(KRESP)) KRESP=IRESP
!
END SUBROUTINE IO_FILE_CLOSE_ll

SUBROUTINE FMCLOS_ll(TPFILE,KRESP,OPARALLELIO,HPROGRAM_ORIG)
!
!!    MODIFICATIONS
!!    -------------
!
!!      J.Escobar   18/10/10   bug with PGI compiler on ADJUSTL
!-------------------------------------------------------------------------------
USE MODD_CONF,  ONLY : CPROGRAM
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODE_IO_ll, ONLY : CLOSE_ll,UPCASE
#if !defined(MNH_SGI)
USE MODI_SYSTEM_MNH
#endif
  use mode_io_file_lfi,  only: io_close_file_lfi
#if defined(MNH_IOCDF4)
  use mode_io_file_nc4,  only: io_close_file_nc4
  use mode_io_write_nc4, only: io_write_coordvar_nc4
#endif
TYPE(TFILEDATA),      INTENT(INOUT)         :: TPFILE ! File structure
INTEGER,              INTENT(OUT), OPTIONAL :: KRESP   ! return-code if problems araised
LOGICAL,              INTENT(IN),  OPTIONAL :: OPARALLELIO
CHARACTER(LEN=*),     INTENT(IN),  OPTIONAL :: HPROGRAM_ORIG !To emulate a file coming from this program

INTEGER                 :: IRESP,IROWF
CHARACTER(LEN=28)       :: YFILEM  ! name of the file
CHARACTER(LEN=8)        :: YRESP
CHARACTER(LEN=10)       :: YCPIO
CHARACTER(LEN=14)       :: YTRANS
CHARACTER(LEN=100)      :: YCOMMAND
INTEGER                 :: IERR, IFITYP
INTEGER, SAVE           :: ICPT=0
LOGICAL                 :: GPARALLELIO

YFILEM  = TPFILE%CNAME

CALL PRINT_MSG(NVERB_DEBUG,'IO','FMCLOS_ll','closing '//TRIM(YFILEM))

IF ( PRESENT(OPARALLELIO) ) THEN
  GPARALLELIO = OPARALLELIO
ELSE
  GPARALLELIO = .TRUE.  !par defaut on active les IO paralleles en Z si possible
ENDIF

IRESP  = 0
IROWF  = 0

IROWF=LEN_TRIM(YFILEM)

IF (IROWF.EQ.0) THEN
  IRESP=-59
  GOTO 1000
ENDIF

#if defined(MNH_IOCDF4)
!Write coordinates variables in NetCDF file
IF (TPFILE%CMODE == 'WRITE' .AND. (TPFILE%CFORMAT=='NETCDF4' .OR. TPFILE%CFORMAT=='LFICDF4')) THEN
  CALL IO_WRITE_COORDVAR_NC4(TPFILE,HPROGRAM_ORIG=HPROGRAM_ORIG)
END IF
#endif

IF (TPFILE%LMASTER) THEN
  if (tpfile%cformat == 'LFI'     .or. tpfile%cformat == 'LFICDF4') call io_close_file_lfi(tpfile,iresp)
#if defined(MNH_IOCDF4)
  if (tpfile%cformat == 'NETCDF4' .or. tpfile%cformat == 'LFICDF4') call io_close_file_nc4(tpfile,iresp)
#endif
  IF (IRESP == 0 .AND. CPROGRAM/='LFICDF') THEN
    !! Write in pipe
#if defined(MNH_LINUX) || defined(MNH_SP4)
    YTRANS='xtransfer.x'
#elif defined(MNH_SX5)
    YTRANS='nectransfer.x'
#else
    YTRANS='fujitransfer.x'
#endif
    IFITYP = TPFILE%NLFITYPE
    
    SELECT CASE (IFITYP)
    CASE(:-1)
      IRESP=-66
      GOTO 500
    CASE(0)
      YCPIO='NIL'
    CASE(1)
      YCPIO='MESONH'
    CASE(2)
      PRINT *,'FILE ',YFILEM,' NOT TRANSFERED'
      GOTO 500
    CASE(3:)
      IRESP=-66
      GOTO 500
    END SELECT
!   WRITE (YCOMMAND,*) YTRANS,' ',YCPIO,' ',YFILEM
#if defined(MNH_LINUX) || defined(MNH_VPP) || defined(MNH_SX5) ||  defined(MNH_SP4)
    ICPT=ICPT+1
    WRITE (YCOMMAND,'(A," ",A," ",A," >> OUTPUT_TRANSFER",I3.3,"  2>&1 &")') TRIM(YTRANS),TRIM(YCPIO),TRIM(YFILEM),ICPT
!JUAN jusqu'a MASDEV4_4    WRITE (YCOMMAND,'(A," ",A," ",A,"  ")') TRIM(YTRANS),TRIM(YCPIO),TRIM(YFILEM)
#endif
#if defined(MNH_SGI)
    WRITE (YCOMMAND,'(A," ",A," ",A," &")') TRIM(YTRANS),TRIM(YCPIO),TRIM(YFILEM)
#endif

    PRINT *,'YCOMMAND =',YCOMMAND
#if !defined(MNH_SGI)
    CALL SYSTEM_MNH(YCOMMAND)
#endif
  END IF
END IF

500 CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
IF (IRESP /= 0) GOTO 1000

CALL CLOSE_ll(TPFILE,IOSTAT=IRESP,OPARALLELIO=GPARALLELIO)

1000 CONTINUE

IF (IRESP.NE.0)  THEN
  WRITE(YRESP,"( I0 )") IRESP
  CALL PRINT_MSG(NVERB_ERROR,'IO','FMCLOS_ll',TRIM(YFILEM)//': exit with IRESP='//TRIM(YRESP))
END IF

IF (PRESENT(KRESP)) KRESP=IRESP

! format: 14c for fujitransfer.x and mesonh/nil
!         32c for file name
! if you have to change this format one day, don't forget the blank after 1H
! 20 FORMAT(A14,1H ,A10,1H ,A32,1H ,A1)
!
END SUBROUTINE FMCLOS_ll

END MODULE MODE_FM
