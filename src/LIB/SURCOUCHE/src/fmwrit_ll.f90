!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  Philippe Wautelet: 10/01/2019: do not write scalars in Z-split files
!  Philippe Wautelet: 10/01/2019: write header also for Z-split files
!-----------------------------------------------------------------

#ifdef MNH_MPI_DOUBLE_PRECISION
#define MPI_FLOAT MPI_DOUBLE_PRECISION
#else
#define MPI_FLOAT MPI_REAL
#endif

#define MNH_SCALARS_IN_SPLITFILES 0

MODULE MODE_FMWRIT

  USE MODD_MPIF
  USE MODD_IO_ll, ONLY: TFILEDATA

  USE MODE_FIELD
  USE MODE_IO_WRITE_LFI
#if defined(MNH_IOCDF4)
  USE MODE_IO_WRITE_NC4
#endif

  IMPLICIT NONE 

  PRIVATE

  INTERFACE IO_WRITE_FIELD
     MODULE PROCEDURE IO_WRITE_FIELD_BYNAME_X0, IO_WRITE_FIELD_BYNAME_X1,  &
                      IO_WRITE_FIELD_BYNAME_X2, IO_WRITE_FIELD_BYNAME_X3,  &
                      IO_WRITE_FIELD_BYNAME_X4, IO_WRITE_FIELD_BYNAME_X5,  &
                      IO_WRITE_FIELD_BYNAME_X6,                            &
                      IO_WRITE_FIELD_BYNAME_N0, IO_WRITE_FIELD_BYNAME_N1,  &
                      IO_WRITE_FIELD_BYNAME_N2, IO_WRITE_FIELD_BYNAME_N3,  &
                      IO_WRITE_FIELD_BYNAME_L0, IO_WRITE_FIELD_BYNAME_L1,  &
                      IO_WRITE_FIELD_BYNAME_C0, IO_WRITE_FIELD_BYNAME_C1,  &
                      IO_WRITE_FIELD_BYNAME_T0,                            &
                      IO_WRITE_FIELD_BYFIELD_X0,IO_WRITE_FIELD_BYFIELD_X1, &
                      IO_WRITE_FIELD_BYFIELD_X2,IO_WRITE_FIELD_BYFIELD_X3, &
                      IO_WRITE_FIELD_BYFIELD_X4,IO_WRITE_FIELD_BYFIELD_X5, &
                      IO_WRITE_FIELD_BYFIELD_X6,                           &
                      IO_WRITE_FIELD_BYFIELD_N0,IO_WRITE_FIELD_BYFIELD_N1, &
                      IO_WRITE_FIELD_BYFIELD_N2,IO_WRITE_FIELD_BYFIELD_N3, &
                      IO_WRITE_FIELD_BYFIELD_L0,IO_WRITE_FIELD_BYFIELD_L1, &
                      IO_WRITE_FIELD_BYFIELD_C0,IO_WRITE_FIELD_BYFIELD_C1, &
                      IO_WRITE_FIELD_BYFIELD_T0
  END INTERFACE

  INTERFACE IO_WRITE_FIELD_BOX
     MODULE PROCEDURE IO_WRITE_FIELD_BOX_BYFIELD_X5
  END INTERFACE

  INTERFACE IO_WRITE_FIELD_LB
     MODULE PROCEDURE IO_WRITE_FIELD_BYNAME_LB, IO_WRITE_FIELD_BYFIELD_LB
  END INTERFACE

  PUBLIC IO_WRITE_FIELD, IO_WRITE_FIELD_BOX, IO_WRITE_FIELD_LB
  PUBLIC IO_WRITE_HEADER

CONTAINS 

  SUBROUTINE FIELD_METADATA_CHECK(TPFIELD,KTYPE,KDIMS,HCALLER)
    TYPE(TFIELDDATA), INTENT(IN) :: TPFIELD ! Field to check
    INTEGER,          INTENT(IN) :: KTYPE   ! Expected datatype
    INTEGER,          INTENT(IN) :: KDIMS   ! Expected number of dimensions
    CHARACTER(LEN=*), INTENT(IN) :: HCALLER ! name of the calling subroutine
    !
    CHARACTER(LEN=2) :: YDIMOK,YDIMKO
    CHARACTER(LEN=8) :: YTYPEOK,YTYPEKO
    !
    IF (TPFIELD%NGRID<0 .OR. TPFIELD%NGRID>8) THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO',HCALLER,'TPFIELD%NGRID is invalid for '//TRIM(TPFIELD%CMNHNAME))
    END IF
    IF (TPFIELD%NTYPE/=KTYPE) THEN
      CALL TYPE_WRITE(KTYPE,YTYPEOK)
      CALL TYPE_WRITE(TPFIELD%NTYPE,YTYPEKO)
      CALL PRINT_MSG(NVERB_WARNING,'IO',HCALLER,&
                     'TPFIELD%NTYPE should be '//YTYPEOK//' instead of '//YTYPEKO//' for '//TRIM(TPFIELD%CMNHNAME))
    END IF
    IF (TPFIELD%NDIMS/=KDIMS) THEN
      WRITE (YDIMOK,'(I2)') KDIMS
      WRITE (YDIMKO,'(I2)') TPFIELD%NDIMS
      CALL PRINT_MSG(NVERB_WARNING,'IO',HCALLER,&
                     'TPFIELD%NDIMS should be '//YDIMOK//' instead of '//YDIMKO//' for '//TRIM(TPFIELD%CMNHNAME))
    END IF
    !
    CONTAINS
    SUBROUTINE TYPE_WRITE(KTYPEINT,HTYPE)
      INTEGER,         INTENT(IN)  :: KTYPEINT
      CHARACTER(LEN=8),INTENT(OUT) :: HTYPE
      !
      SELECT CASE(KTYPEINT)
        CASE(TYPEINT)
          HTYPE = 'TYPEINT'
        CASE(TYPELOG)
          HTYPE = 'TYPELOG'
        CASE(TYPEREAL)
          HTYPE = 'TYPEREAL'
        CASE(TYPECHAR)
          HTYPE = 'TYPECHAR'
        CASE(TYPEDATE)
          HTYPE = 'TYPEDATE'
        CASE DEFAULT
          HTYPE = 'UNKNOWN'
      END SELECT
      !
    END SUBROUTINE TYPE_WRITE
  END SUBROUTINE FIELD_METADATA_CHECK


  SUBROUTINE IO_FILE_WRITE_CHECK(TPFILE,HSUBR,KRESP)
    TYPE(TFILEDATA),  INTENT(IN)  :: TPFILE
    CHARACTER(LEN=*), INTENT(IN)  :: HSUBR
    INTEGER,          INTENT(OUT) :: KRESP
    !
    KRESP = 0
    !
    !Check if file is opened
    IF (.NOT.TPFILE%LOPENED) THEN
      CALL PRINT_MSG(NVERB_ERROR,'IO',HSUBR,TRIM(TPFILE%CNAME)//' is not opened')
      KRESP = -201
      RETURN
    END IF
    !
    !Check if file is in the right opening mode
    IF (TPFILE%CMODE/='WRITE') THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO',HSUBR,&
                    TRIM(TPFILE%CNAME)//': writing in a file opened in '//TRIM(TPFILE%CMODE)//' mode')
    END IF
    !
    !Check fileformat
    IF (TPFILE%CFORMAT/='NETCDF4' .AND. TPFILE%CFORMAT/='LFI' .AND. TPFILE%CFORMAT/='LFICDF4') THEN
      CALL PRINT_MSG(NVERB_FATAL,'IO',HSUBR,&
                    TRIM(TPFILE%CNAME)//': invalid fileformat ('//TRIM(TPFILE%CFORMAT)//')')
      KRESP = -202
      RETURN
    END IF
    !
  END SUBROUTINE IO_FILE_WRITE_CHECK


  SUBROUTINE IO_WRITE_SELECT_FORMAT(TPFILE,OLFI,ONC4)
    TYPE(TFILEDATA), INTENT(IN)  :: TPFILE ! File structure
    LOGICAL,         INTENT(OUT) :: OLFI   ! Write in LFI format?
    LOGICAL,         INTENT(OUT) :: ONC4   ! Write in netCDF format?

    OLFI = .FALSE.
    ONC4 = .FALSE.
    IF (TPFILE%CFORMAT=='LFI'     .OR. TPFILE%CFORMAT=='LFICDF4') OLFI = .TRUE.
    IF (TPFILE%CFORMAT=='NETCDF4' .OR. TPFILE%CFORMAT=='LFICDF4') ONC4 = .TRUE.
  END SUBROUTINE IO_WRITE_SELECT_FORMAT


  SUBROUTINE IO_WRITE_HEADER(TPFILE,HDAD_NAME)
    TYPE(TFILEDATA),          INTENT(IN)  :: TPFILE   ! File structure
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: HDAD_NAME

    integer :: ifile

    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_HEADER_FILE','called for file '//TRIM(TPFILE%CNAME))

    CALL IO_WRITE_HEADER_ONEFILE(TPFILE,HDAD_NAME)

    !Write header also for the Z-split files
    DO IFILE=1,TPFILE%NSUBFILES_IOZ
      CALL IO_WRITE_HEADER_ONEFILE(TPFILE%TFILES_IOZ(IFILE)%TFILE,HDAD_NAME)
    END DO
  END SUBROUTINE IO_WRITE_HEADER


  SUBROUTINE IO_WRITE_HEADER_ONEFILE(TPFILE,HDAD_NAME)
    !
    USE MODD_CONF
    USE MODD_CONF_n,     ONLY: CSTORAGE_TYPE
    USE MODD_PARAMETERS, ONLY: NFILENAMELGTMAXLFI
    !
    TYPE(TFILEDATA),          INTENT(IN)  :: TPFILE   ! File structure
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: HDAD_NAME
    !
    CHARACTER(LEN=:),ALLOCATABLE :: YDAD_NAME
    INTEGER                      :: ILEN,ILEN2
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_HEADER_ONEFILE','called for file '//TRIM(TPFILE%CNAME))
    !
    IF ( ASSOCIATED(TPFILE%TDADFILE) .AND. PRESENT(HDAD_NAME) ) THEN
      IF ( TRIM(TPFILE%TDADFILE%CNAME) /= TRIM(HDAD_NAME) ) THEN
        CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_HEADER_ONEFILE','TPFILE%TDADFILE%CNAME /= HDAD_NAME')
      END IF
    END IF
    !
    CALL IO_WRITE_HEADER_NC4(TPFILE)
    !
    CALL IO_WRITE_FIELD(TPFILE,'MNHVERSION',  NMNHVERSION)
    CALL IO_WRITE_FIELD(TPFILE,'MASDEV',      NMASDEV)
    CALL IO_WRITE_FIELD(TPFILE,'BUGFIX',      NBUGFIX)
    CALL IO_WRITE_FIELD(TPFILE,'BIBUSER',     CBIBUSER)
    CALL IO_WRITE_FIELD(TPFILE,'PROGRAM',     CPROGRAM)
    CALL IO_WRITE_FIELD(TPFILE,'STORAGE_TYPE',CSTORAGE_TYPE)
    CALL IO_WRITE_FIELD(TPFILE,'MY_NAME',     TPFILE%CNAME)
    !
    IF ( ASSOCIATED(TPFILE%TDADFILE) ) THEN
      ILEN  = LEN_TRIM(TPFILE%TDADFILE%CNAME)
      ILEN2 = MAX(NFILENAMELGTMAXLFI,ILEN)
      ALLOCATE(CHARACTER(LEN=ILEN2) :: YDAD_NAME)
      IF(ILEN>0) THEN
        YDAD_NAME(1:ILEN) = TPFILE%TDADFILE%CNAME(1:ILEN)
        YDAD_NAME(ILEN+1:ILEN2) = ' '
      ELSE
        YDAD_NAME(:) = ' '
      END IF
    ELSE IF (PRESENT(HDAD_NAME)) THEN
      ILEN  = LEN_TRIM(HDAD_NAME)
      ILEN2 = MAX(NFILENAMELGTMAXLFI,ILEN)
      ALLOCATE(CHARACTER(LEN=ILEN2) :: YDAD_NAME)
      IF(ILEN>0) THEN
        YDAD_NAME(1:ILEN) = HDAD_NAME(1:ILEN)
        YDAD_NAME(ILEN+1:ILEN2) = ' '
      ELSE
        YDAD_NAME(:) = ' '
      END IF
    ELSE
      CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_HEADER_ONEFILE',TRIM(TPFILE%CNAME)// &
                     ': TPFILE%TDADFILE not associated and HDAD_NAME not provided')
      ALLOCATE(CHARACTER(LEN=NFILENAMELGTMAXLFI) :: YDAD_NAME)
      YDAD_NAME(:) = ' '
    ENDIF
    CALL IO_WRITE_FIELD(TPFILE,'DAD_NAME',YDAD_NAME)
    DEALLOCATE(YDAD_NAME)
    !
  END SUBROUTINE IO_WRITE_HEADER_ONEFILE


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X0(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),           INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),          INTENT(IN) :: HNAME    ! name of the field to write
    REAL,                      INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,          INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X0',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X0


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X0(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll, ONLY: GSMONOPROC,ISP
    !
    USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    REAL,TARGET,                 INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: IRESP
    !
    INTEGER                                  :: IK_FILE
    TYPE(TFILEDATA),POINTER                  :: TZFILE
    LOGICAL                                  :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    TZFILE => NULL()
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X0',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,0,'IO_WRITE_FIELD_BYFIELD_X0')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X0',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
       ELSE ! multiprocesses execution
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocesses execution
#if MNH_SCALARS_IN_SPLITFILES
       IF (TPFILE%NSUBFILES_IOZ>0) THEN
          ! write the data in all Z files
          DO IK_FILE=1,TPFILE%NSUBFILES_IOZ
             TZFILE => TPFILE%TFILES_IOZ(IK_FILE)%TFILE
             IF ( ISP == TZFILE%NMASTER_RANK )  THEN
                IF (GLFI) CALL IO_WRITE_FIELD_LFI(TZFILE,TPFIELD,PFIELD,IRESP)
                IF (GNC4) CALL IO_WRITE_FIELD_NC4(TZFILE,TPFIELD,PFIELD,IRESP)
             END IF
          END DO
       ENDIF
#endif
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X0',YMSG)
    END IF
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X0


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X1(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),           INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),          INTENT(IN) :: HNAME    ! name of the field to write
    REAL,DIMENSION(:),         INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,          INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return-code 
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X1',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X1


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X1(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll, ONLY: GSMONOPROC,ISP
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    REAL,DIMENSION(:),TARGET,    INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: IRESP
    INTEGER                                  :: ISIZEMAX
    REAL,DIMENSION(:),POINTER                :: ZFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    !
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    IRESP = 0
    GALLOC = .FALSE.
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X1',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,1,'IO_WRITE_FIELD_BYFIELD_X1')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X1',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
       ELSE ! multiprocesses execution
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X1','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(ZFIELDP(0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocesses execution
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X1',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(ZFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X1


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X2(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),           INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),          INTENT(IN) :: HNAME    ! name of the field to write
    REAL,DIMENSION(:,:),       INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,          INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return-code 
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X2',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X2


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X2(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_TIMEZ,         ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
#ifdef MNH_GA
    USE MODE_GA
#endif 
    USE MODE_GATHER_ll
    USE MODE_MNH_TIMING, ONLY : SECOND_MNH2
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    REAL,DIMENSION(:,:),TARGET,  INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    REAL,DIMENSION(:,:),POINTER              :: ZFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    !
    REAL*8,DIMENSION(2) :: T0,T1,T2
    REAL*8,DIMENSION(2) :: T11,T22
#ifdef MNH_GA
    REAL,DIMENSION(:,:),POINTER            :: ZFIELD_GA
#endif
    INTEGER                                :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE           :: YMSG
    CHARACTER(LEN=6)                       :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC = .FALSE.
    IHEXTOT = 2*JPHEXT+1
    !
    CALL SECOND_MNH2(T11)
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X2',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,2,'IO_WRITE_FIELD_BYFIELD_X2')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X2',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          !    IF (LPACK .AND. L1D .AND. YDIR=='XY') THEN 
          IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==IHEXTOT .AND. SIZE(PFIELD,2)==IHEXTOT) THEN 
             ZFIELDP=>PFIELD(JPHEXT+1:JPHEXT+1,JPHEXT+1:JPHEXT+1)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
             !    ELSE IF (LPACK .AND. L2D .AND. YDIR=='XY') THEN
          ELSEIF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==IHEXTOT) THEN
             ZFIELDP=>PFIELD(:,JPHEXT+1:JPHEXT+1)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
          END IF
       ELSE ! multiprocesses execution
          CALL SECOND_MNH2(T0)
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X2','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             ! I/O process case
             CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(ZFIELDP(0,0))
             GALLOC = .TRUE.
          END IF
          !   
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             IF (LPACK .AND. L2D) THEN
                CALL GATHER_XXFIELD('XX',PFIELD(:,JPHEXT+1),ZFIELDP(:,1),TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             ELSE
#ifdef MNH_GA
          !
          ! init/create the ga , dim3 = 1
          !
          CALL MNH_INIT_GA(SIZE(PFIELD,1),SIZE(PFIELD,2),1,YRECFM,"WRITE")
         !
         !   copy columun data to global arrays g_a 
         !
         ALLOCATE (ZFIELD_GA (SIZE(PFIELD,1),SIZE(PFIELD,2)))
         ZFIELD_GA = PFIELD
         call nga_put(g_a, lo_col, hi_col,ZFIELD_GA(NIXO_L,NIYO_L) , ld_col)  
         call ga_sync
         DEALLOCATE (ZFIELD_GA)
         IF (ISP == TPFILE%NMASTER_RANK) THEN
            !
            ! this proc get the  Z slide to write
            !
            lo_zplan(JPIZ) = 1
            hi_zplan(JPIZ) = 1
            call nga_get(g_a, lo_zplan, hi_zplan,ZFIELDP, ld_zplan)
         END IF
#else
         CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
#endif
             END IF
          END IF
          CALL SECOND_MNH2(T1)
          TIMEZ%T_WRIT2D_GATH=TIMEZ%T_WRIT2D_GATH + T1 - T0
          !
          IF (ISP == TPFILE%NMASTER_RANK) THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
#ifdef MNH_GA
         call ga_sync
#endif     
          CALL SECOND_MNH2(T2)
          TIMEZ%T_WRIT2D_WRIT=TIMEZ%T_WRIT2D_WRIT + T2 - T1
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X2',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(ZFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
    CALL SECOND_MNH2(T22)
    TIMEZ%T_WRIT2D_ALL=TIMEZ%T_WRIT2D_ALL + T22 - T11
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X2


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X3(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    REAL,DIMENSION(:,:,:),       INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X3',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X3


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X3(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll,            ONLY : GSMONOPROC,ISNPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll,    ONLY : JPHEXT
    USE MODD_TIMEZ,            ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_IO_TOOLS,         ONLY : IO_FILE
    USE MODE_IO_MANAGE_STRUCT, ONLY : IO_FILE_FIND_BYNAME
    USE MODE_MNH_TIMING,       ONLY : SECOND_MNH2
#ifdef MNH_GA
    USE MODE_GA
#endif
    USE MODD_VAR_ll, ONLY : MNH_STATUSES_IGNORE
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),TARGET,      INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    REAL,DIMENSION(:,:,:),TARGET,INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    REAL,DIMENSION(:,:,:),POINTER            :: ZFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    INTEGER                                  :: JK,JKK
    REAL,DIMENSION(:,:),POINTER              :: ZSLICE_ll,ZSLICE
    INTEGER                                  :: IK_FILE,IK_RANK,INB_PROC_REAL,JK_MAX
    INTEGER                                  :: JI,IXO,IXE,IYO,IYE
    REAL,DIMENSION(:,:),POINTER              :: TX2DP
    INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: STATUS
    LOGICAL                                  :: GALLOC_ll
    INTEGER,ALLOCATABLE,DIMENSION(:)         :: REQ_TAB
    INTEGER                                  :: NB_REQ
    TYPE TX_2DP
       REAL,DIMENSION(:,:), POINTER    :: X
    END TYPE TX_2DP
    TYPE(TX_2DP),ALLOCATABLE,DIMENSION(:) :: T_TX2DP
    REAL*8,DIMENSION(2) :: T0,T1,T2
    REAL*8,DIMENSION(2) :: T11,T22
#ifdef MNH_GA
    REAL,DIMENSION(:,:,:),POINTER          :: ZFIELD_GA
#endif
    INTEGER                                  :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    TYPE(TFILEDATA),POINTER                  :: TZFILE
    !
    TZFILE => NULL()
    !
    ZSLICE    => NULL()
    ZSLICE_ll => NULL()
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC    = .FALSE.
    GALLOC_ll = .FALSE.
    IHEXTOT = 2*JPHEXT+1
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X3',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL SECOND_MNH2(T11)
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,3,'IO_WRITE_FIELD_BYFIELD_X3')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X3',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC .AND. TPFILE%NSUBFILES_IOZ==0 ) THEN ! sequential execution
          !    IF (LPACK .AND. L1D .AND. YDIR=='XY') THEN 
          IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==IHEXTOT .AND. SIZE(PFIELD,2)==IHEXTOT) THEN 
             ZFIELDP=>PFIELD(JPHEXT+1:JPHEXT+1,JPHEXT+1:JPHEXT+1,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
             !    ELSE IF (LPACK .AND. L2D .AND. YDIR=='XY') THEN
          ELSEIF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==IHEXTOT) THEN
             ZFIELDP=>PFIELD(:,JPHEXT+1:JPHEXT+1,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
          END IF
       ELSEIF ( TPFILE%NSUBFILES_IOZ==0 .OR. YDIR=='--' ) THEN  ! multiprocesses execution & 1 proc IO
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X3','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          ! write 3D field in 1 time = output for graphique
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(ZFIELDP(0,0,0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             IF (LPACK .AND. L2D) THEN
                CALL GATHER_XXFIELD('XX',PFIELD(:,JPHEXT+1,:),ZFIELDP(:,1,:),TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             ELSE
                CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             END IF
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
          !
       ELSE ! multiprocesses execution & // IO
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X3','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF
          !
          !JUAN BG Z SLICE
          !
          !
#ifdef MNH_GA
          !
          ! init/create the ga
          !
          CALL SECOND_MNH2(T0)
          CALL MNH_INIT_GA(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3),YRECFM,"WRITE")
         !
         !   copy columun data to global arrays g_a 
         !
         ALLOCATE (ZFIELD_GA (SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3)))
         ZFIELD_GA = PFIELD
         call nga_put(g_a, lo_col, hi_col,ZFIELD_GA(NIXO_L,NIYO_L,1) , ld_col)  
         DEALLOCATE(ZFIELD_GA)
         call ga_sync
         CALL SECOND_MNH2(T1)
         TIMEZ%T_WRIT3D_SEND=TIMEZ%T_WRIT3D_SEND + T1 - T0
         !
         ! write the data
         !
         ALLOCATE(ZSLICE_ll(0,0)) ! to avoid bug on test of size
         GALLOC_ll = .TRUE.
         !
         DO JKK=1,IKU_ll
            !
            IK_FILE = IO_FILE(JKK,TPFILE%NSUBFILES_IOZ)
            TZFILE => TPFILE%TFILES_IOZ(IK_FILE+1)%TFILE
            !
            IK_RANK = TZFILE%NMASTER_RANK
            !
            IF (ISP == IK_RANK )  THEN 
               CALL SECOND_MNH2(T0)
               !
               IF ( SIZE(ZSLICE_ll) .EQ. 0 ) THEN
                  DEALLOCATE(ZSLICE_ll)
                  CALL ALLOCBUFFER_ll(ZSLICE_ll,ZSLICE,YDIR,GALLOC_ll)
               END IF
               !
               ! this proc get this JKK slide
               !
               lo_zplan(JPIZ) = JKK
               hi_zplan(JPIZ) = JKK
               call nga_get(g_a, lo_zplan, hi_zplan,ZSLICE_ll, ld_zplan)
               CALL SECOND_MNH2(T1)
               TIMEZ%T_WRIT3D_RECV=TIMEZ%T_WRIT3D_RECV + T1 - T0
               !
               IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZSLICE_ll,IRESP,KVERTLEVEL=JKK,KZFILE=IK_FILE+1)
               IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZSLICE_ll,IRESP,KVERTLEVEL=JKK,KZFILE=IK_FILE+1)
               CALL SECOND_MNH2(T2)
               TIMEZ%T_WRIT3D_WRIT=TIMEZ%T_WRIT3D_WRIT + T2 - T1
            END IF
         END DO
         !
         CALL SECOND_MNH2(T0) 
         call ga_sync
         CALL SECOND_MNH2(T1) 
         TIMEZ%T_WRIT3D_WAIT=TIMEZ%T_WRIT3D_WAIT + T1 - T0     
#else
          !
          ALLOCATE(ZSLICE_ll(0,0))
          GALLOC_ll = .TRUE.
          INB_PROC_REAL = MIN(TPFILE%NSUBFILES_IOZ,ISNPROC)
          Z_SLICE: DO JK=1,SIZE(PFIELD,3),INB_PROC_REAL
             !
             ! collect the data
             !
             JK_MAX=MIN(SIZE(PFIELD,3),JK+INB_PROC_REAL-1)
             !
             NB_REQ=0
             ALLOCATE(REQ_TAB(INB_PROC_REAL))
             ALLOCATE(T_TX2DP(INB_PROC_REAL))
             DO JKK=JK,JK_MAX
                !
                ! get the file & rank to write this level
                !
                IF (TPFILE%NSUBFILES_IOZ .GT. 1 ) THEN
                   IK_FILE = IO_FILE(JKK,TPFILE%NSUBFILES_IOZ)
                   TZFILE => TPFILE%TFILES_IOZ(IK_FILE+1)%TFILE
                ELSE
                   TZFILE => TPFILE
                END IF
                !
                IK_RANK = TZFILE%NMASTER_RANK
                !
                IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
                   STOP " XX NON PREVU SUR BG POUR LE MOMENT "
                   CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
                ELSEIF (YDIR == 'XY') THEN
                   IF (LPACK .AND. L2D) THEN
                      STOP " L2D NON PREVU SUR BG POUR LE MOMENT "
                      CALL GATHER_XXFIELD('XX',PFIELD(:,JPHEXT+1,:),ZFIELDP(:,1,:),TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
                   ELSE
                      CALL SECOND_MNH2(T0)
                      IF ( ISP /= IK_RANK )  THEN
                         ! Other processes
                         CALL GET_DOMWRITE_ll(ISP,'local',IXO,IXE,IYO,IYE)
                         IF (IXO /= 0) THEN ! intersection is not empty
                            NB_REQ = NB_REQ + 1
                            ALLOCATE(T_TX2DP(NB_REQ)%X(IXO:IXE,IYO:IYE))
                            ZSLICE => PFIELD(:,:,JKK)
                            TX2DP=>ZSLICE(IXO:IXE,IYO:IYE)
                            T_TX2DP(NB_REQ)%X=ZSLICE(IXO:IXE,IYO:IYE)
                            CALL MPI_ISEND(T_TX2DP(NB_REQ)%X,SIZE(TX2DP),MPI_FLOAT,IK_RANK-1,99+IK_RANK &
                                          & ,TZFILE%NMPICOMM,REQ_TAB(NB_REQ),IERR)
                            !CALL MPI_BSEND(TX2DP,SIZE(TX2DP),MPI_FLOAT,IK_RANK-1,99+IK_RANK,TZFILE%NMPICOMM,IERR)
                         END IF
                      END IF
                      CALL SECOND_MNH2(T1)
                      TIMEZ%T_WRIT3D_SEND=TIMEZ%T_WRIT3D_SEND + T1 - T0
                   END IF
                END IF
             END DO
             !
             ! write the data
             !
             DO JKK=JK,JK_MAX
                IF (TPFILE%NSUBFILES_IOZ .GT. 1 ) THEN
                   IK_FILE = IO_FILE(JKK,TPFILE%NSUBFILES_IOZ)
                   TZFILE => TPFILE%TFILES_IOZ(IK_FILE+1)%TFILE
                ELSE
                   TZFILE => TPFILE
                ENDIF
                IK_RANK = TZFILE%NMASTER_RANK
                !
                IF (ISP == IK_RANK )  THEN
                   CALL SECOND_MNH2(T0)
                   ! I/O proc case
                   IF ( SIZE(ZSLICE_ll) .EQ. 0 ) THEN
                      DEALLOCATE(ZSLICE_ll)
                      CALL ALLOCBUFFER_ll(ZSLICE_ll,ZSLICE,YDIR,GALLOC_ll)
                   END IF
                   DO JI=1,ISNPROC
                      CALL GET_DOMWRITE_ll(JI,'global',IXO,IXE,IYO,IYE)
                      IF (IXO /= 0) THEN ! intersection is not empty
                         TX2DP=>ZSLICE_ll(IXO:IXE,IYO:IYE)
                         IF (ISP == JI) THEN 
                            CALL GET_DOMWRITE_ll(JI,'local',IXO,IXE,IYO,IYE)
                            ZSLICE => PFIELD(:,:,JKK)
                            TX2DP = ZSLICE(IXO:IXE,IYO:IYE)
                         ELSE 
                            CALL MPI_RECV(TX2DP,SIZE(TX2DP),MPI_FLOAT,JI-1,99+IK_RANK,TZFILE%NMPICOMM,STATUS,IERR)
                         END IF
                      END IF
                   END DO
                   CALL SECOND_MNH2(T1)
                   TIMEZ%T_WRIT3D_RECV=TIMEZ%T_WRIT3D_RECV + T1 - T0
                   IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZSLICE_ll,IRESP,KVERTLEVEL=JKK,KZFILE=IK_FILE+1)
                   IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZSLICE_ll,IRESP,KVERTLEVEL=JKK,KZFILE=IK_FILE+1)
                   CALL SECOND_MNH2(T2)
                   TIMEZ%T_WRIT3D_WRIT=TIMEZ%T_WRIT3D_WRIT + T2 - T1
                END IF
             END DO
             !
             CALL SECOND_MNH2(T0) 
             IF (NB_REQ .GT.0 ) THEN
                CALL MPI_WAITALL(NB_REQ,REQ_TAB,MNH_STATUSES_IGNORE,IERR)
                DO JI=1,NB_REQ ;  DEALLOCATE(T_TX2DP(JI)%X) ; ENDDO
             END IF
             DEALLOCATE(T_TX2DP)
             DEALLOCATE(REQ_TAB)
             CALL SECOND_MNH2(T1) 
             TIMEZ%T_WRIT3D_WAIT=TIMEZ%T_WRIT3D_WAIT + T1 - T0
          END DO Z_SLICE
          !JUAN BG Z SLICE
! end of MNH_GA
#endif
       END IF ! multiprocesses execution
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X3',YMSG)
    END IF
    IF (GALLOC)    DEALLOCATE(ZFIELDP)
    IF (GALLOC_ll) DEALLOCATE(ZSLICE_ll)
    IF (PRESENT(KRESP)) KRESP = IRESP
    CALL SECOND_MNH2(T22)
    TIMEZ%T_WRIT3D_ALL=TIMEZ%T_WRIT3D_ALL + T22 - T11
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X3


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X4(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    REAL,DIMENSION(:,:,:,:),     INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X4',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X4


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X4(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_TIMEZ,         ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_IO_TOOLS,      ONLY : IO_FILE,IO_RANK
    USE MODE_MNH_TIMING,    ONLY : SECOND_MNH2
    USE MODD_VAR_ll,        ONLY : MNH_STATUSES_IGNORE
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),                 INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),                INTENT(IN) :: TPFIELD
    REAL,DIMENSION(:,:,:,:),TARGET,  INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,                INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    REAL,DIMENSION(:,:,:,:),POINTER          :: ZFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    INTEGER                                  :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC    = .FALSE.
    !
    IHEXTOT = 2*JPHEXT+1
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X4',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,4,'IO_WRITE_FIELD_BYFIELD_X4')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X4',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          !    IF (LPACK .AND. L1D .AND. YDIR=='XY') THEN 
          IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==IHEXTOT .AND. SIZE(PFIELD,2)==IHEXTOT) THEN 
             ZFIELDP=>PFIELD(JPHEXT+1:JPHEXT+1,JPHEXT+1:JPHEXT+1,:,:)
             !    ELSE IF (LPACK .AND. L2D .AND. YDIR=='XY') THEN
          ELSEIF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==IHEXTOT) THEN
             ZFIELDP=>PFIELD(:,JPHEXT+1:JPHEXT+1,:,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
          END IF
       ELSE
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X4','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(ZFIELDP(0,0,0,0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             IF (LPACK .AND. L2D) THEN
                CALL GATHER_XXFIELD('XX',PFIELD(:,JPHEXT+1,:,:),ZFIELDP(:,1,:,:),TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             ELSE
                CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             END IF
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocess execution
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X4',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(ZFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X4


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X5(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    REAL,DIMENSION(:,:,:,:,:),   INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X5',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X5


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X5(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_TIMEZ,         ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_IO_TOOLS,      ONLY : IO_FILE,IO_RANK
    USE MODE_MNH_TIMING,    ONLY : SECOND_MNH2
    USE MODD_VAR_ll,        ONLY : MNH_STATUSES_IGNORE
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),                 INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),                INTENT(IN) :: TPFIELD
    REAL,DIMENSION(:,:,:,:,:),TARGET,INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,                INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    REAL,DIMENSION(:,:,:,:,:),POINTER        :: ZFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    INTEGER                                  :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC    = .FALSE.
    !
    IHEXTOT = 2*JPHEXT+1
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X5',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,5,'IO_WRITE_FIELD_BYFIELD_X5')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X5',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          !    IF (LPACK .AND. L1D .AND. YDIR=='XY') THEN 
          IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==IHEXTOT .AND. SIZE(PFIELD,2)==IHEXTOT) THEN 
             ZFIELDP=>PFIELD(JPHEXT+1:JPHEXT+1,JPHEXT+1:JPHEXT+1,:,:,:)
             !    ELSE IF (LPACK .AND. L2D .AND. YDIR=='XY') THEN
          ELSEIF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==IHEXTOT) THEN
             ZFIELDP=>PFIELD(:,JPHEXT+1:JPHEXT+1,:,:,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
          END IF
       ELSE
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X5','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(ZFIELDP(0,0,0,0,0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             IF (LPACK .AND. L2D) THEN
                CALL GATHER_XXFIELD('XX',PFIELD(:,JPHEXT+1,:,:,:),ZFIELDP(:,1,:,:,:),&
                     & TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             ELSE
                CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             END IF
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocess execution
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X5',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(ZFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X5


  SUBROUTINE IO_WRITE_FIELD_BYNAME_X6(TPFILE,HNAME,PFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_X6',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),PFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_X6

  SUBROUTINE IO_WRITE_FIELD_BYFIELD_X6(TPFILE,TPFIELD,PFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC, ISP
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_TIMEZ,         ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_IO_TOOLS,      ONLY : IO_FILE,IO_RANK
    USE MODE_MNH_TIMING,    ONLY : SECOND_MNH2
    USE MODD_VAR_ll,        ONLY : MNH_STATUSES_IGNORE
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),                   INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),                  INTENT(IN) :: TPFIELD
    REAL,DIMENSION(:,:,:,:,:,:),TARGET,INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,OPTIONAL,                  INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    REAL,DIMENSION(:,:,:,:,:,:),POINTER      :: ZFIELDP
    LOGICAL                                  :: GLFI, GNC4
    LOGICAL                                  :: GALLOC
    INTEGER                                  :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC    = .FALSE.
    !
    IHEXTOT = 2*JPHEXT+1
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_X6',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEREAL,6,'IO_WRITE_FIELD_BYFIELD_X6')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_X6',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PFIELD,IRESP)
       ELSE
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(PFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_X6','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(ZFIELDP(0,0,0,0,0,0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocess execution
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_X6',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(ZFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_X6


  SUBROUTINE IO_WRITE_FIELD_BYNAME_N0(TPFILE,HNAME,KFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    INTEGER,                     INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_N0',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),KFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_N0


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_N0(TPFILE,TPFIELD,KFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC, ISP
    !*      0.    DECLARATIONS
    !             ------------
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    INTEGER,                     INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER                      :: IERR
    INTEGER                      :: IRESP
    INTEGER                      :: IK_FILE
    TYPE(TFILEDATA),POINTER      :: TZFILE
    LOGICAL                      :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE :: YMSG
    CHARACTER(LEN=6)             :: YRESP
    !
    IRESP = 0
    TZFILE => NULL()
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_N0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEINT,0,'IO_WRITE_FIELD_BYFIELD_N0')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_N0',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,KFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,KFIELD,IRESP)
       ELSE 
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,KFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,KFIELD,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocess execution
#if MNH_SCALARS_IN_SPLITFILES
       IF (TPFILE%NSUBFILES_IOZ>0) THEN
          ! write the data in all Z files
          DO IK_FILE=1,TPFILE%NSUBFILES_IOZ
             TZFILE => TPFILE%TFILES_IOZ(IK_FILE)%TFILE
             IF ( ISP == TZFILE%NMASTER_RANK )  THEN
                IF (GLFI) CALL IO_WRITE_FIELD_LFI(TZFILE,TPFIELD,KFIELD,IRESP)
                IF (GNC4) CALL IO_WRITE_FIELD_NC4(TZFILE,TPFIELD,KFIELD,IRESP)
             END IF
          END DO
       ENDIF
#endif
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(TPFIELD%CMNHNAME)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_N0',YMSG)
    END IF
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_N0


  SUBROUTINE IO_WRITE_FIELD_BYNAME_N1(TPFILE,HNAME,KFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    INTEGER,DIMENSION(:),        INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_N1',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),KFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_N1


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_N1(TPFILE,TPFIELD,KFIELD,KRESP)
    !
    USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),              INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),             INTENT(IN) :: TPFIELD
    INTEGER,DIMENSION(:),TARGET,  INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,             INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    INTEGER,DIMENSION(:),POINTER             :: IFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC = .FALSE.
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_N1',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEINT,1,'IO_WRITE_FIELD_BYFIELD_N1')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_N1',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,KFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,KFIELD,IRESP)
       ELSE ! multiprocesses execution
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(KFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(KFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_N1','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(IFIELDP(0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,KFIELD,IFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_N1',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(IFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_N1

  
  SUBROUTINE IO_WRITE_FIELD_BYNAME_N2(TPFILE,HNAME,KFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    INTEGER,DIMENSION(:,:),      INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_N2',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),KFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_N2


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_N2(TPFILE,TPFIELD,KFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_TIMEZ,         ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_MNH_TIMING,    ONLY : SECOND_MNH2
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),              INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),             INTENT(IN) :: TPFIELD
    INTEGER,DIMENSION(:,:),TARGET,INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,             INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    INTEGER,DIMENSION(:,:),POINTER           :: IFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    !
    REAL*8,DIMENSION(2) :: T0,T1,T2
    REAL*8,DIMENSION(2) :: T11,T22
    INTEGER                                  :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC = .FALSE.
    !
    IHEXTOT = 2*JPHEXT+1
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_N2',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL SECOND_MNH2(T11)
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEINT,2,'IO_WRITE_FIELD_BYFIELD_N2')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_N2',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (LPACK .AND. L1D .AND. SIZE(KFIELD,1)==IHEXTOT .AND. SIZE(KFIELD,2)==IHEXTOT) THEN 
             IFIELDP=>KFIELD(JPHEXT+1:JPHEXT+1,JPHEXT+1:JPHEXT+1)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
             !    ELSE IF (LPACK .AND. L2D .AND. YDIR=='XY') THEN
          ELSEIF (LPACK .AND. L2D .AND. SIZE(KFIELD,2)==IHEXTOT) THEN
             IFIELDP=>KFIELD(:,JPHEXT+1:JPHEXT+1)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,KFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,KFIELD,IRESP)
          END IF
       ELSE ! multiprocesses execution
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(KFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(KFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_N2','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          CALL SECOND_MNH2(T0)
          IF (ISP == TPFILE%NMASTER_RANK) THEN
             ! I/O process case
             CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(IFIELDP(0,0))
             GALLOC = .TRUE.
          END IF
          !   
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,KFIELD,IFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             IF (LPACK .AND. L2D) THEN
                CALL GATHER_XXFIELD('XX',KFIELD(:,JPHEXT+1),IFIELDP(:,1),TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             ELSE
                CALL GATHER_XYFIELD(KFIELD,IFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             END IF
          END IF
          CALL SECOND_MNH2(T1)
          TIMEZ%T_WRIT2D_GATH=TIMEZ%T_WRIT2D_GATH + T1 - T0
          !
          IF (ISP == TPFILE%NMASTER_RANK) THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
          END IF
          CALL SECOND_MNH2(T2)
          TIMEZ%T_WRIT2D_WRIT=TIMEZ%T_WRIT2D_WRIT + T2 - T1
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_N2',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(IFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
    CALL SECOND_MNH2(T22)
    TIMEZ%T_WRIT2D_ALL=TIMEZ%T_WRIT2D_ALL + T22 - T11
    !
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_N2


  SUBROUTINE IO_WRITE_FIELD_BYNAME_N3(TPFILE,HNAME,KFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    INTEGER,DIMENSION(:,:,:),    INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_N3',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),KFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_N3

  SUBROUTINE IO_WRITE_FIELD_BYFIELD_N3(TPFILE,TPFIELD,KFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_TIMEZ,         ONLY : TIMEZ
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    USE MODE_MNH_TIMING,    ONLY : SECOND_MNH2
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),                INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),               INTENT(IN) :: TPFIELD
    INTEGER,DIMENSION(:,:,:),TARGET,INTENT(IN) :: KFIELD   ! array containing the data field
    INTEGER,OPTIONAL,               INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    INTEGER,DIMENSION(:,:,:),POINTER         :: IFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    !
    REAL*8,DIMENSION(2) :: T11,T22
    INTEGER                                  :: IHEXTOT
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC = .FALSE.
    !
    IHEXTOT = 2*JPHEXT+1
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_N3',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL SECOND_MNH2(T11)
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEINT,3,'IO_WRITE_FIELD_BYFIELD_N3')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_N3',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (LPACK .AND. L1D .AND. SIZE(KFIELD,1)==IHEXTOT .AND. SIZE(KFIELD,2)==IHEXTOT) THEN 
             IFIELDP=>KFIELD(JPHEXT+1:JPHEXT+1,JPHEXT+1:JPHEXT+1,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
             !    ELSE IF (LPACK .AND. L2D .AND. YDIR=='XY') THEN
          ELSEIF (LPACK .AND. L2D .AND. SIZE(KFIELD,2)==IHEXTOT) THEN
             IFIELDP=>KFIELD(:,JPHEXT+1:JPHEXT+1,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,KFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,KFIELD,IRESP)
          END IF
       ELSE ! multiprocesses execution
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(KFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(KFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_N3','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK) THEN
             ! I/O process case
             CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(IFIELDP(0,0,0))
             GALLOC = .TRUE.
          END IF
          !   
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,KFIELD,IFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          ELSEIF (YDIR == 'XY') THEN
             IF (LPACK .AND. L2D) THEN
                CALL GATHER_XXFIELD('XX',KFIELD(:,JPHEXT+1,:),IFIELDP(:,1,:),TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             ELSE
                CALL GATHER_XYFIELD(KFIELD,IFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
             END IF
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK) THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,IFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_N3',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(IFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
    CALL SECOND_MNH2(T22)
    TIMEZ%T_WRIT3D_ALL=TIMEZ%T_WRIT3D_ALL + T22 - T11
    !
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_N3


  SUBROUTINE IO_WRITE_FIELD_BYNAME_L0(TPFILE,HNAME,OFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    LOGICAL,                     INTENT(IN) :: OFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_L0',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),OFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_L0

  SUBROUTINE IO_WRITE_FIELD_BYFIELD_L0(TPFILE,TPFIELD,OFIELD,KRESP)
    USE MODD_IO_ll,            ONLY : GSMONOPROC, ISP
    !
    USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
    !*      0.    DECLARATIONS
    !             ------------
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    LOGICAL,                     INTENT(IN) :: OFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER                      :: IERR
    INTEGER                      :: IRESP
    INTEGER                      :: IK_FILE
    LOGICAL                      :: GLFI, GNC4
    TYPE(TFILEDATA),POINTER      :: TZFILE
    CHARACTER(LEN=:),ALLOCATABLE :: YMSG
    CHARACTER(LEN=6)             :: YRESP
    !
    IRESP = 0
    TZFILE => NULL()
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_L0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPELOG,0,'IO_WRITE_FIELD_BYFIELD_L0')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_L0',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,OFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,OFIELD,IRESP)
       ELSE
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,OFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,OFIELD,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocesses execution
#if MNH_SCALARS_IN_SPLITFILES
       IF (TPFILE%NSUBFILES_IOZ>0) THEN
          ! write the data in all Z files
          DO IK_FILE=1,TPFILE%NSUBFILES_IOZ
             TZFILE => TPFILE%TFILES_IOZ(IK_FILE)%TFILE
             IF ( ISP == TZFILE%NMASTER_RANK )  THEN
                IF (GLFI) CALL IO_WRITE_FIELD_LFI(TZFILE,TPFIELD,OFIELD,IRESP)
                IF (GNC4) CALL IO_WRITE_FIELD_NC4(TZFILE,TPFIELD,OFIELD,IRESP)
             END IF
          END DO
       ENDIF
#endif
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(TPFIELD%CMNHNAME)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_L0',YMSG)
    END IF
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_L0


  SUBROUTINE IO_WRITE_FIELD_BYNAME_L1(TPFILE,HNAME,OFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    LOGICAL,DIMENSION(:),        INTENT(IN) :: OFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_L1',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),OFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_L1


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_L1(TPFILE,TPFIELD,OFIELD,KRESP)
    !
    USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
    !
    USE MODE_ALLOCBUFFER_ll
    USE MODE_GATHER_ll
    !
    IMPLICIT NONE
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),              INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),             INTENT(IN) :: TPFIELD
    LOGICAL,DIMENSION(:),TARGET,  INTENT(IN) :: OFIELD   ! array containing the data field
    INTEGER,OPTIONAL,             INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=2)                         :: YDIR     ! field form
    INTEGER                                  :: IERR
    INTEGER                                  :: ISIZEMAX
    INTEGER                                  :: IRESP
    LOGICAL,DIMENSION(:),POINTER             :: GFIELDP
    LOGICAL                                  :: GALLOC
    LOGICAL                                  :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YDIR     = TPFIELD%CDIR
    !
    IRESP = 0
    GALLOC = .FALSE.
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_L1',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPELOG,1,'IO_WRITE_FIELD_BYFIELD_L1')
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_L1',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,OFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,OFIELD,IRESP)
       ELSE ! multiprocesses execution
#if ( MNH_INT == 4 )
          CALL MPI_ALLREDUCE(SIZE(OFIELD),ISIZEMAX,1,MPI_INTEGER,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#else
          CALL MPI_ALLREDUCE(SIZE(OFIELD),ISIZEMAX,1,MPI_INTEGER8,MPI_MAX,TPFILE%NMPICOMM,IRESP)
#endif
          IF (ISIZEMAX==0) THEN
             CALL PRINT_MSG(NVERB_INFO,'IO','IO_WRITE_FIELD_BYFIELD_L1','ignoring variable with a zero size ('//TRIM(YRECFM)//')')
             IF (PRESENT(KRESP)) KRESP=0
             RETURN
          END IF

          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             CALL ALLOCBUFFER_ll(GFIELDP,OFIELD,YDIR,GALLOC)
          ELSE
             ALLOCATE(GFIELDP(0))
             GALLOC = .TRUE.
          END IF
          !
          IF (YDIR == 'XX' .OR. YDIR =='YY') THEN
             CALL GATHER_XXFIELD(YDIR,OFIELD,GFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM)
          END IF
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,GFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,GFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_L1',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(GFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_L1


  SUBROUTINE IO_WRITE_FIELD_BYNAME_C0(TPFILE,HNAME,HFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    CHARACTER(LEN=*),            INTENT(IN) :: HFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_C0',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),HFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_C0


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_C0(TPFILE,TPFIELD,HFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC, ISP
    !
    !*      0.    DECLARATIONS
    !             ------------
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    CHARACTER(LEN=*),            INTENT(IN) :: HFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER                      :: IERR
    INTEGER                      :: IRESP
    LOGICAL                      :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE :: YMSG
    CHARACTER(LEN=6)             :: YRESP
    !
    IRESP = 0
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_C0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPECHAR,0,'IO_WRITE_FIELD_BYFIELD_C0')
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (LEN(HFIELD)==0 .AND. GLFI) THEN
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_C0',&
                     'zero-size string not allowed if LFI output for '//TRIM(TPFIELD%CMNHNAME))
    END IF
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_C0',IRESP)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,HFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,HFIELD,IRESP)
       ELSE 
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,HFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,HFIELD,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(TPFIELD%CMNHNAME)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_C0',YMSG)
    END IF
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_C0


  SUBROUTINE IO_WRITE_FIELD_BYNAME_C1(TPFILE,HNAME,HFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),              INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),             INTENT(IN) :: HNAME    ! name of the field to write
    CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: HFIELD   ! array containing the data field
    INTEGER,OPTIONAL,             INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_C1',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),HFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_C1


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_C1(TPFILE,TPFIELD,HFIELD,KRESP)
    USE MODD_IO_ll,         ONLY : GSMONOPROC, ISP
    !
    !*      0.    DECLARATIONS
    !             ------------
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),              INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),             INTENT(IN) :: TPFIELD
    CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: HFIELD   ! array containing the data field
    INTEGER,OPTIONAL,             INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER                          :: IERR
    INTEGER                          :: IRESP
    INTEGER                          :: J,JJ
    INTEGER                          :: ILE, IP
    INTEGER,DIMENSION(:),ALLOCATABLE :: IFIELD
    INTEGER                          :: ILENG
    LOGICAL                          :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE     :: YMSG
    CHARACTER(LEN=6)                 :: YRESP
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_C1',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPECHAR,1,'IO_WRITE_FIELD_BYFIELD_C1')
    !
    IRESP = 0
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF(GLFI) THEN
      ILE=LEN(HFIELD)
      IP=SIZE(HFIELD)
      ILENG=ILE*IP
      !
      IF (ILENG==0) THEN
        IP=1
        ILE=1
        ILENG=1
        ALLOCATE(IFIELD(1))
        IFIELD(1)=IACHAR(' ')
      ELSE
        ALLOCATE(IFIELD(ILENG))
        DO JJ=1,IP
          DO J=1,ILE
            IFIELD(ILE*(JJ-1)+J)=IACHAR(HFIELD(JJ)(J:J))
          END DO
        END DO
      END IF
    END IF
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_C1',IRESP)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,HFIELD,IRESP)
       ELSE 
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,IFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,HFIELD,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(TPFIELD%CMNHNAME)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_C1',YMSG)
    END IF
    IF (ALLOCATED(IFIELD)) DEALLOCATE(IFIELD)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_C1


  SUBROUTINE IO_WRITE_FIELD_BYNAME_T0(TPFILE,HNAME,TFIELD,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    TYPE (DATE_TIME),            INTENT(IN) :: TFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_T0',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD(TPFILE,TFIELDLIST(ID),TFIELD,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_T0


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_T0(TPFILE,TPFIELD,TFIELD,KRESP)
    USE MODD_IO_ll, ONLY : GSMONOPROC, ISP
    USE MODD_TYPE_DATE
    !
    !*      0.    DECLARATIONS
    !             ------------
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),            INTENT(IN) :: TPFIELD
    TYPE (DATE_TIME),            INTENT(IN) :: TFIELD   ! array containing the data field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER                      :: IERR
    INTEGER                      :: IRESP
    LOGICAL                      :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE :: YMSG
    CHARACTER(LEN=6)             :: YRESP
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_T0',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
    !
    CALL FIELD_METADATA_CHECK(TPFIELD,TYPEDATE,0,'IO_WRITE_FIELD_BYFIELD_T0')
    !
    IRESP = 0
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_T0',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,TFIELD,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,TFIELD,IRESP)
       ELSE 
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,TFIELD,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,TFIELD,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(TPFIELD%CMNHNAME)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_T0',YMSG)
    END IF
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_T0


  SUBROUTINE IO_WRITE_FIELD_BYNAME_LB(TPFILE,HNAME,KL3D,PLB,KRESP)
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN) :: TPFILE
    CHARACTER(LEN=*),            INTENT(IN) :: HNAME    ! name of the field to write
    INTEGER,                     INTENT(IN) :: KL3D     ! size of the LB array in FM
    REAL,DIMENSION(:,:,:),       INTENT(IN) :: PLB      ! array containing the LB field
    INTEGER,OPTIONAL,            INTENT(OUT):: KRESP    ! return-code
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER :: ID ! Index of the field
    INTEGER :: IRESP ! return_code
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYNAME_LB',TRIM(TPFILE%CNAME)//': writing '//TRIM(HNAME))
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME(HNAME,ID,IRESP)
    !
    IF(IRESP==0) CALL IO_WRITE_FIELD_LB(TPFILE,TFIELDLIST(ID),KL3D,PLB,IRESP)
    !
    IF (PRESENT(KRESP)) KRESP = IRESP
    !
  END SUBROUTINE IO_WRITE_FIELD_BYNAME_LB


  SUBROUTINE IO_WRITE_FIELD_BYFIELD_LB(TPFILE,TPFIELD,KL3D,PLB,KRESP)
    !
    USE MODD_IO_ll,         ONLY : GSMONOPROC,ISNPROC,ISP,L1D,L2D,LPACK
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_VAR_ll,        ONLY : MNH_STATUSES_IGNORE
    !
    USE MODE_DISTRIB_LB
    USE MODE_TOOLS_ll,      ONLY : GET_GLOBALDIMS_ll
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),             INTENT(IN)    :: TPFILE
    TYPE(TFIELDDATA),            INTENT(INOUT) :: TPFIELD
    INTEGER,                     INTENT(IN)    :: KL3D   ! size of the LB array in FM
    REAL,DIMENSION(:,:,:),TARGET,INTENT(IN)    :: PLB    ! array containing the LB field
    INTEGER,OPTIONAL,            INTENT(OUT)   :: KRESP  ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    CHARACTER(LEN=28)                        :: YFILEM   ! FM-file name
    CHARACTER(LEN=NMNHNAMELGTMAX)            :: YRECFM   ! name of the article to write
    CHARACTER(LEN=4)                         :: YLBTYPE  ! 'LBX','LBXU','LBY' or 'LBYV'
    INTEGER                                  :: IRIM     ! size of the LB area
    INTEGER                                  :: IERR
    INTEGER                                  :: IRESP
    REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: Z3D
    REAL,DIMENSION(:,:,:), POINTER           :: TX3DP
    INTEGER                                  :: IIMAX_ll,IJMAX_ll
    INTEGER                                  :: JI
    INTEGER                                  :: IIB,IIE,IJB,IJE
    INTEGER, DIMENSION(MPI_STATUS_SIZE)      :: STATUS
    INTEGER,ALLOCATABLE,DIMENSION(:)         :: REQ_TAB
    INTEGER                                  :: NB_REQ,IKU
    LOGICAL                                  :: GLFI, GNC4
    TYPE TX_3DP
       REAL,DIMENSION(:,:,:), POINTER    :: X
    END TYPE TX_3DP
    TYPE(TX_3DP),ALLOCATABLE,DIMENSION(:)    :: T_TX3DP
    CHARACTER(LEN=:),ALLOCATABLE             :: YMSG
    CHARACTER(LEN=6)                         :: YRESP
    !
    YFILEM   = TPFILE%CNAME
    YRECFM   = TPFIELD%CMNHNAME
    YLBTYPE  = TPFIELD%CLBTYPE
    !
    IRESP = 0
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BYFIELD_LB',TRIM(YFILEM)//': writing '//TRIM(YRECFM))
    !
    IF (YLBTYPE/='LBX' .AND. YLBTYPE/='LBXU' .AND. YLBTYPE/='LBY' .AND. YLBTYPE/='LBYV') THEN
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_LB','unknown LBTYPE ('//YLBTYPE//')')
      RETURN
    END IF
    !
    IF (TPFIELD%CDIR/='') THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO','IO_WRITE_FIELD_BYFIELD_LB','CDIR was set for '//TRIM(YRECFM))
      TPFIELD%CDIR=''
    END IF
    !
    IRIM = (KL3D-2*JPHEXT)/2
    IF (KL3D /= 2*(IRIM+JPHEXT)) THEN
       IRESP = -30
       GOTO 1000
    END IF
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BYFIELD_LB',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN  ! sequential execution
          IF (LPACK .AND. L2D) THEN
             TX3DP=>PLB(:,JPHEXT+1:JPHEXT+1,:)
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,TX3DP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,TX3DP,IRESP)
          ELSE
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,PLB,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,PLB,IRESP)
          END IF
       ELSE
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             ! I/O proc case
             CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
             IF (YLBTYPE == 'LBX' .OR. YLBTYPE == 'LBXU') THEN 
                ALLOCATE(Z3D((IRIM+JPHEXT)*2,IJMAX_ll+2*JPHEXT,SIZE(PLB,3)))
             ELSE ! YLBTYPE == 'LBY' .OR. YLBTYPE == 'LBYV' 
                ALLOCATE(Z3D(IIMAX_ll+2*JPHEXT,(IRIM+JPHEXT)*2,SIZE(PLB,3)))
             END IF
             DO JI = 1,ISNPROC
                CALL GET_DISTRIB_LB(YLBTYPE,JI,'FM','WRITE',IRIM,IIB,IIE,IJB,IJE)
                IF (IIB /= 0) THEN
                   TX3DP=>Z3D(IIB:IIE,IJB:IJE,:)
                   IF (ISP /= JI) THEN
                      CALL MPI_RECV(TX3DP,SIZE(TX3DP),MPI_FLOAT,JI-1,99,TPFILE%NMPICOMM,STATUS,IERR)
                   ELSE
                      CALL GET_DISTRIB_LB(YLBTYPE,JI,'LOC','WRITE',IRIM,IIB,IIE,IJB,IJE)
                      TX3DP = PLB(IIB:IIE,IJB:IJE,:)
                   END IF
                END IF
             END DO
             IF (LPACK .AND. L2D) THEN
                TX3DP=>Z3D(:,JPHEXT+1:JPHEXT+1,:)
             ELSE
                TX3DP=>Z3D
             END IF
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,TX3DP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,TX3DP,IRESP)
          ELSE
             NB_REQ=0
             ALLOCATE(REQ_TAB(1))
             ALLOCATE(T_TX3DP(1))
             IKU = SIZE(PLB,3)
             ! Other processes
             CALL GET_DISTRIB_LB(YLBTYPE,ISP,'LOC','WRITE',IRIM,IIB,IIE,IJB,IJE)
             IF (IIB /= 0) THEN
                TX3DP=>PLB(IIB:IIE,IJB:IJE,:)
                NB_REQ = NB_REQ + 1
                ALLOCATE(T_TX3DP(NB_REQ)%X(IIB:IIE,IJB:IJE,IKU))  
                T_TX3DP(NB_REQ)%X=PLB(IIB:IIE,IJB:IJE,:)
                CALL MPI_ISEND(T_TX3DP(NB_REQ)%X,SIZE(TX3DP),MPI_FLOAT,TPFILE%NMASTER_RANK-1,99, &
                               TPFILE%NMPICOMM,REQ_TAB(NB_REQ),IERR)
                !CALL MPI_BSEND(TX3DP,SIZE(TX3DP),MPI_FLOAT,TPFILE%NMASTER_RANK-1,99,TPFILE%NMPICOMM,IERR)
             END IF
             IF (NB_REQ .GT.0 ) THEN
                CALL MPI_WAITALL(NB_REQ,REQ_TAB,MNH_STATUSES_IGNORE,IERR)
                DEALLOCATE(T_TX3DP(1)%X) 
             END IF
             DEALLOCATE(T_TX3DP,REQ_TAB)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF
    END IF
    !
1000 CONTINUE
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(YRECFM)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BYFIELD_LB',YMSG)
    END IF
    !
    IF (ALLOCATED(Z3D)) DEALLOCATE(Z3D)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BYFIELD_LB


  SUBROUTINE IO_WRITE_FIELD_BOX_BYFIELD_X5(TPFILE,TPFIELD,HBUDGET,PFIELD,KXOBOX,KXEBOX,KYOBOX,KYEBOX,KRESP)
    !
    USE MODD_IO_ll, ONLY : GSMONOPROC, ISP
    !
    USE MODE_GATHER_ll
    !
    !
    !*      0.1   Declarations of arguments
    !
    TYPE(TFILEDATA),                 INTENT(IN) :: TPFILE
    TYPE(TFIELDDATA),                INTENT(IN) :: TPFIELD
    CHARACTER(LEN=*),                INTENT(IN) :: HBUDGET  ! 'BUDGET' (budget)  or 'OTHER' (MesoNH field)
    REAL,DIMENSION(:,:,:,:,:),TARGET,INTENT(IN) :: PFIELD   ! array containing the data field
    INTEGER,                         INTENT(IN) :: KXOBOX   ! 
    INTEGER,                         INTENT(IN) :: KXEBOX   ! Global coordinates of the box
    INTEGER,                         INTENT(IN) :: KYOBOX   ! 
    INTEGER,                         INTENT(IN) :: KYEBOX   ! 
    INTEGER,OPTIONAL,                INTENT(OUT):: KRESP    ! return-code 
    !
    !*      0.2   Declarations of local variables
    !
    INTEGER                             :: IERR
    INTEGER                             :: IRESP
    REAL,DIMENSION(:,:,:,:,:),POINTER   :: ZFIELDP
    LOGICAL                             :: GALLOC
    LOGICAL                             :: GLFI, GNC4
    CHARACTER(LEN=:),ALLOCATABLE        :: YMSG
    CHARACTER(LEN=6)                    :: YRESP
    !
    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_WRITE_FIELD_BOX_BYFIELD_X5',TRIM(TPFILE%CNAME)//': writing '//TRIM(TPFIELD%CMNHNAME))
    !
    IRESP = 0
    GALLOC = .FALSE.
    !
    CALL IO_FILE_WRITE_CHECK(TPFILE,'IO_WRITE_FIELD_BOX_BYFIELD_X5',IRESP)
    !
    CALL IO_WRITE_SELECT_FORMAT(TPFILE,GLFI,GNC4)
    !
    IF (IRESP==0) THEN
       IF (GSMONOPROC) THEN ! sequential execution
          IF (HBUDGET /= 'BUDGET') THEN
             ! take the sub-section of PFIELD defined by the box
             ZFIELDP=>PFIELD(KXOBOX:KXEBOX,KYOBOX:KYEBOX,:,:,:)
          ELSE
             ! take the field as a budget
             ZFIELDP=>PFIELD
          END IF
          IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
          IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
       ELSE ! multiprocesses execution
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             ! Allocate the box
             ALLOCATE(ZFIELDP(KXEBOX-KXOBOX+1,KYEBOX-KYOBOX+1,SIZE(PFIELD,3),&
                  & SIZE(PFIELD,4),SIZE(PFIELD,5)))
             GALLOC = .TRUE.
          ELSE
             ALLOCATE(ZFIELDP(0,0,0,0,0))
             GALLOC = .TRUE.
          END IF
          !
          CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TPFILE%NMASTER_RANK,TPFILE%NMPICOMM,&
               & KXOBOX,KXEBOX,KYOBOX,KYEBOX,HBUDGET)
          !
          IF (ISP == TPFILE%NMASTER_RANK)  THEN
             IF (GLFI) CALL IO_WRITE_FIELD_LFI(TPFILE,TPFIELD,ZFIELDP,IRESP)
             IF (GNC4) CALL IO_WRITE_FIELD_NC4(TPFILE,TPFIELD,ZFIELDP,IRESP)
          END IF
          !
          CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TPFILE%NMASTER_RANK-1,TPFILE%NMPICOMM,IERR)
       END IF ! multiprocesses execution
    END IF
    !
    IF (IRESP.NE.0) THEN
      WRITE(YRESP, '( I6 )') IRESP
      YMSG = 'RESP='//YRESP//' when writing '//TRIM(TPFIELD%CMNHNAME)//' in '//TRIM(TPFILE%CNAME)
      CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELD_BOX_BYFIELD_X5',YMSG)
    END IF
    IF (GALLOC) DEALLOCATE(ZFIELDP)
    IF (PRESENT(KRESP)) KRESP = IRESP
  END SUBROUTINE IO_WRITE_FIELD_BOX_BYFIELD_X5

END MODULE MODE_FMWRIT
