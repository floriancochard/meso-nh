!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE mode_util
  USE MODD_IO_ll,  ONLY: TFILE_ELT
  USE MODD_NETCDF, ONLY: DIMCDF, IDCDF_KIND
  USE MODD_PARAMETERS, ONLY: NLFIMAXCOMMENTLENGTH, NMNHNAMELGTMAX

  USE MODE_FIELD
  USE MODE_FMREAD
  USE MODE_FMWRIT

  USE mode_options

  USE NETCDF

  IMPLICIT NONE 

  INTEGER,PARAMETER :: MAXRAW=10
  INTEGER,PARAMETER :: MAXFILES=100
  INTEGER,parameter :: MAXDATES=100

  INTEGER,PARAMETER :: FM_FIELD_SIZE = 32

  TYPE workfield
     CHARACTER(LEN=NMNHNAMELGTMAX)         :: name   ! nom du champ
     LOGICAL                               :: found  ! T if found in the input file
     LOGICAL                               :: calc   ! T if computed from other variables
     LOGICAL                               :: tbw    ! to be written or not
     LOGICAL                               :: tbr    ! to be read or not
     LOGICAL                               :: LSPLIT = .FALSE. ! TRUE if variable is split by vertical level
     INTEGER                               :: NSIZE = 0 ! Size of the variable (in number of elements)
     INTEGER                               :: NSRC = 0 ! Number of variables used to compute the variable (needed only if calc=.true.)
     INTEGER                               :: NDIMS_FILE  ! Number of dims (as present in input file)
     INTEGER,DIMENSION(:),ALLOCATABLE      :: NDIMSIZES_FILE  ! Dimensions sizes (as present in input file)
     CHARACTER(LEN=NF90_MAX_NAME),DIMENSION(:),ALLOCATABLE :: CDIMNAMES_FILE  ! Dimensions names (as present in input file)
     CHARACTER(LEN=40)                     :: CUNITS_FILE = '' ! Units (as present in input file)
     INTEGER                               :: NGRID_FILE  ! Grid  number (as present in input file)
     INTEGER                               :: NTYPE_FILE  ! netCDF datatype (NF90_CHAR, NF90_INT...) (as present in input file)
     INTEGER,DIMENSION(MAXRAW)             :: src    ! List of variables used to compute the variable (needed only if calc=.true.)
     INTEGER                               :: tgt    ! Target: id of the variable that use it (calc variable)
     TYPE(TFIELDDATA)                      :: TFIELD ! Metadata about the field
     TYPE(DIMCDF),DIMENSION(:),ALLOCATABLE :: TDIMS  ! Dimensions of the field
  END TYPE workfield

  LOGICAL(KIND=LFI_INT), PARAMETER :: ltrue  = .TRUE.
  LOGICAL(KIND=LFI_INT), PARAMETER :: lfalse = .FALSE.

  CHARACTER(LEN=6) :: CPROGRAM_ORIG

CONTAINS 
  SUBROUTINE parse_infiles(infiles, outfiles, KNFILES_OUT, nbvar_infile, nbvar_tbr, nbvar_calc, nbvar_tbw, &
                           tpreclist, options, runmode)
    USE MODD_DIM_n,      ONLY: NIMAX_ll, NJMAX_ll, NKMAX
    USE MODD_PARAMETERS, ONLY: JPHEXT, JPVEXT, NGRIDUNKNOWN

    use mode_io_tools_nc4, only: io_guess_dimids_nc4

    TYPE(TFILE_ELT),DIMENSION(:),         INTENT(IN)  :: infiles
    TYPE(TFILE_ELT),DIMENSION(:),         INTENT(IN)  :: outfiles
    INTEGER,                              INTENT(IN)  :: KNFILES_OUT
    INTEGER,                              INTENT(IN)  :: nbvar_infile, nbvar_tbr, nbvar_calc, nbvar_tbw
    TYPE(workfield), DIMENSION(:),POINTER,INTENT(OUT) :: tpreclist
    TYPE(option),DIMENSION(:),            INTENT(IN)  :: options
    INTEGER,                              INTENT(IN)  :: runmode

    TYPE TLFIDATE
      CHARACTER(LEN=FM_FIELD_SIZE) :: CNAME = ''     !Name of the date variable
      INTEGER                      :: NIDX_DATE = -1 !Index of the date part
      INTEGER                      :: NIDX_TIME = -1 !Index of the time part
    END TYPE TLFIDATE

    CHARACTER(LEN=FM_FIELD_SIZE)             :: yrecfm, YDATENAME
    CHARACTER(LEN=FM_FIELD_SIZE)             :: var_calc
    CHARACTER(LEN=FM_FIELD_SIZE),dimension(MAXRAW) :: var_raw
    CHARACTER(LEN=1)                         :: YNDIMS
    CHARACTER(LEN=32)                        :: YTYPE
    INTEGER                                  :: ji,jj
    INTEGER                                  :: ndb, nde, ndey, idx, idx_out, idx_var, maxvar
    INTEGER                                  :: leng
    INTEGER                                  :: IID, IRESP, IDATES, ICURDATE
    INTEGER                                  :: IDXDATE, IDXTIME
    INTEGER(KIND=LFI_INT)                    :: iresp2,ilu,ileng,ipos
    INTEGER(KIND=IDCDF_KIND)                 :: kcdf_id, kcdf_id2, var_id
    INTEGER(KIND=IDCDF_KIND)                 :: status
    LOGICAL                                  :: ladvan
    LOGICAL                                  :: GOK
    TYPE(TLFIDATE),DIMENSION(MAXDATES)       :: TLFIDATES

    CALL PRINT_MSG(NVERB_DEBUG,'IO','parse_infiles','called')

    IF (options(OPTSPLIT)%set) THEN
      idx_out = 0
    ELSE
      idx_out = 1
    END IF

    IDATES = 0

    IF (runmode==MODECDF2LFI) THEN
      !This file is a dummy one to manage netCDF dims
      idx_out = KNFILES_OUT
    END IF

    IF (INFILES(1)%TFILE%CFORMAT == 'LFI') THEN
      ilu = INFILES(1)%TFILE%NLFIFLU
    ELSE IF (INFILES(1)%TFILE%CFORMAT == 'NETCDF4') THEN
      kcdf_id = INFILES(1)%TFILE%NNCID
    END IF

    PRINT *,'MESONH 3D, 2D articles DIMENSIONS used :'
    PRINT *,'DIMX =',NIMAX_ll+2*JPHEXT
    PRINT *,'DIMY =',NJMAX_ll+2*JPHEXT
    PRINT *,'DIMZ =',NKMAX   +2*JPVEXT

    ! Phase 1 : build articles list to convert.
    !
    !    Pour l'instant tous les articles du fichier LFI sont
    !    convertis. On peut modifier cette phase pour prendre en
    !    compte un sous-ensemble d'article (liste definie par
    !    l'utilisateur par exemple)
    !
    IF (options(OPTVAR)%set) THEN
      ALLOCATE(tpreclist(nbvar_tbr+nbvar_calc))
      DO ji=1,nbvar_tbr+nbvar_calc
        tpreclist(ji)%found  = .FALSE.
        tpreclist(ji)%calc   = .FALSE. !By default variables are not computed from others
        tpreclist(ji)%tbw    = .TRUE.  !By default variables are written
        tpreclist(ji)%tbr    = .TRUE.  !By default variables are read
        tpreclist(ji)%src(:) = -1
        tpreclist(ji)%tgt    = -1
      END DO

       ! A variable list is provided with -v var1,...
       ndb  = 1
       idx_var = 1
       DO ji=1,nbvar_tbw
          !crash compiler GCC 4.2.0: nde = INDEX(TRIM(options(OPTVAR)%cvalue(ndb:)),',')
          nde = INDEX(TRIM(options(OPTVAR)%cvalue(ndb:len(trim(options(OPTVAR)%cvalue)))),',')
          IF (nde == 0) nde = LEN( TRIM(options(OPTVAR)%cvalue(ndb:len(trim(options(OPTVAR)%cvalue)))) ) + 1
          yrecfm = options(OPTVAR)%cvalue(ndb:ndb+nde-2)
          !Detect operations on variables (only + is supported now)
          ndey = INDEX(TRIM(yrecfm),'=')
          idx = 1
          IF (ndey /= 0) THEN
            var_calc = yrecfm(1:ndey-1)
            DO WHILE (ndey /= 0)
              IF (idx>MAXRAW) THEN
                CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','MAXRAW exceeded (too many raw variables for 1 computed one)')
              END IF
              yrecfm = yrecfm(ndey+1:)
              ndey = INDEX(TRIM(yrecfm),'+')
              IF (ndey /= 0) THEN
                var_raw(idx) = yrecfm(1:ndey-1)
              ELSE
                var_raw(idx) = TRIM(yrecfm)
              END IF
              idx = idx + 1
            END DO

            tpreclist(idx_var)%name = trim(var_calc)
            tpreclist(idx_var)%calc = .TRUE.
            tpreclist(idx_var)%tbw  = .TRUE.
            tpreclist(idx_var)%tbr  = .FALSE.
            tpreclist(idx_var)%NSRC = idx-1
            idx_var=idx_var+1
            DO jj = 1, idx-1
              tpreclist(idx_var-jj)%src(jj) = idx_var
              tpreclist(idx_var)%name = trim(var_raw(jj))
              tpreclist(idx_var)%calc = .FALSE.
              tpreclist(idx_var)%tbw  = .FALSE.
              tpreclist(idx_var)%tbr  = .TRUE.
              tpreclist(idx_var)%tgt  = idx_var-jj
              idx_var=idx_var+1
            END DO

          ELSE
            tpreclist(idx_var)%name = trim(yrecfm)
            tpreclist(idx_var)%calc = .FALSE.
            tpreclist(idx_var)%tbw  = .TRUE.
            idx_var=idx_var+1

          END IF

          ndb = nde+ndb
       END DO

       DO ji=1,nbvar_tbr+nbvar_calc
          IF (tpreclist(ji)%calc) CYCLE

          yrecfm = TRIM(tpreclist(ji)%name)
          IF (INFILES(1)%TFILE%CFORMAT == 'LFI') THEN
            CALL LFINFO(iresp2,ilu,trim(yrecfm),ileng,ipos)
            IF (iresp2 == 0 .AND. ileng /= 0) THEN
              tpreclist(ji)%found = .true.
              tpreclist(ji)%NSIZE = ileng - 2 - NLFIMAXCOMMENTLENGTH
            END IF

            IF (iresp2==0 .AND. ileng == 0 .AND. ipos==0 .AND. INFILES(1)%TFILE%NSUBFILES_IOZ>0) THEN
              !Variable not found with no error (iresp2==0 .AND. ileng == 0 .AND. ipos==0)
              !If we are merging, maybe it is one of the split variable
              !In that case, the 1st part of the variable is in the 1st split file with a 0001 suffix
              CALL LFINFO(iresp2,INFILES(1)%TFILE%TFILES_IOZ(1)%TFILE%NLFIFLU,trim(yrecfm)//'0001',ileng,ipos)
              IF (iresp2 == 0 .AND. ileng /= 0) THEN
                tpreclist(ji)%found  = .true.
                tpreclist(ji)%LSPLIT = .true.
                IF (tpreclist(ji)%tgt > 0) THEN !If this variable is used for a calculated one
                  tpreclist(tpreclist(ji)%tgt)%LSPLIT = .true.
                END IF
              END IF
              tpreclist(ji)%NSIZE = (ileng - 2 - NLFIMAXCOMMENTLENGTH) * (NKMAX+2*JPVEXT)
              ileng = tpreclist(ji)%NSIZE + 2 + NLFIMAXCOMMENTLENGTH
            END IF

            leng = ileng
          ELSE IF (INFILES(1)%TFILE%CFORMAT == 'NETCDF4') THEN
            status = NF90_INQ_VARID(kcdf_id,trim(yrecfm),var_id)
            IF (status /= NF90_NOERR .AND. INFILES(1)%TFILE%NSUBFILES_IOZ>0) THEN
              !Variable probably not found (other error possible...)
              !If we are merging, maybe it is one of the split variable
              !In that case, the 1st part of the variable is in the 1st split file with a 0001 suffix
              kcdf_id2 = INFILES(1)%TFILE%TFILES_IOZ(1)%TFILE%NNCID
              status = NF90_INQ_VARID(kcdf_id2,trim(yrecfm)//'0001',var_id)
              IF (status == NF90_NOERR) THEN
                tpreclist(ji)%LSPLIT = .true.
                IF (tpreclist(ji)%tgt > 0) THEN !If this variable is used for a calculated one
                  tpreclist(tpreclist(ji)%tgt)%LSPLIT = .true.
                END IF
              ELSE
                CALL HANDLE_ERR(status,__LINE__)
              END IF
            ELSE IF (status /= NF90_NOERR) THEN
              CALL HANDLE_ERR(status,__LINE__)
            ELSE
              kcdf_id2 = kcdf_id
            ENDIF
            !
            IF (status == NF90_NOERR) THEN
              tpreclist(ji)%found = .true.
              CALL IO_GET_METADATA_NC4(kcdf_id2,var_id,tpreclist(ji))
            END IF
          END IF

          IF (.NOT.tpreclist(ji)%found) THEN
            CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','variable '//TRIM(yrecfm)//' not found => ignored')
             tpreclist(ji)%tbw   = .FAlSE.
             tpreclist(ji)%tbr   = .FAlSE.
          END IF
       END DO

       maxvar = nbvar_tbr+nbvar_calc
DO ji=1,nbvar_tbr+nbvar_calc
  print *,ji,'name=',trim(tpreclist(ji)%name),' calc=',tpreclist(ji)%calc,' tbw=',tpreclist(ji)%tbw,&
          ' tbr=',tpreclist(ji)%tbr,' found=',tpreclist(ji)%found
END DO

    ELSE
       ! Entire file is converted
       ALLOCATE(tpreclist(nbvar_infile))
       DO ji=1,nbvar_infile
         tpreclist(ji)%calc   = .FALSE. !By default variables are not computed from others
         tpreclist(ji)%tbw    = .TRUE.  !By default variables are written
         tpreclist(ji)%src(:) = -1
       END DO

       IF (INFILES(1)%TFILE%CFORMAT == 'LFI') THEN
         CALL LFIPOS(iresp2,ilu)
         ladvan = .TRUE.

         DO ji=1,nbvar_infile
           CALL LFICAS(iresp2,ilu,yrecfm,ileng,ipos,ladvan)
           tpreclist(ji)%name   = trim(yrecfm)
           tpreclist(ji)%found  = .TRUE.
           tpreclist(ji)%NSIZE  = ileng - 2 - NLFIMAXCOMMENTLENGTH

           !Detect if date variable
           IDXDATE = INDEX(trim(yrecfm),"%TDATE",.TRUE.)
           IDXTIME = INDEX(trim(yrecfm),"%TIME", .TRUE.)
           IF (IDXDATE/=0 .AND. IDXTIME/=0) &
             CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','field in LFI file with %TDATE and %TIME in name '//TRIM(YRECFM))
           IDX = MAX(IDXDATE,IDXTIME)
           IF (IDX>0) THEN
             YDATENAME = YRECFM(1:IDX-1)
             !Look if datename is already known
             ICURDATE = 0
             DO JJ=1,IDATES
               IF (TRIM(YDATENAME)==TRIM(TLFIDATES(JJ)%CNAME)) THEN
                 ICURDATE = JJ
                 EXIT
               END IF
             END DO
             !
             IF (ICURDATE == 0) THEN
               !New date name detected
               IDATES = IDATES + 1
               IF (IDATES>MAXDATES) CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','too many dates, increase MAXDATES')
               ICURDATE = IDATES
             END IF
             TLFIDATES(ICURDATE)%CNAME = TRIM(YDATENAME)
             IF (IDXTIME>0) THEN
               IF (TLFIDATES(ICURDATE)%NIDX_TIME /= -1) &
                 CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','NIDX_TIME already set for '//TRIM(YDATENAME))
               TLFIDATES(ICURDATE)%NIDX_TIME = JI
               !Set variable name to truncated name (necessary to correctly identify the variable when read)
               tpreclist(ji)%name = TRIM(YDATENAME)
             END IF
             IF (IDXDATE>0) THEN
               IF (TLFIDATES(ICURDATE)%NIDX_DATE /= -1) &
                 CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','NIDX_DATE already set for '//TRIM(YDATENAME))
               TLFIDATES(ICURDATE)%NIDX_DATE = JI
               !Do not treat this variable (the date part will be read with the time part)
               tpreclist(ji)%name  = 'removed_date'
               tpreclist(ji)%tbw   = .FALSE.
               tpreclist(ji)%tbr   = .FALSE.
               tpreclist(ji)%found = .FALSE.
             END IF
           END IF
         END DO
         !
         DO JI=1,IDATES
           IF (TLFIDATES(JI)%NIDX_DATE == -1 .OR. TLFIDATES(JI)%NIDX_TIME == -1) &
             CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','incomplete DATE/TIME fields for '//TRIM(TLFIDATES(JI)%CNAME))
         END DO

         !
       ELSE IF (INFILES(1)%TFILE%CFORMAT == 'NETCDF4') THEN
         DO ji=1,nbvar_infile
           var_id = ji
           status = NF90_INQUIRE_VARIABLE(kcdf_id,var_id, name = tpreclist(ji)%name)
           IF (status /= NF90_NOERR) CALL HANDLE_ERR(status,__LINE__)
           tpreclist(ji)%found  = .TRUE.
           CALL IO_GET_METADATA_NC4(kcdf_id,var_id,tpreclist(ji))
         END DO
       END IF

       maxvar = nbvar_infile
    END IF

    ! Check if variable is in TFIELDLIST and populate corresponding metadata
    DO ji=1,maxvar
      IF (runmode/=MODECDF2LFI .AND. options(OPTSPLIT)%set .AND. tpreclist(ji)%tbw) idx_out = idx_out + 1

      IF (.NOT.tpreclist(ji)%found .OR. tpreclist(ji)%calc ) CYCLE
      !
      !Do not treat dimension variables (they are automatically added when creating netCDF file)
      IF (      tpreclist(ji)%name == 'ni'          &
           .OR. tpreclist(ji)%name == 'nj'          &
           .OR. tpreclist(ji)%name == 'ni_u'        &
           .OR. tpreclist(ji)%name == 'nj_u'        &
           .OR. tpreclist(ji)%name == 'ni_v'        &
           .OR. tpreclist(ji)%name == 'nj_v'        &
           .OR. tpreclist(ji)%name == 'latitude'    &
           .OR. tpreclist(ji)%name == 'longitude'   &
           .OR. tpreclist(ji)%name == 'latitude_u'  &
           .OR. tpreclist(ji)%name == 'longitude_u' &
           .OR. tpreclist(ji)%name == 'latitude_v'  &
           .OR. tpreclist(ji)%name == 'longitude_v' &
           .OR. tpreclist(ji)%name == 'latitude_f'  &
           .OR. tpreclist(ji)%name == 'longitude_f' &
           .OR. tpreclist(ji)%name == 'level'       &
           .OR. tpreclist(ji)%name == 'level_w'     &
           .OR. tpreclist(ji)%name == 'time'        ) THEN
        tpreclist(ji)%tbw   = .FALSE.
        tpreclist(ji)%tbr   = .FALSE.
        tpreclist(ji)%found = .FALSE.
      ELSE
        CALL FIND_FIELD_ID_FROM_MNHNAME(tpreclist(ji)%name,IID,IRESP,ONOWARNING=.TRUE.)
        IF (IRESP==0) THEN
          tpreclist(ji)%TFIELD = TFIELDLIST(IID)
          ! Determine TDIMS
          IF (runmode==MODELFI2CDF) THEN
            ALLOCATE(tpreclist(ji)%TDIMS(tpreclist(ji)%TFIELD%NDIMS))
            CALL IO_GUESS_DIMIDS_NC4(outfiles(idx_out)%TFILE,tpreclist(ji)%TFIELD,&
                                     tpreclist(ji)%NSIZE,tpreclist(ji)%TDIMS,IRESP)
          ELSE !If we read netCDF4, we already have all necessary data
            !Special case for EMIS (only the first band is read/written) -> NDIMS reduced to 2
            if(tpreclist(ji)%TFIELD%CMNHNAME=="EMIS") tpreclist(ji)%TFIELD%NDIMS = 2

            CALL IO_FILL_DIMS_NC4(outfiles(idx_out)%TFILE,tpreclist(ji),IRESP)
          ENDIF
          IF (IRESP/=0) THEN
            CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','can not guess dimensions for '//tpreclist(ji)%TFIELD%CMNHNAME// &
                                         ' => ignored')
            tpreclist(ji)%tbw   = .FALSE.
            tpreclist(ji)%tbr   = .FALSE.
            tpreclist(ji)%found = .FALSE.
            CYCLE
          END IF
        ELSE !Field not found in list, try to determine characteristics
          tpreclist(ji)%TFIELD%CMNHNAME = TRIM(tpreclist(ji)%name)
          tpreclist(ji)%TFIELD%CSTDNAME = ''
          tpreclist(ji)%TFIELD%CLONGNAME = TRIM(tpreclist(ji)%name)
          tpreclist(ji)%TFIELD%CUNITS = ''
          tpreclist(ji)%TFIELD%CDIR = 'XY' !Assumption...
          tpreclist(ji)%TFIELD%CLBTYPE = 'NONE'
          tpreclist(ji)%TFIELD%CCOMMENT = ''
          !
          IF (runmode==MODELFI2CDF) THEN
            tpreclist(ji)%TFIELD%NGRID = 1 !Assumption
            tpreclist(ji)%TFIELD%NTYPE = TYPEREAL !Assumption
            WRITE(YTYPE,'( A )') 'REAL (forced)'
            IF (tpreclist(ji)%NSIZE>1) THEN
              ALLOCATE(tpreclist(ji)%TDIMS(3))
              ! Determine TDIMS
              CALL PRINT_MSG(NVERB_DEBUG,'IO','parse_infiles',tpreclist(ji)%TFIELD%CMNHNAME//': try 3D')
              tpreclist(ji)%TFIELD%NDIMS = 3 !Try with 3D
              CALL IO_GUESS_DIMIDS_NC4(outfiles(idx_out)%TFILE,tpreclist(ji)%TFIELD,&
                                      tpreclist(ji)%NSIZE,tpreclist(ji)%TDIMS,IRESP)
              !
              IF (IRESP/=0 .OR. tpreclist(ji)%TDIMS(3)%LEN==1) THEN
                CALL PRINT_MSG(NVERB_DEBUG,'IO','parse_infiles',tpreclist(ji)%TFIELD%CMNHNAME//': try 2D')
                !Try again with 2D
                tpreclist(ji)%TFIELD%NDIMS = 2
                CALL IO_GUESS_DIMIDS_NC4(outfiles(idx_out)%TFILE,tpreclist(ji)%TFIELD,&
                                        tpreclist(ji)%NSIZE,tpreclist(ji)%TDIMS,IRESP)
              END IF
              !
              IF (IRESP/=0 .OR. tpreclist(ji)%TDIMS(2)%LEN==1) THEN
                CALL PRINT_MSG(NVERB_DEBUG,'IO','parse_infiles',tpreclist(ji)%TFIELD%CMNHNAME//': try 1D')
                !Try again with 1D
                tpreclist(ji)%TFIELD%NDIMS = 1
                tpreclist(ji)%TFIELD%CDIR = '--' !Assumption...
                CALL IO_GUESS_DIMIDS_NC4(outfiles(idx_out)%TFILE,tpreclist(ji)%TFIELD,&
                                        tpreclist(ji)%NSIZE,tpreclist(ji)%TDIMS,IRESP)
              END IF
              !
              IF (IRESP/=0) THEN !Could not find valid characteristics
                tpreclist(ji)%tbw   = .FALSE.
                tpreclist(ji)%tbr   = .FALSE.
                tpreclist(ji)%found = .FALSE.
                CYCLE
              END IF
            ELSE !NSIZE==0
              tpreclist(ji)%TFIELD%CDIR = '--'
              tpreclist(ji)%TFIELD%NDIMS = 0
              tpreclist(ji)%TFIELD%NGRID = 0
            END IF
            tpreclist(ji)%TFIELD%LTIMEDEP = .FALSE. !Assumption
          ELSE ! Input file is netCDF
            tpreclist(ji)%TFIELD%NGRID = tpreclist(ji)%NGRID_FILE
            SELECT CASE(tpreclist(ji)%NTYPE_FILE)
              CASE (NF90_INT1) !NF90_INT1=NF90_BYTE
                tpreclist(ji)%TFIELD%NTYPE = TYPELOG
                tpreclist(ji)%TFIELD%NDIMS = tpreclist(ji)%NDIMS_FILE
                WRITE(YTYPE,'( A )') 'LOGICAL'
              CASE (NF90_CHAR)
                tpreclist(ji)%TFIELD%NTYPE = TYPECHAR
                tpreclist(ji)%TFIELD%NDIMS = tpreclist(ji)%NDIMS_FILE-1
                WRITE(YTYPE,'( A )') 'CHARACTER'
              CASE (NF90_INT,NF90_INT64)
                tpreclist(ji)%TFIELD%NTYPE = TYPEINT
                tpreclist(ji)%TFIELD%NDIMS = tpreclist(ji)%NDIMS_FILE
                WRITE(YTYPE,'( A )') 'INTEGER'
              CASE (NF90_FLOAT,NF90_DOUBLE)
                tpreclist(ji)%TFIELD%NTYPE = TYPEREAL
                tpreclist(ji)%TFIELD%NDIMS = tpreclist(ji)%NDIMS_FILE
                WRITE(YTYPE,'( A )') 'REAL'
              CASE DEFAULT
                tpreclist(ji)%TFIELD%NTYPE = TYPEUNDEF
                tpreclist(ji)%TFIELD%NDIMS = tpreclist(ji)%NDIMS_FILE
                WRITE(YTYPE,'( A )') 'UNKNOWN'
            END SELECT

            tpreclist(ji)%TFIELD%CUNITS = tpreclist(ji)%CUNITS_FILE

            IF (tpreclist(ji)%TFIELD%NDIMS<2) THEN
              tpreclist(ji)%TFIELD%CDIR  = '--' !Assumption
            ELSE
              tpreclist(ji)%TFIELD%CDIR  = 'XY' !Assumption
            END IF

            CALL IO_FILL_DIMS_NC4(outfiles(idx_out)%TFILE,tpreclist(ji),IRESP)

            IF (tpreclist(ji)%NDIMS_FILE>0) THEN
              IF (tpreclist(ji)%CDIMNAMES_FILE(tpreclist(ji)%NDIMS_FILE)=='time') THEN
                tpreclist(ji)%TFIELD%NDIMS = tpreclist(ji)%TFIELD%NDIMS - 1
              END IF
            END IF
            !
            IF (IRESP/=0) THEN
              tpreclist(ji)%tbw   = .FALSE.
              tpreclist(ji)%tbr   = .FALSE.
              tpreclist(ji)%found = .FALSE.
            END IF
          END IF
          !
          IF (runmode==MODELFI2CDF) THEN
            tpreclist(ji)%TFIELD%NGRID = NGRIDUNKNOWN !Assumption
            IF(tpreclist(ji)%TFIELD%NDIMS == 0 .OR. tpreclist(ji)%TFIELD%NTYPE == TYPECHAR) THEN
              tpreclist(ji)%TFIELD%NGRID = 0
            END IF
          END IF
          !
          IF (.NOT.tpreclist(ji)%found) THEN
            CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','can not guess dimensions for '// &
                           TRIM(tpreclist(ji)%TFIELD%CMNHNAME)//' => ignored')
          ELSE
            IF (tpreclist(ji)%TFIELD%LTIMEDEP) THEN
              WRITE(YNDIMS,'( I1 )') tpreclist(ji)%TFIELD%NDIMS-1
              CALL PRINT_MSG(NVERB_WARNING,'IO','unknown field',tpreclist(ji)%TFIELD%CMNHNAME//' seems to be '// &
                             YNDIMS//'D of type '//TRIM(YTYPE)//' (time dependent)')
            ELSE
              WRITE(YNDIMS,'( I1 )') tpreclist(ji)%TFIELD%NDIMS
              CALL PRINT_MSG(NVERB_WARNING,'IO','unknown field',tpreclist(ji)%TFIELD%CMNHNAME//' seems to be '// &
                             YNDIMS//'D of type '//TRIM(YTYPE))
            END IF
          END IF
        END IF
      END IF
    END DO

    IF (nbvar_calc>0) THEN
    !Calculated variables
    !Done after previous loop to reuse metadata from component variables
    !Derive metadata from its components
    !If same value for all components => take it
    !If not => nothing or default value
    !Check sizes: must be the same for all
    DO ji=1,maxvar
      IF (.NOT.tpreclist(ji)%calc ) CYCLE
      !
      tpreclist(ji)%TFIELD%CMNHNAME  = tpreclist(ji)%name
      tpreclist(ji)%TFIELD%CSTDNAME  = ''
      tpreclist(ji)%TFIELD%CLONGNAME = tpreclist(ji)%name
      !
      GOK = .TRUE.
      DO jj=1,tpreclist(ji)%NSRC
        idx_var = tpreclist(ji)%src(jj)
        IF(.NOT.tpreclist(idx_var)%found) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','some components for calculated variable ' &
                         //TRIM(tpreclist(ji)%name)//' are not known => ignored')
          tpreclist(ji)%tbw   = .FALSE.
          tpreclist(ji)%tbr   = .FALSE.
          tpreclist(ji)%found = .FALSE.
          GOK = .FALSE.
          EXIT
        END IF
      END DO
      !
      IF (GOK) THEN
        idx_var = tpreclist(ji)%src(1)
        tpreclist(ji)%TFIELD%CUNITS   = tpreclist(idx_var)%TFIELD%CUNITS
        tpreclist(ji)%TFIELD%CDIR     = tpreclist(idx_var)%TFIELD%CDIR
        tpreclist(ji)%TFIELD%CLBTYPE  = tpreclist(idx_var)%TFIELD%CLBTYPE
        tpreclist(ji)%TFIELD%CCOMMENT = TRIM(tpreclist(ji)%name)//'='//TRIM(tpreclist(idx_var)%name)
        IF (tpreclist(ji)%NSRC>1) tpreclist(ji)%TFIELD%CCOMMENT = TRIM(tpreclist(ji)%TFIELD%CCOMMENT)//'+'
        tpreclist(ji)%TFIELD%NGRID    = tpreclist(idx_var)%TFIELD%NGRID
        tpreclist(ji)%TFIELD%NTYPE    = tpreclist(idx_var)%TFIELD%NTYPE
        tpreclist(ji)%TFIELD%NDIMS    = tpreclist(idx_var)%TFIELD%NDIMS
#if 0
!PW: TODO?
        tpreclist(ji)%TFIELD%NFILLVALUE
        tpreclist(ji)%TFIELD%XFILLVALUE
        tpreclist(ji)%TFIELD%NVALIDMIN
        tpreclist(ji)%TFIELD%NVALIDMAX
        tpreclist(ji)%TFIELD%XVALIDMIN
        tpreclist(ji)%TFIELD%XVALIDMAX
#endif
        DO jj=2,tpreclist(ji)%NSRC
          idx_var = tpreclist(ji)%src(jj)
          !
          IF (tpreclist(ji)%TFIELD%CUNITS /= tpreclist(idx_var)%TFIELD%CUNITS) THEN
            CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','CUNITS is not uniform between components of calculated variable '&
                           //TRIM(tpreclist(ji)%name)//' => CUNITS not set')
            tpreclist(ji)%TFIELD%CUNITS = ''
          END IF
          !
          IF (tpreclist(ji)%TFIELD%CDIR /= tpreclist(idx_var)%TFIELD%CDIR) THEN
            CALL PRINT_MSG(NVERB_ERROR,'IO','parse_infiles','CDIR is not uniform between components of calculated variable '&
                           //TRIM(tpreclist(ji)%name)//' => CDIR=--')
            tpreclist(ji)%TFIELD%CDIR = '--'
          END IF
          !
          IF (tpreclist(ji)%TFIELD%CLBTYPE /= tpreclist(idx_var)%TFIELD%CLBTYPE) THEN
            CALL PRINT_MSG(NVERB_ERROR,'IO','parse_infiles','CLBTYPE is not uniform between components of calculated variable '&
                           //TRIM(tpreclist(ji)%name)//' => CLBTYPE=NONE')
            tpreclist(ji)%TFIELD%CLBTYPE = 'NONE'
          END IF
          !
          tpreclist(ji)%TFIELD%CCOMMENT = TRIM(tpreclist(ji)%TFIELD%CCOMMENT)//TRIM(tpreclist(idx_var)%name)
          IF (jj<tpreclist(ji)%NSRC) tpreclist(ji)%TFIELD%CCOMMENT = TRIM(tpreclist(ji)%TFIELD%CCOMMENT)//'+'
          !
          IF (tpreclist(ji)%TFIELD%NGRID /= tpreclist(idx_var)%TFIELD%NGRID) THEN
            CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','NGRID is not uniform between components of calculated variable '&
                           //TRIM(tpreclist(ji)%name)//' => NGRID=1')
            tpreclist(ji)%TFIELD%NGRID = 1
          END IF
          !
          IF (tpreclist(ji)%TFIELD%NTYPE /= tpreclist(idx_var)%TFIELD%NTYPE) THEN
            CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','NTYPE is not uniform between components of calculated variable '&
                           //TRIM(tpreclist(ji)%name))
            tpreclist(ji)%TFIELD%NTYPE = TYPEUNDEF
          END IF
          !
          IF (tpreclist(ji)%TFIELD%NDIMS /= tpreclist(idx_var)%TFIELD%NDIMS) THEN
            CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','NDIMS is not uniform between components of calculated variable '&
                           //TRIM(tpreclist(ji)%name))
          END IF
        END DO
        !
        ALLOCATE(tpreclist(ji)%TDIMS(tpreclist(ji)%TFIELD%NDIMS))
        tpreclist(ji)%TDIMS = tpreclist(idx_var)%TDIMS
        !
      END IF
    END DO !ji=1,maxvar
    END IF !nbvar_calc>0

  END SUBROUTINE parse_infiles
  
  SUBROUTINE HANDLE_ERR(status,line)
    INTEGER :: status,line

    IF (status /= NF90_NOERR) THEN
       PRINT *, 'line ',line,': ',NF90_STRERROR(status)
           STOP
    END IF
  END SUBROUTINE HANDLE_ERR

  SUBROUTINE def_ncdf(infiles,outfiles,KNFILES_OUT)
    USE MODD_CONF,   ONLY: NMNHVERSION
    use mode_io_write_nc4, only: io_write_header_nc4

    TYPE(TFILE_ELT),DIMENSION(:),INTENT(IN) :: infiles
    TYPE(TFILE_ELT),DIMENSION(:),INTENT(IN) :: outfiles
    INTEGER,                     INTENT(IN) :: KNFILES_OUT

    CHARACTER(LEN=*),PARAMETER :: YUNKNOWNHIST = 'Previous history is unknown'

    CHARACTER(LEN=16) :: YMNHVERSION
    CHARACTER(LEN=:),ALLOCATABLE :: YHISTORY
    INTEGER :: ilen, ji
    INTEGER(KIND=IDCDF_KIND) :: status
    INTEGER(KIND=IDCDF_KIND) :: kcdf_id


    CALL PRINT_MSG(NVERB_DEBUG,'IO','def_ncdf','called')

    !Copy history attribute for netCDF files
    IF( outfiles(1)%TFILE%CFORMAT=='NETCDF4' ) THEN
      IF( infiles(1)%TFILE%CFORMAT=='NETCDF4' ) THEN
        status = NF90_INQUIRE_ATTRIBUTE(infiles(1)%TFILE%NNCID, NF90_GLOBAL, 'history', LEN=ilen)
        IF (status == NF90_NOERR) THEN
          ALLOCATE(CHARACTER(LEN=ilen) :: YHISTORY)
          status = NF90_GET_ATT(infiles(1)%TFILE%NNCID, NF90_GLOBAL, 'history', YHISTORY)
        ELSE
          YHISTORY = YUNKNOWNHIST
        END IF
      ELSE
        YHISTORY = YUNKNOWNHIST
      END IF

      DO ji = 1,KNFILES_OUT
        kcdf_id = outfiles(ji)%TFILE%NNCID
        status = NF90_PUT_ATT(kcdf_id,NF90_GLOBAL,'history',YHISTORY)
        IF (status /= NF90_NOERR) CALL HANDLE_ERR(status,__LINE__)
      END DO
    END IF

    !Write header for netCDF files
    DO ji = 1,KNFILES_OUT
      kcdf_id = outfiles(ji)%TFILE%NNCID

      ! global attributes
      CALL IO_WRITE_HEADER_NC4(outfiles(ji)%TFILE)
      !
      WRITE(YMNHVERSION,"( I0,'.',I0,'.',I0 )" ) NMNHVERSION(1),NMNHVERSION(2),NMNHVERSION(3)
      status = NF90_PUT_ATT(kcdf_id,NF90_GLOBAL,'lfi2cdf_version',TRIM(YMNHVERSION))
      IF (status /= NF90_NOERR) CALL HANDLE_ERR(status,__LINE__)
    END DO

  END SUBROUTINE def_ncdf

  SUBROUTINE fill_files(infiles,outfiles,tpreclist,knaf,options)
    USE MODD_TYPE_DATE

    TYPE(TFILE_ELT),DIMENSION(:), INTENT(IN)    :: infiles
    TYPE(TFILE_ELT),DIMENSION(:), INTENT(IN)    :: outfiles
    TYPE(workfield), DIMENSION(:),INTENT(INOUT) :: tpreclist
    INTEGER,                      INTENT(IN)    :: knaf
    TYPE(option),DIMENSION(:),    INTENT(IN)    :: options

    INTEGER                                  :: idx, ji, jj
    INTEGER                                  :: IDIMS
    INTEGER                                  :: INSRC
    INTEGER                                  :: ISRC
    INTEGER(KIND=IDCDF_KIND),DIMENSION(NF90_MAX_VAR_DIMS) :: IDIMLEN
    logical,dimension(knaf)                  :: gtimedep_in, gtimedep_out

    CHARACTER(LEN=:),       ALLOCATABLE :: YTAB0D
    INTEGER,DIMENSION(:),   ALLOCATABLE :: ITAB1D, ITAB1D2
    INTEGER,DIMENSION(:,:), ALLOCATABLE :: ITAB2D, ITAB2D2
    LOGICAL,DIMENSION(:),   ALLOCATABLE :: GTAB1D
    REAL,DIMENSION(:),      ALLOCATABLE :: XTAB1D, XTAB1D2
    REAL,DIMENSION(:,:),    ALLOCATABLE :: XTAB2D, XTAB2D2
    REAL,DIMENSION(:,:,:),  ALLOCATABLE :: XTAB3D, XTAB3D2
    REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: XTAB4D, XTAB4D2
    TYPE(DATE_TIME)                     :: TZDATE


    CALL PRINT_MSG(NVERB_DEBUG,'IO','fill_files','called')

    ! For versions of MesoNH <5.4.0, fields were not stored with a time dimension
    ! ->necessary to remove it when reading and to restore to the correct one when writing
    if( infiles(1)%TFILE%NMNHVERSION(1)<5 .OR. &
       (infiles(1)%TFILE%NMNHVERSION(1)==5 .AND. infiles(1)%TFILE%NMNHVERSION(2)<4) ) then
      gtimedep_in(:)  = .false.
    else
      gtimedep_in(:)  = tpreclist(:)%TFIELD%LTIMEDEP
    end if
    gtimedep_out(:)   = tpreclist(:)%TFIELD%LTIMEDEP

    idx = 1
    DO ji=1,knaf
      IF (.NOT.tpreclist(ji)%tbw) CYCLE

      IDIMS = tpreclist(ji)%TFIELD%NDIMS

      SELECT CASE(tpreclist(ji)%TFIELD%NTYPE)
      CASE (TYPEINT)
        IDIMLEN(1:IDIMS) = tpreclist(ji)%TDIMS(1:IDIMS)%LEN

        IF (.NOT.tpreclist(ji)%calc) THEN
          INSRC = 1
          ISRC = ji
        ELSE
          INSRC = tpreclist(ji)%NSRC
          ISRC  = tpreclist(ji)%src(1)
        END IF

        tpreclist(ISRC)%TFIELD%LTIMEDEP = gtimedep_in(ISRC)
        SELECT CASE(IDIMS)
        CASE (0)
          ALLOCATE(ITAB1D(1))
          IF (tpreclist(ji)%calc) ALLOCATE(ITAB1D2(1))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,ITAB1D(1))
        CASE (1)
          ALLOCATE(ITAB1D(IDIMLEN(1)))
          IF (tpreclist(ji)%calc) ALLOCATE(ITAB1D2(IDIMLEN(1)))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,ITAB1D)
        CASE (2)
          ALLOCATE(ITAB2D(IDIMLEN(1),IDIMLEN(2)))
          IF (tpreclist(ji)%calc) ALLOCATE(ITAB2D2(IDIMLEN(1),IDIMLEN(2)))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,ITAB2D)
        CASE DEFAULT
          CALL PRINT_MSG(NVERB_WARNING,'IO','fill_files','too many dimensions for ' &
                         //TRIM(tpreclist(ISRC)%name)//' => ignored')
          CYCLE
        END SELECT

        DO JJ=2,INSRC
          ISRC = tpreclist(ji)%src(jj)
          tpreclist(ISRC)%TFIELD%LTIMEDEP = gtimedep_in(ISRC)

          SELECT CASE(IDIMS)
          CASE (0)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,ITAB1D2(1))
            ITAB1D(1) = ITAB1D(1) + ITAB1D2(1)
          CASE (1)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,ITAB1D2)
            ITAB1D(:) = ITAB1D(:) + ITAB1D2(:)
          CASE (2)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,ITAB2D2)
            ITAB2D(:,:) = ITAB2D(:,:) + ITAB2D2(:,:)
          END SELECT
        END DO

        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_out(ji)
        SELECT CASE(IDIMS)
        CASE (0)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,ITAB1D(1))
          DEALLOCATE(ITAB1D)
          IF (tpreclist(ji)%calc) DEALLOCATE(ITAB1D2)
        CASE (1)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,ITAB1D)
          DEALLOCATE(ITAB1D)
          IF (tpreclist(ji)%calc) DEALLOCATE(ITAB1D2)
        CASE (2)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,ITAB2D)
          DEALLOCATE(ITAB2D)
          IF (tpreclist(ji)%calc) DEALLOCATE(ITAB2D2)
        END SELECT


      CASE (TYPELOG)
        IDIMLEN(1:IDIMS) = tpreclist(ji)%TDIMS(1:IDIMS)%LEN

        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_in(ji)
        SELECT CASE(IDIMS)
        CASE (0)
          ALLOCATE(GTAB1D(1))
          CALL IO_READ_FIELD (INFILES(1)%TFILE,   tpreclist(ji)%TFIELD,GTAB1D(1))
          tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_out(ji)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,GTAB1D(1))
          DEALLOCATE(GTAB1D)
        CASE (1)
          ALLOCATE(GTAB1D(IDIMLEN(1)))
          CALL IO_READ_FIELD (INFILES(1)%TFILE,   tpreclist(ji)%TFIELD,GTAB1D)
          tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_out(ji)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,GTAB1D)
          DEALLOCATE(GTAB1D)
        CASE DEFAULT
          CALL PRINT_MSG(NVERB_WARNING,'IO','fill_files','too many dimensions for ' &
                         //TRIM(tpreclist(ji)%name)//' => ignored')
          CYCLE
        END SELECT


      CASE (TYPEREAL)
        IDIMLEN(1:IDIMS) = tpreclist(ji)%TDIMS(1:IDIMS)%LEN

        IF (.NOT.tpreclist(ji)%calc) THEN
          INSRC = 1
          ISRC = ji
        ELSE
          INSRC = tpreclist(ji)%NSRC
          ISRC  = tpreclist(ji)%src(1)
        END IF

        tpreclist(ISRC)%TFIELD%LTIMEDEP = gtimedep_in(ISRC)
        SELECT CASE(IDIMS)
        CASE (0)
          ALLOCATE(XTAB1D(1))
          IF (tpreclist(ji)%calc) ALLOCATE(XTAB1D2(1))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB1D(1))
        CASE (1)
          ALLOCATE(XTAB1D(IDIMLEN(1)))
          IF (tpreclist(ji)%calc) ALLOCATE(XTAB1D2(IDIMLEN(1)))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB1D)
        CASE (2)
          ALLOCATE(XTAB2D(IDIMLEN(1),IDIMLEN(2)))
          IF (tpreclist(ji)%calc) ALLOCATE(XTAB2D2(IDIMLEN(1),IDIMLEN(2)))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB2D)
        CASE (3)
          ALLOCATE(XTAB3D(IDIMLEN(1),IDIMLEN(2),IDIMLEN(3)))
          IF (tpreclist(ji)%calc) ALLOCATE(XTAB3D2(IDIMLEN(1),IDIMLEN(2),IDIMLEN(3)))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB3D)
        CASE (4)
          ALLOCATE(XTAB4D(IDIMLEN(1),IDIMLEN(2),IDIMLEN(3),IDIMLEN(4)))
          IF (tpreclist(ji)%calc) ALLOCATE(XTAB4D2(IDIMLEN(1),IDIMLEN(2),IDIMLEN(3),IDIMLEN(4)))
          CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB4D)
        CASE DEFAULT
          CALL PRINT_MSG(NVERB_WARNING,'IO','fill_files','too many dimensions for ' &
                         //TRIM(tpreclist(ISRC)%name)//' => ignored')
          CYCLE
        END SELECT

        DO JJ=2,INSRC
          ISRC = tpreclist(ji)%src(jj)
          tpreclist(ISRC)%TFIELD%LTIMEDEP = gtimedep_in(ISRC)

          SELECT CASE(IDIMS)
          CASE (0)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB1D2(1))
            XTAB1D(1) = XTAB1D(1) + XTAB1D2(1)
          CASE (1)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB1D2)
            XTAB1D(:) = XTAB1D(:) + XTAB1D2(:)
          CASE (2)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB2D2)
            XTAB2D(:,:) = XTAB2D(:,:) + XTAB2D2(:,:)
          CASE (3)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB3D2)
            XTAB3D(:,:,:) = XTAB3D(:,:,:) + XTAB3D2(:,:,:)
          CASE (4)
            CALL IO_READ_FIELD(INFILES(1)%TFILE,tpreclist(ISRC)%TFIELD,XTAB4D2)
            XTAB4D(:,:,:,:) = XTAB4D(:,:,:,:) + XTAB4D2(:,:,:,:)
          END SELECT
        END DO

        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_out(ji)
        SELECT CASE(IDIMS)
        CASE (0)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,XTAB1D(1))
          DEALLOCATE(XTAB1D)
          IF (tpreclist(ji)%calc) DEALLOCATE(XTAB1D2)
        CASE (1)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,XTAB1D)
          DEALLOCATE(XTAB1D)
          IF (tpreclist(ji)%calc) DEALLOCATE(XTAB1D2)
        CASE (2)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,XTAB2D)
          DEALLOCATE(XTAB2D)
          IF (tpreclist(ji)%calc) DEALLOCATE(XTAB2D2)
        CASE (3)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,XTAB3D)
          DEALLOCATE(XTAB3D)
          IF (tpreclist(ji)%calc) DEALLOCATE(XTAB3D2)
        CASE (4)
          CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,XTAB4D)
          DEALLOCATE(XTAB4D)
          IF (tpreclist(ji)%calc) DEALLOCATE(XTAB4D2)
        END SELECT


      CASE (TYPECHAR)
        ISRC = ji

        IF (IDIMS/=0) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','fill_files','too many dimensions for ' &
                         //TRIM(tpreclist(ISRC)%name)//' => ignored')
          CYCLE
        END IF

        ALLOCATE(CHARACTER(LEN=tpreclist(ji)%NSIZE)::YTAB0D)
        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_in(ji)
        CALL IO_READ_FIELD (INFILES(1)%TFILE,   tpreclist(ji)%TFIELD,YTAB0D)
        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_out(ji)
        CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,YTAB0D)
        DEALLOCATE(YTAB0D)


      CASE (TYPEDATE)
        ISRC = ji

        IF (IDIMS/=0) THEN
          CALL PRINT_MSG(NVERB_WARNING,'IO','fill_files','too many dimensions for ' &
                         //TRIM(tpreclist(ISRC)%name)//' => ignored')
          CYCLE
        END IF
        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_in(ji)
        CALL IO_READ_FIELD (INFILES(1)%TFILE,   tpreclist(ji)%TFIELD%CMNHNAME,TZDATE)
        tpreclist(ji)%TFIELD%LTIMEDEP = gtimedep_out(ji)
        CALL IO_WRITE_FIELD(outfiles(idx)%TFILE,tpreclist(ji)%TFIELD,TZDATE)


      CASE default
        ISRC = ji

        CALL PRINT_MSG(NVERB_WARNING,'IO','fill_files','invalid datatype for ' &
                       //TRIM(tpreclist(ISRC)%name)//' => ignored')

      END SELECT

      if (options(OPTSPLIT)%set) idx = idx + 1
    END DO
  END SUBROUTINE fill_files


  SUBROUTINE OPEN_FILES(infiles,outfiles,KNFILES_OUT,hinfile,houtfile,nbvar_infile,options,runmode)
    USE MODD_CONF,          ONLY: LCARTESIAN
    USE MODD_CONF_n,        ONLY: CSTORAGE_TYPE
    USE MODD_DIM_n,         ONLY: NIMAX_ll, NJMAX_ll, NKMAX
    USE MODD_GRID,          ONLY: XBETA, XRPK, XLAT0, XLON0, XLATORI, XLONORI
    USE MODD_GRID_n,        ONLY: LSLEVE, XXHAT, XYHAT, XZHAT
    USE MODD_IO_ll,         ONLY: LIOCDF4
    USE MODD_PARAMETERS,    ONLY: JPHEXT
    USE MODD_PARAMETERS_ll, ONLY: JPHEXT_ll=>JPHEXT, JPVEXT_ll=>JPVEXT
    USE MODD_TIME_n,        ONLY: TDTCUR, TDTMOD

    USE MODE_FM,               ONLY: IO_FILE_OPEN_ll, IO_FILE_CLOSE_ll
    USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST

    TYPE(TFILE_ELT),DIMENSION(:),INTENT(OUT) :: infiles
    TYPE(TFILE_ELT),DIMENSION(:),INTENT(OUT) :: outfiles
    INTEGER,                     INTENT(OUT) :: KNFILES_OUT
    CHARACTER(LEN=*),            INTENT(IN)  :: hinfile
    CHARACTER(LEN=*),            INTENT(IN)  :: houtfile
    INTEGER,                     INTENT(OUT) :: nbvar_infile
    TYPE(option),DIMENSION(:),   INTENT(IN)  :: options
    INTEGER,                     INTENT(IN)  :: runmode

    INTEGER                     :: idx, IRESP2
    INTEGER(KIND=IDCDF_KIND)    :: omode
    INTEGER(KIND=IDCDF_KIND)    :: status
    INTEGER(KIND=LFI_INT)       :: ilu,iresp


    CALL PRINT_MSG(NVERB_DEBUG,'IO','OPEN_FILES','called')

    KNFILES_OUT = 0
    !
    ! Infiles
    !
    IF (runmode == MODECDF2CDF .OR. runmode == MODECDF2LFI) THEN
       !
       ! NetCDF
       !
       CALL IO_FILE_ADD2LIST(INFILES(1)%TFILE,HINFILE,'UNKNOWN','READ',HFORMAT='NETCDF4')
       CALL IO_FILE_OPEN_ll(INFILES(1)%TFILE)

       nbvar_infile = INFILES(1)%TFILE%NNCNAR
   ELSE
       !
       ! LFI
       !
       CALL IO_FILE_ADD2LIST(INFILES(1)%TFILE,HINFILE,'UNKNOWN','READ', &
                             HFORMAT='LFI',KLFIVERB=0)
       CALL IO_FILE_OPEN_ll(INFILES(1)%TFILE)

       ilu = INFILES(1)%TFILE%NLFIFLU

       nbvar_infile = INFILES(1)%TFILE%NLFININAR

       IF (options(OPTLIST)%set) THEN
          CALL LFILAF(iresp,ilu,lfalse)
          CALL IO_FILE_CLOSE_ll(INFILES(1)%TFILE)
          return
       END IF
   END IF
   !
   !Read problem dimensions and some grid variables (needed to determine domain size and also by IO_FILE_OPEN_ll to create netCDF files)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'JPHEXT',JPHEXT)
   JPHEXT_ll = JPHEXT
   JPVEXT_ll = JPVEXT
   !
   ALLOCATE(NIMAX_ll,NJMAX_ll,NKMAX)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'IMAX',NIMAX_ll)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'JMAX',NJMAX_ll)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'KMAX',NKMAX,IRESP2)
   IF (IRESP2/=0) NKMAX = 0
   !
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'PROGRAM',CPROGRAM_ORIG)
   !
   ALLOCATE(CSTORAGE_TYPE)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'STORAGE_TYPE',CSTORAGE_TYPE)
   !
   ALLOCATE(XXHAT(NIMAX_ll+2*JPHEXT))
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'XHAT',XXHAT)
   ALLOCATE(XYHAT(NJMAX_ll+2*JPHEXT))
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'YHAT',XYHAT)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'CARTESIAN',LCARTESIAN)
   !
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'LAT0',XLAT0)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'LON0',XLON0)
   CALL IO_READ_FIELD(INFILES(1)%TFILE,'BETA',XBETA)
   !
   IF (.NOT.LCARTESIAN) THEN
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'RPK',   XRPK)
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'LATORI',XLATORI)
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'LONORI',XLONORI)
   ENDIF
   !
   IF (TRIM(CPROGRAM_ORIG)/='PGD' .AND. TRIM(CPROGRAM_ORIG)/='NESPGD' .AND. TRIM(CPROGRAM_ORIG)/='ZOOMPG' &
       .AND. .NOT.(TRIM(CPROGRAM_ORIG)=='REAL' .AND. CSTORAGE_TYPE=='SU') ) THEN !condition to detect PREP_SURFEX
     ALLOCATE(XZHAT(NKMAX+2*JPVEXT))
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'ZHAT',XZHAT)
     ALLOCATE(LSLEVE)
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'SLEVE',LSLEVE)
     ALLOCATE(TDTMOD)
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'DTMOD',TDTMOD,IRESP2)
     IF(IRESP2/=0) DEALLOCATE(TDTMOD)
     ALLOCATE(TDTCUR)
     CALL IO_READ_FIELD(INFILES(1)%TFILE,'DTCUR',TDTCUR,IRESP2)
     IF(IRESP2/=0) DEALLOCATE(TDTCUR)
   END IF
   !
   ! Outfiles
   !
   IF (runmode == MODELFI2CDF .OR. runmode == MODECDF2CDF) THEN
       !
       ! NetCDF
       !
       IF (.NOT.options(OPTSPLIT)%set) THEN
         KNFILES_OUT = KNFILES_OUT + 1

         idx = KNFILES_OUT
         CALL IO_FILE_ADD2LIST(outfiles(idx)%TFILE,HOUTFILE,'UNKNOWN','WRITE', &
                               HFORMAT='NETCDF4',OOLD=.TRUE.)
         CALL IO_FILE_OPEN_ll(outfiles(idx)%TFILE,HPROGRAM_ORIG=CPROGRAM_ORIG)

         IF (options(OPTCOMPRESS)%set) THEN
           outfiles(idx)%tfile%LNCCOMPRESS       = .TRUE.
           outfiles(idx)%tfile%NNCCOMPRESS_LEVEL = options(OPTCOMPRESS)%ivalue
         END IF

         IF (options(OPTREDUCE)%set) THEN
           outfiles(idx)%tfile%LNCREDUCE_FLOAT_PRECISION = .TRUE.
         END IF

         status = NF90_SET_FILL(outfiles(idx)%TFILE%NNCID,NF90_NOFILL,omode)
         IF (status /= NF90_NOERR) CALL HANDLE_ERR(status,__LINE__)
       END IF ! .NOT.osplit
    ELSE
       !
       ! LFI
       !
       KNFILES_OUT = KNFILES_OUT + 1
       idx = KNFILES_OUT
       CALL IO_FILE_ADD2LIST(outfiles(idx)%TFILE,houtfile,'UNKNOWN','WRITE', &
                             HFORMAT='LFI',KLFIVERB=0,OOLD=.TRUE.)
       LIOCDF4 = .FALSE. !Necessary to open correctly the LFI file
       CALL IO_FILE_OPEN_ll(outfiles(idx)%TFILE,HPROGRAM_ORIG=CPROGRAM_ORIG)
       LIOCDF4 = .TRUE.
   END IF
   !
   ! Create a dummy netCDF file necessary to manage correctly the netCDF dims
   IF (runmode == MODECDF2LFI) THEN
     KNFILES_OUT = KNFILES_OUT + 1

     idx = KNFILES_OUT
     CALL IO_FILE_ADD2LIST(outfiles(idx)%TFILE,'dummy_file','UNKNOWN','WRITE', &
                           HFORMAT='NETCDF4',OOLD=.TRUE.)
     CALL IO_FILE_OPEN_ll(outfiles(idx)%TFILE,HPROGRAM_ORIG=CPROGRAM_ORIG)
   END IF

   PRINT *,'--> Converted to file: ', TRIM(houtfile)

  END SUBROUTINE OPEN_FILES

  SUBROUTINE OPEN_SPLIT_NCFILES_OUT(outfiles,KNFILES_OUT,houtfile,nbvar,options)
    USE MODE_FM,               ONLY: IO_FILE_OPEN_ll
    USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST

    TYPE(TFILE_ELT),DIMENSION(:),  INTENT(INOUT) :: outfiles
    INTEGER,                       INTENT(OUT)   :: KNFILES_OUT
    CHARACTER(LEN=*),              INTENT(IN)    :: houtfile
    INTEGER,                       INTENT(IN)    :: nbvar
    TYPE(option),DIMENSION(:),     INTENT(IN)    :: options

    CHARACTER(LEN=:),ALLOCATABLE                   :: filename
    CHARACTER(LEN=:),ALLOCATABLE                   :: YLIST
    CHARACTER(LEN=NMNHNAMELGTMAX),DIMENSION(nbvar) :: YVARS
    INTEGER                  :: ji
    INTEGER                  :: idx1, idx2
    INTEGER(KIND=IDCDF_KIND) :: status
    INTEGER(KIND=IDCDF_KIND) :: omode

    CALL PRINT_MSG(NVERB_DEBUG,'IO','OPEN_SPLIT_NCFILES_OUT','called')

    KNFILES_OUT = nbvar
    YLIST = TRIM(options(OPTVAR)%cvalue)

    DO ji = 1,nbvar-1
      idx1 = INDEX(YLIST,',')
      idx2 = INDEX(YLIST,'=')
      IF (idx1/=0) THEN
        IF (idx2/=0 .AND. idx2<idx1) THEN
          YVARS(ji) = YLIST(1:idx2-1)
        ELSE
          YVARS(ji) = YLIST(1:idx1-1)
        END IF
        YLIST = YLIST(idx1+1:)
      ELSE
        CALL PRINT_MSG(NVERB_FATAL,'IO','OPEN_SPLIT_NCFILES_OUT','problem separating variable names')
      END IF
    END DO
    idx2 = INDEX(YLIST,'=')
    IF (idx2>0) THEN
      YVARS(nbvar) = YLIST(1:idx2-1)
    ELSE
      YVARS(nbvar) = YLIST
    END IF

    DO ji = 1,nbvar
      filename = trim(houtfile)//'.'//TRIM(YVARS(ji))
      CALL IO_FILE_ADD2LIST(outfiles(ji)%TFILE,filename,'UNKNOWN','WRITE', &
                            HFORMAT='NETCDF4')
      CALL IO_FILE_OPEN_ll(outfiles(ji)%TFILE,HPROGRAM_ORIG=CPROGRAM_ORIG)

      IF (options(OPTCOMPRESS)%set) THEN
        outfiles(ji)%tfile%LNCCOMPRESS       = .TRUE.
        outfiles(ji)%tfile%NNCCOMPRESS_LEVEL = options(OPTCOMPRESS)%ivalue
      END IF

      IF (options(OPTREDUCE)%set) THEN
        outfiles(ji)%tfile%LNCREDUCE_FLOAT_PRECISION = .TRUE.
      END IF

      status = NF90_SET_FILL(outfiles(ji)%TFILE%NNCID,NF90_NOFILL,omode)
      IF (status /= NF90_NOERR) CALL HANDLE_ERR(status,__LINE__)
    END DO

  END SUBROUTINE OPEN_SPLIT_NCFILES_OUT
  
  SUBROUTINE CLOSE_FILES(filelist,KNFILES)
    USE MODE_FM,    ONLY: IO_FILE_CLOSE_ll

    TYPE(TFILE_ELT),DIMENSION(:),INTENT(INOUT) :: filelist
    INTEGER,                     INTENT(IN)    :: KNFILES
    
    INTEGER :: ji


    CALL PRINT_MSG(NVERB_DEBUG,'IO','CLOSE_FILES','called')

    DO ji=1,KNFILES
      IF (filelist(ji)%TFILE%LOPENED) CALL IO_FILE_CLOSE_ll(filelist(ji)%TFILE,HPROGRAM_ORIG=CPROGRAM_ORIG)
    END DO

  END SUBROUTINE CLOSE_FILES


  SUBROUTINE IO_GET_METADATA_NC4(KFILE_ID,KVAR_ID,TPREC)
    USE MODD_DIM_n,      ONLY: NKMAX
    USE MODD_PARAMETERS, ONLY: JPVEXT

    INTEGER,        INTENT(IN)    :: KFILE_ID
    INTEGER,        INTENT(IN)    :: KVAR_ID
    TYPE(workfield),INTENT(INOUT) :: TPREC

    INTEGER                                  :: ILENG
    INTEGER                                  :: JDIM
    INTEGER(KIND=IDCDF_KIND)                 :: ISTATUS
    INTEGER(KIND=IDCDF_KIND),DIMENSION(NF90_MAX_VAR_DIMS) :: IDIMS_ID

    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_GET_METADATA_NC4','called')

    ISTATUS = NF90_INQUIRE_VARIABLE(KFILE_ID,KVAR_ID,NDIMS = TPREC%NDIMS_FILE, &
                                    XTYPE = TPREC%NTYPE_FILE, DIMIDS = IDIMS_ID)
    IF (ISTATUS /= NF90_NOERR) CALL HANDLE_ERR(ISTATUS,__LINE__)

    IF (.NOT.TPREC%LSPLIT) THEN
      ALLOCATE(TPREC%NDIMSIZES_FILE(TPREC%NDIMS_FILE))
      ALLOCATE(TPREC%CDIMNAMES_FILE(TPREC%NDIMS_FILE))
    ELSE
      ALLOCATE(TPREC%NDIMSIZES_FILE(TPREC%NDIMS_FILE+1))
      ALLOCATE(TPREC%CDIMNAMES_FILE(TPREC%NDIMS_FILE+1))
    END IF

    IF (TPREC%NDIMS_FILE == 0) THEN
      ! Scalar variable
      ILENG = 1
    ELSE
      ! Fill dimensions info
      ILENG = 1
      DO JDIM=1,TPREC%NDIMS_FILE
        ISTATUS = NF90_INQUIRE_DIMENSION(KFILE_ID,IDIMS_ID(JDIM),                    &
                                                   len =  TPREC%NDIMSIZES_FILE(JDIM), &
                                                   name = TPREC%CDIMNAMES_FILE(JDIM)  )
        IF (ISTATUS /= NF90_NOERR) CALL HANDLE_ERR(ISTATUS,__LINE__)
        ILENG = ILENG*TPREC%NDIMSIZES_FILE(JDIM)
      END DO

      IF (TPREC%NDIMS_FILE>0) THEN
        IF (TPREC%CDIMNAMES_FILE(TPREC%NDIMS_FILE)=='time') THEN
          TPREC%TFIELD%LTIMEDEP = .TRUE.
        ELSE
          TPREC%TFIELD%LTIMEDEP = .FALSE.
        END IF
      ELSE
        TPREC%TFIELD%LTIMEDEP = .FALSE.
      END IF

      IF (TPREC%LSPLIT) THEN
        IF(     (.NOT.TPREC%TFIELD%LTIMEDEP .AND.  TPREC%NDIMS_FILE/=2)   &
            .OR. (     TPREC%TFIELD%LTIMEDEP .AND.  TPREC%NDIMS_FILE/=3) ) &
          CALL PRINT_MSG(NVERB_FATAL,'IO','parse_infiles','split variables can only be 3D')
          !Split variables are Z-split
          ILENG = ILENG * (NKMAX+2*JPVEXT)
          !Move time dimension to last (4th) position
          IF (TPREC%TFIELD%LTIMEDEP) THEN
            TPREC%NDIMSIZES_FILE(4) = TPREC%NDIMSIZES_FILE(3)
            TPREC%CDIMNAMES_FILE(4) = TPREC%CDIMNAMES_FILE(3)
          END IF
          !Add vertical dimension
          TPREC%NDIMSIZES_FILE(3) = NKMAX+2*JPVEXT
          TPREC%CDIMNAMES_FILE(3) = 'level' !Could also be 'level_w'
        END IF
      END IF

      TPREC%NSIZE = ILENG

      ISTATUS = NF90_GET_ATT(KFILE_ID,KVAR_ID,'grid',TPREC%NGRID_FILE)
      !On MesoNH versions < 5.4.0, the grid number was stored in 'GRID' instead of 'grid'
      IF (ISTATUS /= NF90_NOERR) ISTATUS = NF90_GET_ATT(KFILE_ID,KVAR_ID,'GRID',TPREC%NGRID_FILE)
      IF (ISTATUS /= NF90_NOERR) TPREC%NGRID_FILE = 0

      ISTATUS = NF90_GET_ATT(KFILE_ID,KVAR_ID,'units',TPREC%CUNITS_FILE)
      IF (ISTATUS /= NF90_NOERR) TPREC%CUNITS_FILE = ''
  END SUBROUTINE IO_GET_METADATA_NC4


  SUBROUTINE IO_FILL_DIMS_NC4(TPFILE,TPREC,KRESP)
    USE MODD_IO_ll,        ONLY: TFILEDATA
    use mode_io_tools_nc4, only: getdimcdf, io_find_dim_byname_nc4

    TYPE(TFILEDATA),INTENT(IN)    :: TPFILE
    TYPE(workfield),INTENT(INOUT) :: TPREC
    INTEGER,        INTENT(OUT)   :: KRESP

    INTEGER              :: JJ
    TYPE(DIMCDF),POINTER :: TZDIMPTR

    CALL PRINT_MSG(NVERB_DEBUG,'IO','IO_FILL_DIMS_NC4','called')

    KRESP = 0

    IF (TPREC%NDIMS_FILE<TPREC%TFIELD%NDIMS) THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO','IO_FILL_DIMS_NC4','less dimensions than expected for '//TRIM(TPREC%TFIELD%CMNHNAME)// &
                                        ' => ignored')
      TPREC%tbw   = .FALSE.
      TPREC%tbr   = .FALSE.
      TPREC%found = .FALSE.
      RETURN
    END IF

    ALLOCATE(TPREC%TDIMS(TPREC%TFIELD%NDIMS))

    DO JJ=1,TPREC%TFIELD%NDIMS
    !DO JJ=1,TPREC%NDIMS_FILE !NDIMS_FILE can be bigger than NDIMS due to time dimension (it can be ignored here)
      CALL IO_FIND_DIM_BYNAME_NC4(TPFILE,TPREC%CDIMNAMES_FILE(JJ),TPREC%TDIMS(JJ),KRESP)
      !If dimension not found => create it
      IF (KRESP/=0)  THEN
        TZDIMPTR => GETDIMCDF(TPFILE,TPREC%NDIMSIZES_FILE(JJ))
        TPREC%TDIMS(JJ) = TZDIMPTR
        KRESP = 0
      END IF
      IF (TRIM(TPREC%TDIMS(JJ)%name)/='time' .AND. &
        TPREC%TDIMS(JJ)%len /= TPREC%NDIMSIZES_FILE(JJ)) THEN
        CALL PRINT_MSG(NVERB_WARNING,'IO','parse_infiles','problem with dimensions for '//TPREC%TFIELD%CMNHNAME)
        KRESP = -3
        EXIT
      END IF
    END DO

  END SUBROUTINE IO_FILL_DIMS_NC4

END MODULE mode_util
