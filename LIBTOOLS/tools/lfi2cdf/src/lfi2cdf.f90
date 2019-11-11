!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
program LFI2CDF
  USE MODD_CONF,          ONLY: CPROGRAM
  USE MODD_CONFZ,         ONLY: NB_PROCIO_R
  USE MODD_DIM_n,         ONLY: NIMAX_ll, NJMAX_ll, NKMAX
  USE MODD_IO_ll,         ONLY: LVERB_OUTLST, LVERB_STDOUT, NIO_ABORT_LEVEL, NIO_VERB, NGEN_ABORT_LEVEL, NGEN_VERB
  USE MODD_PARAMETERS,    ONLY: JPHEXT, JPVEXT
  USE MODD_TIMEZ,         ONLY: TIMEZ

  USE MODE_IO_ll,         ONLY: INITIO_ll, SET_CONFIO_ll
  USE MODE_FIELD,         ONLY: INI_FIELD_LIST
  USE mode_options
  USE MODE_SPLITTINGZ_ll, ONLY: INI_PARAZ_ll
  USE mode_util
  USE MODI_VERSION

  USE MODN_CONFIO, ONLY: LCDF4, LLFIOUT, LLFIREAD

  IMPLICIT NONE 

  INTEGER :: ji
  INTEGER :: nbvar_infile = 0 ! number of variables available in the input file
  INTEGER :: nbvar_tbr    = 0 ! number of variables to be read
  INTEGER :: nbvar_calc   = 0 ! number of variables to be computed from others
  INTEGER :: nbvar_tbw    = 0 ! number of variables to be written
  INTEGER :: nbvar        = 0 ! number of defined variables
  INTEGER :: IINFO_ll         ! return code of // routines
  INTEGER :: nfiles_out   = 0 ! number of output files
  CHARACTER(LEN=:),allocatable :: hvarlist
  TYPE(TFILE_ELT),DIMENSION(1)        :: infiles
  TYPE(TFILE_ELT),DIMENSION(MAXFILES) :: outfiles

  TYPE(workfield), DIMENSION(:), POINTER :: tzreclist

  type(option),dimension(:),allocatable :: options
  character(len=:),allocatable :: hinfile, houtfile
  integer                      :: runmode


  CPROGRAM = 'LFICDF'

  CALL INITIO_ll()
  CALL VERSION
  CALL INI_CST

  ALLOCATE(TIMEZ) !Used by IO_WRITE_FIELD

  NIO_VERB  = NVERB_WARNING
  NGEN_VERB = NVERB_WARNING
  NIO_ABORT_LEVEL = NVERB_FATAL
  NGEN_ABORT_LEVEL = NVERB_FATAL
  LVERB_OUTLST = .FALSE.
  LVERB_STDOUT = .TRUE.

  call read_commandline(options,hinfile,houtfile,runmode)

  if (options(OPTVERBOSE)%set) then
    NIO_VERB  = NVERB_DEBUG
    NGEN_VERB = NVERB_DEBUG
  end if

  IF (options(OPTMERGE)%set) THEN
    NB_PROCIO_R = options(OPTMERGE)%ivalue
  ELSE
    NB_PROCIO_R = 1
  END IF

  IF (runmode == MODELFI2CDF) THEN
     LCDF4    = .TRUE.
     LLFIOUT  = .FALSE.
     LLFIREAD = .TRUE.
     CALL SET_CONFIO_ll()
  ELSE IF (runmode == MODECDF2CDF) THEN
     LCDF4    = .TRUE.
     LLFIOUT  = .FALSE.
     LLFIREAD = .FALSE.
     CALL SET_CONFIO_ll()
  ELSE
     LCDF4    = .TRUE.
     LLFIOUT  = .TRUE.
     LLFIREAD = .FALSE.
     CALL SET_CONFIO_ll()
  END IF

  CALL INI_FIELD_LIST(1)

  CALL OPEN_FILES(infiles, outfiles, nfiles_out, hinfile, houtfile, nbvar_infile, options, runmode)
  IF (options(OPTLIST)%set) STOP

  !Set and initialize parallel variables (necessary to read splitted files)
  CALL SET_JP_ll(1,JPHEXT,JPVEXT,JPHEXT)
  CALL SET_DAD0_ll()
  CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
  CALL SET_XRATIO_ll(1, 1)
  CALL SET_YRATIO_ll(1, 1)
  CALL SET_XOR_ll(1, 1)
  CALL SET_XEND_ll(NIMAX_ll+2*JPHEXT, 1)
  CALL SET_YOR_ll(1, 1)
  CALL SET_YEND_ll(NJMAX_ll+2*JPHEXT, 1)
  CALL INI_PARAZ_ll(IINFO_ll)

  IF (runmode == MODELFI2CDF .OR. runmode == MODECDF2CDF) THEN
     IF (options(OPTVAR)%set) THEN
        ! nbvar_tbr is computed from number of requested variables
        ! by counting commas, = and +
        nbvar_tbr  = 1
        nbvar_calc = 0
        nbvar_tbw = 1
        hvarlist = options(OPTVAR)%cvalue
        DO ji=1,len(hvarlist)
           IF (hvarlist(ji:ji) == ',' .OR.hvarlist(ji:ji) == '+') THEN
              nbvar_tbr = nbvar_tbr+1
           END IF
           IF (hvarlist(ji:ji) == ',') THEN
              nbvar_tbw = nbvar_tbw+1
           END IF
           IF (hvarlist(ji:ji) == '=') THEN
              nbvar_calc = nbvar_calc+1
           END IF
        END DO
        nbvar = nbvar_calc + nbvar_tbr
     ELSE
        nbvar = nbvar_infile
     END IF
  ELSE
    nbvar = nbvar_infile
  END IF

  IF (runmode == MODELFI2CDF) THEN
     ! Conversion LFI -> NetCDF
     IF (options(OPTSPLIT)%set) call open_split_ncfiles_out(outfiles,nfiles_out,houtfile,nbvar_tbw,options)
     CALL parse_infiles(infiles,outfiles,nfiles_out,nbvar_infile,nbvar_tbr,nbvar_calc,nbvar_tbw,tzreclist,options,runmode)
     CALL def_ncdf(infiles,outfiles,nfiles_out)
     CALL fill_files(infiles,outfiles,tzreclist,nbvar,options)

  ELSE IF (runmode == MODECDF2CDF) THEN
     ! Conversion netCDF -> netCDF
     IF (options(OPTSPLIT)%set) call open_split_ncfiles_out(outfiles,nfiles_out,houtfile,nbvar_tbw,options)
     CALL parse_infiles(infiles,outfiles,nfiles_out,nbvar_infile,nbvar_tbr,nbvar_calc,nbvar_tbw,tzreclist,options,runmode)
     CALL def_ncdf(infiles,outfiles,nfiles_out)
     CALL fill_files(infiles,outfiles,tzreclist,nbvar,options)

  ELSE
     ! Conversion NetCDF -> LFI
     CALL parse_infiles(infiles,outfiles,nfiles_out,nbvar_infile,nbvar_tbr,nbvar_calc,nbvar_tbw,tzreclist,options,runmode)
     CALL fill_files(infiles,outfiles,tzreclist,nbvar,options)
  END IF

  CALL CLOSE_FILES(infiles, 1)
  CALL CLOSE_FILES(outfiles,nfiles_out)
  
end program LFI2CDF
