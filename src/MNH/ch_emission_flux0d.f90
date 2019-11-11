!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ############################
      MODULE MODI_CH_EMISSION_FLUX0D
!!    ############################
!!    
INTERFACE
SUBROUTINE CH_EMISSION_FLUX0D(PTIME, PFLUX, HINPUTFILE, KLUOUT, KVERB)
USE MODD_CH_M9_n,      ONLY: NEQ
IMPLICIT NONE
REAL,                 INTENT(IN)  :: PTIME      ! time of simulation in sec UTC
                                                ! (counting from midnight)
REAL, DIMENSION(NEQ), INTENT(OUT) :: PFLUX      ! emission flux in ppp*m/s
CHARACTER*(*),        INTENT(IN)  :: HINPUTFILE ! name of the input file
INTEGER,              INTENT(IN)  :: KLUOUT     ! output listing channel
INTEGER,              INTENT(IN)  :: KVERB      ! verbosity level
END SUBROUTINE CH_EMISSION_FLUX0D
END INTERFACE
!!
END MODULE MODI_CH_EMISSION_FLUX0D
!!
!!    #################################################################### 
      SUBROUTINE CH_EMISSION_FLUX0D(PTIME, PFLUX, HINPUTFILE, KLUOUT, KVERB)
!!    ####################################################################
!!
!!***  *CH_EMISSION_FLUX0D*
!!
!!    PURPOSE
!!    -------
!!      Return a time-dependant emission flux based on tabulated values
!!
!!**  METHOD
!!    ------
!!      From HINPUTFILE a set of emission fluxes for an arbitrary number of
!!    variables is read. The fluxes are interpolated in time. Data is read
!!    after the keyword EMISDATA in the following format:
!!
!!    -------------------------------------------------
!!    comment line
!!    unit identifier [MIX|CON|MOL]
!!    number of species (N)
!!    number of records (M)
!!    name of species 1
!!     ..
!!    name of species N
!!    input format for the data
!!    time1 emission(spec1,t1) .. emission(specN,t1)
!!     ..    ..                    ..
!!    timeM emission(spec1,tM) .. emission(specN,tM)
!!    -------------------------------------------------
!!
!!    where the unit identifier [MIX|CON|MOL] indicates whether
!!    the flux is given as 
!!    CON: molecules/cm2/s 
!!    MIX: ppp*m/s
!!    MOL: microMol/m2/day
!!    (assuming standard pressure and temperature in the conversion)
!!    The returned flux is given in ppp*m/s, that is standard MesoNH
!!    units so that no conversion is to be applied when introducing 
!!    the emission flux in the 3-D model.
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 26/07/1999
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    EXTERNAL
!!    --------
USE MODD_IO_ll,         ONLY: TFILEDATA
USE MODE_FM,            ONLY: IO_FILE_CLOSE_ll
USE MODE_IO_ll
!
USE MODI_CH_OPEN_INPUT
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_n,      ONLY: NEQ, CNAMES
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!

REAL,                 INTENT(IN)  :: PTIME      ! time of simulation in sec UTC
                                                ! (counting from midnight)
REAL, DIMENSION(NEQ), INTENT(OUT) :: PFLUX      ! emission flux in ppp*m/s
CHARACTER*(*),        INTENT(IN)  :: HINPUTFILE ! name of the input file
INTEGER,              INTENT(IN)  :: KLUOUT     ! output listing channel
INTEGER,              INTENT(IN)  :: KVERB      ! verbosity level
!
!*       0.2  declaration of local variables
!
INTEGER       :: JI, JJ      ! loop control
CHARACTER*80  :: YCOMMENT    ! comment line in the input file
CHARACTER*80  :: YFORMAT     ! format of the input data
INTEGER       :: ICHEMIS     ! number of variables for which a flux is given
                             ! in the input file
INTEGER       :: IIO         ! I/O channel
INTEGER       :: IFAIL       ! return code from CLOSE_ll
REAL          :: ZALPHA      ! interpolation weight
!
CHARACTER(LEN=3)                              :: YUNIT
                             ! unit of the flux
REAL                                          :: ZCONVERSION
                             ! conversion factor to convert from MOL/CON to MIX
CHARACTER(LEN=32), ALLOCATABLE, DIMENSION(:)  :: YNAMEIN 
                             ! name of species in file
REAL,             ALLOCATABLE, DIMENSION(:,:) :: ZFLUXIN 
                             ! fluxes as given in file:
                             ! first index is the species,second index the time
TYPE(TFILEDATA),POINTER :: TZFILE
!
!*       0.3  declaration of saved local variables
!
INTEGER, SAVE                              :: ISTIMES      
                                           ! number of time steps in input file
REAL,    SAVE, ALLOCATABLE, DIMENSION(:)   :: ZSFTIME 
                                           ! time steps for which 
                                           ! fluxes are defined
REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:) :: ZSFLUX  
                                           ! fluxes as fct. of species
                                           ! first index is the species
                                           ! second index is the time
INTEGER, SAVE                              :: ISACT = 1 
                                           ! actual record (for interpolation)
LOGICAL, SAVE                              :: LSFIRSTCALL = .TRUE.  
                                           ! flag to identify first call
!
!------------------------------------------------------------------------------
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
TZFILE => NULL()
!
IF (LSFIRSTCALL) THEN
!
!*       1.   READ DATA 
!        --------------
!
  CALL CH_OPEN_INPUT(HINPUTFILE, "EMISDATA", TZFILE, KLUOUT, KVERB)
  IIO = TZFILE%NLU
!
! read unit identifier
  READ(IIO,'(A3)') YUNIT
!
! read number of variables 
  READ(IIO,*) ICHEMIS
!
! read number of records
  READ(IIO,*) ISTIMES
!
! allocate arrays
!
  ALLOCATE(YNAMEIN(ICHEMIS))
  ALLOCATE(ZSFTIME(ISTIMES))
  ALLOCATE(ZFLUXIN(ICHEMIS,ISTIMES))
  ALLOCATE(ZSFLUX(NEQ,ISTIMES))
!
! read variable names for which fluxes are given
  DO JI = 1, ICHEMIS
    READ(IIO,*) YNAMEIN(JI)
  END DO
!
! read the input format
  READ(IIO,'(A)') YFORMAT
!
! print names of input variables
  IF (KVERB >= 5) THEN
    WRITE(KLUOUT,*) 'CH_EMISSION_FLUX0D: the following ', ICHEMIS, &
	     ' emission fluxes will be read'
    DO JI = 1, ICHEMIS
      WRITE(KLUOUT,*) YNAMEIN(JI)
    END DO
    WRITE(KLUOUT,*) ISTIMES, ' record(s) will be read now ...'
  END IF
!
  DO JJ = 1, ISTIMES
    READ(IIO,YFORMAT) ZSFTIME(JJ), (ZFLUXIN(JI,JJ), JI = 1, ICHEMIS)
  END DO
!
  IF (KVERB >= 10) THEN
    WRITE(KLUOUT,*) 'CH_EMISSION_FLUX0D: the following fluxes were read:'
    DO JJ = 1, ISTIMES
      WRITE(KLUOUT,YFORMAT) ZSFTIME(JJ), (ZFLUXIN(JI,JJ), JI = 1, ICHEMIS)
    END DO
  END IF
!
! close file
!
  CALL IO_FILE_CLOSE_ll(TZFILE)
!
!*       2.   MAP DATA ONTO PROGNOSTIC VARIABLES
!        ---------------------------------------
!
! determine the conversion factor
  SELECT CASE (YUNIT)
  CASE ('MIX') ! flux given ppp*m/s, no conversion required
    ZCONVERSION = 1.0
  CASE ('CON') ! flux given in molecules/cm2/s
               ! where 1 molecule/cm2/s = (224.14/6.022136E23) ppp*m/s
    ZCONVERSION = (224.14/6.022136E23)
  CASE ('MOL') ! flux given in microMol/m2/day
               ! where 1 microMol/m2/day = (22.414/86.400)*1E-12 ppp*m/s
    ZCONVERSION = (22.414/86.400)*1E-12
  CASE DEFAULT
    WRITE(KLUOUT,*) 'CH_EMISSION_FLUX0D: unknow conversion factor: ', YUNIT
!callabortstop
    CALL ABORT
    STOP 'CH_EMISSION_FLUX0D: unknow conversion factor' 
  END SELECT
!
! set all fluxes to zero
  ZSFLUX(:,:) = 0.0
!
! loop over all species in the file and match with those in the
! reaction system
  DO JI = 1, ICHEMIS
    inner: DO JJ = 1, NEQ
      IF (YNAMEIN(JI) .EQ. CNAMES(JJ)) THEN
        ZSFLUX(JJ,:) = ZFLUXIN(JI,:)
        IF (KVERB >= 5) &
          WRITE(KLUOUT,*) 'emission fluxes initialized for ', CNAMES(JJ)
      END IF
    END DO inner
  END DO
!
! convert units
  ZSFLUX(:,:) = ZSFLUX(:,:) * ZCONVERSION
!
! free unused memory
  DEALLOCATE(YNAMEIN)
  DEALLOCATE(ZFLUXIN)
!
  LSFIRSTCALL = .FALSE.
  ISACT = 1
!
END IF ! firstcall
!
!
!------------------------------------------------------------------------------
!
!*    3.  INTERPOLATE SURFACE FLUXES IN TIME
!     ---------------------------------------
!
IF (PTIME .LE. ZSFTIME(1)) THEN
!
! take first record
  PFLUX(:) = ZSFLUX(:,1) 
!
ELSE IF (PTIME .GE. ZSFTIME(ISTIMES)) THEN
!
! take last record
  PFLUX(:) = ZSFLUX(:,ISTIMES) 
!
ELSE
!
! interpolate meteo variables in time
!
  DO WHILE (PTIME .GE. ZSFTIME(ISACT+1))
    ISACT = ISACT+1
  END DO
!
  ZALPHA = (PTIME            - ZSFTIME(ISACT)) &
	 / (ZSFTIME(ISACT+1) - ZSFTIME(ISACT))
!
  PFLUX(:) = ZALPHA      * ZSFLUX(:,ISACT+1) &
  	   + (1.-ZALPHA) * ZSFLUX(:,ISACT)
!
END IF
!
END SUBROUTINE CH_EMISSION_FLUX0D
