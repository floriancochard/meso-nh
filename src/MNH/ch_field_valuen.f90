!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ############################ 
      MODULE MODI_CH_FIELD_VALUE_n
!!    ############################ 
!!
INTERFACE
!!
FUNCTION CH_FIELD_VALUE_n(PZ, HGRIDTYPE, HNAME, &
 			  HUNIT, KLUOUT, KVERB)
IMPLICIT NONE
REAL                         :: CH_FIELD_VALUE_n! the function name
REAL,             INTENT(IN) :: PZ              ! x-y-z coo. to initialize
CHARACTER(LEN=*), INTENT(IN) :: HGRIDTYPE       ! grid type passed as x-y-z coo.
CHARACTER(LEN=*), INTENT(IN) :: HNAME           ! name of species to initialize
CHARACTER(LEN=*), INTENT(OUT):: HUNIT           ! unit ("CON" or "MIX")
INTEGER,          INTENT(IN) :: KLUOUT          ! output listing channel
INTEGER,          INTENT(IN) :: KVERB           ! verbosity level
!!
END FUNCTION CH_FIELD_VALUE_n
!!
END INTERFACE
!!
END MODULE MODI_CH_FIELD_VALUE_n
!!
!!    ########################################################## 
      FUNCTION CH_FIELD_VALUE_n(PZ, HGRIDTYPE, HNAME, &
				HUNIT, KLUOUT, KVERB)
!!    ##########################################################
!!
!!*** *CH_FIELD_VALUE_n*
!!
!!    PURPOSE
!!    -------
!        initialize a given species HNAME at a given point (PZ)
!!
!!**  METHOD
!!    ------
!!       This subroutine may be adapted to the users needs.
!!    The current version initializes horizontally homogeneous fields
!!    From the general purpose ACSII input file, several form-profiles
!!    will be read. Then each species will be associated to one of those
!!    given profile.
!!    Finally one norm-initial value will be read for each species.
!!    The initial value is calculated by linear interpolation of the 
!!    associated form-profile onto the requested PZ coordinate and
!!    multiplication with the given norm-initial value for the requested 
!!    species.
!!       The user may modify this subroutine for more complicated
!!    initializations by passing additional information by the general
!!    purpose input file. The namelist variable CCH_INIT_FIELD_OPT may be
!!    used as a switch between different options. Presently this parameter
!!    is not used.
!!
!!    REFERENCE
!!    ---------
!!    MesoNH-chemistry book 1, 2, 3
!!
!!    AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 03/11/95
!!    05/08/95 (K. Suhre) restructered
!!    11/08/98 (N. Asencio) add parallel code
!!    28/07/99 (V. Crassier & K. Suhre) modify initialization scheme (1-D)
!!    Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_OPEN_INPUT  ! open general purpose ASCII input file
USE MODD_IO_ll,         ONLY: TFILEDATA
USE MODE_FM,            ONLY: IO_FILE_CLOSE_ll
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MNHC_n, ONLY: CCHEM_INPUT_FILE     ! general purpose input file
USE MODD_SUB_CH_FIELD_VALUE_n
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL                         :: CH_FIELD_VALUE_n! the function name
REAL,             INTENT(IN) :: PZ              ! x-y-z coo. to initialize
CHARACTER(LEN=*), INTENT(IN) :: HGRIDTYPE       ! grid type passed as x-y-z coo.
CHARACTER(LEN=*), INTENT(IN) :: HNAME           ! name of species to initialize
CHARACTER(LEN=*), INTENT(OUT):: HUNIT           ! unit ("CON" or "MIX")
INTEGER,          INTENT(IN) :: KLUOUT          ! output listing channel
INTEGER,          INTENT(IN) :: KVERB           ! verbosity level
!
!*      0.2    declarations local variables
!
INTEGER       :: JI, JJ              ! loop control variables
INTEGER       :: ICHANNEL            ! I/O channel for file input
CHARACTER(LEN=40) :: YFORMAT         ! format for input
INTEGER       :: IFAIL               ! return code from CLOSE
!
!
!
INTEGER :: IASSOACT  ! actual index of associated profile
INTEGER :: IINITACT  ! actual index of initial value
REAL    :: IINF      ! inferieur Z-level index use in monotony check for Z-profile
REAL    :: ZINF      ! inferieur Z-level value use in monotony check for Z-profile
INTEGER :: IZINDEX   ! lower index for interpolation
!
REAL    :: ZRETURN   ! return value for function
!
TYPE(TFILEDATA),POINTER :: TZFILE
!-------------------------------------------------------------------------------
TZFILE => NULL()
!
!*       1.    READ DATA FROM GENERAL PURPOSE INPUT FILE
!              -----------------------------------------
!
! on first call, data has to be read from CCHEM_INPUT_FILE
firstcall: IF (GSFIRSTCALL) THEN
  GSFIRSTCALL = .FALSE.
!
  ! print options
  IF (KVERB >= 5) THEN
    WRITE(KLUOUT,*) "CH_FIELD_VALUE_n: grid type is ", HGRIDTYPE
  END IF
!
!
!*       1.1   read form-profiles
!
! open input file
  IF (KVERB >= 5) WRITE(KLUOUT,*) "CH_FIELD_VALUE_n: reading form-profiles"
  CALL CH_OPEN_INPUT(CCHEM_INPUT_FILE, "FORMPROF", TZFILE, KLUOUT, KVERB)
  ICHANNEL = TZFILE%NLU
!
! read number of profiles ISPROF and number of levels ISLEVEL
  READ(ICHANNEL,*) ISPROF, ISLEVEL
  IF (KVERB >= 5) THEN
    WRITE(KLUOUT,*) ISPROF, " profiles and ", ISLEVEL, " levels will be read"
  END IF
!
! read data input format
  READ(ICHANNEL,"(A)") YFORMAT
  IF (KVERB >= 5) THEN
    WRITE(KLUOUT,*) "input format is: ", YFORMAT
  END IF
!
! allocate fields for data input
  ALLOCATE(ZSNORMPROF(ISLEVEL,ISPROF))
  ALLOCATE(ZSZPROF(ISLEVEL))
!
! read data
  DO JI = 1, ISLEVEL
    IF (KVERB >= 5) WRITE(KLUOUT,*) "reading data for level ", JI
    READ(ICHANNEL,FMT=YFORMAT) ZSZPROF(JI), (ZSNORMPROF(JI,JJ),JJ=1,ISPROF)
  ENDDO
  IF (KVERB >= 10) THEN
    WRITE(KLUOUT,*) "The form profiles that have been read are as follows:"
    DO JI = 1, ISLEVEL
      WRITE(KLUOUT,YFORMAT) ZSZPROF(JI), (ZSNORMPROF(JI,JJ),JJ=1,ISPROF)
      WRITE(KLUOUT,*) " level: ",JI," ZSZPROF=",ZSZPROF(JI)
    END DO
  END IF
!
! close file
  CALL IO_FILE_CLOSE_ll(TZFILE)
  TZFILE => NULL()
!
! check if Z-profile is given in increasing order, otherwise stop
  ZINF = ZSZPROF(1) ; IINF=1
  DO JI = 2, ISLEVEL
    IF (ZINF .GE. ZSZPROF(JI)) THEN
      WRITE(KLUOUT,*) &
   	   "CH_FIELD_VALUE_n-Error: Z-profile must be in increasing order!"
      WRITE(KLUOUT,*) " minimum value: ",ZINF," at level ",IINF
      WRITE(KLUOUT,*) " current value: ",ZSZPROF(JI)," at level ",JI
      ! callabortstop
      CALL ABORT
      STOP "Program stopped by CH_FIELD_VALUE_n"
    ENDIF
    ZINF = ZSZPROF(JI) ; IINF=JI
  ENDDO
!
!*       1.2   read species <--> form-profile association
!
! open input file
  IF (KVERB >= 5) WRITE(KLUOUT,*) &
     "CH_FIELD_VALUE_n: reading species <--> form-profile association"
  CALL CH_OPEN_INPUT(CCHEM_INPUT_FILE, "PROFASSO", TZFILE, KLUOUT, KVERB)
  ICHANNEL = TZFILE%NLU
!
! read number of associations ISASSO
  READ(ICHANNEL, *) ISASSO
  IF (KVERB >= 5) WRITE(KLUOUT,*) "number of associations: ", ISASSO
!
! read data input format
  READ(ICHANNEL,"(A)") YFORMAT
  IF (KVERB >= 5) WRITE(KLUOUT,*) "input format is: ", YFORMAT
!
! allocate fields
  ALLOCATE(YSASSONAME(ISASSO))
  ALLOCATE(ISASSOPROF(ISASSO))
!
! read associations
  DO JI = 1, ISASSO
    READ(ICHANNEL,FMT=YFORMAT) YSASSONAME(JI), ISASSOPROF(JI)
    IF (KVERB >= 5) WRITE(KLUOUT,FMT=YFORMAT) YSASSONAME(JI), ISASSOPROF(JI)
  ENDDO
!
! close file
  CALL IO_FILE_CLOSE_ll(TZFILE)
  TZFILE => NULL()
!
!
!*       1.3   read norm-initial values
!
! open input file
  IF (KVERB >= 5) WRITE(KLUOUT,*) "CH_FIELD_VALUE_n:reading norm-initial values"
  CALL CH_OPEN_INPUT(CCHEM_INPUT_FILE, "NORMINIT", TZFILE, KLUOUT, KVERB)
  ICHANNEL = TZFILE%NLU
!
! read units for initial data (may be "CON" for molec./cm3 or "MIX" for ppp)
  READ(ICHANNEL,"(A)") HUNIT
  IF (HUNIT .EQ. "CON") THEN
    IF (KVERB >= 5) THEN
      WRITE(KLUOUT,*) "initial data is given as number density (molec./cm3)"
    END IF
  ELSE IF (HUNIT .EQ. "MIX") THEN
    IF (KVERB >= 5) THEN
      WRITE(KLUOUT,*) "initial data is given as mixing ratio (part per par)"
    END IF
  ELSE
    WRITE(KLUOUT,*) "CH_FIELD_VALUE_n ERROR: unit type unknown: ", HUNIT
    ! callabortstop
    CALL ABORT
    STOP
  ENDIF
!
! read number of initial values ISINIT
  READ(ICHANNEL, *) ISINIT
  IF (KVERB >= 5) WRITE(KLUOUT,*) "number of initial values: ", ISINIT
!
! read data input format
  READ(ICHANNEL,"(A)") YFORMAT
  IF (KVERB >= 5) WRITE(KLUOUT,*) "input format is: ", YFORMAT
!
! allocate fields
  ALLOCATE(YSINITNAME(ISINIT))
  ALLOCATE(ZSINITVAL(ISINIT))
!
! read associations
  DO JI = 1, ISINIT
    READ(ICHANNEL,FMT=YFORMAT) YSINITNAME(JI), ZSINITVAL(JI)
    IF (KVERB >= 5) WRITE(KLUOUT,FMT=YFORMAT) YSINITNAME(JI), ZSINITVAL(JI)
  ENDDO
!
! close file
  CALL IO_FILE_CLOSE_ll(TZFILE)
  TZFILE => NULL()
!
ENDIF firstcall
!
!-------------------------------------------------------------------------------
!
!*       2     INTERPOLATE INITIAL VALUE ON Z 
!              -------------------------------
!
!*       2.1   find associated profile number and initial value for HNAME
!
IASSOACT = 0
find_asso : DO JI = 1, ISASSO
  IF (HNAME .EQ. YSASSONAME(JI)) THEN
    IASSOACT = ISASSOPROF(JI)
    EXIT find_asso
  ENDIF
ENDDO find_asso
!
IINITACT = 0
find_init : DO JI = 1, ISINIT
  IF (HNAME .EQ. YSINITNAME(JI)) THEN
    IINITACT = JI
    EXIT find_init
  ENDIF
ENDDO find_init
!
!*       2.2   if no profile is associated return the norm-initial value 
!
IF ((IASSOACT .EQ. 0) .AND. (IINITACT .NE. 0)) THEN
  CH_FIELD_VALUE_n = ZSINITVAL(IINITACT)
  IF (KVERB >= 15) THEN
    WRITE(KLUOUT,'(A30,F10.2,A,E15.6,A)') &
      TRIM(HNAME) // "(", PZ, ") = ", ZSINITVAL(IINITACT), &
      " (no profile / norm-init)"
  END IF
  RETURN
ENDIF
!
!*       2.3   if no initial value is given return zero
!
IF (IINITACT .EQ. 0) THEN
  CH_FIELD_VALUE_n = 0.0
  IF (KVERB >= 15) THEN
    WRITE(KLUOUT,'(A30,F10.2,A,E15.6,A)') &
      TRIM(HNAME) // "(", PZ, ") = ", 0.0, &
      " (no initial value)"
  END IF
  RETURN
ENDIF
!
!*       2.4   if norm-initial value and a form-profile are associated,
!              interpolate
!
IZINDEX = 1
search_loop : DO JI = ISLEVEL-1, 1, -1
  IF (PZ .GE. ZSZPROF(JI)) THEN
    IZINDEX = JI
    EXIT search_loop
  ENDIF
ENDDO search_loop
!
!*       2.5   check boundaries of IASSOACT and IINITACT
!
IF ((IASSOACT.LE.0).OR.(IASSOACT.GT.ISPROF)) THEN
  WRITE(KLUOUT,*) &
    "CH_FIELD_VALUE_n-ERROR: unproper associated profile value:", IASSOACT
    ! callabortstop
    CALL ABORT
  STOP
ENDIF
!
IF ((IINITACT.LE.0).OR.(IINITACT.GT.ISINIT)) THEN
  WRITE(KLUOUT,*) &
    "CH_FIELD_VALUE_n-ERROR: unproper associated initial value:", IINITACT
    ! callabortstop
    CALL ABORT
  STOP
ENDIF
!
!*       2.6   linear interpolation between IZINDEX and IZINDEX+1, 
!              but return zero if extrapolation generates negative values
!
ZRETURN = MAX(0., ZSINITVAL(IINITACT)                                        &
   *(                                                                        &
      (PZ                 - ZSZPROF(IZINDEX))*ZSNORMPROF(IZINDEX+1,IASSOACT) &
    + (ZSZPROF(IZINDEX+1) - PZ              )*ZSNORMPROF(IZINDEX  ,IASSOACT) &
    )/(ZSZPROF(IZINDEX+1) - ZSZPROF(IZINDEX)))
!
IF (KVERB >= 15) THEN
  WRITE(KLUOUT,'(A30,F10.2,A,E15.6,A)') &
    TRIM(HNAME) // "(", PZ, ") = ", ZRETURN, &
    " (based on profile and norm-init)"
END IF
!
CH_FIELD_VALUE_n = ZRETURN
!
RETURN
!
!-------------------------------------------------------------------------------
!
END FUNCTION CH_FIELD_VALUE_n
