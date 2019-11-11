!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!!    Authors
!!    -------
!
!     D. Gazen
!     Juan 19/08/2005: bug argument optinonel ACCESS --> YACCESS 
!     Juan 22/05/2008: bug mode SPECIFIC in OPEN_ll 
!     Juan 05/11/2009: allow JPMAX_UNIT=48 open files 
!
MODULE MODE_IO_ll

USE MODD_ERRCODES
USE MODD_IO_ll
USE MODE_FD_ll
USE MODD_MPIF

IMPLICIT NONE 

PRIVATE

!INCLUDE 'mpif.h'

INTEGER, PARAMETER :: JPFNULL = 9       !! /dev/null fortran unit
INTEGER, PARAMETER :: JPRESERVED_UNIT   = 11 
INTEGER, PARAMETER :: JPMAX_UNIT        = 48
INTEGER, PARAMETER :: JPMAX_UNIT_NUMBER = JPRESERVED_UNIT+JPMAX_UNIT
! 
CHARACTER(LEN=*),PARAMETER      :: CFILENULL="/dev/null"
!
!! Provisoire
CHARACTER(LEN=*),PARAMETER :: GLOBAL='GLOBAL'
CHARACTER(LEN=*),PARAMETER :: SPECIFIC='SPECIFIC'
!! Provisoire
PUBLIC IONEWFLU,UPCASE,INITIO_ll,OPEN_ll,CLOSE_ll,FLUSH_ll,GLOBAL,SPECIFIC
 
CONTAINS 

FUNCTION IONEWFLU()

INTEGER :: IONEWFLU

INTEGER :: JI
INTEGER :: IOS
LOGICAL :: GEXISTS, GOPENED, GFOUND

GFOUND = .FALSE.

DO JI=JPRESERVED_UNIT, JPMAX_UNIT_NUMBER
  INQUIRE(UNIT=JI, EXIST=GEXISTS, OPENED=GOPENED, IOSTAT=IOS)
  IF (GEXISTS .AND. .NOT. GOPENED .AND. IOS == 0) THEN
    IONEWFLU = JI
    GFOUND   = .TRUE.
    EXIT
  END IF
END DO

IF (.NOT. GFOUND) IONEWFLU = NOSLOTLEFT

END FUNCTION IONEWFLU

FUNCTION UPCASE(HSTRING)
CHARACTER(LEN=*)            :: HSTRING
CHARACTER(LEN=LEN(HSTRING)) :: UPCASE

INTEGER :: JC
INTEGER, PARAMETER :: IAMIN = IACHAR("a")
INTEGER, PARAMETER :: IAMAJ = IACHAR("A")

DO JC=1,LEN(HSTRING)
  IF (HSTRING(JC:JC) >= "a" .AND. HSTRING(JC:JC) <= "z") THEN
    UPCASE(JC:JC) = ACHAR(IACHAR(HSTRING(JC:JC)) - IAMIN + IAMAJ)
  ELSE
    UPCASE(JC:JC) = HSTRING(JC:JC)
  END IF
END DO

END FUNCTION UPCASE


SUBROUTINE INITIO_ll()

INTEGER :: IERR, IOS
LOGICAL :: GISINIT

ISTDERR = 0

CALL MPI_INITIALIZED(GISINIT, IERR)
IF (.NOT. GISINIT) THEN
  CALL MPI_INIT(IERR)
END IF
!! Now MPI is initialized for sure

CALL INITFD()

!! Default number for Processor I/O
ISIOP = 1

!! Get number of allocated processors
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ISNPROC,IERR)
IF (ISNPROC==1) GSMONOPROC = .TRUE.

!! Store proc number
CALL MPI_COMM_RANK(MPI_COMM_WORLD, ISP, IERR)
ISP = ISP + 1

!! Open /dev/null for GLOBAL mode
#if defined(DEV_NULL)
OPEN(UNIT=JPFNULL,FILE=CFILENULL  ,ACTION='WRITE',IOSTAT=IOS)
#else
OPEN(UNIT=JPFNULL,STATUS='SCRATCH',ACTION='WRITE',IOSTAT=IOS)
#endif
IF (IOS > 0) THEN
  WRITE(ISTDERR,*) 'Error OPENING /dev/null...'
  CALL MPI_ABORT(MPI_COMM_WORLD, IOS, IERR)
END IF

!! Init STDOUT and PIPE
IF (ISP == ISIOP) THEN
  ISTDOUT = 6
ELSE
  ISTDOUT = JPFNULL
END IF
    
END SUBROUTINE INITIO_ll
  
SUBROUTINE OPEN_ll(UNIT,    &
                   FILE,    &
                   MODE,    &
                   LFIPAR,  &
                   COMM,    &
                   STATUS,  &
                   ACCESS,  &
                   IOSTAT,  &
                   FORM,    &
                   RECL,    &
                   BLANK,   &
                   POSITION,&
                   ACTION,  &
                   DELIM,   &
                   PAD)
    
INTEGER,         INTENT(OUT)           :: UNIT  !! Different from fortran OPEN
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: FILE
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: MODE
TYPE(LFIPARAM),  POINTER,     OPTIONAL :: LFIPAR
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: STATUS
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: ACCESS
INTEGER,         INTENT(OUT)           :: IOSTAT
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: FORM
INTEGER,         INTENT(IN),  OPTIONAL :: RECL
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: BLANK
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: POSITION
CHARACTER(len=*),INTENT(IN)            :: ACTION
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: DELIM
CHARACTER(len=*),INTENT(IN),  OPTIONAL :: PAD
INTEGER,         INTENT(IN),  OPTIONAL :: COMM

#if defined(MNH_SX5) || defined(MNH_SP4) || defined(NAGf95) || defined(MNH_LINUX)
CHARACTER(len=20)    :: YSTATUS
CHARACTER(len=20)    :: YACCESS
CHARACTER(len=20)    :: YFORM
INTEGER              :: YRECL
INTEGER ,PARAMETER   :: RECL_DEF = 10000
CHARACTER(len=20)    :: YBLANK
CHARACTER(len=20)    :: YPOSITION
CHARACTER(len=20)    :: YDELIM
CHARACTER(len=20)    :: YPAD
!JUAN
#endif
CHARACTER(len=20)    :: YACTION
CHARACTER(len=20)    :: YMODE
INTEGER              :: IOS,IERR
INTEGER              :: ICOMM
INTEGER              :: ICMPRES
TYPE(FD_ll), POINTER :: TZFD, TZFDTEMP
! didier
LOGICAL :: GEXISTS,GOPENED
INTEGER :: IUNIT
! didier
!JUAN SX5 : probleme function retournant un pointer
TYPE(FD_ll), POINTER :: TZJUAN

#ifdef MNH_VPP
!! BUG Fuji avec RECL non fourni en argument de MYOPEN
INTEGER :: IRECSIZE     
IF (PRESENT(RECL)) THEN
  IRECSIZE = RECL
ELSE 
  IRECSIZE = 2147483647  ! Default value for FUJI RECL
END IF
#endif
     
IOS = 0
IF (PRESENT(COMM)) THEN 
  ICOMM = COMM
ELSE
  ICOMM = MPI_COMM_WORLD ! Default communicator
END IF

IF (PRESENT(MODE)) THEN 
  YMODE = UPCASE(TRIM(ADJUSTL(MODE)))
ELSE 
  YMODE = 'GLOBAL'         ! Default Mode
END IF

YACTION = UPCASE(TRIM(ADJUSTL(ACTION)))
IF (YACTION /= "READ" .AND. YACTION /= "WRITE") THEN
  IOSTAT = 99
  UNIT = -1
  WRITE(ISTDERR,*) 'Erreur OPEN_ll : ACTION=',YACTION,' non supportee'
  RETURN
END IF

IF (.NOT. ANY(YMODE == (/'GLOBAL     ','SPECIFIC   ','DISTRIBUTED'/))) THEN
  IOSTAT = 99
  UNIT = -1
  WRITE(ISTDERR,*) 'OPEN_ll error : MODE UNKNOWN'
  RETURN
END IF

!JUAN SX5 : probleme function retournant un pointer
!IF (.NOT. ASSOCIATED(GETFD(FILE))) THEN
TZJUAN=>GETFD(FILE)
IF (.NOT. ASSOCIATED(TZJUAN)) THEN 
!JUAN SX5 : probleme function retournant un pointer
  !! File is not already opened : GOOD
  !! Add a new FD element
  TZFD=>NEWFD()
ELSE 
  !! Error : File already opened
  IOSTAT = 99
  UNIT = -1
  WRITE(ISTDERR,*) 'OPEN_ll error : File', FILE, 'already opened'
  RETURN
END IF

!!$    CALL MPI_ALLREDUCE(ILOCALERR, IGLOBALERR, 1, MPI_INTEGER, MPI_BOR,&
!!$         & ICOMM, IERR)
!!$    IF (IGLOBALERR /= NOERROR) THEN 
!!$       IOSTAT = GLOBALERR
!!$       UNIT = -1
!!$       RETURN 
!!$    END IF



TZFD%NAME = FILE
TZFD%MODE = YMODE
NULLIFY(TZFD%PARAM)

#if defined(MNH_SX5) || defined(MNH_SP4) || defined(NAGf95) || defined(MNH_LINUX)
!JUAN
IF (PRESENT(STATUS)) THEN
YSTATUS=STATUS
ELSE
YSTATUS='UNKNOWN'
ENDIF
IF (PRESENT(ACCESS)) THEN
YACCESS=ACCESS
ELSE
YACCESS='SEQUENTIAL'
ENDIF
IF (PRESENT(FORM)) THEN
YFORM=FORM
ELSE
YFORM='FORMATTED'
ENDIF
IF (PRESENT(RECL)) THEN
YRECL=RECL
ELSE
YRECL=RECL_DEF
ENDIF
IF (PRESENT(BLANK)) THEN
YBLANK=BLANK
ELSE
YBLANK='NULL'
ENDIF
IF (PRESENT(POSITION)) THEN
YPOSITION=POSITION
ELSE
YPOSITION='ASIS'
ENDIF
IF (PRESENT(DELIM)) THEN
YDELIM=DELIM
ELSE
YDELIM='NONE'
ENDIF
IF (PRESENT(PAD)) THEN
YPAD=PAD
ELSE
YPAD='YES'
ENDIF
#endif

SELECT CASE(YMODE)

CASE('GLOBAL')
  IF (YACTION == 'READ') THEN
    TZFD%OWNER = ISP
  ELSE 
    TZFD%OWNER = ISIOP
  END IF
  
  IF (ISP == TZFD%OWNER) THEN 
    !! I/O processor case
    
    TZFD%FLU = IONEWFLU()
#ifdef MNH_VPP
    OPEN(UNIT=TZFD%FLU,       &
         FILE=TRIM(TZFD%NAME),&
         STATUS=STATUS,       &
         ACCESS=ACCESS,       &
         IOSTAT=IOS,          &
         FORM=FORM,           &
         RECL=IRECSIZE,       &
         BLANK=BLANK,         &
         POSITION=POSITION,   &
         ACTION=YACTION,      &
         DELIM=DELIM,         &
         PAD=PAD)
    
#else
#if defined(MNH_SX5) || defined(MNH_SP4) || defined(NAGf95) || defined(MNH_LINUX)
!JUAN : 31/03/2000 modif pour acces direct
IF (YACCESS=='DIRECT') THEN
    OPEN(UNIT=TZFD%FLU,       &
         FILE=TRIM(TZFD%NAME),&
         STATUS=YSTATUS,       &
         ACCESS=YACCESS,       &
         IOSTAT=IOS,          &
         FORM=YFORM,           &
         RECL=YRECL,           &
         ACTION=YACTION)
ELSE
    IF (YFORM=="FORMATTED") THEN
    OPEN(UNIT=TZFD%FLU,       &
         FILE=TRIM(TZFD%NAME),&
         STATUS=YSTATUS,       &
         ACCESS=YACCESS,       &
         IOSTAT=IOS,          &
         FORM=YFORM,           &
         RECL=YRECL,           &
         BLANK=YBLANK,         &
         POSITION=YPOSITION,   &
         ACTION=YACTION,      &
         DELIM=YDELIM,         &
         PAD=YPAD)
    ELSE
    OPEN(UNIT=TZFD%FLU,       &
         FILE=TRIM(TZFD%NAME),&
         STATUS=YSTATUS,       &
         ACCESS=YACCESS,       &
         IOSTAT=IOS,          &
         FORM=YFORM,           &
         RECL=YRECL,           &
         POSITION=YPOSITION,   &
         ACTION=YACTION)
    ENDIF
ENDIF


!print*,' OPEN_ll'
!print*,' OPEN(UNIT=',TZFD%FLU       
!print*,' FILE=',TRIM(TZFD%NAME)
!print*,' STATUS=',YSTATUS       
!print*,' ACCESS=',YACCESS
!print*,' IOSTAT=',IOS
!print*,' FORM=',YFORM
!print*,' RECL=',YRECL
!print*,' BLANK=',YBLANK
!print*,' POSITION=',YPOSITION
!print*,' ACTION=',YACTION
!print*,' DELIM=',YDELIM
!print*,' PAD=',YPAD
#else
    OPEN(UNIT=TZFD%FLU,       &
         FILE=TRIM(TZFD%NAME),&
         STATUS=STATUS,       &
         ACCESS=ACCESS,       &
         IOSTAT=IOS,          &
         FORM=FORM,           &
         RECL=RECL,           &
         BLANK=BLANK,         &
         POSITION=POSITION,   &
         ACTION=YACTION,      &
         DELIM=DELIM,         &
         PAD=PAD)
#endif

#endif

  ELSE 
    !! NON I/O processors case
    IOS = 0
    TZFD%FLU = JPFNULL 
  END IF
  
CASE('SPECIFIC')
  TZFD%OWNER = ISP
  TZFD%FLU = IONEWFLU()
  
#ifdef MNH_VPP
  OPEN(UNIT=TZFD%FLU,                      &
       FILE=TRIM(TZFD%NAME)//SUFFIX(".P"), &
       STATUS=STATUS,                         &
       ACCESS=ACCESS,                         &
       IOSTAT=IOS,                            &
       FORM=FORM,                             &
       RECL=IRECSIZE,                         &
       BLANK=BLANK,                           &
       POSITION=POSITION,                     &
       ACTION=YACTION,                        &
       DELIM=DELIM,                           &
       PAD=PAD)
  
#else
#if defined(MNH_SX5) || defined(MNH_SP4) || defined(NAGf95) || defined(MNH_LINUX)
IF (ACCESS=='DIRECT') THEN
    OPEN(UNIT=TZFD%FLU,       &
         FILE=TRIM(TZFD%NAME)//SUFFIX(".P"), &
         STATUS=YSTATUS,       &
         ACCESS=YACCESS,       &
         IOSTAT=IOS,          &
         FORM=YFORM,           &
         RECL=YRECL,           &
         ACTION=YACTION)
ELSE
  OPEN(UNIT=TZFD%FLU,                      &
       FILE=TRIM(TZFD%NAME)//SUFFIX(".P"), &
       STATUS=YSTATUS,                         &
       ACCESS=YACCESS,                         &
       IOSTAT=IOS,                             &
       FORM=YFORM,                             &
       RECL=YRECL,                             &
       BLANK=YBLANK,                           &
       POSITION=YPOSITION,                     &
       ACTION=YACTION,                         &
       DELIM=YDELIM,                           &
       PAD=YPAD)
ENDIF
#else
  OPEN(UNIT=TZFD%FLU,                      &
       FILE=TRIM(TZFD%NAME)//SUFFIX(".P"), &
       STATUS=STATUS,                         &
       ACCESS=ACCESS,                         &
       IOSTAT=IOS,                            &
       FORM=FORM,                             &
       RECL=RECL,                             &
       BLANK=BLANK,                           &
       POSITION=POSITION,                     &
       ACTION=YACTION,                        &
       DELIM=DELIM,                           &
       PAD=PAD)
#endif
  
#endif

CASE('DISTRIBUTED')
  TZFD%OWNER = ISIOP
  IF (.NOT. PRESENT(LFIPAR)) THEN
    PRINT *,"ERROR OPEN_ll : LFI non present"
    RETURN
  END IF
  TZFD%PARAM=>LFIPAR
  
  IF (ISP == TZFD%OWNER) THEN 
    TZFD%FLU = IONEWFLU()
  ELSE 
    !! NON I/O processors case
    IOS = 0
    TZFD%FLU = -1
  END IF
  
  
  
END SELECT

! Recherche d'un communicateur a reutiliser
! TZFD is the first element

TZFD%COMM = MPI_COMM_NULL

TZFDTEMP=>TZFD%NEXT
DO WHILE(ASSOCIATED(TZFDTEMP))
  CALL MPI_COMM_COMPARE(ICOMM,TZFDTEMP%COMM,ICMPRES,IERR)
  IF (ICMPRES == MPI_CONGRUENT) THEN
    TZFD%COMM = TZFDTEMP%COMM
    EXIT
  END IF
  TZFDTEMP=>TZFDTEMP%NEXT
END DO

IF (TZFD%COMM == MPI_COMM_NULL) THEN
  ! Pas de communicateur equivalent -> on duplique
  !
  CALL MPI_COMM_DUP(ICOMM, TZFD%COMM, IERR)
  !       WRITE(ISTDOUT,*) 'FILE = ',TZFD%NAME,', comm ',TZFD%COMM&
  !            & , ' cree par duplication de comm ', ICOMM
END IF

IOSTAT = IOS
UNIT = TZFD%FLU 

CONTAINS
FUNCTION SUFFIX(HEXT)

CHARACTER(len=*)             :: HEXT
CHARACTER(len=LEN(HEXT)+3)   :: SUFFIX

WRITE(SUFFIX,'(A,i3.3)') TRIM(HEXT), ISP

END FUNCTION SUFFIX

END SUBROUTINE OPEN_ll
  
SUBROUTINE CLOSE_ll(HFILE,IOSTAT,STATUS)

CHARACTER(LEN=*), INTENT(IN)            :: HFILE
INTEGER,          INTENT(OUT), OPTIONAL :: IOSTAT
CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: STATUS

TYPE(FD_ll), POINTER :: TZFD
INTEGER :: OLDCOMM

INTEGER :: IERR, IGLOBALERR, IRESP

CHARACTER(LEN=100)                      :: STATUSL

TZFD=>GETFD(HFILE)

IF (.NOT. ASSOCIATED(TZFD)) THEN
  WRITE(ISTDOUT,*) 'Erreur CLOSE_ll : Fichier : ', HFILE, ' non&
       & present...'
  IF (PRESENT(IOSTAT)) IOSTAT = BADVALUE
  RETURN
END IF

IRESP      = 0
IGLOBALERR = 0
IF (PRESENT(STATUS))  THEN
STATUSL = STATUS
ELSE
STATUSL = "KEEP"
ENDIF

SELECT CASE(TZFD%MODE)
CASE('GLOBAL','SPECIFIC')
  IF (TZFD%OWNER == ISP) THEN
    CLOSE(UNIT=TZFD%FLU, IOSTAT=IRESP,STATUS=STATUSL)
  END IF
  CALL MPI_ALLREDUCE(IRESP,IGLOBALERR,1,MPI_INTEGER,MPI_BOR,TZFD&
       & %COMM,IERR)
CASE('DISTRIBUTED')
  ! nothing to close with FM-Files
END SELECT

OLDCOMM = TZFD%COMM   !! Recopie dans var. temporaire

CALL DELFD(TZFD)

IF (IRESP == IGLOBALERR) THEN
  
  ! liberation du communicateur
  !
  TZFD=>GETFD(OLDCOMM)
  
  IF (.NOT. ASSOCIATED(TZFD)) THEN
    CALL MPI_COMM_FREE(OLDCOMM, IERR)
  END IF
END IF

IF (PRESENT(IOSTAT)) IOSTAT = IGLOBALERR
 
END SUBROUTINE CLOSE_ll

SUBROUTINE FLUSH_ll(HFILE,IRESP)
#if defined(NAGf95)
USE F90_UNIX
#endif
CHARACTER(LEN=*), INTENT(IN)            :: HFILE
INTEGER,          INTENT(OUT), OPTIONAL :: IRESP

TYPE(FD_ll), POINTER :: TZFD
INTEGER              :: IUNIT

IRESP=0
TZFD=>GETFD(HFILE)
IF (.NOT. ASSOCIATED(TZFD)) THEN
  WRITE(ISTDOUT,*) 'Error in FLUSH_ll : file ',TRIM(HFILE),&
       &' not present !'
  IF (PRESENT(IRESP)) IRESP = BADVALUE
  RETURN
END IF

IUNIT=TZFD%FLU
IF (TZFD%OWNER == ISP .AND. TZFD%MODE /= 'DISTRIBUTED') THEN
#if defined(MNH_SP4)
  CALL FLUSH(IUNIT)
#else
  CALL FLUSH(IUNIT)
#endif
END IF

END SUBROUTINE FLUSH_ll

END MODULE MODE_IO_ll
