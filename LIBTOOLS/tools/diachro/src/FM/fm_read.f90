!     ######spl
      SUBROUTINE FM_READ(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                        KLENCH,HCOMMENT,KRESP)
!     ###########################################################
!
!!****  *FM_READ* - routine to read a single data article in a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREAD is to read one single article of data in
!     a Meso-nh file. This routine only holds for LFI-files (not namelists)
!
!!**  METHOD
!!    ------
!!
!!      The unformatted fortran read operation is actually executed in the
!!    routine LFILEC. You just need to indicate the name of the file
!!    without the ".lfi" suffix,
!!    and the name of the article you want to read, as well as the length of
!!    the field. LFILEC then knows how
!!    to get the record number of the desired field by referring to an intern
!!    table of association.
!!      In FMREAD, the data is first stored in IWORK and then split in KGRID
!!    (IWORK(1)=C-grid indicator) and KFIELD (integer or real data field)
!!    which are both stored on the same LFI logical article.
!!
!!    EXTERNAL
!!    --------
!!
!!      FMLOOK,LFINFO,LFILEC,CHAR
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      MODULE: MODD_FMDECLAR contains management parameters and
!!              storage arrays to move information around at the
!!              level of all "FM"-routines.
!!
!!    REFERENCE
!!    ---------
!!
!!      see the Technical Specifications Report for the Meso-nh project
!!      (in French)
!!
!!    AUTHOR
!!    ------
!!
!!      C. FISCHER      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                        06/94
!!      modified by V. Masson               16/09/96 (prints if error occurs)
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR

IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM ! file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM ! name of the desired article

CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field

INTEGER(KIND=8),DIMENSION(1:KLENG),INTENT(OUT)::KFIELD ! array containing 
                                                        ! the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=JPXKRK),     INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems occured

!
!*      0.2   Declarations of local variables
!
INTEGER::IRESP,ILENGA,IPOSEX,ITOTAL,INUMBR,J,IROW,IFMFNL,ILUPRI
INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE::IWORK,IWORKNEW
INTEGER,DIMENSION(1:JPXKRK)::ICOMMENT
CHARACTER(LEN=JPFINL)::YFNLFI
CHARACTER(LEN=LEN(HFILEM))::YINTFN
INTEGER :: DATASIZE,ITYPCOD,NEWSIZE
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKREAD/ILUPRI,INUMBR,IRESP
!DIR$ TASKCOMMON TASKREAD
!
!----------------------------------------------------------------------------
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0 ; IROW = 0 ; ILUPRI = 6
IFMFNL=JPFINL-4

IROW=LEN(HFILEM)

IF (IROW.EQ.0) THEN
   IRESP=-61
   GOTO 1000
ELSEIF (IROW.GT.IFMFNL) THEN
   IRESP=-62
   GOTO 1000
ENDIF
YINTFN=ADJUSTR(HFILEM)
YFNLFI=YINTFN//'.lfi'
YFNLFI=ADJUSTL(YFNLFI)

!
!*      1.2   WE LOOK FOR THE FILE'S LOGICAL UNIT
!
CALL FMLOOK(YFNLFI,HFIPRI,INUMBR,IRESP)
IF (IRESP.NE.0) GOTO 1000

!
!*      2.a   LET'S GET SOME INFORMATION ON THE DESIRED ARTICLE
!
!ILENGA=0
!print *,' ***FM_READ ILENGA mis a 0 avant CALL LFINFO'
CALL LFINFO(IRESP,INUMBR,HRECFM,ILENGA,IPOSEX)
!print *,' ***FM_READ ILENGA,IRESP AP LFINFO ',ILENGA,IRESP
IF (IRESP.NE.0) THEN
        GOTO 1000
ELSEIF (ILENGA.EQ.0) THEN
!print *,' ***FM_READ passage IRESP=-47 GOTO 1000'
        IRESP=-47
        GOTO 1000
ELSEIF (ILENGA.GT.JPXFIE) THEN
        IRESP=-48
        GOTO 1000
ENDIF

!
!*      2.b   UNFORMATTED DIRECT ACCESS READ OPERATION
!
ITOTAL=ILENGA
IF(ALLOCATED(IWORK)) DEALLOCATE(IWORK)
ALLOCATE(IWORK(ITOTAL))

CALL LFILEC(IRESP,INUMBR,HRECFM,IWORK,ITOTAL)
IF (IRESP.NE.0) GOTO 1000
!
!*      2.c   THE GRID INDICATOR AND THE COMMENT STRING
!*            ARE SEPARATED FROM THE DATA
!
KGRID=IWORK(1)
KLENCH=IWORK(2)
IF (KLENCH < 0 .OR. KLENCH > JPXKRK) THEN
  IRESP=-58
  GOTO 1000
END IF
!
DATASIZE=ITOTAL-KLENCH-2
!
CALL GET_COMPHEADER(IWORK(3+KLENCH),DATASIZE,NEWSIZE,ITYPCOD)
IF (NEWSIZE >= 0) THEN
  ! compressed field found
  WRITE (ILUPRI,*) TRIM(HRECFM),' is compressed (old/new/kleng SIZE):',DATASIZE,NEWSIZE,KLENG 
  IF (KLENG /= NEWSIZE) THEN
    IRESP=-63
    GOTO 1000
  ENDIF

  ALLOCATE(IWORKNEW(NEWSIZE))
  CALL DECOMPRESS_FIELD(IWORKNEW,NEWSIZE,IWORK(3+KLENCH),DATASIZE,ITYPCOD)
  KFIELD(1:KLENG) = IWORKNEW(1:KLENG)
  DEALLOCATE(IWORKNEW)
ELSE
  IF (KLENG /= DATASIZE) THEN
    IRESP=-63
    GOTO 1000
  END IF
  KFIELD(1:KLENG)=IWORK(KLENCH+3:ITOTAL)
END IF
!
SELECT CASE (KLENCH)
CASE(-10:-1)
       IRESP=-58
       GOTO 1000
CASE(0)
       KFIELD(1:KLENG)=IWORK(3:ITOTAL)
CASE(1:JPXKRK)
       ICOMMENT(1:KLENCH)=IWORK(3:KLENCH+2)
       DO J=1,KLENCH
          HCOMMENT(J:J)=CHAR(ICOMMENT(J))
       ENDDO
CASE(JPXKRK+1:)
       IRESP=-56
       GOTO 1000
END SELECT
!
DEALLOCATE(IWORK)
!
!  this is a pure binary field: no uncompressing of any kind
!
!*      3.    MESSAGE PRINTING WHATEVER THE ISSUE WAS
!
1000    CONTINUE

IF (IRESP.NE.0) THEN
  YFNLFI=ADJUSTL(HFIPRI)
  DO J=1,JPNXLU
    IF (CNAMFI(J).EQ.YFNLFI) THEN
      ILUPRI=J
      EXIT
    ENDIF
  ENDDO
  WRITE (ILUPRI,*) ' exit from FMREAD with IRESP:',IRESP
  !WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
  WRITE (ILUPRI,*) '   | HRECFM = ',HRECFM
  !WRITE (ILUPRI,*) '   | KLENG  = ',KLENG
  !WRITE (ILUPRI,*) '   | KGRID  = ',KGRID
  !WRITE (ILUPRI,*) '   | KLENCH  = ',KLENCH
  ! Suppression OBLIGATOIRE de l'impression suivante car pb qd IWORK non alloue
  ! (IRESP=-47)
  !WRITE (ILUPRI,*) '   | KLENCH  = ',IWORK(23)
ENDIF
KRESP=IRESP

RETURN
      END SUBROUTINE FM_READ
