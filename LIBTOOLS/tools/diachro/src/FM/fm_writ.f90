!     ###########################################################
      SUBROUTINE FM_WRIT(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                        KLENCH,HCOMMENT,KRESP)
!     ###########################################################
!
!!****  *FM_WRIT* - routine to write a single data article into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRIT is to write one article into a Meso-nh data file.
!     This routine only holds for a LFI-file (not namelist).
!
!!**  METHOD
!!    ------
!!
!!      The unformatted write operation is actually performed by the routine
!!    LFIECR. You need to indicate the file name without the ".lfi"
!!    suffix, the data array and the
!!    length of this array. Furthermore, you have to give a name for the article
!!    you are writing (string) which you better choose by convention.
!!      FMWRIT also appends the grid-indicator (KGRID) at the beginning of
!!    the LFI logical article (IWORK(1)) ; then the length of the comment
!!    string (KLENCH) ; then the comment string itself which is first
!!    converted into integer type using ICHAR.
!!    Finally, it writes the data (integer or
!!    real) itself (rest of array IWORK). We stress that the length KLENG
!!    that the user has to indicate is the length of the real data array
!!    WITHOUT taking the other fields into account.
!!
!!    EXTERNAL
!!    --------
!!
!!      FMLOOK,LFIECR,ICHAR
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
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER(KIND=8),DIMENSION(1:KLENG),INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=KLENCH),     INTENT(IN) ::HCOMMENT ! comment string)

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised

!
!*      0.2   Declarations of local variables
!
INTEGER::IRESP,ITOTAL,INUMBR,J,IROW,IFMFNL,ILUPRI
INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE::IWORK
INTEGER,DIMENSION(1:JPXKRK)::ICOMMENT
CHARACTER(LEN=JPFINL)::YFNLFI
CHARACTER(LEN=LEN(HFILEM))::YINTFN
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKWRIT/ILUPRI,INUMBR,IRESP
!DIR$ TASKCOMMON TASKWRIT
!
!----------------------------------------------------------------------------
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0 ; IROW = 0 ; ILUPRI = 6
IFMFNL=JPFINL-4

IROW=LEN(HFILEM)

IF (IROW.EQ.0) THEN
   IRESP=-64
   GOTO 1000
ELSEIF (IROW.GT.IFMFNL) THEN
   IRESP=-65
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
!*      2.    GRID INDICATOR, COMMENT AND DATA ARE PUT TOGETHER
!
IF (KLENG.LE.0) THEN
    IRESP=-40
    GOTO 1000
ELSEIF (KLENG.GT.JPXFIE) THEN
    IRESP=-43
    GOTO 1000
ELSEIF ((KGRID.LT.0).OR.(KGRID.GT.8)) THEN
    IRESP=-46
    GOTO 1000
ENDIF

ITOTAL=KLENG+1+KLENCH+1
IF(ALLOCATED(IWORK)) DEALLOCATE(IWORK)
ALLOCATE(IWORK(ITOTAL))

IWORK(1)=KGRID

SELECT CASE (KLENCH)
CASE(:-1)
    IRESP=-55
    GOTO 1000
CASE(0)
    IWORK(2)=KLENCH
    IWORK(3:KLENG+2)=KFIELD(1:KLENG)
CASE(1:JPXKRK)
    DO J=1,KLENCH
        ICOMMENT(J)=ICHAR(HCOMMENT(J:J))
    ENDDO
    IWORK(2)=KLENCH
    IWORK(3:KLENCH+2)=ICOMMENT(1:KLENCH)
    IWORK(KLENCH+3:ITOTAL)=KFIELD(1:KLENG)
CASE(JPXKRK+1:)
    IRESP=-57
    GOTO 1000
END SELECT

!
!  no compressing of any kind: the data is pure binary
!
!*      3.    UNFORMATTED, DIRECT ACCESS WRITE OPERATION
!
CALL LFIECR(IRESP,INUMBR,HRECFM,IWORK,ITOTAL)
IF (IRESP.NE.0) GOTO 1000

DEALLOCATE(IWORK)
!
!*      4.    MESSAGE PRINTING WHATEVER THE ISSUE WAS
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
WRITE (ILUPRI,*) ' exit from FMWRIT with IRESP:',IRESP
WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
WRITE (ILUPRI,*) '   | HRECFM = ',HRECFM
WRITE (ILUPRI,*) '   | KLENG  = ',KLENG
WRITE (ILUPRI,*) '   | KGRID  = ',KGRID
WRITE (ILUPRI,*) '   | KLENCH = ',KLENCH
ENDIF
KRESP=IRESP

RETURN
      END SUBROUTINE FM_WRIT
