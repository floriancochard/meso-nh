!     #############################################
      SUBROUTINE FMCLOS(HFILEM,HSTATU,HFIPRI,KRESP)
!     #############################################
!
!!****  *FMCLOS* - routine to close a meso-nh file opened with the "FM"-routines
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMCLOS is to close a mesonh file composed of the DESFM
!     and the LFIFM part. The LFIFM file is closed
!     using the LFI-package for direct access Fortran files. The DESFM file is
!     closed using a classical CLOSE statement.
!
!!**  METHOD
!!    ------
!!
!!      The closure is proceeded in 4 steps:
!!        1. close DESFM
!!        2. close LFIFM by calling LFIFER
!!        3. erase the file from the management arrays (FMFREE)
!!        4. the cpio and storage command is loaded into the pipe
!!           the pipe has the special fortran unit 10
!!
!!    EXTERNAL
!!    --------
!!
!!      FMLOOK,FMFREE,LFIFER,CLOSE,FLUSH,LOCKON,LOCKOFF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      MODULE: MODD_FMDECLAR contains management parameters and
!!              storage arrays to move information around at the
!!              level of all "FM"-routines.
!!              MODD_FMMULTI contains variables for multitasking
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
!!      modified by C. Fischer                    4/11/94 (write in the pipe)
!!      modified by C. Fischer                5/7/95 (locks for multitasking)
!!      modified by P. Jabouille                  26/06/96 (case NFITYP=2 :
!!                                     file is not sent to the remote machine)
!!      modified by V. Masson               16/09/96 (prints if error occurs)
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_FMMULTI

IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),     INTENT(IN) ::HFILEM  ! file name
CHARACTER(LEN=*),     INTENT(IN) ::HSTATU  ! status for the closed file

CHARACTER(LEN=*),     INTENT(IN) ::HFIPRI  ! file for prints in FM

INTEGER,              INTENT(OUT)::KRESP   ! return-code if problems araised

!
!*      0.2   Declarations of local variables
!
INTEGER::IRESP,IROWF,IPOSNU,J,INUMBR,IFMFNL,ILUPRI,IERR
CHARACTER(LEN=7)::YSTATU
CHARACTER(LEN=JPFINL)::YFNDES,YFNLFI
CHARACTER(LEN=LEN(HFILEM))::YINTFN
CHARACTER(LEN=10)::YTRANS,YCPIO
CHARACTER(LEN=100)::YCOMMAND
LOGICAL::GSTATU
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKCLOS/ILUPRI,INUMBR,IRESP,YFNDES,YFNLFI,YSTATU
!DIR$ TASKCOMMON TASKCLOS
!
!----------------------------------------------------------------------------
!
!*      1.1   THE NAME OF DESFM=HFILEM.des
!
IRESP = 0 ; IROWF = 0 ; IPOSNU = 0 ; ILUPRI = 6 ; IERR = 0
IFMFNL=JPFINL-4
YTRANS='transfer.x'

IROWF=LEN(HFILEM)

IF (IROWF.EQ.0) THEN
   IRESP=-59
   GOTO 1000
ELSEIF (IROWF.GT.IFMFNL) THEN
   IRESP=-60
   GOTO 1000
ENDIF
YINTFN=ADJUSTR(HFILEM)
YFNDES=YINTFN//'.des'
YFNDES=ADJUSTL(YFNDES)
!
!*      1.2   TEST FOR FILE EXISTENCE AND SEARCH OF ITS LOGICAL UNIT
!
CALL FMLOOK(YFNDES,HFIPRI,INUMBR,IRESP)
IF (IRESP.NE.0) THEN
        GOTO 1000
ELSEIF (LEN(HSTATU).LE.0) THEN
        IRESP=-41
        GOTO 1000
ELSE
        GSTATU=HSTATU.EQ.'KEEP'.OR.HSTATU.EQ.'DELETE'
        IF (GSTATU) THEN
        YSTATU=HSTATU(1:MIN0(LEN(HSTATU),LEN(YSTATU)))
        ELSE
        YSTATU='DEFAULT'
        ENDIF
ENDIF
!
!*      1.3   THE LOGICAL UNIT OF DESFM IS RELEASED FOR "FM"
!
CALL FMFREE(YFNDES,HFIPRI,IRESP)
IF (IRESP.NE.0) GOTO 1000
!
!*      2.    CLOSURE OF DESFM
!
!  case of a namelist
!
CLOSE (UNIT=INUMBR,IOSTAT=IRESP,STATUS=YSTATU)
IF (IRESP.NE.0) GOTO 1000
!
!*      3.1   THE NAME OF LFIFM=HFILEM.lfi
!
YFNLFI=YINTFN//'.lfi'
YFNLFI=ADJUSTL(YFNLFI)
!
!*      3.2   TEST FOR FILE EXISTENCE AND SEARCH OF ITS LOGICAL UNIT
!
CALL FMLOOK(YFNLFI,HFIPRI,INUMBR,IRESP)
IF (IRESP.NE.0) GOTO 1000
!
!*      3.3   THE LOGICAL UNIT FOR LFIFM IS RELEASED FOR "FM"
!
CALL FMFREE(YFNLFI,HFIPRI,IRESP)
IF (IRESP.NE.0) GOTO 1000
!
!*      4.    CLOSURE OF LFI
!
!  case of a LFI file
!
CALL LFIFER(IRESP,INUMBR,YSTATU)
IF (IRESP.NE.0) GOTO 1000
!
!*      5.    INPUT FOR THE UNIX SYSTEM TO SAVE AND SEND THE FILE
!
PRINT*,'KTYPE=',NFITYP(INUMBR)
SELECT CASE (NFITYP(INUMBR))
CASE(:-1)
  IRESP=-66
  GOTO 1000
CASE(0)
  YCPIO='NIL'
CASE(1)
  YCPIO='MESONH'
CASE(2)
  PRINT*,'FILE ',HFILEM,' NOT TRANSFERED'
  GOTO 1000
CASE(3:)
  IRESP=-66
  GOTO 1000
END SELECT
WRITE (YCOMMAND,20) YTRANS,YCPIO,HFILEM
!
! write into the pipe : the "flush" forces instanteneous buffer transfer
! which is necessary for parallel treatment
!
PRINT*,'YCOMMAND=',YCOMMAND
WRITE (10,'(A100)') YCOMMAND
!CALL FLUSH(10,IERR)
!
!*      6.    UPDATING OF ARRAY NFITYP
!
IF (LFMMUL) CALL LOCKON(NFMLOC)
NFITYP(INUMBR)=JPNIIL
IF (LFMMUL) CALL LOCKOFF(NFMLOC)
!
!*      7.    MESSAGE PRINTING WHATEVER THE ISSUE WAS
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
WRITE (ILUPRI,*) ' exit from FMCLOS with IRESP:',IRESP
WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
WRITE (ILUPRI,*) '   | HSTATU = ',HSTATU
ENDIF
KRESP=IRESP

! format: 10c for transfer.x and mesonh/nil
!         32c for file name
! if you have to change this format one day, don't forget the blank after 1H
20    FORMAT(A10,1H ,A10,1H ,A32)

RETURN
      END SUBROUTINE FMCLOS
