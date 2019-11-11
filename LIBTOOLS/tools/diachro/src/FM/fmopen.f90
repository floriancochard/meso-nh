!     ######spl
      SUBROUTINE FMOPEN(HFILEM,HSTATU,HFIPRI,KNPRAR,KFTYPE,KVERB,&
                        KNINAR,KRESP)
!     ############################################################
!
!!****  *FMOPEN* - routine to open a meso-nh file (DESFM+LFIFM)
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMOPEN is to open a meso-nh file for the "FM"-routines.
!     It is composed of two distinct fortran files: DESFM and LFIFM. DESFM is
!     a namelist formatted file. LFIFM is a LFI file, managed by the LFI-package.
!     LFIFM is a fortran unformatted, direct access file which is
!     manipulated by the FM-routines FMREAD and FMWRIT. 
!     The namelist file is a fortran 90 standard formatted file.
!
!!**  METHOD
!!    ------
!!
!!      The opening is performed in 4 main steps:
!!            1. a logical unit is reserved for DESFM (first call to FMATTR)
!!            2. the DESFM file is created by a
!!               formatted, fortran open. The name of the file is obtained by
!!               appending ".des" to HFILEM.
!!            3. a logical unit is reserved for LFIFM (second call to FMATTR)
!!            4. the LFIFM file is opened in the LFIOUV routine to
!!               which most of the explicit input arguments of FMOPEN are passed.
!!               The name of that file is obtained by appending ".lfi"
!!               to HFILEM.
!!
!!    EXTERNAL
!!    --------
!!
!!      FMATTR,LFIOUV,OPEN,LOCKON,LOCKOFF
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
!!      modified by C. Fischer                5/7/95 (locks for multitasking)
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
CHARACTER(LEN=*),     INTENT(IN) ::HFILEM  ! name of the file
CHARACTER(LEN=*),     INTENT(IN) ::HSTATU  ! status of the file at opening
CHARACTER(LEN=*),     INTENT(IN) ::HFIPRI  ! file for prints in FM

INTEGER,          INTENT(IN) ::KNPRAR  ! number of predicted articles (not vital)
INTEGER,          INTENT(IN) ::KFTYPE  ! type of FM-file
INTEGER,          INTENT(IN) ::KVERB   ! level of verbose

INTEGER,          INTENT(OUT)::KNINAR  ! number of articles initially present in the file
INTEGER,          INTENT(OUT)::KRESP   ! return-code if a problem araised

!
!*      0.2   Declarations of local variables
!
INTEGER::IRESOU,INPRAR,IROWF,IRESP,J,INUMBR,IFMFNL,IMELEV,ILUPRI
CHARACTER(LEN=JPFINL)::YFNDES,YFNLFI
CHARACTER(LEN=LEN(HFILEM))::YINTFN
LOGICAL::GNEWFI,GNAMFI=.TRUE.,GFATER=.TRUE.,GSTATS
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKOPEN/ILUPRI,INUMBR,IRESP,YFNDES,YFNLFI
!DIR$ TASKCOMMON TASKOPEN
!
!----------------------------------------------------------------------------
!
!*      1.    INITIALIZATION
!
INPRAR=KNPRAR+0;KNINAR=0
IRESOU = 0 ; IROWF = 0 ; IRESP = 0 ; ILUPRI = 6
!
!* the model's verbose level is connected to the LFI verbose
!
SELECT CASE (KVERB)
CASE(:2)
   GSTATS=.FALSE. ; IMELEV=0
CASE(3:6)
   GSTATS=.FALSE. ; IMELEV=1
CASE(7:9)
   GSTATS=.FALSE. ; IMELEV=2
CASE(10:)
   GSTATS=.TRUE. ; IMELEV=2
END SELECT

IF (NOPEFI.GE.JPNXFM) THEN
        IRESP=-44
        GOTO 1000
ENDIF
!
!*      2.    LOGICAL UNIT FOR DESFM
!
!  the fortran name for DESFM
!
IFMFNL=JPFINL-4

IROWF=LEN(HFILEM)

IF (IROWF.EQ.0) THEN
   IRESP=-45
   GOTO 1000
ELSEIF (IROWF.GT.IFMFNL) THEN
   IRESP=-49
   GOTO 1000
ENDIF
YINTFN=ADJUSTR(HFILEM)
YFNDES=YINTFN//'.des'
YFNDES=ADJUSTL(YFNDES)

CALL FMATTR(YFNDES,HFIPRI,INUMBR,IRESP)
IF (IRESP.NE.0) GOTO 1000

!
!*      3.    FILE OPENING FOR DESFM
!
!  case of a namelist: sequential, formatted fortran open
!
OPEN(UNIT=INUMBR,FILE=YFNDES,FORM='FORMATTED',DELIM='QUOTE',IOSTAT=IRESP)
IF (IRESP.NE.0) GOTO 1000
!
!*      4.    LOGICAL UNIT FOR LFIFM
!
!  the fortran name for LFIFM
!
YFNLFI=YINTFN//'.lfi'
YFNLFI=ADJUSTL(YFNLFI)

CALL FMATTR(YFNLFI,HFIPRI,INUMBR,IRESP)
IF (IRESP.NE.0) GOTO 1000
!
!*      5.    FILE OPENING FOR LFIFM
!
!  case of a LFI-file: direct access, unformatted open via LFIOUV
!
CALL LFIOUV(IRESOU,INUMBR,GNAMFI,YFNLFI,HSTATU,GFATER,GSTATS,IMELEV,INPRAR,&
            KNINAR)
IF (IRESOU.NE.0.AND.IRESOU.NE.-11) THEN
        IRESP=IRESOU
        GOTO 1000
ENDIF

!
!*      6.    TEST IF FILE IS NEWLY DEFINED
!

GNEWFI=(KNINAR.EQ.0).OR.(KVERB.LT.7)
IF (.NOT.GNEWFI) THEN
YFNLFI=ADJUSTL(HFIPRI)
DO J=1,JPNXLU
    IF (CNAMFI(J).EQ.YFNLFI) THEN
       ILUPRI=J
       EXIT
    ENDIF
ENDDO
WRITE (ILUPRI,*) ' file ',INUMBR,'previously created with LFI'
ENDIF
!
!*      7.    UPDATE OF THE FILE TYPE ARRAY
!
!dino IF (LFMMUL) CALL LOCKON(NFMLOC)
NFITYP(INUMBR)=KFTYPE
!dino IF (LFMMUL) CALL LOCKOFF(NFMLOC)
!
!*      8.    MESSAGE PRINTING WHATEVER THE ISSUE WAS
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
WRITE (ILUPRI,*) ' exit from FMOPEN with IRESP:',IRESP
WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
WRITE (ILUPRI,*) '   | HSTATU = ',HSTATU
WRITE (ILUPRI,*) '   | KNPRAR = ',KNPRAR
WRITE (ILUPRI,*) '   | KFTYPE = ',KFTYPE
ENDIF
KRESP=IRESP

RETURN
      END SUBROUTINE FMOPEN
