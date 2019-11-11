!     ######spl
      SUBROUTINE FMLOOK(HFILEM,HFIPRI,KNUMBR,KRESP)
!     #############################################
!
!!****  *FMLOOK* - routine to look for the logical unit attributed to a file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMLOOK is to look for the logical unit (Fortran)
!     that is associated to the file named HFILEM. This unit was attributed
!     previously to HFILEM by FMATTR.
!
!!**  METHOD
!!    ------
!!
!!      The string HFILEM is searched in array CNAMFI which contains the
!!    names of all files that have been opened for the FM-routines.
!!    The place in array CNAMFI of HFILEM corresponds exactly to
!!    its logical unit.
!!
!!    EXTERNAL
!!    --------
!!
!!      NONE
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
!!      original                                                        04/94
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
CHARACTER(LEN=*),     INTENT(IN) ::HFILEM  ! file name

CHARACTER(LEN=*),     INTENT(IN) ::HFIPRI  ! file for prints in FM

INTEGER,              INTENT(OUT)::KNUMBR  ! logical unit number
INTEGER,              INTENT(OUT)::KRESP   ! return-code if problems araised

!
!*      0.2   Declarations of local variables
!
INTEGER::J,ILOGIQ=0,IRESP=0,ILUPRI
CHARACTER(LEN=JPFINL)::YLOCFN
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKLOOK/ILUPRI
!DIR$ TASKCOMMON TASKLOOK
!
!----------------------------------------------------------------------------
!
!*      1.    WE LOOK FOR THE FILE NAME IN ARRAY CNAMFI
!
ILOGIQ = 0 ; IRESP = 0 ; ILUPRI = 6
IF (NOPEFI.LT.1) THEN
     IRESP=-53
     GOTO 1000
ENDIF
YLOCFN=HFILEM ; YLOCFN=ADJUSTL(YLOCFN)
DO J=1,JPNXLU
     IF (YLOCFN.EQ.CNAMFI(J)) THEN
        ILOGIQ=J
        EXIT
     ENDIF
ENDDO
IF (ILOGIQ.EQ.0) THEN
     IRESP=-54
     GOTO 1000
ENDIF

KNUMBR=ILOGIQ
!
!*      2.     MESSAGE PRINTING WHATEVER THE ISSUE WAS
!
1000    CONTINUE

IF (IRESP.NE.0) THEN
YLOCFN=ADJUSTL(HFIPRI)
DO J=1,JPNXLU
    IF (CNAMFI(J).EQ.YLOCFN) THEN
       ILUPRI=J
       EXIT
    ENDIF
ENDDO
WRITE (ILUPRI,*) ' exit from FMLOOK with IRESP:',IRESP
WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
ENDIF
KRESP=IRESP

RETURN
      END SUBROUTINE FMLOOK
