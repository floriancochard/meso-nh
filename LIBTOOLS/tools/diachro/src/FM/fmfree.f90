!     ######################################
      SUBROUTINE FMFREE(HFILEM,HFIPRI,KRESP)
!     ######################################
!
!!****  *FMFREE* - routine to release a logical unit for FM
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMFREE is to free the logical unit attributed to
!     the file named HFILEM.
!
!!**  METHOD
!!    ------
!!
!!      The association between the file named HFILEM and its logical unit
!!    (ILOGIQ, say) was performed by a previous call to FMATTR. This link
!!    is broken by setting the value CNAMFI(ILOGIQ) back to CPUDFN, so that
!!    HFILEM does not appear anymore in CNAMFI.
!!
!!    EXTERNAL
!!    --------
!!
!!      LOCKON,LOCKOFF
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
CHARACTER(LEN=*),     INTENT(IN) ::HFILEM  ! file name

CHARACTER(LEN=*),     INTENT(IN) ::HFIPRI  ! file for prints in FM

INTEGER,              INTENT(OUT)::KRESP   ! return-code if problems araised

!
!*      0.2   Declarations of local variables
!
INTEGER::IRESP=0,J,ILOGIQ=0,ILUPRI
CHARACTER(LEN=JPFINL)::YLOCFN,YLOCFN2
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKFREE/ILUPRI
!DIR$ TASKCOMMON TASKFREE
!
!----------------------------------------------------------------------------
!
!*      1.    THE NAME IS SEARCHED IN CNAMFI AND ERASED
!
IRESP = 0 ; ILOGIQ = 0 ; ILUPRI = 6
YLOCFN=HFILEM ; YLOCFN=ADJUSTL(YLOCFN)

IF (LFMMUL) CALL LOCKON(NFMLOC)

DO J=1,JPNXLU
   IF (YLOCFN.EQ.CNAMFI(J)) THEN
      ILOGIQ=J
      CNAMFI(J)=CPUDFN
      EXIT
   ENDIF
ENDDO
IF (ILOGIQ.EQ.0) THEN
   IRESP=-42
   GOTO 1000
ENDIF

NOPEFI=NOPEFI-1

IF (LFMMUL) CALL LOCKOFF(NFMLOC)

!
!*      2.    MESSAGE PRINTING WHATEVER THE ISSUE WAS
!
1000    CONTINUE

IF (IRESP.NE.0) THEN
   YLOCFN2=ADJUSTL(HFIPRI)
   IF (YLOCFN2.EQ.YLOCFN) THEN
! special case where HFILEM is the output listing itself: no print in this case
! because we do not know whether this file has already been closed or not
      ILUPRI=ILOGIQ
   ELSE
! most common case is this one
      DO J=1,JPNXLU
         IF (CNAMFI(J).EQ.YLOCFN2) THEN
            ILUPRI=J
            EXIT
         ENDIF
      ENDDO
   WRITE (ILUPRI,*) ' exit from FMFREE with IRESP:',IRESP
   WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
   ENDIF
ENDIF
KRESP=IRESP

RETURN
      END SUBROUTINE FMFREE
