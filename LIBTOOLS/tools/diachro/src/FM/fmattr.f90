!     ######spl
      SUBROUTINE FMATTR(HFILEM,HFIPRI,KNUMBR,KRESP)
!     #############################################
!
!!****  *FMATTR* - routine to attribute a logical unit to a file name
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMATTR is to attribute to the file named HFILEM
!     the logical unit number KNUMBR chosen among the free logical units
!
!!**  METHOD
!!    ------
!!
!!      If FMATTR is called for the very first time, then all the management
!!    arrays used by the FM-routines are initialized in FMINIT. 
!!    Otherwise, the name HFILEM is searched in the array CNAMFI, where
!!    it should not exist ! Finally, a logical unit number is searched
!!    in array CNAMFI. As soon as a free place is found (CNAMFI=CPUDFN),
!!    this place becomes the logical unit number for HFILEM and CNAMFI is
!!    set to HFILEM.
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
!!      original                                                        04/94
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

INTEGER,              INTENT(OUT)::KNUMBR  ! logical unit number
INTEGER,              INTENT(OUT)::KRESP   ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER::IRESP=0,J,ILOGIQ=0,ILUPRI
CHARACTER(LEN=JPFINL)::YLOCFN,YLOCFN2
!
!*      0.3   Taskcommon for logical units
!
COMMON/TASKATTR/ILUPRI
!DIR$ TASKCOMMON TASKATTR
!
!----------------------------------------------------------------------------
!
!*      1.    INITIALISATION AND TEST THAT FILE DOES NOT ALREADY EXIST
!
IRESP = 0 ; ILOGIQ = 0 ; ILUPRI = 6
YLOCFN=HFILEM ; YLOCFN=ADJUSTL(YLOCFN)

!dino IF (LFMMUL) CALL LOCKON(NFMLOC)

IF (LFCATT) THEN
    CALL FMINIT
    LFCATT=.FALSE.
ELSE
    IF (NOPEFI.LT.0) THEN
       IRESP=-50
       GOTO 1000
    ELSE
       DO J=1,JPNXLU
         IF (YLOCFN.EQ.CNAMFI(J)) THEN
           IRESP=-51
           GOTO 1000
         ENDIF
       ENDDO
    ENDIF
ENDIF
!
!*     2.     WE LOOK FOR A FREE PLACE IN ARRAY CNAMFI
!
!   That place will become the number for the logical unit attributed to HFILEM
!

DO J=1,JPNXLU
    IF (CNAMFI(J).EQ.CPUDFN) THEN
       ILOGIQ=J
       CNAMFI(J)=YLOCFN
       EXIT
    ENDIF
ENDDO
IF (ILOGIQ.EQ.0) THEN
    IRESP=-52
    GOTO 1000
ENDIF

KNUMBR=ILOGIQ ; NOPEFI=NOPEFI+1

!dino IF (LFMMUL) CALL LOCKOFF(NFMLOC)

!
!*     3.     MESSAGE PRINTING WHATEVER THE ISSUE WAS
!
1000   CONTINUE

IF (IRESP.NE.0) THEN
   YLOCFN2=ADJUSTL(HFIPRI)
!
! in the special case where FMATTR is called to reserve a logical unit
! for the output file itself (i.e. HFILEM=HFIPRI),
! no print is performed because we do not know
! whether this file was actually opened or not.
!
   IF (YLOCFN2.EQ.YLOCFN) THEN
      ILUPRI=ILOGIQ
   ELSE
      DO J=1,JPNXLU
         IF (CNAMFI(J).EQ.YLOCFN2) THEN
            ILUPRI=J
            EXIT
         ENDIF
      ENDDO
      WRITE (ILUPRI,*) ' exit from FMATTR with IRESP:',IRESP
      WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
   ENDIF
ENDIF
KRESP=IRESP

RETURN
      END SUBROUTINE FMATTR
