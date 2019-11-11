!     ######spl
      SUBROUTINE FMINIT
!     #################
!
!!****  *FMINIT* - routine to initialize the management arrays used by the FM-routines
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMINIT is to initialize the management arrays used
!     by the other FM-routines. These arrays allow to associate each logical
!     unit number with the given file name.
!     FMINIT is only called when FMATTR is called for the very
!     first time.
!       Furthermore, FMINIT opens unit 10 which is dedicated to the pipe
!     in which the transfer orders are written (in FMCLOS). Thus, unit 10
!     is specific and unavailable for common file management.
!
!!**  METHOD
!!    ------
!!
!!      Array intrinsics of fortran 90 are used
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
!!      original                                                        06/94
!!      modified by C. Fischer                    22/11/94  (open unit 10)
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR

IMPLICIT NONE
!----------------------------------------------------------------------------

NOPEFI=0

NFITYP=JPNIIL

CNAMFI=CPUDFN ; CNAMFI(1:10)=CPUNLU

!OPEN(UNIT=10,FILE='pipe_name',FORM='FORMATTED')

RETURN
      END SUBROUTINE FMINIT
