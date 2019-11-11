!     ######spl
      MODULE MODD_FMDECLAR
!     ####################
!
!!****  *MODD_FMDECLAR* - declaration of global variables of the FM-routines
!!
!!    PURPOSE
!!    -------
!
!       The purpose of MODD_FMDECLAR is to declare all the global variables that
!     are needed by the FM-routines. It includes specific FM-software parameters
!     as well as storage arrays. These arrays allow the FM-routines to keep
!     in mind which logical unit is associated with which file name
!     and to state whether a file was actually opened.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      NONE
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
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
IMPLICIT NONE

INTEGER,PARAMETER::JPNXLU=99   ! maximum number of logical units for Fortran
INTEGER,PARAMETER::JPNXFM=JPNXLU-10
                               ! maximum number of files opened at the same time
INTEGER,PARAMETER::JPXFIE=1.5E8! maximum record length for the FM-software
INTEGER,PARAMETER::JPNIIL=-999 ! default value in integer arrays
INTEGER,PARAMETER::JPFINL=32   ! length of the file name strings in FM
INTEGER,PARAMETER::JPXKRK=100  ! maximum length for the comment string

CHARACTER(LEN=JPFINL),PARAMETER::CPUDFN='UNDEFINED_FILE_NAME'
CHARACTER(LEN=JPFINL),PARAMETER::CPUNLU='UNAUTHORIZED_LOGICAL_UNIT'
!
!----------------------------------------------------------------------------
INTEGER::NOPEFI                ! number of opened files

INTEGER,DIMENSION(1:JPNXLU)::NFITYP
                               ! NFITYP contains the type of the FM file which
                               ! will be used in FMCLOS for the Unix save.

CHARACTER(LEN=JPFINL),DIMENSION(1:JPNXLU)::CNAMFI
                               ! management array containing the names of all
                               ! opened files

LOGICAL::LFCATT=.TRUE.         ! This logical is true at the very first call
                               ! to FMATTR and is then set to false.

END MODULE MODD_FMDECLAR
