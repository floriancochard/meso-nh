!     ######spl
      MODULE MODD_FMMULTI
!     ####################
!
!!****  *MODD_FMMULTI* - declaration of global variables for multitasked FM
!!
!!    PURPOSE
!!    -------
!
!       This module contains variables global to all models (tasks).
!     They are used to switch on (off) the multitasked mode in the File Manager.
!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!
!!      see "File structure and content in the Meso-NH model"
!!
!!    AUTHOR
!!    ------
!!
!!      C. FISCHER      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                        07/95
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
IMPLICIT NONE

INTEGER::NFMLOC                   ! identification of the lock
LOGICAL::LFMMUL=.FALSE.           ! becomes TRUE if multitasking is asked

      END MODULE MODD_FMMULTI
