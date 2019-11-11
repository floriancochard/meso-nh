!     ######spl
      MODULE  MODD_DIMGRID_FORDIACHRO
!     ###############################
!
!!****  *MODD_DIMGRID_FORDIACHRO* - 
!!
!!    PURPOSE
!!    -------
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        01/02/96
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE

INTEGER,SAVE  :: NNB, NNBF

INTEGER,DIMENSION(:,:), ALLOCATABLE,SAVE  :: NNUMT, NSIZT, NLENC

CHARACTER*16,DIMENSION(:,:), ALLOCATABLE,SAVE :: CRECFM2T    
!
!
END MODULE MODD_DIMGRID_FORDIACHRO
