!     ######spl
      MODULE  MODD_CONVIJ2XY
!     ######################
!
!!****  *MODD_CONVIJ2XY* - 
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
!!     original        01/04/99
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVIJ
REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVI
REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVJ
REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVX
REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVY
REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVLAT
REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: XCONVLON

!
!
END MODULE MODD_CONVIJ2XY
