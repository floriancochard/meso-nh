!     ######spl
      MODULE  MODD_ALLOC2_FORDIACHRO
!     ##############################
!
!!****  *MODD_ALLOC2_FORDIACHRO* - 
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



REAL,DIMENSION(:,:,:,:,:,:), ALLOCATABLE,SAVE  :: XMASK2
INTEGER,DIMENSION(:), ALLOCATABLE,SAVE  :: NGRIDIA2

REAL,DIMENSION(:,:,:,:,:,:), ALLOCATABLE,SAVE  :: XVAR2
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE          :: XTRAJT2
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE        :: XTRAJX2, XTRAJY2, XTRAJZ2
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE            :: XDATIME2

CHARACTER*100,DIMENSION(:), ALLOCATABLE,SAVE   :: CTITRE2, CUNITE2, CCOMMENT2

!
!
END MODULE MODD_ALLOC2_FORDIACHRO
