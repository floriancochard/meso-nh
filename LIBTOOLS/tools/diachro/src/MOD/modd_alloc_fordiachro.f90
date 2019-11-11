!     ######spl
      MODULE  MODD_ALLOC_FORDIACHRO
!     #############################
!
!!****  *MODD_ALLOC_FORDIACHRO* - 
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


INTEGER,SAVE :: NGRID, NGRIDIAM

REAL,DIMENSION(:,:,:,:,:,:), ALLOCATABLE,SAVE  :: XMASK
INTEGER,DIMENSION(:), ALLOCATABLE,SAVE  :: NGRIDIA

REAL,DIMENSION(:,:,:,:,:,:), ALLOCATABLE,SAVE  :: XVAR 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE          :: XTRAJT
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE        :: XTRAJX, XTRAJY, XTRAJZ
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE            :: XDATIME

CHARACTER*100,DIMENSION(:), ALLOCATABLE,SAVE   :: CTITRE, CUNITE, CCOMMENT

LOGICAL :: LPBREAD=.FALSE.
!
!
END MODULE MODD_ALLOC_FORDIACHRO
