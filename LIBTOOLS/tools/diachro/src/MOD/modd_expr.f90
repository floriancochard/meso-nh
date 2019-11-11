!     ######spl
      MODULE  MODD_EXPR
!     #############################
!
!!****  *MODD_EXPR* - 
!!
!!    PURPOSE
!!    -------
!       Declaration des variables et tableaux intervenant dans la
!       multiplication (ou division) d'un processus par un autre
!       processus
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
!!     original        02/07/01
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


INTEGER,SAVE :: NMD

REAL,DIMENSION(:,:,:,:,:,:), ALLOCATABLE,SAVE  :: XEXPR1,XEXPR2,XEXPR3,&
XEXPR4,XEXPR5,XEXPR6,XEXPR7,XEXPR8,XEXPR9
REAL,DIMENSION(:,:,:,:,:,:), ALLOCATABLE,SAVE  :: XDEXPR1,XDEXPR2,XDEXPR3,&
XDEXPR4,XDEXPR5,XDEXPR6,XDEXPR7,XDEXPR8,XDEXPR9
INTEGER,DIMENSION(100), SAVE  :: NMULTDIV
CHARACTER(LEN=6),DIMENSION(50),SAVE :: CMULTDIV

!
END MODULE MODD_EXPR
