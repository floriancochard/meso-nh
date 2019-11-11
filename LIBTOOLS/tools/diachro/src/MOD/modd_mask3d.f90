!     ######spl
      MODULE  MODD_MASK3D
!     ####################
!
!!****  *MODD_MASK3D* - 
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
!!     original        10/11/96
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


CHARACTER(LEN=16),SAVE :: CGROUPSV3
LOGICAL,SAVE :: LXYZ=.FALSE., LMASK3D=.FALSE., LMSKTOP=.FALSE.
LOGICAL,SAVE :: LSV3=.FALSE., LMARKER=.FALSE.
LOGICAL,SAVE :: LXYZ00=.FALSE.
LOGICAL,SAVE :: LMASK3D_XY=.FALSE.,LMASK3D_XZ=.FALSE.,LMASK3D_YZ=.FALSE.
LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: LXYZT
!
! Masque 3D
LOGICAL,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE    :: LMASK3
!
! limites fournies par l'utilisateur (x,y conformes z altitudes) pour
! definir un masque
REAL,SAVE   :: XXL=0.,XXH=0.,XYL=0.,XYH=0.,XZL=0.,XZH=0.
!
END MODULE MODD_MASK3D
