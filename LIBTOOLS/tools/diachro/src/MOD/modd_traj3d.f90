!     ######spl
      MODULE  MODD_TRAJ3D
!     ####################
!
!!****  *MODD_TRAJ3D* - 
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
!!      JS    "MF"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        10/04/00
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE

INTEGER, PARAMETER :: NPART_MAX=100
REAL, DIMENSION(NPART_MAX),SAVE :: XXPART 
REAL, DIMENSION(NPART_MAX),SAVE :: XYPART
REAL, DIMENSION(NPART_MAX),SAVE :: XZPART
LOGICAL,SAVE :: LTRAJ3D=.FALSE.
LOGICAL,SAVE :: LFLUX3D=.FALSE.
INTEGER,SAVE :: NPART
LOGICAL,SAVE :: LTRAJ_GROUP=.FALSE.
CHARACTER (LEN=16), SAVE :: CTRAJ_GROUP
!
END MODULE MODD_TRAJ3D
