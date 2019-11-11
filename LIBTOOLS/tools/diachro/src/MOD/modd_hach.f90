!     ######spl
      MODULE  MODD_HACH
!     #################
!
!!****  *MODD_HACH* - 
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


! Surfaces en hachures
LOGICAL :: LHACH1=.FALSE., LHACH2=.FALSE., LHACH3=.FALSE., LHACH4=.FALSE., LHACHSEL=.FALSE.
! Surfaces en grises
LOGICAL :: LGREY=.FALSE.
! Label sur la 1ere isoligne
LOGICAL :: LABEL1=.TRUE.
LOGICAL :: LBLUSER1=.FALSE., LBLUSER2=.FALSE., LBLUSER3=.FALSE., LBLUSER4=.FALSE.
INTEGER,SAVE :: NLBL1, NLBL2, NLBL3, NLBL4
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: XLBLUSER1,XLBLUSER2,XLBLUSER3,XLBLUSER4
!
!
END MODULE MODD_HACH
