!     ######spl
      MODULE MODD_TYPE_AND_LH
!     #######################
!
!!****  *MODD_TYPE_AND_LH* - 
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	J. Duron    *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/11/96      
!!              
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!
CHARACTER (LEN=4), SAVE :: CTYPE         ! type of data
!
INTEGER, SAVE :: NKL, NKH                ! lowest and highest K indice values 

LOGICAL, SAVE :: LKCP                    ! switch for compression in K
                                         ! direction
INTEGER, SAVE :: NIL, NIH                ! lowest and highest I indice values 

INTEGER, SAVE :: NJL, NJH                ! lowest and highest J indice values 

LOGICAL, SAVE :: LICP                    ! switch for compression in I
                                         ! direction
LOGICAL, SAVE :: LJCP                    ! switch for comppression in J
                                         ! direction
END MODULE MODD_TYPE_AND_LH
