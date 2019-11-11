!     ######spl
      MODULE MODD_READLH
!     #######################
!
!!****  *MODD_READLH* - 
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
!!	N. Asencio    *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/03/05      
!!              
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, SAVE :: NREADKL, NREADKH        ! lowest and highest K indice values 

INTEGER, SAVE :: NREADIL, NREADIH        ! lowest and highest I indice values 

INTEGER, SAVE :: NREADJL, NREADJH        ! lowest and highest J indice values 

END MODULE MODD_READLH
