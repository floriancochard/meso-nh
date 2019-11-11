!     ######spl
      MODULE MODI_CARINT
!     ##################
!
INTERFACE
!
SUBROUTINE CARINT(HCAR,KOUT)
CHARACTER(LEN=*) :: HCAR
INTEGER          :: KOUT
END SUBROUTINE CARINT
!
END INTERFACE
!
END MODULE MODI_CARINT
!     ######spl
      SUBROUTINE CARINT(HCAR,KOUT)
!     ############################
!
!!****  *CARINT* - 
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       24/11/95
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCAR
INTEGER          :: KOUT
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=LEN(HCAR)) :: YCAR
!------------------------------------------------------------------------------
!
YCAR=HCAR
READ(YCAR,*)KOUT

!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE CARINT
