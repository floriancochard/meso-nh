!     ######spl
      MODULE MODI_CAREAL
!     ##################
!
INTERFACE
!
SUBROUTINE CAREAL(HCAR,POUT)
CHARACTER(LEN=*) :: HCAR
REAL             :: POUT
END SUBROUTINE CAREAL
!
END INTERFACE
!
END MODULE MODI_CAREAL
!     ######spl
      SUBROUTINE CAREAL(HCAR,POUT)
!     ############################
!
!!****  *CAREAL* - 
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
REAL             :: POUT
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=LEN(HCAR)) :: YCAR
!------------------------------------------------------------------------------
!
YCAR=HCAR
READ(YCAR,*)POUT

!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE CAREAL
