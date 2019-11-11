!     ######spl
      MODULE MODI_LOW2UP
!     #####################
!
INTERFACE
!
SUBROUTINE LOW2UP(HCARIN)
CHARACTER(LEN=*), INTENT(INOUT)  :: HCARIN
END SUBROUTINE LOW2UP
!
END INTERFACE
!
END MODULE MODI_LOW2UP
!     ######spl
      SUBROUTINE LOW2UP(HCARIN)
!     ############################
!
!!****  *LOW2UP* - 
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
!!      Original       06/06/94
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------
CHARACTER(LEN=*), INTENT(INOUT)  :: HCARIN
!
!*       0.1   Local variables
!              ---------------
!
CHARACTER(LEN=1),DIMENSION(26),SAVE  :: YLO=(/'a','b','c','d','e','f','g', &
 'h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/)
CHARACTER(LEN=1),DIMENSION(26),SAVE  :: YUP=(/'A','B','C','D','E','F','G', & 
 'H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/)
INTEGER   ::   ILENC
INTEGER   ::   J, JA
!------------------------------------------------------------------------------
!
ILENC = LEN(HCARIN)
!print *,' HCARIN ',LEN(HCARIN)
!print *,HCARIN
DO J=1,ILENC
  DO JA=1,26
    IF(HCARIN(J:J) == YLO(JA))HCARIN(J:J)=YUP(JA)
  ENDDO
ENDDO
!
END SUBROUTINE LOW2UP
