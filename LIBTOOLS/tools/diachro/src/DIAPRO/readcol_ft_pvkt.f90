!     ######spl
      MODULE MODI_READCOL_FT_PVKT
!     ############################
!
INTERFACE
!
SUBROUTINE READCOL_FT_PVKT(HCARIN,KCOLI)
CHARACTER(LEN=*)  :: HCARIN
INTEGER           :: KCOLI
END SUBROUTINE READCOL_FT_PVKT
!
END INTERFACE
END MODULE MODI_READCOL_FT_PVKT
!     ######spl
      SUBROUTINE READCOL_FT_PVKT(HCARIN,KCOLI)
!     ########################################
!
!!****  *READCOL_FT_PVKT* - 
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
!!      Original       2/09/96
!!      Updated   PM   
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

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KCOLI
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: J,JM
CHARACTER(LEN=LEN(HCARIN)) :: YCARIN, YCARIN2

!
!------------------------------------------------------------------------------
KCOLI=0
IF(NBCOLI == 0)THEN
  RETURN
ELSE
  YCARIN(1:LEN(YCARIN))=' '
  YCARIN=ADJUSTL(HCARIN)
  JM=0
  DO J=1,LEN(YCARIN)
    IF(YCARIN(J:J) == ' ')THEN
      JM=J-1
      EXIT
    ENDIF
  ENDDO
  IF(JM /= 0)THEN
    YCARIN2(1:LEN(YCARIN2))=' '
    YCARIN2=YCARIN(1:JM)
    YCARIN(1:LEN(YCARIN))=' '
    YCARIN=ADJUSTL(YCARIN2)
  ENDIF
  DO J=1,NBCOLI
    IF(YCARIN == CCOLI(J))THEN
      KCOLI=NCOLI(J)
      EXIT
    ENDIF
  ENDDO
  RETURN
ENDIF
END SUBROUTINE READCOL_FT_PVKT
