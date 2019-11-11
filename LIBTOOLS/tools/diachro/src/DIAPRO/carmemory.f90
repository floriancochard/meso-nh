!     ######spl
      MODULE MODI_CARMEMORY
!     #####################
!
INTERFACE
!
SUBROUTINE CARMEMORY(HCARIN,KOP)
CHARACTER(LEN=*),INTENT(INOUT)  :: HCARIN
!CHARACTER(LEN=2400),INTENT(INOUT) :: HCARIN
INTEGER          :: KOP 
END SUBROUTINE CARMEMORY
!
END INTERFACE
!
END MODULE MODI_CARMEMORY
!     ######spl
      SUBROUTINE CARMEMORY(HCARIN,KOP)
!     ################################
!
!!****  *CARMEMORY* - 
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
!
!CHARACTER(LEN=2400),INTENT(INOUT)  :: HCARIN
CHARACTER(LEN=*),INTENT(INOUT)  :: HCARIN
INTEGER          :: KOP
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=2400),SAVE :: YCAR
INTEGER,SAVE   ::   ILENC, ILENGP1
!------------------------------------------------------------------------------
!
IF(KOP == 1)THEN
!fuji  HCARIN=ADJUSTL(HCARIN)   !introduit des caracteres genre {Á€W×?Ã
   HCARIN=TRIM(HCARIN)
  YCAR(1:LEN(YCAR))=' '
  YCAR=ADJUSTL(HCARIN)
  ILENC = LEN(YCAR)
if (nverbia > 0)then
!print *, ' *** CARMEMORY 1 ILENC YCAR ',ILENC,YCAR(1:80)
print *, ' *** CARMEMORY 1 ILENC YCAR ',ILENC,YCAR(1:LEN_TRIM(YCAR))
endif
ELSE IF(KOP == 2)THEN
  HCARIN(1:LEN(HCARIN))=' '
  HCARIN=ADJUSTL(YCAR(ILENGP1+1:LEN_TRIM(YCAR)))
  HCARIN=ADJUSTL(HCARIN)
ELSE IF(KOP == 3)THEN
  CGROUPS(1)=ADJUSTL(CGROUPS(1))
  ILENGP1=LEN_TRIM(CGROUPS(1))
ENDIF

!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE CARMEMORY
