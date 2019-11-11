!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_CHANGE_A_GRID
!#################################
!
INTERFACE
      SUBROUTINE CHANGE_A_GRID(PFIELD,KGRID,PFIELDA,KLUOUT)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD      ! values of the field
INTEGER,                INTENT(INOUT) :: KGRID       ! Mesonh grid indicator
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELDA    ! values of the field on the A-grid
INTEGER, INTENT(IN), OPTIONAL       :: KLUOUT      ! unit number of listing
!
END SUBROUTINE CHANGE_A_GRID
END INTERFACE
END MODULE MODI_CHANGE_A_GRID
!     ######spl
      SUBROUTINE CHANGE_A_GRID(PFIELD,KGRID,PFIELDA,KLUOUT)
!     #####################
!
!!****  *CHANGE_A_GRID* - change flux point variables to mass points
!!                         
!!
!!    PURPOSE
!!    -------
!!    
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!      Functions MXF, MYF, MZF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0  : name of output-listing
!!      Module MODD_FIELD1    : contains prognostics  variables 
!!         XUT
!!         XVT
!!         XWT
!!      Module MODD_GRID1
!!         XZZ
!!      Module MODD_DIAG_FIELD1
!!         XUAT
!!         XVAT
!!         XWAT
!!         XZA
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Ducrocq  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/03/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
! 
USE MODD_CONF, ONLY : NVERB
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD      ! values of the field
INTEGER,                INTENT(INOUT) :: KGRID       ! Mesonh grid indicator
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELDA    ! values of the field on the A-grid
INTEGER, INTENT(IN), OPTIONAL       :: KLUOUT      ! unit number of listing
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IIU,IJU,IKU      ! End of arrays
!-------------------------------------------------------------------------------
!
!*         1.     GENERAL CASE
!                 ------------
IKU= SIZE(PFIELD,3)
IIU= SIZE(PFIELD,1)
IJU= SIZE(PFIELD,2)
!
SELECT CASE(KGRID)
  CASE(1)
    IF(PRESENT(KLUOUT)) THEN
      WRITE(KLUOUT,*) ' CHANGE_A_GRID: case 1'
    ELSE
      PRINT*,' CHANGE_A_GRID: case 1'
    ENDIF
    PFIELDA(:,:,:) = PFIELD(:,:,:)
  CASE(2)
    IF(PRESENT(KLUOUT)) THEN
      WRITE(KLUOUT,*) ' CHANGE_A_GRID: case 2'
    ELSE
      PRINT*,' CHANGE_A_GRID: case 2'
    ENDIF
    PFIELDA(:,:,:) = MXF(PFIELD(:,:,:)) 
    PFIELDA(IIU,:,:)=2.*PFIELD(IIU,:,:)-PFIELD(IIU-1,:,:)
    KGRID=1
  CASE(3)
    IF(PRESENT(KLUOUT)) THEN
      WRITE(KLUOUT,*) ' CHANGE_A_GRID: case 3'
    ELSE
      PRINT*,' CHANGE_A_GRID: case 3'
    ENDIF
    PFIELDA(:,:,:) = MYF(PFIELD(:,:,:)) 
    PFIELDA(:,IJU,:)=2*PFIELD(:,IJU,:)-PFIELD(:,IJU-1,:)
    KGRID=1
  CASE(4)
    IF(PRESENT(KLUOUT)) THEN
      WRITE(KLUOUT,*) ' CHANGE_A_GRID: case 4'
    ELSE
      PRINT*,' CHANGE_A_GRID: case 4'
    ENDIF
    PFIELDA(:,:,:) = MZF(PFIELD(:,:,:)) 
    PFIELDA(:,:,IKU)=2*PFIELD(:,:,IKU)-PFIELD(:,:,IKU-1)
    KGRID=1
END SELECT
!
!-------------------------------------------------------------------------------
!
IF (NVERB>=10 .AND. PRESENT(KLUOUT)) &
  WRITE(KLUOUT,*) 'routine CHANGE_A_GRID completed'
!
END SUBROUTINE CHANGE_A_GRID
