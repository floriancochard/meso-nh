!     ######spl
      SUBROUTINE COMPLAT(PLAT)
!     ############################
!
!!****  *COMPLAT* - 
!!****           
!!
!!    PURPOSE
!!    -------
!   
!
!!**  METHOD
!!    ------
!!  
!! 
!!
!!    EXTERNAL
!!    --------
!!      COS  ! trigonometric functions
!!      SIN  !
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
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
!!    MODIFICATIONS
!!    -------------
!!      Original       22/02/2000
!!      Updated   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_COORD
USE MODD_NMGRID
USE MODD_GRID1
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODE_GRIDPROJ 
!
IMPLICIT NONE
!
!*       0.1  Dummy arguments and results
!
REAL, DIMENSION(:,:),  INTENT(OUT) :: PLAT  
!

!*       0.2  Local variables
!
INTEGER             :: II, IJ
INTEGER             :: JILOOP, JJLOOP
!
REAL,DIMENSION(:), ALLOCATABLE,SAVE :: ZY
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: ZLA, ZLO, ZYY, ZX
!
!-------------------------------------------------------------------------------
!
!*        1.   COMPUTING THE LONGITUDINAL AND TRANSVERSE COMPONENTS
!              ----------------------------------------------------
!
!*        1.1  Array sizes calculations
!
II=SIZE(PLAT,1)
IJ=SIZE(PLAT,2)
!
!*        1.2  Array allocations
!
IF (ALLOCATED(ZX))THEN
  DEALLOCATE(ZX)
ENDIF
IF (ALLOCATED(ZY))THEN
  DEALLOCATE(ZY)
ENDIF
IF (ALLOCATED(ZYY))THEN
  DEALLOCATE(ZYY)
ENDIF
IF (ALLOCATED(ZLA))THEN
  DEALLOCATE(ZLA)
ENDIF
IF (ALLOCATED(ZLO))THEN
  DEALLOCATE(ZLO)
ENDIF

ALLOCATE(ZX(II,1),ZY(IJ))
ALLOCATE(ZYY(II,1),ZLA(II,1),ZLO(II,1))
!
ZX(:,1)=XXX(:,NMGRID)
ZY(:)=XXY(:,NMGRID)
DO JJLOOP=1,IJ
  DO JILOOP=1,II
    ZYY(JILOOP,1)=ZY(JJLOOP)
  ENDDO
  CALL SM_LATLON_A(XLATORI,XLONORI,ZX,ZYY,ZLA,ZLO)
  PLAT(:,JJLOOP)=ZLA(:,1)
ENDDO
!------------------------------------------------------------------------------
!
!*        2.     EXIT
!                ----
!
RETURN
END SUBROUTINE COMPLAT
