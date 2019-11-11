!
!      #######################
          MODULE MODD_DRAG_n
!      #######################
!
!!****   *MODD_DRAG*  - declaration of drag (no-slip) parameters
!!
!!     PURPOSE
!!     -------
!!       The purpose of this declarative module is to declare the drag booleans
!!       and parameters
!!
!!**   IMPLICIT ARGUMENTS
!!     ------------------
!!       NONE
!!
!!    
!!     AUTHOR
!!     ------
!1       Jeanne Colin         * Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original            28/10/11
!----------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

!
TYPE DRAG_t

LOGICAL :: LDRAG          ! Logical switch to activate the no-slip condition 
LOGICAL :: LMOUNT          ! Logical switch to activate the no-slip condition only on a moutain
INTEGER :: NSTART         ! Grid point number (in the X-direction) from which the no-slip condition is applied, when LMOUNT = .FALSE.  
REAL  :: XHSTART         ! Height abobe  which the no-slip condition is applied, when LMOUNT = .TRUE.  
REAL, DIMENSION(:,:),POINTER :: XDRAG =>NULL() ! Array defining where the no-slip condition is applied (-1/1)
!where XDRAG = 1 => free slip
! where XDRAG = -1 => no slip
!
END TYPE DRAG_t

TYPE(DRAG_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: DRAG_MODEL

LOGICAL, POINTER :: LDRAG=>NULL()
LOGICAL, POINTER :: LMOUNT=>NULL()
INTEGER, POINTER :: NSTART=>NULL()
REAL, POINTER :: XHSTART=>NULL()
REAL,DIMENSION(:,:), POINTER :: XDRAG =>NULL()

CONTAINS

SUBROUTINE DRAG_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
DRAG_MODEL(KFROM)%XDRAG=>XDRAG
!
! Current model is set to model KTO
LDRAG=>DRAG_MODEL(KTO)%LDRAG
LMOUNT=>DRAG_MODEL(KTO)%LMOUNT
NSTART=>DRAG_MODEL(KTO)%NSTART
XHSTART=>DRAG_MODEL(KTO)%XHSTART
XDRAG=>DRAG_MODEL(KTO)%XDRAG

END SUBROUTINE DRAG_GOTO_MODEL

END MODULE MODD_DRAG_n
