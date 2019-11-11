!     ###################
      MODULE MODN_DRAG_n
!     ###################
!
!!****  *MODN_DRAG* - declaration of namelist NAM_DRAG
!!
!!    PURPOSE
!!    -------
!!      The purpose of this module is to specify the namelist NAM_DRAG
!!    (no-slip condition)
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	    J. Colin                  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    October 2011
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DRAG_n, ONLY: &
     LDRAG_n =>LDRAG, &
     LMOUNT_n =>LMOUNT, &
     NSTART_n =>NSTART, &
     XHSTART_n =>XHSTART
!
IMPLICIT NONE
!
LOGICAL,SAVE  :: LDRAG
LOGICAL,SAVE  :: LMOUNT
INTEGER,SAVE  :: NSTART
REAL,SAVE     :: XHSTART
!
NAMELIST/NAM_DRAGn/LDRAG,LMOUNT,XHSTART,NSTART
!
CONTAINS
!
SUBROUTINE INIT_NAM_DRAGn
  LDRAG = LDRAG_n
  LMOUNT = LMOUNT_n
  NSTART = NSTART_n
  XHSTART = XHSTART_n
END SUBROUTINE INIT_NAM_DRAGn

SUBROUTINE UPDATE_NAM_DRAGn
  LDRAG_n = LDRAG
  LMOUNT_n = LMOUNT
  NSTART_n = NSTART
  XHSTART_n = XHSTART
END SUBROUTINE UPDATE_NAM_DRAGn

END MODULE MODN_DRAG_n
