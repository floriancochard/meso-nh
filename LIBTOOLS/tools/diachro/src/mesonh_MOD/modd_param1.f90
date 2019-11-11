!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_PARAM1
!     ####################
!
!!****  *MODD_PARAM$n* - declaration of parameterization and cloud physics variables 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the
!     parameterization and cloud physics variables.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAMn)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/06/94    
!!      E. Richard  01/06/95  add the selctor for the microphysical scheme
!!      P. Bechtold 26/03/96  add the selector for the deep convection
!!      M. Tomasini 11/12/00  add the selector for the fluxes algorithm over water
!!      JP. Pinty   26/11/02  add the selector for the atmospheric electricity scheme
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
CHARACTER (LEN=4),SAVE :: CTURB    ! Kind of turbulence parameterization
                                   ! 'NONE' if no parameterization 
CHARACTER (LEN=4),SAVE :: CRAD     ! Kind of radiation parameterization 
                                   ! 'NONE' if no parameterization 
CHARACTER (LEN=4),SAVE :: CDRAG    ! Kind of drag parameterization 
                                   ! 'NONE' if no parameterization 
CHARACTER (LEN=4),SAVE :: CCLOUD   ! Kind of cloud parameterization 
                                   ! 'NONE' if no parameterization 
CHARACTER (LEN=4),SAVE :: CDCONV   ! Kind of deep convection
                                   ! 'NONE' if no parameterization
CHARACTER (LEN=4),SAVE :: CELEC    ! Kind of  atmospheric electricity scheme
CHARACTER (LEN=4),SAVE :: CSURF    ! Kind of surface processes parameterization
!
END MODULE MODD_PARAM1
