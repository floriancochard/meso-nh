!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_PARAM_n
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
!!      V. Masson      01/04  externalised surface: remove CGROUND, CSEA_FLUX, add CSURF
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE PARAM_t
!
  CHARACTER (LEN=4) :: CTURB    ! Kind of turbulence parameterization
                                   ! 'NONE' if no parameterization 
  CHARACTER (LEN=4) :: CRAD     ! Kind of radiation parameterization 
                                   ! 'NONE' if no parameterization 
  CHARACTER (LEN=4) :: CCLOUD   ! Kind of cloud parameterization 
                                   ! 'NONE' if no parameterization 
  CHARACTER (LEN=4) :: CDCONV   ! Kind of deep convection
                                   ! 'NONE' if no parameterization
  CHARACTER (LEN=4) :: CSCONV   ! Kind of shallow convection
                                   ! 'NONE' if no parameterization
!  CHARACTER (LEN=4) :: CSURF    ! Kind of surface processes parameterization
!                                   ! 'NONE' if no parameterization
  CHARACTER (LEN=4) :: CELEC    ! Kind of  atmospheric electricity scheme
  CHARACTER (LEN=4) :: CACTCCN  ! Kind of CCN activation scheme
!
END TYPE PARAM_t

TYPE(PARAM_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PARAM_MODEL

CHARACTER (LEN=4), POINTER :: CTURB=>NULL()
CHARACTER (LEN=4), POINTER :: CRAD=>NULL()
CHARACTER (LEN=4), POINTER :: CCLOUD=>NULL()
CHARACTER (LEN=4), POINTER :: CDCONV=>NULL()
CHARACTER (LEN=4), POINTER :: CSCONV=>NULL()
CHARACTER (LEN=4), POINTER :: CSURF=>NULL()
CHARACTER (LEN=4), POINTER :: CELEC=>NULL()
CHARACTER (LEN=4), POINTER :: CACTCCN=>NULL()

CONTAINS

SUBROUTINE PARAM_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
CTURB=>PARAM_MODEL(KTO)%CTURB
CRAD=>PARAM_MODEL(KTO)%CRAD
CCLOUD=>PARAM_MODEL(KTO)%CCLOUD
CDCONV=>PARAM_MODEL(KTO)%CDCONV
CSCONV=>PARAM_MODEL(KTO)%CSCONV
!CSURF=>PARAM_MODEL(KTO)%CSURF !Done in FIELDLIST_GOTO_MODEL
CELEC=>PARAM_MODEL(KTO)%CELEC
CACTCCN=>PARAM_MODEL(KTO)%CACTCCN

END SUBROUTINE PARAM_GOTO_MODEL

END MODULE MODD_PARAM_n
