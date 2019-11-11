!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 14:04:12
!-----------------------------------------------------------------
!     ######################
      MODULE MODD_CURVCOR_n
!     ######################
!
!!****  *MODD_CURVCOR$n* - declaration of curvature and coriolis variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     involved in curvature and coriolis terms. 
!    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_CURVCORn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    16/05/94                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CURVCOR_t
  REAL, DIMENSION(:,:), POINTER :: XCORIOX=>NULL(),XCORIOY=>NULL(),XCORIOZ=>NULL()
                 ! -f(*) sin Gamma, f(*) cos Gamma, f
                 !
  REAL, DIMENSION(:,:), POINTER :: XCURVX=>NULL(),XCURVY=>NULL()     
                 ! curvature variables
                 !
END TYPE CURVCOR_t

TYPE(CURVCOR_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CURVCOR_MODEL

REAL, DIMENSION(:,:), POINTER :: XCORIOX=>NULL(),XCORIOY=>NULL(),XCORIOZ=>NULL()
REAL, DIMENSION(:,:), POINTER :: XCURVX=>NULL(),XCURVY=>NULL()

CONTAINS

SUBROUTINE CURVCOR_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CURVCOR_MODEL(KFROM)%XCORIOX=>XCORIOX
CURVCOR_MODEL(KFROM)%XCORIOY=>XCORIOY
CURVCOR_MODEL(KFROM)%XCORIOZ=>XCORIOZ
CURVCOR_MODEL(KFROM)%XCURVX=>XCURVX
CURVCOR_MODEL(KFROM)%XCURVY=>XCURVY
!
! Current model is set to model KTO
XCORIOX=>CURVCOR_MODEL(KTO)%XCORIOX
XCORIOY=>CURVCOR_MODEL(KTO)%XCORIOY
XCORIOZ=>CURVCOR_MODEL(KTO)%XCORIOZ
XCURVX=>CURVCOR_MODEL(KTO)%XCURVX
XCURVY=>CURVCOR_MODEL(KTO)%XCURVY

END SUBROUTINE CURVCOR_GOTO_MODEL

END MODULE MODD_CURVCOR_n
