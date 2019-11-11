!MNH_LIC Copyright 2016-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODD_CH_BUDGET_n
!     #######################
!
!!****  *MODD_CH_BUDGET_n* - declaration of parameters for chemical budget diagnostic
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the
!     variables that may be used to compute the detailed chemical budget 
!     as diagnostic (DIAG step)
!     Each production and destruction terms of a given species are computed
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!
!!    AUTHOR
!!    ------
!!	    F. Brosse *Laboratoire d'Aerologie UPS-CNRS*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    October 2016
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!      
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE TREAC_BUDGET
  REAL   ,           DIMENSION(:,:,:,:), POINTER :: XB_REAC => NULL()
  INTEGER,           DIMENSION(:)      , POINTER :: NB_REAC => NULL()
END TYPE TREAC_BUDGET

TYPE TCH_BUDGET_t

  CHARACTER(LEN=32),  DIMENSION(:), POINTER :: CNAMES_BUDGET => NULL()
  INTEGER,            DIMENSION(:), POINTER :: NSPEC_BUDGET => NULL()
  INTEGER                                   :: NEQ_BUDGET
  TYPE(TREAC_BUDGET), DIMENSION(:), POINTER :: XTCHEM => NULL()

END TYPE TCH_BUDGET_t

TYPE(TCH_BUDGET_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_BUDGET_MODEL

CHARACTER(LEN=32),  DIMENSION(:), POINTER:: CNAMES_BUDGET=>NULL()
INTEGER,            DIMENSION(:), POINTER:: NSPEC_BUDGET=>NULL()
INTEGER, POINTER :: NEQ_BUDGET=>NULL()
TYPE(TREAC_BUDGET), DIMENSION(:), POINTER :: XTCHEM=>NULL()

CONTAINS

SUBROUTINE CH_BUDGET_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CH_BUDGET_MODEL(KFROM)%CNAMES_BUDGET=>CNAMES_BUDGET
CH_BUDGET_MODEL(KFROM)%NSPEC_BUDGET=>NSPEC_BUDGET
CH_BUDGET_MODEL(KFROM)%XTCHEM=>XTCHEM
!
! Current model is set to model KTO
CNAMES_BUDGET=>CH_BUDGET_MODEL(KTO)%CNAMES_BUDGET
NSPEC_BUDGET=>CH_BUDGET_MODEL(KTO)%NSPEC_BUDGET
NEQ_BUDGET=>CH_BUDGET_MODEL(KTO)%NEQ_BUDGET
XTCHEM=>CH_BUDGET_MODEL(KTO)%XTCHEM

END SUBROUTINE CH_BUDGET_GOTO_MODEL

END MODULE MODD_CH_BUDGET_n
