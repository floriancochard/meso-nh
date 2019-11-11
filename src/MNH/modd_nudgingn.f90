!MNH_LIC Copyright 2006-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_NUDGING_n
!     ###################
!
!!****  *MODD_NUDGING$n* - Variables for nudging towards Large Scale fields
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       15/05/06
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE NUDGING_t
  LOGICAL            :: LNUDGING   ! Logical for nudging term
  REAL               :: XTNUDGING  ! Time scale for nudging
!
END TYPE NUDGING_t

TYPE(NUDGING_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: NUDGING_MODEL

LOGICAL, POINTER :: LNUDGING  => NULL()
REAL,    POINTER :: XTNUDGING => NULL()

CONTAINS

SUBROUTINE NUDGING_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
LNUDGING=>NUDGING_MODEL(KTO)%LNUDGING
XTNUDGING=>NUDGING_MODEL(KTO)%XTNUDGING

END SUBROUTINE NUDGING_GOTO_MODEL

END MODULE MODD_NUDGING_n
