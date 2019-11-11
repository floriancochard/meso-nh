!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 12:31:15
!-----------------------------------------------------------------
!    ########################
     MODULE MODD_SUB_PHYS_PARAM_n  
!    ########################
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE SUB_PHYS_PARAM_t
!
  REAL, DIMENSION(:,:), POINTER :: ZCOSZEN=>NULL()     ! cosine of zenithal angle
  REAL, DIMENSION(:,:), POINTER :: ZSINZEN=>NULL()     ! sine of zenithal angle
  REAL, DIMENSION(:,:), POINTER :: ZAZIMSOL=>NULL()    ! azimuthal solar angle
!
END TYPE SUB_PHYS_PARAM_t

TYPE(SUB_PHYS_PARAM_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: SUB_PHYS_PARAM_MODEL

REAL, DIMENSION(:,:), POINTER :: ZCOSZEN=>NULL()
REAL, DIMENSION(:,:), POINTER :: ZSINZEN=>NULL()
REAL, DIMENSION(:,:), POINTER :: ZAZIMSOL=>NULL()

CONTAINS

SUBROUTINE SUB_PHYS_PARAM_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
SUB_PHYS_PARAM_MODEL(KFROM)%ZCOSZEN=>ZCOSZEN
SUB_PHYS_PARAM_MODEL(KFROM)%ZSINZEN=>ZSINZEN
SUB_PHYS_PARAM_MODEL(KFROM)%ZAZIMSOL=>ZAZIMSOL
!
! Current model is set to model KTO
ZCOSZEN=>SUB_PHYS_PARAM_MODEL(KTO)%ZCOSZEN
ZSINZEN=>SUB_PHYS_PARAM_MODEL(KTO)%ZSINZEN
ZAZIMSOL=>SUB_PHYS_PARAM_MODEL(KTO)%ZAZIMSOL

END SUBROUTINE SUB_PHYS_PARAM_GOTO_MODEL

END MODULE MODD_SUB_PHYS_PARAM_n
