!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 12:42:44
!-----------------------------------------------------------------
!!    ########################
      MODULE MODD_SUB_CH_MONITOR_n
!!    ########################
!!
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE SUB_CH_MONITOR_t
!!
  INTEGER :: ISVECNMASK,ISVECNPT
  INTEGER, DIMENSION(:,:), POINTER :: ISVECMASK=>NULL()
  INTEGER          :: ISTCOUNT = 0  ! counter for JVALUE update
END TYPE SUB_CH_MONITOR_t

TYPE(SUB_CH_MONITOR_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: SUB_CH_MONITOR_MODEL

INTEGER, POINTER :: ISVECNMASK=>NULL(),ISVECNPT=>NULL()
INTEGER, DIMENSION(:,:), POINTER :: ISVECMASK=>NULL()
INTEGER, POINTER :: ISTCOUNT=>NULL()

CONTAINS

SUBROUTINE SUB_CH_MONITOR_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
SUB_CH_MONITOR_MODEL(KFROM)%ISVECMASK=>ISVECMASK
!
! Current model is set to model KTO
ISVECNMASK=>SUB_CH_MONITOR_MODEL(KTO)%ISVECNMASK
ISVECNPT=>SUB_CH_MONITOR_MODEL(KTO)%ISVECNPT
ISVECMASK=>SUB_CH_MONITOR_MODEL(KTO)%ISVECMASK
ISTCOUNT=>SUB_CH_MONITOR_MODEL(KTO)%ISTCOUNT

END SUBROUTINE SUB_CH_MONITOR_GOTO_MODEL

END MODULE MODD_SUB_CH_MONITOR_n
