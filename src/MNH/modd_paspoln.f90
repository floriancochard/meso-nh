!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
! $Source$
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_PASPOL_n
!     ####################
!-------------------------------------------------------------------------------
!***	MODD_PASPOL_n  Declaration of passive pollutant global variables  
!
!!    AUTHOR
!!    ------
!
!	           : Michel Bouzom, DP/SERV/ENV
!	Creation   : 09.10.2001
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE PASPOL_t
! 
REAL, DIMENSION(:,:,:,:), POINTER :: XATC=>NULL() ! 
     ! Atmospheric Transfer Coefficient
! 
END TYPE PASPOL_t

TYPE(PASPOL_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PASPOL_MODEL
REAL, DIMENSION(:,:,:,:), POINTER :: XATC=>NULL() ! 

CONTAINS

SUBROUTINE PASPOL_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
PASPOL_MODEL(KFROM)%XATC=>XATC          
!
! Current model is set to model KTO
XATC=>PASPOL_MODEL(KTO)%XATC 
!
END SUBROUTINE PASPOL_GOTO_MODEL

END MODULE MODD_PASPOL_n
