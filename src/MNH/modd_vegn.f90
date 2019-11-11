!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 13:54:55
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_VEG_n
!     ##################
!
!!****  *MODD_VEG$n* - declaration of veg variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     describing the veg. 
!    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_VEGn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	P. Aumond   *Meteo France*
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX

IMPLICIT NONE

TYPE VEG_t
!
  REAL, DIMENSION(:,:), POINTER :: XLAI_PGD=>NULL()   ! orography
  REAL, DIMENSION(:,:), POINTER :: XH_TREE_PGD=>NULL()   ! height z 

END TYPE VEG_t

TYPE(VEG_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: VEG_MODEL


REAL, DIMENSION(:,:), POINTER :: XLAI_PGD=>NULL()
REAL, DIMENSION(:,:), POINTER :: XH_TREE_PGD=>NULL()

CONTAINS

SUBROUTINE VEG_GOTO_MODEL(KFROM, KTO)

INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
VEG_MODEL(KFROM)%XLAI_PGD=>XLAI_PGD
VEG_MODEL(KFROM)%XH_TREE_PGD=>XH_TREE_PGD
!
! Current model is set to model KTO
XLAI_PGD=>VEG_MODEL(KTO)%XLAI_PGD
XH_TREE_PGD=>VEG_MODEL(KTO)%XH_TREE_PGD

END SUBROUTINE VEG_GOTO_MODEL

END MODULE MODD_VEG_n
