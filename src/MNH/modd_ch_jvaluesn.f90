!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 14:16:58
!-----------------------------------------------------------------
!!    ########################
      MODULE MODD_CH_JVALUES_n
!!    ########################
!!
!!*** *MODD_CH_JVALUES$n*
!!
!!    PURPOSE
!!    -------
!       This module contains the photolysis rates (JVALUES) calculated
!     by the radiative transfer code TUV39
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre      *Laboratoire d'Aerologie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 05/03/97
!!    14/02/99 (K. Suhre) add surface albedo and column dobson to parameter list
!!    16/02/99 (K. Suhre) add offline option
!!    20/01/01 (C. Mari) 3D interpolation of J values 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_JVALUES_t
!
  REAL, DIMENSION(:,:,:,:), POINTER :: XJVALUES=>NULL()  
                              ! photolysis rates after time interpolation
  LOGICAL :: GSFIRSTCALL = .TRUE. ! flag for initialization on first call
!
!-----------------------------------------------------------------------------
END TYPE CH_JVALUES_t

TYPE(CH_JVALUES_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_JVALUES_MODEL

REAL, POINTER, DIMENSION(:,:,:,:) :: XJVALUES=>NULL()
LOGICAL, POINTER :: GSFIRSTCALL=>NULL()

CONTAINS

SUBROUTINE CH_JVALUES_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CH_JVALUES_MODEL(KFROM)%XJVALUES=>XJVALUES
!
! Current model is set to model KTO
XJVALUES=>CH_JVALUES_MODEL(KTO)%XJVALUES
GSFIRSTCALL=>CH_JVALUES_MODEL(KTO)%GSFIRSTCALL

END SUBROUTINE CH_JVALUES_GOTO_MODEL

END MODULE MODD_CH_JVALUES_n
