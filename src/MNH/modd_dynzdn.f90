!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODD_DYNZD_n
!     #################
!
!!****  *MODD_DYN$n* - declaration of dynamic control variables 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the dynamic
!     control variables.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DYNn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/05/94                      
!!      Modifications 16/11/94   (Lafore+Pinty)  For NUM_DIFF
!!      Modifications 06/01/95   (Lafore)        For LSTEADY_DMASS
!!      Modifications 28/07/96   (Masson)        Supress LSTEADY_DMASS
!!      Modifications 15/03/98   (Stein)         Add LHO_RELAX for each variables 
!!      Modifications 22/01/01   (Gazen)         Add LHORELAX_SVC2R2, _SVCHEM, _SVLG
!!      Modifications 29/11/02   (Pinty)         Add  LHORELAX_SVC1R3, _SVELEC
!!      Modifications 03/11/04   (Zängl)         Add fields for truly horizontal diffusion
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
USE MODD_DYNZD
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!
IMPLICIT NONE

TYPE DYNZD_t
!
! Additional variables needed for truly horizontal diffusion (G. Zängl)
! Interpolation coefficients
  TYPE(TYPE_ZDIFFU_HALO2) :: XZDIFFU_HALO2
!
END TYPE DYNZD_t

TYPE(DYNZD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: DYNZD_MODEL

TYPE(TYPE_ZDIFFU_HALO2), POINTER :: XZDIFFU_HALO2=>NULL()

CONTAINS

SUBROUTINE DYNZD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
XZDIFFU_HALO2=>DYNZD_MODEL(KTO)%XZDIFFU_HALO2

END SUBROUTINE DYNZD_GOTO_MODEL

END MODULE MODD_DYNZD_n
