!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 13:37:10
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_NEST_PGD_n
!     ####################
!
!!****  *MODD_NEST_PGD$n* - declaration of mask for gridnesting of pgds
!!
!!    PURPOSE
!!    -------
!!
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_NEST_PGD$n)
!!       
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/09/95
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE NEST_PGD_t
!
  INTEGER, DIMENSION(:,:,:), POINTER :: NNESTMASK=>NULL()
                                                ! mask relatively to model $n
                                                ! grid, to locate where its
                                                ! sons are located:
                                                ! 1st and 2nd dimensions: i, j
                                                ! 3rd dimension: model index
  REAL, DIMENSION(:,:,:,:), POINTER :: XNESTFIELD=>NULL()
                                                ! idem except for physical field
                                                ! (not mask)
                                                ! 4th dimension: equal to 3rd
                                                ! dimension of pgd field if any
                                                ! (COVER fields),
                                                ! equal to 1 in other case (ZS)
!
  INTEGER, DIMENSION(:),    POINTER :: NSON=>NULL()  ! model index corresponding to
                                                ! 3rd dimension of NNESTMASK
!
END TYPE NEST_PGD_t

TYPE(NEST_PGD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: NEST_PGD_MODEL

INTEGER, DIMENSION(:,:,:), POINTER :: NNESTMASK=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XNESTFIELD=>NULL()
INTEGER, DIMENSION(:),     POINTER :: NSON=>NULL()

CONTAINS

SUBROUTINE NEST_PGD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
NEST_PGD_MODEL(KFROM)%NNESTMASK=>NNESTMASK
NEST_PGD_MODEL(KFROM)%XNESTFIELD=>XNESTFIELD
NEST_PGD_MODEL(KFROM)%NSON=>NSON
!
! Current model is set to model KTO
NNESTMASK=>NEST_PGD_MODEL(KTO)%NNESTMASK
XNESTFIELD=>NEST_PGD_MODEL(KTO)%XNESTFIELD
NSON=>NEST_PGD_MODEL(KTO)%NSON

END SUBROUTINE NEST_PGD_GOTO_MODEL

END MODULE MODD_NEST_PGD_n
