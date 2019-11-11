!MNH_LIC Copyright 2007-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ########################
      MODULE MODD_CH_PH_n
!!    ########################
!!
!!*** *MODD_CH_PH$n*
!!
!!    PURPOSE
!!    -------
!       This module contains the pH field of the cloud water of the rainwater
!!
!!**  AUTHOR
!!    ------
!!    M. Leriche      *Laboratoire d'Aerologie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 01/06/07
!!       P. Tulet      Nov 2014 accumulated moles of aqueous species that fall at the surface   
!!       P. Tulet & M. Leriche Nov 2015 add pH in rain at the surface
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
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

TYPE CH_PH_t
!

!  REAL, POINTER, DIMENSION(:,:,:) :: XPHC => NULL() ! cloud
!  REAL, POINTER, DIMENSION(:,:,:) :: XPHR => NULL() ! rain
  REAL, POINTER, DIMENSION(:,:,:) :: XACPRAQ => NULL() ! sum of aqueous chemical species fall at the surface by rain
                                                       ! in moles i / m2 (ratio with XACPRR for concentration
  REAL, POINTER, DIMENSION(:,:) :: XACPHR => NULL()    !  mean PH in accumulated surface rain
!
!-----------------------------------------------------------------------------
END TYPE CH_PH_t

TYPE(CH_PH_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_PH_MODEL

REAL, POINTER, DIMENSION(:,:,:) :: XPHC=>NULL()
REAL, POINTER, DIMENSION(:,:,:) :: XPHR=>NULL()
REAL, POINTER, DIMENSION(:,:) :: XACPHR=>NULL()
REAL, POINTER, DIMENSION(:,:,:) :: XACPRAQ=>NULL()

CONTAINS

SUBROUTINE CH_PH_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!CH_PH_MODEL(KFROM)%XPHC=>XPHC !Done in FIELDLIST_GOTO_MODEL
!CH_PH_MODEL(KFROM)%XPHR=>XPHR !Done in FIELDLIST_GOTO_MODEL
CH_PH_MODEL(KFROM)%XACPHR=>XACPHR
CH_PH_MODEL(KFROM)%XACPRAQ=>XACPRAQ
!
! Current model is set to model KTO
!XPHC=>CH_PH_MODEL(KTO)%XPHC !Done in FIELDLIST_GOTO_MODEL
!XPHR=>CH_PH_MODEL(KTO)%XPHR !Done in FIELDLIST_GOTO_MODEL
XACPHR=>CH_PH_MODEL(KTO)%XACPHR
XACPRAQ=>CH_PH_MODEL(KTO)%XACPRAQ

END SUBROUTINE CH_PH_GOTO_MODEL

END MODULE MODD_CH_PH_n
