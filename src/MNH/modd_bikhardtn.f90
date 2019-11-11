!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 14:09:22
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_BIKHARDT_n
!     ####################
!
!!****  *MODD_BIKHARDT$n* - declaration of coefficients for Bikhardt interpolation
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the coefficients
!      for Bikhardt interpolation. $n refers to the INNER model.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_BIKHARDTn)
!!          
!!    AUTHOR
!!    ------
!!	Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/06/96
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE BIKHARDT_t
!
  REAL, DIMENSION(:), POINTER :: XBMX1=>NULL(),XBMX2=>NULL(), &
                                 XBMX3=>NULL(),XBMX4=>NULL() ! Mass points in X-direc.
  REAL, DIMENSION(:), POINTER :: XBMY1=>NULL(),XBMY2=>NULL(), &
                                 XBMY3=>NULL(),XBMY4=>NULL() ! Mass points in Y-direc.
  REAL, DIMENSION(:), POINTER :: XBFX1=>NULL(),XBFX2=>NULL(), &
                                 XBFX3=>NULL(),XBFX4=>NULL() ! Flux points in X-direc.
  REAL, DIMENSION(:), POINTER :: XBFY1=>NULL(),XBFY2=>NULL(), &
                                 XBFY3=>NULL(),XBFY4=>NULL() ! Flux points in Y-direc.
!
END TYPE BIKHARDT_t

TYPE(BIKHARDT_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: BIKHARDT_MODEL

REAL, DIMENSION(:), POINTER :: XBMX1=>NULL(),XBMX2=>NULL(),XBMX3=>NULL(),XBMX4=>NULL()
REAL, DIMENSION(:), POINTER :: XBMY1=>NULL(),XBMY2=>NULL(),XBMY3=>NULL(),XBMY4=>NULL()
REAL, DIMENSION(:), POINTER :: XBFX1=>NULL(),XBFX2=>NULL(),XBFX3=>NULL(),XBFX4=>NULL()
REAL, DIMENSION(:), POINTER :: XBFY1=>NULL(),XBFY2=>NULL(),XBFY3=>NULL(),XBFY4=>NULL()

CONTAINS

SUBROUTINE BIKHARDT_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
BIKHARDT_MODEL(KFROM)%XBMX1=>XBMX1
BIKHARDT_MODEL(KFROM)%XBMX2=>XBMX2
BIKHARDT_MODEL(KFROM)%XBMX3=>XBMX3
BIKHARDT_MODEL(KFROM)%XBMX4=>XBMX4
BIKHARDT_MODEL(KFROM)%XBMY1=>XBMY1
BIKHARDT_MODEL(KFROM)%XBMY2=>XBMY2
BIKHARDT_MODEL(KFROM)%XBMY3=>XBMY3
BIKHARDT_MODEL(KFROM)%XBMY4=>XBMY4
BIKHARDT_MODEL(KFROM)%XBFX1=>XBFX1
BIKHARDT_MODEL(KFROM)%XBFX2=>XBFX2
BIKHARDT_MODEL(KFROM)%XBFX3=>XBFX3
BIKHARDT_MODEL(KFROM)%XBFX4=>XBFX4
BIKHARDT_MODEL(KFROM)%XBFY1=>XBFY1
BIKHARDT_MODEL(KFROM)%XBFY2=>XBFY2
BIKHARDT_MODEL(KFROM)%XBFY3=>XBFY3
BIKHARDT_MODEL(KFROM)%XBFY4=>XBFY4
!
! Current model is set to model KTO
XBMX1=>BIKHARDT_MODEL(KTO)%XBMX1
XBMX2=>BIKHARDT_MODEL(KTO)%XBMX2
XBMX3=>BIKHARDT_MODEL(KTO)%XBMX3
XBMX4=>BIKHARDT_MODEL(KTO)%XBMX4
XBMY1=>BIKHARDT_MODEL(KTO)%XBMY1
XBMY2=>BIKHARDT_MODEL(KTO)%XBMY2
XBMY3=>BIKHARDT_MODEL(KTO)%XBMY3
XBMY4=>BIKHARDT_MODEL(KTO)%XBMY4
XBFX1=>BIKHARDT_MODEL(KTO)%XBFX1
XBFX2=>BIKHARDT_MODEL(KTO)%XBFX2
XBFX3=>BIKHARDT_MODEL(KTO)%XBFX3
XBFX4=>BIKHARDT_MODEL(KTO)%XBFX4
XBFY1=>BIKHARDT_MODEL(KTO)%XBFY1
XBFY2=>BIKHARDT_MODEL(KTO)%XBFY2
XBFY3=>BIKHARDT_MODEL(KTO)%XBFY3
XBFY4=>BIKHARDT_MODEL(KTO)%XBFY4

END SUBROUTINE BIKHARDT_GOTO_MODEL

END MODULE MODD_BIKHARDT_n
