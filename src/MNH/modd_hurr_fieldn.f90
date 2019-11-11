!MNH_LIC Copyright 2001-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############################
      MODULE MODD_HURR_FIELD_n
!     ############################
!
!!****  *MODD_HURR_FIELD$n* - declaration of filtered fields 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     fields obtained during the different steps of the filtering. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!       This routine is not available in Meso-NH documentation yet.     
!!
!!    AUTHOR
!!    ------
!!
!!       O. Nuissier         * L.A. *
!!       R. Rogers           * NOAA/AOML/HRD (Hurricane Research Division) *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/12/01                      
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

REAL, DIMENSION(:,:,:), POINTER :: XUTOT=>NULL(),XVTOT=>NULL(),XTTOT=>NULL()       ! Total fields
REAL, DIMENSION(:,:),   POINTER :: XPTOT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XUENV=>NULL(),XVENV=>NULL(),XTENV=>NULL()       !  Environmental fields
REAL, DIMENSION(:,:),   POINTER :: XPENV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XUBASIC=>NULL(),XVBASIC=>NULL(),XTBASIC=>NULL() ! Basic fields
REAL, DIMENSION(:,:,:), POINTER :: XPBASIC=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XVTDIS=>NULL()                                  ! Total disturbance tangential wind
REAL, DIMENSION(:,:,:), POINTER :: XQTOT=>NULL(),XQENV=>NULL(),XQBASIC=>NULL()

CONTAINS

SUBROUTINE HURR_FIELD_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
END SUBROUTINE HURR_FIELD_GOTO_MODEL

END MODULE MODD_HURR_FIELD_n

