!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_PAST_FIELD_n
!     ###################
!
!!****  *MODD_PAST_FIELD$n* - declaration of prognostic variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     prognostic variables. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PAST_FIELDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2013                       
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

REAL, DIMENSION(:,:,:), POINTER :: XUM=>NULL(),XVM=>NULL(),XWM=>NULL() ! U,V,W  at time t-dt
REAL, DIMENSION(:,:,:), POINTER :: XDUM=>NULL(),XDVM=>NULL(),XDWM=>NULL()

CONTAINS

SUBROUTINE PAST_FIELD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! ! Save current state for allocated arrays
! !FIELD_MODEL(KFROM)%XUM=>XUM !Done in FIELDLIST_GOTO_MODEL
! !FIELD_MODEL(KFROM)%XVM=>XVM !Done in FIELDLIST_GOTO_MODEL
! !FIELD_MODEL(KFROM)%XWM=>XWM !Done in FIELDLIST_GOTO_MODEL
! !FIELD_MODEL(KFROM)%XDUM=>XDUM !Done in FIELDLIST_GOTO_MODEL
! !FIELD_MODEL(KFROM)%XDVM=>XDVM !Done in FIELDLIST_GOTO_MODEL
! !FIELD_MODEL(KFROM)%XDWM=>XDWM !Done in FIELDLIST_GOTO_MODEL
! !
! ! Current model is set to model KTO
! !XUM=>FIELD_MODEL(KTO)%XUM !Done in FIELDLIST_GOTO_MODEL
! !XVM=>FIELD_MODEL(KTO)%XVM !Done in FIELDLIST_GOTO_MODEL
! !XWM=>FIELD_MODEL(KTO)%XWM !Done in FIELDLIST_GOTO_MODEL
! !XDUM=>FIELD_MODEL(KTO)%XDUM !Done in FIELDLIST_GOTO_MODEL
! !XDVM=>FIELD_MODEL(KTO)%XDVM !Done in FIELDLIST_GOTO_MODEL
! !XDWM=>FIELD_MODEL(KTO)%XDWM !Done in FIELDLIST_GOTO_MODEL

END SUBROUTINE PAST_FIELD_GOTO_MODEL

END MODULE MODD_PAST_FIELD_n
