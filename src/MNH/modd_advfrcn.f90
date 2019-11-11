!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_ADVFRC_n
!     #####################
!
!!***  *MODD_ADVFRC -  Declarative module for the forcing fields of the 2D model
!!
!!    PURPOSE
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_FRC)
!!      
!!    AUTHOR
!!    ------
!!	    P. Peyrille from modd_frc
!!
!!    MODIFICATIONS
!!    -------------
!!      16/10/10 M.Tomasini Grid-nesting
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!
IMPLICIT NONE
!
!*            fields for FORCING
!             ------------------
!
TYPE ADVFRC_t
!
REAL, DIMENSION(:,:,:,:), POINTER   :: XDTHFRC=>NULL() ! Dth frc (I,J,K,TIME)
REAL, DIMENSION(:,:,:,:), POINTER   :: XDRVFRC=>NULL() ! Drv frc 
TYPE (DATE_TIME), DIMENSION(:), POINTER :: TDTADVFRC=>NULL() ! Date and time of the advecting and/or relaxation forcing

!
END TYPE ADVFRC_t
!
TYPE(ADVFRC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: ADVFRC_MODEL
!
INTEGER,          SAVE              :: NADVFRC  ! number of forcing profiles
INTEGER,          SAVE                  :: NPRESSLEV_ADV   ! number of frc levels
REAL, DIMENSION(:,:,:,:), POINTER   :: XDTHFRC=>NULL() 
REAL, DIMENSION(:,:,:,:), POINTER   :: XDRVFRC=>NULL()
TYPE (DATE_TIME), DIMENSION(:), POINTER :: TDTADVFRC=>NULL()
!
!
CONTAINS
!
SUBROUTINE ADVFRC_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
ADVFRC_MODEL(KFROM)%XDTHFRC=>XDTHFRC
ADVFRC_MODEL(KFROM)%XDRVFRC=>XDRVFRC
ADVFRC_MODEL(KFROM)%TDTADVFRC=>TDTADVFRC
!
! Current model is set to model KTO
XDTHFRC=>ADVFRC_MODEL(KTO)%XDTHFRC
XDRVFRC=>ADVFRC_MODEL(KTO)%XDRVFRC
TDTADVFRC=>ADVFRC_MODEL(KTO)%TDTADVFRC

!
END SUBROUTINE ADVFRC_GOTO_MODEL
!
!
END MODULE MODD_ADVFRC_n
