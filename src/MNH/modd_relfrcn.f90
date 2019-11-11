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
      MODULE MODD_RELFRC_n
!     #####################
!
!!***  *MODD_RELFRC -  Declarative module for the forcing fields of the 2D model
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
TYPE RELFRC_t
!
TYPE (DATE_TIME), DIMENSION(:), POINTER :: TDTRELFRC=>NULL() ! Date and time of the advecting and/or relaxation forcing
REAL, DIMENSION(:,:,:,:), POINTER   :: XTHREL=>NULL() ! Dth frc 
REAL, DIMENSION(:,:,:,:), POINTER   :: XRVREL=>NULL() ! Drv frc 

!
END TYPE RELFRC_t
!
TYPE(RELFRC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: RELFRC_MODEL
!
INTEGER,          SAVE              :: NRELFRC  ! number of forcing profiles
INTEGER,          SAVE              :: NPRESSLEV_REL   ! number of rel levels
REAL, DIMENSION(:,:,:,:), POINTER   :: XTHREL=>NULL() 
REAL, DIMENSION(:,:,:,:), POINTER   :: XRVREL=>NULL()
TYPE (DATE_TIME), DIMENSION(:), POINTER :: TDTRELFRC=>NULL()
!
!
CONTAINS
!
SUBROUTINE RELFRC_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
RELFRC_MODEL(KFROM)%TDTRELFRC=>TDTRELFRC
RELFRC_MODEL(KFROM)%XTHREL=>XTHREL
RELFRC_MODEL(KFROM)%XRVREL=>XRVREL
!
! Current model is set to model KTO
TDTRELFRC=>RELFRC_MODEL(KTO)%TDTRELFRC
XTHREL=>RELFRC_MODEL(KTO)%XTHREL
XRVREL=>RELFRC_MODEL(KTO)%XRVREL

!
END SUBROUTINE RELFRC_GOTO_MODEL
!
!
END MODULE MODD_RELFRC_n
