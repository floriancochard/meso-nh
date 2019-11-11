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
!     ############################
      MODULE MODD_STATION_n
!     ############################
!
!!****  *MODD_STATION* - declaration of stations
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the different stations types.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!       
!!    AUTHOR
!!    ------
!!	P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/01/02
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_TYPE_STATION
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE STATION_t
!
!-------------------------------------------------------------------------------------------
!
  LOGICAL                          :: LSTATION    ! flag to use stations
  INTEGER                          :: NUMBSTAT    ! number of stations
  LOGICAL                          :: LSTATLAT    ! positioning in lat/lon
!
  TYPE(STATION) :: TSTATION ! characteristics and records of a station
!
END TYPE STATION_t

TYPE(STATION_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: STATION_MODEL

LOGICAL, POINTER :: LSTATION=>NULL()
INTEGER, POINTER :: NUMBSTAT=>NULL()
LOGICAL, POINTER :: LSTATLAT=>NULL()
TYPE(STATION), POINTER :: TSTATION=>NULL()

CONTAINS

SUBROUTINE STATION_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
LSTATION=>STATION_MODEL(KTO)%LSTATION
NUMBSTAT=>STATION_MODEL(KTO)%NUMBSTAT
LSTATLAT=>STATION_MODEL(KTO)%LSTATLAT
TSTATION=>STATION_MODEL(KTO)%TSTATION

END SUBROUTINE STATION_GOTO_MODEL

END MODULE MODD_STATION_n
