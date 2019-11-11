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
      MODULE MODD_PROFILER_n
!     ############################
!
!!****  *MODD_PROFILER* - declaration of stations
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
USE MODD_TYPE_PROFILER
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE PROFILER_t
!
!-------------------------------------------------------------------------------------------
!
  LOGICAL                          :: LPROFILER    ! flag to use stations
  INTEGER                          :: NUMBPROFILER    ! number of stations
!
  TYPE(PROFILER) :: TPROFILER ! characteristics and records of an aircraft
!
END TYPE PROFILER_t

TYPE(PROFILER_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PROFILER_MODEL

LOGICAL, POINTER :: LPROFILER=>NULL()
INTEGER, POINTER :: NUMBPROFILER=>NULL()
TYPE(PROFILER), POINTER :: TPROFILER=>NULL()

CONTAINS

SUBROUTINE PROFILER_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
LPROFILER=>PROFILER_MODEL(KTO)%LPROFILER
NUMBPROFILER=>PROFILER_MODEL(KTO)%NUMBPROFILER
TPROFILER=>PROFILER_MODEL(KTO)%TPROFILER

END SUBROUTINE PROFILER_GOTO_MODEL

END MODULE MODD_PROFILER_n
