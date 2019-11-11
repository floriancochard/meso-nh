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
!     ###########################
      MODULE MODD_PARAM_KAFR_n
!     ###########################
!
!!****  *MODN_PARAM_KAFR* - Declaration of convection control constants 
!!
!!    PURPOSE
!!    -------
!      The purpose of this declarative module is to declare some control
!      constants (call interval, ice) in the deep convection parameterization.  
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_PARAM_KAFRn)
!!          
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96       
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE PARAM_KAFR_t
!
  REAL    :: XDTCONV   ! Interval of time between two calls of
                           ! convection scheme
  LOGICAL :: LINIDCONV ! logical switch to initialize or read in the FM
                           ! file the deep convection tendencies, surface 
                           ! precipitation flux, cumulated precipitation during
                           ! a convective period, large scale compensating 
                           ! subsidence and temporal informations
  LOGICAL :: LREFRESH_ALL ! logical to refresh all convective columns at
                              ! each call of the convection scheme.
  LOGICAL :: LCHTRANS  ! to activate tranport of passive tracers by the
                           ! convective scheme.
  INTEGER :: NICE      ! flag for ice ( 1 = yes, 0 = no ice ) 
  LOGICAL :: LDOWN     ! activate or not convective downdrafts
  LOGICAL :: LSETTADJ  ! logical to specify user defined convective
                           ! adjustment time
  REAL :: XTADJD    ! deep convective adjustemnt time if LSETTADJ=TRUE
  REAL :: XTADJS    ! shallow convective adjustemnt time if LSETTADJ=TRUE
  LOGICAL :: LDIAGCONV ! to activate convection diagnostics 
                           ! (mass flux, precipitation flux, CAPE, cloud levels)
  INTEGER :: NENSM     ! number of additional convective ensemble members for deep
                           ! standard   : NENSM=0 ( 1 deep + 1 shallow)
                           !            : NENSM=3 ( 3 deep + 1 shallow - for SCM studies
                           !                        or very smooth fields - but is slow ! )
!
END TYPE PARAM_KAFR_t

TYPE(PARAM_KAFR_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PARAM_KAFR_MODEL

REAL, POINTER :: XDTCONV=>NULL()
LOGICAL, POINTER :: LINIDCONV=>NULL()
LOGICAL, POINTER :: LREFRESH_ALL=>NULL()
LOGICAL, POINTER :: LCHTRANS=>NULL()
INTEGER, POINTER :: NICE=>NULL()
LOGICAL, POINTER :: LDOWN=>NULL()
LOGICAL, POINTER :: LSETTADJ=>NULL()
REAL, POINTER :: XTADJD=>NULL()
REAL, POINTER :: XTADJS=>NULL()
LOGICAL, POINTER :: LDIAGCONV=>NULL()
INTEGER, POINTER :: NENSM=>NULL()

CONTAINS

SUBROUTINE PARAM_KAFR_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
XDTCONV=>PARAM_KAFR_MODEL(KTO)%XDTCONV
LINIDCONV=>PARAM_KAFR_MODEL(KTO)%LINIDCONV
LREFRESH_ALL=>PARAM_KAFR_MODEL(KTO)%LREFRESH_ALL
LCHTRANS=>PARAM_KAFR_MODEL(KTO)%LCHTRANS
NICE=>PARAM_KAFR_MODEL(KTO)%NICE
LDOWN=>PARAM_KAFR_MODEL(KTO)%LDOWN
LSETTADJ=>PARAM_KAFR_MODEL(KTO)%LSETTADJ
XTADJD=>PARAM_KAFR_MODEL(KTO)%XTADJD
XTADJS=>PARAM_KAFR_MODEL(KTO)%XTADJS
LDIAGCONV=>PARAM_KAFR_MODEL(KTO)%LDIAGCONV
NENSM=>PARAM_KAFR_MODEL(KTO)%NENSM

END SUBROUTINE PARAM_KAFR_GOTO_MODEL

END MODULE MODD_PARAM_KAFR_n
