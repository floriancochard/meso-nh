!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 12:49:50
!-----------------------------------------------------------------
!     ##########################
      MODULE MODD_SECPGD_FIELD_n
!     ##########################
!
!!****  *MODD_SECPGD_FIELD$n
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHOR
!!    ------
!!	V. Masson
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE SECPGD_FIELD_t
  REAL, DIMENSION(:,:), POINTER :: XAOSIP=>NULL(),XAOSIM=>NULL(),XAOSJP=>NULL(),XAOSJM=>NULL()
! Vertical cross section of subgrid mountains in a direction 
! over surface of the grid mesh (used to recompute Z0EFF)
!
  REAL, DIMENSION(:,:), POINTER :: XHO2IP=>NULL(),XHO2IM=>NULL(),XHO2JP=>NULL(),XHO2JM=>NULL()
! half height of the subgrid mountain in a direction (used to recompute Z0EFF)
!
  REAL, DIMENSION(:,:), POINTER :: XAVG_ZS=>NULL()    ! averaged orography
  REAL, DIMENSION(:,:), POINTER :: XSIL_ZS=>NULL()    ! silhouette orography
!
END TYPE SECPGD_FIELD_t

TYPE(SECPGD_FIELD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: SECPGD_FIELD_MODEL

REAL, DIMENSION(:,:), POINTER :: XAOSIP=>NULL(),XAOSIM=>NULL(),XAOSJP=>NULL(),XAOSJM=>NULL()
REAL, DIMENSION(:,:), POINTER :: XHO2IP=>NULL(),XHO2IM=>NULL(),XHO2JP=>NULL(),XHO2JM=>NULL()
REAL, DIMENSION(:,:), POINTER :: XAVG_ZS=>NULL()
REAL, DIMENSION(:,:), POINTER :: XSIL_ZS=>NULL()

CONTAINS

SUBROUTINE SECPGD_FIELD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
SECPGD_FIELD_MODEL(KFROM)%XAOSIP=>XAOSIP
SECPGD_FIELD_MODEL(KFROM)%XAOSIM=>XAOSIM
SECPGD_FIELD_MODEL(KFROM)%XAOSJP=>XAOSJP
SECPGD_FIELD_MODEL(KFROM)%XAOSJM=>XAOSJM
SECPGD_FIELD_MODEL(KFROM)%XHO2IP=>XHO2IP
SECPGD_FIELD_MODEL(KFROM)%XHO2IM=>XHO2IM
SECPGD_FIELD_MODEL(KFROM)%XHO2JP=>XHO2JP
SECPGD_FIELD_MODEL(KFROM)%XHO2JM=>XHO2JM
SECPGD_FIELD_MODEL(KFROM)%XAVG_ZS=>XAVG_ZS
SECPGD_FIELD_MODEL(KFROM)%XSIL_ZS=>XSIL_ZS
!
! Current model is set to model KTO
XAOSIP=>SECPGD_FIELD_MODEL(KTO)%XAOSIP
XAOSIM=>SECPGD_FIELD_MODEL(KTO)%XAOSIM
XAOSJP=>SECPGD_FIELD_MODEL(KTO)%XAOSJP
XAOSJM=>SECPGD_FIELD_MODEL(KTO)%XAOSJM
XHO2IP=>SECPGD_FIELD_MODEL(KTO)%XHO2IP
XHO2IM=>SECPGD_FIELD_MODEL(KTO)%XHO2IM
XHO2JP=>SECPGD_FIELD_MODEL(KTO)%XHO2JP
XHO2JM=>SECPGD_FIELD_MODEL(KTO)%XHO2JM
XAVG_ZS=>SECPGD_FIELD_MODEL(KTO)%XAVG_ZS
XSIL_ZS=>SECPGD_FIELD_MODEL(KTO)%XSIL_ZS

END SUBROUTINE SECPGD_FIELD_GOTO_MODEL

END MODULE MODD_SECPGD_FIELD_n
