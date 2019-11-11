!MNH_LIC Copyright 2018-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!     ######################
       MODULE MODD_BLOWSNOW_n
!!     ######################
!!
!!     PURPOSE
!!     -------
!!
!!     declaration of variables and types for the BLOWSNOW scheme
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     none
!!
!!
!!     AUTHOR
!!     ------
!!     Vincent Vionnet (CNRM)
!!
!!
!!     MODIFICATIONS
!!     -------------
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE BLOWSNOW_t
!
LOGICAL      :: LSNOWSUBL    ! switch to activate blowing snow sublimation
!
REAL, DIMENSION(:,:,:), POINTER :: XSNWSUBL3D => NULL() ! Drifting snow instataneous
!                              sublimation rate (kg/m3/s)
REAL, DIMENSION(:,:,:), POINTER :: XSNWCANO => NULL() ! Total mass in Canopy at time t
!     (:,:,1) : equivalent number concentration in Canopy (#/kg)
!     (:,:,2) : equivalent mass concentration in Canopy (kg/kg)
!     (:,:,3) : equivalent mass concentration in saltation   (kg/kg)  
REAL, DIMENSION(:,:,:), POINTER :: XRSNWCANOS => NULL() ! Source of (rho*canopy mass) at time t



END TYPE BLOWSNOW_t

TYPE(BLOWSNOW_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: BLOWSNOW_MODEL

REAL, DIMENSION(:,:,:), POINTER :: XSNWSUBL3D=> NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSNWCANO=> NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRSNWCANOS=> NULL()


LOGICAL, POINTER :: LSNOWSUBL=>NULL()

CONTAINS

SUBROUTINE BLOWSNOW_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
BLOWSNOW_MODEL(KFROM)%XSNWSUBL3D=>XSNWSUBL3D
BLOWSNOW_MODEL(KFROM)%XSNWCANO=>XSNWCANO
BLOWSNOW_MODEL(KFROM)%XRSNWCANOS=>XRSNWCANOS

!
! Current model is set to model KTO
XSNWSUBL3D=>BLOWSNOW_MODEL(KTO)%XSNWSUBL3D
XSNWCANO=>BLOWSNOW_MODEL(KTO)%XSNWCANO
XRSNWCANOS=>BLOWSNOW_MODEL(KTO)%XRSNWCANOS

LSNOWSUBL=>BLOWSNOW_MODEL(KTO)%LSNOWSUBL

END SUBROUTINE BLOWSNOW_GOTO_MODEL
END MODULE MODD_BLOWSNOW_n
