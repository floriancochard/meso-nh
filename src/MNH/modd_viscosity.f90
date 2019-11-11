!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #######################
          MODULE MODD_VISCOSITY
!      #######################
!
!!****   *MODD_VISCOSITY*  - declaration of viscosity forces constants
!!
!!     PURPOSE
!!     -------
!        The purpose of this declarative module is to declare the 
!     viscosity forces constants
!!
!!**   IMPLICIT ARGUMENTS
!!     ------------------
!!       NONE
!!
!!    
!!     AUTHOR
!!     ------
!1       Jeanne Colin         * Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original            13/04/11
!----------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE
!
LOGICAL :: LVISC          ! Logical switch to activate viscosity 

LOGICAL :: LVISC_UVW          ! Logical switch to activate viscosity for the
!                               momentum
LOGICAL :: LVISC_TH          ! Logical switch to activate viscosity for the
!                               potential temperature
LOGICAL :: LVISC_SV          ! Logical switch to activate viscosity for the
!                               scalar tracer
LOGICAL :: LVISC_R          ! Logical switch to activate viscosity for the
!                                moisture
REAL, SAVE :: XMU_V        ! Molecular (cinematic) viscosity
REAL, SAVE :: XPRANDTL     ! Prandtl number

END MODULE MODD_VISCOSITY
