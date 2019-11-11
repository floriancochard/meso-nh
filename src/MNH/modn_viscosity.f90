!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###################
      MODULE MODN_VISCOSITY
!     ###################
!
!!****  *MODN_VISCOSITY* - declaration of namelist NAM_VISC
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_VISC
!     which concern the parameters of the viscosity forces for all models
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	    J. Colin                  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    April 2011
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_VISCOSITY
!
IMPLICIT NONE
!
NAMELIST/NAM_VISC/LVISC,LVISC_UVW,LVISC_TH,LVISC_SV,LVISC_R,XMU_V,XPRANDTL
!
END MODULE MODN_VISCOSITY
