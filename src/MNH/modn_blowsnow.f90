!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    #####################
      MODULE MODN_BLOWSNOW
!!    #####################
!!
!!*** *MODN_BLOWSNOW*
!!
!!    PURPOSE
!!    -------
!       Namelist for drifting snow scheme parameters 
!!
!!**  AUTHOR
!!    ------
!!    V.Vionnet      *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 07/04/08
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_BLOWSNOW
IMPLICIT NONE
!
NAMELIST /NAM_BLOWSNOW/ LBLOWSNOW, CSNOWSEDIM,XALPHA_SNOW, XRSNOW
!
END MODULE MODN_BLOWSNOW
