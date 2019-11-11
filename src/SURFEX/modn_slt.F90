!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!!
!!    #####################
      MODULE MODN_SLT
!!    #####################
!!
!!*** *MODN_DUST*
!!
!!    PURPOSE
!!    -------
!       Namelist for SEA SALT EMISSION SCHEME aerosol scheme parameters 
!!
!!**  AUTHOR
!!    ------
!!    A. Grini / P. Tulet     *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/02/05
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_SLT_SURF
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
NAMELIST /NAM_SURF_SLT/  &
       CEMISPARAM_SLT, LVARSIG_SLT, LRGFIX_SLT, JPMODE_SLT            !Parameterization type
!
END MODULE MODN_SLT
