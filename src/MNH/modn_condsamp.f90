!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODN_CONDSAMP
!     ##################
!-------------------------------------------------------------------------------
!***	MODD_CONDSAMP  Declaration of namelist NAM_CONDSAMP
!
!!    AUTHOR
!!    ------
!	           : C.Lac                            
!	Creation   : 05.06.2011
!-------------------------------------------------------------------------------
!
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CONDSAMP
!
IMPLICIT NONE
!
NAMELIST /NAM_CONDSAMP/ &
     LCONDSAMP,NCONDSAMP,XRADIO,XSCAL,XHEIGHT_BASE,XDEPTH_BASE, &
     XHEIGHT_TOP,XDEPTH_TOP
!
END MODULE MODN_CONDSAMP
