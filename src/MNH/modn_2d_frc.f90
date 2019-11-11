!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/11/23 17:22:54
!-----------------------------------------------------------------
!     ########################  
      MODULE MODN_2D_FRC
!     ########################
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_2D_FRC
!
IMPLICIT NONE
!
!
NAMELIST/NAM_2D_FRC/L2D_ADV_FRC,L2D_REL_FRC,XRELAX_HEIGHT_BOT, XRELAX_HEIGHT_TOP,XRELAX_TIME
!

END MODULE MODN_2D_FRC
