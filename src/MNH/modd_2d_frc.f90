!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/11/23 17:28:26
!-----------------------------------------------------------------
!     ######################## 
      MODULE MODD_2D_FRC 
!     ########################
!
!!****  *MODD_2D_FRC* - declaration of the control parameters for
!!                           2d forcong : advection and relaxation
!!
!!
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
LOGICAL, SAVE :: L2D_ADV_FRC 
LOGICAL, SAVE :: L2D_REL_FRC
REAL, SAVE    :: XRELAX_HEIGHT_BOT
REAL, SAVE    :: XRELAX_HEIGHT_TOP
REAL, SAVE    :: XRELAX_TIME
!
END MODULE MODD_2D_FRC
