!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################  
      MODULE MODN_PARAM_C1R3
!     ######################
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAM_C1R3
!
IMPLICIT NONE
!
NAMELIST/NAM_PARAM_C1R3/XALPHAI,XNUI,XALPHAS,XNUS,XALPHAG,XNUG, &
                        XFACTNUC_DEP,XFACTNUC_CON,              &
                        LSEDI,LHHONI,                           &
                        CPRISTINE_ICE_C1R3,CHEVRIMED_ICE_C1R3
!
END MODULE MODN_PARAM_C1R3
