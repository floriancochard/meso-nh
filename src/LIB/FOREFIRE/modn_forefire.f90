!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODN_FOREFIRE
!     ##################
!-------------------------------------------------------------------------------
!***	MODN_FOREFIRE  Declaration of namelist NAM_FOREFIRE
!
!!    AUTHOR
!!    ------
!	           : P. Tulet (LACy / CNRM)
!!              X. Pialat (SPE)
!	Creation   : 09.10.2010
!-------------------------------------------------------------------------------
!
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_FOREFIRE
!
IMPLICIT NONE
!
NAMELIST /NAM_FOREFIRE/LFOREFIRE, COUPLINGRES, &
		NFFSCALARS, FFSVNAMES, &
		LFFCHEM, NFFCHEMVAR, FFCVNAMES, &
		NFFCHEMVAROUT, FFCONAMES, &
		FFOUTUPS, FLOWOUT, PHYSOUT, CHEMOUT
! 
END MODULE MODN_FOREFIRE
