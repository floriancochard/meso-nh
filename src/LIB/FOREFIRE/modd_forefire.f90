!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_FOREFIRE
!     ##################
!-------------------------------------------------------------------------------
!***	MODD_LAVA  Declaration of lava module
!
!!    AUTHOR
!!    ------
!	           : P. Tulet,  LACy / CNRM
!!            : X. Pialat, SPE
!	Creation   : 15.02.2012
!
!-------------------------------------------------------------------------------
!
!
!*    0. DECLARATIONS
!        ------------
! 
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
LOGICAL     :: LFOREFIRE = .FALSE.  				! Switch to activate ForeFire coupling
REAL			:: COUPLINGRES = 100.					! Coupling resolution (above the model is not coupled with ForeFire)
INTEGER		:: NFFSCALARS = 0							! Number of passive tracers in ForeFire
CHARACTER(LEN=6), DIMENSION(10) 	:: FFSVNAMES	! Names of the scalar variables
LOGICAL     :: LFFCHEM = .FALSE.  					! Switch to activate chemistry in ForeFire coupling
INTEGER		:: NFFCHEMVAR = 0							! Number of chemical variables in ForeFire coupling
CHARACTER(LEN=6), DIMENSION(10) 	:: FFCVNAMES	! Names of the chemical variables
INTEGER		:: NFFCHEMVAROUT = 0						! Number of chemical variables for outputs only
CHARACTER(LEN=6), DIMENSION(10) 	:: FFCONAMES	! Names of the chemical variables
REAL, DIMENSION(6)					:: FFOUTUPS		! Outputs updates periods for the MNH models
INTEGER, DIMENSION(6)				:: FLOWOUT		! Booleans for flow variables outputs
INTEGER, DIMENSION(6)				:: PHYSOUT		! Booleans for physical variables outputs
INTEGER, DIMENSION(6)				:: CHEMOUT		! Booleans for physical variables outputs

!
END MODULE MODD_FOREFIRE
