!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_TYPE_ALLVAR
!     ###################
!
!!****  *MODD_TYPE_ALLVAR* - Declaration des types de variables 3D, 2D, 1D,
!!
!!    PURPOSE
!!    -------
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!       
!!    AUTHOR
!!    ------
!!      P Jabouille
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       11/08/97                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
TYPE X_Y_Z_
  CHARACTER(LEN=16)     :: NAME
  INTEGER               :: IGRID
  CHARACTER(LEN=16)     :: UNITS
END TYPE X_Y_Z_
!
TYPE X_Y_
  CHARACTER(LEN=16)     :: NAME
  INTEGER               :: IGRID
  CHARACTER(LEN=16)     :: UNITS
END TYPE X_Y_
!
TYPE VX_VY_VZ_
  CHARACTER(LEN=16),DIMENSION(3)     :: NAME
  INTEGER,DIMENSION(3)               :: IGRID
  CHARACTER(LEN=16),DIMENSION(3)     :: UNITS
END TYPE VX_VY_VZ_
!
TYPE VX_VY_
  CHARACTER(LEN=16),DIMENSION(3)     :: NAME
  INTEGER,DIMENSION(3)               :: IGRID
  CHARACTER(LEN=16),DIMENSION(3)     :: UNITS
END TYPE VX_VY_
!
TYPE Z_
  CHARACTER(LEN=16)     :: NAME
  INTEGER               :: IGRID
  CHARACTER(LEN=16)     :: UNITS
END TYPE Z_
!
END MODULE MODD_TYPE_ALLVAR
