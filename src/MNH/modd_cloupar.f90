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
!     ####################
      MODULE MODD_CLOUDPAR
!     ####################
!
!!****  *MODD_CLOUDPAR* - declaration of Microphysics constants 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the 
!     Microhysics constants.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CLOUDPAR)
!!          
!!    AUTHOR
!!    ------
!!	E. Richard   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/12/95                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
REAL,SAVE :: XCEXVT               ! constant in the rain drop fall velocity
!
REAL,SAVE :: XC1RC, XC2RC         ! constants for autoconversion
!
REAL,SAVE :: XCEXRA, XCRA         ! constants for accretion     
!
REAL,SAVE :: XCEXRE, XC1RE, XC2RE ! constants for rain evaporation 
!
REAL,SAVE :: XCEXRS, XCRS         ! constants for rain sedimentation
!
REAL,SAVE :: XDIVA                ! vapor diffusivity in air
REAL,SAVE :: XTHCO                ! thermal conductivity 
!
END MODULE MODD_CLOUDPAR 
