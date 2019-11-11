!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODN_PASPOL
!     ##################
!-------------------------------------------------------------------------------
!***	MODD_PASPOL  Declaration of namelist NAM_PASPOL
!
!!    AUTHOR
!!    ------
!	           : Michel Bouzom, DP/SERV/ENV
!	Creation   : 09.10.2001
!-------------------------------------------------------------------------------
!
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PASPOL
!
IMPLICIT NONE
!
NAMELIST /NAM_PASPOL/ &
     LPASPOL,NRELEASE,CPPINIT,XPPLAT,XPPLON,XPPMASS , &
     XPPBOT,XPPTOP,CPPT1,CPPT2,CPPT3,CPPT4
!
END MODULE MODN_PASPOL
