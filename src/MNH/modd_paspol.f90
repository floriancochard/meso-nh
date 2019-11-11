!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_PASPOL
!     ##################
!-------------------------------------------------------------------------------
!***	MODD_PASPOL  Declaration of passive pollutants           
!
!!    AUTHOR
!!    ------
!	           : Michel Bouzom, DP/SERV/ENV
!	Creation   : 09.10.2001
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
LOGICAL      :: LPASPOL = .FALSE.  ! Switch to active passive pollutants
!
INTEGER, PARAMETER :: JPRELEASEMAX = 100
!
INTEGER                               :: NRELEASE     ! Number of releases
CHARACTER*3,  DIMENSION(JPRELEASEMAX) :: CPPINIT        ! Type of initialiZation.
REAL,         DIMENSION(JPRELEASEMAX) :: XPPLAT         ! Latitude  of the release
REAL,         DIMENSION(JPRELEASEMAX) :: XPPLON         ! Longitude of the release
REAL,         DIMENSION(JPRELEASEMAX) :: XPPMASS        ! Released mass     
REAL,         DIMENSION(JPRELEASEMAX) :: XPPBOT         ! Bottom of release
REAL,         DIMENSION(JPRELEASEMAX) :: XPPTOP         ! Top of release   
CHARACTER*14, DIMENSION(JPRELEASEMAX) :: CPPT1          ! Begin of release
CHARACTER*14, DIMENSION(JPRELEASEMAX) :: CPPT2          ! Begin of constant
                                                             ! release
CHARACTER*14, DIMENSION(JPRELEASEMAX) :: CPPT3          ! End of constant 
                                                             ! release
CHARACTER*14, DIMENSION(JPRELEASEMAX) :: CPPT4          ! End of release 
!
!
END MODULE MODD_PASPOL
