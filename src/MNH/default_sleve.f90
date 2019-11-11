!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_DEFAULT_SLEVE
!     ######################
!
INTERFACE 
!
      SUBROUTINE DEFAULT_SLEVE(OSLEVE,PLEN1,PLEN2)
!
LOGICAL,                INTENT(OUT)  :: OSLEVE            ! flag for SLEVE coordinate
REAL,                   INTENT(OUT)  :: PLEN1             ! Decay scale for smooth topography
REAL,                   INTENT(OUT)  :: PLEN2             ! Decay scale for small-scale topography deviation
!
END SUBROUTINE DEFAULT_SLEVE
!
END INTERFACE
!
END MODULE MODI_DEFAULT_SLEVE
!
!
!
!     #############################
      SUBROUTINE DEFAULT_SLEVE(OSLEVE,PLEN1,PLEN2)
!     #############################
!
!!****  *DEFAULT_SLEVE* defaults for SLEVE vertical coordinate
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	G. Zangler      * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        nov 2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
LOGICAL,                INTENT(OUT)  :: OSLEVE            ! flag for SLEVE coordinate
REAL,                   INTENT(OUT)  :: PLEN1             ! Decay scale for smooth topography
REAL,                   INTENT(OUT)  :: PLEN2             ! Decay scale for small-scale topography deviation
!
!
!*       0.2   declarations of local variables
!
!
!-------------------------------------------------------------------------------
OSLEVE = .FALSE.
PLEN1 = 7500.
PLEN2 = 2500.
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_SLEVE
