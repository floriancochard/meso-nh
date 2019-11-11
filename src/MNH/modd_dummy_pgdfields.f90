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
      MODULE MODD_DUMMY_PGDFIELDS
!     ####################
!
!!****  *MODD_DUMMY_PGDFIELDS* - declaration of dummy physiographic data arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     dummy physiographic data arrays.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DUMMY_PGDFIELDS)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/01/95                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER                                          :: NDUMMY_PGD_NBR
!                          ! number of dummy pgd fields chosen by user
CHARACTER(LEN=3) , DIMENSION(1000)               :: CDUMMY_PGD_AREA
!                          ! areas where dummy pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
!                          !
CHARACTER(LEN=20), DIMENSION(1000)               :: CDUMMY_PGD_NAME
!                          ! name of the dummy pgd fields (for information)
REAL,              DIMENSION(:,:,:), ALLOCATABLE :: XDUMMY_PGD_FIELDS
!                          ! dummy pgd fields themselves
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_DUMMY_PGDFIELDS
