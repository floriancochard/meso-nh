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
!     ###########################
      MODULE MODD_EMIS_PGDFIELDS
!     ###########################
!
!!****  *MODD_EMIS_PGDFIELDS* - declaration of chemical emission data arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     chemical emission data arrays.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_EMIS_PGDFIELDS)
!!      
!!
!!    AUTHOR
!!    ------
!!	D. Gazen   *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/03/2001                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER,PARAMETER :: JPEMISMAX = 10000
INTEGER                                     :: NEMIS_PGD_NBR
!                          ! number of chemical pgd fields chosen by user
CHARACTER(LEN=3) , DIMENSION(JPEMISMAX) :: CEMIS_PGD_AREA
!                          ! areas where chemical pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
!                          !
CHARACTER(LEN=40), DIMENSION(JPEMISMAX) :: CEMIS_PGD_NAME
!                          ! name of the chemical pgd fields (emitted species)
CHARACTER(LEN=40), DIMENSION(JPEMISMAX) :: CEMIS_PGD_COMMENT ! comment
INTEGER,           DIMENSION(JPEMISMAX) :: NEMIS_PGD_TIME ! emission time
REAL,     DIMENSION(:,:,:), ALLOCATABLE :: XEMIS_PGD_FIELDS
!                          ! emission pgd fields values
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_EMIS_PGDFIELDS
