!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! NEC0 masdev4_7 2007/06/16 01:41:59
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_DUMMY_GR_FIELD_n
!     ####################
!
!!****  *MODD_DUMMY_GR_FIELD$n* - declaration of dummy physiographic data arrays
!!                          for model n
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     dummy physiographic data arrays for model n.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DUMMY_GR_FIELD)
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
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE DUMMY_GR_FIELD_t
!
  INTEGER                                          :: NDUMMY_GR_NBR
!                          ! number of dummy pgd fields chosen by user
!JUAN
  CHARACTER(LEN=3) , DIMENSION(:), POINTER         :: CDUMMY_GR_AREA => NULL()
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
  CHARACTER(LEN=20), DIMENSION(:), POINTER         :: CDUMMY_GR_NAME => NULL()
!                          ! name of the dummy pgd fields (for information)
!JUAN
  REAL,              DIMENSION(:,:,:), POINTER     :: XDUMMY_GR_FIELDS=>NULL()
!                          ! dummy pgd fields themselves
!
!-------------------------------------------------------------------------------
!
END TYPE DUMMY_GR_FIELD_t

TYPE(DUMMY_GR_FIELD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: DUMMY_GR_FIELD_MODEL
!JUAN
LOGICAL    , DIMENSION(JPMODELMAX),         SAVE :: DUMMY_GR_FIELD_FIRST_CALL = .TRUE.
!JUAN

INTEGER, POINTER :: NDUMMY_GR_NBR=>NULL()
CHARACTER(LEN=3) , DIMENSION(:), POINTER :: CDUMMY_GR_AREA=>NULL()
CHARACTER(LEN=20), DIMENSION(:), POINTER :: CDUMMY_GR_NAME=>NULL()
REAL,              DIMENSION(:,:,:), POINTER :: XDUMMY_GR_FIELDS=>NULL()

CONTAINS

SUBROUTINE DUMMY_GR_FIELD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
!JUAN
IF (DUMMY_GR_FIELD_FIRST_CALL(KTO)) THEN
ALLOCATE (DUMMY_GR_FIELD_MODEL(KTO)%CDUMMY_GR_AREA(1000))
ALLOCATE (DUMMY_GR_FIELD_MODEL(KTO)%CDUMMY_GR_NAME(1000))
DUMMY_GR_FIELD_FIRST_CALL(KTO) = .FALSE.
ENDIF
!JUAN
! Save current state for allocated arrays
DUMMY_GR_FIELD_MODEL(KFROM)%XDUMMY_GR_FIELDS=>XDUMMY_GR_FIELDS
!
! Current model is set to model KTO
NDUMMY_GR_NBR=>DUMMY_GR_FIELD_MODEL(KTO)%NDUMMY_GR_NBR
CDUMMY_GR_AREA=>DUMMY_GR_FIELD_MODEL(KTO)%CDUMMY_GR_AREA
CDUMMY_GR_NAME=>DUMMY_GR_FIELD_MODEL(KTO)%CDUMMY_GR_NAME
XDUMMY_GR_FIELDS=>DUMMY_GR_FIELD_MODEL(KTO)%XDUMMY_GR_FIELDS

END SUBROUTINE DUMMY_GR_FIELD_GOTO_MODEL

END MODULE MODD_DUMMY_GR_FIELD_n
