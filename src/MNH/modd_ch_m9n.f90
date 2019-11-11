!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!!    #################
      MODULE MODD_CH_M9_n
!!    #################
!!
!! This code is the MESONH interface to constant and variables defined in module
!! MODD_CH_M9_SCHEME that is proper to one chemical scheme. This interface should
!! reduce the source dependances and then improve the compilation time.
!!
!!*** *MODD_CH_M9*
!!
!!    PURPOSE
!!    -------
!     definition of variables and constant for the chemical core system
!!
!!**  METHOD
!!    ------
!!    The constants NEQ and NREAC are duplicated here in order to avoid
!!    decouple the CCS from the other modules of MNHC.
!!
!!    BEWARE : you must call the procedure 'CH_INIT_SCHEME' before using any
!!             variables from this module.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Didier Gazen (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 19/10/03
!!    12/04/07 (M. Leriche) add NEQAQ for aqueous chemistry
!!    15/07/10 (M. Leriche) add CICNAME array for ice phase species
!!
!!----------------------------------------------------------------------
!!    DECLARATIONS
!!    ------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE METEOTRANSTYPE ! variables from the meteorological part
  REAL,              DIMENSION(20) :: XMETEOVAR  ! the meteorological variables
  CHARACTER(LEN=32), DIMENSION(20) :: CMETEOVAR  ! their names
END TYPE METEOTRANSTYPE

TYPE CH_M9_t
!
  INTEGER :: NEQ           ! number of prognostic chemical species
  INTEGER :: NEQAQ         ! number of prognostic chemical species in aqueous phase
  INTEGER :: NREAC         ! number of chemical reactions
  INTEGER :: NMETEOVARS    ! number of meteorological variables
  INTEGER :: NNONZEROTERMS ! number of non-zero terms returned by CH_TERMS
!
  CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CNAMES=>NULL() ! names of the species
  CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CREACS=>NULL() ! the reaction rate names
  CHARACTER(LEN=256), DIMENSION(:), POINTER  :: CFULLREACS=>NULL() ! the full reactions
!
  CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CICNAMES=>NULL() ! names of ice species
END TYPE CH_M9_t

TYPE(CH_M9_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_M9_MODEL

INTEGER, POINTER :: NEQ=>NULL()
INTEGER, POINTER :: NEQAQ=>NULL()
INTEGER, POINTER :: NREAC=>NULL()
INTEGER, POINTER :: NMETEOVARS=>NULL()
INTEGER, POINTER :: NNONZEROTERMS=>NULL()
CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CNAMES=>NULL()
CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CREACS=>NULL()
CHARACTER(LEN=256), DIMENSION(:), POINTER  :: CFULLREACS=>NULL()
!
CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CICNAMES=>NULL()

CONTAINS

SUBROUTINE CH_M9_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CH_M9_MODEL(KFROM)%CNAMES=>CNAMES
CH_M9_MODEL(KFROM)%CREACS=>CREACS
CH_M9_MODEL(KFROM)%CFULLREACS=>CFULLREACS
CH_M9_MODEL(KFROM)%CICNAMES=>CICNAMES
!
! Current model is set to model KTO
NEQ=>CH_M9_MODEL(KTO)%NEQ
NEQAQ=>CH_M9_MODEL(KTO)%NEQAQ
NREAC=>CH_M9_MODEL(KTO)%NREAC
NMETEOVARS=>CH_M9_MODEL(KTO)%NMETEOVARS
NNONZEROTERMS=>CH_M9_MODEL(KTO)%NNONZEROTERMS
CNAMES=>CH_M9_MODEL(KTO)%CNAMES
CREACS=>CH_M9_MODEL(KTO)%CREACS
CFULLREACS=>CH_M9_MODEL(KTO)%CFULLREACS
CICNAMES=>CH_M9_MODEL(KTO)%CICNAMES

END SUBROUTINE CH_M9_GOTO_MODEL

END MODULE MODD_CH_M9_n
