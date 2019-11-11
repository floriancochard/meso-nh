!MNH_LIC Copyright 2016-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ########################
      MODULE MODD_CH_PRODLOSSTOT_n
!!    ########################
!
!!****  *MODD_CH_PRODLOSSTOT_n* - declaration of parameters for chemical prod/loss diagnostic
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the
!     variables that may be used to compute the integrated chemical production 
!     and destruction terms as diagnostic (DIAG step)
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!          
!!    AUTHOR
!!    ------
!!	    F. Brosse *Laboratoire d'Aerologie UPS-CNRS*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    October 2016
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_PRODLOSSTOT_t
!
    CHARACTER(LEN=32),  DIMENSION(:), POINTER :: CNAMES_PRODLOSST => NULL()
    INTEGER,            DIMENSION(:), POINTER :: NIND_SPEC => NULL()
    INTEGER                                   :: NEQ_PLT
    REAL,         DIMENSION(:,:,:,:), POINTER:: XPROD => NULL()
    REAL,         DIMENSION(:,:,:,:), POINTER:: XLOSS => NULL()
!
!-----------------------------------------------------------------------------

END TYPE CH_PRODLOSSTOT_t

TYPE(CH_PRODLOSSTOT_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_PRODLOSSTOT_MODEL

CHARACTER(LEN=32),  DIMENSION(:), POINTER:: CNAMES_PRODLOSST=>NULL()
INTEGER,            DIMENSION(:), POINTER:: NIND_SPEC=>NULL()
INTEGER, POINTER :: NEQ_PLT=>NULL()
REAL,         DIMENSION(:,:,:,:), POINTER:: XPROD=>NULL()
REAL,         DIMENSION(:,:,:,:), POINTER:: XLOSS=>NULL()

CONTAINS

SUBROUTINE CH_PRODLOSSTOT_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CH_PRODLOSSTOT_MODEL(KFROM)%CNAMES_PRODLOSST=>CNAMES_PRODLOSST
CH_PRODLOSSTOT_MODEL(KFROM)%NIND_SPEC=>NIND_SPEC
CH_PRODLOSSTOT_MODEL(KFROM)%XPROD=>XPROD
CH_PRODLOSSTOT_MODEL(KFROM)%XLOSS=>XLOSS
!
! Current model is set to model KTO
CNAMES_PRODLOSST=>CH_PRODLOSSTOT_MODEL(KTO)%CNAMES_PRODLOSST
NIND_SPEC=>CH_PRODLOSSTOT_MODEL(KTO)%NIND_SPEC
NEQ_PLT=>CH_PRODLOSSTOT_MODEL(KTO)%NEQ_PLT
XPROD=>CH_PRODLOSSTOT_MODEL(KTO)%XPROD
XLOSS=>CH_PRODLOSSTOT_MODEL(KTO)%XLOSS

END SUBROUTINE CH_PRODLOSSTOT_GOTO_MODEL

END MODULE MODD_CH_PRODLOSSTOT_n
