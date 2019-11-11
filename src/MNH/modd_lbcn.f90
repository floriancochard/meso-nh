!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_LBC_n
!     #################
!
!!****  *MODD_LBC$n* - declaration of lateral boundary conditions
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the variables
!     concerning the lateral boundary conditions.  
!        
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_LBCn)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq and J-P. Lafore    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/09/94     
!!                  15/03/95  (J.Stein) remove R from the historical variables
!!                  15/06/95  (J.Stein) add EPS related variables                 
!!                  29/04/02  (P.Jabouille) remove useless variables
!!                  26/06/13  (C.Lac) Introduction of CPHASE_PBL
!!                  01/04/14  (C.Lac) Introduction of CARPKMAX different to
!!                                    RIMKMAX
!!                  08/2015  (Redelsperger & Pianezze) add XPOND coefficient for LBC
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE LBC_t
!
  CHARACTER(LEN=4), DIMENSION(2) :: CLBCX = '' ! X-direction LBC type at left(1)
                                             ! and right(2) boundaries
  CHARACTER(LEN=4), DIMENSION(2) :: CLBCY = '' ! Y-direction LBC type at left(1)
                                             ! and right(2) boundaries
  INTEGER, DIMENSION(2)          :: NLBLX ! X-direction characteristic large
                                             ! scale length at left(1) and 
                                             ! right(2) boundaries ( number of
                                             ! delta x)
  INTEGER, DIMENSION(2)          :: NLBLY ! Y-direction characteristic large
                                             ! scale length at left(1) and 
                                             ! right(2) boundaries ( number of
                                             ! delta y)
  REAL                         :: XCPHASE ! prescribed value of the phase
                                             ! velocity if constant
  REAL                         :: XCPHASE_PBL ! prescribed value of the phase
                                             ! velocity in the PBL if constant
  REAL                         :: XCARPKMAX   ! relaxation coefficient for
                                             !Carpenter, different to XRIMKMAX
  REAL                         :: XPOND      ! relaxation coefficient for LBC
END TYPE LBC_t

TYPE(LBC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: LBC_MODEL

CHARACTER(LEN=4), DIMENSION(:), POINTER :: CLBCX=>NULL()
CHARACTER(LEN=4), DIMENSION(:), POINTER :: CLBCY=>NULL()
INTEGER, DIMENSION(:), POINTER :: NLBLX=>NULL()
INTEGER, DIMENSION(:), POINTER :: NLBLY=>NULL()
REAL, POINTER :: XCPHASE=>NULL()
REAL, POINTER :: XCPHASE_PBL=>NULL()
REAL, POINTER :: XCARPKMAX=>NULL()
REAL, POINTER :: XPOND=>NULL()

CONTAINS

SUBROUTINE LBC_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
CLBCX=>LBC_MODEL(KTO)%CLBCX
CLBCY=>LBC_MODEL(KTO)%CLBCY
NLBLX=>LBC_MODEL(KTO)%NLBLX
NLBLY=>LBC_MODEL(KTO)%NLBLY
XCPHASE=>LBC_MODEL(KTO)%XCPHASE
XCPHASE_PBL=>LBC_MODEL(KTO)%XCPHASE_PBL
XCARPKMAX=>LBC_MODEL(KTO)%XCARPKMAX  
XPOND=>LBC_MODEL(KTO)%XPOND

END SUBROUTINE LBC_GOTO_MODEL

END MODULE MODD_LBC_n
