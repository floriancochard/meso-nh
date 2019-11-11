!MNH_LIC Copyright 2001-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_CH_CONST_n
!     ######################
!
!!
!!    PURPOSE
!!    -------
!     
!   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!
!!    AUTHOR
!!    ------
!!  P. Tulet  (16/01/01) *Meteo France*
!!  M. Leriche (9/12/09) passage en $n
!!
!!    MODIFICATIONS
!!    -------------
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_CONST_t
!

  REAL, DIMENSION(:), POINTER :: XSREALMASSMOLVAL => NULL() ! final molecular
                                                          ! diffusivity value
  REAL, DIMENSION(:), POINTER :: XSREALREACTVAL => NULL() ! final chemical
                                                        ! reactivity factor
                                                        ! with biologie
  REAL, DIMENSION(:,:), POINTER :: XSREALHENRYVAL => NULL() ! chemical Henry
                                                          ! constant value
  REAL                            :: XCONVERSION ! emission unit
                                                     ! conversion factor
!


END TYPE CH_CONST_t

TYPE(CH_CONST_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_CONST_MODEL

REAL, DIMENSION(:), POINTER :: XSREALMASSMOLVAL=>NULL()
REAL, DIMENSION(:), POINTER :: XSREALREACTVAL=>NULL()
REAL, DIMENSION(:,:), POINTER :: XSREALHENRYVAL=>NULL()
REAL, POINTER :: XCONVERSION=>NULL()

CONTAINS

SUBROUTINE CH_CONST_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CH_CONST_MODEL(KFROM)%XSREALMASSMOLVAL=>XSREALMASSMOLVAL
CH_CONST_MODEL(KFROM)%XSREALREACTVAL=>XSREALREACTVAL
CH_CONST_MODEL(KFROM)%XSREALHENRYVAL=>XSREALHENRYVAL
!
! Current model is set to model KTO
XSREALMASSMOLVAL=>CH_CONST_MODEL(KTO)%XSREALMASSMOLVAL
XSREALREACTVAL=>CH_CONST_MODEL(KTO)%XSREALREACTVAL
XSREALHENRYVAL=>CH_CONST_MODEL(KTO)%XSREALHENRYVAL
XCONVERSION=>CH_CONST_MODEL(KTO)%XCONVERSION

END SUBROUTINE CH_CONST_GOTO_MODEL

END MODULE MODD_CH_CONST_n
