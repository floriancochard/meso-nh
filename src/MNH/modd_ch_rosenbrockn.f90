!MNH_LIC Copyright 2007-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!!    #########################
      MODULE MODD_CH_ROSENBROCK_n
!!    #########################
!!
!!*** *MODD_CH_ROSENBROCK*
!!
!!    PURPOSE
!!    -------
!!      This module contains the settings of some of the parameters used in the 
!!    rosenbrock integrator of the kpp distribution
!!
!!**  AUTHOR
!!    ------
!!    J.-P. Pinty    *Laboratoire d'Aerologie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 05/06/07
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_ROSENBROCK_t
!
  INTEGER :: NSPARSEDIM ! size of vectors NSPARSE_IROW and NSPARSE_ICOL
!
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_IROW => NULL() ! row index
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_ICOL => NULL() ! col index
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_CROW => NULL() ! first row element index
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_DIAG => NULL() ! diag index
                                                ! of the sparse JACobian matrix
!
  INTEGER :: NEQ_NAQ    ! number of Non-AQueous species
  INTEGER :: NSPARSEDIM_NAQ !size of vectors NSPARSE_IROW_NAQ and NSPARSE_ICOL_NAQ
!
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_IROW_NAQ => NULL() ! row index
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_ICOL_NAQ => NULL() ! col index
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_CROW_NAQ => NULL() ! first row element index
  INTEGER, POINTER, DIMENSION(:) :: NSPARSE_DIAG_NAQ => NULL() ! diag index
                           ! of the sparse JACobian matrix of NonAQueous species
!
!
!
  INTEGER, DIMENSION(4) :: ICNTRL_ROS = (/ 1,1,4,100 /) !Default
                       ! ICNTRL_ROS(1) =  1   ! autonomous case
                       ! ICNTRL_ROS(2) =  1   ! scalar Tolerance
                       ! ICNTRL_ROS(3) =  4   ! Rodas3 scheme is selected
                       ! ICNTRL_ROS(4) =  100 ! max number of integration steps
  INTEGER, DIMENSION(7) :: RCNTRL_ROS = (/ 0.0,0.0,0.0,0.0,0.0,0.0,0.0 /) !Default
                       ! RCNTRL_ROS(1) =  0.0   ! minimum of the integration 
                                                ! step size
                       ! RCNTRL_ROS(2) =  0.0   ! maximum of the integration 
                                                ! step size
                       ! RCNTRL_ROS(3) =  0.0   ! starting step size
                       ! RCNTRL_ROS(4) =  0.0   ! minimum of step size change
                       ! RCNTRL_ROS(5) =  0.0   ! maximum of step size change
                       ! RCNTRL_ROS(6) =  0.0   ! factor to decrease step after
                                                ! 2 succesive rejections
                       ! RCNTRL_ROS(7) =  0.0   ! safety factor in the 
                                                ! computation of new step size
!
!-----------------------------------------------------------------------------
END TYPE CH_ROSENBROCK_t

TYPE(CH_ROSENBROCK_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_ROSENBROCK_MODEL

INTEGER, POINTER :: NSPARSEDIM=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_IROW=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_ICOL=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_CROW=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_DIAG=>NULL()
INTEGER, POINTER :: NEQ_NAQ=>NULL()
INTEGER, POINTER :: NSPARSEDIM_NAQ=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_IROW_NAQ=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_ICOL_NAQ=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_CROW_NAQ=>NULL()
INTEGER, POINTER, DIMENSION(:) :: NSPARSE_DIAG_NAQ=>NULL()
INTEGER, DIMENSION(:), POINTER :: ICNTRL_ROS=>NULL()
INTEGER, DIMENSION(:), POINTER :: RCNTRL_ROS=>NULL()

CONTAINS

SUBROUTINE CH_ROSENBROCK_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_IROW=>NSPARSE_IROW
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_ICOL=>NSPARSE_ICOL
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_CROW=>NSPARSE_CROW
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_DIAG=>NSPARSE_DIAG
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_IROW_NAQ=>NSPARSE_IROW_NAQ
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_ICOL_NAQ=>NSPARSE_ICOL_NAQ
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_CROW_NAQ=>NSPARSE_CROW_NAQ
CH_ROSENBROCK_MODEL(KFROM)%NSPARSE_DIAG_NAQ=>NSPARSE_DIAG_NAQ
!
! Current model is set to model KTO
NSPARSEDIM=>CH_ROSENBROCK_MODEL(KTO)%NSPARSEDIM
NSPARSE_IROW=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_IROW
NSPARSE_ICOL=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_ICOL
NSPARSE_CROW=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_CROW
NSPARSE_DIAG=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_DIAG
NEQ_NAQ=>CH_ROSENBROCK_MODEL(KTO)%NEQ_NAQ
NSPARSEDIM_NAQ=>CH_ROSENBROCK_MODEL(KTO)%NSPARSEDIM_NAQ
NSPARSE_IROW_NAQ=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_IROW_NAQ
NSPARSE_ICOL_NAQ=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_ICOL_NAQ
NSPARSE_CROW_NAQ=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_CROW_NAQ
NSPARSE_DIAG_NAQ=>CH_ROSENBROCK_MODEL(KTO)%NSPARSE_DIAG_NAQ
ICNTRL_ROS=>CH_ROSENBROCK_MODEL(KTO)%ICNTRL_ROS
RCNTRL_ROS=>CH_ROSENBROCK_MODEL(KTO)%RCNTRL_ROS

END SUBROUTINE CH_ROSENBROCK_GOTO_MODEL

END MODULE MODD_CH_ROSENBROCK_n
