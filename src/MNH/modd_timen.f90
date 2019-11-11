!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_TIME_n
!     ##################
!
!!****  *MODD_TIME$n* - declaration of temporal grid variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables
!     which concern the time for one nested model.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TIME : contains the definition of the types for time 
!!                              variables and time variables for all model
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_TIME$n)
!!       
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/07/94       
!!      J.Stein     27/10/95   add the radiation call's instants               
!!      P.Bechtold  26/03/96   add the last deep convection call
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE (DATE_TIME), POINTER :: TDTMOD=>NULL()        ! Time and date of model beginning
TYPE (DATE_TIME), POINTER :: TDTCUR=>NULL()        ! Current time and date
TYPE (DATE_TIME), POINTER :: TDTRAD_FULL=>NULL()   ! Time and date of last full radiation call
TYPE (DATE_TIME), POINTER :: TDTRAD_CLONLY=>NULL() ! Time and date of last radiation call for only cloudy verticals
TYPE (DATE_TIME), POINTER :: TDTDCONV=>NULL()      ! Time and date of the last deep convection call

CONTAINS

SUBROUTINE TIME_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
END SUBROUTINE TIME_GOTO_MODEL

END MODULE MODD_TIME_n
