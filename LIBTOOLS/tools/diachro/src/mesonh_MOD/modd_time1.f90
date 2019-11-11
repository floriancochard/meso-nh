!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_timen.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_TIME1
!     ##################
!
!!****  *MODD_TIME1* - declaration of temporal grid variables
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
!!      Book2 of documentation of Meso-NH (module MODD_TIME1)
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
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE
!
IMPLICIT NONE
!
TYPE (DATE_TIME),SAVE :: TDTMOD        ! Time and Date of the  model beginning 
TYPE (DATE_TIME),SAVE :: TDTCUR        ! Current Time and Date  
TYPE (DATE_TIME),SAVE :: TDTRAD_FULL   ! Time and Date of the last full
                                       ! radiation call
TYPE (DATE_TIME),SAVE :: TDTRAD_CLONLY ! Time and Date of the last radiation 
                                       ! call for only the cloudy verticals
TYPE (DATE_TIME),SAVE :: TDTDCONV ! Time and Date of the last deep convection
                                  ! call
!
END MODULE MODD_TIME1
