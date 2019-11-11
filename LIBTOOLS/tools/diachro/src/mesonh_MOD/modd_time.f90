!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_time.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE MODD_TIME
!     #################
!
!!****  *MODD_TIME* - declaration of temporal grid variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables
!     which concern the time for all models
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!      Book2 of documentation of Meso-NH (module MODD_TIME)
!!       
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/07/94                      
!!      Modification 10/03/95 (I.Mallet)   add the coupling times
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
TYPE (DATE_TIME), SAVE :: TDTEXP      ! Time and Date of Experiment beginning 
TYPE (DATE_TIME), SAVE :: TDTSEG      ! Time and Date of the  segment beginning 
!
TYPE (DATE_TIME), SAVE, DIMENSION(JPCPLFILEMAX) :: TDTCPL ! Time and Date of 
                                                          ! the CouPLing files
END MODULE MODD_TIME
