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
