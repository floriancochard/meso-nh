!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_LUNIT
!     #################
!
!!****  *MODD_LUNIT* - declaration of names and logical unit numbers of files 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the 
!     logical unit numbers  of  output file for all models.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_LUNIT)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94                      
!!      V. Masson   01/2004 add file names for use in externalized surface
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
IMPLICIT NONE 
!
CHARACTER(LEN=16),SAVE :: CLUOUT0    ! Name of output_listing file
TYPE(TFILEDATA),POINTER :: TLUOUT0 => NULL() ! output_listing file
TYPE(TFILEDATA),POINTER :: TOUTDATAFILE => NULL() ! output data file being written
TYPE(TFILEDATA),POINTER :: TPGDFILE     => NULL() ! PGD file
!
END MODULE MODD_LUNIT
