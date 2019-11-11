!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODD_IO_NAM
!     ######################
!
!!****  *MODD_IO_NAM* Keep in memory the namelist file names
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
IMPLICIT NONE
!------------------------------------------------------------------------------
!
TYPE(TFILEDATA),POINTER :: TNAM  => NULL() ! namelist file
TYPE(TFILEDATA),POINTER :: TFILE => NULL() ! file
!
!------------------------------------------------------------------------------
!
END MODULE MODD_IO_NAM

