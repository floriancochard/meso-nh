!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      MODULE MODN_DIAG_BLANK
!     ######################
!
!!****  *MODN_DIAG_BLANK* -  Namelist module for MesoNH developpers namelist in programm diag
!!
!!    PURPOSE
!!    -------
!!
!!       The purpose of this module is to specify the namelist NAM_DIAG_BLANK
!!      which offer dummy real, integer, logical and character variables for
!!      test and debugging purposes.
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DIAG_BLANK : contains declaration of dummy variables
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!	K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!! 
!!    Original 25/04/96
!!    Modification 14/12/00 (P.Jabouille) add dummy arrays
!!    Modification 29/05/01 (G.Jaubert)  add _DIAG at the end of the DUMMY names
!!                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DIAG_BLANK
!
IMPLICIT NONE
!
NAMELIST /NAM_DIAG_BLANK/ XDUMMY_DIAG,NDUMMY_DIAG,LDUMMY_DIAG,CDUMMY_DIAG
!
END MODULE MODN_DIAG_BLANK
