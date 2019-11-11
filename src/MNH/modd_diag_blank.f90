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
!     ######################
      MODULE MODD_DIAG_BLANK
!     ######################
!
!!****  *MODD_DIAG_BLANK* -  Declarative module for MesoNH developpers namelist used in diag programm
!!
!!    PURPOSE
!!    -------
!!
!!      Offer dummy real, integer, logical and character variables for
!!    test and debugging purposes.
!!
!!**  METHOD
!!    ------
!!
!!      Eight dummy real, integer, logical and character*80 variables are
!!    defined and passed through the namelist read operations. None of the
!!    MesoNH routines uses any of those variables. When a developper choses
!!    to introduce temporarily a parameter to some subroutine, he has to
!!    introduce a USE MODD_BLANK statement into that subroutine. Then he
!!    can use any of the variables defined here and change them easily via
!!    the namelist input.
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
!!      updated     17/11/00  (P Jabouille) Use dummy array
!!      updated     29/05/01  (G.Jaubert) add _DIAG at the end of the DUMMY names
!!                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
REAL,    SAVE, DIMENSION(JPDUMMY) :: XDUMMY_DIAG
INTEGER, SAVE, DIMENSION(JPDUMMY) :: NDUMMY_DIAG
LOGICAL, SAVE, DIMENSION(JPDUMMY) :: LDUMMY_DIAG
CHARACTER*80, SAVE, DIMENSION(JPDUMMY) :: CDUMMY_DIAG
!
END MODULE MODD_DIAG_BLANK
