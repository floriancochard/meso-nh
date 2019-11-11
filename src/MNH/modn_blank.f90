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
!     #################
      MODULE MODN_BLANK
!     #################
!
!!****  *MODN_BLANK* -  Namelist module for MesoNH developpers namelist
!!
!!    PURPOSE
!!    -------
!!
!!       The purpose of this module is to specify the namelist NAM_BLANK
!!      which offer dummy real, integer, logical and character variables for
!!      test and debugging purposes.
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_BLANK : contains declaration of dummy variables
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
!!                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_BLANK
!
IMPLICIT NONE
!
NAMELIST /NAM_BLANK/ XDUMMY1, XDUMMY2, XDUMMY3, XDUMMY4, &
                     XDUMMY5, XDUMMY6, XDUMMY7, XDUMMY8, &
                     NDUMMY1, NDUMMY2, NDUMMY3, NDUMMY4, &
                     NDUMMY5, NDUMMY6, NDUMMY7, NDUMMY8, &
                     LDUMMY1, LDUMMY2, LDUMMY3, LDUMMY4, &
                     LDUMMY5, LDUMMY6, LDUMMY7, LDUMMY8, &
                     CDUMMY1, CDUMMY2, CDUMMY3, CDUMMY4, &
                     CDUMMY5, CDUMMY6, CDUMMY7, CDUMMY8, &
                     XDUMMY,NDUMMY,LDUMMY,CDUMMY
!
END MODULE MODN_BLANK
