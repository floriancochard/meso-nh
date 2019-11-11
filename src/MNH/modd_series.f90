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
!     ##################
      MODULE MODD_SERIES
!     ##################
!
!!****  *MODD_SERIES*- declaration of variables related with the diagnostics
!!                     for diachro files
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify
!   the conditions of realization of the diagnostics (box and slice definition)
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
!!      V. Ducrocq              *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      29/01/98
!!                Oct. 10,1998 (Lafore) adaptation of Diagnostics 
!!                                      to the sequential nesting version
!!                Oct. 2011 : (P.Le Moigne) Surface series
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
LOGICAL, SAVE :: LSERIES ! switch for temporal series 
LOGICAL, SAVE :: LMASKLANDSEA ! logical for additional diagnostics that separate  sea and land points
LOGICAL, SAVE :: LWMINMAX ! switch to compute max and min of W
LOGICAL, SAVE :: LSURF    ! switch to compute surface diagnostics 
!
END MODULE MODD_SERIES
