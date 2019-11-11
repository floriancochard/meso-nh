!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 14:17:24
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_MEAN_FIELD
!     ###################
!
!!****  *MODD_MEAN_FIELD* - declaration of prognostic variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     prognostic variables. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_MEAN_FIELDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!    P. Aumond   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       12/2009                     
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS

IMPLICIT NONE

LOGICAL :: LMEAN_FIELD          

END MODULE MODD_MEAN_FIELD
