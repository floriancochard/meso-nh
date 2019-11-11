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
!      #####################
       MODULE MODD_STAND_ATM
!      #####################
!
!!****  *MODD_STAND_ATM * - declaration of 5 standard atmospheres
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to allocate 5 arrays for
!!    the tropical, summer/winter mid-latitudes and summer/winter polar
!!    standard atmospheres
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_RADIATIONS$n)
!!
!!    AUTHOR
!!    ------
!!     J.-P. Pinty   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      26/02/95
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XSTROATM ! Standard tropical atmosphere
REAL, DIMENSION(:,:), ALLOCATABLE :: XSMLSATM ! Standard summer atmosphere
					       !    in the mid-latitudes
REAL, DIMENSION(:,:), ALLOCATABLE :: XSMLWATM ! Standard winter atmosphere
					       !    in the mid-latitudes
REAL, DIMENSION(:,:), ALLOCATABLE :: XSPOSATM ! Standard summer atmosphere
         	        		       !   above the polar circles
REAL, DIMENSION(:,:), ALLOCATABLE :: XSPOWATM ! Standard winter atmosphere
					       !   above the polar circles
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_STAND_ATM
