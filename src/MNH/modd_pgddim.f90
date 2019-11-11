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
      MODULE MODD_PGDDIM
!     ##################
!
!!****  *MODD_PGDDIM* - declaration of dimensions
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the dimensions 
!     of the physiographic data arrays in physiographic data file.   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PGDDIM)
!!          
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/01/95                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
INTEGER :: NPGDIMAX,NPGDJMAX  !  Dimensions respectively  in x , y  directions.
!
!
!-------------------------------------------------------------------------------
END MODULE MODD_PGDDIM
