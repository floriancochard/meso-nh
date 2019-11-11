!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_DIM_ll
!     ##################
!
!!****  *MODD_DIM_ll* - declaration of dimensions 
!!                      and grid-nesting informations
!
!!    Purpose
!!    -------
!       The purpose of this declarative module is to specify  the dimensions
!       of the array, the parameters and the lateral boundaries conditions 
!
!!    Reference
!!    ---------
!          
!!    Author
!!    ------
!
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!     N. Gicquel                * CNRM *
!
!!    Implicit Arguments
!!    ------------------
!
!      None
!
!!    Modifications
!!    -------------
!     Original 04/05/98
!     Philippe Wautelet: 12/01/2018: renamed dimension variables NIMAX_TMP_ll,NJMAX_TMP_ll, NKMAX_TMP_ll
!                                    to prevent mix-up with modd_dimn
!
!-------------------------------------------------------------------------------
!
!     DIMENSIONS
! 
! Dimensions respectively  in x and y directions of the physical domain
!
  INTEGER :: NIMAX_TMP_ll,NJMAX_TMP_ll, NKMAX_TMP_ll
!
! Lower bound and upper bound of the arrays in x direction
!
  INTEGER :: NIINF, NISUP
!
! Lower bound and upper bound of the arrays in y direction
!
  INTEGER :: NJINF, NJSUP
!
!-------------------------------------------------------------------------------
!
!     GRID-NESTING
!
! Model number of the father of each model m
!
  INTEGER, DIMENSION(:), ALLOCATABLE :: NDAD
!
! Ratio for all models
!
  INTEGER, DIMENSION(:), ALLOCATABLE :: NDXRATIO_ALL ! in x-direction 
  INTEGER, DIMENSION(:), ALLOCATABLE :: NDYRATIO_ALL ! in y-direction
!
! Location of models
!   ( horizontal position (i,j) of the ORigin and END of a model m
!     relative to its father NDAD(m) )
!
  INTEGER, DIMENSION(:), ALLOCATABLE :: NXOR_ALL, NYOR_ALL
  INTEGER, DIMENSION(:), ALLOCATABLE :: NXEND_ALL,NYEND_ALL
!
! LBC type
!
! X-direction LBC type at left(1) and right(2) boundaries
!
  CHARACTER(LEN=4), DIMENSION(:, :), ALLOCATABLE :: CLBCX
!
! Y-direction LBC type at left(1) and right(2) boundaries
!
  CHARACTER(LEN=4), DIMENSION(:, :), ALLOCATABLE :: CLBCY
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_DIM_ll
