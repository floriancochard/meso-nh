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
!     ##########################
      MODULE MODD_VER_INTERP_LIN
!     ##########################
!
!!****  *MODD_VER_INTERP_LIN* - declaration of linear vertical interpolation
!!                              coefficients
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_VER_INTERP_LIN)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/07/97
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NKLIN      ! level in grid1 just below
                                                     ! level in grid2
                                                     ! (level for which XCOEFLIN
                                                     !  is computed)
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: XCOEFLIN   ! interpolating coefficient
                                                     ! 0< <1 : interpolation
                                                     !    <0 : extrapolation up
                                                     ! 1<    : extrapolation down
!-------------------------------------------------------------------------------
!
END MODULE MODD_VER_INTERP_LIN
