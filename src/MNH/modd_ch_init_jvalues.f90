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
!!    ###########################
      MODULE MODD_CH_INIT_JVALUES
!!    ###########################
!!
!!*** *MODD_CH_INIT_JVALUES*
!!
!!    PURPOSE 
!!    -------
!!    Store J values variables common to all models 
!!     (i.e. before spatial and temporal interpolation)
!!    XJDATA is calculated at the first call of model 1. 
!!    XJDATA is calculated for a discrete number 
!!     of solar zenith angle, altitude and albedo.
!
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
!
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: XJDATA 
INTEGER                               :: NSZA_INCR = 99 + 1
REAL, ALLOCATABLE, DIMENSION(:)       :: XSZA_JVAL
INTEGER, PARAMETER                    :: NZZ_JVAL = 30 + 1
REAL, ALLOCATABLE, DIMENSION(:)       :: XZZ_JVAL
INTEGER, PARAMETER                    :: JPJVMAX = 42    
INTEGER                               :: NBALB = 10   
!
END MODULE MODD_CH_INIT_JVALUES
