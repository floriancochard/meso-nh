!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######################
      MODULE MODD_ELEC_FLASH
!     ######################
!
!!****  *MODD_ELEC_FLASH* - declaration of flash map arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define the storage of the
!     flash maps
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!       
!!    AUTHOR
!!    ------
!!	J.-P. Pinty *Laboratoire Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/11/13
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!* data records: accumulated flash locations
!
INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: NMAP_TRIG_IC   ! X-Y trig point
INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: NMAP_IMPACT_CG ! X-Y ground impact
INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: NMAP_2DAREA_IC ! IC: 2D max extent
INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: NMAP_2DAREA_CG ! CG: 2D max extent
INTEGER, DIMENSION(:,:,:), SAVE, ALLOCATABLE :: NMAP_3DIC ! reached gridpoints
INTEGER, DIMENSION(:,:,:), SAVE, ALLOCATABLE :: NMAP_3DCG ! reached gridpoints
!
!------------------------------------------------------------------------------
!
END MODULE MODD_ELEC_FLASH
