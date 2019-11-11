!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_real 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #####################
      MODULE MODI_WATER_SUM
!     #####################
INTERFACE
      FUNCTION WATER_SUM(PR) RESULT (PSUM_R)
!
REAL,DIMENSION(:,:,:,:),              INTENT(IN) :: PR     ! water species
REAL,DIMENSION(SIZE(PR,1),SIZE(PR,2),SIZE(PR,3)) :: PSUM_R ! sum of water species
!
END FUNCTION WATER_SUM
END INTERFACE
END MODULE MODI_WATER_SUM
!     ######################################
      FUNCTION WATER_SUM(PR) RESULT (PSUM_R)
!     ######################################
!
!!****  *WATER_SUM* - summation on the water species (fourth array index)
!! 
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/10/98
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL,DIMENSION(:,:,:,:),              INTENT(IN) :: PR     ! water species
REAL,DIMENSION(SIZE(PR,1),SIZE(PR,2),SIZE(PR,3)) :: PSUM_R ! sum of water species
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IRR      ! number of water species
INTEGER :: JRR      ! loop on water species
!-------------------------------------------------------------------------------
!
IRR = SIZE(PR,4)
!
PSUM_R(:,:,:) = 0.
!
DO JRR=1,IRR
  PSUM_R(:,:,:) = PSUM_R(:,:,:) + PR(:,:,:,JRR)
END DO
!-------------------------------------------------------------------------------
!
END FUNCTION WATER_SUM
