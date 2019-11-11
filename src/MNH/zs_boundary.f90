!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2006/05/18 13:07:25
!-----------------------------------------------------------------
!###################
MODULE MODI_ZS_BOUNDARY
!###################
INTERFACE
      SUBROUTINE ZS_BOUNDARY(PZS,PZS_LS)
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PZS     ! fine orography
REAL, DIMENSION(:,:), INTENT(IN)   :: PZS_LS  ! coarse orography
!
END SUBROUTINE ZS_BOUNDARY
END INTERFACE
END MODULE MODI_ZS_BOUNDARY
!     ##################################
      SUBROUTINE ZS_BOUNDARY(PZS,PZS_LS)
!     ##################################
!
!!****  *ZS_BOUNDARY* - fill zs on boundary
!! 
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!     This routine called from 
!!    READ_ALL_DATA_MESONH_CASE
!!    SPAWN_ZS (called from SPAWN_GRID2 and FILL_ZSMTn
!!     ensures that fine orography is set to coarse orography on non-physical points.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS
!!         JPHEXT
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
!!      Original    12/01/05
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PZS     ! fine orography
REAL, DIMENSION(:,:), INTENT(IN)   :: PZS_LS  ! coarse orography
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IIB,IJB,IIE,IJE
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION
!
IIB = 1+JPHEXT
IJB = 1+JPHEXT
IIE = SIZE(PZS,1) - JPHEXT
IJE = SIZE(PZS,2) - JPHEXT 
!
!-------------------------------------------------------------------------------
!
!*            2. 
!
PZS(IIB-1,:) = PZS_LS(IIB-1,:)
PZS(IIE+1,:) = PZS_LS(IIE+1,:)
PZS(:,IJB-1) = PZS_LS(:,IJB-1)
PZS(:,IJE+1) = PZS_LS(:,IJE+1)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ZS_BOUNDARY
