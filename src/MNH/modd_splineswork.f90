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
!     ###################
      MODULE MODD_SPLINESWORK
!     ###################
!
!!****  *MODD_SPLINESWORK* - declaration of work arrays
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
!!          
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/05/96               
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*        0.1    local variables
!                ---------------
!
REAL,    DIMENSION(:),   ALLOCATABLE   :: W
INTEGER, DIMENSION(:,:), ALLOCATABLE   :: IW
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_SPLINESWORK
