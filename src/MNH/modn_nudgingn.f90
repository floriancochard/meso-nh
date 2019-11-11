!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/29 12:00:43
!-----------------------------------------------------------------
!     ###################
      MODULE MODN_NUDGING_n
!     ###################
!
!!****  *MODN_NUDGING_n* - declaration of namelist NAM_NUDGINGn
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_NUDGINGn 
!!           LNUDGING   ! Logical for nudging term
!!           XTNUDGING  ! Time scale for nudging
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
!!      Original       15/05/06
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_NUDGING_n, ONLY: &
         LNUDGING_n => LNUDGING, &
         XTNUDGING_n => XTNUDGING
!
IMPLICIT NONE
!
LOGICAL, SAVE  :: LNUDGING
REAL,    SAVE  :: XTNUDGING
!
NAMELIST/NAM_NUDGINGn/LNUDGING,XTNUDGING
!
CONTAINS
!
SUBROUTINE INIT_NAM_NUDGINGn
  LNUDGING = LNUDGING_n
  XTNUDGING = XTNUDGING_n
END SUBROUTINE INIT_NAM_NUDGINGn
!
SUBROUTINE UPDATE_NAM_NUDGINGn
  LNUDGING_n = LNUDGING
  XTNUDGING_n = XTNUDGING
END SUBROUTINE UPDATE_NAM_NUDGINGn
!
END MODULE MODN_NUDGING_n
