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
!     ###################
      SUBROUTINE MKFLAG_n
!     ###################
!
!!****  *MKFLAG_n * - subroutine to flag surface variables.
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!! where data has no signification, it is replaced by XUNDEF
!!
!!
!!    EXTERNAL
!!    --------
!!
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson    * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     15/03/99
!       F.solmon      06/00 adaptation for patch approach
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!*       0.2    Declarations of local variables for print on FM file
!
!-------------------------------------------------------------------------------
!
      RETURN
!-------------------------------------------------------------------------------
!
END SUBROUTINE MKFLAG_n
