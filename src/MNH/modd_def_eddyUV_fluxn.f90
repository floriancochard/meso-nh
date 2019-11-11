!MNH_LIC Copyright 2004-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##############################
      MODULE MODD_DEF_EDDYUV_FLUX_n
!     ##############################
!
!!****  *MODDB_DEF_EDDYUV_FLUX_n* - declaration FLUX V'U'
!!
!!    PURPOSE
!!    -------
!! To write fluxes in FM FILE
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    AUTHOR
!!    ------
!!	P.Peyrille 18/02/04
!!
!!    MODIFICATIONS
!!    -------------
!!    05/05/09 M.Tomasini Grid-nesting
!!    25/06/11 M.Tomasini Add a source term for the one_way grid-nesting
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), POINTER :: XVU_FLUX_M=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRVS_EDDY_FLUX=>NULL()
!
CONTAINS
!
SUBROUTINE EDDYUV_FLUX_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
END SUBROUTINE EDDYUV_FLUX_GOTO_MODEL
!
END MODULE MODD_DEF_EDDYUV_FLUX_n
