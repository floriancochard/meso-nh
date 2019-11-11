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
!     ######################
      MODULE MODD_DEF_GR_FIELD
!     ######################
!
!!****  *MODD_DEF_GR_FIELD* - declaration of default values of surface pgd fields
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/08/97
!!                     15/03/99 (Masson) add new defaults
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL, SAVE                              :: TS, T2, WG, W2
!                                          initial values of ground prognostic
!                                          variables
!                                          TS = surface temperature
!                                          T2 = deep-soil temperature
!                                          WG = superficial soil-moisture
!                                                  volumetric water content
!                                          W2 = deep-soil moisutre
!                                                  volumetric water content
!
REAL, SAVE                              :: CLAY, SAND, LAND, D2
REAL, SAVE                              :: Z0VEG, Z0HVEG, Z0REL, ALBVIS, ALBNIR 
REAL, SAVE                              :: EMIS, VEG
REAL, SAVE                              :: LAI, RSMIN, GAMMA, RGL, CV, SST
!                                          default values of the corresponding
REAL, SAVE                              :: WSNOW, WGI
REAL, SAVE                              :: TS_TOWN, WS_TOWN, WSNOW_TOWN
REAL, SAVE                              :: D3
REAL, SAVE                              :: RUNOFFB
REAL, SAVE                              :: WG1 ! superficial soil reservoir (kg/kg)
REAL, SAVE                              :: WG2 ! root soil reservoir (kg/kg)
REAL, SAVE                              :: WG3 ! deep soil reservoir (kg/kg)
!
END MODULE MODD_DEF_GR_FIELD
