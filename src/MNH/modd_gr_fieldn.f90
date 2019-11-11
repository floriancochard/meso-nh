!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_GR_FIELD_n
!     ######################
!
!!****  *MODD_GR_FIELD$n* - declaration of prognostic variables related
!!                          to the ground parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     prognostic variables related to the ground. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_GR_FIELDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!	S. Belair   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       29/04/95      
!!      (J.Stein)      15/11/95  add the surface fluxes                
!!      (V.Masson)     10/06/96  add lake fraction
!!       Modification 09/07/97 (M.Georgelin) add directional z0
!!       Modification 09/07/97 (V.Masson) add Subgrid-Scale-Orography
!!       Modification 22/12/97 (V.Masson) add A/S, h/2, town and water variables
!!       Modification 15/03/99 (V.Masson) *** MAJOR MODIFICATION ***
!!                                        add cover types, additional ISBA and town variables.
!!       Modification 01/06/00 (F.Solmon) adaptation for patch approach 
!!                                        +1D for vegetation variable
!!      (V.Masson)        01/04  surface externalisation (cleaning)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE


REAL, DIMENSION(:,:), POINTER :: XSSO_STDEV=>NULL()      ! standard deviation of
                                                         ! Subgrid Scale Orography
REAL, DIMENSION(:,:), POINTER :: XSSO_ANISOTROPY=>NULL() ! anisotropy of S.S.O.
REAL, DIMENSION(:,:), POINTER :: XSSO_DIRECTION=>NULL()  ! direction of S.S.O.
REAL, DIMENSION(:,:), POINTER :: XSSO_SLOPE=>NULL()      ! slope of S.S.O.
REAL, DIMENSION(:,:), POINTER :: XAVG_ZS=>NULL()         ! averaged orography
REAL, DIMENSION(:,:), POINTER :: XSIL_ZS=>NULL()         ! silhouette orography
REAL, DIMENSION(:,:), POINTER :: XMIN_ZS=>NULL()         ! minimum orography
REAL, DIMENSION(:,:), POINTER :: XMAX_ZS=>NULL()         ! maximum orography

CONTAINS

SUBROUTINE GR_FIELD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!GR_FIELD_MODEL(KFROM)%XSSO_STDEV=>XSSO_STDEV !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XSSO_ANISOTROPY=>XSSO_ANISOTROPY !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XSSO_DIRECTION=>XSSO_DIRECTION !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XSSO_SLOPE=>XSSO_SLOPE !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XAVG_ZS=>XAVG_ZS !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XSIL_ZS=>XSIL_ZS !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XMIN_ZS=>XMIN_ZS !Done in FIELDLIST_GOTO_MODEL
!GR_FIELD_MODEL(KFROM)%XMAX_ZS=>XMAX_ZS !Done in FIELDLIST_GOTO_MODEL
!
! Current model is set to model KTO
!XSSO_STDEV=>GR_FIELD_MODEL(KTO)%XSSO_STDEV !Done in FIELDLIST_GOTO_MODEL
!XSSO_ANISOTROPY=>GR_FIELD_MODEL(KTO)%XSSO_ANISOTROPY !Done in FIELDLIST_GOTO_MODEL
!XSSO_DIRECTION=>GR_FIELD_MODEL(KTO)%XSSO_DIRECTION !Done in FIELDLIST_GOTO_MODEL
!XSSO_SLOPE=>GR_FIELD_MODEL(KTO)%XSSO_SLOPE !Done in FIELDLIST_GOTO_MODEL
!XAVG_ZS=>GR_FIELD_MODEL(KTO)%XAVG_ZS !Done in FIELDLIST_GOTO_MODEL
!XSIL_ZS=>GR_FIELD_MODEL(KTO)%XSIL_ZS !Done in FIELDLIST_GOTO_MODEL
!XMIN_ZS=>GR_FIELD_MODEL(KTO)%XMIN_ZS !Done in FIELDLIST_GOTO_MODEL
!XMAX_ZS=>GR_FIELD_MODEL(KTO)%XMAX_ZS !Done in FIELDLIST_GOTO_MODEL

END SUBROUTINE GR_FIELD_GOTO_MODEL

END MODULE MODD_GR_FIELD_n
