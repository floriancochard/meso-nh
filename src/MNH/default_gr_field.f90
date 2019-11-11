!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_ideal 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_DEFAULT_GR_FIELD
!     #########################
!
INTERFACE
!
SUBROUTINE DEFAULT_GR_FIELD
END SUBROUTINE DEFAULT_GR_FIELD
!
END INTERFACE
!
END MODULE MODI_DEFAULT_GR_FIELD
!
!
!
!     ########################
      SUBROUTINE DEFAULT_GR_FIELD
!     ########################
!
!!****  *DEFAULT_GR_FIELD * - set default values to variables in EXPRE file
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set default values for the ground
!     variables by filling the corresponding variables stored in modules.
!     These variables  will be stored in LFIFM file, and therefore they have
!     not been set to default values by DEFAULT_DESFM.
!       
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson        * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                                           27/11/96
!!               04/08/97 change the declarative module used (Masson) 
!!               15/03/99 most defaults are removed (because of new
!!                        PGD fields treatment) (Masson)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_DEF_GR_FIELD
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*      1.    SET DEFAULT VALUES FOR MODD_DEF_GR_FIELD :
!             ----------------------------------------
!
TS_TOWN = XUNDEF
WS_TOWN = 0.
WSNOW_TOWN = 0.

TS      = XUNDEF
T2      = XUNDEF
WG1     = XUNDEF 
WG2     = XUNDEF
WG3     = XUNDEF
WGI     = 0.
WSNOW   = 0.

CLAY    = 0.33
SAND    = 0.33

LAND    = 1.0

Z0REL   = 0.001

Z0VEG   = XUNDEF
Z0HVEG  = XUNDEF

ALBVIS  = XUNDEF
ALBNIR  = XUNDEF
EMIS    = XUNDEF

VEG     = XUNDEF
LAI     = XUNDEF
RSMIN   = XUNDEF
GAMMA   = XUNDEF
RGL     = XUNDEF
CV      = XUNDEF
D2      = XUNDEF
D3      = XUNDEF
!
SST     = XUNDEF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_GR_FIELD
