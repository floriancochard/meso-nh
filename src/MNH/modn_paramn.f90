!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/18 13:07:53
!-----------------------------------------------------------------
!     ###################
      MODULE MODN_PARAM_n
!     ###################
!
!!****  *MODN_PARAM$n* - declaration of namelist NAM_PARAMn
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_PARAMn
!     which concern the parameterization and cloud physics variables
!     of one nested  model.   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAM$n : contains declaration of parameterization and
!!    cloud physics variables
!! 
!!         CTURB   :  Kind of turbulence parameterization
!!                     'NONE' if no parameterization 
!!         CRAD    :  Kind of radiation parameterization 
!!                     'NONE' if no parameterization 
!!         CGROUND :  Kind of surface processes parameterization 
!!                      'NONE' if no parameterization 
!!         CCLOUD  :  Kind of microphysical scheme
!!                      'NONE' if no parameterization
!!         CDCONV  :  Kind of deep convection scheme
!!                      'NONE' if no parameterization
!!         CSCONV  :  Kind of shallow convection scheme
!!                      'NONE' if no parameterization
!!        CSEA_FLUX:  Kind of algrithm for the fluxes over water
!!                      'DIRECT' for the old algorithm
!!                      'ITERAT' for the new algorithm (cf Fairall et al 96)
!!                               adapted for weak winds
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAMn)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/06/94            
!!      E. Richard  01/06/95   add CCLOUD          
!!      P. Bechtold 26/03/96   add CDCONV
!!      M. Tomasini 11/12/00   add CSEA_FLUX
!!      JP. Pinty   26/11/02   add CELEC
!!      V. Masson   01/2004    removes CGROUND and CSEA_FLUX
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAM_n, ONLY: &
         CTURB_n => CTURB, &
         CRAD_n => CRAD, &
         CCLOUD_n => CCLOUD, &
         CDCONV_n => CDCONV, &
         CSCONV_n => CSCONV, &
         CELEC_n => CELEC,  &
         CACTCCN_n => CACTCCN
!
IMPLICIT NONE
!
CHARACTER (LEN=4),SAVE  :: CTURB
CHARACTER (LEN=4),SAVE  :: CRAD
CHARACTER (LEN=4),SAVE  :: CCLOUD
CHARACTER (LEN=4),SAVE  :: CDCONV
CHARACTER (LEN=4),SAVE  :: CSCONV
CHARACTER (LEN=4),SAVE  :: CELEC
CHARACTER (LEN=4),SAVE  :: CACTCCN
!
NAMELIST/NAM_PARAMn/CTURB,CRAD,CCLOUD,CDCONV,CSCONV,CELEC,CACTCCN
!
CONTAINS
!
SUBROUTINE INIT_NAM_PARAMn
  CTURB = CTURB_n
  CRAD = CRAD_n
  CCLOUD = CCLOUD_n
  CDCONV = CDCONV_n
  CSCONV = CSCONV_n
  CELEC = CELEC_n
  CACTCCN = CACTCCN_n
END SUBROUTINE INIT_NAM_PARAMn

SUBROUTINE UPDATE_NAM_PARAMn
  CTURB_n = CTURB
  CRAD_n = CRAD
  CCLOUD_n = CCLOUD
  CDCONV_n = CDCONV
  CSCONV_n = CSCONV
  CELEC_n = CELEC
  CACTCCN_n = CACTCCN
END SUBROUTINE UPDATE_NAM_PARAMn

END MODULE MODN_PARAM_n
