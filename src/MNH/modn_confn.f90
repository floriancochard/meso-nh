!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##################
      MODULE MODN_CONF_n
!     ##################
!
!!****  *MODN_CONF$n* - declaration of namelist NAM_CONFn
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_CONFn
!     which concern the configuration of one nested  model.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONF$n : contains declaration of configuration variables
!!
!!         LUSERV   : Logical to use rv
!!         LUSERC   : Logical to use rc
!!         LUSERR   : Logical to use rr
!!         LUSERI   : Logical to use ri
!!         LUSERS   : Logical to use rs
!!         LUSERG   : Logical to use rg
!!         LUSERH   : Logical to use rh
!!         LUSECI   : Logical to use ci
!!         NSV_USER : number of user scalar variables
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_CONFn)
!!
!!       
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/06/94                      
!!      J.-P. Pinty 11/04/96  include the ice concentration
!!      D. Gazen    22/01/01  replace NSV by NSV_USER
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_NSV, ONLY: &
         NSV_USER_n => NSV_USER
USE MODD_CONF_n, ONLY: &
         LUSERV_n => LUSERV, &
         LUSERC_n => LUSERC, &
         LUSERR_n => LUSERR, &
         LUSERI_n => LUSERI, &
         LUSERS_n => LUSERS, &
         LUSERG_n => LUSERG, &
         LUSERH_n => LUSERH, &
         LUSECI_n => LUSECI
!
IMPLICIT NONE
!
LOGICAL,SAVE  :: LUSERV
LOGICAL,SAVE  :: LUSERC
LOGICAL,SAVE  :: LUSERR
LOGICAL,SAVE  :: LUSERI
LOGICAL,SAVE  :: LUSERS
LOGICAL,SAVE  :: LUSERG
LOGICAL,SAVE  :: LUSERH
LOGICAL,SAVE  :: LUSECI
INTEGER,SAVE  :: NSV_USER
!
NAMELIST/NAM_CONFn/LUSERV,LUSERC,LUSERR,LUSERI,LUSERS,LUSERG,LUSERH,LUSECI,NSV_USER
!
CONTAINS
!
SUBROUTINE INIT_NAM_CONFn
  LUSERV = LUSERV_n
  LUSERC = LUSERC_n
  LUSERR = LUSERR_n
  LUSERI = LUSERI_n
  LUSERS = LUSERS_n
  LUSERG = LUSERG_n
  LUSERH = LUSERH_n
  LUSECI = LUSECI_n
  NSV_USER = NSV_USER_n
END SUBROUTINE INIT_NAM_CONFn

SUBROUTINE UPDATE_NAM_CONFn
  LUSERV_n = LUSERV
  LUSERC_n = LUSERC
  LUSERR_n = LUSERR
  LUSERI_n = LUSERI
  LUSERS_n = LUSERS
  LUSERG_n = LUSERG
  LUSERH_n = LUSERH
  LUSECI_n = LUSECI
  NSV_USER_n = NSV_USER
END SUBROUTINE UPDATE_NAM_CONFn

END MODULE MODN_CONF_n
