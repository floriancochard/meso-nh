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
!     ####################
      MODULE MODN_SERIES_n
!     ####################
!****  *MODN_SERIES$n* - declaration of namelist NAM_SERIESn
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_SERIESn
!     which concern the diagnostics for diachro files.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_SERIES$n : contains declaration of some parameters for 
!! diagnostics to store in  diachro files
!!         NKCLA,NKCLS   : K level respectively in CLS and in CLA
!!         NKLOW,NKMID,NKUP : K levels  in the mid troposphere
!!                           ( average are done between NKLOW and NKUP)
!!         NBJSLICE : Number of y-slices for (x,t) series
!!         NIBOXL,NJBOXL : Lower indices of the horizontal box 
!!         NIBOXH,NJBOXH : Higher indices of the horizontal box  
!!         NJSLICEL,NJSLICEH :  Lower and  Higher index along y-axe 
!!         of the y-slice
!!         NFREQSERIES   : time freqency    of diagnostic writting       
!!
!!    REFERENCE
!!    ---------
!!     none
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/02/98
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SERIES_n, ONLY: &
         NFREQSERIES_n => NFREQSERIES, &
         NIBOXL_n => NIBOXL, &
         NJBOXL_n => NJBOXL, &
         NIBOXH_n => NIBOXH, &
         NJBOXH_n => NJBOXH, &
         NKCLA_n => NKCLA, &
         NKCLS_n => NKCLS, &
         NKLOW_n => NKLOW, &
         NKMID_n => NKMID, &
         NKUP_n => NKUP, &
         NBJSLICE_n => NBJSLICE, &
         NJSLICEL_n => NJSLICEL, &
         NJSLICEH_n => NJSLICEH
!
IMPLICIT NONE
!
INTEGER, SAVE  :: NFREQSERIES
INTEGER, SAVE  :: NIBOXL
INTEGER, SAVE  :: NJBOXL
INTEGER, SAVE  :: NIBOXH
INTEGER, SAVE  :: NJBOXH
INTEGER, SAVE  :: NKCLA
INTEGER, SAVE  :: NKCLS
INTEGER, SAVE  :: NKLOW
INTEGER, SAVE  :: NKMID
INTEGER, SAVE  :: NKUP
INTEGER, SAVE  :: NBJSLICE
INTEGER, SAVE, DIMENSION(20)  :: NJSLICEL
INTEGER, SAVE, DIMENSION(20)  :: NJSLICEH
!
NAMELIST/NAM_SERIESn/NFREQSERIES,NIBOXL,NJBOXL,NIBOXH,NJBOXH,&
                     NKCLA,NKCLS,NKLOW,NKMID,NKUP,NBJSLICE,&
                     NJSLICEL,NJSLICEH
!
CONTAINS
!
SUBROUTINE INIT_NAM_SERIESn
  NFREQSERIES = NFREQSERIES_n
  NIBOXL = NIBOXL_n
  NJBOXL = NJBOXL_n
  NIBOXH = NIBOXH_n
  NJBOXH = NJBOXH_n
  NKCLA = NKCLA_n
  NKCLS = NKCLS_n
  NKLOW = NKLOW_n
  NKMID = NKMID_n
  NKUP = NKUP_n
  NBJSLICE = NBJSLICE_n
  NJSLICEL = NJSLICEL_n
  NJSLICEH = NJSLICEH_n
END SUBROUTINE INIT_NAM_SERIESn

SUBROUTINE UPDATE_NAM_SERIESn
  NFREQSERIES_n = NFREQSERIES
  NIBOXL_n = NIBOXL
  NJBOXL_n = NJBOXL
  NIBOXH_n = NIBOXH
  NJBOXH_n = NJBOXH
  NKCLA_n = NKCLA
  NKCLS_n = NKCLS
  NKLOW_n = NKLOW
  NKMID_n = NKMID
  NKUP_n = NKUP
  NBJSLICE_n = NBJSLICE
  NJSLICEL_n = NJSLICEL
  NJSLICEH_n = NJSLICEH
END SUBROUTINE UPDATE_NAM_SERIESn

END MODULE MODN_SERIES_n
