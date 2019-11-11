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
!     ###########################
      MODULE MODN_PARAM_KAFR_n
!     ###########################
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAM_KAFR_n, ONLY: &
         XDTCONV_n => XDTCONV, &
         LREFRESH_ALL_n => LREFRESH_ALL, &
         LCHTRANS_n => LCHTRANS, &
         LDOWN_n => LDOWN, &
         NICE_n => NICE, &
         LSETTADJ_n => LSETTADJ, &
         XTADJD_n => XTADJD, &
         XTADJS_n => XTADJS, &
         LDIAGCONV_n => LDIAGCONV, &
         NENSM_n => NENSM
!
IMPLICIT NONE
!
REAL, SAVE  :: XDTCONV
LOGICAL, SAVE  :: LREFRESH_ALL
LOGICAL, SAVE  :: LCHTRANS
LOGICAL, SAVE  :: LDOWN
INTEGER, SAVE  :: NICE
LOGICAL, SAVE  :: LSETTADJ
REAL, SAVE  :: XTADJD
REAL, SAVE  :: XTADJS
LOGICAL, SAVE  :: LDIAGCONV
INTEGER, SAVE  :: NENSM
!
NAMELIST/NAM_PARAM_KAFRn/XDTCONV,LREFRESH_ALL,LCHTRANS,&
                            LDOWN,NICE,LSETTADJ,XTADJD,XTADJS,&
                            LDIAGCONV,NENSM
!
CONTAINS
!
SUBROUTINE INIT_NAM_PARAM_KAFRn
  XDTCONV = XDTCONV_n
  LREFRESH_ALL = LREFRESH_ALL_n
  LCHTRANS = LCHTRANS_n
  LDOWN = LDOWN_n
  NICE = NICE_n
  LSETTADJ = LSETTADJ_n
  XTADJD = XTADJD_n
  XTADJS = XTADJS_n
  LDIAGCONV = LDIAGCONV_n
  NENSM = NENSM_n
END SUBROUTINE INIT_NAM_PARAM_KAFRn

SUBROUTINE UPDATE_NAM_PARAM_KAFRn
  XDTCONV_n = XDTCONV
  LREFRESH_ALL_n = LREFRESH_ALL
  LCHTRANS_n = LCHTRANS
  LDOWN_n = LDOWN
  NICE_n = NICE
  LSETTADJ_n = LSETTADJ
  XTADJD_n = XTADJD
  XTADJS_n = XTADJS
  LDIAGCONV_n = LDIAGCONV
  NENSM_n = NENSM
END SUBROUTINE UPDATE_NAM_PARAM_KAFRn

END MODULE MODN_PARAM_KAFR_n
