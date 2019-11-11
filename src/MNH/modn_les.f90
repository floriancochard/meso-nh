!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/10/17 13:46:43
!-----------------------------------------------------------------
!     ####################
      MODULE MODN_LES
!     ###################
!
!!****  *MODN_LES* - declaration of namelist NAM_LES
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_LES
!     which controls the averages used in the LES mode.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_LES : 
!! 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_TURBn)
!!          
!!    AUTHOR
!!    ------
!!	    J. Cuxart and J. Stein     * I.N.M. and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Sept  25,1995
!!     (V. Masson)  Feb   03,2000  LES & //
! ------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_LES
!
IMPLICIT NONE
!
NAMELIST/NAM_LES/LLES_MEAN, LLES_RESOLVED, LLES_SUBGRID,                      &
                 LLES_UPDRAFT, LLES_DOWNDRAFT, LLES_SPECTRA,                  &
                 NLES_LEVELS, XLES_ALTITUDES,                                 &
                 NSPECTRA_LEVELS, XSPECTRA_ALTITUDES,                         &
                 NLES_TEMP_SERIE_I, NLES_TEMP_SERIE_J, NLES_TEMP_SERIE_Z,     &
                 CLES_NORM_TYPE, CBL_HEIGHT_DEF,                              &
                 XLES_TEMP_SAMPLING, XLES_TEMP_MEAN_START, XLES_TEMP_MEAN_END,&
                 XLES_TEMP_MEAN_STEP,                                         &
                 LLES_CART_MASK,                                              &
                 NLES_IINF, NLES_ISUP, NLES_JINF, NLES_JSUP,                  &
                 LLES_NEB_MASK, LLES_CORE_MASK, LLES_MY_MASK, LLES_CS_MASK,   &
                 NLES_MASKS_USER
!
NAMELIST/NAM_PDF/LLES_PDF, NPDF,                                  &
                 XTH_PDF_MIN, XTH_PDF_MAX, XW_PDF_MIN,            &
                 XW_PDF_MAX, XTHV_PDF_MIN,                        &
                 XTHV_PDF_MAX, XRV_PDF_MIN,                       &
                 XRV_PDF_MAX, XRC_PDF_MIN, XRC_PDF_MAX,           &
                 XRR_PDF_MIN, XRR_PDF_MAX,                        &
                 XRI_PDF_MIN, XRI_PDF_MAX, XRS_PDF_MIN,           &
                 XRS_PDF_MAX,XRG_PDF_MIN,                         &
                 XRG_PDF_MAX, XRT_PDF_MIN, XRT_PDF_MAX,           &
                 XTHL_PDF_MIN, XTHL_PDF_MAX
!
END MODULE MODN_LES
