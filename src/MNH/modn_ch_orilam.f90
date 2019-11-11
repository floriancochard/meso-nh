!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/10/18 12:10:16
!-----------------------------------------------------------------
!!    #####################
      MODULE MODN_CH_ORILAM
!!    #####################
!!
!!*** *MODN_CH_ORILAM*
!!
!!    PURPOSE
!!    -------
!       Namelist for ORILAM aerosol scheme parameters 
!!
!!**  AUTHOR
!!    ------
!!    P. Tulet      *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/02/05
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_AEROSOL, ONLY: LORILAM, XN0IMIN, XN0JMIN, LSEDIMAERO, LAERINIT,&
                          LHETEROSO4, CNUCLEATION, XINISIGI, XINISIGJ,  & 
                          XINIRADIUSI, XINIRADIUSJ, LVARSIGI,&
                          LVARSIGJ, CMINERAL, CORGANIC,&
                          XSIGIMIN, XSIGIMAX,XSIGJMIN, XSIGJMAX,  & 
                          XCOEFRADIMAX, XCOEFRADIMIN, XCOEFRADJMAX, XCOEFRADJMIN,&
                          CRGUNIT, LRGFIX, LDEPOS_AER
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
NAMELIST /NAM_CH_ORILAM/  LORILAM, XN0IMIN, XN0JMIN, LSEDIMAERO, LAERINIT, &
                          LHETEROSO4, CNUCLEATION, XINISIGI, XINISIGJ,  & 
                          XINIRADIUSI, XINIRADIUSJ, LVARSIGI,&
                          LVARSIGJ, CMINERAL, CORGANIC, &
                          XSIGIMIN, XSIGIMAX,XSIGJMIN, XSIGJMAX,  & 
                          XCOEFRADIMAX, XCOEFRADIMIN, XCOEFRADJMAX, XCOEFRADJMIN,&
                          CRGUNIT, LRGFIX, LDEPOS_AER

!
END MODULE MODN_CH_ORILAM
