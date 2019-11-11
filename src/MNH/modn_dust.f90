!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!     ######spl
      MODULE MODN_DUST
!!    #####################
!!
!!*** *MODN_DUST*
!!
!!    PURPOSE
!!    -------
!       Namelist for desertic dust scheme parameters 
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
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_DUST
IMPLICIT NONE
!
NAMELIST /NAM_DUST/ LDUST, CRGUNITD, LVARSIG, LSEDIMDUST, XN0MIN, XINIRADIUS, &
                    XINISIG, XSIGMIN, XSIGMAX, XCOEFRADMAX, XCOEFRADMIN,      &
                    NMODE_DST, LRGFIX_DST, LDEPOS_DST
!
END MODULE MODN_DUST
