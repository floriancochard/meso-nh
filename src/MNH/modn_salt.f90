!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/12/21 14:53:53
!-----------------------------------------------------------------
!!    #####################
      MODULE MODN_SALT
!!    #####################
!!
!!*** *MODN_SALT*
!!
!!    PURPOSE
!!    -------
!       Namelist for marine sea salt scheme parameters 
!!
!!**  AUTHOR
!!    ------
!!    P. Tulet      *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/02/05
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_SALT
IMPLICIT NONE
!
NAMELIST /NAM_SALT/ LSALT, CRGUNITS, LVARSIG_SLT,LSEDIMSALT,XN0MIN_SLT, XINIRADIUS_SLT, &
               XINISIG_SLT, XSIGMIN_SLT, XSIGMAX_SLT, XCOEFRADMAX_SLT, XCOEFRADMIN_SLT, &
               NMODE_SLT, LRGFIX_SLT, LDEPOS_SLT, LONLY
!
END MODULE MODN_SALT
