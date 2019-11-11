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
!!    ############################# 
      MODULE MODN_CH_MODEL0D
!!    #############################
!!
!!*** *MODN_CH_MODEL0D*
!!
!!    PURPOSE
!!    -------
!     contains namelist parameters for the box model
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D, ONLY: XTBEGIN, XTEND, XDTACT,         &
			  XDTOUT, XDTDIAG,                 &
			  CRUNID,                          &
			  CINITFILE, COUTFILE, CMETEOFILE, &
			  CRESULTFILE, CRESULTFORMAT,      &
			  CDIAGFILE, CDIAGFORMAT,          &
			  NVERB
!
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!     
! variables to be put into the namelist
NAMELIST /NAM_CH_MODEL0D/ XTBEGIN, XTEND, XDTACT,          &
			  XDTOUT, XDTDIAG,                 &
			  CRUNID,                          &
			  CINITFILE, COUTFILE, CMETEOFILE, &
			  CRESULTFILE, CRESULTFORMAT,      &
			  CDIAGFILE, CDIAGFORMAT,          &
			  NVERB
!
END MODULE MODN_CH_MODEL0D
