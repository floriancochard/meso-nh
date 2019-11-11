!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ################
      MODULE MODN_DYN
!     ################
!
!!****  *MODN_DYN* - declaration of namelist NAM_DYN
!!
!!    PURPOSE
!!    -------
!       The purpose of this  module is to specify the namelist NAM_DYN 
!    which contains the  dynamic control variables for all models.    
!     
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DYN : contains the  declaration of dynamic control 
!!                        variables for all models
!!         XSEGLEN  : Duration of segment (in seconds)
!!         XASSELIN : Asselin coefficient
!!         LCORIO   : Switch to take into account the Earth rotation
!!         LNUMDIFF : Switch to apply the numerical diffusion
!!         XALKTOP  : Damping coef. at the top of the absorbing layer
!!         XALZBOT  : Height of the absorbing layer base
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DYN)
!!      Asencio N. et al., 1994, "Le projet de modele non-hydrostatique 
!!    commun CNRM-LA, specifications techniques", Note CNRM/GMME, 26, 139p,
!!    (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      01/06/94
!!      Modifications 17/10/94  (Stein)   For LCORIO                 
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DYN
USE MODD_DYNZD
!
IMPLICIT NONE
! 
NAMELIST/NAM_DYN/XSEGLEN,XASSELIN,XASSELIN_SV,LCORIO,LNUMDIFU,LNUMDIFTH, &
                 LNUMDIFSV,XALKTOP,XALZBOT,LZDIFFU,XALKGRD,XALZBAS
!
END MODULE MODN_DYN
