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
      MODULE MODN_SERIES
!     ##################
!****  *MODN_SERIES* - declaration of namelist NAM_SERIES
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_SERIES
!     which concern the diagnostics for diachro files.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_SERIES : contains declaration of some parameters for 
!! diagnostics to store in  diachro files
!!         LSERIES
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
!!                Oct. 2011 : (P.Le Moigne) Surface series
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SERIES
!
IMPLICIT NONE
!
NAMELIST/NAM_SERIES/LSERIES,LMASKLANDSEA,LWMINMAX,LSURF
!
END MODULE MODN_SERIES
