!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODN_MEAN
!     ###################
!
!!****  *MODN_MEAN* - declaration of namelist NAM_MEAN
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_MEAN
!     which controls the averages used in the MEAN mode.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_MEAN : 
!! 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	    P.Aumond * CNRM *                                           
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    July, 2011         
! ------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_MEAN_FIELD
!
IMPLICIT NONE
!
NAMELIST/NAM_MEAN/LMEAN_FIELD
!
!
END MODULE MODN_MEAN
