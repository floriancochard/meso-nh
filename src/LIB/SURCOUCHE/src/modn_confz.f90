!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ################
      MODULE MODN_CONFZ
!     ################
!
!!****  *MODN_CONFZ* - declaration of namelist NAM_CONFZ
!!
!!    PURPOSE
!!    -------
!    configuration of ZSPLITTING 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONFZ : contains declaration of configuration variables
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	J. Escobar L.A.
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/02/2009    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CONFZ
!
IMPLICIT NONE
!
NAMELIST/NAM_CONFZ/ NZ_VERB,NZ_PROC,NB_PROCIO_R,NB_PROCIO_W,MPI_BUFFER_SIZE,LMNH_MPI_BSEND & 
                   ,LMNH_MPI_ALLTOALLV_REMAP,NZ_SPLITTING !JUAN Z_SPLITTING
!
END MODULE MODN_CONFZ
