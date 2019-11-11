!MNH_LIC Copyright 1996-2017 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODN_OUTPUT
!     ##################
!
!!****  *MODN_OUTPUT* - declaration of namelist NAM_OUTPUT
!!
!!    PURPOSE
!!    -------
!       The purpose of this  module is to specify the namelist  NAM_OUTPUT
!       which concerns the instants and some parameters (compression and precision reduction)
!       of the outputs realized by all models.
!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_BAKOUT : contains declaration of the variables describing
!!                           the instants and some parameters (compression and
!!                           precision reduction) of the outputs
!!
!!    REFERENCE
!!    ---------
!!      Book2 of Meso-NH documentation (module MODD_BAKOUT)
!!          
!!    AUTHOR
!!    ------
!!	J.P. Lafore      *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/07/96                      
!!      Ph. Wautelet 2016       new structures for outputs/backups
!!      Ph. Wautelet 02/10/2017 split NAM_OUTPUT in NAM_BACKUP and NAM_OUTPUT
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_BAKOUT
!
IMPLICIT NONE
!
NAMELIST/NAM_OUTPUT/LOUT_BEG,LOUT_END,&
                   XOUT_TIME,NOUT_STEP,&
                   NOUT_STEP_FREQ,NOUT_STEP_FREQ_FIRST,&
                   XOUT_TIME_FREQ,XOUT_TIME_FREQ_FIRST, &
                   COUT_VAR, &
                   LOUT_REDUCE_FLOAT_PRECISION, &
                   LOUT_COMPRESS, NOUT_COMPRESS_LEVEL,&
                   COUT_DIR
!
END MODULE MODN_OUTPUT
