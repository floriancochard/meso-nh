!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ################
      MODULE MODN_CONF
!     ################
!
!!****  *MODN_CONF* - declaration of namelist NAM_CONF
!!
!!    PURPOSE
!!    -------
!       The purpose of this  module is to specify  the namelist NAMCONF
!     which concerns the configuration of all models. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONF : contains declaration of configuration variables
!!
!!         CCONF      : configuration of models
!!                          'START' for start configuration
!!                          'RESTA' for restart configuration
!!         LTHINSHELL : Logical for thinshell approximation
!!                          .TRUE.  = thinshell approximation
!!                          .FALSE. = no thinshell approximation
!!         LFLAT      :  Logical for zero ororography
!!                          .TRUE.  = no orography (zs=0.)
!!                          .FALSE. = orography  
!!         NMODEL     : Number of nested models
!!         CEQNSYS    : EQuatioN SYStem resolved by the MESONH model
                                 ! 'LHE' Lipps and HEmler anelastic system
                                 ! 'DUR' approximated form of the DURran version
                                 ! of the anelastic sytem
                                 ! 'MAE' classical Modified Anelastic Equations
                                 ! but with not any approximation in the
                                 ! momentum equation
!!         NVERB      : Level of informations on output-listing
!!                          0 for minimum of prints
!!                          5 for intermediate level of prints
!!                         10 for maximum of prints 
!!         CEXP       : Experiment name
!!         CSEG       : Name of segment
!!         LFORCING   : Logical for forcing
!!                          .TRUE.  = forcing fields available
!!                          .FALSE. = no forcing fields
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_CONF)
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
!!      Original    01/06/94                      
!!      Jabouille P (03/03/97)     suppress LTHINSHELL
!!      Pinty J-P   (03/03/97)     add LFORCING
!!      J. Stein    (25/07/97)     add the equation system switch
!!      Jabouille P (10/05/98)     add LPACK flag
!!      P Jabouille (21/07/99)     add NHALO and CSPLIT
!!      P Jabouille (26/06/01)     lagrangian variable management
!!      V Masson    (03/01/05)     suppress L1D,L2D,LPACK
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CONF
USE MODD_PARAMETERS, ONLY : JPHEXT
!
IMPLICIT NONE
!
NAMELIST/NAM_CONF/CCONF,LFLAT,NMODEL,CEQNSYS,NVERB,CEXP,CSEG,LFORCING, &
                  NHALO,CSPLIT,LLG,LINIT_LG,CINIT_LG,LNOMIXLG,LCHECK, &
                  JPHEXT
!
END MODULE MODN_CONF
