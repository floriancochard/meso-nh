!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!     ##############
      MODULE MODE_ll
!     ##############
!
!!     Purpose
!!     -------
!
!      The purpose of this module is to provide subroutines and functions
!      of the user interface
!
!------------------------------------------------------------------------------
! 
  USE MODD_ARGSLIST_ll
!
  USE MODI_INIT_ll
!
  USE MODI_ADDnDFIELD_ll
  USE MODI_DELnDFIELD_ll
!
  USE MODI_ADDDELFIELD2_ll
!
  USE MODI_UPDATE_ll
  USE MODI_REMAP_ll
!
  USE MODI_SUM_ll 
!
  USE MODI_GET_ll 
  USE MODI_LOCATION_ll 
!
  USE MODI_NEST_ll 
! 
END MODULE MODE_ll
