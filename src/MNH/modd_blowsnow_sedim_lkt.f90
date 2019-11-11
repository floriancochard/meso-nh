!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!!     ########################
       MODULE MODD_BLOWSNOW_SEDIM_LKT
!!     ########################
!!
!!     PURPOSE
!!     -------
!!
!! Purpose: Contains look up tables for settling velocity of drifitng snow particles
!! The parameters to be looked up are: 
!! 1) Number-averaged settling velocity
!! 2) Mass-averaged settling velocity
!! 
!! All values are pre-calculated using matlab. 
!! 
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     Based on MODD_DUST_OPT_LKT (Pierre Tulet)
!!
!!
!!     AUTHOR
!!     ------
!!     Vincent VIONNET (CNRM)
!!
!!
!!     MODIFICATIONS
!!     -------------
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER    :: NMAX_RADIUS_LKT=196 !Max number of radii in look up tables ()
  INTEGER, PARAMETER    :: NMAX_PRESSURE_LKT=4   !Max number of pressure in lkt

  !Declaration of the look up tables 
  REAL, DIMENSION(NMAX_RADIUS_LKT,NMAX_PRESSURE_LKT)             :: XNUMB_SPEED_LKT
  REAL, DIMENSION(NMAX_RADIUS_LKT,NMAX_PRESSURE_LKT)             :: XMASS_SPEED_LKT

  !Declaration of the max and min values taken into account in the tables
  REAL, PARAMETER      :: XRADIUS_LKT_MIN = 5         ![um] smallest number median radius taken into account
  REAL, PARAMETER      :: XRADIUS_LKT_MAX = 200       ![um] largest number median radius taken into account
  REAL, PARAMETER      :: XPRESSURE_LKT_MIN = 45000   ![Pa] smallest pressure coefficient taken into account
  REAL, PARAMETER      :: XPRESSURE_LKT_MAX = 105000  ![Pa] largest pressure coefficient taken into account

END MODULE MODD_BLOWSNOW_SEDIM_LKT
