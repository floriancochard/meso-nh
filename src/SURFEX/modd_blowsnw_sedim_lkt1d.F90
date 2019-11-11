
!!     ########################
       MODULE MODD_BLOWSNW_SEDIM_LKT1D
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

  INTEGER, PARAMETER    :: NMAX_RADIUS_LKT1D=196 !Max number of radii in look up tables ()
  INTEGER, PARAMETER    :: NMAX_PRESSURE_LKT1D=4   !Max number of pressure in lkt

  !Declaration of the look up tables 
  REAL, DIMENSION(NMAX_RADIUS_LKT1D,NMAX_PRESSURE_LKT1D)             :: XNUMB_SPEED_LKT1D
  REAL, DIMENSION(NMAX_RADIUS_LKT1D,NMAX_PRESSURE_LKT1D)             :: XMASS_SPEED_LKT1D

  !Declaration of the max and min values taken into account in the tables
  REAL, PARAMETER      :: XRADIUS_LKT1D_MIN = 5         ![um] smallest number median radius taken into account
  REAL, PARAMETER      :: XRADIUS_LKT1D_MAX = 200       ![um] largest number median radius taken into account
  REAL, PARAMETER      :: XPRESSURE_LKT1D_MIN = 45000   ![Pa] smallest pressure coefficient taken into account
  REAL, PARAMETER      :: XPRESSURE_LKT1D_MAX = 105000  ![Pa] largest pressure coefficient taken into account

END MODULE MODD_BLOWSNW_SEDIM_LKT1D
