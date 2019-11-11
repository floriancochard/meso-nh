!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2007/02/19 11:21:57
!-----------------------------------------------------------------
!     ######################
       MODULE MODI_INI_CLOUD 
!     ######################
!
INTERFACE
      SUBROUTINE INI_CLOUD ( PTSTEP,  PDZMIN, KSPLITR )
!
INTEGER,                 INTENT(OUT):: KSPLITR   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
!
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
END SUBROUTINE INI_CLOUD
!
END INTERFACE
!
END MODULE MODI_INI_CLOUD
!
!
!
!     ###############################################
      SUBROUTINE INI_CLOUD ( PTSTEP, PDZMIN, KSPLITR )
!     ###############################################
!
!!****  *INI_CLOUD * - initialize the constants necessary for the resolved
!!                         microphysics
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constant necessary 
!!    for the resolved cloud and precipitation computation. The number of small
!!    time steps leading to stable scheme for the rain sedimentation is also
!!    computed (time-splitting technic).     
!!
!!**  METHOD
!!    ------
!!      
!!      The constants are set to their numerical values and the number of small
!!    time is computed by dividing the 2* Deltat time interval of the Leap-frog
!!    scheme so that the stability criterion for the rain sedimentation is
!!    fulfilled for a Raindrop maximal fall velocity equal VTRMAX.  
!!
!!    EXTERNAL
!!    --------
!!      GAMMA    :  gamma function
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XP00                 ! Reference pressure
!!        XRD                  ! Gaz constant for dry air
!!        XRHOLW               ! Liquid water density
!!      Module MODD_REF
!!        XTHVREFZ             ! Reference virtual pot.temp. without orography
!!      Module MODD_CLOUDPAR
!!        XCEXVT               ! constant in the rain drop fall velocity
!!        XC1RC, XC2RC         ! constants for autoconversion
!!        XCEXRA, XCRA         ! constants for accretion     
!!        XCEXRE, XC1RE, XC2RE ! constants for rain evaporation 
!!        XCEXRS, XCRS         ! constants for rain sedimentation
!!        XDIVA                ! vapor diffusivity in air
!!        XTHCO                ! thermal conductivity 
!!      Module MODD_PARAMETERS
!!        JPVEXT               !
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_CLOUD )
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/12/94 
!!      (J.Stein)   30/06/95  use 2*PTSTEP to compute the number of small 
!!                            timesteps for the rain sedimentation
!!      (N. Asencio) 11/08/98 parallel code: PDZMIN is computed outside in ini_modeln
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_REF
USE MODD_CLOUDPAR
USE MODD_PARAMETERS
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                 INTENT(OUT):: KSPLITR   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
!
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
!
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB                ! Coordinates of the first  physical
                              ! points along z
REAL :: ZN0                   ! Raindrop spectrum parameter
REAL :: ZA, ZB                ! Constants in the raindrop fall velocity formula 
REAL :: ZGBP3, ZGBP5O2, ZGBP4 ! Values of gamma function
REAL :: ZF                    ! Ventilation factor coefficient 
REAL :: ZNU                   ! Kinematic viscosity
REAL :: ZVTRMAX               ! Raindrop maximal fall velocity
REAL :: ZZ, ZZR, ZT           ! Work variables
!
!  
!-------------------------------------------------------------------------------
!
!*       1.     RAINDROP SPECTRUM          
!	        -----------------
!
ZN0 = 1.E7
ZZ = XPI * XRHOLW * ZN0
!
!
!-------------------------------------------------------------------------------
!
!*       2.     RAINDROP FALL VELOCITY                    
!	        ----------------------
!
XCEXVT = 0.4
IKB = 1 + JPVEXT
ZZR = (XP00/(XRD*XTHVREFZ(IKB))) ** XCEXVT
ZA = 842.
ZB = 0.8
ZGBP3 = GAMMA(ZB+3.)
ZGBP5O2 = GAMMA((ZB+5.)/2.)
ZGBP4 = GAMMA(ZB+4.)
!
!
!-------------------------------------------------------------------------------
!
!*       3.     AUTOCONVERSION                   
!	        --------------
!
XC1RC = 1.E-3
XC2RC = 0.5E-3
!
!
!-------------------------------------------------------------------------------
!
!*       4.     ACCRETION
!	        ---------
!
XCEXRA = (ZB+3.)/4.
XCRA = 0.25 * XPI * ZA * ZN0 *ZGBP3 * ZZR / (ZZ ** XCEXRA)
!
!-------------------------------------------------------------------------------
!
!*       5.     RAIN EVAPORATION
!	        ----------------
!
ZF = 0.22
ZNU = 0.15E-4
XDIVA = 226.E-7
XTHCO = 24.3E-3
XCEXRE = (ZB+5.)/8.
XC1RE = 2. * XPI * ZN0 / SQRT(ZZ)
XC2RE = 2. * XPI * ZN0 * ZF * SQRT(ZA/ZNU) * SQRT(ZZR) * ZGBP5O2  / (ZZ**XCEXRE)
!
!-------------------------------------------------------------------------------
!
!*       6.     RAIN SEDIMENTATION
!	        ------------------
!
XCEXRS = 1. + ZB /4.
XCRS = (ZA/6.) * ZGBP4 * ZZR / (ZZ ** (ZB/4.))
!
!*       6.1    Set the raindrop maximum fall velocity
!               --------------------------------------
ZVTRMAX = 10.                          
!
!
!*       6.2    Compute the number of small time step integration
!               -------------------------------------------------
KSPLITR = 1
SPLIT : DO
  ZT = PTSTEP / FLOAT(KSPLITR)
  IF ( ZT * ZVTRMAX / PDZMIN .LT. 1.) EXIT SPLIT
  KSPLITR = KSPLITR + 1
END DO SPLIT
!
!
END SUBROUTINE INI_CLOUD
