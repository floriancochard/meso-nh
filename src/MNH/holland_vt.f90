!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_HOLLAND_VT
!     ######################
INTERFACE
!
      SUBROUTINE HOLLAND_VT(PZMAX,PZ,PR_M,PCORIO,PVT)
!
REAL, INTENT(IN)        :: PZMAX! Altitude for which the wind vanishes (m)
REAL, DIMENSION(:), INTENT(IN) :: PZ ! Altitude of Meso-NH's point (m)
REAL, INTENT(IN)        :: PR_M ! Distance from the considered point to the 
                                !center of circulation (meters)
REAL, INTENT(IN)        :: PCORIO ! Coriolis Parameter at the considered point
REAL, DIMENSION(:), INTENT(OUT):: PVT  ! Value of tangential wind Vt 
                                !from the Holland's formulation (m/s)
!
END SUBROUTINE HOLLAND_VT
END INTERFACE
!
END MODULE MODI_HOLLAND_VT
!
!
!
!     #############################################
      SUBROUTINE HOLLAND_VT(PZMAX,PZ,PR_M,PCORIO,PVT)
!     #############################################
!
!!****  *POLAR_MEAN* - This subroutine calculates tangential wind VT_aj(r,z)
!!                          from the Holland (1980) analytical formulation 
!!
!!    PURPOSE
!!    -------
!         This subroutine calculates the tangential wind field VT_aj(r,z)
!           from the Holland's formulation
!
!
!!**  METHOD
!!    ------
!!        The Holland's law describing VT_aj with respect to r and z is
!!       assumed to have the following expression including gradient wind :
!!
!!         VT(r,z) = SQRT ( VTmax(z)**2 * exp {1. - 1./[r/RMW(z)]^B(z)} /
!!                          ([r/RMW(z)]^B(z)) + (fr/2)**2 ) - |fr/2|
!!
!!         with f :Coriolis parameter
!!
!!         The vertical evolution of the maximum tangential wind,
!!                                   the radius of maximum wind,
!!                                   the B parameter.
!!
!!          VTmax(z)= VTmax_0 * ( 1 - z / Zmax)^ C
!!            with : VTmax_0: the maximum tangential wind 
!!                           near the surface or about 500 m altitude
!!                   Zmax: the level for which the tangential wind vanish
!!                   C: a standard coefficient 
!!          RMW(z) = RMW_0 * (1. + rho_z*(z/Zmax) + rho_zz*(z/Zmax)^2)
!!            with : RMW_0 : the radius of maximum wind near the surface 
!!                          or about 500 m altitude
!!                   Zmax: the level for which the tangential wind vanish
!!                   rho_zz and rho_zz: standard coefficients 
!!          B(z) = B_0 * (1. + beta_z*(z/Zmax) + beta_zz*(z/Zmax)^2)
!!            with : B_0 : the B parameter value near the surface
!!                        or about 500 m altitude
!!                   Zmax: the level for which the tangential wind vanish
!!                   beta_zz and beta_zz: standard coefficients  
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      F. Roux               * L.A  *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              10/12/03
!!      BARBARY D (01/2008) : Uses namelists for Holland parameters
!!      BARBARY D (02/2008) : Add Coriolis parameter in Holland formulation 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_HURR_CONF,  ONLY: XVTMAXSURF ,&  ! maximum tangential wind (m/s)
                           XRADWINDSURF,& ! radius of maximum wind (km)
			   XC,&           ! standard coefficient for VTmax(z)
			   XRHO_Z, XRHO_ZZ,& ! standard coefficient for RMW(z)
			   XB_0, XBETA_Z, XBETA_ZZ ! standard coefficient for B(z)
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!
REAL, INTENT(IN)        :: PZMAX! Altitude for which the wind vanishes
REAL, DIMENSION(:), INTENT(IN) :: PZ ! Altitude of Meso-NH's point (m)
REAL, INTENT(IN)        :: PR_M ! Distance from the considered point to the 
                                !center of circulation (meters)
REAL, INTENT(IN)        :: PCORIO ! Coriolis Parameter at the considered point
REAL, DIMENSION(:), INTENT(OUT):: PVT  ! Value of tangential wind Vt 
                                !from the Holland's formulation (m/s)
!
!*       0.2   Local variables
!
REAL, DIMENSION(SIZE(PVT)):: ZVTMAX, ZRMW, ZB
REAL, DIMENSION(SIZE(PVT)):: ZRADIM, ZRADIM_IB
!
!-------------------------------------------------------------------------------
!
!*       1. INITIALIZATIONS
!           ---------------
! Use values in Namelist NAM_HURR_CONF for XC, XHRO_Z, XRHO_ZZ, XBETA_Z, XBETA_ZZ

!
!-------------------------------------------------------------------------------
!
!*       2. COMPUTE THE VALUE OF VT_aj AT A GIVEN RADIUS AND FOR ALL ALTITUDE
!           ---------------
!
ZVTMAX(:) =  XVTMAXSURF * ( 1 - (PZ(:) / PZMAX)) ** XC
ZRMW  (:) = XRADWINDSURF*1000. *(1. + XRHO_Z*(PZ(:)/PZMAX)      &
                                    + XRHO_ZZ*((PZ(:)/PZMAX)**2.) )
ZB    (:) = XB_0 *(1. + XBETA_Z*(PZ(:)/PZMAX) + XBETA_ZZ*((PZ(:)/PZMAX)**2.))
!
ZRADIM(:) = PR_M / ZRMW(:)
ZRADIM_IB(:) = ZRADIM(:) ** ZB(:)
!
!
PVT(:) = SQRT ( ZVTMAX(:)**2 * EXP (1. - 1. / ZRADIM_IB(:) ) / ZRADIM_IB(:) + &
		(PCORIO * PR_M / 2.)**2 ) - ABS(PCORIO * PR_M / 2.)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE HOLLAND_VT
