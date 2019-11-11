!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_RADAR
!     ################
!
!!****  *MODD_RADAR* - declaration of flags related to the radar diagnostics
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      NONE
!!
!!    AUTHOR
!!    ------
!!	O. Caumont  & V. Ducrocq   * Meteo France *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/11/03 
!!      O. Caumont     14/09/09 Removal of XAZIM
!!      O. Caumont     14/09/09 Possibility to use polar coordinates
!!      C. Augros      2013 Simulator RADAR : add LSNRT XSNRMIN 
!!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL, DIMENSION(30)                       :: XLAT_RAD  ! latitude of radars
REAL, DIMENSION(SIZE(XLAT_RAD))           :: XLON_RAD  ! longitude of radars
REAL, DIMENSION(SIZE(XLAT_RAD))           :: XALT_RAD  ! altitude of radars
CHARACTER(LEN=5),DIMENSION(SIZE(XLAT_RAD)):: CNAME_RAD ! name of radars
REAL, DIMENSION(SIZE(XLAT_RAD))           :: XLAM_RAD  ! radar wavelengths
REAL, DIMENSION(SIZE(XLAT_RAD))           :: XDT_RAD   ! angles at -3 dB of radar (deg)
!
REAL, DIMENSION(SIZE(XLAT_RAD),25)        :: XELEV     !  elevation of source synthetic signal (in deg)
                                                       !  for each radar
REAL, DIMENSION(SIZE(XLAT_RAD))   :: XX_INI    ! x positions of the ground source signal   
REAL, DIMENSION(SIZE(XLAT_RAD))   :: XY_INI    ! y positions of the ground source signal   
REAL, DIMENSION(SIZE(XLAT_RAD))   :: XZ_INI    ! z positions of the ground source signal 
REAL                              :: XSTEP_RAD ! step of integration along the ray beam = gate length
INTEGER                           :: NBRAD     ! number of radars (size of XX_INI)
INTEGER,DIMENSION(SIZE(XLAT_RAD)) :: NBELEV    ! number of elevations computed (sizeof XELEV)
INTEGER                           :: NBAZIM    ! number of azimuths computed
INTEGER                           :: NBSTEPMAX ! maximum number of integration step

INTEGER :: NCURV_INTERPOL ! 0 if effective 4/3 earth radius parametrization 
                          ! 1 if interpolation of P, T, e for computation of curvature of the radar beam
INTEGER :: NDIFF ! 0 : RAYLEIGH ; 1 : MIE ; 2 : T-MATRIX
LOGICAL :: LATT ! FALSE : no attenuation; TRUE: attenuation
LOGICAL :: LCART_RAD ! true if  interpolation of  reflectivity on a cartesian grid ; false if polar
INTEGER :: NPTS_GAULAG ! Only used for NDIFF==1 ; number of points for Gauss-Laguerre quadrature (default )
INTEGER :: NPTS_H
INTEGER :: NPTS_V
LOGICAL :: LQUAD ! false Gauss-Hermite; true Gauss-Legendre
REAL :: XGRID,XVALGROUND=-99999.
INTEGER :: NMAX ! =INT(NBSTEPMAX*XSTEP_RAD/XGRID) (defined in module write_lfifm1_for_diag.f90)
CHARACTER(LEN=5):: CARF ! axis ratio function (PB70=Pruppacher & Beard 1970; AND99=Andsager et al. 1999)
LOGICAL :: LREFR,LDNDZ ! refractivity
INTEGER :: NDGS        !NDGS: number of division points in computing integrals over the surface particles (default=2)
LOGICAL :: LFALL ! take fall speeds into account when simulating Doppler winds
LOGICAL :: LWBSCS ! weighting by backscattering cross sections
LOGICAL :: LWREFL ! weighting by reflectivities
REAL :: XREFLMIN ! min val for reflectivities
REAL :: XREFLVDOPMIN ! min val for Doppler velocities
LOGICAL :: LSNRT !if .TRUE. the threshold on Z and V is function of SNR
REAL :: XSNRMIN !SNR threshold under which reflectivity is set to -XUNDEF if (LSNRT=.TRUE.)



END MODULE MODD_RADAR

