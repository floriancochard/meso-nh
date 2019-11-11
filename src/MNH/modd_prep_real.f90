!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###############
      MODULE MODD_PREP_REAL
!     ###############
!
!!****  *MODD_PREP_REAL* - declaration of work arrays in PREP_REAL_CASE
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
!!
!!    AUTHOR
!!    ------
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   05/05
!!                 05/06 (I.Mallet) add *_SV_* variables to allow chemical
!!                                 initialization from HCHEMFILE
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!* variables allocated in case of Gribex input file
!
REAL                               :: XP00_LS  ! reference pressure in eta
REAL,DIMENSION(:),     ALLOCATABLE :: XA_LS    ! function A in definition of eta
REAL,DIMENSION(:),     ALLOCATABLE :: XB_LS    ! function B in definition of eta
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XT_LS    ! temperature
REAL,DIMENSION(:,:,:,:),ALLOCATABLE:: XQ_LS    ! specific ratio of humidity and
                                               !other hydrometeors
REAL                               :: XP00_SV_LS  ! reference pressure in eta
REAL,DIMENSION(:),     ALLOCATABLE :: XA_SV_LS    ! function A in definition of eta
REAL,DIMENSION(:),     ALLOCATABLE :: XB_SV_LS    ! function B in definition of eta
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XT_SV_LS    ! temperature
REAL,DIMENSION(:,:,:,:),ALLOCATABLE:: XQ_SV_LS    ! specific humidity
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XDUMMY_2D   ! 2D dummy fields read in 
CHARACTER(LEN=16),DIMENSION(:), ALLOCATABLE :: CDUMMY_2D   !GRIBex file
!
!* variables allocated in case of Mesonh input file
!
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XTH_LS   ! potential temperature
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSU_LS  ! large scale pseudo zonal wind component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSV_LS  ! large scale pseudo meridian wind component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSW_LS  ! large scale vertical wind speed
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSTH_LS ! large scale potential temperature
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSRV_LS ! large scale vapor mixing ratios
REAL,DIMENSION(:),     ALLOCATABLE    :: XZHAT_LS ! altitude coordinate
!
LOGICAL                               :: LSLEVE_LS! flag for sleve coordinate
REAL                                  :: XLEN1_LS ! Decay scale for smooth topography
REAL                                  :: XLEN2_LS ! Decay scale for small-scale topography deviation
!
!* variables allocated in both cases
!
REAL,DIMENSION(:,:),   ALLOCATABLE :: XPS_LS   ! surface pressure
REAL,DIMENSION(:,:),   ALLOCATABLE :: XZS_LS   ! orography
REAL,DIMENSION(:,:),   ALLOCATABLE :: XZSMT_LS ! smooth orography
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XZFLUX_LS! altitude of pressure points
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XZMASS_LS! altitude of mass points
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XPMHP_LS ! pressure minus hyd. pressure
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XTHV_LS  ! virtual potential temperature
REAL,DIMENSION(:,:,:,:), ALLOCATABLE :: XR_LS  ! water mixing ratios
REAL,DIMENSION(:,:,:,:), ALLOCATABLE :: XSV_LS ! scalar mixing ratios
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XHU_LS   ! relative humidity
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XTKE_LS  ! turbulence kinetic energy
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XU_LS    ! pseudo zonal wind component
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XV_LS    ! pseudo meridian wind component
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XW_LS    ! vertical wind speed
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XPMASS_LS! pressure of mass points
REAL,DIMENSION(:,:),   ALLOCATABLE :: XPS_SV_LS   ! surface pressure
REAL,DIMENSION(:,:),   ALLOCATABLE :: XZS_SV_LS   ! orography
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XZFLUX_SV_LS! altitude of pressure points
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XZMASS_SV_LS! altitude of mass points
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XPMASS_SV_LS! pressure of mass points
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XTHV_SV_LS  ! virtual potential temperature
REAL,DIMENSION(:,:,:,:), ALLOCATABLE :: XR_SV_LS  ! water mixing ratios
REAL,DIMENSION(:,:,:), ALLOCATABLE :: XHU_SV_LS   ! relative humidity
!
!* variables on MiXed grid
!
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XTHV_MX  ! thetav
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: XR_MX    ! water mixing ratios
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: XSV_MX   ! scalar mixing ratios
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XTKE_MX  ! turbulence kinetic energy
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XPMHP_MX ! pressure minus hyd. pressure
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XU_MX    ! pseudo zonal wind component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XV_MX    ! pseudo meridian wind component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XW_MX    ! vertical wind speed
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSTH_MX ! Large scale pot. temperature
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSRV_MX ! Large scale vapor mixing ratio
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSU_MX  ! Large scale U component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSV_MX  ! Large scale V component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XLSW_MX  ! Large scale W component
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XZFLUX_MX! altitude of the pressure levels
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XZMASS_MX! altitude of the mass levels
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: XRHOD_MX ! rho dry
REAL,DIMENSION(:,:),   ALLOCATABLE    :: XEXNTOP2D  ! local top Exner function
REAL,DIMENSION(:,:),   ALLOCATABLE    :: XPSURF     ! Surface pressure
!
END MODULE MODD_PREP_REAL
