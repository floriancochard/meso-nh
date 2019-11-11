!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_PGDFIELDS
!     #####################
!
!!****  *MODD_PGDFIELDS* - declaration of physiographic data arrays
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     physiographic data arrays.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PGDFIELDS)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/01/95                      
!!      Modified    22/11/96 (V. Masson) add some primitive orographic fields,
!!                                       used in ADVANCED_PREP_PGD, but not in
!!                                       the model.
!!      Modified    18/02/97 (M. Georgelin) add directional z0
!!      Modified    19/02/97 (V. Masson) add map factor
!!      Modified    15/03/99 (V. Masson) *** MAJOR MODIFICATION ***
!!                                       add cover types, additional ISBA and town fields
!!      Modified    01/06/00 (F.solmon) patch approach +1D for vegetation isba variable 
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSEA       ! subgrid sea fraction
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDLAND      ! subgrid land fraction
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDZS        ! orography
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDZ0REL     ! subgrid-orography roughness length 
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDCLAY      ! clay percentage of the soil
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSAND      ! sand percentage of the soil
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDALBVIS    ! albedo for visible radiation
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDALBNIR    ! albedo for infrared radiation
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDEMIS      ! emissivity
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDZ0VEG     ! roughness length for momentum exchange
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDZ0HVEG    ! roughness length for heat exchange
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDVEG       ! fraction of vegetation
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDLAI       ! leaf area index
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDRSMIN     ! minimum stomatal resistance
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDGAMMA     ! parameter in the calculation of RS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDRGL       ! parameter in the calculation of RS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDCV        ! inverse of the heat capacity of the vegetation
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSST       ! Sea Surface Temperature
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDAVG_ZS    ! averaged orography
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSIL_ZS    ! silhouette orography
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSSO_STDEV     ! standard deviation of
!                                                      ! Subgrid Scale Orography
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDZ0EFFIP,XPGDZ0EFFIM,XPGDZ0EFFJP,XPGDZ0EFFJM
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDAOSIP,XPGDAOSIM,XPGDAOSJP,XPGDAOSJM
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDHO2IP,XPGDHO2IM,XPGDHO2JP,XPGDHO2JM
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSSO_ANISOTROPY ! anisotropy of S.S.O.
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSSO_DIRECTION  ! direction of S.S.O.
REAL, DIMENSION(:,:), ALLOCATABLE :: XPGDSSO_SLOPE      ! slope of S.S.O.
!
LOGICAL, DIMENSION(100)             :: LNOCLASS_PGD0 ! flag to setermine if a
!                                                    ! field was separately
!                                                    ! specified by user in
!                                                    ! PREP_PGD step
!
INTEGER                             :: NPT_USER      ! patch number where
                                                     ! default data found
                                                     ! from cover classes are
                                                     ! replaced by
                                                     ! the fields given by user in prep_pgd.
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDMIN_ZS    ! minimum orography
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDMAX_ZS    ! maximum orography
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDCOVER     ! nature of coverage
                                                     ! 3rd dim = number of classes
!
!* primary surface type parameters
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDTOWN    ! artificial cover
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDNATURE  ! natural and cultivated cover
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDWATER   ! inland waters
!
!* secondary surface type parameters
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDALBNIR_VEG  ! near-infrared albedo
!                                                       ! of artificial surfaces
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDALBVIS_VEG  ! visible albedo
!                                                       ! of artificial surfaces
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_DRY  ! near infra-red and
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_DRY  ! visible albedo
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_WET  ! for dry and wet bare
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_WET  ! soils
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDALBNIR_ECO  !  near infra-red and visible albedo
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDALBVIS_ECO  ! for (prog.)snow-free ecosystems
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDEMIS_ECO    ! emissivity
!                                                       ! of natural surfaces
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_NAT  !  near infra-red and visible albedo
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_NAT  ! for averaged for continental natural surfaces
!
!                                                       ! of artificial surfaces
INTEGER                             :: NPGD_GROUND_LAYER ! number of ground layers
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: XPGDDG          ! layer soil depths
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDRUNOFFB     ! runoff parameter
INTEGER                             :: NPGDVEGTYPE     ! number of veg. types
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDVEGTYPE     ! fraction of each veg.
!                                                      ! type
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDGMES        ! mesophyll conductance (mm s-1)
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDBSLAI       ! ratio d(biomass)/d(lai)
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDLAIMIN      ! minimum LAI
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDSEFOLD      ! e-folding time for senescence (days)
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XPGDH_TREE       ! height of vegetation
!
!
INTEGER                             :: NPGD_ROOF_LAYER  ! number of roof layers
INTEGER                             :: NPGD_ROAD_LAYER  ! number of road layers
INTEGER                             :: NPGD_WALL_LAYER  ! number of wall layers
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDZ0_TOWN  ! town roughness length
                                                    ! for momentum
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDZ0H_TOWN ! town roughness length
                                                    ! for heat
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_TOWN ! near-infrared albedo
                                                       ! of artificial surfaces
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_TOWN ! visible albedo
                                                       ! of artificial surfaces
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDEMIS_TOWN   ! emissivity
                                                       ! of artificial surfaces
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_ROOF ! near-infrared albedo
                                                       ! of roofs
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_ROOF ! visible albedo
                                                       ! of roofs
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDEMIS_ROOF   ! emissivity
                                                       ! of roofs
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDHC_ROOF     ! heat capacity
                                                       ! for roof layers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDTC_ROOF     ! thermal conductivity
                                                       ! for roof layers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDD_ROOF      ! width of roof layers

REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_ROAD ! near-infrared albedo
                                                       ! of roads
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_ROAD ! visible albedo
                                                       ! of roads
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDEMIS_ROAD   ! emissivity
                                                       ! of roads
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDHC_ROAD     ! heat capacity
                                                       ! for road layers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDTC_ROAD     ! thermal conductivity
                                                       ! for road layers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDD_ROAD      ! width of road layers
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDSVF_ROAD    ! sky-view-factor
                                                       ! of roads

REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBNIR_WALL ! near-infrared albedo
                                                       ! of walls
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDALBVIS_WALL ! visible albedo
                                                       ! of walls
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDEMIS_WALL   ! emissivity
                                                       ! of walls
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDHC_WALL     ! heat capacity
                                                       ! for wall layers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDTC_WALL     ! thermal conductivity
                                                       ! for wall layers
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XPGDD_WALL      ! width of wall layers
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDSVF_WALL    ! sky-view-factor
                                                       ! of walls

REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDBLD_HEIGHT   ! building height h
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDBLD_HL_RATIO ! building h/L ratios
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDCAN_HW_RATIO ! canyon h/W ratios
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDBLD          ! fraction of buildings
                                                        ! in artificial areas
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDH_TRAFFIC   ! anthropogenic sensible
!                                                      ! heat fluxes due to traffic
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDLE_TRAFFIC  ! anthropogenic latent
!                                                      ! heat fluxes due to traffic
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDH_INDUSTRY  ! anthropogenic sensible
!                                                      ! heat fluxes due to factories
REAL, DIMENSION(:,:),   ALLOCATABLE :: XPGDLE_INDUSTRY ! anthropogenic latent
!                                                      ! heat fluxes due to factories
!-------------------------------------------------------------------------------
!
END MODULE MODD_PGDFIELDS
