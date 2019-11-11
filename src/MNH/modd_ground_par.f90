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
!     ######################
      MODULE MODD_GROUND_PAR
!     ######################
!
!!****  *MODD_GROUND_PAR* - declaration of prognostic variables related
!!                          to the ground parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization ISBA. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_GR_FIELDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!	S. Belair   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       29/04/95                      
!!      (V.Masson)     05/10/98 add XCDZ0EFF, XZ0_O_Z0H, XRHOSMIN, XRHOSMAX
!!      (V.Masson)     15/03/99 add number of layers
!       (F.solmon)     01/06/00 add number of patch
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, PARAMETER    :: JPCOVER = 255
! Maximum number of cover classes
!
REAL, PARAMETER       :: XBETAS = 0.408
!
!
REAL, PARAMETER       :: XLAMI = 2.22
!                        conductivity of ice
!
REAL, PARAMETER       :: XWCRN = 10.0
!                        critical value of the equivalent water content
!                        of the snow reservoir
!
REAL, PARAMETER       :: XZ0SN = 1.E-3
REAL, PARAMETER       :: XEMISSN = 1.0
!                        roughness length and emissivity of snow
REAL, PARAMETER       :: XALBNIRWAT =  0.20
!                        water near-infra-red albedo
!
REAL, PARAMETER       :: XALBVISWAT =  0.07
!                        water visible albedo
!
REAL, PARAMETER       :: XEMISWAT =  0.98
!                        water emissivity
!
!
REAL, PARAMETER       :: XEMISSOIL = 0.94
!                        bare soil emissivity
!
REAL, PARAMETER       :: XEMISVEG = 0.97
!                        vegetation emissivity
!
!
REAL, PARAMETER       :: XANSMIN = 0.50
REAL, PARAMETER       :: XANSMAX = 0.85
!                        minimum and maximum values of the albedo of snow
!
REAL, PARAMETER       :: XTODRY = 0.008
!

REAL, PARAMETER       :: XCDZ0EFF = 0.8
! drag coefficient in z0eff computation
!
REAL, PARAMETER       :: XZ0_O_Z0H = 10.
! ratio between roughness lengths for momentum and heat
!
REAL                  :: XZ0REL_MIN = 0.001
!
REAL, PARAMETER       :: XRHOSMIN = 100.
REAL, PARAMETER       :: XRHOSMAX = 300.
! minimum and maximum values of the density of snow
!
REAL, PARAMETER       :: XLAISH = 6.
! LAI over vegetation height ratio for general C3 type and C3 grassland
!
REAL, PARAMETER       :: XWGMIN   = 0.001   ! (m3/m3)
! minimum allowable volumetric liquid water content of soil
!
!
! Peters-Lidard et al. (JAS, 1998) from method of Johanssen (1975)
! thermal conductivity (option) parameters:
!
REAL, PARAMETER       :: XSPHSOIL  = 733.   ! J/(kg K) Soil specific heat
REAL, PARAMETER       :: XDRYWGHT  = 2700.0 ! kg/m3    Soil solids dry weight
REAL, PARAMETER       :: XCONDQRTZ = 7.7    ! W/(m K)  Quartz thermal conductivity
REAL, PARAMETER       :: XCONDOTH1 = 2.0    ! W/(m K)  Other thermal conductivity
REAL, PARAMETER       :: XCONDOTH2 = 3.0    ! W/(m K)  Other thermal conductivity
REAL, PARAMETER       :: XCONDWTR  = 0.57   ! W/(m K)  Water thermal conductivity
!
INTEGER                                :: NVEGTYPE
! number of vegetation types
!
INTEGER                                :: NVT_C4
INTEGER                                :: NVT_TREE
INTEGER                                :: NVT_CONI
INTEGER                                :: NVT_EVER
INTEGER                                :: NVT_GRAS
INTEGER                                :: NVT_TROG
INTEGER                                :: NVT_C3
INTEGER                                :: NVT_NO
INTEGER                                :: NVT_ROCK
INTEGER                                :: NVT_SNOW
INTEGER                                :: NVT_IRR
INTEGER                                :: NVT_PARK
! indexes of these types (C4 cultures types, Broadleaf TREEs, CONIferous trees,
!                         EVERgreen broadleaf trees, GRASsland, TROpical Grass,
!                         C3 CULTures types, NO vegetation, ROCKs,
!                         IRRigated crops, irrigated PARKs gardens or peat bogs)
!
!
INTEGER                                :: NPATCH
! number of different vegetation patches  considered in mosaic approach, patches
! are based on VEGTYPE classification.
!
!
INTEGER                                :: NGROUND_LAYER
! number of ground layers

INTEGER                                :: NROOT_LAYER
! index or root base layer
!
!
INTEGER                                :: NROOF_LAYER
! number of roof layers
!
INTEGER                                :: NROAD_LAYER
! number of road layers
!
INTEGER                                :: NWALL_LAYER
! number of wall layers
!
!

END MODULE MODD_GROUND_PAR
