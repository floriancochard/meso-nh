!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_RESOLVED_CLOUD
!     ##########################
INTERFACE
      SUBROUTINE RESOLVED_CLOUD ( HCLOUD, HACTCCN, HSCONV, HMF_CLOUD,                  &
                                  KRR, KSPLITR, KSPLITG, KMI, KTCOUNT,                 &
                                  HLBCX, HLBCY, TPFILE, HRAD, HTURBDIM,                &
                                  OCLOSE_OUT, OSUBG_COND, OSIGMAS, HSUBG_AUCV,         &
                                  PTSTEP, PZZ, PRHODJ, PRHODREF, PEXNREF,              &
                                  PPABST, PTHT, PRT, PSIGS, PSIGQSAT, PMFCONV,         &
                                  PTHM, PRCM, PPABSM,                                  &
                                  PW_ACT,PDTHRAD, PTHS, PRS, PSVT, PSVS, PSRCS, PCLDFR,&
                                  PCIT, OSEDIC, OACTIT, OSEDC, OSEDI,                  &
                                  ORAIN, OWARM, OHHONI, OCONVHG,                       &
                                  PCF_MF,PRC_MF, PRI_MF,                               &
                                  PINPRC,PINPRC3D,PINPRR,PINPRR3D, PEVAP3D,            &
                                  PINPRS,PINPRS3D,PINPRG,PINPRG3D,PINPRH,PINPRH3D,     &
                                  PSOLORG,PMI,                                         &
                                  PSPEEDC, PSPEEDR, PSPEEDS, PSPEEDG, PSPEEDH,         &
                                  PINDEP, PSUPSAT,  PNACT, PNPRO,PSSPRO, PRAINFR,      &
                                  PSEA,PTOWN          )   
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud
CHARACTER(LEN=4),         INTENT(IN)   :: HACTCCN  ! kind of CCN activation scheme
                                                   ! paramerization
CHARACTER(LEN=4),         INTENT(IN)   :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HMF_CLOUD! Type of statistical cloud
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
INTEGER,                  INTENT(IN)   :: KSPLITG  ! Number of small time step
                                       ! integrations for  ice  sedimendation
INTEGER,                  INTENT(IN)   :: KMI      ! Model index
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE   ! Output file
CHARACTER*4,              INTENT(IN)   :: HRAD     ! Radiation scheme name
CHARACTER*4,              INTENT(IN)   :: HTURBDIM ! Dimensionality of the
                                                   ! turbulence scheme
LOGICAL,                  INTENT(IN)   :: OCLOSE_OUT ! Conditional closure of
                                                   ! the OUTPUT FM-file
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)   :: PTSTEP ! Time step :XTSTEP in namelist
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL,                     INTENT(IN)   :: PSIGQSAT! coeff applied to qsat variance contribution
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHM    ! Theta at time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABSM   ! Pressure time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRCM    ! Cloud water m.r. at time t-Dt
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_ACT ! W for CCN activation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD ! THeta RADiative Tendancy
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT  ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS  ! Scalar variable sources
!
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PCLDFR! Cloud fraction
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PCIT  ! Pristine ice number
                                                 ! concentration at time t
LOGICAL,                  INTENT(IN)    :: OSEDIC! Switch to activate the
                                                 ! cloud droplet sedimentation
                                                 ! for ICE3            
LOGICAL,                  INTENT(IN)    :: OACTIT ! Switch to activate the
                                                 ! activation through temp.
                                                 ! evolution in C2R2 and KHKO
LOGICAL,                  INTENT(IN)    :: OSEDC ! Switch to activate the
                                                 ! cloud droplet sedimentation
                                                 ! for C2R2 or KHKO
LOGICAL,                  INTENT(IN)    :: OSEDI ! Switch to activate the
                                                 ! cloud crystal sedimentation
LOGICAL,                  INTENT(IN)    :: ORAIN ! Switch to activate the
                                                 ! raindrop formation
LOGICAL,                  INTENT(IN)    :: OWARM ! Control of the rain formation
                                                 !  by slow warm microphysical
                                                 !         processes
LOGICAL,                  INTENT(IN)    :: OHHONI! enable haze freezing
LOGICAL,                  INTENT(IN)    :: OCONVHG! Switch for conversion from
                                                  ! hail to graupel
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux solid mixing ratio
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRR3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PEVAP3D  ! evap profile
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRH! Hail instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRC3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRS3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRG3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRH3D ! sed flux of precip
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMI !
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDC ! Cloud sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDR ! Rain sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDS ! Snow sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDG ! Graupel sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDH ! Hail sedimentation speed
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINDEP! Cloud instant deposition
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSUPSAT  !sursat
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNACT  !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNPRO  !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSPRO   !sursat
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRAINFR ! Rain fraction                
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA      ! Land Sea mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN      ! Town fraction
!
END SUBROUTINE RESOLVED_CLOUD
END INTERFACE
END MODULE MODI_RESOLVED_CLOUD
!
!     ##########################################################################
      SUBROUTINE RESOLVED_CLOUD ( HCLOUD, HACTCCN, HSCONV, HMF_CLOUD,                  &
                                  KRR, KSPLITR, KSPLITG, KMI, KTCOUNT,                 &
                                  HLBCX, HLBCY, TPFILE, HRAD, HTURBDIM,                &
                                  OCLOSE_OUT, OSUBG_COND, OSIGMAS, HSUBG_AUCV,         &
                                  PTSTEP, PZZ, PRHODJ, PRHODREF, PEXNREF,              &
                                  PPABST, PTHT, PRT, PSIGS, PSIGQSAT, PMFCONV,         &
                                  PTHM, PRCM, PPABSM,                                  &
                                  PW_ACT,PDTHRAD, PTHS, PRS, PSVT, PSVS, PSRCS, PCLDFR,&
                                  PCIT, OSEDIC, OACTIT, OSEDC, OSEDI,                  &
                                  ORAIN, OWARM, OHHONI, OCONVHG,                       &
                                  PCF_MF,PRC_MF, PRI_MF,                               &
                                  PINPRC,PINPRC3D,PINPRR,PINPRR3D, PEVAP3D,            &
                                  PINPRS,PINPRS3D,PINPRG,PINPRG3D,PINPRH,PINPRH3D,     &
                                  PSOLORG,PMI,                                         &
                                  PSPEEDC, PSPEEDR, PSPEEDS, PSPEEDG, PSPEEDH,         &
                                  PINDEP, PSUPSAT,  PNACT, PNPRO,PSSPRO, PRAINFR,      &
                                  PSEA,PTOWN          )   
!     ##########################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation
!!
!!
!!**  METHOD
!!    ------
!!      The main actions of this routine is to call the routines computing the
!!    microphysical sources. Before that:
!!        - it computes the real absolute pressure,
!!        - negative values of the current guess of all mixing ratio are removed.
!!          This is done by a global filling algorithm based on a multiplicative
!!          method (Rood, 1987), in order to conserved the total mass in the
!!          simulation domain.
!!        - Sources are transformed in physical tendencies, by removing the
!!          multiplicative term Rhod*J.
!!        - External points values are filled owing to the use of cyclic
!!          l.b.c., in order to performe computations on the full domain.
!!      After calling to microphysical routines, the physical tendencies are
!!    switched back to prognostic variables.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine SLOW_TERMS: Computes the explicit microphysical sources
!!      Subroutine FAST_TERMS: Performs the saturation adjustment for l
!!      Subroutine RAIN_ICE  : Computes the explicit microphysical sources for i
!!      Subroutine ICE_ADJUST: Performs the saturation adjustment for i+l
!!      MIN_ll,SUM3D_ll : distributed functions equivalent to MIN and SUM
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains declarations of parameter variables
!!         JPHEXT       : Horizontal external points number
!!         JPVEXT       : Vertical external points number
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD                ! Gaz  constant for dry air
!!          XCPD               ! Cpd (dry air)
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 and book2 of documentation ( routine RESOLVED_CLOUD )
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/12/94
!!      Modifications: June 8, 1995 ( J.Stein )
!!                                   Cleaning to improve efficienty and clarity
!!                                  in agreement with the MESO-NH coding norm
!!                     March 1, 1996 ( J.Stein )
!!                                   store the cloud fraction
!!                     March 18, 1996 ( J.Stein )
!!                                   check that ZMASSPOS /= 0
!!                     Oct.  12, 1996 ( J.Stein )
!!                                   remove the negative values correction
!!                                   for the KES2 case
!!      Modifications: Dec 14, 1995 (J.-P. Pinty)
!!                                   Add the mixed-phase option
!!      Modifications: Jul 01, 1996 (J.-P. Pinty)
!!                                   Change arg. list in routine FAST_TERMS
!!      Modifications: Jan 27, 1997 (J.-P. Pinty)
!!                                   add W and SV in arg. list
!!      Modifications: March 23, 98 (E.Richard)
!!                                   correction of negative value based on
!!                                  rv+rc+ri and thetal or thetail conservation
!!      Modifications: April 08, 98 (J.-P. Lafore and V. Ducrocq )
!!                                  modify the  correction of negative values
!!      Modifications: June 08, 00  (J.-P. Pinty and J.-M. Cohard)
!!                                  add the C2R2 scheme
!!      Modifications: April 08, 01  (J.-P. Pinty)
!!                                  add the C3R5 scheme
!!      Modifications: July  21, 01  (J.-P. Pinty)
!!                                  Add OHHONI and PW_ACT (for haze freezing)
!!      Modifications: Sept 21, 01  (J.-P. Pinty)
!!                                  Add XCONC_CCN limitation
!!      Modifications: Nov  21, 02  (J.-P. Pinty)
!!                                  Add ICE4 and C3R5 options
!!                     June, 2005   (V. Masson)
!!                                  Technical change in interface for scalar arguments
!!      Modifications : March, 2006 (O.Geoffroy)
!!                                  Add KHKO scheme
!!      Modifications : March 2013  (O.Thouron)
!!                                  Add prognostic supersaturation
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      M.Mazoyer : 04/2016 : Temperature radiative tendency used for  
!!                            activation by cooling (OACTIT)
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!                     10/2016 M.Mazoyer New KHKO output fields
!!                    10/2016 (C.Lac) Add droplet deposition
!!      S.Riette  : 11/2016 : ice_adjust before and after rain_ice
!!                            ICE3/ICE4 modified, old version under LRED=F   
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      P. Wautelet: 01/02/2019: ZRSMIN is now allocatable (instead of size of XRTMIN which was sometimes not allocated)
!!                   02/2019 C.Lac add rain fraction as an output field
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_BUDGET,         ONLY: LBUDGET_TH, LBUDGET_RC, LBUDGET_RG, LBUDGET_RH, LBUDGET_RI, LBUDGET_RR, LBUDGET_RS, LBUDGET_RV, &
                               LBUDGET_SV
USE MODD_CH_AEROSOL,     ONLY: LORILAM
USE MODD_DUST,           ONLY: LDUST
USE MODD_CST,            ONLY: XCI, XCL, XCPD, XCPV, XLSTT, XLVTT, XMNH_TINY, XP00, XRD, XRHOLW, XTT
USE MODD_DUST ,          ONLY: LDUST
USE MODD_IO_ll,          ONLY: TFILEDATA
USE MODD_NSV,            ONLY: NSV_C1R3END, NSV_C2R2BEG, NSV_C2R2END,                            &
                               NSV_LIMA_BEG, NSV_LIMA_END, NSV_LIMA_CCN_FREE, NSV_LIMA_IFN_FREE, &
                               NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_NR
USE MODD_PARAM_C2R2,     ONLY: LSUPSAT
USE MODD_PARAMETERS,     ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_ICE,      ONLY: CSEDIM, LADJ_BEFORE, LADJ_AFTER, CFRAC_ICE_ADJUST, LRED
USE MODD_PARAM_LIMA,     ONLY: LCOLD, XCONC_CCN_TOT, NMOD_CCN, NMOD_IFN, NMOD_IMM, LPTSPLIT, &
                               YRTMIN=>XRTMIN, YCTMIN=>XCTMIN
USE MODD_RAIN_ICE_DESCR, ONLY: XRTMIN
USE MODD_SALT,           ONLY: LSALT
!
USE MODE_ll
!
USE MODI_BUDGET
USE MODI_C2R2_ADJUST
USE MODI_C3R5_ADJUST
USE MODI_FAST_TERMS
USE MODI_GET_HALO
USE MODI_ICE_ADJUST
USE MODI_ICE_C1R3
USE MODI_KHKO_NOTADJUST
USE MODI_LIMA
USE MODI_LIMA_WARM
USE MODI_LIMA_COLD
USE MODI_LIMA_MIXED
USE MODI_LIMA_ADJUST
USE MODI_RAIN_C2R2_KHKO
USE MODI_RAIN_ICE
USE MODI_RAIN_ICE_RED
USE MODI_SHUMAN
USE MODI_SLOW_TERMS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud
                                                   ! paramerization
CHARACTER(LEN=4),         INTENT(IN)   :: HACTCCN  ! kind of CCN activation
CHARACTER(LEN=4),         INTENT(IN)   :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HMF_CLOUD! Type of statistical cloud
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
INTEGER,                  INTENT(IN)   :: KSPLITG  ! Number of small time step
                                       ! integrations for  ice  sedimendation
INTEGER,                  INTENT(IN)   :: KMI      ! Model index
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE   ! Output file
CHARACTER*4,              INTENT(IN)   :: HRAD     ! Radiation scheme name
CHARACTER*4,              INTENT(IN)   :: HTURBDIM ! Dimensionality of the
                                                   ! turbulence scheme
LOGICAL,                  INTENT(IN)   :: OCLOSE_OUT ! Conditional closure of
                                                   ! the OUTPUT FM-file
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)   :: PTSTEP ! Time step :XTSTEP in namelist
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL,                     INTENT(IN)   :: PSIGQSAT! coeff applied to qsat variance contribution
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHM    ! Theta at time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABSM   ! Pressure time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRCM    ! Cloud water m.r. at time t-Dt
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_ACT ! W for CCN activation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD ! THeta RADiative Tendancy
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT  ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS  ! Scalar variable sources
!
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PCLDFR! Cloud fraction
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PCIT  ! Pristine ice number
                                                 ! concentration at time t
LOGICAL,                  INTENT(IN)    :: OSEDIC! Switch to activate the
                                                 ! cloud droplet sedimentation
                                                 ! for ICE3            
LOGICAL,                  INTENT(IN)    :: OACTIT ! Switch to activate the
                                                 ! activation through temp.
                                                 ! evolution in C2R2 and KHKO
LOGICAL,                  INTENT(IN)    :: OSEDC ! Switch to activate the
                                                 ! cloud droplet sedimentation
LOGICAL,                  INTENT(IN)    :: OSEDI ! Switch to activate the
                                                 ! cloud crystal sedimentation
LOGICAL,                  INTENT(IN)    :: ORAIN ! Switch to activate the
                                                 ! raindrop formation
LOGICAL,                  INTENT(IN)    :: OWARM ! Control of the rain formation
                                                 !  by slow warm microphysical
                                                 !         processes
LOGICAL,                  INTENT(IN)    :: OHHONI! enable haze freezing
LOGICAL,                  INTENT(IN)    :: OCONVHG! Switch for conversion from
                                                  ! hail to graupel
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux solid mixing ratio
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRR3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PEVAP3D  ! evap profile
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRH! Hail instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRC3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRS3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRG3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRH3D ! sed flux of precip
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMI !
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDC ! Cloud sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDR ! Rain sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDS ! Snow sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDG ! Graupel sedimentation speed
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSPEEDH ! Hail sedimentation speed
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP! Cloud instant deposition
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSUPSAT  !sursat
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNACT  !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNPRO  !concentrtaion d'aérosols activés au temps t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSPRO   !sursat
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRAINFR ! Rain fraction                
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA      ! Land Sea mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN      ! Town fraction
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR,JSV       ! Loop index for the moist and scalar variables
INTEGER :: IIB           !  Define the physical domain
INTEGER :: IIE           !
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !
INTEGER :: IKU
INTEGER :: IINFO_ll      ! return code of parallel routine
INTEGER :: JK,JI,JL
!
!
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZDZZ
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZT,ZEXN,ZLV,ZLS,ZCPH
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZCOR
                                    ! for the correction of negative rv
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZZZ
                                    ! model layer height
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR
!
INTEGER                               :: ISVBEG ! first scalar index for microphysics
INTEGER                               :: ISVEND ! last  scalar index for microphysics
REAL, DIMENSION(:),       ALLOCATABLE :: ZRSMIN ! Minimum value for tendencies
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSVT   ! scalar variable for microphysics only
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSVS   ! scalar tendency for microphysics only
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: LLMICRO ! mask to limit computation
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3), KRR) :: ZFPR
!
INTEGER                               :: JMOD, JMOD_IFN
! BVIE work array waiting for PINPRI
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)):: ZINPRI
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
IKU=SIZE(PZZ,3)
!
IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO') THEN
  ISVBEG = NSV_C2R2BEG
  ISVEND = NSV_C2R2END
ELSE IF (HCLOUD == 'C3R5') THEN
  ISVBEG = NSV_C2R2BEG
  ISVEND = NSV_C1R3END
ELSE IF (HCLOUD == 'LIMA') THEN
  ISVBEG = NSV_LIMA_BEG
  ISVEND = NSV_LIMA_END
ELSE
  ISVBEG = 0
  ISVEND = 0
END IF
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
  ALLOCATE(ZSVT(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),ISVEND - ISVBEG + 1))
  ALLOCATE(ZSVS(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),ISVEND - ISVBEG + 1))
  ZSVT(:,:,:,:) = PSVT(:,:,:,ISVBEG:ISVEND)
  ZSVS(:,:,:,:) = PSVS(:,:,:,ISVBEG:ISVEND)
END IF
IF (HCLOUD(1:3)=='ICE') THEN
  ALLOCATE(ZRSMIN(SIZE(XRTMIN)))
  ZRSMIN(:) = XRTMIN(:) / PTSTEP
END IF
!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
PTHS(:,:,:) = PTHS(:,:,:) / PRHODJ(:,:,:)
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) / PRHODJ(:,:,:)
END DO
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
  DO JSV = 1,SIZE(ZSVS,4)
    ZSVS(:,:,:,JSV) = ZSVS(:,:,:,JSV) / PRHODJ(:,:,:)
  ENDDO
ENDIF
!
!  complete the lateral boundaries to avoid possible problems
!
DO JI=1,JPHEXT
 PTHS(JI,:,:) = PTHS(IIB,:,:)
 PTHS(IIE+JI,:,:) = PTHS(IIE,:,:)
 PTHS(:,JI,:) = PTHS(:,IJB,:)
 PTHS(:,IJE+JI,:) = PTHS(:,IJE,:)
!
 PRS(JI,:,:,:) = PRS(IIB,:,:,:)
 PRS(IIE+JI,:,:,:) = PRS(IIE,:,:,:)
 PRS(:,JI,:,:) = PRS(:,IJB,:,:)
 PRS(:,IJE+JI,:,:) = PRS(:,IJE,:,:)
END DO
!
!  complete the physical boundaries to avoid some computations
!
IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL')  PRT(:IIB-1,:,:,2:) = 0.0
IF(LEAST_ll()  .AND. HLBCX(2) /= 'CYCL')  PRT(IIE+1:,:,:,2:) = 0.0
IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL')  PRT(:,:IJB-1,:,2:) = 0.0
IF(LNORTH_ll() .AND. HLBCY(2) /= 'CYCL')  PRT(:,IJE+1:,:,2:) = 0.0
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
DO JI=1,JPHEXT
  ZSVS(JI,:,:,:) = ZSVS(IIB,:,:,:)
  ZSVS(IIE+JI,:,:,:) = ZSVS(IIE,:,:,:)
  ZSVS(:,JI,:,:) = ZSVS(:,IJB,:,:)
  ZSVS(:,IJE+JI,:,:) = ZSVS(:,IJE,:,:)
END DO
 !
!  complete the physical boundaries to avoid some computations
!
  IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL')  ZSVT(:IIB-1,:,:,:) = 0.0
  IF(LEAST_ll()  .AND. HLBCX(2) /= 'CYCL')  ZSVT(IIE+1:,:,:,:) = 0.0
  IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL')  ZSVT(:,:IJB-1,:,:) = 0.0
  IF(LNORTH_ll() .AND. HLBCY(2) /= 'CYCL')  ZSVT(:,IJE+1:,:,:) = 0.0
ENDIF
!
!  complete the vertical boundaries
!
PTHS(:,:,IKB-1) = PTHS(:,:,IKB)
PTHS(:,:,IKE+1) = PTHS(:,:,IKE)
!
PRS(:,:,IKB-1,:) = PRS(:,:,IKB,:)
PRS(:,:,IKE+1,:) = PRS(:,:,IKE,:)
!
PRT(:,:,IKB-1,:) = PRT(:,:,IKB,:)
PRT(:,:,IKE+1,:) = PRT(:,:,IKE,:)
!
IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'C3R5' .OR. HCLOUD == 'KHKO' &
                                           .OR. HCLOUD == 'LIMA') THEN
  ZSVS(:,:,IKB-1,:) = ZSVS(:,:,IKB,:)
  ZSVS(:,:,IKE+1,:) = ZSVS(:,:,IKE,:)
  ZSVT(:,:,IKB-1,:) = ZSVT(:,:,IKB,:)
  ZSVT(:,:,IKE+1,:) = ZSVT(:,:,IKE,:)
ENDIF
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
ZEXN(:,:,:)= (PPABST(:,:,:)/XP00)**(XRD/XCPD)
ZT(:,:,:)= PTHT(:,:,:)*ZEXN(:,:,:)
ZLV(:,:,:)=XLVTT +(XCPV-XCL) *(ZT(:,:,:)-XTT)
ZLS(:,:,:)=XLSTT +(XCPV-XCI) *(ZT(:,:,:)-XTT)
ZCPH(:,:,:)=XCPD +XCPV*PRT(:,:,:,1)
!
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
IF (HCLOUD == 'KESS' .OR. HCLOUD == 'ICE3'                                &
    .OR. HCLOUD == 'C2R2' .OR. &
    HCLOUD == 'C3R5' .OR. HCLOUD == 'KHKO' .OR. HCLOUD=='LIMA' ) THEN
!
  DO JRR = 3,KRR
    SELECT CASE (JRR)
      CASE(3,5,6,7) ! rain, snow, graupel and hail

        IF ( MIN_ll( PRS(:,:,:,JRR), IINFO_ll) < 0.0 ) THEN
!
! compute the total water mass computation
!
          ZMASSTOT = MAX( 0. , SUM3D_ll( PRS(:,:,:,JRR), IINFO_ll ) )
!
! remove the negative values
!
          PRS(:,:,:,JRR) = MAX( 0., PRS(:,:,:,JRR) )
!
! compute the new total mass
!
          ZMASSPOS = MAX(XMNH_TINY,SUM3D_ll( PRS(:,:,:,JRR), IINFO_ll ) )
!
! correct again in such a way to conserve the total mass
!
          ZRATIO = ZMASSTOT / ZMASSPOS
          PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * ZRATIO
!
        END IF
    END SELECT
  END DO
END IF
!
!*       3.2    Adjustement for liquid and solid cloud
!
SELECT CASE ( HCLOUD )
  CASE('KESS')
    WHERE (PRS(:,:,:,2) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2) = 0.0
    END WHERE
!
!
! CASE('C2R2','KHKO')                                 
!   CALL GET_HALO(PRS(:,:,:,2))
!   CALL GET_HALO(ZSVS(:,:,:,2))
!   WHERE (PRS(:,:,:,2) < 0. .OR. ZSVS(:,:,:,2) < 0.)
!     ZSVS(:,:,:,1) = 0.0
!   END WHERE
!   DO JSV = 2, 3
!     WHERE (PRS(:,:,:,JSV) < 0. .OR. ZSVS(:,:,:,JSV) < 0.)
!       PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,JSV)
!       PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,JSV) * ZLV(:,:,:) /  &
!            ZCPH(:,:,:) / ZEXN(:,:,:)
!       PRS(:,:,:,JSV)  = 0.0
!       ZSVS(:,:,:,JSV) = 0.0
!     END WHERE
!   ENDDO
! Commented 03/2013 O.Thouron 
! (at least necessary to be commented for supersaturation variable)
!  ZSVS(:,:,:,:) = MAX( 0.0,ZSVS(:,:,:,:) )
!
!
  CASE('ICE3','ICE4')
    WHERE (PRS(:,:,:,4) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLS(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,4) = 0.
    END WHERE
!
!   cloud
    WHERE (PRS(:,:,:,2) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2) = 0.
    END WHERE
!
! if rc or ri are positive, we can correct negative rv
!   cloud
    WHERE ((PRS(:,:,:,1) <0.) .AND. (PRS(:,:,:,2)> 0.) )
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2) = 0.
    END WHERE
!   ice
    IF(KRR > 3) THEN
      WHERE ((PRS(:,:,:,1) < 0.).AND.(PRS(:,:,:,4) > 0.))
        ZCOR(:,:,:)=MIN(-PRS(:,:,:,1),PRS(:,:,:,4))
        PRS(:,:,:,1) = PRS(:,:,:,1) + ZCOR(:,:,:)
        PTHS(:,:,:) = PTHS(:,:,:) - ZCOR(:,:,:) * ZLS(:,:,:) /  &
             ZCPH(:,:,:) / PEXNREF(:,:,:)
        PRS(:,:,:,4) = PRS(:,:,:,4) -ZCOR(:,:,:)
      END WHERE
    END IF
!
   CASE('C3R5')
    WHERE (PRS(:,:,:,2) < 0. .OR. ZSVS(:,:,:,2) < 0.)
      ZSVS(:,:,:,1) = 0.0
    END WHERE
    DO JSV = 2, 3
      WHERE (PRS(:,:,:,JSV) < 0. .OR. ZSVS(:,:,:,JSV) < 0.)
        PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,JSV)
        PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,JSV) * ZLV(:,:,:) /  &
             ZCPH(:,:,:) / PEXNREF(:,:,:)
        PRS(:,:,:,JSV)  = 0.0
        ZSVS(:,:,:,JSV) = 0.0
      END WHERE
    ENDDO
    ZSVS(:,:,:,:) = MAX( 0.0,ZSVS(:,:,:,:) )
!   ice
    WHERE (PRS(:,:,:,4) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,4)  = 0.0
      PSVS(:,:,:,4) = 0.0
    END WHERE
!   cloud
    WHERE (PRS(:,:,:,2) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2)  = 0.0
      PSVS(:,:,:,2) = 0.0
    END WHERE
    PSVS(:,:,:,:) = MAX( 0.0,PSVS(:,:,:,:) )
!
   CASE('LIMA')   
! Correction where rc<0 or Nc<0
      IF (OWARM) THEN
         WHERE (PRS(:,:,:,2) < 0. .OR. ZSVS(:,:,:,NSV_LIMA_NC) < 0.)
            PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
            PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
                 ZCPH(:,:,:) / ZEXN(:,:,:)
            PRS(:,:,:,2)  = 0.0
            ZSVS(:,:,:,NSV_LIMA_NC) = 0.0
         END WHERE
      END IF
! Correction where rr<0 or Nr<0
      IF (OWARM .AND. ORAIN) THEN
         WHERE (PRS(:,:,:,3) < YRTMIN(3)/PTSTEP .OR. ZSVS(:,:,:,NSV_LIMA_NR) < YCTMIN(3)/PTSTEP)
            PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,3)
            PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,3) * ZLV(:,:,:) /  &
                 ZCPH(:,:,:) / ZEXN(:,:,:)
            PRS(:,:,:,3)  = 0.0
            ZSVS(:,:,:,NSV_LIMA_NR) = 0.0
         END WHERE
      END IF
! Correction where ri<0 or Ni<0
      IF (LCOLD) THEN
         WHERE (PRS(:,:,:,4) < YRTMIN(4)/PTSTEP .OR. ZSVS(:,:,:,NSV_LIMA_NI) < YCTMIN(4)/PTSTEP)
            PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
            PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLS(:,:,:) /  &
                 ZCPH(:,:,:) / ZEXN(:,:,:)
            PRS(:,:,:,4)  = 0.0
            ZSVS(:,:,:,NSV_LIMA_NI) = 0.0
         END WHERE
      END IF
!
     ZSVS(:,:,:,:) = MAX( 0.0,ZSVS(:,:,:,:) )
     PRS(:,:,:,:) = MAX( 0.0,PRS(:,:,:,:) )
!
END SELECT
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------
!
IF ((HCLOUD /= 'KHKO') .AND. (HCLOUD /= 'C2R2') ) THEN
 IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:)  * PRHODJ(:,:,:), 4,'NEGA_BU_RTH')
 IF (LBUDGET_RV) CALL BUDGET (PRS(:,:,:,1) * PRHODJ(:,:,:), 6,'NEGA_BU_RRV')
 IF (LBUDGET_RC) CALL BUDGET (PRS(:,:,:,2) * PRHODJ(:,:,:), 7,'NEGA_BU_RRC')
END IF
IF (LBUDGET_RR) CALL BUDGET (PRS(:,:,:,3) * PRHODJ(:,:,:), 8,'NEGA_BU_RRR')
IF (LBUDGET_RI) CALL BUDGET (PRS(:,:,:,4) * PRHODJ(:,:,:) ,9,'NEGA_BU_RRI')
IF (LBUDGET_RS) CALL BUDGET (PRS(:,:,:,5) * PRHODJ(:,:,:),10,'NEGA_BU_RRS')
IF (LBUDGET_RG) CALL BUDGET (PRS(:,:,:,6) * PRHODJ(:,:,:),11,'NEGA_BU_RRG')
IF (LBUDGET_RH) CALL BUDGET (PRS(:,:,:,7) * PRHODJ(:,:,:),12,'NEGA_BU_RRH')
IF (LBUDGET_SV .AND. (HCLOUD == 'LIMA')) THEN
   IF (OWARM) CALL BUDGET (ZSVS(:,:,:,NSV_LIMA_NC) * PRHODJ(:,:,:),12+NSV_LIMA_NC,'NEGA_BU_RSV')
   IF (OWARM.AND.ORAIN) CALL BUDGET (ZSVS(:,:,:,NSV_LIMA_NR) * PRHODJ(:,:,:),12+NSV_LIMA_NR,'NEGA_BU_RSV')
   IF (LCOLD) CALL BUDGET (ZSVS(:,:,:,NSV_LIMA_NI) * PRHODJ(:,:,:),12+NSV_LIMA_NI,'NEGA_BU_RSV')
   IF (NMOD_CCN.GE.1) THEN
      DO JL=1, NMOD_CCN
         CALL BUDGET ( ZSVS(:,:,:,NSV_LIMA_CCN_FREE+JL-1)* &
              PRHODJ(:,:,:),12+NSV_LIMA_CCN_FREE+JL-1,'NEGA_BU_RSV') 
      END DO
   END IF
   IF (NMOD_IFN.GE.1) THEN
      DO JL=1, NMOD_IFN
         CALL BUDGET ( ZSVS(:,:,:,NSV_LIMA_IFN_FREE+JL-1)* &
              PRHODJ(:,:,:),12+NSV_LIMA_IFN_FREE+JL-1,'NEGA_BU_RSV') 
      END DO
   END IF
END IF
!

!*       3.4    Limitations of Na and Nc to the CCN max number concentration
!
! Commented by O.Thouron 03/2013
!IF ((HCLOUD == 'C2R2' .OR. HCLOUD == 'C3R5' .OR. HCLOUD == 'KHKO') &
!     .AND.(XCONC_CCN > 0)) THEN
!  IF ((HACTCCN /= 'ABRK')) THEN
!  ZSVT(:,:,:,1) = MIN( ZSVT(:,:,:,1),XCONC_CCN )
!  ZSVT(:,:,:,2) = MIN( ZSVT(:,:,:,2),XCONC_CCN )
!  ZSVS(:,:,:,1) = MIN( ZSVS(:,:,:,1),XCONC_CCN )
!  ZSVS(:,:,:,2) = MIN( ZSVS(:,:,:,2),XCONC_CCN )
!  END IF
!END IF
!
!
!-------------------------------------------------------------------------------
!
SELECT CASE ( HCLOUD )
  CASE ('REVE')
!
!*       4.     REVERSIBLE MICROPHYSICAL SCHEME
!               -------------------------------
!
    CALL FAST_TERMS ( KRR, KMI, HRAD, HTURBDIM,                                &
                      HSCONV, HMF_CLOUD, OSUBG_COND, PTSTEP,                   &
                      PRHODJ, PSIGS, PPABST,                                   &
                      PCF_MF,PRC_MF,                                           &
                      PRVT=PRT(:,:,:,1), PRCT=PRT(:,:,:,2),                    &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                      PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR                    )
!
  CASE ('KESS')
!
!*       5.     KESSLER MICROPHYSICAL SCHEME
!               ----------------------------
!
!
!*       5.1    Compute the explicit microphysical sources
!
    CALL SLOW_TERMS ( KSPLITR, PTSTEP, KMI, HSUBG_AUCV,                       &
                      PZZ, PRHODJ, PRHODREF, PCLDFR,                          &
                      PTHT, PRT(:,:,:,1), PRT(:,:,:,2), PRT(:,:,:,3), PPABST, &
                      PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),         &
                      PINPRR, PINPRR3D, PEVAP3D                         )
!
!*       5.2    Perform the saturation adjustment
!
    CALL FAST_TERMS ( KRR, KMI, HRAD, HTURBDIM,                                &
                      HSCONV, HMF_CLOUD, OSUBG_COND, PTSTEP,                   &
                      PRHODJ, PSIGS, PPABST,                                   &
                      PCF_MF,PRC_MF,                                           &
                      PRVT=PRT(:,:,:,1), PRCT=PRT(:,:,:,2),                    &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2), PRRS=PRS(:,:,:,3), &
                      PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR                    )
!
!
  CASE ('C2R2','KHKO')
!
!*       7.     2-MOMENT WARM MICROPHYSICAL SCHEME C2R2 or KHKO
!               ---------------------------------------
!
!
!*       7.1    Compute the explicit microphysical sources
!
!
    CALL RAIN_C2R2_KHKO ( HCLOUD, OACTIT, OSEDC, ORAIN, KSPLITR, PTSTEP, KMI,     &
                     TPFILE, OCLOSE_OUT, PZZ, PRHODJ, PRHODREF, PEXNREF,          &
                     PPABST, PTHT, PRT(:,:,:,1), PRT(:,:,:,2),  PRT(:,:,:,3),     &
                     PTHM, PRCM, PPABSM,                                          &
                     PW_ACT,PDTHRAD,PTHS, PRS(:,:,:,1),PRS(:,:,:,2),PRS(:,:,:,3), &
                     ZSVT(:,:,:,1), ZSVT(:,:,:,2), ZSVT(:,:,:,3),                 &
                     ZSVS(:,:,:,1), ZSVS(:,:,:,2), ZSVS(:,:,:,3),                 &
                     PINPRC, PINPRR, PINPRR3D, PEVAP3D ,                          &
                     PSVT(:,:,:,:), PSOLORG, PMI, HACTCCN,                        &
                     PINDEP, PSUPSAT, PNACT                                       )
!
!
!*       7.2    Perform the saturation adjustment
!
   IF (LSUPSAT) THEN
    CALL KHKO_NOTADJUST (KRR, KTCOUNT,TPFILE, HRAD, OCLOSE_OUT,                  &
                         PTSTEP, PRHODJ, PPABSM, PPABST, PRHODREF, PZZ,          &
                         PTHT,PRT(:,:,:,1),PRT(:,:,:,2),PRT(:,:,:,3),            &
                         PTHS,PRS(:,:,:,1),PRS(:,:,:,2),PRS(:,:,:,3),            &
                         ZSVS(:,:,:,2),ZSVS(:,:,:,1),                            &
                         ZSVS(:,:,:,4), PCLDFR, PSRCS , PNPRO,PSSPRO             )
!
   ELSE
    CALL C2R2_ADJUST ( KRR,TPFILE, HRAD,                              &
                       HTURBDIM, OCLOSE_OUT, OSUBG_COND, PTSTEP,               &
                       PRHODJ, PSIGS, PPABST,                                  &
                       PTHS=PTHS, PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),        &
                       PCNUCS=ZSVS(:,:,:,1), PCCS=ZSVS(:,:,:,2),               &
                       PSRCS=PSRCS, PCLDFR=PCLDFR, PRRS=PRS(:,:,:,3)           )
!
   END IF
!
  CASE ('ICE3')
!
!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!
!*       9.1    Compute the explicit microphysical sources
!
!
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZZZ = MZF(1,IKU,1, PZZ )
    IF(LRED .AND. LADJ_BEFORE) THEN
      CALL ICE_ADJUST (1,IKU,1, KRR, CFRAC_ICE_ADJUST, 'ADJU',                 &
                      OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,                    &
                      PRHODJ, PEXNREF,  PSIGS, PMFCONV, PPABST, ZZZ,           &
                      ZEXN, PCF_MF,PRC_MF,PRI_MF,                              &   
                      PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,        &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                      PTH=PTHS*PTSTEP, PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR,  &
                      PRR=PRS(:,:,:,3)*PTSTEP,                                 &
                      PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),              &
                      PRS=PRS(:,:,:,5)*PTSTEP,                                 &
                      PRG=PRS(:,:,:,6)*PTSTEP                                  )
    ENDIF
    IF (LRED) THEN
      LLMICRO(:,:,:)=PRT(:,:,:,2)>XRTMIN(2) .OR. &
                   PRT(:,:,:,3)>XRTMIN(3) .OR. &
                   PRT(:,:,:,4)>XRTMIN(4) .OR. &
                   PRT(:,:,:,5)>XRTMIN(5) .OR. &
                   PRT(:,:,:,6)>XRTMIN(6)
      LLMICRO(:,:,:)=LLMICRO(:,:,:) .OR. &
                   PRS(:,:,:,2)>ZRSMIN(2) .OR. &
                   PRS(:,:,:,3)>ZRSMIN(3) .OR. &
                   PRS(:,:,:,4)>ZRSMIN(4) .OR. &
                   PRS(:,:,:,5)>ZRSMIN(5) .OR. &
                   PRS(:,:,:,6)>ZRSMIN(6)
      CALL RAIN_ICE_RED ( OSEDIC,CSEDIM, HSUBG_AUCV, OWARM,1,IKU,1,      &
                    PTSTEP, KRR, LLMICRO, ZEXN,                          &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT,PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                    &
                    PRT(:,:,:,3), PRT(:,:,:,4),                          &
                    PRT(:,:,:,5), PRT(:,:,:,6),                          &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),      &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),            &
                    PINPRC,PINPRR, PINPRR3D, PEVAP3D,                    &
                    PINPRS, PINPRG, PSIGS, PINDEP, PRAINFR, PSEA,PTOWN, PFPR=ZFPR)
    ELSE 
      CALL RAIN_ICE ( OSEDIC,CSEDIM, HSUBG_AUCV, OWARM,1,IKU,1,          &
                    KSPLITR, PTSTEP, KRR,                                &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT,PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                    &
                    PRT(:,:,:,3), PRT(:,:,:,4),                          &
                    PRT(:,:,:,5), PRT(:,:,:,6),                          &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),      &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),            &
                    PINPRC,PINPRR, PINPRR3D, PEVAP3D,                    &
                    PINPRS, PINPRG, PSIGS,PINDEP, PRAINFR,               &
                    PSEA,PTOWN, PFPR=ZFPR)
    END IF
!
!*       9.2    Perform the saturation adjustment over cloud ice and cloud water
!
!
    IF (.NOT. LRED .OR. (LRED .AND. LADJ_AFTER) ) THEN
      CALL ICE_ADJUST (1,IKU,1, KRR, CFRAC_ICE_ADJUST, 'DEPI',    &
                    OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,                    &
                    PRHODJ, PEXNREF,  PSIGS, PMFCONV, PPABST, ZZZ,           &
                    ZEXN, PCF_MF,PRC_MF,PRI_MF,                                    &   
                    PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,                    &
                    PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                    PTH=PTHS*PTSTEP, PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR,                   &
                    PRR=PRS(:,:,:,3)*PTSTEP,                                 &
                    PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),              &
                    PRS=PRS(:,:,:,5)*PTSTEP,                                 &
                    PRG=PRS(:,:,:,6)*PTSTEP                                  )
    END IF
!
  CASE ('ICE4')
!
!*       10.    MIXED-PHASE MICROPHYSICAL SCHEME (WITH 4 ICE SPECIES)
!               -----------------------------------------------------
!
!
!*       10.1   Compute the explicit microphysical sources
!
!
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZZZ = MZF(1,IKU,1, PZZ )
    IF(LRED .AND. LADJ_BEFORE) THEN
            CALL ICE_ADJUST (1,IKU,1, KRR, CFRAC_ICE_ADJUST, 'ADJU',                 &
                      OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,                    &
                      PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, ZZZ,            &
                      ZEXN, PCF_MF,PRC_MF,PRI_MF,                              & 
                      PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,        &
                      PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                      PTH=PTHS*PTSTEP, PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR,  &
                      PRR=PRS(:,:,:,3)*PTSTEP,                                 &
                      PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),              &
                      PRS=PRS(:,:,:,5)*PTSTEP,                                 &
                      PRG=PRS(:,:,:,6)*PTSTEP,                                 &
                      PRH=PRS(:,:,:,7)*PTSTEP                                  )
    ENDIF
    IF  (LRED) THEN
      LLMICRO(:,:,:)=PRT(:,:,:,2)>XRTMIN(2) .OR. &
                   PRT(:,:,:,3)>XRTMIN(3) .OR. &
                   PRT(:,:,:,4)>XRTMIN(4) .OR. &
                   PRT(:,:,:,5)>XRTMIN(5) .OR. &
                   PRT(:,:,:,6)>XRTMIN(6) .OR. &
                   PRT(:,:,:,7)>XRTMIN(7)
      LLMICRO(:,:,:)=LLMICRO(:,:,:) .OR. &
                   PRS(:,:,:,2)>ZRSMIN(2) .OR. &
                   PRS(:,:,:,3)>ZRSMIN(3) .OR. &
                   PRS(:,:,:,4)>ZRSMIN(4) .OR. &
                   PRS(:,:,:,5)>ZRSMIN(5) .OR. &
                   PRS(:,:,:,6)>ZRSMIN(6) .OR. &
                   PRS(:,:,:,7)>ZRSMIN(7)
      CALL RAIN_ICE_RED ( OSEDIC,CSEDIM, HSUBG_AUCV, OWARM,1,IKU,1,             &
                    PTSTEP, KRR, LLMICRO, ZEXN,             &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                     &
                    PRT(:,:,:,3), PRT(:,:,:,4),                           &
                    PRT(:,:,:,5), PRT(:,:,:,6),                           &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                    PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                    PINPRS, PINPRG, PSIGS, PINDEP, PRAINFR, PSEA, PTOWN,  &
                    PRT(:,:,:,7), PRS(:,:,:,7), PINPRH, PFPR=ZFPR         )
    ELSE
      CALL RAIN_ICE ( OSEDIC,CSEDIM, HSUBG_AUCV, OWARM,1,IKU,1,         &
                    KSPLITR, PTSTEP, KRR,                                 &
                    ZDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                    PTHT, PRT(:,:,:,1), PRT(:,:,:,2),                     &
                    PRT(:,:,:,3), PRT(:,:,:,4),                           &
                    PRT(:,:,:,5), PRT(:,:,:,6),                           &
                    PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                    PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                    PINPRC, PINPRR, PINPRR3D, PEVAP3D,           &
                    PINPRS, PINPRG, PSIGS,PINDEP, PRAINFR,                &
                    PSEA, PTOWN,                                          &
                    PRT(:,:,:,7),  PRS(:,:,:,7), PINPRH,PFPR=ZFPR )
     END IF


!
!*       10.2   Perform the saturation adjustment over cloud ice and cloud water
!
    IF (.NOT. LRED .OR. (LRED .AND. LADJ_AFTER) ) THEN
     CALL ICE_ADJUST (1,IKU,1, KRR, CFRAC_ICE_ADJUST, 'DEPI',                 &
                    OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,                    &
                    PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, ZZZ,            &
                    ZEXN, PCF_MF,PRC_MF,PRI_MF,                              &                     
                    PRV=PRS(:,:,:,1)*PTSTEP, PRC=PRS(:,:,:,2)*PTSTEP,        &
                    PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),                    &
                    PTH=PTHS*PTSTEP, PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR,  &
                    PRR=PRS(:,:,:,3)*PTSTEP,                                 &
                    PRI=PRS(:,:,:,4)*PTSTEP, PRIS=PRS(:,:,:,4),              &
                    PRS=PRS(:,:,:,5)*PTSTEP,                                 &
                    PRG=PRS(:,:,:,6)*PTSTEP,                                 &
                    PRH=PRS(:,:,:,7)*PTSTEP                                  )
    END IF
!           
!
!*       12.    2-MOMENT MIXED-PHASE MICROPHYSICAL SCHEME LIMA
!               --------------------------------------------------------------
!
!
!*       12.1   Compute the explicit microphysical sources
!
  CASE ('LIMA')
     !
    DO JK=IKB,IKE
      ZDZZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)    
    ENDDO
    ZZZ = MZF(1,IKU,1, PZZ )
     IF (LPTSPLIT) THEN
        CALL LIMA (1, IKU, 1,                                              &
                   PTSTEP, TPFILE, OCLOSE_OUT,                             &
                   PRHODREF, PEXNREF, ZDZZ,                                &
                   PRHODJ, PPABSM, PPABST,                                 &
                   NMOD_CCN, NMOD_IFN, NMOD_IMM,                           &
                   PTHM, PTHT, PRT, ZSVT, PW_ACT,                          &
                   PTHS, PRS, ZSVS,                                        &
                   PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, PINPRH, &
                   PEVAP3D                                         )
     ELSE

        IF (OWARM) CALL LIMA_WARM(OACTIT, OSEDC, ORAIN, KSPLITR, PTSTEP, KMI,   &
                                  TPFILE, OCLOSE_OUT, KRR, PZZ, PRHODJ,         &
                                  PRHODREF, PEXNREF, PW_ACT, PPABSM, PPABST,    &
                                  PTHM, PRCM,                                   &
                                  PTHT, PRT, ZSVT,                              &
                                  PTHS, PRS, ZSVS,                              &
                                  PINPRC, PINPRR, PINDEP, PINPRR3D, PEVAP3D     )
!
        IF (LCOLD) CALL LIMA_COLD(OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,               &
                                  KRR, PZZ, PRHODJ,                                  &
                                  PRHODREF, PEXNREF, PPABST, PW_ACT,                 &
                                  PTHM, PPABSM,                                      &
                                  PTHT, PRT, ZSVT,                                   &
                                  PTHS, PRS, ZSVS,                                   &
                                  PINPRS, PINPRG, PINPRH)
!
        IF (OWARM .AND. LCOLD) CALL LIMA_MIXED(OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,    &
                                               KRR, PZZ, PRHODJ,                       &
                                               PRHODREF, PEXNREF, PPABST, PW_ACT,      &
                                               PTHM, PPABSM,                           &
                                               PTHT, PRT, ZSVT,                        &
                                               PTHS, PRS, ZSVS                         )
     ENDIF
!
!*       12.2   Perform the saturation adjustment
!
     CALL LIMA_ADJUST(KRR, KMI, TPFILE, HRAD,                           &
                      HTURBDIM, OCLOSE_OUT, OSUBG_COND, PTSTEP,         &
                      PRHODREF, PRHODJ, PEXNREF, PPABST, PSIGS, PPABST, &
                      PRT, PRS, ZSVT, ZSVS,                             &
                      PTHS, PSRCS, PCLDFR                               )
!
END SELECT
!
IF(HCLOUD=='ICE3' .OR. HCLOUD=='ICE4' ) THEN
  PINPRC3D=ZFPR(:,:,:,2) / XRHOLW
  PINPRR3D=ZFPR(:,:,:,3) / XRHOLW
  PINPRS3D=ZFPR(:,:,:,5) / XRHOLW
  PINPRG3D=ZFPR(:,:,:,6) / XRHOLW
  IF(KRR==7) PINPRH3D=ZFPR(:,:,:,7) / XRHOLW
  WHERE (PRT(:,:,:,2) > 1.E-04 )
    PSPEEDC=ZFPR(:,:,:,2) / (PRT(:,:,:,2) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,3) > 1.E-04 )
    PSPEEDR=ZFPR(:,:,:,3) / (PRT(:,:,:,3) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,5) > 1.E-04 )
    PSPEEDS=ZFPR(:,:,:,5) / (PRT(:,:,:,5) * PRHODREF(:,:,:))
  ENDWHERE
  WHERE (PRT(:,:,:,6) > 1.E-04 )
    PSPEEDG=ZFPR(:,:,:,6) / (PRT(:,:,:,6) * PRHODREF(:,:,:))
  ENDWHERE
  IF(KRR==7) THEN
    WHERE (PRT(:,:,:,7) > 1.E-04 )
      PSPEEDH=ZFPR(:,:,:,7) / (PRT(:,:,:,7) * PRHODREF(:,:,:))
    ENDWHERE
  ENDIF
ENDIF
!
IF ( (HCLOUD == 'KHKO') .OR. (HCLOUD == 'C2R2') ) THEN
!    CALL GET_HALO(PRS(:,:,:,2))
!    CALL GET_HALO(ZSVS(:,:,:,2))
!    CALL GET_HALO(ZSVS(:,:,:,3))
    WHERE (PRS(:,:,:,2) < 0. .OR. ZSVS(:,:,:,2) < 0.)
      ZSVS(:,:,:,1) = 0.0
    END WHERE
    DO JSV = 2, 3
      WHERE (PRS(:,:,:,JSV) < 0. .OR. ZSVS(:,:,:,JSV) < 0.)
        PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,JSV)
        PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,JSV) * ZLV(:,:,:) /  &
             ZCPH(:,:,:) / PEXNREF(:,:,:)
        PRS(:,:,:,JSV)  = 0.0
        ZSVS(:,:,:,JSV) = 0.0
      END WHERE
    ENDDO
 IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:)  * PRHODJ(:,:,:), 4,'NECON_BU_RTH')
 IF (LBUDGET_RV) CALL BUDGET (PRS(:,:,:,1) * PRHODJ(:,:,:), 6,'NECON_BU_RRV')
 IF (LBUDGET_RC) CALL BUDGET (PRS(:,:,:,2) * PRHODJ(:,:,:), 7,'NECON_BU_RRC')
END IF
!-------------------------------------------------------------------------------
!
!
!*      13.     SWITCH BACK TO THE PROGNOSTIC VARIABLES
!               ---------------------------------------
!
PTHS(:,:,:) = PTHS(:,:,:) * PRHODJ(:,:,:)
!
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) * PRHODJ(:,:,:)
END DO
!
IF (HCLOUD=='C2R2' .OR. HCLOUD=='C3R5' .OR. HCLOUD=='KHKO' .OR. HCLOUD=='LIMA') THEN
  DO JSV = 1,SIZE(ZSVS,4)
    PSVS(:,:,:,JSV+ISVBEG-1) = ZSVS(:,:,:,JSV) * PRHODJ(:,:,:)
  ENDDO
  DO JSV = 1,SIZE(ZSVT,4)
    PSVT(:,:,:,JSV+ISVBEG-1) = ZSVT(:,:,:,JSV)
  ENDDO
  DEALLOCATE(ZSVS)
  DEALLOCATE(ZSVT)
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RESOLVED_CLOUD
