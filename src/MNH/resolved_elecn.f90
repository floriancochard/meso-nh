!MNH_LIC Copyright 2009-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!     ###########################
      MODULE MODI_RESOLVED_ELEC_n
!     ###########################
!
INTERFACE
      SUBROUTINE RESOLVED_ELEC_n (HCLOUD, HSCONV, HMF_CLOUD,                              &
                                  KRR, KSPLITR, KMI, KTCOUNT, OEXIT,                      &
                                  HLBCX, HLBCY, HRAD, HTURBDIM,                           &
                                  OSUBG_COND, OSIGMAS,PSIGQSAT, HSUBG_AUCV,               &
                                  PTSTEP, PZZ, PRHODJ, PRHODREF, PEXNREF,                 &
                                  PPABST, PTHT, PTHS, PWT,                                & 
                                  PRT, PRS, PSVT, PSVS, PCIT,                             & 
                                  PSIGS, PSRCS, PCLDFR, PMFCONV, PCF_MF, PRC_MF,          &
                                  PRI_MF, OSEDIC, OWARM,                                  &
                                  PINPRC, PINPRR, PINPRR3D, PEVAP3D,                      &
                                  PINPRS, PINPRG, PINPRH,                                 &
                                  PSEA, PTOWN                                             )   
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud
CHARACTER(LEN=4),         INTENT(IN)   :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HMF_CLOUD! Type of statistical cloud
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                                   ! integrations for  rain sedimendation
INTEGER,                  INTENT(IN)   :: KMI      ! Model index
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
LOGICAL,                  INTENT(IN)   :: OEXIT    ! switch for the end of the temporal loop
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
CHARACTER*4,              INTENT(IN)   :: HRAD     ! Radiation scheme name
CHARACTER*4,              INTENT(IN)   :: HTURBDIM ! Dimensionality of the
                                                   ! turbulence scheme
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                                   ! use values computed in CONDENSATION
                                                   ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                                   ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)   :: PTSTEP ! Double Time step
                                                 ! (single if cold start)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PMFCONV ! convective mass flux
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT  ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS  ! Scalar variable sources
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT  ! Pristine ice number
                                                 ! concentration at time t
!
LOGICAL,                  INTENT(IN) :: OSEDIC! Switch to activate the
                                              ! cloud droplet sedimentation
LOGICAL,                  INTENT(IN) :: OWARM ! Control of the rain formation
                                              !  by slow warm microphysical
                                              !         processes
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRI_MF! Convective Mass Flux solid mixing ratio
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRC   ! Cloud instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRR   ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PINPRR3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D  ! evap profile
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRS   ! Snow instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRG   ! Graupel instant precip
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PINPRH   ! Hail instant precip
!
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA   ! Land Sea mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN  ! Town fraction
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWT  ! vertical velocity at time t-dt
!
END SUBROUTINE RESOLVED_ELEC_n
END INTERFACE
END MODULE MODI_RESOLVED_ELEC_n
!
!     #####################################################################################
      SUBROUTINE RESOLVED_ELEC_n (HCLOUD, HSCONV, HMF_CLOUD,                              &
                                  KRR, KSPLITR, KMI, KTCOUNT, OEXIT,                      &
                                  HLBCX, HLBCY, HRAD, HTURBDIM,                           &
                                  OSUBG_COND, OSIGMAS,PSIGQSAT, HSUBG_AUCV,               &
                                  PTSTEP, PZZ, PRHODJ, PRHODREF, PEXNREF,                 &
                                  PPABST, PTHT, PTHS, PWT,                                & 
                                  PRT, PRS, PSVT, PSVS, PCIT,                             & 
                                  PSIGS, PSRCS, PCLDFR, PMFCONV, PCF_MF, PRC_MF,          &
                                  PRI_MF, OSEDIC, OWARM,                                  &
                                  PINPRC, PINPRR, PINPRR3D, PEVAP3D,                      &
                                  PINPRS, PINPRG, PINPRH,                                 &
                                  PSEA, PTOWN                                             )   
!     #####################################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the resolved clouds and 
!!    precipitation, the associated cloud electrification, and the charge 
!!    neutralization associated to lightning flashes
!!
!!
!!    METHOD
!!    ------
!!      The main action of this routine is to call the routines computing the
!!    microphysical and electrical sources. Before that:
!!        - it computes the real absolute pressure,
!!        - negative values of the current guess of all mixing ratio are removed.
!!          This is done by a global filling algorithm based on a multiplicative
!!          method (Rood, 1987), in order to conserved the total mass in the
!!          simulation domain.
!!        - Sources are transformed in physical tendencies, by removing the
!!          multiplicative term Rhod*J.
!!        - External points values are filled owing to the use of cyclic
!!          l.b.c., in order to performe computations on the full domain.
!!      After calling to microphysical and electrical routines, the physical 
!!    tendencies are switched back to prognostic variables.
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
!!    AUTHOR
!!    ------
!!      C. Barthe       * LACy *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/11/09
!!      Modifications: 
!!      M. Chong      26/01/10  Add Small ions parameters
!!      M. Chong      31/07/14  Add explicit LiNOx
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 10/01/2019: use NEWUNIT argument of OPEN
!  P. Wautelet 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!  P. Wautelet 14/03/2019: bugfix: correct management of files
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ELEC_ll
USE MODE_FM,               ONLY: IO_FILE_CLOSE_ll,IO_FILE_OPEN_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST,IO_FILE_FIND_BYNAME
USE MODE_ll
!
USE MODD_METRICS_n, ONLY : XDXX, XDYY, XDZX, XDZY, XDZZ 
USE MODD_FIELD_n, ONLY : XRSVS
USE MODD_CONF, ONLY : L1D, L2D, CEXP
USE MODD_CST
USE MODD_IO_ll,            ONLY: TFILEDATA, TFILE_DUMMY
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_ELEC_DESCR
USE MODD_ELEC_n
USE MODD_BUDGET
USE MODD_NSV
USE MODD_CH_MNHC_n,    ONLY: LUSECHEM,LCH_CONV_LINOX
USE MODD_DYN_n, ONLY: NSTOP, XTSTEP
USE MODD_ARGSLIST_ll, ONLY : LIST_ll

USE MODD_TIME_n
USE MODD_LMA_SIMULATOR
!
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODI_RAIN_ICE_ELEC
USE MODI_ICE_ADJUST_ELEC
USE MODI_TO_ELEC_FIELD_n
USE MODI_FLASH_GEOM_ELEC_n
USE MODI_SHUMAN
USE MODI_BUDGET
USE MODI_ION_ATTACH_ELEC
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODI_ION_DRIFT
USE MODI_SERIES_CLOUD_ELEC
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=4),         INTENT(IN)   :: HCLOUD   ! kind of cloud
                                                   ! paramerization
CHARACTER(LEN=4),         INTENT(IN)   :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)   :: HMF_CLOUD! Type of statistical cloud
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
                                       ! integrations for  ice  sedimendation
INTEGER,                  INTENT(IN)   :: KMI      ! Model index
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
LOGICAL,                  INTENT(IN)   :: OEXIT    ! switch for the end of the temporal loop
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
CHARACTER*4,              INTENT(IN)   :: HRAD     ! Radiation scheme name
CHARACTER*4,              INTENT(IN)   :: HTURBDIM ! Dimensionality of the
                                                   ! turbulence scheme
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)   :: PTSTEP   ! Double Time step
                                                   ! (single if cold start)
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
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PMFCONV ! convective mass flux
!
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT  ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS  ! Scalar variable sources
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PSRCS ! Second-order flux
                                               ! s'rc'/2Sigma_s2 at time t+1
                                               ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCLDFR! Cloud fraction
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCIT  ! Pristine ice number
                                               ! concentration at time t
LOGICAL,                  INTENT(IN) :: OSEDIC! Switch to activate the
                                              ! cloud droplet sedimentation
                                              ! for ICE3            
LOGICAL,                  INTENT(IN) :: OWARM ! Control of the rain formation
                                              !  by slow warm microphysical
                                              !         processes
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRI_MF! Convective Mass Flux solid mixing ratio
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PINPRC ! Cloud instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PINPRR ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PINPRR3D ! sed flux of precip
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D  ! evap profile
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PINPRS ! Snow instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PINPRG ! Graupel instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PINPRH ! Hail instant precip
!
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA   ! Land Sea mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN  ! Town fraction
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWT  ! vertical velocity at time t-dt
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
INTEGER :: IPROC         ! my proc number
INTEGER :: IERR          ! error status
INTEGER :: ILU           ! unit number for IO
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZT,   &
                                                       ZEXN, &
                                                       ZLV,  &
                                                       ZLS,  &
                                                       ZCPH
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZCOR
                                    ! for the correction of negative rv
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZZZ
                                    ! model layer height
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZQTOT
                                    ! total charge source term
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZION_NUMBER  !nearly Nb
                         !  of elementary charge in hydrometeor charge
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZADD         ! ratio (0
                         !  or 1) of ZION_NUMBER to add to positive
                         !              or negative ion number
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZIONTOT
!
REAL :: ZMASSTOT         ! total mass  for one water category
                         ! including the negative values
REAL :: ZMASSPOS         ! total mass  for one water category
                         ! after removing the negative values
REAL :: ZRATIO           ! ZMASSTOT / ZMASSCOR
!
INTEGER, DIMENSION(3) :: IMINLOC, IMAXLOC
!
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: GMASSCOR ! mask for
                                                               ! mass correction
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: GATTACH  ! mask for
                                     !ion recombination and attachment
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange

INTEGER, DIMENSION(3) :: IM_LOC
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZDRIFT
INTEGER :: IPROCMIN, IK
INTEGER :: IXOR, IYOR  ! origin of the extended subdomain
CHARACTER (LEN=32) :: YASCFILE
!
REAL               :: ZTEMP_DIST
CHARACTER (LEN=18) :: YNAME
LOGICAL            :: GLMA_FILE
LOGICAL,                 SAVE :: GFIRST_CALL = .TRUE.
TYPE(TFILEDATA),POINTER, SAVE :: TZFILE_FGEOM_COORD       => NULL()
TYPE(TFILEDATA),POINTER, SAVE :: TZFILE_FGEOM_DIAG        => NULL()
TYPE(TFILEDATA),POINTER, SAVE :: TZFILE_LMA               => NULL()
TYPE(TFILEDATA),POINTER, SAVE :: TZFILE_SERIES_CLOUD_ELEC => NULL()
!
NULLIFY(TZFIELDS_ll)
!
IF ( GFIRST_CALL ) THEN
  GFIRST_CALL = .FALSE.
  TZFILE_FGEOM_COORD       => TFILE_DUMMY
  TZFILE_FGEOM_DIAG        => TFILE_DUMMY
  TZFILE_LMA               => TFILE_DUMMY
  TZFILE_SERIES_CLOUD_ELEC => TFILE_DUMMY
END IF
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PZZ,3) - JPVEXT
IKU = SIZE(PZZ,3)
!
!
!------------------------------------------------------------------------------
!
!*       2.     MICROPHYSICS AND CLOUD ELECTRIFICATION
!               --------------------------------------
!
!*       2.1    Transformation into physical tendencies
!
! X-Component per m3.s into X-Component per kg.s
PTHS(:,:,:) = PTHS(:,:,:) / PRHODJ(:,:,:)
DO JRR = 1, KRR
  PRS(:,:,:,JRR) = PRS(:,:,:,JRR) / PRHODJ(:,:,:)
END DO
!
DO JSV = NSV_ELECBEG, NSV_ELECEND
  PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) / PRHODJ(:,:,:)
ENDDO
!
!  complete the lateral boundaries to avoid possible problems
!
PTHS(IIB-1,:,:) = PTHS(IIB,:,:)
PTHS(IIE+1,:,:) = PTHS(IIE,:,:)
PTHS(:,IJB-1,:) = PTHS(:,IJB,:)
PTHS(:,IJE+1,:) = PTHS(:,IJE,:)
!
PRS(IIB-1,:,:,1) = PRS(IIB,:,:,1)
PRS(IIE+1,:,:,1) = PRS(IIE,:,:,1)
PRS(:,IJB-1,:,1) = PRS(:,IJB,:,1)
PRS(:,IJE+1,:,1) = PRS(:,IJE,:,1)
!
PRS(IIB-1,:,:,2:) = 0.0
PRS(IIE+1,:,:,2:) = 0.0
PRS(:,IJB-1,:,2:) = 0.0
PRS(:,IJE+1,:,2:) = 0.0
!
! positive ion source
PSVS(IIB-1,:,:,NSV_ELECBEG) = PSVS(IIB,:,:,NSV_ELECBEG)
PSVS(IIE+1,:,:,NSV_ELECBEG) = PSVS(IIE,:,:,NSV_ELECBEG) 
PSVS(:,IJB-1,:,NSV_ELECBEG) = PSVS(:,IJB,:,NSV_ELECBEG)
PSVS(:,IJE+1,:,NSV_ELECBEG) = PSVS(:,IJE,:,NSV_ELECBEG)
! source of hydrometeor charge
PSVS(IIB-1,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0    
PSVS(IIE+1,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0   
PSVS(:,IJB-1,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
PSVS(:,IJE+1,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
! negative ion source
PSVS(IIB-1,:,:,NSV_ELECEND) = PSVS(IIB,:,:,NSV_ELECEND) 
PSVS(IIE+1,:,:,NSV_ELECEND) = PSVS(IIE,:,:,NSV_ELECEND)
PSVS(:,IJB-1,:,NSV_ELECEND) = PSVS(:,IJB,:,NSV_ELECEND)
PSVS(:,IJE+1,:,NSV_ELECEND) = PSVS(:,IJE,:,NSV_ELECEND)
!
!  complete the physical boundaries to avoid some computations
!
IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL')  PRT(IIB-1,:,:,2:) = 0.0
IF(LEAST_ll()  .AND. HLBCX(2) /= 'CYCL')  PRT(IIE+1,:,:,2:) = 0.0
IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL')  PRT(:,IJB-1,:,2:) = 0.0
IF(LNORTH_ll() .AND. HLBCY(2) /= 'CYCL')  PRT(:,IJE+1,:,2:) = 0.0
!
IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL')  &
                        PSVT(IIB-1,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
IF(LEAST_ll()  .AND. HLBCX(2) /= 'CYCL')  &
                        PSVT(IIE+1,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL')  &
                        PSVT(:,IJB-1,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
IF(LNORTH_ll() .AND. HLBCY(2) /= 'CYCL')  &
                        PSVT(:,IJE+1,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
!
!  complete the vertical boundaries
!
PTHS(:,:,IKB-1) = PTHS(:,:,IKB)
PTHS(:,:,IKE+1) = PTHS(:,:,IKE)
!
PRS(:,:,IKB-1,1) = PRS(:,:,IKB,1)
PRS(:,:,IKE+1,1) = PRS(:,:,IKE,1)
PRS(:,:,IKB-1,2:) = 0.0
PRS(:,:,IKE+1,2:) = 0.0
!
PRT(:,:,IKB-1,1) = PRT(:,:,IKB,1)
PRT(:,:,IKE+1,1) = PRT(:,:,IKE,1)
PRT(:,:,IKB-1,2:) = 0.0
PRT(:,:,IKE+1,2:) = 0.0
!
PSVS(:,:,IKB-1,NSV_ELECBEG) = PSVS(:,:,IKB,NSV_ELECBEG)    ! Positive ion
PSVT(:,:,IKB-1,NSV_ELECBEG) = PSVT(:,:,IKB,NSV_ELECBEG)
PSVS(:,:,IKB-1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0          ! Hydrometeor charge
PSVS(:,:,IKE+1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
PSVT(:,:,IKB-1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
PSVT(:,:,IKE+1,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
PSVS(:,:,IKB-1,NSV_ELECEND) = PSVS(:,:,IKB,NSV_ELECEND)    ! Negative ion
PSVT(:,:,IKB-1,NSV_ELECEND) = PSVT(:,:,IKB,NSV_ELECEND)
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
ZEXN(:,:,:) = (PPABST(:,:,:) / XP00)**(XRD / XCPD)
ZT(:,:,:)   = PTHT(:,:,:) * ZEXN(:,:,:)
ZLV(:,:,:)  = XLVTT + (XCPV - XCL) * (ZT(:,:,:) - XTT)
ZLS(:,:,:)  = XLSTT + (XCPV - XCI) * (ZT(:,:,:) - XTT)
ZCPH(:,:,:) = XCPD + XCPV * PTSTEP * PRS(:,:,:,1)
!
!
!------------------------------------------------------------------------------
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
DO JRR = 3, KRR
  SELECT CASE (JRR)
    CASE(3,5,6,7) ! rain, snow, graupel and hail
!
      IF (MIN_ll(PRS(:,:,:,JRR), IINFO_ll) < 0.0) THEN
        GMASSCOR = PRS(:,:,:,JRR) < 0.
!
! compute the total water mass computation
        ZMASSTOT = MAX( 0. , SUM3D_ll( PRS(:,:,:,JRR), IINFO_ll ) )
!
! remove the negative values
        PRS(:,:,:,JRR) = MAX( 0., PRS(:,:,:,JRR) )
!
! compute the new total mass
        ZMASSPOS = MAX(XMNH_TINY,SUM3D_ll( PRS(:,:,:,JRR), IINFO_ll ) )
!
! correct again in such a way to conserve the total mass
        ZRATIO = ZMASSTOT / ZMASSPOS
        PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * ZRATIO
!
! No electric charge without hydrometeors
        WHERE( GMASSCOR )
          PSVS(:,:,:,NSV_ELECBEG-1+JRR) = 0.
        ENDWHERE
      END IF
  END SELECT
END DO
!
!
!*       3.2    Adjustement for liquid and solid cloud
!
WHERE (PRS(:,:,:,4) < 0.)                          ! ice particles
  PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
  PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLS(:,:,:) /  &
                ZCPH(:,:,:) / ZEXN(:,:,:)
  PRS(:,:,:,4) = 0.
!
  ZION_NUMBER(:,:,:) = ABS(PSVS(:,:,:,NSV_ELECBEG+3)) / XECHARGE
  ZADD(:,:,:) = 0.5 + SIGN(0.5, PSVS(:,:,:,NSV_ELECBEG+3))
  PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) +  &
                            ZADD(:,:,:) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) +  &
                            (1. - ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECBEG+3) = 0.0
END WHERE
!
! cloud
WHERE (PRS(:,:,:,2) < 0.)
  PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
  PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
       ZCPH(:,:,:) / ZEXN(:,:,:)
  PRS(:,:,:,2) = 0.
!
  ZION_NUMBER(:,:,:) = ABS(PSVS(:,:,:,NSV_ELECBEG+1)) / XECHARGE
  ZADD(:,:,:) = 0.5 + SIGN(0.5, PSVS(:,:,:,NSV_ELECBEG+1))
  PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) +  &
                            ZADD(:,:,:) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) +  &
                            (1. - ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECBEG+1) = 0.0
END WHERE
!
! if rc or ri are positive, we can correct negative rv
! cloud
WHERE ((PRS(:,:,:,1) < 0.) .AND. (PRS(:,:,:,2) > 0.) )
  PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
  PTHS(:,:,:)  = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
                 ZCPH(:,:,:) / ZEXN(:,:,:)
  PRS(:,:,:,2) = 0.
!
  ZION_NUMBER(:,:,:) = ABS(PSVS(:,:,:,NSV_ELECBEG+1)) / XECHARGE
  ZADD(:,:,:) = 0.5 + SIGN(0.5, PSVS(:,:,:,NSV_ELECBEG+1))
  PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) +  &
                            ZADD(:,:,:) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) +  &
                            (1.-ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECBEG+1) = 0.0
END WHERE
!
! ice
IF(KRR > 3) THEN
  WHERE ((PRS(:,:,:,1) < 0.).AND.(PRS(:,:,:,4) > 0.))
    ZCOR(:,:,:)  = MIN(-PRS(:,:,:,1),PRS(:,:,:,4))
    PRS(:,:,:,1) = PRS(:,:,:,1) + ZCOR(:,:,:)
    PTHS(:,:,:)  = PTHS(:,:,:) - ZCOR(:,:,:) * ZLS(:,:,:) /  &
                   ZCPH(:,:,:) / ZEXN(:,:,:)
    PRS(:,:,:,4) = PRS(:,:,:,4) -ZCOR(:,:,:)
!
  ZION_NUMBER(:,:,:) = ABS(PSVS(:,:,:,NSV_ELECBEG+3)) / XECHARGE
  ZADD(:,:,:) = 0.5 + SIGN(0.5, PSVS(:,:,:,NSV_ELECBEG+3))
  PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) +  &
                            ZADD(:,:,:) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) +  &
                            (1. - ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
  PSVS(:,:,:,NSV_ELECBEG+3) = 0.0
  END WHERE
END IF
!
!
!*       3.3     cascade the electric charges in absence of hydrometeor
!
DO JRR = KRR, 5, -1
  WHERE(PRS(:,:,:,JRR) < XRTMIN_ELEC(JRR))
    PSVS(:,:,:,NSV_ELECBEG-2+JRR) = PSVS(:,:,:,NSV_ELECBEG-2+JRR) + &
                                    PSVS(:,:,:,NSV_ELECBEG-1+JRR)
    PSVS(:,:,:,NSV_ELECBEG-1+JRR) = 0.0
  END WHERE
END DO
JRR = 3
WHERE(PRS(:,:,:,JRR) < XRTMIN_ELEC(JRR))
  PSVS(:,:,:,NSV_ELECBEG-2+JRR) = PSVS(:,:,:,NSV_ELECBEG-2+JRR) + &
                                  PSVS(:,:,:,NSV_ELECBEG-1+JRR)
  PSVS(:,:,:,NSV_ELECBEG-1+JRR) = 0.0
END WHERE
DO JRR = 4, 2, -2
  WHERE(PRS(:,:,:,JRR) < XRTMIN_ELEC(JRR))
!
    ZION_NUMBER(:,:,:) = ABS(PSVS(:,:,:,NSV_ELECBEG-1+JRR)) / XECHARGE
    ZADD(:,:,:) = 0.5 + SIGN(0.5, PSVS(:,:,:,NSV_ELECBEG-1+JRR))
    PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) +  &
                              ZADD(:,:,:) * ZION_NUMBER(:,:,:)
    PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) +  &
                              (1. - ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
    PSVS(:,:,:,NSV_ELECBEG-1+JRR) = 0.0
  END WHERE
END DO
!
!
!*       3.4     store the budget terms
!
IF (LBUDGET_RV) CALL BUDGET (PRS(:,:,:,1) * PRHODJ(:,:,:), 6,'NEGA_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRS(:,:,:,2) * PRHODJ(:,:,:), 7,'NEGA_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRS(:,:,:,3) * PRHODJ(:,:,:), 8,'NEGA_BU_RRR')
IF (LBUDGET_RI) CALL BUDGET (PRS(:,:,:,4) * PRHODJ(:,:,:) ,9,'NEGA_BU_RRI')
IF (LBUDGET_RS) CALL BUDGET (PRS(:,:,:,5) * PRHODJ(:,:,:),10,'NEGA_BU_RRS')
IF (LBUDGET_RG) CALL BUDGET (PRS(:,:,:,6) * PRHODJ(:,:,:),11,'NEGA_BU_RRG')
IF (LBUDGET_RH) CALL BUDGET (PRS(:,:,:,7) * PRHODJ(:,:,:),12,'NEGA_BU_RRH')
IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:)  * PRHODJ(:,:,:), 4,'NEGA_BU_RTH')
!
IF (LBUDGET_SV) THEN
  DO JSV = NSV_ELECBEG, NSV_ELECEND
    CALL BUDGET (PSVS(:,:,:,JSV) * PRHODJ(:,:,:), 12+JSV, 'NEGA_BU_RSV')
  END DO
END IF
!
!
!------------------------------------------------------------------------------
!
!*       4.     ION SOURCE FROM DRIFT MOTION AND COSMIC RAYS
!               ---------------------------------------------
!
!* 	 4.1	Compute the electric field at mass points
!
PSVT(:,:,:,NSV_ELECBEG) = XECHARGE*PSVT(:,:,:,NSV_ELECBEG)    ! 1/kg --> C/kg
PSVT(:,:,:,NSV_ELECEND) =-XECHARGE*PSVT(:,:,:,NSV_ELECEND)
!
CALL TO_ELEC_FIELD_n (PRT, PSVT(:,:,:,NSV_ELECBEG:NSV_ELECEND), PRHODJ, &
                      KTCOUNT, KRR,                                     &
                      XEFIELDU, XEFIELDV, XEFIELDW                      )
!
PSVT(:,:,:,NSV_ELECBEG) = PSVT(:,:,:,NSV_ELECBEG)/XECHARGE    ! back to 1/kg 
PSVT(:,:,:,NSV_ELECEND) =-PSVT(:,:,:,NSV_ELECEND)/XECHARGE
!
!
!*       4.2    Compute source term from -/+(Div (N.mu E)) at mass points, 
!               N positive or negative ion number per kg of air (= PSVT)
!               This is a contribution of drift motion to Source PSVS for ions
!               in 1/(kg.s)
!
CALL MYPROC_ELEC_ll (IPROC)
!
!     Hereafter, ZCPH and ZCOR are used temporarily to store the drift sources 
!          of the positive and negative ions, respectively
!
CALL ION_DRIFT(ZCPH, ZCOR, PSVT, HLBCX, HLBCY)


PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) + ZCPH(:,:,:)
PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) + ZCOR(:,:,:)
!
!*       4.3    Add Cosmic Ray source
!
PSVS(:,:,:,NSV_ELECBEG) = PSVS(:,:,:,NSV_ELECBEG) +           &
                              XIONSOURCEFW(:,:,:) / PRHODREF(:,:,:)
PSVS(:,:,:,NSV_ELECEND) = PSVS(:,:,:,NSV_ELECEND) +           &
                              XIONSOURCEFW(:,:,:) / PRHODREF(:,:,:)
!
!-------------------------------------------------------------------------------
!
SELECT CASE (HCLOUD)
!
  CASE ('ICE3')
!
!*       5.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!*       5.1    Compute the explicit microphysical sources and
!*              the explicit charging rates
!
    CALL RAIN_ICE_ELEC (OSEDIC, HSUBG_AUCV, OWARM,                            &
                        KSPLITR, PTSTEP, KMI, KRR,                            &
                        PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                        PTHT, PRT(:,:,:,1), PRT(:,:,:,2), PRT(:,:,:,3),       &
                        PRT(:,:,:,4), PRT(:,:,:,5), PRT(:,:,:,6),             &
                        PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                        PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                        PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                        PINPRS, PINPRG, PSIGS,                                &
                        PSVT(:,:,:,NSV_ELECBEG),   PSVT(:,:,:,NSV_ELECBEG+1), &
                        PSVT(:,:,:,NSV_ELECBEG+2), PSVT(:,:,:,NSV_ELECBEG+3), &
                        PSVT(:,:,:,NSV_ELECBEG+4), PSVT(:,:,:,NSV_ELECBEG+5), &
                        PSVT(:,:,:,NSV_ELECEND),                              &
                        PSVS(:,:,:,NSV_ELECBEG),   PSVS(:,:,:,NSV_ELECBEG+1), &
                        PSVS(:,:,:,NSV_ELECBEG+2), PSVS(:,:,:,NSV_ELECBEG+3), &
                        PSVS(:,:,:,NSV_ELECBEG+4), PSVS(:,:,:,NSV_ELECBEG+5), &
                        PSVS(:,:,:,NSV_ELECEND),                              &
                        PSEA, PTOWN                                           )

!
!*       5.2    Perform the saturation adjustment over cloud ice and cloud water
!
    ZZZ = MZF(1,IKU,1, PZZ )
    CALL ICE_ADJUST_ELEC (KRR, KMI, HRAD, HTURBDIM,                           &
                          HSCONV, HMF_CLOUD,                                  &
                          OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,               &
                          PRHODJ, PEXNREF, PSIGS, PPABST, ZZZ,                &
                          PMFCONV, PCF_MF, PRC_MF, PRI_MF,                    &
                          PRVT=PRT(:,:,:,1), PRCT=PRT(:,:,:,2),               &
                          PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),               &
                          PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR,              &
                          PRRT=PRT(:,:,:,3), PRRS=PRS(:,:,:,3),               &
                          PRIT=PRT(:,:,:,4), PRIS=PRS(:,:,:,4),               &
                          PRST=PRT(:,:,:,5), PRSS=PRS(:,:,:,5),               &
                          PRGT=PRT(:,:,:,6), PRGS=PRS(:,:,:,6),               &
                          PQPIT=PSVT(:,:,:,NSV_ELECBEG),  & !..PI.. Positive
                          PQPIS=PSVS(:,:,:,NSV_ELECBEG),  & !  Ion Mixing Ratio
                          PQCT=PSVT(:,:,:,NSV_ELECBEG+1), &
                          PQCS=PSVS(:,:,:,NSV_ELECBEG+1), &
                          PQRT=PSVT(:,:,:,NSV_ELECBEG+2), &
                          PQRS=PSVS(:,:,:,NSV_ELECBEG+2), &
                          PQIT=PSVT(:,:,:,NSV_ELECBEG+3), &
                          PQIS=PSVS(:,:,:,NSV_ELECBEG+3), &
                          PQST=PSVT(:,:,:,NSV_ELECBEG+4), &
                          PQSS=PSVS(:,:,:,NSV_ELECBEG+4), &
                          PQGT=PSVT(:,:,:,NSV_ELECBEG+5), &
                          PQGS=PSVS(:,:,:,NSV_ELECBEG+5), &
                          PQNIT=PSVT(:,:,:,NSV_ELECEND),  & !..NI.. Negative
                          PQNIS=PSVS(:,:,:,NSV_ELECEND))    !  Ion Mixing Ratio
!
!
!-------------------------------------------------------------------------------
!
!*       6.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 4 ICE SPECIES)
!               -----------------------------------------------------
!
!*       6.1    Compute the explicit microphysical sources and
!*              the explicit charging rates
!
  CASE ('ICE4')
!
    CALL RAIN_ICE_ELEC (OSEDIC, HSUBG_AUCV, OWARM,                            &
                        KSPLITR, PTSTEP, KMI, KRR,                            &
                        PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                        PTHT, PRT(:,:,:,1), PRT(:,:,:,2), PRT(:,:,:,3),       &
                        PRT(:,:,:,4), PRT(:,:,:,5), PRT(:,:,:,6),             &
                        PTHS, PRS(:,:,:,1), PRS(:,:,:,2), PRS(:,:,:,3),       &
                        PRS(:,:,:,4), PRS(:,:,:,5), PRS(:,:,:,6),             &
                        PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                        PINPRS, PINPRG, PSIGS,                   &
                        PSVT(:,:,:,NSV_ELECBEG),   PSVT(:,:,:,NSV_ELECBEG+1), &
                        PSVT(:,:,:,NSV_ELECBEG+2), PSVT(:,:,:,NSV_ELECBEG+3), &
                        PSVT(:,:,:,NSV_ELECBEG+4), PSVT(:,:,:,NSV_ELECBEG+5), &
                        PSVT(:,:,:,NSV_ELECEND),          &
                        PSVS(:,:,:,NSV_ELECBEG),   PSVS(:,:,:,NSV_ELECBEG+1), &
                        PSVS(:,:,:,NSV_ELECBEG+2), PSVS(:,:,:,NSV_ELECBEG+3), &
                        PSVS(:,:,:,NSV_ELECBEG+4), PSVS(:,:,:,NSV_ELECBEG+5), &
                        PSVS(:,:,:,NSV_ELECEND),          &
                        PSEA, PTOWN,                   &
                        PRT(:,:,:,7), PRS(:,:,:,7), PINPRH,                   &
                        PSVT(:,:,:,NSV_ELECBEG+6), PSVS(:,:,:,NSV_ELECBEG+6)  )
! Index NSV_ELECBEG: Positive ion , NSV_ELECEND: Negative ion
!
!
!*       6.2    Perform the saturation adjustment over cloud ice and cloud water
!
    ZZZ = MZF(1,IKU,1, PZZ )
    CALL ICE_ADJUST_ELEC (KRR, KMI, HRAD,                                     &
                          HTURBDIM, HSCONV, HMF_CLOUD,                        &
                          OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,               &
                          PRHODJ, PEXNREF, PSIGS, PPABST, ZZZ,                &
                          PMFCONV, PCF_MF, PRC_MF, PRI_MF,                    & 
                          PRVT=PRT(:,:,:,1), PRCT=PRT(:,:,:,2),               &
                          PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2),               &
                          PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR,              &
                          PRRT=PRT(:,:,:,3), PRRS=PRS(:,:,:,3),               &
                          PRIT=PRT(:,:,:,4), PRIS=PRS(:,:,:,4),               &
                          PRST=PRT(:,:,:,5), PRSS=PRS(:,:,:,5),               &
                          PRGT=PRT(:,:,:,6), PRGS=PRS(:,:,:,6),               &
                          PQPIT=PSVT(:,:,:,NSV_ELECBEG), &  !..PI.. Positive
                          PQPIS=PSVS(:,:,:,NSV_ELECBEG),  & !  Ion Mixing Ratio
                          PQCT=PSVT(:,:,:,NSV_ELECBEG+1), &
                          PQCS=PSVS(:,:,:,NSV_ELECBEG+1), &
                          PQRT=PSVT(:,:,:,NSV_ELECBEG+2), &
                          PQRS=PSVS(:,:,:,NSV_ELECBEG+2), &
                          PQIT=PSVT(:,:,:,NSV_ELECBEG+3), &
                          PQIS=PSVS(:,:,:,NSV_ELECBEG+3), &
                          PQST=PSVT(:,:,:,NSV_ELECBEG+4), &
                          PQSS=PSVS(:,:,:,NSV_ELECBEG+4), &
                          PQGT=PSVT(:,:,:,NSV_ELECBEG+5), &
                          PQGS=PSVS(:,:,:,NSV_ELECBEG+5), &
                          PQNIT=PSVT(:,:,:,NSV_ELECEND),  & !..NI.. Negative
                          PQNIS=PSVS(:,:,:,NSV_ELECEND),  & !  Ion Mixing Ratio
                          PRHT=PRT(:,:,:,7), PRHS=PRS(:,:,:,7),               &
                          PQHT=PSVT(:,:,:,NSV_ELECBEG+6), &
                          PQHS=PSVS(:,:,:,NSV_ELECBEG+6) )
!
END SELECT
!
IF(KTCOUNT .EQ. 1 .AND. IPROC .EQ. 0) PRINT *,'KSPLITR=', KSPLITR
!
!-------------------------------------------------------------------------------
!
!*      7.      SWITCH BACK TO THE PROGNOSTIC VARIABLES
!               ---------------------------------------
!
! Convert source into component per m3 of air and sec., i.e. volumetric source
!
PTHS(:,:,:) = PTHS(:,:,:) * PRHODJ(:,:,:)
!
DO JRR = 1,KRR
  PRS(:,:,:,JRR)  = PRS(:,:,:,JRR) * PRHODJ(:,:,:)
END DO
!
DO JSV = NSV_ELECBEG, NSV_ELECEND
  PSVS(:,:,:,JSV) = PSVS(:,:,:,JSV) * PRHODJ(:,:,:)
ENDDO
!
! Note that the LiNOx Conc. (in mol/mol) is PSVS (:,::,NSV_LNOXBEG)
! but there is no need to *PRHODJ(:,:,:) as it is done implicitly
! during unit conversion in flash_geom.
!
PSVS(:,:,:,NSV_ELECBEG) = MAX(0., PSVS(:,:,:,NSV_ELECBEG))
PSVS(:,:,:,NSV_ELECEND) = MAX(0., PSVS(:,:,:,NSV_ELECEND))
!
!-------------------------------------------------------------------------------
!
!*      8.      ION RECOMBINATION AND ATTACHMENT
!               --------------------------------
!
GATTACH(:,:,:) = .FALSE.
GATTACH(IIB:IIE, IJB:IJE, IKB:IKE) = .TRUE.
!
IF (PRESENT(PSEA)) THEN
  CALL ION_ATTACH_ELEC(KTCOUNT, KRR, PTSTEP, PRHODREF,                   &
                       PRHODJ, PSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),      & 
                       PRS, PTHT, PCIT, PPABST, XEFIELDU,                &
                       XEFIELDV, XEFIELDW, GATTACH, PTOWN, PSEA          )
ELSE
  CALL ION_ATTACH_ELEC(KTCOUNT, KRR, PTSTEP, PRHODREF,                   &
                       PRHODJ, PSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND),      &
                       PRS, PTHT, PCIT, PPABST, XEFIELDU,                &
                       XEFIELDV, XEFIELDW, GATTACH                       )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      9.      OPEN THE OUTPUT ASCII FILES
!               ---------------------------
!
IF (KTCOUNT==1 .AND. IPROC==0) THEN
  IF (LFLASH_GEOM) THEN
    YASCFILE = CEXP//"_fgeom_diag.asc"
    TZFILE_FGEOM_DIAG => NULL()
    CALL IO_FILE_ADD2LIST(TZFILE_FGEOM_DIAG,YASCFILE,'TXT','WRITE')
    CALL IO_FILE_OPEN_ll(TZFILE_FGEOM_DIAG,HPOSITION='APPEND',HSTATUS='NEW')
    ILU = TZFILE_FGEOM_DIAG%NLU
    WRITE (UNIT=ILU, FMT='(A)') '--------------------------------------------------------'
    WRITE (UNIT=ILU, FMT='(A)') '*FLASH CHARACTERISTICS FROM FLASH_GEOM_ELEC*'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 1 : total flash number          --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 2 : time (s)                    --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 3 : cell number                 --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 4 : flash number/cell/time step --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 5 : flash type 1=IC, 2=CGN, 3=CGP '
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 6 : number of segments          --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 7 : trig electric field (kV/m)  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 8 : x coord. trig. point        --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 9 : y coord. trig. point        --'
    WRITE (UNIT=ILU, FMT='(A)') '--         --> x,y in km if lcartesian=t, --'
    WRITE (UNIT=ILU, FMT='(A)') '--                    deg otherwise       --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 10 : z coord. trig. point (km)  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 11: neutr. positive charge (C)  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 12: neutr. negative charge (C)  --'
    WRITE (UNIT=ILU, FMT='(A)') '--------------------------------------------'
    FLUSH(UNIT=ILU)
!
    IF (LSAVE_COORD) THEN
      YASCFILE = CEXP//"_fgeom_coord.asc"
      TZFILE_FGEOM_COORD => NULL()
      CALL IO_FILE_ADD2LIST(TZFILE_FGEOM_COORD,YASCFILE,'TXT','WRITE')
      CALL IO_FILE_OPEN_ll(TZFILE_FGEOM_COORD,HPOSITION='APPEND',HSTATUS='NEW')
      ILU = TZFILE_FGEOM_COORD%NLU
      WRITE (UNIT=ILU,FMT='(A)') '------------------------------------------'
      WRITE (UNIT=ILU,FMT='(A)') '*****FLASH COORD. FROM FLASH_GEOM_ELEC****'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 1 : flash number             --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 2 : time (s)                 --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 3 : type                     --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 4 : coordinate along X (km)  --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 5 : coordinate along Y (km)  --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 6 : coordinate along Z (km)  --'
      WRITE (UNIT=ILU,FMT='(A)') '------------------------------------------'
      FLUSH(UNIT=ILU)
    END IF
  END IF
!
  IF (LSERIES_ELEC) THEN
    YASCFILE = CEXP//"_series_cloud_elec.asc"
    TZFILE_SERIES_CLOUD_ELEC => NULL()
    CALL IO_FILE_ADD2LIST(TZFILE_SERIES_CLOUD_ELEC,YASCFILE,'TXT','WRITE')
    CALL IO_FILE_OPEN_ll(TZFILE_SERIES_CLOUD_ELEC,HPOSITION='APPEND',HSTATUS='NEW')
    ILU = TZFILE_SERIES_CLOUD_ELEC%NLU
    WRITE (UNIT=ILU, FMT='(A)') '----------------------------------------------------'
    WRITE (UNIT=ILU, FMT='(A)') '********* RESULTS FROM of LSERIES_ELEC *************'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 1 : Time (s)                            --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 2 : Cloud top height / Z > 20 dBZ (m)   --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 3 : Cloud top height / m.r. > 1.e-4 (m) --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 4 : Maximum radar reflectivity (dBZ)    --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 5 : Maximum vertical velocity (m/s)     --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 6 : Updraft volume for W > 5 m/s        --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 7 : Updraft volume for W > 10 m/s       --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 8 : Cloud water mass (kg)               --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 9 : Rain water mass (kg)                --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 10 : Ice crystal mass (kg)              --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 11 : Snow mass (kg)                     --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 12 : Graupel mass (kg)                  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 13 : Precipitation ice mass (kg)        --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 14 : Ice mass flux product (kg2 m2/s2)  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 15 : Precip. ice mass flux (kg m/s)     --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 16 : Non-precip. ice mass flux (kg m/s) --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 17 : Ice water path (kg/m2)             --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 18 : Cloud volume (m3)                  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 19 : Maximum rain inst. precip. (mm/H)  --'
    WRITE (UNIT=ILU, FMT='(A)') '-- Column 20 : Rain instant. precip. (mm/H)       --'
    WRITE (UNIT=ILU, FMT='(A)') '----------------------------------------------------'
    FLUSH(UNIT=ILU)
  END IF
END IF
!
IF (LFLASH_GEOM .AND. LLMA) THEN
!
! test to see if a new LMA file should be created
!
  GLMA_FILE = .FALSE.
!
  IF (TDTCUR%TIME >= TDTLMA%TIME-PTSTEP .AND. CLMA_FILE(1:5) /= "BEGIN") THEN
    LWRITE_LMA  = .TRUE.
  END IF
  IF (TDTCUR%TIME >= TDTLMA%TIME) THEN
    TDTLMA%TIME = TDTLMA%TIME + XDTLMA
    GLMA_FILE   = .TRUE.
    LWRITE_LMA  = .FALSE.
  END IF
!
  IF (GLMA_FILE) THEN
    IF(CLMA_FILE(1:5) /= "BEGIN") THEN ! close previous file if exists
      CALL IO_FILE_FIND_BYNAME(CLMA_FILE,TZFILE_LMA,IERR)
      CALL IO_FILE_CLOSE_ll(TZFILE_LMA)
      TZFILE_LMA => NULL()
    ENDIF
!
    TDTLMA%TIME = TDTLMA%TIME - XDTLMA
    WRITE (YNAME,FMT='(3I2.2,A1,3I2.2,A1,I4.4)')                               &
          ABS(TDTCUR%TDATE%YEAR-2000),TDTCUR%TDATE%MONTH,TDTCUR%TDATE%DAY,'_', &
            INT(TDTLMA%TIME/3600.),INT(FLOAT(MOD(INT(TDTLMA%TIME),3600))/60.), &
                                      MOD(INT(TDTLMA%TIME),60), '_', INT(XDTLMA)
    TDTLMA%TIME = MOD(TDTLMA%TIME + XDTLMA,86400.)
    CLMA_FILE = CEXP//"_SIMLMA_"//YNAME//".dat"
!
    IF ( IPROC .EQ. 0 ) THEN
      TZFILE_LMA => NULL()
      CALL IO_FILE_ADD2LIST(TZFILE_LMA,CLMA_FILE,'TXT','WRITE')
      CALL IO_FILE_OPEN_ll(TZFILE_LMA,HPOSITION='APPEND',HSTATUS='NEW')
      ILU = TZFILE_LMA%NLU
      WRITE (UNIT=ILU,FMT='(A)') '----------------------------------------'
      WRITE (UNIT=ILU,FMT='(A)') '*** FLASH COORD. FROM LMA SIMULATOR ****'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 1  : flash number           --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 2  : time (s)               --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 3  : type                   --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 4  : coordinate along X (km)--'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 5  : coordinate along Y (km)--'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 6  : coordinate along Z (km)--'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 7  : cld drop. mixing ratio --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 8  : rain mixing ratio      --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 9  : ice cryst mixing ratio --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 10 : snow mixing ratio      --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 11 : graupel mixing ratio   --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 12 : rain charge neut       --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 13 : ice cryst. charge neut --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 14 : snow charge neut       --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 15 : graupel charge neut    --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 16 : positive ions neut     --'
      WRITE (UNIT=ILU,FMT='(A)') '-- Column 17 : negative ions neut     --'
      WRITE (UNIT=ILU,FMT='(A)') '----------------------------------------'
      FLUSH(UNIT=ILU)
    END IF
  END IF
END IF
!
!
!-------------------------------------------------------------------------------
!
!*      10.     LIGHTNING FLASHES AND NOX PRODUCTION
!               ------------------------------------
!
! the lightning scheme is now called at each time step
! but only if there's electric charge in the domain
!
ZQTOT(:,:,:) = XECHARGE * (PSVT(:,:,:,NSV_ELECBEG) - PSVT(:,:,:,NSV_ELECEND))
DO JSV = NSV_ELECBEG+1, NSV_ELECEND-1
  ZQTOT(:,:,:) = ZQTOT(:,:,:) + PSVT(:,:,:,JSV)
END DO
!
IF ((.NOT. LOCG) .AND. LELEC_FIELD .AND.  MAX_ll(ABS(ZQTOT),IINFO_ll)>0.) THEN
  IF (LFLASH_GEOM) THEN
    CALL FLASH_GEOM_ELEC_n (KTCOUNT, KMI, KRR, PTSTEP, OEXIT,                                 &
                            PRHODJ, PRHODREF, PRT, PCIT, PSVS(:,:,:,NSV_ELECBEG:NSV_ELECEND), &
                            PRS, PTHT, PPABST, XEFIELDU, XEFIELDV, XEFIELDW,                  &
                            PZZ, PSVS(:,:,:,NSV_LNOXBEG),                                     &
                            TZFILE_FGEOM_DIAG, TZFILE_FGEOM_COORD, TZFILE_LMA,                &
                            PTOWN, PSEA)
  END IF
!
  PSVS(:,:,:,NSV_ELECBEG) = MAX(0., PSVS(:,:,:,NSV_ELECBEG))
  PSVS(:,:,:,NSV_ELECEND) = MAX(0., PSVS(:,:,:,NSV_ELECEND))
!
END IF
!
!
!-------------------------------------------------------------------------------
!
!*      11.     LOOK FOR FLASH RATE PROXIES
!               ---------------------------
!
IF (LSERIES_ELEC) THEN
  CALL SERIES_CLOUD_ELEC (KTCOUNT, PTSTEP,                &
                          PZZ, PRHODJ, PRHODREF, PEXNREF, &
                          PRT, PRS, PSVT,                 &
                          PTHT, PWT, PPABST, PCIT,        &
                          TZFILE_SERIES_CLOUD_ELEC,       &
                          PINPRR                          )
END IF
!
!
!-------------------------------------------------------------------------------
!
!   Close Ascii Files if KTCOUNT = NSTOP
IF (OEXIT .AND. IPROC==0) THEN
  IF (LFLASH_GEOM) THEN
    CALL IO_FILE_CLOSE_ll(TZFILE_FGEOM_DIAG)
    TZFILE_FGEOM_DIAG => TFILE_DUMMY
  END IF
  IF (LFLASH_GEOM .AND. LSAVE_COORD) THEN
    CALL IO_FILE_CLOSE_ll(TZFILE_FGEOM_COORD)
    TZFILE_FGEOM_COORD => TFILE_DUMMY
  END IF
  IF (LSERIES_ELEC) THEN
    CALL IO_FILE_CLOSE_ll(TZFILE_SERIES_CLOUD_ELEC)
    TZFILE_SERIES_CLOUD_ELEC => TFILE_DUMMY
  END IF
  IF (LFLASH_GEOM .AND. LLMA) THEN
    CALL IO_FILE_CLOSE_ll(TZFILE_LMA)
    TZFILE_LMA => TFILE_DUMMY
  END IF
ENDIF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RESOLVED_ELEC_n
