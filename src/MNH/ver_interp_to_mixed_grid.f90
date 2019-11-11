!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################################
      MODULE MODI_VER_INTERP_TO_MIXED_GRID
!     ####################################
INTERFACE
      SUBROUTINE VER_INTERP_TO_MIXED_GRID(HFILE,OEXTR,PZS_LS,PZSMT_LS,   &
                          PZMASS_LS,PSV_LS,                              &
                          PZFLUX_LS,PPS_LS,PPMHP_LS,                     &
                          PTHV_LS,PR_LS,PHU_LS,PTKE_LS,                  &
                          PU_LS,PV_LS,PW_LS,HWLOC_LS,                    &
                          PLSU_LS,PLSV_LS,PLSW_LS,PLSTH_LS,PLSRV_LS      )

!
CHARACTER(LEN=4),      INTENT(IN) :: HFILE    ! which file ('ATM ' or 'CHEM')
LOGICAL,               INTENT(IN) :: OEXTR    ! T: extremum vertical values are used
!                                             ! F: point under the soil is not used
!
REAL,DIMENSION(:,:),   INTENT(IN) :: PZS_LS   ! orography
REAL,DIMENSION(:,:),   INTENT(IN) :: PZSMT_LS ! smooth orography
REAL,DIMENSION(:,:,:), INTENT(IN) :: PZMASS_LS! altitude of mass points
REAL,DIMENSION(:,:,:,:), INTENT(IN):: PSV_LS  ! scalar mixing ratios
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PZFLUX_LS! altitude of pressure points
REAL,DIMENSION(:,:),   INTENT(IN), OPTIONAL :: PPS_LS   ! surface pressure
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PPMHP_LS ! pressure minus hyd. pressure
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PTHV_LS  ! virtual potential temperature
REAL,DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL :: PR_LS  ! water mixing ratios
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PHU_LS   ! relative humidity
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PTKE_LS  ! turbulence kinetic energy
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PU_LS    ! pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PV_LS    ! pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PW_LS    ! vertical wind speed
CHARACTER(LEN=4),      INTENT(IN), OPTIONAL :: HWLOC_LS ! localisation of vertical wind speed
!
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSU_LS  ! large scale pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSV_LS  ! large scale pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSW_LS  ! large scale vertical wind speed
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSTH_LS ! large scale potential temperature
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSRV_LS ! large scale vapor mixing ratios
!
END SUBROUTINE VER_INTERP_TO_MIXED_GRID
END INTERFACE
END MODULE MODI_VER_INTERP_TO_MIXED_GRID
!     ####################################################################
      SUBROUTINE VER_INTERP_TO_MIXED_GRID(HFILE,OEXTR,PZS_LS,PZSMT_LS,   &
                          PZMASS_LS,PSV_LS,                              &
                          PZFLUX_LS,PPS_LS,PPMHP_LS,                     &
                          PTHV_LS,PR_LS,PHU_LS,PTKE_LS,                  &
                          PU_LS,PV_LS,PW_LS,HWLOC_LS,                    &
                          PLSU_LS,PLSV_LS,PLSW_LS,PLSTH_LS,PLSRV_LS      )
!     ####################################################################
!
!!****  *VER_INTERP_TO_MIXED_GRID* - interpolates onto the mixed grid.
!!
!!    PURPOSE
!!    -------
!!    This routine computes:
!!  1 the interpolation of thetav, rv and covariant wind components to the
!!     mixed grid  (change of discretization)
!!  2 the value of the Exner function at the top of the MESO-NH model
!!  3 the reference state Exner function at MESO-NH top
!!  4 the value of rho dry on the mixed grid
!!
!!**  METHOD
!!    ------
!!
!!
!!  1 thetav, rv, u and v are interpolated on the Gal-Chen Sommerville grid
!!    defined by XZHAT and the grib orography. This is only a change of
!!    discretization.
!!
!!    CAUTION: The upper grigrib mass point is not used, since its Exner
!!             function, and therefore its altitude, could not be correctly
!!             computed.
!!
!!
!!  2 the local top Exner function is integrated from the hydrostatic relation
!!    from bottom to top.
!!
!!  3 the reference anelastic state top Exner function is computed as the
!!    mean value of the local top Exner function.
!!
!!  4 the density of dry air on the mixed grid, according to the perfect gas law
!!
!!                  ~ (Cpd/Rd -1)
!!            p00 (PI)
!!    rhod=  ____________________
!!
!!            Rd thetav (1+rv)
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    routine  VER_INTERP    : to vertically interpolate a field
!!    routine COMPUTE_EXNER_FROM_GROUND : to compute hydrostatic Exner function
!!    Module MODI_VER_INTERP : interface for function VER_INTERP
!!    Module MODI_SHUMAN     : interface for Shuman operators
!!    Module MODI_COMPUTE_EXNER_FROM_GROUND
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB      : verbosity level for output-listing
!!      Module MODD_CONF1
!!         NRR
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_CST       : contains physical constants
!!         XRD : gas constant for dry air
!!         XRV : gas constant for vapor
!!         XP00: reference pressure
!!         XCPD: specific heat for dry air
!!         XG  : gravity constant
!!         XRADIUS : earth radius
!!      Module MODD_REF       : contains anelastic reference state variables.
!!         XEXNTOP
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!         JPVEXT
!!      Module MODD_GRID1
!!         XXHAT
!!         XYHAT
!!         XZHAT
!!         XMAP
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/12/94
!!                  20/08/96 (V. Masson) definition of density of dry air
!!                                         (rw instead of rv)
!!                  29/10/96 (V. Masson) add positivity control on rv
!!                  12/12/96 (V. Masson) add vertical wind velocity
!!                  10/01/96 (V. Masson) S T O P if model top to high
!!                  10/06/97 (V. Masson) add non-hydrostatic pressure
!!                  12/06/97 (V. Masson) relative humidity interpolation
!!                  10/07/97 (V. Masson) add epsilon
!!                  11/07/97 (V. Masson) add scalar variables
!!                  26/08/97 (V. Masson) call to new linear vertical
!!                                       interpolation routine
!!                  21/11/97 (V. Masson) bug in rhod computation (since masdev3_1)
!!                  22/01/01 (D. Gazen)  add MODD_NSV for NSV access
!!                  20/05/06             Remove EPS
!!                  10/04/2014 (J.Escobar &  M.Faivre ) add reprod_sum on XEXNTOP
!!                  24/04/2014 (J.escobar) bypass CRAY internal compiler error on IIJ computation
!!                      2014 (M.Faivre)
!!                  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_THERMO
USE MODE_FM
USE MODE_IO_ll
!
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
USE MODI_SHUMAN
USE MODI_COMPUTE_EXNER_FROM_GROUND
USE MODI_WATER_SUM
USE MODI_VERT_COORD
!
USE MODD_CONF           ! declaration modules
USE MODD_NSV, ONLY : NSV
USE MODD_CONF_n
USE MODD_LUNIT
USE MODD_CST
USE MODD_REF
USE MODD_GRID_n
USE MODD_PARAMETERS
USE MODD_PARAM_n
USE MODD_VER_INTERP_LIN
USE MODD_PREP_REAL
!20131028 add MODD_DIMn to use NIMAX,JMAX
USE MODD_DIM_n
USE MODD_PGDDIM
!20131028 add REPRO_SUM
USE MODE_REPRO_SUM
!JUAN REALZ
USE MODE_ll
USE MODE_EXTRAPOL
!JUAN REALZ
USE MODE_MSG
USE MODE_REPRO_SUM
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
CHARACTER(LEN=4),      INTENT(IN) :: HFILE    ! which file ('ATM ' or 'CHEM')
LOGICAL,               INTENT(IN) :: OEXTR    ! T: extremum vertical values are used
!                                             ! F: point under the soil is not used
!
REAL,DIMENSION(:,:),   INTENT(IN) :: PZS_LS   ! orography
REAL,DIMENSION(:,:),   INTENT(IN) :: PZSMT_LS ! smooth orography
REAL,DIMENSION(:,:,:), INTENT(IN) :: PZMASS_LS! altitude of mass points
REAL,DIMENSION(:,:,:,:), INTENT(IN):: PSV_LS  ! scalar mixing ratios
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PZFLUX_LS! altitude of pressure points
REAL,DIMENSION(:,:),   INTENT(IN), OPTIONAL :: PPS_LS   ! surface pressure
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PPMHP_LS ! pressure minus hyd. pressure
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PTHV_LS  ! virtual potential temperature
REAL,DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL :: PR_LS  ! water mixing ratios
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PHU_LS   ! relative humidity
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PTKE_LS  ! turbulence kinetic energy
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PU_LS    ! pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PV_LS    ! pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PW_LS    ! vertical wind speed
CHARACTER(LEN=4),      INTENT(IN), OPTIONAL :: HWLOC_LS ! localisation of vertical wind speed
!
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSU_LS  ! large scale pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSV_LS  ! large scale pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSW_LS  ! large scale vertical wind speed
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSTH_LS ! large scale potential temperature
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSRV_LS ! large scale vapor mixing ratios
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER                                                  :: IRESP, ILUOUT0
INTEGER                                                  :: IIU,IJU,IKU,ILU
INTEGER                                                  :: IIB,IIE,IJB,IJE,IKE
INTEGER                                                  :: JK,JRR,JSV
INTEGER, DIMENSION(2)                                    :: IIJ
INTEGER                                                  :: IK4000
!                        ! upper level where interpolation occured
INTEGER                                                  :: ISV
REAL,DIMENSION(SIZE(PZMASS_LS,1),SIZE(PZMASS_LS,2),SIZE(XZHAT)):: ZHEXNFLUX_MX,    &
                                                              ZHEXNMASS_MX,    &
                                                              ZHU_MX,ZPMASS_MX
!                        ! local hyd. Exner function at pressure and mass points
!                        ! relative hymidity and pressure
!                        !at mass points on the mixed grid
REAL,DIMENSION(SIZE(PZMASS_LS,1),SIZE(PZMASS_LS,2)):: ZEXNSURF2D_MX
!                        ! local Exner function at ground
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZZFLUX_MX,ZZMASS_MX
!JUAN REALZ
REAL                                                     :: ZCOUNT
INTEGER                                                  :: IINFO_ll
!JUAN REALZ

!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIU=SIZE(PZMASS_LS,1)
IJU=SIZE(PZMASS_LS,2)
IKU=SIZE(XZHAT)
IKE=IKU-JPVEXT
ILU=SIZE(PZMASS_LS,3)
!
!-------------------------------------------------------------------------------
!
!*       1.    INTERPOLATION TO THE MIXED GRID
!              -------------------------------
!
!*       1.1   Grid definition
!              ---------------
!
IF (MINVAL(PZMASS_LS (:,:,ILU))<0.5*(XZHAT(IKE)+XZHAT(IKE+1))) THEN
  WRITE(ILUOUT0,*)
  WRITE(ILUOUT0,*) '+-----------------------------------------------------+'
  WRITE(ILUOUT0,*) '| MESONH highest mass level above highest input level |'
  WRITE(ILUOUT0,*) '+-----------------------------------------------------+'
  WRITE(ILUOUT0,*) 'MESONH top (highest W level)          : ', XZHAT(IKE+1)
  WRITE(ILUOUT0,*) 'MESONH highest physical mass level    : ', 0.5*(XZHAT(IKE)+XZHAT(IKE+1))
  WRITE(ILUOUT0,*) 'input  highest mass level             : ', MINVAL(PZMASS_LS (:,:,ILU))
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','VER_INTERP_TO_MIXED_GRID','MESONH highest mass level above highest input level')
ENDIF
!
IF (HFILE=='ATM ') THEN
  ALLOCATE(XZFLUX_MX(IIU,IJU,IKU))
  ALLOCATE(XZMASS_MX(IIU,IJU,IKU))
  CALL VERT_COORD(LSLEVE,PZS_LS,PZSMT_LS,XLEN1,XLEN2,XZHAT,XZFLUX_MX)
  XZMASS_MX(:,:,:)=MZF(1,IKU,1,XZFLUX_MX)
  XZMASS_MX(:,:,IKU)=1.5*XZFLUX_MX(:,:,IKU)-0.5*XZFLUX_MX(:,:,IKU-1)
ELSE IF (HFILE=='CHEM') THEN
  ALLOCATE(ZZFLUX_MX(IIU,IJU,IKU))
  ALLOCATE(ZZMASS_MX(IIU,IJU,IKU))
  CALL VERT_COORD(LSLEVE,PZS_LS,PZSMT_LS,XLEN1,XLEN2,XZHAT,ZZFLUX_MX)
  ZZMASS_MX(:,:,:)=MZF(1,IKU,1,ZZFLUX_MX)
  ZZMASS_MX(:,:,IKU)=1.5*ZZFLUX_MX(:,:,IKU)-0.5*ZZFLUX_MX(:,:,IKU-1)
END IF
!
!
!*       1.2   Interpolation onto this grid
!              ----------------------------
!
IF (HFILE=='ATM ') THEN
  ALLOCATE(XTHV_MX(IIU,IJU,IKU))
  ALLOCATE(XR_MX(IIU,IJU,IKU,NRR))
  ALLOCATE(XU_MX(IIU,IJU,IKU))
  ALLOCATE(XV_MX(IIU,IJU,IKU))
  ALLOCATE(XW_MX(IIU,IJU,IKU))
  ALLOCATE(XPMHP_MX(IIU,IJU,IKU))
  IF ( PRESENT(PLSTH_LS) ) THEN
    ALLOCATE(XLSTH_MX(IIU,IJU,IKU))
    ALLOCATE(XLSRV_MX(IIU,IJU,IKU))
    ALLOCATE(XLSU_MX(IIU,IJU,IKU))
    ALLOCATE(XLSV_MX(IIU,IJU,IKU))
    ALLOCATE(XLSW_MX(IIU,IJU,IKU))
  END IF
  CALL COEF_VER_INTERP_LIN(PZMASS_LS  (:,:,:),XZMASS_MX  (:,:,:))
ELSE IF (HFILE=='CHEM') THEN
  CALL COEF_VER_INTERP_LIN(PZMASS_LS  (:,:,:),ZZMASS_MX  (:,:,:))
END IF
!
IF (ALLOCATED(XSV_MX)) DEALLOCATE(XSV_MX)
ISV=SIZE(PSV_LS,4)
ALLOCATE(XSV_MX(IIU,IJU,IKU,ISV))
DO JSV=1,ISV
  XSV_MX(:,:,:,JSV)=VER_INTERP_LIN(PSV_LS(:,:,:,JSV),NKLIN(:,:,:),XCOEFLIN(:,:,:))
END DO
!
!
IF (HFILE=='ATM ') THEN
  !
  XPMHP_MX(:,:,:)=VER_INTERP_LIN(PPMHP_LS   (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  XTHV_MX (:,:,:)=VER_INTERP_LIN(PTHV_LS    (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  ZHU_MX  (:,:,:)=VER_INTERP_LIN(PHU_LS     (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  IF (PRESENT(PLSTH_LS)) THEN
    XLSTH_MX(:,:,:)=VER_INTERP_LIN(PLSTH_LS   (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    XLSRV_MX(:,:,:)=VER_INTERP_LIN(PLSRV_LS   (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  END IF
  !
  WHERE( XZMASS_MX(:,:,:)<SPREAD(PZMASS_LS(:,:,1),3,IKU) )
    ZHU_MX(:,:,:) = MIN(ZHU_MX(:,:,:),100.)
  END WHERE
  !
  DO JRR=2,NRR
    XR_MX(:,:,:,JRR)  =VER_INTERP_LIN(PR_LS      (:,:,:,JRR),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    XR_MX(:,:,:,JRR)  =MAX(XR_MX(:,:,:,JRR),0.)
  END DO
  !
  IF (CTURB=='TKEL' ) THEN
    ALLOCATE(XTKE_MX(IIU,IJU,IKU))
    XTKE_MX(:,:,:)=VER_INTERP_LIN(PTKE_LS    (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    XTKE_MX(:,:,:)=MAX(XTKE_MX(:,:,:),0.)
  ELSE
    ALLOCATE(XTKE_MX(0,0,0))
  END IF
  !
  !
  XU_MX(:,:,:)  =VER_INTERP_LIN(PU_LS      (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  XV_MX(:,:,:)  =VER_INTERP_LIN(PV_LS      (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  IF (PRESENT(PLSTH_LS)) THEN
    XLSU_MX(:,:,:)=VER_INTERP_LIN(PLSU_LS   (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    XLSV_MX(:,:,:)=VER_INTERP_LIN(PLSV_LS   (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
  END IF
  !
  IF (HWLOC_LS=='FLUX') THEN
    CALL COEF_VER_INTERP_LIN(PZFLUX_LS  (:,:,:),XZFLUX_MX  (:,:,:))
    XW_MX(:,:,:)  =VER_INTERP_LIN(PW_LS      (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    IF  (PRESENT(PLSTH_LS)) THEN
      XLSW_MX(:,:,:)  =VER_INTERP_LIN(PLSW_LS(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    END IF
  ELSE IF (HWLOC_LS=='MASS') THEN
    CALL COEF_VER_INTERP_LIN(PZMASS_LS  (:,:,:),XZFLUX_MX  (:,:,:))
    XW_MX(:,:,:)  =VER_INTERP_LIN(PW_LS      (:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    IF  (PRESENT(PLSTH_LS)) THEN
      XLSW_MX(:,:,:)  =VER_INTERP_LIN(PLSW_LS(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
    END IF
  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'GEN','VER_INTERP_TO_MIXED_GRID','Bad HWLOC_LS input argument')
  END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION OF THE LOCAL EXNER FUNCTION
!              ---------------------------------------
!
  ZEXNSURF2D_MX(:,:)=(PPS_LS(:,:)/XP00)**(XRD/XCPD)
  CALL COMPUTE_EXNER_FROM_GROUND(XTHV_MX,XZFLUX_MX,ZEXNSURF2D_MX,ZHEXNFLUX_MX,ZHEXNMASS_MX)
  !
  ALLOCATE(XEXNTOP2D(IIU,IJU))
  XEXNTOP2D(:,:)=ZHEXNFLUX_MX(:,:,IKE+1)
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF THE WATER VAPOR MIXING RATIO
!              -------------------------------------------
!
  ZPMASS_MX(:,:,:)=XP00*(ZHEXNMASS_MX(:,:,:))**(XCPD/XRD) + XPMHP_MX(:,:,:)
  !
  XR_MX(:,:,:,1)=SM_PMR_HU(ZPMASS_MX(:,:,:),                                   &
                           XTHV_MX(:,:,:)*(ZPMASS_MX(:,:,:)/XP00)**(XRD/XCPD), &
                           ZHU_MX(:,:,:),XR_MX(:,:,:,:),KITERMAX=100)
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTATION OF THE REFERENCE STATE TOP EXNER FUNCTION
!              -----------------------------------------------------
!
!!$  XEXNTOP=SUM(ZHEXNFLUX_MX(IIB:IIE,IJB:IJE,IKE+1))/FLOAT((IIE-IIB+1)*(IJE-IJB+1))
!JUAN REALZ
!!!  XEXNTOP  = SUM(ZHEXNFLUX_MX(IIB:IIE,IJB:IJE,IKE+1))
!20131028 in Mymodif --> 20131129 in MNHorig
XEXNTOP=SUM_DD_R2_ll(ZHEXNFLUX_MX(IIB:IIE,IJB:IJE,IKE+1))
ZCOUNT   = FLOAT((IIE-IIB+1)*(IJE-IJB+1))
!$20140227 disable reduce no xexntop !!
!$ CALL REDUCESUM_ll(XEXNTOP,IINFO_ll)
  CALL REDUCESUM_ll(ZCOUNT,IINFO_ll)
  XEXNTOP = XEXNTOP / ZCOUNT

!JUAN REALZ
!
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTATION OF THE DENSITY OF DRY AIR
!              -------------------------------------
!
  ALLOCATE(XRHOD_MX(IIU,IJU,IKU))
!
  XRHOD_MX(:,:,:)=ZPMASS_MX(:,:,:)/(ZPMASS_MX(:,:,:)/XP00)**(XRD/XCPD) &
                   /(XRD*XTHV_MX(:,:,:)*(1.+WATER_SUM(XR_MX(:,:,:,:))))
!
 CALL EXTRAPOL('E',XRHOD_MX)
!
!-------------------------------------------------------------------------------
!
!*       6.    PRINTS
!              ------
!
  IF (NVERB>=1 .AND. ANY(XZHAT>=5000.) ) THEN
    IK4000 = COUNT(XZHAT(:)<4000.)
    IK4000 = COUNT(XZHAT(:)<4000.)
    IIJ = MAXLOC(        SUM(ZHU_MX(IIB:IIE,IJB:IJE,JPVEXT+1:IK4000),3),                  &
                  MASK=COUNT(ZHU_MX(IIB:IIE,IJB:IJE,JPVEXT+1:IKE)                         &
                             >=MAXVAL(ZHU_MX(IIB:IIE,IJB:IJE,JPVEXT+1:IKE))-0.01,DIM=3 )  &
                        >=1                                                   )           &
          + JPHEXT
    WRITE(ILUOUT0,*) ' '
    WRITE(ILUOUT0,*) 'Altitude and humidity on large-scale grid (I=',IIJ(1),';J=',IIJ(2),')'
    DO JK=1,ILU
      WRITE(ILUOUT0,'(6Hlevel ,F6.0,5H m : ,F6.2,2H %)') PZMASS_LS(IIJ(1),IIJ(2),JK),PHU_LS(IIJ(1),IIJ(2),JK)
    END DO
    !
    WRITE(ILUOUT0,*) ' '
    WRITE(ILUOUT0,*) 'Altitude and humidity on mixed grid       (I=',IIJ(1),';J=',IIJ(2),')'
    DO JK=JPVEXT+1,IKE
      WRITE(ILUOUT0,'(6Hlevel ,F6.0,5H m : ,F6.2,2H %)') XZMASS_MX(IIJ(1),IIJ(2),JK),ZHU_MX(IIJ(1),IIJ(2),JK)
    END DO
  END IF
!
END IF
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'Routine VER_INTERP_TO_MIXED_GRID for ',HFILE,' completed'
!
END SUBROUTINE VER_INTERP_TO_MIXED_GRID
