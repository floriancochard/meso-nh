!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision
! MASDEV4_7 prep_real 2006/05/23 14:49:51
!-----------------------------------------------------------------
!     ################################
      MODULE MODI_VER_PREP_GRIBEX_CASE
!     ################################
INTERFACE
      SUBROUTINE VER_PREP_GRIBEX_CASE(HFILE,PDIAG)
!
CHARACTER(LEN=4),      INTENT(IN) :: HFILE    ! which file ('ATM0','ATM1' or 'CHEM')
REAL, INTENT(OUT)                 :: PDIAG    ! diagnostics computing time
!
END SUBROUTINE VER_PREP_GRIBEX_CASE
END INTERFACE
END MODULE MODI_VER_PREP_GRIBEX_CASE
!     ####################################################################
      SUBROUTINE VER_PREP_GRIBEX_CASE(HFILE,PDIAG)
!     ####################################################################
!
!!****  *VER_PREP_GRIBEX_CASE* - monitors the preparation to orographic change
!!
!!    PURPOSE
!!    -------
!!    This routine monitors the preparation of variables to future change
!!    of orography, according to the type of input file.
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    function MZF
!!    routine VER_INTERP_TO_MIXED_GRID
!!    routine CHANGE_GRIBEX_VAR
!!
!!    module MODI_SHUMAN
!!    module MODI_VER_INTERP_TO_MIXED_GRID
!!    module MODI_CHANGE_GRIBEX_VAR
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF1     : contains configuration variables for all models.
!!         NVERB      : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_CST       : contains physical constants
!!         XRD : gas constant for dry air
!!         XRV : gas constant for vapor
!!         XP00: reference pressure
!!         XCPD: specific heat for dry air
!!         XG  : gravity constant
!!         XRADIUS : earth radius
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
!!                  Jan, 31 1996 (V. Masson) duplication of the routine
!!                               to accept different input fields
!!                  May, 25 1996 (V. Masson) take into account the upper level
!!                  Aug, 20 1996 (V. Masson) correction on theta
!!                  Oct, 20 1996 (V. Masson) add deallocations
!!                  Dec, 06 1996 (V. Masson) add air temperature at ground
!!                  Dec, 12 1996 (V. Masson) add vertical wind velocity
!!                  May, 07 1997 (V. Masson) add null tke
!!                  Jun, 10 1997 (V. Masson) add null difference between
!!                                           pressure and hydrostatic pressure
!!                  Jul, 11 1997 (V. Masson) add null scalar variables
!!                  Nov, 22 2000 (I. Mallet) add scalar variables
!!                  Nov, 22 2000 (P. Jabouille) change routine name
!!                  May 2006                 Remove EPS
!!                  Apr, 09 2018 (J.-P. Chaboureau) add isobaric surface 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS, ONLY : JPVEXT, XUNDEF
USE MODD_PREP_REAL
!
USE MODE_THERMO
!
USE MODI_CHANGE_GRIBEX_VAR
USE MODI_COMPUTE_EXNER_FROM_TOP
USE MODI_RMS_AT_Z
USE MODI_SECOND_MNH
USE MODI_SHUMAN
USE MODI_VER_INTERP_TO_MIXED_GRID
USE MODI_WATER_SUM
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
CHARACTER(LEN=4),      INTENT(IN) :: HFILE    ! which file ('ATM0','ATM1' or 'CHEM')
REAL, INTENT(OUT)                 :: PDIAG    ! diagnostics computing time
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER                            :: ILUOUT0
INTEGER                            :: IIU,IJU,ILU
REAL                               :: ZTIME1, ZTIME2
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZTH_LS    ! potential temperature
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZTH_MX    ! potential temperature
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZPMASS_MX    ! pressure
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHEXNFLUX_MX ! pressure function
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHEXNMASS_MX ! pressure function
!
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZZFLUX_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZZMASS_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZPMHP_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZTHV_LS
REAL,DIMENSION(:,:,:,:),ALLOCATABLE:: ZR_LS
REAL,DIMENSION(:,:,:,:),ALLOCATABLE:: ZSV_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHU_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZU_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZV_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZW_LS
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZTKE_LS
INTEGER                            :: JRR     ! loop counter
INTEGER                            :: JSV     ! loop counter
INTEGER                            :: JK      ! loop counter
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!*       1.    CHANGING OF VARIABLES
!              ---------------------
!
IF (HFILE(1:3)=='ATM') THEN
  IIU=SIZE(XT_LS,1)
  IJU=SIZE(XT_LS,2)
  ILU=SIZE(XT_LS,3)
ELSE IF (HFILE=='CHEM') THEN
  IIU=SIZE(XT_SV_LS,1)
  IJU=SIZE(XT_SV_LS,2)
  ILU=SIZE(XT_SV_LS,3)
ELSE
  WRITE (ILUOUT0,'(A)') ' -> Bad input argument in ver_prep_gribex_case - abort'
END IF
!
!
IF (HFILE(1:3)=='ATM') THEN
  ALLOCATE(XPMASS_LS(IIU,IJU,ILU))
  ALLOCATE(XZMASS_LS(IIU,IJU,ILU),XZFLUX_LS(IIU,IJU,ILU))
  ALLOCATE(XTHV_LS(IIU,IJU,ILU),XR_LS(IIU,IJU,ILU,NRR),XHU_LS(IIU,IJU,ILU))
  ALLOCATE(XW_LS(IIU,IJU,ILU))
  CALL CHANGE_GRIBEX_VAR(XA_LS,XB_LS,XP00_LS,XPS_LS,XZS_LS,         &
                         XT_LS,XQ_LS,XPMASS_LS,XZFLUX_LS,XZMASS_LS, &
                         XTHV_LS,XR_LS,XHU_LS,XU_LS,XV_LS,XW_LS     )
ELSE IF (HFILE=='CHEM') THEN
  ALLOCATE(XPMASS_SV_LS(IIU,IJU,ILU))
  ALLOCATE(XZMASS_SV_LS(IIU,IJU,ILU),XZFLUX_SV_LS(IIU,IJU,ILU))
  ALLOCATE(XTHV_SV_LS(IIU,IJU,ILU),XR_SV_LS(IIU,IJU,ILU,NRR),XHU_SV_LS(IIU,IJU,ILU))
  CALL CHANGE_GRIBEX_VAR(XA_SV_LS,XB_SV_LS,XP00_SV_LS,XPS_SV_LS,XZS_SV_LS, &
                         XT_SV_LS,XQ_SV_LS,XPMASS_SV_LS,XZFLUX_SV_LS,XZMASS_SV_LS, &
                         XTHV_SV_LS,XR_SV_LS,XHU_SV_LS                      )
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    INTERPOLATION TO MIXED GRID AND DIAGNOSTIC VARIABLES
!              ----------------------------------------------------
!
IF (HFILE(1:3)=='ATM') THEN
  ALLOCATE(XPMHP_LS(IIU,IJU,ILU))
  XPMHP_LS(:,:,:)=0.
END IF
!
!* Add extra points below and above grids, in order to use MESONH linear
! vertical interpolation programs with all ILU physical points
!
ALLOCATE(ZZMASS_LS(IIU,IJU,ILU+2*JPVEXT))
ALLOCATE(ZSV_LS(IIU,IJU,ILU+2*JPVEXT,SIZE(XSV_LS,4)))
IF (HFILE(1:3)=='ATM') THEN
  ALLOCATE(ZZFLUX_LS(IIU,IJU,ILU+2*JPVEXT))
  ALLOCATE(ZPMHP_LS(IIU,IJU,ILU+2*JPVEXT))
  ALLOCATE(ZTHV_LS(IIU,IJU,ILU+2*JPVEXT))
  ALLOCATE(ZR_LS(IIU,IJU,ILU+2*JPVEXT,NRR))
  ALLOCATE(ZHU_LS(IIU,IJU,ILU+2*JPVEXT))
  ALLOCATE(ZU_LS(IIU,IJU,ILU+2*JPVEXT))
  ALLOCATE(ZV_LS(IIU,IJU,ILU+2*JPVEXT))
  ALLOCATE(ZW_LS(IIU,IJU,ILU+2*JPVEXT))
  IF (SIZE(XTKE_LS)>0) THEN
    ALLOCATE(ZTKE_LS(IIU,IJU,ILU+2*JPVEXT))
  ELSE
    ALLOCATE(ZTKE_LS(0,0,0))
  END IF
END IF
!
IF (HFILE(1:3)=='ATM') THEN
  ZZMASS_LS (:,:,JPVEXT+1:JPVEXT+ILU) = XZMASS_LS(:,:,:)
  DO JK=1,JPVEXT
    ZZMASS_LS(:,:,           JK) = XZMASS_LS(:,:,1)   - (XZMASS_LS(:,:,2)  -XZMASS_LS(:,:,1)    )*(JPVEXT+1-JK)
    ZZMASS_LS(:,:,ILU+JPVEXT+JK) = XZMASS_LS(:,:,ILU) + (XZMASS_LS(:,:,ILU)-XZMASS_LS(:,:,ILU-1))*          JK
  END DO
ELSE IF (HFILE=='CHEM') THEN
  ZZMASS_LS (:,:,JPVEXT+1:JPVEXT+ILU) = XZMASS_SV_LS(:,:,:)
  DO JK=1,JPVEXT
    ZZMASS_LS(:,:,           JK) = XZMASS_SV_LS(:,:,1)   - (XZMASS_SV_LS(:,:,2)  -XZMASS_SV_LS(:,:,1)    )*(JPVEXT+1-JK)
    ZZMASS_LS(:,:,ILU+JPVEXT+JK) = XZMASS_SV_LS(:,:,ILU) + (XZMASS_SV_LS(:,:,ILU)-XZMASS_SV_LS(:,:,ILU-1))*          JK
  END DO
END IF
!
ZSV_LS    = XUNDEF
IF (HFILE(1:3)=='ATM') THEN
  ZSV_LS(:,:,:,:) = 0.
ELSE IF (HFILE=='CHEM') THEN
DO JSV=1,SIZE(XSV_LS,4)
  ZSV_LS    (:,:,JPVEXT+1:JPVEXT+ILU,JSV) = XSV_LS (:,:,:,JSV)
END DO
END IF
!
IF (HFILE(1:3)=='ATM') THEN
  ZZFLUX_LS (:,:,JPVEXT+1:JPVEXT+ILU) = XZFLUX_LS(:,:,:)
  DO JK=1,JPVEXT
    ZZFLUX_LS(:,:,           JK) = XZFLUX_LS(:,:,1)   - (XZFLUX_LS(:,:,2)  -XZFLUX_LS(:,:,1)    )*(JPVEXT+1-JK)
    ZZFLUX_LS(:,:,ILU+JPVEXT+JK) = XZFLUX_LS(:,:,ILU) + (XZFLUX_LS(:,:,ILU)-XZFLUX_LS(:,:,ILU-1))*          JK
  END DO
  !
  ZPMHP_LS  = XUNDEF
  ZTHV_LS   = XUNDEF
  ZR_LS     = XUNDEF
  ZHU_LS    = XUNDEF
  ZU_LS     = XUNDEF
  ZV_LS     = XUNDEF
  ZW_LS     = XUNDEF
  IF (SIZE(ZTKE_LS)>0) ZTKE_LS = XUNDEF
  !
  ZPMHP_LS  (:,:,JPVEXT+1:JPVEXT+ILU) = XPMHP_LS (:,:,:)
  ZTHV_LS   (:,:,JPVEXT+1:JPVEXT+ILU) = XTHV_LS  (:,:,:)
  ZHU_LS    (:,:,JPVEXT+1:JPVEXT+ILU) = XHU_LS   (:,:,:)
  ZU_LS     (:,:,JPVEXT+1:JPVEXT+ILU) = XU_LS    (:,:,:)
  ZV_LS     (:,:,JPVEXT+1:JPVEXT+ILU) = XV_LS    (:,:,:)
  ZW_LS     (:,:,JPVEXT+1:JPVEXT+ILU) = XW_LS    (:,:,:)
  IF (SIZE(ZTKE_LS)>0) ZTKE_LS(:,:,JPVEXT+1:JPVEXT+ILU) = XTKE_LS (:,:,:)
  DO JRR=1,NRR
    ZR_LS     (:,:,JPVEXT+1:JPVEXT+ILU,JRR) = XR_LS  (:,:,:,JRR)
  END DO
END IF
!
IF (HFILE(1:3)=='ATM') THEN

  IF (SIZE(XB_LS)/=0) THEN   ! hybrid level (w at flux points)
  CALL VER_INTERP_TO_MIXED_GRID('ATM ',.TRUE.,XZS_LS,XZSMT_LS,    &
                                ZZMASS_LS,ZSV_LS,                &
                                ZZFLUX_LS,XPS_LS,ZPMHP_LS,       &
                                ZTHV_LS,ZR_LS,                   &
                                ZHU_LS,                          &
                                ZTKE_LS,                         &
                                ZU_LS,ZV_LS,                     &
                                ZW_LS,'FLUX'                     )
  ELSE                      ! isobaric surfaces (w at mass points)
    CALL VER_INTERP_TO_MIXED_GRID('ATM ',.TRUE.,XZS_LS,XZSMT_LS,    &
                                  ZZMASS_LS,ZSV_LS,                &
                                  ZZFLUX_LS,XPS_LS,ZPMHP_LS,       &
                                  ZTHV_LS,ZR_LS,                   &
                                  ZHU_LS,                          &
                                  ZTKE_LS,                         &
                                  ZU_LS,ZV_LS,                     &
                                  ZW_LS,'MASS'                     )
  END IF
ELSE IF (HFILE=='CHEM') THEN
  CALL VER_INTERP_TO_MIXED_GRID(HFILE,.TRUE.,XZS_SV_LS,XZS_SV_LS,&
                                ZZMASS_LS,ZSV_LS                 )
END IF
!
DEALLOCATE(ZZMASS_LS)
DEALLOCATE(ZSV_LS)
IF (HFILE(1:3)=='ATM') THEN
  DEALLOCATE(ZZFLUX_LS)
  DEALLOCATE(ZPMHP_LS)
  DEALLOCATE(ZTHV_LS)
  DEALLOCATE(ZR_LS)
  DEALLOCATE(ZHU_LS)
  DEALLOCATE(ZU_LS)
  DEALLOCATE(ZV_LS)
  DEALLOCATE(ZW_LS)
  DEALLOCATE(ZTKE_LS)
  !
  DEALLOCATE(XPMHP_LS)
END IF
!-------------------------------------------------------------------------------
!
!*       3.    ERROR CONTROL
!              -------------
!
CALL SECOND_MNH(ZTIME1)
IF (HFILE(1:3)=='ATM' .AND. NVERB>=5) THEN
!
  ALLOCATE(ZTH_LS  (SIZE(XTHV_LS  ,1),SIZE(XTHV_LS  ,2),SIZE(XTHV_LS  ,3)))
  ALLOCATE(ZTH_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  ZTH_LS(:,:,:)  =XTHV_LS  (:,:,:)/(1.+XRV/XRD*XR_LS  (:,:,:,1))*(1.+WATER_SUM(XR_LS  (:,:,:,:)))
  ZTH_MX(:,:,:)=XTHV_MX(:,:,:)/(1.+XRV/XRD*XR_MX(:,:,:,1))*(1.+WATER_SUM(XR_MX(:,:,:,:)))
!
  CALL RMS_AT_Z(ZTH_LS,XZS_LS,XZMASS_LS,ZTH_MX,XZS_LS,XZMASS_MX, &
                'RMS on theta between GRIB file grid and mixed grid (K):                          ')
!
  ALLOCATE(ZPMASS_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  ALLOCATE(ZHEXNMASS_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  ALLOCATE(ZHEXNFLUX_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  CALL COMPUTE_EXNER_FROM_TOP(XTHV_MX,XZFLUX_MX,XEXNTOP2D,ZHEXNFLUX_MX,ZHEXNMASS_MX)
  ZPMASS_MX(:,:,:)=XP00*(ZHEXNMASS_MX(:,:,:))**(XCPD/XRD)
!
  CALL RMS_AT_Z(XPMASS_LS,XZS_LS,XZMASS_LS,ZPMASS_MX,XZS_LS,XZMASS_MX, &
                'RMS on pressure between GRIB file grid and mixed grid (Pa):                      ')
!
  DEALLOCATE(ZTH_LS)
  DEALLOCATE(ZTH_MX)
  DEALLOCATE(ZPMASS_MX)
  DEALLOCATE(ZHEXNFLUX_MX)
  DEALLOCATE(ZHEXNMASS_MX)
!
END IF
CALL SECOND_MNH(ZTIME2)
PDIAG = ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       4.    DEALLOCATIONS
!              -------------
!
IF (HFILE(1:3)=='ATM') THEN
  DEALLOCATE(XA_LS)
  DEALLOCATE(XB_LS)
  DEALLOCATE(XQ_LS)
  DEALLOCATE(XU_LS)
  DEALLOCATE(XV_LS)
  DEALLOCATE(XZMASS_LS)
  DEALLOCATE(XZFLUX_LS)
  DEALLOCATE(XTHV_LS)
  DEALLOCATE(XR_LS)
  DEALLOCATE(XHU_LS)
  IF (HFILE=='ATM0') DEALLOCATE(XSV_LS)
ELSE IF (HFILE=='CHEM') THEN
  DEALLOCATE(XA_SV_LS)
  DEALLOCATE(XB_SV_LS)
  DEALLOCATE(XT_SV_LS)
  DEALLOCATE(XQ_SV_LS)
  DEALLOCATE(XZMASS_SV_LS)
  DEALLOCATE(XZFLUX_SV_LS)
  DEALLOCATE(XTHV_SV_LS)
  DEALLOCATE(XR_SV_LS)
  DEALLOCATE(XHU_SV_LS)
  DEALLOCATE(XSV_LS)
END IF
!
!
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'Routine VER_PREP_GRIBEX_CASE for ',HFILE,' file completed'
!
END SUBROUTINE VER_PREP_GRIBEX_CASE
