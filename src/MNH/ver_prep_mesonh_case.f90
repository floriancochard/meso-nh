!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################################
      MODULE MODI_VER_PREP_MESONH_CASE
!     ################################
INTERFACE
      SUBROUTINE VER_PREP_MESONH_CASE(PDIAG)
!
REAL, INTENT(OUT)                 :: PDIAG    ! diagnostics computing time
!
END SUBROUTINE VER_PREP_MESONH_CASE
END INTERFACE
END MODULE MODI_VER_PREP_MESONH_CASE
!     ####################################################################
      SUBROUTINE VER_PREP_MESONH_CASE(PDIAG)
!     ####################################################################
!
!!****  *VER_PREP_MESONH_CASE* - monitors the preparation to orographic change
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
!!
!!    module MODI_SHUMAN
!!    module MODI_VER_INTERP_TO_MIXED_GRID
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB      : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         TLUOUT0 : output-listing file
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
!!                  Aug, 20 1996 (V. Masson) correction in virtual temperature
!!                               computation
!!                  Oct, 20 1996 (V. Masson) add deallocations
!!                  Dec, 06 1996 (V. Masson) add air temperature at ground
!!                  Dec, 12 1996 (V. Masson) add vertical velocity
!!                  May, 07 1997 (V. Masson) add tke
!!                  Jun, 10 1997 (V. Masson) add non-hydrostatic pressure
!!                  Jul, 10 1997 (V. Masson) add epsilon
!!                  Jul, 11 1997 (V. Masson) add scalar variables
!!                          2014 (M.Faivre)  parallelization
!!                  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ARGSLIST_ll, ONLY: LIST_ll
USE MODD_CONF
USE MODD_CST
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_PREP_REAL
!
USE MODE_ll
USE MODE_MPPDB
USE MODE_THERMO
!
USE MODI_COMPUTE_EXNER_FROM_GROUND
USE MODI_COMPUTE_EXNER_FROM_TOP
USE MODI_RMS_AT_Z
USE MODI_SECOND_MNH
USE MODI_SHUMAN
USE MODI_VER_INTERP_TO_MIXED_GRID
USE MODI_VERT_COORD
USE MODI_WATER_SUM
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL, INTENT(OUT)                 :: PDIAG    ! diagnostics computing time
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER                            :: IIU,IJU,ILU
REAL                               :: ZTIME1, ZTIME2
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHPMASS_LS  ! hyd. pressure
REAL,DIMENSION(:,:),   ALLOCATABLE :: ZHEXNSURF2D ! surface hyd. pressure function
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHEXNFLUX_LS! hyd. pressure function
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHEXNMASS_LS! hyd. pressure function
!
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZES_LS      ! vapor saturation pressure
!
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZTH_MX      ! potential temperature
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHEXNFLUX_MX! hyd. pressure function
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZHEXNMASS_MX! hyd. pressure function
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZPMASS_MX   ! pressure
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZWORK       ! work array
!20131105 add vars related to ADD3DFIELD and UPDATE_HALO
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
!
INTEGER                           :: IIB,IJB,IIE,IJE
!-------------------------------------------------------------------------------
!
!*       1.    CHANGING OF VARIABLES
!              ---------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
IIU=SIZE(XZS_LS,1)
IJU=SIZE(XZS_LS,2)
ILU=SIZE(XZHAT_LS)
!
ALLOCATE(XTHV_LS(IIU,IJU,ILU))
ALLOCATE(XZFLUX_LS(IIU,IJU,ILU))
ALLOCATE(XZMASS_LS(IIU,IJU,ILU))
!
!*       1.1   Virtual potential temperature
!              -----------------------------
!
XTHV_LS(:,:,:)=XTH_LS(:,:,:)*(1.+XRV/XRD*XR_LS(:,:,:,1))/(1.+WATER_SUM(XR_LS(:,:,:,:)))
!
!20131112 add update_halo for this type of calculation
CALL ADD3DFIELD_ll(TZFIELDS_ll, XTHV_LS)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(XTHV_LS,"ver_prep_mesonh_case1.1::XTHV_LS",PRECISION)
!
!*       1.2   Altitudes
!              ---------
!
CALL VERT_COORD(LSLEVE_LS,XZS_LS,XZSMT_LS,XLEN1_LS,XLEN2_LS,XZHAT_LS,XZFLUX_LS)
!
!20140217 upgrade MZF
!$XZMASS_LS(:,:,:)=MZF(XZFLUX_LS(:,:,:))
XZMASS_LS(:,:,:)=MZF(1,ILU,1,XZFLUX_LS(:,:,:))
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(XZMASS_LS,"ver_prep_mesonh_case1.2a::XZMASS_LS",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll, XZMASS_LS)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(XZMASS_LS,"ver_prep_mesonh_case1.2b::XZMASS_LS",PRECISION)
!
XZMASS_LS(:,:,SIZE(XZMASS_LS,3))=1.5*XZFLUX_LS(:,:,SIZE(XZFLUX_LS,3)  )    &
                                -0.5*XZFLUX_LS(:,:,SIZE(XZFLUX_LS,3)-1)
!
!*       1.3   Winds on Arakawa A grid
!              -----------------------
!
!20131104 check XU_LS calculated in read_prc_fmfile before any changes in
!ver_prep_mesonh_case
CALL MPPDB_CHECK3D(XU_LS,"ver_prep_mesonh_case1.3::XU_LS",PRECISION)
!
ALLOCATE(ZWORK(IIU,IJU,ILU))
ZWORK = XU_LS
XU_LS(1:IIU-1,:,:)=0.5*ZWORK(2:IIU,:,:)+0.5*ZWORK(1:IIU-1,:,:)
XU_LS(IIU    ,:,:)=1.5*ZWORK(IIU  ,:,:)-0.5*ZWORK(IIU-1  ,:,:)
XU_LS(IIE+1    ,:,:)=1.5*ZWORK(IIE+1  ,:,:)-0.5*ZWORK(IIE  ,:,:) ! for JPHEXT <> 1 

ZWORK = XV_LS
XV_LS(:,1:IJU-1,:)=0.5*ZWORK(:,2:IJU,:)+0.5*ZWORK(:,1:IJU-1,:)
XV_LS(:,IJU    ,:)=1.5*ZWORK(:,IJU  ,:)-0.5*ZWORK(:,IJU-1  ,:)
XV_LS(:,IJE+1    ,:)=1.5*ZWORK(:,IJE+1  ,:)-0.5*ZWORK(:,IJE  ,:) ! for JPHEXT <> 1 

ZWORK = XLSU_LS
XLSU_LS(1:IIU-1,:,:)=0.5*ZWORK(2:IIU,:,:)+0.5*ZWORK(1:IIU-1,:,:)
XLSU_LS(IIU    ,:,:)=1.5*ZWORK(IIU  ,:,:)-0.5*ZWORK(IIU-1  ,:,:)
XLSU_LS(IIE+1    ,:,:)=1.5*ZWORK(IIE+1  ,:,:)-0.5*ZWORK(IIE  ,:,:) ! for JPHEXT <> 1 

ZWORK = XLSV_LS
XLSV_LS(:,1:IJU-1,:)=0.5*ZWORK(:,2:IJU,:)+0.5*ZWORK(:,1:IJU-1,:)
XLSV_LS(:,IJU    ,:)=1.5*ZWORK(:,IJU  ,:)-0.5*ZWORK(:,IJU-1  ,:)
XLSV_LS(:,IJE+1    ,:)=1.5*ZWORK(:,IJE+1  ,:)-0.5*ZWORK(:,IJE  ,:) ! for JPHEXT <> 1 

DEALLOCATE(ZWORK)
!
!*       1.5   Difference between pressure and hydrostatic pressure
!              ----------------------------------------------------
!
ALLOCATE(ZHEXNSURF2D(IIU,IJU))
ALLOCATE(ZHPMASS_LS(IIU,IJU,ILU))
ALLOCATE(ZHEXNMASS_LS(IIU,IJU,ILU))
ALLOCATE(ZHEXNFLUX_LS(IIU,IJU,ILU))
ALLOCATE(XPMHP_LS(IIU,IJU,ILU))
ZHEXNSURF2D(:,:)=(XPS_LS(:,:)/XP00)**(XRD/XCPD)
!
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK2D(ZHEXNSURF2D,"ver_prep_mesonh_case1.5::ZHEXNSURF2D",PRECISION)
CALL ADD2DFIELD_ll(TZFIELDS_ll, ZHEXNSURF2D)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check2d
CALL MPPDB_CHECK2D(ZHEXNSURF2D,"ver_prep_mesonh_case1.5::ZHEXNSURF2D",PRECISION)
!
CALL COMPUTE_EXNER_FROM_GROUND(XTHV_LS,XZFLUX_LS,ZHEXNSURF2D,ZHEXNFLUX_LS,ZHEXNMASS_LS)
ZHPMASS_LS(:,:,:)=XP00*(ZHEXNMASS_LS(:,:,:))**(XCPD/XRD)
!
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(ZHPMASS_LS,"ver_prep_mesonh_case1.5a::ZHPMASS_LS",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZHPMASS_LS)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(ZHPMASS_LS,"ver_prep_mesonh_case1.5b::ZHPMASS_LS",PRECISION)
!
XPMHP_LS(:,:,:)=XPMASS_LS(:,:,:)-ZHPMASS_LS(:,:,:)
!
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(XPMHP_LS,"ver_prep_mesonh_case1.5c::XPMHP_LS",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XPMHP_LS)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(XPMHP_LS,"ver_prep_mesonh_case1.5d::XPMHP_LS",PRECISION)
!
DEALLOCATE(ZHPMASS_LS)
DEALLOCATE(ZHEXNSURF2D)
DEALLOCATE(ZHEXNFLUX_LS)
DEALLOCATE(ZHEXNMASS_LS)
!
!*       1.6   Relative humidity
!              -----------------
!
ALLOCATE(ZES_LS(IIU,IJU,ILU))
ALLOCATE(XHU_LS(IIU,IJU,ILU))
!
ZES_LS(:,:,:)=SM_FOES( XTHV_LS(:,:,:)                 &
                      *(1.+WATER_SUM(XR_LS(:,:,:,:))) &
                      /(1.+XRV/XRD*XR_LS(:,:,:,1))    &
                      *(XPMASS_LS(:,:,:)/XP00)**(XRD/XCPD)  )
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(ZES_LS,"ver_prep_mesonh_case1.6a::ZES_LS",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZES_LS)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(ZES_LS,"ver_prep_mesonh_case1.6b::ZES_LS",PRECISION)
!
XHU_LS(:,:,:)=100.*XPMASS_LS(:,:,:)/(XRD/XRV/MAX(XR_LS(:,:,:,1),1.E-12)+1.)/ZES_LS(:,:,:)
!
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(XHU_LS,"ver_prep_mesonh_case1.6c::XHU_LS",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XHU_LS)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(XHU_LS,"ver_prep_mesonh_case1.6d::XHU_LS",PRECISION)
!
DEALLOCATE(ZES_LS)
!-------------------------------------------------------------------------------
!
!*       2.    INTERPOLATION TO MIXED GRID AND DIAGNOSTIC VARIABLES
!              ----------------------------------------------------
!
!20131104 add check on XU_LS : IN arg of ver_interp_to_mixed in calculation of
!XU_MX
CALL MPPDB_CHECK3D(XU_LS,"ver_prep_mesonh_case2::XU_LS",PRECISION)
!
CALL VER_INTERP_TO_MIXED_GRID('ATM ',.FALSE.,XZS_LS,XZSMT_LS,                 &
                        XZMASS_LS,XSV_LS,XZFLUX_LS,XPS_LS,XPMHP_LS,           &
                        XTHV_LS,XR_LS,XHU_LS,XTKE_LS,                         &
                        XU_LS,XV_LS,XW_LS,'FLUX',                             &
                        XLSU_LS,XLSV_LS,XLSW_LS,XLSTH_LS,XLSRV_LS             )
!
!-------------------------------------------------------------------------------
!
!*       3.    ERROR CONTROL
!              -------------
!
CALL SECOND_MNH(ZTIME1)
IF (NVERB>=5) THEN
!
  ALLOCATE(ZTH_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  ZTH_MX(:,:,:)=XTHV_MX(:,:,:)/(1.+XRV/XRD*XR_MX(:,:,:,1))*(1.+WATER_SUM(XR_MX(:,:,:,:)))
!
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(ZTH_MX,"ver_prep_mesonh_case3a::ZTH_MX",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZTH_MX)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(ZTH_MX,"ver_prep_mesonh_case3b::ZTH_MX",PRECISION)
!
  CALL RMS_AT_Z(XTH_LS,XZS_LS,XZMASS_LS,ZTH_MX,XZS_LS,XZMASS_MX, &
                'RMS on theta between input Mesonh grid and mixed grid (K):                      ')
!
  ALLOCATE(ZPMASS_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  ALLOCATE(ZHEXNMASS_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  ALLOCATE(ZHEXNFLUX_MX(SIZE(XTHV_MX,1),SIZE(XTHV_MX,2),SIZE(XTHV_MX,3)))
  CALL COMPUTE_EXNER_FROM_TOP(XTHV_MX,XZFLUX_MX,XEXNTOP2D,ZHEXNFLUX_MX,ZHEXNMASS_MX)
  ZPMASS_MX(:,:,:)=XP00*(ZHEXNMASS_MX(:,:,:))**(XCPD/XRD) + XPMHP_MX(:,:,:)
!
!20131112 add update_halo for this type of calculation
CALL MPPDB_CHECK3D(ZPMASS_MX,"ver_prep_mesonh_case3c::ZPMASS_MX",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZPMASS_MX)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!20131112 check3d
CALL MPPDB_CHECK3D(ZPMASS_MX,"ver_prep_mesonh_case3d::ZPMASS_MX",PRECISION)
!
  CALL RMS_AT_Z(XPMASS_LS,XZS_LS,XZMASS_LS,ZPMASS_MX,XZS_LS,XZMASS_MX, &
                'RMS on pressure between input Mesonh grid and mixed grid (Pa):                  ')
!
  DEALLOCATE(ZTH_MX)
  DEALLOCATE(ZPMASS_MX)
  DEALLOCATE(ZHEXNFLUX_MX)
  DEALLOCATE(ZHEXNMASS_MX)
END IF
CALL SECOND_MNH(ZTIME2)
PDIAG = ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       4.    DEALLOCATIONS
!              -------------
!
DEALLOCATE(XPS_LS)
DEALLOCATE(XTH_LS)
DEALLOCATE(XR_LS)
DEALLOCATE(XTKE_LS)
DEALLOCATE(XU_LS)
DEALLOCATE(XV_LS)
DEALLOCATE(XW_LS)
DEALLOCATE(XPMASS_LS)
DEALLOCATE(XZHAT_LS)
DEALLOCATE(XLSU_LS)
DEALLOCATE(XLSV_LS)
DEALLOCATE(XLSW_LS)
DEALLOCATE(XLSTH_LS)
DEALLOCATE(XLSRV_LS)
!
DEALLOCATE(XZFLUX_LS)
DEALLOCATE(XZMASS_LS)
DEALLOCATE(XPMHP_LS)
DEALLOCATE(XHU_LS)
!
DEALLOCATE(XTHV_LS)
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT0%NLU,*) 'Routine VER_PREP_MESONH_CASE completed'
!
END SUBROUTINE VER_PREP_MESONH_CASE
