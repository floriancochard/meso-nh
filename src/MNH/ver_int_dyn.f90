!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_VER_INT_DYN
!     #######################
INTERFACE
      SUBROUTINE VER_INT_DYN(OSHIFT,PRHODU_MX,PRHODV_MX,PZFLUX_MX,PZMASS_MX,PZS_LS,&
                             PRHODUA,PRHODVA)
!
LOGICAL,                  INTENT(IN)  :: OSHIFT     ! T: vertical shift of BL (used for GRIB file data)
!                                                   ! F: no vertical shift (used for MESONH data)
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODU_MX ! rhoU on the mixed A-grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODV_MX ! rhoV on the mixed A-grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX_MX ! altitude of the pressure
!                                                  ! points of the mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! altitude of the mass points
!                                                  ! of the mixed grid
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PRHODUA   ! rhoU on MESO-NH A-grid
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PRHODVA   ! rhoV on MESO-NH A-grid
!
END SUBROUTINE VER_INT_DYN
END INTERFACE
END MODULE MODI_VER_INT_DYN
!     ######spl
      SUBROUTINE VER_INT_DYN(OSHIFT,PRHODU_MX,PRHODV_MX,PZFLUX_MX,PZMASS_MX,PZS_LS,&
                             PRHODUA,PRHODVA)
!     ###########################################################################
!
!!****  *VER_INT_DYN* - Vertical shift and interpolation of rhoU and rhoV.
!!
!!    PURPOSE
!!    -------
!!    This routine computes the 3D fields of horizontal momentum components
!!    on the MESO-NH grid with the MESO-NH orography
!!    from the corresponding fields on the mixed grid with large scale
!!    orography. The fields are kept on an Arakawa A-grid.
!!
!!**  METHOD
!!    ------
!!  * The change of orography is performed by a altitude shifting method:
!!      - the shifting of the grid PZMASS_MX occurs in VER_SHIFT --> ZZMASS_SH
!!      - the shifting of the grid PZFLUX_MX occurs in VER_SHIFT --> ZZFLUX_SH
!!      - the values on the shifted grid ZZMASS_SH are computed as:
!!
!!   PRHODU_MX(PZMASS_MX(k))*(PZFLUX_MX(k+1)-PZFLUX_MX(k))/(ZZFLUX_SH(k+1)-ZZFLUX_SH(k))
!!   PRHODV_MX(PZMASS_MX(k))*(PZFLUX_MX(k+1)-PZFLUX_MX(k))/(ZZFLUX_SH(k+1)-ZZFLUX_SH(k))
!!
!!  * The interpolation from the hybrid shifted grid to the MESO-NH grid is
!!    performed by the function VER_INTERP. Only the inner points are passed
!!    as arguments in this function (particularly along the vertical direction).
!!
!!    EXTERNAL
!!    --------
!!
!!    function VER_SHIFT    : to shift a array of altitudes
!!    function VER_INTERP   : to interpolate one field from one grid to another
!!    MZF                   : Shuman operator
!!
!!    Module MODI_SHUMAN    : interface for Shuman operators
!!    module MODI_VER_SHIFT : interface module for function VER_SHIFT
!!    module MODI_VER_INTERP: interface module for function VER_INTERP
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1     : contains grid variables for model1
!!         XZS   : orography of MESO-NH
!!         XZZ   : altitude of the w points in the MESO-NH grid.
!!      Module MODD_PARAMETERS
!!         JPVEXT
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
!!      Original    15/12/94
!!                  26/08/97 (V. Masson) call to new linear vertical
!!                                       interpolation routine
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_GRID_n
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_VER_INTERP_LIN
!
USE MODI_COEF_VER_INTERP_LIN
USE MODI_SHUMAN
USE MODI_VER_INTERP_LIN
USE MODI_VER_SHIFT
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
LOGICAL,                  INTENT(IN)  :: OSHIFT     ! T: vertical shift of BL (used for GRIB file data)
!                                                   ! F: no vertical shift (used for MESONH data)
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODU_MX ! rhodU on the mixed A-grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODV_MX ! rhodV on the mixed A-grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX_MX ! altitude of the pressure
!                                                  ! points of the mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! altitude of the mass points
!                                                  ! of the mixed grid
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PRHODUA   ! rhodU on MESO-NH A-grid
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PRHODVA   ! rhodV on MESO-NH A-grid
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER                    ::IKU
!
REAL,DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2),SIZE(PZMASS_MX,3))  &
                           ::ZZMASS_SH      ! altitude of the mass points
!                                           ! of the shifted grid
REAL,DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2),SIZE(PZMASS_MX,3))  &
                           ::ZZFLUX_SH      ! altitude of the pressure points
!                                           ! of the shifted grid
REAL,DIMENSION(SIZE(PRHODU_MX,1),SIZE(PRHODU_MX,2),SIZE(PRHODU_MX,3))  &
                           ::ZRHODU_SH,&   ! rhodU on the shifted grid
                             ZRHODV_SH     ! rhodV on the shifted grid
REAL,DIMENSION(SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3))                          &
                           ::ZZMASS        ! altitude of the mass points in
!                                          ! the MESO-NH grid.
!-------------------------------------------------------------------------------
!
!
IKU=SIZE(XZZ,3)
!
!-------------------------------------------------------------------------------
!
!*       1.    SHIFT OF THE HYBRID GRID
!              ------------------------
!
IF (OSHIFT) THEN
  ZZMASS_SH(:,:,:)=VER_SHIFT(PZMASS_MX,PZS_LS,XZS)
  ZZFLUX_SH(:,:,:)=VER_SHIFT(PZFLUX_MX,PZS_LS,XZS)
ELSE
  ZZFLUX_SH(:,:,:)=PZFLUX_MX
  ZZMASS_SH(:,:,:)=PZMASS_MX
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    SHIFT OF RHODU AND RHODV
!              ------------------------
!
ZRHODU_SH(:,:,1:IKU-1) = PRHODU_MX(:,:,1:IKU-1)               &
                 * (PZFLUX_MX(:,:,2:IKU)-PZFLUX_MX(:,:,1:IKU-1))      &
                 / (ZZFLUX_SH(:,:,2:IKU)-ZZFLUX_SH(:,:,1:IKU-1))
ZRHODU_SH(:,:,IKU)=2.*ZRHODU_SH(:,:,IKU-1)-ZRHODU_SH(:,:,IKU-2)
!
ZRHODV_SH(:,:,1:IKU-1) = PRHODV_MX(:,:,1:IKU-1)               &
                 * (PZFLUX_MX(:,:,2:IKU)-PZFLUX_MX(:,:,1:IKU-1))      &
                 / (ZZFLUX_SH(:,:,2:IKU)-ZZFLUX_SH(:,:,1:IKU-1))
ZRHODV_SH(:,:,IKU)=2.*ZRHODV_SH(:,:,IKU-1)-ZRHODV_SH(:,:,IKU-2)
!
!* no extrapolation below. Same value set to non-physical point
ZRHODU_SH(:,:,1) = ZRHODU_SH(:,:,2)
ZRHODV_SH(:,:,1) = ZRHODV_SH(:,:,2)
!-------------------------------------------------------------------------------
!
!*       3.    INTERPOLATION OF RHODU AND RHODV ON THE MESO-NH A-GRID
!              ------------------------------------------------------
!
!*       3.1   Altitude of the mass points on the MESO-NH grid
!              -----------------------------------------------
!
ZZMASS(:,:,:)=MZF(1,IKU,1,XZZ(:,:,:))
ZZMASS(:,:,SIZE(XZZ,3))=1.5*XZZ(:,:,SIZE(XZZ,3))-0.5*XZZ(:,:,SIZE(XZZ,3)-1)
!
!*       3.2   Interpolation on the MESO-NH grid
!              ---------------------------------
!
CALL COEF_VER_INTERP_LIN(ZZMASS_SH(:,:,:),ZZMASS(:,:,:))
PRHODUA(:,:,:)=VER_INTERP_LIN(ZRHODU_SH(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
PRHODVA(:,:,:)=VER_INTERP_LIN(ZRHODV_SH(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
!
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT0%NLU,*) 'Routine VER_INT_DYN completed'
!
END SUBROUTINE VER_INT_DYN
