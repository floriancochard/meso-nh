!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_VER_DYN
!     ###################
INTERFACE
      SUBROUTINE VER_DYN(OSHIFT,                                                &
                         PU_MX,PV_MX,PW_MX,PRHOD_MX,PZFLUX_MX,PZMASS_MX,PZS_LS, &
                         PDXX,PDYY,PDZZ,PDZX,PDZY,PJ,HATMFILETYPE,              &
                         PLSU_MX,PLSV_MX,PLSW_MX                                )
!
LOGICAL,                  INTENT(IN)  :: OSHIFT     ! T: vertical shift of BL (used for GRIB file data)
!                                                   ! F: no vertical shift (used for MESONH data)
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PU_MX     ! U on LS grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PV_MX     ! V on LS grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PW_MX     ! W on LS grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PRHOD_MX  ! local rhod on mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX_MX ! altitude of pressure
!                                                  ! points on mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! altitude of mass
!                                                  ! points on mixed grid
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! metric coefficient dyy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZZ      ! metric coefficient dzz
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZX      ! metric coefficient dzx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZY      ! metric coefficient dzy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! jacobian
CHARACTER(LEN=6)                  :: HATMFILETYPE  ! type of the Atmospheric file
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSU_MX ! large scale U component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSV_MX ! large scale V component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSW_MX ! large scale W component
!
END SUBROUTINE VER_DYN
END INTERFACE
END MODULE MODI_VER_DYN

!     ######spl
      SUBROUTINE VER_DYN(OSHIFT,                                               &
                         PU_MX,PV_MX,PW_MX,PRHOD_MX,PZFLUX_MX,PZMASS_MX,PZS_LS,&
                         PDXX,PDYY,PDZZ,PDZX,PDZY,PJ,HATMFILETYPE,             &
                         PLSU_MX,PLSV_MX,PLSW_MX                               )
!     ##########################################################################
!
!!****  *VER_DYN* -  initializes dynamical fields in MESO-NH for a real case
!!                   from Aladin fields.
!!
!!    PURPOSE
!!    -------
!!    This routine initializes the three components of the momentum
!!    on the MESO-NH Arakawa C-grid from the horizontal fields on the mixed
!!    Arakawa A-grid defined by the altitudes of its pressure and
!!    mass points.
!!
!!**  METHOD
!!    ------
!!
!!  1 The values of wind components are multiplied by the density to obtain
!!    the horizontal momentum components (ZRHODU_MX,ZRHODV_MX).
!!
!!  2 The first guess of horizontal momentum is initialized in VER_INT_DYN on
!!    the Arakawa A-grid (ZRHODUA,ZRHODVA).
!!
!!  3 The values on the Arakawa C-grid are deduced (ZRHODU,ZRHODV)
!!
!!  4 The first guess of vertical momentum is the large scale field, which is
!!    absurd near the ground.
!!    At the time being, there is a forcing to zero value in the upper quarter
!!    of the domain
!!
!!  5 The interpolated fields ZRHODJU, ZRHODJV are placed in the
!!    module MODD_FIELD1 in place of the prognostic variables XUM, XVM.
!!
!!
!!    EXTERNAL
!!    --------
!!    subroutine VER_INT_DYN    : to initialize the horizontal momentum
!!    subroutine WGUESS         : to initialize vertical momentum
!!    subroutine ANEL_BALANCE1  : to apply the anelastic correction
!!    functions MXM ,MYM ,MZM   : Shuman operators
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_FIELD1    : contains prognostics  variables
!!         XUM : U (:,:,:)         at t-dt
!!         XVM : V (:,:,:,:)       at t-dt
!!         XWM : w (:,:,:,:)       at t-dt
!!      Module MODD_LBC1      : contains lateral boundary conditions
!!         CLBCX   : X-direction LBC type at left(1) and right(2) boundaries
!!         CLBCY   : Y-direction LBC type at left(1) and right(2) boundaries
!!      Module MODD_REF1      : contains 3D reference state variables for model1
!!         XRHODJ  : rhod * Jacobian
!!         XTHVREF_STAR : THvref * (1 + Rvref)
!!      Module MODD_PARAMETERS:
!!         JPVEXT
!!      Module MODD_GRID1    :
!!         XZZ   : height of the w-points
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
!!      J.Stein and J.P. Lafore   18/04/96  change the anel_balance CALL
!!      J.Stein                   15/05/96  change the wguess CALL
!!      V.Masson                  11/10/96  L1D and L2D configurations
!!      V.Masson                  12/12/96  add LS vertical wind velocity
!!      Stein,Lafore              15/01/97  Durran anelastic equation
!!      V.Masson                  26/08/97  call to new linear vertical
!!                                          interpolation routine
!!      V.Masson                  24/11/97  use of the 3D dry density
!!      J.Stein                   20:01/98  add the LS field interpolation
!!      M.Faivre                  2014
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ARGSLIST_ll, ONLY: LIST_ll
USE MODD_CONF
USE MODD_CST
USE MODD_DIM_n
USE MODD_DYN_n
USE MODD_FIELD_n,     ONLY: XUT,XVT,XWT,XPABST,XTHT,XRT
USE MODD_GRID_n
USE MODD_LBC_n
USE MODD_LSFIELD_n
USE MODD_LUNIT,       ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_REF_n
USE MODD_VER_INTERP_LIN
!
USE MODE_EXTRAPOL
USE MODE_ll
USE MODE_MPPDB
!
USE MODI_ANEL_BALANCE_n
USE MODI_COEF_VER_INTERP_LIN
USE MODI_SHUMAN
USE MODI_VER_INT_DYN
USE MODI_VER_INTERP_LIN
USE MODI_VER_SHIFT
USE MODI_WGUESS
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
LOGICAL,                  INTENT(IN)  :: OSHIFT     ! T: vertical shift of BL (used for GRIB file data)
!                                                   ! F: no vertical shift (used for MESONH data)
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PU_MX     ! U on LS grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PV_MX     ! V on LS grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PW_MX     ! W on LS grid
REAL,   DIMENSION(:,:,:), INTENT(IN)     :: PRHOD_MX  ! local rhod on mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX_MX ! altitude of pressure
!                                                  ! points on mixed grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! altitude of mass
!                                                  ! points on mixed grid
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! metric coefficient dyy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZZ      ! metric coefficient dzz
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZX      ! metric coefficient dzx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZY      ! metric coefficient dzy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! jacobian
CHARACTER(LEN=6)                  :: HATMFILETYPE  ! type of the Atmospheric file
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSU_MX ! large scale U component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSV_MX ! large scale V component
REAL,DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PLSW_MX ! large scale W component
!
!*       0.2   Declaration of local variables
!              ------------------------------
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))    :: ZRHODU_MX, ZRHODV_MX
!                  ! horizontal momentum components on the mixed grid
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))    :: ZRHODUA, ZRHODVA
!                  ! horizontal momentum components on the MESONH Arakawa A grid
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))    :: ZRHODJU, ZRHODJV
!                  ! momentum components on the MESONH Arakawa C grid
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))    :: ZRHOD
!                  ! dry density on MESO-NH grid
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))    :: ZZFLUX_SH
                   ! shifted flux grid
REAL,DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3))    ::  ZCOEF
                  ! coefficient  for weight function
!
INTEGER :: IIB,IIE,IIU
INTEGER :: IJB,IJE,IJU
INTEGER :: IKB,IKE,IKU
INTEGER :: ILBX,ILBY
!
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
!
!20131105 declare vars related to add3dfield and update_halo_ll
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()  ! list of fields to exchange
!
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIU=SIZE(PJ,1)
IJU=SIZE(PJ,2)
IKB=JPVEXT+1
IKE=SIZE(PJ,3)-JPVEXT
IKU=SIZE(PJ,3)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF MOMENTUM ON THE MIXED GRID
!              -----------------------------------------
!
ZRHODU_MX=PU_MX*PRHOD_MX
ZRHODV_MX=PV_MX*PRHOD_MX
!
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZATION OF HORIZONTAL MOMENTUM ON MESO-NH ARAKAWA A-GRID
!              ---------------------------------------------------------------
!
CALL VER_INT_DYN(OSHIFT,ZRHODU_MX,ZRHODV_MX,PZFLUX_MX,PZMASS_MX,PZS_LS,ZRHODUA,ZRHODVA)
!
CALL EXTRAPOL('E',ZRHODUA,ZRHODVA)

CALL MPPDB_CHECK3D(ZRHODUA,"VERDYN::ZRHODUA",PRECISION)
CALL MPPDB_CHECK3D(ZRHODVA,"VERDYN::ZRHODVA",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       3.    CHANGE TO ARAKAWA C-GRID
!              ------------------------
!
!20131105 add UPDATE_HALO on ZRHODUA, ZRHODVA and PJ(needed) : ok ZRHODJU,V
!20131112 impact of PJ update_halo =>ZRHODJU,V error => XUM,XVM error
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODUA)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODVA)
CALL ADD3DFIELD_ll(TZFIELDS_ll,PJ)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!
ZRHODJU(:,:,:)=MXM(ZRHODUA(:,:,:)*PJ(:,:,:))
ZRHODJV(:,:,:)=MYM(ZRHODVA(:,:,:)*PJ(:,:,:))
!
!20131112 add update_halo_ll though checking on vars is correct
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODJU)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODJV)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!
!20131105 add check
CALL MPPDB_CHECK3D(ZRHODJU,"VERDYN3-before extrapol::ZRHODJU",PRECISION)
CALL MPPDB_CHECK3D(ZRHODJV,"VERDYN3-before extrapol::ZRHODJV",PRECISION)
!
CALL EXTRAPOL('W',ZRHODJU)
CALL EXTRAPOL('S',ZRHODJV)
!
!CALL MPPDB_CHECK3D(ZRHODJU,"VERDYN::ZRHODJU",PRECISION)
!CALL MPPDB_CHECK3D(ZRHODJV,"VERDYN::ZRHODJV",PRECISION)
!
!20131104 add check
CALL MPPDB_CHECK3D(ZRHODJU,"VERDYN3-after extrapol::ZRHODJU",PRECISION)
CALL MPPDB_CHECK3D(ZRHODJV,"VERDYN3-after extrapol::ZRHODJV",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       4.    STORAGE IN MODD_FIELD1
!              ----------------------
!
ALLOCATE(XUT(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
ALLOCATE(XVT(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
ALLOCATE(XWT(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
!
!20131104 add check on xpabsm
CALL MPPDB_CHECK3D(XPABST,"VER_DYN4::XPABST",PRECISION)
CALL MPPDB_CHECK3D(XTHT,"VER_DYN4::XTHT",PRECISION)
!
ZRHOD(:,:,:)=XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**(XRD/XCPD) &
            /(XRD*XTHT(:,:,:)*(1.+XRV/XRD*XRT(:,:,:,1)))
!
CALL MPPDB_CHECK3D(ZRHOD,"VER_DYN4::ZRHOD",PRECISION)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHOD)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!
XUT(:,:,:)=ZRHODJU(:,:,:)/MXM(ZRHOD(:,:,:)*PJ(:,:,:))
XVT(:,:,:)=ZRHODJV(:,:,:)/MYM(ZRHOD(:,:,:)*PJ(:,:,:))

CALL ADD3DFIELD_ll(TZFIELDS_ll,XUT)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XVT)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
CALL EXTRAPOL('W',XUT)
CALL EXTRAPOL('S',XVT)

CALL MPPDB_CHECK3D(XUT,"VER_DYN4-after extrapol::XUT",PRECISION)
CALL MPPDB_CHECK3D(XVT,"VER_DYN4-after extrapol::XVT",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       5.    INITIALIZATION OF LS HORIZONTAL MOMENTUM ON MESO-NH ARAKAWA A-GRID
!              ------------------------------------------------------------------
!
IF( HATMFILETYPE == 'MESONH' ) THEN
  !
  ZRHODU_MX=PLSU_MX*PRHOD_MX
  ZRHODV_MX=PLSV_MX*PRHOD_MX
  CALL VER_INT_DYN(OSHIFT,ZRHODU_MX,ZRHODV_MX,PZFLUX_MX,PZMASS_MX,PZS_LS,ZRHODUA,ZRHODVA)
  !
  ALLOCATE(XLSUM(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
  ALLOCATE(XLSVM(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
  !
  !20131104 add check on zrhodua, zrhodju
  CALL MPPDB_CHECK3D(ZRHODUA,"VER_DYN5::ZRHODUA",PRECISION)
  CALL MPPDB_CHECK3D(ZRHODJU,"VER_DYN5::ZRHODJU",PRECISION)
  !
  !20131105 add UPDATE_HALO on ZRHODUA, ZRHODVA and PJ(needed): ok ZRHODJU,V
  CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODUA)
  CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODVA)
  CALL ADD3DFIELD_ll(TZFIELDS_ll,PJ)
     CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
        CALL CLEANLIST_ll(TZFIELDS_ll)  
  !
  ZRHODJU(:,:,:)=MXM(ZRHODUA(:,:,:)*PJ(:,:,:))
  ZRHODJV(:,:,:)=MYM(ZRHODVA(:,:,:)*PJ(:,:,:))
  !
  !20131112 add update_halo_ll though checking on vars is correct
  CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODJU)
  CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHODJV)
     CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
        CALL CLEANLIST_ll(TZFIELDS_ll)
  !
  !20131105 add UPDATE_HALO on ZRHOD
    CALL ADD3DFIELD_ll(TZFIELDS_ll,ZRHOD)
       CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
          CALL CLEANLIST_ll(TZFIELDS_ll)
  !
  XLSUM(:,:,:)=ZRHODJU(:,:,:)/MXM(ZRHOD(:,:,:)*PJ(:,:,:))
  XLSVM(:,:,:)=ZRHODJV(:,:,:)/MYM(ZRHOD(:,:,:)*PJ(:,:,:))
  !
  !20131112 add update_halo_ll though checking on vars is correct
  CALL ADD3DFIELD_ll(TZFIELDS_ll,XLSUM)
  CALL ADD3DFIELD_ll(TZFIELDS_ll,XLSVM)
     CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
        CALL CLEANLIST_ll(TZFIELDS_ll)
  !20131126 add check on XLSUN, XLSVM
  CALL MPPDB_CHECK3D(XLSUM,"VER_DYN5-beforeextrapol::XLSUM",PRECISION)
  CALL MPPDB_CHECK3D(XLSVM,"VER_DYN5-beforeextrapol::XLSVM",PRECISION)
  !
END IF
!
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTATION OF FIRST GUESS OF W
!              -------------------------------
!
ZZFLUX_SH(:,:,:)=VER_SHIFT(PZFLUX_MX,PZS_LS,XZS)
CALL COEF_VER_INTERP_LIN(ZZFLUX_SH(:,:,:),XZZ(:,:,:))
XWT(:,:,:)=VER_INTERP_LIN(PW_MX(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
!
IF ( HATMFILETYPE == 'MESONH' ) THEN
  ALLOCATE(XLSWM(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)))
  XLSWM(:,:,:)=VER_INTERP_LIN(PLSW_MX(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
END IF
!
!20131126 add check on XWM,XLSWM
CALL MPPDB_CHECK3D(XWT,"VER_DYN5::XWT",PRECISION)
! CALL MPPDB_CHECK3D(XLSWM,"VER_DYN5::XLSWM",PRECISION)
!
DEALLOCATE(NKLIN)
DEALLOCATE(XCOEFLIN)
!
!*       5.2   forcing to zero value at top (to be removed when solver allows other values)
!
ZCOEF(:,:,:)=(       XZZ(:,:,:)           -SPREAD(XZZ(:,:,IKB),3,IKU)) &
            /(SPREAD(XZZ(:,:,IKE+1),3,IKU)-SPREAD(XZZ(:,:,IKB),3,IKU))
XWT(:,:,:)=XWT(:,:,:)*MAX(MIN( (4.-4.*ZCOEF(:,:,:)) ,1.),0.)
!-------------------------------------------------------------------------------
!
!20131112 add update_halo_ll though checking on vars is correct
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZCOEF)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XWT)
   CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
      CALL CLEANLIST_ll(TZFIELDS_ll)
!
!------------------------------------------------------------------------------
!
!*       6.    STORAGE OF LARGE SCALE FIELDS
!              ------------------------------
!
IF ( HATMFILETYPE == 'GRIBEX' ) THEN
  ALLOCATE(XLSUM(SIZE(XUT,1),SIZE(XUT,2),SIZE(XUT,3)))
  ALLOCATE(XLSVM(SIZE(XVT,1),SIZE(XVT,2),SIZE(XVT,3)))
  ALLOCATE(XLSWM(SIZE(XWT,1),SIZE(XWT,2),SIZE(XWT,3)))
  XLSUM(:,:,:)=XUT(:,:,:)
  XLSVM(:,:,:)=XVT(:,:,:)
  XLSWM(:,:,:)=XWT(:,:,:)
END IF
! enforce zero gradient along the vertical under and above the vertical
! boundaries
XLSUM(:,:,IKB-1)=XLSUM(:,:,IKB)
XLSUM(:,:,IKE+1)=XLSUM(:,:,IKE)
XLSVM(:,:,IKB-1)=XLSVM(:,:,IKB)
XLSVM(:,:,IKE+1)=XLSVM(:,:,IKE)
XLSWM(:,:,IKB-1)=XLSWM(:,:,IKB)
XLSWM(:,:,IKE+1)=XLSWM(:,:,IKE)

CALL EXTRAPOL('W',XLSUM)
CALL EXTRAPOL('E',XLSUM)
CALL EXTRAPOL('S',XLSVM)
CALL EXTRAPOL('E',XLSVM)
!
!20131126 add check on XLSUN, XLSVM
CALL MPPDB_CHECK3D(XLSUM,"VER_DYN5-afterextrapol::XLSUM",PRECISION)
CALL MPPDB_CHECK3D(XLSVM,"VER_DYN5-afterextrapol::XLSVM",PRECISION)
!
! FROM PREP_IDEAL_CASE
!
! 3D case
!
  CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,               &
                    IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
                    IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  CALL GET_SIZEY_LB(NIMAX_ll,NJMAX_ll,NRIMY,               &
                    IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV, &
                    IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)

  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2*JPHEXT
    NSIZELBXU_ll=2*NRIMX+2*JPHEXT
    NSIZELBY_ll=2*NRIMY+2*JPHEXT
    NSIZELBYV_ll=2*NRIMY+2*JPHEXT
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,IKU))
    ALLOCATE(XLBYUM(IISIZEYF,IJSIZEYF,IKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYVM(IISIZEYFV,IJSIZEYFV,IKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYWM(IISIZEYF,IJSIZEYF,IKU))
    !ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,IKU))
    !ALLOCATE(XLBYTHM(IISIZEYF,IJSIZEYF,IKU))
  ELSE
    NSIZELBX_ll= 2*JPHEXT  ! 2
    NSIZELBXU_ll= 2*(JPHEXT+1) ! 4
    NSIZELBY_ll=2*JPHEXT  !  2
    NSIZELBYV_ll= 2*(JPHEXT+1) ! 4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,IKU))
    ALLOCATE(XLBYUM(IISIZEY2,IJSIZEY2,IKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYVM(IISIZEY4,IJSIZEY4,IKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYWM(IISIZEY2,IJSIZEY2,IKU))
    !ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,IKU))
    !ALLOCATE(XLBYTHM(IISIZEY2,IJSIZEY2,IKU))
  END IF  

ILBX=SIZE(XLBXUM,1)
ILBY=SIZE(XLBYUM,2)
IF(LWEST_ll() .AND. .NOT. L1D) THEN
  XLBXUM(1:NRIMX+JPHEXT,        :,:)     = XUT(2:NRIMX+JPHEXT+1,        :,:)
  XLBXVM(1:NRIMX+JPHEXT,        :,:)     = XVT(1:NRIMX+JPHEXT,        :,:)
  XLBXWM(1:NRIMX+JPHEXT,        :,:)     = XWT(1:NRIMX+JPHEXT,        :,:)

ENDIF
IF(LEAST_ll() .AND. .NOT. L1D) THEN
  XLBXUM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)     = XUT(IIU-NRIMX-JPHEXT+1:IIU,    :,:)
  XLBXVM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)     = XVT(IIU-NRIMX-JPHEXT+1:IIU,    :,:)
  XLBXWM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)     = XWT(IIU-NRIMX-JPHEXT+1:IIU,    :,:)

ENDIF
IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYUM(:,1:NRIMY+JPHEXT,        :)     = XUT(:,1:NRIMY+JPHEXT,      :)
  XLBYVM(:,1:NRIMY+JPHEXT,        :)     = XVT(:,2:NRIMY+JPHEXT+1,      :)
  XLBYWM(:,1:NRIMY+JPHEXT,        :)     = XWT(:,1:NRIMY+JPHEXT,  :)

ENDIF
IF(LNORTH_ll().AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYUM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)     = XUT(:,IJU-NRIMY-JPHEXT+1:IJU,  :)
  XLBYVM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)     = XVT(:,IJU-NRIMY-JPHEXT+1:IJU,  :)
  XLBYWM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)     = XWT(:,IJU-NRIMY-JPHEXT+1:IJU,  :)

ENDIF

!

CALL  MPPDB_CHECKLB(XLBXUM,"ver_dyn::XLBXUM::",PRECISION,'LBXU',NRIMX)
CALL  MPPDB_CHECKLB(XLBXVM,"ver_dyn::XLBXVM::",PRECISION,'LBXU',NRIMX)
CALL  MPPDB_CHECKLB(XLBXWM,"ver_dyn::XLBXWM::",PRECISION,'LBXU',NRIMX)


CALL  MPPDB_CHECKLB(XLBYUM,"ver_dyn::XLBYUM::",PRECISION,'LBYV',NRIMY)
CALL  MPPDB_CHECKLB(XLBYVM,"ver_dyn::XLBYVM::",PRECISION,'LBYV',NRIMY)
CALL  MPPDB_CHECKLB(XLBYWM,"ver_dyn::XLBYWM::",PRECISION,'LBYV',NRIMY)
!
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT0%NLU,*) 'Routine VER_DYN completed.'
!
END SUBROUTINE VER_DYN
