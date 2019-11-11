!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_ANEL_BALANCE_n
!     ##########################
!
INTERFACE
!
SUBROUTINE ANEL_BALANCE_n(PRESIDUAL)
!
REAL, OPTIONAL                 :: PRESIDUAL
END SUBROUTINE ANEL_BALANCE_n
!
END INTERFACE
!
END MODULE MODI_ANEL_BALANCE_n
!
!
!
!     ################################
      SUBROUTINE ANEL_BALANCE_n(PRESIDUAL)
!
!     ################################
!
!
!!****  *ANEL_BALANCE_n* - routine to apply an anelastic correction
!!                   
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to  fulfill the anelastic balance
!    in case of non-vanishing orography
!
!
!
!!**  METHOD
!!    ------
!!      The coefficients for the flat operator are first computed. Then the
!!    pressure equation is solved and the pressure gradient is added to the wind
!!    components in order to render this wind field non-divergent.
!!
!!    EXTERNAL
!!    --------
!!    TRID     : to compute coefficients for the flat operator
!!    PRESSURE : to solve the pressure equation and add the pressure term to
!!               wind
!!    MXM,MYM,MZM : to average a field at mass point in the x,y,z directions
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF    : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!
!!      Module MODD_LUNIT   :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!
!!      Module MODD_GRID_n  : contains grid variables
!!         XMAP,XXHAT,XYHAT,XZZ
!!
!!      Module MODD_REF_n   : contains reference state variables 
!!         XRHODJ,XTHVREF,XEXNREF,XRHODREF,XRVREF
!!
!!      Module MODD_REF_n   : contains reference state variables 
!!         XLINMASS : lineic mass along the lateral boundaries
!!
!!      Module MODD_FIELD_n : contains prognostic variables
!!         XUT,XVT,XWT,XTHT,XRT
!! 
!!      Module MODD_DYN_n : contains parameters for the dynamics
!!         CPRESOPT : option for the pressure solver
!!         NITR     : number of iterations for the solver
!!         XRELAX   : relaxation coefficient used in the Richardson method
!! 
!!      Module MODD_LBC_n : contains parameters relative to the boundaries
!!         CLBCX    : choice of lateral boundary condition along x
!!         CLBCY    : choice of lateral boundary condition along y
!!
!!    REFERENCE
!!    ---------
!!      NONE
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     6/09/94
!!      J. Stein     4/11/94  put the pressure solver parameters in namelist
!!      J. Stein     2/12/94  source cleaning
!!      J.P. Lafore 03/01/95  call to PRESSURE to account for absolute pressure
!!      J. Stein    17/01/95  bug in the pressure call
!!      J. Stein    15/03/95  remove R from the historical variables
!!      J.Stein and J.P. lafore 17/04/96 new version including the way to choose
!!            the model number and the instant where the projection is performed
!!      Stein,Lafore 14/01/97 new anelastic equations
!!      M.Faivre    2014
!!      M.Moge      08/2015   removing UPDATE_HALO_ll(XRHODJ) + EXTRAPOL on ZRU and ZRV in part 3.1
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_MODELN_HANDLER
!
USE MODD_CONF    ! declarative modules
USE MODD_LUNIT
USE MODD_PARAMETERS
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_METRICS_n
USE MODD_REF_n
USE MODD_FIELD_n
USE MODD_DYN_n
USE MODD_LBC_n
!
USE MODI_TRIDZ    ! interface modules
USE MODI_PRESSUREZ
USE MODE_SPLITTINGZ_ll
USE MODI_SHUMAN
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODE_MPPDB
USE MODE_EXTRAPOL
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments :
!
REAL, OPTIONAL                 :: PRESIDUAL
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IRESP                  ! return code
INTEGER :: IIY,IJY                ! same variable for Y decomposition
INTEGER :: ITCOUNT                ! counter value of temporal loop set to 1 ( this
                                  ! means that no guess of the pressure is available for
                                  ! the pressure solver
REAL    :: ZDXHATM                   ! mean grid increment in the x direction
REAL    :: ZDYHATM                   ! mean grid increment in the y direction
REAL, DIMENSION (SIZE(XRHODJ,3)) :: ZRHOM   !  mean of XRHODJ on the plane x y
                                            !  localized at a mass level
!
REAL, DIMENSION(SIZE(XRHODJ,3))    :: ZAF,ZCF  ! vector giving the nonvanishing
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBFY    ! elements of the tri-diag matrix
                                               ! in the pressure equation
REAL, DIMENSION(:), ALLOCATABLE   :: ZTRIGSX   ! arrays of sin or cos values for
REAL, DIMENSION(:), ALLOCATABLE   :: ZTRIGSY   ! the FFT in x and y directions
INTEGER, DIMENSION(19)            :: IIFAXX    ! decomposition in prime numbers
INTEGER, DIMENSION(19)            :: IIFAXY    ! for the FFT in x and y
                                               ! directions
REAL, DIMENSION(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)) :: ZPABST
                                               ! Potential at time t
REAL, DIMENSION(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)) :: ZRU,ZRV,ZRW
                                               ! Rhod * (U,V,W)
REAL, DIMENSION(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)) :: ZTH
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZRR
!
INTEGER                       ::       IRR     ! Total number of water variables
INTEGER                       ::       IRRL    ! Number of liquid water variables
INTEGER                       ::       IRRI    ! Number of solid water variables
REAL                          ::  ZDRYMASST    ! Mass of dry air Md
REAL                          ::  ZREFMASS     ! Total mass of the ref. atmosphere
REAL                          ::  ZMASS_O_PHI0 ! Mass / Phi0
LOGICAL                       ::  GCLOSE_OUT   ! switch for the LFI writing
CHARACTER (LEN= 28)           ::  YFMFILE      ! virtual FM file
INTEGER                       ::  IMI          ! model index
!JUAN
INTEGER                               ::  IIU_B,IJU_B,IKU
INTEGER                               ::  IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll
REAL, DIMENSION(:,:,:), ALLOCATABLE   ::  ZBFB,ZBF_SXP2_YP1_Z                         
!JUAN
!
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
!
!-------------------------------------------------------------------------------
!
!*       1.     PROLOGUE  :
!               --------
!
CALL GET_DIM_EXT_ll('Y',IIY,IJY)
IF (L2D) THEN
  ALLOCATE(ZBFY(IIY,IJY,SIZE(XRHODJ,3)))
ELSE
  ALLOCATE(ZBFY(IJY,IIY,SIZE(XRHODJ,3)))
ENDIF
ALLOCATE(ZTRIGSX(3*(NIMAX_ll+2*JPHEXT)))
ALLOCATE(ZTRIGSY(3*(NJMAX_ll+2*JPHEXT)))
!JUAN Z_SPLITING
IKU=SIZE(XRHODJ,3)
CALL GET_DIM_EXT_ll('B',IIU_B,IJU_B)
ALLOCATE(ZBFB(IIU_B,IJU_B,IKU))
CALL GET_DIM_EXTZ_ll('SXP2_YP1_Z',IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll)
ALLOCATE(ZBF_SXP2_YP1_Z(IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll))
!JUAN Z_SPLITING
CALL MPPDB_CHECK3D(XRHODJ,"anel_balancen1-::XRHODJ",PRECISION)
CALL MPPDB_CHECK3D(XUT,"anel_balancen1-::XUT",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       2.     PRESSURE SOLVER INITIALIZATION :
!               -------------------------------
!
!
CALL TRIDZ(CLBCX,CLBCY,XMAP,XDXHAT,XDYHAT,ZDXHATM,ZDYHATM,ZRHOM,        &
          ZAF,ZCF,ZTRIGSX,ZTRIGSY,IIFAXX,IIFAXY,XRHODJ,XTHVREF,XZZ,ZBFY,&
          ZBFB,ZBF_SXP2_YP1_Z) 
CALL MPPDB_CHECK3D(XRHODJ,"anel_balancen1-after TRIDZ::XRHODJ",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       3.      ANELASTIC CORRECTION :
!                ---------------------
!
!
!*       3.1     multiplication by RHODJ
!
!$20140710 UPHALO on XRHODJ
!CALL ADD3DFIELD_ll(TZFIELDS_ll,XRHODJ)
!CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
!CALL CLEANLIST_ll(TZFIELDS_ll)
CALL MPPDB_CHECK3D(XRHODJ,"anel_balancen3.1-after update halo::XRHODJ",PRECISION)
CALL MPPDB_CHECK3D(XUT,"anel_balancen3.1-after update halo::XUT",PRECISION)
CALL MPPDB_CHECK3D(XWT,"anel_balancen3.1-after update halo::XWT",PRECISION)
!
ZRU(:,:,:) = MXM(XRHODJ) * XUT(:,:,:)
ZRV(:,:,:) = MYM(XRHODJ) * XVT(:,:,:)
ZRW(:,:,:) = MZM(1,IKU,1,XRHODJ) * XWT(:,:,:)
ZTH(:,:,:) = XTHT(:,:,:)
ALLOCATE(ZRR(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),SIZE(XRT,4)))
ZRR(:,:,:,:) = XRT(:,:,:,:)
!20131112 appli update_halo_ll
CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRU)
CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRV)
CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRW)
CALL ADD3DFIELD_ll(TZFIELDS_ll, ZTH)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
CALL MPPDB_CHECK3D(ZRU,"anel_balancen3.1-after1stupdhalo::ZRU",PRECISION)
!$20131125 add extrapol on ZRU to have correct boundaries
!CALL EXTRAPOL('W',ZRU)  ! ZRU boundaries now correct
CALL MPPDB_CHECK3D(ZRU,"anel_balancen3.1-afterextrapol W::ZRU",PRECISION)
!20131126 add extrapol on ZRV to have correct boundaries
!CALL EXTRAPOL('S',ZRV)  ! ZRV boundaries now correct
CALL MPPDB_CHECK3D(ZRV,"anel_balancen3.1-afterextrapol S::ZRV",PRECISION)
CALL MPPDB_CHECK3D(ZRW,"anel_balancen3.1-afterextrapol S::ZRW",PRECISION)
!
!
!
!
!*       3.2     satisfy the anelastic constraint
!
ITCOUNT      =-1     ! no first guess of the pressure is available
ZPABST(:,:,:)= 0.    !       ==================CAUTION=====================
ZDRYMASST    = 0.    !      |   Initialization necessary for the           |
ZREFMASS     = 0.    !      |  computation of the absolute pressure,       |
ZMASS_O_PHI0 = 1.    !      |  which is here not needed                    |
IRR  = 0             !      |                                              |
IRRL = 0             !      |                                              |
IRRI = 0             !       ==============================================
GCLOSE_OUT=.FALSE.
YFMFILE='UNUSED'
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL PRESSUREZ(CLBCX,CLBCY,CPRESOPT,NITR,LITRADJ,ITCOUNT,XRELAX,IMI,  &
               XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY,ZDXHATM,ZDYHATM,ZRHOM, &
               ZAF,ZBFY,ZCF,ZTRIGSX,ZTRIGSY,IIFAXX,IIFAXY,            &
               IRR,IRRL,IRRI,ZDRYMASST,ZREFMASS,ZMASS_O_PHI0,         &
               ZTH,ZRR,XRHODREF,XTHVREF,XRVREF,XEXNREF, XLINMASS,     &
               ZRU,ZRV,ZRW,ZPABST,                                    &
               ZBFB,ZBF_SXP2_YP1_Z,PRESIDUAL                          )
!
CALL MPPDB_CHECK3D(XRHODJ,"anel_balancen3.2-after pressurez halo::XRHODJ",PRECISION)
CALL MPPDB_CHECK3D(ZRU,"anel_balancen3.2-after pressurez::ZRU",PRECISION)
CALL MPPDB_CHECK3D(ZRV,"anel_balancen3.2-after pressurez::ZRV",PRECISION)
!
DEALLOCATE(ZBFY,ZTRIGSX,ZTRIGSY,ZRR,ZBF_SXP2_YP1_Z)
!*       3.2     return to the historical variables
!
!20131112 appli update_halo_ll and associated operations
XUT(:,:,:) = ZRU(:,:,:) / MXM(XRHODJ)
XVT(:,:,:) = ZRV(:,:,:) / MYM(XRHODJ)
XWT(:,:,:) = ZRW(:,:,:) / MZM(1,IKU,1,XRHODJ)
!20131112 appli update_halo_ll to XUT,XVT,XWT
CALL ADD3DFIELD_ll(TZFIELDS_ll, XUT)
CALL ADD3DFIELD_ll(TZFIELDS_ll, XVT)
CALL ADD3DFIELD_ll(TZFIELDS_ll, XWT)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
CALL MPPDB_CHECK3D(XUT,"anel_balancen3.2-afterupdhalo::XUT",PRECISION)
CALL MPPDB_CHECK3D(XVT,"anel_balancen3.2-afterupdhalo::XVT",PRECISION)
!20131125 apply extrapol to fix boundary issue in //
CALL EXTRAPOL('W',XUT)
CALL EXTRAPOL('S',XVT)
CALL MPPDB_CHECK3D(XUT,"anel_balancen3.2-after extrapolW::XUT",PRECISION)
CALL MPPDB_CHECK3D(XVT,"anel_balancen3.2-after extrapolS::XVT",PRECISION)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ANEL_BALANCE_n
