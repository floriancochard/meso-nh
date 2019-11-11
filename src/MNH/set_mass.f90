!MNH_LIC Copyright 2010-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ########################
      MODULE MODI_SET_MASS
!     ########################
!
INTERFACE
!
SUBROUTINE SET_MASS(TPFILE,OPROFILE_IN_PROC, PZFLUX_PROFILE,                           &
                    KILOC,KJLOC,PZS_MX,PZMASS_MX,PZFLUX_MX,PPGROUND,                   &
                    PTHVM,PMRM,PUW,PVW,OSHIFT,OBOUSS,PJ,HFUNU,HFUNV,PMRCM,PMRIM,PCORIOZ)
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),        INTENT(IN) :: TPFILE    ! File characteristics
LOGICAL,                INTENT(IN) :: OPROFILE_IN_PROC ! initialization profile in current processor
REAL, DIMENSION(:),     INTENT(IN) :: PZFLUX_PROFILE   ! Z at flux points on the initialization point column
INTEGER,                INTENT(IN) :: KILOC     ! I Localisation of vertical profile
INTEGER,                INTENT(IN) :: KJLOC     ! J Localisation of vertical profile
REAL, DIMENSION(:,:),   INTENT(IN) :: PZS_MX    ! zs on the mixed grid
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZFLUX_MX ! Z of mixed grid at flux points
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZMASS_MX ! Z of mixed grid at mass points
REAL,                   INTENT(IN) :: PPGROUND  ! pressure at PV ground
REAL, DIMENSION(:),     INTENT(IN) :: PUW,PVW   ! Wind at w model grid levels
REAL, DIMENSION(:),     INTENT(IN) :: PMRM      ! vapor mixing ratio at mass model
                                                ! grid levels 
REAL, DIMENSION(:),     INTENT(IN) :: PTHVM     ! Temperature at mass model grid levels
LOGICAL,                INTENT(IN) :: OSHIFT    ! logical switch for vertical shift
LOGICAL,                INTENT(IN) :: OBOUSS    ! logical switch for Boussinesq version
REAL, DIMENSION(:,:,:), INTENT(IN) :: PJ ! jacobien
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
REAL, DIMENSION(:),     INTENT(IN),OPTIONAL  :: PMRCM
REAL, DIMENSION(:),     INTENT(IN),OPTIONAL  :: PMRIM
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ ! Coriolis parameternn
                                               ! (exceptionnaly 3D array)
!
END SUBROUTINE SET_MASS
!
END INTERFACE
!

END MODULE MODI_SET_MASS
!
!
!     ##########################################################################
SUBROUTINE SET_MASS(TPFILE,OPROFILE_IN_PROC, PZFLUX_PROFILE,                           &
                    KILOC,KJLOC,PZS_MX,PZMASS_MX,PZFLUX_MX,PPGROUND,                   &
                    PTHVM,PMRM,PUW,PVW,OSHIFT,OBOUSS,PJ,HFUNU,HFUNV,PMRCM,PMRIM,PCORIOZ)
!     ##########################################################################
!
!!****  *SET_MASS * - routine to initialize mass and wind fields on MESONH grid
!!
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize mass and wind fields on the 3D
!     model grid (with orography) from a vertical profile of (u,v,thetav,r) on
!     the mixed grid (PZMASS_MX for mass points, PZFLUX_MX for flux points)
!
!       In the Boussinesq case, THVREF is taken constant and equal to the
!     value of THV in the first mass level, and RHOREF is also constant and
!     equal to P00/(Rd*THVREF)
!
!     The variation of U in y and z-directions is possible:
!        * HFUNU  ="ZZZ"   purely vertical profile u(y,z)=U(z)
!        * HFUNU  ="Y*Z"   separable fonction of y and z  u(y,z)=f(y)*U(z)
!        * HFUNU  ="Y,Z"   non-separable fonction of y and z  u(y,z)=f(y,z)
!     The variation of V in x and z-directions is possible:
!        * HFUNV  ="ZZZ"   purely vertical profile v(x,z)=V(z)
!        * HFUNV  ="X*Z"   separable fonction of x and z  v(x,z)=f(x)*V(z)
!        * HFUNV  ="X,Z"   non-separable fonction of x and z  v(x,z)=f(x,z)
!
!
!!**  METHOD
!!    ------
!!    
!!      - First we compute homogenous 3D fields from the vertical profil on 
!!        the mixed grid for theta and r
!!      - Then the vertical wind profile (on mixed grid) is considered :
!!          Following the type of geometry (cartesian or conformal projection),
!!      the wind component along x and y model axes are computed from the zonal
!!      and meridian components (PUW,PVW) at w-grid levels.
!!    The variation of U in y and z-directions is selected according to:
!!       * HFUNU  ="ZZZ"   purely vertical profile u(y,z)=PUW(z)
!!       * HFUNU  ="Y*Z"   separable fonction of y and z  u(y,z)=f(y)*PUW(z)
!!       * HFUNU  ="Y,Z"   non-separable fonction of y and z  u(y,z)=f(y,z)
!!                         (PUW information is ignored)
!!    The variation of V in x and z-directions is possible:
!!       * HFUNV  ="ZZZ"   purely vertical profile v(x,z)=PVW(z)
!!       * HFUNV  ="X*Z"   separable fonction of x and z  v(x,z)=f(x)*PVW(z)
!!       * HFUNV  ="X,Z"   non-separable fonction of x and z  v(x,z)=f(x,z)
!!                         (PVW information is ignored)
!!      - Interpolate fields on the MESONH grid by using a shifting function
!!        method ( as in PREP_REAL_CASE) with VER_INT_THERMO/VER_INT_DYN
!!      - In LGEOSBAL=.TRUE. case (PCORIOZ is present in input of this
!!        subroutine), LSHIFT=.FALSE. and the mixed grid as no orography,
!!        SET_GEOSBAL is called to interpolate fields with the thermal wind
!!        balance
!!      - Compute anelastic reference profil by CALL SET_REFZ
!!
!!
!!
!!    AUTHOR
!!    ------
!!     G.TANGUY        * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original : oct 2010 
!!    Tout a été modifié pour se rapprocher de PREP_REAL_CASE
!!    J. Escobar  27/03/2012 modif for reprod sum
!!    V.Masson    12/08/13  Parallelization of the initilization profile
!!    M.Moge      08/2015   add UPDATE_HALO_ll on XTHT, ZTHV3D, XRT(:,:,1,:) after computation
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!    
!-------------------------------------------------------------------------------
!!
! use des modules
USE MODD_GRID_n ! declarative modules
USE MODD_GRID
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_CONF
USE MODD_CONF_n
USE MODD_FIELD_n
USE MODD_CST
USE MODD_REF
USE MODD_PARAMETERS
USE MODD_DIM_n
!
USE MODE_GATHER_ll
USE MODE_ll
!
USE MODI_VER_INT_THERMO ! interface modules
USE MODI_WATER_SUM
USE MODI_SET_REFZ
USE MODI_FUN
USE MODI_VER_INT_DYN
USE MODI_SHUMAN
USE MODI_COMPUTE_EXNER_FROM_GROUND
USE MODI_COMPUTE_EXNER_FROM_TOP
USE MODI_SET_GEOSBAL
USE MODE_REPRO_SUM
USE MODE_MPPDB
USE MODE_SUM_ll, ONLY : SUM_DIM1_DD_ll

IMPLICIT NONE
!
!*       0.1   Declarations of arguments :
!
TYPE(TFILEDATA),        INTENT(IN)             :: TPFILE    ! File characteristics
LOGICAL,                INTENT(IN)             :: OPROFILE_IN_PROC ! initialization profile in current processor
REAL, DIMENSION(:),     INTENT(IN)             :: PZFLUX_PROFILE   ! Z at flux points on the initialization point column
INTEGER,                INTENT(IN)             :: KILOC     ! I Localisation of vertical profile
INTEGER,                INTENT(IN)             :: KJLOC     ! J Localisation of vertical profile
REAL, DIMENSION(:,:),   INTENT(IN)             :: PZS_MX    ! zs on the mixed grid
REAL, DIMENSION(:,:,:), INTENT(IN)             :: PZFLUX_MX ! Z of mixed grid at flux points
REAL, DIMENSION(:,:,:), INTENT(IN)             :: PZMASS_MX ! Z of mixed grid at mass points
REAL,                   INTENT(IN)             :: PPGROUND  ! pressure at ground
REAL, DIMENSION(:),     INTENT(IN)             :: PUW,PVW   ! Wind at w model grid levels
REAL, DIMENSION(:),     INTENT(IN)             :: PMRM      ! vapor mixing ratio at mass model grid levels 
REAL, DIMENSION(:),     INTENT(IN)             :: PTHVM     ! Temperature at mass model grid levels
LOGICAL,                INTENT(IN)             :: OSHIFT    ! logical switch for vertical shift
LOGICAL,                INTENT(IN)             :: OBOUSS    ! logical switch for Boussinesq version
REAL, DIMENSION(:,:,:), INTENT(IN)             :: PJ        ! jacobien
CHARACTER(LEN=*),       INTENT(IN)             :: HFUNU     ! type of variation of U in y direction
CHARACTER(LEN=*),       INTENT(IN)             :: HFUNV     ! type of variation of V in x direction
REAL, DIMENSION(:),     INTENT(IN),  OPTIONAL  :: PMRCM
REAL, DIMENSION(:),     INTENT(IN),  OPTIONAL  :: PMRIM
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  :: PCORIOZ ! Coriolis parameter (exceptionnaly 3D array)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIU,IJU,IIB,IIE,IJB,IJE,IKE,IKU
INTEGER :: JI,JK,JJ
!* 0.2.1  fields on the mixed grid : 
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZTHV3D_MX     ! virtual potential temperature (mass level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZTHVREF3D     ! virtual potential temperature (mass level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT),NRR):: ZMR3D_MX      ! vapor mixing ratio (mass level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZUW3D_FL      ! zonal wind component (flux level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZVW3D_FL      ! meridian wind component (flux level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZPMHP_MX      ! pressure minus hyd. pressure (mass level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHOD_MX      ! local rhod (mass level)
REAL,DIMENSION(SIZE(XZHAT))                            :: ZRHOD_PROFILE ! local rhod (mass level) at initialization profile column
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZPMASS_MX     ! pressure (mass level)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT))                :: ZEXNSURF2D_MX ! local Exner function at ground
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZHEXNFLUX_MX  ! local hyd. Exner function at flux points on the mixed grid
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZHEXNMASS_MX  ! local hyd. Exner function at mass points on the mixed grid
!
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT))                :: ZEXNTOP2D     ! Exner function at top
!* 0.2.2 for VER_INT_THERMO call
REAL                                                   :: ZDIAG         ! diagnostics computing time
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZTHV3D        ! virtual potential temperature on MESONH grid
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZPMHP         ! pressure minus hyd. pressure on MNH grid with orography (mass level)
!* 0.2.3 for wind (application of HFUNU/HFUNV)
!!$REAL, DIMENSION(SIZE(XXHAT),1,1)                       :: ZNFLX_TOT     ! total normalized mass flux
REAL, DIMENSION(SIZE(XYHAT),SIZE(XZHAT))               :: ZUYZ          ! vertical variations for U
!!$REAL, DIMENSION(1,SIZE(XYHAT),1)                       :: ZNFLY_TOT     ! total normalized mass flux
REAL, DIMENSION(SIZE(XXHAT),SIZE(XZHAT))               :: ZVXZ          ! vertical variations for V
REAL, DIMENSION(:,:), ALLOCATABLE                      :: ZNFLX_TOT_ll,ZNFLY_TOT_ll
INTEGER                                                :: IXOR_ll,IYOR_ll! origin's coordinates of extended subdomain
REAL                                                   :: ZD1           ! DELTA1 (switch 0/1) for thinshell approximation
REAL, DIMENSION(SIZE(XZHAT))                           :: ZDZZFLUX_MX   ! deltaz for mixed grid at flux level
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))   :: ZDXX,ZDYY     ! metric coefficients dxx,dyy for the MNH grid
                                                                        ! without orography
REAL,DIMENSION(:,:),ALLOCATABLE                        :: ZGAMMA,ZGAMMA_ll ! K(lambda-lambda0) -beta
REAL                                                   :: ZGAMMALOC
INTEGER                                                :: IINFO_ll  
REAL                                                   :: ZRADSDG        !Pi/180
INTEGER                                                :: IIU_ll,IJU_ll
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHODU_MX      ! horizontal momentum components on
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHODV_MX      ! the mixed grid
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHODUA        ! horizontal momentum components on
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHODVA        ! the MESONH Arakawa A grid
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHODJU        ! horizontal momentum components on
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHODJV        ! the MESONH Arakawa C grid
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZHEXNFLUX      ! local hyd. Exner function at flux points (MNH grid)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZHEXNMASS      ! local hyd. Exner function at mass points (MNH grid)
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))    :: ZRHOD          ! dry density on MESO-NH grid
!
!!$INTEGER                                                :: IIBP,IIEP,IJBP,IJEP
REAL, DIMENSION(:,:), ALLOCATABLE                      :: ZNFLXZ_TOT,ZNFLYZ_TOT
REAL, DIMENSION(:)  , ALLOCATABLE                      :: ZNFLXZ_TOT_ll,ZNFLYZ_TOT_ll ! total normalized mass flux
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!*	 1.     PROLOGUE : INITIALIZE SOME CONSTANTS
!-------------------------------------------------------------------------------
!
IIU=SIZE(XXHAT)
IJU=SIZE(XYHAT)
IKU=SIZE(XZHAT)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!!$CALL GET_PHYSICAL_ll(IIBP,IJBP,IIEP,IJEP)
IKE=IKU-JPVEXT
ZRADSDG = XPI/180.
!
!-------------------------------------------------------------------------------
!*        2. COMPUTE FIELDS ON MIXED GRID
!-------------------------------------------------------------------------------
!
!* 2.1 Compute fields from vertical profile on mixed grid
!
DO JI=1,IIU
  DO JJ=1,IJU
    ZTHV3D_MX(JI,JJ,:)=PTHVM(:)
    ZMR3D_MX(JI,JJ,:,1)=PMRM(:)
  ENDDO
ENDDO
ZPMHP_MX(:,:,:)=0.
ZMR3D_MX(:,:,:,2:)=0.
IF(PRESENT(PMRCM)) THEN
  DO JI=1,IIU
    DO JJ=1,IJU      
      ZMR3D_MX(JI,JJ,:,2)= PMRCM(:)
    ENDDO
  ENDDO
ENDIF
IF(PRESENT(PMRIM)) THEN
  ZMR3D_MX(:,:,:,3:)= 0      
  DO JI=1,IIU
    DO JJ=1,IJU             
      ZMR3D_MX(JI,JJ,:,4)= PMRIM(:)
        ENDDO
  ENDDO
  
ENDIF
!------------------------------
!* 2.2 compute exner function on mixed grid
!
ZEXNSURF2D_MX(:,:)=(PPGROUND/XP00)**(XRD/XCPD)    
CALL COMPUTE_EXNER_FROM_GROUND(ZTHV3D_MX,PZFLUX_MX,&
            ZEXNSURF2D_MX,ZHEXNFLUX_MX,ZHEXNMASS_MX)
ZEXNTOP2D(:,:)=ZHEXNFLUX_MX(:,:,IKE+1)
ZPMASS_MX(:,:,:)=XP00*(ZHEXNMASS_MX(:,:,:))**(XCPD/XRD)
ZRHOD_MX(:,:,:)=ZPMASS_MX(:,:,:)/(ZPMASS_MX(:,:,:)/XP00)**(XRD/XCPD) &
                 /(XRD*ZTHV3D_MX(:,:,:)*(1.+WATER_SUM(ZMR3D_MX(:,:,:,:))))

XEXNTOP=SUM_DD_R2_ll(ZHEXNFLUX_MX(IIB:IIE,IJB:IJE,IKE+1))/FLOAT(NIMAX_ll*NJMAX_ll)


!------------------------------
!*  2.3 Rotate wind in model axis and take into account variations in x,y
!      directions on the mixed grid
!
CALL GET_OR_ll('B',IXOR_ll,IYOR_ll)
!
IF (OPROFILE_IN_PROC) THEN
  ZRHOD_PROFILE(:) = ZRHOD_MX(KILOC-IXOR_ll+1,KJLOC-IYOR_ll+1,:)
ELSE
  ZRHOD_PROFILE(:) = 0.
END IF
DO JK = 1,IKU
  CALL REDUCESUM_ll(ZRHOD_PROFILE(JK), IINFO_ll)
END DO

IF (LTHINSHELL .OR. LCARTESIAN) THEN
  ZD1=0.
ELSE
  ZD1=1.
END IF
IIU_ll=NIMAX_ll+2*JPHEXT
IJU_ll=NJMAX_ll+2*JPHEXT
IF (LCARTESIAN) THEN                 ! cartesian geometry
  ZGAMMALOC = -XBETA *ZRADSDG
  DO JK=1,IKU ;  DO JJ=1,IJU ;  DO JI=1,IIU 
    ZDXX(JI,JJ,JK) = XDXHAT(JI)
    ZDYY(JI,JJ,JK) = XDYHAT(JJ)
  ENDDO ; ENDDO ; ENDDO ;
  ZDXX  = MXM( ZDXX )
  ZDYY  = MYM( ZDYY )
ELSE                                 ! conformal projection
  ALLOCATE(ZGAMMA(IIU,IJU))
  ALLOCATE(ZGAMMA_ll(IIU_ll,IJU_ll))
  ZGAMMA(:,:) = XRPK * (XLON(:,:) -XLON0) * ZRADSDG -(XBETA *ZRADSDG)
  CALL GATHERALL_FIELD_ll('XY',ZGAMMA,ZGAMMA_ll,IINFO_ll)
  ZGAMMALOC = ZGAMMA_ll(KILOC,KJLOC)
  DEALLOCATE(ZGAMMA,ZGAMMA_ll)
  DO JK=1,IKU ;  DO JJ=1,IJU ;  DO JI=1,IIU 
    ZDXX(JI,JJ,JK) = ( 1.+ ZD1*XZHAT(JK)/XRADIUS ) * ( XDXHAT(JI) /XMAP(JI,JJ) ) ! XDXHAT(JI)
    ZDYY(JI,JJ,JK) = ( 1.+ ZD1*XZHAT(JK)/XRADIUS ) * ( XDYHAT(JJ) /XMAP(JI,JJ) ) ! XDYHAT(JJ)
  ENDDO ; ENDDO ; ENDDO ;
  ZDXX  = MXM(MZF(1,IKU,1,ZDXX))
  ZDYY  = MYM(MZF(1,IKU,1,ZDYY))
END IF
!
SELECT CASE(HFUNU)
  CASE('ZZZ')
    DO JK=1,IKU
      ZUYZ(:,JK)=PUW(JK)*COS(ZGAMMALOC)-PVW(JK)*SIN(ZGAMMALOC)
    END DO
    !
  CASE('Y*Z')
    DO JK=1,IKU
      ZUYZ(:,JK)=FUNUY(IJU)*(PUW(JK)*COS(ZGAMMALOC) -   &
                             PVW(JK)*SIN(ZGAMMALOC)   )
    END DO
  CASE('Y,Z')
    ZUYZ(:,:)=FUNUYZ(IJU,IKU,PZFLUX_PROFILE(:))
END SELECT
!
ZDZZFLUX_MX(2:IKU) = PZFLUX_PROFILE(2:IKU)-PZFLUX_PROFILE(1:IKU-1)
ZDZZFLUX_MX(1) = ZDZZFLUX_MX(2)
!
!!$ZNFLX_TOT  = 0.
!!$DO JK = 2,IKU-1
!!$  DO JJ=IJB,IJE
!!$    ZNFLX_TOT(:,1,1)=ZNFLX_TOT(:,1,1)+ZDYY(:,JJ,JK)*ZUYZ(JJ,JK)*ZDZZFLUX_MX(JK)* &
!!$                   ZRHOD_PROFILE(JK)
!!$  END DO
!!$END DO
!JUAN
ALLOCATE(ZNFLXZ_TOT(IIU,IJU))
ZNFLXZ_TOT  = 0.
DO JK = 2,IKU-1
  DO JJ=IJB,IJE
    ZNFLXZ_TOT(:,JJ)=ZNFLXZ_TOT(:,JJ)+ZDYY(:,JJ,JK)*ZUYZ(JJ,JK)*ZDZZFLUX_MX(JK)* &
                   ZRHOD_PROFILE(JK)
  END DO
END DO
!
!!$ALLOCATE(ZNFLX_TOT_ll(IIU_ll,1))
!!$CALL SUM_DIM1_ll(ZNFLX_TOT,ZNFLX_TOT_ll,IINFO_ll)
!!$ZNFLX_TOT_ll=SIGN(1.,ZNFLX_TOT_ll)*MAX(ABS(ZNFLX_TOT_ll),TINY(ZNFLX_TOT_ll))
!
ALLOCATE(ZNFLXZ_TOT_ll(IIU_ll))
CALL SUM_DIM1_DD_ll(ZNFLXZ_TOT,ZNFLXZ_TOT_ll,KDIM=2,KINFO=IINFO_ll)
ZNFLXZ_TOT_ll=SIGN(1.,ZNFLXZ_TOT_ll)*MAX(ABS(ZNFLXZ_TOT_ll),TINY(ZNFLXZ_TOT_ll))

DO JI=1,IIU
!!$  ZUW3D_MX(JI,:,:)=ZUYZ(:,:)* ( ZNFLX_TOT_ll(KILOC,1)/ZNFLX_TOT_ll(IXOR_ll-1+JI,1) ) ! add () for reproductibility
  ZUW3D_FL(JI,:,:)=ZUYZ(:,:)* ( ZNFLXZ_TOT_ll(KILOC)/ZNFLXZ_TOT_ll(IXOR_ll-1+JI) ) 
END DO

!!$DEALLOCATE(ZNFLX_TOT_ll)
DEALLOCATE(ZNFLXZ_TOT,ZNFLXZ_TOT_ll)
!
!
SELECT CASE(HFUNV)
  CASE('ZZZ')
    DO JK=1,IKU
      ZVXZ(:,JK)=PUW(JK)*SIN(ZGAMMALOC)+PVW(JK)*COS(ZGAMMALOC)
    END DO
    !
  CASE('X*Z')
    DO JK=1,IKU
      ZVXZ(:,JK)=FUNVX(IIU)*(PUW(JK)*SIN(ZGAMMALOC) +   &
                             PVW(JK)*COS(ZGAMMALOC)   )
    END DO
  CASE('X,Z')
    ZVXZ(:,:)=FUNVXZ(IIU,IKU,PZFLUX_PROFILE(:))
END SELECT
!
!!$ZNFLY_TOT  = 0.
!!$DO JK = 2,IKU-1
!!$  DO JI=IIB,IIE
!!$    ZNFLY_TOT(1,:,1) = ZNFLY_TOT(1,:,1) + ZDXX(JI,:,JK)*ZVXZ(JI,JK)* ZDZZFLUX_MX(JK)* &
!!$                   ZRHOD_PROFILE(JK)
!!$  END DO
!!$END DO
!
ALLOCATE(ZNFLYZ_TOT(IIU,IJU))
ZNFLYZ_TOT  = 0.
DO JK = 2,IKU-1
  DO JI=IIB,IIE
    ZNFLYZ_TOT(JI,:) = ZNFLYZ_TOT(JI,:) + ZDXX(JI,:,JK)*ZVXZ(JI,JK)* ZDZZFLUX_MX(JK)* &
                   ZRHOD_PROFILE(JK)
  END DO
END DO
!
!!$ALLOCATE(ZNFLY_TOT_ll(IJU_ll,1))
!!$CALL SUM_DIM1_ll(ZNFLY_TOT,ZNFLY_TOT_ll,IINFO_ll)
!!$ZNFLY_TOT_ll=SIGN(1.,ZNFLY_TOT_ll)*MAX(ABS(ZNFLY_TOT_ll),TINY(ZNFLY_TOT_ll))
!
ALLOCATE(ZNFLYZ_TOT_ll(IJU_ll))
CALL SUM_DIM1_DD_ll(ZNFLYZ_TOT,ZNFLYZ_TOT_ll,KDIM=1,KINFO=IINFO_ll)
ZNFLYZ_TOT_ll=SIGN(1.,ZNFLYZ_TOT_ll)*MAX(ABS(ZNFLYZ_TOT_ll),TINY(ZNFLYZ_TOT_ll))
!
!
DO JJ=1,IJU
!!$    ZVW3D_FL(:,JJ,:)= ZVXZ(:,:) * ( ZNFLY_TOT_ll(KJLOC,1)/ZNFLY_TOT_ll(IYOR_ll-1+JJ,1) ) ! add () for reproductibility
    ZVW3D_FL(:,JJ,:)= ZVXZ(:,:) * ( ZNFLYZ_TOT_ll(KJLOC)/ZNFLYZ_TOT_ll(IYOR_ll-1+JJ) )
END DO

  CALL MPPDB_CHECK3DM("SET_MASS:ZUW3D_FL,ZVW3D_FL",PRECISION,&
                   & ZUW3D_FL,ZVW3D_FL   )

!!$DEALLOCATE(ZNFLY_TOT_ll)
DEALLOCATE(ZNFLYZ_TOT,ZNFLYZ_TOT_ll)
!
!-------------------------------------------------------------------------------
!*                     3. INTERPOLATION ON MNH GRID
!-------------------------------------------------------------------------------
!

IF (PRESENT(PCORIOZ)) THEN
  CALL SET_GEOSBAL(ZUW3D_FL,ZVW3D_FL,PTHVM,PMRM, &
                    KILOC,KJLOC,OBOUSS,ZTHV3D,PCORIOZ)
  CALL COMPUTE_EXNER_FROM_TOP(ZTHV3D,XZZ,ZEXNTOP2D,ZHEXNFLUX,ZHEXNMASS)
  XPABSM(:,:,:)=XP00*ZHEXNMASS(:,:,:) ** (XCPD/XRD)
ELSE
! 
! Interpolation of theta and r
!
 IF (SIZE(ZTHV3D_MX,3) > 3) THEN
  CALL VER_INT_THERMO(TPFILE,OSHIFT,ZTHV3D_MX,ZMR3D_MX,PZS_MX,PZS_MX,PZMASS_MX,&
                      PZFLUX_MX,ZPMHP_MX,ZEXNTOP2D, &
                      ZTHV3D,XRT,ZPMHP,ZDIAG)
 ELSE
   ZTHV3D = ZTHV3D_MX
   XRT    = ZMR3D_MX
   ZDIAG  = 0.
 END IF
  XTHT(:,:,:)=ZTHV3D(:,:,:)*(1.+WATER_SUM(XRT(:,:,:,:)))/(1.+XRV/XRD*XRT(:,:,:,1))
  ZTHV3D(:,:,1)=ZTHV3D(:,:,2)
  XTHT(:,:,1)=XTHT(:,:,2)
  XRT(:,:,1,:)=XRT(:,:,2,:)
NULLIFY( TZFIELDS_ll )
CALL ADD3DFIELD_ll(TZFIELDS_ll,XTHT)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZTHV3D)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XRT(:,:,1,:))
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)

!
  IF (NRR>=3) THEN
    WHERE  (XRT(:,:,:,3)<1.E-20)
      XRT(:,:,:,3)=0.
    END WHERE
  END IF
!
! Interpolation of the wind      
!
  ZRHODU_MX=MZF(1,IKU,1,ZUW3D_FL)*ZRHOD_MX
  ZRHODV_MX=MZF(1,IKU,1,ZVW3D_FL)*ZRHOD_MX
  CALL MPPDB_CHECK3DM("SET_MASS:ZRHODU_MX,ZRHODV_MX,PZFLUX_MX,PZMASS_MX",PRECISION,&
                   &  ZRHODU_MX,ZRHODV_MX,PZFLUX_MX,PZMASS_MX  )
  CALL VER_INT_DYN(OSHIFT,ZRHODU_MX,ZRHODV_MX,PZFLUX_MX,PZMASS_MX,PZS_MX,ZRHODUA,ZRHODVA)
  ZRHODJU(:,:,:)=MXM(ZRHODUA(:,:,:)*PJ(:,:,:))
  ZRHODJV(:,:,:)=MYM(ZRHODVA(:,:,:)*PJ(:,:,:))
  CALL COMPUTE_EXNER_FROM_TOP(ZTHV3D,XZZ,ZEXNTOP2D,ZHEXNFLUX,ZHEXNMASS)
  XPABST(:,:,:)=ZPMHP(:,:,:) + XP00*ZHEXNMASS(:,:,:) ** (XCPD/XRD)
  ZRHOD(:,:,:)=XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**(XRD/XCPD) &
            /(XRD*XTHT(:,:,:)*(1.+XRV/XRD*XRT(:,:,:,1)))
  XUT(:,:,:)=ZRHODJU(:,:,:)/MXM(ZRHOD(:,:,:)*PJ(:,:,:))
  XVT(:,:,:)=ZRHODJV(:,:,:)/MYM(ZRHOD(:,:,:)*PJ(:,:,:))
  XWT(:,:,:)=0
  CALL MPPDB_CHECK3DM("SET_MASS:XVT,ZRHODJV,PJ,ZRHODVA",PRECISION,&
                   &    XVT,ZRHODJV,PJ,ZRHODVA )
  ENDIF

!
!-------------------------------------------------------------------------------
!*                   4. COMPUTE ANELASTIC REFERENCE (PV)
!-------------------------------------------------------------------------------
!
IF (.NOT. OBOUSS) THEN
  DEALLOCATE(XTHVREFZ)
  DEALLOCATE(XRHODREFZ)
  CALL SET_REFZ(ZTHV3D,XRT(:,:,:,1))
ELSE 
  IF (OPROFILE_IN_PROC) THEN
    XTHVREFZ(:) = ZTHV3D(KILOC-IXOR_ll+1,KJLOC-IYOR_ll+1,2)
  ELSE
    XTHVREFZ(:) = 0.
  END IF
  DO JK = 1,IKU
    CALL REDUCESUM_ll(XTHVREFZ(JK), IINFO_ll)
  END DO

  XRHODREFZ(:) = XP00/ (XRD* XTHVREFZ(:))
  ZTHVREF3D(:,:,:)=XTHVREFZ(2)
  CALL COMPUTE_EXNER_FROM_GROUND(ZTHVREF3D,PZFLUX_MX,&
          ZEXNSURF2D_MX,ZHEXNFLUX,ZHEXNMASS)

  XEXNTOP=SUM_DD_R2_ll(ZHEXNFLUX(IIB:IIE,IJB:IJE,IKE+1))/FLOAT(NIMAX_ll*NJMAX_ll)

  ZEXNTOP2D=ZHEXNFLUX(:,:,IKE+1)
  CALL COMPUTE_EXNER_FROM_TOP(ZTHVREF3D,XZZ,ZEXNTOP2D,ZHEXNFLUX,ZHEXNMASS)
  XPABST(:,:,:)=ZPMHP(:,:,:) + XP00*ZHEXNMASS(:,:,:) ** (XCPD/XRD)   
ENDIF
!---------------------------------------------------------------------------------
END SUBROUTINE SET_MASS     
