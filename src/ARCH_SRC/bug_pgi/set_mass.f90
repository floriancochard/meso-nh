!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_ideal 2006/07/06 15:17:49
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_MASS
!     ####################
!
INTERFACE
!
SUBROUTINE SET_MASS(PUW,PVW,PTHVM,PMRM,PZHATM,PPGROUND,PTHVGROUND,PZGROUND, &
                    HFUNU,HFUNV,KILOC,KJLOC,OBOUSS,OPV_PERT,ORMV_BL,PCORIOZ)
!
REAL, DIMENSION(:), INTENT(IN)    :: PUW,PVW ! Wind at w model grid levels
                                             ! (without orography)  : zonal and
                                             ! meridian components
REAL, DIMENSION(:), INTENT(IN)    :: PTHVM   ! Temperature at mass model grid
                                             ! levels (without orography)
REAL, DIMENSION(:), INTENT(IN)    :: PMRM    ! vapor mixing ratio at mass model
                                             ! grid levels (without orography)
REAL, DIMENSION(:), INTENT(IN)    :: PZHATM  ! Height of mass model grid levels
                                             ! (without orography)
REAL              , INTENT(IN)    :: PPGROUND! pressure at the ground heigth of
                                             ! the RS
REAL              , INTENT(IN)    :: PTHVGROUND! potential temperature at the
                                             ! ground height of the RS
REAL              , INTENT(IN)    :: PZGROUND! ground height of the RS
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
LOGICAL,                INTENT(IN)  :: OPV_PERT! logical switch for PV inversion
LOGICAL,                INTENT(IN)  :: ORMV_BL! logical switch for remouve boundary layer
!
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ  ! Coriolis parameter
                                                ! (exceptionally 3D array)
!
END SUBROUTINE SET_MASS
!
END INTERFACE
!
END MODULE MODI_SET_MASS
!
!
!     ##########################################################################
      SUBROUTINE SET_MASS(PUW,PVW,PTHVM,PMRM,PZHATM,PPGROUND,PTHVGROUND,PZGROUND, &
                    HFUNU,HFUNV,KILOC,KJLOC,OBOUSS,OPV_PERT,ORMV_BL,PCORIOZ)
!     ##########################################################################
!
!!****  *SET_MASS * - routine to initialize mass and wind fields using thermal wind
!!                    balance
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize mass and wind fields on the 3D
!     model grid (with orography) from a vertical profile of (U,V,thetav,r).
!     The vertical profile of wind is given at w-grid levels  (without
!     orography), i.e. at height XZHAT, with u along latitude circles and v
!     along meridians.
!     The vertical profile of mass variable is given at mass-grid levels
!     (without orography), i.e. at height PZHATM.
!     The variation of U in y and z-directions is possible:
!        * HFUNU  ="ZZZ"   purely vertical profile u(y,z)=U(z)
!        * HFUNU  ="Y*Z"   separable fonction of y and z  u(y,z)=f(y)*U(z)
!        * HFUNU  ="Y,Z"   non-separable fonction of y and z  u(y,z)=f(y,z)
!     The variation of V in x and z-directions is possible:
!        * HFUNV  ="ZZZ"   purely vertical profile v(x,z)=V(z)
!        * HFUNV  ="X*Z"   separable fonction of x and z  v(x,z)=f(x)*V(z)
!        * HFUNV  ="X,Z"   non-separable fonction of x and z  v(x,z)=f(x,z)
!       The thermal wind balance is used to compute the potential virtual
!     temperature from wind field and the vertical profile of thetav.
!       The vapor mixing ratio is taken  uniform in horizontal plane, i.e.
!     its vertical profile is only interpolated on the mass-grid levels with
!     orography, i.e. at height XZZ.
!
!       In the Boussinesq case, THVREF is taken constant and equal to the
!     value of THV in the first mass level, and RHOREF is also constant and
!     equal to P00/(Rd*THVREF)
!
!       All the fields are interpolated at the end of the subroutine on
!     the 3D-model grid with orography.
!
!       The Coriolis parameter (f= 2 omega sin(lat)) is stored as side-products
!
!
!!**  METHOD
!!    ------
!!
!!      - The vertical wind profile is first considered :
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
!!      As side-products, the Coriolis parameter is also computed.  Note
!!    that the wind is computed from the RS or function informations in order to
!!    approximatively satisfy the anelastic constrain in any geometrical domain.
!!
!!      - The virtual potential  temperature is then computed by integration of
!!      the thermal wind balance from the wind field and from the vertical
!!      profile of thetav at point (KILOC,KJLOC). The thermal wind balance has
!!      the following discretized form :
!!      Lipps and Hemler case:
!!                                                ------x
!!       d thetav          f                 1     d V
!!       --------  = dxx  ---  [ thetavref ----- -------  ]
!!        d xbar           g                ---z  d zhat
!!                                          dzzhat
!!
!!                                                   ------y
!!       d thetav          f                    1      d U
!!       --------  = dyy  ---  [ - thetavref  -----  -------  ]
!!        d ybar          g                   ---z   d zhat
!!                                             dzzhat
!!      Durran and Modified Anelastic Equations' case:
!!                                                  ------x
!!       d 1/thetav          f        1        1     d V
!!       ----------  =-dxx  ---  [ --------- ----- -------  ]
!!        d xbar             g     thetavref  ---z  d zhat
!!                                            dzzhat
!!
!!                                                      ------y
!!       d 1/thetav          f          1         1      d U
!!       ----------  =-dyy  ---  [ -  --------- -----  -------  ]
!!        d ybar             g         thetavref ---z   d zhat
!!                                               dzzhat
!!
!!        with the variables localized as follows (in I,K) plane :
!!
!!        ---+-------------*------------------   zhatm(k)
!!           dxx     thetav,f,thetavref
!!
!!        -----------------#------------------   zhat(k)
!!                         U,V
!!
!!       The integration is done in four steps :
!!
!!            1st step :             integration
!!                                   ------>
!!  KJLOC                *---------*------*-----------------*---->
!!                       KILOC                                IKU   I
!!
!!            2nd step :
!!
!!         integration
!!             <-----
!!  KJLOC *---*------*---*---------*------*--------*---------*---->
!!        1              KILOC                                IKU   I
!!
!!            3rd step :
!!  KJU   ^
!!        |
!!        |*---*------*---*---------*------*-----------------------> I
!!        |^   ^      ^   ^         ^      ^
!!        ||   |      |   |         |      |integration
!!        |*---*------*---*---------*------*--------*---------*----> I
!!        |
!!        |
!!  KJLOC |*---*------*---*---------*------*--------*---------*----> I
!!         1              KILOC                                IKU
!!
!!            4th step :
!!  KJU   ^
!!        |
!!        |*---*------*---*---------*------*--------*---------*----> I
!!        |
!!        |
!!        |*---*------*---*---------*------*--------*---------*----> I
!!        |
!!        |
!!  KJLOC |*---*------*---*---------*------*--------*---------*----> I
!!        |
!!        |*---*------*---*---------*------*-----------------------> I
!!        ||   |      |   |         |      |
!!        |!   !      !   !         !      ! integration
!!     1  |*---*------*---*---------*------*--------*---------*----> I
!!         1              KILOC                               IKU
!!
!!      The anelastic refernce state is unknown at the beginning of the
!!      determination of the Thetav field. It is set to the vertical profile
!!      at the point (KILOC,KJLOC). In order to obtain weaker deviation from
!!      this reference state, we average the obtained thetav field on the whole
!!      horizontal domain and restart the thetav computation with this new
!!      reference state ( thetavref). This procedure is iterated iitr_ref times.
!!      One iteration seems nevertheless suffisant.
!!
!!      If the Coriolis parameter is not present (i.e. LGEOSBAL =.FALSE.),
!!      virtual potential temperature is set equal to the input vertical profile
!!      at (KILOC,KJLOC).
!!
!!
!!      - The anelastic reference Exner function (ZEXNREFZ)
!!     is deduced by integration of the hydrostatic relation
!!     from ZHAT=0 and the potential virtual temperature (XTHVREFZ).:
!!
!!        d Exnrefz      g        1
!!        -------- = - ---  -------------
!!        d z          Cpd    ----------z
!!                            thetavrefz
!!
!!        Then the anelastic reference dry density (XRHODREFZ) is determined
!!     as :
!!                                   Cpd/Rd -1
!!                          (Exnrefz)
!!         xrhodrefz = P00  ---------------------
!!                          Rd thetavrefz (1. + r)
!!
!!     So that, the 1D anelastic reference state is completely defined.
!!
!!     - The orography is then taken into account :
!!         The potential virtual temperature, the vapor mixing ratio and the wind
!!      are interpolated on the model grid (with orography).
!!
!!         From this interpolated potential virtual temperature and vapor mixing
!!     ratio, the potential dry temperature is computed.
!!     So that, the prognostic mass variables (thetad and r) are completely
!!     defined.
!!
!!      If the Coriolis parameter is not present (i.e. LGEOSBAL =.FALSE.),
!!    the wind components are simply interpolated (no geostraphic balance)
!!
!!      In the Boussinesq case, THVREF is taken constant and equal to the
!!    value of THV in the first mass level, and RHODREF is also constant and
!!    equal to P00/(Rd*THVREF)
!!
!!    Parallelisation of the geostrophic balance computation is achieved by
!!    using a Y slice domain decomposition and the GET_SLICE routine for
!!    X direction
!!    EXTERNAL
!!    --------
!!      FMLOOK : to retrieve logical unit number
!!      FUNUY  :  to compute a y-variation
!!      FUNVX  :  to compute a x-variation
!!      FUNUYZ :  to compute a y and z-variations
!!      FUNVXZ :  to compute a x and z-variations
!!
!!
!!      Module MODI_SHUMAN: interface for shuman operators
!!      Module MODI_FUNUY : interface for function FUNUY
!!      Module MODI_FUNVX : interface for function FUNVX
!!      Module MODI_FUNUYZ : interface for function FUNUYZ
!!      Module MODI_FUNVXZ : interface for function FUNVXZ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains parameters
!!        JPVEXT : Vertical External points number
!!
!!      Module MODD_CST        : contains physical constants
!!        XRADIUS : earth radius
!!        XPI   : Pi
!!        XRV   : Gas constant for vapor
!!        XRD   : Gas constant for dry air
!!        XCPD  : Specific heat for dry air at constant pressure
!!        XG    : Gravity constant
!!        XP00  : reference pressure
!!        XOMEGA : earth rotation
!!
!!      Module MODD_LUNIT1  : contains logical unit names
!!        CLUOUT : name of output-listing
!!
!!      Module MODD_CONF    : contains configuration variables for all models.
!!
!!        L2D        : logical for 2D model version
!!        LTHINSHELL : Logical for thinshell approximation
!!                     .TRUE. = Thinshell approximation done
!!        LCARTESIAN : Logical for cartesian geometry :
!!                    .TRUE.  = cartesian geometry
!!                    .FALSE. = conformal projection
!!        NVERB      : verbosity level for output-listing
!!
!!      Module MODD_GRID  : contains grid variables
!!        XLON0,XLAT0 : Reference longitude and latitude
!!                      for the conformal projection
!!        XBETA,XRPK  : Rotation angle and projection parameter
!!                      for the conformal projection
!!
!!      Module MODD_GRID1  : contains grid variables
!!        XXHAT  : Position x in the conformal
!!                 plane or on the cartesian plane
!!        XYHAT  : Position y in the conformal
!!                 plane or on the cartesian plane
!!        XDXHAT        : horizontal stretching in x
!!        XDYHAT        : horizontal stretching in y
!!        XLONOR,XLATOR : Longitude and latitude of the Origine point
!!                        for the conformal projection
!!        XLON,XLAT     : Longitude and latitude
!!        XZHAT   : height level (without orography)
!!        XZZ     : height z
!!
!!      Module MODD_REF     : contains reference state variables for all models
!!        XRHODREFZ : rhod(z) for reference state
!!                    without orography
!!        XTHVREFZ : Thetav(z) for reference state
!!                   without orography
!!        XEXNTOP  : Exner function at model top
!!
!!      Module MODD_FIELD1 : contains prognostic variables
!!        XTHM :  theta at time t
!!        XRM  :  Moist variables  at time t
!!        XUM  :  x -component of the wind
!!        XVM  :  y -component of the wind
!!        XWM  :  z -component of the wind
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine SET_MASS)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/08/94
!!      J.Stein     15/11/94  generate fields in geostrophic equilibrium in only
!!                            one step
!!      J.Stein     16/11/94  bug correction in the case LCORIO = false
!!      J.Stein     17/11/94  initialize w
!!      J.Stein     06/12/94  add an average for the coriolis parameter +
!!                            change the options for HFUNU and HFUNV
!!      J.Cuxart    Jan 23,1995 add the Boussinesq values for THVREF and
!!                             RHOREF
!!      J.Stein     16/03/95  remove R from the historical variables
!!      J.Stein     16/04/95  put the same names of the declarative modules
!!                            in the descriptive part
!!      J.Stein     27/09/95  compute the wind in  anel. balance
!!      J.Stein     30/01/96  use the RS ground pressure to initialize the
!!                            hydrostatic pressure computation
!!      P.Jabouille 28/11/97  bug in the U,V extrapolation
!!      J.Stein     10/07/97  add the other equation systems
!!      P.Jabouille 28/08/00  parallelization (to be improved latter)
!!      J.Stein     12/03/01  correct the Durran and MAE case
!!      P.Jabouille 03/12/02   add no thinshell condition
!!      J.Escobar   09/02/05  bug XTHVREFZ(JK-1<=>0)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_FM
USE MODE_ll
!
USE MODD_PARAMETERS  ! declarative modules
USE MODD_CONF
USE MODD_CST 
USE MODD_LUNIT_n
USE MODD_CONF
USE MODD_GRID
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_REF
USE MODD_FIELD_n
!
USE MODI_FUN       ! interface modules
USE MODI_SHUMAN
USE MODI_GATHER_ll
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments :
!
!
REAL, DIMENSION(:), INTENT(IN)    :: PUW,PVW ! Wind at w model grid levels
                                             ! (without orography)  : zonal and
                                             ! meridian components
REAL, DIMENSION(:), INTENT(IN)    :: PTHVM   ! Temperature at mass model grid
                                             ! levels (without orography)
REAL, DIMENSION(:), INTENT(IN)    :: PMRM    ! vapor mixing ratio at mass model
                                             ! grid levels (without orography)
REAL, DIMENSION(:), INTENT(IN)    :: PZHATM  ! Height of mass model grid levels
                                             ! (without orography)
REAL              , INTENT(IN)    :: PPGROUND! pressure at the ground heigth of
                                             ! the RS
REAL              , INTENT(IN)    :: PTHVGROUND! potential temperature at the
                                             ! ground height of the RS
REAL              , INTENT(IN)    :: PZGROUND! ground height of the RS
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
LOGICAL,                INTENT(IN)  :: OPV_PERT! logical switch for PV inversion
LOGICAL,                INTENT(IN)  :: ORMV_BL! logical switch for remouve boundary layer
!
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ  ! Coriolis parameter
                                                ! (exceptionally 3D array)
!
!
!*       0.2   Declarations of local variables :
!
INTEGER             :: IITRREF  !  iteration number for the determination of
                                 !  the thetav field
INTEGER             :: IRESP     !  FM return code
INTEGER             :: ILUOUT    ! Logical unit number for
                                             ! output-listing
REAL,DIMENSION(:,:),ALLOCATABLE :: ZGAMMA,ZGAMMA_ll ! K(lambda-lambda0) -beta
REAL                            :: ZGAMMALOC
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZTHV3D ! Potential
                                       ! virtual temperature on mass level grid
                                       ! (without orography)
                                       ! and deduced from thermal wind balance
REAL,DIMENSION(:,:),ALLOCATABLE   :: ZDTHVDXB,ZVJLOC,ZTHVJLOC,ZCORIOZJLOC,ZDXXJLOC
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZDTHVDYB
                                       ! x-gradient and y-gradient of thetav
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZTHV   ! potential
                                       ! virtual temperature on 3D model grid
                                       ! (mass levels, with orography)
REAL, DIMENSION(SIZE(XZHAT))                         :: ZEXNREFZ ! 1D anelastic
                                       ! reference pressure
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))  :: ZZM,ZZUM,ZZVM
                                       ! height at mass gridpoints, at u
                                       ! gridpoints  and at v gridpoints
                                       ! of grid with orogrpahy
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZDXX,ZDYY ! metric
                                       ! coefficients dxx,dyy for the grid
                                       ! without orography
REAL, DIMENSION(SIZE(XZHAT))                         :: ZDZHATM,ZDZHAT  ! deltaz
                                       ! for PZHATM and XZHAT levels
REAL :: ZGSCPD,ZRADSDG,         &      ! g/Cpd, Pi/180.,
        ZCPDSRD,ZRVSRD,ZRDSCPD         ! Cpd/Rd,Rv/Rd,Rd/Cpd
INTEGER :: IKB,IKE                 ! useful area in z direction
INTEGER :: IIU,IJU,IKU             ! Upper bounds in x,y,z directions
INTEGER :: IIU_ll,IJU_ll
INTEGER :: IXOR,IYOR,IILOC,IJLOC,IBEG,IEND
INTEGER :: IINFO_ll                ! return code of // routines
REAL                :: ZD1      ! DELTA1 (switch 0/1) for thinshell
                                ! approximation
REAL    :: ZDZ1SDH,ZDZ2SDH      ! working arrays
INTEGER :: JI,JJ,JK,JKS,JITR    ! Loop indexes
INTEGER :: IXOR_ll,IYOR_ll      ! origin's coordinates of extended subdomain
LOGICAL :: GGEOSBAL             ! logical to take into account the geostrophic
                                ! .TRUE.  => geostrophic balance
                                ! .FALSE. => unbalanced wind and mass fields
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZU3D,ZV3D ! Wind components
                                    !  in model axes at the w-point localization
REAL, DIMENSION(SIZE(XXHAT),1,1)         :: ZNFLX_TOT ! total normalized mass flux
REAL, DIMENSION(SIZE(XYHAT),SIZE(XZHAT)) ::  ZUYZ     ! vertical variations for U
REAL, DIMENSION(1,SIZE(XYHAT),1)         :: ZNFLY_TOT ! total normalized mass flux
REAL, DIMENSION(SIZE(XXHAT),SIZE(XZHAT)) ::  ZVXZ     ! vertical variations for V
REAL, DIMENSION(:,:), ALLOCATABLE        :: ZNFLX_TOT_ll,ZNFLY_TOT_ll,ZTHVREF2D
!
!-------------------------------------------------------------------------------
!
!*       1.     PROLOGUE : RETRIEVE LOGICAL UNIT NUMBER, INITIALIZE SOME
!                          CONSTANTS
!               ------------------------------------------------------------
!
CALL FMLOOK_ll(CLUOUT,CLUOUT,ILUOUT,IRESP)
!
ZRVSRD  = XRV / XRD
ZGSCPD   = XG / XCPD
ZRADSDG = XPI/180.
ZCPDSRD = XCPD/XRD
ZRDSCPD = XRD/XCPD
!
IIU = SIZE(XXHAT)
IJU = SIZE(XYHAT)
IKU = SIZE(XZHAT)
IKE = IKU - JPVEXT
IKB = 1 + JPVEXT
IIU_ll=NIMAX_ll+2*JPHEXT
IJU_ll=NJMAX_ll+2*JPHEXT
CALL GET_OR_ll('B',IXOR_ll,IYOR_ll)
!
!
!
!-------------------------------------------------------------------------------
!
!*       2.     ROTATE  WIND IN MODEL AXES, TAKE INTO ACCOUNT VARIATIONS
!               IN X,Y DIRECTION AND COMPUTE CORIOLIS PARAMETER, AND METRIC
!               COEFFICIENTS :
!               -------------------------------------------------------------
!
!*       2.1    Compute a first guess of the dry density
!
IF (LTHINSHELL .OR. LCARTESIAN) THEN
  ZD1=0.
ELSE
  ZD1=1.
END IF
!
! Exnrefz at IKB mass level
ZEXNREFZ(IKB)=((PPGROUND/XP00)**ZRDSCPD*(1.+ZD1*2./7.*(PZGROUND-PZHATM(IKB))/  &
                                     (XRADIUS+(PZHATM(IKB)+PZGROUND)/2.))  &
   + 2.*ZGSCPD/(PTHVM(IKB)+PTHVGROUND)*(PZGROUND-PZHATM(IKB)))/ &
(1.-ZD1*2./7.*(PZGROUND-PZHATM(IKB))/(XRADIUS+(PZHATM(IKB)+PZGROUND)/2.))
!
! integration of hydrostatic relation from IKB mass level to IKU mass level
DO JK =IKB+1,IKU
  ZEXNREFZ(JK)=(ZEXNREFZ(JK-1)*(1.+ZD1*2./7.*(PZHATM(JK-1)-PZHATM(JK))/  &
                                           (XRADIUS+XZHAT(JK)))+ &
     2.*ZGSCPD/(PTHVM(JK)+PTHVM(JK-1))*(PZHATM(JK-1)-PZHATM(JK)))/  &
  (1.-ZD1*2./7.*(PZHATM(JK-1)-PZHATM(JK))/(XRADIUS+XZHAT(JK)))
END DO
!
! integration of hydrostatic relation from IKB mass level to lowest mass level
DO JK = IKB-1, 1, -1
  ZEXNREFZ(JK)=(ZEXNREFZ(JK+1)*(1.+ZD1*2./7.*(PZHATM(JK+1)-PZHATM(JK))/  &
                                           (XRADIUS+XZHAT(JK+1)))+ &
     2.*ZGSCPD/(PTHVM(JK)+PTHVM(JK+1))*(PZHATM(JK+1)-PZHATM(JK)))/&
  (1.-ZD1*2./7.*(PZHATM(JK+1)-PZHATM(JK))/(XRADIUS+XZHAT(JK+1)))
END DO
!
! compute Exner function at model top
XEXNTOP=(ZEXNREFZ(IKE)*(1.+ZD1*2./7.*(PZHATM(IKE)-XZHAT(IKE+1))/  &
                                   (XRADIUS+(PZHATM(IKE)+XZHAT(IKE+1))/2.))  &
  + ZGSCPD/PTHVM(IKE)*(PZHATM(IKE)-XZHAT(IKE+1)))/ &
(1.-ZD1*2./7.*(PZHATM(IKE)-XZHAT(IKE+1))/(XRADIUS+(PZHATM(IKE)+XZHAT(IKE+1))/2.))
!
IF (OBOUSS) THEN
  XRHODREFZ(:) = XP00/ (XRD* PTHVM(:))
ELSE
  XRHODREFZ(:) = XP00* (ZEXNREFZ(:)** (ZCPDSRD -1.)) & ! rhodrefz is the density
           / (XRD* PTHVM(:)*(1.+PMRM(:)))           ! of dry air
ENDIF
!
!
!*       2.2    Compute the metric coefficients without orography and the Coriolis parameter
!
!
IF (LCARTESIAN) THEN                 ! cartesian geometry
  ZGAMMALOC = -XBETA *ZRADSDG
  IF (PRESENT(PCORIOZ)) THEN
    GGEOSBAL=.TRUE.
    PCORIOZ(:,:,:) =  2. * XOMEGA * SIN(XLAT0*ZRADSDG) ! f=2 omega sin(lat0)
  ELSE
    GGEOSBAL=.FALSE.
  END IF
  ZDXX(:,:,:) = MXM( SPREAD( SPREAD(XDXHAT(1:IIU),2,IJU),3,IKU) ) ! dxx (without orog.)
  ZDYY(:,:,:) = MYM( SPREAD( SPREAD(XDYHAT(1:IJU),1,IIU),3,IKU) ) ! dyy (without orog.)
ELSE                                 ! conformal projection
  ALLOCATE(ZGAMMA(IIU,IJU))
  ALLOCATE(ZGAMMA_ll(IIU_ll,IJU_ll))
  ZGAMMA(:,:) = XRPK * (XLON(:,:) -XLON0) * ZRADSDG -(XBETA *ZRADSDG)
  CALL GATHERALL_FIELD_ll('XY',ZGAMMA,ZGAMMA_ll,IINFO_ll)
  ZGAMMALOC = ZGAMMA_ll(KILOC,KJLOC)
  DEALLOCATE(ZGAMMA,ZGAMMA_ll)
!
  IF (PRESENT(PCORIOZ)) THEN
    GGEOSBAL=.TRUE.
    PCORIOZ(:,:,:)   = SPREAD(  2. * XOMEGA * SIN(XLAT(:,:)*ZRADSDG),3,IKU)
                                                    ! f=2 omega sin(lat)
  ELSE
    GGEOSBAL=.FALSE.
  END IF
  ZDXX(:,:,:) = MXM(                                                        &
             MZF(SPREAD(SPREAD( 1.+ ZD1*XZHAT(:)/XRADIUS ,1,IIU),2,IJU ))  &
                    * SPREAD( SPREAD(XDXHAT(1:IIU),2,IJU) /XMAP(:,:),3,IKU)   )
                                                   ! dxx (without orography)
  ZDYY(:,:,:) = MYM(                                                        &
              MZF(SPREAD(SPREAD( 1.+ ZD1*XZHAT(:)/XRADIUS,1,IIU),2,IJU ))  &
                    * SPREAD( SPREAD(XDYHAT(1:IJU),1,IIU) /XMAP(:,:),3,IKU)   )
                                                   ! dyy (without orography)
END IF
!
!*       2.3    Set the horizontal wind in approximate anelastic balance
!
!
ZDZHATM(2:IKU) = PZHATM(2:IKU)-PZHATM(1:IKU-1)
ZDZHATM(1) = ZDZHATM(2)
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
    ZUYZ(:,:)=FUNUYZ(IJU,IKU)
END SELECT
!
ZNFLX_TOT  = 0.
DO JK = 2,IKU-1
  DO JJ=2,IJU-1
    ZNFLX_TOT(:,1,1)=ZNFLX_TOT(:,1,1)+ZDYY(:,JJ,JK)*ZUYZ(JJ,JK)*ZDZHATM(JK)* &
                   XRHODREFZ(JK)
  END DO
END DO
!
ALLOCATE(ZNFLX_TOT_ll(IIU_ll,1))
CALL SUM_DIM1_ll(ZNFLX_TOT,ZNFLX_TOT_ll,IINFO_ll)
ZNFLX_TOT_ll=SIGN(1.,ZNFLX_TOT_ll)*MAX(ABS(ZNFLX_TOT_ll),1.E-40)
!
DO JI=1,IIU
  ZU3D(JI,:,:)=ZUYZ(:,:)*ZNFLX_TOT_ll(KILOC,1)/ZNFLX_TOT_ll(IXOR_ll-1+JI,1)
END DO
DEALLOCATE(ZNFLX_TOT_ll)
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
    ZVXZ(:,:)=FUNVXZ(IIU,IKU)
END SELECT
!
ZNFLY_TOT  = 0.
DO JK = 2,IKU-1
  DO JI=2,IIU-1
    ZNFLY_TOT(1,:,1) = ZNFLY_TOT(1,:,1) + ZDXX(JI,:,JK)*ZVXZ(JI,JK)* ZDZHATM(JK)* &
                   XRHODREFZ(JK)
  END DO
END DO
!
ALLOCATE(ZNFLY_TOT_ll(IJU_ll,1))
CALL SUM_DIM1_ll(ZNFLY_TOT,ZNFLY_TOT_ll,IINFO_ll)
ZNFLY_TOT_ll=SIGN(1.,ZNFLY_TOT_ll)*MAX(ABS(ZNFLY_TOT_ll),1.E-40)
!
!
DO JJ=1,IJU
    ZV3D(:,JJ,:)= ZVXZ(:,:) * ZNFLY_TOT_ll(KJLOC,1)/ZNFLY_TOT_ll(IYOR_ll-1+JJ,1)
END DO
DEALLOCATE(ZNFLY_TOT_ll)
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTE THETAV GRADIENT AND DEDUCE THETAV
!               -----------------------------------------
!
ZTHV3D(:,:,:) = SPREAD(SPREAD(PTHVM(:),1,IIU),2,IJU)  ! initialize  with
                                              ! (KILOC,KJLOC) vertical profile
IF (GGEOSBAL) THEN     ! thermal wind balance
!
  ZDZHATM(2:IKU) = PZHATM(2:IKU)-PZHATM(1:IKU-1)
  ZDZHAT(1:IKU-1) = XZHAT(2:IKU) - XZHAT(1:IKU-1)
!
  IF (OBOUSS) THEN
    IITRREF = 1              ! in this case, no iteration is useful  because
    XTHVREFZ(:) = PTHVM(2)   ! XTHVREFZ is uniform
  ELSE
    IITRREF = 2
  ENDIF
!
  CALL GET_OR_ll('B',IXOR,IYOR)
  IJLOC=KJLOC+1-IYOR
  IILOC=KILOC+1-IXOR
  ALLOCATE(ZDTHVDXB(IIU_ll,IKU))
  ALLOCATE(ZDXXJLOC(IIU_ll,IKU))
  ALLOCATE(ZCORIOZJLOC(IIU_ll,IKU))
  ALLOCATE(ZVJLOC(IIU_ll,IKU))
  ALLOCATE(ZTHVJLOC(IIU_ll,IKU))
  ALLOCATE(ZTHVREF2D(IIU,IKU))
  CALL GET_GLOBALSLICE_ll (ZDXX,'X',KJLOC,ZDXXJLOC,1,IIU_ll,1,IKU,IINFO_ll)
  CALL GET_GLOBALSLICE_ll (PCORIOZ,'X',KJLOC,ZCORIOZJLOC,1,IIU_ll,1,IKU,IINFO_ll)
  CALL GET_GLOBALSLICE_ll (ZV3D,'X',KJLOC,ZVJLOC,1,IIU_ll,1,IKU,IINFO_ll)
  DO JITR = 1, IITRREF
!
    CALL GET_GLOBALSLICE_ll (ZTHV3D,'X',KJLOC,ZTHVJLOC,1,IIU_ll,1,IKU,IINFO_ll)
    IF ( .NOT. OBOUSS ) THEN
!  compute the anelastic reference state when the geostrophic equilibrium is
!  taken into account
      XTHVREFZ(:)= SUM2D_ll(ZTHV3D,1,2,IINFO_ll,1,1,1,IIU_ll,IJU_ll,IKU) &
                                  /FLOAT(IIU_ll*IJU_ll)
    END IF
!
!*       3.1 Integration from I=ILOC to I=IIU
!
    IF ( CEQNSYS == 'LHE' ) THEN
      DO JI = KILOC , IIU_ll-1                             ! 1st step of integration
        ZDTHVDXB(JI,1:IKU-1)=XTHVREFZ(1:IKU-1)*0.25*             &
             (  ZVJLOC(JI,2:IKU)-ZVJLOC(JI,1:IKU-1)              &
               +ZVJLOC(JI+1,2:IKU)-ZVJLOC(JI+1,1:IKU-1)          &
             )/ ZDZHAT(1:IKU-1)*ZDXXJLOC(JI+1,1:IKU-1)*          &
             ( ZCORIOZJLOC(JI,1:IKU-1)+ZCORIOZJLOC(JI+1,1:IKU-1))/XG
        !
        ZDTHVDXB(JI,IKU) =  ZDTHVDXB(JI,IKU-1)
        !
        ZTHVJLOC(JI+1,:) = ZTHVJLOC(JI,:) + ZDTHVDXB(JI,:)
      END DO
    ELSEIF( CEQNSYS == 'MAE' ) THEN
      DO JI = KILOC , IIU_ll-1                             ! 1st step of integration
        ZDTHVDXB(JI,2:IKU-1)= 0.5  * XTHVREFZ(2:IKU-1)**2  *               &
             (  ZVJLOC(JI,3:IKU)    /(XTHVREFZ(2:IKU-1)+XTHVREFZ(3:IKU)  ) &
               -ZVJLOC(JI,2:IKU-1)  /(XTHVREFZ(1:IKU-2)+XTHVREFZ(2:IKU-1)) &
               +ZVJLOC(JI+1,3:IKU)  /(XTHVREFZ(2:IKU-1)+XTHVREFZ(3:IKU)  ) &
               -ZVJLOC(JI+1,2:IKU-1)/(XTHVREFZ(1:IKU-2)+XTHVREFZ(2:IKU-1)) &
             ) / ZDZHAT(2:IKU-1) * ZDXXJLOC(JI+1,2:IKU-1)*                 &
             ( ZCORIOZJLOC(JI,2:IKU-1)+ZCORIOZJLOC(JI+1,2:IKU-1) ) / XG
        !
        ZDTHVDXB(JI,1)   =  ZDTHVDXB(JI,2)
        ZDTHVDXB(JI,IKU) =  ZDTHVDXB(JI,IKU-1)
        !
        ZTHVJLOC(JI+1,:) = ZTHVJLOC(JI,:)+ZDTHVDXB(JI,:)
      END DO
    ELSEIF( CEQNSYS == 'DUR' ) THEN
      DO JI = KILOC , IIU_ll-1                             ! 1st step of integration
        ZDTHVDXB(JI,2:IKU-1) =                                              &
        ( ZVJLOC(JI,2:IKU-1)/ ( ZTHVJLOC(JI,2:IKU-1)+ZTHVJLOC(JI,1:IKU-2) ) &
         -ZVJLOC(JI,3:IKU)  / ( ZTHVJLOC(JI,3:IKU)  +ZTHVJLOC(JI,2:IKU-1) ) &
        ) / ZDZHAT(2:IKU-1)  * ZDXXJLOC(JI+1,2:IKU-1) *                     &
         ( ZCORIOZJLOC(JI,2:IKU-1) +ZCORIOZJLOC(JI+1,2:IKU-1) ) /XG
        !
        ZDTHVDXB(JI,1) =  ZDTHVDXB(JI,2)
        ZDTHVDXB(JI,IKU) =  ZDTHVDXB(JI,IKU-1)
        !
        ZTHVJLOC(JI+1,:) = 1./(1./ZTHVJLOC(JI,:)+ZDTHVDXB(JI,:))
      END DO
    END IF
!
    IF(IJLOC>=1 .AND. IJLOC<=IJU .AND. IILOC+1<=IIU) THEN
      IBEG=MAX(IILOC+1,1)
      ZTHV3D(IBEG:IIU,IJLOC,:) = ZTHVJLOC(IBEG+IXOR-1:IIU+IXOR-1,:)
    ENDIF
!
!*       3.2 Integration from I=ILOC to I=1
!
    IF ( CEQNSYS == 'LHE' ) THEN
      DO JI = KILOC, 2, -1                            ! 2nd step of integration
        ZDTHVDXB(JI,1:IKU-1) =  XTHVREFZ(1:IKU-1) * 0.25 *               &
             (  ZVJLOC(JI,2:IKU)   -ZVJLOC(JI,1:IKU-1)                   &
               +ZVJLOC(JI-1,2:IKU) -ZVJLOC(JI,1:IKU-1)                   &
             )/ ZDZHAT(1:IKU-1)  *  ZDXXJLOC(JI,1:IKU-1) *               &
             ( ZCORIOZJLOC(JI,1:IKU-1) + ZCORIOZJLOC(JI-1,1:IKU-1) ) /XG
        !
        ZDTHVDXB(JI,IKU) =  ZDTHVDXB(JI,IKU-1)
        !
        ZTHVJLOC(JI-1,:) = ZTHVJLOC(JI,:) - ZDTHVDXB(JI,:)
      END DO
    ELSEIF( CEQNSYS == 'MAE' ) THEN
      DO JI = KILOC, 2, -1                            ! 2nd step of integration
        ZDTHVDXB(JI,2:IKU-1) = 0.5  * XTHVREFZ(2:IKU-1)**2  *              &
             (  ZVJLOC(JI,3:IKU)    /(XTHVREFZ(2:IKU-1)+XTHVREFZ(3:IKU)  ) &
               -ZVJLOC(JI,2:IKU-1)  /(XTHVREFZ(1:IKU-2)+XTHVREFZ(2:IKU-1)) &
               +ZVJLOC(JI-1,3:IKU)  /(XTHVREFZ(2:IKU-1)+XTHVREFZ(3:IKU)  ) &
               -ZVJLOC(JI-1,2:IKU-1)/(XTHVREFZ(1:IKU-2)+XTHVREFZ(2:IKU-1)) &
             ) / ZDZHAT(2:IKU-1) * ZDXXJLOC(JI  ,2:IKU-1)*                 &
             ( ZCORIOZJLOC(JI,2:IKU-1)+ZCORIOZJLOC(JI-1,2:IKU-1) ) / XG
        !
        ZDTHVDXB(JI,1)   =  ZDTHVDXB(JI,2)
        ZDTHVDXB(JI,IKU) =  ZDTHVDXB(JI,IKU-1)
        !
        ZTHVJLOC(JI-1,:) = ZTHVJLOC(JI,:)-ZDTHVDXB(JI,:)
      END DO
    ELSEIF( CEQNSYS == 'DUR' ) THEN
      DO JI = KILOC, 2, -1                            ! 2nd step of integration
        ZDTHVDXB(JI,2:IKU-1) =                                                    &
        ( ZVJLOC(JI,2:IKU-1)  / ( ZTHVJLOC(JI,2:IKU-1)  +ZTHVJLOC(JI,1:IKU-2))     &
         -ZVJLOC(JI,3:IKU)    / ( ZTHVJLOC(JI,3:IKU)    +ZTHVJLOC(JI,2:IKU-1))     &
        ) / ZDZHAT(2:IKU-1)  * ZDXXJLOC(JI  ,2:IKU-1) *                           &
         ( ZCORIOZJLOC(JI,2:IKU-1) +ZCORIOZJLOC(JI-1,2:IKU-1) ) / XG
        !
        ZDTHVDXB(JI,1) =  ZDTHVDXB(JI,2)
        ZDTHVDXB(JI,IKU) =  ZDTHVDXB(JI,IKU-1)
        !
        ZTHVJLOC(JI-1,:) = 1./(1./ZTHVJLOC(JI,:)-ZDTHVDXB(JI,:))
        !
      END DO
    END IF
!
    IF(IJLOC>=1 .AND. IJLOC<=IJU .AND. IILOC>=2) THEN
      IEND=MIN(IILOC-1,IIU)
      ZTHV3D(1:IEND,IJLOC,:) = ZTHVJLOC(IXOR:IEND+IXOR-1,:)
    ENDIF
!
!
!*       3.3 Integration from J=JLOC to J=JU
!
    IF (.NOT. L2D ) THEN
      IF ( CEQNSYS == 'LHE' ) THEN
        DO JJ = IJLOC, IJU-1                 ! 3rd step of integration
          ZDTHVDYB(:,JJ,1:IKU-1) = 0.25 *                               &
                ( ZU3D(:,JJ,1:IKU-1)   -ZU3D(:,JJ,2:IKU)                &
                 +ZU3D(:,JJ+1,1:IKU-1) -ZU3D(:,JJ+1,2:IKU)              &
                ) * SPREAD(XTHVREFZ(1:IKU-1)/ZDZHAT(1:IKU-1),1,IIU)     &
                  * ZDYY(:,JJ+1,1:IKU-1)   *                            &
                ( PCORIOZ(:,JJ,1:IKU-1) + PCORIOZ(:,JJ+1,1:IKU-1) ) /XG
!
          ZDTHVDYB(:,JJ,IKU) =  ZDTHVDYB(:,JJ,IKU-1)
!
          ZTHV3D(:,JJ+1,:) = ZTHV3D(:,JJ,:) + ZDTHVDYB(:,JJ,:)
       END DO
      ELSEIF( CEQNSYS == 'MAE' ) THEN
        ZTHVREF2D=SPREAD(XTHVREFZ,1,IIU)
        DO JJ = IJLOC, IJU-1                 ! 3rd step of integration
          ZDTHVDYB(:,JJ,2:IKU-1) = -0.5 * ZTHVREF2D(:,2:IKU-1)**2          *        &
                ( ZU3D(:,JJ,3:IKU)    /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -ZU3D(:,JJ,2:IKU-1)  /(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
                 +ZU3D(:,JJ+1,3:IKU)  /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -ZU3D(:,JJ+1,2:IKU-1)/(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
                ) / SPREAD(ZDZHAT(2:IKU-1),1,IIU) * ZDYY(:,JJ+1,2:IKU-1)   *        &
                ( PCORIOZ(:,JJ,2:IKU-1) + PCORIOZ(:,JJ+1,2:IKU-1) ) /XG
!
          ZDTHVDYB(:,JJ,1)   =  ZDTHVDYB(:,JJ,2)
          ZDTHVDYB(:,JJ,IKU) =  ZDTHVDYB(:,JJ,IKU-1)
!
          ZTHV3D(:,JJ+1,:) = ZTHV3D(:,JJ,:) + ZDTHVDYB(:,JJ,:)
        END DO
      ELSEIF( CEQNSYS == 'DUR' ) THEN
        DO JJ = IJLOC, IJU-1                 ! 3rd step of integration
          ZDTHVDYB(:,JJ,2:IKU-1) =                                                    &
                ( ZU3D(:,JJ,3:IKU)     /(ZTHV3D(:,JJ,2:IKU-1)+ZTHV3D(:,JJ,3:IKU)    ) &
                 -ZU3D(:,JJ,2:IKU-1)   /(ZTHV3D(:,JJ,1:IKU-2)+ZTHV3D(:,JJ,2:IKU-1)  ) &
                ) / SPREAD(ZDZHAT(2:IKU-1),1,IIU) * ZDYY(:,JJ+1,2:IKU-1)   *          &
                ( PCORIOZ(:,JJ,2:IKU-1) + PCORIOZ(:,JJ+1,2:IKU-1) ) /XG

!
          ZDTHVDYB(:,JJ,1) =  ZDTHVDYB(:,JJ,2)
          ZDTHVDYB(:,JJ,IKU) =  ZDTHVDYB(:,JJ,IKU-1)
!
          ZTHV3D(:,JJ+1,:) = 1. /        &
          (1./ ZTHV3D(:,JJ,:) + ZDTHVDYB(:,JJ,:))
       END DO
     ENDIF
!
    ELSE
      DO JJ = KJLOC, IJU-1                 ! uniformity along the y direction
        ZTHV3D(:,JJ+1,:) = ZTHV3D(:,KJLOC,:)
      END DO
    END IF
!
!*       3.4 Integration from J=JLOC to J=1
!
    IF (.NOT. L2D ) THEN
      IF ( CEQNSYS == 'LHE' ) THEN
        DO JJ = KJLOC, 2, -1                       ! 4th step of integration
          ZDTHVDYB(:,JJ,1:IKU-1) = 0.25 *                              &
                ( ZU3D(:,JJ,1:IKU-1)   - ZU3D(:,JJ,2:IKU)              &
                 +ZU3D(:,JJ-1,1:IKU-1) - ZU3D(:,JJ-1,2:IKU)            &
                ) * SPREAD(XTHVREFZ(1:IKU-1)/ZDZHAT(1:IKU-1),1,IIU)    &
                  * ZDYY(:,JJ,1:IKU-1) *                               &
                ( PCORIOZ(:,JJ,1:IKU-1) + PCORIOZ(:,JJ-1,1:IKU-1) ) /XG
!
          ZDTHVDYB(:,JJ,IKU) =  ZDTHVDYB(:,JJ,IKU-1)
!
          ZTHV3D(:,JJ-1,:) = ZTHV3D(:,JJ,:) - ZDTHVDYB(:,JJ,:)
        END DO
      ELSEIF( CEQNSYS == 'MAE' ) THEN
        DO JJ = KJLOC, 2, -1                       ! 4th step of integration
          ZDTHVDYB(:,JJ,2:IKU-1) =  (-0.5) * ZTHVREF2D(:,2:IKU-1)**2          *     &
                ( ZU3D(:,JJ,3:IKU)    /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -ZU3D(:,JJ,2:IKU-1)  /(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
                 +ZU3D(:,JJ-1,3:IKU)  /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -ZU3D(:,JJ-1,2:IKU-1)/(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
                ) / SPREAD(ZDZHAT(2:IKU-1),1,IIU) * ZDYY(:,JJ  ,2:IKU-1)   *        &
                ( PCORIOZ(:,JJ,2:IKU-1) + PCORIOZ(:,JJ-1,2:IKU-1) ) /XG
!
          ZDTHVDYB(:,JJ,1)   =  ZDTHVDYB(:,JJ,2)
          ZDTHVDYB(:,JJ,IKU) =  ZDTHVDYB(:,JJ,IKU-1)
!
          ZTHV3D(:,JJ-1,:) = ZTHV3D(:,JJ,:) - ZDTHVDYB(:,JJ,:)
        END DO
      ELSEIF( CEQNSYS == 'DUR' ) THEN
        DO JJ = KJLOC, 2, -1                       ! 4th step of integration
          ZDTHVDYB(:,JJ,2:IKU-1) =                                                   &
                ( ZU3D(:,JJ,3:IKU)    /(ZTHV3D(:,JJ,2:IKU-1)+ZTHV3D(:,JJ,3:IKU)  )   &
                 -ZU3D(:,JJ,2:IKU-1)  /(ZTHV3D(:,JJ,1:IKU-2)+ZTHV3D(:,JJ,2:IKU-1))   &
                ) / SPREAD(ZDZHAT(2:IKU-1),1,IIU) * ZDYY(:,JJ  ,2:IKU-1)   *         &
                ( PCORIOZ(:,JJ,2:IKU-1) + PCORIOZ(:,JJ-1,2:IKU-1) ) /XG

!
          ZDTHVDYB(:,JJ,1) =  ZDTHVDYB(:,JJ,2)
          ZDTHVDYB(:,JJ,IKU) =  ZDTHVDYB(:,JJ,IKU-1)
!
          ZTHV3D(:,JJ-1,:) =  1. /        &
          (1./ ZTHV3D(:,JJ,:) - ZDTHVDYB(:,JJ,:))
        END DO
      END IF
    ELSE
      DO JJ = KJLOC-1, 1, -1                 ! uniformity along the y direction
        ZTHV3D(:,JJ,:) = ZTHV3D(:,KJLOC,:)
      END DO
    END IF
!
  END DO                         ! end of the iterative loop for ZTHV3D
ELSE
!  compute the anelastic reference state when the geostrophic equilibrium is
!  not taken into account
  IF (OBOUSS) THEN
    XTHVREFZ(:)=PTHVM(IKB)
  ELSE
    XTHVREFZ(:)=PTHVM(:)
  ENDIF
!
END IF                           ! end of the if block structure for Coriolis
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPLETE 1D ANELASTIC REFERENCE STATE
!               -------------------------------------
!
!
! Exnrefz at IKB mass level
ZEXNREFZ(IKB)=((PPGROUND/XP00)**ZRDSCPD*(1.+ZD1*2./7.*(PZGROUND-PZHATM(IKB))/  &
                                     (XRADIUS+(PZHATM(IKB)+PZGROUND)/2.))+  &
   2.*ZGSCPD/(PTHVM(IKB)+PTHVGROUND)*(PZGROUND-PZHATM(IKB)))/&
(1.-ZD1*2./7.*(PZGROUND-PZHATM(IKB))/(XRADIUS+(PZHATM(IKB)+PZGROUND)/2.))
!
! integration of hydrostatic relation from IKB mass level to IKU mass level
DO JK =IKB+1,IKU
  ZEXNREFZ(JK)=(ZEXNREFZ(JK-1)*(1.+ZD1*2./7.*(PZHATM(JK-1)-PZHATM(JK))/  &
                                           (XRADIUS+XZHAT(JK)))+ &
     2.*ZGSCPD/(XTHVREFZ(JK)+XTHVREFZ(JK-1))*(PZHATM(JK-1)-PZHATM(JK)))/  &
  (1.-ZD1*2./7.*(PZHATM(JK-1)-PZHATM(JK))/(XRADIUS+XZHAT(JK)))
END DO
!
! integration of hydrostatic relation from IKB mass level to lowest mass level
DO JK = IKB-1, 1, -1
  ZEXNREFZ(JK)=(ZEXNREFZ(JK+1)*(1.+ZD1*2./7.*(PZHATM(JK+1)-PZHATM(JK))/  &
                                           (XRADIUS+XZHAT(JK+1)))+ &
     2.*ZGSCPD/(XTHVREFZ(JK+1)+XTHVREFZ(JK))*(PZHATM(JK+1)-PZHATM(JK)))/&
  (1.-ZD1*2./7.*(PZHATM(JK+1)-PZHATM(JK))/(XRADIUS+XZHAT(JK+1)))
END DO
!
! compute Exner function at model top
XEXNTOP=(ZEXNREFZ(IKE)*(1.+ZD1*2./7.*(PZHATM(IKE)-XZHAT(IKE+1))/  &
                                   (XRADIUS+(PZHATM(IKE)+XZHAT(IKE+1))/2.))  &
  + ZGSCPD/XTHVREFZ(IKE)*(PZHATM(IKE)-XZHAT(IKE+1)))/ &
(1.-ZD1*2./7.*(PZHATM(IKE)-XZHAT(IKE+1))/(XRADIUS+(PZHATM(IKE)+XZHAT(IKE+1))/2.))
!
IF (OBOUSS) THEN
  XRHODREFZ(:) = XP00/ (XRD* XTHVREFZ(:))
ELSE
  XRHODREFZ(:) = XP00* (ZEXNREFZ(:)** (ZCPDSRD -1.)) & ! rhodrefz is the density
           / (XRD* XTHVREFZ(:)*(1.+PMRM(:)))           ! of dry air
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     INTERPOLATE THETAV, MR ON MODEL GRID (WITH OROGRAPHY)
!               ------------------------------------------------------------
!
ZZM(:,:,:)   = MZF(XZZ)                       ! compute height at mass level
                                              ! of grid with orography
!
ZZM(:,:,IKU) = 2. * XZZ(:,:,IKU) - ZZM(:,:,IKU-1) ! extrapolate on IKU mass level
!  ZZM(:,:,1) is always greater than or equal to PZHATM(1)
!  ZZM(:,:,IKU) is always smaller than or equal to PZHATM(IKU)
!
DO JI = 1,IIU
  DO JJ = 1,IJU
!
    DO JK = 1,IKU            ! loop on vertical levels of grid with orography
!
      IF (ZZM(JI,JJ,JK) >= PZHATM(IKU)) THEN     ! copy out when
        ZTHV(JI,JJ,JK)  =  ZTHV3D (JI,JJ,IKU)    ! ZZM(IKU)= PZHATM(IKU)
        XRM(JI,JJ,JK,1) = PMRM(IKU)              ! (in case zs=0.)
!
      ELSEIF (ZZM(JI,JJ,JK) < PZHATM(1)) THEN    ! copy out when  
        ZTHV(JI,JJ,JK)  =  ZTHV3D (JI,JJ,1)      ! ZZM(1)< PZHATM(1)  
        XRM(JI,JJ,JK,1) = PMRM(1)                ! (in case zs=0.)
!
      ELSE                  ! search levels on the mass grid without orography
        DO JKS = 2,IKU      ! that surrounded JK
          IF((ZZM(JI,JJ,JK) >= PZHATM(JKS-1)).AND.(ZZM(JI,JJ,JK) < PZHATM(JKS))) &
          THEN              ! interpolation with the values on the grid without
                            ! orography
            ZDZ1SDH = (ZZM(JI,JJ,JK)-PZHATM(JKS-1))/ (PZHATM(JKS)-PZHATM(JKS-1))
            ZDZ2SDH = 1. - ZDZ1SDH
            ZTHV(JI,JJ,JK)  = (ZDZ1SDH * ZTHV3D(JI,JJ,JKS)   )       &
                            + (ZDZ2SDH * ZTHV3D(JI,JJ,JKS-1) )
            XRM(JI,JJ,JK,1) = (ZDZ1SDH  * PMRM(JKS)         )        &
                            + (ZDZ2SDH  * PMRM(JKS-1)       )
          END IF
        END DO
      END IF
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       5.     DEDUCE THETA FROM THETAV AND MR ON MODEL GRID
!               ---------------------------------------------
!
XTHM(:,:,:) = ZTHV(:,:,:) * (1.+SUM(XRM(:,:,:,:),DIM=4)) /(1. + ZRVSRD *XRM(:,:,:,1))
!
!
!-------------------------------------------------------------------------------
!
!*       6.     COMPUTE (U,V,W) ON MODEL GRID
!               ----------------------------------
!
ZZUM(:,:,:) = MXM(ZZM)   ! Compute height at u grid-point
ZZVM(:,:,:) = MYM(ZZM)   ! compute height at v gridpoint
ZZUM(1,:,:) = 2. * ZZM(1,:,:) - ZZUM(2,:,:)
ZZVM(:,1,:) = 2. * ZZM(:,1,:) - ZZUM(:,2,:)
!
!
!*       6.1 U interpolation
!
DO JI = 2,IIU
  DO JJ =1, IJU
!
    DO JK = 1,IKU            ! loop on vertical levels of grid with orography
!
      IF (ZZUM(JI,JJ,JK) >= XZHAT(IKU)) THEN     ! extrapolation
        ZDZ1SDH = (ZZUM(JI,JJ,JK)-XZHAT(IKU))/ (XZHAT(IKU)-XZHAT(IKU-1))
        XUM(JI,JJ,JK) =  0.5*( ZU3D(JI,JJ,IKU) + ZU3D(JI-1,JJ,IKU) )  &
                         * (1.+ ZDZ1SDH)                              &
                        -0.5*(ZU3D(JI,JJ,IKU-1)+ ZU3D(JI-1,JJ,IKU-1)) &
                         * ZDZ1SDH
!
      ELSE IF (ZZUM(JI,JJ,JK) < XZHAT(1)) THEN     ! extrapolation
        ZDZ1SDH = (ZZUM(JI,JJ,JK)-XZHAT(1))/ (XZHAT(2)-XZHAT(1))
        XUM(JI,JJ,JK) = 0.5*( ZU3D(JI,JJ,1) + ZU3D(JI-1,JJ,1) )  &
                         * (1.- ZDZ1SDH)                         &
                        +0.5*(ZU3D(JI,JJ,2)+ ZU3D(JI-1,JJ,2) )   &
                         * ZDZ1SDH
!
      ELSE                  ! search levels on the mass grid without orography
        DO JKS = 2,IKU      ! that surrounded JK
          IF((ZZUM(JI,JJ,JK) >= XZHAT(JKS-1)).AND.(ZZUM(JI,JJ,JK) < XZHAT(JKS))) &
          THEN              ! interpolation with the values on the grid without
                          ! orography
            ZDZ1SDH = (ZZUM(JI,JJ,JK)-XZHAT(JKS-1))/ (XZHAT(JKS)-XZHAT(JKS-1))
            ZDZ2SDH = 1. - ZDZ1SDH
            XUM(JI,JJ,JK) = ZDZ1SDH *.5 *( ZU3D(JI,JJ,JKS)  + ZU3D(JI-1,JJ,JKS)) &
                           +ZDZ2SDH *.5 *( ZU3D(JI,JJ,JKS-1)+ ZU3D(JI-1,JJ,JKS-1))
          END IF
        END DO
      END IF
    END DO
  END DO
END DO
!
XUM(1,:,:)=-999.
!
!
!*       6.2 V interpolation
!
DO JI = 1,IIU
  DO JJ = 2,IJU
!
    DO JK = 1,IKU            ! loop on vertical levels of grid with orography
!
      IF (ZZVM(JI,JJ,JK) >= XZHAT(IKU)) THEN     ! extrapolation
        ZDZ1SDH = (ZZVM(JI,JJ,JK)-XZHAT(IKU))/ (XZHAT(IKU)-XZHAT(IKU-1))
        XVM(JI,JJ,JK) =   0.5*( ZV3D(JI,JJ,IKU) + ZV3D(JI,JJ-1,IKU) )  &
                         * (1.+ ZDZ1SDH)                               &
                        -0.5*(ZV3D(JI,JJ,IKU-1)+ ZV3D(JI,JJ-1,IKU-1))  &
                         * ZDZ1SDH
!
      ELSE IF (ZZVM(JI,JJ,JK) < XZHAT(1)) THEN     ! extrapolation
        ZDZ1SDH = (ZZVM(JI,JJ,JK)-XZHAT(1))/ (XZHAT(2)-XZHAT(1))
        XVM(JI,JJ,JK) =  0.5*( ZV3D(JI,JJ,1) + ZV3D(JI,JJ-1,1) )  &
                         * (1.- ZDZ1SDH)                         &
                        +0.5*(ZV3D(JI,JJ,2)+ ZV3D(JI,JJ-1,2) )   &
                         * ZDZ1SDH
!
      ELSE                  ! search levels on the mass grid without orography
        DO JKS = 2,IKU      ! that surrounded JK
          IF((ZZVM(JI,JJ,JK) >= XZHAT(JKS-1)).AND.(ZZVM(JI,JJ,JK) < XZHAT(JKS))) &
          THEN              ! interpolation with the values on the grid without
                          ! orography
            ZDZ1SDH = (ZZVM(JI,JJ,JK)-XZHAT(JKS-1))/ (XZHAT(JKS)-XZHAT(JKS-1))
            ZDZ2SDH = 1. - ZDZ1SDH
            XVM(JI,JJ,JK) = ZDZ1SDH *.5 *( ZV3D(JI,JJ,JKS)  + ZV3D(JI,JJ-1,JKS)) &
                           +ZDZ2SDH *.5 *( ZV3D(JI,JJ,JKS-1)+ ZV3D(JI,JJ-1,JKS-1))
          END IF
        END DO
      END IF
!
    END DO
  END DO
END DO
!
XVM(:,1,:)=-999.
!
!*       6.3 W initialization
!
XWM(:,:,:)=0.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_MASS
