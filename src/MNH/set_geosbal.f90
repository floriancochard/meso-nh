!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_GEOSBAL
!     ####################
!
INTERFACE
!
SUBROUTINE SET_GEOSBAL(PU3D,PV3D,PTHVM,PMRM, &
                    KILOC,KJLOC,OBOUSS,PTHV,PCORIOZ)
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PU3D,PV3D ! Wind at w model grid levels
                                             ! (without orography)  : zonal and
                                             ! meridian components
REAL, DIMENSION(:), INTENT(IN)    :: PTHVM   ! Temperature at mass model grid
                                             ! levels (without orography)
REAL, DIMENSION(:), INTENT(IN)    :: PMRM    ! vapor mixing ratio at mass model
                                             ! grid levels (without orography)
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PCORIOZ  ! Coriolis parameter
                                                ! (exceptionally 3D array)
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PTHV  ! temperature at mass level on mesonh grid
!
END SUBROUTINE SET_GEOSBAL
!
END INTERFACE
!
END MODULE MODI_SET_GEOSBAL
!
!
!     ##########################################################################
      SUBROUTINE SET_GEOSBAL(PU3D,PV3D,PTHVM,PMRM, &
                    KILOC,KJLOC,OBOUSS,PTHV,PCORIOZ)
!     ##########################################################################
!
!!****  *SET_GEOSBAL * - routine to initialize mass and wind fields using thermal wind
!!                    balance
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize mass and wind fields on the 3D
!     model grid (with orography) from a vertical profile of (thetav,r) and wind
!     fields (at w-grid levels without ororgraphy i.e. at height XZHAT,
!     with u along latitude circles and v along meridians.)
!     The vertical profile of mass variable is given at mass-grid levels
!     (without orography), i.e. at height ZZHATM.
!       The thermal wind balance is used to compute the potential virtual
!     temperature from wind field and the vertical profile of thetav.
!       The vapor mixing ratio is taken  uniform in horizontal plane, i.e.
!     its vertical profile is only interpolated on the mass-grid levels with
!     orography, i.e. at height XZZ.
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
!!     - The orography is then taken into account :
!!         The potential virtual temperature, the vapor mixing ratio and the wind
!!      are interpolated on the model grid (with orography).
!!
!!         From this interpolated potential virtual temperature and vapor mixing
!!     ratio, the potential dry temperature is computed.
!!     So that, the prognostic mass variables (thetad and r) are completely
!!     defined.
!!
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
!!
!!      Module MODI_SHUMAN: interface for shuman operators
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
!!
!!      Module MODD_GRID  : contains grid variables
!!       XLON0,XLAT0 : Reference longitude and latitude
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
!!        XLON,XLAT     : Longitude and latitude
!!        XZHAT   : height level (without orography)
!!        XZZ     : height z
!!
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
!!     G.TANGUY        * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original : oct 2010 
!!    crée à partir de l'ancienne routine set_mass.f90 en prenant la partie
!!    concernant la balance geostrophique uniquement
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST 
USE MODD_DIM_n
USE MODD_FIELD_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_PARAMETERS
USE MODD_REF
!
USE MODE_GATHER_ll
USE MODE_ll
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments :
!
!
REAL, DIMENSION(:,:,:), INTENT(IN):: PU3D    ! Wind at w model grid levels
REAL, DIMENSION(:,:,:), INTENT(IN):: PV3D    ! (without orography)  : zonal and
                                             ! meridian components
REAL, DIMENSION(:), INTENT(IN)    :: PTHVM   ! vertical profil of Temperature at mass model grid
                                             ! levels (without orography)
REAL, DIMENSION(:), INTENT(IN)    :: PMRM    ! vertical profil of  vapor mixing ratio at mass model
                                             ! grid levels (without orography)
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PCORIOZ  ! Coriolis parameter
                                                ! (exceptionally 3D array)
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PTHV  !  potential virtual temperature 
                                             !  at mass level on mesonh grid                                                
!
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(XZHAT))    :: ZZHATM  ! Height of mass model grid levels
                                             ! (without orography)

INTEGER             :: IITRREF  !  iteration number for the determination of
                                 !  the thetav field
REAL,DIMENSION(:,:),ALLOCATABLE :: ZGAMMA,ZGAMMA_ll ! K(lambda-lambda0) -beta
REAL                            :: ZGAMMALOC
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZTHV3D ! Potential
                                       ! virtual temperature on mass level grid
                                       ! (without orography)
                                       ! and deduced from thermal wind balance
REAL,DIMENSION(:,:),ALLOCATABLE   :: ZDTHVDXB,ZVJLOC,ZTHVJLOC,ZCORIOZJLOC,ZDXXJLOC
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZDTHVDYB
                                       ! x-gradient and y-gradient of thetav
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT))  :: ZZM,ZZUM,ZZVM
                                       ! height at mass gridpoints, at u
                                       ! gridpoints  and at v gridpoints
                                       ! of grid with orogrpahy
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZDXX,ZDYY ! metric
                                       ! coefficients dxx,dyy for the grid
                                       ! without orography
REAL, DIMENSION(SIZE(XZHAT))                         :: ZDZHATM,ZDZHAT  ! deltaz
                                       ! for ZZHATM and XZHAT levels
REAL    :: ZRADSDG,ZRVSRD         !Pi/180,Rv/Rd
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
REAL, DIMENSION(:,:), ALLOCATABLE        ::ZTHVREF2D
!
!-------------------------------------------------------------------------------
!
!*       1.     PROLOGUE : RETRIEVE LOGICAL UNIT NUMBER, INITIALIZE SOME
!                          CONSTANTS
!               ------------------------------------------------------------
!
ZRVSRD  = XRV / XRD
ZRADSDG = XPI/180.
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
  ZZHATM(1:IKU-1) = 0.5 * (XZHAT(2:IKU)+XZHAT(1:IKU-1))  ! ZHATm(k)=
                                                     ! 0.5(ZHAT(k+1) +ZHAT(k)
  ZZHATM(IKU) = 2.* XZHAT(IKU) - ZZHATM(IKU-1)    ! extrapolation for IKU 
                                                  ! based on deltazhat(iku+1) =
                                                  ! deltazhat(iku) and  Zhatm(k)
                                                  ! is the middle point between 
                                                  ! Zhat(k) and Zhat(k+1)  
!                                                  
!*       2.1    Compute a first guess of the dry density
!
IF (LTHINSHELL .OR. LCARTESIAN) THEN
  ZD1=0.
ELSE
  ZD1=1.
END IF
!
!*       2.2    Compute the metric coefficients without orography and the Coriolis parameter
!
!
IF (LCARTESIAN) THEN                 ! cartesian geometry
  ZGAMMALOC = -XBETA *ZRADSDG
  PCORIOZ(:,:,:) =  2. * XOMEGA * SIN(XLAT0*ZRADSDG) ! f=2 omega sin(lat0)
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
  PCORIOZ(:,:,:)   = SPREAD(  2. * XOMEGA * SIN(XLAT(:,:)*ZRADSDG),3,IKU)
  ZDXX(:,:,:) = MXM(                                                        &
             MZF(1,IKU,1,SPREAD(SPREAD( 1.+ ZD1*XZHAT(:)/XRADIUS ,1,IIU),2,IJU ))  &
                    * SPREAD( SPREAD(XDXHAT(1:IIU),2,IJU) /XMAP(:,:),3,IKU)   )
                                                   ! dxx (without orography)
  ZDYY(:,:,:) = MYM(                                                        &
              MZF(1,IKU,1,SPREAD(SPREAD( 1.+ ZD1*XZHAT(:)/XRADIUS,1,IIU),2,IJU ))  &
                    * SPREAD( SPREAD(XDYHAT(1:IJU),1,IIU) /XMAP(:,:),3,IKU)   )
                                                   ! dyy (without orography)
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTE THETAV GRADIENT AND DEDUCE THETAV
!               -----------------------------------------
!
ZTHV3D(:,:,:) = SPREAD(SPREAD(PTHVM(:),1,IIU),2,IJU)  ! initialize  with
                                              ! (KILOC,KJLOC) vertical profile
  ZDZHATM(2:IKU) = ZZHATM(2:IKU)-ZZHATM(1:IKU-1)
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
  CALL GET_GLOBALSLICE_ll (PV3D,'X',KJLOC,ZVJLOC,1,IIU_ll,1,IKU,IINFO_ll)
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
                ( PU3D(:,JJ,1:IKU-1)   -PU3D(:,JJ,2:IKU)                &
                 +PU3D(:,JJ+1,1:IKU-1) -PU3D(:,JJ+1,2:IKU)              &
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
                ( PU3D(:,JJ,3:IKU)    /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -PU3D(:,JJ,2:IKU-1)  /(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
                 +PU3D(:,JJ+1,3:IKU)  /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -PU3D(:,JJ+1,2:IKU-1)/(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
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
                ( PU3D(:,JJ,3:IKU)     /(ZTHV3D(:,JJ,2:IKU-1)+ZTHV3D(:,JJ,3:IKU)    ) &
                 -PU3D(:,JJ,2:IKU-1)   /(ZTHV3D(:,JJ,1:IKU-2)+ZTHV3D(:,JJ,2:IKU-1)  ) &
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
                ( PU3D(:,JJ,1:IKU-1)   - PU3D(:,JJ,2:IKU)              &
                 +PU3D(:,JJ-1,1:IKU-1) - PU3D(:,JJ-1,2:IKU)            &
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
                ( PU3D(:,JJ,3:IKU)    /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -PU3D(:,JJ,2:IKU-1)  /(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
                 +PU3D(:,JJ-1,3:IKU)  /(ZTHVREF2D(:,2:IKU-1)+ZTHVREF2D(:,3:IKU)  )  &
                 -PU3D(:,JJ-1,2:IKU-1)/(ZTHVREF2D(:,1:IKU-2)+ZTHVREF2D(:,2:IKU-1))  &
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
                ( PU3D(:,JJ,3:IKU)    /(ZTHV3D(:,JJ,2:IKU-1)+ZTHV3D(:,JJ,3:IKU)  )   &
                 -PU3D(:,JJ,2:IKU-1)  /(ZTHV3D(:,JJ,1:IKU-2)+ZTHV3D(:,JJ,2:IKU-1))   &
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
!-------------------------------------------------------------------------------
!
!*       4.     INTERPOLATE THETAV, MR ON MODEL GRID (WITH OROGRAPHY)
!               ------------------------------------------------------------
!
ZZM(:,:,:)   = MZF(1,IKU,1,XZZ)                       ! compute height at mass level
                                              ! of grid with orography
!
ZZM(:,:,IKU) = 2. * XZZ(:,:,IKU) - ZZM(:,:,IKU-1) ! extrapolate on IKU mass level
!  ZZM(:,:,1) is always greater than or equal to ZZHATM(1)
!  ZZM(:,:,IKU) is always smaller than or equal to ZZHATM(IKU)
!
DO JI = 1,IIU
  DO JJ = 1,IJU
!
    DO JK = 1,IKU            ! loop on vertical levels of grid with orography
!
      IF (ZZM(JI,JJ,JK) >= ZZHATM(IKU)) THEN     ! copy out when
        PTHV(JI,JJ,JK)  =  ZTHV3D (JI,JJ,IKU)    ! ZZM(IKU)= ZZHATM(IKU)
        XRT(JI,JJ,JK,1) = PMRM(IKU)              ! (in case zs=0.)
!
      ELSEIF (ZZM(JI,JJ,JK) < ZZHATM(1)) THEN    ! copy out when  
        PTHV(JI,JJ,JK)  =  ZTHV3D (JI,JJ,1)      ! ZZM(1)< ZZHATM(1)  
        XRT(JI,JJ,JK,1) = PMRM(1)                ! (in case zs=0.)
!
      ELSE                  ! search levels on the mass grid without orography
        DO JKS = 2,IKU      ! that surrounded JK
          IF((ZZM(JI,JJ,JK) >= ZZHATM(JKS-1)).AND.(ZZM(JI,JJ,JK) < ZZHATM(JKS))) &
          THEN              ! interpolation with the values on the grid without
                            ! orography
            ZDZ1SDH = (ZZM(JI,JJ,JK)-ZZHATM(JKS-1))/ (ZZHATM(JKS)-ZZHATM(JKS-1))
            ZDZ2SDH = 1. - ZDZ1SDH
            PTHV(JI,JJ,JK)  = (ZDZ1SDH * ZTHV3D(JI,JJ,JKS)   )       &
                            + (ZDZ2SDH * ZTHV3D(JI,JJ,JKS-1) )
            XRT(JI,JJ,JK,1) = (ZDZ1SDH  * PMRM(JKS)         )        &
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
XTHT(:,:,:) = PTHV(:,:,:) * (1.+SUM(XRT(:,:,:,:),DIM=4)) /(1. + ZRVSRD *XRT(:,:,:,1))
!
!
!-------------------------------------------------------------------------------
!
!*       6.     INTERPOLATE (U,V,W) ON MODEL GRID (WITH OROGRAPHY)
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
        XUT(JI,JJ,JK) =  0.5*( PU3D(JI,JJ,IKU) + PU3D(JI-1,JJ,IKU) )  &
                         * (1.+ ZDZ1SDH)                              &
                        -0.5*(PU3D(JI,JJ,IKU-1)+ PU3D(JI-1,JJ,IKU-1)) &
                         * ZDZ1SDH
!
      ELSE IF (ZZUM(JI,JJ,JK) < XZHAT(1)) THEN     ! extrapolation
        ZDZ1SDH = (ZZUM(JI,JJ,JK)-XZHAT(1))/ (XZHAT(2)-XZHAT(1))
        XUT(JI,JJ,JK) = 0.5*( PU3D(JI,JJ,1) + PU3D(JI-1,JJ,1) )  &
                         * (1.- ZDZ1SDH)                         &
                        +0.5*(PU3D(JI,JJ,2)+ PU3D(JI-1,JJ,2) )   &
                         * ZDZ1SDH
!
      ELSE                  ! search levels on the mass grid without orography
        DO JKS = 2,IKU      ! that surrounded JK
          IF((ZZUM(JI,JJ,JK) >= XZHAT(JKS-1)).AND.(ZZUM(JI,JJ,JK) < XZHAT(JKS))) &
          THEN              ! interpolation with the values on the grid without
                          ! orography
            ZDZ1SDH = (ZZUM(JI,JJ,JK)-XZHAT(JKS-1))/ (XZHAT(JKS)-XZHAT(JKS-1))
            ZDZ2SDH = 1. - ZDZ1SDH
            XUT(JI,JJ,JK) = ZDZ1SDH *.5 *( PU3D(JI,JJ,JKS)  + PU3D(JI-1,JJ,JKS)) &
                           +ZDZ2SDH *.5 *( PU3D(JI,JJ,JKS-1)+ PU3D(JI-1,JJ,JKS-1))
          END IF
        END DO
      END IF
    END DO
  END DO
END DO
!
XUT(1,:,:)=-999.
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
        XVT(JI,JJ,JK) =   0.5*( PV3D(JI,JJ,IKU) + PV3D(JI,JJ-1,IKU) )  &
                         * (1.+ ZDZ1SDH)                               &
                        -0.5*(PV3D(JI,JJ,IKU-1)+ PV3D(JI,JJ-1,IKU-1))  &
                         * ZDZ1SDH
!
      ELSE IF (ZZVM(JI,JJ,JK) < XZHAT(1)) THEN     ! extrapolation
        ZDZ1SDH = (ZZVM(JI,JJ,JK)-XZHAT(1))/ (XZHAT(2)-XZHAT(1))
        XVT(JI,JJ,JK) =  0.5*( PV3D(JI,JJ,1) + PV3D(JI,JJ-1,1) )  &
                         * (1.- ZDZ1SDH)                         &
                        +0.5*(PV3D(JI,JJ,2)+ PV3D(JI,JJ-1,2) )   &
                         * ZDZ1SDH
!
      ELSE                  ! search levels on the mass grid without orography
        DO JKS = 2,IKU      ! that surrounded JK
          IF((ZZVM(JI,JJ,JK) >= XZHAT(JKS-1)).AND.(ZZVM(JI,JJ,JK) < XZHAT(JKS))) &
          THEN              ! interpolation with the values on the grid without
                          ! orography
            ZDZ1SDH = (ZZVM(JI,JJ,JK)-XZHAT(JKS-1))/ (XZHAT(JKS)-XZHAT(JKS-1))
            ZDZ2SDH = 1. - ZDZ1SDH
            XVT(JI,JJ,JK) = ZDZ1SDH *.5 *( PV3D(JI,JJ,JKS)  + PV3D(JI,JJ-1,JKS)) &
                           +ZDZ2SDH *.5 *( PV3D(JI,JJ,JKS-1)+ PV3D(JI,JJ-1,JKS-1))
          END IF
        END DO
      END IF
!
    END DO
  END DO
END DO
!
XVT(:,1,:)=-999.
!
!*       6.3 W initialization
!
XWT(:,:,:)=0.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_GEOSBAL
