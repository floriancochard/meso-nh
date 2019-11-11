! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_ANEL_BALANCE_n
!     ##########################
!
INTERFACE
!
SUBROUTINE ANEL_BALANCE_n(OINST,PRESIDUAL)
!
CHARACTER (LEN=1), INTENT(IN)  :: OINST      ! selected instant to enforce the
                                             ! anelastic constraint
!JUAN
REAL, OPTIONAL                 :: PRESIDUAL
!JUAN
END SUBROUTINE ANEL_BALANCE_n
!
END INTERFACE
!
END MODULE MODI_ANEL_BALANCE_n
!
!
!
!     ################################
      SUBROUTINE ANEL_BALANCE_n(OINST,PRESIDUAL)
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
!!    FMLOOK   : to retrieve a logical unit number associated with a file
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODE_ll
USE MODE_FM
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
USE MODD_LUNIT_n
!
USE MODI_TRID    ! interface modules
USE MODI_PRESSURE
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments :
!
!
CHARACTER (LEN=1), INTENT(IN)  :: OINST      ! selected instant to enforce the
                                             ! anelastic constraint
!JUAN
REAL, OPTIONAL                 :: PRESIDUAL
!JUAN
!
!*       0.2   Declarations of local variables :
!
INTEGER :: ILUOUT,IRESP           ! Logical unit number for output listing and
                                  ! return code
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
REAL, DIMENSION(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)) :: ZPABSM,ZPABST
                                               ! Potential at time t-dt  and t
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
!
!-------------------------------------------------------------------------------
!
!*       1.     PROLOGUE  :
!               --------
!
CALL FMLOOK_ll(CLUOUT0,CLUOUT0,ILUOUT,IRESP)
CALL GET_DIM_EXT_ll('Y',IIY,IJY)
IF (L2D) THEN
  ALLOCATE(ZBFY(IIY,IJY,SIZE(XRHODJ,3)))
ELSE
  ALLOCATE(ZBFY(IJY,IIY,SIZE(XRHODJ,3)))
ENDIF
ALLOCATE(ZTRIGSX(3*(NIMAX_ll+2*JPHEXT)))
ALLOCATE(ZTRIGSY(3*(NJMAX_ll+2*JPHEXT)))
!
!-------------------------------------------------------------------------------
!
!*       2.     PRESSURE SOLVER INITIALIZATION :
!               -------------------------------
!

!
CALL TRID(CLUOUT0,CLBCX,CLBCY,XMAP,XDXHAT,XDYHAT,ZDXHATM,ZDYHATM,ZRHOM,  &
          ZAF,ZCF,ZTRIGSX,ZTRIGSY,IIFAXX,IIFAXY,XRHODJ,XTHVREF,XZZ,ZBFY )
!
!-------------------------------------------------------------------------------
!
!*       3.      ANELASTIC CORRECTION :
!                ---------------------
!
!
!*       3.1     multiplication by RHODJ
!
IF (OINST == 'T') THEN
  ZRU(:,:,:) = MXM(XRHODJ) * XUT(:,:,:)
  ZRV(:,:,:) = MYM(XRHODJ) * XVT(:,:,:)
  ZRW(:,:,:) = MZM(XRHODJ) * XWT(:,:,:)
  ZTH(:,:,:) = XTHT(:,:,:)
  ALLOCATE(ZRR(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),SIZE(XRT,4)))
  ZRR(:,:,:,:) = XRT(:,:,:,:)
ELSEIF (OINST == 'M') THEN
  ZRU(:,:,:) = MXM(XRHODJ) * XUM(:,:,:)
  ZRV(:,:,:) = MYM(XRHODJ) * XVM(:,:,:)
  ZRW(:,:,:) = MZM(XRHODJ) * XWM(:,:,:)
  ZTH(:,:,:) = XTHM(:,:,:)
  ALLOCATE(ZRR(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),SIZE(XRM,4)))
  ZRR(:,:,:,:) = XRM(:,:,:,:)
ELSE
  STOP
END IF
!
!
!
!
!*       3.2     satisfy the anelastic constraint
!
ITCOUNT      =-1     ! no first guess of the pressure is available
ZPABSM(:,:,:)= 0.    !       ==================CAUTION=====================
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
CALL PRESSURE(GCLOSE_OUT,YFMFILE,CLUOUT,                             &
              CLBCX,CLBCY,CPRESOPT,NITR,LITRADJ,ITCOUNT,XRELAX,IMI,  &
              XRHODJ,XDXX,XDYY,XDZZ,XDZX,XDZY,ZDXHATM,ZDYHATM,ZRHOM, &
              ZAF,ZBFY,ZCF,ZTRIGSX,ZTRIGSY,IIFAXX,IIFAXY,ZPABSM,     &
              IRR,IRRL,IRRI,ZDRYMASST,ZREFMASS,ZMASS_O_PHI0,         &
              ZTH,ZRR,XRHODREF,XTHVREF,XRVREF,XEXNREF, XLINMASS,     &
              ZRU,ZRV,ZRW,ZPABST,PRESIDUAL                   )
!
DEALLOCATE(ZBFY,ZTRIGSX,ZTRIGSY,ZRR)
!*       3.2     return to the historical variables
!
IF (OINST == 'T') THEN
  XUT(:,:,:) = ZRU(:,:,:) / MXM(XRHODJ)
  XVT(:,:,:) = ZRV(:,:,:) / MYM(XRHODJ)
  XWT(:,:,:) = ZRW(:,:,:) / MZM(XRHODJ)
ELSEIF (OINST == 'M') THEN
  XUM(:,:,:) = ZRU(:,:,:) / MXM(XRHODJ)
  XVM(:,:,:) = ZRV(:,:,:) / MYM(XRHODJ)
  XWM(:,:,:) = ZRW(:,:,:) / MZM(XRHODJ)
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ANEL_BALANCE_n
