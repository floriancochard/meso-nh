!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_FORCING
!     ###################
!
INTERFACE
!
      SUBROUTINE FORCING ( PTSTEP, OUSERV, PRHODJ, PCORIOZ,                &
                           PZHAT,  PZZ,  TPDTCUR,                          &
                           PUFRC_PAST, PVFRC_PAST,                         &
                           PUT,  PVT,  PWT,  PTHT,  PTKET,  PRT,  PSVT,    &
                           PRUS, PRVS, PRWS, PRTHS, PRTKES, PRRS, PRSVS,   &
                           KMI,PJ)
!
USE MODD_TIME, ONLY: DATE_TIME
!
REAL,                   INTENT(IN) :: PTSTEP  ! time-step
LOGICAL               , INTENT(IN) :: OUSERV  ! Logical to use rv
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! ( rhod J ) = dry density
              ! for reference state * Jacobian of the GCS transformation.
REAL, DIMENSION(:,:),   INTENT(IN) :: PCORIOZ ! f: Coriolis parameter
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT   ! height level without orography
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ     ! height z
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR ! current date and time
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PUFRC_PAST, PVFRC_PAST 
!                                             ! forcing at previous time-step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT,PVT,PWT,PTHT,PTKET 
                                          ! wind, potential temperature and
                                          ! TKE at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT !  moist variables at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT!  scalar variables at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS,PRVS,PRWS,PRTHS,PRTKES
                                          ! wind, potential temperature and
                                          ! TKE tendencies at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS ! moist variables at time t+1
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS! scalar variables at time t+1
!
INTEGER,                  INTENT(IN)    :: KMI      ! Model index
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PJ

!
END SUBROUTINE FORCING
!
END INTERFACE
!
END MODULE MODI_FORCING
!
!     ######################################################################
      SUBROUTINE FORCING ( PTSTEP, OUSERV, PRHODJ, PCORIOZ,                &
                           PZHAT,  PZZ,  TPDTCUR,                          &
                           PUFRC_PAST, PVFRC_PAST,                         &
                           PUT,  PVT,  PWT,  PTHT,  PTKET,  PRT,  PSVT,    &
                           PRUS, PRVS, PRWS, PRTHS, PRTKES, PRRS, PRSVS,   &
                           KMI,PJ)
!     ######################################################################
!
!!***  *FORCING* - routine to compute the forced terms 
!!
!!    PURPOSE
!!    -------
!!      The routine prepares (linear interpolations) and integrates each
!!    specified forcing terms. The forcing terms can be
!!      - a geostrophic wind (u_frc, v_frc)
!!      - a thermal wind     (Dtheta/Dx, Dtheta/Dy)
!!      - a tendency forcing (dth/dt, drv/dt )
!!      - a vertical transport (w_frc)
!!      - a newtonian relaxation (u_frc, v_frc, theta_frc, r_vfrc)
!!   
!!**  METHOD
!!    ------
!!       For its first call, the routine looks for a starting sounding with
!!    a date_and_time immediately lower or close to that the current
!!    date_and_time of the model. Then the temporal interpolation or extension
!!    is performed according to the position of the current date_and_time
!!    as compared to that of the forcing soundings. In case of non-flat
!!    terrain, vertical interpolations are necessary because the forcing terms
!!    are horizontally homogeneous.
!!      All the necessary interpolations are linear. Integration of each forcing
!!    term is enabled by a dedicated switch. The forced vertical motion is
!!    computed for each prognostic variable with an upstream scheme applied to
!!    the lagged fields. The forced advection of water vapor and the integration
!!    of the thermal and geostrophic wind is time centered. If a forced
!!    relaxation is enabled, a mask is defined to restrict the application of
!!    the forcing.
!!
!!    EXTERNAL
!!    --------
!!      Shuman functions       (finite differences operators)
!!      Upstream_z function    (upstream selection of the vertical advective
!!                              tendency)
!!      Temporal_lt function   (compare 2 TYPEd date_and_time data)
!!      Temporal_dist function (compute the number of seconds between
!!                              2 TYPEd date_and_time data)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONF: declaration of configuration variables
!!        LFLAT: tells if orography is present
!!        L1D  : tells if 1D model version is running
!!        NVERB: level of printed information on the output listing
!!      Module MODD_DYN: declaration of dynamic control variables
!!        LCORIO: Coriolis parameter flag
!!      Module MODD_FRC: declaration of the forcing variables
!!        NFRC  : number of forcing variables
!!        TDTFRC: date of each forcing profile
!!        XUFRC,XVFRC,XWFRC,XTHFRC,XRVFRC: large scale forcing variables
!!        XGXTHFRC,XGYTHFRC: large scale gradient of Theta
!!        XTENDTHFRC,XTENDRVFRC: large scale tendencies for Theta and Rv
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         TLUOUT0 : output-listing file
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPVEXT: define the number of marginal points out of the 
!!        physical domain along the vertical direction.    
!!      Module MODD_TIME: contains the structure of the TYPEd date_and_time
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	M. Georgelin      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/12/95 
!!      30/07/96 (Mari, Pinty, Suhre) restructuring
!!      18/11/96 J.-P. Pinty          remove $n
!!      10/12/96 J.-P. Pinty          reverse the loop order for vertical
!!                                    interpolation in case of orography
!!      24/01/97 J.-P. Pinty          add budget call
!!      24/01/98 P. Bechtold          use tendency forcing instead of rv_advect
!!                                    store large scale w
!!                                    add SST and surface pressure forcing
!!
!!      01/2004  V. Masson            surface externalization, removes SST
!!                                    forcing
!!      06/2012  V. Masson            Adds tendency of geostrophic wind itself to wind tendency
!!      01/2014  J. escobar           correction for // initialisation geostrophic ZUF,ZVF,ZWF 
!!      09/2017 Q.Rodier add LTEND_UV_FRC
!!      28/03/2018 P. Wautelet        Replace TEMPORAL_DIST by DATETIME_DISTANCE
!!                                    use overloaded comparison operator for date_time
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_DATETIME
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
!
USE MODD_CONF
USE MODD_DYN
USE MODD_FRC
USE MODD_LUNIT
USE MODD_PARAMETERS
USE MODD_TIME
USE MODD_BUDGET
USE MODD_CST
!
USE MODI_SHUMAN
USE MODI_UPSTREAM_Z
USE MODI_BUDGET
!
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                   INTENT(IN) :: PTSTEP  ! time-step
LOGICAL               , INTENT(IN) :: OUSERV  ! Logical to use rv
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! ( rhod J ) = dry density
              ! for reference state * Jacobian of the GCS transformation.
REAL, DIMENSION(:,:),   INTENT(IN) :: PCORIOZ ! f: Coriolis parameter
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT   ! height level without orography
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ     ! height z
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR ! current date and time
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PUFRC_PAST, PVFRC_PAST 
!                                             ! forcing at previous time-step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT,PVT,PWT,PTHT,PTKET 
                                          ! wind, potential temperature and
                                          ! TKE at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT !  moist variables at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT!  scalar variables at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS,PRVS,PRWS,PRTHS,PRTKES
                                          ! wind, potential temperature and
                                          ! TKE tendencies at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS ! moist variables at time t+1
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS! scalar variables at time t+1
!
INTEGER,                  INTENT(IN)    :: KMI      ! Model index
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PJ
!
!*       0.2   Declarations of local variables
!
INTEGER                         :: IIU, IJU, IKU      ! dimensions
INTEGER, SAVE                   :: JSX                ! saved loop index
INTEGER                         :: JI, JJ, JK, JL, JXP! loop indexes 
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWF, ZUF, ZVF  ! 3D forcing fields on
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTHF, ZRVF     !  the model grid mesh
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZGXTHF, ZGYTHF !          at
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTENDTHF, ZTENDRVF   !        time t
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTENDVF, ZTENDUF
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDUF, ZDVF     ! evolution of geostrophic wind
!                                                     ! during the time step
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCOEF          ! coefficient to take into
!                                                     ! account geostrophic wind evolution 
!                                                     ! in wind tendencies
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZDUT, ZDVT ! variation of Wind components
!                                                                    ! due to geostrophic wind evolution
!
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZZF, ZZA       ! altitudes
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDZZ, ZRWCF    ! dzz and ZWF contravar.
!
REAL, DIMENSION(SIZE(PUT,3)) :: ZXWFRC, ZXUFRC, ZXVFRC! 1D forcing fields
REAL, DIMENSION(SIZE(PUT,3)) :: ZXTHFRC, ZXRVFRC      !       after
REAL, DIMENSION(SIZE(PUT,3)) :: ZXGXTHFRC, ZXGYTHFRC  !       time
REAL, DIMENSION(SIZE(PUT,3)) :: ZXTENDTHFRC, ZXTENDRVFRC  !   interpolation 
REAL, DIMENSION(SIZE(PUT,3)) :: ZXTENDUFRC, ZXTENDVFRC
REAL                         :: ZXPGROUNDFRC          ! ground fields interpol.
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZOMEGA ! vertical velocity forcing (Pa/s)
!
LOGICAL, SAVE :: GSFIRSTCALL = .TRUE. ! control switch for the first call
!
REAL :: ZDZ, ZDT, ZALPHA ! height and time rate
REAL, SAVE :: ZSDTJX
REAL :: ZDZHAT_INV, ZDZHAT_INV_IKU ! Inverse of the vertical mesh size
!          
INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
INTEGER  :: IRESP   ! Return code of FM-routines
!
LOGICAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: GRELAX_MASK_FRC
!
!----------------------------------------------------------------------------
!
IIU=SIZE(PUT,1) 
IJU=SIZE(PUT,2) 
IKU=SIZE(PUT,3) 
!
ILUOUT0 = TLUOUT0%NLU
!
!*        1.   PREPARATION OF FORCING
!              ----------------------
!
IF (GSFIRSTCALL) THEN
!
  GSFIRSTCALL = .FALSE.
!
!*        1.1  printout number of forcing profiles
!
  WRITE(UNIT=ILUOUT0,FMT='(" THERE ARE ",I2," FORCING SOUNDINGS AT:")') NFRC
  DO JSX = 1 , NFRC
    WRITE(UNIT=ILUOUT0,FMT='(F9.0, "s, date:", I3, "/", I3, "/", I5)') &
      TDTFRC(JSX)%TIME,        &
      TDTFRC(JSX)%TDATE%DAY,   &
      TDTFRC(JSX)%TDATE%MONTH, &
      TDTFRC(JSX)%TDATE%YEAR
  END DO
!
!*        1.2  find first sounding to be used 
!
  JSX = 0
  IF( TPDTCUR < TDTFRC(1) ) THEN
    WRITE(UNIT=ILUOUT0,FMT='(" THE INITIAL FORCING FIELDS ARE NULL ")') 
    ELSE IF( TPDTCUR >= TDTFRC(NFRC) ) THEN
      WRITE(UNIT=ILUOUT0,FMT='(" THE FORCING FIELDS WILL REMAIN STATIONARY ")')
      ELSE
!
    TIM_FOR:    DO JI = NFRC-1, 1, -1
                  JSX = JI
                  IF( TPDTCUR >= TDTFRC(JSX) ) EXIT TIM_FOR
                END DO TIM_FOR
!
    WRITE(UNIT=ILUOUT0,FMT='(" THE INITIAL FORCING FIELDS ARE INTERPOLATED" , &
                       & " IN TIME STARTING FROM THE SOUNDING NUMBER ",I2)') JSX
    JSX = JSX - 1
  END IF
!
!*        1.3  printout all forcing field
!
  IF (NVERB >= 5) THEN
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "FORCING: the following forcing fields will be used:"
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XUFRC: geostrophic wind in X"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8F10.2))') &
           JK, (XUFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
           "XVFRC: geostrophic wind in Y"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8F10.2))') &
           JK, (XVFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XTHFRC: THETA for relaxation"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8F10.2))') &
           JK, (XTHFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XRVFRC: RV for relaxation"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8F10.6))') &
           JK, (XRVFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XWFRC: vertical transport velocity"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8F10.7))') &
           JK, (XWFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XGXTHFRC: thermal wind advection in X"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8E10.3))') &
           JK, (XGXTHFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XGYTHFRC: thermal wind advection in Y"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8E10.3))') &
           JK, (XGYTHFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XTENDTHFRC: Theta tendency"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8E10.3))') &
           JK, (XTENDTHFRC(JK,JL), JL=1, NFRC)       
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XTENDRVFRC: humidity tendency"
    DO JK = 1, IKU
      WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8E10.3))') &
           JK, (XTENDRVFRC(JK,JL), JL=1, NFRC)       
    END DO
!    
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
       "XTENDUFRC : wind advection tendency in X"
    DO JK = 1, IKU
     WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8E10.3))') &
         JK, (XTENDUFRC(JK,JL), JL=1, NFRC)
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
       "XTENDVFRC : wind advection tendency in Y"
    DO JK = 1, IKU
     WRITE(UNIT=ILUOUT0,FMT='(I10,99(/8E10.3))') &
         JK, (XTENDVFRC(JK,JL), JL=1, NFRC)
    END DO
!
    WRITE(UNIT=ILUOUT0,FMT='(A)') &
         "XPGROUNDFRC: SURF PRESSURE FORCING"
    WRITE(UNIT=ILUOUT0,FMT='(10X,99(/8E10.3))') &
         (XPGROUNDFRC(JL), JL=1, NFRC)       
!
  END IF
!
END IF
!
!*       2.     INTERPOLATE IN TIME
!    	        -------------------
!
IF( TPDTCUR < TDTFRC(1) ) THEN
  ZXUFRC(:)    = XUFRC(:,1)
  ZXVFRC(:)    = XVFRC(:,1)
  ZXWFRC(:)    = XWFRC(:,1)
  ZXTHFRC(:)   = XTHFRC(:,1)
  ZXRVFRC(:)   = XRVFRC(:,1)
  ZXTENDTHFRC(:) = XTENDTHFRC(:,1)
  ZXTENDRVFRC(:) = XTENDRVFRC(:,1)
  ZXGXTHFRC(:) = XGXTHFRC(:,1)
  ZXGYTHFRC(:) = XGYTHFRC(:,1)
  ZXTENDUFRC(:) = XTENDUFRC(:,1)
  ZXTENDVFRC(:) = XTENDVFRC(:,1)
  ZXPGROUNDFRC  = XPGROUNDFRC(1)
ELSE IF ( TPDTCUR >= TDTFRC(NFRC) ) THEN
  ZXUFRC(:)    = XUFRC(:,NFRC)
  ZXVFRC(:)    = XVFRC(:,NFRC)
  ZXWFRC(:)    = XWFRC(:,NFRC)
  ZXTHFRC(:)   = XTHFRC(:,NFRC)
  ZXRVFRC(:)   = XRVFRC(:,NFRC)
  ZXTENDTHFRC(:) = XTENDTHFRC(:,NFRC)
  ZXTENDRVFRC(:) = XTENDRVFRC(:,NFRC)
  ZXGXTHFRC(:) = XGXTHFRC(:,NFRC)
  ZXGYTHFRC(:) = XGYTHFRC(:,NFRC)
  ZXTENDUFRC(:) = XTENDUFRC(:,NFRC)
  ZXTENDVFRC(:) = XTENDVFRC(:,NFRC)
  ZXPGROUNDFRC  = XPGROUNDFRC(NFRC)
ELSE
  JXP = JSX + 1
  IF( TPDTCUR >= TDTFRC(JXP) ) THEN
    JSX = JSX +1
    JXP= JSX +1
    WRITE(UNIT=ILUOUT0,FMT='(" THE FORCING FIELDS ARE INTERPOLATED NOW" ,&
    & " BETWEEN SOUNDING NUMBER ",I2," AND SOUNDING NUMBER ",I2)') JSX,JXP
    CALL DATETIME_DISTANCE(TDTFRC(JSX),TDTFRC(JXP),ZSDTJX)
  END IF
!
  CALL DATETIME_DISTANCE(TDTFRC(JSX),TPDTCUR,ZDT)
!
  ZALPHA = ZDT / ZSDTJX
!
  ZXUFRC(:)    = XUFRC(:,JSX)   +(XUFRC(:,JXP)-XUFRC(:,JSX))*ZALPHA
  ZXVFRC(:)    = XVFRC(:,JSX)   +(XVFRC(:,JXP)-XVFRC(:,JSX))*ZALPHA
  ZXWFRC(:)    = XWFRC(:,JSX)   +(XWFRC(:,JXP)-XWFRC(:,JSX))*ZALPHA
  ZXTHFRC(:)   = XTHFRC(:,JSX)  +(XTHFRC(:,JXP)-XTHFRC(:,JSX))*ZALPHA
  ZXRVFRC(:)   = XRVFRC(:,JSX)  +(XRVFRC(:,JXP)-XRVFRC(:,JSX))*ZALPHA
  ZXTENDTHFRC(:) = XTENDTHFRC(:,JSX)+(XTENDTHFRC(:,JXP)-XTENDTHFRC(:,JSX))*ZALPHA
  ZXTENDRVFRC(:) = XTENDRVFRC(:,JSX)+(XTENDRVFRC(:,JXP)-XTENDRVFRC(:,JSX))*ZALPHA
  ZXTENDUFRC(:) = XTENDUFRC(:,JSX)+(XTENDUFRC(:,JXP)-XTENDUFRC(:,JSX))*ZALPHA
  ZXTENDVFRC(:) = XTENDVFRC(:,JSX)+(XTENDVFRC(:,JXP)-XTENDVFRC(:,JSX))*ZALPHA
  ZXGXTHFRC(:) = XGXTHFRC(:,JSX)+(XGXTHFRC(:,JXP)-XGXTHFRC(:,JSX))*ZALPHA
  ZXGYTHFRC(:) = XGYTHFRC(:,JSX)+(XGYTHFRC(:,JXP)-XGYTHFRC(:,JSX))*ZALPHA
  ZXPGROUNDFRC  = XPGROUNDFRC(JSX) +(XPGROUNDFRC(JXP)-XPGROUNDFRC(JSX))*ZALPHA
  !
END IF
!
!
!*       3.     INTERPOLATE IN SPACE
!   	        --------------------
!
ALLOCATE(ZUF(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)))
ALLOCATE(ZVF(SIZE(PVT,1),SIZE(PVT,2),SIZE(PVT,3)))
ALLOCATE(ZWF(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)))
ALLOCATE(ZTHF(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
ALLOCATE(ZRVF(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
ALLOCATE(ZGXTHF(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)))
ALLOCATE(ZGYTHF(SIZE(PVT,1),SIZE(PVT,2),SIZE(PVT,3)))
ALLOCATE(ZTENDTHF(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
ALLOCATE(ZTENDRVF(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)))
ALLOCATE(ZDUF(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)))
ALLOCATE(ZDVF(SIZE(PVT,1),SIZE(PVT,2),SIZE(PVT,3)))
ALLOCATE(ZTENDUF(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)))
ALLOCATE(ZTENDVF(SIZE(PVT,1),SIZE(PVT,2),SIZE(PVT,3)))
!
IF (LFLAT) THEN
!
!*       3.1    flat terrain case
!
  ZUF(:,:,:)    = SPREAD( SPREAD( ZXUFRC(:),1,IIU )  ,2,IJU )
  ZVF(:,:,:)    = SPREAD( SPREAD( ZXVFRC(:),1,IIU )  ,2,IJU )
  ZWF(:,:,:)    = SPREAD( SPREAD( ZXWFRC(:),1,IIU )  ,2,IJU )
  ZTHF(:,:,:)   = SPREAD( SPREAD( ZXTHFRC(:),1,IIU ) ,2,IJU )
  ZRVF(:,:,:)   = SPREAD( SPREAD( ZXRVFRC(:),1,IIU ) ,2,IJU )
  ZTENDTHF(:,:,:)  = SPREAD( SPREAD( ZXTENDTHFRC(:),1,IIU ),2,IJU )
  ZTENDRVF(:,:,:)  = SPREAD( SPREAD( ZXTENDRVFRC(:),1,IIU ),2,IJU )
  ZTENDUF(:,:,:)  = SPREAD( SPREAD( ZXTENDUFRC(:),1,IIU ),2,IJU )
  ZTENDVF(:,:,:)  = SPREAD( SPREAD( ZXTENDVFRC(:),1,IIU ),2,IJU )
  ZGXTHF(:,:,:) = SPREAD( SPREAD( ZXGXTHFRC(:),1,IIU ),2,IJU )
  ZGYTHF(:,:,:) = SPREAD( SPREAD( ZXGYTHFRC(:),1,IIU ),2,IJU )
ELSE
!
!*       3.2    case with orography
!
  ALLOCATE(ZZA(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)))
  ALLOCATE(ZZF(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)))
  ZZA(:,:,:)   = MZF(1,IKU,1, PZZ(:,:,:) )
  ZZA(:,:,IKU) = 2.0*PZZ(:,:,IKU) - ZZA(:,:,IKU-1)
  ZDZHAT_INV_IKU = 1.0 / ( PZHAT(IKU)-PZHAT(IKU-1) )
!
  ZZF(:,:,:)   = MXM( ZZA(:,:,:) )
  ZZF(1,:,:)   = 2.0*ZZA(1,:,:) - ZZF(2,:,:)
!
  ZUF = 0.0
!
  DO JL=1,IKU-1
    ZDZHAT_INV = 1.0 / ( PZHAT(JL+1)-PZHAT(JL) )
    DO JK=1,IKU
      DO JJ=1,IJU
        DO JI=1,IIU
          IF( ZZF(JI,JJ,JK) >= PZHAT(JL).AND.ZZF(JI,JJ,JK) <= PZHAT(JL+1) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(JL))  * ZDZHAT_INV
            ZUF(JI,JJ,JK) = ZXUFRC(JL+1)*ZDZ + ZXUFRC(JL)*(1-ZDZ)
          ELSE IF( ZZF(JI,JJ,JK) > PZHAT(IKU) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(IKU)) * ZDZHAT_INV_IKU
            ZUF(JI,JJ,JK) = ZXUFRC(IKU)*ZDZ  + ZXUFRC(IKU-1)*(1-ZDZ)
          END IF
        END DO
      END DO
    END DO
  END DO
  CALL GET_HALO(ZUF)
!
  ZZF(:,:,:) = MYM( ZZA(:,:,:) )
  ZZF(:,1,:) = 2.0*ZZA(:,1,:)-ZZF(:,2,:)
!
  ZVF = 0.0
!
  DO JL=1,IKU-1
    ZDZHAT_INV = 1.0 / ( PZHAT(JL+1)-PZHAT(JL) )
    DO JK=1,IKU
      DO JJ=1,IJU
        DO JI=1,IIU
          IF( ZZF(JI,JJ,JK) >= PZHAT(JL).AND.ZZF(JI,JJ,JK) <= PZHAT(JL+1) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(JL))  * ZDZHAT_INV
            ZVF(JI,JJ,JK) = ZXVFRC(JL+1)*ZDZ + ZXVFRC(JL)*(1-ZDZ)
          ELSE IF( ZZF(JI,JJ,JK) > PZHAT(IKU) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(IKU)) * ZDZHAT_INV_IKU
            ZVF(JI,JJ,JK) = ZXVFRC(IKU)*ZDZ  + ZXVFRC(IKU-1)*(1-ZDZ)
          END IF
        END DO
      END DO
    END DO
  END DO
  CALL GET_HALO(ZVF)
!
  ZWF = 0.0
!
  DO JL=1,IKU-1
    ZDZHAT_INV = 1.0 / ( PZHAT(JL+1)-PZHAT(JL) )
    DO JK=1,IKU
      DO JJ=1,IJU
        DO JI=1,IIU
          IF( ZZF(JI,JJ,JK) >= PZHAT(JL).AND.ZZF(JI,JJ,JK) <= PZHAT(JL+1) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(JL))  * ZDZHAT_INV
            ZWF(JI,JJ,JK) = ZXWFRC(JL+1)*ZDZ + ZXWFRC(JL)*(1-ZDZ)
          ELSE IF( ZZF(JI,JJ,JK) > PZHAT(IKU) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(IKU)) * ZDZHAT_INV_IKU
            ZWF(JI,JJ,JK) = ZXWFRC(IKU)*ZDZ  + ZXWFRC(IKU-1)*(1-ZDZ)
          END IF
        END DO
      END DO
    END DO
  END DO
  CALL GET_HALO(ZWF)
!
  ZZF(:,:,:)   = MZF(1,IKU,1, PZZ(:,:,:) )
  ZZF(:,:,IKU) = 2.0*PZZ(:,:,IKU)-ZZF(:,:,IKU-1)
!
  DO JL=1,IKU-1
    ZDZHAT_INV = 1.0 / ( PZHAT(JL+1)-PZHAT(JL) )
    DO JK=1,IKU
      DO JJ=1,IJU
        DO JI=1,IIU
          IF( ZZF(JI,JJ,JK) >= PZHAT(JL).AND.ZZF(JI,JJ,JK) <= PZHAT(JL+1) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(JL))  * ZDZHAT_INV
            ZTHF(JI,JJ,JK)   = ZXTHFRC(JL+1)*ZDZ   + ZXTHFRC(JL)*(1-ZDZ)
            ZRVF(JI,JJ,JK)   = ZXRVFRC(JL+1)*ZDZ   + ZXRVFRC(JL)*(1-ZDZ)
            ZGXTHF(JI,JJ,JK) = ZXGXTHFRC(JL+1)*ZDZ + ZXGXTHFRC(JL)*(1-ZDZ)
            ZGYTHF(JI,JJ,JK) = ZXGYTHFRC(JL+1)*ZDZ + ZXGYTHFRC(JL)*(1-ZDZ)
            ZTENDTHF(JI,JJ,JK)  = ZXTENDTHFRC(JL+1)*ZDZ + ZXTENDTHFRC(JL)*(1-ZDZ)
            ZTENDRVF(JI,JJ,JK)  = ZXTENDRVFRC(JL+1)*ZDZ + ZXTENDRVFRC(JL)*(1-ZDZ)
            ZTENDUF(JI,JJ,JK)  = ZXTENDUFRC(JL+1)*ZDZ + ZXTENDUFRC(JL)*(1-ZDZ)
            ZTENDVF(JI,JJ,JK)  = ZXTENDVFRC(JL+1)*ZDZ + ZXTENDVFRC(JL)*(1-ZDZ)
          ELSE IF( ZZF(JI,JJ,JK) > PZHAT(IKU) ) THEN
            ZDZ = (ZZF(JI,JJ,JK)-PZHAT(IKU)) * ZDZHAT_INV_IKU
            ZTHF(JI,JJ,JK)   = ZXTHFRC(IKU)*ZDZ   + ZXTHFRC(IKU-1)*(1-ZDZ)
            ZRVF(JI,JJ,JK)   = ZXRVFRC(IKU)*ZDZ   + ZXRVFRC(IKU-1)*(1-ZDZ)
            ZGXTHF(JI,JJ,JK) = ZXGXTHFRC(IKU)*ZDZ + ZXGXTHFRC(IKU-1)*(1-ZDZ)
            ZGYTHF(JI,JJ,JK) = ZXGYTHFRC(IKU)*ZDZ + ZXGYTHFRC(IKU-1)*(1-ZDZ)
            ZTENDTHF(JI,JJ,JK)  = ZXTENDTHFRC(IKU)*ZDZ + ZXTENDTHFRC(IKU-1)*(1-ZDZ)
            ZTENDRVF(JI,JJ,JK)  = ZXTENDRVFRC(IKU)*ZDZ + ZXTENDRVFRC(IKU-1)*(1-ZDZ)
            ZTENDUF(JI,JJ,JK)  = ZXTENDUFRC(IKU)*ZDZ + ZXTENDUFRC(IKU-1)*(1-ZDZ)
            ZTENDVF(JI,JJ,JK)  = ZXTENDVFRC(IKU)*ZDZ + ZXTENDVFRC(IKU-1)*(1-ZDZ)                  
          END IF
        END DO
      END DO
    END DO
  END DO
END IF
!
!!============================
!!
!! Ligne to add if you want W in Pa/s in namelist instead of m/s (omega = - w/(rho*g))
!!
!ZWF(:,:,:) = - ZWF(:,:,:)/(XG*MZM(1,IKU,1,(PRHODJ(:,:,:)/PJ(:,:,:))))
!
!!============================
!
!
!
! under the ground, forcings do not exist.
!
DO JK=1,JPVEXT
  ZUF(:,:,JK)    = 0.
  ZVF(:,:,JK)    = 0.
  ZWF(:,:,JK)    = 0.
  ZTHF(:,:,JK)   = 0.
  ZRVF(:,:,JK)   = 0.
  ZGXTHF(:,:,JK) = 0.
  ZGYTHF(:,:,JK) = 0.
  ZTENDTHF(:,:,JK)  = 0.
  ZTENDRVF(:,:,JK)  = 0.
  ZTENDUF(:,:,JK)  = 0.
  ZTENDVF(:,:,JK)  = 0.
END DO
!
! store large scale w in module to be used later
! in convection scheme
XWTFRC(:,:,:) = ZWF(:,:,:) 
!
!* computes evolution of forcing wind
WHERE(PUFRC_PAST==XUNDEF) PUFRC_PAST = ZUF(:,:,:)
WHERE(PVFRC_PAST==XUNDEF) PVFRC_PAST = ZVF(:,:,:)
!
ZDUF(:,:,:) = ZUF(:,:,:) - PUFRC_PAST(:,:,:)
ZDVF(:,:,:) = ZVF(:,:,:) - PVFRC_PAST(:,:,:)
!
!
IF (.NOT.LFLAT) THEN
!
  DEALLOCATE(ZZA)
  DEALLOCATE(ZZF)
!
END IF
!
!
!*       4.     INTEGRATION OF THE FORCINGS IN THE SOURCES
!   	        ------------------------------------------
!
ALLOCATE(ZDZZ(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)))
ALLOCATE(ZRWCF(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)))
!
!*       4.1    integration of vertical motion (upstream scheme)
!
IF (LVERT_MOTION_FRC) THEN
  ZDZZ(:,:,:) = DZM(1,IKU,1,MZF(1,IKU,1,PZZ(:,:,:)))
  ZDZZ(:,:,IKU) = PZZ(:,:,IKU) - PZZ(:,:,IKU-1) ! same delta z in IKU and IKU -1
!
  ZRWCF(:,:,:) =   ZWF(:,:,:) * MZM(1,IKU,1,PRHODJ(:,:,:)) / ZDZZ(:,:,:)
  ZRWCF(:,:,1) = - ZRWCF(:,:,3)     ! Mirror hypothesis
!
! forced vertical transport of U and V
!
  ZDZZ(:,:,:) = MXF(ZRWCF(:,:,:)) *DZM(1,IKU,1,PUT(:,:,:))
  PRUS(:,:,:) = PRUS(:,:,:) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
  ZDZZ(:,:,:) = MYF(ZRWCF(:,:,:)) *DZM(1,IKU,1,PVT(:,:,:))
  PRVS(:,:,:) = PRVS(:,:,:) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
!
! forced vertical transport of W
!
  IF( .NOT.L1D ) THEN
    ZDZZ(:,:,:) = MZF(1,IKU,1,ZRWCF(:,:,:)) *DZF(1,IKU,1,PWT(:,:,:))
    PRWS(:,:,:) = PRWS(:,:,:) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
  END IF
!
! forced vertical transport of THETA
!
  ZDZZ(:,:,:) = ZRWCF(:,:,:)  *DZM(1,IKU,1,PTHT(:,:,:))
  PRTHS(:,:,:) = PRTHS(:,:,:) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
!
! forced vertical transport of TKE (if allocated)
!
  IF( SIZE(PTKET) == SIZE(ZDZZ) ) THEN
    ZDZZ(:,:,:) = ZRWCF(:,:,:)  *DZM(1,IKU,1,PTKET(:,:,:))
    PRTKES(:,:,:) = PRTKES(:,:,:) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
  END IF
!
! forced vertical transport of water variables
!
  DO JL = 1 , SIZE(PRRS,4)
    ZDZZ(:,:,:) = ZRWCF(:,:,:) *DZM(1,IKU,1,PRT(:,:,:,JL))
    PRRS(:,:,:,JL) = PRRS(:,:,:,JL) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
  END DO
!
! forced vertical transport of scalar variables
!
  DO JL = 1 , SIZE(PRSVS,4)
    ZDZZ(:,:,:) = ZRWCF(:,:,:) *DZM(1,IKU,1,PSVT(:,:,:,JL))
    PRSVS(:,:,:,JL) = PRSVS(:,:,:,JL) - UPSTREAM_Z(ZDZZ(:,:,:),ZRWCF(:,:,:))
  END DO
!
END IF
!
!*       4.2    integration of the tendency forcing for th and rv
!
IF ( LTEND_THRV_FRC ) THEN
  PRTHS(:,:,:) = PRTHS(:,:,:) + PRHODJ(:,:,:) * ZTENDTHF(:,:,:)
  IF ( OUSERV ) THEN
      PRRS(:,:,:,1) = PRRS(:,:,:,1) + PRHODJ(:,:,:) * ZTENDRVF(:,:,:)
  END IF
END IF
!
!*       4.2.1    integration of the tendency forcing for uv
!
IF ( LTEND_UV_FRC ) THEN
 PRUS(:,:,:) = PRUS(:,:,:) +  MXM(PRHODJ) * ZTENDUF(:,:,:)
 PRVS(:,:,:) = PRVS(:,:,:) +  MYM(PRHODJ) * ZTENDVF(:,:,:)
END IF
!
!*       4.3    integration of the thermal and geostrophic wind
!
IF( LCORIO ) THEN
!
! thermal wind advection
!
  IF (LGEOST_TH_FRC) THEN
    PRTHS(:,:,:) = PRTHS(:,:,:) - PRHODJ(:,:,:)*(MXF(PUT(:,:,:))*ZGXTHF(:,:,:) &
                                               + MYF(PVT(:,:,:))*ZGYTHF(:,:,:) )
  END IF
!
! geostrophic wind for U and V due to large scale pressure gradients
!
  IF (LGEOST_UV_FRC) THEN
    ! Adds pressure force (in the form of geostrophic wind) to U component
    PRUS(:,:,:) = PRUS(:,:,:)                                                    &
                - MXM( MYF(ZVF(:,:,:))*PRHODJ(:,:,:)*SPREAD(PCORIOZ(:,:),3,IKU))
    ! Adds pressure force (in the form of geostrophic wind) to V component
    PRVS(:,:,:) = PRVS(:,:,:)                                                    &
                + MYM( MXF(ZUF(:,:,:))*PRHODJ(:,:,:)*SPREAD(PCORIOZ(:,:),3,IKU))
    ! adds tendency of geostrophic wind to force wind in the free troposphere to
    ! follow the geostrophic wind when the latter changes. 
    ! When winds differs from the geotrophic wind, the impact of this tendency is reduced.
    IF ( .NOT. LTEND_UV_FRC ) THEN
      ALLOCATE(ZCOEF(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)))
      ZCOEF(:,:,:) = (MXF(PUT       **2)+MYF(PVT       **2))        &
                 /MAX(MXF(PUFRC_PAST**2)+MYF(PVFRC_PAST**2), 1.E-3)
      !
      ZCOEF(:,:,:) = MIN(1.,SQRT(ZCOEF))
      !
      ZDUT(:,:,:) = ZDUF(:,:,:) * MXM(ZCOEF) 
      ZDVT(:,:,:) = ZDVF(:,:,:) * MYM(ZCOEF) 
      !
      PRUS(:,:,:) = PRUS(:,:,:) + ZDUT(:,:,:) * MXM(PRHODJ) / PTSTEP
      !
      PRVS(:,:,:) = PRVS(:,:,:) + ZDVT(:,:,:) * MYM(PRHODJ) / PTSTEP
      !
      !
      ! Takes into acount the Coriolis force due to this evolution
      PRUS(:,:,:) = PRUS(:,:,:)                                                    &
                  + MXM( MYF(ZDVT(:,:,:))*PRHODJ(:,:,:)*SPREAD(PCORIOZ(:,:),3,IKU))
      PRVS(:,:,:) = PRVS(:,:,:)                                                    &
                  - MYM( MXF(ZDUT(:,:,:))*PRHODJ(:,:,:)*SPREAD(PCORIOZ(:,:),3,IKU))
    !
      DEALLOCATE(ZCOEF)
    END IF
  END IF
!
END IF
!
! stores new forcing wind
PUFRC_PAST(:,:,:) = ZUF(:,:,:)
PVFRC_PAST(:,:,:) = ZVF(:,:,:)
!
!
!*       4.4    integration of the thermal, moisture and wind relaxation
!
IF( LRELAX_THRV_FRC .OR. LRELAX_UV_FRC ) THEN
!
  ZDZZ(:,:,:) = DZM(1,IKU,1,MZF(1,IKU,1,PZZ(:,:,:)))
  ZDZZ(:,:,IKU) = PZZ(:,:,IKU) - PZZ(:,:,IKU-1)
!
! define the mask where the relaxation is to be applied
!
  SELECT CASE (CRELAX_HEIGHT_TYPE)
    CASE ("FIXE")
      GRELAX_MASK_FRC(:,:,:) = .TRUE.
    CASE ("THGR")
      CALL DEFINE_RELAX_FORCING(GRELAX_MASK_FRC,PTHT,ZDZZ,JPHEXT,JPVEXT)
    CASE DEFAULT
      ! the following error should not occur, since tests are made earlier
!callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','FORCING','wrong CRELAX_HEIGHT_TYPE option')
  END SELECT
  WHERE ( MZF(1,IKU,1,PZZ(:,:,:)) .LE. XRELAX_HEIGHT_FRC )
    GRELAX_MASK_FRC = .FALSE.
  END WHERE
!
  IF( LRELAX_THRV_FRC ) THEN
!
!  apply THETA relaxation
!
    WHERE( GRELAX_MASK_FRC )
      PRTHS(:,:,:) = PRTHS(:,:,:) - PRHODJ(:,:,:)*(PTHT(:,:,:)-ZTHF(:,:,:)) &
                                                 / XRELAX_TIME_FRC
    END WHERE
!
!   apply humidity relaxation
!
    IF( OUSERV ) THEN
      WHERE( GRELAX_MASK_FRC )
        PRRS(:,:,:,1) = PRRS(:,:,:,1) &
		      - PRHODJ(:,:,:)*(PRT(:,:,:,1)-ZRVF(:,:,:)) &
                                                 / XRELAX_TIME_FRC
      END WHERE
!
    END IF
!
  END IF
!
  IF( LRELAX_UV_FRC ) THEN
!
!   apply UV relaxation
!
    WHERE( GRELAX_MASK_FRC )
      PRUS(:,:,:) = PRUS(:,:,:) - MXM(PRHODJ(:,:,:))*(PUT(:,:,:)-ZUF(:,:,:)) &
                                                 / XRELAX_TIME_FRC
      PRVS(:,:,:) = PRVS(:,:,:) - MYM(PRHODJ(:,:,:))*(PVT(:,:,:)-ZVF(:,:,:)) &
                                                 / XRELAX_TIME_FRC
    END WHERE
!
  END IF
!
END IF
!
!
!*       4.6    ground pressure forcing
!
IF( LPGROUND_FRC ) THEN
!
! to be implemented as a function of ZXPGROUNFRC
!
END IF
!
!
!*       5.     BUDGET CALLS
!   	        ------------
!
!
IF (LBUDGET_U)   CALL BUDGET (PRUS,1,'FRC_BU_RU')
IF (LBUDGET_V)   CALL BUDGET (PRVS,2,'FRC_BU_RV')
IF (LBUDGET_W)   CALL BUDGET (PRWS,3,'FRC_BU_RW')
IF (LBUDGET_TH)  CALL BUDGET (PRTHS,4,'FRC_BU_RTH')
IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'FRC_BU_RTKE')
IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'FRC_BU_RRV')
IF (LBUDGET_RC)  CALL BUDGET (PRRS(:,:,:,2),7,'FRC_BU_RRC')
IF (LBUDGET_RR)  CALL BUDGET (PRRS(:,:,:,3),8,'FRC_BU_RRR')
IF (LBUDGET_RI)  CALL BUDGET (PRRS(:,:,:,4),9,'FRC_BU_RRI')
IF (LBUDGET_RS)  CALL BUDGET (PRRS(:,:,:,5),10,'FRC_BU_RRS')
IF (LBUDGET_RG)  CALL BUDGET (PRRS(:,:,:,6),11,'FRC_BU_RRG')
IF (LBUDGET_RH)  CALL BUDGET (PRRS(:,:,:,7),12,'FRC_BU_RRH')
IF (LBUDGET_SV) THEN
  DO JL = 1 , SIZE(PRSVS,4)
    CALL BUDGET (PRSVS(:,:,:,JL),JL+12,'FRC_BU_RSV')
  END DO
END IF
!
!----------------------------------------------------------------------------
!
! deallocate work arrays
!
DEALLOCATE(ZUF)
DEALLOCATE(ZVF)
DEALLOCATE(ZWF)
DEALLOCATE(ZTHF)
DEALLOCATE(ZRVF)
DEALLOCATE(ZGXTHF)
DEALLOCATE(ZGYTHF)
DEALLOCATE(ZTENDTHF)
DEALLOCATE(ZTENDRVF)
DEALLOCATE(ZTENDUF)
DEALLOCATE(ZTENDVF)
DEALLOCATE(ZDZZ)
DEALLOCATE(ZRWCF)
DEALLOCATE(ZDUF)
DEALLOCATE(ZDVF)
!
!----------------------------------------------------------------------------
!
CONTAINS
  SUBROUTINE DEFINE_RELAX_FORCING( OMASK_RELAX_FRC,PTHT,PDZZ, &
                                                KPHEXT,KPVEXT )
!
  LOGICAL, DIMENSION(:,:,:), INTENT(OUT) :: OMASK_RELAX_FRC
  REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PTHT
  REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PDZZ
  INTEGER,                   INTENT(IN)  :: KPHEXT, KPVEXT
!
  INTEGER :: JI, JJ, JK, IIBEG, IIEND, IJBEG, IJEND, IKBEG, IKEND
  INTEGER :: IKGRAD_TH_MAX
  REAL    :: ZGRAD_TH, ZGRAD_TH_MAX
!
  IIBEG = KPHEXT + 1
  IIEND = SIZE(PTHT,1) - KPHEXT
  IJBEG = KPHEXT + 1
  IJEND = SIZE(PTHT,2) - KPHEXT
  IKBEG = KPHEXT + 1
  IKEND = SIZE(PTHT,3) - KPVEXT
!
  OMASK_RELAX_FRC(:,:,:) = .TRUE.
!
  DO JI = IIBEG , IIEND
    DO JJ = IJBEG , IJEND
      ZGRAD_TH_MAX = -1.E30
      DO JK = IKBEG+1 , IKEND-1
        ZGRAD_TH = (PTHT(JI,JJ,JK+1)-PTHT(JI,JJ,JK))/PDZZ(JI,JJ,JK+1)
        IF( ZGRAD_TH > ZGRAD_TH_MAX ) THEN
          IKGRAD_TH_MAX = JK
          ZGRAD_TH_MAX = ZGRAD_TH
        END IF
      END DO
      IKGRAD_TH_MAX = MIN( SIZE(PTHT,3),IKGRAD_TH_MAX+2 )
      OMASK_RELAX_FRC(:,:,1:IKGRAD_TH_MAX) = .FALSE.
    END DO
  END DO
!
  IF (NVERB >= 10) THEN
    WRITE(ILUOUT0,*) 'DEFINE_RELAX_FORCING: IKGRAD_TH_MAX = ',IKGRAD_TH_MAX
  END IF
!
  END SUBROUTINE DEFINE_RELAX_FORCING 
!
!----------------------------------------------------------------------------
!
END SUBROUTINE FORCING
