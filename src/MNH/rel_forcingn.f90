!MNH_LIC Copyright 2010-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #####################
      MODULE MODI_REL_FORCING_n
!     #####################
!
INTERFACE
!
      SUBROUTINE REL_FORCING_n ( PRHODJ, TPDTCUR,PTHM,PRM, PZZ,PRTHS, PRRS)
!
USE MODD_TIME, ONLY: DATE_TIME
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! ( rhod J ) = dry density
              ! for reference state * Jacobian of the GCS transformation.
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR ! current date and time
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::PRTHS ! potential temperature tendencies at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS ! moist variables tendencies at time t+1
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ     ! height z
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM !  moist variables at time t-dt
!
END SUBROUTINE REL_FORCING_n
!
END INTERFACE
!
END MODULE MODI_REL_FORCING_n
!
!     ######################################################################
      SUBROUTINE REL_FORCING_n ( PRHODJ, TPDTCUR, PTHM,PRM, PZZ,PRTHS, PRRS)
!     ######################################################################
!
!!***  *REL_FORCING* - routine to compute the relaxtaion-forced terms for 2D runs 
!!
!!    PURPOSE
!!    -------
!!      The routine prepares (linear interpolations) and integrates 
!!      non homogenous relaxation on theta and rv.
!!   
!!**  METHOD
!!    ------
!!      For its first call, the routine looks for a starting advecting-forcing 
!!    with a date_and_time immediately lower or close to that the current
!!    date_and_time of the model. Then the temporal interpolation or extension
!!    is performed according to the position of the current date_and_time
!!    as compared to that of the advecting-forcing. In case of non-flat
!!    terrain, no interpolation is anticipated.
!!      All the necessary interpolations are linear. 
!!
!!   NB:   For relaxation forcing, only mask=FIXE has been implemented for simplicity
!!
!!   DUMMIES: LDUMMY(2)=T allows ADV forcing
!!            LDUMMY3=T ------- REL -------
!!                   with XDUMMY1=lower limit of relaxation (m)
!!                        XDUMMY2=top limit of relxation (m)
!!                        XDUMMY3=relaxation timsescale (s)
!!
!!    EXTERNAL
!!    --------
!!      Temporal_lt function   (compare 2 TYPEd date_and_time data)
!!      Temporal_dist function (compute the number of seconds between
!!                              2 TYPEd date_and_time data)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODDB_ADVFRCn: declaration of the advecting-forcing variables
!!        NADVFRC  : number of advecting-forcing variables
!!        TDTADVFRC: date of each advecting-forcing profile
!!        XUFRC,XVFRC,XWFRC,XTHFRC,XRVFRC: advecting-forcing variables
!!      Module MODD_LUNIT :  contains logical unit names for all models
!!        CLUOUT0 : name of output-listing
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPVEXT: define the number of marginal points out of the 
!!        physical domain along the vertical direction.    
!!      Module MODD_TIME: contains the structure of the TYPEd date_and_time
!!      Module MODD_BLANK: Uses LDUMMY(2)=T to activate the time varying adv frc 
!!
!!    REFERENCE
!!    ---------
!!      Peyrille&Lafore JAS 2007 vol64 n°8 An idealized two-dimensional framework to study 
!!      the west african monsoon. Part II Large-scale advection and the diurnal cycle
!!
!!    AUTHOR
!!    ------
!!	    M. Tomasini (CNRM) from forcing.f90 
!!      and P.Peyrille (CNRM)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/11/10
!!     28/03/2018 P. Wautelet: replace TEMPORAL_DIST by DATETIME_DISTANCE
!!                             use overloaded comparison operator for date_time
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_2D_FRC
USE MODD_BUDGET
USE MODD_DYN
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_RELFRC_n     ! Modules for time evolving advfrc
USE MODD_TIME
!
USE MODE_DATETIME
USE MODE_IO_ll
!
USE MODI_BUDGET
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! ( rhod J ) = dry density
              ! for reference state * Jacobian of the GCS transformation.
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR ! current date and time
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::PRTHS ! potential temperature tendencies at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS ! moist variables tendencies at time t+1
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ     ! height z
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM !  moist variables at time t-dt
!
!*       0.2   Declarations of local variables
!
INTEGER                         :: JN, JK, JXP
INTEGER, SAVE                   :: JSX_REL                ! saved loop index

LOGICAL, SAVE :: GSFIRSTCALL = .TRUE. ! control switch for the first call
!
REAL :: ZDT, ZALPHA ! height and time rate
REAL, SAVE :: ZSDTJX
!          
INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
!
REAL, DIMENSION(SIZE(PRTHS,1),SIZE(PRTHS,2),SIZE(PRTHS,3)) :: ZXADVTHFRC,ZXADVRVFRC
REAL, DIMENSION(SIZE(PRTHS,1),SIZE(PRTHS,2),SIZE(PRTHS,3)) :: ZTHREL,ZRVREL
LOGICAL,DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) :: GRELAX_MASK_FRC ! MAsk for relaxation
REAL :: ZRELAX_HEIGHT_TOP,ZRELAX_HEIGHT_BOT, ZRELAX_TIME
INTEGER :: IKU

!----------------------------------------------------------------------------
!
IKU = SIZE(PTHM,3)
!*        1.   PREPARATION OF FORCING
!              ----------------------
!

IF (GSFIRSTCALL) THEN
!
  GSFIRSTCALL = .FALSE.
!
  ILUOUT0 = TLUOUT0%NLU
!!*        1.1  printout number of forcing profiles
!
  WRITE(UNIT=ILUOUT0,FMT='(" THERE ARE ",I2," REL FORCING FIELDs  AT:")') NRELFRC
  DO JSX_REL = 1 , NRELFRC
    WRITE(UNIT=ILUOUT0,FMT='(F9.0, "s, date:", I3, "/", I3, "/", I5)') &
      TDTRELFRC(JSX_REL)%TIME,        &
      TDTRELFRC(JSX_REL)%TDATE%DAY,   &
      TDTRELFRC(JSX_REL)%TDATE%MONTH, &
      TDTRELFRC(JSX_REL)%TDATE%YEAR
  END DO

!*        1.2  find first sounding to be used 
  JSX_REL = 0
  IF( TPDTCUR < TDTRELFRC(1) ) THEN
    WRITE(UNIT=ILUOUT0,FMT='(" THE INITIAL REL FORCING FIELDS ARE NULL ")') 
  ELSE IF( TPDTCUR >= TDTRELFRC(nRELFRC) ) THEN
      WRITE(UNIT=ILUOUT0,FMT='(" THE REL FORCING FIELDS WILL REMAIN STATIONARY ")')
  ELSE
    TIM1_FOR:  DO JN = NRELFRC-1, 1, -1
                  JSX_REL = JN
                  IF( TPDTCUR >= TDTRELFRC(JSX_REL) ) EXIT TIM1_FOR
               END DO TIM1_FOR
!
    WRITE(UNIT=ILUOUT0,FMT='(" THE INITIAL FORCING FIELDS ARE INTERPOLATED" , &
                       & " IN TIME STARTING FROM THE SOUNDING NUMBER ",I2)') JSX_REL
    JSX_REL = JSX_REL - 1
  END IF
END IF
!
!*       2.     INTEGRATION OF TH and RV ADVECTING FORCINGS TENDANCY IN THE SOURCES
!   	        ---------------------------------------------------------------------
!
!    2.1 Temporal interpolation of each term
!   ------------------------------------------
IF( TPDTCUR < TDTRELFRC(1) ) THEN
  ZTHREL(:,:,:)   = 0.
  ZRVREL(:,:,:)   = 0.
ELSE IF ( TPDTCUR >= TDTRELFRC(NRELFRC) ) THEN
   ZTHREL(:,:,:)   = XTHREL(:,:,:,NRELFRC)
   ZRVREL(:,:,:)   = XRVREL(:,:,:,NRELFRC)
ELSE
  JXP = JSX_REL + 1

  IF( TPDTCUR >= TDTRELFRC(JXP) ) THEN
    JSX_REL = JSX_REL +1
    JXP= JSX_REL +1
    ILUOUT0 = TLUOUT0%NLU
    WRITE(UNIT=ILUOUT0,FMT='(" THE REL FORCING FIELDS ARE INTERPOLATED NOW" ,&
    & " BETWEEN SOUNDING NUMBER ",I2," AND SOUNDING NUMBER ",I2)') JSX_REL,JXP
    CALL DATETIME_DISTANCE(TDTRELFRC(JSX_REL),TDTRELFRC(JXP),ZSDTJX)
  END IF
!
  CALL DATETIME_DISTANCE(TDTRELFRC(JSX_REL),TPDTCUR,ZDT)
!
  ZALPHA = ZDT / ZSDTJX
!
!  temperature and humidity 2D reference profiles depending on time
!
  ZTHREL(:,:,:)   = XTHREL(:,:,:,JSX_REL)  +(XTHREL(:,:,:,JXP)-XTHREL(:,:,:,JSX_REL))*ZALPHA
  ZRVREL(:,:,:)   = XRVREL(:,:,:,JSX_REL)  +(XRVREL(:,:,:,JXP)-XRVREL(:,:,:,JSX_REL))*ZALPHA
END IF
!
!    2.2 Integration of the  forcing in the source
!   ------------------------------------------


!
!    2.2.2 relaxation forcing in 
! ---------------------------------------------------
  ZRELAX_HEIGHT_BOT=XRELAX_HEIGHT_BOT
  ZRELAX_HEIGHT_TOP=XRELAX_HEIGHT_TOP
  ZRELAX_TIME=XRELAX_TIME
  !
 ! define the mask where the relaxation is to be applied
 ! Corresponds to CASE=FIXE of forcing.f90
 !
  GRELAX_MASK_FRC(:,:,:) = .TRUE.
  WHERE ((MZF(1,IKU,1,PZZ).LT.ZRELAX_HEIGHT_BOT).OR.(MZF(1,IKU,1,PZZ).GT.ZRELAX_HEIGHT_TOP)) 
    GRELAX_MASK_FRC = .FALSE.
  END WHERE
!
  WHERE( GRELAX_MASK_FRC )
    PRTHS(:,:,:) = PRTHS(:,:,:) - PRHODJ(:,:,:)*(PTHM-ZTHREL) / ZRELAX_TIME
  END WHERE
!!
  WHERE( GRELAX_MASK_FRC )
    PRRS(:,:,:,1) = PRRS(:,:,:,1) - & 
    PRHODJ(:,:,:) * (PRM(:,:,:,1)-ZRVREL(:,:,:))/ ZRELAX_TIME
  END WHERE
!
!*       3.     BUDGET CALLS
!   	        ------------
IF (LBUDGET_TH)  CALL BUDGET (PRTHS,4,'2DREL_BU_RTH')
IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'2DREL_BU_RRV')
!----------------------------------------------------------------------------
!
END SUBROUTINE REL_FORCING_n
