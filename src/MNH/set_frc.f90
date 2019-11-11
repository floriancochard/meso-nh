!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_SET_FRC
!     ###################
!
INTERFACE
!
SUBROUTINE SET_FRC(TPEXPREFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
END SUBROUTINE SET_FRC
!
END INTERFACE
!
END MODULE MODI_SET_FRC
!
!     ##########################
      SUBROUTINE SET_FRC(TPEXPREFILE)
!     ##########################
!
!!*** *SET_FRC * -  to initialize forcing fields from successive soundings
!!
!!    PURPOSE
!!    -------
!!      The routine reads the free-format soundings and performs vertical
!!    interpolations to get the forcing fields on the model levels.
!!
!!**  METHOD
!!    ------
!!      The data are read in a prescribed format. If the data are given on
!!    pressure levels, a sounding is read to compute the corresponding heights
!!    with accuracy. The data are then vertically interpolated on the model
!!    levels and the horizontal gradient of theta is computed from the
!!    thermal wind balance.
!!
!!    EXTERNAL
!!    --------
!!      Module MODE_THERMO : contains thermodynamic routines
!!      Module  MODI_HEIGHT_PRESS: to compute height from pressure and
!!                                 potential virtual temperature
!!      Module  MODI_THETAVPU_THETAVPM: to interpolate thetav on wind levels
!!                                      from thetav on mass levels
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST       : contains Physical constants
!!        XPI   : Pi
!!        XOMEGA: Earth rotation
!!        XRV   : Gas constant for vapor
!!        XRD   : Gas constant for dry air
!!      Module MODD_LUNIT1  : contains logical unit names
!!        CLUOUT : name of output-listing
!!      Module MODD_GRID1: declaration of grid variables
!!        XZHAT: height levels without orography
!!      Module  MODD_CONF   : contains configuration variables for all models
!!        L1D  : 1D version flag
!!        NVERB: level of printed information on the output listing
!!      Module MODD_FRC: declaration of the forcing variables
!!        NFRC  : number of forcing variables
!!        TDTFRC: date of each forcing profile
!!        XUFRC,XVFRC,XWFRC,XTHFRC,XRVFRC: large scale forcing variables
!!        XGXTHFRC,XGYTHFRC: large scale gradients of THETA
!!        XTENDTHFRC,XTENDRVFRC: large scale tendencies of THETA, RV
!!      Module  MODD_GRID: declaration of grid variables
!!        XLAT0: Reference latitude to compute the Coriolis parameter
!!      Module  MODD_REF: declaration of reference state profile
!!        XTHVREFZ: Thetav(z) for reference state without orography
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine SET_FRC)
!!
!!    AUTHOR
!!    ------
!!      M. Georgelin       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/12/95
!!      29/07/96    (Pinty, Mari, Suhre) restructured
!!      10/12/96    (Pinty)              add the effect of humidity to the
!!                                       thermal wind balance
!!      24/01/98    (Bechtold) substitute humidity advection by tendency forc.
!!                             add SST and ground pressure forcing
!!         06/12    (Masson)   Removes extrapolations below or above forcing
!!                             data. Reproduces the same data instead.
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!      27/11/17    (Chaboureau) fix bug in allocation relative to LTEND_UV_FRC 
!!      28/03/2018  (P.Wautelet) use overloaded comparison operator for date_time
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST 
USE MODD_LUNIT_n
USE MODD_GRID_n
USE MODD_CONF
USE MODD_FRC
USE MODD_GRID
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_REF
USE MODD_PARAMETERS
!
USE MODE_DATETIME
USE MODE_THERMO
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_HEIGHT_PRESS  ! interface modules
USE MODI_PRESS_HEIGHT
USE MODI_THETAVPU_THETAVPM
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments :
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
!*       0.2   Declarations of local variables :
!
INTEGER :: ILUPRE,IRESP ! logical unit number of the  EXPRE and FM return code
INTEGER :: ILUOUT       ! Logical unit number for output-listing
INTEGER :: JK,JKT,JKLEV,JKU,JKM,JL ! Loop control
INTEGER :: IKU       ! Upper bound in z direction of model arrays
INTEGER :: IKE       ! Lower bound in z direction of model arrays
INTEGER :: ILEVELF   ! number of wind levels of the forcing sounding
INTEGER :: ILEVELMF  ! number of mass levels of the forcing sounding
!
REAL    :: ZF        ! Coriolis parameter (1D)
REAL    :: ZRVSRD    ! Rv/Rd
REAL    :: ZDZSDH,ZDZ1SDH,ZDZ2SDH ! for interpolation
REAL    :: ZZGROUNDF,ZPGROUNDF    ! height and Pressure at ground
REAL    :: ZTHDGROUNDF,ZMRGROUNDF ! temp. and moist forcing variables at ground
REAL, DIMENSION(:), ALLOCATABLE :: ZWF,ZUF,ZVF ! Local variables for
REAL, DIMENSION(:), ALLOCATABLE :: ZTHF,ZRVF   ! the data reading
REAL, DIMENSION(:), ALLOCATABLE :: ZGXRF,ZGYRF !          "
REAL, DIMENSION(:), ALLOCATABLE :: ZHEIGHTF    !          "
REAL, DIMENSION(:), ALLOCATABLE :: ZTUF, ZTVF
REAL, DIMENSION(:), ALLOCATABLE :: ZPRESSUF    !          "
REAL, DIMENSION(:), ALLOCATABLE :: ZPRESSMF    !          "
REAL, DIMENSION(:), ALLOCATABLE :: ZTHVUF      ! Thetav at wind levels
REAL, DIMENSION(:), ALLOCATABLE :: ZHEIGHTMF   ! Height at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTHVF       ! Thetav at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTHDF       ! Theta (dry) at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZMRF        ! Vapor mixing ratio at mass lev.
REAL, DIMENSION(SIZE(XZHAT))    :: ZZHATM      ! Height of mass model grid
                                               ! levels  without orography
REAL, DIMENSION(SIZE(XZHAT))    :: ZSHEAR      ! vertical wind shear
CHARACTER(LEN=4)                :: YZP         ! choice of zfrc or pfrc
CHARACTER(LEN=100)              :: YMSG
!
REAL, DIMENSION(SIZE(XZHAT))    :: ZEXREFZ,  & ! Pi_ref
                                   ZRVREFZ,  & ! r_vref
                                   ZTHREFZ     ! Th_ref
REAL                            :: ZEPS,     & ! XRV/XRD-1.0
                                   ZG_O_CPD    ! XG/CPD
!
!-------------------------------------------------------------------------------
!
!*       1.     PROLOGUE : RETRIEVE LOGICAL UNIT NUMBERS
!               ----------------------------------------
!
ILUPRE = TPEXPREFILE%NLU
ILUOUT = TLUOUT%NLU
!
ZRVSRD  = XRV/XRD
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE FORCING FIELDS PROFILES
!               -------------------------------
!
BACKSPACE(ILUPRE)
READ(ILUPRE,*) YZP   ! Z-altitude or P-altitude
!
IF( YZP/='PFRC' .AND. YZP/='ZFRC' ) THEN
 !callabortstop
  WRITE(YMSG,*) 'undefined type of forcing: ',TRIM(YZP),'. It should be PFRC or ZFRC'
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_FRC',YMSG)
END IF
!
READ(ILUPRE,*) NFRC  ! Number of time-dependent forcing soundings
!
! CAUTION: number of forcing times is limited by the WRITE format 99(8E10.3)
!          and also by the name of forcing variables (format I3.3)
!          You have to modify those if you need more forcing times :-(
IF (NFRC > 99*8) THEN
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_FRC','maximum forcing times NFRC is 99*8')
 !callabortstop
END IF
!
!* Allocate the MODD_FRC forcing arrays
!
IKU = SIZE(XZHAT)
IKE = IKU - JPVEXT
ALLOCATE(TDTFRC(NFRC))
ALLOCATE(XUFRC(IKU,NFRC))
ALLOCATE(XVFRC(IKU,NFRC))
ALLOCATE(XWFRC(IKU,NFRC))
ALLOCATE(XTHFRC(IKU,NFRC))
ALLOCATE(XRVFRC(IKU,NFRC))
ALLOCATE(XTENDTHFRC(IKU,NFRC))
ALLOCATE(XTENDRVFRC(IKU,NFRC))
ALLOCATE(XGXTHFRC(IKU,NFRC))
ALLOCATE(XGYTHFRC(IKU,NFRC))
ALLOCATE(XPGROUNDFRC(NFRC))
ALLOCATE(XTENDUFRC(IKU,NFRC))
ALLOCATE(XTENDVFRC(IKU,NFRC))
!
! Reading the forcing sounding written in prep_idea1.nam
!
DO JKT = 1,NFRC
!
  WRITE(UNIT=ILUOUT,FMT='(A, I4)') "SET_FRC: start reading forcing field ", JKT
!
  IF( YZP=='ZFRC' ) THEN
    READ(ILUPRE,*) TDTFRC(JKT)%TDATE%YEAR,  &
                   TDTFRC(JKT)%TDATE%MONTH, &
                   TDTFRC(JKT)%TDATE%DAY,   &
                   TDTFRC(JKT)%TIME
    READ(ILUPRE,*) ZZGROUNDF  ! data at ground level
    READ(ILUPRE,*) ZPGROUNDF
    READ(ILUPRE,*) ZTHDGROUNDF
    READ(ILUPRE,*) ZMRGROUNDF
!
    READ(ILUPRE,*) ILEVELF    ! Number of level
    IF (JKT.GT.1) THEN
      DEALLOCATE(ZHEIGHTF)
      DEALLOCATE(ZUF)
      DEALLOCATE(ZVF)
      DEALLOCATE(ZWF)
      DEALLOCATE(ZTHF)
      DEALLOCATE(ZRVF)
      DEALLOCATE(ZGXRF)
      DEALLOCATE(ZGYRF)
      DEALLOCATE(ZTUF)
      DEALLOCATE(ZTVF)
    ENDIF
    ALLOCATE(ZHEIGHTF(ILEVELF))
    ALLOCATE(ZUF(ILEVELF))
    ALLOCATE(ZVF(ILEVELF))
    ALLOCATE(ZWF(ILEVELF))
    ALLOCATE(ZTHF(ILEVELF))
    ALLOCATE(ZRVF(ILEVELF))
    ALLOCATE(ZGXRF(ILEVELF))
    ALLOCATE(ZGYRF(ILEVELF))
    ALLOCATE(ZTUF(ILEVELF))
    ALLOCATE(ZTVF(ILEVELF))
!
    DO JKU =1,ILEVELF
      READ(ILUPRE,*) ZHEIGHTF(JKU)                         &
                    ,ZUF(JKU),ZVF(JKU),ZTHF(JKU),ZRVF(JKU) &
                    ,ZWF(JKU),ZGXRF(JKU),ZGYRF(JKU),ZTUF(JKU)&
                    ,ZTVF(JKU)

    END DO
  END IF
!
  IF( YZP=='PFRC' ) THEN
    READ(ILUPRE,*) TDTFRC(JKT)%TDATE%YEAR,  &
                   TDTFRC(JKT)%TDATE%MONTH, &
                   TDTFRC(JKT)%TDATE%DAY,   &
                   TDTFRC(JKT)%TIME
!
    READ(ILUPRE,*) ZZGROUNDF  ! data at ground level
    READ(ILUPRE,*) ZPGROUNDF
    READ(ILUPRE,*) ZTHDGROUNDF
    READ(ILUPRE,*) ZMRGROUNDF
!
!
    READ(ILUPRE,*) ILEVELF    ! Number of level
    IF (JKT.GT.1) THEN
      DEALLOCATE(ZPRESSUF)
      DEALLOCATE(ZTHVUF)
      DEALLOCATE(ZHEIGHTF)
      DEALLOCATE(ZUF)
      DEALLOCATE(ZVF)
      DEALLOCATE(ZWF)
      DEALLOCATE(ZTHF)
      DEALLOCATE(ZRVF)
      DEALLOCATE(ZGXRF)
      DEALLOCATE(ZGYRF)
      DEALLOCATE(ZTUF)
      DEALLOCATE(ZTVF)
    END IF
    ALLOCATE(ZPRESSUF(ILEVELF))
    ALLOCATE(ZTHVUF(ILEVELF))
    ALLOCATE(ZHEIGHTF(ILEVELF))
    ALLOCATE(ZUF(ILEVELF))
    ALLOCATE(ZVF(ILEVELF))
    ALLOCATE(ZWF(ILEVELF))
    ALLOCATE(ZTHF(ILEVELF))
    ALLOCATE(ZRVF(ILEVELF))
    ALLOCATE(ZGXRF(ILEVELF))
    ALLOCATE(ZGYRF(ILEVELF))
    ALLOCATE(ZTUF(ILEVELF))
    ALLOCATE(ZTVF(ILEVELF))
!
    DO JKU =1,ILEVELF
      READ(ILUPRE,*) ZPRESSUF(JKU)                         &
                    ,ZUF(JKU),ZVF(JKU),ZTHF(JKU),ZRVF(JKU) &
                    ,ZWF(JKU),ZGXRF(JKU),ZGYRF(JKU),ZTUF(JKU)&
                    ,ZTVF(JKU)
    END DO
!
!  read sounding
!
    READ(ILUPRE,*) ILEVELMF
    IF (JKT.GT.1) THEN
      DEALLOCATE(ZPRESSMF)
      DEALLOCATE(ZTHDF)
      DEALLOCATE(ZMRF)
      DEALLOCATE(ZTHVF)
      DEALLOCATE(ZHEIGHTMF)
    END IF
    ALLOCATE(ZPRESSMF(ILEVELMF))
    ALLOCATE(ZTHDF(ILEVELMF))
    ALLOCATE(ZMRF(ILEVELMF))
    ALLOCATE(ZHEIGHTMF(ILEVELMF))
    ALLOCATE(ZTHVF(ILEVELMF))
!
    DO JKM=2,ILEVELMF
      READ(ILUPRE,*) ZPRESSMF(JKM),ZTHDF(JKM),ZMRF(JKM)
    END DO
!
! Complete the mass arrays with the ground informations read in EXPRE file
!
    ZPRESSMF(1) = ZPGROUNDF   ! Mass level 1 is at the ground level
    ZTHDF(1)    = ZTHDGROUNDF
    ZMRF(1)     = ZMRGROUNDF
!
! Compute thetav at the mass levels of the RS
!
    ZTHVF(:) = ZTHDF(:) * (1. + ZRVSRD *ZMRF(:))/(1.+ZMRF(:))
!
! Compute the heights at the mass levels of the RS
!
    ZHEIGHTMF(:) = HEIGHT_PRESS(ZPRESSMF,ZTHVF,ZPGROUNDF,ZTHVF(1),ZZGROUNDF)
!
! Compute thetav and heights of the wind levels
!
    ZTHVUF(:) = THETAVPU_THETAVPM(ZPRESSMF,ZPRESSUF,ZTHVF)
    ZHEIGHTF(:) = HEIGHT_PRESS(ZPRESSUF,ZTHVUF,ZPGROUNDF,ZTHVF(1),ZZGROUNDF)
!
  END IF
!
! Interpolate and extrapolate Wfrc on w-vertical-grid levels :
!
  DO JK = 1,IKU
    IF (XZHAT(JK) <= ZHEIGHTF(1)) THEN
!
! copy below the first level
!
      XWFRC(JK,JKT) = ZWF(1)
    ELSE IF (XZHAT(JK) > ZHEIGHTF(ILEVELF) ) THEN
!
! copy above the last level
!
      XWFRC(JK,JKT) = ZWF(ILEVELF)
    ELSE
!
! interpolation between first and last levels
!
      DO JKLEV = 1,ILEVELF-1
        IF ( (XZHAT(JK) > ZHEIGHTF(JKLEV)).AND. &
             (XZHAT(JK) <= ZHEIGHTF(JKLEV+1))   ) THEN
          ZDZ1SDH = (XZHAT(JK) - ZHEIGHTF(JKLEV)) /   &
                    (ZHEIGHTF(JKLEV+1)-ZHEIGHTF(JKLEV))
          ZDZ2SDH = 1.- ZDZ1SDH
          XWFRC(JK,JKT) = ZWF(JKLEV) * ZDZ2SDH + ZWF(JKLEV+1) *ZDZ1SDH
        END IF
      END DO
    END IF
  END DO
!
! Interpolate and extrapolate Ufrc on u-vertical-grid levels
!           the other forcing variables.
!
  ZZHATM(1:IKU-1) = 0.5*(XZHAT(2:IKU)+XZHAT(1:IKU-1))
  ZZHATM(IKU) = 2.*XZHAT(IKU)-ZZHATM(IKU-1)
!
  DO JK = 1,IKU
    IF (ZZHATM(JK) <= ZHEIGHTF(1)) THEN
!
! copy below the first level
!
      XUFRC(JK,JKT) = ZUF(1)
      XVFRC(JK,JKT) = ZVF(1)
      XTHFRC(JK,JKT) = ZTHF(1)
      XRVFRC(JK,JKT) = ZRVF(1)
      XTENDTHFRC(JK,JKT) = ZGXRF(1)
      XTENDRVFRC(JK,JKT) = ZGYRF(1)
      XTENDUFRC(JK,JKT) = ZTUF(1)
      XTENDVFRC(JK,JKT) = ZTVF(1)        
    ELSE IF (ZZHATM(JK) > ZHEIGHTF(ILEVELF) ) THEN
!
! copy above the last level
!
      XUFRC(JK,JKT)  = ZUF(ILEVELF)
      XVFRC(JK,JKT)  = ZVF(ILEVELF)
      XTHFRC(JK,JKT) = ZTHF(ILEVELF)
      XRVFRC(JK,JKT) = ZRVF(ILEVELF)
      XTENDTHFRC(JK,JKT)=ZGXRF(ILEVELF)
      XTENDRVFRC(JK,JKT)=ZGYRF(ILEVELF)
      XTENDUFRC(JK,JKT)=ZTUF(ILEVELF)
      XTENDVFRC(JK,JKT)=ZTVF(ILEVELF)
    ELSE
!
! interpolation between first and last levels
!
      DO JKLEV = 1,ILEVELF-1
        IF ( (ZZHATM(JK) > ZHEIGHTF(JKLEV)).AND. &
             (ZZHATM(JK) <= ZHEIGHTF(JKLEV+1))   ) THEN
          ZDZ1SDH = (ZZHATM(JK) - ZHEIGHTF(JKLEV)) /  &
                    (ZHEIGHTF(JKLEV+1)-ZHEIGHTF(JKLEV))
          ZDZ2SDH = 1.- ZDZ1SDH
          XUFRC(JK,JKT)  = ZUF(JKLEV)*ZDZ2SDH  + ZUF(JKLEV+1)*ZDZ1SDH
          XVFRC(JK,JKT ) = ZVF(JKLEV)*ZDZ2SDH  + ZVF(JKLEV+1)*ZDZ1SDH
          XTHFRC(JK,JKT) = ZTHF(JKLEV)*ZDZ2SDH + ZTHF(JKLEV+1)*ZDZ1SDH
          XRVFRC(JK,JKT) = ZRVF(JKLEV)*ZDZ2SDH + ZRVF(JKLEV+1)*ZDZ1SDH
          XTENDTHFRC(JK,JKT)=ZGXRF(JKLEV)*ZDZ2SDH + ZGXRF(JKLEV+1)*ZDZ1SDH
          XTENDRVFRC(JK,JKT)=ZGYRF(JKLEV)*ZDZ2SDH + ZGYRF(JKLEV+1)*ZDZ1SDH
          XTENDUFRC(JK,JKT)=ZTUF(JKLEV)*ZDZ2SDH + ZTUF(JKLEV+1)*ZDZ1SDH    
          XTENDVFRC(JK,JKT)=ZTVF(JKLEV)*ZDZ2SDH + ZTVF(JKLEV+1)*ZDZ1SDH
        END IF
      END DO
    END IF
  END DO
!
! store ground pressure for forcing
!
  XPGROUNDFRC(JKT) = ZPGROUNDF
!
END DO ! End of loop in time
!
! check whether soundings have been entered in increasing temporal order
!

DO JKT = 2,NFRC-1
  IF ( TDTFRC(JKT-1) >= TDTFRC(JKT) ) THEN
    WRITE(ILUOUT,*) &
      "SET_FRC ERROR: sounding", JKT-1, " is given for a later time than", JKT
    WRITE(ILUOUT,*) &
      "               soundings have to be entered in increasing temporal order"
    WRITE(ILUOUT,*) "SOUNDING TIME ", JKT-1, " IS: "
    WRITE(ILUOUT,*) TDTFRC(JKT-1)%TDATE%YEAR,  &
                    TDTFRC(JKT-1)%TDATE%MONTH, &
                    TDTFRC(JKT-1)%TDATE%DAY,   &
                    TDTFRC(JKT-1)%TIME
    WRITE(ILUOUT,*) "SOUNDING TIME ", JKT, " IS: "
    WRITE(ILUOUT,*) TDTFRC(JKT)%TDATE%YEAR,  &
                    TDTFRC(JKT)%TDATE%MONTH, &
                    TDTFRC(JKT)%TDATE%DAY,   &
                    TDTFRC(JKT)%TIME
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_FRC','')
  END IF
END DO
!
! deduce gradient of theta from the thermal wind balance
!
ZF=2*XOMEGA*SIN(XLAT0*XPI/180.)
!
! Coriolis parameter is constant (Cartesian case)
! Humidity correction is performed
!
ZEPS      = XRV/XRD-1.0
ZG_O_CPD  = XG/XCPD
ZEXREFZ(IKE) = XEXNTOP + 0.5 * ZG_O_CPD * (XZHAT(IKE+1) - XZHAT(IKE)) &
                           / XTHVREFZ(IKE)
DO JK = IKE-1, 1, -1
  ZEXREFZ(JK)= ZEXREFZ(JK+1) + ZG_O_CPD * (XZHAT(JK+2) - XZHAT(JK))  &
                         / (XTHVREFZ(JK+1)+XTHVREFZ(JK))
END DO
DO JK = IKE+1, IKU-1
  ZEXREFZ(JK)= ZEXREFZ(JK-1) + ZG_O_CPD * (XZHAT(JK-1) - XZHAT(JK+1)) &
                         / (XTHVREFZ(JK-1)+XTHVREFZ(JK))
END DO
ZEXREFZ(IKU)= ZEXREFZ(IKU-1) + ZG_O_CPD * 2.*(XZHAT(IKU-1) - XZHAT(IKU)) &
                       / (XTHVREFZ(IKU-1)+XTHVREFZ(IKU))

ZRVREFZ(:) = ( (ZEXREFZ(:)**(XCPD/XRD-1.)*XP00)/      &
                  (XRD*XTHVREFZ(:)*XRHODREFZ(:)) ) - 1.0 ! r_vref
ZTHREFZ(:) = XTHVREFZ(:)/(1.+ZRVREFZ(:)*(XRV/XRD))       ! Th_ref
DO JKT=1,NFRC
  DO JK=1,IKU-1
    ZSHEAR(JK)=(XVFRC(JK+1,JKT)-XVFRC(JK,JKT)) /(XZHAT(JK+1)-XZHAT(JK))
  END DO
  DO JK=2,IKU-1
    XGXTHFRC(JK,JKT)=-(ZF*XTHVREFZ(JK)/XG)*0.5*(ZSHEAR(JK)+ZSHEAR(JK-1))
    XGXTHFRC(JK,JKT)=( XGXTHFRC(JK,JKT) - ZTHREFZ(JK)*ZEPS*XTENDTHFRC(JK,JKT) ) &
                               /( 1.0 + ZEPS*ZRVREFZ(JK) )
  END DO
  XGXTHFRC(1,JKT)= 0.
  XGXTHFRC(IKU,JKT)=XGXTHFRC(IKU-1,JKT)
!
  DO JK=1,IKU-1
    ZSHEAR(JK)=(XUFRC(JK+1,JKT)-XUFRC(JK,JKT)) /(XZHAT(JK+1)-XZHAT(JK))
  END DO
  DO JK=2,IKU-1
    XGYTHFRC(JK,JKT)=+(ZF*XTHVREFZ(JK)/XG)*0.5*(ZSHEAR(JK)+ZSHEAR(JK-1))
    XGYTHFRC(JK,JKT)=( XGYTHFRC(JK,JKT) - ZTHREFZ(JK)*ZEPS*XTENDRVFRC(JK,JKT) )&
                               /( 1.0 + ZEPS*ZRVREFZ(JK) )
  END DO
  XGYTHFRC(1,JKT)= 0.
  XGYTHFRC(IKU,JKT)=XGYTHFRC(IKU-1,JKT)
END DO
!
!-------------------------------------------------------------------------------
!
!*        3.    PRINT FORCING FIELDS
!         --------------------------
!
WRITE(UNIT=ILUOUT,FMT='(" THERE ARE ",I2," FORCING SOUNDINGS AT:")') NFRC
DO JL = 1 , NFRC
  WRITE(UNIT=ILUOUT,FMT='(F9.0, "s, date:", I3, "/", I3, "/", I5)') &
    TDTFRC(JL)%TIME,        &
    TDTFRC(JL)%TDATE%DAY,   &
    TDTFRC(JL)%TDATE%MONTH, &
    TDTFRC(JL)%TDATE%YEAR
END DO
!
IF (NVERB >= 10) THEN
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "FORCING: the following forcing fields will be used:"
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XUFRC: geostrophic wind in X"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8F10.2))') &
         JK, (XUFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
         "XVFRC: geostrophic wind in Y"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8F10.2))') &
         JK, (XVFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XTHFRC: THETA for relaxation"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8F10.2))') &
         JK, (XTHFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XRVFRC: RV for relaxation"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8F10.6))') &
         JK, (XRVFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XWFRC: vertical transport velocity"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8F10.7))') &
         JK, (XWFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XGXTHFRC: thermal wind advection in X"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8E10.3))') &
         JK, (XGXTHFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XGYTHFRC: thermal wind advection in Y"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8E10.3))') &
         JK, (XGYTHFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XTENDTHFRC: Theta tendency"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8E10.3))') &
         JK, (XTENDTHFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XTENDRVFRC: humidity tendency"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8E10.3))') &
         JK, (XTENDRVFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XTENDUFRC : wind advection tendency in X"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8E10.3))') &
         JK, (XTENDUFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XTENDVFRC : wind advection tendency in Y"
  DO JK = 1, IKU
    WRITE(UNIT=ILUOUT,FMT='(I10,99(/8E10.3))') &
         JK, (XTENDVFRC(JK,JL), JL=1, NFRC)
  END DO
!
  WRITE(UNIT=ILUOUT,FMT='(A)') &
       "XPGROUNDFRC: SURF PRESSURE FORCING"
  WRITE(UNIT=ILUOUT,FMT='(10X,99(/8E10.3))') &
       (XPGROUNDFRC(JL), JL=1, NFRC)
!
END IF
!
!-------------------------------------------------------------------------------
END SUBROUTINE SET_FRC
