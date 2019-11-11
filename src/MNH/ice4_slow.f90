!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_SLOW
INTERFACE
SUBROUTINE ICE4_SLOW(KSIZE, LDSOFT, LDCOMPUTE, PRHODREF, PT,&
                     &PSSI, PLVFACT, PLSFACT, &
                     &PRVT, PRCT, PRIT, PRST, PRGT,&
                     &PLBDAS, PLBDAG,&
                     &PAI, PCJ,&
                     &PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG, &
                     &PA_TH, PA_RV, PA_RC, PA_RI, PA_RS, PA_RG)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel/hail m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCHONI  ! Homogeneous nucleation
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRVDEPS  ! Deposition on r_s
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIAGGS  ! Aggregation on r_s
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIAUTS  ! Autoconversion of r_i for r_s production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRVDEPG  ! Deposition on r_g
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RV
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
END SUBROUTINE ICE4_SLOW
END INTERFACE
END MODULE MODI_ICE4_SLOW
SUBROUTINE ICE4_SLOW(KSIZE, LDSOFT, LDCOMPUTE, PRHODREF, PT, &
                     &PSSI, PLVFACT, PLSFACT, &
                     &PRVT, PRCT, PRIT, PRST, PRGT, &
                     &PLBDAS, PLBDAG, &
                     &PAI, PCJ, &
                     &PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG, &
                     &PA_TH, PA_RV, PA_RC, PA_RI, PA_RS, PA_RG)
!!
!!**  PURPOSE
!!    -------
!!      Computes the slow process
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel/hail m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCHONI  ! Homogeneous nucleation
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRVDEPS  ! Deposition on r_s
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIAGGS  ! Aggregation on r_s
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIAUTS  ! Autoconversion of r_i for r_s production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRVDEPG  ! Deposition on r_g
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RV
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PRHODREF)) :: ZCRIAUTI
REAL            :: ZTIMAUTIC
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMASK
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!
GMASK(:)=PT(:)<XTT-35.0 .AND. PRCT(:)>XRTMIN(2) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRCHONI(:) = 0.
  END WHERE
ELSE
  PRCHONI(:) = 0.
  WHERE(GMASK(:))
    PRCHONI(:) = XHON*PRHODREF(:)*PRCT(:)       &
                                 *EXP( XALPHA3*(PT(:)-XTT)-XBETA3 )
  ENDWHERE
ENDIF
PA_RI(:) = PA_RI(:) + PRCHONI(:)
PA_RC(:) = PA_RC(:) - PRCHONI(:)
PA_TH(:) = PA_TH(:) + PRCHONI(:)*(PLSFACT(:)-PLVFACT(:))
!
!*       3.4    compute the deposition, aggregation and autoconversion sources
!
!
!*       3.4.2  compute the riming-conversion of r_c for r_i production: RCAUTI
!
!  ZZW(:) = 0.0
!  ZTIMAUTIC = SQRT( XTIMAUTI*XTIMAUTC )
!  WHERE ( (PRCT(:)>0.0) .AND. (PRIT(:)>0.0) .AND. (PRCS(:)>0.0) )
!    ZZW(:) = MIN( PRCS(:),ZTIMAUTIC * MAX( SQRT( PRIT(:)*PRCT(:) ),0.0 ) )
!    PRIS(:) = PRIS(:) + ZZW(:)
!    PRCS(:) = PRCS(:) - ZZW(:)
!    PTHS(:) = PTHS(:) + ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RCAUTI))
!  END WHERE
!
!*       3.4.3  compute the deposition on r_s: RVDEPS
!
GMASK(:)=PRVT(:)>XRTMIN(1) .AND. PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRVDEPS(:) = 0.
  END WHERE
ELSE
  PRVDEPS(:) = 0.
  WHERE(GMASK(:))
    PRVDEPS(:) = ( PSSI(:)/(PRHODREF(:)*PAI(:)) ) *                               &
                 ( X0DEPS*PLBDAS(:)**XEX0DEPS + X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS )
  END WHERE
ENDIF
PA_RS(:) = PA_RS(:) + PRVDEPS(:)
PA_RV(:) = PA_RV(:) - PRVDEPS(:)
PA_TH(:) = PA_TH(:) + PRVDEPS(:)*PLSFACT(:)
!
!*       3.4.4  compute the aggregation on r_s: RIAGGS
!
GMASK(:)=PRIT(:)>XRTMIN(4) .AND. PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRIAGGS(:) = 0.
  END WHERE
ELSE
  PRIAGGS(:) = 0.
  WHERE(GMASK(:))
    PRIAGGS(:) = XFIAGGS * EXP( XCOLEXIS*(PT(:)-XTT) ) &
                         * PRIT(:)                      &
                         * PLBDAS(:)**XEXIAGGS          &
                         * PRHODREF(:)**(-XCEXVT)
  END WHERE
ENDIF
PA_RS(:) = PA_RS(:) + PRIAGGS(:)
PA_RI(:) = PA_RI(:) - PRIAGGS(:)
!
!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS
!
GMASK(:)=PRIT(:)>XRTMIN(4) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRIAUTS(:) = 0.
  END WHERE
ELSE
  PRIAUTS(:) = 0.
  !ZCRIAUTI(:)=MIN(XCRIAUTI,10**(0.06*(PT(:)-XTT)-3.5))
  ZCRIAUTI(:)=MIN(XCRIAUTI,10**(XACRIAUTI*(PT(:)-XTT)+XBCRIAUTI))
  WHERE(GMASK(:))
    PRIAUTS(:) = XTIMAUTI * EXP( XTEXAUTI*(PT(:)-XTT) ) &
                          * MAX( PRIT(:)-ZCRIAUTI(:),0.0 )
  END WHERE
ENDIF
PA_RS(:) = PA_RS(:) + PRIAUTS(:)
PA_RI(:) = PA_RI(:) - PRIAUTS(:)
!
!*       3.4.6  compute the deposition on r_g: RVDEPG
!
!
GMASK(:)=PRVT(:)>XRTMIN(1) .AND. PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRVDEPG(:) = 0.
  END WHERE
ELSE
  PRVDEPG(:) = 0.
  WHERE(GMASK(:))
    PRVDEPG(:) = ( PSSI(:)/(PRHODREF(:)*PAI(:)) ) *                               &
                 ( X0DEPG*PLBDAG(:)**XEX0DEPG + X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG )
  END WHERE
ENDIF
PA_RG(:) = PA_RG(:) + PRVDEPG(:)
PA_RV(:) = PA_RV(:) - PRVDEPG(:)
PA_TH(:) = PA_TH(:) + PRVDEPG(:)*PLSFACT(:)
!
!
END SUBROUTINE ICE4_SLOW
