!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_FAST_RI
INTERFACE
SUBROUTINE ICE4_FAST_RI(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, &
                       &PAI, PCJ, PCIT, &
                       &PSSI, &
                       &PRCT, PRIT, &
                       &PRCBERI, PA_TH, PA_RC, PA_RI)
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
USE MODI_BUDGET
USE MODD_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCBERI  ! Bergeron-Findeisen effect
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
END SUBROUTINE ICE4_FAST_RI
END INTERFACE
END MODULE MODI_ICE4_FAST_RI
SUBROUTINE ICE4_FAST_RI(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, &
                       &PAI, PCJ, PCIT, &
                       &PSSI, &
                       &PRCT, PRIT, &
                       &PRCBERI, PA_TH, PA_RC, PA_RI)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast ri process
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
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCBERI  ! Bergeron-Findeisen effect
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PRHODREF)) :: ZZW
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMASK
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!
GMASK(:)=PSSI(:)>0. .AND. PRCT(:)>XRTMIN(2) .AND. PRIT(:)>XRTMIN(4) .AND. &
        &PCIT(:)>0. .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRCBERI(:) = 0.
  END WHERE
ELSE
  PRCBERI(:) = 0.
  WHERE(GMASK(:))
    PRCBERI(:) = MIN(1.E8, XLBI*(PRHODREF(:)*PRIT(:)/PCIT(:))**XLBEXI) ! Lbda_i
    PRCBERI(:) = ( PSSI(:) / (PRHODREF(:)*PAI(:)) ) * PCIT(:) * &
                 ( X0DEPI/PRCBERI(:) + X2DEPI*PCJ(:)*PCJ(:)/PRCBERI(:)**(XDI+2.0) )
  END WHERE
ENDIF
PA_RC(:) = PA_RC(:) - PRCBERI(:)
PA_RI(:) = PA_RI(:) + PRCBERI(:)
PA_TH(:) = PA_TH(:) + PRCBERI(:)*(PLSFACT(:)-PLVFACT(:))
!
!
END SUBROUTINE ICE4_FAST_RI
