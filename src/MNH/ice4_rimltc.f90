!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_RIMLTC
INTERFACE
SUBROUTINE ICE4_RIMLTC(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT, &
                       &PTHT, PRIT, &
                       &PRIMLTC_MR, PB_TH, PB_RC, PB_RI)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                  INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Cloud ice at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIMLTC_MR ! Mixing ratio change due to cloud ice melting
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RI
END SUBROUTINE ICE4_RIMLTC
END INTERFACE
END MODULE MODI_ICE4_RIMLTC
SUBROUTINE ICE4_RIMLTC(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT, &
                       &PTHT, PRIT, &
                       &PRIMLTC_MR, PB_TH, PB_RC, PB_RI)
!!
!!**  PURPOSE
!!    -------
!!      Computes the RIMLTC process
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
USE MODD_PARAM_ICE, ONLY : LFEEDBACKT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Cloud ice at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIMLTC_MR ! Mixing ratio change due to cloud ice melting
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RI
!
!*       0.2  declaration of local variables
!
LOGICAL, DIMENSION(KSIZE) :: GMASK
!
!-------------------------------------------------------------------------------
!
!*       7.1    cloud ice melting
!
PRIMLTC_MR(:)=0.
IF(.NOT. LDSOFT) THEN
  GMASK(:)=PRIT(:)>0. .AND. PT(:)>XTT .AND. LDCOMPUTE(:)
  WHERE(GMASK(:))
    PRIMLTC_MR(:)=PRIT(:)
  END WHERE

  IF(LFEEDBACKT) THEN
    !Limitation due to 0 crossing of temperature
    WHERE(GMASK(:))
      PRIMLTC_MR(:)=MIN(PRIMLTC_MR(:), MAX(0., (PTHT(:)-XTT/PEXN(:)) / (PLSFACT(:)-PLVFACT(:))))
    END WHERE
  ENDIF
ENDIF
PB_RC(:) = PB_RC(:) + PRIMLTC_MR(:)
PB_RI(:) = PB_RI(:) - PRIMLTC_MR(:)
PB_TH(:) = PB_TH(:) - PRIMLTC_MR(:)*(PLSFACT(:)-PLVFACT(:))
!
!
END SUBROUTINE ICE4_RIMLTC
