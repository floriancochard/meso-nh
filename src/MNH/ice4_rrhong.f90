!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_RRHONG
INTERFACE
SUBROUTINE ICE4_RRHONG(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT,   PRRT, &
                       &PTHT, &
                       &PRRHONG_MR, PB_TH, PB_RR, PB_RG)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                  INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRHONG_MR ! Mixing ratio change due to spontaneous freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
END SUBROUTINE ICE4_RRHONG
END INTERFACE
END MODULE MODI_ICE4_RRHONG
SUBROUTINE ICE4_RRHONG(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT,   PRRT, &
                       &PTHT, &
                       &PRRHONG_MR, PB_TH, PB_RR, PB_RG)
!!
!!**  PURPOSE
!!    -------
!!      Computes the RRHONG process
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
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                  INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRHONG_MR ! Mixing ratio change due to spontaneous freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
!
!*       0.2  declaration of local variables
!
LOGICAL, DIMENSION(SIZE(PRRT)) :: GMASK
!
!-------------------------------------------------------------------------------
!
!*       3.3     compute the spontaneous freezing source: RRHONG
!
PRRHONG_MR(:) = 0.
IF(.NOT. LDSOFT) THEN
  GMASK(:)=PT(:)<XTT-35.0 .AND. PRRT(:)>XRTMIN(3) .AND. LDCOMPUTE(:)
  WHERE(GMASK(:))
    PRRHONG_MR(:) = PRRT(:)
  ENDWHERE
  IF(LFEEDBACKT) THEN
    !Limitation due to -35 crossing of temperature
    WHERE(GMASK(:))
      PRRHONG_MR(:)=MIN(PRRHONG_MR(:), MAX(0., ((XTT-35.)/PEXN(:)-PTHT)/(PLSFACT(:)-PLVFACT(:))))
    ENDWHERE
  ENDIF
ENDIF
PB_RG(:) = PB_RG(:) + PRRHONG_MR(:)
PB_RR(:) = PB_RR(:) - PRRHONG_MR(:)
PB_TH(:) = PB_TH(:) + PRRHONG_MR(:)*(PLSFACT(:)-PLVFACT(:))
!
!
END SUBROUTINE ICE4_RRHONG
