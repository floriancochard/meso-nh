!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ###############################
       MODULE MODI_LIMA_NUCLEATION_PROCS
!      ###############################
!
INTERFACE
   SUBROUTINE LIMA_NUCLEATION_PROCS (PTSTEP, TPFILE, OCLOSE_OUT, PRHODJ,          &
                                     PRHODREF, PEXNREF, PPABST, PT, PTM, PW_NU,    &
                                     PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,     &
                                     PCCT, PCRT, PCIT,                             &
                                     PNFT, PNAT, PIFT, PINT, PNIT, PNHT            )
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE     ! Output file
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTM        ! Temperature at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU      ! updraft velocity used for
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT       ! Theta at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT       ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST       ! Snow m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT       ! Graupel m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT       ! Prinstine ice conc. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFT       ! CCN C. available at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAT       ! CCN C. activated at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFT       ! IFN C. available at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT       ! IFN C. activated at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIT       ! Coated IFN activated at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHT       ! CCN hom freezing
!
END SUBROUTINE LIMA_NUCLEATION_PROCS
END INTERFACE
END MODULE MODI_LIMA_NUCLEATION_PROCS
!     #############################################################################
SUBROUTINE LIMA_NUCLEATION_PROCS (PTSTEP, TPFILE, OCLOSE_OUT, PRHODJ,                &
                                  PRHODREF, PEXNREF, PPABST, PT, PTM, PW_NU,          &
                                  PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,           &
                                  PCCT, PCRT, PCIT,                                   &
                                  PNFT, PNAT, PIFT, PINT, PNIT, PNHT                  )
!     #############################################################################
!
!!    PURPOSE
!!    -------
!!      Compute nucleation processes for the time-splitted version of LIMA
!!
!!    AUTHOR
!!    ------
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!-------------------------------------------------------------------------------
!
USE MODD_PARAM_LIMA, ONLY : LCOLD, LNUCL, LMEYERS, LSNOW, LWARM, LACTI, LRAIN, LHHONI, &
                            NMOD_CCN, NMOD_IFN, NMOD_IMM
USE MODD_BUDGET,     ONLY : LBU_ENABLE, LBUDGET_TH, LBUDGET_RV, LBUDGET_RC, LBUDGET_RR,&
                            LBUDGET_RI, LBUDGET_RS, LBUDGET_RG, LBUDGET_RH, LBUDGET_SV
USE MODD_NSV,        ONLY : NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_CCN_FREE,               &
                            NSV_LIMA_NI, NSV_LIMA_IFN_FREE
!
USE MODD_IO_ll,   ONLY: TFILEDATA
USE MODI_BUDGET
USE MODI_LIMA_CCN_ACTIVATION
USE MODI_LIMA_PHILLIPS_IFN_NUCLEATION
USE MODI_LIMA_MEYERS_NUCLEATION
USE MODI_LIMA_CCN_HOM_FREEZING
!
!-------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE     ! Output file
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTM        ! Temperature at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU      ! updraft velocity used for
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT       ! Theta at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT       ! Rain water m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT       ! Prinstine ice conc. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFT       ! CCN C. available at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAT       ! CCN C. activated at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFT       ! IFN C. available at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT       ! IFN C. activated at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIT       ! Coated IFN activated at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHT       ! CCN hom. freezing
!
!-------------------------------------------------------------------------------
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))          :: Z_TH_HIND, Z_RI_HIND, Z_CI_HIND, Z_RC_HINC, Z_CC_HINC
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))          :: ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))          :: ZCCT, ZCRT, ZCIT
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3),NMOD_CCN) :: ZNFT, ZNAT
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3),NMOD_IFN) :: ZIFT, ZINT
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3),NMOD_IMM) :: ZNIT
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))          :: ZNHT
!
INTEGER :: JL
!-------------------------------------------------------------------------------
!
ZTHT(:,:,:) = PTHT(:,:,:)
ZRVT(:,:,:) = PRVT(:,:,:)
ZRCT(:,:,:) = PRCT(:,:,:)
ZCCT(:,:,:) = PCCT(:,:,:)
ZRRT(:,:,:) = PRRT(:,:,:)
ZCRT(:,:,:) = PCRT(:,:,:)
ZRIT(:,:,:) = PRIT(:,:,:)
ZCIT(:,:,:) = PCIT(:,:,:)
ZRST(:,:,:) = PRST(:,:,:)
ZRGT(:,:,:) = PRGT(:,:,:)
ZNFT(:,:,:,:) = PNFT(:,:,:,:)
ZNAT(:,:,:,:) = PNAT(:,:,:,:)
ZIFT(:,:,:,:) = PIFT(:,:,:,:)
ZINT(:,:,:,:) = PINT(:,:,:,:)
ZNIT(:,:,:,:) = PNIT(:,:,:,:)
ZNHT(:,:,:) = PNHT(:,:,:)
!
!-------------------------------------------------------------------------------
!
IF (LWARM .AND. LACTI .AND. NMOD_CCN.GE.1) THEN
   CALL LIMA_CCN_ACTIVATION (PTSTEP, TPFILE, OCLOSE_OUT,                 &
                             PRHODREF, PEXNREF, PPABST, PT, PTM, PW_NU,   &
                             ZTHT, ZRVT, ZRCT, ZCCT, ZRRT, ZNFT, ZNAT)
   PTHT(:,:,:) = ZTHT(:,:,:)
   PRVT(:,:,:) = ZRVT(:,:,:)
   PRCT(:,:,:) = ZRCT(:,:,:)
   PCCT(:,:,:) = ZCCT(:,:,:)
   PNFT(:,:,:,:) = ZNFT(:,:,:,:)
   PNAT(:,:,:,:) = ZNAT(:,:,:,:)
!
! Call budgets
!
   IF (LBU_ENABLE) THEN
      IF (LBUDGET_TH) CALL BUDGET (PTHT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   4,                        'HENU_BU_RTH')
      IF (LBUDGET_RV) CALL BUDGET (PRVT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   6,                        'HENU_BU_RRV')
      IF (LBUDGET_RC) CALL BUDGET (PRCT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   7,                        'HENU_BU_RRC')
      IF (LBUDGET_SV) THEN
                      CALL BUDGET (PCCT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   12+NSV_LIMA_NC,           'HENU_BU_RSV')
            DO JL=1, NMOD_CCN
                      CALL BUDGET (PNFT(:,:,:,JL)*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_CCN_FREE+JL-1,'HENU_BU_RSV') 
            END DO
      END IF
   END IF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD .AND. LNUCL .AND. .NOT.LMEYERS .AND. NMOD_IFN.GE.1) THEN
   CALL LIMA_PHILLIPS_IFN_NUCLEATION (PTSTEP,                                           &
                                      PRHODREF, PEXNREF, PPABST,                        &
                                      ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,         &
                                      ZCCT, ZCIT, ZNAT, ZIFT, ZINT, ZNIT,               &
                                      Z_TH_HIND, Z_RI_HIND, Z_CI_HIND,                  &
                                      Z_RC_HINC, Z_CC_HINC                              )
!
! Call budgets
!
   IF (LBU_ENABLE) THEN
      IF (LBUDGET_TH) CALL BUDGET ((PTHT(:,:,:)+Z_TH_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,4,                        'HIND_BU_RTH')
      IF (LBUDGET_RV) CALL BUDGET ((PRVT(:,:,:)-Z_RI_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,6,                        'HIND_BU_RRV')
      IF (LBUDGET_RI) CALL BUDGET ((PRIT(:,:,:)+Z_RI_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,9,                        'HIND_BU_RRI')
      IF (LBUDGET_SV) THEN
                      CALL BUDGET ((PCIT(:,:,:)+Z_CI_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_NI,           'HIND_BU_RSV')
         IF (NMOD_IFN.GE.1) THEN
            DO JL=1, NMOD_IFN
                      CALL BUDGET ((ZIFT(:,:,:,JL))*PRHODJ(:,:,:)/PTSTEP,              12+NSV_LIMA_IFN_FREE+JL-1,'HIND_BU_RSV') 
            END DO
         END IF
      END IF
!
      IF (LBUDGET_TH) CALL BUDGET (ZTHT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,4,'HINC_BU_RTH')
      IF (LBUDGET_RC) CALL BUDGET (ZRCT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,7,'HINC_BU_RRC')
      IF (LBUDGET_RI) CALL BUDGET (ZRIT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,9,'HINC_BU_RRI')
      IF (LBUDGET_SV) THEN
         CALL BUDGET (ZCCT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_NC,'HINC_BU_RSV')
         CALL BUDGET (ZCIT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_NI,'HINC_BU_RSV')
      END IF
   END IF
!
   PTHT(:,:,:) = ZTHT(:,:,:)
   PRVT(:,:,:) = ZRVT(:,:,:)
   PRCT(:,:,:) = ZRCT(:,:,:)
   PCCT(:,:,:) = ZCCT(:,:,:)
   PRIT(:,:,:) = ZRIT(:,:,:)
   PCIT(:,:,:) = ZCIT(:,:,:)
   PNAT(:,:,:,:) = ZNAT(:,:,:,:)
   PIFT(:,:,:,:) = ZIFT(:,:,:,:)
   PINT(:,:,:,:) = ZINT(:,:,:,:)
   PNIT(:,:,:,:) = ZNIT(:,:,:,:)
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD .AND. LNUCL .AND. LMEYERS) THEN
   CALL LIMA_MEYERS_NUCLEATION (PTSTEP,                                     &
                                PRHODREF, PEXNREF, PPABST,                  &
                                ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,   &
                                ZCCT, ZCIT, ZINT,                           &
                                Z_TH_HIND, Z_RI_HIND, Z_CI_HIND,            &
                                Z_RC_HINC, Z_CC_HINC                        )
!
! Call budgets
!
   IF (LBU_ENABLE) THEN
      IF (LBUDGET_TH) CALL BUDGET ((PTHT(:,:,:)+Z_TH_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,4,             'HIND_BU_RTH')
      IF (LBUDGET_RV) CALL BUDGET ((PRVT(:,:,:)-Z_RI_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,6,             'HIND_BU_RRV')
      IF (LBUDGET_RI) CALL BUDGET ((PRIT(:,:,:)+Z_RI_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,9,             'HIND_BU_RRI')
      IF (LBUDGET_SV) CALL BUDGET ((PCIT(:,:,:)+Z_CI_HIND(:,:,:))*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_NI,'HIND_BU_RSV')
!
      IF (LBUDGET_TH) CALL BUDGET (ZTHT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,4,'HINC_BU_RTH')
      IF (LBUDGET_RC) CALL BUDGET (ZRCT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,7,'HINC_BU_RRC')
      IF (LBUDGET_RI) CALL BUDGET (ZRIT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,9,'HINC_BU_RRI')
      IF (LBUDGET_SV) THEN
         CALL BUDGET (ZCCT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_NC,'HINC_BU_RSV')
         CALL BUDGET (ZCIT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_NI,'HINC_BU_RSV')
      END IF
   END IF
!
PTHT(:,:,:) = ZTHT(:,:,:)
PRVT(:,:,:) = ZRVT(:,:,:)
PRCT(:,:,:) = ZRCT(:,:,:)
PCCT(:,:,:) = ZCCT(:,:,:)
PRIT(:,:,:) = ZRIT(:,:,:)
PCIT(:,:,:) = ZCIT(:,:,:)
PINT(:,:,:,:) = ZINT(:,:,:,:)
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD .AND. LNUCL .AND. LHHONI .AND. NMOD_CCN.GE.1) THEN
   CALL LIMA_CCN_HOM_FREEZING (PRHODREF, PEXNREF, PPABST, PW_NU,            &
                               ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,    &
                               ZCCT, ZCRT, ZCIT, ZNFT, ZNHT                 )
!
! Call budgets
!
   IF (LBU_ENABLE) THEN
     IF (LBUDGET_TH) CALL BUDGET (ZTHT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   4,                        'HONH_BU_RTH')
     IF (LBUDGET_RV) CALL BUDGET (ZRVT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   6,                        'HONH_BU_RRV')
     IF (LBUDGET_RI) CALL BUDGET (ZRIT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   9,                        'HONH_BU_RRI')
     IF (LBUDGET_SV) THEN
                     CALL BUDGET (ZCIT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,   12+NSV_LIMA_NI,           'HONH_BU_RSV')
          DO JL=1, NMOD_CCN
                     CALL BUDGET (ZNFT(:,:,:,JL)*PRHODJ(:,:,:)/PTSTEP,12+NSV_LIMA_CCN_FREE+JL-1,'HONH_BU_RSV') 
          END DO
     END IF
  END IF
!
PTHT(:,:,:) = ZTHT(:,:,:)
PRVT(:,:,:) = ZRVT(:,:,:)
PRIT(:,:,:) = ZRIT(:,:,:)
PCIT(:,:,:) = ZCIT(:,:,:)
PNHT(:,:,:) = ZNHT(:,:,:)
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_NUCLEATION_PROCS
