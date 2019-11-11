!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_SEDIMENTATION_SPLIT_OLD
INTERFACE
SUBROUTINE ICE4_SEDIMENTATION_SPLIT_OLD(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, KSPLITR, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PINPRH, PRHT, PRHS, PFPR)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                      INTENT(IN)              :: KKL     !vert. levels type 1=MNH -1=ARO
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
INTEGER,                      INTENT(IN)              :: KSPLITR ! Number of small time step integration for  rain sedimendation
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT_OLD
END INTERFACE
END MODULE MODI_ICE4_SEDIMENTATION_SPLIT_OLD
SUBROUTINE ICE4_SEDIMENTATION_SPLIT_OLD(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, KSPLITR, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PINPRH, PRHT, PRHS, PFPR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the sedimentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
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
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODI_BUDGET
USE MODD_BUDGET
USE MODI_GAMMA
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                      INTENT(IN)              :: KKL     !vert. levels type 1=MNH -1=ARO
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
INTEGER,                      INTENT(IN)              :: KSPLITR ! Number of small time step integration for  rain sedimendation
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
    :: GSEDIM ! Test where to compute the SED processes
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT

REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D, & !  droplet condensation
                                                                     & ZRAY,   & ! Cloud Mean radius
                                                                     & ZLBC,   & ! XLBC weighted by sea fraction
                                                                     & ZFSEDC, &
                                                                     & ZPRCS,ZPRRS,ZPRIS,ZPRSS,ZPRGS,ZPRHS, &   ! Mixing ratios created during the time step
                                                                     & ZW, & ! work array
                                                                     & ZRCT, &
                                                                     & ZRRT, &
                                                                     & ZRIT, &
                                                                     & ZRST, &
                                                                     & ZRGT, &
                                                                     & ZRHT
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1) :: ZWSED        ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZCONC_TMP    ! Weighted concentration
REAL :: ZINVTSTEP
INTEGER :: ISEDIM ! ! Case number of sedimentation
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
INTEGER :: JJ, JK, JN, JL
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!
!        O. Initialization of for sedimentation
!
ZINVTSTEP=1./PTSTEP
ZTSPLITR=PTSTEP/FLOAT(KSPLITR)
IF (OSEDIC) PINPRC (:,:) = 0.
PINPRR (:,:) = 0.
PINPRI (:,:) = 0.
PINPRS (:,:) = 0.
PINPRG (:,:) = 0.
IF ( KRR == 7 ) PINPRH (:,:) = 0.
!
!*       1. Parameters for cloud sedimentation
!
IF (OSEDIC) THEN
  ZRAY(:,:,:)   = 0.
  ZCONC_TMP(:,:)=PSEA(:,:)*XCONC_SEA+(1.-PSEA(:,:))*XCONC_LAND

  DO JK=KKTB, KKTE
    ZLBC(:,:,JK)   = PSEA(:,:)*XLBC(2)+(1.-PSEA(:,:))*XLBC(1)
    ZFSEDC(:,:,JK) = (PSEA(:,:)*XFSEDC(2)+(1.-PSEA(:,:))*XFSEDC(1))
    ZFSEDC(:,:,JK) = MAX(MIN(XFSEDC(1),XFSEDC(2)),ZFSEDC(:,:,JK))
    ZCONC3D(:,:,JK)= (1.-PTOWN(:,:))*ZCONC_TMP(:,:)+PTOWN(:,:)*XCONC_URBAN
    ZRAY(:,:,JK)   = 0.5*((1.-PSEA(:,:))*GAMMA(XNUC+1.0/XALPHAC)/(GAMMA(XNUC)) + &
            PSEA(:,:)*GAMMA(XNUC2+1.0/XALPHAC2)/(GAMMA(XNUC2)))
  END DO
  ZRAY(:,:,:)      = MAX(1.,ZRAY(:,:,:))
  ZLBC(:,:,:)      = MAX(MIN(XLBC(1),XLBC(2)),ZLBC(:,:,:))
ENDIF
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately
!
! External tendecies
IF (OSEDIC) ZPRCS(:,:,:) = PRCS(:,:,:)-PRCT(:,:,:)* ZINVTSTEP
ZPRRS(:,:,:) = PRRS(:,:,:)-PRRT(:,:,:)* ZINVTSTEP
ZPRIS(:,:,:) = PRIS(:,:,:)-PRIT(:,:,:)* ZINVTSTEP
ZPRSS(:,:,:) = PRSS(:,:,:)-PRST(:,:,:)* ZINVTSTEP
ZPRGS(:,:,:) = PRGS(:,:,:)-PRGT(:,:,:)* ZINVTSTEP
IF ( KRR == 7 ) ZPRHS(:,:,:) = PRHS(:,:,:)-PRHT(:,:,:)* ZINVTSTEP
!
! mr values inside the time-splitting loop
ZRCT(:,:,:) = PRCT(:,:,:)
ZRRT(:,:,:) = PRRT(:,:,:)
ZRIT(:,:,:) = PRIT(:,:,:)
ZRST(:,:,:) = PRST(:,:,:)
ZRGT(:,:,:) = PRGT(:,:,:)
IF (KRR==7) ZRHT(:,:,:) = PRHT(:,:,:)
!
DO JK = KKTB , KKTE
  ZW(:,:,JK) =ZTSPLITR/(PRHODREF(:,:,JK)* PDZZ(:,:,JK))
END DO
!
DO JN = 1 , KSPLITR
  !We add part of the external tendencies
  IF (OSEDIC) ZRCT(:,:,:) = ZRCT(:,:,:) + ZPRCS(:,:,:)*ZTSPLITR
  ZRRT(:,:,:) = ZRRT(:,:,:) + ZPRRS(:,:,:)*ZTSPLITR
  ZRIT(:,:,:) = ZRIT(:,:,:) + ZPRIS(:,:,:)*ZTSPLITR
  ZRST(:,:,:) = ZRST(:,:,:) + ZPRSS(:,:,:)*ZTSPLITR
  ZRGT(:,:,:) = ZRGT(:,:,:) + ZPRGS(:,:,:)*ZTSPLITR
  IF (KRR==7) ZRHT(:,:,:) = ZRHT(:,:,:) + ZPRHS(:,:,:)*ZTSPLITR
  !
  !
  !*       2.1   for cloud
  !
  IF (OSEDIC) THEN
    GSEDIM(:,:,:)=.FALSE.
    GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                    ZRCT(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(2)
    ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                             &SIZE(I1),I1(:),I2(:),I3(:))
    CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                          &ISEDIM, GSEDIM, I1, I2, I3, &
                          &PRHODREF, ZW, PPABST, PTHT, PSEA, PTOWN, ZTSPLITR, PTSTEP, &
                          &2, &
                          &ZRCT, PRCS, ZWSED, &
                          &ZRAY, ZLBC, ZFSEDC, ZCONC3D)
    IF (PRESENT(PFPR)) THEN
      DO JK = KKTB , KKTE
        PFPR(:,:,JK,2)=ZWSED(:,:,JK)
      ENDDO
    ENDIF
    PINPRC(:,:) = PINPRC(:,:) + ZWSED(:,:,KKB) / XRHOLW / KSPLITR
  END IF
  !
  !*       2.2   for rain
  !
  GSEDIM(:,:,:)=.FALSE.
  GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  ZRRT(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(3)
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                          &ISEDIM, GSEDIM, I1, I2, I3, &
                          &PRHODREF, ZW, PPABST, PTHT, PSEA, PTOWN, ZTSPLITR, PTSTEP, &
                          &3, &
                          &ZRRT, PRRS, ZWSED)
  IF (PRESENT(PFPR)) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,3)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRR(:,:) = PINPRR(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  !
  !*       2.3   for pristine ice
  !
  GSEDIM(:,:,:)=.FALSE.
  GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  ZRIT(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(4)
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                          &ISEDIM, GSEDIM, I1, I2, I3, &
                          &PRHODREF, ZW, PPABST, PTHT, PSEA, PTOWN, ZTSPLITR, PTSTEP, &
                          &4, &
                          &ZRIT, PRIS, ZWSED)
  IF (PRESENT(PFPR)) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,4)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRI(:,:) = PINPRI(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  !
  !*       2.4   for aggregates/snow
  !
  GSEDIM(:,:,:)=.FALSE.
  GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  ZRST(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(5)
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                          &ISEDIM, GSEDIM, I1, I2, I3, &
                          &PRHODREF, ZW, PPABST, PTHT, PSEA, PTOWN, ZTSPLITR, PTSTEP, &
                          &5, &
                          &ZRST, PRSS, ZWSED)
  IF (PRESENT(PFPR)) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,5)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRS(:,:) = PINPRS(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  !
  !*       2.5   for graupeln
  !
  GSEDIM(:,:,:)=.FALSE.
  GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  ZRGT(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(6)
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                          &ISEDIM, GSEDIM, I1, I2, I3, &
                          &PRHODREF, ZW, PPABST, PTHT, PSEA, PTOWN, ZTSPLITR, PTSTEP, &
                          &6, &
                          &ZRGT, PRGS, ZWSED)
  IF (PRESENT(PFPR)) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,6)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRG(:,:) = PINPRG(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  !
  !*       2.6   for hail
  !
  IF ( KRR == 7 ) THEN
    GSEDIM(:,:,:)=.FALSE.
    GSEDIM(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                    ZRHT(KIB:KIE,KJB:KJE,KKTB:KKTE)>XRTMIN(7)
    ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                             &SIZE(I1),I1(:),I2(:),I3(:))
    CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                            &ISEDIM, GSEDIM, I1, I2, I3, &
                            &PRHODREF, ZW, PPABST, PTHT, PSEA, PTOWN, ZTSPLITR, PTSTEP, &
                            &7, &
                            &ZRHT, PRHS, ZWSED)
    IF (PRESENT(PFPR)) THEN
      DO JK = KKTB , KKTE
        PFPR(:,:,JK,7)=ZWSED(:,:,JK)
      ENDDO
    ENDIF
    PINPRH(:,:) = PINPRH(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  END IF
  !
END DO
!
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE INTERNAL_SEDIM_SPLI(KIT, KJT, KKT, KKL, &
                                &KSEDIM, LDSEDIM, I1, I2, I3, &
                                &PRHODREF, PTSORHODZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, PTOTAL_TSTEP, &
                                &KSPE, &
                                &PRXT, PRXS, PWSED, &
                                &PRAY, PLBC, PFSEDC, PCONC3D)
    !
    !*      0. DECLARATIONS
    !          ------------
    !
    USE MODD_RAIN_ICE_DESCR
    USE MODD_RAIN_ICE_PARAM
    !
    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments :
    !
    INTEGER, INTENT(IN) :: KIT, KJT, KKT, KKL
    INTEGER, INTENT(IN) :: KSEDIM
    LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: LDSEDIM
    INTEGER, DIMENSION(KSEDIM), INTENT(IN) :: I1, I2, I3
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF ! Reference density
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTSORHODZ ! TimeStep Over (Rhodref time delta Z)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST
    REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PSEA    ! Sea Mask
    REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PTOWN   ! Fraction that is town
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT
    REAL,                         INTENT(IN)              :: PTSTEP  ! small timestep
    REAL,                         INTENT(IN)              :: PTOTAL_TSTEP ! total timestep
    INTEGER,                      INTENT(IN)              :: KSPE ! 1 for rc, 2 for rr...
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXT ! mr of specy X
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXS !Tendency of the specy KSPE
    REAL, DIMENSION(KIT,KJT,0:KKT+1), INTENT(OUT)         :: PWSED ! sedimentation flux
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN), OPTIONAL :: PRAY, PLBC, PFSEDC, PCONC3D
    !
    !*       0.2  declaration of local variables
    !
    !
    INTEGER :: JK, JL, JI, JJ
    REAL :: ZINVTOTAL_TSTEP
    REAL :: ZZWLBDC, ZRAY, ZZT, ZZWLBDA, ZZCC
    REAL :: ZFSED, ZEXSED
    REAL, DIMENSION(KIT, KJT) :: ZMRCHANGE
    !
    !-------------------------------------------------------------------------------
    !
    !
    !*       1. Parameters for cloud sedimentation
    !
    !
    !*       2.    compute the fluxes
    !
    !
    ZINVTOTAL_TSTEP = 1./PTOTAL_TSTEP
    PWSED(:,:,:) = 0.
    IF(KSPE==2) THEN
      !******* for cloud
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        ZZWLBDC = PLBC(JI,JJ,JK) * PCONC3D(JI,JJ,JK) / &
                 (PRHODREF(JI,JJ,JK) * PRXT(JI,JJ,JK))
        ZZWLBDC = ZZWLBDC**XLBEXC
        ZRAY = PRAY(JI,JJ,JK) / ZZWLBDC
        ZZT = PTHT(JI,JJ,JK) * (PPABST(JI,JJ,JK)/XP00)**(XRD/XCPD)
        ZZWLBDA = 6.6E-8*(101325./PPABST(JI,JJ,JK))*(ZZT/293.15)
        ZZCC = XCC*(1.+1.26*ZZWLBDA/ZRAY)
        PWSED(JI, JJ, JK) = PRHODREF(JI,JJ,JK)**(-XCEXVT +1 ) *   &
                 ZZWLBDC**(-XDC)*ZZCC*PFSEDC(JI,JJ,JK) * PRXT(JI,JJ,JK)
      ENDDO
    ELSEIF(KSPE==4) THEN
      ! ******* for pristine ice
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        IF(PRXT(JI, JJ, JK) .GT. MAX(XRTMIN(4), 1.0E-7)) THEN
          PWSED(JI, JJ, JK) =  XFSEDI * PRXT(JI, JJ, JK) *  &
                              & PRHODREF(JI,JJ,JK)**(1.-XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ,JK)) )**XEXCSEDI
        ENDIF
      ENDDO
    ELSE
      ! ******* for other species
      IF(KSPE==3) THEN
        ZFSED=XFSEDR
        ZEXSED=XEXSEDR
      ELSEIF(KSPE==5) THEN
        ZFSED=XFSEDS
        ZEXSED=XEXSEDS
      ELSEIF(KSPE==6) THEN
        ZFSED=XFSEDG
        ZEXSED=XEXSEDG
      ELSEIF(KSPE==7) THEN
        ZFSED=XFSEDH
        ZEXSED=XEXSEDH
      ELSE
        WRITE(*,*) ' STOP'
        WRITE(*,*) ' NO SEDIMENTATION PARAMETER FOR KSPE==', KSPE
        CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_SEDIMENTATION_SPLIT_OLD','')
      ENDIF
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        PWSED(JI, JJ, JK) = ZFSED  * PRXT(JI, JJ, JK)**ZEXSED *   &
                                     PRHODREF(JI, JJ, JK)**(ZEXSED-XCEXVT)
      ENDDO
    ENDIF
    ZMRCHANGE(:,:) = 0.
    DO JK = KKTB , KKTE
       ZMRCHANGE(:,:) = PTSORHODZ(:,:,JK)*(PWSED(:,:,JK+KKL)-PWSED(:,:,JK))
       PRXT(:,:,JK) = PRXT(:,:,JK) + ZMRCHANGE(:,:)
       PRXS(:,:,JK) = PRXS(:,:,JK) + ZMRCHANGE(:,:) * ZINVTOTAL_TSTEP
    ENDDO
  END SUBROUTINE INTERNAL_SEDIM_SPLI
  !
  FUNCTION ICE4_SEDIMENTATION_SPLIT_COUNTJV(LTAB,KIT,KJT,KKT,KSIZE,I1,I2,I3) RESULT(IC)
  !
  !*      0. DECLARATIONS
  !          ------------
  !
  IMPLICIT NONE
  !
  !*       0.2  declaration of local variables
  !
  INTEGER, INTENT(IN) :: KIT,KJT,KKT,KSIZE
  LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)  :: LTAB ! Mask
  INTEGER, DIMENSION(KSIZE), INTENT(OUT) :: I1,I2,I3 ! Used to replace the COUNT and PACK
  INTEGER :: JI,JJ,JK,IC
  !
  !-------------------------------------------------------------------------------
  !
  IC = 0
  DO JK = 1,SIZE(LTAB,3)
    DO JJ = 1,SIZE(LTAB,2)
      DO JI = 1,SIZE(LTAB,1)
        IF( LTAB(JI,JJ,JK) ) THEN
          IC = IC +1
          I1(IC) = JI
          I2(IC) = JJ
          I3(IC) = JK
        END IF
      END DO
    END DO
  END DO
  !
  END FUNCTION ICE4_SEDIMENTATION_SPLIT_COUNTJV
  !
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT_OLD
