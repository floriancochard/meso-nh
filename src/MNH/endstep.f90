!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /srv/cvsroot/MNH-VX-Y-Z/src/MNH/endstep.f90,v $ $Revision: 1.2.2.2.2.2.16.1.2.5 $ $Date: 2014/04/22 14:31:38 $
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_ENDSTEP
!     ###################
!
INTERFACE
!
      SUBROUTINE ENDSTEP        (PTSTEP,KRR,KSV,KTCOUNT,KMI,               &
                                 HUVW_ADV_SCHEME,HTEMP_SCHEME, PRHODJ,     &
                                 PUS,PVS,PWS,PDRYMASSS,                    &
                                 PTHS,PRS,PTKES,PSVS,                      &
                                 PLSUS,PLSVS,PLSWS,                        &
                                 PLSTHS,PLSRVS,PLSZWSS,                    &
                                 PLBXUS,PLBXVS,PLBXWS,                     &
                                 PLBXTHS,PLBXRS,PLBXTKES,PLBXSVS,          &
                                 PLBYUS,PLBYVS,PLBYWS,                     &
                                 PLBYTHS,PLBYRS,PLBYTKES,PLBYSVS,          &
                                 PUM,PVM,PWM,PZWS,                         &
                                 PUT,PVT,PWT,PPABST,PDRYMASST,             &
                                 PTHT,PRT,PTHM,PRCM,PPABSM,PTKET,PSVT,     &
                                 PLSUM,PLSVM,PLSWM,                        &
                                 PLSTHM,PLSRVM,PLSZWSM,                    &
                                 PLBXUM,PLBXVM,PLBXWM,                     &
                                 PLBXTHM,PLBXRM,PLBXTKEM,PLBXSVM,          &
                                 PLBYUM,PLBYVM,PLBYWM,                     &
                                 PLBYTHM,PLBYRM,PLBYTKEM,PLBYSVM           )
!
REAL,                     INTENT(IN) :: PTSTEP        !  Time step
INTEGER,                  INTENT(IN) :: KRR           !  Number of water var.
INTEGER,                  INTENT(IN) :: KSV           !  Number of scal. var.
INTEGER,                  INTENT(IN) :: KTCOUNT       !  Temporal loop COUNTer
INTEGER,                  INTENT(IN) :: KMI           !  Model index
CHARACTER(LEN=6),         INTENT(IN) :: HUVW_ADV_SCHEME ! advection scheme for wind
CHARACTER(LEN=4),         INTENT(IN) :: HTEMP_SCHEME  ! Temporal scheme
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ        ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUS,PVS,PWS,   & ! 
                                        PTHS,PTKES       ! variables at 
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRS,PSVS         !    t+dt
!
REAL,                     INTENT(IN) :: PDRYMASSS           !   Md source
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLSUS,PLSVS,PLSWS,& ! Large Scale 
                                        PLSTHS,PLSRVS       ! fields tendencies
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PLSZWSS               ! Large Scale fields tendencies
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBXUS,PLBXVS,PLBXWS,  &  !
                                        PLBXTHS,PLBXTKES          ! LBX tendancy 
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PLBXRS,PLBXSVS            ! 
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBYUS,PLBYVS,PLBYWS,&    !     
                                        PLBYTHS,PLBYTKES          ! LBY tendancy 
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PLBYRS,PLBYSVS            !  
!  
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PUM,PVM,PWM! Variables at t-dt
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PUT,PVT,PWT,PPABST,PTHT,&!
                                         PTKET              ! Variables at
REAL, DIMENSION(:,:,:,:),INTENT(INOUT):: PRT,PSVT                 !     t
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PTHM, PRCM,PPABSM  ! Variables at t-Dt
REAL,                    INTENT(INOUT):: PDRYMASST                !
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLSUM,PLSVM,PLSWM,& ! Large Scale fields
                                         PLSTHM,PLSRVM       !     at t-dt
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PLSZWSM  ! Large Scale fields at t-dt
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PLBXUM,PLBXVM,PLBXWM,   & ! 
                                          PLBXTHM,PLBXTKEM          ! LBX fields
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PLBXRM,PLBXSVM            !
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PLBYUM,PLBYVM,PLBYWM,   & ! 
                                          PLBYTHM,PLBYTKEM          ! LBY fields
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PLBYRM,PLBYSVM            ! 
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PZWS                  ! significant wave height
!
END SUBROUTINE ENDSTEP
!
END INTERFACE
!
END MODULE MODI_ENDSTEP
!
!
!
!     ######################################################################
      SUBROUTINE ENDSTEP        (PTSTEP,KRR,KSV,KTCOUNT,KMI,               &
                                 HUVW_ADV_SCHEME,HTEMP_SCHEME, PRHODJ,     &
                                 PUS,PVS,PWS,PDRYMASSS,                    &
                                 PTHS,PRS,PTKES,PSVS,                      &
                                 PLSUS,PLSVS,PLSWS,                        &
                                 PLSTHS,PLSRVS,PLSZWSS,                    &
                                 PLBXUS,PLBXVS,PLBXWS,                     &
                                 PLBXTHS,PLBXRS,PLBXTKES,PLBXSVS,          &
                                 PLBYUS,PLBYVS,PLBYWS,                     &
                                 PLBYTHS,PLBYRS,PLBYTKES,PLBYSVS,          &
                                 PUM,PVM,PWM,PZWS,                         &
                                 PUT,PVT,PWT,PPABST,PDRYMASST,             &
                                 PTHT,PRT,PTHM,PRCM,PPABSM,PTKET,PSVT,     &
                                 PLSUM,PLSVM,PLSWM,                        &
                                 PLSTHM,PLSRVM,PLSZWSM,                    &
                                 PLBXUM,PLBXVM,PLBXWM,                     &
                                 PLBXTHM,PLBXRM,PLBXTKEM,PLBXSVM,          &
                                 PLBYUM,PLBYVM,PLBYWM,                     &
                                 PLBYTHM,PLBYRM,PLBYTKEM,PLBYSVM           )
!     ######################################################################
!
!!****  *ENDSTEP* - temporal advance and asselin filter for all variables
!!        (replaces the previous endstep_dyn and endstep_scalar subroutines)
!!
!!    PURPOSE
!!    -------
!!
!!    The purpose of ENDSTEP is to apply the asselin filter, perform
!!    the time advance and thereby finalize the time step.
!
!
!!**  METHOD
!!    ------
!!    
!!    The filtered values of the prognostic variables at t is obtained
!!    by linear combination of variables at t-dt, t, and t+dt.
!!    This value is put into the array containing the t-dt value.
!!    To perform the time swapping, the t+dt values are put into the arrays 
!!    containing the t values.
!!
!!    In case of cold start (first time step), indicated by the value 'START'
!!    of CCONF in module MODD_CONF, a simple time advance is performed.
!!
!!    The swapping for the absolute pressure function is only a copy of time t in
!!    time (t-dt).
!!
!!    Temporal advances of large scale, lateral boundarie and SST fields
!!    are also made in this subroutine.
!!
!!    The different sources terms are stored for the budget computations.
!!    
!!    EXTERNAL
!!    --------
!!      BUDGET      : Stores the different budget components
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    MODULE MODD_DYN containing XASSELIN
!!    MODULE MODD_CONF containing CCONF
!!    MODULE MODD_CTURB containing XTKEMIN, XEPSMIN
!!    MODULE MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         NBUPROCCTR   : process counter used for each budget variable
!!         NBUTSHIFT    : temporal shift for budgets writing
!!
!!    REFERENCE
!!    ---------
!!    Book2 of documentation
!!
!!    AUTHOR
!!    ------
!!    P. Bougeault  Meteo France
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    original     22/06/94
!!    corrections  01/09/94 (J. P. Lafore)
!!     "           07/11/94 (J.Stein)   pressure function swapping
!!    update       03/01/94 (J. P. Lafore) Total mass of dry air Md evolution
!!                 20/03/95 (J.Stein )   remove R from the historical variables
!!                                      + switch for TKE unused
!!                 01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                 30/08/95 (J.Stein)    remove the positivity control and 
!!                        correct the bug for PRM and PSVM for the cold start
!!                 16/10/95 (J. Stein)     change the budget calls 
!!                 12/10/96 (J. Stein)     add the SRC temporal evolution
!!                 20/12/96 (J.-P. Pinty)  update the CALL BUDGET
!!                 03/09/96 (J. P. Lafore) temporal advance of LS scalar fields
!!                 22/06/97 (J. Stein)     add the absolute pressure
!!                 13/03/97 (J. P. Lafore) add "surfacic" LS fields
!!                 24/09/97 (V. Masson)    positive values for ls fields
!!                 10/01/98 (J. Stein)     use the LB fields
!!                 20/04/98 (P. Josse)     temporal evolution of SST
!!                 18/09/98 (P. Jabouille) merge endstep_dyn and endstep_scalar
!!                 08/12/00 (P. Jabouille) minimum values for hydrometeors
!!                 22/06/01 (P. Jabouille) use XSVMIN
!!                 06/11/02 (V. Masson)    update the budget calls
!!                 01/2004  (V. Masson)  surface externalization
!!                 05/2006                Remove KEPS
!!                 10/2006  (Maric, Lac)  modification for PPM schemes
!!                 10/2009  (C.Lac)       Correction on FIT temporal scheme for variables
!!                                         advected with PPM
!!                 04/2013  (C.Lac)       FIT for all the variables     
!!                 04/2014  (C.Lac)       Check on the positivity of PSVT
!!                 J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!                 02/2019  (S. Bielli)  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODD_DYN 
USE MODD_CONF
USE MODD_CTURB
USE MODD_GRID_n
USE MODD_BUDGET
USE MODD_NSV, ONLY : XSVMIN, NSV_CHEMBEG, NSV_CHEMEND, &
                     NSV_AERBEG, NSV_AEREND,&
                     NSV_DSTBEG, NSV_DSTEND,&
                     NSV_SNWBEG, NSV_SNWEND
USE MODD_CH_AEROSOL, ONLY : LORILAM
USE MODD_DUST,       ONLY : LDUST
USE MODD_PARAM_C2R2, ONLY : LACTIT
USE MODD_PARAM_LIMA, ONLY : LACTIT_LIMA=>LACTIT
USE MODD_LBC_n, ONLY : CLBCX, CLBCY
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODI_BUDGET
USE MODI_SHUMAN
!
USE MODE_ll
!
IMPLICIT NONE
!
!*      0.1  DECLARATIONS OF ARGUMENTS
!
!
REAL,                     INTENT(IN) :: PTSTEP        !  Time step
INTEGER,                  INTENT(IN) :: KRR           !  Number of water var.
INTEGER,                  INTENT(IN) :: KSV           !  Number of scal. var.
INTEGER,                  INTENT(IN) :: KTCOUNT       !  Temporal loop COUNTer
INTEGER,                  INTENT(IN) :: KMI           !  Model index
CHARACTER(LEN=6),         INTENT(IN) :: HUVW_ADV_SCHEME ! advection scheme for wind
CHARACTER(LEN=4),         INTENT(IN) :: HTEMP_SCHEME  ! Temporal scheme
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ        ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUS,PVS,PWS,   & ! 
                                        PTHS,PTKES       ! variables at 
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRS,PSVS         !    t+dt
!
REAL,                     INTENT(IN) :: PDRYMASSS           !   Md source
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLSUS,PLSVS,PLSWS,& !    Large Scale 
                                        PLSTHS,PLSRVS       ! fields tendencies
REAL, DIMENSION(:,:),     INTENT(IN) :: PLSZWSS               ! Large Scale fields tendencies
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBXUS,PLBXVS,PLBXWS,  &  !
                                        PLBXTHS,PLBXTKES          ! LBX tendancy 
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PLBXRS,PLBXSVS            ! 
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBYUS,PLBYVS,PLBYWS,&    !     
                                        PLBYTHS,PLBYTKES          ! LBY tendancy 
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PLBYRS,PLBYSVS            !  
!  
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PUM,PVM,PWM! Variables at t-dt
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PUT,PVT,PWT,PPABST,PTHT,&!
                                         PTKET              ! Variables at
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PTHM, PRCM, PPABSM  ! Variables at t-Dt
REAL, DIMENSION(:,:,:,:),INTENT(INOUT):: PRT,PSVT                 !     t
REAL,                    INTENT(INOUT):: PDRYMASST                ! 
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLSUM,PLSVM,PLSWM,& ! Large Scale fields
                                         PLSTHM,PLSRVM       !     at t-dt
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PLSZWSM  ! Large Scale fields at t-dt
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PLBXUM,PLBXVM,PLBXWM,   & ! 
                                          PLBXTHM,PLBXTKEM          ! LBX fields
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PLBXRM,PLBXSVM            !
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PLBYUM,PLBYVM,PLBYWM,   & ! 
                                          PLBYTHM,PLBYTKEM          ! LBY fields
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PLBYRM,PLBYSVM            !   
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PZWS                  ! significant wave height
!
!*      0.2  DECLARATIONS OF LOCAL VARIABLES
!
INTEGER:: JSV                  ! loop counters
INTEGER :: IKU
INTEGER :: IIB, IIE  ! index of first and last inner mass points along x
INTEGER :: IJB, IJE  ! index of first and last inner mass points along y
!
!------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKU=SIZE(XZHAT)
!
!*      1.   ASSELIN FILTER
!
IF ((HUVW_ADV_SCHEME(1:3)=='CEN').AND. (HTEMP_SCHEME == 'LEFR')) THEN
  IF( KTCOUNT /= 1 .OR. CCONF /= 'START' ) THEN
     PUM(:,:,:)=(1.-XASSELIN)*PUT(:,:,:)+0.5*XASSELIN*(PUM(:,:,:)+PUS(:,:,:))
     PVM(:,:,:)=(1.-XASSELIN)*PVT(:,:,:)+0.5*XASSELIN*(PVM(:,:,:)+PVS(:,:,:))
     PWM(:,:,:)=(1.-XASSELIN)*PWT(:,:,:)+0.5*XASSELIN*(PWM(:,:,:)+PWS(:,:,:))
  END IF
END IF

!*      1.   TEMPORAL ADVANCE OF PROGNOSTIC VARIABLES
!
PPABSM(:,:,:) = PPABST(:,:,:)
!
IF (LACTIT .OR. LACTIT_LIMA) THEN
   PTHM(:,:,:)   = PTHT(:,:,:)
   PRCM(:,:,:)   = PRT(:,:,:,2)
END IF

PUT(:,:,:)=PUS(:,:,:)
PVT(:,:,:)=PVS(:,:,:)
PWT(:,:,:)=PWS(:,:,:)
!
PDRYMASST = PDRYMASST + PTSTEP * PDRYMASSS
!
PTHT(:,:,:)=PTHS(:,:,:)
!
! Moisture
!
PRT(:,:,:,1:KRR)=PRS(:,:,:,1:KRR)
!
! Turbulence
!
IF (SIZE(PTKET,1) /= 0) PTKET(:,:,:)=PTKES(:,:,:)
!
! Other scalars
!
PSVT(:,:,:,1:KSV)=PSVS(:,:,:,1:KSV)
!
IF(LBLOWSNOW) THEN
   DO JSV=1,(NBLOWSNOW_2D)
       XSNWCANO(:,:,JSV) = XRSNWCANOS(:,:,JSV)
   END DO
!*         MINIMUM VALUE FOR BLOWING SNOW
!
   WHERE(XSNWCANO(:,:,:)<1.E-20)
       XSNWCANO(:,:,:)=0.
   END WHERE

   IF (SIZE(PSVT,4) > 1) THEN
     WHERE(PSVT(:,:,:,NSV_SNWBEG:NSV_SNWEND)<1.E-20)
       PSVT(:,:,:,NSV_SNWBEG:NSV_SNWEND)=0.
     END WHERE
   END IF
!
END IF
!
IF (LWEST_ll( ) .AND. CLBCX(1)=='OPEN') THEN
 DO JSV=1,KSV
   PSVT(IIB,:,:,JSV)=MAX(PSVT(IIB,:,:,JSV),XSVMIN(JSV))
   PSVT(IIB-1,:,:,JSV)=MAX(PSVT(IIB-1,:,:,JSV),XSVMIN(JSV))
 END DO
END IF
!
IF (LEAST_ll( ) .AND. CLBCX(2)=='OPEN') THEN
 DO JSV=1,KSV
   PSVT(IIE,:,:,JSV)=MAX(PSVT(IIE,:,:,JSV),XSVMIN(JSV))
   PSVT(IIE+1,:,:,JSV)=MAX(PSVT(IIE+1,:,:,JSV),XSVMIN(JSV))
 END DO
END IF
!
IF (LSOUTH_ll( ) .AND. CLBCY(1)=='OPEN') THEN
 DO JSV=1,KSV
   PSVT(:,IJB,:,JSV)=MAX(PSVT(:,IJB,:,JSV),XSVMIN(JSV))
   PSVT(:,IJB-1,:,JSV)=MAX(PSVT(:,IJB-1,:,JSV),XSVMIN(JSV))
 END DO
END IF
!
IF (LNORTH_ll( ) .AND. CLBCY(2)=='OPEN') THEN
 DO JSV=1,KSV
   PSVT(:,IJE,:,JSV)=MAX(PSVT(:,IJE,:,JSV),XSVMIN(JSV))
   PSVT(:,IJE+1,:,JSV)=MAX(PSVT(:,IJE+1,:,JSV),XSVMIN(JSV))
 END DO
END IF
!------------------------------------------------------------------------------
!
!*      4.   TEMPORAL ADVANCE OF THE LARGE SCALE FIELDS
!
!
IF (SIZE(PLSUS,1) /= 0) THEN
  PLSUM(:,:,:)  = PLSUM(:,:,:)  + PTSTEP * PLSUS(:,:,:)
  PLSVM(:,:,:)  = PLSVM(:,:,:)  + PTSTEP * PLSVS(:,:,:)
  PLSWM(:,:,:)  = PLSWM(:,:,:)  + PTSTEP * PLSWS(:,:,:)
END IF
!
IF (SIZE(PLSTHS,1) /= 0) THEN
  PLSTHM(:,:,:) = PLSTHM(:,:,:) + PTSTEP * PLSTHS(:,:,:)
ENDIF
!
IF (SIZE(PLSRVS,1) /= 0) THEN
  PLSRVM(:,:,:) = MAX( PLSRVM(:,:,:) + PTSTEP * PLSRVS(:,:,:) , 0.)
ENDIF

IF (SIZE(PLSZWSS,1) /= 0) THEN
  PLSZWSM(:,:) = MAX( PLSZWSM(:,:) + PTSTEP * PLSZWSS(:,:) , 0.)
  PZWS(:,:) = PLSZWSM(:,:)
ENDIF
!
!------------------------------------------------------------------------------
!
!*      5.   TEMPORAL ADVANCE OF THE LATERAL BOUNDARIES FIELDS
!
IF (SIZE(PLBXUS,1) /= 0) THEN
  PLBXUM(:,:,:)  = PLBXUM(:,:,:)  + PTSTEP * PLBXUS(:,:,:)
  PLBXVM(:,:,:)  = PLBXVM(:,:,:)  + PTSTEP * PLBXVS(:,:,:)
  PLBXWM(:,:,:)  = PLBXWM(:,:,:)  + PTSTEP * PLBXWS(:,:,:)
ENDIF
IF (SIZE(PLBYUS,1) /= 0) THEN
  PLBYUM(:,:,:)  = PLBYUM(:,:,:)  + PTSTEP * PLBYUS(:,:,:)
  PLBYVM(:,:,:)  = PLBYVM(:,:,:)  + PTSTEP * PLBYVS(:,:,:)
  PLBYWM(:,:,:)  = PLBYWM(:,:,:)  + PTSTEP * PLBYWS(:,:,:)
ENDIF
!
IF (SIZE(PLBXTHS,1) /= 0) THEN
  PLBXTHM(:,:,:) = PLBXTHM(:,:,:) + PTSTEP * PLBXTHS(:,:,:)
END IF
IF (SIZE(PLBYTHS,1) /= 0) THEN
  PLBYTHM(:,:,:) = PLBYTHM(:,:,:) + PTSTEP * PLBYTHS(:,:,:)
END IF
!
IF (SIZE(PLBXTKES,1) /= 0) THEN
  PLBXTKEM(:,:,:) = MAX( PLBXTKEM(:,:,:) + PTSTEP * PLBXTKES(:,:,:), XTKEMIN)
END IF
IF (SIZE(PLBYTKES,1) /= 0) THEN
  PLBYTKEM(:,:,:) = MAX( PLBYTKEM(:,:,:) + PTSTEP * PLBYTKES(:,:,:), XTKEMIN)
END IF
!
IF (SIZE(PLBXRS,1) /= 0) THEN
  PLBXRM(:,:,:,:) = MAX( PLBXRM(:,:,:,:) + PTSTEP * PLBXRS(:,:,:,:), 0.)
END IF
IF (SIZE(PLBYRS,1) /= 0) THEN
  PLBYRM(:,:,:,:) = MAX( PLBYRM(:,:,:,:) + PTSTEP * PLBYRS(:,:,:,:), 0.)
END IF
!
IF (SIZE(PLBXSVS,1) /= 0) THEN
  DO JSV = 1,KSV
    PLBXSVM(:,:,:,JSV) = MAX( PLBXSVM(:,:,:,JSV) + PTSTEP * PLBXSVS(:,:,:,JSV),XSVMIN(JSV))
  ENDDO
ENDIF
IF (SIZE(PLBYSVS,1) /= 0) THEN
  DO JSV = 1,KSV
    PLBYSVM(:,:,:,JSV) = MAX( PLBYSVM(:,:,:,JSV) + PTSTEP * PLBYSVS(:,:,:,JSV),XSVMIN(JSV))
  ENDDO
END IF
!
!------------------------------------------------------------------------------
!
!*      6.   MINIMUM VALUE FOR HYDROMETEORS
!
IF (SIZE(PRT,4) > 1) THEN
  WHERE(PRT(:,:,:,2:)<1.E-20)
    PRT(:,:,:,2:)=0.
  END WHERE
END IF
IF (SIZE(PLBXRM,4) > 1) THEN
  WHERE(PLBXRM(:,:,:,2:)<1.E-20)
    PLBXRM(:,:,:,2:)=0.
  END WHERE
END IF
IF (SIZE(PLBYRM,4) > 1) THEN
  WHERE(PLBYRM(:,:,:,2:)<1.E-20)
    PLBYRM(:,:,:,2:)=0.
  END WHERE
END IF
!
!------------------------------------------------------------------------------
!
!*      7.   MINIMUM VALUE FOR CHEMISTRY
!
IF ((SIZE(PLBXSVM,4) > NSV_CHEMEND-1).AND.(SIZE(PLBXSVM,1) /= 0))  THEN
  DO JSV=NSV_CHEMBEG, NSV_CHEMEND
    PLBXSVM(:,:,:,JSV) = MAX(PLBXSVM(:,:,:,JSV), XSVMIN(JSV))
  END DO
END IF
IF ((SIZE(PLBYSVM,4) > NSV_CHEMEND-1).AND.(SIZE(PLBYSVM,1) /= 0)) THEN
  DO JSV=NSV_CHEMBEG, NSV_CHEMEND
    PLBYSVM(:,:,:,JSV) = MAX(PLBYSVM(:,:,:,JSV), XSVMIN(JSV))
  END DO
END IF
!
!------------------------------------------------------------------------------
!
!*      8.   MINIMUM VALUE FOR AEROSOLS
!
IF (LORILAM) THEN
  IF ((SIZE(PLBXSVM,4) > NSV_AEREND-1).AND.(SIZE(PLBXSVM,1) /= 0))  THEN
    DO JSV=NSV_AERBEG, NSV_AEREND
      PLBXSVM(:,:,:,JSV) = MAX(PLBXSVM(:,:,:,JSV), XSVMIN(JSV))
    END DO
  END IF
  IF ((SIZE(PLBYSVM,4) > NSV_AEREND-1).AND.(SIZE(PLBYSVM,1) /= 0)) THEN
    DO JSV=NSV_AERBEG, NSV_AEREND
      PLBYSVM(:,:,:,JSV) = MAX(PLBYSVM(:,:,:,JSV), XSVMIN(JSV))
    END DO
  END IF
END IF
!
!------------------------------------------------------------------------------
!
!*      9.   MINIMUM VALUE FOR DUSTS
!
IF (LDUST) THEN
  IF ((SIZE(PLBXSVM,4) > NSV_DSTEND-1).AND.(SIZE(PLBXSVM,1) /= 0))  THEN
    DO JSV=NSV_DSTBEG, NSV_DSTEND
      PLBXSVM(:,:,:,JSV) = MAX(PLBXSVM(:,:,:,JSV), XSVMIN(JSV))
    END DO
  END IF
  IF ((SIZE(PLBYSVM,4) > NSV_DSTEND-1).AND.(SIZE(PLBYSVM,1) /= 0)) THEN
    DO JSV=NSV_DSTBEG, NSV_DSTEND
      PLBYSVM(:,:,:,JSV) = MAX(PLBYSVM(:,:,:,JSV), XSVMIN(JSV))
    END DO
  END IF
END IF
!
!------------------------------------------------------------------------------
!
!*      10.   STORAGE IN BUDGET ARRAYS
!
IF (LBU_ENABLE) THEN
  NBUPROCCTR(1:12+KSV)=3
  NBUCTR_ACTV(1:12+KSV)=3
!
  IF (LBUDGET_U)   CALL BUDGET (PUT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,1,'AVEF_BU_RU')
  IF (LBUDGET_V)   CALL BUDGET (PVT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,2,'AVEF_BU_RV')
  IF (LBUDGET_W)   CALL BUDGET (PWT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,3,'AVEF_BU_RW')
  IF (LBUDGET_TH)  CALL BUDGET (PTHT(:,:,:)*PRHODJ(:,:,:)/PTSTEP,4,'AVEF_BU_RTH')
  IF (LBUDGET_TKE) CALL BUDGET (PTKET(:,:,:)*PRHODJ(:,:,:)/PTSTEP,5,'AVEF_BU_RTKE')
  IF (LBUDGET_RV)  CALL BUDGET (PRT(:,:,:,1)*PRHODJ(:,:,:)/PTSTEP,6,'AVEF_BU_RRV')
  IF (LBUDGET_RC)  CALL BUDGET (PRT(:,:,:,2)*PRHODJ(:,:,:)/PTSTEP,7,'AVEF_BU_RRC')
  IF (LBUDGET_RR)  CALL BUDGET (PRT(:,:,:,3)*PRHODJ(:,:,:)/PTSTEP,8,'AVEF_BU_RRR')
  IF (LBUDGET_RI)  CALL BUDGET (PRT(:,:,:,4)*PRHODJ(:,:,:)/PTSTEP,9,'AVEF_BU_RRI')
  IF (LBUDGET_RS)  CALL BUDGET (PRT(:,:,:,5)*PRHODJ(:,:,:)/PTSTEP,10,'AVEF_BU_RRS')
  IF (LBUDGET_RG)  CALL BUDGET (PRT(:,:,:,6)*PRHODJ(:,:,:)/PTSTEP,11,'AVEF_BU_RRG')
  IF (LBUDGET_RH)  CALL BUDGET (PRT(:,:,:,7)*PRHODJ(:,:,:)/PTSTEP,12,'AVEF_BU_RRH')
  IF (LBUDGET_SV) THEN
    DO JSV=1,KSV
      CALL BUDGET (PSVT(:,:,:,JSV)*PRHODJ(:,:,:)/PTSTEP,12+JSV,'AVEF_BU_RSV')
    END DO
  END IF
!
  NBUPROCCTR(1:12+KSV)=2
  NBUCTR_ACTV(1:12+KSV)=2
!
  IF (LBUDGET_U)   CALL BUDGET (PUS*MXM(PRHODJ)/PTSTEP,1,'ENDF_BU_RU')
  IF (LBUDGET_V)   CALL BUDGET (PVS*MYM(PRHODJ)/PTSTEP,2,'ENDF_BU_RV')
  IF (LBUDGET_W)   CALL BUDGET (PWS*MZM(1,IKU,1,PRHODJ)/PTSTEP,3,'ENDF_BU_RW')
  IF (LBUDGET_TH)  CALL BUDGET (PTHS*PRHODJ/PTSTEP,4,'ENDF_BU_RTH')
  IF (LBUDGET_TKE) CALL BUDGET (PTKES*PRHODJ/PTSTEP,5,'ENDF_BU_RTKE')
  IF (LBUDGET_RV)  CALL BUDGET (PRS(:,:,:,1)*PRHODJ/PTSTEP,6,'ENDF_BU_RRV')
  IF (LBUDGET_RC)  CALL BUDGET (PRS(:,:,:,2)*PRHODJ/PTSTEP,7,'ENDF_BU_RRC')
  IF (LBUDGET_RR)  CALL BUDGET (PRS(:,:,:,3)*PRHODJ/PTSTEP,8,'ENDF_BU_RRR')
  IF (LBUDGET_RI)  CALL BUDGET (PRS(:,:,:,4)*PRHODJ/PTSTEP,9,'ENDF_BU_RRI')
  IF (LBUDGET_RS)  CALL BUDGET (PRS(:,:,:,5)*PRHODJ/PTSTEP,10,'ENDF_BU_RRS')
  IF (LBUDGET_RG)  CALL BUDGET (PRS(:,:,:,6)*PRHODJ/PTSTEP,11,'ENDF_BU_RRG')
  IF (LBUDGET_RH)  CALL BUDGET (PRS(:,:,:,7)*PRHODJ/PTSTEP,12,'ENDF_BU_RRH')
  IF (LBUDGET_SV) THEN
    DO JSV=1,KSV
      CALL BUDGET (PSVS(:,:,:,JSV)*PRHODJ/PTSTEP,JSV+12,'ENDF_BU_RSV')
    END DO
  END IF
END IF
!
!------------------------------------------------------------------------------
!
!*      11.  COMPUTATION OF PHASE VELOCITY
!            -----------------------------
!
!   It is temporarily set to a constant value
!
!------------------------------------------------------------------------------
!
!
END SUBROUTINE ENDSTEP
