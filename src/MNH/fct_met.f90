!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/06/06 12:00:47
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_FCT_MET
!     ######################
!
INTERFACE
  SUBROUTINE FCT_MET  (HLBCX, HLBCY, KRR,                             &
                       PTSTEP, PRHODJ, PTHM, PRM, PTKEM,              &
                       PTHT, PRT, PTKET,                              &
                       PRUCT, PRVCT, PRWCT,                           &
                       PRTHS, PRRS, PRTKES                            )
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Double Time step 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM 
                                                  ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT
                                                  ! Contravariant component of
                                                  ! momentum 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
                                                  ! Sources terms 
!
END SUBROUTINE FCT_MET 
!
END INTERFACE
!
END MODULE MODI_FCT_MET 
!
!
!
! #####################################################################
  SUBROUTINE FCT_MET  (HLBCX, HLBCY, KRR,                         &
                       PTSTEP, PRHODJ, PTHM, PRM, PTKEM,          &
                       PTHT, PRT, PTKET,                          &
                       PRUCT, PRVCT, PRWCT,                       &
                       PRTHS, PRRS, PRTKES                        )
! #####################################################################
!
!!****  *FCT_MET * - routine to call the Flux-Corrected Transport for
!!                           meteorological scalars
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to call the FLUX-CORRected routine
!!    for meteorological scalars variables.
!!
!!**  METHOD
!!    ------
!!     The Flux-Corrected Transport method correct the fluxes using a limiting
!!    factor. This method ensures that the advection scheme is definite 
!!    positive. A minimum value of the scalar (MIN) equal to 0 is used for
!!    the positiveness of the scheme.
!!
!!    EXTERNAL
!!    --------
!!     Functions MXM,MYM,MZM : computes the averages along three directins
!!     Functions DXF,DYF,DZF  : computes the finite differences  
!!     Subroutine FLUX_CORR : corrects the advective fluxes 
!!     Subroutine BUDGET    : stores the sources for budget purposes
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_BUDGET : LBU_R* ( individual budget switches) 
!!                    CBUTYPE, NBUMOD
!!    REFERENCE
!!    ---------
!!      Book1 and book2 ( routine ADVECTION )
!!
!!    AUTHOR
!!    ------
!!      J. Vila & JP. Lafore    *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/10/94 
!!      Stein       27/06/96  add the budgets
!!      Pinty       20/12/96  update the budgets
!!      Lafore      27/03/98  call to DFLUX_CORR
!!      Lafore      01/04/98  FCT only on total water (rv+rc+ri) and
!!                            precipitating hydrometeors,
!!                            remove 4D flux local arrays
!!      Stein       20/04/99  remove KMI from the list of argument of DFLUX_CORR
!!      Masson      07/11/02  update the budgets
!!      Lac         24/04/06  split meteorological and passive tracer routines
!!                     05/06  Remove KEPS
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODI_SHUMAN
USE MODI_DFLUX_CORR
USE MODI_BUDGET
USE MODD_GRID_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Double Time step 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM 
                                                  ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT
                                                  ! Contravariant component of
                                                  ! momentum 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
                                                  ! Sources terms 
!
!
!*       0.2   declarations of local variables
!
INTEGER                                  :: JRR
                                                  ! Loop index 
!
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3))  &
                     :: ZFX  ,ZFY  ,ZFZ    ! Advective flux components for each
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) ::  &  
              ZRTFX,ZRTFY,ZRTFZ  ! variables and for total water (rv+rc+ri)
!
REAL                                                    :: ZMINR,ZMINTKE
                                                  ! Absolute minimum values of
                                                  ! water substances, TKE
INTEGER :: IKU                                                  
!-------------------------------------------------------------------------------
!
IKU=SIZE(XZHAT)
!*       1.   FLUX-CORRECTED TRANSPORT ADVECTION SCHEME for the HMET group
!
!
!*       1.1 Temperature:        ---> advected by a CEN scheme
                                       ! NB! It is not necessary to make the
                                       ! flux correction since temperature is
                                       ! always positive.
!
!
  PRTHS(:,:,:) = PRTHS(:,:,:)                            & 
                - DXF(PRUCT(:,:,:)*MXM (PTHT(:,:,:)))     
  IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVX_BU_RTH')
!
  PRTHS(:,:,:) = PRTHS(:,:,:)                            &
                - DYF(PRVCT(:,:,:)*MYM (PTHT(:,:,:)))     
  IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVY_BU_RTH')
!
  PRTHS(:,:,:) = PRTHS(:,:,:)                            &
                - DZF(1,IKU,1,PRWCT(:,:,:)*MZM (1,IKU,1,PTHT(:,:,:)))     
  IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVZ_BU_RTH')
!
!*       1.2 No condensation case: Vapor ---> advected by a FCT scheme
!
  ZMINR=0. ! Absolute minimum water substances 
!
  IF (KRR == 1) THEN
    CALL DFLUX_CORR (HLBCX, HLBCY, PTSTEP, ZMINR, PRHODJ,                &
                      PRM(:,:,:,1), PRT(:,:,:,1),                        &
                      PRUCT, PRVCT, PRWCT,                               &
                      ZFX(:,:,:), ZFY(:,:,:), ZFZ(:,:,:)                 )
!
    PRRS(:,:,:,1) = PRRS(:,:,:,1) - DXF(ZFX(:,:,:))
    IF (LBUDGET_RV)                          &
                              CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVX_BU_RRV')
!
    PRRS(:,:,:,1) = PRRS(:,:,:,1) - DYF(ZFY(:,:,:))
    IF (LBUDGET_RV)                          & 
                              CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVY_BU_RRV')
!
    PRRS(:,:,:,1) = PRRS(:,:,:,1) - DZF(1,IKU,1,ZFZ(:,:,:))
    IF (LBUDGET_RV)                          &
                              CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVZ_BU_RRV')
  END IF 
!
!*       1.3 No ice case:          rv+rc ---> advected by the FCT scheme
!                                     rc ---> advected by the CEN scheme
!
  IF (KRR == 2 .OR. KRR == 3 ) THEN
    CALL DFLUX_CORR (HLBCX,HLBCY,PTSTEP,ZMINR,PRHODJ,                      &
                      PRM(:,:,:,1)+PRM(:,:,:,2),PRT(:,:,:,1)+PRT(:,:,:,2), &
                      PRUCT, PRVCT, PRWCT,                                 &
                      ZRTFX(:,:,:),ZRTFY(:,:,:),ZRTFZ(:,:,:)               )
!
  ZFX(:,:,:) = PRUCT(:,:,:) * MXM (PRT(:,:,:,2))    !
  ZFY(:,:,:) = PRVCT(:,:,:) * MYM (PRT(:,:,:,2))    ! CENtred scheme for rc
  ZFZ(:,:,:) = PRWCT(:,:,:) * MZM (1,IKU,1,PRT(:,:,:,2))    !
!
  ZRTFX(:,:,:) = ZRTFX(:,:,:) - ZFX(:,:,:)        !
  ZRTFY(:,:,:) = ZRTFY(:,:,:) - ZFY(:,:,:)        !  rv fluxes deduction
  ZRTFZ(:,:,:) = ZRTFZ(:,:,:) - ZFZ(:,:,:)        !
!
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - DXF(ZRTFX(:,:,:))
  PRRS(:,:,:,2) = PRRS(:,:,:,2) - DXF(  ZFX(:,:,:)) 
  IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVX_BU_RRV')
  IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVX_BU_RRC')
!
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - DYF(ZRTFY(:,:,:))
  PRRS(:,:,:,2) = PRRS(:,:,:,2) - DYF(  ZFY(:,:,:))
  IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVY_BU_RRV')
  IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVY_BU_RRC')
!
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - DZF(1,IKU,1,ZRTFZ(:,:,:))
  PRRS(:,:,:,2) = PRRS(:,:,:,2) - DZF(1,IKU,1,  ZFZ(:,:,:))
  IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVZ_BU_RRV')
  IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVZ_BU_RRC')
!
  END IF
!
!*       1.4 Ice case:          rv+rc+ri ---> advected by the FCT scheme
!                                     rc ---> advected by the CEN scheme
!                                     ri ---> advected by the CEN scheme
!
  IF ( KRR >= 4 ) THEN
    CALL DFLUX_CORR (HLBCX,HLBCY,PTSTEP,ZMINR,PRHODJ,                      &
                      PRM(:,:,:,1)+PRM(:,:,:,2)+PRM(:,:,:,4),              &
                      PRT(:,:,:,1)+PRT(:,:,:,2)+PRT(:,:,:,4),              &
                      PRUCT, PRVCT, PRWCT,                                 &
                      ZRTFX(:,:,:),ZRTFY(:,:,:),ZRTFZ(:,:,:)               )
!
  ZFX(:,:,:) = PRUCT(:,:,:) * MXM (PRT(:,:,:,2))    !
  ZFY(:,:,:) = PRVCT(:,:,:) * MYM (PRT(:,:,:,2))    ! CENtred scheme for rc
  ZFZ(:,:,:) = PRWCT(:,:,:) * MZM (1,IKU,1,PRT(:,:,:,2))    !
!
  ZRTFX(:,:,:) = ZRTFX(:,:,:) - ZFX(:,:,:)        !
  ZRTFY(:,:,:) = ZRTFY(:,:,:) - ZFY(:,:,:)        ! rv+ri fluxes deduction
  ZRTFZ(:,:,:) = ZRTFZ(:,:,:) - ZFZ(:,:,:)        !
!
  PRRS(:,:,:,2) = PRRS(:,:,:,2) - DXF(  ZFX(:,:,:))
  IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVX_BU_RRC')
!
  PRRS(:,:,:,2) = PRRS(:,:,:,2) - DYF(  ZFY(:,:,:))
  IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVY_BU_RRC')
!
  PRRS(:,:,:,2) = PRRS(:,:,:,2) - DZF(1,IKU,1,  ZFZ(:,:,:))
  IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVZ_BU_RRC')
!
!
  ZFX(:,:,:) = PRUCT(:,:,:) * MXM (PRT(:,:,:,4))    !
  ZFY(:,:,:) = PRVCT(:,:,:) * MYM (PRT(:,:,:,4))    ! CENtred scheme for ri
  ZFZ(:,:,:) = PRWCT(:,:,:) * MZM (1,IKU,1,PRT(:,:,:,4))    !
!
  ZRTFX(:,:,:) = ZRTFX(:,:,:) - ZFX(:,:,:)        !
  ZRTFY(:,:,:) = ZRTFY(:,:,:) - ZFY(:,:,:)        !  rv fluxes deduction
  ZRTFZ(:,:,:) = ZRTFZ(:,:,:) - ZFZ(:,:,:)        !
!
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - DXF(ZRTFX(:,:,:))
  PRRS(:,:,:,4) = PRRS(:,:,:,4) - DXF(  ZFX(:,:,:)) 
  IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVX_BU_RRV')
  IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVX_BU_RRI')
!
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - DYF(ZRTFY(:,:,:))
  PRRS(:,:,:,4) = PRRS(:,:,:,4) - DYF(  ZFY(:,:,:)) 
  IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVY_BU_RRV')
  IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVY_BU_RRI')
!
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - DZF(1,IKU,1,ZRTFZ(:,:,:))
  PRRS(:,:,:,4) = PRRS(:,:,:,4) - DZF(1,IKU,1,  ZFZ(:,:,:))
  IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVZ_BU_RRV')
  IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVZ_BU_RRI')
!
  END IF
!
!*       1.5 Rain case:                rr ---> advected by the FCT scheme
!
  IF ( KRR >= 3 ) THEN
    CALL DFLUX_CORR (HLBCX, HLBCY,                                       &
                      PTSTEP, ZMINR, PRHODJ, PRM(:,:,:,3), PRT(:,:,:,3), &
                      PRUCT, PRVCT, PRWCT,                               &
                      ZFX(:,:,:), ZFY(:,:,:), ZFZ(:,:,:)                 )
!
  PRRS(:,:,:,3) = PRRS(:,:,:,3) - DXF(  ZFX(:,:,:)) 
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVX_BU_RRR')
!
  PRRS(:,:,:,3) = PRRS(:,:,:,3) - DYF(  ZFY(:,:,:)) 
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVY_BU_RRR')
!
  PRRS(:,:,:,3) = PRRS(:,:,:,3) - DZF(1,IKU,1,  ZFZ(:,:,:)) 
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVZ_BU_RRR')
!
  END IF
!
!*       1.6 Other hydrometeors: rs,rg,rh ---> advected by the FCT scheme
!
  DO JRR = 5, KRR
    CALL DFLUX_CORR (HLBCX, HLBCY,                                           &
                      PTSTEP, ZMINR, PRHODJ, PRM(:,:,:,JRR), PRT(:,:,:,JRR), &
                      PRUCT, PRVCT, PRWCT,                                   &
                      ZFX(:,:,:), ZFY(:,:,:), ZFZ(:,:,:)                     )
!
    PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) - DXF(ZFX(:,:,:)) 
    IF (JRR==5.AND.LBUDGET_RS) &
                                    CALL BUDGET (PRRS(:,:,:,5),10,'ADVX_BU_RRS')
    IF (JRR==6.AND.LBUDGET_RG) &
                                    CALL BUDGET (PRRS(:,:,:,6),11,'ADVX_BU_RRG')
    IF (JRR==7.AND.LBUDGET_RH) &
                                    CALL BUDGET (PRRS(:,:,:,7),12,'ADVX_BU_RRH')
!
    PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) - DYF(ZFY(:,:,:))
    IF (JRR==5.AND.LBUDGET_RS) &
                                    CALL BUDGET (PRRS(:,:,:,5),10,'ADVY_BU_RRS')
    IF (JRR==6.AND.LBUDGET_RG) &
                                    CALL BUDGET (PRRS(:,:,:,6),11,'ADVY_BU_RRG')
    IF (JRR==7.AND.LBUDGET_RH) &
                                    CALL BUDGET (PRRS(:,:,:,7),12,'ADVY_BU_RRH')
!
    PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) - DZF(1,IKU,1,ZFZ(:,:,:))
    IF (JRR==5.AND.LBUDGET_RS) &
                                    CALL BUDGET (PRRS(:,:,:,5),10,'ADVZ_BU_RRS')
    IF (JRR==6.AND.LBUDGET_RG) &
                                    CALL BUDGET (PRRS(:,:,:,6),11,'ADVZ_BU_RRG')
    IF (JRR==7.AND.LBUDGET_RH) &
                                    CALL BUDGET (PRRS(:,:,:,7),12,'ADVZ_BU_RRH')
!
  END DO
!
!*       1.6 TKE                          ---> advected by the FCT scheme
!
  IF (SIZE(PTKET,1) /= 0) THEN
!
    ZMINTKE=0.                    ! Absolute minimum TKE 
!
    CALL DFLUX_CORR    (HLBCX, HLBCY,                            &
                         PTSTEP, ZMINTKE, PRHODJ, PTKEM,PTKET,   &
                         PRUCT, PRVCT, PRWCT,                    &
                         ZFX(:,:,:), ZFY(:,:,:), ZFZ(:,:,:)      )
!
    PRTKES(:,:,:) = PRTKES(:,:,:) - DXF(ZFX(:,:,:))     
    IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVX_BU_RTKE')
!
    PRTKES(:,:,:) = PRTKES(:,:,:) - DYF(ZFY(:,:,:))   
    IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVY_BU_RTKE')
!
    PRTKES(:,:,:) = PRTKES(:,:,:) - DZF(1,IKU,1,ZFZ(:,:,:)) 
    IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVZ_BU_RTKE')
!
  END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FCT_MET
