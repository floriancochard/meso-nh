!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/05/22 19:02:00
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_ADVECMET 
!     #######################
INTERFACE
      SUBROUTINE ADVECMET ( KRR, PTHT, PRT, PTKET,     &
                            PRUCT, PRVCT, PRWCT,       &
                            PRTHS, PRRS, PRTKES        )
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
                                                  ! Sources terms 
END SUBROUTINE ADVECMET
!
END INTERFACE
!
END MODULE MODI_ADVECMET 
!
!
!
!     ######################################################################
      SUBROUTINE ADVECMET ( KRR, PTHT, PRT, PTKET,        &
                            PRUCT, PRVCT, PRWCT,          &
                            PRTHS, PRRS, PRTKES           )
!     ######################################################################
!
!!****  *ADVECMET * - routine to compute the advection tendancies of the
!!                       meterological scalar fields.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the total advection 
!!    tendancies of the meteorological scalar fields, written in flux form.
!!      The advection velocity is taken as the contravariant form of 
!!    the momentum for extension to non-cartesian geometry and 
!!    conformal projection cases.
!!
!!
!!**  METHOD
!!    ------
!!      The left and right lateral EXTernal zones, have been previously
!!    prepared in routine LBC_S, to avoid particular cases close to the
!!    Lateral Boundaries in this routine.
!!      The Shuman functions are used to write the mean and finite 
!!    differences operators.
!!      The different sources terms are stored for the budget
!!    computations.
!!
!!    EXTERNAL
!!    --------
!!      MXM,MYM,MZM : Shuman functions (mean operators)
!!      DXM,DYM,DZM : Shuman functions (finite differences operators)
!!      BUDGET      : Stores the different budget components
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH        
!!                        .FALSE. = no budget of RTH
!!         LBU_RTKE     : logical for budget of RTKE (turbulent kinetic energy)
!!                        .TRUE. = budget of RTKE       
!!                        .FALSE. = no budget of RTKE
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV 
!!                        .FALSE. = no budget of RRV 
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC 
!!                        .FALSE. = no budget of RRC 
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR 
!!                        .FALSE. = no budget of RRR 
!!         LBU_RRI      : logical for budget of RRI (ice)
!!                        .TRUE. = budget of RRI 
!!                        .FALSE. = no budget of RRI 
!!         LBU_RRS      : logical for budget of RRS (snow)
!!                        .TRUE. = budget of RRS 
!!                        .FALSE. = no budget of RRS 
!!         LBU_RRG      : logical for budget of RRG (graupel)
!!                        .TRUE. = budget of RRG 
!!                        .FALSE. = no budget of RRG 
!!         LBU_RRH      : logical for budget of RRH (hail)
!!                        .TRUE. = budget of RRH 
!!                        .FALSE. = no budget of RRH 
!!        
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine ADVECMET )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!  	J.-P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/07/94 
!!      Corrections 06/09/94 (J.-P. Lafore)
!!                  16/03/95 (J. Stein)     remove R from the historical var.
!!                  01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                  16/10/95 (J. Stein)     change the budget calls 
!!                  19/12/96 (J.-P. Pinty)  update the budget calls 
!!                  07/11/02 (V. Masson)    update the budget calls
!!                  24/04/06 (C.Lac)        Split meteorological scalar and passive
!!                                          tracer routines
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_GRID_n
!
USE MODI_SHUMAN
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
                                                  ! Sources terms 
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR           ! Loop index for  moist variables
INTEGER :: IKU
!
!  
!-------------------------------------------------------------------------------
!
IKU=SIZE(XZHAT)
!*       1.     COMPUTES THE ADVECTIVE TENDENCIES
!     	        ---------------------------------
!
                                        ! Thermodynamical variable
PRTHS(:,:,:) = PRTHS(:,:,:)                            &
              -DXF( PRUCT(:,:,:) * MXM (PTHT(:,:,:)) ) 
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVX_BU_RTH')
!
PRTHS(:,:,:) = PRTHS(:,:,:)                            &
              -DYF( PRVCT(:,:,:) * MYM (PTHT(:,:,:)) ) 
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVY_BU_RTH')
!
PRTHS(:,:,:) = PRTHS(:,:,:)                            &
              -DZF(1,IKU,1, PRWCT(:,:,:) * MZM (1,IKU,1,PTHT(:,:,:)) )
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVZ_BU_RTH')
!
                                        ! Case with KRR moist variables 
DO JRR=1,KRR
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)                            &
                   -DXF( PRUCT(:,:,:) * MXM (PRT(:,:,:,JRR)) ) 
END DO
!
IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVX_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVX_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVX_BU_RRR')
IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVX_BU_RRI')
IF (LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'ADVX_BU_RRS')
IF (LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'ADVX_BU_RRG')
IF (LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'ADVX_BU_RRH')
!
DO JRR=1,KRR
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)                            &
                   -DYF( PRVCT(:,:,:) * MYM (PRT(:,:,:,JRR)) ) 
END DO
!
IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVY_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVY_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVY_BU_RRR')
IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVY_BU_RRI')
IF (LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'ADVY_BU_RRS')
IF (LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'ADVY_BU_RRG')
IF (LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'ADVY_BU_RRH')
!
DO JRR=1,KRR
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)                            &
                   -DZF(1,IKU,1, PRWCT(:,:,:) * MZM (1,IKU,1,PRT(:,:,:,JRR)) )
END DO
!
IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVZ_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVZ_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVZ_BU_RRR')
IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVZ_BU_RRI')
IF (LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'ADVZ_BU_RRS')
IF (LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'ADVZ_BU_RRG')
IF (LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'ADVZ_BU_RRH')
!
                                        ! TKE variable
IF (SIZE(PTKET,1) /= 0) THEN
  PRTKES(:,:,:) = PRTKES(:,:,:)                            &
                 -DXF( PRUCT(:,:,:) * MXM (PTKET(:,:,:)) ) 
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVX_BU_RTKE')
!
  PRTKES(:,:,:) = PRTKES(:,:,:)                            &
                 -DYF( PRVCT(:,:,:) * MYM (PTKET(:,:,:)) ) 
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVY_BU_RTKE')
!
   PRTKES(:,:,:) = PRTKES(:,:,:)                           &
                 -DZF(1,IKU,1, PRWCT(:,:,:) * MZM (1,IKU,1,PTKET(:,:,:)) )
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVZ_BU_RTKE')
END IF
!
! 
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECMET
