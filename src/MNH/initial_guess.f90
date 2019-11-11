!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_INITIAL_GUESS
!     #########################
!
INTERFACE
!
      SUBROUTINE INITIAL_GUESS ( KRR, KSV, KTCOUNT,PRHODJ, KMI, PTSTEP,        &
                         PRUS, PRVS, PRWS, PRTHS, PRRS, PRTKES, PRSVS,         &
                         PUT, PVT, PWT, PTHT, PRT, PTKET, PSVT )
!
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of Scalar Variables
INTEGER,                  INTENT(IN)  :: KTCOUNT ! Temporal loop COUNTer
                                                 ! (=1 at the segment beginning)
INTEGER,                  INTENT(IN)  :: KMI     ! Model index
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL,                     INTENT(IN)  :: PTSTEP !  timestep 
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PRRS, PRSVS              !  terms
!
! variables at time t (needed for PPM schemes)
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PUT, PVT, PWT
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHT, PTKET
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRT, PSVT
!
END SUBROUTINE INITIAL_GUESS
!
END INTERFACE
!
END MODULE MODI_INITIAL_GUESS 
!
!     #########################################################################
      SUBROUTINE INITIAL_GUESS ( KRR, KSV, KTCOUNT,PRHODJ, KMI, PTSTEP,        &
                         PRUS, PRVS, PRWS, PRTHS, PRRS, PRTKES, PRSVS,         &
                         PUT, PVT, PWT, PTHT, PRT, PTKET, PSVT )
!     #########################################################################
!
!!****  *INITIAL_GUESS * - routine to initialize the source terms
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to integrate the prognostic variables
!!    at t-dt into their respective source terms.
!!
!!**  METHOD
!!    ------
!!      The fields at t-dt divided by 2*TSTEP (1*TSTEP for the first time step
!!    in case of START configuration) are initializing the source term arrays.
!!      The different sources terms are initialized for the budget computations.
!!     
!!
!!    EXTERNAL
!!    --------
!!      MXM,MYM,MZM : Mean Shuman operators in the x,y,z directions
!!      BUDGET      : Stores the different budget components
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF   : contains configuration variable 
!!           CCONF :  Configuration of models
!!    MODULE MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         LBU_BEG      : logical for budget begnning
!!                       .TRUE. = budget begining
!!                       .FALSE. = no budget begining
!!         NBUPROCCTR   : process counter used for each budget variable
!!         Switches for budgets activations:
!!         
!!         LBU_RU       : logical for budget of RU (wind component along x)
!!                        .TRUE. = budget of RU         
!!                        .FALSE. = no budget of RU 
!!         LBU_RV       : logical for budget of RV (wind component along y)
!!                        .TRUE. = budget of RV         
!!                        .FALSE. = no budget of RV 
!!         LBU_RW        : logical for budget of RW (wind component along z)
!!                        .TRUE. = budget of RW         
!!                        .FALSE. = no budget of RW 
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
!!         LBU_RSV      : logical for budget of RSVx (scalar variable)
!!                        .TRUE. = budget of RSVx 
!!                        .FALSE. = no budget of RSVx
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INITIAL_GUESS )
!!
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/07/94 
!!                  20/03/95 (J.Stein) : remove R from the historical variables
!!                  01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                  16/10/95 (J. Stein)     change the budget calls 
!!                  19/12/96 (J.-P. Pinty)  update the budget calls 
!!                  06/11/02 (V. Masson)    update the budget calls 
!!                  20/05/06                Remove KEPS
!!                  10/09    (C.Lac)        FIT for variables advected with PPM
!!                  04/13    (C.Lac)        FIT for all variables 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF 
USE MODD_GRID_n
USE MODD_BUDGET
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
!
USE MODI_SHUMAN
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of Scalar Variables
INTEGER,                  INTENT(IN)  :: KTCOUNT ! Temporal loop COUNTer
                                                 ! (=1 at the segment beginning)
INTEGER,                  INTENT(IN)  :: KMI     ! Model index
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL,                     INTENT(IN)  :: PTSTEP !  timestep 
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PRRS, PRSVS  !  terms
!
!
! variables at time t (needed for PPM schemes)
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PUT, PVT, PWT
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHT, PTKET
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRT, PSVT
!
!*       0.2   declarations of local variables
!
INTEGER                               :: JRR, JSV
INTEGER                               :: IKU
REAL                                  :: ZINVTSTEP
!
!-------------------------------------------------------------------------------
!
IKU=SIZE(XZHAT)
!*       1.     COMPUTES THE INVERSE OF THE APPLICABLE TIMESTEP
!   	        -----------------------------------------------
!
ZINVTSTEP = 1./PTSTEP                          
!
!
!*       2.     COMPUTES THE FIRST SOURCE TERMS
!   	        -------------------------------
! 
! *** momentum
! forward-in-time time-marching scheme
PRUS = PUT * ZINVTSTEP * MXM(PRHODJ)
PRVS = PVT * ZINVTSTEP * MYM(PRHODJ)
PRWS = PWT * ZINVTSTEP * MZM(1,IKU,1,PRHODJ)
!
! *** meteorological variables
!
PRTHS(:,:,:) = PTHT(:,:,:) * ZINVTSTEP * PRHODJ(:,:,:)
IF (SIZE(PTKET,1) /= 0) THEN 
  PRTKES(:,:,:) = PTKET(:,:,:) * ZINVTSTEP * PRHODJ(:,:,:)
END IF
!
! Case with KRR moist variables 
DO JRR=1,KRR
  PRRS(:,:,:,JRR) = PRT(:,:,:,JRR) * ZINVTSTEP * PRHODJ(:,:,:) 
END DO
!
! *** passive tracers
!
! Case with KSV Scalar Variables
DO JSV=1,KSV
  PRSVS(:,:,:,JSV) = PSVT(:,:,:,JSV) * ZINVTSTEP * PRHODJ(:,:,:)
END DO
!
IF(LBLOWSNOW) THEN
  DO JSV=1,(NBLOWSNOW_2D)
    XRSNWCANOS(:,:,JSV) = XSNWCANO(:,:,JSV) * ZINVTSTEP * PRHODJ(:,:,1)
  END DO
END IF
!
IF (LBU_ENABLE) THEN
  IF (LBU_BEG) THEN
    NBUPROCCTR(:)=1
    NBUCTR_ACTV(:)=1
!
    IF (LBUDGET_U)   CALL BUDGET (PRUS,1,'INIF_BU_RU')
    IF (LBUDGET_V)   CALL BUDGET (PRVS,2,'INIF_BU_RV')
    IF (LBUDGET_W)   CALL BUDGET (PRWS,3,'INIF_BU_RW')
    IF (LBUDGET_TH)  CALL BUDGET (PRTHS,4,'INIF_BU_RTH')
    IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'INIF_BU_RTKE')
    IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'INIF_BU_RRV')
    IF (LBUDGET_RC)  CALL BUDGET (PRRS(:,:,:,2),7,'INIF_BU_RRC')
    IF (LBUDGET_RR)  CALL BUDGET (PRRS(:,:,:,3),8,'INIF_BU_RRR')
    IF (LBUDGET_RI)  CALL BUDGET (PRRS(:,:,:,4),9,'INIF_BU_RRI')
    IF (LBUDGET_RS)  CALL BUDGET (PRRS(:,:,:,5),10,'INIF_BU_RRS')
    IF (LBUDGET_RG)  CALL BUDGET (PRRS(:,:,:,6),11,'INIF_BU_RRG')
    IF (LBUDGET_RH)  CALL BUDGET (PRRS(:,:,:,7),12,'INIF_BU_RRH')
    DO JSV=1,KSV
      IF (LBUDGET_SV)  CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'INIF_BU_RSV')
    END DO
!
    NBUPROCCTR(:)=2
    NBUCTR_ACTV(:)=2
!
    IF (LBUDGET_U)   CALL BUDGET (PRUS,1,'ENDF_BU_RU')
    IF (LBUDGET_V)   CALL BUDGET (PRVS,2,'ENDF_BU_RV')
    IF (LBUDGET_W)   CALL BUDGET (PRWS,3,'ENDF_BU_RW')
    IF (LBUDGET_TH)  CALL BUDGET (PRTHS,4,'ENDF_BU_RTH')
    IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ENDF_BU_RTKE')
    IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'ENDF_BU_RRV')
    IF (LBUDGET_RC)  CALL BUDGET (PRRS(:,:,:,2),7,'ENDF_BU_RRC')
    IF (LBUDGET_RR)  CALL BUDGET (PRRS(:,:,:,3),8,'ENDF_BU_RRR')
    IF (LBUDGET_RI)  CALL BUDGET (PRRS(:,:,:,4),9,'ENDF_BU_RRI')
    IF (LBUDGET_RS)  CALL BUDGET (PRRS(:,:,:,5),10,'ENDF_BU_RRS')
    IF (LBUDGET_RG)  CALL BUDGET (PRRS(:,:,:,6),11,'ENDF_BU_RRG')
    IF (LBUDGET_RH)  CALL BUDGET (PRRS(:,:,:,7),12,'ENDF_BU_RRH')
    DO JSV=1,KSV
      IF (LBUDGET_SV)  CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ENDF_BU_RSV')
    END DO
!
    LBU_BEG=.FALSE.
  END IF    
!
  NBUPROCCTR(:)=4
  NBUCTR_ACTV(:)=4
!
!  stores the Asselin source term
!
  IF (LBUDGET_U)   CALL BUDGET (PRUS,1,'ASSE_BU_RU')
  IF (LBUDGET_V)   CALL BUDGET (PRVS,2,'ASSE_BU_RV')
  IF (LBUDGET_W)   CALL BUDGET (PRWS,3,'ASSE_BU_RW')
  IF (LBUDGET_TH)  CALL BUDGET (PRTHS,4,'ASSE_BU_RTH')
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ASSE_BU_RTKE')
  IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'ASSE_BU_RRV')
  IF (LBUDGET_RC)  CALL BUDGET (PRRS(:,:,:,2),7,'ASSE_BU_RRC')
  IF (LBUDGET_RR)  CALL BUDGET (PRRS(:,:,:,3),8,'ASSE_BU_RRR')
  IF (LBUDGET_RI)  CALL BUDGET (PRRS(:,:,:,4),9,'ASSE_BU_RRI')
  IF (LBUDGET_RS)  CALL BUDGET (PRRS(:,:,:,5),10,'ASSE_BU_RRS')
  IF (LBUDGET_RG)  CALL BUDGET (PRRS(:,:,:,6),11,'ASSE_BU_RRG')
  IF (LBUDGET_RH)  CALL BUDGET (PRRS(:,:,:,7),12,'ASSE_BU_RRH')
  DO JSV=1,KSV
    IF (LBUDGET_SV)  CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ASSE_BU_RSV')
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INITIAL_GUESS
