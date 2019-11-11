!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_ICE_ADJUST
!     ######################
!
INTERFACE
!
      SUBROUTINE ICE_ADJUST (KKA, KKU, KKL, KRR, HFRAC_ICE,                    &
                             HBUNAME, OSUBG_COND, OSIGMAS,                     &
                             PTSTEP, PSIGQSAT,                                 &
                             PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, PZZ,     &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR ,  &
                             PRR, PRI, PRIS, PRS, PRG,                         &
                             PRH, POUT_RV, POUT_RC, POUT_RI, POUT_TH           )
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER*1,              INTENT(IN)    :: HFRAC_ICE
CHARACTER*4,              INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL                                 :: OSIGMAS  ! Switch for Sigma_s: 
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ     ! height of model layer
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXN    ! Exner function
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRIS   ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRR    ! Rain water m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRI    ! Cloud ice  m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRS    ! Aggregate  m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRG    ! Graupel    m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
!
!
END SUBROUTINE ICE_ADJUST
!
END INTERFACE
!
END MODULE MODI_ICE_ADJUST

!     ##########################################################################
      SUBROUTINE ICE_ADJUST (KKA, KKU, KKL, KRR, HFRAC_ICE,                    &
                             HBUNAME, OSUBG_COND, OSIGMAS,                     &
                             PTSTEP, PSIGQSAT,                                 &
                             PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, PZZ,     &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR ,  &
                             PRR, PRI, PRIS, PRS, PRG, PRH,                    &
                             POUT_RV, POUT_RC, POUT_RI, POUT_TH                )
!     #########################################################################
!
!!****  *ICE_ADJUST* -  compute the ajustment of water vapor in mixed-phase 
!!                      clouds
!!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute the fast microphysical sources
!!    through a saturation ajustement procedure in case of mixed-phase clouds.
!!
!!
!!**  METHOD
!!    ------
!!    Langlois, Tellus, 1973 for the cloudless version.
!!    When cloud water is taken into account, refer to book 1 of the
!!    documentation.
!!
!!     
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!         XP00               ! Reference pressure
!!         XMD,XMV            ! Molar mass of dry air and molar mass of vapor
!!         XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
!!         XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
!!         XCL                ! Cl (liquid)
!!         XCI                ! Ci (ice)
!!         XTT                ! Triple point temperature
!!         XLVTT              ! Vaporization heat constant
!!         XLSTT              ! Sublimation  heat constant
!!         XALPW,XBETAW,XGAMW ! Constants for saturation vapor over liquid
!!                            !  pressure  function 
!!         XALPI,XBETAI,XGAMI ! Constants for saturation vapor over ice
!!                            !  pressure  function 
!!      Module  MODD_CONF 
!!         CCONF
!!      Module MODD_BUDGET:
!!         NBUMOD 
!!         CBUTYPE
!!         NBUPROCCTR 
!!         LBU_RTH    
!!         LBU_RRV  
!!         LBU_RRC  
!!         LBU_RRI  
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 and Book2 of documentation ( routine ICE_ADJUST )
!!      Langlois, Tellus, 1973
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty    * Laboratoire d'Aerologie*
!!   
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/12/96 
!!      M. Tomasini 27/11/00 Change CND and DEP fct of the T instead of rc and ri
!!                           Avoid the sub- and super-saturation before the ajustment
!!                           Avoid rc>0 if T<T00 before the ajustment
!!      P Bechtold 12/02/02  change subgrid condensation
!!      JP Pinty   29/11/02  add ICE2 and IC4 cases
!!      (P. Jabouille) 27/05/04 safety test for case where esw/i(T)> pabs (~Z>40km)
!!      J.Pergaud and S.Malardel Add EDKF case
!!      S. Riette ice for EDKF
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      2016-07 S. Riette: adjustement is now realized on state variables (PRV, PRC, PRI, PTH)
!!                         whereas tendencies are still applied on S variables.
!!                         This modification allows to call ice_adjust on T variable
!!                         or to call it on S variables
!!      2016-11 S. Riette: all-or-nothing adjustment now uses condensation
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CONF
USE MODD_BUDGET
!
USE MODI_CONDENSATION
USE MODI_BUDGET
USE MODE_FMWRIT
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KKA  !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER*1,              INTENT(IN)    :: HFRAC_ICE
CHARACTER*4,              INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL                                 :: OSIGMAS  ! Switch for Sigma_s: 
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ     ! height of model layer
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXN    ! Exner function
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
!
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRR  ! Rain water m.r. to adjust
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRI  ! Cloud ice  m.r. to adjust
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRS  ! Aggregate  m.r. to adjust
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRG  ! Graupel    m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
!
!*       0.2   Declarations of local variables :
!
!
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                         :: ZT,   &  ! adjusted temperature
                   ZRV, ZRC, ZRI, &  ! adjusted state
                            ZCPH, &  ! guess of the CPh for the mixing
                            ZLV,  &  ! guess of the Lv at t+1
                            ZLS,  &  ! guess of the Ls at t+1
                          ZW1,ZW2    ! Work arrays for intermediate fields
!
INTEGER             :: IIU,IJU,IKU! dimensions of dummy arrays
INTEGER             :: IIB,IJB    ! Horz index values of the first inner mass points
INTEGER             :: IIE,IJE    ! Horz index values of the last inner mass points
INTEGER             :: IKB        ! K index value of the first inner mass point
INTEGER             :: IKE        ! K index value of the last inner mass point
INTEGER             :: JITER,ITERMAX ! iterative loop for first order adjustment
!
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) :: ZSIGS,ZSRCS
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
IIU = SIZE(PEXNREF,1)
IJU = SIZE(PEXNREF,2)
IKU = SIZE(PEXNREF,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
!
ITERMAX=1
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!
!    beginning of the iterative loop (to compute the adjusted state)
ZRV(:,:,:)=PRV(:,:,:)
ZRC(:,:,:)=PRC(:,:,:)
ZRI(:,:,:)=PRI(:,:,:)
ZT(:,:,:)=PTH(:,:,:) * PEXN(:,:,:)
!
DO JITER =1,ITERMAX
  !
  !*       2.3    compute the latent heat of vaporization Lv(T*) at t+1
  !                   and the latent heat of sublimation  Ls(T*) at t+1
  !
  ZLV(:,:,:) = XLVTT + ( XCPV - XCL ) * ( ZT(:,:,:) -XTT )
  ZLS(:,:,:) = XLSTT + ( XCPV - XCI ) * ( ZT(:,:,:) -XTT )
  !
  !*       2.4    compute the specific heat for moist air (Cph) at t+1
  !
  IF     ( KRR == 7 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV * ZRV(:,:,:)                             &
                       + XCL  * (ZRC(:,:,:) + PRR(:,:,:))             &
                       + XCI  * (ZRI(:,:,:) + PRS(:,:,:) + PRG(:,:,:) + PRH(:,:,:))
  ELSE IF( KRR == 6 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV * ZRV(:,:,:)                             &
                       + XCL  * (ZRC(:,:,:) + PRR(:,:,:))             &
                       + XCI  * (ZRI(:,:,:) + PRS(:,:,:) + PRG(:,:,:))
  ELSE IF( KRR == 5 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV * ZRV(:,:,:)                             &
                       + XCL  * (ZRC(:,:,:) + PRR(:,:,:))             &
                       + XCI  * (ZRI(:,:,:) + PRS(:,:,:))
  ELSE IF( KRR == 3 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV * ZRV(:,:,:)               &
                       + XCL  * (ZRC(:,:,:) + PRR(:,:,:))
  ELSE IF( KRR == 2 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV * ZRV(:,:,:) &
                       + XCL  * ZRC(:,:,:)
  END IF
  !
  IF ( OSUBG_COND ) THEN
    !
    !*       3.     SUBGRID CONDENSATION SCHEME
    !               ---------------------------
    !
    !   PSRC= s'rci'/Sigma_s^2
    !   ZT, ZRV, ZRC and ZRI are INOUT
    CALL CONDENSATION(IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, KKL,                  &
         HFRAC_ICE,                                                                      &
         PPABST, PZZ, ZT, ZRV, ZRC, ZRI, PRS, PRG, PSIGS, PMFCONV, PCLDFR, PSRCS, .TRUE., OSIGMAS, &
         PSIGQSAT, PLV=ZLV, PLS=ZLS, PCPH=ZCPH )
  ELSE
    !
    !*       4.     ALL OR NOTHING CONDENSATION SCHEME
    !                            FOR MIXED-PHASE CLOUD
    !               -----------------------------------------------
    !
    !
    !   ZT, ZRV, ZRC and ZRI are INOUT
    !
    !CALL ADJUST_LANGLOIS(IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, KKL,                &
    !                     PPABST, ZT, ZRV, ZRC, ZRI, ZLV, ZLS, ZCPH) HFRAC_ICE must be implemented in Langlois before using it again
    ZSIGS=0.
    ZSRCS=0.
    CALL CONDENSATION(IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, KKL,        &
         HFRAC_ICE,                                                            &
         PPABST, PZZ, ZT, ZRV, ZRC, ZRI, PRS, PRG, ZSIGS, PMFCONV, PCLDFR,     &
         ZSRCS, .TRUE., OSIGMAS=.TRUE.,                                        &
         PSIGQSAT=0., PLV=ZLV, PLS=ZLS, PCPH=ZCPH )
  END IF
ENDDO         ! end of the iterative loop
!
!*       5.     COMPUTE THE SOURCES AND STORES THE CLOUD FRACTION
!               -------------------------------------------------
!
!
!*       5.0    compute the variation of mixing ratio
!
                                                      !         Rc - Rc*
ZW1(:,:,:) = (ZRC(:,:,:) - PRC(:,:,:)) / PTSTEP       ! Pcon = ----------
                                                      !         2 Delta t

ZW2(:,:,:) = (ZRI(:,:,:) - PRI(:,:,:)) / PTSTEP       ! idem ZW1 but for Ri
!
!*       5.1    compute the sources
!
WHERE( ZW1(:,:,:) < 0.0 )
  ZW1(:,:,:) = MAX ( ZW1(:,:,:), -PRCS(:,:,:) )
ELSEWHERE
  ZW1(:,:,:) = MIN ( ZW1(:,:,:),  PRVS(:,:,:) )
END WHERE
PRVS(:,:,:) = PRVS(:,:,:) - ZW1(:,:,:)
PRCS(:,:,:) = PRCS(:,:,:) + ZW1(:,:,:)
PTHS(:,:,:) = PTHS(:,:,:) +        &
                ZW1(:,:,:) * ZLV(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
!
WHERE( ZW2(:,:,:) < 0.0 )
  ZW2(:,:,:) = MAX ( ZW2(:,:,:), -PRIS(:,:,:) )
ELSEWHERE
  ZW2(:,:,:) = MIN ( ZW2(:,:,:),  PRVS(:,:,:) )
END WHERE
PRVS(:,:,:) = PRVS(:,:,:) - ZW2(:,:,:)
PRIS(:,:,:) = PRIS(:,:,:) + ZW2(:,:,:)
PTHS(:,:,:) = PTHS(:,:,:) +        &
              ZW2(:,:,:) * ZLS(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
!
!
!*       5.2    compute the cloud fraction PCLDFR
!
IF ( .NOT. OSUBG_COND ) THEN
  WHERE (PRCS(:,:,:) + PRIS(:,:,:) > 1.E-12 / PTSTEP)
    PCLDFR(:,:,:)  = 1.
  ELSEWHERE
    PCLDFR(:,:,:)  = 0. 
  ENDWHERE 
  IF ( SIZE(PSRCS,3) /= 0 ) THEN
    PSRCS(:,:,:) = PCLDFR(:,:,:) 
  END IF
ELSE
  !We limit PRC_MF+PRI_MF to PRVS*PTSTEP to avoid negative humidity
  ZW1(:,:,:)=PRC_MF(:,:,:)/PTSTEP
  ZW2(:,:,:)=PRI_MF(:,:,:)/PTSTEP
  WHERE(ZW1(:,:,:)+ZW2(:,:,:)>PRVS(:,:,:))
    ZW1(:,:,:)=ZW1(:,:,:)*PRVS(:,:,:)/(ZW1(:,:,:)+ZW2(:,:,:))
    ZW2(:,:,:)=PRVS(:,:,:)-ZW1(:,:,:)
  ENDWHERE
  PCLDFR(:,:,:)=MIN(1.,PCLDFR(:,:,:)+PCF_MF(:,:,:))
  PRCS(:,:,:)=PRCS(:,:,:)+ZW1(:,:,:)
  PRIS(:,:,:)=PRIS(:,:,:)+ZW2(:,:,:)
  PRVS(:,:,:)=PRVS(:,:,:)-(ZW1(:,:,:)+ZW2(:,:,:))
  PTHS(:,:,:) = PTHS(:,:,:) + &
                (ZW1 * ZLV(:,:,:) + ZW2 * ZLS(:,:,:)) / ZCPH(:,:,:)     &
                /  PEXNREF(:,:,:)
  IF(PRESENT(POUT_RV) .OR. PRESENT(POUT_RC) .OR. &
    &PRESENT(POUT_RI) .OR. PRESENT(POUT_TH)) THEN
    ZW1(:,:,:)=PRC_MF(:,:,:)
    ZW2(:,:,:)=PRI_MF(:,:,:)
    WHERE(ZW1(:,:,:)+ZW2(:,:,:)>ZRV(:,:,:))
      ZW1(:,:,:)=ZW1(:,:,:)*ZRV(:,:,:)/(ZW1(:,:,:)+ZW2(:,:,:))
      ZW2(:,:,:)=ZRV(:,:,:)-ZW1(:,:,:)
    ENDWHERE
    ZRC(:,:,:)=ZRC(:,:,:)+ZW1(:,:,:)
    ZRI(:,:,:)=ZRI(:,:,:)+ZW2(:,:,:)
    ZRV(:,:,:)=ZRV(:,:,:)-(ZW1(:,:,:)+ZW2(:,:,:))
    ZT(:,:,:) = ZT(:,:,:) + &
                (ZW1 * ZLV(:,:,:) + ZW2 * ZLS(:,:,:)) / ZCPH(:,:,:)
    IF(PRESENT(POUT_RV)) POUT_RV=ZRV
    IF(PRESENT(POUT_RC)) POUT_RC=ZRC
    IF(PRESENT(POUT_RI)) POUT_RI=ZRI
    IF(PRESENT(POUT_TH)) POUT_TH=ZT / PEXN(:,:,:)
  ENDIF
ENDIF
!
!
!
!*       6.  STORE THE BUDGET TERMS
!            ----------------------
!
IF (LBUDGET_RV) CALL BUDGET (PRVS(:,:,:) * PRHODJ(:,:,:),6,HBUNAME//'_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRCS(:,:,:) * PRHODJ(:,:,:),7,HBUNAME//'_BU_RRC')
IF (LBUDGET_RI) CALL BUDGET (PRIS(:,:,:) * PRHODJ(:,:,:),9,HBUNAME//'_BU_RRI')
IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:) * PRHODJ(:,:,:),4,HBUNAME//'_BU_RTH')
!
!------------------------------------------------------------------------------
!
!
END SUBROUTINE ICE_ADJUST 
