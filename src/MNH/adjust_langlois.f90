!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######################
      MODULE MODI_ADJUST_LANGLOIS
!     ######################
!
INTERFACE
!
      SUBROUTINE ADJUST_LANGLOIS(KIU, KJU, KKU, KIB, KIE, KJB, KJE, KKB, KKE, KKL, &
                                 PPABS, PT, PRV, PRC, PRI, &
                                 PLV, PLS, PCPH )
INTEGER,                      INTENT(IN)    :: KIU    ! horizontal dimension in x
INTEGER,                      INTENT(IN)    :: KJU    ! horizontal dimension in y
INTEGER,                      INTENT(IN)    :: KKU    ! vertical dimension
INTEGER,                      INTENT(IN)    :: KIB    ! value of the first point in x
INTEGER,                      INTENT(IN)    :: KIE    ! value of the last  point in x
INTEGER,                      INTENT(IN)    :: KJB    ! value of the first point in y
INTEGER,                      INTENT(IN)    :: KJE    ! value of the last  point in y
INTEGER,                      INTENT(IN)    :: KKB    ! value of the first point in z
INTEGER,                      INTENT(IN)    :: KKE    ! value of the last  point in z
INTEGER,                      INTENT(IN)    :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PLV
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PLS
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PCPH
!
!
END SUBROUTINE ADJUST_LANGLOIS
!
END INTERFACE
!
END MODULE MODI_ADJUST_LANGLOIS
!     ##########################################################################
      SUBROUTINE ADJUST_LANGLOIS(KIU, KJU, KKU, KIB, KIE, KJB, KJE, KKB, KKE, KKL, &
                                 PPABS, PT, PRV, PRC, PRI, &
                                 PLV, PLS, PCPH )
!     #########################################################################
!
!!****  *ADJUST_LABGLOIS* -  compute the ajustment of water vapor in mixed-phase
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
!!      S. Riette from the source code from J.-P. Pinty    * Laboratoire d'Aerologie*
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/07/2016 from the splitting of ice_adjust
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
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                      INTENT(IN)    :: KIU    ! horizontal dimension in x
INTEGER,                      INTENT(IN)    :: KJU    ! horizontal dimension in y
INTEGER,                      INTENT(IN)    :: KKU    ! vertical dimension
INTEGER,                      INTENT(IN)    :: KIB    ! value of the first point in x
INTEGER,                      INTENT(IN)    :: KIE    ! value of the last  point in x
INTEGER,                      INTENT(IN)    :: KJB    ! value of the first point in y
INTEGER,                      INTENT(IN)    :: KJE    ! value of the last  point in y
INTEGER,                      INTENT(IN)    :: KKB    ! value of the first point in z
INTEGER,                      INTENT(IN)    :: KKE    ! value of the last  point in z
INTEGER,                      INTENT(IN)    :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PLV
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PLS
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PCPH
!
!*       0.2   Declarations of local variables :
!
!
REAL  :: ZEPS  ! Mv/Md
REAL  :: ZT00,ZT0   ! Min and max temperature for the mixed phase liquid and solid water
                    ! for the coeff CND of the barycentric mixing ratio
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: &
      ZW1,ZW2,ZW3,ZW4,ZW5,ZW6,ZW7,&  ! Work arrays for intermediate fields
                            ZCND     ! CND=(T-T00)/(T0-T00) cf sc doc and TAO etal (89)
!
LOGICAL             :: LPRETREATMENT, LNEW_ADJUST
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
!
ZEPS= XMV / XMD
!
LPRETREATMENT=.TRUE.     ! FALSE to retreive the previous MASDEV4_1 version
LNEW_ADJUST  =.TRUE.     ! FALSE to retreive the previous MASDEV4_1 version
ZT0  = XTT               ! Usefull if LPRETREATMENT=T or LNEW_ADJUST=T
ZT00 = XTT-40.           ! Usefull if LPRETREATMENT=T or LNEW_ADJUST=T
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!
!
!*       4.     SECOND ORDER ALL OR NOTHING CONDENSATION SCHEME
!                            FOR MIXED-PHASE CLOUD
!               -----------------------------------------------
!
!
!*       4.1    Eventually pretreatment
!
IF (LPRETREATMENT) THEN
  !
  !*    compute the saturation vapor pressures at t+1
  !
  ZW1(:,:,:) = EXP( XALPW - XBETAW/PT(:,:,:) - XGAMW*ALOG(PT(:,:,:)) ) ! e_sw
  ZW2(:,:,:) = EXP( XALPI - XBETAI/PT(:,:,:) - XGAMI*ALOG(PT(:,:,:)) ) ! e_si
  ZW1(:,:,:) =  MIN(PPABS(:,:,:)/2.,ZW1(:,:,:))   ! safety limitation
  ZW2(:,:,:) =  MIN(PPABS(:,:,:)/2.,ZW2(:,:,:))   ! safety limitation
  !
  !*    compute the saturation mixing ratios at t+1
  !
  ZW3(:,:,:) =  ZW1(:,:,:) * ZEPS     &
             / ( PPABS(:,:,:) - ZW1(:,:,:) )     ! r_sw
  ZW4(:,:,:) =  ZW2(:,:,:) * ZEPS     &
             / ( PPABS(:,:,:) - ZW2(:,:,:) )     ! r_si
  !
  WHERE(PRV(:,:,:).LT.ZW4(:,:,:).AND. PRC(:,:,:).GT.0..AND. PT(:,:,:).LT.XTT)
    !
    !       Subsaturation case with respect to rsi(T,P) (and case rv<0):
    !       Evaporation of rc>0 (while enough) to decrease the lack of vapor
    !
    ZW5 (:,:,:)= MIN( PRC , ZW4(:,:,:) - PRV(:,:,:) )                ! RVCNDC
    PRV(:,:,:)= PRV(:,:,:) + ZW5(:,:,:)
    PRC(:,:,:)= PRC(:,:,:) - ZW5(:,:,:)
    PT(:,:,:)= PT(:,:,:) - ZW5(:,:,:) * PLV(:,:,:) /PCPH(:,:,:)
    !
  END WHERE
  !
  WHERE ( PRV(:,:,:) .GT. ZW3(:,:,:) )
    !
    !       Supersaturation case with respect to rsw(T,P):
    !       Condensation of the vapor that is left
    !
    ZW5 (:,:,:)= PRV(:,:,:) - ZW3(:,:,:)
    PRV(:,:,:)= PRV(:,:,:) - ZW5(:,:,:)                                  ! RVCNDC
    PRC(:,:,:)= PRC(:,:,:) + ZW5(:,:,:)
    PT(:,:,:)= PT(:,:,:) + ZW5(:,:,:) * PLV(:,:,:) /PCPH(:,:,:)
    !
  END WHERE
  !
  WHERE ( PRC(:,:,:).GT.0. .AND. PT(:,:,:).LT.ZT00 )
    !
    !       Treatment of rc>0 if T<T00:
    !
    PRI(:,:,:)= PRI(:,:,:) + PRC(:,:,:)
    PT(:,:,:)= PT(:,:,:) +                                                     &
                 PRC(:,:,:) * (PLS(:,:,:)-PLV(:,:,:)) /PCPH(:,:,:)
    PRC(:,:,:)= 0.
    !
  END WHERE
  !
END IF  !end PRETREATMENT
!
!*       4.3    compute the saturation vapor pressures at t+1
!
ZW1(:,:,:) = EXP( XALPW - XBETAW/PT(:,:,:) - XGAMW*ALOG(PT(:,:,:)) ) ! e_sw
ZW2(:,:,:) = EXP( XALPI - XBETAI/PT(:,:,:) - XGAMI*ALOG(PT(:,:,:)) ) ! e_si
ZW1(:,:,:) =  MIN(PPABS(:,:,:)/2.,ZW1(:,:,:))   ! safety limitation
ZW2(:,:,:) =  MIN(PPABS(:,:,:)/2.,ZW2(:,:,:))   ! safety limitation
!
!*       4.4    compute the saturation mixing ratios at t+1
!
ZW3(:,:,:) =  ZW1(:,:,:) * ZEPS     &
           / ( PPABS(:,:,:) - ZW1(:,:,:) ) ! r_sw
ZW4(:,:,:) =  ZW2(:,:,:) * ZEPS     &
           / ( PPABS(:,:,:) - ZW2(:,:,:) ) ! r_si
!
!*       4.5    compute the saturation mixing ratio derivatives (r'_vs)
!
ZW1(:,:,:) = (( XBETAW/PT(:,:,:) - XGAMW ) / PT(:,:,:)) & ! r'_sw
                * ZW3(:,:,:) * ( 1. + ZW3(:,:,:)/ZEPS )
ZW2(:,:,:) = (( XBETAI/PT(:,:,:) - XGAMI ) / PT(:,:,:)) & ! r'_si
                * ZW4(:,:,:) * ( 1. + ZW4(:,:,:)/ZEPS )
!
IF (LNEW_ADJUST) THEN
  ZCND(:,:,:)= ( PT(:,:,:) - ZT00 ) / (ZT0-ZT00)       ! Like Tao et al 89
  ZCND(:,:,:)= MAX ( MIN(ZCND(:,:,:),1.) , 0. )
ELSE
  ZCND(:,:,:) = 1.
  WHERE( (PRC(:,:,:)+PRI(:,:,:)) .GT. 1.0E-20 )      &
  ZCND(:,:,:)= PRC(:,:,:) / (PRC(:,:,:)+PRI(:,:,:)) ! Like the original version
END IF
!
!*       4.5    compute L_v CND + L_s DEP and F'(T)
!
WHERE( (PRC(:,:,:)+PRI(:,:,:)) .GT. 1.0E-20 )
  !
  ZW5(:,:,:) = PLS(:,:,:) + ( PLV(:,:,:)-PLS(:,:,:) ) * ZCND(:,:,:)
  ZW6(:,:,:) = PCPH(:,:,:) * ( PRC(:,:,:) + PRI(:,:,:) ) +               &
               ZW5(:,:,:)  * ( PRC(:,:,:)*ZW1(:,:,:)                      &
                                           + PRI(:,:,:)*ZW2(:,:,:) )
  !
  !*       4.6    compute Delta 2
  !
  ZW7(:,:,:) = (ZW5(:,:,:)/(ZW6(:,:,:)*PT(:,:,:))) *                       &
               ( PRC(:,:,:)*ZW1(:,:,:) *                                  &
        ((-2.*XBETAW+XGAMW*PT(:,:,:)) / (XBETAW-XGAMW*PT(:,:,:))           &
           + (XBETAW-XGAMW*PT(:,:,:))*(1.0+2.0*ZW3(:,:,:)/ZEPS)/PT(:,:,:)) &
               + PRI(:,:,:)*ZW2(:,:,:) *                                  &
        ((-2.*XBETAI+XGAMI*PT(:,:,:)) / (XBETAI-XGAMI*PT(:,:,:))           &
           + (XBETAI-XGAMI*PT(:,:,:))*(1.0+2.0*ZW4(:,:,:)/ZEPS)/PT(:,:,:)) )
  !
  !*       4.7    compute Delta 1
  !
  ZW6(:,:,:) = ZW5(:,:,:)*( PRC(:,:,:)*ZW3(:,:,:) + PRI(:,:,:)*ZW4(:,:,:)&
                          - PRV(:,:,:)*(PRC(:,:,:) + PRI(:,:,:))    &
                          )/ZW6(:,:,:)
  !
  !*       4.8    compute the sources
  !
  ZW3(:,:,:) = (PCPH(:,:,:)/ZW5(:,:,:)) *                               &
               (-ZW6(:,:,:) * ( 1.0 + 0.5*ZW6(:,:,:)*ZW7(:,:,:) ))
  ! With ZCND computed as in the original version, the following limitations were not needed
  ! because ZCND is in agrrement with the actual PRC and PRI (if PRC is 0, ZCND is also 0
  ! for instance).
  ! With ZCND computed following Tao et al 89, we can have, for instance, PRC=0 with ZCND!=0
  ! and ZW3 < 0, in this case we need to limit the tendencies.
  ! The way these limitations are implemented here is not correct because the resulting PRC
  ! and PRI are not correctly adjusted when limitation accurs.
  ZW3(:,:,:) = MIN(ZW3(:,:,:),  PRV(:,:,:)) ! Do not know if we need this limitation
  ZW6(:,:,:) = -MIN(-ZW3(:,:,:)*ZCND(:,:,:), PRC(:,:,:))
  ZW7(:,:,:) = -MIN(-ZW3(:,:,:)*( 1.0 - ZCND(:,:,:) ), PRI(:,:,:))

  PRV(:,:,:) = PRV(:,:,:) - ZW3(:,:,:)
  PRC(:,:,:) = PRC(:,:,:) + ZW6(:,:,:) !ZW3(:,:,:)*ZCND(:,:,:)                               ! RVCNDC
  PRI(:,:,:) = PRI(:,:,:) + ZW7(:,:,:) !ZW3(:,:,:)*( 1.0 - ZCND(:,:,:) )                     ! RVDEPI
  !PT(:,:,:)= PT(:,:,:) + (ZCND(:,:,:) * PLV(:,:,:) + ( 1.0 - ZCND(:,:,:) ) * PLS(:,:,:)) * &
  !                           ZW3(:,:,:) /PCPH(:,:,:)
  PT(:,:,:)= PT(:,:,:) + (ZW6(:,:,:) * PLV(:,:,:) + ZW7(:,:,:) * PLS(:,:,:)) / PCPH(:,:,:)
  !
ELSEWHERE
  !
  !*       4.9    special case when both r_c and r_i are zero
  !
  ZW6(:,:,:) = PCPH(:,:,:) + PLV(:,:,:) * ZW1(:,:,:)               ! F'(T)
  ZW7(:,:,:) = (PLV(:,:,:)/(ZW6(:,:,:)*PT(:,:,:))) * &             ! Delta 2
               ( ZW1(:,:,:) *                                              &
      ( (-2.*XBETAW+XGAMW*PT(:,:,:)) / (XBETAW-XGAMW*PT(:,:,:))            &
          + (XBETAW-XGAMW*PT(:,:,:))*(1.0+2.0*ZW3(:,:,:)/ZEPS)/PT(:,:,:) ) )
  ZW6(:,:,:) = PLV(:,:,:)*( ZW3(:,:,:)-PRV(:,:,:) )/ZW6(:,:,:)! Delta 1
  ZW3(:,:,:) = (PCPH(:,:,:)/PLV(:,:,:)) *                          &! RVCNDC
               (-ZW6(:,:,:) * ( 1.0 + 0.5*ZW6(:,:,:)*ZW7(:,:,:) ))
  ZW3(:,:,:) = MAX(ZW3(:,:,:), 0.) ! This line was not here before sept 2016
  PRV(:,:,:) = PRV(:,:,:) - ZW3(:,:,:)
  PRC(:,:,:) = PRC(:,:,:) + ZW3(:,:,:)
  PT(:,:,:)= PT(:,:,:) + ZW3(:,:,:) * PLV(:,:,:) /PCPH(:,:,:)
  !
END WHERE
print*, "Langlois must be modified: limitations are not correct in first where case and ZCND is not used in second where case"
!
!
!------------------------------------------------------------------------------
!
!
END SUBROUTINE ADJUST_LANGLOIS
