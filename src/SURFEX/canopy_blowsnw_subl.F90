!     #########################################
      SUBROUTINE CANOPY_BLOWSNW_SUBL(KI,KLVL,KSNW,PTSTEP,PBLOWSNW,PT,PQ,PRHOA,PP,         &
                         PBLOWSNWS,PDZ,PSFTQ,PSFTH,PSNW_SUBL )

USE MODD_CSTS
USE MODD_BLOWSNW_SURF
!
USE MODE_BLOWSNW_SURF
!
USE MODI_BLOWSNW_VELGRAV1D

IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!      
INTEGER,                  INTENT(IN)    :: KI        ! number of horizontal points
INTEGER,                  INTENT(IN)    :: KLVL      ! number of levels in canopy
INTEGER,                  INTENT(IN)    :: KSNW      ! number of snow variables in canopy
REAL,                     INTENT(IN)    :: PTSTEP    ! time-step    
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PT        ! Temperature at canopy levels          (K) 
REAL, DIMENSION(KI,KLVL),INTENT(IN)     :: PP        ! air presseure at canopy levels (Pa)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PQ        ! humidity at canopy levels          (kg/m3)
REAL, DIMENSION(KI,KLVL,KSNW), INTENT(IN)  :: PBLOWSNW  ! blowing snow variables at canopy levels   (1: #_s/kg, 2: kg_s/kg)
REAL, DIMENSION(KI,KLVL), INTENT(IN)    :: PDZ       ! depth   of canopy levels              (m)

REAL, DIMENSION(KI), INTENT(IN)                 :: PRHOA      !  air density (kg/m3)

REAL, DIMENSION(KI), INTENT(INOUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(INOUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)

REAL, DIMENSION(KI,KLVL,KSNW), INTENT(INOUT)  :: PBLOWSNWS   ! snow variables tendency due to sublimation at Canopy levels
!                                                           (1: #_s/kg/s, 2: kg_s/kg/s)
REAL, DIMENSION(KI,KSNW+1), INTENT(INOUT)  :: PSNW_SUBL   ! Diagnostic !Sublimation rate
                                                    ! 1: Instantaneous number 2: Instantaneous mass 3: Accumulated Mass 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER JK,JI,JSV

REAL, DIMENSION(KI,KLVL)      :: ZBET   ! Scale parameter of the gamma distribution (m)
REAL, DIMENSION(KI,KLVL,KSNW) :: ZVGK   ! Terminal fallspeed (m/s)
REAL, DIMENSION(KI,KLVL)      :: ZRG    ! Mean radius  (m)
REAL, DIMENSION(KI,KLVL)      :: ZEXN   ! Exner function        at full levels
REAL, DIMENSION(KI,KLVL)      :: ZW     ! working array
REAL, DIMENSION(KI)           :: ZT     ! average temprature of Canopy
REAL, DIMENSION(KI,KLVL)      :: ZLSFACT! L_s/(C_ph)

REAL, DIMENSION(KI,KLVL,KSNW) :: ZSNW      ! work variable : pot. temp at futur instant 
!                                       ! (or past at the end of the routine)
REAL, DIMENSION(KI,KLVL-1)      :: ZMU    ! air kinematic viscosity
REAL, DIMENSION(KI,KLVL-1)      :: ZNU    ! Nusselt number  (-)
REAL, DIMENSION(KI,KLVL-1)      :: ZKA    ! air thermal conductivity 
REAL, DIMENSION(KI,KLVL-1)      :: ZDV    ! diffusivity of water vapor in the air. 
REAL, DIMENSION(KI,KLVL-1)      :: ZAI    ! denominator in Thorpe and Masson (66) equation
REAL, DIMENSION(KI,KLVL-1)      :: ZUSI   ! undersaturation over ice 
REAL, DIMENSION(KI,KLVL-1)      :: ZSUBL_NUM,ZSUBL_MASS  ! sublimation rate (number and mass) 
REAL, DIMENSION(KI)             :: ZCOL_SUBL ! integrated sublimation rate (kg/m2/s)
!
LOGICAL                         :: LSUBL_ALPINE3D
REAL                            :: ZR_EFF   ! Effective radius for the computation of sublimation
REAL, DIMENSION(KI,KLVL)        :: ZAA   ! Coeff used in the computation of terminal fallspeed (m/s)
REAL, DIMENSION(KI,KLVL)        :: ZBB   ! Coeff used in the computation of terminal fallspeed (m/s)
!
!-------------------------------------------------------------------------------
!
!
!*    1. initializations
! 
!-------------------------------------------------------------------------------
DO JK=1,(KLVL-1)
! Sutherland's equation for kinematic viscosity
   ZMU(:,JK)=1.8325d-5*416.16/(PT(:,JK)+120)*(PT(:,JK)/296.16)*SQRT(PT(:,JK)/296.16)/PRHOA(:)
! Thermal conductivity of the air
   ZKA(:,JK) = 2.38E-2 + 0.0071E-2 * ( PT(:,JK) - XTT )          ! k_a
! Diffusivity of water vapor in the air. 
   ZDV(:,JK) = 0.211E-4 * (PT(:,JK)/XTT)**1.94 * (XP00/PP(:,JK)) ! D_v
END DO

!*       Compute the denominator in the Thorpe and Masson (66) equation
! 
ZAI(:,:) = EXP( XALPI - XBETAI/PT(:,1:(KLVL-1)) - XGAMI*ALOG(PT(:,1:(KLVL-1))))  ! es_i
!
DO JK=1,(KLVL-1)
  ! Undersaturation over ice
   ZUSI(:,JK) = PQ(:,JK)/PRHOA(:)*( PP(:,JK)-ZAI(:,JK) ) / ( (XMV/XMD) * ZAI(:,JK) ) - 1.0
  ! denominator in the Thorpe and Masson (66) equation
   ZAI(:,JK) = ( XLSTT + (XCPV-XCI)*(PT(:,JK)-XTT) ) / (ZKA(:,JK)*PT(:,JK))                  &
             *( ( XLSTT + (XCPV-XCI)*(PT(:,JK)-XTT) ) / (XRV*PT(:,JK)) - 1.)  &
             + (XRV*PT(:,JK)) / (ZDV(:,JK)*ZAI(:,JK))
END DO

ZUSI(:,:)  =MIN(ZUSI(:,:),0.)

! 
DO JK=1,(KLVL)
     ZW(:,JK)  = XCPD+XCPV*PQ(:,JK)/PRHOA(:)
END DO
!
ZLSFACT(:,:) = (XLSTT+(XCPV-XCI)*(PT(:,:)-XTT))/ZW(:,:) ! L_s/(*C_ph)
!
ZEXN = (PP/XP00)**(XRD/XCPD)
!
ZSUBL_NUM = 0.
ZSUBL_MASS = 0.
ZCOL_SUBL=0.
!
LSUBL_ALPINE3D    = .FALSE.   ! Compute sublimation using the method of reprsentative
                              ! radius implemented in Alpine 3D (Groot and al, 2011)
!

IF(LSUBL_ALPINE3D) THEN
!
!-------------------------------------------------------------------------------
!
!
!*    2. Compute settling velocity based on the effective radius
!        ---------------
!
!-------------------------------------------------------------------------------
         ZR_EFF = 73.5e-6  ! Effective radius computed following the Swiss
                          ! method. This effective radius give the same total 
                          ! sublimation for a equal concentration an ensemble of
                          ! gamma distributed particles with rm = 35e-6 m and
                          ! alpha=3
          ZRG(:,:)  = ZR_EFF
          ZBET(:,:) = ZR_EFF/XEMIALPHA_SNW
!
!
!*         compute gravitational velocities
!
          CALL BLOWSNW_VELGRAV1D(ZBET, ZRG, PT, PRHOA,PP, ZVGK)
!
!-------------------------------------------------------------------------------
!
!
!*    3. Compute sublimation rate based on formulation of Thorpe and Masson (66)
!        Note that no integration is computed on the particle spectra
!        ---------------
!
!-------------------------------------------------------------------------------
!
!       Nusselt Number using effective radius
ZNU(:,:)    =    NUSSELT(ZRG(:,1:(KLVL-1)),ZMU(:,1:(KLVL-1)),ZVGK(:,1:(KLVL-1),2))

ZSUBL_MASS(:,:) = 3*PBLOWSNW(:,1:(KLVL-1),2)*ZNU(:,:)*ZUSI(:,:)  &
                       /(2*ZAI(:,:)*XRHOLI*ZR_EFF**2)


ELSE


!-------------------------------------------------------------------------------
!
!
!*    2. Compute settling velocity based on the size distribution of the previous time
!        step: used as ventilation velocity
!        ---------------
!
!-------------------------------------------------------------------------------

!
!           Convert canopy variables in _/m3
!
DO JK=1,KLVL
  DO JSV=1,KSNW
     ZSNW(:,JK,JSV)=PBLOWSNW(:,JK,JSV)*PRHOA(:)
  ENDDO
ENDDO
!
!*         compute BETA, RG and moments
!
CALL SNOWMOMENT2SIZE(ZSNW, PBETA1D=ZBET, PRG1D=ZRG )
!
!*         compute gravitational velocities
!
CALL BLOWSNW_VELGRAV1D(ZBET, ZRG, PT, PRHOA,PP, ZVGK)
!
!-------------------------------------------------------------------------------
!
!
!*    3. Compute sublimation rate based on formulation of Thorpe and Masson (66)
!        ---------------
!
!-------------------------------------------------------------------------------
!
!       Nusselt Number using mean radius of particle size distribution
ZNU(:,:)    =    NUSSELT(ZRG(:,1:(KLVL-1)),ZMU(:,1:(KLVL-1)),ZVGK(:,1:(KLVL-1),2))


! mass averaged sublimation rate follows Dery and Yau (1999) and avoids
! numerical integration over the particle spectrum
ZSUBL_MASS(:,:) = PBLOWSNW(:,1:(KLVL-1),2)*ZNU(:,:)*ZUSI(:,:)/   &
                        (ZAI(:,:)*2*XRHOLI*ZRG(:,1:(KLVL-1))**2)
!
!
ENDIF                        
!   Restriction of ZSUM_SUBL to insure coherence between vapor and blowing snow
!   mixing ratio
! 
DO JI=1,KI
    DO JK=1,(KLVL-1)
          ZSUBL_MASS(JI,JK) =  MIN( PQ(JI,JK)/(PRHOA(JI)*PTSTEP),ZSUBL_MASS(JI,JK))*               &
                               (0.5+SIGN(0.5,ZSUBL_MASS(JI,JK)))-                                  &
                               MIN(PBLOWSNW(JI,JK,2)/PTSTEP,ABS(ZSUBL_MASS(JI,JK)))*                   &
                               (0.5-SIGN(0.5,ZSUBL_MASS(JI,JK)))
          ZSUBL_MASS(JI,JK)=MIN(0.,ZSUBL_MASS(JI,JK)) ! Sink of snow
!
!         number-averaged sublimation rate    
!  Change in concentration rate  Sn = Sb*N/qb (Dery and Yau,2000)
!
          IF(ZSUBL_MASS(JI,JK)<0.) THEN
                ZSUBL_NUM(JI,JK) = ZSUBL_MASS(JI,JK)*PBLOWSNW(JI,JK,1)/PBLOWSNW(JI,JK,2)
          END IF
          ZCOL_SUBL(JI) =  ZCOL_SUBL(JI)+ZSUBL_MASS(JI,JK)*PRHOA(JI)*PDZ(JI,JK)
!*       Compute modification of heat and vapor fluxes due to
!        sublimation of blowing snow particles in Canopy
!        These fluxes are then sent to MNH
         PSFTQ(JI) = PSFTQ(JI) - ZSUBL_MASS(JI,JK)*PRHOA(JI)*PDZ(JI,JK)
         PSFTH(JI) = PSFTH(JI) + ZSUBL_MASS(JI,JK)*PRHOA(JI)*PDZ(JI,JK)* &
                                (XLSTT+(XCPV-XCI)*(PT(JI,JK)-XTT))/ZEXN(JI,JK)
    END DO
!    ZT(JI) = SUM(PT(JI,1:(KLVL-1)))/(KLVL-1)
END DO   
!
!       Store number and mass sublimation rate
PBLOWSNWS(:,1:(KLVL-1),1) =  ZSUBL_NUM(:,:)
PBLOWSNWS(:,1:(KLVL-1),2) =  ZSUBL_MASS(:,:)

!  Store integrated sublimation rate in mmSWE/day
PSNW_SUBL(:,2) = ZCOL_SUBL(:)*3600*24/XRHOLW*1000
!  Store accumulated sublimation rate in mmSWE
PSNW_SUBL(:,3) = PSNW_SUBL(:,3)+ZCOL_SUBL(:)*PTSTEP

CONTAINS

FUNCTION NUSSELT(PR,PMU,PVEL_VENT) RESULT(PNU)
!
!!    PURPOSE
!!    -------
!     Calculate the Nusselt number for a given particle radius 
!     Formulation based on Lee (1975)
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)              :: PR ! (m)
REAL, DIMENSION(:,:), INTENT(IN)              :: PMU ! (m2/s)
REAL, DIMENSION(:,:), INTENT(IN)              :: PVEL_VENT ! (m/s)
!
REAL, DIMENSION(SIZE(PMU,1),SIZE(PMU,2))      :: PNU ! (m/s)
!
!
!*      0.2    declaration of local variables
!
REAL,  DIMENSION(SIZE(PMU,1),SIZE(PMU,2))    :: ZRE
!
!
!*      1    Calculate Reynolds number 
!
ZRE(:,:) = 2*PR(:,:)*PVEL_VENT(:,:)/PMU(:,:)
!
!*      2    Calculate Nusselt number
!
WHERE(ZRE(:,:)<10)
        PNU(:,:) = 1.79+0.606*ZRE(:,:)**0.5
ELSEWHERE
        PNU(:,:) = 1.88+0.580*ZRE(:,:)**0.5
END WHERE

END FUNCTION NUSSELT

END SUBROUTINE CANOPY_BLOWSNW_SUBL
