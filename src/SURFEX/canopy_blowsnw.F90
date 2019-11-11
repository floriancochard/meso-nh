!     #########################################
      SUBROUTINE CANOPY_BLOWSNW(SB,KI,KSNW,PTSTEP,BLOWSNW,  &
      				PZZ,PRHOA,PBLOWSNWA,PSFLUX_BLOWSNW,				&
      				PSFTH,PSFTQ)
!     #########################################
!
!!****  *CANOPY_BLOWSNW* - driver for evolution of blowing snow variables in Canopy
!!
!!
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	  V. Vionnet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2014
!!      Modif  
!!          J.Escobar 27/04/2018 : BUG?! => uncomment USE MODI_CANOPY_BLOWSNW_SUBL
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CANOPY_n, ONLY : CANOPY_t
!
USE MODD_BLOWSNW_n,  ONLY : BLOWSNW_t
USE MODD_BLOWSNW_SURF
USE MODD_CANOPY_TURB, ONLY : XCMFS, XCSHF

USE MODE_BLOWSNW_SURF

USE MODI_RMC01_SURF
USE MODI_CANOPY_EVOL_BLOWSNW
USE MODI_CANOPY_BLOWSNW_SUBL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(CANOPY_t), INTENT(INOUT) :: SB
TYPE(BLOWSNW_t), INTENT(INOUT)                :: BLOWSNW
!
INTEGER,                  INTENT(IN)    :: KI        ! number of horizontal points
INTEGER,                  INTENT(IN)    :: KSNW      ! number of blowing snow variables
!
REAL,                     INTENT(IN)    :: PTSTEP    ! atmospheric time-step
!
REAL, DIMENSION(KI,SB%NLVL), INTENT(IN)    :: PZZ       ! Mixing length generic profile at mid levels (-)
REAL, DIMENSION(KI),      INTENT(IN)    :: PRHOA     ! Air density at forcing level          (kg/m3)
REAL, DIMENSION(KI,KSNW), INTENT(IN)    :: PBLOWSNWA ! Blowing snow variables at forcing level
						! 1: #/kg  2: kg/kg
!
REAL, DIMENSION(KI,KSNW), INTENT(INOUT)    :: PSFLUX_BLOWSNW! surface flux of blowing snow
!                                 IN : surface turbulent flux 1:(# m-2 s-1) 2:(kg(snow) m2 s-1)
!                                 OUT : net flux (turb+sedim) at Top of Canopy 1: (#/m2/s) 2: (kg(snow)/m2/s)
REAL, DIMENSION(KI), INTENT(INOUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(INOUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JLAYER,JSV                              ! loop counter

REAL, DIMENSION(KI,SB%NLVL) :: ZK                 ! mixing coefficient
REAL, DIMENSION(KI,SB%NLVL)   :: ZLM       ! mixing length                         (m)
REAL, DIMENSION(KI,SB%NLVL)   :: ZLEPS     ! dissipative length      (m)

REAL, DIMENSION(KI,SB%NLVL,KSNW)  :: ZBLOWSNWS   ! snow variables tendency due to sublimation at Canopy levels
!                                                           1: #_s/kg/s, 2: kg_s/kg/s

REAL, DIMENSION(KI,KSNW)       :: ZSFLUX_BLOWSNW! surface flux of blowing snow variables
!                                      1:(#/kg(air) m s-1) 2:(kg(snow)/kg(air) m s-1)


REAL, DIMENSION(KI,KSNW)  :: ZBLOWSNW_DEP   ! sedimentation flux of blowing snow variables at Canopy bottom
!                                                           1: (#/m2/s) 2: (kg(snow)/m2/s
REAL, DIMENSION(KI,SB%NLVL,KSNW) :: ZSNW      ! Blowing snow variables at Canopy levels
					    ! 1: #/m3  2: kg/m3
REAL, DIMENSION(KI,SB%NLVL)      :: ZRG    ! Mean radius  (m)                                            


REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CANOPY_BLOWSNW',0,ZHOOK_HANDLE)

!
! convert flux from kg/m2/s to [kg(snow)/kg(dry air).m.s-1]
DO JSV= 1,KSNW
     ZSFLUX_BLOWSNW(:,JSV) = PSFLUX_BLOWSNW(:,JSV)/ PRHOA(:)
END DO

ZBLOWSNW_DEP = 0.
ZBLOWSNWS    = 0.

!-------------------------------------------------------------------------------
!
!*    1. mixing and dissipative lengths (at full levels)
!        ------------------------------
!
CALL RMC01_SURF(PZZ,SB%XLMO,ZLM,ZLEPS)
!
!*    2. mixing coefficient for scalars at mid layers (at half levels)
!         ---------------------------
!
ZK = -999.
DO JLAYER=2,SB%NLVL
   ZK(:,JLAYER) = 0.5 * XCSHF * ZLM(:,JLAYER)   * SQRT(SB%XTKE(:,JLAYER)  ) &
                + 0.5 * XCSHF * ZLM(:,JLAYER-1) * SQRT(SB%XTKE(:,JLAYER-1))
END DO

!
!*    3. Compute the impact of blowing snow sublimation  on air temperature,
! humidity and blowing snow number and mass.
!         ---------------------------
!
IF(LBLOWSNW_CANOSUBL) THEN

    CALL  CANOPY_BLOWSNW_SUBL(KI,SB%NLVL,KSNW,PTSTEP,SB%XBLOWSNW,SB%XT,SB%XQ,PRHOA,SB%XP,  &
                      ZBLOWSNWS,SB%XDZ,PSFTQ,PSFTH,BLOWSNW%XSNW_SUBL)

ENDIF

!
!*    4. effects of turbulent diffusion and sedimentation on blowing snow variables in Canopy
!         ---------------------------
!

 CALL CANOPY_EVOL_BLOWSNW(KI,SB%NLVL,KSNW,PTSTEP,PBLOWSNWA,ZK,ZSFLUX_BLOWSNW   ,&
                     SB%XZ,SB%XDZ,SB%XDZF,SB%XBLOWSNW,PRHOA,ZBLOWSNW_DEP,SB%XT, &
                     SB%XP,ZBLOWSNWS,SB%XTKE)

 !
 !*   5. Security at atmospheric forcing level
 !        -------------------------------------
 !
DO JSV=1,KSNW
  	SB%XBLOWSNW(:,SB%NLVL,JSV) = PBLOWSNWA(:,JSV)
END DO

!
 !*   6. Store and update fluxes and diagnostics
 !        -------------------------------------
 !
!
!Convert canopy variables in _/m3
DO JLAYER=1,SB%NLVL
  DO JSV=1,KSNW
     ZSNW(:,JLAYER,JSV)=SB%XBLOWSNW(:,JLAYER,JSV)*PRHOA(:)
  ENDDO
ENDDO
!
! Compute mean radius
!
CALL SNOWMOMENT2SIZE(ZSNW, PRG1D=ZRG )
! 
! Store diagnostic variables
!
BLOWSNW%XSNW_CANO_RGA(:,:)   = ZRG(:,:)
BLOWSNW%XSNW_CANO_VAR(:,:,:) = ZSNW(:,:,:)
!
! Store turbulent flux at the bottom of Canopy as a diagnostic
BLOWSNW%XSNW_FTURB(:,1:2) =  -PSFLUX_BLOWSNW(:,1:2)   ! Instantaneous fluxes (number and mass)
BLOWSNW%XSNW_FTURB(:,3) =  BLOWSNW%XSNW_FTURB(:,3)+BLOWSNW%XSNW_FTURB(:,2)*PTSTEP ! Accumulated flux (mass)

! Store total net flux (saltation+suspension) at the top of the snowpack
BLOWSNW%XSNW_FNET(:,1:2) =  BLOWSNW%XSNW_FTURB(:,1:2)+BLOWSNW%XSNW_FSED(:,1:2)  ! Instantaneous flux (1: number; 2: mass)
BLOWSNW%XSNW_FNET(:,3) =  BLOWSNW%XSNW_FNET(:,3)+BLOWSNW%XSNW_FNET(:,2)*PTSTEP ! Accumulated flux (mass)

! Store sedimentation flux at the bottom of Canopy to be sent to Crocus
BLOWSNW%XSFSNW(:,1:2) = ZBLOWSNW_DEP(:,1:2)   ! Instantaneous fluxes (number and mass)

! Update flux : net flux (turb+sedim) at Top of Canopy 1: (#/m2/s) 2: (kg(snow)/m2/s)
DO JSV=1,KSNW
  PSFLUX_BLOWSNW(:,JSV) = ZSFLUX_BLOWSNW(:,JSV)
END DO


IF (LHOOK) CALL DR_HOOK('CANOPY_BLOWSNW',1,ZHOOK_HANDLE)

!
!-------------------------------------------------------------------------------
END SUBROUTINE CANOPY_BLOWSNW
