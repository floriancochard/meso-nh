!     ######spl
MODULE MODE_BLOWSNW_CANOPY
!     ####################
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
CONTAINS

SUBROUTINE INIT_BLOWSNW_SBL(KI,KLVL,KSV_SNWBEG,KSV_SNWEND,KSNWEQ,K2D_SNWBEG,  &
               K2D_SNWEND ,K2DSNWEQ, PTSTEP,BLOWSNW,PDZ,PU,PUREF,           &
               PRHOA,PSV,PBLOWSNWA,PBLOWSNW)
!
USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
!
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
USE MODD_BLOWSNW_n
USE MODD_BLOWSNW_SURF


IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(BLOWSNW_t), INTENT(INOUT)                :: BLOWSNW
INTEGER, INTENT(IN)				 :: KI        ! number of grid points
INTEGER, INTENT(IN)				 :: KLVL      ! number of levels in Canopy
INTEGER, INTENT(IN)    :: KSV_SNWBEG, KSV_SNWEND    ! index of first and last blowing snow related scalar variable
INTEGER, INTENT(IN)    :: KSNWEQ                    ! number of blowing snow related species in scalar variables list
INTEGER, INTENT(IN)    :: K2D_SNWBEG, K2D_SNWEND    ! index of first and last blowing snow 2D variable sent to MNH
INTEGER, INTENT(IN)    :: K2DSNWEQ                  ! number of blowing snow 2D related species in scalar variables list
REAL, INTENT(IN)  				 :: PTSTEP    ! atmospheric time-step                 (s)


REAL, DIMENSION(:,:), INTENT(IN) :: PDZ
REAL, DIMENSION(:,:), INTENT(IN) :: PU


REAL, DIMENSION(:), INTENT(IN)   :: PUREF
REAL, DIMENSION(:), INTENT(IN)   :: PRHOA

REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSV
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PBLOWSNW
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PBLOWSNWA   ! Blowing snow variables at forcing levels
										              !  1:#/kg(air)   2:kg/kg(air)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JLAYER,JI,JSV   ! loop counter

REAL, DIMENSION(KI)        :: ZZCANO ! Total thickness of Canopy layer
REAL, DIMENSION(KI,KSNWEQ) :: ZSNWC ! Total number and mass in Canopy (previous time step)
										     !  1:#/m2   2:kg/m2
REAL, DIMENSION(KI,KSNWEQ) :: ZSNWFC ! Transport rate in Canopy (over the whole layer) (previous time step)
										     !  1:#/m/s   2:kg/m/s
REAL, DIMENSION(KI,KSNWEQ) :: ZSNWC_NEW  ! Total number and mass in Canopy (after advection in MNH)
										     !  1:#/m2   2:kg/m2
REAL, DIMENSION(KI)        :: ZFSALT_NEW ! Updated streamise saltation flux
REAL, DIMENSION(KI)        :: ZFSALT_DIV ! Divergence of vector transport in saltation

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_CANOPY:INIT_BLOWSNW_SBL',0,ZHOOK_HANDLE)
!
!*       1.   Initialize blowing snow variables in Canopy
!              -------------------------------
!
! Initialize profiles of blowing snow variables in Canopy at 1st time step
IF(ANY(PBLOWSNW(:,:,:) == XUNDEF)) THEN
     PBLOWSNW(:,:,:)=0.
ENDIF

! Insure coherence between mass and number at Canopy levels
DO JLAYER=1,(KLVL-1)
	 WHERE( PBLOWSNW(:,JLAYER,1) == 0. .AND. PBLOWSNW(:,JLAYER,2) >0  )
	     PBLOWSNW(:,JLAYER,2)=0.
	 END WHERE

	 WHERE( PBLOWSNW(:,JLAYER,1) > 0. .AND. PBLOWSNW(:,JLAYER,2) ==0  )
	     PBLOWSNW(:,JLAYER,1)=0.
	 END WHERE
END DO

! Initialization of blowing snow variables at forcing level
DO JSV=KSV_SNWBEG,KSV_SNWEND
     PBLOWSNWA(:,JSV)= PSV(:,JSV)/PRHOA(:)     
END DO


! Put blowing snow variables at top of Canopy equals to blowing snow variables at forcing level
DO JSV=1,KSNWEQ
     PBLOWSNW(:,KLVL,JSV) = PBLOWSNWA(:,JSV)
ENDDO

! Insure coherence at top of Canopy
WHERE(PBLOWSNWA(:,1)==0. .AND. PBLOWSNWA(:,2)>0.)
     PBLOWSNW(:,KLVL,2) = 0.
END WHERE
WHERE(PBLOWSNWA(:,2)==0. .AND. PBLOWSNWA(:,1)>0.)
     PBLOWSNW(:,KLVL,1) = 0.
END WHERE

! Compute thickness of Canopy layer
DO JI=1,KI
     ZZCANO(JI) = SUM(PDZ(JI,1:(KLVL-1)))
ENDDO

!
!*       2.   Store variables necessary for compuation of surface fluxes
!              -------------------------------
!
! Store sedimentation flux at the bottom of Canopy as a diagnostic
BLOWSNW%XSNW_FSED(:,1:2) =  BLOWSNW%XSFSNW(:,1:2)   ! Instantaneous fluxes (number and mass)
BLOWSNW%XSNW_FSED(:,3) =  BLOWSNW%XSNW_FSED(:,3)+BLOWSNW%XSNW_FSED(:,2)*PTSTEP  ! Accumulated flux (mass)

!
!*       3.   Update variables if transport in saltation is computed
!              -------------------------------
!

IF(LBLOWSNW_ADV) THEN  ! Account for snow redistribution in saltation
!
! Updated value of streamwise saltation flux accounting for advection computed in MNH
     ZFSALT_NEW(:)= PSV(:,K2D_SNWEND)*2*PUREF(:)*PU(:,KLVL)
!
!
! Store contribution of saltation transport to snowpack mass balance:
!    Divergence of vector transport in saltation
!    (to be sent to snowpack_evol.f90)
!
     WHERE(PU(:,1)>0)
         ZFSALT_DIV(:) = (ZFSALT_NEW(:)-BLOWSNW%XSFSNW(:,3))/(PU(:,KLVL)*PTSTEP)
     END WHERE
!
! Store contribution of divergence of saltation flux to net surface flux
!
     WHERE(PU(:,1)>0)
        BLOWSNW%XSNW_FSALT(:,1) =  BLOWSNW%XSFSNW(:,3) ! Instantaneous streamwise saltation flux
        BLOWSNW%XSNW_FSALT(:,2) = (ZFSALT_NEW(:)-BLOWSNW%XSFSNW(:,3))/  &
                     (PU(:,KLVL)*PTSTEP)                      ! Instantaneous flux
        BLOWSNW%XSNW_FSALT(:,3) = BLOWSNW%XSNW_FSALT(:,3) + BLOWSNW%XSNW_FSALT(:,2)*PTSTEP     ! Accumulated flux
     END WHERE

! Store contribution of saltation transport to snowpack mass balance
     BLOWSNW%XSFSNW(:,3)=ZFSALT_DIV(:)

ELSE
     BLOWSNW%XSFSNW(:,3) = 0.
ENDIF 

!
!*       4.   Update blowing snow variables in Canopy
!				to account for advection effects computed in MNH
!              -------------------------------
!

ZSNWC =0.
ZSNWFC=0.
ZSNWC_NEW=0.

! Compute number and mass content and fluxes of Canopy layer.
DO JSV=1,KSNWEQ
     DO JLAYER=1,(KLVL-1)
        ZSNWC(:,JSV)  = ZSNWC(:,JSV)+PBLOWSNW(:,JLAYER,JSV)*PDZ(:,JLAYER)*PRHOA(:)
        ZSNWFC(:,JSV) = ZSNWFC(:,JSV)+PBLOWSNW(:,JLAYER,JSV)*PU(:,JLAYER)*PDZ(:,JLAYER)*PRHOA(:)
     END DO
END DO


! Insure coherence between number and mass profile
WHERE(ZSNWC(:,1)==0. .AND. ZSNWC(:,2)>0.)
     ZSNWC(:,2)=0.
END WHERE
WHERE(ZSNWC(:,2)==0. .AND. ZSNWC(:,1)>0.)
     ZSNWC(:,1)=0.
END WHERE

! Update total mass in Canopy to account for advection effect :
!      M'_Cano = M_Cano + M_Adv - M_Eq
! where : M'_Cano : total mass accounting for advection
!         M_Cano  : initial total mass
!         M_Adv   : equivalent mass after advection
!         M_Eq    : equivelent mass before advection
IF(LBLOWSNW_ADV) THEN
DO JSV=1,KSNWEQ
     ZSNWC_NEW(:,JSV)=ZSNWC(:,JSV)+(PSV(:,JSV+K2D_SNWBEG-1)*2*PUREF(:)*PU(:,KLVL)-ZSNWFC(:,JSV))/PU(:,KLVL)
END DO
ELSE
DO JSV=1,KSNWEQ
     ZSNWC_NEW(:,JSV)=ZSNWC(:,JSV)
END DO
ENDIF

! Apply advection contribution to canopy variables at each level: number, mass
DO JI=1,KI
  DO JSV=1,KSNWEQ
    DO JLAYER=1,(KLVL-1)
      IF(ZSNWC(JI,JSV)>0) THEN
!
! Blowing snow already in Canopy: update profile PBLOWSNW to account for advection
!
           PBLOWSNW(JI,JLAYER,JSV) = PBLOWSNW(JI,JLAYER,JSV)*ZSNWC_NEW(JI,JSV)/ZSNWC(JI,JSV)
           PBLOWSNW(JI,JLAYER,JSV) = MAX(PBLOWSNW(JI,JLAYER,JSV),0.)
      ELSE
           IF(ZSNWC_NEW(JI,JSV)>0)  THEN
!
! Snow transported by advection in a initially snow-free Canopy layer:
! Use a vertically uniform profile of number and mass for this Canopy layer
!
                PBLOWSNW(JI,JLAYER,JSV)=ZSNWC_NEW(JI,JSV)*ZZCANO(JI)/PRHOA(JI)
           ELSE
                PBLOWSNW(JI,JLAYER,JSV) = 0.
           END IF
      END IF
    END DO
  END DO
END DO

! Update blowing snow variables to be sent to snowpack_evol.f90 for computation of mass exchanges
! between the snow surface and the atmosphere.

DO JSV=KSV_SNWBEG,KSV_SNWEND
    PSV(:,JSV)=  PBLOWSNW(:,1,JSV)*PRHOA(:)
END DO

!
!*       5.   Store blowing snow sedimentation flux and net contribution of saltation transport
!               to be sent to snowpack_evol.f90
!              -------------------------------
!
PSV(:,K2D_SNWBEG:K2D_SNWEND) = BLOWSNW%XSFSNW(:,:)


IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_CANOPY:INIT_BLOWSNW_SBL',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_BLOWSNW_SBL


SUBROUTINE UPDATE_BLOWSNW_SBL(KI,KLVL,KSV_SNWBEG,KSV_SNWEND,KSNWEQ,K2D_SNWBEG,  &
               K2D_SNWEND ,K2DSNWEQ,BLOWSNW,PDZ,PU,PRHOA,PUREF,PBLOWSNW,PSFTS)

USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
!
USE MODD_SURF_PAR,      ONLY : XUNDEF
USE MODD_BLOWSNW_n
USE MODD_BLOWSNW_SURF

IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(BLOWSNW_t), INTENT(INOUT)                :: BLOWSNW
INTEGER, INTENT(IN)				 :: KI        ! number of grid points
INTEGER, INTENT(IN)				 :: KLVL      ! number of levels in Canopy
INTEGER, INTENT(IN)    :: KSV_SNWBEG, KSV_SNWEND    ! index of first and last blowing snow related scalar variable
INTEGER, INTENT(IN)    :: KSNWEQ                    ! number of blowing snow related species in scalar variables list
INTEGER, INTENT(IN)    :: K2D_SNWBEG, K2D_SNWEND    ! index of first and last blowing snow 2D variable sent to MNH
INTEGER, INTENT(IN)    :: K2DSNWEQ                  ! number of blowing snow 2D related species in scalar variables list

REAL, DIMENSION(:,:), INTENT(IN) :: PDZ
REAL, DIMENSION(:,:), INTENT(IN) :: PU

REAL, DIMENSION(:), INTENT(IN)   :: PUREF
REAL, DIMENSION(:), INTENT(IN)   :: PRHOA

REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PBLOWSNW

REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSFTS

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JLAYER,JSV   ! loop counter


REAL, DIMENSION(KI,KSNWEQ) :: ZSNWFC ! Transport rate in Canopy (over the whole layer) (previous time step)
										     !  1:#/m/s   2:kg/m/s

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_CANOPY:UPDATE_BLOWSNW_SBL',0,ZHOOK_HANDLE)
!
!*       1.   Update variables sent to MNH for computation of saltation divergence
!              -------------------------------

IF(LBLOWSNW_ADV) THEN ! Account for snow redistribution in saltation
!
! Store streamwise saltation flux
!
	BLOWSNW%XSFSNW(:,3) = PSFTS(:,K2D_SNWEND)
!
! Adapted streamwise saltation flux sent to MNH to compute advection
!
	PSFTS(:,K2D_SNWEND) =   PSFTS(:,K2D_SNWEND)/(2*PUREF(:)*PU(:,KLVL))
END IF

!
!*       2.   Update variables sent to MNH for computation of saltation divergence
!              -------------------------------

! Insure coherence between mass and number
DO JLAYER=1,(KLVL-1)

	WHERE( PBLOWSNW(:,JLAYER,1) == 0. .AND. PBLOWSNW(:,JLAYER,2) >0  )
	    PBLOWSNW(:,JLAYER,2)=0.
	END WHERE
	WHERE( PBLOWSNW(:,JLAYER,1) > 0. .AND. PBLOWSNW(:,JLAYER,2) ==0  )
	    PBLOWSNW(:,JLAYER,1)=0.
	END WHERE

END DO


ZSNWFC =0.
! Compute number and mass fluxes in Canopy layer.
DO JSV=1,KSNWEQ
    DO JLAYER=1,(KLVL-1)
        ZSNWFC(:,JSV) = ZSNWFC(:,JSV)+PBLOWSNW(:,JLAYER,JSV)*PU(:,JLAYER)*PDZ(:,JLAYER)*PRHOA(:)
    END DO
END DO


! Equivalent concentration sent to MNH to compute the contribution of advection
! to Canopy snow variables
DO JSV=1,KSNWEQ
     PSFTS(:,K2D_SNWBEG+JSV-1) = ZSNWFC(:,JSV)/(2*PUREF(:)*PU(:,KLVL))
END DO



!Insure coherence between number and mass
! Remove number where mass is not present
WHERE(PSFTS(:,K2D_SNWBEG+1)==0.)
	PSFTS(:,K2D_SNWBEG)=0.
END WHERE


IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_CANOPY:UPDATE_BLOWSNW_SBL',1,ZHOOK_HANDLE)

END SUBROUTINE UPDATE_BLOWSNW_SBL

END MODULE MODE_BLOWSNW_CANOPY
