!     #########################################
      SUBROUTINE CANOPY_EVOL_BLOWSNW(KI,KLVL,KSNW,PTSTEP,PSNWA,PK,         &
                                  PSFLUX_SNW,PZ,PDZ,PDZF,PSNW,PRHOA,PSNWDEP,          &
                                  PTHT,PPABST,PSNW_SUBL,PTKE)
!     #########################################
!
!!****  *CANOPY_EVOL_BLOWSNW* - evolution of blowing snow variables in canopy
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
!!	V. Vionnet   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BLOWSNW_n
USE MODD_BLOWSNW_SURF

USE MODE_BLOWSNW_SURF

USE MODI_TRIDIAG_GROUND
USE MODI_BLOWSNW_VELGRAV1D
USE MODI_BLOWSNW_DIFFK

!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                  INTENT(IN)    :: KI        ! number of horizontal points
INTEGER,                  INTENT(IN)    :: KLVL      ! number of levels in canopy
INTEGER,                  INTENT(IN)    :: KSNW      ! number of drifting snow variables
REAL,                     INTENT(IN)    :: PTSTEP    ! time-step                             (s)
REAL, DIMENSION(KI,KSNW),      INTENT(IN)    :: PSNWA ! Blowing snow variables at forcing levels  (__ /kg)
REAL, DIMENSION(KI,KLVL), INTENT(IN)         :: PK        ! mixing exchange coefficient           (m2/s)
REAL, DIMENSION(KI,KSNW),      INTENT(INOUT)    :: PSFLUX_SNW  ! surface flux w'snw'             (mkg/kg/s)
REAL, DIMENSION(KI,KLVL), INTENT(IN)            :: PZ       ! heights of canopy levels    (m)
REAL, DIMENSION(KI,KLVL), INTENT(IN)            :: PDZ       ! deltaZ between canopy half levels,
!                                                            ! located at full levels                (m)
REAL, DIMENSION(KI,KLVL), INTENT(IN)            :: PDZF      ! deltaZ between canopy (full) levels,
!                                                             ! located at half levels                (m)
REAL, DIMENSION(KI,KLVL,KSNW), INTENT(INOUT)    :: PSNW  ! drifting snow at canopy levels           (kg/kg)
REAL, DIMENSION(KI,KLVL,KSNW), INTENT(IN)    :: PSNW_SUBL  ! drifting snow sublimation rate at canopy levels   (kg/kg/s)
                                                           ! zero when sublimation is not activated
REAL, DIMENSION(KI), INTENT(IN)                 :: PRHOA      !  air density (kg/m3)
REAL,  DIMENSION(KI,KLVL),    INTENT(IN)        :: PTHT       !  air temperature (K)
REAL,  DIMENSION(KI,KLVL),    INTENT(IN)        :: PPABST     !  air pressure (Pa)
REAL,  DIMENSION(KI,KLVL),    INTENT(IN)        :: PTKE     !  air pressure (Pa)
REAL,  DIMENSION(KI,KSNW), INTENT(OUT)          :: PSNWDEP  ! sedimentation flux at the bottom of Canopy

!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                     :: JLAYER   ! loop counter on layers
INTEGER                     :: JSV      ! loop counter on blowing snow variables
INTEGER                     :: JI       ! loop counter on grid points
!
REAL, DIMENSION(KI)           :: ZZREF  ! height of forcing level
REAL, DIMENSION(KI,KSNW)      :: ZSNWC  ! Blowing Snow content in canopy layers     (__ /m3)

REAL, DIMENSION(KI,KLVL)      :: ZBET   ! Scale parameter of the gamma distribution (m)
REAL, DIMENSION(KI,KLVL,KSNW) :: ZVGK  ! Terminal fallspeed (m/s)
REAL, DIMENSION(KI,KLVL,KSNW) :: ZKSNW  ! particle eddy diffusivity 
REAL, DIMENSION(KI,KLVL)      :: ZRG    ! Mean radius  (m)
REAL, DIMENSION(KI,KLVL-1,KSNW) :: ZA,ZB,ZC,ZRHS   ! Term of the tridiagonal matrix  
REAL, DIMENSION(KI,KLVL,KSNW) :: ZSNW      ! work variable : pot. temp at futur instant 
!                                       ! (or past at the end of the routine) 
LOGICAL, DIMENSION(KI)        :: GEMIS  ! Logical=TRUE if snow is emitted in the atmosphere for the 
                                        ! first time 
INTEGER    ::    IPRINT
!
!-------------------------------------------------------------------------------
!
!
!*    1. initializations
!    ---------------
!
IPRINT = 1
!
GEMIS(:)= .FALSE.
!
!-------------------------------------------------------------------------------
!
!
!*    2. Compute settling velocity based on the size distribution of the previous time
!        step
!        ---------------
!
!-------------------------------------------------------------------------------
!
!           Convert canopy variables in _/m3
!
DO JLAYER=1,KLVL
  DO JSV=1,KSNW
     ZSNW(:,JLAYER,JSV)=PSNW(:,JLAYER,JSV)*PRHOA(:)
  ENDDO
ENDDO

!
!*         compute BETA, RG and moments
!
CALL SNOWMOMENT2SIZE(ZSNW, PBETA1D=ZBET, PRG1D=ZRG )
!
!*          Initialize profile of radius if no snow is present in the atmosphere
!           We use Pomeroy's theoretical profile.  
!

DO JI=1,KI
   IF(PSFLUX_SNW(JI,2)>0 .AND. ALL(PSNW(JI,:,2)==0.)) THEN
           GEMIS(JI) = .TRUE.
           DO JLAYER=1,KLVL
               ZRG(JI,JLAYER)=XEMIRADIUS_SNW*(PZ(JI,JLAYER)/0.05)**(-0.258)
               ZBET(JI,JLAYER) = ZRG(JI,JLAYER)/XEMIALPHA_SNW
           END DO    
   END IF
END DO
!*         compute gravitational velocities
!
!
CALL BLOWSNW_VELGRAV1D(ZBET, ZRG, PTHT, PRHOA,PPABST, ZVGK)

!
! Mettre en place cas spécial avec flux positif et pas de neige dans l'atmo 
!  initialisation du profil de rayon et de vitesse de chute suivant la méthode
!  de Déry
!
!*         compute particle eddy diffusivity
!
!
CALL BLOWSNW_DIFFK(PK, PTKE, ZVGK,KSNW, ZKSNW)
!
!
!-------------------------------------------------------------------------------
!
!
!*    3. Initialize coefficients of the tridiagonal matrix
!        ---------------
!
!-------------------------------------------------------------------------------
!
!       Lower boundary: turbulent flux PSFLUX_SNW at the surface
!
DO JSV = 1,KSNW
    ZA(:,1,JSV)  = 0.
    ZB(:,1,JSV)  = 1+(ZKSNW(:,2,JSV)/PDZF(:,2)+ZVGK(:,1,JSV))*PTSTEP/PDZ(:,1)
    ZC(:,1,JSV)  = -(ZKSNW(:,2,JSV)/PDZF(:,2)+ZVGK(:,2,JSV))*PTSTEP/PDZ(:,1)
    ZRHS(:,1,JSV)=  PSNW(:,1,JSV)+PSFLUX_SNW(:,JSV)*PTSTEP/PDZ(:,1)+  &
                   PSNW_SUBL(:,1,JSV)*PTSTEP
END DO
!
!       Upper boundary: imposed concentration at level KLVL
!
DO JSV = 1,KSNW
    ZC(:,KLVL-1,JSV)  = 0.
    ZRHS(:,KLVL-1,JSV)= PSNW(:,KLVL-1,JSV)+PSNW_SUBL(:,KLVL-1,JSV)*PTSTEP+     &
                       PSNW(:,KLVL,JSV)*PTSTEP/                                &
                       PDZ(:,KLVL-1)*(ZVGK(:,KLVL,JSV)+ZKSNW(:,KLVL,JSV)/PDZF(:,KLVL))
END DO
!
!       Values at other levels
!
DO JSV = 1,KSNW
  DO JLAYER=2,(KLVL-1)
     ZA(:,JLAYER,JSV)= -(ZKSNW(:,JLAYER,JSV)/PDZF(:,JLAYER))*PTSTEP/PDZ(:,JLAYER)
     ZB(:,JLAYER,JSV)=1+PTSTEP/PDZ(:,JLAYER)*(ZKSNW(:,JLAYER+1,JSV)/PDZF(:,JLAYER+1)+  &
                      ZKSNW(:,JLAYER,JSV)/PDZF(:,JLAYER)+ZVGK(:,JLAYER,JSV))
     IF(JLAYER<(KLVL-1)) THEN
             ZC(:,JLAYER,JSV) = -(ZKSNW(:,JLAYER+1,JSV)/PDZF(:,JLAYER+1)+ZVGK(:,JLAYER+1,JSV)) &
                                  *PTSTEP/PDZ(:,JLAYER)
             ZRHS(:,JLAYER,JSV) = PSNW(:,JLAYER,JSV)+PSNW_SUBL(:,JLAYER,JSV)*PTSTEP
     END IF   
   END DO     
END DO
!
!-------------------------------------------------------------------------------
!
!
!*    4. Solve tridiagonal system
!        ---------------
!
!-------------------------------------------------------------------------------
!
DO JSV = 1,KSNW

   CALL TRIDIAG_GROUND(ZA(:,:,JSV),ZB(:,:,JSV),ZC(:,:,JSV),ZRHS(:,:,JSV), &
                      ZSNW(:,1:KLVL-1,JSV))
END DO
!
!
!-------------------------------------------------------------------------------
!
!*    5. updated fluxes 
!        ----------------------
!
! Compute net flux sent to Meso-NH: use top of canopy instead of surface
!
DO JI=1,KI
  DO JSV=1,KSNW
       IF(GEMIS(JI)) THEN 
              PSFLUX_SNW(JI,JSV) =   MAX(-ZKSNW (JI,KLVL,JSV) * (PSNW(JI,KLVL,JSV)-ZSNW(JI,KLVL-1,JSV))/PDZF(JI,KLVL)* &
                     (1-ZVGK(JI,KLVL,JSV)*PTSTEP/(2*PZ(JI,KLVL))),0.)
       ELSE
              PSFLUX_SNW(JI,JSV) =   -ZKSNW (JI,KLVL,JSV) * (PSNW(JI,KLVL,JSV)-ZSNW(JI,KLVL-1,JSV))/PDZF(JI,KLVL)- &
                                   ZVGK(JI,KLVL,JSV)*PSNW(JI,KLVL,JSV)
       ENDIF
       PSFLUX_SNW(JI,JSV) =  PSFLUX_SNW(JI,JSV)*PRHOA(JI)  ! Convert flux in kg_{snow}/m2/s
  END DO
END DO
!
! Update sedimentation flux send to Crocus at the next time step
!

DO JSV=1,KSNW
  PSNWDEP(:,JSV) = ZVGK(:,1,JSV)*ZSNW(:,1,JSV)*PRHOA(:)
END DO
!
! Store sedimentation flux at the bottom of Canopy as a diagnostic
!XSNSED(:,1:2) =  -PSNWDEP(:,1:2)   ! Instantaneous fluxes (number and mass) 
!XSNSED(:,3) =  XSNSED(:,3)+XSNSED(:,2)*PTSTEP  ! Accumulated flux (mass)


!
!-------------------------------------------------------------------------------
!
!*    6. New value of blowing snow variables (at full levels)
!        ----------------------------------
!
PSNW(:,1:(KLVL-1),:) = ZSNW(:,1:(KLVL-1),:)
!
!-------------------------------------------------------------------------------
!
!*    7. Update total Canopy variables
!        ----------------------------------
!
!ZSNWC =0.
! Compute number and mass content of Canopy layer.  
!DO JI=1,KI
!  ZZREF(JI) = SUM(PDZ(JI,1:(KLVL-1)))
!  DO JSV=1,KSNW
!     DO JLAYER=1,(KLVL-1)
!        ZSNWC(JI,JSV) = ZSNWC(JI,JSV)+PSNW(JI,JLAYER,JSV)*PDZ(JI,JLAYER)*PRHOA(JI)
!     END DO
!     PSNOW_CANO(JI,JSV)=ZSNWC(JI,JSV)/ZZREF(JI)
!  END DO
!END DO
!

END SUBROUTINE CANOPY_EVOL_BLOWSNW
