SUBROUTINE SNOWPACK_EVOL(HSNOWRES,PBLOWSNW_FLUX,PSNOWHEAT,PSNOWSWE,PSNOWRHO,   &
                        PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE,              &
                        PTSTEP,PRHOA,PTA,PBLOWSNW_CONC,                        &
                        PVMOD,PQA,PPS,PUREF,PEXNS,PDIRCOSZW,PZREF,             &
                        PZ0EFF,PZ0H,PBLOWSNW_DEP,PTG)
!     ######################################################################
!     ##########################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to simulate the mass exchange between the
!!      snowpack and the atmosphere when coupled to Meso-NH. A net mass flux flux is computed.
!!           - If erosion occurs, it is simumated by this routine.
!!           - If deposition occurs, the deposited mass flux is computed, sent to Crocus
!!                and is added to the snowfall
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    Vionnet et al, The Cryosphere, 2014
!!
!!    AUTHOR
!!    ------
!!      V. Vionnet      * CNRM/GAME*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/10/11
!
!*       0.    DECLARATIONS

USE MODD_CSTS,        ONLY : XLMTT, XTT, XCI,XRHOLW
USE MODD_SNOW_PAR,    ONLY : X_RI_MAX
USE MODD_BLOWSNW_SURF
USE MODD_BLOWSNW_n
USE MODD_SNOW_METAMO

USE MODE_SNOW3L
USE MODE_THERMOS

USE MODI_SURFACE_RI
USE MODI_SURFACE_CD

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE
! ===============================================================
!
!       0.1 declarations of arguments        
!  
! ===============================================================
! 
    CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES
!                                      HSNOWRES  = ISBA-SNOW3L turbulant exchange option
!                                      'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                                      'RIL' = Limit Richarson number under very stable 
!                                              conditions

    REAL, DIMENSION(:), INTENT(IN) :: PTA,PRHOA,PVMOD,PZ0EFF,PQA,PEXNS, &
                                      PUREF,PZREF,PDIRCOSZW,PPS,PZ0H,PTG
!                                      PTA     = atmospheric temperature at level za (K)
!                                      PRHOA   = air density (kg/m3)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PZ0EFF    = roughness length for momentum 
!                                      PQA     = atmospheric specific humidity
!                                                at level za
!                                      PEXNS     = Exner function
!                                      PUREF     = reference height of the wind
!                                      PZREF     = reference height of the first
!                                      PDIRCOSZW = Cosinus of the angle between the 
!                                                  normal to the surface and the vertical
!                                                  atmospheric level
!                                      PPS     = surface pressure
!                                      PZ0H   = grid box average roughness length for heat
!                                      PTG       = Surface soil temperature (effective 
!                                                  temperature of the layer lying below snowpack)

    REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWHEAT,PSNOWSWE,PSNOWRHO
!                                      PSNOWHEAT = Snow layer(s) heat content (J/m3)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)

    REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1,PSNOWGRAN2,       &
                                         PSNOWHIST,PSNOWAGE
!                                      PSNOWGRAN1 = Snow layer(s) grain parameter 1
!                                      PSNOWGRAN2 = Snow layer(s) grain parameter 2
!                                      PSNOWHIST  = Snow layer(s) grain historical parameter
!                                      PSNOWAGE   = Snow layer(s) age

    REAL, DIMENSION(:,:)   , INTENT(IN)   :: PBLOWSNW_CONC
!                                      PBLOWSNW_CONC = blowing snow particle concentration
!                                                   (1: number; 2: mass)

    REAL, DIMENSION(:,:)   , INTENT(INOUT)   :: PBLOWSNW_FLUX
!                                      PBLOWSNW_FLUX  = Blowing snow particles flux
!                                      IN : sedimentation flux from previous time step (1,2)
!											contribution of saltation transport to snowpack mass balance (3)
!                                                   1: number (#/m2/s); 2: mass (kg/m2/s), 3: (kg/m2/s)
!                                      OUT: emitted turbulent flux + streamwise saltation flux
!                                                   1: number (#/m2/s); 2: mass (kg/m2/s), 3: kg/m/s
!
    REAL, DIMENSION(:)   , INTENT(OUT)   :: PBLOWSNW_DEP
!                                      PBLOWSNW_DEP    = deposion flux of blowing snow particles (sent to Crocus)
!
    REAL, INTENT(IN)               :: PTSTEP
!                                      PTSTEP    = time step of the integration

! ===============================================================
!       
!       0.2 declarations of local variables
!
! ===============================================================
    INTEGER JWRK, JJ   ! Loop counters
    INTEGER INLVLS     ! maximum number of snow layers
    REAL    SWE_MAX    ! Maximum SWE that can be removed from a snow layer
    REAL    SWE_DEP    ! SWE which is deposited
    REAL    ZMOB       ! Mobility index of the eroded layer


    REAL, DIMENSION(SIZE(PRHOA)) ::    T_EROSION  ! Erosion time step

    INTEGER, DIMENSION(SIZE(PRHOA)) :: INLVLS_USE ! decrit nbre de couches effectives initiales
    INTEGER, DIMENSION(SIZE(PRHOA)) :: INLVLS_USE_TMP ! decrit nbre de couches effectives en cours d'�volution

    LOGICAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: LAYER_USE !
    LOGICAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: LAYER_USE_TMP !

    REAL   , DIMENSION(SIZE(PRHOA),2) :: ZTURBFLUX  ! blowing snow turbulent flux
    REAL   , DIMENSION(SIZE(PRHOA),2) :: ZSEDFLUX   ! blowing snow sedimentation flux
    REAL   , DIMENSION(SIZE(PRHOA)) :: ZSALT_CONTR  ! contribution of saltation transport to snowpack mass balance
    REAL   , DIMENSION(SIZE(PRHOA),2) :: ZTOTFLUX_TURB ! total turbulent blowing snow flux
    REAL   , DIMENSION(SIZE(PRHOA))   :: ZTOTFLUX_SALT ! total saltation mass flux
    REAL   , DIMENSION(SIZE(PRHOA))   :: ZSALTFLUX  ! streamwise saltation mass flux
    REAL   , DIMENSION(SIZE(PRHOA))  :: ZEMIMOB     ! emitted mobility index
    REAL   , DIMENSION(SIZE(PRHOA))  :: ZDEPMOB     ! deposited mobility index
    REAL   , DIMENSION(SIZE(PRHOA)) :: ZUSTAR     ! friction velocity
    REAL   , DIMENSION(SIZE(PRHOA)) :: ZRISNOW     ! Richardson number over snow
    REAL   , DIMENSION(SIZE(PRHOA)) :: ZCDSNOW     ! Drag coefficient for momentum over snow


    REAL, DIMENSION(SIZE(PRHOA)) :: ZEXNA, ZQSAT,  ZCDN,ZTEMP_RI

    ! Temporary snow variables
    REAL   , DIMENSION(SIZE(PRHOA,1)) :: ZSNOW     ! total snow depth

    REAL,DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))  :: ZSNOWGRAN1,ZSNOWGRAN2, &
                                         ZSNOWHIST,ZSNOWAGE

    REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEAT,ZSNOWDZ,      &
                                            ZSNOWSWE,ZSNOWRHO,ZSNOWTEMP,ZSNOWLIQ, &
                                            ZSCAP

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
! ===============================================================
!
!       1. Initialization of local variables
!
! ===============================================================
!
IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL',0,ZHOOK_HANDLE)

INLVLS       = SIZE(PSNOWSWE(:,:),2)



ZSEDFLUX(:,:)= -PBLOWSNW_FLUX(:,1:2)
ZSALT_CONTR(:) = -PBLOWSNW_FLUX(:,3)

ZTOTFLUX_TURB(:,:)=0.

ZTOTFLUX_SALT(:)=0.
ZSALTFLUX(:)=0.

ZEMIMOB(:) =0.
ZDEPMOB(:) =0.

!       Initialize snowpack characteristics
ZSNOWGRAN1(:,:)=PSNOWGRAN1(:,:)
ZSNOWGRAN2(:,:)=PSNOWGRAN2(:,:)
ZSNOWHIST(:,:) =PSNOWHIST(:,:)
ZSNOWAGE(:,:)  =PSNOWAGE(:,:)
ZSNOWHEAT(:,:) =PSNOWHEAT(:,:)
ZSNOWRHO(:,:)  =PSNOWRHO(:,:)
ZSNOWSWE(:,:)  =PSNOWSWE(:,:)


! Initially : no deposition of blowing snow particles
PBLOWSNW_DEP(:) = 0.


WHERE (PSNOWSWE(:,:)>0.) 
        LAYER_USE(:,:)=.TRUE.
        ZSNOWDZ(:,:) = PSNOWSWE(:,:)/PSNOWRHO(:,:)
ELSEWHERE
        LAYER_USE(:,:)=.FALSE.
        ZSNOWDZ(:,:) =0.
ENDWHERE

ZSNOWTEMP(:,:)= 0.

! Initialize temperature where snow is present
WHERE (LAYER_USE(:,:))
ZSCAP(:,:)     = SNOW3LSCAP(PSNOWRHO)
!
ZSNOWTEMP(:,:) = XTT + ( ((PSNOWHEAT(:,:)/ZSNOWDZ(:,:))                   &
                 + XLMTT*PSNOWRHO(:,:))/ZSCAP(:,:) )
!
ZSNOWLIQ(:,:)  = MAX(0.0,ZSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*                  &
                 ZSNOWDZ(:,:)/(XLMTT*XRHOLW)
!
ZSNOWTEMP(:,:) = MIN(XTT,ZSNOWTEMP(:,:))
END WHERE

! Temperature used to compute Richardson number depending if the soil is
! covered by snow or not. 
ZTEMP_RI(:) = 0.
WHERE(ZSNOWTEMP(:,1)>0.)
     ZTEMP_RI(:) = ZSNOWTEMP(:,1)
ELSEWHERE
     ZTEMP_RI(:) = PTG(:)
END WHERE


ZQSAT(:)    = QSATI(ZTEMP_RI(:),PPS(:))
ZEXNA(:)  = PEXNS(:)
!
!      Compute Richardson number
!
CALL SURFACE_RI(ZTEMP_RI(:), ZQSAT, PEXNS, ZEXNA, PTA, PQA,                  &
                PZREF, PUREF, PDIRCOSZW, PVMOD, ZRISNOW                  )
!
! Simple adaptation of method by Martin and Lejeune (1998)
! to apply a lower limit to turbulent transfer coef
! by defining a maximum Richarson number for stable
! conditions:
!
IF(HSNOWRES=='RIL') ZRISNOW(:) = MIN(X_RI_MAX, ZRISNOW(:))
!
! Compute drag coefficient for momentum
!
CALL SURFACE_CD(ZRISNOW, PZREF, PUREF, PZ0EFF, PZ0H, ZCDSNOW, ZCDN)
!
! Compute friction velocity
!
ZUSTAR(:) = SQRT(ZCDSNOW(:)*PVMOD(:)**2)
!
!
! Initialize snow layering
!


INLVLS_USE(:)=0
DO JJ=1,INLVLS
    WHERE(LAYER_USE(:,JJ)) INLVLS_USE(:)=JJ
ENDDO   
INLVLS_USE_TMP(:)=INLVLS_USE(:)
LAYER_USE_TMP(:,:)=LAYER_USE(:,:)

! ===============================================================
!
!       2. Compute net flux of blowing snow particles
!         and update of snowpack accordingly
!
! ===============================================================
!       Loop on snow covered points
         DO JWRK=1,SIZE(PSNOWSWE,1)
           T_EROSION(JWRK)=PTSTEP
           DO JJ=1,INLVLS_USE(JWRK)
              ZSNOW(JWRK) = SUM(ZSNOWDZ(JWRK,1:INLVLS_USE_TMP(JWRK)))

!             Turbulent flux              
              IF(LSNOW_WIND) THEN
              ! Buoyancy effects in presence of suspended particles are taken into account.
              ! Iterations are necessary to reach a balance between u* and
              ! concentration in the saltation layer. 
              CALL SNOWFLUX_GET_TURB(ZSNOWGRAN1(JWRK,1),ZSNOWGRAN2(JWRK,1),               &
                        ZSNOWLIQ(JWRK,1),ZSNOWHIST(JWRK,1),PRHOA(JWRK),ZTURBFLUX(JWRK,:), & 
                        PBLOWSNW_CONC(JWRK,:),PVMOD(JWRK),ZCDSNOW(JWRK),PZ0EFF(JWRK),ZMOB,   &
                        ZRISNOW(JWRK),PZREF(JWRK),PZ0H(JWRK),PUREF(JWRK), HSNOWRES,       &
                        PDIRCOSZW(JWRK),ZSALTFLUX(JWRK))
              ELSE
              ! Buoyancy effects in presence of suspended particles are not
              ! taken into account. Only thermal effect are included. Snow turbulent flux is directly computed.
              CALL SNOWFLUX_GET(ZSNOWGRAN1(JWRK,1),ZSNOWGRAN2(JWRK,1),                                       &
                        ZSNOWLIQ(JWRK,1),ZSNOWHIST(JWRK,1),ZUSTAR(JWRK),PRHOA(JWRK),ZTURBFLUX(JWRK,:),       & 
                        PBLOWSNW_CONC(JWRK,:),PVMOD(JWRK),ZCDSNOW(JWRK),PZ0EFF(JWRK),ZMOB,ZSALTFLUX(JWRK))
              END IF
              ! Negative net flux : snow deposition
              IF(ZTURBFLUX(JWRK,2)+ZSEDFLUX(JWRK,2)+ZSALT_CONTR(JWRK)<=0.) THEN
                  SWE_DEP=-(ZTURBFLUX(JWRK,2)+ZSEDFLUX(JWRK,2)+ZSALT_CONTR(JWRK))*T_EROSION(JWRK)
                  PBLOWSNW_DEP(JWRK) = SWE_DEP / PTSTEP
                  ZTOTFLUX_TURB(JWRK,:)=ZTOTFLUX_TURB(JWRK,:)+T_EROSION(JWRK)/PTSTEP*ZTURBFLUX(JWRK,:)
                  ZTOTFLUX_SALT(JWRK)=ZTOTFLUX_SALT(JWRK)+T_EROSION(JWRK)/PTSTEP*ZSALTFLUX(JWRK)
                  ZEMIMOB(JWRK) = ZEMIMOB(JWRK)+T_EROSION(JWRK)/PTSTEP*ZMOB
                  T_EROSION(JWRK)=0.
                  EXIT

              ! Positive net flux : snow erosion                 
              ELSE IF(ZTURBFLUX(JWRK,2)+ZSEDFLUX(JWRK,2)+ZSALT_CONTR(JWRK)>0.) THEN
                  ! SWE potentially eroded during the T_EROSION
                  SWE_MAX=(ZTURBFLUX(JWRK,2)+ZSEDFLUX(JWRK,2)+ZSALT_CONTR(JWRK))*T_EROSION(JWRK)

                  ! Only layer JJ is partially eroded
                  IF(SWE_MAX<=ZSNOWSWE(JWRK,1)) THEN   ! Only layer JJ is partially eroded
                      CALL SNOWEROSION(SWE_MAX,ZSNOWHEAT(JWRK,:),ZSNOWDZ(JWRK,:),ZSNOWSWE(JWRK,:), &
                                ZSNOWRHO(JWRK,:),ZSNOWGRAN1(JWRK,:),ZSNOWGRAN2(JWRK,:)             &
                                ,ZSNOWHIST(JWRK,:),ZSNOWAGE(JWRK,:),ZSNOWLIQ(JWRK,:),              &
                                INLVLS_USE_TMP(JWRK),LAYER_USE_TMP(JWRK,:),PVMOD(JWRK),PTSTEP)
                      ZTOTFLUX_TURB(JWRK,:)=ZTOTFLUX_TURB(JWRK,:)+T_EROSION(JWRK)/PTSTEP*ZTURBFLUX(JWRK,:)
                      ZTOTFLUX_SALT(JWRK)=ZTOTFLUX_SALT(JWRK)+T_EROSION(JWRK)/PTSTEP*ZSALTFLUX(JWRK)
                      ZEMIMOB(JWRK) = ZEMIMOB(JWRK)+T_EROSION(JWRK)/PTSTEP*ZMOB
                      T_EROSION(JWRK)=0.
                      EXIT

                  ! Layer JJ is totally eroded and removed from the snowpack, next layer is then checked
                  ! to see if blowing snow occurs
                  ELSE                            
                      T_EROSION(JWRK)=T_EROSION(JWRK)*(1-ZSNOWSWE(JWRK,1)/SWE_MAX) 
                      ZTOTFLUX_TURB(JWRK,:)=ZTOTFLUX_TURB(JWRK,:)+(1-T_EROSION(JWRK)/PTSTEP)*ZTURBFLUX(JWRK,:)
                      ZTOTFLUX_SALT(JWRK)=ZTOTFLUX_SALT(JWRK)+(1-T_EROSION(JWRK)/PTSTEP)*ZSALTFLUX(JWRK)  
                      CALL SNOWEROSION(ZSNOWSWE(JWRK,1),ZSNOWHEAT(JWRK,:),ZSNOWDZ(JWRK,:),         &
                           ZSNOWSWE(JWRK,:),  ZSNOWRHO(JWRK,:),ZSNOWGRAN1(JWRK,:),                 &
                           ZSNOWGRAN2(JWRK,:),ZSNOWHIST(JWRK,:),ZSNOWAGE(JWRK,:),ZSNOWLIQ(JWRK,:), &
                           INLVLS_USE_TMP(JWRK),LAYER_USE_TMP(JWRK,:),PVMOD(JWRK),PTSTEP)
                      ZEMIMOB(JWRK) = ZEMIMOB(JWRK)+(1-T_EROSION(JWRK)/PTSTEP)*ZMOB    
                  END IF
              END IF 
           END DO


           IF(T_EROSION(JWRK)>0.AND.(ZSEDFLUX(JWRK,2)+ZSALT_CONTR(JWRK))<-XUEPSI) THEN
                  ! Take into account two cases : 
                  !        -   Snow pack has been totally eroded during less
                  !              than one time step and potential deposition is
                  !              computed for the remaining part of the time step
                  !        -   Deposition of blown snow particles on bare ground
              !write(*,*) 'Deposition en excedent'   

             SWE_DEP=-(ZSEDFLUX(JWRK,2)+ZSALT_CONTR(JWRK))*T_EROSION(JWRK)
             PBLOWSNW_DEP(JWRK) = SWE_DEP / PTSTEP
                 write(*,*)'pb T_EROSION',T_EROSION(JWRK),JWRK
                 write(*,*)'pb INLVLS_OLD',SWE_DEP
                 write(*,*) ZTURBFLUX(JWRK,2)
!                 STOP 'deposition en excedent'
            ENDIF

             ! Compute final snow depth that may have changed because of erosion.
           IF(INLVLS_USE_TMP(JWRK)>0) THEN
                ZSNOW(JWRK) = SUM(ZSNOWDZ(JWRK,1:INLVLS_USE_TMP(JWRK)))
           ELSE
                ZSNOW(JWRK) = 0.
           ENDIF

        END DO


! ===============================================================
!
!       3. Final computation 
!
! ===============================================================

! Update final snowpack
        PSNOWGRAN1(:,:)=ZSNOWGRAN1(:,:)
        PSNOWGRAN2(:,:)=ZSNOWGRAN2(:,:)
        PSNOWHIST (:,:)=ZSNOWHIST (:,:)
        PSNOWAGE  (:,:)=ZSNOWAGE  (:,:)
        PSNOWHEAT (:,:)=ZSNOWHEAT (:,:)
        PSNOWRHO  (:,:)=ZSNOWRHO  (:,:)
        PSNOWSWE  (:,:)=ZSNOWDZ(:,:)*PSNOWRHO(:,:)

! Snow fluxes sent to Canopy or directly to Meso-NH if Canopy is activated

        PBLOWSNW_FLUX(:,1)=ZTOTFLUX_TURB(:,1)          ! Save total number turbulent flux
        PBLOWSNW_FLUX(:,2)=ZTOTFLUX_TURB(:,2)          ! Save total mass turbulent flux
        PBLOWSNW_FLUX(:,3)=ZTOTFLUX_SALT(:)            ! Save streamwise saltation flux

IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL',1,ZHOOK_HANDLE)

CONTAINS

        SUBROUTINE SNOWFLUX_GET(PSNOWGRAN1,PSNOWGRAN2,                      &
                        PSNOWLIQ,PSNOWHIST,PUSTAR,PRHOA,PBLOWSNW_FLUX,      &
                        PBLOWSNW_CONC,PVMOD,PCDSNOW,PZ0EFF,PMOB,PQSALT)
!
!!    PURPOSE
!!    -------
!     Calculates turbulent blowing snow flux and streamwise saltation flux as a function
!     of wind speed and snowpack properties
!
!       Returns :
!          - emitted turbulent flux
!          - mobility index of emitted snow
!          - streamwise saltation flux
!
!
     USE MODE_BLOWSNW_MBLUTL

     IMPLICIT NONE
!
!       0.1 declarations of arguments        
!             
!                                                    
        REAL,  INTENT(IN)    :: PSNOWGRAN1,PSNOWGRAN2,PSNOWLIQ,PSNOWHIST
!                                                                       
        REAL,  INTENT(IN)  :: PUSTAR,PRHOA
        REAL,  INTENT(IN)  :: PVMOD,PCDSNOW,PZ0EFF
        REAL, DIMENSION(:)   , INTENT(IN)   :: PBLOWSNW_CONC
       
        REAL,  INTENT(OUT) :: PMOB
        REAL, DIMENSION(:)   , INTENT(OUT)   :: PBLOWSNW_FLUX
        REAL, INTENT(OUT)   :: PQSALT ! Streamwise saltation mass flux (kg/m/s)
!
!       0.2 declaration of local variables
!        
        REAL ::  ZMBLINDEX
        REAL ::  ZWIND_RFR_THR,ZWIND_FRC,ZWIND_FRC_THR,            &
                                              ZWIND_FRC_SALT
                                              
        REAL ::  ZCONC_SALT

        REAL(KIND=JPRB) :: ZHOOK_HANDLE
!       
!      0.3     Initialization
!
IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL:SNOWFLUX_GET',0,ZHOOK_HANDLE)

!                            
        ZWIND_FRC = PUSTAR
        PMOB = SNOW_MBL_INDEX(PSNOWGRAN1, PSNOWGRAN2)
!
!       1.     Threshold Wind at 5m
!        
        CALL WIND_RFR_THR(ZWIND_RFR_THR,PSNOWGRAN1,PSNOWGRAN2)
!
!       2.     Threshold Friction velocity
!        
        CALL WIND_FRC_THR(ZWIND_RFR_THR,PZ0EFF,ZWIND_FRC_THR)
!                      
!       3. Blown snow particles concentration (kg/m3) in the saltation layer
        
!          
        CALL CONC_SALT(ZCONC_SALT,ZWIND_FRC,ZWIND_FRC_THR,PRHOA,PQSALT, &
                     PSNOWLIQ,PSNOWHIST)
!
!       4. Turbulent flux of blowing snow particles     
!                                
        CALL SNOW_FLUX(ZCONC_SALT,PVMOD,PBLOWSNW_CONC,PCDSNOW,PBLOWSNW_FLUX)
!    
!       Snowdrfit does not occur when: 
!          - surface layer is wet
!          - historical variable is higher than 1 : crust of non-transportable snow
        IF(PSNOWLIQ>0 .OR. PSNOWHIST>0) THEN 
           PBLOWSNW_FLUX(:)=0
        ENDIF          
!

IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL:SNOWFLUX_GET',1,ZHOOK_HANDLE)

       END SUBROUTINE SNOWFLUX_GET 


        SUBROUTINE SNOWFLUX_GET_TURB(PSNOWGRAN1,PSNOWGRAN2,                  &
                        PSNOWLIQ,PSNOWHIST,PRHOA,PBLOWSNW_FLUX,                  &
                        PBLOWSNW_CONC,PVMOD,PCDSNOW,PZ0EFF,PMOB,PRI,PZREF,PZ0H, &
                        PUREF, HSNOWRES, PDIRCOSZW,PQSALT )

     USE MODE_BLOWSNW_MBLUTL
     USE MODD_BLOWSNW_SURF

IMPLICIT NONE
!
!
!!    PURPOSE
!!    -------
!       Routine that compute emitted turbulent flux accounting for influence of 
!       blown snow particles on the wind profile. 
!       Return : 
!          - emitted turbulent flux
!          - mobility index of emitted snow   
!          - streamwise saltation flux
!          - Richardson number including particle-induced buoyancy
!          - modified drag coefficient for momentum
! 
!      Routine to be tested!!
!
!       0.1 declarations of arguments        
!             
!                                                    
        REAL,  INTENT(IN)    :: PSNOWGRAN1,PSNOWGRAN2,PSNOWLIQ,PSNOWHIST
!                                                                       
        REAL,  INTENT(IN)    :: PRHOA
        REAL,  INTENT(IN)    :: PVMOD,PZ0EFF
        REAL,  INTENT(IN)    :: PZREF, PZ0H, PUREF, PDIRCOSZW
        REAL, DIMENSION(:)   , INTENT(IN)   :: PBLOWSNW_CONC

        CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES          
!
        REAL,  INTENT(INOUT) :: PCDSNOW, PRI
!
        REAL,  INTENT(OUT) :: PMOB
        REAL, DIMENSION(:)   , INTENT(OUT)   :: PBLOWSNW_FLUX
        REAL, INTENT(OUT)   :: PQSALT ! Streamwise saltation mass flux (kg/m/s)       
!
!       0.2 declaration of local variables
!        
        REAL ::  ZMBLINDEX
        REAL ::  ZWIND_RFR_THR,ZWIND_FRC,ZWIND_FRC_THR
        REAL ::  ZCD_TEMP, ZCDN 
        REAL ::  ZCONC_SALT
        REAL ::  ZRI_TOT       ! Total richardson number (thermal + particle
!        buoyancy included) 
        REAL ::  ZRI_PART      ! Particle richardson number ( particle
!        buoyancy included)
        REAL :: ZB,ZDELTA,ZY,ZY_DRIV,ZY_TEMP
        INTEGER :: II,JJ
        INTEGER :: NITER

        REAL(KIND=JPRB) :: ZHOOK_HANDLE


!       
!      0.3     Initialization             
!
IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL:SNOWFLUX_GET_TURB',0,ZHOOK_HANDLE)
!
        ZWIND_FRC = SQRT(PCDSNOW*PVMOD**2)
        ZCD_TEMP = PCDSNOW  ! Drag coefficient without particle-induced buoyancy
!        and increase in aerodynamic rougnhess due to saltation
        ZRI_TOT = PRI       ! Richardson number without particle-induced buoyancy
!       Parameters for Newton's method
        ZDELTA = 0.01
        JJ=0
        ZB = ZWIND_FRC-999
        NITER = 5  ! Number of iterations to compute u* and Conc_salt
        !
        ! Compute mobility index based on surface layer properties
        PMOB = SNOW_MBL_INDEX(PSNOWGRAN1, PSNOWGRAN2);
!
!       1.     Threshold Wind at 5m
!        
        CALL WIND_RFR_THR(ZWIND_RFR_THR,PSNOWGRAN1,PSNOWGRAN2)
!
!       2.     Threshold Friction velocity
!
        CALL WIND_FRC_THR(ZWIND_RFR_THR,PZ0EFF,ZWIND_FRC_THR)
!
!       3. Newton method to find u*, Ri and CD     
!
DO WHILE(ABS(ZWIND_FRC-ZB)>0.001)
      JJ=JJ+1
      ZB=ZWIND_FRC
      IF(JJ>NITER) EXIT
      CALL SOLVE_CD(ZWIND_FRC+ZDELTA,PVMOD,PRHOA,PBLOWSNW_CONC(2),PZREF,PUREF,PDIRCOSZW,       &
                   ZWIND_FRC_THR,PRI,PZ0EFF,HSNOWRES,PZ0H,ZY_TEMP,ZRI_TOT,ZCD_TEMP,         &
                               PSNOWLIQ,PSNOWHIST)

      CALL SOLVE_CD(ZWIND_FRC,PVMOD,PRHOA,PBLOWSNW_CONC(2),PZREF,PUREF,PDIRCOSZW,              &
                   ZWIND_FRC_THR,PRI,PZ0EFF,HSNOWRES,PZ0H,ZY,ZRI_TOT,ZCD_TEMP,              &
                             PSNOWLIQ,PSNOWHIST)

      ZY_DRIV = (ZY_TEMP-ZY)/ZDELTA 
      ZWIND_FRC = ZB- ZY/ZY_DRIV
END DO      
!
!       4. Blown snow particles concentration in saltation (kg/m3)     
!
        CALL CONC_SALT(ZCONC_SALT,ZWIND_FRC,ZWIND_FRC_THR,PRHOA,PQSALT,  &
                             PSNOWLIQ,PSNOWHIST)
!          
!       5. Turbulent flux of blown snow particles     
!                                
        CALL SNOW_FLUX(ZCONC_SALT,PVMOD,PBLOWSNW_CONC,ZCD_TEMP,PBLOWSNW_FLUX)
!    
!       Snowdrfit does not occur when : 
!          - surface layer is wet
!          - historical variable is higher than 1 : crust of non-transportable snow
!
        IF(PSNOWLIQ>0 .OR. PSNOWHIST>0) THEN
           PBLOWSNW_FLUX(:)=0
        ENDIF 

!       Update Richardson number and drag coefficient for momentum
        PRI     = ZRI_TOT
        PCDSNOW = ZCD_TEMP
       
       IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL:SNOWFLUX_GET_TURB',1,ZHOOK_HANDLE)

        END SUBROUTINE SNOWFLUX_GET_TURB 
       

       SUBROUTINE SNOWEROSION(PSWE_REM,PSNOWHEAT,PSNOWDZ,PSNOWSWE,       &
                                PSNOWRHO,PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,&
                                PSNOWAGE,PSNOWLIQ,INLVLS_USE,LAYER_USE,  &
                                PVMOD, PTSTEP)
!
USE MODE_SNOW3L
USE MODD_SNOW_PAR                                
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
USE MODE_BLOWSNW_MBLUTL
!
IMPLICIT NONE
!
!!    PURPOSE
!!    -------
!
!      Erosion of the snowpack surface layer :
!           reduces heat content, SWE and thickness
!
!      ASSUMPTION: the snowpack surface layer is never totally eroded.
!        Total erosion is insured earlier.
!                        
!       0.1 declarations of arguments        
!             
!       
        REAL, INTENT(IN) :: PSWE_REM   ! Snow mass to be removed
!        from the snow top layer (never exceed top layer thickness)
        REAL, INTENT(IN) :: PVMOD, PTSTEP
!
!       Top layers properties to be updated
!
        REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWHEAT,PSNOWDZ,PSNOWSWE,PSNOWRHO, &
                      PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE,PSNOWLIQ

        LOGICAL, DIMENSION(:), INTENT(INOUT) :: LAYER_USE 
        INTEGER,INTENT(INOUT)                :: INLVLS_USE

!       0.2 declaration of local variables
!        
        REAL :: ZSNOW_REM
        REAL ::ZSNOWHMASS,ZSNOWTEMP
        REAL :: ZRMOB,ZRDRIFT,ZPROFEQ,ZRT,ZSNOWRHO2
        REAL :: ZDRO, ZDGR1,ZDGR2
        REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSCAP      
        INTEGER :: JLAYER,KNLVLS  
        LOGICAL  :: OEVOL_GRAIN
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!       0.3 initialization
!         
IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL:SNOWEROSION',0,ZHOOK_HANDLE)
!
ZSNOW_REM=0.0 
ZPROFEQ=0.0
KNLVLS = SIZE(PSNOWSWE,1)
!
! Option to simulated evolution of snow properties in areas where 
! erosion occurs 
!
OEVOL_GRAIN = .FALSE.

                                    
!       1. Remove snow from surface layer

        IF(PSWE_REM==PSNOWSWE(1)) THEN ! Surface layer is totally removed.
               DO JLAYER=1,INLVLS_USE-1
                PSNOWSWE(JLAYER)=PSNOWSWE(JLAYER+1)
                PSNOWHEAT(JLAYER)=PSNOWHEAT(JLAYER+1)
                PSNOWDZ(JLAYER)=PSNOWDZ(JLAYER+1)
                PSNOWRHO(JLAYER)=PSNOWRHO(JLAYER+1)
                PSNOWLIQ(JLAYER)=PSNOWLIQ(JLAYER+1)
                PSNOWGRAN1(JLAYER)=PSNOWGRAN1(JLAYER+1)
                PSNOWGRAN2(JLAYER)=PSNOWGRAN2(JLAYER+1)
                PSNOWHIST(JLAYER)=PSNOWHIST(JLAYER+1)
                PSNOWAGE(JLAYER)=PSNOWAGE(JLAYER+1)
              ENDDO 
                PSNOWSWE(INLVLS_USE)=0.0
                PSNOWRHO(INLVLS_USE)=999.
                PSNOWDZ(INLVLS_USE)=0.
                PSNOWGRAN1(INLVLS_USE)=0.
                PSNOWGRAN2(INLVLS_USE)=0.
                PSNOWHIST(INLVLS_USE)=0.
                PSNOWAGE(INLVLS_USE)=0.
                PSNOWHEAT(INLVLS_USE)=0.
                PSNOWLIQ(INLVLS_USE)=0.
                LAYER_USE(INLVLS_USE)=.FALSE.           
                INLVLS_USE=INLVLS_USE-1 
        ELSE
                ZSCAP      = SNOW3LSCAP(PSNOWRHO)
                ZSNOWTEMP  = XTT + (PSNOWHEAT(1) +                &
                  XLMTT*PSNOWRHO(1)*PSNOWDZ(1))/                    &
                  (ZSCAP(1)*MAX(XSNOWDMIN/KNLVLS,PSNOWDZ(1)))
                ZSNOWTEMP  = MIN(XTT, ZSNOWTEMP)
                ZSNOWHMASS = PSWE_REM*(XCI*(ZSNOWTEMP-XTT)-XLMTT)
!
!       2.Reduce total pack depth                
                ZSNOW_REM=PSWE_REM/PSNOWRHO(1)
                PSNOWDZ(1)=PSNOWDZ(1)-ZSNOW_REM
!       3.Reduce heat content
                PSNOWHEAT(1)=PSNOWHEAT(1)-ZSNOWHMASS
                PSNOWSWE(1)=PSNOWDZ(1)*PSNOWRHO(1)
!       4.Evolution of the remaining layer to simulate subgrid density increase
!         and grain shape and size evolution due snowdrift.
                IF(OEVOL_GRAIN) THEN
                DO JLAYER=1,INLVLS_USE
                   ZSNOWRHO2  = PSNOWRHO(JLAYER)
                   ZRMOB = SNOW_MBL_INDEX(PSNOWGRAN1(JLAYER),PSNOWGRAN2(JLAYER))
                   IF(PSNOWHIST(JLAYER).ge.2) ZRMOB = MIN(ZRMOB,-0.0583)
                   ZRDRIFT = ZRMOB-(2.868*EXP(-0.085*PVMOD)-1.)
                   ZPROFEQ = ZPROFEQ + 0.5* PSNOWDZ(JLAYER) * (3.25-ZRDRIFT)
                   ZRT=MAX(0.,ZRDRIFT*EXP(-ZPROFEQ/0.1))

                   IF(PSNOWRHO(JLAYER).lt.350.) THEN
                      ZDRO=(350.-PSNOWRHO(JLAYER))*PTSTEP/12./3600.*ZRT
                      PSNOWRHO(JLAYER)=MIN(350.,PSNOWRHO(JLAYER)+ZDRO)
                      PSNOWDZ(JLAYER)=PSNOWDZ(JLAYER)*ZSNOWRHO2/   &
                            PSNOWRHO(JLAYER)
                   END IF
!                  Dendritic snow
                   IF (PSNOWGRAN1(JLAYER).lt.0.) THEN
                      ZDGR1=-PSNOWGRAN1(JLAYER)*PTSTEP/12./3600.*ZRT*0.5
                      ZDGR1=MIN(ZDGR1,-0.99*PSNOWGRAN1(JLAYER))
                      ZDGR2=ZRT*PTSTEP/12./3600.*(99.-PSNOWGRAN2(JLAYER))
                      PSNOWGRAN1(JLAYER)=PSNOWGRAN1(JLAYER)+ZDGR1
                      PSNOWGRAN2(JLAYER)=MIN(99.,PSNOWGRAN2(JLAYER)+zdgr2)                      
                  ELSE    ! Non dendritic snow    
                     ZDGR1=PTSTEP/12./3600.*ZRT*(99.-PSNOWGRAN1(JLAYER))
                     ZDGR2=PTSTEP/12./3600.*ZRT*5./10000.
                     PSNOWGRAN1(JLAYER)=MIN(99.,PSNOWGRAN1(JLAYER)+ZDGR1)
                     PSNOWGRAN2(JLAYER)=MAX(0.0003,PSNOWGRAN2(JLAYER)-ZDGR2)
                  END IF
                END DO
                END IF

            END IF

     IF (LHOOK) CALL DR_HOOK('SNOWPACK_EVOL:SNOWEROSION',1,ZHOOK_HANDLE)

     END SUBROUTINE SNOWEROSION

END SUBROUTINE SNOWPACK_EVOL
