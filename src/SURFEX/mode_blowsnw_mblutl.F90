!-----------------------------------------------------------------
!     ##################
      MODULE MODE_BLOWSNW_MBLUTL
!     ##################
!
!!****  *MODE_BLOWSNW_MBLUTL * - contains subroutines to determine:
!                  - threshold velocity for snow erosion as a function of snow
!                       pack properties,
!                  - blown snow particles concentration in the saltation layer,
!                  - surface turbulent flux of blown snow particles
!!
!!    PURPOSE
!!    -------
!
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
!!    AUTHOR
!!    ------
!!	V. Vionnet       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        21/03/08
!!      
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
!-------------------------------------------------------------------------------
CONTAINS


       FUNCTION SNOW_MBL_INDEX(PSNOWGRAN1, PSNOWGRAN2)                  &
RESULT(PMBLINDEX)
!
!!    PURPOSE
!!    -------      
!     Calculate the mobility index as a function of snow grain type 
!      Based on G&M,1998 : PROTEON, Centre d'Etude de la Neige 
!
USE MODD_SNOW_METAMO, ONLY : XGRAN,XEPSI
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
!*      0.1    Declarations of arguments
!
REAL, INTENT(IN)              :: PSNOWGRAN1, PSNOWGRAN2
!
REAL                          :: PMBLINDEX
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!       0.3    Initialization
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SNOW_MBL_INDEX',0,ZHOOK_HANDLE)
!
PMBLINDEX=0.
!
!      1.       Computation       
!  
IF(PSNOWGRAN1<XEPSI)THEN
        PMBLINDEX=-0.75*PSNOWGRAN1/XGRAN-0.5*PSNOWGRAN2/        &
                        XGRAN+0.5
ELSE IF(PSNOWGRAN1>=0 .AND. PSNOWGRAN2>0) THEN
        PMBLINDEX=-0.583*PSNOWGRAN2*1000-0.833*PSNOWGRAN1/       &
                        XGRAN+0.833
ELSE
        PMBLINDEX=0.                        
ENDIF

PMBLINDEX=MAX(MIN(PMBLINDEX,1.),-0.9999)

IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SNOW_MBL_INDEX',1,ZHOOK_HANDLE)
                                  
END FUNCTION SNOW_MBL_INDEX   
      
      
        SUBROUTINE WIND_RFR_THR(PWIND_RFR_THR,PSNOWGRAN1,               &
                                PSNOWGRAN2)
!
!!    PURPOSE
!!    -------      
!                                
!       Computation of threshold wind speed at level of reference (5
!       meter from G&M)
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
        IMPLICIT NONE
!
!       0.1 declarations of arguments        
!             
!                                                    
        REAL, INTENT(IN)      :: PSNOWGRAN1, PSNOWGRAN2
                                                                       
        REAL, INTENT(INOUT)   ::PWIND_RFR_THR                                                                 
!
!       0.2 declaration of local variables
!        
        REAL  :: ZMBLINDEX
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!       1. mobility index                                                
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:WIND_RFR_THR',0,ZHOOK_HANDLE)
!
       ZMBLINDEX= SNOW_MBL_INDEX(PSNOWGRAN1,PSNOWGRAN2)
!
!       2. threshold 5-m wind speed
!
       PWIND_RFR_THR=(log(2.868)-log(1+ZMBLINDEX))/0.085
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:WIND_RFR_THR',1,ZHOOK_HANDLE)
!
        END SUBROUTINE WIND_RFR_THR  

         SUBROUTINE WIND_FRC_THR(PWIND_RFR_THR,PZ0EFF,    &
                                PWIND_FRC_THR)
!
!!    PURPOSE
!!    -------      
!                                 
!       Compute threshold friction velocity from the threshold
!       velocity at a reference height of 5 m. We assume a log profile 
!       in the SBL with z0 the same as used to compute CDSNOW
!
!
        USE MODD_CSTS,ONLY : XKARMAN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
        IMPLICIT NONE
!
!       0.1 declarations of arguments        
!                                                                 
        REAL,  INTENT(IN)   :: PZ0EFF    
        REAL, INTENT(IN)    :: PWIND_RFR_THR
!
        REAL, INTENT(OUT)   :: PWIND_FRC_THR         
!
!       0.2 declaration of local variables
        REAL           :: ZCD
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

!            !
!      0.3     Initialization 
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:WIND_FRC_THR',0,ZHOOK_HANDLE)
!
       ZCD = (XKARMAN/LOG(5/PZ0EFF))**2

!       1. Threshold friction velocity
!
       PWIND_FRC_THR=SQRT(ZCD)*PWIND_RFR_THR   
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:WIND_FRC_THR',1,ZHOOK_HANDLE)
!
      END SUBROUTINE WIND_FRC_THR
      
      SUBROUTINE CONC_SALT(PCONC_SALT,PWIND_FRC_SALT,PWIND_FRC_THR,   &
                                        PRHOA,PQSALT,PSNOWLIQ,PSNOWHIST)
!
!!    PURPOSE
!!    -------      
!              
!        Saltation layer parameterization : computes
!           - snow particle concentration (kg/m3) 
!           - streamwise saltation flux (kg/m/s)
!            
! 
USE MODD_CSTS, ONLY : XG,XRHOLI  
USE MODD_BLOWSNW_SURF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

        IMPLICIT NONE
!
!       0.1 declarations of arguments        
!                                                                 
        REAL, INTENT(IN)      :: PWIND_FRC_SALT,PRHOA,    &
                                              PWIND_FRC_THR
                                    
        REAL, INTENT(IN)       :: PSNOWLIQ,PSNOWHIST
        REAL,  INTENT(INOUT)   :: PCONC_SALT  
        REAL,  INTENT(OUT)     :: PQSALT
!
!       0.2 declaration of local variables
!        
       REAL   :: ZSALT_FLUX
       REAL   :: ZHSALT
       REAL   :: ZUPART
       REAL(KIND=JPRB) :: ZHOOK_HANDLE

!       
!      0.3     Initialization             
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:CONC_SALT',0,ZHOOK_HANDLE)
!
!       1. Streamiwise saltation flux : horizontal saltation, transport rate under
!       equilibrium conditions

!            
       IF(PWIND_FRC_SALT>=PWIND_FRC_THR) THEN ! Transport in saltation
               IF(CSNOW_SALT=='POME' .OR. CSNOW_SALT=='TPOM') THEN   
!  
!               Formulation of Pomeroy and Gray, 1990 
!
                ZSALT_FLUX=0.68*PRHOA*PWIND_FRC_THR*              &
                        (PWIND_FRC_SALT**2-PWIND_FRC_THR**2)/     &
                        (PWIND_FRC_SALT*XG)

               ELSE IF(CSNOW_SALT=='SORE'.OR. CSNOW_SALT=='TSOR') THEN
!                       
!              Formulation of Sorensen (2004) (initially developped for sand) 
!              with coefficients adapted for snow based on the wind-tunnel
!              measurement of Nishimura and Hunt (2000). More details in the PhD
!              thesis of Vionnet (2012)
!
                ZSALT_FLUX=PRHOA/XG*PWIND_FRC_SALT**3.*               &
                         (1.-PWIND_FRC_THR**2./PWIND_FRC_SALT**2.)*   &
                         (2.6+2.*PWIND_FRC_THR/PWIND_FRC_SALT+2.5*    &         
                         PWIND_FRC_THR**2./PWIND_FRC_SALT**2.) 
               ENDIF 
       ELSE ! No transport
         ZSALT_FLUX=0.
       ENDIF  
!    
!       Saltation does not occur when : 
!          - surface layer is wet
!          - historical variable is higher than 1 : crust of non-transportable snow
        IF(PSNOWLIQ>0 .OR. PSNOWHIST>0) THEN 
           ZSALT_FLUX=0
        ENDIF    
!    
!       Store streamise saltation flux for future advection
        PQSALT=ZSALT_FLUX   
!
!       2. Saltation height (Pomeroy and Male, 1992)                        
!   
         ZHSALT = 0.08436*PWIND_FRC_SALT**1.27
!
!       3. Horizontalparticle velocity Up=2,8.U*t 
!          according to Pomeroy and Gray (1900) 
!
        ZUPART = 2.8*PWIND_FRC_THR
!
!       4. Saltating particle concentration  
!
        IF(ZSALT_FLUX>0.) THEN
          PCONC_SALT=ZSALT_FLUX/(ZUPART*ZHSALT)
        ELSE
          PCONC_SALT=0.
        ENDIF
            
        IF(CSNOW_SALT=='TPOM' .OR. CSNOW_SALT=='TSOR') THEN
!        Compute concentration at the top of the saltation layer assuming 
!        an exponential flux profile in the saltation layer (cf Doorschot
!        et al., 2004 and Nishimura and Hunt, 2000).
             IF(ZSALT_FLUX>0.) THEN        
                PCONC_SALT = ZSALT_FLUX*0.45*XG/PWIND_FRC_SALT**2    &
                            *EXP(-0.45*ZHSALT*XG/PWIND_FRC_SALT**2)  &
                             /(ZUPART)
             ENDIF                    
        ENDIF

        IF(LSNOW_SALT) THEN  
!       Imposed concentration in the saltation layer (only for sensitivity
!       studies)
          PCONC_SALT= XCONC_SALT
          PQSALT    = XCONC_SALT*ZUPART*ZHSALT
        END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:CONC_SALT',1,ZHOOK_HANDLE)
!
      END SUBROUTINE CONC_SALT


        SUBROUTINE SNOW_FLUX(PCONC_SALT,PVMOD,PCONC_SURF,               &
                                PCDSNOW,PSNOW_FLUX)
!       Compute snow turbulent vertical flux as in Gallee et al (2001)
!       Flux=u*q*=Cd.V.(Q(Surf)-Q(Salt))    (kg/(m2.s))
!       Cd : drag coefficient of momentum   
!       V wind speed in the lowest level of the atmo. model (m/s)
!       Q(Surf) blown snow concentration in the lowest level of the atmo. model (kg/m3)
!       Q(Salt) blown snow concentration in the saltation layer (kg/m3)
!       Only positive flux : from the saltation layer to the model lowest level. Deposition is computed later
!
!       NOTE: the lowest level of the atmo. model can be Canopy or directly Meso-NH if Canopy is
!             not activated
!
USE MODD_CSTS, ONLY : XG,XRHOLI, XPI 
USE MODD_BLOWSNW_SURF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

        IMPLICIT NONE
!
!       0.1 declarations of arguments        
!
        REAL, INTENT(IN)      :: PCONC_SALT,PCDSNOW,PVMOD
        REAL, DIMENSION(:),INTENT(IN)   :: PCONC_SURF                                                 
!       
        REAL, DIMENSION(:),INTENT(OUT)   :: PSNOW_FLUX         
!
!
!       0.2 declaration of local variables
!        
       REAL   :: ZNUM_SALT ! Number particle concentration in the saltation
!                          layer (#.m^{-3}) 
       REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!       1. Initialisation following the gamma distribution
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SNOW_FLUX',0,ZHOOK_HANDLE)
!
        ZNUM_SALT = 3*PCONC_SALT/(4*XPI*XRHOLI*(XEMIALPHA_SNW+2.)*(XEMIALPHA_SNW+1.)*   &
                XEMIALPHA_SNW*(XEMIRADIUS_SNW/XEMIALPHA_SNW)**3)
!                
!       2. Turbulent fluxes         
!
!       Number flux
!
        PSNOW_FLUX(1)=MAX(0.,-PCDSNOW*PVMOD*(PCONC_SURF(1)-           &
                                ZNUM_SALT))
!       Mass flux                                
        PSNOW_FLUX(2)=MAX(0.,-PCDSNOW*PVMOD*(PCONC_SURF(2)-           &
                                PCONC_SALT))
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SNOW_FLUX',1,ZHOOK_HANDLE)
!      
END SUBROUTINE SNOW_FLUX

      
    SUBROUTINE SURFACE_RI_PART(PCONC_SALT, PCONC_SURF,PRHOA,          &
                             PZREF, PUREF, PDIRCOSZW, PVMOD, PRI_PART )
!   ######################################################################
!
!!****  *SURFACE_RI*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the component of richardson number near the ground associated
!     with particle-induced buoyancy.
!         
!     
!!**  METHOD
!!    ------
!
!      Based on surface_ri.f90 
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!    MODD_GROUND_PAR
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Vionnet           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/10/10 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XRV, XRD, XG
USE MODI_WIND_THRESHOLD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, INTENT(IN)    :: PCONC_SALT ! Concentration of blowing snow particle
!                                                 in the saltation layer (kg/m3) 
REAL,   INTENT(IN)    :: PCONC_SURF ! Concentration of blowing snow particle
!                                                 at the 1st atmospheric level (kg/m3) 
REAL,   INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL,   INTENT(IN)    :: PRHOA    ! air density (kg/m3)
!
REAL,   INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL,   INTENT(IN)    :: PUREF    ! reference height of the wind
!                                             ! NOTE this is different from ZZREF
!                                             ! ONLY in stand-alone/forced mode,
!                                             ! NOT when coupled to a model (MesoNH)
REAL,   INTENT(IN)    :: PDIRCOSZW! Cosine of the angle between
!                                             ! the normal to the surface and
!                                             ! the vertical
!
REAL,   INTENT(OUT)   :: PRI_PART ! "Particle" Richardson number
!
!*      0.2    declarations of local variables
!
REAL                  :: ZVMOD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       1.     Richardson number
!              -----------------
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SURFACE_RI_PART',0,ZHOOK_HANDLE)
!
ZVMOD = PVMOD
!
!                                                 
                                                ! Richardson's number
PRI_PART = XG * PDIRCOSZW/PRHOA * PUREF * PUREF              &
        * (PCONC_SALT-PCONC_SURF)                            &
        / (ZVMOD*ZVMOD) /PZREF
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SURFACE_RI_PART',1,ZHOOK_HANDLE)
!
END SUBROUTINE SURFACE_RI_PART

!   #################################################################
    SUBROUTINE SURFACE_CD_PART(PRI, PZREF, PUREF, PZ0EFF, PZ0H,   &
                             PCD, PCDN)
!   #################################################################
!
!!****  *SURFACE_CD_PART*
!!
!!    PURPOSE
!!    -------
!
!     Computes the drag coefficients for momentum near the ground
!         including particle-induced buoyancy.
!     
!!**  METHOD
!!    ------
!
!    1 and 2 : computation of relative humidity near the ground
!
!    3 : richardson number including particle-induced buoyancy
!
!    4 : the aerodynamical resistance for heat transfers is deduced
!
!    5 : the drag coefficient for momentum ZCD is computed
!           including particle-induced buoyancy
!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Vionnet           * Meteo-France *
!!    Based on surface_cd.f90 by V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/10/10

!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XKARMAN
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL,  INTENT(IN)    :: PRI      ! Richardson number
REAL,  INTENT(IN)    :: PZREF    ! reference height of the first
                                             ! atmospheric level
REAL,  INTENT(IN)    :: PUREF    ! reference height of the wind
!                                             ! NOTE this is different from ZZREF
!                                             ! ONLY in stand-alone/forced mode,
!                                             ! NOT when coupled to a model (MesoNH)
REAL,  INTENT(IN)    :: PZ0EFF   ! roughness length for momentum
                                              ! with subgrid-scale orography
REAL,  INTENT(IN)    :: PZ0H     ! roughness length for heat
!
REAL,  INTENT(OUT)   :: PCD      ! drag coefficient for momentum
REAL,  INTENT(OUT)   :: PCDN     ! neutral drag coefficient for momentum
!
!*      0.2    declarations of local variables
!
!
REAL                       :: ZZ0EFF, ZZ0H, ZMU,     &
                              ZCMSTAR, ZPM, ZCM, ZFM
INTEGER                    :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------
!
!*       1.     Drag coefficient for momentum transfers
!               ---------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SURFACE_CD_PART',0,ZHOOK_HANDLE)
!
  ZZ0EFF = MIN(PZ0EFF,PUREF*0.5)
  ZZ0H   = MIN(ZZ0EFF,PZ0H)
!
  ZMU = LOG( MIN(ZZ0EFF/ZZ0H,200.) )
!
  PCDN = (XKARMAN/LOG(PUREF/ZZ0EFF))**2

  ZCMSTAR = CMSTAR(ZMU)
  ZPM     = PM(ZMU)
!
  ZCM = 10.*ZCMSTAR*PCDN*( PUREF/ZZ0EFF )**ZPM
!
  IF ( PRI > 0.0 ) THEN
    ZFM = 1. + 10.*PRI / SQRT( 1.+5.*PRI )
    ZFM = 1. / ZFM
  ELSE
    ZFM = 1. - 10.*PRI / ( 1.+ZCM*SQRT(-PRI) )
  ENDIF
!
  PCD = PCDN*ZFM
!
!
!-------------------------------------------------------------------------------
!
!
!
IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SURFACE_CD_PART',1,ZHOOK_HANDLE)

!
CONTAINS
!
!                                              functions used in the calculation
!                                              the terms Cm
!  
  FUNCTION CMSTAR(X)
  REAL, INTENT(IN)     :: X
  REAL                :: CMSTAR
  !
  CMSTAR = 6.8741 + 2.6933*X - 0.3601*X*X + 0.0154*X*X*X
  !
  END FUNCTION CMSTAR
  !
  !
  FUNCTION PM(X)
  REAL, INTENT(IN)     :: X
  REAL                 :: PM
  !
  PM = 0.5233 - 0.0815*X + 0.0135*X*X - 0.0010*X*X*X
  !
  END FUNCTION PM
                               
!
!-------------------------------------------------------------------------------


!
END SUBROUTINE SURFACE_CD_PART

SUBROUTINE SOLVE_CD(PUSTAR,PVMOD,PRHOA,PCONC_SURF,PZREF,PUREF,PDIRCOSZW,   &
                   PWIND_FRC_THR,PRI,PZ0EFF,HSNOWRES,PZ0H,PRES,PRI_TOT,    & 
                   PCDSNOW,PSNOWLIQ,PSNOWHIST)

!   
!    Compute the value of the function used in the Newton's iterative algorithm. 
!


USE MODD_BLOWSNW_SURF
USE MODD_SNOW_PAR,    ONLY : X_RI_MAX
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,     INTENT(IN)    :: PUSTAR
REAL,     INTENT(IN)    :: PSNOWLIQ,PSNOWHIST
REAL,     INTENT(IN)    :: PVMOD,PRHOA,PCONC_SURF,PZREF,PUREF
REAL,     INTENT(IN)    :: PDIRCOSZW,PWIND_FRC_THR,PRI,PZ0EFF,PZ0H
REAL,     INTENT(OUT)   :: PCDSNOW, PRI_TOT,PRES
CHARACTER(LEN=*), INTENT(IN)    :: HSNOWRES 
!
!*      0.2    declarations of local variables
!
REAL :: ZCONC_SALT,ZRI_PART,ZCD, ZCDN, ZZ0_TEMP,ZQSALT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

PRI_TOT = PRI
ZZ0_TEMP = PZ0EFF 

IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SOLVE_CD',0,ZHOOK_HANDLE)

IF(XSNOW_BUOYANCY==1) THEN
!                        
!       Blown snow particles concentration (kg/m3) in the saltation
!       layer    
!                 
CALL CONC_SALT(ZCONC_SALT,PUSTAR,PWIND_FRC_THR,PRHOA,ZQSALT  ,&
               PSNOWLIQ,PSNOWHIST)

CALL SURFACE_RI_PART(ZCONC_SALT,PCONC_SURF, PRHOA,      &
        PZREF, PUREF, PDIRCOSZW, PVMOD, ZRI_PART        )

!                
PRI_TOT = PRI + ZRI_PART 
!
! Simple adaptation of method by Martin and Lejeune (1998)
! to apply a lower limit to turbulent transfer coef
! by defining a maximum Richarson number for stable
! conditions:
!
IF(HSNOWRES=='RIL') PRI_TOT = MIN(X_RI_MAX, PRI_TOT)
END IF

IF(XSNOW_ROUGHNESS==1 .AND. PUSTAR > PWIND_FRC_THR  ) THEN
      ! Increase in surface rougnhess due to saltation following
      !  Dover (1993) z0_s = z0_ns*(u*/uth*)�    
      ! Limit increase to value observed over snow surface      
      ZZ0_TEMP = MIN(PZ0EFF * (PUSTAR/PWIND_FRC_THR)**2,0.01)
END IF

CALL SURFACE_CD_PART(PRI_TOT, PZREF, PUREF, ZZ0_TEMP, PZ0H, PCDSNOW, ZCDN)

PRES = PUSTAR - SQRT(PCDSNOW)* PVMOD

IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_MBLUTL:SOLVE_CD',1,ZHOOK_HANDLE)


END SUBROUTINE SOLVE_CD

END MODULE MODE_BLOWSNW_MBLUTL
