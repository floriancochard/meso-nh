!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.

!-----------------------------------------------------------------
!!   ########################
     MODULE MODI_ABDULRAZZAK
!!   ########################
!!
INTERFACE
!
SUBROUTINE ABDULRAZZAK(KAP, PRHOP, PNUE, PEPS, PMI, PRHODREF, PTEMP, PW, PTRAD, &
                        PABST, PRG, PSIG, PN0, PCTOTA, PNCN, PMCN, PSMAX, PT0, PFRACT_ACT1, &
                        PFRACT_ACT2, PKELVIN1, PKELVIN2, PRAOULT1,PRAOULT2, PR0AP1, PR0AP2)

IMPLICIT NONE

INTEGER, INTENT(IN)                   :: KAP 
REAL,                   INTENT(IN)    :: PRHOP
REAL,  DIMENSION(:,:),INTENT(IN)    :: PNUE
REAL,  DIMENSION(:,:,:),  INTENT(IN)    :: PEPS
REAL,  DIMENSION(:,:,:),  INTENT(IN)    :: PMI
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PN0
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PRG
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PSIG
REAL,  DIMENSION(:,:,:),INTENT(IN)    :: PCTOTA
REAL,  DIMENSION(:),    INTENT(IN)    :: PTEMP
REAL,  DIMENSION(:),    INTENT(IN)    :: PABST
REAL,  DIMENSION(:),    INTENT(IN)    :: PW
REAL,  DIMENSION(:),    INTENT(IN)    :: PTRAD
REAL,  DIMENSION(:),    INTENT(IN)    :: PRHODREF
REAL,  DIMENSION(:,:),  INTENT(OUT)   :: PNCN    ! OUT : Number of actived aerosols (#/m3)
REAL,  DIMENSION(:,:),  INTENT(OUT)   :: PMCN    ! OUT : Mass of actived aerosols (ug/m3)
REAL,  DIMENSION(:),  INTENT(OUT)   :: PSMAX   ! OUT : Maximum supersaturation
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PT0     ! OUT : Surface tension
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)     :: PFRACT_ACT1  !OUT : Activated fraction mode 1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)     :: PFRACT_ACT2  !OUT : Activated fraction mode 2
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PKELVIN1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PKELVIN2
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PRAOULT1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PRAOULT2
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PR0AP1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PR0AP2

END SUBROUTINE ABDULRAZZAK
!!
END INTERFACE
!!
END MODULE MODI_ABDULRAZZAK
!!
!     ######spl
SUBROUTINE ABDULRAZZAK( KAP, PRHOP, PNUE, PEPS, PMI, PRHODREF, PTEMP, PW, PTRAD,    &
                        PABST, PRG, PSIG, PN0, PCTOTA, PNCN, PMCN, PSMAX, PT0, PFRACT_ACT1, &
                        PFRACT_ACT2, PKELVIN1, PKELVIN2, PRAOULT1,PRAOULT2, PR0AP1, PR0AP2  )
!!   #######################################
!!
!   PURPOSE
!!   -------
!!   Activation of aerosols from Abdul-Razzak parameterization.
!!   Source of activation came from updrafts velocities and/or 
!!   cooling radiative rate.
!!   Aerosols size distribution, number of modes and chemical composition
!!   (number of inorganic and organic salt, surfactants) came from ORILAM
!!   aerosol scheme (Tulet et al., 2005 and 2006 of J. Geo. Res.), and 
!!   computed in MesoNH routine ch_aer_activation.
!!   Raoult term is computed from inorganic and organic salt concentration.
!!   Kelvin term is computed from organic surfactant concentration that changes
!!   the surface tension. Micelle formation can limit the surfactant concentration. 
!!   

!!
!!   REFERENCE
!!   ---------
!!   Abdul-Razzak, 1998, 2000 and 2004, J. Geo. Res. (Parameterization)
!!   Li et al., 1998, J. Atm. Sci. (Micelle formation)
!!
!!   AUTHOR
!!    ------
!!     J. Rangognio (CNRM/GMME) &  P. Tulet (CNRM/GMEI) 
!!     R.L. Curier (LAMP) & M. Leriche (LAMP) 
!!
!!   MODIFICATIONS
!!    -------------
!!
!!   IMPLICIT ARGUMENTS

USE MODD_CH_AEROSOL

IMPLICIT NONE

!! Argument variables
INTEGER, INTENT(IN)                   :: KAP     ! Number of aerosol modes
REAL,                   INTENT(IN)    :: PRHOP   ! Aerosol density (kg/m3)
REAL,  DIMENSION(:,:),INTENT(IN)    :: PNUE    ! Number of dissociative ions
REAL,  DIMENSION(:,:,:),  INTENT(IN)    :: PEPS    ! Soluble fraction of aerosol
REAL,  DIMENSION(:,:,:),  INTENT(IN)    :: PMI     ! Molecular weight (g/mol)
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PN0     ! Aerosol number concentration (#/m3)
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PRG     ! Aerosol mean radius (um)
REAL,  DIMENSION(:,:),  INTENT(IN)    :: PSIG    ! Aerosol standard deviation
REAL,  DIMENSION(:,:,:),INTENT(IN)    :: PCTOTA  ! Aerosol mass (ug/m3)
REAL,  DIMENSION(:),    INTENT(IN)    :: PTEMP   ! Air temperature (K)
REAL,  DIMENSION(:),    INTENT(IN)    :: PABST   ! Air pressure (Pa)
REAL,  DIMENSION(:),    INTENT(IN)    :: PW      ! Activation vertical velocity (m/s)
REAL,  DIMENSION(:),    INTENT(IN)    :: PTRAD   ! Activation cooling radiative tendency (K/s)
REAL,  DIMENSION(:),    INTENT(IN)    :: PRHODREF ! Air density (kg/m3)
REAL,  DIMENSION(:,:),  INTENT(OUT)   :: PNCN    ! OUT : Number of actived aerosols (#/m3)
REAL,  DIMENSION(:,:),  INTENT(OUT)   :: PMCN    ! OUT : Mass of actived aerosols (ug/m3)
REAL,  DIMENSION(:),    INTENT(OUT)   :: PSMAX   ! OUT : Maximum supersaturation
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PT0     ! OUT : Surface tension
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)     :: PFRACT_ACT1  !OUT : Activated fraction mode 1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)     :: PFRACT_ACT2  !OUT : Activated fraction mode 2
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PKELVIN1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PKELVIN2
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PRAOULT1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PRAOULT2
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PR0AP1
REAL,  DIMENSION(:), OPTIONAL, INTENT(OUT)   :: PR0AP2

!! Local Variables

!parameters for lognormal distributions of ap
REAL,  DIMENSION(SIZE(PCTOTA,1),SIZE(PCTOTA,2),SIZE(PCTOTA,3)) :: ZCCTOT
REAL,  DIMENSION(SIZE(PCTOTA,1))  :: ZSUM
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: ZRCI ! radius of the smalest activated particle (in um)
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: R0AP ! mean geometric radius for AP spectrum (m)
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: XNAP            ! amplitude of ap spectrum (m-3)
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: VAP             ! mixing ratio in volume (m3/m3 of air)

!physico-chemical properties of AP
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: rhoap           !volumetric mass of ap (g/cm3)
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: xms             !ap molecular weight (g/mole)

!parameters for resolution of activation using Abdul-Razzak param.
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: eta,f,g1        !fi, gi, eta in Smax formula
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: u,um            !factor u in erf
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: B               !Raoult term total
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: SM              !supersaturation for each ap mode
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: ratio           !factor in Smax formula for each mode
REAL,  DIMENSION(SIZE(PN0,1))     :: sum1, sum2, sum3, sum4
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: xmap            !ap mass, 
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: vapac,xmapac    !activated volume, activated mass
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: xD0AP           !use for new diameter calculation
REAL,  DIMENSION(SIZE(PN0,1),KAP) :: XNNUCS          !number of activated ap for each ap mode (m-3)
REAL, DIMENSION(SIZE(PN0,1))   ::  ESW               ! saturated vapor pressure in Pa
REAL, DIMENSION(SIZE(PN0,1))   ::  EDIFU             ! diffusivity of water vapor
REAL, DIMENSION(SIZE(PN0,1))   ::  EK                ! thermal conductivity of air
REAL, DIMENSION(SIZE(PN0,1))   ::  COND1,COND2,COND  ! term used to calculate 1/A3
REAL, DIMENSION(SIZE(PN0,1))   ::  SIGMAWV           ! water surface tension (kg/s2)
REAL, DIMENSION(SIZE(PN0,1))   ::  A0                ! Kelvin term for pure water
REAL, DIMENSION(SIZE(PN0,1),KAP)   ::  A             ! Final Kelvin term
REAL, DIMENSION(SIZE(PN0,1))   ::  ZOMEGAP           ! Abdul-Razzak parameter
REAL, DIMENSION(SIZE(PN0,1))   ::  ZOMEGAG           ! Abdul-Razzak parameter
REAL, DIMENSION(SIZE(PN0,1))   ::  ZX                ! Abdul-Razzak parameter
REAL, DIMENSION(SIZE(PN0,1))   ::  ZBETA             ! Abdul-Razzak parameter
REAL, DIMENSION(SIZE(PN0,1))   ::  ZQUI              ! Abdul-Razzak parameter
REAL, DIMENSION(SIZE(PN0,1))   ::  ZPHI              ! Abdul-Razzak parameter
REAL, DIMENSION(SIZE(PN0,1))   ::  ZTO               ! Surface tension (dyn/m2)
REAL, DIMENSION(SIZE(PN0,1))   ::  ZR, ZR1           ! Critical radius (um)
REAL, DIMENSION(SIZE(PN0,1))   ::  ZCMC              ! Limit of surfactant concentration
                                                     ! (micelle formation)
REAL, DIMENSION(SIZE(PN0,1))   ::  ZTM               ! Maximum surface excess concentration
REAL, DIMENSION(SIZE(PN0,1))   ::  ZCSAL             ! Salt concentration (ug/m3)
REAL, DIMENSION(SIZE(PN0,1))   ::  ZCSFTM            ! Surfactant concentration (ug/m3)
REAL, DIMENSION(SIZE(PN0,1))   ::  ZFRACSEL          ! Fraction de sels dissous / masse seche aerosols 
REAL, DIMENSION(SIZE(PN0,1),KAP) ::  ZEPSM           ! Surfactant mass fraction 
REAL, DIMENSION(SIZE(PN0,1),KAP) ::  ZMSFT           ! Surfactant mean molecular weight
REAL, DIMENSION(SIZE(PN0,1),KAP) ::  ZMSAL           ! Salt mean molecular weight
REAL, DIMENSION(SIZE(PN0,1))   ::  alpha, gama       ! terms in prognistic equation of supersaturation
REAL, DIMENSION(SIZE(PN0,1))   ::  cc0, dzeta        ! factor used in eta, dzeta for Smax calculation
REAL, DIMENSION(SIZE(PN0,1))   ::  somratio          ! sum of ratio for Smax calculation
REAL, DIMENSION(SIZE(PN0,1))   ::  XNUCNC            ! total number of activated ap (m-3)
REAL, DIMENSION(SIZE(PN0,1))   ::  XNNUC             ! total number of activated ap (m-3)
REAL, DIMENSION(SIZE(PN0,1))   ::  CHLAV             ! vaporisation latent heat calculation

! thermodynamics and physics constants  
REAL   ::  ALAT1,ALAT2 ! vaporisation latent heat calculation
REAL   ::  PI        ! pi value
REAL   ::  XMW, XMA  ! water and air molecular weight
REAL   ::  CP        ! calorific capacity at constant pressure
REAL   ::  TZ0       ! temperature 0°C in Kelvin
REAL   ::  RR        ! universal gas constant
REAL   ::  RD        ! gas constant of dry air
REAL   ::  RV        ! gas constant of water vapor
REAL   ::  GAPS      ! acceleration of gravity
REAL   ::  RHOW      ! volumetric mass of water

! dynamics variables
REAL, DIMENSION(SIZE(PN0,1))   ::  TEMP      ! temperature (K)
REAL, DIMENSION(SIZE(PN0,1))   ::  P         ! pressure (hPa)
REAL, DIMENSION(SIZE(PN0,1))   ::  w         ! vertical velocity (m/s)	
REAL, DIMENSION(SIZE(PN0,1))   ::  PSIRAD    ! Latent condensation heat 


! Local index
INTEGER ::  i, JSV, JJ, JK ! index    

!  0. Initialization of variables
!     ===========================
TEMP(:)= PTEMP(:)
w(:)   = PW(:)
P(:)   = PABST(:)*1E-2  ! pression in hPa

!constants and initialization
PI = 3.141592654
XMW = 18.016 !water molecular weight
XMA = 28.9644 !air molecular weight
CP = 287.04*3.5
TZ0 = 273.15
RR= 8.3144 
RD = 287.04
RV = 461.6 ! gas constant of water vapor
GAPS = 9.807
RHOW = 1000. !in kg/m3
XNNUC = 0.

! Aerosol number concentration 
XNAP(:,:) = PN0(:,:) ! particles/m3 air
! Aerosol mean radius  (of lognormal fonction) 
R0AP(:,:) = PRG(:,:)*1.E-6 ! in meter
! Volume particle in m3:
VAP(:,:) = 4.*PI/3. *  PN0(:,:) *&
               (R0AP(:,:)**3)*EXP(4.5 * LOG(PSIG(:,:))**2)


!calculate parameters for nucleation
ESW(:) = 610.7*exp(17.15*(TEMP(:)-TZ0)/(TEMP(:)-38.25)) !TEMP=temperature 
EDIFU(:) = 0.211*1.e-4*1013.25*((TEMP(:)/TZ0)**1.94)/P(:) !P=pression in hPa
COND1(:) = RV*TEMP(:)/(ESW(:)*EDIFU(:))
EK(:) = 4.18684E-3*(5.69+0.017*(TEMP(:)-TZ0))
ALAT1 = 4186.84*(597.3+TZ0*0.566)/CP
ALAT2 = 4186.84*0.566/CP
CHLAV(:) = CP*(ALAT1-ALAT2*TEMP(:))
COND2(:) = CHLAV*(CHLAV(:)/(RV*TEMP(:))-1.)/(EK(:)*TEMP(:))
COND(:) = RHOW*(COND1(:)+COND2(:)) ! equal to term 1/A3 in Pruppacher
PSIRAD(:) = -1*XMW*CHLAV(:)/(XMA*RD*(TEMP(:)**2)) 


!  1. Raoult term computation from Abdul-Razzak et al. 1998 and 2000
!     ==============================================================

somratio(:) = 0.
! Loop on every aerosol modes
do i = 1, kap 
  ratio(:,i) = 0.
    sum1(:) = 0. 
    sum2(:) = 0. 
    sum3(:) = 0. 
    sum4(:) = 0. 
    ZMSFT (:,i) = 0.
    ZMSAL (:,i) = 0.
    ZEPSM (:,i) = 0.

    DO JSV=1,NSP+NCARB+NSOA
      IF (JSV /= JP_AER_H2O) THEN
      ! sum1 en moles.kg-1 air
      sum1(:) = sum1(:) + PCTOTA(:,JSV,i)*1E-9 / &  ! kg.m-3 air
                PRHODREF(:)*PNUE(JSV,i)* &        ! kg.kg-1 air
                PEPS(:,JSV,i) / (PMI(:,JSV,i)*1E-3)     ! moles.kg-1 air  
                                  
      sum2(:) = sum2(:) + PCTOTA(:,JSV,i) * 1E-9 / PRHODREF(:) / PRHOP
      
      sum3(:) = sum3(:) + PCTOTA(:,JSV,i) * 1E-9    ! en kg.m-3  (aerosol dry mass)
      sum4(:) = sum4(:) + PEPS(:,JSV,i) * PCTOTA(:,JSV,i) * 1E-9 ! kg/m3 (aerosol acqueous mass)
      END IF
    ENDDO 
   
! Final  Raoult term  
    B(:,i) = XMW*1.E-3*sum1(:)/ (RHOW*sum2(:)) 

!  2. Compute Kelvin term of Abdul-Razzak et al. 2004 parameterization
!     ================================================================

! Aerosol salt mass fraction 
ZFRACSEL(:) = sum4(:) / sum3(:)

! Salt inorganic mean molecular weight (kg.mol-1)
   sum3(:) = 0.
   DO JSV= 1,NSP
       IF ((JSV /= JP_AER_H2O)) THEN
                ZMSAL(:,i) = ZMSAL(:,i) + PMI(:,JSV,i) * 1E-3 *&  
                             PEPS(:,JSV,i) * PCTOTA(:,JSV,i)  !!en kg.mol-1
        sum3(:) = sum3(:) + PEPS(:,JSV,i) * PCTOTA(:,JSV,i) 
       ENDIF
   ENDDO
   WHERE (sum3(:) .GT. 1E-10)
     ZMSAL(:,i) = ZMSAL(:,i) / sum3(:)
   ELSEWHERE
     ZMSAL(:,i) = 64.E-3
   ENDWHERE
   
! Fraction of surfactant in the soluble aerosol mass:  
! surfactant mass over the total acqueous mass 
   DO JSV=NSP,NSP+NCARB+NSOA
         ZEPSM(:,i) = ZEPSM(:,i) + PEPS(:,JSV,i) * PCTOTA(:,JSV,i) * 1E-9
   ENDDO
   WHERE (sum3(:) .GT. 0.)
     ZEPSM(:,i) = ZEPSM(:,i) / sum4(:)
   ELSEWHERE
     ZEPSM(:,i) = 0.
   ENDWHERE
     

! Surfactant mean molecular weight (kg.mol-1)
   sum3(:) = 0.
   DO JSV=NSP+1,NSP+NCARB+NSOA
        ZMSFT(:,i) = ZMSFT(:,i) + PMI(:,JSV,i) * 1E-3 * &
                     PEPS(:,JSV,i) * PCTOTA(:,JSV,i) 
        sum3(:) = sum3(:) + PEPS(:,JSV,i) * PCTOTA(:,JSV,i)
   ENDDO
   WHERE (sum3(:) .GT. 1E-10)
     ZMSFT(:,i) = ZMSFT(:,i) / sum3(:)
   ELSEWHERE
     ZMSFT(:,i) = 64.E-3
   ENDWHERE

! Surface tension of pure water
SIGMAWV(:) = 1.e-3*(76.1-(0.155*(TEMP(:)-TZ0))) ! N.m-1 ou kg/s2

! Kelvin term of pure water
A0(:) = 2.*XMW*1.e-3*SIGMAWV(:)/(RR*RHOW*TEMP(:))   

ZR(:) = (3.*B(:,i)*(R0AP(:,i))**3/A0(:))**(1./2.)   

! initial condition of the iteration loop
ZQUI(:) = 1.
ZOMEGAG(:) = 0.
ZOMEGAP(:) = 2.
ZR(:) = MAX(ZR(:), 1E-10)
ZR1(:) = ZR(:)

! Iteration loop 
do JJ=1,10
ZR(:) = SQRT(ZR(:)*ZR1(:)) 

! Salt concentration (moles per liter)
ZCSAL(:) = (PRHOP * (1. - ZEPSM(:,i)) * ZFRACSEL(:) * &
           (R0AP(:,i) / ZR(:))**3. / ZMSAL(:,i)) ! mole.m-3
ZCSAL(:) = ZCSAL(:) * 1E-3 !moles.liter
ZCSAL(:) = MAX(ZCSAL(:), 1E-10)

! Surfactant concentration (moles per liter)
ZCSFTM(:) = PRHOP * ZEPSM(:,i) * ZFRACSEL(:) * & 
            (R0AP(:,i)/ZR(:))**3. / ZMSFT(:,i)
ZCSFTM(:) = ZCSFTM(:) * 1E-3
   
! Limitation of Surfactant concentration by micelle formation 
WHERE (ZCSFTM(:) .GE. 0.01)  
ZCMC(:) = EXP(-0.5389*LOG(ZCSAL(:)) - 3.3366)
ELSEWHERE
ZCMC(:) = EXP(0.008 - 0.26*ZCSAL(:)) 
ENDWHERE
WHERE (ZCSFTM(:) .GT. ZCMC(:))
 ZCSFTM(:) = ZCMC(:) 
ENDWHERE
ZCSFTM(:) = MAX(ZCSFTM(:), 1.E-20)

! Maximum excess surface concentration (moles.dm-2) (* T * R)
ZTM(:) = 0.
WHERE (ZCSAL(:) .LT. 0.001)
  ZTM(:) = (- 2288.28 * ZCSAL(:) + 18.9907) / (RR * TEMP(:))  
ENDWHERE
WHERE ((ZCSAL(:) .GE. 0.001).AND.(ZCSAL(:) .LT. 0.01))
  ZTM(:) = (- 697.651 * ZCSAL(:) + 18.237) / (RR * TEMP(:))   
ENDWHERE
WHERE ((ZCSAL(:) .GE. 0.01).AND.(ZCSAL(:) .LT. 0.2))  
  ZTM(:) = (11.4271 - 9.2084 * ZCSAL(:) + 94.2226 * ZCSAL(:)**2.) / (RR * TEMP(:))  
ENDWHERE
WHERE (ZCSAL(:) .GE. 0.2)
  ZTM(:) = (1.2859 * ZCSAL(:) + 13.097) / (RR * TEMP(:))   
ENDWHERE

!Szyskowski constant
ZBETA(:) = 0.
WHERE (ZCSAL(:) .LT. 0.001)
ZBETA(:) = -0.365 * ZCSAL(:) + 0.001658   
ELSEWHERE
ZBETA(:) = 1.2E-5 * ZCSAL(:)**(-0.6808)  
ENDWHERE

! Final surface tension ZTO
ZTO(:) =  1E3*SIGMAWV(:) - &    !g/s2
          RR * TEMP(:) * ZTM(:)*  LOG(1+ZCSFTM(:) / (ZBETA(:)*(1+ZOMEGAG(:)))) 
ZTO(:) =  MAX(ZTO(:), 10.)

ZPHI(:) = (3.*ZTM(:))/(2.*ZOMEGAP(:)*ZR(:)*1E6*ZBETA(:)) &  
          -(0.5*(1.+ZCSFTM(:)/ZBETA(:)))

ZOMEGAG(:) = ZPHI(:) + (ZPHI(:)**2. + & 
             (3.*ZTM(:))/(ZOMEGAP(:)*ZR(:)*1E6*ZBETA(:)*1E3))**0.5

ZOMEGAP(:) = 1. + (ZEPSM(:,i)/ZMSFT(:,i)) / &
            (ZEPSM(:,i)/ZMSFT(:,i) + ((1+ZOMEGAG(:))*(1-ZEPSM(:,i))) / ZMSAL(:,i))
                       
ZQUI(:) =  (3.*ZCSFTM(:)* RR * TEMP(:) * ZTM(:)) / &  
           (ZTO(:) *(ZCSFTM(:) + ZBETA(:)*(1.+ZOMEGAG(:)) ))  
            
ZQUI(:) = MIN(ZQUI(:), 0.9999)
ZQUI(:) =  (ZTO(:)/(SIGMAWV(:)*1E3)) * (1. - ZQUI(:)) 
ZQUI(:) = ZQUI(:)**0.5  

! Final Critical radius of dropplet
ZR(:) = ((3.*B(:,i)*(R0AP(:,i))**3./A0(:))**(1./2.)) / ZQUI(:)
ZR(:) = MAX(ZR(:), 1E-10)
ZR1(:) = ZR(:)

enddo

ZX(:) = 1.5 * ((ZTO(:)*ZQUI(:))/(SIGMAWV(:)*1E3) - (ZQUI(:)**3.) / 3.) 

! Environmental supersaturation 
WHERE (B(:,i) .GT. 0.)
  SM(:,i) =  ZX(:) * ((4.*A0(:)**3.)/(27.*B(:,i)*R0AP(:,i)**3.))**0.5  
ELSEWHERE
  SM(:,i) = 0.
ENDWHERE

! Kelvin final term
A(:,i) = ((27.* R0AP(:,i)**3. * SM(:,i)**2. * B(:,i)) / 4.)**(1./3.)  

alpha(:) = (chlav(:)/(RV*TEMP(:)*TEMP(:)*CP) &
                -(1./(RD*TEMP(:))))*GAPS
gama(:) = (RR*TEMP(:))/(esw(:)*xmw*1.e-3)+&
              (chlav(:)*chlav(:)*xmw/(CP*P(:)*100.*TEMP(:)*xma))
cc0(:) = (alpha(:)*abs(w(:))+PSIRAD(:)*PTRAD(:))*COND(:)

dzeta(:) = 2./3.*SQRT(cc0(:))*A(:,i)

eta(:,i) = (cc0(:)**1.5/(2.*RHOW*PI*gama(:)))/XNAP(:,i)
f(:,i) = 0.5*exp(2.5*LOG(PSIG(:,i))**2.)
g1(:,i) = 1.+0.25*LOG(PSIG(:,i))

WHERE (B(:,i) .GT. 0.) 
ratio(:,i) = (f(:,i)*(dzeta/eta(:,i))**1.5+&
                  (g1(:,i)*(SM(:,i)*SM(:,i)/  &
                  (eta(:,i)+(3.*dzeta)))**0.75))/(SM(:,i)*SM(:,i))
ENDWHERE

somratio(:) = somratio(:) + ratio(:,i)

ENDDO

! Maximum supersaturation
WHERE (somratio(:).GT.0.) 
  PSMAX(:) = 1./SQRT(somratio(:))
ELSEWHERE
  PSMAX(:) = 0.
ENDWHERE

!  3. Calculation of cloud droplet number and mass concentration
!     ==========================================================
!
XNNUCS(:,:) = 0.
xmapac(:,:) = 0.
DO i = 1, kap
 
! number of activated aerosol pp/m3
  U(:,i) = 2.*LOG(SM(:,i)/PSMAX(:))/(3.*SQRT(2.)*LOG(PSIG(:,i)))
  XNNUCS(:,i) = 0.5*XNAP(:,i)*(1.-erfa(U(:,i)))
  XNNUC(:) = XNNUC(:) + XNNUCS(:,i)

! calculation of activated aerosols mass (kg/m3)
  xmap(:,i) = PRHOP*VAP(:,i) 
  um(:,i) =u(:,i)-3./sqrt(2.)*LOG(PSIG(:,i)) 
  xmapac(:,i) = 0.5*xmap(:,i)*(1.-erfa(um(:,i)))
  
  ZRCI(:,i) = PRG(:,i) * (SM(:,i)/PSMAX(:))**(2./3.)
ENDDO

PNCN(:,:) =  MAX(XNNUCS(:,:), 0.)
PMCN(:,:) =  MAX(xmapac(:,:), 0.)
IF(PRESENT(PT0)) PT0(:)= ZTO(:) !! en g/s² Variable for diaprog synchronic

!!Calcul de la fraction activée par mode
IF(PRESENT(PFRACT_ACT1)) PFRACT_ACT1(:) = (PNCN(:,1) / PN0(:,1))*100 !! en % Variable for diaprog synchronic
IF(PRESENT(PFRACT_ACT2)) PFRACT_ACT2(:) = (PNCN(:,2) / PN0(:,2))*100 !! en % Variable for diaprog synchronic


IF(PRESENT(PRAOULT1)) PRAOULT1(:) = B(:,1) !!Raoult term mode 1
IF(PRESENT(PRAOULT2)) PRAOULT2(:) = B(:,2) !!Raoult term mode 2
IF(PRESENT(PKELVIN1)) PKELVIN1(:) = A(:,1) !!Kelvin term mode 1
IF(PRESENT(PKELVIN2)) PKELVIN2(:) = A(:,2) !!Kelvin term mode 2

IF(PRESENT(PR0AP1)) PR0AP1(:) = R0AP(:,1) !!Aerosol mean radius (m) mode 1
IF(PRESENT(PR0AP2)) PR0AP2(:) = R0AP(:,2) !!Aerosol mean radius (m) mode 2


CONTAINS

FUNCTION ERFA(X) RESULT(PERF)
USE MODI_ERF
REAL, DIMENSION(:), INTENT(IN) :: X
REAL, DIMENSION(SIZE(X,1)) :: PERF
INTEGER :: II
DO II=1,SIZE(X,1)
PERF(II) = erf(X(II))
ENDDO
END FUNCTION ERFA

END SUBROUTINE ABDULRAZZAK
