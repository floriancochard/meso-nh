!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/Attic/mode_salt_psd.f90,v $ $Revision: 1.1.2.1.2.1.2.1.2.1 $ $Date: 2013/07/12 13:55:08 $
!-----------------------------------------------------------------
!!   ########################
     MODULE MODE_SALT_PSD_WET
!!   ########################
!!
!!    PURPOSE
!!    -------
!! MODULE SALT PSD (Particle Size Distribution)
!! Purpose: Contains subroutines to convert from transported variables (ppp)
!! to understandable aerosol variables, e.g. #/m3, kg/m3, sigma, R_{n}
!!
!!    AUTHOR
!!    ------
!!      Alf Grini (CNRM/GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!      M. Claeys - (CNRM-GMEI) 2015
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
USE MODD_CSTS_SALT         !Constants which are important for sea salt calculations
USE MODD_SALT              !Dust module which contains even more constants
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XAVOGADRO        & ![molec/mol] avogadros number
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD             & !  Cpd (dry air)
      ,XRHOLW           & ! Densité de l'eau
      ,XMV              & ! Molar weight of water
      ,XALPI            &
      ,XBETAI           &
      ,XGAMI            &
      ,XTT             
USE MODD_CST, ONLY : XMNH_TINY
USE MODE_THERMO             ! Pour calcul de la pression de vapeur saturante
USE MODD_PARAM_n,    ONLY : CCLOUD


!
IMPLICIT NONE
!
CONTAINS
!
!!   ############################################################
  SUBROUTINE PPP2SALT_WET(          &
       PSVT                         & !I [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PPABST                     & !I Pression
       , PTHT                       & !I Potential temperature
       , PRT                        & !I Large scale vapor mixing ratio
       , PSIG3D                     & !O [-] standard deviation of aerosol distribution
       , PRG3D                      & !O [um] number median wet radius of aerosol distribution
       , PN3D                       & !O [#/m3] number concentration of aerosols
       , PMASS3D                    & !O [kg/m3]wet mass concentration of aerosol
       , PM3D                       & !O aerosols moments 0, 3 and 6
       , PDENSITY_WET               & !O [g/m2] density of wet aerosol (water + salt)
       )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the three moments M0, M3 and M6 given in ppp into
!!    Values which can be understood more easily (R, sigma, N, M)
!!    
!!    Calcul the wet radius of the particles, using RH and Gerber (1985) relation
!!    The mass of the aerosols is calculated using the new radius and the
!density of water and salt
!!
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AEROS(PSVT, PRHODREF, PSIG3D=SIGVAR,  &
!!       PRG3D=RVAR, PN3D=NVAR, PM3D=MASSVAR)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    2005 Alf Grini (CNRM)
!!    2006 Jean-Pierre Chaboureau (LA)
!!    2015 Marine Claeys (CNRM)
!!    EXTERNAL
!!    --------
!!    None
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:,:),           INTENT(INOUT) :: PSVT      !I [ppp] first moment
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRHODREF  !I [kg/m3] density of air
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PPABST    !I Pression
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PTHT      !I Potential temperature
REAL, DIMENSION(:,:,:,:),           INTENT(IN)    :: PRT       !I Large scale vapor mixing ratio

REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PSIG3D   !O [-] standard deviation
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PRG3D    !O [um] number median radius 
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PN3D     !O [#/m3] number concentration
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PMASS3D  !O [kg_{aer}/m3] wet mass concentration
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PM3D     !O aerosols moments 
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PDENSITY_WET  !O Density of wet aerosol (water + salt) 
!
!
!*      0.2    declarations local variables
!
REAL                                  :: ZRHOI    ! [kg/m3] density of aerosol
REAL                                  :: ZRHOLW   ! [kg/m3] density of water
REAL                                  :: ZMI      ! [kg/mol] molar weight of aerosol
REAL                                  :: ZMV      ! [kg/mol] molar weight of water
REAL                                  :: ZRGMIN   ! [um] minimum radius accepted
REAL                                  :: ZSIGMIN  ! minimum standard deviation accepted

REAL, PARAMETER                       :: C1 = 0.7674
REAL, PARAMETER                       :: C2 = 3.079
REAL, PARAMETER                       :: C3 = 2.572E-11
REAL, PARAMETER                       :: C4 = -1.424

REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM       ! [aerosol units] local array which goes to output later
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSV      ! [sea salts moment concentration]
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN    ! [aerosol units] minimum values for N, sigma, M
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM0      ! [idx] index for Mode 0 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM3      ! [idx] indexes for Mode 3 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM6      ! [idx] indexes for Mode 6 in passed variables
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS      ! initial mean radius
INTEGER                               :: JN,IMODEIDX,JJ  ! [idx] loop counters

REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3), NMODE_SLT)  :: ZMASS3D, &
                                                                       ZMASS3D_SLT         ![kg/m3] mass of one sea salt mode

REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3))  :: ZTEMP, ZREHU, ZREHU_tmp
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3), NMODE_SLT)  :: ZDENSITY_WET      ! [g/m2] Aerosol wet density (salt + water) 
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3))  :: ZSIGMA            ! [-] standard deviation 
REAL,DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3))  :: ZRG, ZRG_WET      ! [um] number median radius, and nuber median wet radius 
!
!-------------------------------------------------------------------------------
!

!+ Marine
!  Calcul de RH
! Pris dans write_lfi_for_diag pour le calcul de RH
ZTEMP(:,:,:) = PTHT(:,:,:) * (PPABST(:,:,:) / XP00)**(XRD/XCPD)

ZREHU_tmp(:,:,:) = SM_FOES(ZTEMP(:,:,:))  ! SM_FOES = to compute saturation vapor pressure
ZREHU_tmp(:,:,:) = (XMV / XMD) * ZREHU_tmp(:,:,:) / (PPABST(:,:,:) - ZREHU_tmp(:,:,:))
!XMD,XMV      ! Molar mass of dry air and molar mass of vapor, PPABST: pression

ZREHU(:,:,:) = PRT(:,:,:,1) / ZREHU_tmp(:,:,:)

IF (CCLOUD(1:3) =='ICE' .OR. CCLOUD =='C3R5')  THEN
  WHERE ( ZTEMP(:,:,:) <  XTT) ! XTT : Triple point temperature
    ZREHU_tmp(:,:,:) = EXP( XALPI - XBETAI/ZTEMP(:,:,:) & ! XALPI,XBETAI,XGAMI ! Constants for saturation vapor pressure 
                       - XGAMI*ALOG(ZTEMP(:,:,:)) ) !saturation over ice
    ZREHU_tmp(:,:,:) = (XMV / XMD) * ZREHU_tmp(:,:,:) / (PPABST(:,:,:) - ZREHU_tmp(:,:,:))
    ZREHU(:,:,:) = PRT(:,:,:,1) / ZREHU_tmp(:,:,:)
  END WHERE
END IF

ZREHU(:,:,:) = MIN(MAX(ZREHU(:,:,:), 0.02),0.95)


!        1.1    initialisation 
!
!Calculations here are for one mode only
!
ALLOCATE (NM0(NMODE_SLT))
ALLOCATE (NM3(NMODE_SLT))
ALLOCATE (NM6(NMODE_SLT))
ALLOCATE (ZM(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), NMODE_SLT*3))
ALLOCATE (ZMMIN(NMODE_SLT*3))
ALLOCATE (ZSV(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), SIZE(PSVT,4)))
ALLOCATE (ZINIRADIUS(NMODE_SLT))
    
ZSV(:,:,:,:) = MAX(PSVT(:,:,:,:), XMNH_TINY)

DO JN = 1, NMODE_SLT
  IMODEIDX = JPSALTORDER(JN)
  !Calculations here are for one mode only
  IF (CRGUNITS == "MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF

  !Set counter for number, M3 and M6
  NM0(JN) = 1 + (JN - 1) * 3
  NM3(JN) = 2 + (JN - 1) * 3
  NM6(JN) = 3 + (JN - 1) * 3
  !Get minimum values possible
  ZMMIN(NM0(JN)) = XN0MIN_SLT(IMODEIDX)
  ZRGMIN         = ZINIRADIUS(JN)
  IF (LVARSIG_SLT) THEN
    ZSIGMIN = XSIGMIN_SLT
  ELSE
    ZSIGMIN = XINISIG_SLT(IMODEIDX)
  ENDIF
  ZMMIN(NM3(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2) 
  ZMMIN(NM6(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
END DO
!
!Set density of aerosol, here sea salt (kg/m3) and water
ZRHOI = XDENSITY_SALT
ZRHOLW = XRHOLW
!Set molecular weight of sea salt and water!NOTE THAT THIS IS NOW IN KG
ZMI   = XMOLARWEIGHT_SALT
ZMV   = XMV
!
!
DO JN = 1, NMODE_SLT
  !
  IF (LVARSIG_SLT) THEN ! give M6 (case of variable standard deviation)
  !
  !Get number concentration (#/molec_{air}==>#/m3)
    ZM(:,:,:,NM0(JN))=                         &
         ZSV(:,:,:,1+(JN-1)*3)                 & !#/molec_{air}
         * XAVOGADRO                           & !==>#/mole
         / XMD                                 & !==>#/kg_{air}
         * PRHODREF(:,:,:)                       !==>#/m3
  ! 
  !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
    ZM(:,:,:,NM3(JN)) =                        &
         ZSV(:,:,:,2+(JN-1)*3)                 & !molec_{aer}/molec_{aer}
         * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
         * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
         * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
         * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
         / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
     !Limit mass concentration to minimum value
    ZM(:,:,:,NM3(JN)) = MAX(ZM(:,:,:,NM3(JN)), ZMMIN(NM3(JN)))
  ! 
    ZM(:,:,:,NM6(JN)) = ZSV(:,:,:,3+(JN-1)*3)  & !um6/molec_{air}*(cm3/m3)
         * 1.d-6                               & !==> um6/molec_{air}
         * XAVOGADRO                           & !==> um6/mole_{air}
         / XMD                                 & !==> um6/kg_{air}
         * PRHODREF(:,:,:)                       !==> um6/m3_{air}
     !Limit m6 concentration to minimum value
    ZM(:,:,:,NM6(JN)) =  MAX(ZM(:,:,:,NM6(JN)), ZMMIN(NM6(JN)))
  !
  !Get sigma (only if sigma is allowed to vary)
    !Get intermediate values for sigma M3^2/(M0*M6) (ORILAM paper, eqn 8)
    ZSIGMA(:,:,:)=ZM(:,:,:,NM3(JN))**2/(ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM6(JN)))
    !Limit the intermediate value, can not be larger than 1
    ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
    !Limit the value for intermediate, can not be smaller than 0
    ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
    !Calculate log(sigma)
    ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
    !Finally get the real sigma the negative sign is because of 
    !The way the equation is written (M3^2/(M0*M6)) instead of (M0*M6)/M3^3
    ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
    !Limit the value to reasonable ones
    ZSIGMA(:,:,:) =  MAX( XSIGMIN_SLT, MIN( XSIGMAX_SLT, ZSIGMA(:,:,:) ) )

  !
    !Put back M6 so that it fits the sigma which is possibly modified above
    !The following makes M6 consistent with N, R, SIGMA
    ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN)) &
         * ( (ZM(:,:,:,NM3(JN))/ZM(:,:,:,NM0(JN)))**(1./3.) &
         * exp(-(3./2.)*log(ZSIGMA(:,:,:))**2))**6 &
         * exp(18.*log(ZSIGMA(:,:,:))**2)

  ELSE ! compute M6 from M0, M3 and SIGMA
    ! 
    ZSIGMA(:,:,:) = XINISIG_SLT(JPSALTORDER(JN))
    IF (LRGFIX_SLT) THEN

    !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
    ZM(:,:,:,NM3(JN)) =                        &
         ZSV(:,:,:,JN)                         & !molec_{aer}/molec_{aer}
         * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
         * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
         * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
         * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
         / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)

    ZM(:,:,:,NM3(JN)) = MAX(ZM(:,:,:,NM3(JN)), ZMMIN(NM3(JN)))
       PSVT(:,:,:,JN) = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)

    ZM(:,:,:,NM0(JN))=  ZM(:,:,:,NM3(JN))/&
       ((ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(XINISIG_SLT(JPSALTORDER(JN)))**2))



    ELSE

    !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
      ZM(:,:,:,NM3(JN)) =                        &
           ZSV(:,:,:,2+(JN-1)*2)                 & !molec_{aer}/molec_{aer}
           * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
           * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
           * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
           * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
           / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)

 
    !Get number concentration (#/molec_{air}==>#/m3)
       ZM(:,:,:,NM0(JN))=                         &
           ZSV(:,:,:,1+(JN-1)*2)                 & !#/molec_{air}
           * XAVOGADRO                           & !==>#/mole
           / XMD                                 & !==>#/kg_{air}
           * PRHODREF(:,:,:)                       !==>#/m3

    ! Limit concentration to minimum values
      WHERE ((ZM(:,:,:,NM0(JN)) < ZMMIN(NM0(JN)) ).OR. &
             (ZM(:,:,:,NM3(JN)) < ZMMIN(NM3(JN)) )) 
         ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
         ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
         PSVT(:,:,:,1+(JN-1)*2) = ZM(:,:,:,NM0(JN)) * XMD / &
         (XAVOGADRO * PRHODREF(:,:,:) )
         PSVT(:,:,:,2+(JN-1)*2) = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                                (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
      ENDWHERE
    END IF

    ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN))                   &
         * ( (ZM(:,:,:,NM3(JN))/ZM(:,:,:,NM0(JN)))**(1./3.) &
         * exp(-(3./2.)*log(ZSIGMA(:,:,:))**2))**6          &
         * exp(18.*log(ZSIGMA(:,:,:))**2)
    
  !
  END IF
  !   
  !Get number median radius (eqn. 7 in Orilam manuscript)
  ZRG(:,:,:)=    &
          (      &
          ZM(:,:,:,NM3(JN))*ZM(:,:,:,NM3(JN))*ZM(:,:,:,NM3(JN))*ZM(:,:,:,NM3(JN))    &
          /(ZM(:,:,:,NM6(JN))*ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM0(JN))*ZM(:,:,:,NM0(JN))) &
          )                                                                          &
          ** XSIXTH_SALT 


  !ZRG(:,:,:)=MIN(ZRG(:,:,:),ZINIRADIUS(JN))
  !Give the sigma-values to the passed array
  IF(PRESENT(PSIG3D)) PSIG3D(:,:,:,JN) = ZSIGMA(:,:,:)
  !
  !Set the number concentrations in the passed array
  IF(PRESENT(PN3D)) PN3D(:,:,:,JN) = ZM(:,:,:,NM0(JN))
  !
!  !Get the number median radius
!  IF(PRESENT(PRG3D)) PRG3D(:,:,:,JN)= ZRG(:,:,:)
  !
  !
  !+ Marine
!!!!!!!!!!! Wet radius calculus

!Test pour Marine

  ZRG_WET(:,:,:) =  C1 * (ZRG(:,:,:)*1.d-4)**C2 ! Pour le calcul, ZRG en cm! (d'où 1.d-4)

!+ Marine test

  ZRG_WET(:,:,:) =  ZRG_WET(:,:,:) / (C3 * ((ZRG(:,:,:)*1.d-4)**C4) - LOG10(ZREHU(:,:,:)))
  ZRG_WET(:,:,:) =  ZRG_WET(:,:,:) + (ZRG(:,:,:)*1.d-4)**3 
  ZRG_WET(:,:,:) = ( ZRG_WET(:,:,:)**(1./3) )*1.d4 ! *1.d4 pour repasser de cm à micromètres

  !Get the number median radius
  IF(PRESENT(PRG3D)) PRG3D(:,:,:,JN) = ZRG_WET(:,:,:)



! Wet density calcul
  ZDENSITY_WET(:,:,:,JN)=(ZRHOI * ZRG(:,:,:) * ZRG(:,:,:) * ZRG(:,:,:) + &
          ZRHOLW * (ZRG_WET(:,:,:) * ZRG_WET(:,:,:) * ZRG_WET(:,:,:)- &
          ZRG(:,:,:) * ZRG(:,:,:) * ZRG(:,:,:))) &
         / (ZRG_WET(:,:,:) * ZRG_WET(:,:,:) * ZRG_WET(:,:,:)) 

!Wet mass
  ZMASS3D(:,:,:,JN)=    &
          ZM(:,:,:,NM0(JN)) & !#/m^3_{air}
           * XPI*4./3.       &
           * ZDENSITY_WET(:,:,:,JN)           &    !==>kg/m^3_{aeros}/m^3_{air}
           * ZRG_WET(:,:,:) * ZRG_WET(:,:,:) * ZRG_WET(:,:,:) &
           * XUM3TOM3_SALT        &    !==>kg/m^3_{air}
           * exp(4.5*log(ZSIGMA(:,:,:))*log(ZSIGMA(:,:,:)))
 
! Salt Mass
  ZMASS3D_SLT(:,:,:,JN)=     &
              ZM(:,:,:,NM0(JN)) &    !#/m^3_{air}
              * XPI*4./3.       &
              * ZRHOI           &    !==>kg/m^3_{aeros}/m^3_{air}
              * ZRG(:,:,:) * ZRG(:,:,:) * ZRG(:,:,:) &
              * XUM3TOM3_SALT        &    !==>kg/m^3_{air}
              * exp(4.5*log(ZSIGMA(:,:,:))*log(ZSIGMA(:,:,:)))

  IF(PRESENT(PMASS3D)) THEN 
    PMASS3D(:,:,:,JN)= ZMASS3D(:,:,:,JN)
  ENDIF

  IF(PRESENT(PDENSITY_WET)) THEN
    PDENSITY_WET(:,:,:,JN) = ZDENSITY_WET(:,:,:,JN)
  ENDIF
!
END DO  !Loop on modes
!
IF(PRESENT(PM3D)) PM3D(:,:,:,:) = ZM(:,:,:,:)
!
DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZSV)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZM)
DEALLOCATE(NM6)
DEALLOCATE(NM3)
DEALLOCATE(NM0)

!+ Marine

END SUBROUTINE PPP2SALT_WET

!!   ############################################################
  SUBROUTINE SALT2PPP(             &
       PSVT                         & !IO [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG3D                     & !I [-] standard deviation of aerosol distribution
       , PRG3D                      & !I [um] number median diameter of aerosol distribution
       )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the sea salt Mass, RG and SIGMA in the  three moments M0, M3 and M6 given in ppp 
!! 
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AEROS(PSVT, PRHODREF, PSIG3D=SIGVAR,  &
!!       PRG3D=RVAR, PN3D=NVAR)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Alf Grini (CNRM)
!!
!!    EXTERNAL
!!    --------
!!    None
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
    !INPUT
    REAL,       DIMENSION(:,:,:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air
    REAL,       DIMENSION(:,:,:,:),  INTENT(IN)     :: PSIG3D   !O [-] standard deviation
    REAL,       DIMENSION(:,:,:,:),  INTENT(IN)     :: PRG3D    !O [um] number median diameter

    !OUTPUT
    REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)  :: PSVT  !IO [#/molec_{air}] first moment
                                                                !IO [molec_{aer}/molec_{air} 3rd moment
                                                                !IO [um6/molec_{air}*(cm3/m3)] 6th moment
!
!
!*      0.2    declarations local variables
!
    REAL                                  :: ZRHOI               ! [kg/m3] density of aerosol
    REAL                                  :: ZMI                 ! [kg/mol] molar weight of aerosol
    REAL                                  :: ZRGMIN              ! [um] minimum radius accepted
    REAL                                  :: ZSIGMIN             ! minimum standard deviation accepted
    REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM                  ! [aerosol units] local array which goes to output later
    REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZSIGMA              ! aersol standard deviation
    REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN               ! [aerosol units] minimum values for N, sigma, M
    REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS          ! initial mean radius
    INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM0                 ! [idx] index for Mode 0 in passed variables
    INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM3                 ! [idx] indexes for Mode 3 in passed variables
    INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM6                 ! [idx] indexes for Mode 6 in passed variables
    INTEGER                               :: JJ, JN              ! [idx] loop counters
    INTEGER                               :: IMODEIDX
!
!-------------------------------------------------------------------------------
!
!        1.1    initialisation 


    ALLOCATE (NM0(NMODE_SLT))
    ALLOCATE (NM3(NMODE_SLT))
    ALLOCATE (NM6(NMODE_SLT))
    ALLOCATE (ZINIRADIUS(NMODE_SLT))
    ALLOCATE (ZMMIN(NMODE_SLT*3))
    ALLOCATE (ZM(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), NMODE_SLT*3))
    ALLOCATE (ZSIGMA(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)))

    !Set density of aerosol, here sea salt (kg/m3)
    ZRHOI = XDENSITY_SALT
    !Set molecular weight of sea salt !NOTE THAT THIS IS NOW IN KG
    ZMI   = XMOLARWEIGHT_SALT
!

    ! PSVT need to be positive
!    PSVT(:,:,:,:) = MAX(PSVT(:,:,:,:), XMNH_TINY)
    
    DO JN=1,NMODE_SLT
      IMODEIDX = JPSALTORDER(JN)
    !Calculations here are for one mode only
      IF (CRGUNITS=="MASS") THEN
        ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
      ELSE
        ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
      END IF

    !Set counter for number, M3 and M6
      NM0(JN) = 1+(JN-1)*3
      NM3(JN) = 2+(JN-1)*3
      NM6(JN) = 3+(JN-1)*3

    !Get minimum values possible
      ZMMIN(NM0(JN)) = XN0MIN_SLT(IMODEIDX)
      ZRGMIN     =  ZINIRADIUS(JN)
      IF (LVARSIG_SLT) THEN
        ZSIGMIN = XSIGMIN_SLT
      ELSE
        ZSIGMIN = XINISIG_SLT(IMODEIDX)
      ENDIF
      ZMMIN(NM3(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2) 
      ZMMIN(NM6(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
    END DO

    !Set density of aerosol, here sea salt (kg/m3)
    ZRHOI = XDENSITY_SALT
    !Set molecular weight of sea salt !NOTE THAT THIS IS NOW IN KG
    ZMI   = XMOLARWEIGHT_SALT
!
    DO JN=1,NMODE_SLT
     !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
     IF (LVARSIG_SLT) THEN
       ZM(:,:,:,NM3(JN)) =                        &
            PSVT(:,:,:,2+(JN-1)*3)                & !molec_{aer}/molec_{aer}
            * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
            * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
            * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
            * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
            / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
    ELSE 
      IF ((LRGFIX_SLT)) THEN
        ZM(:,:,:,NM3(JN)) =                   &
            PSVT(:,:,:,JN)                        & !molec_{aer}/molec_{aer}
            * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
            * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
            * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
            * XM3TOUM3_SALT                       & !==>um3_{aer}/m3_{air}
            / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
           ZM(:,:,:,NM3(JN)) = MAX(ZM(:,:,:,NM3(JN)), ZMMIN(NM3(JN)))
      ELSE
        ZM(:,:,:,NM3(JN)) =                        &
            PSVT(:,:,:,2+(JN-1)*2)                & !molec_{aer}/molec_{aer}
            * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
            * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
            * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
            * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
            / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
      END IF
    END IF
! calculate moment 0 from dispersion and mean radius
    ZM(:,:,:,NM0(JN))=  ZM(:,:,:,NM3(JN))/&
       ((PRG3D(:,:,:,JN)**3)*EXP(4.5 * LOG(PSIG3D(:,:,:,JN))**2))

! calculate moment 6 from dispersion and mean radius
    ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN)) * (PRG3D(:,:,:,JN)**6) * &
               EXP(18 *(LOG(PSIG3D(:,:,:,JN)))**2)

    IF (LVARSIG_SLT) THEN
      WHERE ((ZM(:,:,:,NM0(JN)) .LT. ZMMIN(NM0(JN))).OR.&
             (ZM(:,:,:,NM3(JN)) .LT. ZMMIN(NM3(JN))).OR.&
             (ZM(:,:,:,NM6(JN)) .LT. ZMMIN(NM6(JN))))
        ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
        ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
        ZM(:,:,:,NM6(JN)) = ZMMIN(NM6(JN))
      END WHERE
    ELSE  IF (.NOT.(LRGFIX_SLT)) THEN

      WHERE ((ZM(:,:,:,NM0(JN)) .LT. ZMMIN(NM0(JN))).OR.&
             (ZM(:,:,:,NM3(JN)) .LT. ZMMIN(NM3(JN))))
        ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
        ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
      END WHERE
    ENDIF

     
     ! return to concentration #/m3 =>  (#/molec_{air}
    IF (LVARSIG_SLT) THEN
      PSVT(:,:,:,1+(JN-1)*3) = ZM(:,:,:,NM0(JN)) * XMD / &
                               (XAVOGADRO*PRHODREF(:,:,:))

      PSVT(:,:,:,2+(JN-1)*3) = ZM(:,:,:,NM3(JN)) * XMD  * XPI * 4./3 * ZRHOI / &
                               (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)

      PSVT(:,:,:,3+(JN-1)*3) = ZM(:,:,:,NM6(JN)) * XMD  / &
                               ( XAVOGADRO*PRHODREF(:,:,:) * 1.d-6) 
    ELSE IF (LRGFIX_SLT) THEN
      PSVT(:,:,:,JN)         = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                               (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
    ELSE
      PSVT(:,:,:,1+(JN-1)*2) = ZM(:,:,:,NM0(JN)) * XMD / &
                               (XAVOGADRO*PRHODREF(:,:,:))

      PSVT(:,:,:,2+(JN-1)*2) = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                               (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
    END IF
!
END DO  !Loop on modes

DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZM)
DEALLOCATE(NM6)
DEALLOCATE(NM3)
DEALLOCATE(NM0)
!
END SUBROUTINE SALT2PPP
!
!!   ############################################################
  SUBROUTINE PPP2SALT1D(             &
       PSVT                         & !I [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG1D                     & !O [-] standard deviation of aerosol distribution
       , PRG1D                      & !O [um] number median diameter of aerosol distribution
       , PN1D                       & !O [#/m3] number concentration of aerosols
       , PMASS1D                    & !O [kg/m3] mass concentration of aerosol
       , PM1D                       & !O aerosols moments 0, 3 and 6
       )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the three moments M0, M3 and M6 given in ppp into
!!    Values which can be understood more easily (R, sigma, N, M)
!! 
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AEROS(PSVT, PRHODREF, PSIG3D=SIGVAR,  &
!!       PRG3D=RVAR, PN3D=NVAR, PM3D=MASSVAR)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    2005 Alf Grini (CNRM)
!!    2006 Jean-Pierre Chaboureau (LA)
!!
!!    EXTERNAL
!!    --------
!!    None
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,       DIMENSION(:,:),  INTENT(INOUT)  :: PSVT      !I [ppp] first moment
REAL,       DIMENSION(:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air

REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PSIG1D   !O [-] standard deviation
REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PRG1D    !O [um] number median diameter
REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PN1D     !O [#/m3] number concentration
REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PMASS1D  !O [kg_{aer}/m3] mass concentration
REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PM1D     !O aerosols moments 
!
!
!*      0.2    declarations local variables
!
REAL                                  :: ZRHOI               ! [kg/m3] density of aerosol
REAL                                  :: ZMI                 ! [kg/mol] molar weight of aerosol
REAL                                  :: ZRGMIN              ! [um] minimum radius accepted
REAL                                  :: ZSIGMIN             ! minimum standard deviation accepted
REAL,DIMENSION(:,:), ALLOCATABLE  :: ZM                  ! [aerosol units] local array which goes to output later
REAL,DIMENSION(:,:), ALLOCATABLE  :: ZSV                 ! [sea salts moment concentration]
REAL,DIMENSION(:),   ALLOCATABLE  :: ZSIGMA              ! [-] standard deviation
REAL,DIMENSION(:),   ALLOCATABLE  :: ZRG                 ! [um] number median diameter
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN               ! [aerosol units] minimum values for N, sigma, M
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM0                 ! [idx] index for Mode 0 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM3                 ! [idx] indexes for Mode 3 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM6                 ! [idx] indexes for Mode 6 in passed variables
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS          ! initial mean radius
INTEGER                               :: JN,IMODEIDX,JJ      ! [idx] loop counters
!
!-------------------------------------------------------------------------------
!
!        1.1    initialisation 
!
!Calculations here are for one mode only
!
ALLOCATE (NM0(NMODE_SLT))
ALLOCATE (NM3(NMODE_SLT))
ALLOCATE (NM6(NMODE_SLT))
ALLOCATE (ZM(SIZE(PSVT,1), NMODE_SLT*3))
ALLOCATE (ZMMIN(NMODE_SLT*3))
ALLOCATE (ZSIGMA(SIZE(PSVT,1)))
ALLOCATE (ZRG(SIZE(PSVT,1)))
ALLOCATE (ZSV(SIZE(PSVT,1), SIZE(PSVT,2)))
ALLOCATE (ZINIRADIUS(NMODE_SLT))
    
!ZSV(:,:) = MAX(PSVT(:,:), XMNH_TINY)

DO JN=1,NMODE_SLT
  IMODEIDX = JPSALTORDER(JN)
  !Calculations here are for one mode only
  IF (CRGUNITS=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF

  !Set counter for number, M3 and M6
  NM0(JN) = 1+(JN-1)*3
  NM3(JN) = 2+(JN-1)*3
  NM6(JN) = 3+(JN-1)*3
  !Get minimum values possible
  ZMMIN(NM0(JN)) = XN0MIN_SLT(IMODEIDX)
  ZRGMIN         = ZINIRADIUS(JN)
  IF (LVARSIG_SLT) THEN
    ZSIGMIN = XSIGMIN_SLT
  ELSE
    ZSIGMIN = XINISIG_SLT(IMODEIDX)
  ENDIF
  ZMMIN(NM3(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2) 
  ZMMIN(NM6(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
END DO
!
!Set density of aerosol, here sea salt (kg/m3)
ZRHOI = XDENSITY_SALT
!Set molecular weight of sea salt !NOTE THAT THIS IS NOW IN KG
ZMI   = XMOLARWEIGHT_SALT
!
!
DO JN=1,NMODE_SLT
  !
  IF (LVARSIG_SLT) THEN ! give M6 (case of variable standard deviation)
  !
  !Get number concentration (#/molec_{air}==>#/m3)
    ZM(:,NM0(JN))=                         &
         ZSV(:,1+(JN-1)*3)                 & !#/molec_{air}
         * XAVOGADRO                           & !==>#/mole
         / XMD                                 & !==>#/kg_{air}
         * PRHODREF(:)                       !==>#/m3
  ! 
  !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
    ZM(:,NM3(JN)) =                        &
         ZSV(:,2+(JN-1)*3)                 & !molec_{aer}/molec_{aer}
         * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
         * PRHODREF(:)                     & !==>kg_{aer}/m3_{air}
         * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
         * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
         / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
     !Limit mass concentration to minimum value
    ZM(:,NM3(JN)) = MAX(ZM(:,NM3(JN)), ZMMIN(NM3(JN)))
  ! 
    ZM(:,NM6(JN)) = ZSV(:,3+(JN-1)*3)  & !um6/molec_{air}*(cm3/m3)
         * 1.d-6                               & !==> um6/molec_{air}
         * XAVOGADRO                           & !==> um6/mole_{air}
         / XMD                                 & !==> um6/kg_{air}
         * PRHODREF(:)                       !==> um6/m3_{air}
     !Limit m6 concentration to minimum value
    ZM(:,NM6(JN)) =  MAX(ZM(:,NM6(JN)), ZMMIN(NM6(JN)))
  !
  !Get sigma (only if sigma is allowed to vary)
    !Get intermediate values for sigma M3^2/(M0*M6) (ORILAM paper, eqn 8)
    ZSIGMA(:)=ZM(:,NM3(JN))**2/(ZM(:,NM0(JN))*ZM(:,NM6(JN)))
    !Limit the intermediate value, can not be larger than 1
    ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
    !Limit the value for intermediate, can not be smaller than 0
    ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
    !Calculate log(sigma)
    ZSIGMA(:)= LOG(ZSIGMA(:))
    !Finally get the real sigma the negative sign is because of 
    !The way the equation is written (M3^2/(M0*M6)) instead of (M0*M6)/M3^3
    ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
    !Limit the value to reasonable ones
    ZSIGMA(:) =  MAX( XSIGMIN_SLT, MIN( XSIGMAX_SLT, ZSIGMA(:) ) )

  !
    !Put back M6 so that it fits the sigma which is possibly modified above
    !The following makes M6 consistent with N, R, SIGMA
    ZM(:,NM6(JN)) = ZM(:,NM0(JN)) &
         * ( (ZM(:,NM3(JN))/ZM(:,NM0(JN)))**(1./3.) &
         * exp(-(3./2.)*log(ZSIGMA(:))**2))**6 &
         * exp(18.*log(ZSIGMA(:))**2)

  ELSE ! compute M6 from M0, M3 and SIGMA
    ! 
    ZSIGMA(:) = XINISIG_SLT(JPSALTORDER(JN))
    IF (LRGFIX_SLT) THEN

    !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
      ZM(:,NM3(JN)) =                        &
           ZSV(:,JN)                         & !molec_{aer}/molec_{aer}
           * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
           * PRHODREF(:)                     & !==>kg_{aer}/m3_{air}
           * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
           * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
           / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
      ZM(:,NM3(JN)) = MAX(ZM(:,NM3(JN)), ZMMIN(NM3(JN)))

      ZM(:,NM0(JN))=  ZM(:,NM3(JN))/&
         ((ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(XINISIG_SLT(JPSALTORDER(JN)))**2))

    ELSE

    !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
      ZM(:,NM3(JN)) =                        &
           ZSV(:,2+(JN-1)*2)                 & !molec_{aer}/molec_{aer}
           * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
           * PRHODREF(:)                     & !==>kg_{aer}/m3_{air}
           * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
           * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
           / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)

 
    !Get number concentration (#/molec_{air}==>#/m3)
      ZM(:,NM0(JN))=                         &
          ZSV(:,1+(JN-1)*2)                 & !#/molec_{air}
          * XAVOGADRO                           & !==>#/mole
          / XMD                                 & !==>#/kg_{air}
          * PRHODREF(:)                       !==>#/m3

    ! Limit concentration to minimum values
      WHERE ((ZM(:,NM0(JN)) < ZMMIN(NM0(JN)) ).OR. &
             (ZM(:,NM3(JN)) < ZMMIN(NM3(JN)) )) 
        ZM(:,NM0(JN)) = ZMMIN(NM0(JN))
        ZM(:,NM3(JN)) = ZMMIN(NM3(JN))
        PSVT(:,1+(JN-1)*2) = ZM(:,NM0(JN)) * XMD / &
                            (XAVOGADRO * PRHODREF(:) )
        PSVT(:,2+(JN-1)*2) = ZM(:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                             (ZMI*PRHODREF(:)*XM3TOUM3_SALT)
      ENDWHERE

    END IF

    ZM(:,NM6(JN)) = ZM(:,NM0(JN))                   &
         * ( (ZM(:,NM3(JN))/ZM(:,NM0(JN)))**(1./3.) &
         * exp(-(3./2.)*log(ZSIGMA(:))**2))**6          &
         * exp(18.*log(ZSIGMA(:))**2)
    
  !
  END IF
  !   
  !Get number median radius (eqn. 7 in Orilam manuscript)
  ZRG(:)=    &
          (      &
          ZM(:,NM3(JN))*ZM(:,NM3(JN))*ZM(:,NM3(JN))*ZM(:,NM3(JN))    &
          /(ZM(:,NM6(JN))*ZM(:,NM0(JN))*ZM(:,NM0(JN))*ZM(:,NM0(JN))) &
          )                                                                          &
          ** XSIXTH_SALT 
  !ZRG(:)=MIN(ZRG(:),ZINIRADIUS(JN))
  !Give the sigma-values to the passed array
  IF(PRESENT(PSIG1D)) PSIG1D(:,JN) = ZSIGMA(:)
  !
  !Set the number concentrations in the passed array
  IF(PRESENT(PN1D)) PN1D(:,JN) = ZM(:,NM0(JN))
  !
  !Get the number median radius
  IF(PRESENT(PRG1D)) PRG1D(:,JN)= ZRG(:)
  !
  IF(PRESENT(PMASS1D))THEN
       PMASS1D(:,JN)=     &
            ZM(:,NM0(JN)) &    !#/m^3_{air}
            * XPI*4./3.       &    
            * ZRHOI           &    !==>kg/m^3_{aeros}/m^3_{air}
            * ZRG(:) * ZRG(:) * ZRG(:) &
            * XUM3TOM3_SALT        &    !==>kg/m^3_{air}
            * exp(4.5*log(ZSIGMA(:))*log(ZSIGMA(:)))
  ENDIF
!
END DO  !Loop on modes
!
IF(PRESENT(PM1D)) PM1D(:,:) = ZM(:,:)
!
DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZSV)
DEALLOCATE(ZRG)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZM)
DEALLOCATE(NM6)
DEALLOCATE(NM3)
DEALLOCATE(NM0)
!
END SUBROUTINE PPP2SALT1D

!!   ############################################################
END MODULE MODE_SALT_PSD_WET
