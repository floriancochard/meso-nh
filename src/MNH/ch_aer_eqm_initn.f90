!ORILAM_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!!   ########################
     MODULE MODI_CH_AER_EQM_INIT_n
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_EQM_INIT_n(PCHEM, PAERO, PM3D, PRHOP3D, PSIG3D, PRG3D, &
                             PN3D, PRHODREF, PCTOTA)
IMPLICIT NONE
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)   :: PCHEM, PAERO
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)   :: PM3D, PRHOP3D, PSIG3D, PRG3D, PN3D
REAL,       DIMENSION(:,:,:,:,:),INTENT(INOUT)   :: PCTOTA
REAL,       DIMENSION(:,:,:),    INTENT(IN)      :: PRHODREF
END SUBROUTINE CH_AER_EQM_INIT_n
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_EQM_INIT_n
!!
!!
!!   ############################################################
     SUBROUTINE CH_AER_EQM_INIT_n(PCHEM,PAERO, PM3D, PRHOP3D, PSIG3D, PRG3D, &
                             PN3D, PRHODREF, PCTOTA)
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Realise l'equilibre entre les moments via la masse contenue 
!!    dans les aerosols, les diametres moyens et la dispersion.
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
!!    M.Leriche 2015 : masse molaire Black carbon à 12 g/mol
!  P. Wautelet 05/03/2019: modify allocation procedure for XMI and XSOLORG
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL
USE MODD_CSTS_DUST, ONLY : XDENSITY_DUST
USE MODD_CH_AERO_n
USE MODD_CH_M9_n,     ONLY : CNAMES
USE MODD_CH_MNHC_n, ONLY : LCH_INIT_FIELD
USE MODD_NSV
USE MODD_CONF
USE MODE_ll
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_BLANK , ONLY : CDUMMY1
USE MODD_CST, ONLY :       & 
          XMNH_TINY        &
         ,XAVOGADRO        & ![molec/mol] avogadros number
         ,XMD               ![kg/mol] molar weight of air
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
REAL,   DIMENSION(:,:,:,:),    INTENT(INOUT)   :: PCHEM, PAERO
REAL,   DIMENSION(:,:,:,:),    INTENT(INOUT)   :: PM3D, PRHOP3D, PSIG3D, PRG3D, PN3D
REAL,   DIMENSION(:,:,:,:,:),  INTENT(INOUT)   :: PCTOTA
REAL,   DIMENSION(:,:,:),      INTENT(IN)      :: PRHODREF


!
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),NSP+NCARB+NSOA,JPMODE) :: ZCCTOT
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3)) :: ZSUM
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3)) :: ZSIGMA
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3)) :: ZBCMINI, ZBCMINJ, ZOCMINI, ZOCMINJ, ZDSTMINI, ZDSTMINJ
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),JPMODE*3) :: ZPM
REAL,DIMENSION(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3)) :: ZPOIDS, ZWORK
REAL    ::  ZMASS, ZM6MIN
INTEGER :: JN, JJ,  JK  ! loop counter
INTEGER :: IINFO_ll
REAL    :: ZDEN2MOL, ZRHODREFMIN, ZCOEFAEROBC, ZCOEFAEROOC, ZCOEFAERODST
REAL    :: ZVALBC, ZVALOC, ZMINRGI, ZMINRGJ, ZVALDST
REAL    :: ZSUMAEROCO, ZSUMRHOD
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
!
!-------------------------------------------------------------------------------
!

!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!        1.1    initialisation 
! Index gas scheme <=> Index Orilam

DO JJ=1,SIZE(CNAMES)
IF (CNAMES(JJ) == "CO")  JP_CH_CO   = JJ
END DO

ZDEN2MOL = 1E-6 * XAVOGADRO / XMD

IF ( ASSOCIATED(XMI) ) THEN
  IF ( SIZE(XMI) == 0 ) THEN
    DEALLOCATE( XMI )
    XMI => NULL()
  END IF
END IF
IF (.NOT.(ASSOCIATED(XMI))) THEN
  ALLOCATE(XMI(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),NSP+NCARB+NSOA))
END IF

IF ( ASSOCIATED(XSOLORG) ) THEN
  IF ( SIZE(XSOLORG) == 0 ) THEN
    DEALLOCATE( XSOLORG )
    XSOLORG => NULL()
  END IF
END IF
IF (.NOT.(ASSOCIATED(XSOLORG))) THEN
  ALLOCATE(XSOLORG(SIZE(PCHEM,1),SIZE(PCHEM,2),SIZE(PCHEM,3),10))
  XSOLORG(:,:,:,:) = 0.
END IF
!
! Default values of molar mass

XMI(:,:,:,:) = 250.
XMI(:,:,:,JP_AER_SO4)  = 98.
XMI(:,:,:,JP_AER_NO3)  = 63.
XMI(:,:,:,JP_AER_NH3)  = 17.
XMI(:,:,:,JP_AER_H2O)  = 18.
XMI(:,:,:,JP_AER_BC) = 12.
XMI(:,:,:,JP_AER_DST)  = 100.
IF (NSOA .EQ. 10) THEN
XMI(:,:,:,JP_AER_SOA1) = 88. 
XMI(:,:,:,JP_AER_SOA2) = 180.
XMI(:,:,:,JP_AER_SOA3) = 1.5374857E+02
XMI(:,:,:,JP_AER_SOA4) = 1.9586780E+02
XMI(:,:,:,JP_AER_SOA5) = 195.
XMI(:,:,:,JP_AER_SOA6) = 195.
XMI(:,:,:,JP_AER_SOA7) = 165.
XMI(:,:,:,JP_AER_SOA8) = 195.
XMI(:,:,:,JP_AER_SOA9) = 270.
XMI(:,:,:,JP_AER_SOA10) = 210.
END IF


! Moments index
NM0(1) = 1 
NM3(1) = 2
NM6(1) = 3 
NM0(2) = 4 
NM3(2) = 5 
NM6(2) = 6 

IF (CRGUNIT=="MASS") THEN
  ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
  ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
  ZINIRADIUSI = XINIRADIUSI 
  ZINIRADIUSJ = XINIRADIUSJ
END IF
ZMINRGI = ZINIRADIUSI * XCOEFRADIMIN
ZMINRGJ = ZINIRADIUSJ * XCOEFRADJMIN


! Aerosol Density
! Cf Ackermann (all to black carbon except water)
XRHOI(:) = 1.8e3
XRHOI(JP_AER_H2O) = 1.0e3   ! water
XRHOI(JP_AER_DST) = XDENSITY_DUST  ! dusts

PCHEM(:,:,:,:) = MAX(PCHEM(:,:,:,:), 0.)
PAERO(:,:,:,:) = MAX(PAERO(:,:,:,:), XMNH_TINY )
!

DO JJ=1,NSP+NCARB+NSOA
  XFAC(JJ)=(4./3.)*3.14292654*XRHOI(JJ)*1.e-9
ENDDO
!
!
!*       1.n    transfer aerosol mass from gas to aerosol variables
!               (and conversion of part/part --> microgram/m3)
!
DO JJ=1,NSV_AER
  PAERO(:,:,:,JJ) =  PAERO(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
ENDDO

!
PCTOTA(:,:,:,:,:) =0.
! mineral phase
  PCTOTA(:,:,:,JP_AER_SO4,1) = PAERO(:,:,:,JP_CH_SO4i)*XMI(:,:,:,JP_AER_SO4)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SO4,2) = PAERO(:,:,:,JP_CH_SO4j)*XMI(:,:,:,JP_AER_SO4)/6.0221367E+11

  PCTOTA(:,:,:,JP_AER_NO3,1) = PAERO(:,:,:,JP_CH_NO3i)*XMI(:,:,:,JP_AER_NO3)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_NO3,2) = PAERO(:,:,:,JP_CH_NO3j)*XMI(:,:,:,JP_AER_NO3)/6.0221367E+11

  PCTOTA(:,:,:,JP_AER_NH3,1) = PAERO(:,:,:,JP_CH_NH3i)*XMI(:,:,:,JP_AER_NH3)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_NH3,2) = PAERO(:,:,:,JP_CH_NH3j)*XMI(:,:,:,JP_AER_NH3)/6.0221367E+11

! water
  PCTOTA(:,:,:,JP_AER_H2O,1) = PAERO(:,:,:,JP_CH_H2Oi)*XMI(:,:,:,JP_AER_H2O)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_H2O,2) = PAERO(:,:,:,JP_CH_H2Oj)*XMI(:,:,:,JP_AER_H2O)/6.0221367E+11

!
! primary organic carbon
  PCTOTA(:,:,:,JP_AER_OC,1) = PAERO(:,:,:,JP_CH_OCi)*XMI(:,:,:,JP_AER_OC)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_OC,2) = PAERO(:,:,:,JP_CH_OCj)*XMI(:,:,:,JP_AER_OC)/6.0221367E+11

! primary black carbon
  PCTOTA(:,:,:,JP_AER_BC,1) = PAERO(:,:,:,JP_CH_BCi)*XMI(:,:,:,JP_AER_BC)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_BC,2) = PAERO(:,:,:,JP_CH_BCj)*XMI(:,:,:,JP_AER_BC)/6.0221367E+11

!dust
  PCTOTA(:,:,:,JP_AER_DST,1) = PAERO(:,:,:,JP_CH_DSTi)*XMI(:,:,:,JP_AER_DST)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_DST,2) = PAERO(:,:,:,JP_CH_DSTj)*XMI(:,:,:,JP_AER_DST)/6.0221367E+11


!
IF (NSOA .EQ. 10) THEN
  PCTOTA(:,:,:,JP_AER_SOA1,1) = PAERO(:,:,:,JP_CH_SOA1i)*XMI(:,:,:,JP_AER_SOA1)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA1,2) = PAERO(:,:,:,JP_CH_SOA1j)*XMI(:,:,:,JP_AER_SOA1)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA2,1) = PAERO(:,:,:,JP_CH_SOA2i)*XMI(:,:,:,JP_AER_SOA2)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA2,2) = PAERO(:,:,:,JP_CH_SOA2j)*XMI(:,:,:,JP_AER_SOA2)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA3,1) = PAERO(:,:,:,JP_CH_SOA3i)*XMI(:,:,:,JP_AER_SOA3)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA3,2) = PAERO(:,:,:,JP_CH_SOA3j)*XMI(:,:,:,JP_AER_SOA3)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA4,1) = PAERO(:,:,:,JP_CH_SOA4i)*XMI(:,:,:,JP_AER_SOA4)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA4,2) = PAERO(:,:,:,JP_CH_SOA4j)*XMI(:,:,:,JP_AER_SOA4)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA5,1) = PAERO(:,:,:,JP_CH_SOA5i)*XMI(:,:,:,JP_AER_SOA5)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA5,2) = PAERO(:,:,:,JP_CH_SOA5j)*XMI(:,:,:,JP_AER_SOA5)/6.0221367E+11

  PCTOTA(:,:,:,JP_AER_SOA6,1) = PAERO(:,:,:,JP_CH_SOA6i)*XMI(:,:,:,JP_AER_SOA6)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA6,2) = PAERO(:,:,:,JP_CH_SOA6j)*XMI(:,:,:,JP_AER_SOA6)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA7,1) = PAERO(:,:,:,JP_CH_SOA7i)*XMI(:,:,:,JP_AER_SOA7)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA7,2) = PAERO(:,:,:,JP_CH_SOA7j)*XMI(:,:,:,JP_AER_SOA7)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA8,1) = PAERO(:,:,:,JP_CH_SOA8i)*XMI(:,:,:,JP_AER_SOA8)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA8,2) = PAERO(:,:,:,JP_CH_SOA8j)*XMI(:,:,:,JP_AER_SOA8)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA9,1) = PAERO(:,:,:,JP_CH_SOA9i)*XMI(:,:,:,JP_AER_SOA9)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA9,2) = PAERO(:,:,:,JP_CH_SOA9j)*XMI(:,:,:,JP_AER_SOA9)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA10,1) = PAERO(:,:,:,JP_CH_SOA10i)*XMI(:,:,:,JP_AER_SOA10)/6.0221367E+11
  PCTOTA(:,:,:,JP_AER_SOA10,2) = PAERO(:,:,:,JP_CH_SOA10j)*XMI(:,:,:,JP_AER_SOA10)/6.0221367E+11
END IF


!*       1.1    calculate moment 3 from mass
    
PM3D(:,:,:,2) = 0.
PM3D(:,:,:,5) = 0.
PCTOTA(:,:,:,:,:) = MAX(PCTOTA(:,:,:,:,:), 0.)
DO JJ = 1,NSP+NCARB+NSOA
  PM3D(:,:,:,2) = PM3D(:,:,:,2)+PCTOTA(:,:,:,JJ,1)/XFAC(JJ)
  PM3D(:,:,:,5) = PM3D(:,:,:,5)+PCTOTA(:,:,:,JJ,2)/XFAC(JJ)
ENDDO
!
!
!
IF ((CCONF=="START").AND.(CPROGRAM/='DIAG  ')) THEN
!*       1.2    calculate moment 0 from dispersion and mean radius
    PM3D(:,:,:,1)= PM3D(:,:,:,2) / &
               ((ZINIRADIUSI**3)*EXP(4.5 * (LOG(XINISIGI))**2))

    PM3D(:,:,:,4)= PM3D(:,:,:,5) / &
               ((ZINIRADIUSJ**3)*EXP(4.5 * (LOG(XINISIGJ))**2))

!*       1.3    calculate moment 6 from dispersion and mean radius
   PM3D(:,:,:,3) = PM3D(:,:,:,1) * (ZINIRADIUSI**6) *EXP(18 *(LOG(XINISIGI))**2)
   PM3D(:,:,:,6) = PM3D(:,:,:,4) * (ZINIRADIUSJ**6) *EXP(18 *(LOG(XINISIGJ))**2)

ELSE
!*       1.2    give  moment 0 
   PM3D(:,:,:,1)= MAX(PAERO(:,:,:,JP_CH_M0i) * 1E+6 , 0.)
   PM3D(:,:,:,4)= MAX(PAERO(:,:,:,JP_CH_M0j) * 1E+6 , 0.)
!
!*       1.3    give  moment 6 
   
IF (LVARSIGI) THEN ! set M6 variable standard deviation
 PM3D(:,:,:,3) = MAX(PAERO(:,:,:,JP_CH_M6i), XMNH_TINY)
ELSE ! fixed standard deviation
 PM3D(:,:,:,3) = PM3D(:,:,:,1) &
          * ( (PM3D(:,:,:,2)/PM3D(:,:,:,1))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGI)**2))**6 &
          * exp(18.*log(XINISIGI)**2)
END IF

IF (LVARSIGJ) THEN ! set M6 variable standard deviation
 PM3D(:,:,:,6) = MAX(PAERO(:,:,:,JP_CH_M6j), XMNH_TINY)
ELSE ! fixed standard deviation
 PM3D(:,:,:,6) = PM3D(:,:,:,4) &
          * ( (PM3D(:,:,:,5)/PM3D(:,:,:,4))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGJ)**2))**6 &
          * exp(18.*log(XINISIGJ)**2)
END IF

!
!
ENDIF
!
!**********************************************
! Calcul de XRHOP3D
!**********************************************

PRHOP3D(:,:,:,:)=0.
DO JN=1,JPMODE
  ZSUM(:,:,:)=0.
  DO JJ=1,NSP+NCARB+NSOA
   ZSUM(:,:,:)=ZSUM(:,:,:)+PCTOTA(:,:,:,JJ,JN)/XRHOI(JJ)
  ENDDO
  DO JJ=1,NSP+NCARB+NSOA
  ZCCTOT(:,:,:,JJ,JN)=PCTOTA(:,:,:,JJ,JN)/XRHOI(JJ)/ZSUM(:,:,:)
  PRHOP3D(:,:,:,JN)=PRHOP3D(:,:,:,JN)+ZCCTOT(:,:,:,JJ,JN)*XRHOI(JJ)
  ENDDO
ENDDO


DO JN=1,JPMODE
  IF (JN .EQ. 1) THEN

    IF (LVARSIGI) THEN ! variable dispersion for mode 1

      ZSIGMA(:,:,:)=PM3D(:,:,:,NM3(JN))**2/(PM3D(:,:,:,NM0(JN))*PM3D(:,:,:,NM6(JN)))
      ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
      WHERE (ZSIGMA(:,:,:) > XSIGIMAX)
      ZSIGMA(:,:,:) =  XSIGIMAX
      END WHERE
      WHERE (ZSIGMA(:,:,:) < XSIGIMIN)
      ZSIGMA(:,:,:) =  XSIGIMIN
      END WHERE

    ELSE ! fixed dispersion for mode 1
      ZSIGMA(:,:,:) = XINISIGI
    END IF
  END IF

!
  IF (JN .EQ. 2) THEN

    IF (LVARSIGJ) THEN ! variable dispersion for mode 2

      ZSIGMA(:,:,:)=PM3D(:,:,:,NM3(JN))**2/(PM3D(:,:,:,NM0(JN))*PM3D(:,:,:,NM6(JN)))
      ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
      ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
      WHERE (ZSIGMA(:,:,:) > XSIGJMAX)
      ZSIGMA(:,:,:) =  XSIGJMAX
      END WHERE
      WHERE (ZSIGMA(:,:,:) < XSIGJMIN)
      ZSIGMA(:,:,:) =  XSIGJMIN
      END WHERE

    ELSE ! fixed dispersion for mode 2
      ZSIGMA(:,:,:) = XINISIGJ
    END IF
  END IF
!


!*       1.4    calculate modal parameters from moments
PSIG3D(:,:,:,JN) = ZSIGMA(:,:,:)
PN3D(:,:,:,JN) = PM3D(:,:,:,NM0(JN))

ZSIGMA(:,:,:)=LOG(PSIG3D(:,:,:,JN))**2

PRG3D(:,:,:,JN)=(PM3D(:,:,:,NM3(JN))/PN3D(:,:,:,JN))**(1./3.)*EXP(-1.5*ZSIGMA(:,:,:))

PM3D(:,:,:,NM6(JN))=PN3D(:,:,:,JN)*PRG3D(:,:,:,JN)**6*EXP(18.*ZSIGMA(:,:,:))
!
ENDDO
!
!
PAERO(:,:,:,JP_CH_M0i) = PM3D(:,:,:,1) * 1E-6 
PAERO(:,:,:,JP_CH_M0j) = PM3D(:,:,:,4) * 1E-6
IF (LVARSIGI) PAERO(:,:,:,JP_CH_M6i) = PM3D(:,:,:,3) 
IF (LVARSIGJ) PAERO(:,:,:,JP_CH_M6j) = PM3D(:,:,:,6)

!
  DO JJ=1,NSV_AER
  PAERO(:,:,:,JJ) =  PAERO(:,:,:,JJ) / (ZDEN2MOL * PRHODREF(:,:,:))
  ENDDO

XSEDA(:,:,:,:)=0. ! no sedimentation for the first time step
  
!*      0.3    définition of minimum values
!Minimum values for gaseous (interact with aerosol phase) 
XSVMIN(NSV_AERBEG:NSV_AEREND) = XMNH_TINY
XSVMIN(NSV_CHEMBEG-1+JP_CH_CO) =  1E-10
! For i mode
ZRHODREFMIN = MAX_ll( PRHODREF(:,:,:), IINFO_ll)
ZMASS  = XN0IMIN *  ((ZMINRGI**3)*EXP(4.5 * (LOG(XSIGIMIN))**2))
ZM6MIN = XN0IMIN *  ((ZMINRGI**6)*EXP(18. * (LOG(XSIGIMIN))**2))
XSVMIN(NSV_AERBEG-1+JP_CH_BCi) = ZMASS * XFAC(JP_AER_BC) * 6.0221367E+11/(ZDEN2MOL*12.*ZRHODREFMIN)
XSVMIN(NSV_AERBEG-1+JP_CH_DSTi) = ZMASS * XFAC(JP_AER_DST) * 6.0221367E+11/(ZDEN2MOL*12.*ZRHODREFMIN)
XSVMIN(NSV_AERBEG-1+JP_CH_M0i) = XN0IMIN * 1E-6 / (ZDEN2MOL*ZRHODREFMIN)
IF (LVARSIGI) XSVMIN(NSV_AERBEG-1+JP_CH_M6i) = ZM6MIN  / (ZDEN2MOL*ZRHODREFMIN)
!
! For j mode
ZMASS  = XN0JMIN *  ((ZMINRGJ**3)*EXP(4.5 * (LOG(XSIGJMIN))**2))
ZM6MIN = XN0JMIN *  ((ZMINRGJ**6)*EXP(18. * (LOG(XSIGJMIN))**2))
XSVMIN(NSV_AERBEG-1+JP_CH_BCj) = ZMASS * XFAC(JP_AER_BC) * 6.0221367E+11/(ZDEN2MOL*12.*ZRHODREFMIN)
XSVMIN(NSV_AERBEG-1+JP_CH_DSTj) = ZMASS * XFAC(JP_AER_DST) * 6.0221367E+11/(ZDEN2MOL*12.*ZRHODREFMIN)
XSVMIN(NSV_AERBEG-1+JP_CH_M0j) = XN0JMIN * 1E-6 / (ZDEN2MOL*ZRHODREFMIN)
IF (LVARSIGJ) XSVMIN(NSV_AERBEG-1+JP_CH_M6j) = ZM6MIN  / (ZDEN2MOL*ZRHODREFMIN)
!
XSVMIN(NSV_AERBEG:NSV_AEREND) = 0.
!
END SUBROUTINE CH_AER_EQM_INIT_n
