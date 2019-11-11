!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!    #########################
     MODULE MODI_CH_AER_DRIVER
!    ######################### 
! 
INTERFACE
! 
SUBROUTINE CH_AER_DRIVER(PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT,  PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT, &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PMASK,&
                           PTIME, PSOLORG)
IMPLICIT NONE
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMASK
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
END SUBROUTINE CH_AER_DRIVER
! 
END INTERFACE
! 
END MODULE MODI_CH_AER_DRIVER
! 
!#####################################################################################
SUBROUTINE CH_AER_DRIVER(PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT, PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT,  &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PMASK,&
                           PTIME, PSOLORG)
!#####################################################################################
!!
!!    PURPOSE
!!    -------
!!
!!    compute the right hand side of the moment equations
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!       M.Leriche 2015 Calcul de la fraction massique entre les modes
!!       M.Leriche 08/16 suppress moments index declaration already in modd_aerosol
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_AER_COAG
USE MODI_CH_AER_GROWTH
USE MODI_CH_AER_SOLV
!!
!!  IMPLICIT ARGUMENTS
!!  ------------------
!!
USE MODD_CH_AEROSOL
!
!
IMPLICIT NONE
!  Declaration arguments
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMASK
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
!
!  Declarations variables internes
!
INTEGER                       :: II, JI, JJ

! Variables utilisees pour le tranfert de moment de chaque espece
! pour la condensation 
!----------------------------------------------------------------

REAL, DIMENSION(SIZE(PM,1),(JPMODE)*3)  :: ZDMINTRA,ZDMINTER,ZDMCOND

REAL                          :: ZGASMW       ! Molecular weight of background
                                           ! gas (g/mol) 
REAL, DIMENSION(SIZE(PM,1))   :: ZPGAS        ! background gas pressure (Pa)
REAL, DIMENSION(SIZE(PM,1))   :: ZRH,PSAT            ! Relative humidity
REAL                          :: ZDT          ! Pas de temps
REAL, DIMENSION(SIZE(PM,1))   :: ZPKM, ZPKH2O, ZSUM

!-----------------------------------------------------------------------------

ZDT=PDTACT

!*************************************************************
! Calcul de la fraction massique entre les modes
!*************************************************************
ZSUM (:) = 0.
DO JI=1,JPMODE
  DO JJ=1,NSP+NCARB+NSOA
    ZSUM (:) = ZSUM (:) + PCTOTA(:,JJ,JI)
  ENDDO
ENDDO
POM(:,:) = 0.
DO JI=1,JPMODE
  DO JJ=1,NSP+NCARB+NSOA
    POM(:,JI)  =  POM(:,JI) + PCTOTA(:,JJ,JI) / ZSUM (:)
  ENDDO
ENDDO


!******************************************************
!      Thermodynamic variables initialization
!         from Meso-NHC
!******************************************************
  
ZPKM(:) = 1E-3*PDENAIR(:) * 6.0221367E+23 / 28.9644
ZPKH2O(:) = ZPKM(:)*1.6077*PRV(:)
PSAT(:)=0.611*EXP(17.2694*(PTEMP(:)-273.16)/(PTEMP(:)-35.86))
PSAT(:)=PSAT(:)*1000.
ZRH(:)=(ZPKH2O(:)/(ZPKM(:)*1.6077))*PPRESSURE(:)/&
      &(0.622+(ZPKH2O(:)/(ZPKM(:)*1.6077)))/PSAT(:)

ZPGAS(:)=PPRESSURE(:)
ZGASMW=29.

!******************************************************
!      calculate gas viscosity and mean free path
!******************************************************
PMU(:)=0.003661*PTEMP(:)
PMU(:)=.0066164*PMU(:)*sqrt(PMU(:))/(PTEMP(:)+114.d0)

PLAMBDA(:)=PMU(:)/PDENAIR(:)*sqrt(1.89d-4*ZGASMW/PTEMP(:))*1.e6

CALL CH_AER_COAG(PM, PSIG0, PRG0, PN0,ZDMINTRA,ZDMINTER,&
                 PTEMP,PMU,PLAMBDA,PRHOP0)


CALL CH_AER_GROWTH(PM, PSIG0, PRG0, ZDMCOND,PDENAIR,ZGASMW,&
                     ZPGAS,PTEMP,ZRH,POM,PSO4RAT,PDTACT)

DO II=1,JPMODE
ZDMINTRA(:,NM0(II)) = ZDMINTRA(:,NM0(II)) * PMASK(:,II)
ZDMINTRA(:,NM3(II)) = ZDMINTRA(:,NM3(II)) * PMASK(:,II)
ZDMINTRA(:,NM6(II)) = ZDMINTRA(:,NM6(II)) * PMASK(:,II)
ZDMINTER(:,NM0(II)) = ZDMINTER(:,NM0(II)) * PMASK(:,II)
ZDMINTER(:,NM3(II)) = ZDMINTER(:,NM3(II)) * PMASK(:,II)
ZDMINTER(:,NM6(II)) = ZDMINTER(:,NM6(II)) * PMASK(:,II)
ZDMCOND(:,NM0(II)) = ZDMCOND(:,NM0(II)) * PMASK(:,II)
ZDMCOND(:,NM3(II)) = ZDMCOND(:,NM3(II)) * PMASK(:,II)
ZDMCOND(:,NM6(II)) = ZDMCOND(:,NM6(II)) * PMASK(:,II)
POM(:,II)          = POM(:,II)  * PMASK(:,II)
PSEDA(:,NM0(II))   = PSEDA(:,NM0(II))   * PMASK(:,II)
PSEDA(:,NM3(II))   = PSEDA(:,NM3(II))   * PMASK(:,II)
PSEDA(:,NM6(II))   = PSEDA(:,NM6(II))   * PMASK(:,II)
END DO


CALL CH_AER_SOLV(PM, PSIG0, PRG0, PN0,PCTOTG, PCTOTA, PCCTOT, &
                   ZDMINTRA,ZDMINTER,ZDMCOND,PSEDA,ZDT,POM,&
                   PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PTIME, PSOLORG)


END SUBROUTINE CH_AER_DRIVER
