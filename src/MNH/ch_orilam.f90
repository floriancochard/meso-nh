!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/ch_orilam.f90,v $ $Revision: 1.1.2.1.2.1.18.1 $
! MASDEV4_7 chimie 2007/03/02 13:59:37
!-----------------------------------------------------------------
!!   #########################
     MODULE MODI_CH_ORILAM
!!   ######################### 
!!
INTERFACE
!!
SUBROUTINE CH_ORILAM(PAERO, PCHEM, PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT, PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT,  &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PFRAC, PMI,&
                           PTIME, GSCHEME, PSOLORG)
IMPLICIT NONE
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCHEM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PAERO
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PFRAC
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMI
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
CHARACTER(LEN=10),                      INTENT(IN)    :: GSCHEME


END SUBROUTINE CH_ORILAM
!!
END INTERFACE
!!
END MODULE MODI_CH_ORILAM
!!
!! #####################################################################################
SUBROUTINE CH_ORILAM(PAERO, PCHEM, PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT, PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT,  &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PFRAC, PMI,&
                           PTIME, GSCHEME, PSOLORG)
!! #####################################################################################
!!
!!    PURPOSE
!!    -------
!!
!!    ORILAM aerosol Code
!!    
!!
!!    Inputs:
!!    PCHEM : Chemical (gaseous and aerosol) species (in molec./cm3)  
!!    PSEDA : Moments 
!!
!!    Outputs:
!!
!!
!!
!!    REFERENCE
!!    ---------
!!    P. Tulet, V. Crassier, F. Cousin, K. Suhre, R. Rosset, jgr
!!    ORILAM, A three moment lognormal aerosol scheme for mesoscale atmospheric 
!!    model.
!!    On-line coupling into the Meso-NH-C model and validation  on the Escompte 
!!    campaign.
!!
!!    AUTHOR
!!    ------
!!    Pierre Tulet (GMEI) and Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!    M. Leriche (08/16) add initialization of ZMASK
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_AER_TRANS
USE MODI_CH_AER_DRIVER

!!
!!  IMPLICIT ARGUMENTS
!!  ------------------
!!
USE MODD_CH_AEROSOL
!
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCHEM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PAERO
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PFRAC
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMI
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
CHARACTER(LEN=10),                      INTENT(IN)    :: GSCHEME

REAL, DIMENSION(SIZE(PAERO,1),JPMODE)                 :: ZMASK
!
!-------------------------------------------------------------------------------
!initialize ZMASK
ZMASK(:,:) = 1.
!
! transfer gas phase variables into aerosol variables
CALL CH_AER_TRANS(0, PM, PSIG0, PRG0, PN0, PRHOP0,PAERO, PCHEM, PCTOTG, PCTOTA, PCCTOT,&
                         PFRAC, PMI, ZMASK,GSCHEME)

! integrate aerosol variables
CALL CH_AER_DRIVER(PM,PSIG0, PRG0, PN0,  PCTOTG, PCTOTA, PCCTOT,         &
                      PDTACT, PSEDA, PMU, PLAMBDA, PRHOP0, POM, PSO4RAT, &
                      PRV, PDENAIR, PPRESSURE, PTEMP, PRC, ZMASK, PTIME, &
                      PSOLORG)
!
! transfer aerosol variables back into gas phase variables
 CALL CH_AER_TRANS(1, PM, PSIG0, PRG0, PN0, PRHOP0, PAERO, PCHEM, PCTOTG, PCTOTA, PCCTOT,&
                      PFRAC, PMI, ZMASK,GSCHEME)
!
!
END SUBROUTINE CH_ORILAM
