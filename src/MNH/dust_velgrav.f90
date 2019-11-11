!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_DUST_VELGRAV
!!   ##############################
!!
INTERFACE
!!
SUBROUTINE DUST_VELGRAV(PSIG, PRG, PTHT, PABST,PRHODREF, PRHOP,&
                        PMU, PVGK,PDPK, PVGG, PDPG, PCOR, PTEMP)
IMPLICIT NONE
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PSIG, PRG
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHT, PABST, PRHODREF
REAL,                     INTENT(IN)  :: PRHOP
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PMU
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PVGK,PDPK
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PVGG
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PDPG
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(INOUT) :: PCOR
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PTEMP
END SUBROUTINE DUST_VELGRAV
!!
END INTERFACE
!!
END MODULE MODI_DUST_VELGRAV
!!
!!   #######################################
SUBROUTINE DUST_VELGRAV(PSIG, PRG, PTHT, PABST,PRHODREF, PRHOP,&
                        PMU, PVGK,PDPK, PVGG, PDPG, PCOR, PTEMP)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!   P. Tulet (CNRM/GMEI)
!!
!!   MODIFICATIONS
!!    -------------
!!
! Entry variables:
!
! PM(IN)       -Array of moments
!
!*************************************************************
! Exit variables:
!
! PFSED(IN)  -Array of moment variation due to dry deposition
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! PDPK       -Polydisperse diffusivity (m2/s)
! PVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!-----------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DUST
USE MODD_CSTS_DUST
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XAVOGADRO        & ![molec/mol] avogadros number
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD               !  Cpd (dry air)
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PSIG, PRG
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHT, PABST, PRHODREF
REAL,                     INTENT(IN)  :: PRHOP
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PMU
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PVGK,PDPK
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PVGG
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PDPG
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(INOUT) :: PCOR
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PTEMP
!
!*      0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3)) :: ZLAMBDA, ZTEMP
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3)) :: ZRG,ZLN2S
!
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3)) :: ZKNG
REAL, DIMENSION(SIZE(PSIG,1),SIZE(PSIG,2),SIZE(PSIG,3),SIZE(PSIG,4)) :: ZCOR
!
!
REAL, PARAMETER :: gasmw=28.9644E-3
REAL, PARAMETER :: gasr=8.3143
REAL :: ZK
!
INTEGER :: JI,JJ
!
!-----------------------------------------------------------------
!temperature
ZTEMP(:,:,:)=PTHT(:,:,:)*(PABST(:,:,:)/XP00)**(XRD/XCPD)
!
! Sutherland's equation for viscosity
PMU(:,:,:)=1.8325d-5*416.16/(ZTEMP(:,:,:)+120)*(ZTEMP(:,:,:)/296.16)*SQRT(ZTEMP(:,:,:)/296.16)
!
! Mean free path (Seinfeld and Pandis p455)
ZLAMBDA(:,:,:)=2*PMU(:,:,:)/(PABST(:,:,:)*SQRT(8*gasmw/(XPI*gasr*ZTEMP(:,:,:))))
!
!
DO JI=1,NMODE_DST
  ZRG(:,:,:)=PRG(:,:,:,JI) * 1E-6 
  ZLN2S(:,:,:)=LOG(PSIG(:,:,:,JI))**2 
  !
  ZKNG(:,:,:)=ZLAMBDA(:,:,:) / PRG(:,:,:,JI) 
  !
  ! 
  !Slip Correction Factor
  ZCOR(:,:,:,JI)=1+ZLAMBDA(:,:,:)/ZRG(:,:,:)*(1.257+ &
              0.4*exp(-1.1*ZRG(:,:,:)/ZLAMBDA(:,:,:)))
  ZCOR(:,:,:,JI)=MAX(ZCOR(:,:,:,JI),1.0)

  PVGG(:,:,:,JI)= 2.*XG*PRHOP*ZRG(:,:,:)**2 /(9.*PMU(:,:,:))*ZCOR(:,:,:,JI)
  PDPG(:,:,:,JI)=XBOLTZ*ZTEMP(:,:,:)/ (6.*XPI* ZRG(:,:,:)*PMU(:,:,:))
  IF(PRESENT(PCOR)) PCOR(:,:,:,JI) = ZCOR(:,:,:,JI)

  IF(PRESENT(PTEMP)) PTEMP(:,:,:) = ZTEMP(:,:,:)
  ! 
  DO JJ=0,2
    ZK=REAL(3*JJ)
    PDPK(:,:,:,3*JI+JJ-2)=PDPG(:,:,:,JI)*(exp((-2.*ZK+1.)/2.*ZLN2S(:,:,:))+1.246*ZKNG(:,:,:)*&
              exp((-4.*ZK+4)/2.*ZLN2S(:,:,:)))
    !
    PVGK(:,:,:,3*JI+JJ-2)=PVGG(:,:,:,JI)*&
    (exp((4.*ZK+4.)/2.*ZLN2S(:,:,:)) + 1.246*ZKNG(:,:,:)* exp((2.*ZK+1.)/2.*ZLN2S(:,:,:)))
  ENDDO
  !
ENDDO
!
!
END SUBROUTINE DUST_VELGRAV
