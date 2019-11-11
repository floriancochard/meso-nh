!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_CH_AER_VELGRAV_n
!!   ##############################
!!
INTERFACE
!!
SUBROUTINE CH_AER_VELGRAV_n(PSIG, PRG, PTHT, PABST,PRHODREF, PRHOP, PMU, PVGK,PDPK, PVGG)
IMPLICIT NONE
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGK,PDPK
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PMU
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGG
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRHOP
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PSIG, PRG
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT, PABST, PRHODREF
END SUBROUTINE CH_AER_VELGRAV_n
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_VELGRAV_n
!!
!!   #######################################
SUBROUTINE CH_AER_VELGRAV_n(PSIG, PRG, PTHT, PABST, PRHODREF, PRHOP, PMU, PVGK,PDPK, PVGG)
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
!!   P. Tulet (meteo france)
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
!!   IMPLICIT ARGUMENTS
USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XAVOGADRO        & ![molec/mol] avogadros number
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD               !  Cpd (dry air)
!!
!-------------------------------------------------------------------------------
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGK,PDPK
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PMU
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGG
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRHOP
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSIG, PRG
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PABST, PRHODREF
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(size(PSIG,1),size(PSIG,2),size(PSIG,3)) :: ZTEMP,ZLAMBDA
REAL, DIMENSION(size(PSIG,1),size(PSIG,2),size(PSIG,3)) :: ZRG,ZLN2S
REAL, DIMENSION(size(PSIG,1),size(PSIG,2),size(PSIG,3)) :: ZKNG
REAL, DIMENSION(size(PSIG,1),size(PSIG,2),size(PSIG,3),JPMODE) :: ZDPG
!
REAL, PARAMETER :: gasmw=28.9644d0
REAL :: ZK
!
INTEGER :: JI,JJ
!
!-------------------------------------------------------------------------------
!temperature
ZTEMP(:,:,:)=PTHT(:,:,:)*(PABST(:,:,:)/XP00)**(XRD/XCPD)
!
! Sutherland's equation for viscosity
PMU(:,:,:)=1.8325d-5*416.16/(ZTEMP(:,:,:)+120)*(ZTEMP(:,:,:)/296.16)*SQRT(ZTEMP(:,:,:)/296.16)
!PMU(:,:,:)=0.003661*ZTEMP(:,:,:)
!PMU(:,:,:)=.0066164*PMU(:,:,:)*sqrt(PMU(:,:,:))/(ZTEMP(:,:,:)+114.d0)
!
! Mean free path (Seinfeld and Pandis p455)
ZLAMBDA(:,:,:)=PMU(:,:,:)/PRHODREF(:,:,:)*sqrt(1.89d-4*gasmw/ZTEMP(:,:,:))*1.e6
!
DO JI=1,JPMODE
  !
  ZRG(:,:,:)=PRG(:,:,:,JI) * 1E-6 
  ZLN2S(:,:,:)=LOG(PSIG(:,:,:,JI))**2 
  !
  ZKNG(:,:,:)=ZLAMBDA(:,:,:) / PRG(:,:,:,JI) 
  !
  PVGG(:,:,:,JI)= 2.*XG*PRHOP(:,:,:,JI)*ZRG(:,:,:)**2 /(9.*PMU(:,:,:))
  ZDPG(:,:,:,JI)=XBOLTZ*ZTEMP(:,:,:)/ (6.*XPI* ZRG(:,:,:)*PMU(:,:,:))
  !
  DO JJ=0,2
    !
    ZK=real(3*JJ)
    PDPK(:,:,:,3*JI+JJ-2)=ZDPG(:,:,:,JI)*(exp((-2.*ZK+1.)/2.*ZLN2S(:,:,:))+1.246*ZKNG(:,:,:)*&
              exp((-4.*ZK+4)/2.*ZLN2S(:,:,:)))

    PVGK(:,:,:,3*JI+JJ-2)=PVGG(:,:,:,JI)*&
    (exp((4.*ZK+4.)/2.*ZLN2S(:,:,:)) + 1.246*ZKNG(:,:,:)* exp((2.*ZK+1.)/2.*ZLN2S(:,:,:)))
  !
  ENDDO
  !
ENDDO
!
END SUBROUTINE CH_AER_VELGRAV_n
