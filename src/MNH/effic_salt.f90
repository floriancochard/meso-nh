!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_EFFIC_SALT
!!   ##############################
!!
INTERFACE
!
SUBROUTINE EFFIC_SALT(  &
      PTHT              & !I [K] theta
     ,PRHODREF          & !I [kg/m3] air density
     ,PPABST            & !I [Pa] pressure
     ,PURR              & !I
     ,PSVT              & !I [scalar variable, ppp] sea salt concentration
     ,PEFFIC            & !O [scalar variable, ppp] sea salt concentration
     )

IMPLICIT NONE

REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PURR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT   !scalar variable 
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PEFFIC !scavenging efficiency


END SUBROUTINE EFFIC_SALT
!!
END INTERFACE
!!
END MODULE MODI_EFFIC_SALT
!!
!!   #######################################
     SUBROUTINE EFFIC_SALT(PTHT,PRHODREF,&
                       PPABST,PURR,PSVT,PEFFIC)
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
!!    Pierre TULET & Joris PIANEZZE (LACy)
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
! Entry variables:
!
! PSVTS(INOUT)       -Array of moments included in PSVTS
!
!*************************************************************
! Exit variables:
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! ZVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
!
!
USE MODD_CSTS_SALT
USE MODD_SALT
USE MODI_SALT_VELGRAV
USE MODI_AER_EFFIC3D
USE MODD_CST,         ONLY : XP00, XRD
USE MODD_PARAMETERS , ONLY : JPVEXT
USE MODE_SALT_PSD
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PURR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT   !scalar variable 
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PEFFIC !scavenging efficiency
!
!*      0.2    declarations of local variables
!
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NMODE_SLT)   :: ZRG, ZSIG, ZDENSITY_AER
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4)):: ZSVT

REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZMU,ZMUW, ZTEMP  
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),3*NMODE_SLT) :: ZVGK, ZDPK
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NMODE_SLT)   :: ZVG, ZDPG, ZCOR
REAL ::  ZRHOI
REAL ::  A0, A1, A2, A3     ! Constants for computing viscosity
!
!*       0.3   initialize constant
!
ZRHOI = XDENSITY_SALT
ZSVT(:,:,:,:) = PSVT(:,:,:,:)

CALL PPP2SALT(ZSVT, PRHODREF, PSIG3D=ZSIG, PRG3D=ZRG)
!
CALL SALT_VELGRAV(ZSIG, ZRG, PTHT, PPABST, PRHODREF, ZRHOI, &
                  ZMU, ZVGK, ZDPK, ZVG, ZDPG, PCOR=ZCOR, PTEMP=ZTEMP)
!
!calcul de ZMUW

A0=1.76
A1=-5.5721e-2
A2=-1.3943e-3
A3=-4.3015e-5
ZMUW(:,:,:)=A0*EXP(A1*(ZTEMP(:,:,:)-273.15) &
        +A2*(ZTEMP(:,:,:)-273.15) + A3*(ZTEMP(:,:,:)-273.15))*1.e-3

A1=-3.5254e-2
A2=4.7163e-4
A3=-6.0667e-6

 WHERE(ZTEMP(:,:,:)>273.15)
         ZMUW(:,:,:)=A0*EXP(A1*(ZTEMP(:,:,:)-273.15) &
        +A2*(ZTEMP(:,:,:)-273.15) + A3*(ZTEMP(:,:,:)-273.15))*1.e-3

 END WHERE
 ZMUW(:,:,:)=MAX(ZMUW(:,:,:),1.e-12)

ZDENSITY_AER(:,:,:,:) = XDENSITY_SALT

CALL AER_EFFIC3D(ZRG,ZVG,     &
                 PRHODREF,    &
                 ZMUW,  ZMU,  &
                 ZDPG,        &
                 PURR,        &
                 NMODE_SLT,   &
                 ZTEMP, ZCOR, &
                 ZDENSITY_AER,&
                 PEFFIC       )       

END SUBROUTINE EFFIC_SALT
