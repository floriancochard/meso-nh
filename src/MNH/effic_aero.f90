!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_EFFIC_AERO
!!   ##############################
!!
INTERFACE
!
SUBROUTINE EFFIC_AERO(  &
     PTHT               & !I [K] theta
     ,PRHODREF          & !I [kg/m3] air density
     ,PPABST            & !I [Pa] pressure
     ,PURR              & !I
     ,PSVT              & !I [scalar variable, ppp] sea salt concentration
     ,PEFFIC_AER            & !O [scalar variable, ppp] sea salt concentration
     )

IMPLICIT NONE

REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PURR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT   !scalar variable 
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PEFFIC_AER !scavenging efficiency


END SUBROUTINE EFFIC_AERO
!!
END INTERFACE
!!
END MODULE MODI_EFFIC_AERO
!!
!!   #######################################
     SUBROUTINE EFFIC_AERO(PTHT,PRHODREF,&
                       PPABST,PURR,PSVT,PEFFIC_AER)
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
!!    Pierre TULET (GMEI) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
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
USE MODI_AERO_VELGRAV
USE MODI_AERO_EFFIC3D
USE MODD_CST,         ONLY : XP00, XRD
USE MODD_PARAMETERS , ONLY : JPVEXT
USE MODE_AERO_PSD
USE MODD_CH_AEROSOL,  ONLY: JPMODE, XRHOI

!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PURR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT   !scalar variable 
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PEFFIC_AER !scavenging efficiency
!
!*      0.2    declarations of local variables
!
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZRG, ZSIG, ZDENSITY_AER
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),SIZE(PSVT,4)):: ZSVT
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NSP+NCARB+NSOA,JPMODE):: ZCTOTA

REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZMU,ZMUW, ZTEMP  
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),3*JPMODE) :: ZVGK, ZDPK
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZVG, ZDPG, ZCOR, ZRHOI, ZSUM1, ZSUM2
REAL ::  A0, A1, A2, A3     ! Constants for computing viscosity
INTEGER :: JSV, II
!
!*       0.3   initialize constant
!

ZSVT(:,:,:,:) = PSVT(:,:,:,:)

CALL PPP2AERO(ZSVT, PRHODREF, PSIG3D=ZSIG, PRG3D=ZRG, PCTOTA=ZCTOTA)

!
ZSUM1(:,:,:,:) = 0.
ZSUM2(:,:,:,:) = 0.
DO II=1,JPMODE
 DO JSV=1,NSP+NCARB+NSOA
   ZSUM1(:,:,:,II) = ZSUM1(:,:,:,II) + ZCTOTA(:,:,:,JSV,II)*XRHOI(JSV)
   ZSUM2(:,:,:,II) = ZSUM2(:,:,:,II) + ZCTOTA(:,:,:,JSV,II)
 ENDDO
ZRHOI(:,:,:,II) = ZSUM1(:,:,:,II)/ZSUM2(:,:,:,II) 
ENDDO
CALL AERO_VELGRAV(ZSIG, ZRG, PTHT, PPABST, PRHODREF, ZRHOI, &
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


CALL AERO_EFFIC3D(ZRG,ZVG,     &
                 PRHODREF,    &
                 ZMUW,  ZMU,  &
                 ZDPG,        &
                 PURR,        &
                 JPMODE,   &
                 ZTEMP, ZCOR, &
                 ZRHOI,&
                 PEFFIC_AER    )       

END SUBROUTINE EFFIC_AERO
