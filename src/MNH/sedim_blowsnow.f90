!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!   ##############################
     MODULE MODI_SEDIM_BLOWSNOW
!!   ##############################
!!
INTERFACE
!
SUBROUTINE SEDIM_BLOWSNOW(  &
     PTHT               & !I [K] theta
     ,PDTMONITOR        & !I Time step
     ,PRHODREF          & !I [kg/m3] air density
     ,PPABST            & !I [Pa] pressure
     ,PZZ               & !I [m] height of layers
     ,PSVT              & !IO [scalar variable, ppp] Blowing snow concentration
     ,PSVS              & !IO ! Blowing snow variable source
     ,PVGK              &  !I [m/s] Blowing snow variable settling velocity 
     )

IMPLICIT NONE

REAL,  INTENT(IN) :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSVT   !scalar variable
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSVS   !scalar variable
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT)    :: PVGK   !Settling velocity of blowing snow variable



END SUBROUTINE SEDIM_BLOWSNOW
!!
END INTERFACE
!!
END MODULE MODI_SEDIM_BLOWSNOW

!!   #######################################
     SUBROUTINE SEDIM_BLOWSNOW(PTHT,PDTMONITOR,&
                       PRHODREF,PPABST,PZZ,PSVT,&
                       PSVS,PVGK)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   Based on sedim_dust.f90 from Pierre Tulet
!!
!!   AUTHOR
!!    ------
!!    Vincent Vionnet (GMME) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_BLOWSNOW
USE MODD_CSTS_BLOWSNOW
USE MODI_BLOWSNOW_VELGRAV
USE MODE_BLOWSNOW_PSD
!
!
IMPLICIT NONE
!*      0.1    declarations of arguments
!
REAL,  INTENT(IN) :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT)    :: PSVT  !scalar variable 
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT)    :: PSVS  !scalar variable
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT)    :: PVGK   !Settling velocity of blowing snow variable

!
!*      0.2    declarations of local variables
!
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   :: ZRG, ZBETA!,ZMOB
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NBLOWSNOW3D) :: ZPM, ZPMOLD
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)+1,NBLOWSNOW3D) :: ZFLUXSED, ZFLUXMAX
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZH, ZW, ZVSNUMMAX
REAL,  DIMENSION(2)   :: ZPMIN
REAL                  :: ZRGMIN,ZHMIN,ZTSPLITR,ZVSMAX,ZT
INTEGER               :: ILU  ! indice K End       in z direction
INTEGER               :: JK,JN,JT ! Loop counters
INTEGER               ::  ISPLITA 
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NBLOWSNOW3D) :: ZPSVS
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NBLOWSNOW3D) :: ZSVS
!
!*       0.3   initialize constant
!
ZH(:,:,:)      = 0.
ZW(:,:,:)      = 0.
ZFLUXSED(:,:,:,:) = 0.
ILU = SIZE(PSVT,3)
ZPSVS(:,:,:,:) = 0.0
!
!
!*       1.   compute dimensions of arrays
!
!Get minimum values possible
  ZPMIN(1) = XN0MIN_SNW
  ZRGMIN   = XINIRADIUS_SNW
  ZPMIN(2) = 4*XPI*XRHOLI/3*(ZRGMIN/XALPHA_SNOW)**3.*(XALPHA_SNOW+2)*(XALPHA_SNOW+1)*XALPHA_SNOW*XN0MIN_SNW
!
!
!*       2.   compute BETA, RG and moments using profile at time t (PSVT)
!
!
CALL PPP2SNOW(PSVT, PRHODREF, PBET3D=ZBETA, PRG3D=ZRG, PM3D=ZPM)! ,PMOB3D=ZMOB)

!
!*       3.  Mobility Index evolution
!
! Temporal decrease towards fine grain in 6 hours from fresh snow
!WHERE(PSVS(:,:,:,2)>0)
!        ZMOB(:,:,:) = PSVS(:,:,:,3)/PSVS(:,:,:,2)
!END WHERE
!ZMOB(:,:,:) = MAX(1.2,ZMOB(:,:,:)-PDTMONITOR/(3600*6))

! Update transported mobility index (M3*Mob)
!PSVS(:,:,:,3) = ZMOB(:,:,:)* PSVT(:,:,:,2)  / PDTMONITOR 
!PSVS(:,:,:,3) = ZMOB(:,:,:)* PSVS(:,:,:,2)  


!*       4.   Adjust settling velocity
!
! No sedimentation at fisrt atmospheric level since sedimentation for this level 
! is already included in the surface net flux
PVGK(:,:,1,:) = 0
!
!*       5.   Compute time-splitting condition
!
ZH=9999.
ZVSNUMMAX(:,:,:)   = 0.

DO JK=1,ILU
   ZH(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)
   ! Maximum velocity 
   ZVSNUMMAX(:,:,JK) = MIN(10. *ZH(:,:,JK) / PDTMONITOR,20.)
ENDDO
!
ZHMIN=MINVAL(ZH(:,:,1:ILU))

ZVSMAX = 2.      
ISPLITA = 1
SPLIT : DO
  ZT = PDTMONITOR / FLOAT(ISPLITA)
  IF ( ZT * ZVSMAX / ZHMIN .LT. 1.) EXIT SPLIT
  ISPLITA = ISPLITA + 1
END DO SPLIT
ZTSPLITR  = PDTMONITOR / FLOAT(ISPLITA)

ZFLUXSED(:,:,:,:) = 0.
ZFLUXMAX(:,:,:,:) = 0.


DO JN=1,NBLOWSNOW3D ! Compute sedimentation for both moments (M0 and M3) and mobility index
  !
  ZFLUXSED(:,:,ILU+1,JN) = 0.
!
! ZPSVS = Specie SV source creating during the current time step
! PSVS  = Source of the previous time step 
!  
  ZPSVS(:,:,:,JN) = PSVS(:,:,:,JN)-PSVT(:,:,:,JN)/PDTMONITOR
  PSVS(:,:,:,JN)  = PSVT(:,:,:,JN)/PDTMONITOR  

DO JT = 1 , ISPLITA

  IF( JT==1 ) THEN
    DO JK=1,ILU
       ZW(:,:,JK)  = ZTSPLITR /(PRHODREF(:,:,JK)*(PZZ(:,:,JK+1)-PZZ(:,:,JK)))
    END DO
    PSVS(:,:,:,JN) = MAX(0.,PSVS(:,:,:,JN) + ZPSVS(:,:,:,JN)/ISPLITA)
  ELSE
    PSVS(:,:,:,JN) = MAX(0.,PSVS(:,:,:,JN) + ZPSVS(:,:,:,JN)*ZTSPLITR)
  END IF

  IF( JT==1 ) PSVS(:,:,:,JN) = PSVS(:,:,:,JN) *PDTMONITOR

  ! Compute concentration averaged verticaly within one layer
  ZSVS(:,:,:,JN) =  PSVS(:,:,:,JN)*PRHODREF(:,:,1:ILU) 
  ! Compute sedimentation flux  F = C*V    [kg/m2/s]
  ZFLUXSED(:,:,1:ILU,JN)= PVGK(:,:,1:ILU,JN)* ZSVS(:,:,1:ILU,JN)

  DO JK=1,ILU
      PSVS(:,:,JK,JN)= MAX(0.,PSVS(:,:,JK,JN) + &
                      ZW(:,:,JK)*(ZFLUXSED(:,:,JK+1,JN)- ZFLUXSED(:,:,JK,JN)))                      
  END DO

  IF( JT==ISPLITA ) THEN
        PSVS(:,:,:,JN) = MAX(0.,PSVS(:,:,:,JN) / PDTMONITOR)
  END IF
END DO
END DO

END SUBROUTINE SEDIM_BLOWSNOW
