!SFX_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!!   #######################################
SUBROUTINE BLOWSNW_VELGRAV1D(PBETA, PRG, PTA, PRHODREF,PPABST,PVGK)

!!   -------
!!
!!   Compute number- and mass-averaged settling velocity for
!!   blowing snow particles based on several methods :
!!         - Mitchell (1996) : numerical integration assuming spherical
!!                             particles (expensive)
!!         - Carrier (1953) and Dover (1993) : numerical integration (expensive)
!!         - look-up table based on Carrier (1953) depending on mean radius and
!!               pressure
!!         - None : assume no settling velocity
!!
!!   REFERENCE
!!   ---------
!!
!!       Mitchell (1996) : Use of mass- and area-dimensionla power laws for
!!             determining precipitation particle terminal velocities, JAS,
!!             53(12),1710-1723
!!       Carrier, C. : On Slow Viscous Flow, Tech. rep., Office of Naval Research, Contract Nonr-653(00), Brown
!!       University, Providence, RI, 1953.
!!       Dover, S. : Numerical Modelling of Blowing Snow, Ph.D. thesis, University of Leeds, U.K., 1993.
!!
!!   AUTHOR
!!    ------
!!   V. Vionnet (CNRM/GMME/MOSAYC)
!!
!!   MODIFICATIONS
!!    -------------
!!  Philippe Wautelet 28/05/2018: corrected truncated integer division (1*10**(-6) -> 1E-6)
!!
!!   NB : this routine is similar to the routine implemented in Meso-NH (blowsnow_velgrav.f90)
!
!!

USE MODD_BLOWSNW_SURF
USE MODD_CSTS

USE MODI_GAMMA_INC_LOW
USE MODI_GAMMA_SURF

USE MODE_BLOWSNW_SEDIM_LKT1D

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE
  !
  !*       0.1   Declarations of dummy arguments :
  !
       REAL, DIMENSION(:,:),   INTENT(IN)    :: PBETA        ! distribution scale parameter [m]
       REAL, DIMENSION(:,:),   INTENT(IN)    :: PRG          !  mean radius [m]
       REAL, DIMENSION(:,:),   INTENT(IN)       :: PTA       ! air temperature
       REAL, DIMENSION(:,:),   INTENT(IN)       :: PPABST     ! air pressure
       REAL, DIMENSION(:),     INTENT(IN)       :: PRHODREF   ! air density
       REAL, DIMENSION(:,:,:),     INTENT(OUT)  :: PVGK       !  terminal fallspeed [m/s]
  !
  !
  !*       0.2   Declarations of local variables :
  !   

      REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZMU                !air kinematic viscosity [m2/s] 
      REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZR1                ! first limit radius in integration
      REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZR2                ! second limit radius in integration [m]
      REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZAM1               ! Parameter avr in Mitchell
      REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZAM2               ! 
      REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZAM3               ! 

REAL, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: ZAA 
REAL, DIMENSION(SIZE(PTA,1))             :: ZBB    
INTEGER, DIMENSION(SIZE(PTA,1),SIZE(PTA,2)) :: NMAX
INTEGER                           :: ZNUM_EXP,ZMAS_EXP
REAL                              :: ZGAM,ZVEL_CARRIER,ZR
REAL                              :: ZW_M3,ZW_M0
REAL                              :: ZSUM_VEL_M3,ZSUM_M3,ZSUM_VEL_M0,ZSUM_M0
REAL                              :: ZDELTAR
REAL                              :: ZGAMB,ZGAM_BM3, ZGAM_BM3B
      INTEGER JI,JJ,II,JK,ILAYER ! Loop counter

REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !
  !*       0.3   Initialization of variables :
  !
IF (LHOOK) CALL DR_HOOK('BLOWSNW_VELGRAV1D',0,ZHOOK_HANDLE)

ZDELTAR = 1e-6 
ZGAM      = GAMMA_SURF(XEMIALPHA_SNW)
ILAYER=SIZE(PTA,2)

! Sutherland's equation for kinematic viscosity
DO JK=1,ILAYER
   ZMU(:,JK)=1.8325d-5*416.16/(PTA(:,JK)+120)*(PTA(:,JK)/296.16)*SQRT(PTA(:,JK)/296.16)/PRHODREF(:)
END DO

  !
  !*       1   Compute number- and mass-averaged settling velocity :
  !

IF(CSNOW_SEDIM=='TABC')  THEN
!
! Sedimentation of snow particles is computed according to Carrier's drag coofficient.
! To reduce computational time; look-up tables are used. They depend on the
! average radius and the pressure (interpolation)
!
  CALL BLOWSNW_SEDIM_LKT1D(PRG,PPABST,PVGK)

ELSE IF(CSNOW_SEDIM=='MITC') THEN  ! Sedimentation following Mitchell (1996)

  ZGAMB     = GAMMA_SURF(XEMIALPHA_SNW+3)
  ZGAM_BM3  = GAMMA_SURF(3*XBM3-1+XEMIALPHA_SNW)
  ZGAM_BM3B = GAMMA_SURF(3*XBM3+2+XEMIALPHA_SNW)

  ! Compute limit radius for integration of Mitchell's formulation
  ZR1(:,:)=RLIM(ZMU,PRHODREF,XBESTL_1)
  ZR2(:,:)=RLIM(ZMU,PRHODREF,XBESTL_2)
  ! Compute parameter avr for integration of Mitchell's formulation
  ZAM1(:,:)=AVR(XAM1,XBM1,PRHODREF,ZMU)
  ZAM2(:,:)=AVR(XAM2,XBM2,PRHODREF,ZMU)
  ZAM3(:,:)=AVR(XAM3,XBM3,PRHODREF,ZMU)

DO JI=1,SIZE(PTA,1)
   DO JK=1,ILAYER
!Number weighted terminal velocity
          PVGK(JI,JK,1)=(PBETA(JI,JK)**(3*XBM1-1)*ZAM1(JI,JK)*                          &
                                           GAMMA_INC_LOW(3*XBM1-1+XEMIALPHA_SNW,ZR1(JI,JK)/PBETA(JI,JK)) +     &
          PBETA(JI,JK)**(3*XBM2-1)*ZAM2(JI,JK)*                                                              &
          (GAMMA_INC_LOW(3*XBM2-1+XEMIALPHA_SNW,ZR2(JI,JK)/PBETA(JI,JK))-                                      &
          GAMMA_INC_LOW(3*XBM2-1+XEMIALPHA_SNW,ZR1(JI,JK)/PBETA(JI,JK)))+                                      & 
          PBETA(JI,JK)**(3*XBM3-1)*ZAM3(JI,JK)*                                                              &
          (ZGAM_BM3-GAMMA_INC_LOW(3*XBM3-1+XEMIALPHA_SNW,ZR2(JI,JK)/PBETA(JI,JK))))/ZGAM
!Mass weighted terminal velocity
          PVGK(JI,JK,2)=(PBETA(JI,JK)**(3*XBM1-1)*ZAM1(JI,JK)*                         &
                                           GAMMA_INC_LOW(3*XBM1+2+XEMIALPHA_SNW,ZR1(JI,JK)/PBETA(JI,JK)) +    &
          PBETA(JI,JK)**(3*XBM2-1)*ZAM2(JI,JK)*                                                             &
          (GAMMA_INC_LOW(3*XBM2+2+XEMIALPHA_SNW,ZR2(JI,JK)/PBETA(JI,JK))-                                     &
          GAMMA_INC_LOW(3*XBM2+2+XEMIALPHA_SNW,ZR1(JI,JK)/PBETA(JI,JK)))+                                     & 
          PBETA(JI,JK)**(3*XBM3-1)*ZAM3(JI,JK)*                                                             &
          (ZGAM_BM3B-GAMMA_INC_LOW(3*XBM3+2+XEMIALPHA_SNW,ZR2(JI,JK)/PBETA(JI,JK))))/ZGAMB
    ENDDO
ENDDO


ELSE IF(CSNOW_SEDIM=='CARR') THEN
! Settling velocity is computed according to Carrier's drag coofficient.
! This method is used in other blowing snow model such as PIEKTUK or SNOWSTORM
! We perfom a numerical integration since no analytical solution exists. 

ZAA(:,:) = 6.203*ZMU(:,:)/2
ZBB(:) = 5.516*XRHOLI/(4*PRHODREF(:))*XG
NMAX(:,:)=GET_INDEX(PBETA,ZDELTAR)

! Exponent used to weight the number-averaged falling speed
ZNUM_EXP=0.
! Exponent used to weight the mass-averaged falling speed
ZMAS_EXP=3.


DO JI=1,SIZE(PTA,1)
   DO JK=1,ILAYER
       ZSUM_M3=0.
       ZSUM_M0=0.
       ZSUM_VEL_M0=0.
       ZSUM_VEL_M3=0.
       DO II=1,NMAX(JI,JK)
          ZR = 1E-6+(II-0.5)*ZDELTAR
          ZVEL_CARRIER = - ZAA(JI,JK)/ZR+((ZAA(JI,JK)/ZR)**2.+ZBB(JI)*ZR)**0.5
          ZW_M0=ZR**(XEMIALPHA_SNW-1)*exp(-ZR/PBETA(JI,JK))/(PBETA(JI,JK)**XEMIALPHA_SNW*ZGAM)

          ZW_M3=ZR**ZMAS_EXP*ZW_M0
          ZW_M0=ZR**ZNUM_EXP*ZW_M0
            
          ZSUM_M3 = ZSUM_M3+ZW_M3*ZDELTAR
          ZSUM_M0 = ZSUM_M0+ZW_M0*ZDELTAR


          ZSUM_VEL_M0 = ZSUM_VEL_M0+ ZW_M0*ZVEL_CARRIER*ZDELTAR
          ZSUM_VEL_M3 = ZSUM_VEL_M3+ ZW_M3*ZVEL_CARRIER*ZDELTAR
       ENDDO
         PVGK(JI,JK,1) = ZSUM_VEL_M0/ZSUM_M0
         PVGK(JI,JK,2) = ZSUM_VEL_M3/ZSUM_M3
   ENDDO
ENDDO

ELSE 

! Settling velocity is set equal to 0.
PVGK(:,:,1) = 0.
PVGK(:,:,2) = 0.


END IF

IF (LHOOK) CALL DR_HOOK('BLOWSNW_VELGRAV1D',1,ZHOOK_HANDLE)

CONTAINS


FUNCTION RLIM(PMU,PRHODREF,PBEST_LIM) RESULT(PRLIM)
!
!!    PURPOSE
!!    -------
!     Calculate the radius of a sperical particle for a given Best Number 
!
!
USE MODD_CSTS,     ONLY : XRHOLI,XG
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)                                  :: PRHODREF ! (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)                                  :: PMU      ! (m2/s)
REAL,               INTENT(IN)                                  :: PBEST_LIM! (-)

!
REAL, DIMENSION(SIZE(PMU,1),SIZE(PMU,2)) :: PRLIM ! (m)
INTEGER    JK
!
DO JK=1,SIZE(PMU,2)
   PRLIM(:,JK)=(3./32.*PRHODREF(:)/(XRHOLI*XG)*PMU(:,JK)**2.*PBEST_LIM)**0.333333333
END DO

END FUNCTION RLIM

FUNCTION AVR(PARE,PBRE,PRHODEF,PMU) RESULT(PAVR)
!
!!    PURPOSE
!!    -------
!     Calculate the parameter av_r in Michell formulation 
!
!
USE MODD_CSTS,     ONLY : XRHOLI,XG

!
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                                  :: PARE      ! (-)
REAL,               INTENT(IN)                                  :: PBRE      ! (-)
REAL, DIMENSION(:), INTENT(IN)                                  :: PRHODEF ! (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)                                  :: PMU      ! (m2/s)

!
REAL, DIMENSION(SIZE(PMU,1),SIZE(PMU,2)) :: PAVR ! (-)
!
INTEGER JK

DO JK=1,SIZE(PMU,2)
    PAVR(:,JK)=2.**(3.*PBRE-1.)*PARE*PMU(:,JK)**(1.-2.*PBRE)*(4./3.*XRHOLI/PRHODEF(:)*XG)**PBRE
END DO

END FUNCTION AVR

FUNCTION GET_INDEX(PBETA,PDELTAR) RESULT(KMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the upper index in numerical integration of Carrier's formulation 
!     Index equals to 5* mean radius
!
!
USE MODD_BLOWSNW_SURF

!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                                  :: PDELTAR      ! (-)
REAL, DIMENSION(:,:), INTENT(IN)                                  :: PBETA ! (kg/m3)

!
INTEGER, DIMENSION(SIZE(PBETA,1),SIZE(PBETA,2)) :: KMAX ! (-)
!

KMAX(:,:)=int(PBETA(:,:)*XEMIALPHA_SNW*5/PDELTAR)


END FUNCTION GET_INDEX

END SUBROUTINE BLOWSNW_VELGRAV1D
