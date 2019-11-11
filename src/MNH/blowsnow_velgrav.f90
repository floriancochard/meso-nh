!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_BLOWSNOW_VELGRAV
!!   ##############################
!!
INTERFACE
!!
SUBROUTINE BLOWSNOW_VELGRAV(PSVT, PTHT, PABST,PRHODREF,PVGK)
IMPLICIT NONE
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT) :: PSVT    ! Blowing snow concentration
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT, PABST, PRHODREF
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGK
END SUBROUTINE BLOWSNOW_VELGRAV
!!
END INTERFACE
!!
END MODULE MODI_BLOWSNOW_VELGRAV
!!
!!   #######################################
     SUBROUTINE BLOWSNOW_VELGRAV(PSVT, PTHT, PABST, PRHODREF,PVGK)
!!   #######################################
!!
!!   PURPOSE
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
!!
!-----------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_CSTS_BLOWSNOW
USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODI_GAMMA_INC_LOW
USE MODE_BLOWSNOW_SEDIM_LKT
USE MODE_BLOWSNOW_PSD
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT) :: PSVT    ! Blowing snow concentration
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT, PABST, PRHODREF
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PVGK
!
!*      0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZTEMP,ZMU
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   :: ZRG, ZBETA,ZMOB

REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZR1,ZR2
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZAM1,ZAM2,ZAM3
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZAA, ZBB   
INTEGER, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: NMAX  
INTEGER                           :: ZNUM_EXP,ZMAS_EXP
REAL                                                       :: ZGAM,ZVEL_CARRIER,ZR
REAL                                                       :: ZW_M0,ZW_M3
REAL                                                       :: ZSUM_VEL_M0,ZSUM_VEL_M3,ZSUM_M3,ZSUM_M0
REAL                                                       :: ZDELTAR
REAL                                                       :: ZGAM_BM1,ZGAM_BM2,ZGAM_BM3,ZGAMB
REAL                                                       :: ZGAM_BM1B,ZGAM_BM2B,ZGAM_BM3B
INTEGER                                                 :: JI,JJ,JK,II     !Loop counter
LOGICAL :: LNONEFFICIENT
!

ZDELTAR = 1e-6                   ! Bin size (m)

ZGAM      = GAMMA(XALPHA_SNOW)

!
!-----------------------------------------------------------------
!
!*       2.   compute BETA and RG
!
CALL PPP2SNOW(PSVT, PRHODREF, PBET3D=ZBETA, PRG3D=ZRG)
!
!-----------------------------------------------------------------
!
!*       3.   compute temperature and kinematic viscosity
!
!  Temperature
ZTEMP(:,:,:)=PTHT(:,:,:)*(PABST(:,:,:)/XP00)**(XRD/XCPD)
!
!  Sutherland's equation for kinematic viscosity
ZMU(:,:,:)=1.8325d-5*416.16/(ZTEMP(:,:,:)+120)*(ZTEMP(:,:,:)/296.16)*SQRT(ZTEMP(:,:,:)/296.16)/PRHODREF(:,:,:)
!
!-----------------------------------------------------------------
!
!*       4.   compute number and mass-averaged settling velocity
!
IF(CSNOWSEDIM=='NONE') THEN ! No sedimentation
DO JI=1,SIZE(PSVT,1)
   DO JJ=1,SIZE(PSVT,2)
      DO JK=1,SIZE(PSVT,3)        
         PVGK(JI,JJ,JK,1)= 0.
         PVGK(JI,JJ,JK,2)= 0.
!         PVGK(JI,JJ,JK,3)= 0.
      ENDDO
    ENDDO
ENDDO
END IF

IF(CSNOWSEDIM=='MITC') THEN ! Sedimentation following Mitchell (1996)

LNONEFFICIENT = .FALSE.

IF(LNONEFFICIENT) THEN    

ZGAMB     = GAMMA(XALPHA_SNOW+3)
ZGAM_BM1  = GAMMA(3*XBM1-1+XALPHA_SNOW)
ZGAM_BM2  = GAMMA(3*XBM2-1+XALPHA_SNOW)
ZGAM_BM3  = GAMMA(3*XBM3-1+XALPHA_SNOW)
ZGAM_BM1B = GAMMA(3*XBM1+2+XALPHA_SNOW)
ZGAM_BM2B = GAMMA(3*XBM2+2+XALPHA_SNOW)
ZGAM_BM3B = GAMMA(3*XBM3+2+XALPHA_SNOW)

  ! Compute limit radius for integration of Mitchell's formulation
ZR1(:,:,:)=RLIM(ZMU,PRHODREF,XBESTL_1)
ZR2(:,:,:)=RLIM(ZMU,PRHODREF,XBESTL_2)
  ! Compute parameter avr for integration of Mitchell's formulation
ZAM1(:,:,:)=AVR(XAM1,XBM1,PRHODREF,ZMU)
ZAM2(:,:,:)=AVR(XAM2,XBM2,PRHODREF,ZMU)
ZAM3(:,:,:)=AVR(XAM3,XBM3,PRHODREF,ZMU)
        
DO JI=1,SIZE(PSVT,1)
   DO JJ=1,SIZE(PSVT,2)
      DO JK=1,SIZE(PSVT,3)
!Number weighted terminal velocity
          PVGK(JI,JJ,JK,1)=(ZBETA(JI,JJ,JK)**(3*XBM1-1)*ZAM1(JI,JJ,JK)*ZGAM_BM1*                          &
                                           GAMMA_INC(3*XBM1-1+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)) +                  &
          ZBETA(JI,JJ,JK)**(3*XBM2-1)*ZAM2(JI,JJ,JK)*ZGAM_BM2*                                                              &
          (GAMMA_INC(3*XBM2-1+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))-                                                   &
          GAMMA_INC(3*XBM2-1+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)))+                                                   & 
          ZBETA(JI,JJ,JK)**(3*XBM3-1)*ZAM3(JI,JJ,JK)*ZGAM_BM3*                                                              &
          (1.-GAMMA_INC(3*XBM3-1+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))))/ZGAM
!Mass weighted terminal velocity
          PVGK(JI,JJ,JK,2)=(ZBETA(JI,JJ,JK)**(3*XBM1-1)*ZAM1(JI,JJ,JK)*ZGAM_BM1B*                         &
                                           GAMMA_INC(3*XBM1+2+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)) +                  &
          ZBETA(JI,JJ,JK)**(3*XBM2-1)*ZAM2(JI,JJ,JK)*ZGAM_BM2B*                                                             &
          (GAMMA_INC(3*XBM2+2+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))-                                                   &
          GAMMA_INC(3*XBM2+2+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)))+                                                   & 
          ZBETA(JI,JJ,JK)**(3*XBM3-1)*ZAM3(JI,JJ,JK)*ZGAM_BM3B*                                                             &
          (1.-GAMMA_INC(3*XBM3+2+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))))/ZGAMB
          !Mass weighted terminal velocity for mobility index          
!          PVGK(JI,JJ,JK,3)= PVGK(JI,JJ,JK,2)
       ENDDO
    ENDDO
ENDDO

ELSE
! Fast integration of the incomplete gamma function following Blahak (2010)
! Blahak U., Efficient approximation of the incomplete gamma function for use
!   in cloud model applications, GMD, 3, 329-336, 2010

ZGAMB     = GAMMA(XALPHA_SNOW+3)
ZGAM_BM3  = GAMMA(3*XBM3-1+XALPHA_SNOW)
ZGAM_BM3B = GAMMA(3*XBM3+2+XALPHA_SNOW)

  ! Compute limit radius for integration of Mitchell's formulation
ZR1(:,:,:)=RLIM(ZMU,PRHODREF,XBESTL_1)
ZR2(:,:,:)=RLIM(ZMU,PRHODREF,XBESTL_2)
  ! Compute parameter avr for integration of Mitchell's formulation
ZAM1(:,:,:)=AVR(XAM1,XBM1,PRHODREF,ZMU)
ZAM2(:,:,:)=AVR(XAM2,XBM2,PRHODREF,ZMU)
ZAM3(:,:,:)=AVR(XAM3,XBM3,PRHODREF,ZMU)

DO JI=1,SIZE(PSVT,1)
   DO JJ=1,SIZE(PSVT,2)
      DO JK=1,SIZE(PSVT,3)
!Number weighted terminal velocity
          PVGK(JI,JJ,JK,1)=(ZBETA(JI,JJ,JK)**(3*XBM1-1)*ZAM1(JI,JJ,JK)*                          &
                                           GAMMA_INC_LOW(3*XBM1-1+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)) +     &
          ZBETA(JI,JJ,JK)**(3*XBM2-1)*ZAM2(JI,JJ,JK)*                                                              &
          (GAMMA_INC_LOW(3*XBM2-1+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))-                                      &
          GAMMA_INC_LOW(3*XBM2-1+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)))+                                      & 
          ZBETA(JI,JJ,JK)**(3*XBM3-1)*ZAM3(JI,JJ,JK)*                                                              &
          (ZGAM_BM3-GAMMA_INC_LOW(3*XBM3-1+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))))/ZGAM
!Mass weighted terminal velocity
          PVGK(JI,JJ,JK,2)=(ZBETA(JI,JJ,JK)**(3*XBM1-1)*ZAM1(JI,JJ,JK)*                         &
                                           GAMMA_INC_LOW(3*XBM1+2+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)) +    &
          ZBETA(JI,JJ,JK)**(3*XBM2-1)*ZAM2(JI,JJ,JK)*                                                             &
          (GAMMA_INC_LOW(3*XBM2+2+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))-                                     &
          GAMMA_INC_LOW(3*XBM2+2+XALPHA_SNOW,ZR1(JI,JJ,JK)/ZBETA(JI,JJ,JK)))+                                     & 
          ZBETA(JI,JJ,JK)**(3*XBM3-1)*ZAM3(JI,JJ,JK)*                                                             &
          (ZGAM_BM3B-GAMMA_INC_LOW(3*XBM3+2+XALPHA_SNOW,ZR2(JI,JJ,JK)/ZBETA(JI,JJ,JK))))/ZGAMB
          !Mass weighted terminal velocity for mobility index          
!          PVGK(JI,JJ,JK,3)= PVGK(JI,JJ,JK,2)
       ENDDO
    ENDDO
ENDDO

END IF

END IF

IF(CSNOWSEDIM=='CARR') THEN
! Settling velocity is computed according to Carrier's drag coofficient.
! This method is used in other blowing snow model such as PIEKTUK or SNOWSTORM
! We perfom a numerical integration since no analytical solution exists. 

ZAA(:,:,:) = 6.203*ZMU(:,:,:)/2
ZBB(:,:,:) = 5.516*XRHOLI/(4*PRHODREF(:,:,:))*XG
NMAX(:,:,:)=GET_INDEX(ZBETA,ZDELTAR)


! Exponent used to weight the number-averaged falling speed
ZNUM_EXP=0.
! Exponent used to weight the mass-averaged falling speed
ZMAS_EXP=3.

DO JI=1,SIZE(PSVT,1)
   DO JJ=1,SIZE(PSVT,2)
      DO JK=1,SIZE(PSVT,3)
         ZSUM_M3=0.
         ZSUM_M0=0.
         ZSUM_VEL_M0=0.
         ZSUM_VEL_M3=0.
         DO II=1,NMAX(JI,JJ,JK)
            ZR = 1E-6+(II-0.5)*ZDELTAR
            ZVEL_CARRIER = - ZAA(JI,JJ,JK)/ZR+((ZAA(JI,JJ,JK)/ZR)**2.+ZBB(JI,JJ,JK)*ZR)**0.5
            ZW_M0=ZR**(XALPHA_SNOW-1)*exp(-ZR/ZBETA(JI,JJ,JK))/(ZBETA(JI,JJ,JK))**XALPHA_SNOW*ZGAM

          ZW_M3=ZR**ZMAS_EXP*ZW_M0
          ZW_M0=ZR**ZNUM_EXP*ZW_M0
            
          ZSUM_M3 = ZSUM_M3+ZW_M3*ZDELTAR
          ZSUM_M0 = ZSUM_M0+ZW_M0*ZDELTAR

          ZSUM_VEL_M0 = ZSUM_VEL_M0+ ZW_M0*ZVEL_CARRIER*ZDELTAR
          ZSUM_VEL_M3 = ZSUM_VEL_M3+ ZW_M3*ZVEL_CARRIER*ZDELTAR
         ENDDO
         PVGK(JI,JJ,JK,1) = ZSUM_VEL_M0/ZSUM_M0
         PVGK(JI,JJ,JK,2) = ZSUM_VEL_M3/ZSUM_M3
         !PVGK(JI,JJ,JK,3) = PVGK(JI,JJ,JK,2) 
      ENDDO
   ENDDO
ENDDO

END IF


IF(CSNOWSEDIM=='TABC') THEN
! Sedimentation of snow particles is computed according to Carrier's drag coofficient.
! To reduce computational time; look-up tables are used. They depend on the
! average radius and the pressure (interpolation)
CALL BLOWSNOW_SEDIM_LKT(ZRG,PABST,PVGK)
END IF


CONTAINS

FUNCTION RLIM(PMU,PRHODREF,PBEST_LIM) RESULT(PRLIM)
!
!!    PURPOSE
!!    -------
!     Calculate the radius of a sperical particle for a given Best Number 
!
!
USE MODD_CSTS_BLOWSNOW,     ONLY : XRHOLI,XG
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PRHODREF ! (kg/m3)
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PMU      ! (m2/s)
REAL,                   INTENT(IN)                                  :: PBEST_LIM! (-)

!
REAL, DIMENSION(SIZE(PMU,1),SIZE(PMU,2),SIZE(PMU,3)) :: PRLIM ! (m)
!
PRLIM(:,:,:)=(3./32.*PRHODREF(:,:,:)/(XRHOLI*XG)*PMU(:,:,:)**2.*PBEST_LIM)**0.333333333

END FUNCTION RLIM

FUNCTION AVR(PARE,PBRE,PRHODREF,PMU) RESULT(PAVR)
!
!!    PURPOSE
!!    -------
!     Calculate the parameter av_r in KC02 formulation (Eq. 3.1)
!
!
USE MODD_CSTS_BLOWSNOW,     ONLY : XRHOLI,XG


!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                                  :: PARE      ! (-)
REAL,               INTENT(IN)                                  :: PBRE      ! (-)
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PRHODREF ! (kg/m3)
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PMU      ! (m2/s)

!
REAL, DIMENSION(SIZE(PMU,1),SIZE(PMU,2),SIZE(PMU,3)) :: PAVR ! (-)
!


PAVR(:,:,:)=2.**(3.*PBRE-1.)*PARE*PMU(:,:,:)**(1.-2.*PBRE)*(4./3.*XRHOLI/PRHODREF(:,:,:)*XG)**PBRE

END FUNCTION AVR

FUNCTION GET_INDEX(PBETA,PDELTAR) RESULT(KMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the upper index in numerical integration of Carrier's formulation 
!     Index equals to 5* mean radius
!
!
USE MODD_BLOWSNOW,     ONLY : XALPHA_SNOW


!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                                  :: PDELTAR      ! (-)
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PBETA ! (kg/m3)

!
INTEGER, DIMENSION(SIZE(PBETA,1),SIZE(PBETA,2),SIZE(PBETA,3)) :: KMAX ! (-)
!

KMAX(:,:,:)=int(PBETA(:,:,:)*XALPHA_SNOW*5/PDELTAR)


END FUNCTION GET_INDEX

END SUBROUTINE BLOWSNOW_VELGRAV
