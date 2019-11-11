!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE COUPLING_SLT_n (SLT, &
      KI,                       &!I [nbr] number of sea points 
      KSLT,                     &!I [nbr] number of sea salt variables 
      PWIND,                    &!I Wind velocity
! ++ PIERRE / MARINE SSA - MODIF ++
      PWHEIGHT,                 &! Significant height of wind-generated waves (in ECMWF analyses)
                                 ! local pour l'instant, PWHEIGHT plus tard
      PSST,                     &! Sea water temperature (C) 
      PUSTAR,                   &! Friction velocity (ecmwf?) Calcule dans coupling_seafluxn.F90
! -- PIERRE / MARINE SSA - MODIF --
      PSFSLT                    &!O [kg/m2/sec] production flux of sea salt
      ) 
  
!PURPOSE
!-------
!  Compute sea salt emission  upon Vignatti et al, 2001
! ++ PIERRE / MARINE SSA - MODIF ++
!  Compute sea salt emission  upon Ovadnevaite et al, 2014
! -- PIERRE / MARINE SSA - MODIF --
!
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!AUTHOR
!-------
! P. Tulet
!
!
USE MODD_SLT_n, ONLY : SLT_t
!
USE MODD_CSTS, ONLY : XAVOGADRO, XPI
USE MODD_SLT_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!INPUT
!
TYPE(SLT_t), INTENT(INOUT) :: SLT
!
INTEGER, INTENT(IN)                :: KI             !I Number of sea points
INTEGER, INTENT(IN)                :: KSLT           !I Number of sea salt emission variables
REAL, DIMENSION(KI),      INTENT(IN)  :: PWIND       !I wind velocity
REAL, DIMENSION(KI,KSLT), INTENT(OUT) :: PSFSLT      !Out: kg/m2/s (index #2)
! ++ PIERRE / MARINE SSA - MODIF ++
REAL, DIMENSION(KI),      INTENT(INOUT)  :: PWHEIGHT !Significant height of wind-generated waves (in ECMWF analyses)
REAL, DIMENSION(KI),      INTENT(IN)  :: PUSTAR   !Friction velocity (ecmwf?) : Unite: m.s^(-2)?
REAL, DIMENSION(KI),      INTENT(IN)  :: PSST     ! Sea surface temperature (K)
! -- PIERRE / MARINE SSA - MODIF --

!LOCAL VARIABLES
REAL,DIMENSION(KI,JPMODE_SLT)  :: ZSFSLT_MDE ! sea salt flux from modes
INTEGER                        :: JN, JI, II !Counter for sea salt modes
REAL, DIMENSION(KI)            :: DZSPEED 
INTEGER, DIMENSION(KI)         :: WCL
REAL                           :: ZCONVERTFACM0_SLT  ![kg/mole*mole/molec] conversion factor
                                                     !for moment fluxes and used fluxes
REAL                           :: ZCONVERTFACM3_SLT 
REAL                           :: ZCONVERTFACM6_SLT 
!
! ++ PIERRE / MARINE SSA - MODIF ++

REAL, DIMENSION(5)             :: ZNUWATER   !  Temperature-dependant kinematic viscosity of 
                                             ! sea-water (table of data to interpolate) (m².s-¹)
REAL, DIMENSION(5)             :: ZWT ! Sea water temperature in table
REAL, DIMENSION(KI)            :: ZREYNOLDS ! Reynolds Number
REAL, DIMENSION(KI)            :: ZHVAGUE  ! sea wave height from wind if ZWS is unknown.
REAL, DIMENSION(KI)            :: ZVISCO ! Temperature-dependant kinematic viscosity
                                         ! of sea-water interpolated
! -- PIERRE / MARINE SSA - MODIF --
!
!REAL, PARAMETER :: mass1flux(0:40) = (/  &
!        0.000E+00, 2.483E-15, 2.591E-14, 1.022E-13, 2.707E-13, 5.761E-13,  &
!        1.068E-12, 1.800E-12, 2.829E-12, 4.215E-12, 6.023E-12, 8.317E-12, &
!        1.117E-11, 1.464E-11, 1.882E-11, 2.378E-11, 2.959E-11, 3.633E-11, &
!        4.409E-11, 5.296E-11, 6.301E-11, 7.433E-11, 8.693E-11, 1.012E-10, &
!        1.168E-10, 1.342E-10, 1.532E-10, 1.741E-10, 1.970E-10, 2.219E-10, &
!        2.489E-10, 2.781E-10, 3.097E-10, 3.437E-10, 3.803E-10, 4.195E-10, &
!        4.616E-10, 5.065E-10, 5.544E-10, 6.054E-10, 6.711E-10             /) 

!REAL, PARAMETER :: mass2flux(0:40) = (/  &
!        0.000E+00, 2.319E-13, 2.411E-12, 9.481E-12, 2.505E-11, 5.321E-11,  &
!        9.850E-11, 1.658E-10, 2.602E-10, 3.874E-10, 5.529E-10, 7.628E-10,  &
!       1.023E-09, 1.341E-09, 1.722E-09, 2.175E-09, 2.704E-09, 3.319E-09,  &
!       4.026E-09, 4.832E-09, 5.746E-09, 6.776E-09, 7.925E-09, 9.214E-09,  &
!       1.064E-08, 1.221E-08, 1.394E-08, 1.584E-08, 1.791E-08, 2.016E-08,  &
!        2.261E-08, 2.526E-08, 2.812E-08, 3.120E-08, 3.451E-08, 3.806E-08,  &
!        4.186E-08, 4.592E-08, 5.025E-08, 5.486E-08, 6.014E-08             /) 

!REAL, PARAMETER :: mass3flux(0:40) = (/ 0.0, &
!      1.783E-12, 1.579E-11, 5.852E-11, 1.501E-10, 3.134E-10, 5.740E-10, &
!      9.597E-10, 1.500E-09, 2.227E-09, 3.175E-09, 4.378E-09, 5.872E-09, &
!      7.698E-09, 9.897E-09, 1.250E-08, 1.556E-08, 1.912E-08, 2.323E-08, &
!      2.792E-08, 3.325E-08, 3.927E-08, 4.608E-08, 5.356E-08, 6.194E-08, &
!        7.121E-08, 8.143E-08, 9.266E-08, 1.049E-07, 1.183E-07, 1.329E-07, &
!        1.487E-07, 1.658E-07, 1.843E-07, 2.041E-07, 2.255E-07, 2.484E-07, &
!        2.729E-07, 2.991E-07, 3.270E-07, 3.517E-07 /) 

REAL, PARAMETER :: HVAGUE(1:9) = (/ 0., 0.1, 0.5, 1.25, 2.5, 4., 6., 9., 14. /)
REAL, PARAMETER :: VVENT(1:9) = (/  1., 2.7, 4.1, 6.3, 8.3, 11.1, 13.8, &
                                    16.6, 19.4/)

REAL, PARAMETER :: NUMB1FLUX(0:40) = (/ &
         0.000E+00, 3.004E+01, 3.245E+02, 1.306E+03, 3.505E+03, 7.542E+03,  &
          1.410E+04, 2.394E+04, 3.787E+04, 5.674E+04, 8.147E+04, 1.130E+05,  &
          1.523E+05, 2.005E+05, 2.586E+05, 3.278E+05, 4.091E+05, 5.037E+05,  &
          6.129E+05, 7.379E+05, 8.800E+05, 1.041E+06, 1.220E+06, 1.422E+06,  &
          1.646E+06, 1.893E+06, 2.166E+06, 2.466E+06, 2.794E+06, 3.152E+06,  &
          3.541E+06, 3.962E+06, 4.419E+06, 4.911E+06, 5.441E+06, 6.011E+06,  &
          6.621E+06, 7.274E+06, 7.972E+06, 8.716E+06, 8.801E+06             /) 

REAL, PARAMETER :: NUMB2FLUX(0:40) = (/  &
          0.000E+00, 1.934E+01, 2.068E+02, 8.271E+02, 2.211E+03, 4.741E+03,  &
          8.841E+03, 1.497E+04, 2.363E+04, 3.534E+04, 5.066E+04, 7.017E+04,  &
          9.447E+04, 1.242E+05, 1.600E+05, 2.025E+05, 2.525E+05, 3.106E+05,  &
          3.776E+05, 4.542E+05, 5.413E+05, 6.395E+05, 7.501E+05, 8.726E+05,  &
          1.009E+06, 1.160E+06, 1.327E+06, 1.509E+06, 1.709E+06, 1.927E+06,  &
          2.163E+06, 2.420E+06, 2.697E+06, 2.996E+06, 3.318E+06, 3.664E+06,  &
          4.034E+06, 4.430E+06, 4.852E+06, 5.303E+06, 5.740E+06             /) 

REAL, PARAMETER :: NUMB3FLUX(0:40) = (/ 0.0, &
        4.340E-01, 5.217E+00, 2.241E+01, 6.301E+01, 1.404E+02, 2.703E+02, &
        4.699E+02, 7.584E+02, 1.157E+03, 1.687E+03, 2.373E+03, 3.240E+03, &
        4.314E+03, 5.625E+03, 7.197E+03, 9.063E+03, 1.126E+04, 1.380E+04, &
        1.674E+04, 2.011E+04, 2.393E+04, 2.827E+04, 3.311E+04, 3.853E+04, &
        4.457E+04, 5.126E+04, 5.864E+04, 6.675E+04, 7.564E+04, 8.535E+04, &
        9.592E+04, 1.074E+05, 1.198E+05, 1.333E+05, 1.478E+05, 1.633E+05, &
        1.801E+05, 1.980E+05, 2.172E+05, 2.353E+05 /) 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!!
!!    MESONH carries the following units during transport:
!!    M0=#/molec_{air}
!!    M6=um6/molec_{air}*1.d6
!!    The surface model should have (for sea salt)
!!    M0=#/m3*[kg_{slt}/mole_{slt}/XAVOGADRO]
!!    M3=kg/m3
!!    M6=um6/m3
!!    REFERENCE
!!    ---------
!!    Tulet et al, ORILAM manuscript for transformation of modal parameters
!!    J. Geophys. Res., 110, D18201, doi:10.1029/2004JD005716
!
!Initialize output which is total flux of sea salt (kg/m2/sec). 
IF (LHOOK) CALL DR_HOOK('COUPLING_SLT_N',0,ZHOOK_HANDLE)
!
!Factor which is needed so that all gains normal units when leaving ground paramn
ZCONVERTFACM0_SLT = XMOLARWEIGHT_SLT / XAVOGADRO !(kg_slt/mol_slt)/(molec/mol)
!Factor which is needed for moment 6, there is a factor 1.d6 transported around in M6 in MESONH
ZCONVERTFACM6_SLT = XMOLARWEIGHT_SLT / XAVOGADRO*1.d6
ZCONVERTFACM3_SLT = 4./3.*XPI*XDENSITY_SLT / 1.d18
!
PSFSLT(:,:)=0.d0
!
!+ Marine
IF (CEMISPARAM_SLT .eq. "Ova14") THEN ! Rajouter Ova14 dans fichier initialisation
  ZHVAGUE(:)  = 0.
  DO II = 1, 8
!++cb++19/10/16 modif de la formule : + de vent => vagues + hautes
!    WHERE ((PWIND(:) .GT. VVENT(II)).AND.(PWIND(:) .LT. VVENT(II+1)))
    WHERE ((PWIND(:) .GT. VVENT(II)).AND.(PWIND(:) .LT. VVENT(II+1)))
!      ZHVAGUE(:) = HVAGUE(II) + (VVENT(II+1)  - PWIND(:)) * &
      ZHVAGUE(:) = HVAGUE(II) + (PWIND(:)     - VVENT(II+1)) * &
                                (HVAGUE(II+1) - HVAGUE(II)) / &
                                (VVENT(II+1)  - VVENT(II))
!--cb--
    ENDWHERE
  ENDDO

  WHERE (PWIND(:) .GE. VVENT(9))
    ZHVAGUE(:) = HVAGUE(9)
  END WHERE

  WHERE (PWHEIGHT(:) .EQ. -1.)
    PWHEIGHT(:) = ZHVAGUE(:)
  END WHERE

  ZWT = (/ 273.15, 283.15, 293.15, 303.15, 313.15 /)  ! Unite : K
  ZNUWATER = (/ 1.854E-6, 1.36E-6, 1.051E-6, 0.843E-6, 0.695E-6 /) 
! Unite : m².s^(-1) Pour une salinite = 35g/kg.
! En mer Mediterranee = 38.5g/kg (Lewis and Schwartz)

! Initialisation des valeurs de ZVISCO, ZREYNOLDS
  ZVISCO(:)    = 0.
  ZREYNOLDS(:) = 0.

  ! Tableau d'interpolation pour calculer ZNUWATER en fonction de la SST
  ! Cas ou 0 < SST < 10 C
  WHERE ((PSST(:) >= 273.15).AND.(PSST(:) < 283.15))
    ZVISCO(:) = ZNUWATER(1) + (PSST(:) - ZWT(1)) * (ZNUWATER(2)-ZNUWATER(1)) / &
                (ZWT(2) - ZWT(1))
  ENDWHERE

  ! Cas ou 10 < SST < 20 C
  WHERE ((PSST(:) >= 283.15).AND.(PSST(:) < 293.15))
    ZVISCO(:) = ZNUWATER(2) + (PSST(:) - ZWT(2)) * (ZNUWATER(3)-ZNUWATER(2)) / &
                (ZWT(3) - ZWT(2))
  ENDWHERE

  ! Cas ou 20 < SST < 30 C
  WHERE ((PSST(:) >= 293.15).AND.(PSST(:) < 303.15))
    ZVISCO(:) = ZNUWATER(3) + (PSST(:) - ZWT(3)) * (ZNUWATER(4)-ZNUWATER(3)) / &
                (ZWT(4) - ZWT(3))
  ENDWHERE

  ! Cas ou 30 < SST < 40 C
  WHERE ((PSST(:) >= 303.15).AND.(PSST(:) < 313.15))
    ZVISCO(:) = ZNUWATER(4) + (PSST(:) - ZWT(4)) * (ZNUWATER(5)-ZNUWATER(4)) / &
                (ZWT(5) - ZWT(4))
  ENDWHERE

! Calcul du nombre de Reynolds
  ZREYNOLDS(:) = (PUSTAR(:) * PWHEIGHT(:)) / ZVISCO(:)

! Calcul du flux en nombre pour chaque mode

! Ovadnevaite et al. 2014 
!!!!! Total number flux, Unite ZSDSLT_MDE ne correspond pas au total number
!flux mais au size dependent SSA production flux

! Ecrire equation integration pour chaque mode

!Condition d'emission : ZREYNOLDS > 1E5

  WHERE (ZREYNOLDS(:) > 1.E5)
    ZSFSLT_MDE(:,1) = 104.51 * ( ZREYNOLDS(:) - 1.E5)**0.556
    ZSFSLT_MDE(:,2) = 0.044  * ( ZREYNOLDS(:) - 1.E5)**1.08
    ZSFSLT_MDE(:,3) = 149.64 * ( ZREYNOLDS(:) - 1.E5)**0.545
    ZSFSLT_MDE(:,4) = 2.96   * ( ZREYNOLDS(:) - 1.E5)**0.79
  ENDWHERE
  WHERE (ZREYNOLDS(:) > 2.E5)
    ZSFSLT_MDE(:,5) = 0.52   * ( ZREYNOLDS(:) - 2.E5)**0.87
  ENDWHERE



  WHERE (ZREYNOLDS(:) <= 1.E5)
    ZSFSLT_MDE(:,1) = 1.E-10
    ZSFSLT_MDE(:,2) = 1.E-10
    ZSFSLT_MDE(:,3) = 1.E-10
    ZSFSLT_MDE(:,4) = 1.E-10
  ENDWHERE
  WHERE (ZREYNOLDS(:) <= 2.E5)
    ZSFSLT_MDE(:,5) = 1.E-10
  ENDWHERE

! Controle avec des valeurs limites , Pas besoin de la conversion 1E4 pour Ova
! car deja en m-2
  ZSFSLT_MDE(:,1) = MAX(ZSFSLT_MDE(:,1) , 1.E-10)
  ZSFSLT_MDE(:,2) = MAX(ZSFSLT_MDE(:,2) , 1.E-10)
  ZSFSLT_MDE(:,3) = MAX(ZSFSLT_MDE(:,3) , 1.E-10)
  ZSFSLT_MDE(:,4) = MAX(ZSFSLT_MDE(:,4) , 1.E-10)
  ZSFSLT_MDE(:,5) = MAX(ZSFSLT_MDE(:,5) , 1.E-10)
!- Marine

ELSEIF (CEMISPARAM_SLT .eq. "Vig01") THEN
! Vignatti et al. 2001 (in particles.cm-2.s-1) : en #.cm-3 en fait
  ZSFSLT_MDE(:,1) =  10.**(0.09  *PWIND(:) + 0.283)   ! fine mode
  ZSFSLT_MDE(:,2) =  10.**(0.0422*PWIND(:) + 0.288)   ! median mode
  ZSFSLT_MDE(:,3) =  10.**(0.069 *PWIND(:) - 3.5)     ! coarse mode

! convert into  particles.m-2.s-1)
  ZSFSLT_MDE(:,1) = MAX(ZSFSLT_MDE(:,1) * 1.E4, 1.E-10)
  ZSFSLT_MDE(:,2) = MAX(ZSFSLT_MDE(:,2) * 1.E4, 1.E-10)
  ZSFSLT_MDE(:,3) = MAX(ZSFSLT_MDE(:,3) * 1.E4, 1.E-10)
!
ELSEIF (CEMISPARAM_SLT .eq. "Sch04") THEN! Use Schultz et al., 2004
  WCL(:) = INT(PWIND(:))
  WCL(:) = MAX (0, MIN(WCL(:), 39))
 
  DZSPEED(:) = MAX(0., MIN(PWIND(:) - FLOAT(WCL(:)), 1.))
 !
 ! Flux given  in  particles.m-2 s-1
 !
  DO JI = 1, KI
   !plm-gfortran
    ZSFSLT_MDE(JI,1) = NUMB1FLUX(WCL(JI)) + &
                      (NUMB1FLUX(WCL(JI)+1)-NUMB1FLUX(WCL(JI)))*DZSPEED(JI)
    ZSFSLT_MDE(JI,2) = NUMB2FLUX(WCL(JI)) + &
                      (NUMB2FLUX(WCL(JI)+1)-NUMB2FLUX(WCL(JI)))*DZSPEED(JI)
    ZSFSLT_MDE(JI,3) = NUMB3FLUX(WCL(JI)) + &
                      (NUMB3FLUX(WCL(JI)+1)-NUMB3FLUX(WCL(JI)))*DZSPEED(JI)
   !plm-gfortran
  END DO
END IF
!
DO JN = 1, JPMODE_SLT

! convert  particles.m-2 s-1 into kg.m-2.s-1
! N'est calculé que pour le moment 3 (en masse), la conversion pour les autres
! flux de moments se fait plus tard (mode_dslt_surf.F90 MASSFLUX2MOMENTFLUX)
!+Marine
  !
  IF (LVARSIG_SLT) THEN ! cas 3 moment

    PSFSLT(:,2+(JN-1)*3) = ZSFSLT_MDE(:,JORDER_SLT(JN)) &
                            * ((SLT%XEMISRADIUS_SLT(JORDER_SLT(JN))**3) &
                            * EXP(4.5 * LOG(SLT%XEMISSIG_SLT(JORDER_SLT(JN)))**2)) &
                            * ZCONVERTFACM3_SLT

  ELSEIF (LRGFIX_SLT) THEN ! cas 1 moment
    PSFSLT(:,JN) =  ZSFSLT_MDE(:,JORDER_SLT(JN)) &
                      * (SLT%XEMISRADIUS_SLT(JORDER_SLT(JN))**3) &
                      * EXP(4.5 * LOG(SLT%XEMISSIG_SLT(JORDER_SLT(JN)))**2) &
                      * ZCONVERTFACM3_SLT

  ELSE ! cas 2 moments

    PSFSLT(:,2+(JN-1)*2) = ZSFSLT_MDE(:,JORDER_SLT(JN)) &
                            * ((SLT%XEMISRADIUS_SLT(JORDER_SLT(JN))**3) &
                            * EXP(4.5 * LOG(SLT%XEMISSIG_SLT(JORDER_SLT(JN)))**2)) &
                            * ZCONVERTFACM3_SLT
! -- PIERRE / MARINE SSA - MODIF --
  END IF
END DO


IF (LHOOK) CALL DR_HOOK('COUPLING_SLT_N',1,ZHOOK_HANDLE)
END SUBROUTINE COUPLING_SLT_n
