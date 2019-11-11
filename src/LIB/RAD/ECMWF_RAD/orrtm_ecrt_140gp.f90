!****************** SUBROUTINE RRTM_ECRT_140GP **************************

SUBROUTINE ORRTM_ECRT_140GP &
 &( iplon, klon , klev, kcld &
 &, paer , paph , pap &
 &, pts  , pth  , pt &
 &, zemis, zemiw &
 &, pq   , pcco2, pozn, pcldf, ptaucld, ptclear &
 &, CLDFRAC,TAUCLD,COLDRY,WKL,WX &
 &, TAUAERL,PAVEL,TAVEL,PZ,TZ,TBOUND,NLAYERS,SEMISS,IREFLECT)

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     Read in atmospheric profile from ECMWF radiation code, and prepare it
!     for use in RRTM.  Set other RRTM input parameters.  Values are passed
!     back through existing RRTM arrays and commons.

!- Modifications

!     2000-05-15 Deborah Salmond  Speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,JPLAY   ,&
            &JPINPX
USE OYOERAD   , ONLY : NOVLP
USE OYOERDI   , ONLY : RCARDI   ,RCH4     ,RN2O    ,RCFC11  ,RCFC12
USE OYOESW    , ONLY : RAER

!------------------------------Arguments--------------------------------


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: iplon
INTEGER_M :: kcld

!     DUMMY REAL SCALARS
REAL_B :: ptclear

INTEGER_M :: kidia                   ! First atmosphere index
INTEGER_M :: kfdia                   ! Last atmosphere index
INTEGER_M :: klon                    ! Number of atmospheres (longitudes)
INTEGER_M :: klev                    ! Number of atmospheric layers
REAL_B :: paer(klon,6,klev)          ! Aerosol optical thickness
REAL_B :: pap(klon,klev)             ! Layer pressures (Pa)
REAL_B :: paph(klon,klev+1)          ! Interface pressures (Pa)
REAL_B :: pts(klon)                  ! Surface temperature (K)
REAL_B :: pth(klon,klev+1)           ! Interface temperatures (K)
REAL_B :: pt(klon,klev)              ! Layer temperature (K)
REAL_B :: zemis(klon)                ! Non-window surface emissivity
REAL_B :: zemiw(klon)                ! Window surface emissivity
REAL_B :: pq(klon,klev)              ! H2O specific humidity (mmr)
REAL_B :: pozn(klon,klev)            ! O3 mass mixing ratio
REAL_B :: pcco2                      ! CO2 mass mixing ratio
!      real rch4                       ! CH4 mass mixing ratio
!      real rn2o                       ! N2O mass mixing ratio
!      real rcfc11                     ! CFC11 mass mixing ratio
!      real rcfc12                     ! CFC12 mass mixing ratio
REAL_B :: pcldf(klon,klev)           ! Cloud fraction
REAL_B :: ptaucld(klon,klev,JPBAND)  ! Cloud optical depth
REAL_B :: CLDFRAC(JPLAY)             ! Cloud fraction
REAL_B :: TAUCLD(JPLAY,JPBAND)       ! Spectral optical thickness
REAL_B :: COLDRY(JPLAY)
REAL_B :: WKL(JPINPX,JPLAY)
REAL_B :: WX(JPXSEC,JPLAY)           ! Amount of trace gases

!- from AER
REAL_B :: TAUAERL(JPLAY,JPBAND)

!- from PROFILE             
REAL_B :: PAVEL(JPLAY)
REAL_B :: TAVEL(JPLAY)
REAL_B :: PZ(0:JPLAY)
REAL_B :: TZ(0:JPLAY)
REAL_B :: TBOUND
INTEGER_M :: NLAYERS

!- from SURFACE             
REAL_B :: SEMISS(JPBAND)
INTEGER_M :: IREFLECT

REAL_B :: ztauaer(5)
REAL_B :: zc1j(0:klev)               ! total cloud from top and level k
INTEGER_M ::  IXINDX(JPINPX)         ! Indices of trace gases accounted for

REAL_B :: amd                  ! Effective molecular weight of dry air (g/mol)
REAL_B :: amw                  ! Molecular weight of water vapor (g/mol)
REAL_B :: amco2                ! Molecular weight of carbon dioxide (g/mol)
REAL_B :: amo                  ! Molecular weight of ozone (g/mol)
REAL_B :: amch4                ! Molecular weight of methane (g/mol)
REAL_B :: amn2o                ! Molecular weight of nitrous oxide (g/mol)
REAL_B :: amc11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL_B :: amc12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL_B :: avgdro               ! Avogadro's number (molecules/mole)
REAL_B :: gravit               ! Gravitational acceleration (cm/sec2)

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
data amd   /  28.970_JPRB    /
data amw   /  18.0154_JPRB   /
data amco2 /  44.011_JPRB    /
data amo   /  47.9982_JPRB   /
data amch4 /  16.043_JPRB    /
data amn2o /  44.013_JPRB    /
data amc11 / 137.3686_JPRB   /
data amc12 / 120.9140_JPRB   /
data avgdro/ 6.02214E23_JPRB /
data gravit/ 9.80665E02_JPRB /

!     LOCAL INTEGER SCALARS
INTEGER_M :: IATM, IMOL, IX, IXMAX, J1, J2, JAE, JB, JK, JL, L, JIS
INTEGER_M :: NMOL, NXMOL

!     LOCAL REAL SCALARS
REAL_B :: amm, ZCLDLY, ZCLEAR, ZCLOUD, ZEPSEC

! ***

! *** mji
! Initialize all molecular amounts and aerosol optical depths to zero here, 
! then pass ECRT amounts into RRTM arrays below.

!      DATA ZWKL /MAXPRDW*0.0/
!      DATA ZWX  /MAXPROD*0.0/
!      DATA KREFLECT /0/

! Activate cross section molecules:
!     NXMOL     - number of cross-sections input by user
!     IXINDX(I) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in RRTM
!                 = 1 -- CCL4
!                 = 2 -- CFC11
!                 = 3 -- CFC12
!                 = 4 -- CFC22
!      DATA KXMOL  /2/
!      DATA KXINDX /0,2,3,0,31*0/

!      IREFLECT=KREFLECT
!      NXMOL=KXMOL

!print *,'Just entering RRTM_ECRT_140GP KLEV=',KLEV,' IPLON=',IPLON

IREFLECT=0
NXMOL=2

DO J1=1,35
  IXINDX(J1)=0
  DO J2=1,KLEV
    WKL(J1,J2)=_ZERO_
  ENDDO
ENDDO
IXINDX(2)=2
IXINDX(3)=3
DO J1=1,JPXSEC
  DO J2=1,KLEV
    WX(J1,J2)=_ZERO_
  ENDDO
ENDDO

!     Set parameters needed for RRTM execution:
IATM    = 0
!      IXSECT  = 1
!      NUMANGS = 0
!      IOUT    = -1
IXMAX   = 4

!     Bands 6,7,8 are considered the 'window' and allowed to have a
!     different surface emissivity (as in ECMWF).  Eli wrote this part....
SEMISS(1)  = ZEMIS(IPLON)
SEMISS(2)  = ZEMIS(IPLON)
SEMISS(3)  = ZEMIS(IPLON)
SEMISS(4)  = ZEMIS(IPLON)
SEMISS(5)  = ZEMIS(IPLON)
SEMISS(6)  = ZEMIW(IPLON)
SEMISS(7)  = ZEMIW(IPLON)
SEMISS(8)  = ZEMIW(IPLON)
SEMISS(9)  = ZEMIS(IPLON)
SEMISS(10) = ZEMIS(IPLON)
SEMISS(11) = ZEMIS(IPLON)
SEMISS(12) = ZEMIS(IPLON)
SEMISS(13) = ZEMIS(IPLON)
SEMISS(14) = ZEMIS(IPLON)
SEMISS(15) = ZEMIS(IPLON)
SEMISS(16) = ZEMIS(IPLON)

!print *,'after SEMISS'

!     Set surface temperature.  

TBOUND = pts(iplon)
!print *,'after TBOUND=',TBOUND

!     Install ECRT arrays into RRTM arrays for pressure, temperature,
!     and molecular amounts.  Pressures are converted from Pascals
!     (ECRT) to mb (RRTM).  H2O, CO2, O3 and trace gas amounts are 
!     converted from mass mixing ratio to volume mixing ratio.  CO2
!     converted with same dry air and CO2 molecular weights used in 
!     ECRT to assure correct conversion back to the proper CO2 vmr.
!     The dry air column COLDRY (in molec/cm2) is calculated from 
!     the level pressures PZ (in mb) based on the hydrostatic equation
!     and includes a correction to account for H2O in the layer.  The
!     molecular weight of moist air (amm) is calculated for each layer.
!     Note: RRTM levels count from bottom to top, while the ECRT input
!     variables count from the top down and must be reversed here.

NLAYERS = klev
NMOL = 6
PZ(0) = paph(iplon,klev+1)/100._JPRB
TZ(0) = pth(iplon,klev+1)
DO L = 1, KLEV
  PAVEL(L) = pap(iplon,KLEV-L+1)/100._JPRB
  TAVEL(L) = pt(iplon,KLEV-L+1)
  PZ(L) = paph(iplon,KLEV-L+1)/100._JPRB
  TZ(L) = pth(iplon,KLEV-L+1)
  WKL(1,L) = pq(iplon,KLEV-L+1)*amd/amw
  WKL(2,L) = pcco2*amd/amco2
  WKL(3,L) = pozn(iplon,KLEV-L+1)*amd/amo
  WKL(4,L) = rn2o*amd/amn2o
  WKL(6,L) = rch4*amd/amch4
  amm = (1-WKL(1,L))*amd + WKL(1,L)*amw
  COLDRY(L) = (PZ(L-1)-PZ(L))*1.E3_JPRB*avgdro/(gravit*amm*(1+WKL(1,L)))
ENDDO

!print *,'after WKL'
!print 9001,((RAER(JIS,JAE),JAE=1,6),JIS=1,5)
9001 format(1x,6E12.5)


!- Fill RRTM aerosol arrays with operational ECMWF aerosols,
!  do the mixing and distribute over the 16 spectral intervals

DO L=1,KLEV
  JK=KLEV-L+1
!  print 9002,JK,(PAER(IPLON,JK,JAE),JAE=1,6)
9002 format(1x,I3,6E12.5)
  
  
!       DO JAE=1,5
  JAE=1
  ZTAUAER(JAE) =&
   &(RAER(JAE,1)*PAER(IPLON,1,JK)+RAER(JAE,2)*PAER(IPLON,2,JK)&
   &+RAER(JAE,3)*PAER(IPLON,3,JK)+RAER(JAE,4)*PAER(IPLON,4,JK)&
   &+RAER(JAE,5)*PAER(IPLON,5,JK)+RAER(JAE,6)*PAER(IPLON,6,JK))
!   &/(PAPH(IPLON,JK+1)-PAPH(IPLON,JK))
  TAUAERL(L, 1)=ZTAUAER(1)
  TAUAERL(L, 2)=ZTAUAER(1)
  JAE=2
  ZTAUAER(JAE) =&
   &(RAER(JAE,1)*PAER(IPLON,1,JK)+RAER(JAE,2)*PAER(IPLON,2,JK)&
   &+RAER(JAE,3)*PAER(IPLON,3,JK)+RAER(JAE,4)*PAER(IPLON,4,JK)&
   &+RAER(JAE,5)*PAER(IPLON,5,JK)+RAER(JAE,6)*PAER(IPLON,6,JK))
!   &/(PAPH(IPLON,JK+1)-PAPH(IPLON,JK))
  TAUAERL(L, 3)=ZTAUAER(2)
  TAUAERL(L, 4)=ZTAUAER(2)
  TAUAERL(L, 5)=ZTAUAER(2)
  JAE=3
  ZTAUAER(JAE) =&
   &(RAER(JAE,1)*PAER(IPLON,1,JK)+RAER(JAE,2)*PAER(IPLON,2,JK)&
   &+RAER(JAE,3)*PAER(IPLON,3,JK)+RAER(JAE,4)*PAER(IPLON,4,JK)&
   &+RAER(JAE,5)*PAER(IPLON,5,JK)+RAER(JAE,6)*PAER(IPLON,6,JK))
!   &/(PAPH(IPLON,JK+1)-PAPH(IPLON,JK))
  TAUAERL(L, 6)=ZTAUAER(3)
  TAUAERL(L, 8)=ZTAUAER(3)
  TAUAERL(L, 9)=ZTAUAER(3)
  JAE=4
  ZTAUAER(JAE) =&
   &(RAER(JAE,1)*PAER(IPLON,1,JK)+RAER(JAE,2)*PAER(IPLON,2,JK)&
   &+RAER(JAE,3)*PAER(IPLON,3,JK)+RAER(JAE,4)*PAER(IPLON,4,JK)&
   &+RAER(JAE,5)*PAER(IPLON,5,JK)+RAER(JAE,6)*PAER(IPLON,6,JK))
!   &/(PAPH(IPLON,JK+1)-PAPH(IPLON,JK))
  TAUAERL(L, 7)=ZTAUAER(4)
  JAE=5
  ZTAUAER(JAE) =&
   &(RAER(JAE,1)*PAER(IPLON,1,JK)+RAER(JAE,2)*PAER(IPLON,2,JK)&
   &+RAER(JAE,3)*PAER(IPLON,3,JK)+RAER(JAE,4)*PAER(IPLON,4,JK)&
   &+RAER(JAE,5)*PAER(IPLON,5,JK)+RAER(JAE,6)*PAER(IPLON,6,JK))
!   &/(PAPH(IPLON,JK+1)-PAPH(IPLON,JK))
!       END DO
  TAUAERL(L,10)=ZTAUAER(5)
  TAUAERL(L,11)=ZTAUAER(5)
  TAUAERL(L,12)=ZTAUAER(5)
  TAUAERL(L,13)=ZTAUAER(5)
  TAUAERL(L,14)=ZTAUAER(5)
  TAUAERL(L,15)=ZTAUAER(5)
  TAUAERL(L,16)=ZTAUAER(5)
!  print 9003,L,(ZTAUAER(JAE),JAE=1,5)
9003 format(1x,'rrtm_ecrt ZTAUAER:',I3,5e13.6)  
ENDDO

DO L = 1, KLEV
!- Set cross section molecule amounts from ECRT; convert to vmr
  WX(2,L) = rcfc11*amd/amc11
  WX(3,L) = rcfc12*amd/amc12
!-- DS_000515  
END DO

!- Here, all molecules in WKL and WX are in volume mixing ratio; convert to
!  molec/cm2 based on COLDRY for use in RRTM
DO IMOL = 1, NMOL
  DO L = 1, KLEV
!-- DS_000515  
    WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
  END DO  
ENDDO
  
DO IX = 1,JPXSEC
  IF (IXINDX(IX)  /=  0) THEN
!-- DS_000515  
    DO L=1 , KLEV
      WX(IXINDX(IX),L) = COLDRY(L) * WX(IX,L) * 1.E-20_JPRB
    END DO  
  ENDIF
ENDDO


!- Approximate treatment for various cloud overlaps
ZCLEAR=_ONE_
ZCLOUD=_ZERO_
ZC1J(0)=_ZERO_
ZEPSEC=1.E-03_JPRB
JL=IPLON

!++MODIF_MESONH
IF ((NOVLP == 1).OR.(NOVLP ==6).OR.(NOVLP ==8)) THEN
!--MODIF_MESONH

  DO JK=1,KLEV
    IF (pcldf(JL,JK) > ZEPSEC) THEN
      ZCLDLY=pcldf(JL,JK)
      ZCLEAR=ZCLEAR &
       &*(_ONE_-MAX( ZCLDLY , ZCLOUD ))&
       &/(_ONE_-MIN( ZCLOUD , _ONE_-ZEPSEC ))
      ZCLOUD = ZCLDLY
      ZC1J(JK)= _ONE_ - ZCLEAR
    ELSE
      ZCLDLY=_ZERO_
      ZCLEAR=ZCLEAR &
       &*(_ONE_-MAX( ZCLDLY , ZCLOUD ))&
       &/(_ONE_-MIN( ZCLOUD , _ONE_-ZEPSEC ))
      ZCLOUD = ZCLDLY
      ZC1J(JK)= _ONE_ - ZCLEAR
    ENDIF
  ENDDO

!++MODIF_MESONH
ELSEIF ((NOVLP == 2).OR.(NOVLP ==7)) THEN
!--MODIF_MESONH

  DO JK=1,KLEV
    IF (pcldf(JL,JK) > ZEPSEC) THEN
      ZCLDLY=pcldf(JL,JK)
      ZCLOUD = MAX( ZCLDLY , ZCLOUD )
      ZC1J(JK) = ZCLOUD
    ELSE
      ZCLDLY=_ZERO_
      ZCLOUD = MAX( ZCLDLY , ZCLOUD )
      ZC1J(JK) = ZCLOUD
    ENDIF
  ENDDO

!++MODIF_MESONH
ELSEIF ((NOVLP == 3).OR.(NOVLP ==5)) THEN
!--MODIF_MESONH

  DO JK=1,KLEV
    IF (pcldf(JL,JK) > ZEPSEC) THEN
      ZCLDLY=pcldf(JL,JK)
      ZCLEAR = ZCLEAR * (_ONE_-ZCLDLY)
      ZCLOUD = _ONE_ - ZCLEAR
      ZC1J(JK) = ZCLOUD
    ELSE
      ZCLDLY=_ZERO_
      ZCLEAR = ZCLEAR * (_ONE_-ZCLDLY)
      ZCLOUD = _ONE_ - ZCLEAR
      ZC1J(JK) = ZCLOUD
    ENDIF
  ENDDO

ENDIF
PTCLEAR=_ONE_-ZC1J(KLEV)

! Transfer cloud fraction and cloud optical depth to RRTM arrays; 
! invert array index for pcldf to go from bottom to top for RRTM

!- clear-sky column
IF (PTCLEAR  >  _ONE_-ZEPSEC) THEN
  KCLD=0
  DO L = 1, KLEV
    CLDFRAC(L) = _ZERO_
  ENDDO
  DO JB=1,JPBAND
    DO L=1,KLEV
      TAUCLD(L,JB) = _ZERO_
    ENDDO
  ENDDO

ELSE

!- cloudy column
!   The diffusivity factor (Savijarvi, 1997) on the cloud optical 
!   thickness TAUCLD has already been applied in RADLSW

  KCLD=1
  DO L=1,KLEV
    CLDFRAC(L) = pcldf(iplon,L)
  ENDDO
  DO JB=1,JPBAND
    DO L=1,KLEV
      TAUCLD(L,JB) = ptaucld(iplon,L,JB)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE ORRTM_ECRT_140GP

