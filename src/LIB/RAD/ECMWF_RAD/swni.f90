SUBROUTINE SWNI &
 &( KIDIA , KFDIA , KLON  , KLEV , KAER  , KNU &
 &, PAER  , PAKI  , PALBD , PALBP, PCG   , PCLD, PCLEAR &
 &, PDSIG , POMEGA, POZ   , PRMU , PSEC  , PTAU &
 &, PUD   , PWV   , PQS &
 &, PFDOWN, PFUP  , PCDOWN, PCUP , PSUDU2, PDIFF,PDIRF &
!++MODIF_MESONH
&, ODUST,PPIZA_DST,PCGA_DST,PTAUREL_DST )
!--MODIF_MESONH

!**** *SWNI* - SHORTWAVE RADIATION, NEAR-INFRARED SPECTRAL INTERVALS

!     PURPOSE.
!     --------

!          COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE NEAR-INFRARED 
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SWNI* IS CALLED FROM *SW*.


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
!     CONTINUUM SCATTERING
!          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
!     A GREY MOLECULAR ABSORPTION
!          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
!     OF ABSORBERS
!          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
!          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWDE*, *SWTT*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
!        95-12-07   J.-J. MORCRETTE    NEAR-INFRARED SW
!        990128     JJMorcrette        Sunshine duration
!        99-05-25   JJMorcrette        Revised aerosols

!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE OYOERAD   , ONLY : NOVLP 
USE OYOESW    , ONLY : RRAY     ,RSUN     ,RSWCE    ,RSWCP
USE OYOERAD   , ONLY : NSW
USE YOERDU   , ONLY : REPLOG
!
USE MODI_SWCLR
USE MODI_SWR
USE MODI_SWDE
USE MODI_SWTT1
USE MODI_SWTT
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KNU




!#include "yoeaer.h"
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAER(KLON,6,KLEV)    , PAKI(KLON,2,NSW)&
  &,  PALBD(KLON,NSW)      , PALBP(KLON,NSW)&
  &,  PCG(KLON,NSW,KLEV)   , PCLD(KLON,KLEV)&
  &,  PCLEAR(KLON)         , PDSIG(KLON,KLEV)&
  &,  POMEGA(KLON,NSW,KLEV), POZ(KLON,KLEV)&
  &,  PQS(KLON,KLEV)&
  &,  PRMU(KLON)           , PSEC(KLON)&
  &,  PTAU(KLON,NSW,KLEV)  , PUD(KLON,5,KLEV+1)&
  &,  PWV(KLON,KLEV)

REAL_B :: PFDOWN(KLON,KLEV+1)  , PFUP(KLON,KLEV+1)&
  &,  PCDOWN(KLON,KLEV+1)  , PCUP(KLON,KLEV+1)&
  &,  PSUDU2(KLON)         , PDIFF(KLON,KLEV)&
  &,  PDIRF(KLON,KLEV) 

!Quentin
REAL_B :: ZCLDIR 
REAL_B :: ZTA1(KLON)  
  
!++MODIF_MESONH
LOGICAL           :: ODUST                   ! flag for DUST
REAL_B  :: PPIZA_DST(KLON,KLEV)    !wvl ssa dust for current wavelength
REAL_B  :: PCGA_DST(KLON,KLEV)     !wvl assym.fact dust for current wavelength
REAL_B  :: PTAUREL_DST(KLON,KLEV)  !wvl tau_{wvl}/tau_{550} for current wavelength
!--MODIF_MESONH


!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER_M :: IIND2(2), IIND3(3)
REAL_B :: ZCGAZ(KLON,KLEV)  , ZDIFF(KLON)         , ZDIRF(KLON)&
  &,  ZFD(KLON,KLEV+1)  , ZFU(KLON,KLEV+1) &
  &,  ZG(KLON)          , ZGG(KLON)
REAL_B :: ZPIZAZ(KLON,KLEV)&
  &,  ZRAYL(KLON)       , ZRAY1(KLON,KLEV+1)  , ZRAY2(KLON,KLEV+1)&
  &,  ZREF(KLON)        , ZREFZ(KLON,2,KLEV+1)&
  &,  ZRE1(KLON)        , ZRE2(KLON)&
  &,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
  &,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
  &,  ZRL(KLON,8)&
  &,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)  , ZRMUZ(KLON)&
  &,  ZRNEB(KLON)       , ZRUEF(KLON,8)       , ZR1(KLON) &
  &,  ZR2(KLON,2)       , ZR3(KLON,3)         , ZR4(KLON)&
  &,  ZR21(KLON)        , ZR22(KLON)
REAL_B :: ZS(KLON)&
  &,  ZTAUAZ(KLON,KLEV) , ZTO1(KLON)          , ZTR(KLON,2,KLEV+1)&
  &,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
  &,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
  &,  ZTR1(KLON)        , ZTR2(KLON)&
  &,  ZW(KLON)          , ZW1(KLON)           , ZW2(KLON,2)&
  &,  ZW3(KLON,3)       , ZW4(KLON)           , ZW5(KLON)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IABS, IKL, IKM1, JABS, JAJ, JAJP, JK, JKKI,&
             &JKKP4, JKL, JKLP1, JKM1, JL, JN, JN2J, JREF

!     LOCAL REAL SCALARS
REAL_B :: ZAA, ZBB, ZCNEB, ZRE11, ZRKI, ZRMUM1, ZWH2O
REAL_B :: ZB_ODI(KLON)



!     ------------------------------------------------------------------

!*         1.     NEAR-INFRARED SPECTRAL INTERVAL (0.68-4.00 MICRON)
!                 --------------------------------------------------



!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------


DO JL = KIDIA,KFDIA
  ZRMUM1 = _ONE_ - PRMU(JL)
  ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1 &
   &* (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1 &
   &* (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
  ZRAYL(JL) = MAX(ZRAYL(JL),_ZERO_) 
ENDDO


!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------


!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------


!++MODIF_MESONH
   CALL SWCLR &
        &( KIDIA , KFDIA , KLON ,  KLEV , KAER , KNU &
        &, PAER  , PALBP , PDSIG , ZRAYL, PSEC &
        &, ZCGAZ , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
        &, ZRK0  , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
        &, ODUST,PPIZA_DST,PCGA_DST,PTAUREL_DST &
        &)
!--MODIF_MESONH


!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------


CALL SWR &
  &( KIDIA , KFDIA , KLON , KLEV  , KNU &
  &, PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU &
  &, ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2 , ZREFZ, ZRJ  , ZRK, ZRMUE &
  &, ZTAUAZ, ZTRA1 , ZTRA2, ZTRCLD &
  &)

!     ------------------------------------------------------------------

!*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
!                ------------------------------------------------------


JN = 2

DO JABS=1,2


!*         3.1  SURFACE CONDITIONS
!               ------------------


  DO JL = KIDIA,KFDIA
    ZREFZ(JL,2,1) = PALBD(JL,KNU)
    ZREFZ(JL,1,1) = PALBD(JL,KNU)
  ENDDO


!*         3.2  INTRODUCING CLOUD EFFECTS
!               -------------------------


  DO JK = 2 , KLEV+1
    JKM1 = JK - 1
    IKL=KLEV+1-JKM1
    DO JL = KIDIA,KFDIA
      ZRNEB(JL) = PCLD(JL,JKM1)
      IF (JABS == 1.AND. ZRNEB(JL) > _TWO_*REPLOG) THEN
        ZWH2O=MAX(PWV(JL,IKL),REPLOG)
        ZCNEB=MAX(REPLOG,MIN(ZRNEB(JL),_ONE_-REPLOG))
        ZBB=PUD(JL,JABS,JKM1)*PQS(JL,IKL)/ZWH2O
        ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(_ONE_-ZCNEB),REPLOG)
!-- just to test Box-type computations
!        ZAA=PUD(JL,JABS,JKM1)
!        ZBB=ZAA
      ELSE
        ZAA=PUD(JL,JABS,JKM1)
        ZBB=ZAA
      ENDIF
      ZRKI = PAKI(JL,JABS,KNU)
      ZS(JL) = EXP(-ZRKI * ZAA * 1.66_JPRB)
      ZG(JL) = EXP(-ZRKI * ZAA / ZRMUE(JL,JK))
      ZTR1(JL) = _ZERO_
      ZRE1(JL) = _ZERO_
      ZTR2(JL) = _ZERO_
      ZRE2(JL) = _ZERO_

!++MODIF_MESONH
    IF (NOVLP.GE.5)THEN !MESONH VERSION
       ZW(JL) =PCG(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)
       ZTO1(JL) = PTAU(JL,KNU,JKM1)*(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
       ZW(JL) =POMEGA(JL,KNU,JKM1)*(1-ZW(JL))/(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
       ZGG(JL) =PCG(JL,KNU,JKM1)/(1+PCG(JL,KNU,JKM1))
       ZGG(JL)=ZTO1(JL)*ZW(JL)*ZGG(JL)+ZTAUAZ(JL,JKM1)*ZPIZAZ(JL,JKM1)*ZCGAZ(JL,JKM1)
       ZGG(JL)=ZGG(JL)/(ZTO1(JL)*ZW(JL)+ZTAUAZ(JL,JKM1)*ZPIZAZ(JL,JKM1))
       ZB_ODI(JL)=ZTO1(JL) / ZW(JL)&
         &+ ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
     !if g=0 tau/w=tau'/w'
         &+ ZBB * ZRKI
       ZB_ODI(JL)=(1/( (ZTO1(JL) / ZW(JL))&
         &+ (ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)) ))-(1/ZB_ODI(JL))
       ZB_ODI(JL)=((ZTO1(JL) +  ZTAUAZ(JL,JKM1))**2)*ZB_ODI(JL)
       ZW(JL)=ZTO1(JL)*ZW(JL)+ZTAUAZ(JL,JKM1)*ZPIZAZ(JL,JKM1)-ZB_ODI(JL)
       ZTO1(JL) = ZTO1(JL) +  ZTAUAZ(JL,JKM1)
       ZW(JL)=ZW(JL)/ZTO1(JL)
     ELSE !ECMWF VERSION
       ZW(JL)= POMEGA(JL,KNU,JKM1)
       ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)&
         &+ ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
         &+ ZBB * ZRKI
       ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
       ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
       ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
         &+ (_ONE_ - ZR22(JL)) * ZCGAZ(JL,JKM1)
       ZW(JL) = ZR21(JL) / ZTO1(JL)
    ENDIF
!--MODIF_MESONH
       ZREF(JL) = ZREFZ(JL,1,JKM1)
       ZRMUZ(JL) = ZRMUE(JL,JK)
    ENDDO

    CALL SWDE ( KIDIA, KFDIA, KLON &
     &, ZGG  , ZREF , ZRMUZ, ZTO1, ZW &
     &, ZRE1 , ZRE2 , ZTR1 , ZTR2     )

    DO JL = KIDIA,KFDIA

      ZREFZ(JL,2,JK) = (_ONE_-ZRNEB(JL)) * (ZRAY1(JL,JKM1)&
       &+ ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
       &* ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)&
       &+ ZRNEB(JL) * ZRE1(JL)

      ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)&
       &+ (ZTRA1(JL,JKM1)) * ZG(JL) * (_ONE_-ZRNEB(JL))

      ZREFZ(JL,1,JK)=(_ONE_-ZRNEB(JL))*(ZRAY1(JL,JKM1)&
       &+ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
       &/(_ONE_-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)&
       &+ ZRNEB(JL) * ZRE2(JL)

      ZTR(JL,1,JKM1)= ZRNEB(JL) * ZTR2(JL)&
       &+ (ZTRA1(JL,JKM1)/(_ONE_-ZRAY2(JL,JKM1)&
       &* ZREFZ(JL,1,JKM1)))&
       &* ZG(JL) * (_ONE_ -ZRNEB(JL))

    ENDDO
  ENDDO

!*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!               -------------------------------------------------


  DO JREF=1,2

    JN = JN + 1

    DO JL = KIDIA,KFDIA
      ZRJ(JL,JN,KLEV+1) = _ONE_
      ZRK(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
        ZRJ(JL,JN,JKL) = ZRE11
        ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
      ENDDO
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         4.    INVERT GREY AND CONTINUUM FLUXES
!                --------------------------------



!*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
!                ---------------------------------------------


DO JK = 1 , KLEV+1
  DO JAJ = 1 , 5 , 2
    JAJP = JAJ + 1
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
      ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , REPLOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , REPLOG )
    ENDDO
  ENDDO
ENDDO

DO JK = 1 , KLEV+1
  DO JAJ = 2 , 6 , 2
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , REPLOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , REPLOG )
    ENDDO
  ENDDO
ENDDO

!*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
!                 ---------------------------------------------


DO JK = 1 , KLEV+1
  JKKI = 1
  DO JAJ = 1 , 2
    IIND2(1)=JAJ
    IIND2(2)=JAJ
    DO JN = 1 , 2
      JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4

!*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
!                 --------------------------


      DO JL = KIDIA,KFDIA
        ZW2(JL,1) = LOG( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK))/ PAKI(JL,JAJ,KNU)
        ZW2(JL,2) = LOG( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK))/ PAKI(JL,JAJ,KNU)
      ENDDO

!*         4.2.2  TRANSMISSION FUNCTION
!                 ---------------------


      CALL SWTT1 ( KIDIA,KFDIA,KLON, KNU, 2, IIND2 &
       &, ZW2 &
       &, ZR2                              )

      DO JL = KIDIA,KFDIA
        ZRL(JL,JKKI) = ZR2(JL,1)
        ZRUEF(JL,JKKI) = ZW2(JL,1)
        ZRL(JL,JKKP4) = ZR2(JL,2)
        ZRUEF(JL,JKKP4) = ZW2(JL,2)
      ENDDO

      JKKI=JKKI+1
    ENDDO
  ENDDO

!*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
!                 ------------------------------------------------------


  DO JL = KIDIA,KFDIA
    PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
     &+ ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
    PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
     &+ ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
!                ----------------------------------------



!*         5.1   DOWNWARD FLUXES
!                ---------------


JAJ = 2
IIND3(1)=1
IIND3(2)=2
IIND3(3)=3

DO JL = KIDIA,KFDIA
  ZW3(JL,1)=_ZERO_
  ZW3(JL,2)=_ZERO_
  ZW3(JL,3)=_ZERO_
  ZW4(JL)  =_ZERO_
  ZW5(JL)  =_ZERO_
  ZR4(JL)  =_ONE_
  ZFD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1)
ENDDO
DO JK = 1 , KLEV
  IKL = KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
    ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
    ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
    ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
  ENDDO

  CALL SWTT1 ( KIDIA,KFDIA,KLON, KNU, 3, IIND3 &
   &, ZW3 &
   &, ZR3                              )

  DO JL = KIDIA,KFDIA
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
    ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)* ZRJ0(JL,JAJ,IKL)
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZDIFF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLD(JL)
  ZDIRF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLR(JL)
  PSUDU2(JL) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
   &+PCLEAR(JL) * ZDIRF(JL)) * RSUN(KNU)
ENDDO


!*         5.2   UPWARD FLUXES
!                -------------


DO JL = KIDIA,KFDIA
  ZFU(JL,1) = ZFD(JL,1)*PALBP(JL,KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
    ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
    ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKM1)*1.66_JPRB
    ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKM1)*1.66_JPRB
    ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKM1)*1.66_JPRB
  ENDDO

  CALL SWTT1 ( KIDIA,KFDIA,KLON, KNU, 3, IIND3 &
   &, ZW3 &
   &, ZR3                              )

  DO JL = KIDIA,KFDIA
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
    ZFU(JL,JK) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)* ZRK0(JL,JAJ,JK)
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
!                 --------------------------------------------------

IABS=3

!*         6.1    DOWNWARD FLUXES
!                 ---------------

DO JL = KIDIA,KFDIA
  ZW1(JL)=_ZERO_
  ZW4(JL)=_ZERO_
  ZW5(JL)=_ZERO_
  ZR1(JL)=_ZERO_
  PFDOWN(JL,KLEV+1) = ((_ONE_-PCLEAR(JL))*PFDOWN(JL,KLEV+1)&
   &+ PCLEAR(JL) * ZFD(JL,KLEV+1)) * RSUN(KNU)
  PCDOWN(JL,KLEV+1) = ZFD(JL,KLEV+1) * RSUN(KNU)
ENDDO

! Quentin
DO JL = KIDIA,KFDIA
  ZTA1(JL)=_ZERO_
  ZTO1(JL)=_ZERO_
ENDDO
! Quentin

DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
    ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
    ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
  ENDDO

  CALL SWTT ( KIDIA,KFDIA,KLON, KNU, IABS, ZW1, ZR1 )
  
! Quentin
  DO JL = KIDIA,KFDIA
     PFDOWN(JL,IKL) = ((_ONE_-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,&
     &IKL)&
     &+PCLEAR(JL)*ZFD(JL,IKL)) * RSUN(KNU)
     PCDOWN(JL,IKL) = ZFD(JL,IKL) * RSUN(KNU)
     ZTA1(JL)=ZTA1(JL)+ZTAUAZ(JL,IKL)
     ZTO1(JL) = PTAU(JL,KNU,IKL)*(1.-(POMEGA(JL,KNU,IKL)* &
     &           PCG(JL,KNU,IKL)*PCG(JL,KNU,IKL))) + ZTO1(JL)
     ZCLDIR     = ZFD(JL,IKL)/ZRJ0(JL,JAJ,IKL)*EXP(-ZTA1(JL)/PRMU(JL))
     
     PDIRF(JL,IKL) = ((_ONE_-PCLEAR(JL))*ZCLDIR*EXP(-ZTO1(JL)/PRMU(JL)) + &
     &              PCLEAR(JL)*ZCLDIR) * RSUN(KNU)
     PDIRF(JL,IKL) = MIN(PFDOWN(JL,IKL),PDIRF(JL,IKL))
     PDIFF(JL,IKL) = PFDOWN(JL,IKL) - PDIRF(JL,IKL)
    ! PDIFF(JL,IKL)=ZR1(JL)*ZR4(JL)*PFDOWN(JL,IKL)*RSUN(KNU)*&
    ! &              (_ONE_-PCLEAR(JL))
    ! PDIRF(JL,IKL)=ZFD(JL,IKL)*RSUN(KNU)* PCLEAR(JL)
  ENDDO
ENDDO
! 
! PRINT*,"SWNI,PDIFF,PDIRF"
! print*, KLEV,shape(PDIRF),shape(PFDOWN)
! PRINT*,PDIRF(1,:)-PFDOWN(1,1:KLEV)
! PRINT*,PDIFF(:,1)
! PRINT*,PDIRF(:,1)

!*         6.2    UPWARD FLUXES
!                 -------------

DO JL = KIDIA,KFDIA
  PFUP(JL,1) = ((_ONE_-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)&
   &+PCLEAR(JL)*ZFU(JL,1)) * RSUN(KNU)
  PCUP(JL,1) = ZFU(JL,1) * RSUN(KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66_JPRB
    ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66_JPRB
    ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66_JPRB
    ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
  ENDDO

  CALL SWTT ( KIDIA,KFDIA,KLON, KNU, IABS, ZW1, ZR1 )

  DO JL = KIDIA,KFDIA
    PFUP(JL,JK) = ((_ONE_-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)&
     &+PCLEAR(JL)*ZFU(JL,JK)) * RSUN(KNU)
    PCUP(JL,JK) = ZFU(JL,JK) * RSUN(KNU)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWNI




