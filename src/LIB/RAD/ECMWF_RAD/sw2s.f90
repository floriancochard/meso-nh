!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:46
!-----------------------------------------------------------------
SUBROUTINE SW2S ( KIDIA, KFDIA, KLON , KLEV , KAER, KNU &
  &,  PAER  ,PAKI, PALBD, PALBP, PCG   , PCLD, PCLEAR, PCLDSW &
  &,  PDSIG ,POMEGA,POZ , PRMU , PSEC  , PTAU &
  &,  PUD   ,PWV , PQS &
  &,  PFDOWN,PFUP,PFDOWNC,PFUPC                               )

!**** *SW2* - SHORTWAVE RADIATION, 2ND SPECTRAL INTERVAL

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
!     SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW2S* IS CALLED FROM *SW*.


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
!        96-05-30 Michel Deque (security in EXP()) 
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE OYOESW    , ONLY : RRAY     ,RSUN
USE YOERDU   , ONLY : REPLOG


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

REAL_B :: PAER(KLON,KLEV,5), PAKI(KLON,2)&
  &,  PALBD(KLON,2)     , PALBP(KLON,2)&
  &,  PCG(KLON,2,KLEV) , PCLD(KLON,KLEV)&
  &,  PCLDSW(KLON,KLEV)&
  &,  PCLEAR(KLON)      , PDSIG(KLON,KLEV)&
  &,  POMEGA(KLON,2,KLEV),POZ(KLON,KLEV)&
  &,  PQS(KLON,KLEV)&
  &,  PRMU(KLON)        , PSEC(KLON)&
  &,  PTAU(KLON,2,KLEV), PUD(KLON,5,KLEV+1)&
  &,  PWV(KLON,KLEV)

REAL_B :: PFDOWN(KLON,KLEV+1),PFUP(KLON,KLEV+1)
REAL_B :: PFDOWNC(KLON,KLEV+1),PFUPC(KLON,KLEV+1)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER_M :: IIND2(2), IIND3(3)
REAL_B :: ZCGAZ(KLON,KLEV)&
  &,  ZFD(KLON,KLEV+1), ZFU(KLON,KLEV+1) &
  &,  ZG(KLON), ZGG(KLON)&
  &,  ZPIZAZ(KLON,KLEV)&
  &,  ZRAYL(KLON), ZRAY1(KLON,KLEV+1)&
  &,  ZRAY2(KLON,KLEV+1), ZREF(KLON), ZREFZ(KLON,2,KLEV+1)&
  &,  ZRE1(KLON), ZRE2(KLON)&
  &,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
  &,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
  &,  ZRL(KLON,8)&
  &,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)&
  &,  ZRMUZ(KLON), ZRNEB(KLON),  ZRUEF(KLON,8)&
  &,  ZR1(KLON) , ZR2(KLON,2), ZR3(KLON,3), ZR4(KLON)&
  &,  ZR21(KLON), ZR22(KLON)&
  &,  ZS(KLON)&
  &,  ZTAUAZ(KLON,KLEV), ZTO1(KLON), ZTR(KLON,2,KLEV+1)&
  &,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
  &,  ZTR1(KLON), ZTR2(KLON)&
  &,  ZW(KLON)   , ZW1(KLON), ZW2(KLON,2)&
  &,  ZW3(KLON,3), ZW4(KLON), ZW5(KLON)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IABS, IKL, IKM1, JABS, JAJ, JAJP, JK, JKKI,&
             &JKKP4, JKL, JKLP1, JKM1, JL, JN, JN2J, JREF

!     LOCAL REAL SCALARS
REAL_B :: ZAA, ZBB, ZCNEB, ZRE11, ZRKI, ZRMUM1, ZWH2O



!     ------------------------------------------------------------------

!*         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
!                 -------------------------------------------



!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------


DO JL = KIDIA,KFDIA
  ZRMUM1 = _ONE_ - PRMU(JL)
  ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1 &
   &* (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1 &
   &* (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
ENDDO


!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------


!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------


CALL SWCLR ( KIDIA , KFDIA , KLON, KLEV, KAER, KNU &
  &, PAER   , PALBP , PDSIG , ZRAYL, PSEC &
  &, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
  &, ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2           )


!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------


CALL SWR  ( KIDIA , KFDIA, KLON, KLEV , KAER , KNU &
  &, PALBD , PCG   , PCLD , PDSIG, POMEGA, PSEC , PTAU &
  &, ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2, ZREFZ , ZRJ  , ZRK, ZRMUE &
  &, ZTAUAZ, ZTRA1 , ZTRA2 )


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
      ELSE
        ZAA=PUD(JL,JABS,JKM1)
        ZBB=ZAA
      ENDIF
      ZRKI = PAKI(JL,JABS)
      ZS(JL) = EXP(MIN(200._JPRB,-ZRKI * ZAA * 1.66_JPRB))
      ZG(JL) = EXP(MIN(200._JPRB,-ZRKI * ZAA / ZRMUE(JL,JK)))
      ZTR1(JL) = _ZERO_
      ZRE1(JL) = _ZERO_
      ZTR2(JL) = _ZERO_
      ZRE2(JL) = _ZERO_

      ZW(JL)= POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)&
       &+ ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
       &+ ZBB * ZRKI

      ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
       &+ (_ONE_ - ZR22(JL)) * ZCGAZ(JL,JKM1)
      ZW(JL) = ZR21(JL) / ZTO1(JL)
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
        ZW2(JL,1) = LOG( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK))/ PAKI(JL,JAJ)
        ZW2(JL,2) = LOG( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK))/ PAKI(JL,JAJ)
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
!     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
    ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)* ZRJ0(JL,JAJ,IKL)
  ENDDO
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
!     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
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
  PFDOWNC(JL,KLEV+1) = ZFD(JL,KLEV+1) * RSUN(KNU)
ENDDO

DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
    ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
    ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
!     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
  ENDDO

  CALL SWTT ( KIDIA,KFDIA,KLON, KNU, IABS, ZW1, ZR1 )

  DO JL = KIDIA,KFDIA
    PFDOWN(JL,IKL) = ((_ONE_-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,&
     &IKL)&
     &+PCLEAR(JL)*ZFD(JL,IKL)) * RSUN(KNU)
    PFDOWNC(JL,IKL) = ZFD(JL,IKL) * RSUN(KNU)
  ENDDO
ENDDO


!*         6.2    UPWARD FLUXES
!                 -------------

DO JL = KIDIA,KFDIA
  PFUP(JL,1) = ((_ONE_-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)&
   &+PCLEAR(JL)*ZFU(JL,1)) * RSUN(KNU)
  PFUPC(JL,1) = ZFU(JL,1) * RSUN(KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66_JPRB
    ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66_JPRB
    ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66_JPRB
!     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
  ENDDO

  CALL SWTT ( KIDIA,KFDIA,KLON, KNU, IABS, ZW1, ZR1 )

  DO JL = KIDIA,KFDIA
    PFUP(JL,JK) = ((_ONE_-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)&
     &+PCLEAR(JL)*ZFU(JL,JK)) * RSUN(KNU)
    PFUPC(JL,JK) = ZFU(JL,JK) * RSUN(KNU)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SW2S
