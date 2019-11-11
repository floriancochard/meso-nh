!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE LWU &
  &( KIDIA, KFDIA, KLON, KLEV &
  &, PAER , PCCO2, PDP , PPMB, PQOF , PTAVE, PVIEW, PWV &
  &, PABCU &
  &)

!**** *LWU* - LONGWAVE EFFECTIVE ABSORBER AMOUNTS

!     PURPOSE.
!     --------
!           COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
!           TEMPERATURE EFFECTS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PCCO2  :                   ; CONCENTRATION IN CO2 (PA/PA)
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS (PA)
! PPMB   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
! PTAVE  : (KLON,KLEV)       ; TEMPERATURE
! PWV    : (KLON,KLEV)       ; SPECIFIC HUMIDITY PA/PA
! PVIEW  : (KLON)            ; COSECANT OF VIEWING ANGLE
!     ==== OUTPUTS ===
! PABCU  :(KLON,NUA,3*KLEV+1); EFFECTIVE ABSORBER AMOUNTS

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
!     ABSORBERS.

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        JJ Morcrette 97-04-18 Revised Continuum + Clean-up

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOMCST   , ONLY : RG
USE OYOESW    , ONLY : RAER
USE YOELW    , ONLY : NSIL     ,NUA      ,NG1      ,NG1P1    ,&
            &ALWT     ,BLWT     ,RO3T     ,RT1      ,TREF     ,&
            &RVGCO2   ,RVGH2O   ,RVGO3
USE OYOERDI   , ONLY : RCH4     ,RN2O     ,RCFC11   ,RCFC12
USE YOERDU   , ONLY : R10E     ,REPSCO   ,REPSCQ


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON

!     DUMMY REAL SCALARS
REAL_B :: PCCO2



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAER(KLON,6,KLEV), PDP(KLON,KLEV)&
  &,  PPMB(KLON,KLEV+1), PQOF(KLON,KLEV)&
  &,  PTAVE(KLON,KLEV) , PVIEW(KLON),  PWV(KLON,KLEV)

REAL_B :: PABCU(KLON,NUA,3*KLEV+1)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------
REAL_B :: ZABLY(KLON,7,3*KLEV+1)  , ZDPM(KLON,3*KLEV)&
  &,  ZDUC(KLON, 3*KLEV+1)    , ZFACT(KLON)&
  &,  ZUPM(KLON,3*KLEV)
REAL_B :: ZPHIO(KLON),ZPSC2(KLON) , ZPSC3(KLON), ZPSH1(KLON)&
  &,  ZPSH2(KLON),ZPSH3(KLON) , ZPSH4(KLON), ZPSH5(KLON)&
  &,  ZPSH6(KLON),ZPSIO(KLON) , ZTCON(KLON)&
  &,  ZPHM6(KLON),ZPSM6(KLON) , ZPHN6(KLON), ZPSN6(KLON)
REAL_B :: ZSSIG(KLON,3*KLEV+1)    , ZTAVI(KLON)&
  &,  ZUAER(KLON,NSIL)        , ZXOZ(KLON) , ZXWV(KLON)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IAE1, IAE2, IAE3, IC, ICP1, IG1, IJ, IJPN,&
             &IKIP1, IKJ, IKJP, IKJPN, IKJR, IKL, JA, JAE, &
             &JK, JKI, JKK, JL

!     LOCAL REAL SCALARS
REAL_B :: ZALUP, ZCAC8, ZCAH1, ZCAH2, ZCAH3, ZCAH4,&
          &ZCAH5, ZCAH6, ZCBC8, ZCBH1, ZCBH2, ZCBH3, &
          &ZCBH4, ZCBH5, ZCBH6, ZDIFF, ZDPMG, ZDPMP0, &
          &ZFPPW, ZTX, ZTX2, ZU6, ZUP, ZUPMCO2, ZUPMG, &
          &ZUPMH2O, ZUPMO3, ZZABLY


!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!-----------------------------------------------------------------------


!*         2.    PRESSURE OVER GAUSS SUB-LEVELS
!                ------------------------------

DO JL = KIDIA,KFDIA
  ZSSIG(JL, 1 ) = PPMB(JL,1) * 100._JPRB
ENDDO

DO JK = 1 , KLEV
  IKJ=(JK-1)*NG1P1+1
  IKJR = IKJ
  IKJP = IKJ + NG1P1
  DO JL = KIDIA,KFDIA
    ZSSIG(JL,IKJP)=PPMB(JL,JK+1)* 100._JPRB
  ENDDO
  DO IG1=1,NG1
    IKJ=IKJ+1
    DO JL = KIDIA,KFDIA
      ZSSIG(JL,IKJ)= (ZSSIG(JL,IKJR) + ZSSIG(JL,IKJP)) * _HALF_ &
       &+ RT1(IG1) * (ZSSIG(JL,IKJP) - ZSSIG(JL,IKJR)) * _HALF_
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------


!*         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
!                --------------------------------------------------

DO JKI=1,3*KLEV
  IKIP1=JKI+1
  DO JL = KIDIA,KFDIA
    ZUPM(JL,JKI)=(ZSSIG(JL,JKI)+ZSSIG(JL,IKIP1))*_HALF_
    ZDPM(JL,JKI)=(ZSSIG(JL,JKI)-ZSSIG(JL,IKIP1))/(10._JPRB*RG)
  ENDDO
ENDDO

DO JK = 1 , KLEV
  IKL = KLEV+1 - JK
  DO JL = KIDIA,KFDIA
    ZXWV(JL) = MAX (PWV(JL,IKL) , REPSCQ )
    ZXOZ(JL) = MAX (PQOF(JL,IKL) / PDP(JL,IKL) , REPSCO )
  ENDDO
  IKJ=(JK-1)*NG1P1+1
  IKJPN=IKJ+NG1
  DO JKK=IKJ,IKJPN
    DO JL = KIDIA,KFDIA
      ZDPMG = ZDPM(JL,JKK)
      ZDPMP0 = ZDPMG / 101325._JPRB
      ZUPMG = ZUPM(JL,JKK) * ZDPMP0
      ZUPMCO2 = ( ZUPM(JL,JKK) + RVGCO2 ) * ZDPMP0
      ZUPMH2O = ( ZUPM(JL,JKK) + RVGH2O ) * ZDPMP0
      ZUPMO3  = ( ZUPM(JL,JKK) + RVGO3  ) * ZDPMP0
      ZDUC(JL,JKK) = ZDPMG
      ZABLY(JL,6,JKK) = ZXOZ(JL) * ZDPMG
      ZABLY(JL,7,JKK) = ZXOZ(JL) * ZUPMO3
      ZU6 = ZXWV(JL) * ZUPMG
      ZFPPW = 1.6078_JPRB * ZXWV(JL) / (_ONE_+0.608_JPRB*ZXWV(JL))
      ZABLY(JL,1,JKK)  = ZXWV(JL) * ZUPMH2O
      ZABLY(JL,5,JKK) = ZU6 * ZFPPW
      ZABLY(JL,4,JKK) = ZU6 * (_ONE_-ZFPPW)
      ZABLY(JL,3,JKK)  = PCCO2 * ZUPMCO2
      ZABLY(JL,2,JKK)  = PCCO2 * ZDPMG
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------


!*         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
!                --------------------------------------------------

DO JA = 1, NUA
  DO JL = KIDIA,KFDIA
    PABCU(JL,JA,3*KLEV+1) = _ZERO_
  ENDDO
ENDDO

DO JK = 1 , KLEV
  IJ=(JK-1)*NG1P1+1
  IJPN=IJ+NG1
  IKL=KLEV+1-JK


!*         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
!               --------------------------------------------------
! --            NB: 'PAER' AEROSOLS ARE ENTERED FROM TOP TO BOTTOM

  IAE1=3*KLEV+1-IJ
  IAE2=3*KLEV+1-(IJ+1)
  IAE3=3*KLEV+1-IJPN
  DO JAE=1,6
    DO JL = KIDIA,KFDIA
      ZUAER(JL,JAE) =&
       &(RAER(JAE,1)*PAER(JL,1,JK)+RAER(JAE,2)*PAER(JL,2,JK)&
       &+RAER(JAE,3)*PAER(JL,3,JK)+RAER(JAE,4)*PAER(JL,4,JK)&
       &+RAER(JAE,5)*PAER(JL,5,JK)+RAER(JAE,6)*PAER(JL,5,JK))&
       &/(ZDUC(JL,IAE1)+ZDUC(JL,IAE2)+ZDUC(JL,IAE3))
    ENDDO
  ENDDO



!*         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
!               --------------------------------------------------

  DO JL = KIDIA,KFDIA
    ZTAVI(JL)=PTAVE(JL,IKL)
    ZFACT(JL)=_ONE_-ZTAVI(JL)/296._JPRB
    ZTCON(JL)=EXP(6.08_JPRB*(296._JPRB/ZTAVI(JL)-_ONE_))
!     ZTCON(JL)=EXP(6.08*ZFACT(JL))
    ZTX=ZTAVI(JL)-TREF
    ZTX2=ZTX*ZTX
    ZZABLY = ZABLY(JL,1,IAE1)+ZABLY(JL,1,IAE2)+ZABLY(JL,1,IAE3)
    ZUP=MIN( MAX( _HALF_*R10E*LOG( ZZABLY ) + 5._JPRB, _ZERO_), 6.0_JPRB)
    ZCAH1=ALWT(1,1)+ZUP*(ALWT(1,2)+ZUP*(ALWT(1,3)))
    ZCBH1=BLWT(1,1)+ZUP*(BLWT(1,2)+ZUP*(BLWT(1,3)))
    ZPSH1(JL)=EXP( ZCAH1 * ZTX + ZCBH1 * ZTX2 )
    ZCAH2=ALWT(2,1)+ZUP*(ALWT(2,2)+ZUP*(ALWT(2,3)))
    ZCBH2=BLWT(2,1)+ZUP*(BLWT(2,2)+ZUP*(BLWT(2,3)))
    ZPSH2(JL)=EXP( ZCAH2 * ZTX + ZCBH2 * ZTX2 )
    ZCAH3=ALWT(3,1)+ZUP*(ALWT(3,2)+ZUP*(ALWT(3,3)))
    ZCBH3=BLWT(3,1)+ZUP*(BLWT(3,2)+ZUP*(BLWT(3,3)))
    ZPSH3(JL)=EXP( ZCAH3 * ZTX + ZCBH3 * ZTX2 )
    ZCAH4=ALWT(4,1)+ZUP*(ALWT(4,2)+ZUP*(ALWT(4,3)))
    ZCBH4=BLWT(4,1)+ZUP*(BLWT(4,2)+ZUP*(BLWT(4,3)))
    ZPSH4(JL)=EXP( ZCAH4 * ZTX + ZCBH4 * ZTX2 )
    ZCAH5=ALWT(5,1)+ZUP*(ALWT(5,2)+ZUP*(ALWT(5,3)))
    ZCBH5=BLWT(5,1)+ZUP*(BLWT(5,2)+ZUP*(BLWT(5,3)))
    ZPSH5(JL)=EXP( ZCAH5 * ZTX + ZCBH5 * ZTX2 )
    ZCAH6=ALWT(6,1)+ZUP*(ALWT(6,2)+ZUP*(ALWT(6,3)))
    ZCBH6=BLWT(6,1)+ZUP*(BLWT(6,2)+ZUP*(BLWT(6,3)))
    ZPSH6(JL)=EXP( ZCAH6 * ZTX + ZCBH6 * ZTX2 )
    ZPHM6(JL)=EXP(-5.81E-4_JPRB * ZTX - 1.13E-6_JPRB * ZTX2 )
    ZPSM6(JL)=EXP(-5.57E-4_JPRB * ZTX - 3.30E-6_JPRB * ZTX2 )
    ZPHN6(JL)=EXP(-3.46E-5_JPRB * ZTX + 2.05E-7_JPRB * ZTX2 )
    ZPSN6(JL)=EXP( 3.70E-3_JPRB * ZTX - 2.30E-6_JPRB * ZTX2 )
  ENDDO

  DO JL = KIDIA,KFDIA
    ZTAVI(JL)=PTAVE(JL,IKL)
    ZTX=ZTAVI(JL)-TREF
    ZTX2=ZTX*ZTX
    ZZABLY = ZABLY(JL,3,IAE1)+ZABLY(JL,3,IAE2)+ZABLY(JL,3,IAE3)
    ZALUP = R10E * LOG ( ZZABLY )
    ZUP   = MAX( _ZERO_ , 5.0_JPRB + _HALF_ * ZALUP )
    ZPSC2(JL) = (ZTAVI(JL)/TREF) ** ZUP
    ZCAC8=ALWT(8,1)+ZUP*(ALWT(8,2)+ZUP*(ALWT(8,3)))
    ZCBC8=BLWT(8,1)+ZUP*(BLWT(8,2)+ZUP*(BLWT(8,3)))
    ZPSC3(JL)=EXP( ZCAC8 * ZTX + ZCBC8 * ZTX2 )
    ZPHIO(JL) = EXP( RO3T(1) * ZTX + RO3T(2) * ZTX2)
    ZPSIO(JL) = EXP( _TWO_* (RO3T(3)*ZTX+RO3T(4)*ZTX2))
  ENDDO

  DO JKK=IJ,IJPN
    IC=3*KLEV+1-JKK
    ICP1=IC+1
    DO JL = KIDIA,KFDIA
      ZDIFF = PVIEW(JL)
!- H2O continuum      
      PABCU(JL,10,IC)=PABCU(JL,10,ICP1)+ ZABLY(JL,4,IC)          *ZDIFF
      PABCU(JL,11,IC)=PABCU(JL,11,ICP1)+ ZABLY(JL,5,IC)*ZTCON(JL)*ZDIFF
!- O3      
      PABCU(JL,12,IC)=PABCU(JL,12,ICP1)+ ZABLY(JL,6,IC)*ZPHIO(JL)*ZDIFF
      PABCU(JL,13,IC)=PABCU(JL,13,ICP1)+ ZABLY(JL,7,IC)*ZPSIO(JL)*ZDIFF
!- CO2
      PABCU(JL,7,IC)=PABCU(JL,7,ICP1)+ ZABLY(JL,3,IC)*ZPSC2(JL)*ZDIFF
      PABCU(JL,8,IC)=PABCU(JL,8,ICP1)+ ZABLY(JL,3,IC)*ZPSC3(JL)*ZDIFF
      PABCU(JL,9,IC)=PABCU(JL,9,ICP1)+ ZABLY(JL,3,IC)*ZPSC3(JL)*ZDIFF
!- H2O
      PABCU(JL,1,IC)=PABCU(JL,1,ICP1)+ ZABLY(JL,1,IC)*ZPSH1(JL)
      PABCU(JL,2,IC)=PABCU(JL,2,ICP1)+ ZABLY(JL,1,IC)*ZPSH2(JL)
      PABCU(JL,3,IC)=PABCU(JL,3,ICP1)+ ZABLY(JL,1,IC)*ZPSH5(JL)*ZDIFF
      PABCU(JL,4,IC)=PABCU(JL,4,ICP1)+ ZABLY(JL,1,IC)*ZPSH3(JL)
      PABCU(JL,5,IC)=PABCU(JL,5,ICP1)+ ZABLY(JL,1,IC)*ZPSH4(JL)
      PABCU(JL,6,IC)=PABCU(JL,6,ICP1)+ ZABLY(JL,1,IC)*ZPSH6(JL)*ZDIFF
!- aerosols
      PABCU(JL,14,IC)=PABCU(JL,14,ICP1)+ ZUAER(JL,1)    *ZDUC(JL,IC)*ZDIFF
      PABCU(JL,15,IC)=PABCU(JL,15,ICP1)+ ZUAER(JL,2)    *ZDUC(JL,IC)*ZDIFF
      PABCU(JL,16,IC)=PABCU(JL,16,ICP1)+ ZUAER(JL,3)    *ZDUC(JL,IC)*ZDIFF
      PABCU(JL,17,IC)=PABCU(JL,17,ICP1)+ ZUAER(JL,4)    *ZDUC(JL,IC)*ZDIFF
      PABCU(JL,18,IC)=PABCU(JL,18,ICP1)+ ZUAER(JL,5)    *ZDUC(JL,IC)*ZDIFF
!- CH4
      PABCU(JL,19,IC)=PABCU(JL,19,ICP1)&
       &+ ZABLY(JL,2,IC)*RCH4/PCCO2*ZPHM6(JL)*ZDIFF
      PABCU(JL,20,IC)=PABCU(JL,20,ICP1)&
       &+ ZABLY(JL,3,IC)*RCH4/PCCO2*ZPSM6(JL)*ZDIFF
!- N2O
      PABCU(JL,21,IC)=PABCU(JL,21,ICP1)&
       &+ ZABLY(JL,2,IC)*RN2O/PCCO2*ZPHN6(JL)*ZDIFF
      PABCU(JL,22,IC)=PABCU(JL,22,ICP1)&
       &+ ZABLY(JL,3,IC)*RN2O/PCCO2*ZPSN6(JL)*ZDIFF
!- CFC11
      PABCU(JL,23,IC)=PABCU(JL,23,ICP1)&
       &+ ZABLY(JL,2,IC)*RCFC11/PCCO2        *ZDIFF
!- CFC12
      PABCU(JL,24,IC)=PABCU(JL,24,ICP1)&
       &+ ZABLY(JL,2,IC)*RCFC12/PCCO2        *ZDIFF
    ENDDO
  ENDDO

ENDDO
!      print *,'END OF LWU'

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE LWU
