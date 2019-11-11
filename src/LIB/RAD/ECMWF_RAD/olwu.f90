!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWU ( KIDIA, KFDIA, KLON, KLEV &
     &  ,  PAER, PCCO2, PDP, PAPH, PQOF, PT, PVIEW, PWV &
     &  ,  PABCU     )
!
!**** *LWU* - LONGWAVE EFFECTIVE ABSORBER AMOUNTS
!
!     PURPOSE.
!     --------
!           COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
!           TEMPERATURE EFFECTS
!
!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PCCO2  :                   ; CONCENTRATION IN CO2 (PA/PA)
! PAPH   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PQOF   : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
! PT     : (KLON,KLEV)       ; TEMPERATURE
! PWV    : (KLON,KLEV)       ; SPECIFIC HUMIDITY PA/PA
! PVIEW  : (KLON)            ; COSECANT OF VIEWING ANGLE
!     ==== OUTPUTS ===
! PABCU  :(KLON,NUA,3*KLEV+1); EFFECTIVE ABSORBER AMOUNTS
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
!     ABSORBERS.
!
!     EXTERNALS.
!     ----------
!
!          NONE
!
!     REFERENCE.
!     ----------
!
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
!
!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOMCST   , ONLY : RG
USE OYOESW    , ONLY : RAER
USE YOEOLW   , ONLY : NISP     ,NUA      ,NG1      ,NG1P1    ,&
            &AT       ,BT       ,RT1      ,TREF     ,OCT      ,&
            &RVGCO2   ,RVGH2O   ,RVGO3
USE OYOERDI   , ONLY : RCH4     ,RN2O     ,RCFC11   ,RCFC12
USE YOERDU   , ONLY : R10E     ,REPSCO   ,REPSCQ
USE YOEDBUG  , ONLY : LDEBUG


IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON

!     DUMMY REAL SCALARS
REAL_B :: PCCO2

CHARACTER*3 CVAR


!-----------------------------------------------------------------------
!
!*       0.1   ARGUMENTS
!              ---------
!
REAL_B :: PAER(KLON,6,KLEV) , PDP(KLON,KLEV) &
     &  ,  PAPH(KLON,KLEV+1), PQOF(KLON,KLEV) &
     &  ,  PT(KLON,KLEV), PVIEW(KLON),  PWV(KLON,KLEV)
!
REAL_B :: PABCU(KLON,NUA,3*KLEV+1)
!
!-----------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
REAL_B :: ZABLY(KLON,7, 3*KLEV),  ZDUC(KLON, 3*KLEV) &
     &  , ZDPM(KLON, 3*KLEV), ZUPM(KLON, 3*KLEV)
REAL_B :: ZPHIO(KLON),ZPSC2(KLON),ZPSC3(KLON),ZPSH1(KLON) &
     &  ,  ZPSH2(KLON),ZPSH3(KLON),ZPSH4(KLON),ZPSH5(KLON) &
     &  ,  ZPSH6(KLON),ZPSIO(KLON),ZTCON(KLON) &
     &  ,  ZPHM6(KLON),ZPSM6(KLON),ZPHN6(KLON),ZPSN6(KLON)
REAL_B :: ZSSIG(KLON,3*KLEV+1), ZTAVI(KLON) &
     &  ,  ZUAER(KLON,6), ZXOZ(KLON), ZXWV(KLON)
!

!     LOCAL INTEGER SCALARS
INTEGER_M :: IAE1, IAE2, IAE3, IC, ICCC, ICP1, IG1, IJ, IJPN,&
             &IKIP1, IKJ, IKJP, IKJPN, IKJR, IKL, IPK, JA, JAE, &
             &JK, JKI, JKK, JL, JAE1, JAE2, JAE3, JC, JCP1, JJ, JJPN, &
             & JKJ, JKJR, JKJP, JKIP1, JKP1, JKJPN

!     LOCAL REAL SCALARS
REAL_B :: ZALUP, ZCAC8, ZCAH1, ZCAH2, ZCAH3, ZCAH4,&
          &ZCAH5, ZCAH6, ZCBC8, ZCBH1, ZCBH2, ZCBH3, &
          &ZCBH4, ZCBH5, ZCBH6, ZDIFF, ZDPMG, ZDPMP0, &
          &ZFPPW, ZTX, ZTX2, ZU6, ZUP, ZUPMCO2, ZUPMG, &
          &ZUPMH2O, ZUPMO3, ZZABLY, ZDPMG, ZUPMG, ZDPMP0

!
!-----------------------------------------------------------------------
!
!*         1.    INITIALIZATION
!                --------------
!
!-----------------------------------------------------------------------
!
!
!*         2.    PRESSURE OVER GAUSS SUB-LEVELS
!                ------------------------------
!
DO JL = KIDIA,KFDIA
  ZSSIG(JL, 1 ) = PAPH(JL,KLEV+1)
END DO
!
DO JK = 1 , KLEV
  JKJ=(JK-1)*NG1P1+1
  JKJR = JKJ
  JKJP = JKJ + NG1P1
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZSSIG(JL,JKJP)=PAPH(JL,IKL)
  END DO
  DO IG1=1,NG1
    JKJ=JKJ+1
    DO JL = KIDIA,KFDIA
      ZSSIG(JL,JKJ)= (ZSSIG(JL,JKJR)+ZSSIG(JL,JKJP))*0.5_JPRB &
      &  + RT1(IG1) * (ZSSIG(JL,JKJP) - ZSSIG(JL,JKJR)) * 0.5_JPRB
    END DO
  END DO
END DO
!
!-----------------------------------------------------------------------
!
!
!*         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
!                --------------------------------------------------
!
DO JKI=1,3*KLEV
  JKIP1=JKI+1
  DO JL = KIDIA,KFDIA
    ZUPM(JL,JKI)=(ZSSIG(JL,JKI)+ZSSIG(JL,JKIP1))*0.5_JPRB
    ZDPM(JL,JKI)=(ZSSIG(JL,JKI)-ZSSIG(JL,JKIP1))/(10._JPRB*RG)
  END DO
END DO
!
DO JK = 1 , KLEV
  JKP1=JK+1
  IKL = KLEV+1 - JK
  DO JL = KIDIA,KFDIA
    ZXWV(JL) = MAX (PWV(JL,IKL) , REPSCQ )
    ZXOZ(JL) = MAX (PQOF(JL,IKL) / PDP(JL,IKL) , REPSCO )
  END DO
  JKJ=(JK-1)*NG1P1+1
  JKJPN=JKJ+NG1
  DO JKK=JKJ,JKJPN
    DO JL = KIDIA,KFDIA
      ZDPMG = ZDPM(JL,JKK)
      ZDPMP0 = ZDPMG / 101325._JPRB
      ZUPMG =   ZUPM(JL,JKK) * ZDPMP0
      ZUPMCO2 = ( ZUPM(JL,JKK) + RVGCO2 ) * ZDPMP0
      ZUPMH2O = ( ZUPM(JL,JKK) + RVGH2O ) * ZDPMP0
      ZUPMO3  = ( ZUPM(JL,JKK) + RVGO3  ) * ZDPMP0
      ZDUC(JL,JKK) = ZDPMG
      ZABLY(JL,6,JKK) = ZXOZ(JL) * ZDPMG
      ZABLY(JL,7,JKK) = ZXOZ(JL) * ZUPMO3
      ZU6 = ZXWV(JL) * ZUPMG
      ZFPPW = 1.6078_JPRB * ZXWV(JL) / (1._JPRB+0.608_JPRB*ZXWV(JL))
      ZABLY(JL,1,JKK) = ZXWV(JL) * ZUPMH2O
      ZABLY(JL,5,JKK) = ZU6 * ZFPPW
      ZABLY(JL,4,JKK) = ZU6 * (1._JPRB-ZFPPW)
      ZABLY(JL,3,JKK) = PCCO2 * ZUPMCO2
      ZABLY(JL,2,JKK) = PCCO2 * ZDPMG
    END DO
  END DO
END DO
!
!-----------------------------------------------------------------------
!
!
!*         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
!                --------------------------------------------------
!
DO JA = 1, NUA
  DO JL = KIDIA,KFDIA
    PABCU(JL,JA,3*KLEV+1) = 0.
  END DO
END DO
!

DO JK = 1 , KLEV
  JJ=(JK-1)*NG1P1+1
  JJPN=JJ+NG1
  IKL=KLEV+1-JK
!
!
!*         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
!               --------------------------------------------------
! --            NB: 'PAER' AEROSOLS ARE ENTERED FROM TOP TO BOTTOM
!
  JAE1=3*KLEV+1-JJ
  JAE2=3*KLEV+1-(JJ+1)
  JAE3=3*KLEV+1-JJPN
  DO JAE=1,5
    DO JL = KIDIA,KFDIA
      ZUAER(JL,JAE) = (RAER(JAE,1)*PAER(JL, 1, JK ) &
      &      +RAER(JAE,2)*PAER(JL, 2, JK)+RAER(JAE,3)*PAER(JL, 3, JK) &
      &      +RAER(JAE,4)*PAER(JL, 4, JK)+RAER(JAE,5)*PAER(JL, 5, JK)) &
      &      /(ZDUC(JL,JAE1)+ZDUC(JL,JAE2)+ZDUC(JL,JAE3))
    END DO
  END DO
!
!
!*         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
!               --------------------------------------------------
! --            NB: 'PT' TEMPERATURES ARE ENTERED FROM TOP TO BOTTOM
!
  DO JL = KIDIA,KFDIA
  
    ZTAVI(JL)=PT(JL,JK)
    ZTCON(JL)=EXP(6.08_JPRB*(296._JPRB/ZTAVI(JL)-1._JPRB))
    ZTX=ZTAVI(JL)-TREF
    ZTX2=ZTX*ZTX
    ZZABLY = ZABLY(JL,3,JAE1)+ZABLY(JL,3,JAE2)+ZABLY(JL,3,JAE3)
    ZUP=MIN( MAX( 0.5_JPRB*R10E*LOG( ZZABLY ) + 5._JPRB, 0._JPRB), 6.0_JPRB)
    ZCAH1=AT(1,1)+ZUP*(AT(1,2)+ZUP*(AT(1,3)))
    ZCBH1=BT(1,1)+ZUP*(BT(1,2)+ZUP*(BT(1,3)))
    ZPSH1(JL)=EXP( ZCAH1 * ZTX + ZCBH1 * ZTX2 )
    ZCAH2=AT(2,1)+ZUP*(AT(2,2)+ZUP*(AT(2,3)))
    ZCBH2=BT(2,1)+ZUP*(BT(2,2)+ZUP*(BT(2,3)))
    ZPSH2(JL)=EXP( ZCAH2 * ZTX + ZCBH2 * ZTX2 )
    ZCAH3=AT(3,1)+ZUP*(AT(3,2)+ZUP*(AT(3,3)))
    ZCBH3=BT(3,1)+ZUP*(BT(3,2)+ZUP*(BT(3,3)))
    ZPSH3(JL)=EXP( ZCAH3 * ZTX + ZCBH3 * ZTX2 )
    ZCAH4=AT(4,1)+ZUP*(AT(4,2)+ZUP*(AT(4,3)))
    ZCBH4=BT(4,1)+ZUP*(BT(4,2)+ZUP*(BT(4,3)))
    ZPSH4(JL)=EXP( ZCAH4 * ZTX + ZCBH4 * ZTX2 )
    ZCAH5=AT(5,1)+ZUP*(AT(5,2)+ZUP*(AT(5,3)))
    ZCBH5=BT(5,1)+ZUP*(BT(5,2)+ZUP*(BT(5,3)))
    ZPSH5(JL)=EXP( ZCAH5 * ZTX + ZCBH5 * ZTX2 )
    ZCAH6=AT(6,1)+ZUP*(AT(6,2)+ZUP*(AT(6,3)))
    ZCBH6=BT(6,1)+ZUP*(BT(6,2)+ZUP*(BT(6,3)))
    ZPSH6(JL)=EXP( ZCAH6 * ZTX + ZCBH6 * ZTX2 )
    ZPHM6(JL)=EXP(-5.81E-4_JPRB * ZTX - 1.13E-6 * ZTX2 )
    ZPSM6(JL)=EXP(-5.57E-4_JPRB * ZTX - 3.30E-6 * ZTX2 )
    ZPHN6(JL)=EXP(-3.46E-5_JPRB * ZTX + 2.05E-7 * ZTX2 )
    ZPSN6(JL)=EXP( 3.70E-3_JPRB * ZTX - 2.30E-6 * ZTX2 )
  END DO
  
!
  DO JL = KIDIA,KFDIA
  
    ZTAVI(JL)=PT(JL,JK)
    ZTX=ZTAVI(JL)-TREF
    ZTX2=ZTX*ZTX
    ZZABLY = ZABLY(JL,5,JAE1)+ZABLY(JL,5,JAE2)+ZABLY(JL,5,JAE3)
    ZALUP = R10E * LOG ( ZZABLY )
    ZUP   = MAX( 0.0_JPRB , 5.0_JPRB + 0.5_JPRB * ZALUP )
    ZPSC2(JL) = (ZTAVI(JL)/TREF) ** ZUP
    ZCAC8=AT(8,1)+ZUP*(AT(8,2)+ZUP*(AT(8,3)))
    ZCBC8=BT(8,1)+ZUP*(BT(8,2)+ZUP*(BT(8,3)))
    ZPSC3(JL)=EXP( ZCAC8 * ZTX + ZCBC8 * ZTX2 )
    ZPHIO(JL) = EXP( OCT(1) * ZTX + OCT(2) * ZTX2)
    ZPSIO(JL) = EXP( 2.* (OCT(3)*ZTX+OCT(4)*ZTX2))
  END DO
  

  ICCC=2
  
  DO JKK=JJ,JJPN
    JC=3*KLEV+1-JKK
    JCP1=JC+1
    
    ICCC=ICCC+1
    
    DO JL = KIDIA,KFDIA
      ZDIFF = PVIEW(JL)
      PABCU(JL,10,JC)=PABCU(JL,10,JCP1)+ ZABLY(JL,4,JC)          *ZDIFF
      PABCU(JL,11,JC)=PABCU(JL,11,JCP1)+ ZABLY(JL,5,JC)*ZTCON(JL)*ZDIFF
!
      PABCU(JL,12,JC)=PABCU(JL,12,JCP1)+ ZABLY(JL,6,JC)*ZPHIO(JL)*ZDIFF
      PABCU(JL,13,JC)=PABCU(JL,13,JCP1)+ ZABLY(JL,7,JC)*ZPSIO(JL)*ZDIFF
!
      PABCU(JL, 7,JC)=PABCU(JL, 7,JCP1)+ ZABLY(JL,3,JC)*ZPSC2(JL)*ZDIFF
      PABCU(JL, 8,JC)=PABCU(JL, 8,JCP1)+ ZABLY(JL,3,JC)*ZPSC3(JL)*ZDIFF
      PABCU(JL, 9,JC)=PABCU(JL, 9,JCP1)+ ZABLY(JL,3,JC)*ZPSC3(JL)*ZDIFF
!
      PABCU(JL, 1,JC)=PABCU(JL, 1,JCP1)+ ZABLY(JL,1,JC)*ZPSH1(JL)*ZDIFF
      PABCU(JL, 2,JC)=PABCU(JL, 2,JCP1)+ ZABLY(JL,1,JC)*ZPSH2(JL)*ZDIFF
      PABCU(JL, 3,JC)=PABCU(JL, 3,JCP1)+ ZABLY(JL,1,JC)*ZPSH5(JL)*ZDIFF
      PABCU(JL, 4,JC)=PABCU(JL, 4,JCP1)+ ZABLY(JL,1,JC)*ZPSH3(JL)*ZDIFF
      PABCU(JL, 5,JC)=PABCU(JL, 5,JCP1)+ ZABLY(JL,1,JC)*ZPSH4(JL)*ZDIFF
      PABCU(JL, 6,JC)=PABCU(JL, 6,JCP1)+ ZABLY(JL,1,JC)*ZPSH6(JL)*ZDIFF
!
      PABCU(JL,14,JC)=PABCU(JL,14,JCP1)+ ZUAER(JL,1) *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,15,JC)=PABCU(JL,15,JCP1)+ ZUAER(JL,2) *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,16,JC)=PABCU(JL,16,JCP1)+ ZUAER(JL,3) *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,17,JC)=PABCU(JL,17,JCP1)+ ZUAER(JL,4) *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,18,JC)=PABCU(JL,18,JCP1)+ ZUAER(JL,5) *ZDUC(JL,JC)*ZDIFF
! 
      PABCU(JL,19,JC)=PABCU(JL,19,JCP1) &
      &              +ZABLY(JL,2,JC)*RCH4/PCCO2*ZPHM6(JL)*ZDIFF
      PABCU(JL,20,JC)=PABCU(JL,20,JCP1) &
      &              +ZABLY(JL,3,JC)*RCH4/PCCO2*ZPSM6(JL)*ZDIFF
      PABCU(JL,21,JC)=PABCU(JL,21,JCP1) &
      &              +ZABLY(JL,2,JC)*RN2O/PCCO2*ZPHN6(JL)*ZDIFF
      PABCU(JL,22,JC)=PABCU(JL,22,JCP1) &
      &              +ZABLY(JL,3,JC)*RN2O/PCCO2*ZPSN6(JL)*ZDIFF
!
      PABCU(JL,23,JC)=PABCU(JL,23,JCP1) &
      &               +ZABLY(JL,2,JC)*RCFC11/PCCO2       *ZDIFF
      PABCU(JL,24,JC)=PABCU(JL,24,JCP1) &
      &               +ZABLY(JL,2,JC)*RCFC12/PCCO2       *ZDIFF
      
    END DO
  END DO
!
END DO

!
RETURN
END SUBROUTINE OLWU
