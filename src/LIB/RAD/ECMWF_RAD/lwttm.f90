!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
!OPTIONS XOPT(HSFUN)
SUBROUTINE LWTTM ( KIDIA, KFDIA, KLON, PGA  , PGB, PUU1 , PUU2 , PTT        )

!**** *LWTTM* - LONGWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
!     INTERVALS.

!**   INTERFACE.
!     ----------
!          *LWTTM* IS CALLED FROM *LWVD*


!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PGA, PGB                    ; PADE APPROXIMANTS
! PUU1   : (KLON,NUA)         ; ABSORBER AMOUNTS FROM TOP TO LEVEL 1
! PUU2   : (KLON,NUA)         ; ABSORBER AMOUNTS FROM TOP TO LEVEL 2  
!     ==== OUTPUTS ===
! PTT    : (KLON,NTRA)        ; TRANSMISSION FUNCTIONS

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
!     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
!          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
!          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
!     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.

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
!        ORIGINAL : 88-12-15
!        97-04-18 JJ Morcrette        Revised continuum

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : NTRA     ,NUA      ,RPTYPE   ,RETYPE   ,&
            &RO1H     ,RO2H     ,RPIALF0


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLON



!     ------------------------------------------------------------------

!*        0.1   ARGUMENTS
!               ---------

REAL_B :: PUU1(KLON,NUA), PUU2(KLON,NUA), PTT(KLON,NTRA)&
  &,  PGA(KLON,8,2) , PGB(KLON,8,2)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JA, JL

!     LOCAL REAL SCALARS
REAL_B :: ZA11, ZA12, ZAERCN, ZEU, ZEU10, ZEU11, ZEU12,&
          &ZEU13, ZODH41, ZODH42, ZODN21, ZODN22, ZPU, &
          &ZPU10, ZPU11, ZPU12, ZPU13, ZSQ1, ZSQ2, ZSQH41, &
          &ZSQH42, ZSQN21, ZSQN22, ZTO1, ZTO2, ZTTF11, &
          &ZTTF12, ZUU11, ZUU12, ZUXY, ZVXY, ZX, ZXCH4, &
          &ZXD, ZXN, ZXN2O, ZY, ZYCH4, ZYN2O, ZZ


!     ------------------------------------------------------------------
!DIR$ VFUNCTION SQRTHF


!*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
!                 -----------------------------------------------

DO JA = 1 , 8
  DO JL = KIDIA,KFDIA
    ZZ  = SQRT(PUU1(JL,JA) - PUU2(JL,JA))
    ZXD = PGB( JL,JA,1) + ZZ * (PGB( JL,JA,2) + ZZ )
    ZXN = PGA( JL,JA,1) + ZZ * (PGA( JL,JA,2) )
    PTT(JL,JA) = ZXN / ZXD
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  PTT(JL,3)=MAX(PTT(JL,3),_ZERO_)
ENDDO
!     ------------------------------------------------------------------

!*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
!                 ---------------------------------------------------

DO JL = KIDIA,KFDIA
  PTT(JL, 9) = PTT(JL, 8)

!-  CONTINUUM ABSORPTION: E- AND P-TYPE

  ZPU   = (PUU1(JL,10) - PUU2(JL,10))
  ZPU10 = RPTYPE(1) * ZPU
  ZPU11 = RPTYPE(2) * ZPU
  ZPU12 = RPTYPE(3) * ZPU
  ZPU13 = RPTYPE(4) * ZPU
  ZEU   = (PUU1(JL,11) - PUU2(JL,11))
  ZEU10 = RETYPE(1) * ZEU
  ZEU11 = RETYPE(2) * ZEU
  ZEU12 = RETYPE(3) * ZEU
  ZEU13 = RETYPE(4) * ZEU

!-  OZONE ABSORPTION

  ZX = (PUU1(JL,12) - PUU2(JL,12))
  ZY = (PUU1(JL,13) - PUU2(JL,13))
  ZUXY = 4._JPRB * ZX * ZX / (RPIALF0 * ZY)
  ZSQ1 = SQRT(_ONE_ + RO1H * ZUXY ) - _ONE_
  ZSQ2 = SQRT(_ONE_ + RO2H * ZUXY ) - _ONE_
  ZVXY = RPIALF0 * ZY / (_TWO_ * ZX)
  ZAERCN = (PUU1(JL,17) -PUU2(JL,17)) + ZEU12 + ZPU12
  ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
  ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )

!-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)

!* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1

  ZXCH4 = (PUU1(JL,19) - PUU2(JL,19))
  ZYCH4 = (PUU1(JL,20) - PUU2(JL,20))
  ZUXY = 4._JPRB * ZXCH4*ZXCH4/(0.103_JPRB*ZYCH4)
  ZSQH41 = SQRT(_ONE_ + 33.7_JPRB * ZUXY) - _ONE_
  ZVXY = 0.103_JPRB * ZYCH4 / (_TWO_ * ZXCH4)
  ZODH41 = ZVXY * ZSQH41

!* N2O IN INTERVAL 800-970 + 1110-1250 CM-1

  ZXN2O = (PUU1(JL,21) - PUU2(JL,21))
  ZYN2O = (PUU1(JL,22) - PUU2(JL,22))
  ZUXY = 4._JPRB * ZXN2O*ZXN2O/(0.416_JPRB*ZYN2O)
  ZSQN21 = SQRT(_ONE_ + 21.3_JPRB * ZUXY) - _ONE_
  ZVXY = 0.416_JPRB * ZYN2O / (_TWO_ * ZXN2O)
  ZODN21 = ZVXY * ZSQN21

!* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1

  ZUXY = 4._JPRB * ZXCH4*ZXCH4/(0.113_JPRB*ZYCH4)
  ZSQH42 = SQRT(_ONE_ + 400._JPRB * ZUXY) - _ONE_
  ZVXY = 0.113_JPRB * ZYCH4 / (_TWO_ * ZXCH4)
  ZODH42 = ZVXY * ZSQH42

!* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1

  ZUXY = 4._JPRB * ZXN2O*ZXN2O/(0.197_JPRB*ZYN2O)
  ZSQN22 = SQRT(_ONE_ + 2000._JPRB * ZUXY) - _ONE_
  ZVXY = 0.197_JPRB * ZYN2O / (_TWO_ * ZXN2O)
  ZODN22 = ZVXY * ZSQN22

!* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1

  ZA11 = (PUU1(JL,23) - PUU2(JL,23)) * 4.404E+05_JPRB
  ZTTF11 = _ONE_ - ZA11 * 0.003225_JPRB

!* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1

  ZA12 = (PUU1(JL,24) - PUU2(JL,24)) * 6.7435E+05_JPRB
  ZTTF12 = _ONE_ - ZA12 * 0.003225_JPRB

  ZUU11 = - (PUU1(JL,15) - PUU2(JL,15)) - ZEU10 - ZPU10
  ZUU12 = - (PUU1(JL,16) - PUU2(JL,16)) - ZEU11 - ZPU11 -ZODH41 - ZODN21
  PTT(JL,10) = EXP( - (PUU1(JL,14)- PUU2(JL,14)) )
  PTT(JL,11) = EXP( ZUU11 )
  PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
  PTT(JL,13) = 0.7554_JPRB * ZTO1 + 0.2446_JPRB * ZTO2
  PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
  PTT(JL,15) = EXP( - (PUU1(JL,14) - PUU2(JL,14)) - ZODH42-ZODN22 )

ENDDO

RETURN
END SUBROUTINE LWTTM
