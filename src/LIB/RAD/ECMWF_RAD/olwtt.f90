!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
!OPTIONS XOPT(HSFUN)
SUBROUTINE OLWTT (KIDIA,KFDIA,KLON, PGA,PGB,PUU, PTT)
!
!**** *LWTT* - LONGWAVE TRANSMISSION FUNCTIONS
!
!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
!     INTERVALS.
!
!**   INTERFACE.
!     ----------
!          *LWTT* IS CALLED FROM *LWVN*, *LWVD*, *LWVB*
!
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KND    :                     ; WEIGHTING INDEX
! PUU    : (KLON,NUA)         ; ABSORBER AMOUNTS
!     ==== OUTPUTS ===
! PTT    : (KLON,NTRA)        ; TRANSMISSION FUNCTIONS
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
!     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
!          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
!          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
!     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
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
!        ORIGINAL : 88-12-15
!
!-----------------------------------------------------------------------
!      IMPLICIT LOGICAL (L)
!
!#include "yoelw.h"
!#include "yoerad.h"
!#include "yoerdu.h"

#include "tsmbkind.h"

USE YOEOLW   , ONLY : NTRA     ,NUA      ,&
            & O1H     , O2H     ,RPIALF0


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLON

!     ------------------------------------------------------------------
!
!*        0.1   ARGUMENTS
!               ---------
!
REAL_B :: PUU(KLON,NUA), PTT(KLON,NTRA) &
     &  ,  PGA(KLON,8,2), PGB(KLON,8,2)
!
!     ------------------------------------------------------------------
!
!*        0.2   LOCAL ARRAYS
!               ------------
!
!     REAL ZXN(KLON,8),ZXD(KLON,8),ZZ(KLON,8)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JA, JL

!     LOCAL REAL SCALARS
REAL_B :: ZA11, ZA12, ZAERCN, ZEU, ZEU10, ZEU11, ZEU12,&
          &ZEU13, ZODH41, ZODH42, ZODN21, ZODN22, ZPU, ZPU10, &
          &ZPU11, ZPU12, ZPU13, ZSQ1, ZSQ2, ZSQH41, &
          &ZSQH42, ZSQN21, ZSQN22, ZTO1, ZTO2, ZTTF11, &
          &ZTTF12, ZUU11, ZUU12, ZUXY, ZVXY, ZX, ZXCH4, &
          &ZXD, ZXN, ZXN2O, ZY, ZYCH4, ZYN2O, ZZ

!     ------------------------------------------------------------------
!#!DIR$ VFUNCTION SQRTHF
!
!*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
!                 -----------------------------------------------
!
DO JA = 1 , 8
  DO JL = KIDIA,KFDIA
    ZZ      =SQRT(PUU(JL,JA))
!    ZXD(JL,1)=PGB( JL, 1,1) + ZZ(JL, 1)*(PGB( JL, 1,2) + ZZ(JL, 1))
!    ZXN(JL,1)=PGA( JL, 1,1) + ZZ(JL, 1)*(PGA( JL, 1,2) )
!    PTT(JL,1)=ZXN(JL,1)/ZXD(JL,1)
    ZXD      =PGB( JL,JA,1) + ZZ       *(PGB( JL,JA,2) + ZZ       )
    ZXN      =PGA( JL,JA,1) + ZZ       *(PGA( JL,JA,2) )
    PTT(JL,JA)=ZXN      /ZXD
  END DO
END DO
!
!     ------------------------------------------------------------------
!
!*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
!                 ---------------------------------------------------
!
DO JL = KIDIA,KFDIA
  PTT(JL, 9) = PTT(JL, 8)
!
!-  CONTINUUM ABSORPTION: E- AND P-TYPE
!
!-   13     350-500 cm-1
!-   10     500-800 cm-1
!-   11     800-970 + 1110-1250 cm-1
!-   12     970-1110 cm-1
!
!  IF (INWCONT.EQ.0) THEN
!- original ECMWF 16r1 coefficients      
    ZPU   = 0.002 * PUU(JL,10)
    ZPU10 = 112. * ZPU
    ZPU11 = 6.25 * ZPU
    ZPU12 = 5.00 * ZPU
    ZPU13 = 80.0 * ZPU
    ZEU   =  PUU(JL,11)
    ZEU10 =  12. * ZEU
    ZEU11 = 6.25 * ZEU
    ZEU12 = 5.00 * ZEU
    ZEU13 = 80.0 * ZEU
!  ELSE IF (INWCONT.EQ.1) THEN
!- coefficients proposed by Giorgetta and Wild
!    ZPU10 =  0.8109 * PUU(JL,10)
!    ZPU11 =  0.0208 * PUU(JL,10)
!    ZPU12 =  0.0106 * PUU(JL,10)
!    ZPU13 = 12.331  * PUU(JL,10)
!    ZEU10 =   47.7  * PUU(JL,11)
!    ZEU11 =    8.31 * PUU(JL,11)
!    ZEU12 =    5.87 * PUU(JL,11)
!    ZEU13 =   209.  * PUU(JL,11)
!  ELSE IF (INWCONT.EQ.2) THEN
!- coefficients adjusted from Clough CKD22 
!c   ZPU10 =  0.00488 * 81.63 * PUU(JL,10)
!c   ZPU11 =  0.00022 *  8.43 * PUU(JL,10)
!c   ZPU12 =  0.00001 *  5.08 * PUU(JL,10)
!c   ZPU13 =  0.03638 * 721.8 * PUU(JL,10)
!    ZPU10 =  0.3981  * PUU(JL,10)
!    ZPU10 =  0.18    * PUU(JL,10)
!    ZPU11 =  0.00127 * PUU(JL,10)
!    ZPU12 =  0.00071 * PUU(JL,10)
!    ZPU13 =  26.26   * PUU(JL,10)
!         
!!    ZEU10 =   81.63  * PUU(JL,11)
!    ZEU10 =   18.    * PUU(JL,11)
!    ZEU11 =    8.43  * PUU(JL,11)
!    ZEU12 =    5.08  * PUU(JL,11)
!    ZEU13 =   721.8  * PUU(JL,11)
!  END IF   
!
!  IF (LNOCONT) THEN
!    ZPU10 = 0. 
!    ZPU11 = 0.
!    ZPU12 = 0. 
!    ZPU13 = 0. 
!
!    ZEU10 = 0. 
!    ZEU11 = 0. 
!    ZEU12 = 0. 
!    ZEU13 = 0. 
!  END IF
!
!-  OZONE ABSORPTION
!
  ZX = PUU(JL,12)
  ZY = PUU(JL,13)
  ZUXY = 4. * ZX * ZX / (RPIALF0 * ZY)
  ZSQ1 = SQRT(1. + O1H * ZUXY ) - 1.
  ZSQ2 = SQRT(1. + O2H * ZUXY ) - 1.
  ZVXY = RPIALF0 * ZY / (2. * ZX)
  ZAERCN = PUU(JL,17) + ZEU12 + ZPU12
  ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
  ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )
      
!  IF (LNOOZON) THEN
!    ZTO1 = EXP( - ZAERCN )
!    ZTO2 = EXP( - ZAERCN )
!  END IF
!
!-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)
!
!* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1
!
  ZXCH4 = PUU(JL,19)
  ZYCH4 = PUU(JL,20)
  ZUXY = 4. * ZXCH4*ZXCH4/(0.103*ZYCH4)
  ZSQH41 = SQRT(1. + 33.7 * ZUXY) - 1.
  ZVXY = 0.103 * ZYCH4 / (2. * ZXCH4)
  ZODH41 = ZVXY * ZSQH41
!
!* N2O IN INTERVAL 800-970 + 1110-1250 CM-1
!
  ZXN2O = PUU(JL,21)
  ZYN2O = PUU(JL,22)
  ZUXY = 4. * ZXN2O*ZXN2O/(0.416*ZYN2O)
  ZSQN21 = SQRT(1. + 21.3 * ZUXY) - 1.
  ZVXY = 0.416 * ZYN2O / (2. * ZXN2O)
  ZODN21 = ZVXY * ZSQN21
!  
!* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1
!
  ZUXY = 4. * ZXCH4*ZXCH4/(0.113*ZYCH4)
  ZSQH42 = SQRT(1. + 400. * ZUXY) - 1.
  ZVXY = 0.113 * ZYCH4 / (2. * ZXCH4)
  ZODH42 = ZVXY * ZSQH42
!
!* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1
!
  ZUXY = 4. * ZXN2O*ZXN2O/(0.197*ZYN2O)
  ZSQN22 = SQRT(1. + 2000. * ZUXY) - 1.
  ZVXY = 0.197 * ZYN2O / (2. * ZXN2O)
  ZODN22 = ZVXY * ZSQN22
!
!* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1
!
  ZA11 = 2. * PUU(JL,23) * 4.404E+05
  ZTTF11 = 1. - ZA11 * 0.003225
!
!* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1
!
  ZA12 = 2. * PUU(JL,24) * 6.7435E+05
  ZTTF12 = 1. - ZA12 * 0.003225
!
      
!  IF (LNOUMG) THEN
!    PTT(JL,7) = 1.
!    PTT(JL,8) = 1.
!    PTT(JL,9) = 1.
!    ZODH41 = 0.
!    ZODH42 = 0.
!    ZODN21 = 0.
!    ZODN22 = 0.
!    ZTTF11 = 1.
!    ZTTF12 = 1.
!  END IF  
       
  ZUU11 = - PUU(JL,15) - ZEU10 - ZPU10
  ZUU12 = - PUU(JL,16) - ZEU11 - ZPU11 - ZODH41 - ZODN21
  PTT(JL,10) = EXP( - PUU(JL,14) )
  PTT(JL,11) = EXP( ZUU11 )
  PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
  PTT(JL,13) = 0.7554 * ZTO1 + 0.2446 * ZTO2
  PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
  PTT(JL,15) = EXP ( - PUU(JL,14) - ZODH42 - ZODN22 )
END DO
!
RETURN
END SUBROUTINE OLWTT
