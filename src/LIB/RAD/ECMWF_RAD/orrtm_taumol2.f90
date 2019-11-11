!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL2 (KLEV,TAU,COLDRY,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
  &COLH2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up
!     JJMorcrette 2000-07-14 bugfix


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG2   ,NGS1
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA2 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &FORREF   ,KA     ,KB     ,SELFREF , REFPARAM
USE OYOERRTBG2, ONLY : CORR1  ,CORR2

!  Input
!#include "yoeratm.h"


IMPLICIT NONE

REAL_B :: COLDRY(JPLAY)

!  Output
REAL_B :: TAU   (JPGPT,JPLAY)

!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from AER
REAL_B :: TAUAERL(JPLAY,JPBAND)

!- from INTFAC      
REAL_B :: FAC00(JPLAY)
REAL_B :: FAC01(JPLAY)
REAL_B :: FAC10(JPLAY)
REAL_B :: FAC11(JPLAY)
REAL_B :: FORFAC(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY)
INTEGER_M :: JT(JPLAY)
INTEGER_M :: JT1(JPLAY)

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

REAL_B :: FC00(JPLAY),FC01(JPLAY),FC10(JPLAY),FC11(JPLAY)
!      REAL TAUAER(JPLAY)
REAL_B :: FRACINT(JPLAY)
INTEGER_M :: IND0(JPLAY),IND1(JPLAY),INDS(JPLAY), INDEX(JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IFP, IFRAC, IG,  JFRAC, LAY

!     LOCAL REAL SCALARS
REAL_B :: FP, H2OPARAM, WATER

!      EQUIVALENCE (TAUAERL(1,2),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum is 
!     interpolated (in temperature) separately.

DO LAY = 1, LAYTROP
  WATER = 1.E20_JPRB * COLH2O(LAY) / COLDRY(LAY)
  H2OPARAM = WATER/(WATER +.002_JPRB)
  
!  DO IFRAC = 2, 12
!    IF (H2OPARAM  >=  REFPARAM(IFRAC)) GO TO 1900
!  ENDDO
!  1900 CONTINUE
!  FRACINT(LAY) = (H2OPARAM-REFPARAM(IFRAC))/&
!   &(REFPARAM(IFRAC-1)-REFPARAM(IFRAC))

  IF (H2OPARAM >= REFPARAM(2)) THEN
    INDEX(LAY)=2
  ELSE
    DO JFRAC = 2, 12
      IF (H2OPARAM < REFPARAM(JFRAC)) THEN
        INDEX(LAY)=JFRAC+1
      END IF  
    ENDDO
  ENDIF  
  
!---- JJM_000714
  IFRAC=INDEX(LAY)
  FRACINT(LAY) = (H2OPARAM-REFPARAM(IFRAC))/&
   &(REFPARAM(IFRAC-1)-REFPARAM(IFRAC))
ENDDO

DO LAY = 1, LAYTROP

  FP = FAC11(LAY) + FAC01(LAY)
  IFP = 2.E2_JPRB*FP+_HALF_

!---MI 981104        
!       IF (IFP.LE.0) IFP=0

  IFP=MAX(0,IFP)

  FC00(LAY) = FAC00(LAY) * CORR2(IFP)
  FC10(LAY) = FAC10(LAY) * CORR2(IFP)
  FC01(LAY) = FAC01(LAY) * CORR1(IFP)
  FC11(LAY) = FAC11(LAY) * CORR1(IFP)
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(2) + 1
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(2) + 1
  INDS(LAY) = INDSELF(LAY)
ENDDO

!-- DS_000515  
DO IG = 1, NG2
  DO LAY = 1, LAYTROP
!-- JJM_000714
    IFRAC=INDEX(LAY)  
!-- DS_000515  
    TAU (NGS1+IG,LAY) = COLH2O(LAY) *&
     &(FC00(LAY) * ABSA(IND0(LAY)  ,IG) +&
     & FC10(LAY) * ABSA(IND0(LAY)+1,IG) +&
     & FC01(LAY) * ABSA(IND1(LAY)  ,IG) +&
     & FC11(LAY) * ABSA(IND1(LAY)+1,IG) +&
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG)))&
     &+ FORFAC(LAY) * FORREF(IG) ) &
     &+ TAUAERL(LAY,2)
    PFRAC(NGS1+IG,LAY) = FRACREFA(IG,IFRAC) + FRACINT(LAY) *&
     &(FRACREFA(IG,IFRAC-1)-FRACREFA(IG,IFRAC))
  ENDDO
ENDDO

DO LAY = LAYTROP+1, KLEV
  FP = FAC11(LAY) + FAC01(LAY)
  IFP = 2.E2_JPRB*FP+_HALF_

!---MI 981104        
  IF (IFP <= 0) IFP=0

  FC00(LAY) = FAC00(LAY) * CORR2(IFP)
  FC10(LAY) = FAC10(LAY) * CORR2(IFP)
  FC01(LAY) = FAC01(LAY) * CORR1(IFP)
  FC11(LAY) = FAC11(LAY) * CORR1(IFP)
  IND0(LAY) = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(2) + 1
  IND1(LAY) = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(2) + 1
ENDDO

!-- JJM_000517
DO IG = 1, NG2
  DO LAY = LAYTROP+1, KLEV
!-- JJM_000517
    TAU (NGS1+IG,LAY) = COLH2O(LAY) *&
     &(FC00(LAY) * ABSB(IND0(LAY)  ,IG) +&
     & FC10(LAY) * ABSB(IND0(LAY)+1,IG) +&
     & FC01(LAY) * ABSB(IND1(LAY)  ,IG) +&
     & FC11(LAY) * ABSB(IND1(LAY)+1,IG)&
     &+ FORFAC(LAY) * FORREF(IG) ) &
     &+ TAUAERL(LAY,2)
    PFRAC(NGS1+IG,LAY) = FRACREFB(IG)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL2
