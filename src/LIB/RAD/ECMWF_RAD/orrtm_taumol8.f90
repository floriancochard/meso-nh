!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!*******************************************************************************
SUBROUTINE ORRTM_TAUMOL8 (KLEV,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,COLO3,COLN2O,CO2MULT,LAYSWTCH,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG8   ,NGS7
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA8 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &KA      , KB     ,SELFREF,ABSCO2A , ABSCO2B ,&
            &ABSN2OA , ABSN2OB,CFC12  ,CFC22ADJ, H2OREF  ,&
            &N2OREF  , O3REF

!  Input
!#include "yoeratm.h"


IMPLICIT NONE

REAL_B :: WX(JPXSEC,JPLAY)             ! Amount of trace gases
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

!- from INTIND
INTEGER_M :: JP(JPLAY)
INTEGER_M :: JT(JPLAY)
INTEGER_M :: JT1(JPLAY)

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY)
REAL_B :: COLO3 (JPLAY)
REAL_B :: COLN2O(JPLAY)
REAL_B :: CO2MULT(JPLAY)
INTEGER_M :: LAYSWTCH

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

INTEGER_M :: IND0(JPLAY),IND1(JPLAY),INDS(JPLAY)

!      REAL TAUAER(JPLAY)
REAL_B :: N2OMULT(JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, LAY

!     LOCAL REAL SCALARS
REAL_B :: COLREF1, COLREF2, CURRN2O, FP, RATIO, WCOMB1, WCOMB2

!      EQUIVALENCE (TAUAERL(1,8),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  

DO LAY = 1, LAYSWTCH
  FP = FAC01(LAY) + FAC11(LAY)
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(8) + 1
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(8) + 1
  INDS(LAY) = INDSELF(LAY)
  COLREF1 = N2OREF(JP(LAY))
  COLREF2 = N2OREF(JP(LAY)+1)
  WCOMB1 = _ONE_/H2OREF(JP(LAY))
  WCOMB2 = _ONE_/H2OREF(JP(LAY)+1)
  RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
  CURRN2O = COLH2O(LAY) * RATIO
  N2OMULT(LAY) = COLN2O(LAY) - CURRN2O
ENDDO

!-- DS_000515
DO IG = 1, NG8
  DO LAY = 1, LAYSWTCH
!-- DS_000515
    TAU (NGS7+IG,LAY) = COLH2O(LAY) *&
     &(FAC00(LAY) * ABSA(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSA(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSA(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSA(IND1(LAY)+1,IG) +&
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG))))&
     &+ WX(3,LAY) * CFC12(IG)&
     &+ WX(4,LAY) * CFC22ADJ(IG)&
     &+ CO2MULT(LAY) * ABSCO2A(IG)&
     &+ N2OMULT(LAY) * ABSN2OA(IG)&
     &+ TAUAERL(LAY,8)
    PFRAC(NGS7+IG,LAY) = FRACREFA(IG)
  ENDDO
ENDDO

DO LAY = LAYSWTCH+1, KLEV
  FP = FAC01(LAY) + FAC11(LAY)
  IND0(LAY) = ((JP(LAY)-7)*5+(JT(LAY)-1))*NSPB(8) + 1
  IND1(LAY) = ((JP(LAY)-6)*5+(JT1(LAY)-1))*NSPB(8) + 1
  COLREF1 = N2OREF(JP(LAY))
  COLREF2 = N2OREF(JP(LAY)+1)
  WCOMB1 = _ONE_/O3REF(JP(LAY))
  WCOMB2 = _ONE_/O3REF(JP(LAY)+1)
  RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
  CURRN2O = COLO3(LAY) * RATIO
  N2OMULT(LAY) = COLN2O(LAY) - CURRN2O
ENDDO

!-- JJM_000517
DO IG = 1, NG8
  DO LAY = LAYSWTCH+1, KLEV
!-- JJM_000517
    TAU (NGS7+IG,LAY) = COLO3(LAY) *&
     &(FAC00(LAY) * ABSB(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSB(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSB(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSB(IND1(LAY)+1,IG)) &
     &+ WX(3,LAY) * CFC12(IG)&
     &+ WX(4,LAY) * CFC22ADJ(IG)&
     &+ CO2MULT(LAY) * ABSCO2B(IG)&
     &+ N2OMULT(LAY) * ABSN2OB(IG)&
     &+ TAUAERL(LAY,8)
    PFRAC(NGS7+IG,LAY) = FRACREFB(IG)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL8
