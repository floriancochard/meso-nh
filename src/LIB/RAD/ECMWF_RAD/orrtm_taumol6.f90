!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL6 (KLEV,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,CO2MULT,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG6   ,NGS5
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA6 , ONLY : ABSA   ,ABSCO2 ,CFC11ADJ , CFC12  ,&
            &FRACREFA, KA     ,SELFREF

!  Input
!#include "yoeratm.h"

!      REAL TAUAER(JPLAY)

IMPLICIT NONE

REAL_B :: WX(JPXSEC,JPLAY)              ! Amount of trace gases
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
REAL_B :: CO2MULT(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

INTEGER_M :: IND0(JPLAY),IND1(JPLAY),INDS(JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, LAY

!      EQUIVALENCE (TAUAERL(1,6),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature. The water vapor self-continuum is interpolated
!     (in temperature) separately.  

DO LAY = 1, LAYTROP
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(6) + 1
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(6) + 1
  INDS(LAY) = INDSELF(LAY)
ENDDO

!-- DS_000515  
DO IG = 1, NG6
  DO LAY = 1, LAYTROP
!-- DS_000515  
    TAU (NGS5+IG,LAY) = COLH2O(LAY) *&
     &(FAC00(LAY) * ABSA(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSA(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSA(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSA(IND1(LAY)+1,IG) +&
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY)*&
     &(SELFREF(INDS(LAY)+1,IG)-SELFREF(INDS(LAY),IG))))&
     &+ WX(2,LAY) * CFC11ADJ(IG)&
     &+ WX(3,LAY) * CFC12(IG)&
     &+ CO2MULT(LAY) * ABSCO2(IG)&
     &+ TAUAERL(LAY,6)
    PFRAC(NGS5+IG,LAY) = FRACREFA(IG)
  ENDDO
ENDDO

!     Nothing important goes on above LAYTROP in this band.
!-- JJM_000517
DO IG = 1, NG6
  DO LAY = LAYTROP+1, KLEV
!-- JJM_000517
    TAU (NGS5+IG,LAY) = _ZERO_ &
     &+ WX(2,LAY) * CFC11ADJ(IG)&
     &+ WX(3,LAY) * CFC12(IG)&
     &+ TAUAERL(LAY,6)
    PFRAC(NGS5+IG,LAY) = FRACREFA(IG)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL6
