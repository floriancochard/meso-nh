!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:41
!-----------------------------------------------------------------
!******************************************************************************
SUBROUTINE ORRTM_TAUMOL11 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG11  ,NGS10
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA11, ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &KA      , KB     ,SELFREF

!  Input
!#include "yoeratm.h"

!      REAL TAUAER(JPLAY)

IMPLICIT NONE

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

!      EQUIVALENCE (TAUAERL(1,11),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.
  
DO LAY = 1, LAYTROP
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(11) + 1
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(11) + 1
  INDS(LAY) = INDSELF(LAY)
ENDDO

!-- DS_000515  
DO IG = 1, NG11
  DO LAY = 1, LAYTROP
!-- DS_000515  
    TAU (NGS10+IG,LAY) = COLH2O(LAY) *&
     &(FAC00(LAY) * ABSA(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSA(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSA(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSA(IND1(LAY)+1,IG) +&
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG))))&
     &+ TAUAERL(LAY,11)
    PFRAC(NGS10+IG,LAY) = FRACREFA(IG)
  ENDDO
ENDDO

DO LAY = LAYTROP+1, KLEV
  IND0(LAY) = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(11) + 1
  IND1(LAY) = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(11) + 1
ENDDO

!-- JJM_000517
DO IG = 1, NG11
  DO LAY = LAYTROP+1, KLEV
!-- JJM_000517
    TAU (NGS10+IG,LAY) = COLH2O(LAY) *&
     &(FAC00(LAY) * ABSB(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSB(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSB(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSB(IND1(LAY)+1,IG)) &
     &+ TAUAERL(LAY,11)
    PFRAC(NGS10+IG,LAY) = FRACREFB(IG)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL11
