!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!******************************************************************************
SUBROUTINE ORRTM_TAUMOL14 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLCO2,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

! Modifications
!
!     D Salmond 1999-07-14 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG14  ,NGS13
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA14, ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
              &KA    , KB     ,SELFREF

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
REAL_B :: COLCO2(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, LAY


!  Input
!#include "yoeratm.h"

!      REAL TAUAER(JPLAY)
!      EQUIVALENCE (TAUAERL(1,14),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.  

DO LAY = 1, LAYTROP
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(14) + 1
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(14) + 1
  INDS = INDSELF(LAY)
!-- DS_990714  
!         DO IG = 1, NG14
  IG=1
  TAU (NGS13+IG,LAY) = COLCO2(LAY) *&
   &(FAC00(LAY) * ABSA(IND0  ,IG) +&
   & FAC10(LAY) * ABSA(IND0+1,IG) +&
   & FAC01(LAY) * ABSA(IND1  ,IG) +&
   & FAC11(LAY) * ABSA(IND1+1,IG) +&
   &SELFFAC(LAY) * (SELFREF(INDS,IG) + &
   &SELFFRAC(LAY) *&
   &(SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))&
   &+ TAUAERL(LAY,14)
  PFRAC(NGS13+IG,LAY) = FRACREFA(IG)
  IG=2
  TAU (NGS13+IG,LAY) = COLCO2(LAY) *&
   &(FAC00(LAY) * ABSA(IND0  ,IG) +&
   & FAC10(LAY) * ABSA(IND0+1,IG) +&
   & FAC01(LAY) * ABSA(IND1  ,IG) +&
   & FAC11(LAY) * ABSA(IND1+1,IG) +&
   &SELFFAC(LAY) * (SELFREF(INDS,IG) +&
   &SELFFRAC(LAY) *&
   &(SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))&
   &+ TAUAERL(LAY,14)
  PFRAC(NGS13+IG,LAY) = FRACREFA(IG)
!         END DO
!-- DS_990714  
ENDDO

DO LAY = LAYTROP+1, KLEV
  IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(14) + 1
  IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(14) + 1
!-- DS_990714  
!         DO IG = 1, NG14
  IG=1
  TAU (NGS13+IG,LAY) = COLCO2(LAY) *&
   &(FAC00(LAY) * ABSB(IND0  ,IG) +&
   & FAC10(LAY) * ABSB(IND0+1,IG) +&
   & FAC01(LAY) * ABSB(IND1  ,IG) +&
   & FAC11(LAY) * ABSB(IND1+1,IG)) &
   &+ TAUAERL(LAY,14)
  PFRAC(NGS13+IG,LAY) = FRACREFB(IG)
  IG=2
  TAU (NGS13+IG,LAY) = COLCO2(LAY) *&
   &(FAC00(LAY) * ABSB(IND0  ,IG) +&
   & FAC10(LAY) * ABSB(IND0+1,IG) +&
   & FAC01(LAY) * ABSB(IND1  ,IG) +&
   & FAC11(LAY) * ABSB(IND1+1,IG)) &
   &+ TAUAERL(LAY,14)
  PFRAC(NGS13+IG,LAY) = FRACREFB(IG)
!         END DO
!-- DS_990714  
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL14
