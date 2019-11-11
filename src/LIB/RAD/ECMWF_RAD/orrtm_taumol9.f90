!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL9 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLN2O,COLCH4,LAYTROP,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG9   ,NGS8
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA9 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &KA      , KB     ,SELFREF,ABSN2O  , CH4REF  ,&
            &ETAREF  , H2OREF ,N2OREF ,STRRAT

!  Input
!#include "yoeratm.h"


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

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY)
REAL_B :: COLN2O(JPLAY)
REAL_B :: COLCH4(JPLAY)
INTEGER_M :: LAYTROP
INTEGER_M :: LAYSWTCH
INTEGER_M :: LAYLOW

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

INTEGER_M :: JFRAC(JPLAY)
REAL_B :: FFRAC(JPLAY),ZFS(JPLAY),SPECCOMB(JPLAY)
INTEGER_M :: IND0(JPLAY),IND1(JPLAY),INDS(JPLAY),IIOFF(JPLAY)

!      REAL TAUAER(JPLAY)
REAL_B :: N2OMULT(JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IOFF, JS, LAY, NS

!     LOCAL REAL SCALARS
REAL_B :: COLREF1, COLREF2, CURRN2O, FAC000, FAC001,&
          &FAC010, FAC011, FAC100, FAC101, FAC110, FAC111, &
          &FP, FS, RATIO, SPECMULT, SPECPARM, WCOMB1, &
          &WCOMB2

!      EQUIVALENCE (TAUAERL(1,9),TAUAER)

IOFF = 0

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.
  
DO LAY = 1, LAYTROP
  SPECCOMB(LAY) = COLH2O(LAY) + STRRAT*COLCH4(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB(LAY)
  SPECPARM=MIN(ONEMINUS,SPECPARM)
  SPECMULT = 8._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  JFRAC(LAY) = JS
  FS = MOD(SPECMULT,_ONE_)
  FFRAC(LAY) = FS
  IF (JS  ==  8) THEN
    IF (FS.LE. 0.68_JPRB) THEN
      FS = FS/0.68_JPRB
    ELSEIF (FS  <=  0.92_JPRB) THEN
      JS = JS + 1
      FS = (FS-0.68_JPRB)/0.24_JPRB
    ELSE
      JS = JS + 2
      FS = (FS-0.92_JPRB)/0.08_JPRB
    ENDIF
  ELSEIF (JS  == 9) THEN
    JS = 10
    FS = _ONE_
    JFRAC(LAY) = 8
    FFRAC(LAY) = _ONE_
  ENDIF
  FP = FAC01(LAY) + FAC11(LAY)
  NS = JS + INT(FS + _HALF_)
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(9) + JS
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(9) + JS
  INDS(LAY) = INDSELF(LAY)
  IF (LAY  ==  LAYLOW) IOFF = NG9
  IF (LAY  ==  LAYSWTCH) IOFF = 2*NG9
  COLREF1 = N2OREF(JP(LAY))
  COLREF2 = N2OREF(JP(LAY)+1)
  IF (NS  ==  11) THEN
    WCOMB1 = _ONE_/H2OREF(JP(LAY))
    WCOMB2 = _ONE_/H2OREF(JP(LAY)+1)
  ELSE
    WCOMB1 = (_ONE_-ETAREF(NS))/(STRRAT * CH4REF(JP(LAY)))
    WCOMB2 = (_ONE_-ETAREF(NS))/(STRRAT * CH4REF(JP(LAY)+1))
  ENDIF
  RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
  CURRN2O = SPECCOMB(LAY) * RATIO
  N2OMULT(LAY) = COLN2O(LAY) - CURRN2O

  ZFS(LAY)=FS
  IIOFF(LAY)=IOFF

ENDDO

!-- DS_000515
DO IG = 1, NG9
  DO LAY = 1, LAYTROP
!-- DS_000515

    FS=ZFS(LAY)
    IOFF=IIOFF(LAY)
!---jjm
!    FAC000 = (_ONE_ - FS) * FAC00(LAY)
!    FAC010 = (_ONE_ - FS) * FAC10(LAY)
!    FAC100 = FS * FAC00(LAY)
!    FAC110 = FS * FAC10(LAY)
!    FAC001 = (_ONE_ - FS) * FAC01(LAY)
!    FAC011 = (_ONE_ - FS) * FAC11(LAY)
!    FAC101 = FS * FAC01(LAY)
!    FAC111 = FS * FAC11(LAY)
!------         

    TAU (NGS8+IG,LAY) = SPECCOMB(LAY) *&
!-- DS_000515
!     &(FAC000 * ABSA(IND0(LAY)   ,IG) +&
!     & FAC100 * ABSA(IND0(LAY)+ 1,IG) +&
!     & FAC010 * ABSA(IND0(LAY)+11,IG) +&
!     & FAC110 * ABSA(IND0(LAY)+12,IG) +&
!     & FAC001 * ABSA(IND1(LAY)   ,IG) +&
!     & FAC101 * ABSA(IND1(LAY)+ 1,IG) +&
!     & FAC011 * ABSA(IND1(LAY)+11,IG) +&
!     & FAC111 * ABSA(IND1(LAY)+12,IG))+&
     &( (1. - FS) *(FAC00(LAY) * ABSA(IND0(LAY)   ,IG) +   &
     &              FAC10(LAY) * ABSA(IND0(LAY)+11,IG) +   &
     &              FAC01(LAY) * ABSA(IND1(LAY)   ,IG) +   &
     &              FAC11(LAY) * ABSA(IND1(LAY)+11,IG))+   &
     &     FS     *(FAC00(LAY) * ABSA(IND0(LAY)+ 1,IG) +   &
     &              FAC10(LAY) * ABSA(IND0(LAY)+12,IG) +   &
     &              FAC01(LAY) * ABSA(IND1(LAY)+ 1,IG) +   &
     &              FAC11(LAY) * ABSA(IND1(LAY)+12,IG))) + &
!-- DS_000515
     &COLH2O(LAY) * &
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG)))&
     &+ N2OMULT(LAY) * ABSN2O(IG+IOFF)&
     &+ TAUAERL(LAY,9)
    PFRAC(NGS8+IG,LAY) = FRACREFA(IG,JFRAC(LAY)) + FFRAC(LAY) *&
     &(FRACREFA(IG,JFRAC(LAY)+1) - FRACREFA(IG,JFRAC(LAY)))
  ENDDO
ENDDO

DO LAY = LAYTROP+1, KLEV
  IND0(LAY) = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(9) + 1
  IND1(LAY) = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(9) + 1
ENDDO

!-- JJM_000517
DO IG = 1, NG9
  DO LAY = LAYTROP+1, KLEV
!-- JJM_000517
    TAU (NGS8+IG,LAY) = COLCH4(LAY) *&
     &(FAC00(LAY) * ABSB(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSB(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSB(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSB(IND1(LAY)+1,IG))&
     &+ TAUAERL(LAY,9)
    PFRAC(NGS8+IG,LAY) = FRACREFB(IG)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL9
