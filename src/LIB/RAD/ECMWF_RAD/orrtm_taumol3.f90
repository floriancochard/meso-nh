!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL3 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLN2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG3   ,NGS2
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA3 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &FORREF   ,KA     ,KB     ,SELFREF , ABSN2OA ,&
            &ABSN2OB  ,ETAREF ,H2OREF ,N2OREF  , CO2REF  ,&
            &STRRAT

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
REAL_B :: FORFAC(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY)
INTEGER_M :: JT(JPLAY)
INTEGER_M :: JT1(JPLAY)

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY)
REAL_B :: COLCO2(JPLAY)
REAL_B :: COLN2O(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

INTEGER_M :: IJS(JPLAY)
REAL_B :: ZFS(JPLAY),SPECCOMB(JPLAY)
INTEGER_M :: IND0(JPLAY),IND1(JPLAY),INDS(JPLAY)
REAL_B :: N2OMULT(JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, JS, LAY, NS

!     LOCAL REAL SCALARS
REAL_B :: COLREF1, COLREF2, CURRN2O, FAC000, FAC001,&
          &FAC010, FAC011, FAC100, FAC101, FAC110, FAC111, &
          &FP, FS, RATIO, SPECMULT, SPECPARM, WCOMB1, &
          &WCOMB2

!      EQUIVALENCE (TAUAERL(1,3),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

DO LAY = 1, LAYTROP
  SPECCOMB(LAY) = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB(LAY)
  SPECPARM=MIN(ONEMINUS,SPECPARM)
  SPECMULT = 8._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT,_ONE_)
  IF (JS  ==  8) THEN
    IF (FS  >=  0.9_JPRB) THEN
      JS = 9
      FS = 10._JPRB * (FS - 0.9_JPRB)
    ELSE
      FS = FS/0.9_JPRB
    ENDIF
  ENDIF

  NS = JS + INT(FS + _HALF_)
  FP = FAC01(LAY) + FAC11(LAY)
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(3) + JS
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(3) + JS
  INDS(LAY) = INDSELF(LAY)
  COLREF1 = N2OREF(JP(LAY))
  COLREF2 = N2OREF(JP(LAY)+1)
  IF (NS  ==  10) THEN
    WCOMB1 = _ONE_/H2OREF(JP(LAY))
    WCOMB2 = _ONE_/H2OREF(JP(LAY)+1)
  ELSE
    WCOMB1 = (_ONE_-ETAREF(NS))/(STRRAT * CO2REF(JP(LAY)))
    WCOMB2 = (_ONE_-ETAREF(NS))/(STRRAT * CO2REF(JP(LAY)+1))
  ENDIF
  RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
  CURRN2O = SPECCOMB(LAY) * RATIO
  N2OMULT(LAY) = COLN2O(LAY) - CURRN2O

  ZFS(LAY)=FS
  IJS(LAY)=JS

ENDDO

!-- DS_000515
DO IG = 1, NG3
  DO LAY = 1, LAYTROP
!-- DS_000515

    FS=ZFS(LAY)
    JS=IJS(LAY)

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

    TAU (NGS2+IG,LAY) = SPECCOMB(LAY) *   &
!-- DS_000515
     & ( (1. - FS) *(FAC00(LAY) * ABSA(IND0(LAY)   ,IG) +   &
     &               FAC10(LAY) * ABSA(IND0(LAY)+10,IG) +   &
     &               FAC01(LAY) * ABSA(IND1(LAY)   ,IG) +   &
     &               FAC11(LAY) * ABSA(IND1(LAY)+10,IG))+   &
     &      FS     *(FAC00(LAY) * ABSA(IND0(LAY)+ 1,IG) +   &
     &               FAC10(LAY) * ABSA(IND0(LAY)+11,IG) +   &
     &               FAC01(LAY) * ABSA(IND1(LAY)+ 1,IG) +   &
     &               FAC11(LAY) * ABSA(IND1(LAY)+11,IG))) + &
!     &(FAC000 * ABSA(IND0(LAY)   ,IG) +&
!     & FAC100 * ABSA(IND0(LAY)+ 1,IG) +&
!     & FAC010 * ABSA(IND0(LAY)+10,IG) +&
!     & FAC110 * ABSA(IND0(LAY)+11,IG) +&
!     & FAC001 * ABSA(IND1(LAY),   IG) +&
!     & FAC101 * ABSA(IND1(LAY)+ 1,IG) +&
!     & FAC011 * ABSA(IND1(LAY)+10,IG) +&
!     & FAC111 * ABSA(IND1(LAY)+11,IG))+&
!-- DS_000515
     &COLH2O(LAY) * &
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG))&
     &+ FORFAC(LAY) * FORREF(IG) ) &
     &+ N2OMULT(LAY) * ABSN2OA(IG) &
     &+ TAUAERL(LAY,3)
    PFRAC(NGS2+IG,LAY) = FRACREFA(IG,JS) + FS *&
     &(FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
  ENDDO
ENDDO

DO LAY = LAYTROP+1, KLEV
  SPECCOMB(LAY) = COLH2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB(LAY)
  SPECPARM=MIN(ONEMINUS,SPECPARM)
  SPECMULT = 4._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT,_ONE_)
  NS = JS + INT(FS + _HALF_)
  FP = FAC01(LAY) + FAC11(LAY)
  IND0(LAY) = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(3) + JS
  IND1(LAY) = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(3) + JS
  COLREF1 = N2OREF(JP(LAY))
  COLREF2 = N2OREF(JP(LAY)+1)
  IF (NS  ==  5) THEN
    WCOMB1 = _ONE_/H2OREF(JP(LAY))
    WCOMB2 = _ONE_/H2OREF(JP(LAY)+1)
  ELSE
    WCOMB1 = (_ONE_-ETAREF(NS))/(STRRAT * CO2REF(JP(LAY)))
    WCOMB2 = (_ONE_-ETAREF(NS))/(STRRAT * CO2REF(JP(LAY)+1))
  ENDIF
  RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
  CURRN2O = SPECCOMB(LAY) * RATIO
  N2OMULT(LAY) = COLN2O(LAY) - CURRN2O

  ZFS(LAY)=FS
  IJS(LAY)=JS

ENDDO

DO LAY = LAYTROP+1, KLEV

  FS=ZFS(LAY)
  JS=IJS(LAY)
!---jjm
!  FAC000 = (_ONE_ - FS) * FAC00(LAY)
!  FAC010 = (_ONE_ - FS) * FAC10(LAY)
!  FAC100 = FS * FAC00(LAY)
!  FAC110 = FS * FAC10(LAY)
!  FAC001 = (_ONE_ - FS) * FAC01(LAY)
!  FAC011 = (_ONE_ - FS) * FAC11(LAY)
!  FAC101 = FS * FAC01(LAY)
!  FAC111 = FS * FAC11(LAY)
!---        

  DO IG = 1, NG3
    TAU (NGS2+IG,LAY) = SPECCOMB(LAY) *   &
!-- DS_000515
     & ( (1. - FS) *(FAC00(LAY) * ABSB(IND0(LAY)  ,IG) +   &
     &               FAC10(LAY) * ABSB(IND0(LAY)+5,IG) +   &
     &               FAC01(LAY) * ABSB(IND1(LAY)  ,IG) +   & 
     &               FAC11(LAY) * ABSB(IND1(LAY)+5,IG))+   &
     &      FS     *(FAC00(LAY) * ABSB(IND0(LAY)+1,IG) +   &
     &               FAC10(LAY) * ABSB(IND0(LAY)+6,IG) +   &
     &               FAC01(LAY) * ABSB(IND1(LAY)+1,IG) +   &
     &               FAC11(LAY) * ABSB(IND1(LAY)+6,IG)))   &
!     &(FAC000 * ABSB(IND0(LAY)  ,IG) +&
!     & FAC100 * ABSB(IND0(LAY)+1,IG) +&
!     & FAC010 * ABSB(IND0(LAY)+5,IG) +&
!     & FAC110 * ABSB(IND0(LAY)+6,IG) +&
!     & FAC001 * ABSB(IND1(LAY)  ,IG) +&
!     & FAC101 * ABSB(IND1(LAY)+1,IG) +&
!     & FAC011 * ABSB(IND1(LAY)+5,IG) +&
!     & FAC111 * ABSB(IND1(LAY)+6,IG))&
!-- DS_000515
     &+ COLH2O(LAY)*FORFAC(LAY)*FORREF(IG) &
     &+ N2OMULT(LAY) * ABSN2OB(IG)&
     &+ TAUAERL(LAY,3)
    PFRAC(NGS2+IG,LAY) = FRACREFB(IG,JS) + FS *&
     &(FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL3
