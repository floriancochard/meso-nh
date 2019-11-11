!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL15 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLN2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)

! Modifications
!
!     D Salmond 1999-07-14 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG15  ,NGS14
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA15, ONLY : ABSA   ,FRACREFA,KA    ,SELFREF,STRRAT


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
REAL_B :: COLCO2(JPLAY)
REAL_B :: COLN2O(JPLAY)
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, IND0, IND1, INDS, JS, LAY

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
          &FAC110, FAC111, FS, SPECCOMB, SPECMULT, SPECPARM


!  Input
!#include "yoeratm.h"

!      REAL TAUAER(JPLAY)
!      EQUIVALENCE (TAUAERL(1,15),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
DO LAY = 1, LAYTROP
  SPECCOMB = COLN2O(LAY) + STRRAT*COLCO2(LAY)
  SPECPARM = COLN2O(LAY)/SPECCOMB
  SPECPARM = MIN(SPECPARM,ONEMINUS)
  SPECMULT = 8._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT,_ONE_)
!-----jjm        
  FAC000 = (_ONE_ - FS) * FAC00(LAY)
  FAC010 = (_ONE_ - FS) * FAC10(LAY)
  FAC100 = FS * FAC00(LAY)
  FAC110 = FS * FAC10(LAY)
  FAC001 = (_ONE_ - FS) * FAC01(LAY)
  FAC011 = (_ONE_ - FS) * FAC11(LAY)
  FAC101 = FS * FAC01(LAY)
  FAC111 = FS * FAC11(LAY)
!------         
  IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(15) + JS
  IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(15) + JS
  INDS = INDSELF(LAY)
!-- DS_990714  
!         DO IG = 1, NG15
  IG=1
  TAU (NGS14+IG,LAY) = SPECCOMB *&
!     &((1. - FS)*(FAC00(LAY) * ABSA(IND0,IG) +
!     &            FAC10(LAY) * ABSA(IND0+9,IG) +
!     &            FAC01(LAY) * ABSA(IND1,IG) + 
!     &            FAC11(LAY) * ABSA(IND1+9,IG)) +
!     &     FS *  (FAC01(LAY) * ABSA(IND1+1,IG) +
!     &            FAC10(LAY) * ABSA(IND0+10,IG) +
!     &            FAC00(LAY) * ABSA(IND0+1,IG) +
!     &            FAC11(LAY) * ABSA(IND1+10,IG))) +
   &(FAC000 * ABSA(IND0   ,IG) +&
   & FAC100 * ABSA(IND0+ 1,IG) +&
   & FAC010 * ABSA(IND0+ 9,IG) +&
   & FAC110 * ABSA(IND0+10,IG) +&
   & FAC001 * ABSA(IND1   ,IG) +&
   & FAC101 * ABSA(IND1+ 1,IG) +&
   & FAC011 * ABSA(IND1+ 9,IG) +&
   & FAC111 * ABSA(IND1+10,IG))+&
   &COLH2O(LAY) * &
   &SELFFAC(LAY) * (SELFREF(INDS,IG) + &
   &SELFFRAC(LAY) *&
   &(SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))&
   &+ TAUAERL(LAY,15)
  PFRAC(NGS14+IG,LAY) = FRACREFA(IG,JS) + FS *&
   &(FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
  IG=2
  TAU (NGS14+IG,LAY) = SPECCOMB *&
!     &((1. - FS)*(FAC00(LAY) * ABSA(IND0,IG) +
!     &            FAC10(LAY) * ABSA(IND0+9,IG) +
!     &            FAC01(LAY) * ABSA(IND1,IG) +
!     &            FAC11(LAY) * ABSA(IND1+9,IG)) +
!     &     FS *  (FAC01(LAY) * ABSA(IND1+1,IG) +
!     &            FAC10(LAY) * ABSA(IND0+10,IG) +
!     &            FAC00(LAY) * ABSA(IND0+1,IG) +
!     &            FAC11(LAY) * ABSA(IND1+10,IG))) +
   &(FAC000 * ABSA(IND0   ,IG) +&
   & FAC100 * ABSA(IND0+ 1,IG) +&
   & FAC010 * ABSA(IND0+ 9,IG) +&
   & FAC110 * ABSA(IND0+10,IG) +&
   & FAC001 * ABSA(IND1   ,IG) +&
   & FAC101 * ABSA(IND1+ 1,IG) +&
   & FAC011 * ABSA(IND1+ 9,IG) +&
   & FAC111 * ABSA(IND1+10,IG))+&
   &COLH2O(LAY) *&
   &SELFFAC(LAY) * (SELFREF(INDS,IG) +&
   &SELFFRAC(LAY) *&
   &(SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))&
   &+ TAUAERL(LAY,15)
  PFRAC(NGS14+IG,LAY) = FRACREFA(IG,JS) + FS *&
   &(FRACREFA(IG,JS+1) - FRACREFA(IG,JS))

!         END DO
!-- DS_990714  
ENDDO

DO LAY = LAYTROP+1, KLEV
!         DO IG = 1, NG15
!-- DS_990714  
  IG=1
  TAU (NGS14+IG,LAY) = TAUAERL(LAY,15)
  PFRAC(NGS14+IG,LAY) = _ZERO_
  IG=2
  TAU (NGS14+IG,LAY) = TAUAERL(LAY,15)
  PFRAC(NGS14+IG,LAY) = _ZERO_
!-- DS_990714  
!         END DO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL15
