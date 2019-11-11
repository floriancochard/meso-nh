!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL13 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLN2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG13  ,NGS12
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA13, ONLY : ABSA   ,FRACREFA,KA    ,SELFREF,STRRAT

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

!- from PRECISE             
REAL_B :: ONEMINUS

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY)
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

!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, JS, LAY

!     LOCAL REAL SCALARS
REAL_B :: FAC000, FAC001, FAC010, FAC011, FAC100, FAC101,&
          &FAC110, FAC111, FS, SPECMULT, SPECPARM


!      EQUIVALENCE (TAUAERL(1,13),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
DO  LAY = 1, LAYTROP
  SPECCOMB(LAY) = COLH2O(LAY) + STRRAT*COLN2O(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB(LAY)
  SPECPARM=MIN(ONEMINUS,SPECPARM)
  SPECMULT = 8._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT,_ONE_)
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(13) + JS
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(13) + JS
  INDS(LAY) = INDSELF(LAY)
  
  ZFS(LAY)=FS
  IJS(LAY)=JS

ENDDO

!-- DS_000515
DO IG = 1, NG13
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
!----
    TAU (NGS12+IG,LAY) = SPECCOMB(LAY) *&
!-- DS_000515
     &( (1. - FS) *( FAC00(LAY) * ABSA(IND0(LAY)   ,IG) +   &
     &               FAC10(LAY) * ABSA(IND0(LAY)+ 9,IG) +   &
     &               FAC01(LAY) * ABSA(IND1(LAY)   ,IG) +   & 
     &               FAC11(LAY) * ABSA(IND1(LAY)+ 9,IG))+   &
     &    FS      *( FAC00(LAY) * ABSA(IND0(LAY)+ 1,IG) +   & 
     &               FAC10(LAY) * ABSA(IND0(LAY)+10,IG) +   &
     &               FAC01(LAY) * ABSA(IND1(LAY)+ 1,IG) +   &
     &               FAC11(LAY) * ABSA(IND1(LAY)+10,IG))) + &
!     &(FAC000 * ABSA(IND0(LAY)   ,IG) +&
!     & FAC100 * ABSA(IND0(LAY)+ 1,IG) +&
!     & FAC010 * ABSA(IND0(LAY)+ 9,IG) +&
!     & FAC110 * ABSA(IND0(LAY)+10,IG) +&
!     & FAC001 * ABSA(IND1(LAY)   ,IG) +&
!     & FAC101 * ABSA(IND1(LAY)+ 1,IG) +&
!     & FAC011 * ABSA(IND1(LAY)+ 9,IG) +&
!     & FAC111 * ABSA(IND1(LAY)+10,IG))+&
!-- DS_000515
     &COLH2O(LAY) * &
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG)))&
     &+ TAUAERL(LAY,13)
    PFRAC(NGS12+IG,LAY) = FRACREFA(IG,JS) + FS *&
     &(FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
  ENDDO
ENDDO

!-- JJM_000517
DO IG = 1, NG13
  DO LAY = LAYTROP+1, KLEV
!-- JJM_000517
    TAU (NGS12+IG,LAY) = TAUAERL(LAY,13)
    PFRAC(NGS12+IG,LAY) = _ZERO_
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL13
