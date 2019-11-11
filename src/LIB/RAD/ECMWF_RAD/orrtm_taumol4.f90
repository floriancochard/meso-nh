!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ORRTM_TAUMOL4 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLO3,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG4   ,NGS3
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA4 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &KA       ,KB     ,SELFREF,STRRAT1 , STRRAT2

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
REAL_B :: COLO3 (JPLAY)
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

!      EQUIVALENCE (TAUAERL(1,4),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
DO LAY = 1, LAYTROP
  SPECCOMB(LAY) = COLH2O(LAY) + STRRAT1*COLCO2(LAY)
  SPECPARM = COLH2O(LAY)/SPECCOMB(LAY)
  SPECPARM=MIN(ONEMINUS,SPECPARM)
  SPECMULT = 8._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT,_ONE_)
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(4) + JS
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(4) + JS
  INDS(LAY) = INDSELF(LAY)

  ZFS(LAY)=FS
  IJS(LAY)=JS

ENDDO

!-- DS_000515
DO IG = 1, NG4
  DO LAY = 1, LAYTROP
!-- DS_000515

    FS=ZFS(LAY)
    JS=IJS(LAY)
!--jjm         
!    FAC000 = (_ONE_ - FS) * FAC00(LAY)
!    FAC010 = (_ONE_ - FS) * FAC10(LAY)
!    FAC100 = FS * FAC00(LAY)
!    FAC110 = FS * FAC10(LAY)
!    FAC001 = (_ONE_ - FS) * FAC01(LAY)
!    FAC011 = (_ONE_ - FS) * FAC11(LAY)
!    FAC101 = FS * FAC01(LAY)
!    FAC111 = FS * FAC11(LAY)
!---

    TAU (NGS3+IG,LAY) = SPECCOMB(LAY) *   & 
!-- DS_000515
     &((1. - FS) *(FAC00(LAY) * ABSA(IND0(LAY)   ,IG) +   &
     &             FAC10(LAY) * ABSA(IND0(LAY)+ 9,IG) +   &
     &             FAC01(LAY) * ABSA(IND1(LAY)   ,IG) +   & 
     &             FAC11(LAY) * ABSA(IND1(LAY)+ 9,IG))+   &
     &    FS     *(FAC00(LAY) * ABSA(IND0(LAY)+ 1,IG) +   &
     &             FAC10(LAY) * ABSA(IND0(LAY)+10,IG) +   &
     &             FAC01(LAY) * ABSA(IND1(LAY)+ 1,IG) +   &
     &             FAC11(LAY) * ABSA(IND1(LAY)+10,IG))) + &
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
     &+ TAUAERL(LAY,4)
    PFRAC(NGS3+IG,LAY) = FRACREFA(IG,JS) + FS *&
     &(FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
  ENDDO
ENDDO

DO LAY = LAYTROP+1, KLEV
  SPECCOMB(LAY) = COLO3(LAY) + STRRAT2*COLCO2(LAY)
  SPECPARM = COLO3(LAY)/SPECCOMB(LAY)
  SPECPARM=MIN(ONEMINUS,SPECPARM)
  SPECMULT = 4._JPRB*(SPECPARM)
  JS = 1 + INT(SPECMULT)
  FS = MOD(SPECMULT,_ONE_)
  IF (JS  >  1) THEN
    JS = JS + 1
  ELSEIF (FS  >=  0.0024_JPRB) THEN
    JS = 2
    FS = (FS - 0.0024_JPRB)/0.9976_JPRB
  ELSE
    JS = 1
    FS = FS/0.0024_JPRB
  ENDIF
  IND0(LAY) = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(4) + JS
  IND1(LAY) = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(4) + JS
  ZFS(LAY)=FS
  IJS(LAY)=JS
ENDDO

DO LAY = LAYTROP+1, KLEV
  FS=ZFS(LAY)
  JS=IJS(LAY)
!--- jjm
!  FAC000 = (_ONE_ - FS) * FAC00(LAY)
!  FAC010 = (_ONE_ - FS) * FAC10(LAY)
!  FAC100 = FS * FAC00(LAY)
!  FAC110 = FS * FAC10(LAY)
!  FAC001 = (_ONE_ - FS) * FAC01(LAY)
!  FAC011 = (_ONE_ - FS) * FAC11(LAY)
!  FAC101 = FS * FAC01(LAY)
!  FAC111 = FS * FAC11(LAY)
!------                   
  DO IG = 1, NG4
    TAU (NGS3+IG,LAY) = SPECCOMB(LAY) *   &
!-- DS_000515
     &( (1. - FS) *(FAC00(LAY) * ABSB(IND0(LAY)  ,IG) +   &
     &              FAC10(LAY) * ABSB(IND0(LAY)+6,IG) +   &
     &              FAC01(LAY) * ABSB(IND1(LAY)  ,IG) +   & 
     &              FAC11(LAY) * ABSB(IND1(LAY)+6,IG))+   &
     &     FS     *(FAC00(LAY) * ABSB(IND0(LAY)+1,IG) +   &
     &              FAC10(LAY) * ABSB(IND0(LAY)+7,IG) +   &
     &              FAC01(LAY) * ABSB(IND1(LAY)+1,IG) +   &
     &              FAC11(LAY) * ABSB(IND1(LAY)+7,IG)))   &
!     &(FAC000 * ABSB(IND0(LAY)   ,IG) +&
!     & FAC100 * ABSB(IND0(LAY)+ 1,IG) +&
!     & FAC010 * ABSB(IND0(LAY)+ 6,IG) +&
!     & FAC110 * ABSB(IND0(LAY)+ 7,IG) +&
!     & FAC001 * ABSB(IND1(LAY)   ,IG) +&
!     & FAC101 * ABSB(IND1(LAY)+ 1,IG) +&
!     & FAC011 * ABSB(IND1(LAY)+ 6,IG) +&
!     & FAC111 * ABSB(IND1(LAY)+ 7,IG))&
!-- DS_000515
     &+ TAUAERL(LAY,4)
    PFRAC(NGS3+IG,LAY) = FRACREFB(IG,JS) + FS *&
     &(FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL4
