!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB8
!***************************************************************************

!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG8

USE OYOERRTO8 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, ABSCO2AO,ABSCO2BO,ABSN2OAO   ,ABSN2OBO   ,&
           &CFC12O   , CFC22ADJO
USE OYOERRTA8 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
           &FRACREFB , ABSCO2A ,ABSCO2B ,ABSN2OA    ,ABSN2OB    ,&
           &CFC12    , CFC22ADJ
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JP, JT

!     LOCAL REAL SCALARS
REAL_B :: SUMF1, SUMF2, SUMK, SUMK1, SUMK2, SUMK3, SUMK4, SUMK5, SUMK6


DO JT = 1,5
  DO JP = 1,7
    IPRSM = 0
    DO IGC = 1,NGC(8)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(7)+IGC)
        IPRSM = IPRSM + 1
        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+112)
      ENDDO
      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 7,59
    IPRSM = 0
    DO IGC = 1,NGC(8)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(7)+IGC)
        IPRSM = IPRSM + 1
        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+112)
      ENDDO
      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(8)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(7)+IGC)
      IPRSM = IPRSM + 1
      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+112)
    ENDDO
    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(8)
  SUMF1= _ZERO_
  SUMF2= _ZERO_
  SUMK1= _ZERO_
  SUMK2= _ZERO_
  SUMK3= _ZERO_
  SUMK4= _ZERO_
  SUMK5= _ZERO_
  SUMK6= _ZERO_
  DO IPR = 1, NGN(NGS(7)+IGC)
    IPRSM = IPRSM + 1
    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
    SUMK1= SUMK1+ ABSCO2AO(IPRSM)*RWGT(IPRSM+112)
    SUMK2= SUMK2+ ABSCO2BO(IPRSM)*RWGT(IPRSM+112)
    SUMK3= SUMK3+ ABSN2OAO(IPRSM)*RWGT(IPRSM+112)
    SUMK4= SUMK4+ ABSN2OBO(IPRSM)*RWGT(IPRSM+112)
    SUMK5= SUMK5+ CFC12O(IPRSM)*RWGT(IPRSM+112)
    SUMK6= SUMK6+ CFC22ADJO(IPRSM)*RWGT(IPRSM+112)
  ENDDO
  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
  ABSCO2A(IGC) = SUMK1
  ABSCO2B(IGC) = SUMK2
  ABSN2OA(IGC) = SUMK3
  ABSN2OB(IGC) = SUMK4
  CFC12(IGC) = SUMK5
  CFC22ADJ(IGC) = SUMK6
ENDDO

DO IGC = 1,NGC(8)
  FREFA(NGS(7)+IGC,1) = FRACREFA(IGC)
  FREFB(NGS(7)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB8
