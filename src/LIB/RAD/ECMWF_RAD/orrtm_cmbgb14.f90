!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB14
!***************************************************************************

!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG14

USE OYOERRTO14, ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,FRACREFBO
USE OYOERRTA14, ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,FRACREFB
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JP, JT

!     LOCAL REAL SCALARS
REAL_B :: SUMF1, SUMF2, SUMK


DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(14)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(13)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+208)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(14)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(13)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+208)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(14)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(13)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+208)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(14)
  SUMF1= _ZERO_
  SUMF2= _ZERO_
  DO IPR = 1, NGN(NGS(13)+IGC)
    IPRSM = IPRSM + 1


    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
  ENDDO


  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
ENDDO

DO IGC = 1,NGC(14)


  FREFA(NGS(13)+IGC,1) = FRACREFA(IGC)
  FREFB(NGS(13)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB14
