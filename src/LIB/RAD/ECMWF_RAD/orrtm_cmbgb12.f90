!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB12
!***************************************************************************

!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG12

USE OYOERRTO12, ONLY : KAO     ,SELFREFO   ,FRACREFAO
USE OYOERRTA12, ONLY : KA      ,SELFREF    ,FRACREFA
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JP, JT

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK


DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(12)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(11)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+176)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(12)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(11)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+176)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(12)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(11)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(12)

    FREFA(NGS(11)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(12)


    FREFADF(NGS(11)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB12
