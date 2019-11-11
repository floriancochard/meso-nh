!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB16
!***************************************************************************

!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG16

USE OYOERRTO16, ONLY : KAO     ,SELFREFO   ,FRACREFAO
USE OYOERRTA16, ONLY : KA      ,SELFREF    ,FRACREFA
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
      DO IGC = 1,NGC(16)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(15)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+240)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(16)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(15)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+240)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(16)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(15)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(16)

    FREFA(NGS(15)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO

DO JP = 1,8
  DO IGC = 1,NGC(16)


    FREFADF(NGS(15)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB16
