!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB5
!***************************************************************************

!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG5

USE OYOERRTO5 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, CCL4O
USE OYOERRTA5 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
           &FRACREFB , CCL4
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
      DO IGC = 1,NGC(5)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(4)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+64)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JN = 1,5
  DO JT = 1,5
    DO JP = 13,59
      IPRSM = 0
      DO IGC = 1,NGC(5)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(4)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KBO(JN,JT,JP,IPRSM)*RWGT(IPRSM+64)
        ENDDO

        KB(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(5)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(4)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+64)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(5)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(4)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

DO JP = 1,5
  IPRSM = 0
  DO IGC = 1,NGC(5)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(4)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFBO(IPRSM,JP)
    ENDDO

    FRACREFB(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(5)
  SUMK = _ZERO_
  DO IPR = 1, NGN(NGS(4)+IGC)
    IPRSM = IPRSM + 1

    SUMK = SUMK + CCL4O(IPRSM)*RWGT(IPRSM+64)
  ENDDO

  CCL4(IGC) = SUMK
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(5)

    FREFA(NGS(4)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(5)


    FREFADF(NGS(4)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,5
  DO IGC = 1,NGC(5)

    FREFB(NGS(4)+IGC,JP) = FRACREFB(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,4
  DO IGC = 1,NGC(5)


    FREFBDF(NGS(4)+IGC,JP) = FRACREFB(IGC,JP+1) -FRACREFB(IGC,JP)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB5
