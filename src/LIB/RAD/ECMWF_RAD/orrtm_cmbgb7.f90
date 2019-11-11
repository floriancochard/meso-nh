!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB7
!***************************************************************************

!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG7

USE OYOERRTO7 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, ABSCO2O
USE OYOERRTA7 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
           &FRACREFB , ABSCO2
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
      DO IGC = 1,NGC(7)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(6)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+96)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(7)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(6)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+96)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(7)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(6)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+96)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(7)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(6)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(7)
  SUMF = _ZERO_
  SUMK = _ZERO_
  DO IPR = 1, NGN(NGS(6)+IGC)
    IPRSM = IPRSM + 1


    SUMF = SUMF + FRACREFBO(IPRSM)
    SUMK = SUMK + ABSCO2O(IPRSM)*RWGT(IPRSM+96)
  ENDDO


  FRACREFB(IGC) = SUMF
  ABSCO2(IGC) = SUMK
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(7)

    FREFA(NGS(6)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(7)


    FREFADF(NGS(6)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO IGC = 1,NGC(7)

  FREFB(NGS(6)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB7
