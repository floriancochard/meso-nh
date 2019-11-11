!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB11
!***************************************************************************

!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG11

USE OYOERRTO11, ONLY : KAO     ,KBO     ,SELFREFO    ,FRACREFAO ,FRACREFBO
USE OYOERRTA11, ONLY : KA      ,KB      ,SELFREF     ,FRACREFA  ,FRACREFB
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
    DO IGC = 1,NGC(11)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(10)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+160)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(11)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(10)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+160)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(11)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(10)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+160)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(11)
  SUMF1= _ZERO_
  SUMF2= _ZERO_
  DO IPR = 1, NGN(NGS(10)+IGC)
    IPRSM = IPRSM + 1


    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
  ENDDO


  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
ENDDO

DO IGC = 1,NGC(11)


  FREFA(NGS(10)+IGC,1) = FRACREFA(IGC)
  FREFB(NGS(10)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB11
