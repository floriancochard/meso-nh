!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB10
!***************************************************************************

!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG10

USE OYOERRTO10, ONLY : KAO     ,KBO      ,FRACREFAO   ,FRACREFBO
USE OYOERRTA10, ONLY : KA      ,KB       ,FRACREFA    ,FRACREFB
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
    DO IGC = 1,NGC(10)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(9)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+144)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(10)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(9)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+144)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(10)
  SUMF1= _ZERO_
  SUMF2= _ZERO_
  DO IPR = 1, NGN(NGS(9)+IGC)
    IPRSM = IPRSM + 1


    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
  ENDDO


  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
ENDDO

DO IGC = 1,NGC(10)


  FREFA(NGS(9)+IGC,1) = FRACREFA(IGC)
  FREFB(NGS(9)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB10
