!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB2
!***************************************************************************

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG2

USE OYOERRTO2 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO  ,FORREFO
USE OYOERRTA2 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
           &FRACREFB   ,FORREF     ,REFPARAM 
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JP, JT

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK


DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(2)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(1)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+16)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(2)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(1)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+16)
      ENDDO
!               KBC(JT,JP,IGC) = SUMK
      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(2)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(1)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+16)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,13
  IPRSM = 0
  DO IGC = 1,NGC(2)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(1)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(2)
  SUMK = _ZERO_
  SUMF = _ZERO_
  DO IPR = 1, NGN(NGS(1)+IGC)
    IPRSM = IPRSM + 1


    SUMK = SUMK + FORREFO(IPRSM)*RWGT(IPRSM+16)
    SUMF = SUMF + FRACREFBO(IPRSM)
  ENDDO


  FORREF(IGC) = SUMK
  FRACREFB(IGC) = SUMF
ENDDO

DO JP = 1,13
  DO IGC = 1,NGC(2)

    FREFA(NGS(1)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 2,13
  DO IGC = 1,NGC(2)


    FREFADF(NGS(1)+IGC,JP) = FRACREFA(IGC,JP-1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO IGC = 1,NGC(2)

  FREFB(NGS(1)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB2
