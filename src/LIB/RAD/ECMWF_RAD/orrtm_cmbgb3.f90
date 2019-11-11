!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB3
!***************************************************************************

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG3

USE OYOERRTO3 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO  ,FORREFO    ,ABSN2OAO   ,ABSN2OBO
USE OYOERRTA3 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
           &FRACREFB   ,FORREF    ,ABSN2OA   ,ABSN2OB
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JP, JT

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK, SUMK1, SUMK2, SUMK3


DO JN = 1,10
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(3)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(2)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
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
      DO IGC = 1,NGC(3)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(2)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KBO(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
        ENDDO

        KB(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(3)
    SUMK = _ZERO_
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(2)+IGC)
      IPRSM = IPRSM + 1


      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+32)
      SUMF = SUMF + FRACREFAO(IPRSM,JT)
    ENDDO


    SELFREF(JT,IGC) = SUMK
    FRACREFA(IGC,JT) = SUMF
  ENDDO
ENDDO

DO JP = 1,5
  IPRSM = 0
  DO IGC = 1,NGC(3)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(2)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFBO(IPRSM,JP)
    ENDDO

    FRACREFB(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(3)
  SUMK1= _ZERO_
  SUMK2= _ZERO_
  SUMK3= _ZERO_
  DO IPR = 1, NGN(NGS(2)+IGC)
    IPRSM = IPRSM + 1



    SUMK1= SUMK1+ FORREFO(IPRSM)*RWGT(IPRSM+32)
    SUMK2= SUMK2+ ABSN2OAO(IPRSM)*RWGT(IPRSM+32)
    SUMK3= SUMK3+ ABSN2OBO(IPRSM)*RWGT(IPRSM+32)
  ENDDO



  FORREF(IGC) = SUMK1
  ABSN2OA(IGC) = SUMK2
  ABSN2OB(IGC) = SUMK3
ENDDO

DO JP = 1,10
  DO IGC = 1,NGC(3)

    FREFA(NGS(2)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,9
  DO IGC = 1,NGC(3)


    FREFADF(NGS(2)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,5
  DO IGC = 1,NGC(3)

    FREFB(NGS(2)+IGC,JP) = FRACREFB(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,4
  DO IGC = 1,NGC(3)


    FREFBDF(NGS(2)+IGC,JP) = FRACREFB(IGC,JP+1) -FRACREFB(IGC,JP)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB3
