!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB9
!***************************************************************************

!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG9

USE OYOERRTO9 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, ABSN2OO
USE OYOERRTA9 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA  ,&
           &FRACREFB , ABSN2O
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JND, JNDC, JP, JT

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK


DO JN = 1,11
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(9)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(8)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+128)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(9)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(8)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+128)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(9)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(8)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+128)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JN = 1,3
  IPRSM = 0
  DO IGC = 1,NGC(9)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(8)+IGC)
      IPRSM = IPRSM + 1
      JND = (JN-1)*16

      SUMK = SUMK + ABSN2OO(JND+IPRSM)*RWGT(IPRSM+128)
    ENDDO
    JNDC = (JN-1)*NGC(9)

    ABSN2O(JNDC+IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(9)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(8)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(9)
  SUMF = _ZERO_
  DO IPR = 1, NGN(NGS(8)+IGC)
    IPRSM = IPRSM + 1

    SUMF = SUMF + FRACREFBO(IPRSM)
  ENDDO

  FRACREFB(IGC) = SUMF
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(9)

    FREFA(NGS(8)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(9)


    FREFADF(NGS(8)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO IGC = 1,NGC(9)

  FREFB(NGS(8)+IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB9
