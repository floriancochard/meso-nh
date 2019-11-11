!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_CMBGB1
!***************************************************************************

!  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
!  in arrays FRACREFA and FRACREFB are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTM.

!  BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,NG1

USE OYOERRTO1 , ONLY : KAO, KBO, SELFREFO, FORREFO, FRACREFAO,FRACREFBO
USE OYOERRTA1 , ONLY : KA , KB , SELFREF , FORREF , FRACREFA ,FRACREFB
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
    DO IGC = 1,NGC(1)
      SUMK = _ZERO_
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(1)
      SUMK = _ZERO_
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(1)
    SUMK = _ZERO_
    DO IPR = 1, NGN(IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(1)
  SUMK = _ZERO_
  SUMF1 = _ZERO_
  SUMF2 = _ZERO_
  DO IPR = 1, NGN(IGC)
    IPRSM = IPRSM + 1



    SUMK = SUMK + FORREFO(IPRSM)*RWGT(IPRSM)
    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
  ENDDO



  FORREF(IGC) = SUMK
  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
ENDDO

DO IGC = 1,NGC(1)


  FREFA(IGC,1) = FRACREFA(IGC)
  FREFB(IGC,1) = FRACREFB(IGC)
ENDDO

RETURN
END SUBROUTINE ORRTM_CMBGB1
