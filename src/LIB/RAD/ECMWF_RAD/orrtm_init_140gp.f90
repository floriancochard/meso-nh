!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:37
!-----------------------------------------------------------------
!***************************************************************************
SUBROUTINE ORRTM_INIT_140GP
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

! Parameters
#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT
USE OYOERRTWN , ONLY : NG       ,NSPA     ,NSPB
USE OYOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT
! Output
USE OYOERRTBG2, ONLY : CORR1    ,CORR2
USE OYOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
!
!USE MODI_ORRTM_CMBGB1
!USE MODI_ORRTM_CMBGB2
!USE MODI_ORRTM_CMBGB3
!USE MODI_ORRTM_CMBGB4
!USE MODI_ORRTM_CMBGB5
!USE MODI_ORRTM_CMBGB6
!USE MODI_ORRTM_CMBGB7
!USE MODI_ORRTM_CMBGB8
!USE MODI_ORRTM_CMBGB9
!USE MODI_ORRTM_CMBGB10
!USE MODI_ORRTM_CMBGB11
!USE MODI_ORRTM_CMBGB12
!USE MODI_ORRTM_CMBGB13
!USE MODI_ORRTM_CMBGB14
!USE MODI_ORRTM_CMBGB15
!USE MODI_ORRTM_CMBGB16
! Local

IMPLICIT NONE
REAL_B :: WTSM(JPG)

!     LOCAL INTEGER SCALARS
INTEGER_M :: I, IBND, IG, IGC, IGCSM, IND, IPR, IPRSM, IPT

!     LOCAL REAL SCALARS
REAL_B :: FP, RTFP, WTSUM


!  Calculate lookup tables for functions needed in routine TAUMOL (TAUGB2)
CORR1(0) = _ONE_
CORR1(200) = _ONE_
CORR2(0) = _ONE_
CORR2(200) = _ONE_
DO I = 1,199
  FP = 0.005_JPRB*REAL(I)
  RTFP = SQRT(FP)
  CORR1(I) = RTFP/FP
  CORR2(I) = (_ONE_-RTFP)/(_ONE_-FP)
ENDDO

!  Perform g-point reduction from 16 per band (256 total points) to
!  a band dependant number (140 total points) for all absorption
!  coefficient input data and Planck fraction input data.
!  Compute relative weighting for new g-point combinations.

IGCSM = 0
DO IBND = 1,JPBAND
  IPRSM = 0
  IF (NGC(IBND) < 16) THEN
    DO IGC = 1,NGC(IBND)
      IGCSM = IGCSM + 1
      WTSUM = _ZERO_
      DO IPR = 1, NGN(IGCSM)
        IPRSM = IPRSM + 1
        WTSUM = WTSUM + WT(IPRSM)
      ENDDO
      WTSM(IGC) = WTSUM
    ENDDO
    DO IG = 1,NG(IBND)
      IND = (IBND-1)*16 + IG
      RWGT(IND) = WT(IG)/WTSM(NGM(IND))
    ENDDO
  ELSE
    DO IG = 1,NG(IBND)
      IGCSM = IGCSM + 1
      IND = (IBND-1)*16 + IG
      RWGT(IND) = _ONE_
    ENDDO
  ENDIF
ENDDO

!  Initialize arrays for combined Planck fraction data.

DO IPT = 1,13
  DO IPR = 1, JPGPT
    FREFA(IPR,IPT) = _ZERO_
    FREFADF(IPR,IPT) = _ZERO_
  ENDDO
ENDDO
DO IPT = 1,6
  DO IPR = 1, JPGPT
    FREFB(IPR,IPT) = _ZERO_
    FREFBDF(IPR,IPT) = _ZERO_
  ENDDO
ENDDO

!  Reduce g-points for relevant data in each LW spectral band.

CALL ORRTM_CMBGB1
CALL ORRTM_CMBGB2
CALL ORRTM_CMBGB3
CALL ORRTM_CMBGB4
CALL ORRTM_CMBGB5
CALL ORRTM_CMBGB6
CALL ORRTM_CMBGB7
CALL ORRTM_CMBGB8
CALL ORRTM_CMBGB9
CALL ORRTM_CMBGB10
CALL ORRTM_CMBGB11
CALL ORRTM_CMBGB12
CALL ORRTM_CMBGB13
CALL ORRTM_CMBGB14
CALL ORRTM_CMBGB15
CALL ORRTM_CMBGB16

RETURN
END SUBROUTINE ORRTM_INIT_140GP
