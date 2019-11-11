!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:34
!-----------------------------------------------------------------
SUBROUTINE LEGTRI (PSIN,KCP,KDIM,PALP)

!**** *LEGTRI* - *LEGENDRE FUNCTIONS FOR A TRIANGULAR TRUNCATION.

!     J.F.GELEYN     E.C.M.W.F.     03/06/82.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE VALUES *PALP* FOR THE ARGUMENT
!     *PSIN* OF THE NORMALISED *LEGENDRE ASSOCIATED FUNCTIONS IN THE
!     ORDER ((JN1=JM1,KCP),JM1=1,KCP) FOR JN=JN1-1 AND JM=JM1-1 .

!**   INTERFACE.
!     ----------

!          *LEGTRI* IS CALLED FROM *RADMOD*.
!          THERE ARE THREE DUMMY ARGUMENTS: *PSIN* IS THE SINE OF
!     LATITUDE.
!                                           *KCP* IS ONE PLUS THE LIMIT
!     WAVE NUMBER.
!                                           *PALP* IS THE ARRAY OF THE
!     RESULTS.

!     METHOD.
!     -------

!          SIMPLE RECURENCE FORMULA.

!     EXTERNALS.
!     ----------

!          NONE.

!     REFERENCE.
!     ----------

!          NONE.


#include "tsmbkind.h"

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KCP
INTEGER_M :: KDIM

!     DUMMY REAL SCALARS
REAL_B :: PSIN

REAL_B :: PALP(KDIM)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IC, ICP, II, IM, IM2, JM1, JN

!     LOCAL REAL SCALARS
REAL_B :: Z2M, ZCOS, ZE1, ZE2, ZF1M, ZF2M, ZM, ZN, ZN2, ZRE1, ZSIN


!     ------------------------------------------------------------------

!*         1.     PRELIMINARY SETTING.
!                 ----------- --------


ZSIN=PSIN
ICP=KCP

!     ------------------------------------------------------------------

!*         2.     COMPUTATIONS.
!                 -------------


IC=ICP-1
ZCOS=SQRT(_ONE_-ZSIN**2)
II=2
PALP(1)=_ONE_
ZF1M=SQRT(3._JPRB)
PALP(2)=ZF1M*ZSIN
DO JM1=1,ICP
  IM=JM1-1
  ZM=IM
  Z2M=ZM+ZM
  ZRE1=SQRT(Z2M+3._JPRB)
  ZE1=_ONE_/ZRE1
  IF(IM == 0) GO TO 201
  ZF2M=ZF1M*ZCOS/SQRT(Z2M)
  ZF1M=ZF2M*ZRE1
  II=II+1
  PALP(II)=ZF2M
  IF(IM == IC) GO TO 203
  II=II+1
  PALP(II)=ZF1M*ZSIN
  IF(JM1 == IC) GO TO 203
  201 CONTINUE
  IM2=IM+2
  DO JN=IM2,IC
    ZN=JN
    ZN2=ZN**2
    ZE2=SQRT((4._JPRB*ZN2-_ONE_)/(ZN2-ZM**2))
    II=II+1
    PALP(II)=ZE2*(ZSIN*PALP(II-1)-ZE1*PALP(II-2))
    ZE1=_ONE_/ZE2
  ENDDO
  203 continue
ENDDO

!     ------------------------------------------------------------------

!*         3.     RETURN.
!                 -------


RETURN
END SUBROUTINE LEGTRI



