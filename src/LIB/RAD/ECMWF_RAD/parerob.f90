!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PAREROB


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!*
!     PAREROB   D. VASILJEVIC  91/03/05

!     DIMENSIONS FOR ARRAYS FOR OBS. ERROR CORRELATION


!     VARIABLE     / TYPE /               MEANING
!     --------       ----                 ------
!     JPMXXENT        I     MAX. NO. OF X ENTRIES IN THE CORREL. TABLE
!     JPMXYENT        I     MAX. NO. OF Y ENTRIES IN THE CORREL. TABLE
!     JPRZPAR         I
!     JPCZR           I
!     JPCXR           I
!     JPCZRCOR        I
!*

INTEGER_M, PARAMETER :: JPMXXENT=9
INTEGER_M, PARAMETER :: JPMXYENT=9
INTEGER_M, PARAMETER :: JPRZPAR=10
INTEGER_M, PARAMETER :: JPCZR=8
INTEGER_M, PARAMETER :: JPCXR=8
INTEGER_M, PARAMETER :: JPCZRCOR=2

!     ------------------------------------------------------------------

END MODULE PAREROB
