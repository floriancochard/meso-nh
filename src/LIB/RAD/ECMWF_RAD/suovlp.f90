!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:46
!-----------------------------------------------------------------
SUBROUTINE SUOVLP ( KLEV, PAZ )
  
!***** *SUOVLP*   -INITIALIZE PROFILE OF ALPHA1

!**   INTERFACE.
!     ----------
!        CALL *SUOVLP* FROM *SUECRAD*
!              ------        -------
!
!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 01-02-16
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOERAD   , ONLY : RAOVLP, RBOVLP
USE YOEOVLP  , ONLY : RA1OVLP

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAZ(KLEV)

INTEGER_M :: KLEV

!     ------------------------------------------------------------------

!     LOCAL INTEGER PARAMETERS

INTEGER_M :: JK

DO JK=1,KLEV
  RA1OVLP(JK)=RAOVLP*PAZ(JK)+RBOVLP
END DO  

!     ------------------------------------------------------------------
RETURN
END SUBROUTINE SUOVLP
