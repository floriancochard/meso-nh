!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
SUBROUTINE SWUVO3 &
  &( KIDIA,KFDIA,KLON,KNU,KABS &
  &, PU, PTR &
  & )
  
!**** *SWUVO3* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR OZONE
!     IN THE UV and VISIBLE SPECTRAL INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWUVO3* IS CALLED FROM *SW1S*.


!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING SUMS OF EXPONENTIALS

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        Ph.DUBUISSON/B.BONNEL L.O.A.

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 00-12-18
   
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE OYOESW    , ONLY : NEXPO3, REXPO3


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KABS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLON
INTEGER_M :: KNU

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PU(KLON,KABS)
REAL_B :: PTR(KLON,KABS)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

!     LOCAL INTEGER SCALARS
INTEGER_M :: JA, JL, IEXP, JX


IEXP=NEXPO3(KNU)
!!print *,'KIDIA, KFDIA, KLON, KNU, KABS=',KIDIA, KFDIA, KLON, KNU, KABS
!!print *,'IEXP(',KNU,')=',IEXP
!!print *,(REXPO3(KNU,1,JX),JX=1,IEXP)
!!print *,(REXPO3(KNU,2,JX),JX=1,IEXP)
!!print *,'PU(KLON,KABS)=', PU

DO JA = 1,KABS
  DO JL=KIDIA,KFDIA
    PTR(JL,JA)=_ZERO_
  END DO
    
!
! JP Pinty: a MIN() is added for caution in the exponential
!
  DO JX=1,IEXP
    DO JL = KIDIA,KFDIA
      PTR(JL,JA) = PTR(JL,JA) &
        &+REXPO3(KNU,1,JX)*EXP(-MIN(100.0,REXPO3(KNU,2,JX)*PU(JL,JA)))
    END DO
  END DO    
ENDDO

RETURN
END SUBROUTINE SWUVO3
