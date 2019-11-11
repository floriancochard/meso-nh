!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOEAERD


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOEAERD* SPECTRAL DISTRIBUTION OF AEROSOLS
!     ------------------------------------------------------------------

!REAL_B,ALLOCATABLE:: CVDAES(:)
!REAL_B,ALLOCATABLE:: CVDAEL(:)
!REAL_B,ALLOCATABLE:: CVDAEU(:)
!REAL_B,ALLOCATABLE:: CVDAED(:)

REAL_B :: CVDAES(61)
REAL_B :: CVDAEL(61)
REAL_B :: CVDAEU(61)
REAL_B :: CVDAED(61)

REAL_B :: RAESC(66)
REAL_B :: RAESS(55)
REAL_B :: RAELC(66)
REAL_B :: RAELS(55)
REAL_B :: RAEUC(66)
REAL_B :: RAEUS(55)
REAL_B :: RAEDC(66)
REAL_B :: RAEDS(55)

REAL_B :: RCAEOPS
REAL_B :: RCAEOPL
REAL_B :: RCAEOPU
REAL_B :: RCAEOPD
REAL_B :: RCTRBGA
REAL_B :: RCVOBGA
REAL_B :: RCSTBGA
REAL_B :: RCTRPT
REAL_B :: RCAEADK(3)
REAL_B :: RCAEADM
REAL_B :: RCAEROS


!*     *YOEAERD* SPECTRAL DISTRIBUTION OF AEROSOLS.
!                     (TRIANGULAR *T10* TRUNCATION FOR AEROSOLS).

!     R.G AND M.J        E.C.M.W.F.     29/11/82.
!     J.-J. MORCRETTE    E.C.M.W.F.     92/09/24  Adaptation to IFS

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *CAES_*   REAL      *REFERS TO *SEA AEROSOLS.
!     *CAEL_*   REAL      *REFERS TO *LAND AEROSOLS.
!     *CAEU_*   REAL      *REFERS TO *URBAN AEROSOLS.
!     *CAED_*   REAL      *REFERS TO *DESERT AEROSOLS.
!     *C___C*   REAL      *REFERS TO *COS COMPONENT.
!     *C___S*   REAL      *REFERS TO *SIN COMPONENT.
!     *CAEOP_*  REAL      *CONSTANTS USED FOR AEROSOL COMPUTATIONS.
!     *C___S*   REAL      *REFERS TO *SEA AEROSOLS.
!     *C___L*   REAL      *REFERS TO *LAND AEROSOLS.
!     *C___U*   REAL      *REFERS TO *URBAN AEROSOLS.
!     *C___D*   REAL      *REFERS TO *DESERT AEROSOLS.
!     *C__BGA*  REAL      *CONSTANTS USED FOR AEROSOL COMPUTATIONS.
!     *CVDAE_*  REAL      *CONSTANTS USED FOR AEROSOL COMPUTATIONS.(NFLEV+1)
!     *RCTRPT*   REAL      *CONSTANTS USED FOR AEROSOL COMPUTATIONS.
!     *RCAEADK*  REAL      *CONSTANTS USED FOR AEROSOL COMPUTATIONS.
!     *RCAEADM*  REAL      *CONSTANTS USED FOR AEROSOL COMPUTATIONS.
!     *RCAEROS*  REAL      *BACKGROUND VALUE IN ABSENCE OF AEROSOLS.

!     ------------------------------------------------------------------
END MODULE YOEAERD
