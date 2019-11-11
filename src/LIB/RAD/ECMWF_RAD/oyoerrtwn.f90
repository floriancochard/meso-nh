!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE OYOERRTWN


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER_M , DIMENSION(16) :: NG
INTEGER_M , DIMENSION(16) :: NSPA
INTEGER_M , DIMENSION(16) :: NSPB

REAL_B , DIMENSION(16) :: WAVENUM1
REAL_B , DIMENSION(16) :: WAVENUM2
REAL_B , DIMENSION(16) :: DELWAVE

REAL_B , DIMENSION(181,16) :: TOTPLNK
REAL_B , DIMENSION(181)    :: TOTPLK16

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! TOTPLNK : REAL    :
! TOTPLK16: REAL    :
!     -----------------------------------------------------------------
END MODULE OYOERRTWN
