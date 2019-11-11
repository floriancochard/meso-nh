!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
MODULE MODD_COMPPAR
IMPLICIT NONE 
! Debug mode : set LPDEBUG to .TRUE.
LOGICAL,PARAMETER :: LPDEBUG = .FALSE.


! contains coding parameters for (de)compress routines

INTEGER,PARAMETER :: JPCSTENCOD = 1 ! constant array 
INTEGER,PARAMETER :: JPSOPENCOD = 2 ! second order packing 
INTEGER,PARAMETER :: JPEXTENCOD = 3 ! second order packing with min/max values excluded

! Extended code when JPEXTENCOD enabled
!
! BE CAREFUL : 3 bits are reserved for coding this code => max value is 7
INTEGER,PARAMETER :: JPCONST   = 0 ! constant value array
INTEGER,PARAMETER :: JPNORM    = 1 ! same as JPSOPENCOD
INTEGER,PARAMETER :: JPMINEXCL = 2 ! Min value is isolated
INTEGER,PARAMETER :: JPMAXEXCL = 3 ! Max value is isolated
INTEGER,PARAMETER :: JPMINMAXEXCL = 4 ! Min&Max values are isolated
INTEGER,PARAMETER :: JP2VAL       = 5 ! 2 different values in array
INTEGER,PARAMETER :: JP3VAL       = 6 ! 3 different values in array
INTEGER,PARAMETER :: JPOTHER      = 7 ! for future use
INTEGER,PARAMETER :: JPLOG        = 8
END MODULE MODD_COMPPAR
