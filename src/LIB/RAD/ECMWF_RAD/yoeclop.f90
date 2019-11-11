!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOECLOP


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     * YOECLOP* PARAMETERS FOR CLOUD OPTICAL PROPERTIES
!     ------------------------------------------------------------------

REAL_B :: RYFWCA(4)
REAL_B :: RYFWCB(4)
REAL_B :: RYFWCC(4)
REAL_B :: RYFWCD(4)
REAL_B :: RYFWCE(4)
REAL_B :: RYFWCF(4)
REAL_B :: REBCUA(4)
REAL_B :: REBCUB(4)
REAL_B :: REBCUC(4)
REAL_B :: REBCUD(4)
REAL_B :: REBCUE(4)
REAL_B :: REBCUF(4)
REAL_B :: RASWCA(4)
REAL_B :: RASWCB(4)
REAL_B :: RASWCC(4)
REAL_B :: RASWCD(4)
REAL_B :: RASWCE(4)
REAL_B :: RASWCF(4)
REAL_B :: REBCUG
REAL_B :: REBCUH
REAL_B :: REBCUI(6)
REAL_B :: REBCUJ(6)
REAL_B :: REFFIA
REAL_B :: REFFIB
REAL_B :: RTIW
REAL_B :: RRIW


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!*    FOUQUART (1987) WATER CLOUD OPTICAL PROPERTIES

! RYFWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RYFWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! RYFWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RYFWCF :  REAL   : ASSYMETRY FACTOR

!*    SLINGO (1989) WATER CLOUD OPTICAL PROPERTIES

! RASWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RASWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! RASWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
! RASWCF :  REAL   : ASSYMETRY FACTOR

!*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM EBERT-CURRY (1992)

! REBCUA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! REBCUB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
! REBCUC :  REAL   : 1-C3  IN SINGLE SCATTERING ALBEDO FORMULA
! REBCUD :  REAL   : C4 IN SINGLE SCATTERING ALBEDO FORMULA
! REBCUE :  REAL   : C5 IN ASSYMETRY FACTOR FORMULA
! REBCUF :  REAL   : C6 IN ASSYMETRY FACTOR FORMULA
! REBCUG :  REAL   : C7 IN MASS ABSORPTION COEFFICIENT FORMULA
! REBCUH :  REAL   : C8 IN MASS ABSORPTION COEFFICIENT FORMULA
! REBCUI :  REAL   : C7 IN LW MASS ABSORPTION COEFFICIENT FORMULA
! REBCUJ :  REAL   : C8 IN LW MASS ABSORPTION COEFFICIENT FORMULA

! REFFIA :  REAL   : C9  IN EFFECTIVE RADIUS FORMULA
! REFFIB :  REAL   : C10 IN EFFECTIVE RADIUS FORMULA

!*    TRANSITION BETWEEN LIQUID AND SOLID WATER

! RTIW   :  REAL   : TEMPERATURE THRESHOLD
! RRIW   :  REAL   : TRANSITION RANGE
!     -----------------------------------------------------------------
END MODULE YOECLOP