!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.wsous.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      FUNCTION WSOUS(PT,PP)
!     #####################
!
!!****  *WSOUS* - Computes the saturation mixing ratio
!!            
!!
!!    PURPOSE
!!    -------
!       Computes the saturation mixing ratio for a given temperature
!       used in the emagram routine of TRACE
!
!!**  METHOD
!!    ------ 
!!      Explicit analytical formula
!!
!!    EXTERNAL
!!    --------
!!      ESAT : computes the saturation water vapor pressure at a
!!             given temperature
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!
!!      Among many others, see for instance:
!!       Bluestein H. B., 1992, "Synoptic-Dynamic Meteorology in mid-latitudes"
!!       Volume 1, Priciples of Kinematics and Dynamics, Section 4.3, p. 195,
!!       Oxford University Press.
!!
!!
!!    AUTHOR
!!    ------
!!      - Initial version Peridot TRACE Program, P.Bougeault *Meteo-France*,
!!      modified by R. Benoit (mc2, april 91) for the PYREX Oracle data base.
!!      - Present version J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   10/01/95
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments and result
!
REAL,INTENT(IN)                :: PT, PP
REAL                           :: WSOUS
!
!*       0.2   Declaration of local variables
!
REAL        :: ZX
!
!*       0.3   Declaration of external function interfaces
!
INTERFACE
  FUNCTION ESAT(PT)
  REAL,INTENT(IN)                :: PT
  REAL                           :: ESAT
  END FUNCTION ESAT
END INTERFACE
!-------------------------------------------------------------------------------
!
!*       1.    CALCULATION OF WSOUS
!              --------------------
!
! W (GRAMS WATER VAPOR/KILOGRAM DRY AIR), PP (MILLIBARS)
!
ZX = ESAT(PT)
WSOUS = 622.*ZX/(PP-ZX)
IF(PT.GE.999.)WSOUS = 0.
!
!------------------------------------------------------------------------------
!
!*       2.    EXIT
!              ----
!
RETURN
END FUNCTION WSOUS
