!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.tsa.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      FUNCTION TSA(POS,PP)
!     ####################
!
!!****  *TSA* -  Computation of the wet-bulb potential temperature
!!            
!!
!!    PURPOSE
!!    -------
!      Computation of the wet-bulb potential temperature from given
!     equivalent potential temperature and pressure used in the
!     emagram routine of TRACE
!
!!**  METHOD
!!    ------ 
!!     Iterative formula 
!!
!!    EXTERNAL
!!    --------
!!      WSOUS: computes the saturation miwing ration at given temperature 
!!             and moisture
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
!!    AUTHOR
!!    ------
!!      - Initial version Peridot TRACE Program, P.Bougeault *Meteo-France*,
!!      modified by R. Benoit (mc2, april 91) for the PYREX Oracle data base.
!!      - Present version J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/07/94 
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
REAL,INTENT(IN)                :: POS, PP
REAL                           :: TSA

!
!*       0.2   Declaration of local variables
!
!
INTEGER     :: I
REAL        :: ZA, ZTQ, ZD, ZX
!
!*       0.3   Declaration of external function interface 
!
INTERFACE
  FUNCTION WSOUS(PT,PP)
  REAL,INTENT(IN)                :: PT, PP
  REAL                           :: WSOUS
  END FUNCTION WSOUS
END INTERFACE
!-------------------------------------------------------------------------------
!
!*       1.    CALCULATION OF TSA
!              ------------------
!
! TSA and OS (KELVIN), PP (MILLIBARS)
! SIGN(ZA,ZB) REPLACES THE ALGEBRAIC SIGN OF ZA WITH THE SIGN OF ZB
!
ZA = POS
ZTQ = 253.16
ZD = 120.
! If the temperature difference ZX is small, exit this loop
DO I = 1,12
  ZD=ZD/2.
  ZX=ZA*EXP(-2.6518986*WSOUS(ZTQ,PP)/ZTQ)-ZTQ*((1000./PP)**.286)
  IF(ABS(ZX).LT.0.01)EXIT
  ZTQ = ZTQ + SIGN(ZD,ZX)
ENDDO
TSA = ZTQ
!
!------------------------------------------------------------------------------
!
!*      2.     EXIT
!              ----
!
RETURN
END FUNCTION TSA
