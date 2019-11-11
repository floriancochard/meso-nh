!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.dewp.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      FUNCTION DEWP(PQ,PP)
!     ####################
!
!!****  *DEWP* - Computes the dewpoint temperature
!!
!!    PURPOSE
!!    -------
!        Computes the dewpoint temperature for given mixing ratio and pressure
!      used for the emagram routine of TRACE
! 
!!**  METHOD
!!    ------ 
!!       Analytical formula inverting the Tetens formula
!!
!!    EXTERNAL
!!    --------
!!      NONE 
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
REAL,INTENT(IN)                :: PQ ! Mixing ratio ( g/kg)
REAL,INTENT(IN)                :: PP ! Pressure (millibars)
!
REAL                           :: DEWP ! Dewpoint temperature (Kelvin)

!
!*       0.2   Declaration of local variables
!
REAL        :: ZX, ZY, ZPQ
!
!-------------------------------------------------------------------------------
!
!*       1.    CALCULATION OF DEWP
!              -------------------
!  PQ (G/KG),  PP (MILIBARS), DEWP (KELVIN)
!
!
ZPQ=PQ
IF(PQ <= 0.)ZPQ=1.E-16
ZX = PP*ZPQ/(622.+ZPQ)
!ZX = PP*PQ/(622.+PQ)
ZY = ALOG(ZX/6.1078)
DEWP = ZY*237.3/(17.2693882-ZY)
!
!-----------------------------------------------------------------------------
!
!*       2.    EXIT
!              ----
!
RETURN
END FUNCTION DEWP
