!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.os.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      FUNCTION OS(PT,PP)
!     ##################
!
!!****  *OS* - Computes the equivalent potential temperature
!!
!!    PURPOSE
!!    -------
!       Computes the equivalent potential temperature
!      at a given temperature and pressure, used in the 
!      emagram plotting utility of TRACE.
!
!!**  METHOD
!!    ------ 
!!      Explicit analytical formula. 
!!
!!    EXTERNAL
!!    --------
!!      WSOUS: computes the saturation mixing ratio at a given
!!             temperature and pressure
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
!*       0.1   Declarations of argument and result
!
REAL,INTENT(IN)                :: PT, PP
REAL                           :: OS

!
!*       0.2   Declaration of external function interface
!
INTERFACE
  FUNCTION WSOUS(PT,PP)
  REAL,INTENT(IN)                :: PT, PP
  REAL                           :: WSOUS
  END FUNCTION WSOUS
END INTERFACE
!-------------------------------------------------------------------------------
!
!*       1.    CALCULATION OF OS
!              -----------------
!
! OS and PT (KELVIN), PP  (MILLIBARS)
!
OS = PT*((1000./PP)**.286)/(EXP(-2.6518986*WSOUS(PT,PP)/PT))
!
!------------------------------------------------------------------------------
!
!*       2.    EXIT
!              ----
!
END FUNCTION OS
