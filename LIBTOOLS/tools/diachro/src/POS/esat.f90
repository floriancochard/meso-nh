!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.esat.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      FUNCTION ESAT(PT)
!     #################
!
!!****  *ESAT* - Computes the saturation water vapor pressure
!!            
!!
!!    PURPOSE
!!    -------
!       Computes the saturation water vapor pressure at a given temperature,
!      used in the emagram routine of TRACE
!
!!**  METHOD
!!    ------ 
!!      Analytical formula of Tetens (1930)
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
!*       0.1   Declaration of argument and result
!
REAL,INTENT(IN)                :: PT
REAL                           :: ESAT
!
!*       0.2   Declaration of local variables
!
!
REAL        :: ZABZ, ZTC
!
!-------------------------------------------------------------------------------
!
!*       1.    CALCULATION OF ESAT
!              -------------------
!
! ESAT (MILLIBARS), PT (KELVIN)
!
ZABZ=273.16
ZTC = PT-ZABZ
ESAT = 6.1078*EXP((17.2693882*ZTC)/(ZTC+237.3))
!
!------------------------------------------------------------------------------
!
!*       2.     EXIT
!               ----
!
RETURN
END FUNCTION ESAT
