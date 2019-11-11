!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.fleche.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE FLECHE(PX,PY,PU,PV,KLEN,PHA)
!     #######################################
!
!!****  *FLECHE* - Draws a single arrow for emagram wind display
!!
!!    PURPOSE
!!    -------
!
!    This routine draws an emagram wind vector by invoking the NCAR 
!  "DRWVEC" utility (drawing of a single vector). The wind arrow is
!  drawn in the appropriate direction and location for the emagram
!  environment. KLEN and PHA are input only scaling factors received 
!  from the "ECHELLE" routine.
!
!
!!**  METHOD
!!    ------
!!      A simple call to DRWVEC, which has stand after scaling by
!!  "ECHELLE" to set KLEN and PHA.
!!
!!   NOTICE:  The DRWVEC and the NCAR graphical utilities are NOT written
!!   ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
!!            does not follow the Meso-NH usual rules: communication has
!!            to be made using the /VEC1/ COMMON stack with  static memory
!!            allocation. See the ECHELLE routine for details.
!!
!!    EXTERNAL
!!    --------
!!      FL2INT : Given a coordinate pair in the NCAR user system, returns the
!!               corresponding coordinate pair in the metacode system;
!!      VVSETI : Sets an integer NCAR parameter to select an option in the
!!               NCAR vector environment
!!      DRWVEC : Draws a single vector given by two pairs of metacode
!!               coordinates, CALL  DRWVEC (M1,M2,M3,M4,LABEL,NC), where
!!               (M1,M2) coordinate of arrow base on a 2**15x2**15 grid,
!!               (M3,M4) coordinate of arrow head on a 2**15x2**15 grid,
!!               LABEL   character label to be put above arrow, and
!!               NC      number of character in label. This routine is
!!               given and documented in the VELVECT NCAR sources, but
!!               not really documented elsewhere... Sorry for this!
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!     MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!     NCAR Graphics Technical documentation, UNIX version 3.2,
!!     Scientific computing division, NCAR/UCAR, Boulder, USA.
!!      Volume 1: Fundamentals, Vers. 1, May 1993
!!      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
!!
!!     For the vector utilities not documented in the NCAR package
!!     Version 3 idocumentation, a better reference is:
!!      The NCAR GKS-Compatible Graphics System Version 2,
!!      SPPS an NCAR System Plot Package Simulator.
!!      NCAR Technical note 267+1A, April 1986, NCAR/UCAR, Boulder, USA.
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
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
!*       0.1   Dummy arguments and results
!
INTEGER           :: KLEN            ! Maximum arrow size which can be 
                                     ! plotted (given in metacode units)
REAL              :: PX, PY          ! Arrow tail location, given in NCAR 
                                     ! user coordinate system.
REAL              :: PU, PV          ! Wind components U and V to be plotted,
                                     ! given in m/s.
REAL              :: PHA             ! Maximum wind modulus which can be 
                                     ! plotted (given in m/s). Values of KLEN 
                                     ! and PHA have to be mutually consistent.
!
!*       0.2   Local variables  
!
INTEGER           :: IM1, IM2, IM3, IM4 ! Tail and head locations of the
                                        ! arrow, given in metacode coordinates
CHARACTER(LEN=10) :: YLABEL='AAAAAAAAAA'             ! Arrow label (i.e.: its scale)
!
INTERFACE
  SUBROUTINE DRWVEC (M1,M2,M3,M4,LABEL,NC)
   CHARACTER*10 LABEL
   INTEGER :: M1,M2,M3,M4,NC
  END SUBROUTINE DRWVEC
END INTERFACE
!-------------------------------------------------------------------------------
!
!*       1.    ARROW DRAWING
!              -------------
! 
!*       1.1   Converts tail location from user to metacode coordinates
!*                     (also called fractional) coordinates
!
CALL FL2INT(PX,PY,IM1,IM2)
!
!*       1.2   Computes the head location in metacode coordinates
!
IM3=IM1+INT(PU*FLOAT(KLEN)/PHA)
IM4=IM2+INT(PV*FLOAT(KLEN)/PHA)
!
!*       1.3   Draws the arrow
!
! Setting VPO >0, the tail of the vector arrow is 
! placed at the grid point location
!
CALL VVSETI('VPO',1)
!
! As the last argument for DRWVEC 
! is 0, no label is actually written
!
CALL DRWVEC(IM1,IM2,IM3,IM4,YLABEL,0)
!      CALL PWRITX(PU,PV,6H'KGU'-,6,10,0,0)
!
!------------------------------------------------------------------------------
!
!*       2.    EXIT
!              ----
!
RETURN
!
END SUBROUTINE FLECHE
