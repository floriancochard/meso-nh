!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.echelle.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE ECHELLE(KLEN,PHA)
!     ############################
!
!!****  *ECHELLE* - Sets the arrow scales for the emagram environment
!!
!!    PURPOSE
!!    -------
!
!    This routine initialize the emagram wind vector plotting by invoking
!  the NCAR "DRWVEC" utility (drawing of a single vector). KLEN and PHA
!  are returned to the calling program.
!
!!**  METHOD
!!    ------
!!     The scaling is made is made by converting to the old-fashioned 
!!    NCAR "metacode coordinate", see NCAR documentation volume I, page 345.
!!    A scaling vector is drawn to the page bottom as a visual guidance.
!!    Returned values are: KLEN maximum arrow size which can be plotted 
!!    (given in metacode units), PHA maximum wind modulus which can be 
!!    plotted (given in m/s). Values of KLEN and PHA have to be mutually
!!    consistent.
!!
!!    EXTERNAL
!!    --------
!!      GETSI  : Retrieves the parameters defining the size of the plotter
!!               in the plotter coordinate system. Size assumed between 1 and
!!               2**ISX-1 and 2**ISY-1. This old-fashioned  NCAR routine is
!!               documented in the SSPS reference manual of the Version 2
!!               (not in version 3!) of the NCAR package. We sincerely
!!               apologize for the inconvenience.
!!      GSCLIP : Controls NCAR window clipping.
!!      GETSET : Returns the current mapping of the NCAR user coordinate
!!               onto the current GKS viewport in normalized device coordinate.
!!               See NCAR reference manual volume 1, page 343 for details.
!!      CFUX   : Converts a X  "fractional coordinate" value into its 
!!               X "user coordinate" counterpart. See NCAR manual volume 1, 
!!               page 346 for details.
!!      CFUY   : Converts a Y  "fractional coordinate" value into its
!!               Y "user coordinate" counterpart. See NCAR manual volume 1,
!!               page 346 for details.
!!      FL2INT : Given a coordinate pair in the NCAR user system, returns the 
!!               coresponding coordinate pair in the metacode system;
!!      DRWVEC : Draws a single vector given by two pairs of metacode 
!!               coordinates, CALL  DRWVEC (M1,M2,M3,M4,LABEL,NC), where
!!               (M1,M2) coordinate of arrow base on a 2**15x2**15 grid,
!!               (M3,M4) coordinate of arrow head on a 2**15x2**15 grid,
!!               LABEL   character label to be put above arrow, and
!!               NC      number of character in label. This routine is 
!!               and documented in the VELVECT NCAR sources, but
!!               not really documented elsewhere... Sorry for this!
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
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   11/01/59
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments and results
!
INTEGER, INTENT(OUT) :: KLEN  ! KLEN maximum arrow size which can be plotted
                              ! (given in metacode units)
REAL,    INTENT(OUT) :: PHA   ! PHA maximum wind modulus which can be plotted
                              ! (given in m/s)
!
!*       0.2   Local variables
!
INTEGER            :: ILENGTH, IDUM5, IM1, IM2, IM3, IM4, IPHAS4, IL

CHARACTER(LEN=10)  :: YLABEL
CHARACTER(LEN=1)  :: Y1
REAL               :: ZU, ZV
REAL               :: ZFXMIN,ZFXMAX,ZFYMIN,ZFYMAX,ZUMIN,ZUMAX,ZVMIN,ZVMAX
!
!*       0.3   TRACE interface with the DRWVEC routine of the NCAR package
!
! NOTICE:  The DRWVEC and the NCAR graphical utilities are NOT written
! ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
!          does not follow the Meso-NH usual rules: communication has
!          to be made using the /VEC1/ COMMON stack with  static memory 
!          allocation.  Actually used variables are: 
!          ICTRFG  arrow centering control flag
!          ISX     plotter size along x in plotter units
!          ISY     plotter size along y in plotter units
!          ZMN     plotter size along x in metacode units 
!          ZMX     plotter size along y in metacode units
!
INTEGER           :: ICTRFG, ILAB, IOFFD, IOFFM, ISX, ISY
REAL              :: ASH, EXT, RMN, RMX, SIDE, SIZE, XLT, YBT, ZMN, ZMX
!
COMMON /VEC1/   ASH        ,EXT        ,ICTRFG     ,ILAB       ,  &
IOFFD      ,IOFFM      ,ISX        ,ISY        ,  &
RMN        ,RMX        ,SIDE       ,SIZE       ,  &
XLT        ,YBT        ,ZMN        ,ZMX
!
!*       0.4   Interface declarations
!
INTERFACE
  FUNCTION CFUX (RX)
  REAL  :: RX, CFUX
  END FUNCTION CFUX
END INTERFACE
!
INTERFACE
  FUNCTION CFUY (RY)
  REAL  :: RY, CFUY
  END FUNCTION CFUY
END INTERFACE
INTERFACE
  SUBROUTINE DRWVEC (M1,M2,M3,M4,LABEL,NC)
  INTEGER :: M1,M2,M3,M4,NC
  CHARACTER(LEN=10) LABEL
  END SUBROUTINE DRWVEC
END INTERFACE
!
!---------------------------------------------------------------------------
!
!*      1.      ARROW SCALE CALCULATION
!
!*      1.0     Sets the plotter dimensions in metacode units
!*              and some upper bound wind value
!
ILENGTH=160  ! ILENGTH is the maximum possible arrow length in plotter units
             ! (i.e.: with respect to the 2**10-1 default value)
PHA=80.      ! PHA is the maximum possible wind value corresponding to the
             ! maximum possible arrow size given above. Thes two values have
             ! to be consistent
!
! Retrieves plotter size, first in plotter units
!
CALL GETSI(ISX,ISY)  
ISX=2**(15-ISX)     
ISY=2**(15-ISY)
!
! Converts the maximum possiblble arrow length in metacode units
! (i.e. with respect to 2**15-1)
!
KLEN=ILENGTH*ISX
ZMN=0.
ZMX=FLOAT(KLEN)+0.01
!
!*       1.1    Computes appropriate scale 
!
CALL GSCLIP(0) ! Enables leader writing out of the frame
!
! Prepares header and scale.
! Retrieves current window limits in normalized 
! device coordinate and NCAR user coordinate.
!
CALL GETSET(ZFXMIN,ZFXMAX,ZFYMIN,ZFYMAX,ZUMIN,ZUMAX,ZVMIN,ZVMAX,IDUM5)
!
! Computes the normalized device coordinates of the point located by
! user coordinates (ZFXMAX-0.05,ZFYMIN-0.04)
!
ZU=CFUX(ZFXMAX-0.05)
ZV=CFUY(ZFYMIN-0.04)
!
! Then, convert result to metacode coordinates
!
CALL FL2INT(ZU,ZV,IM1,IM2)
IM3=IM1+KLEN/4
IM4=IM2
IPHAS4=IFIX(PHA/4)
!
!*       1.2    Draws a unit vector under the plot
!
!               
! The unit vector is 1/4 of the maximum possible wind PHA
!
CALL PCGETC('FC',Y1)
!print *,' **echelle Y1',Y1
CALL PCSETC('FC','?')
YLABEL=' '
WRITE(YLABEL,'(I2,'' M/S    '')')IPHAS4
YLABEL=ADJUSTL(YLABEL)
!print *,' ECHELLE AV DRW..',YLABEL
!CALL DRWVEC(IM1,IM2,IM3,IM4,YLABEL(1:LEN_TRIM(YLABEL)),LEN_TRIM(YLABEL))
IL=10
IF(LRS .OR. LRS1)THEN
IL=0
CALL DRWVEC(IM1,IM2,IM3,IM4,YLABEL,IL)
CALL GSLWSC(1.)
CALL PLCHHQ(25.5553226,-1.4807138,YLABEL(1:LEN_TRIM(YLABEL)),7.,0.,0.)
CALL GSLWSC(2.)
ELSE
CALL DRWVEC(IM1,IM2,IM3,IM4,YLABEL,IL)
ENDIF
CALL SFLUSH
!print *,' ECHELLE AP DRW..'
! 
!  Setting the ICTRFG flag controls the arrow centering.
!  Arrow is centered with ICTRFG=0,  and the tail of the 
!  arrow is placed at the grid point location with ICTRFG=1.
!
ICTRFG=1
!
! Window clipping restored after header writing 
!
!CALL GSCLIP(1)
CALL PCSETC('FC',Y1)
!
!----------------------------------------------------------------------------
!
!*       2.      EXIT
!                ----
!
RETURN
!
END SUBROUTINE ECHELLE
