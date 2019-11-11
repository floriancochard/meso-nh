!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.tracexy.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE TRACEXY
!     ##################
!
!!****  *TRACEXY* - Overlays a gridpoint location stencil over a
!!                  horizontal cross-section plot.
!!
!!    PURPOSE
!!    -------
!       When LXY=.T. shows the gridpoint location on horizontal
!    cross-section plots.
!
!!**  METHOD
!!    ------
!!      Draws polylines between gridpoints corresponding to the NMGRID value.
!!
!!    EXTERNAL
!!    --------
!!      GSLN : NCAR routine to set a line type.
!!      GPL  : NCAR routine to draw a polyline.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_NMGRID : declares global variable  NMGRID
!!         NMGRID  : Current MESO-NH grid indicator
!!
!!      Module MODD_OUT    : Defines a log. unit for printing
!!         NIMAXT, NJMAXT:  Size of the displayed window within a 
!!                          MESO-NH field arrays
!!
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!         XXX,XXY : x, y coordinate values for all the MESO-NH grids
!!
!!      Module MODD_DIM1       : Contains dimensions
!!         NIMAX,NJMAX :  x, and y array dimensions
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODD_OUT
USE MODD_COORD
USE MODD_DIM1

IMPLICIT NONE
!
!*       0.1   Variables locales

INTEGER           :: JJLOOP, JILOOP

REAL,DIMENSION(2) :: ZX, ZY
!
!-------------------------------------------------------------------------------
!
!*      1.    GRIDPOINT STENCIL DRAWING
!             -------------------------
!
CALL GSLN(3)
!
!*      1.1    Draws a "w" grid stencil
!
DO JILOOP=1,NIMAXT
  ZX(1)=XXX(NIINF+JILOOP-1,4)
  ZX(2)=XXX(NIINF+JILOOP-1,4)
  ZY(1)=XXY(NJINF,4)
  ZY(2)=XXY(NJSUP,4)
  CALL GPL(2,ZX,ZY)
ENDDO
!
DO JJLOOP=1,NJMAXT
  ZX(1)=XXX(NIINF,4)
  ZX(2)=XXX(NISUP,4)
  ZY(1)=XXY(NJINF+JJLOOP-1,4)
  ZY(2)=XXY(NJINF+JJLOOP-1,4)
  CALL GPL(2,ZX,ZY)
ENDDO
!
!*      1.2   Draws the NMGRID grid stencil
!
IF(NMGRID.EQ.4)CALL GSLN(3)
IF(NMGRID.EQ.2)CALL GSLN(2)
IF(NMGRID.EQ.3)CALL GSLN(4)
IF(NMGRID.EQ.1)CALL GSLN(5)
!
DO JILOOP=1,NIMAXT
  ZX(1)=XXX(NIINF+JILOOP-1,NMGRID)
  ZX(2)=XXX(NIINF+JILOOP-1,NMGRID)
  ZY(1)=XXY(NJINF,NMGRID)
  ZY(2)=XXY(NJSUP,NMGRID)
  CALL GPL(2,ZX,ZY)
ENDDO
!
DO JJLOOP=1,NJMAXT
  ZX(1)=XXX(NIINF,NMGRID)
  ZX(2)=XXX(NISUP,NMGRID)
  ZY(1)=XXY(NJINF+JJLOOP-1,NMGRID)
  ZY(2)=XXY(NJINF+JJLOOP-1,NMGRID)
  CALL GPL(2,ZX,ZY)
ENDDO
!
!*      2.    EXIT
!             ----
!
CALL GSLN(1)
!
RETURN
END SUBROUTINE TRACEXY
