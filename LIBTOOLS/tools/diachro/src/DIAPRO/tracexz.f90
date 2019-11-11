!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.tracexz.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE TRACEXZ
!     ##################
!
!!****  *TRACEXZ* - Overlays a gridpoint location stencil over a
!!                  West-East vertical  cross-section plot.
!!
!!    PURPOSE
!!    -------
!
!       When LXZ=.T., and in the special case of a vertical cross-section 
!    located using the grid index format, shows a model level stencil 
!    overlaid on the plot.
!
!!**  METHOD
!!    ------
!!      Draws polylines between gridpoints corresponding to the NMGRID value.
!!
!!    EXTERNAL
!!    --------
!!      GSLN       : NCAR routine to set a line type.
!!      GPL        : NCAR routine to draw a polyline.
!!      VALNGRID   : loads current grid number in the NMGRID global variable
!!      COMPCOORD  : computes true altitudes for NMGRID grid location
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_NMGRID : declares global variable  NMGRID
!!         NMGRID  : Current MESO-NH grid indicator
!!
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!         XXX : XHAT coordinate values for all the MESO-NH grids
!!
!!      Module MODN_PARA   : defines NAM_DOMAIN_POS namelist 
!!         Module MODD_DIM1  : contains dimension of data array
!!                  NIMAX,NKMAX    :  x, and z array dimensions
!!
!!      Module MODD_GRID1  : declares grid variables (Model module)
!!         XZZ : true z altitude for the current NMGRID grid location
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
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
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   01/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODD_COORD
USE MODN_PARA            !NOTICE: MODN_PARA includes MODD_DIM1
USE MODD_GRID1
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Local variables
!
INTEGER           :: JKLOOP,JILOOP 
INTEGER           :: IIU, IKU, IMGRID
!
REAL,DIMENSION(200) :: ZX, ZY
!
!-------------------------------------------------------------------------------
!
!*      1.     MODEL LEVELS STENCIL DRAWING
!              ----------------------------
!
IIU=NIMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
!
CALL GSLN(3)
!
!*      1.1    Draws the "w" level stencil
!
!print *,' Tracexz NMGRID ',NMGRID
IMGRID=NMGRID
CALL COMPCOORD_FORDIACHRO(4)      ! computes NMGRID grid true altitudes
!print *,' Tracexz IMGRID ',IMGRID
!CALL VALNGRID('WM')
DO JKLOOP=1,IKU
  DO JILOOP=1,IIU
    ZX(JILOOP)=XXX(JILOOP,4)-XXX(NIDEBCOU,IMGRID)
    ZY(JILOOP)=XZZ(JILOOP,NJDEBCOU,JKLOOP)
  ENDDO
  CALL GPL(IIU,ZX,ZY)
ENDDO
!
!*      1.2   Draws the NMGRID model level stencil
!
NMGRID=IMGRID
CALL COMPCOORD_FORDIACHRO(NMGRID)      ! computes NMGRID grid true altitudes
!print *,' Tracexz NMGRID ',NMGRID
!
IF(NMGRID.EQ.4)CALL GSLN(3)
IF(NMGRID.EQ.2)CALL GSLN(2)
IF(NMGRID.EQ.3)CALL GSLN(4)
IF(NMGRID.EQ.1)CALL GSLN(5)
!
DO JKLOOP=1,IKU
  DO JILOOP=1,IIU
    ZX(JILOOP)=XXX(JILOOP,NMGRID)-XXX(NIDEBCOU,NMGRID)
    ZY(JILOOP)=XZZ(JILOOP,NJDEBCOU,JKLOOP)
  ENDDO
  CALL GPL(IIU,ZX,ZY)
ENDDO
!
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
!
CALL GSLN(1)
!
RETURN
END SUBROUTINE TRACEXZ
