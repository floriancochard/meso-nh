!     ######spl
      SUBROUTINE CCOLR(XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C     #################################################
C
C
CC****  *CCOLR* - Performs color filling of the contour intervals
CC
CC    PURPOSE
CC    -------
C       When contour plot is drawn, the successive contour intervals are
C     filled with colors as given by a color index which is a function of
C     contor level.
C
CC**  METHOD
CC    ------
CC
CC      In IMAGE, IMAGEv or IMCOU.., as the contour plots are prepared, the
CC    areas between successive contour levels are identified using "area 
CC    numbers". CCOLR uses these area numbers to select a color in the
CC    current color table, and fills the corresponding area using a GKS
CC    fill area call.  See the NCAR manual to understand how "area numbers"
CC    work, this topic is slightly involved.. (NCAR contouring tutorial, 
CC    Vol. 2, pages 12-19, page 120,  and pages 130-133).
CC
CC      To summarize, all the lines composing a plot are grouped by "edge
CC    groups" which may be individually accessed using "group numbers" to
CC    perform specific tasks. For the present purpose only the lines drawn
CC    by CONPACK are important, and they belong to group number 3.
CC      When the contours are computed, CONPACK  assigns "area numbers" to the
CC    different sub-regions of the plot: typically screen points out of the
CC    model domain are given a negative  area number,  and areas between
CC    isocontours receive area numbers greater than 2, with increasing area
CC    numbers from the lower contour to the higher one.
CC      The coloring is therefore performed by scanning the group and area
CC    numbers to locate the screen locations to be colored, as follows:
CC    - CCOLR is called by CONPACK for each contour polygon, with XWRK-YWRK
CC    containing the NWRK points of the current contour, and IAREA-IGRP
CC    containing the corresponding group and area numbers;
CC    - First, the group number is checked to select CONPACK items only,
CC    - Second, the area numbers are checked to select positive ones, and
CC    a color values are picked up in the ICOL color table.
CC    - If so, the color parameter is set (GSFACI) and the color filling 
CC    routine is called to fill the current contour (XWRK-YWRK) with the 
CC    prescribed color.
CC
CC NOTICE:    CCOLR and the NCAR graphical utilities are NOT written
CC ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
CC          does not follow the Meso-NH usual rules: it has to be directly
CC          called by the NCAR CONPACK utility.
CC
CC    EXTERNAL
CC    --------
CC     None
CC
CC    EXPLICIT ARGUMENTS
CC    ------------------
CC
CC       XWRK : x-coordinates (in NCAR fractional system) of the successive
CC              points forming a given contour enclosing a polygonal area.
CC       YWRK : y-coordinates (in NCAR fractional system) of the successive
CC              points forming a given contour enclosing a polygonal area.
CC       NWRK : Number of points in XWRK-YWRK to build the contour.
CC       IAREA: Area identifiers for the polygon defined by the XWRK-YWRK and
CC              for each of the NGRPS groups of edges in this plot.
CC       IGRP : Group identifiers for the polygon defined by the XWRK-YWRK and
CC              for each of the NGRPS groups of edges in this plot.
CC       NGRPS: Maximum number of edge groups defined in this plot.
CC
CC       NOTICE: All these dummy arguments are required
CC       ------  by the NCAR CALLS
CC
CC    IMPLICIT ARGUMENTS
CC    ------------------
CC
CC       Common COLAREA : color table information
CC         ICOL  : Array of the possible values of the GKS color index. These
CC                 GKS color index values are initialized earlier in the TRACE
CC                 run by reading a user provided color table file.
CC
CC    REFERENCE
CC    ---------
CC
CC      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
CC       + Book1: Concepts and Fundamentals, to appear in 1994;
CC       + Book2: Technical Reference and Flowcharts, to appear in 1994;
CC       + Book3: Tutorial, November 1994.
CC
CC     NCAR Graphics Technical documentation, UNIX version 3.2,
CC     Scientific computing division, NCAR/UCAR, Boulder, USA.
CC      Volume 1: Fundamentals, Vers. 1, May 1993
CC      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
CC
CC    AUTHOR
CC    ------
CC
CC      J. Duron    * Laboratoire d'Aerologie *
CC
CC    MODIFICATIONS
CC    -------------
CC      Original       01/07/94
CC      Updated   PM   24/01/95
C-------------------------------------------------------------------------------
C
C*     0.   DECLARATIONS
C           ------------
C
C>>>>>>>DRAGOON NOTICE: I ENFORCED "IMPLICIT NONE" IT'S WISE CHECKING...
C
      IMPLICIT NONE
C
C*     0.0  Dummy arguments
C
      REAL XWRK(*), YWRK(*)
      INTEGER IAREA(*), IGRP(*)
      INTEGER NWRK,NGRPS
C
C*     0.1  Commons
C
      COMMON/COLAREA/ICOL(300)
      INTEGER ICOL
C
C*     0.2  Local variables
C

      REAL RSCR(10000)
      INTEGER ISCR(10000)
      INTEGER I,IA
C
C-----------------------------------------------------------------------------
C
C*     1.    PERFORMS CONTOUR INTERVAL COLORING
C            ----------------------------------
C
C*     1.1   Select a color index for each area number
C*           when edge group=3 (contour edges) and area number
C*           is positive (within plot limits)
C
      DO I=1,NGRPS
C     print *,' IGRP IAREA',IGRP(I),IAREA(I),'  I',I
      IF(IGRP(I).EQ.3)THEN
        IA=-5
        IF(IAREA(I).GT.0)IA=ICOL(IAREA(I))
C       IF(IAREA(I).GT.0)IA=IAREA(I)+2
      END IF
      ENDDO
C
C*     1.2   Fills the (XWRK,YWRK) polygon with selected color
C
      IF(IA.GT.0)THEN
        CALL GSFACI(IA)
        CALL GFA(NWRK-1,XWRK,YWRK)
      ENDIF
      RETURN
C
C----------------------------------------------------------------------------
C
C*    2.    EXIT
C           ----
C
      END
