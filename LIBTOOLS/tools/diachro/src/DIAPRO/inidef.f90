!     ######spl
      SUBROUTINE INIDEF
!     #################
!
!!****  *INIDEF* - Sets defaults values of TRACE namelists variables
!!
!!    PURPOSE
!!    -------
!      Sets defaults values of TRACE namelists variables
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST  : declares model physical constants
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist 
!!                         (former NCAR common)
!!
!!       NIOFFD     : Label normalisation (=0 none, =/=0 active)
!!       NULBLL     : Nb of contours between 2 labelled contours
!!       NIOFFM     : =0    --> message at picture bottom
!!                    =/= 0 --> no message
!!       NIOFFP     : Special point value detection
!!                    (=0 none, =/=0 active)
!!       NHI        : Extrema detection
!!                    (=0 --> H+L, <0 nothing)
!!       NINITA     : For streamlimes
!!       NINITB     : Not yet implemented
!!       NIGRNC     : Not yet implemented
!!       NDOT       : Line style
!!                    (=0|1|1023|65535 --> solid lines;
!!                    <0 --> solid lines for positive values and
!!                    dotted lines(ABS(NDOT))for negative values;
!!                    >0 --> dotted lines(ABS(NDOT)) )
!!       NIFDC      : Coastline data style (0 none, 1 NCAR, 2 IGN)
!!       NLPCAR     : Number of land-mark points to be plotted
!!       NIMNMX     : Contour selection option
!!                    (=-1 Min, max and inc. automatically set;
!!                    =0 Min, max automatically set; inc. given;
!!                    >0 Min, max, inc. given by user)
!!       NISKIP     : Rate for drawing velocity vectors
!!       CTYPHOR    : Horizontal cross-section type
!!                    (='K' --> model level section;
!!                     ='Z' --> constant-altitude section;
!!                     ='P' --> isobar section (planned)
!!                     ='T' --> isentrope section (planned)
!!       XSPVAL     : Special value
!!       XUINT      : Increment contour value for UM, UT
!!       XVINT      : Increment contour value for VM, VT
!!       XWINT      : Increment contour value for WM, WT
!!       XTHINT     : Increment contour value for THM,THT
!!       XPABSINT   : Increment contour value for PABSM, PABST
!!       XSIZEL     : Label size
!!       XLATCAR, XLONCAR :  Lat. and Long. of land-mark points
!!       LXY        : If =.TRUE., plots  a grid-mesh stencil background
!!       LXZ        : If =.TRUE., plots  a model-level stencil background 
!!
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist 
!!                          (former PARA common)
!!
!!       XIDEBCOU, XJDEBCOU : Origin of a vertical cross-section
!!                            in cartesian (or conformal) real values
!!       XHMIN      : Altitude of the vert. cross-section
!!                    bottom (in meters above sea-level)
!!       XHMAX      : Altitude of the vert. cross-section
!!                    top (in meters above sea-level)
!!
!!      Module MODD_ALLVAR 
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
!!     NCAR Graphics Technical documentation, UNIX version 3.2,
!!     Scientific computing division, NCAR/UCAR, Boulder, USA.
!!      Volume 1: Fundamentals, Vers. 1, May 1993
!!      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
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
USE MODN_NCAR
USE MODN_PARA
USE MODD_CST
USE MODD_ALLVAR
USE MODD_RESOLVCAR, ONLY : XISOLEV, XLW1, XLW2, XLW3, XLW4

IMPLICIT NONE
!
!*       0.1   Local variables
!              ---------------

LOGICAL :: LSUPER   !  TO BE COMPLETED <<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!------------------------------------------------------------------------------
!
!*       1.   SETS DEFAULTS FOR THE NAMELISTS VARIABLES
!             -----------------------------------------
!
NIINF=0; NJINF=0; NISUP=0; NJSUP=0
XISOLEV(:)=9999.; XLATCAR(:)=9999.; XLONCAR(:)=9999.; XICAR(:)=9999.;
XJCAR(:)=9999.
!
NIOFFD=0
NULBLL=1
NIOFFM=1
XSIZEL=0.02
NIOFFP=1
XSPVAL=999.
XHMIN=0.
XHMAX=0.
NHI=-1
NDOT=-21845
NIFDC=1
NIGRNC=1
NLPCAR=0
XLATCAR(1)=44.52
XLONCAR(1)=.3
NINITA=2
NINITB=2
NIMNMX=0
CTYPHOR='Z'
LSUPER=.FALSE.
LCOLAREA=.FALSE.
LSPOT=.FALSE.
LCOLBR=.TRUE.
LTABCOLDEF=.TRUE.
LCOLAREASEL=.FALSE.
LCOLINESEL=.FALSE.
LISOWHI=.FALSE.
LCOLINE=.FALSE.
LSPOT=.FALSE.
NISKIP=1
XIDEBCOU=-999.
XJDEBCOU=-999.
LXY=.FALSE.
LXZ=.FALSE.
NVAR3D=0
NVAR2D=0
X3DINT(:)=0.
X2DINT(:)=0.
XAMX=.2
XVHC=0.
XVRL=0.
LARROVL=.FALSE.
LISO=.TRUE.
LVECTMNMX=.FALSE.
LVPTUSER=.FALSE.
XVPTL=.1
XVPTR=.9
XVPTB=.1
XVPTT=.9
LVPTVUSER=.FALSE.
XVPTVL=.1
XVPTVR=.9
XVPTVB=.1
XVPTVT=.9
LVPTPVUSER=.FALSE.
XVPTPVL=.13
XVPTPVR=.9
XVPTPVB=.1
XVPTPVT=.9
LMINMAX=.FALSE.
LDATFILE=.TRUE.
XLWDEF=1.
XLWVDEF=0.5
!
XLW=-1.; XLW1=-1.; XLW2=-1.; XLW3=-1.; XLW4=-1.; XLWV=-1.
NIMNMX=-1

!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE INIDEF
