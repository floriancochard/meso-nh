!     ######spl
      MODULE  MODN_NCAR
!     #################
!
!!****  *MODN_NCAR* - defines the NAM_DIRTRA_POS namelist (former NCAR common)
!!
!!    PURPOSE
!!    -------
!      This declarative module defines the NAM_DIRTRA_POS namelist, which
!     contains the parameters controlling the NCAR plotting environnement
!     parameters.
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_FIELD1_CV2D), to appear in 1994 
!!
!!     NCAR Graphics Technical documentation, UNIX version 3.2,
!!     Scientific computing division, NCAR/UCAR, Boulder, USA.
!!      Volume 1: Fundamentals, Vers. 1, May 1993
!!      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        02/06/94
!!     updated   PM    19/11/94
!!     JS change the pressure variable 25/07/97
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!
IMPLICIT NONE
!
!*     0.1  Former namelist NAM_DIRTRA_POS
!
INTEGER            :: NIOFFD,  &  ! Label normalisation (=0 none, =/=0 active)
                      NULBLL,  &  ! Nb of contours between 2 labelled contours
                      NIOFFM,  &  ! =0    --> message at picture bottom
                                  ! =/= 0 --> no message
                      NIOFFP,  &  ! Special point value detection
                                  ! (=0 none, =/=0 active)
                      NHI,     &  ! Extrema detection
                                  ! (=0 --> H+L, <0 nothing)
                      NINITA,  &  ! For streamlimes
                      NINITB,  &  ! (Not yet implemented)
                      NIGRNC,  &  ! 
                      NDOT,    &  ! Line style
                                  ! (=0|1|1023|65535 --> solid lines;
                                  !  <0 --> solid lines for positive values and
                                  !  dotted lines(ABS(NDOT))for negative values;
                                  !  >0 --> dotted lines(ABS(NDOT)) )
                      NIFDC,   &  ! Coastline data style (0 none, 1 NCAR, 2 IGN)
                      NLPCAR,  &  ! Number of land-mark points to be plotted
                      NIMNMX,  &  ! Contour selection option
                                  ! (=-1 Min, max and inc. automatically set;
                                  !  =0 Min, max automatically set; inc. given;
                                  !  >0 Min, max, inc. given by user)
		      NISKIP      ! Rate for drawing velocity vectors
! Nov 2000
INTEGER            :: NIJCAR=0    ! Cartes. Equivalent de NLPCAR en proj. cart.
		    
CHARACTER(LEN=8)   :: CTYPHOR     ! Horizontal cross-section type
                                  ! (='K' --> model level section;
                                  !  ='Z' --> constant-altitude section;
                                  !  ='P' --> isobar section (planned)
                                  !  ='T' --> isentrope section (planned)

REAL               :: XSPVAL,  &  ! Special value
                      XSIZEL      ! Label size
REAL               :: XVHC,XVRL,XAMX

REAL,DIMENSION(100) :: X3DINT, X2DINT

! Nov 2000
REAL,DIMENSION(400) :: XICAR, XJCAR !  En cartesien en indices de grilles
! les precedents sont les equivalents des suivants et leur nb=NIJCAR
! Nov 2000
REAL,DIMENSION(400) :: XLATCAR, XLONCAR ! Lat. and Long. of land-mark points

LOGICAL :: LXY,  &  ! If =.T., plots  a grid-mesh stencil background
           LXZ,  &  ! If =.T., plots  a model-level stencil background 
           LCOLAREA, LCOLAREASEL, LTABCOLDEF, &
           LCOLINE, LCOLINESEL, LISOWHI, LCOLBR, LARROVL, LISO,  &
           LDATFILE, LVECTMNMX, LMINMAX,      &
           LSPOT
!
!*     0.2  Former namelist NAM_DIRTRA2_POS
!
! Gestion taille fenetre affichage
! ********************************
! Cas coupes horizontales (isocontours et vecteurs)
REAL               :: XVPTL, XVPTR, XVPTB, XVPTT   
REAL               :: XWINL, XWINR, XWINB, XWINTT   
!
! Cas coupes verticales (isocontours et vecteurs)
REAL               :: XVPTVL, XVPTVR, XVPTVB, XVPTVT   
REAL               :: XWINVL, XWINVR, XWINVB, XWINVT   
!
! Cas profils verticaux
REAL               :: XVPTPVL, XVPTPVR, XVPTPVB, XVPTPVT   
REAL               :: XWINPVL, XWINPVR, XWINPVB, XWINPVT   
!
! Cas TRAXY
REAL               :: XVPTXYL, XVPTXYR, XVPTXYB, XVPTXYT   
REAL               :: XWINXYL, XWINXYR, XWINXYB, XWINXYT   
!
! CH
LOGICAL :: LVPTUSER, LWINUSER
! CV
LOGICAL :: LVPTVUSER, LWINVUSER
! PV
LOGICAL :: LVPTPVUSER, LWINPVUSER
! XY
LOGICAL :: LVPTXYUSER, LWINXYUSER
!
! Gestion epaisseur traits isocontours (CH et CV)
! ***********************************************
REAL                :: XLWDEF, XLW, XLWVDEF, XLWV
REAL                :: XLWIDTH
!
END MODULE MODN_NCAR
