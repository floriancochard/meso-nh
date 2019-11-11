!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modn/s.modn_para.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE  MODN_PARA
!     #################
!
!!****  *MODN_PARA* - defines the NAM_DOMAIN_POS namelist (former PARA common)
!!
!!    PURPOSE
!!    -------
!       This declarative module declares the variables of the NAM_DOMAIN_POS
!      namelist, which specify all the geometrical characteristics of the
!      plotted domain as requested by the user.
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module MODD_DIM1 : contains dimensions of data arrays
!!       NIINF, NISUP   : lower and upper bounds of arrays to be plotted in
!!                        x direction
!!       NJINF, NJSUP   : lower and upper bounds of arrays to be plotted in
!!                        y direction
!!
!!    REFERENCE
!!    ---------
!!     Bougeault et al., 1994, "The MESO-NH user's guide", Chapter 4: Run a
!!     post-processing session, Internal technical note, CNRM/GMME, Toulouse
!!
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODN_PARA), to appear in 1994 
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
!!     updated   PM    21/11/94
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!
USE MODD_DIM1

IMPLICIT NONE

LOGICAL,SAVE     :: LHORIZ,  & ! =.T. to perform horizontal cross-sections
                               ! (LVERTI must be = to .F.)
                    LVERTI,  & ! =.T. to perform vertical cross-sections, in-
                               ! -cluding vert. 1D profiles. (LHORIZ must be
                               ! = to .F.)
                    L3D        ! =.T. to draw 3D perspective plots (LHORIZ and 
                               ! LVERTI must be = to .F.).(Not yet implemented)

INTEGER,SAVE     :: NIDEBCOU,NJDEBCOU,  & ! Origin of a vertical cross-section
                                          ! in grid index integer values
                                          ! (XIDEBCOU and XJDEBCOU must be = to
                                          ! -999.)
                    NLANGLE,            & ! Angle between X Meso-NH axis and 
                                          ! cross-section direction in degrees
                                          ! (Integer value anticlockwise)
                    NLMAX,              & ! Number of points horizontally along
                                          ! the vertical section
                    NIFLAG         

REAL,SAVE        :: XIDEBCOU,XJDEBCOU,  & ! Origin of a vertical cross-section
                                          ! in cartesian (or conformal) real 
                                          ! values
                    XHMIN,              & ! altitude of the vert. cross-section
                                          ! bottom (in meters above sea-level)
                    XHMAX,              & ! altitude of the vert. cross-section
                                          ! top (in meters above sea-level)
                    XDZTRA                ! Not yet used

REAL,DIMENSION(3):: XEYE                  ! Not yet used

!
!*     0.1  Namelist NAM_DOMAIN_POS
!
NAMELIST/NAM_DOMAIN_POS/LHORIZ,NIINF,NISUP,NJINF,NJSUP,LVERTI,NIDEBCOU,NJDEBCOU,  &
XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,XHMIN,XHMAX,XDZTRA,L3D,XEYE,NIFLAG
!
END MODULE MODN_PARA
