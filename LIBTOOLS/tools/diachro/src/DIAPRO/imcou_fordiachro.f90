!     ######spl
      SUBROUTINE IMCOU_FORDIACHRO(PTABV,PINT,HLEGEND,HTEXT)
!     #####################################################
!
!!****  *IMCOU_FORDIACHRO* - Contour plot manager for vertical cross-sections
!!
!!    PURPOSE
!!    -------
!       Draws contour plots in the vertical cross-section case
!
!!**  METHOD
!!    ------
!!       Calls the NCAR contour routines and defines the display environment
!!      for the vertical cross-sections
!!
!!    EXTERNAL
!!    --------
!!      GMNMX    computes min, max and contour increment for current field
!!      TRACEXZ  draws a model-level stencil background if i) the current
!!               plot is a East-West cross-section, ii) the section origin
!!               is directly defined by grid indexes, and iii) if LXZ = .TRUE.
!!
!!      CURVE    draws a curve made by a series of data points      !
!!      SFSETR   sets parameters for NCAR softfill environment      !
!!      SFWRLD   fills the inside of a closed curve as requested by !
!!               the previous SFSETR calls                          !
!!                                                                  !
!!      CPSETI !                                          INTEGER   !
!!      CPSETR ! gives a value to a NCAR variabe, type:   REEL      !
!!      CPSETC !                                          CHARACTER !
!!      CPGETI !                                          INTEGER   !Routines
!!      CPGETI !                                          INTEGER   !
!!      CPGETR ! retrieves a NCAR parmeter value, type    REEL      !
!!      CPGETC !                                          CHARACTER !
!!      CPRECT   initialize contour drawing                         !
!!      CPPKCL   selects the contour values                         !
!!      CPCLDR   draws the contours                                 !
!!      CPLBDR   activates High and Low option                      !
!!      CPRSET   restores NCAR default values                       !
!!                                                                  !
!!      GSLWSC   sets line widths                                   !
!!      SET      defines the display window limits in both          !
!!               normalised and user coordinates                    !
!!      GETSET   retrieves the user and normalized coordinate ranges!
!!               for current window  for the current display window.! 
!!      PLCHHQ   prints high qualty text                            !
!!      GSCLIP   CLIPS the display window                           !
!!
!!      CPMPXY   TRACE provided FORTRAN-77 routine directly called  
!!               within CONPACK to map the array space onto the    
!!               Gal-Chen stretched  space
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist
!!                          (former PARA common)
!!          XIDEBCOU, XJDEBCOU : Origin of a vertical cross-section
!!                               in cartesian (or conformal) real values
!!          XHMIN      : Altitude of the vert. cross-section
!!                       bottom (in meters above sea-level)
!!          XHMAX      : Altitude of the vert. cross-section
!!                       top (in meters above sea-level)
!!          LHORIZ     : Horizontal mode selector
!!                       =.TRUE. to perform horizontal cross-sections
!!                       (LVERTI must be = to .FALSE.)
!!          LVERTI     : Vertical mode selector 
!!                       =.TRUE. to perform vertical cross-sections, 
!!                       including vert. 1D profiles. 
!!                       (LHORIZ must be = to .FALSE.)
!!          NIDEBCOU,  : Origin of a vertical cross-section
!!          NJDEBCOU     in grid index integer values
!!                       (XIDEBCOU and XJDEBCOU must be = to -999.)
!!          NLANGLE    : Angle between X Meso-NH axis and
!!                       cross-section direction in degrees
!!                       (Integer value anticlockwise)
!!          NLMAX      : Number of points horizontally along
!!                       the vertical section
!!          Module MODD_DIM1 : contains dimensions of data arrays
!!              NKMAX       : z array dimension
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!          JPHEXT : Horizontal external points number
!!          JPVEXT : Vertical external points number
!!
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID    : Current MESO-NH grid indicator
!!
!!     Module MODD_CVERT:  Declares work arrays for vertical cross-sections
!!          XWORKZ   : working array for true altitude storage (all grids)
!!          XWZ      : working array for topography (all grids)
!!
!!      Module MODD_COORD  : declares gridpoint coordinates
!!                           (TRACE use only)
!!          XDS      : Abscissa array along the horizontal axis of an oblique
!!                     vertical cross-section (meters), for all grid locations
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist
!!                         (former NCAR common)
!!          NIOFFD     : Label normalisation (=0 none, =/=0 active)
!!          NULBLL     : Nb of contours between 2 labelled contours
!!          NIOFFM     : =0    --> message at picture bottom
!!                       =/= 0 --> no message
!!          NDOT       : Line style
!!                        (=0|1|1023|65535 --> solid lines;
!!                        <0 --> solid lines for positive values and
!!                        dotted lines(ABS(NDOT))for negative values;
!!                        >0 --> dotted lines(ABS(NDOT)) )
!!          NHI        : Extrema detection
!!                       (=0 --> H+L, <0 nothing)
!!          NIMNMX     : Contour selection option
!!                       (=-1 Min, max and inc. automatically set;
!!                       =0 Min, max automatically set; inc. given;
!!                       >0 Min, max, inc. given by user)
!!          XSPVAL     : Special value
!!          XSIZEL     : Label size
!!
!!      Module MODD_SUPER   : defines plot overlay control variables
!!         LSUPER   : =.TRUE. --> plot overlay is active
!!                    =.FALSE. --> plot overlay is not active
!!         NSUPER   : Rank of the current plot in the overlay
!!                    sequence. The initial plot is rank 1.
!!
!!      Module MODD_ALLVAR
!!         >>>>>>>>>>DRAGOON QUERY: Is this one really necessary????
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
!!
!!    AUTHOR
!!    ------
!!
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   19/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef NAGf95
USE F90_UNIX  ! for FLUSH and GETENV
#endif

USE MODN_PARA
USE MODD_PARAMETERS
USE MODD_NMGRID
USE MODD_CVERT
USE MODD_COORD
USE MODD_CONF
USE MODD_GRID
USE MODD_GRID1
USE MODD_DIM1
USE MODD_TYPE_AND_LH
USE MODN_NCAR
USE MODD_SUPER
USE MODD_ALLVAR
USE MODD_TITLE
USE MODD_LUNIT1
USE MODD_OUT
USE MODD_PVT
USE MODD_RSISOCOL
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO
USE MODI_RESOLV_TIT
USE MODI_RESOLV_TITY
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODI_READMNMXINT_ISO
USE MODI_READREFINT_ISO
USE MODI_READXISOLEVP
USE MODD_TIT
USE MODD_HACH
USE MODD_DEFCV
USE MODE_GRIDPROJ
USE MODD_CTL_AXES_AND_STYL
USE MODD_MASK3D
!
USE MODI_CREATLINK
USE MODI_WRITEDIR
!      
IMPLICIT NONE
INTERFACE
SUBROUTINE AXELOGPRES(PHMIN,PHMAX)
REAL :: PHMIN,PHMAX
END SUBROUTINE AXELOGPRES
END INTERFACE
!
!*        0.0   TRACE interface with the "CPMPXY" routine of the NCAR package
!
! NOTICE:  The CPMPXY and the NCAR graphical utilities are NOT written
! ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
!          does not follow the Meso-NH usual rules: it has to be made using
!          a COMMON stack with  static memory allocation of XZZXX and
!          XZZXY arrays.
!
COMMON/TEMV/XZWORKZ,XZZDS,NINX,NINY
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
COMMON/COLAREA/ICOL(300)
COMMON/HACHAREA/IHACH(300)
#include "big.h"
REAL,DIMENSION(N2DVERTX,2500) :: XZWORKZ
!REAL,DIMENSION(1000,400) :: XZWORKZ
!REAL,DIMENSION(200,200) :: XZWORKZ
REAL,DIMENSION(N2DVERTX)     :: XZZDS
!REAL,DIMENSION(1000)     :: XZZDS
!REAL,DIMENSION(200)     :: XZZDS
INTEGER                 :: NINX, NINY
LOGICAL                 :: LVERT, LHOR, LPT, LXABS
INTEGER  :: ICOL
!
!*       0.1    Work arrays for NCAR
!
INTEGER,PARAMETER       :: JPLRWK=50000, JPLIWK=50000
INTEGER,PARAMETER       :: JPRSCR=20000, JPISCR=20000
INTEGER,PARAMETER       :: JPMAP=NPMAP, JPAREAGRP=300, JPWRK=50000
!INTEGER,PARAMETER       :: JPMAP=800000, JPAREAGRP=300, JPWRK=50000
!
REAL,DIMENSION(JPLRWK)      :: ZRWRK
INTEGER,DIMENSION(JPLIWK)      :: IWRK
REAL,DIMENSION(JPRSCR)      :: ZRSCR
INTEGER,DIMENSION(JPISCR)      :: ISCR
INTEGER,DIMENSION(JPMAP)    :: IIMAP
INTEGER,DIMENSION(JPAREAGRP):: IAREA, IGRP
REAL,DIMENSION(JPWRK)       :: ZXWRK, ZYWRK
INTEGER                     :: IHACH
!
!*       0.2   Dummy arguments and results
!
REAL,DIMENSION(:,:) :: PTABV            !  Vertical section data array 
                                        !  to be plotted
REAL                :: PINT             !  Contour increment fo the 
                                        !  current plot
CHARACTER(LEN=*)    :: HTEXT            !  PLot heading with section location 
CHARACTER(LEN=*)    :: HLEGEND          !  PLot heading with variable name
!CHARACTER(LEN=8) :: YDAT8, YTIM8, YTEM8
CHARACTER(LEN=32):: YLBL
CHARACTER(LEN=80)               :: YCAR80 
CHARACTER(LEN=160)               :: YCAR160,YCAR161
!
!*       0.3   Local variables
!
INTEGER :: IA, IB
INTEGER :: IKU, IKB, IKE, JILOOP, JKLOOP, J, JU
INTEGER :: ICL, INCL2, ILMAX
INTEGER :: INCL, I, ICLD, III, IO
INTEGER :: INBC, IDX, INBCT
INTEGER :: JJD, JJF, JI, JJ
INTEGER :: JB, ISTOK
INTEGER,SAVE :: ILUCOL, IRESP, ID, IDD
INTEGER,SAVE :: ISUIT, ISUI, INDISTM
INTEGER :: ILENT, IND, II2,IJ2
INTEGER             :: JLBL, JL
INTEGER             :: ISTA, IER, IWK, INB, INBB
INTEGER,SAVE        :: IH, IHT, IMI, ILE
INTEGER,DIMENSION(32):: INDHACHREF=(/0,54,52,60,14,59,58,1,57,56,55,54,53,52,51,50, &
			1,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35/)
INTEGER :: INUM, ILOOP, JLOOPI, IDEB,IFIN, II, JLOOPJ
INTEGER,DIMENSION(:),ALLOCATABLE       :: ICOL2
INTEGER,DIMENSION(:),ALLOCATABLE       :: IE
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE  :: ISTM
#ifdef RHODES
INTEGER          :: ISTAF
#endif

REAL    :: ZWLC, ZWRC, ZWBC, ZWTC
REAL    :: ZTA, ZTB, ZTD, ZTF, ZTINT,ZINTV
REAL    :: ZINT, ZMIN, ZMAX
REAL    :: ZINTT, ZH, ZJ, ZJJ, ZWBBB
REAL    :: ZISO
REAL    :: ZTEMP
REAL,SAVE :: ZWL, ZWR, ZWB, ZWT
REAL,SAVE :: ZWLL, ZWRR, ZWBB, ZWTT
REAL,SAVE :: ZVL, ZVR, ZVB, ZVT
REAL    :: ZCLV, ZINTERV, ZCLV2
REAL    :: ZCLVD, ZCLVF
REAL    :: RED, GREEN, BLUE
REAL    :: ZMN, ZMX
REAL    :: ZDIXEPS 
REAL    :: ZX, ZY, ZXE, ZYE
REAL    :: ZLAT, ZLON
REAL    :: ZMI, ZMA, ZMIG, ZMAG
REAL    :: ZVLDEF, ZWIDTH
REAL    :: ZSC
REAL    :: ZXPOSTITT1, ZXYPOSTITT1
REAL    :: ZXPOSTITT2, ZXYPOSTITT2
REAL    :: ZXPOSTITT3, ZXYPOSTITT3
REAL    :: ZXPOSTITB1, ZXYPOSTITB1
REAL    :: ZXPOSTITB2, ZXYPOSTITB2
REAL    :: ZXPOSTITB3, ZXYPOSTITB3
REAL,DIMENSION(5)   :: ZX5, ZY5
REAL                :: ZEPX, ZEPYD, ZEPYU
!
REAL,SAVE           :: ZD, ZF, ZVERA, ZINTE
REAL,DIMENSION(SIZE(PTABV,1),SIZE(PTABV,2)):: ZTEMV, ZTEMV2
REAL,DIMENSION(N2DVERTX+20)                        :: ZDS, ZWZ
!REAL,DIMENSION(1020)                        :: ZDS, ZWZ
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZDS2, ZWZ2
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZLA, ZLO
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZZCLV2, ZTDX

!REAL,DIMENSION(300)                        :: ZDS, ZWZ
REAL,DIMENSION(300)                        :: ZLEV, ZISOLEVP
!
CHARACTER(LEN=5)   :: YFORMAT
CHARACTER(LEN=82),SAVE :: YCARCOU, YCAR
CHARACTER(LEN=100) :: YTEM
CHARACTER(LEN=1)   :: YREP
CHARACTER(LEN=2)   :: YC2
CHARACTER(LEN=3)   :: YC3  
CHARACTER(LEN=4)   :: YC4  
CHARACTER(LEN=8),DIMENSION(300) :: YLLBS
CHARACTER(LEN=32),SAVE          :: YNAMTABCOL
CHARACTER(LEN=40)  :: YTEXT
CHARACTER(LEN=45)  :: YTEX  ! 45=40+5
CHARACTER(LEN=8)   :: YC8  
CHARACTER(LEN=20)  :: YXYO 
CHARACTER(LEN=20)  :: YCAR20
CHARACTER(LEN=10) :: FORMAX, FORMAY,FORMA160
!
EXTERNAL SFILL     
EXTERNAL SFILLH     
EXTERNAL CCOLR
!
!-----------------------------------------------------------------------------
!
!*       1.     DISPLAY ENVIRONMENT SETUP
!               -------------------------
!
!-----------------------------------------------------------------------------
if(nverbia > 0)then
  print *,' ENTREE IMCOU'
  print *,'  LEN_TRIM(HTEXT) ',LEN_TRIM(HTEXT),HTEXT(1:LEN_TRIM(HTEXT))
  print *,'  LPRESY,XHMIN,XHMAX CTIMEC ',LPRESY,XHMIN,XHMAX,CTIMEC
  print *,'  CLEGEND2 ',CLEGEND2
endif
ZVLDEF=.1
YTEXT(1:LEN(YTEXT))=' '
YTEX(1:LEN(YTEX))=' '
!HTEXT=ADJUSTL(HTEXT)
JU=0
DO J=1,LEN_TRIM(HTEXT)
  IF(HTEXT(J:J) == ' ')THEN
    JU=JU+1
    YTEXT(1:J-1)=HTEXT(1:J-1)
    IF(YTEXT(1:4) == 'MASK')THEN
      IF(JU == 2)THEN
        IF(YTEXT(1:4) == 'MASK')THEN
          IF(YTEXT(6:6) /= ' ')THEN
            YTEXT(1:6)=' '
          ELSE
            YTEXT(1:5)=' '
          ENDIF
          YTEXT=ADJUSTL(YTEXT)
          EXIT
        ENDIF
      ENDIF
    ELSE
    EXIT
    ENDIF
  ENDIF
  IF(J == LEN_TRIM(HTEXT))THEN
    YTEXT=HTEXT
    YTEXT=ADJUSTL(YTEXT)
    IF(YTEXT(1:4) == 'MASK')THEN
      IF(YTEXT(6:6) /= ' ')THEN
        YTEXT(1:6)=' '
      ELSE
        YTEXT(1:5)=' '
      ENDIF
      YTEXT=ADJUSTL(YTEXT)
    ENDIF
  ENDIF
ENDDO

IF(nverbia > 0)then
  print *,' IMCOU NMGRID YTEXT ',NMGRID,YTEXT
  print *,' PTABV',size(PTABV,1),size(PTABV,2),PTABV(1,1),PTABV(size(PTABV,1),6)
endif
NLUOUT=6

IF(LPRINT)THEN
  
! IF(LDEFCV2CC)THEN                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! print *,' Pour l''instant, cette operation n''est prevue que pour une coupe definie avec :'
! print *,' NIDEBCOU= NJDEBCOU= NLANGLE= NLMAX= '
! print *,' A suivre ........ '
! ELSE                                 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=SIZE(PTABV,1)/5
  IF(ILOOP * 5 < SIZE(PTABV,1))ILOOP=ILOOP+1
  IF(.NOT.LPVT)THEN
    WRITE(INUM,'(''CV  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'',''   (1-NLMAX,1-IKU)'')')CGROUP,&
&   CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
  ELSE
    WRITE(INUM,'(''CV  '',''G:'',A16,'' P:'',A25)')CGROUP,&
&   CTITRE(NLOOPP)(1:25)
  ENDIF
  IF(LMINUS .OR. LPLUS)THEN
    WRITE(INUM,'(A70)')CTITB3
  ELSE
    WRITE(INUM,'(A40)')CTITGAL
  ENDIF
  IF(.NOT.LPVT)THEN
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2)THEN
        WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,SIZE(PTABV,2),ILOOP
      ELSE IF(LDEFCV2LL)THEN
        WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,SIZE(PTABV,2),ILOOP
      ELSE IF(LDEFCV2IND)THEN
        WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,SIZE(PTABV,2),ILOOP
      ENDIF
    ELSE
      IF(XIDEBCOU /= -999.)THEN
        WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,SIZE(PTABV,2),ILOOP
      ELSE
        WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,SIZE(PTABV,2),ILOOP
      ENDIF
    ENDIF
  ELSE
    WRITE(INUM,'(''NBVAL en I (TIME): '',i4, &
&  '' NBVAL en K (Z)'',i4,''    iter'',i3)') &
  & SIZE(PTABV,1),SIZE(PTABV,2),ILOOP
  ENDIF
  DO JLOOPI=1,ILOOP
    IF(JLOOPI == 1)THEN
      IDEB=1; IFIN=5
    ELSE
      IDEB=IFIN+1; IFIN=IFIN+5
    ENDIF
    IF(JLOOPI == ILOOP)THEN
      IFIN=SIZE(PTABV,1)
    ENDIF
    
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'** IMCOU XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
    ELSE
      WRITE(INUM,'(1X,75(1H*))')
      WRITE(INUM,'(1X,''    Dates courante   *     modele      *   experience    *      segment'')')
      WRITE(INUM,'(1X,'' J   An  M  J  Sec.  * An  M  J  Sec.  * An  M  J  Sec.  * An  M  J  Sec.'')')
      WRITE(INUM,'(1X,75(1H*))')
      DO J=1,SIZE(XPRDAT,2)
        WRITE(INUM,'(1X,I3,1X,3(I4,I3,I3,I6,'' *''),I4,I3,I3,I6)')J,INT(XPRDAT(:,J))
      ENDDO
    ENDIF
  ENDIF
    WRITE(INUM,'(1X,79(1H*))')
    WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
    WRITE(INUM,'(''.'',79(1H*))')
    DO JLOOPJ=SIZE(PTABV,2),1,-1
      WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(PTABV(II,JLOOPJ),II=IDEB,IFIN)
!     WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(PTABV(II,JLOOPJ),II=IDEB,IFIN)
    ENDDO
    WRITE(INUM,'(1X,79(1H*))')
  ENDDO
! ENDIF                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENDIF


ZDIXEPS=1.E-11 
!print *,' ZDIXEPS ',ZDIXEPS
IKU=NKMAX+2*JPVEXT
IKB=1+JPVEXT
IKE=IKU-JPVEXT
LVERTI=.TRUE.; LHORIZ=.FALSE.
!IF(.NOT.LCVXZ .AND. .NOT.LCVYZ)THEN
!
!*       1.1    Window definition, NDC and user coordinate setting
!
XLWIDTH=XLWDEF
IF(LSUPER)THEN
  NSUPER=NSUPER+1
  SELECT CASE(NSUPER)
    CASE(1)
      IF(XLW >= 0)THEN
	XLWIDTH=XLW
      ENDIF
      IF(XLW1 >= 0)THEN
	XLWIDTH=XLW1
      ENDIF

      IH=0; IHT=0

      IF(LHACH2 .AND. LHACH3 .AND. LHACH4)THEN
     
        IHT=3
      ELSE IF((LHACH2 .AND. LHACH3 .AND. .NOT.LHACH4) .OR.  &
              (LHACH2 .AND. LHACH4 .AND. .NOT.LHACH3) .OR.  &
	      (LHACH3 .AND. LHACH4 .AND. .NOT.LHACH2))THEN
	      IHT=2
      ELSE IF((LHACH2 .AND. .NOT.LHACH3 .AND. .NOT.LHACH4) .OR.  &
 	      (LHACH3 .AND. .NOT.LHACH2 .AND. .NOT.LHACH4) .OR.  &
              (LHACH4 .AND. .NOT.LHACH2 .AND. .NOT.LHACH3))THEN
	      IHT=1
      ENDIF

    CASE(2)
      IF(XLW2 >= 0)THEN
	XLWIDTH=XLW2
      ENDIF
    CASE(3)
      IF(XLW3 >= 0)THEN
	XLWIDTH=XLW3
      ENDIF
    CASE(4)
      IF(XLW4 >= 0)THEN
	XLWIDTH=XLW4
      ENDIF
  END SELECT
ELSE
  IF(XLW >= 0)THEN
    XLWIDTH=XLW
  ENDIF
  IF(XLW1 >= 0)THEN
    XLWIDTH=XLW1
  ENDIF
  IH=0; IHT=0
END IF

LPT=LPXT
IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT .AND. .NOT.LCVXZ .AND. .NOT.LCVYZ)THEN

  IF((.NOT.LSUPER) .OR. (LSUPER .AND. NSUPER == 1))THEN
    ZWL=XDS(1,NMGRID)
    ZWR=XDS(NLMAX,NMGRID)
! Nov 2000
    IF(LPRESY)THEN
      IF(XHMIN<=XHMAX)THEN
! Bornes en altitude -> besoin de calculer bornes en pression +loin
        !XHMIN=0.
        !XHMAX=XWORKZ(1,IKE,NMGRID)
        IF(XPMIN==XPMAX)THEN
          print*,' ordonnee en Log(P): indiquez XPMIN et XPMAX'
          read(5,*) XPMIN,XPMAX
          CALL WRITEDIR(NDIR,XPMIN)
          CALL WRITEDIR(NDIR,XPMAX)
        ENDIF
        IF(XPMIN<XPMAX) THEN
          ZTEMP=XPMIN
          XPMIN=XPMAX
          XPMAX=ZTEMP
        ENDIF
        XHMIN=XPMIN
        XHMAX=XPMAX
      ENDIF
! Bornes fournies en pression . Verifier qu'elles sont en pascals
! Besoin de calculer les bornes en altitudes +loin
      IF (XHMIN < 1500)THEN
        XHMIN=XHMIN*100
      ENDIF
      IF (XHMAX < 1500)THEN
        XHMAX=XHMAX*100
      ENDIF
    ELSE
      IF((XHMIN==0..AND.XHMAX==0.).OR.(XHMAX<=XHMIN))THEN
! Nov 2000 -> Petite modif a signaler aux utilisateurs
        XHMIN=0.
!       XHMIN=XWORKZ(1,IKB,NMGRID)
        XHMAX=XWORKZ(1,IKE,NMGRID)
      ENDIF
    ENDIF
    ZWB=XHMIN
    ZWT=XHMAX
    IF (.NOT. LPRESY .AND. ZWB==ZWT) THEN
      print *,' min, max identiques pour la 2e direction: ',XHMIN,XHMAX
      print *,'entrez 2 valeurs telles que XHMIN < XHMAX '
      read(5,*) ZWB,ZWT
      CALL WRITEDIR(NDIR,ZWB)
      CALL WRITEDIR(NDIR,ZWT)
    END IF
!
    if(nverbia > 0)then
      print *,' ****** IMCOU_FORDIACHRO ZWL R B T',ZWL,ZWR,ZWB,ZWT
    endif
    LVERT=LVERTI
    LHOR=LHORIZ
!
! Nov 2000
    IF(LPRESY)THEN
      CALL SETUSV('MI',1)
      CALL SETUSV('LS',2)
      IF(LVPTVUSER)THEN
        CALL SET(XVPTVL,XVPTVR,XVPTVB,XVPTVT,ZWL,ZWR,ZWB,ZWT,2)
      ELSE
        CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,2)
      ENDIF
    ELSE
! Nov 2000
      CALL SETUSV('MI',1)
      IF(LVPTVUSER)THEN
        CALL SET(XVPTVL,XVPTVR,XVPTVB,XVPTVT,ZWL,ZWR,ZWB,ZWT,1)
      ELSE
        CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
      ENDIF
! Nov 2000
    ENDIF
! Nov 2000
  END IF

ELSE

  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
    
  IF(LCVXZ .OR. LCVYZ)THEN
    IF(LVPTVUSER)THEN
      CALL SET(XVPTVL,XVPTVR,XVPTVB,XVPTVT,ZWL,ZWR,ZWB,ZWT,1)
    ELSE
! Dans ce cas definition de la fenetre ds OPER avec .1,.9,.1,.9
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
    ENDIF
  ELSE 
!!!!!PROVI
  IF(LPXT .AND. .NOT.LXABSC .AND. LXMINTOP)THEN
    CALL SETUSV('MI',2)
! Attention ici inversion de ZWB et ZWT
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWT,ZWB,ID)
  ELSEIF(LPVT .AND. LPRESY)THEN
    CALL SETUSV('MI',1)
    CALL SETUSV('LS',2)
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,2)
!   CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  ELSE
    CALL SETUSV('MI',1)
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  ENDIF
  ENDIF

ENDIF
CALL GETUSV('MI',IMI)
!
IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT .AND. .NOT.(LCVXZ.AND.LJCP) .AND. .NOT.(LCVYZ.AND.LICP))THEN
  CALL GSCLIP(1)              ! Display clipping activated
!
  CALL CPSETI('SET',0)  ! Compack keeps user's call to set
  CALL CPSETI('MAP',4)  ! Customized vertical z-stretching used in CPMPXY
!
!*      1.2    Topography outline drawing
! 
  
  ZDS(1)=XDS(1,NMGRID)
  ZWZ(1)=XHMIN
  IF(LCVYZ .AND.LICP)THEN
    ZWZ(1)=0.
  ENDIF
  IF(LCVYZ .AND. .NOT.LICP)THEN
    ZWZ(2:NLMAX+1)=XXZS(NIDEBCOU,NJDEBCOU:NJDEBCOU+NLMAX-1,NMGRID)
  ENDIF
  DO JILOOP=2,NLMAX+1
    ZDS(JILOOP)=XDS(JILOOP-1,NMGRID)
    IF(LCVYZ .AND. .NOT.LICP)THEN
    ELSEIF(LCVYZ .AND.LICP)THEN
      ZWZ(JILOOP)=0.
    ELSE
      ZWZ(JILOOP)=XWZ(JILOOP-1,NMGRID)
    ENDIF
  ENDDO
  ZDS(NLMAX+2)=ZDS(NLMAX+1)
  ZWZ(NLMAX+2)=XHMIN
  IF(LCVYZ .AND.LICP)THEN
    ZWZ(NLMAX+2)=0.
  ENDIF
!
  IF(ALLOCATED(ZDS2))THEN
   DEALLOCATE(ZDS2)
  ENDIF
  IF(ALLOCATED(ZWZ2))THEN
   DEALLOCATE(ZWZ2)
  ENDIF
  ALLOCATE(ZDS2(NLMAX+2))
  ALLOCATE(ZWZ2(NLMAX+2))
  ZDS2=ZDS(1:NLMAX+2)
  ZWZ2=ZWZ(1:NLMAX+2)
  if(nverbia > 4)then
print *,' ********IMCOU_FORDIACHRO NLMAX  ZDS',NLMAX
print *,(ZDS(JILOOP),JILOOP=1,NLMAX)
print *,' ********IMCOU_FORDIACHRO ZWZ'
print *,(ZWZ(JILOOP),JILOOP=1,NLMAX)
  endif
!
  IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN
  IF(.NOT. LPRESY) THEN
    CALL CURVE(ZDS2,ZWZ2,NLMAX+2)                          ! draws Topo outline 
!   CALL CURVE(ZDS,ZWZ,NLMAX+2)                            ! draws Topo outline 
    CALL SFSETR('SP',.008)                                 ! Softfill setting
    CALL SFSETR('AN',45.)                                  ! Softfill setting
    CALL SFSETI('DO',0)                                  ! Softfill setting
    CALL SFWRLD(ZDS2,ZWZ2,NLMAX+2,ZRSCR,JPRSCR,ISCR,JPISCR)  ! Hatched under 
!   CALL SFWRLD(ZDS,ZWZ,NLMAX+2,ZRSCR,JPRSCR,ISCR,JPISCR)  ! Hatched under 
!                                                      ! topography
  ENDIF
!
!*     1.3     If required, draws a model-level background
!
    IF(.NOT.LDEFCV2CC)THEN              !%%%%%%%%%%%%%%%%%%%%%%

    IF(NLANGLE.EQ.0.AND.XIDEBCOU.EQ.-999..AND.LXZ)THEN
      CALL GSCLIP(0)
      CALL TRACEXZ
      CALL GSCLIP(1)
    END IF

    ENDIF                               !%%%%%%%%%%%%%%%%%%%%%%

  ENDIF

ENDIF
!
!-----------------------------------------------------------------------------
!
!*     2.        CONTOUR DRAWING 
!                ---------------
!
!*     2.1       Loads abscissa and true-altitudes along
!*               the section in work arrays 

IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT .AND. .NOT.LCVXZ .AND. .NOT.LCVYZ)THEN

  NINX=NLMAX
  NINY=IKU
  DO JILOOP=1,NLMAX
    XZZDS(JILOOP)=XDS(JILOOP,NMGRID)
  ENDDO
!print *,' ********IMCOU_FORDIACHRO NLMAX  XZZDS',NLMAX
!print *,(XZZDS(JILOOP),JILOOP=1,NLMAX)
  DO JILOOP=1,NLMAX
    DO JKLOOP=1,IKU
      XZWORKZ(JILOOP,JKLOOP)=XWORKZ(JILOOP,JKLOOP,NMGRID)
    ENDDO
  ENDDO

ENDIF
!-----------------------------------------------------------------------------
IF(LPRINTXY)THEN
! IF(LDEFCV2CC .OR. XIDEBCOU /= -999.)THEN   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! print *,' Pour l''instant, cette operation n''est prevue que pour une coupe definie avec :'
! print *,' NIDEBCOU= NJDEBCOU= NLANGLE= NLMAX= '
! print *,' A suivre ........ '
! ELSE                                 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=SIZE(PTABV,1)/5
  IF(ILOOP * 5 < SIZE(PTABV,1))ILOOP=ILOOP+1
  IF(.NOT. LPVT)THEN
!!Oct 2002
    IF(LCVYZ)THEN
    WRITE(INUM,'(''CV YZ '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'','' (1-NLMAX,1-IKU)'')')CGROUP, &
&   CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ELSE
!!Oct 2002
    WRITE(INUM,'(''CV XZ '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'','' (1-NLMAX,1-IKU)'')')CGROUP, &
&   CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ENDIF
  ELSE
    WRITE(INUM,'(''CV TIMEZ '',''   G:'',A16,'' P:'',A40)')CGROUP, &
!&   CTITGAL
&   CTITRE(NLOOPP)(1:40)
  ENDIF
  IF(LMINUS .OR. LPLUS)THEN
    WRITE(INUM,'(A70)')CTITB3
  ELSE
    WRITE(INUM,'(A40)')CTITGAL
  ENDIF
  IF(.NOT. LPVT)THEN
    IF(.NOT.LCARTESIAN)THEN
      ALLOCATE(ZLA(NLMAX),ZLO(NLMAX))
      DO J=1,NLMAX
	ZX=XDSX(J,NMGRID)
	ZY=XDSY(J,NMGRID)
	CALL SM_LATLON_S(XLATORI,XLONORI,ZX,ZY,ZLAT,ZLON)
	ZLA(J)=ZLAT
	ZLO(J)=ZLON
      ENDDO
      IF(LDEFCV2LL)THEN
	ZLA(1)=XIDEBCVLL
	ZLO(1)=XJDEBCVLL
      ENDIF
      if(nverbia > 0)then
!     print *,' ZLA'
!     print *,ZLA
!     print *,' ZLO'
!     print *,ZLO
      endif
!     DEALLOCATE(ZLA,ZLO)
    ENDIF
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2)THEN
        WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,SIZE(PTABV,2),ILOOP
      ELSE IF(LDEFCV2LL)THEN
        WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,SIZE(PTABV,2),ILOOP
      ELSE IF(LDEFCV2IND)THEN
        WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,SIZE(PTABV,2),ILOOP
      ENDIF
    ELSE
      IF(XIDEBCOU /= -999.)THEN
        WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,SIZE(PTABV,2),ILOOP
      ELSE
    WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4, &
&  '' iku'',i4,''    iter'',i3)') &
  & NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,SIZE(PTABV,2),ILOOP
      ENDIF
    ENDIF
    IF(LCARTESIAN)THEN
      WRITE(INUM,'(1X,41(1H*))')
      WRITE(INUM,'(18X,''X'',12X,''RELIEF'')')
      WRITE(INUM,'(1X,41(1H*))')
      DO JLOOPI=1,NLMAX
        IF(JLOOPI == 1)THEN
          WRITE(INUM,'(''   1 '',I5,2(1X,E15.8))')JLOOPI,XDS(JLOOPI,NMGRID), &
          XWZ(JLOOPI,NMGRID)
        ELSE IF(JLOOPI == NLMAX)THEN
          WRITE(INUM,'(''NLMAX'',I5,2(1X,E15.8))')JLOOPI,XDS(JLOOPI,NMGRID), &
          XWZ(JLOOPI,NMGRID)
        ELSE
          WRITE(INUM,'(''     '',I5,2(1X,E15.8))')JLOOPI,XDS(JLOOPI,NMGRID), &
          XWZ(JLOOPI,NMGRID)
        ENDIF
      ENDDO
      WRITE(INUM,'(1X,41(1H*))')
    ELSE
      WRITE(INUM,'(1X,66(1H*))')
      WRITE(INUM,'(18X,''X'',12X,''RELIEF'',11X,''LAT'',10X,''LONG'')')
      WRITE(INUM,'(1X,66(1H*))')
      DO JLOOPI=1,NLMAX
        IF(JLOOPI == 1)THEN
          IF(LCVYZ)THEN
            WRITE(INUM,'(''   1 '',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
            ZWZ(JLOOPI+1),ZLA(JLOOPI),ZLO(JLOOPI)
          ELSE
            WRITE(INUM,'(''   1 '',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
            XWZ(JLOOPI,NMGRID),ZLA(JLOOPI),ZLO(JLOOPI)
          END IF
        ELSE IF(JLOOPI == NLMAX)THEN
          IF(LCVYZ)THEN
            WRITE(INUM,'(''NLMAX'',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
            ZWZ(JLOOPI+1),ZLA(JLOOPI),ZLO(JLOOPI)
          ELSE
            WRITE(INUM,'(''NLMAX'',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
            XWZ(JLOOPI,NMGRID),ZLA(JLOOPI),ZLO(JLOOPI)
          END IF
        ELSE
          IF(LCVYZ)THEN
            WRITE(INUM,'(''     '',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
            ZWZ(JLOOPI+1),ZLA(JLOOPI),ZLO(JLOOPI)
          ELSE
            WRITE(INUM,'(''     '',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
            XWZ(JLOOPI,NMGRID),ZLA(JLOOPI),ZLO(JLOOPI)
          END IF
        ENDIF
      ENDDO
      WRITE(INUM,'(1X,66(1H*))')
      DEALLOCATE(ZLA,ZLO)
    ENDIF
  
    DO JLOOPI=1,ILOOP
      IF(JLOOPI == 1)THEN
        IDEB=1; IFIN=5
      ELSE
        IDEB=IFIN+1; IFIN=IFIN+5
      ENDIF
      IF(JLOOPI == ILOOP)THEN
        IFIN=SIZE(PTABV,1)
      ENDIF
      
      WRITE(INUM,'(''ALTITUDES   (1-NLMAX,1-IKU)'')')
      WRITE(INUM,'(1X,79(1H*))')
      WRITE(INUM,'(''  K  X->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',79(1H*))')
      DO JLOOPJ=SIZE(PTABV,2),1,-1
        IF(LCVYZ)THEN
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(XZWORKZ(II,JLOOPJ),II=IDEB,IFIN)
        ELSE
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(XWORKZ(II,JLOOPJ,NMGRID),II=IDEB,IFIN)
!       WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(XWORKZ(II,JLOOPJ,NMGRID),II=IDEB,IFIN)
        ENDIF
      ENDDO
      WRITE(INUM,'(1X,79(1H*))')
    ENDDO

  ELSE

    WRITE(INUM,'(''NBVAL en I (TIME): '',i4, &
&  '' NBVAL en K (Z)'',i4)') &
  & SIZE(PTABV,1),SIZE(PTABV,2)
    ZMIG=MINVAL(XZWORKZ(1:NINX,1:NINY))
    ZMAG=MAXVAL(XZWORKZ(1:NINX,1:NINY))
    ZMI=MINVAL(XZWORKZ(NINX/2,1:NINY))
    ZMA=MAXVAL(XZWORKZ(NINX/2,1:NINY))
!   print *,' ZMIG,ZMAG,ZMI,ZMA ',ZMIG,ZMAG,ZMI,ZMA

    IF(ZMIG == ZMI .AND. ZMAG == ZMA)THEN

      II=MAX(SIZE(PTABV,1),SIZE(PTABV,2))
      WRITE(INUM,'(1X,43(1H*))')
      WRITE(INUM,'(2X,''  I'',7X,''TIME'',10X,''K'',9X,''Z'')')
      WRITE(INUM,'(1X,43(1H*))')
      DO JLOOPJ=1,II
        IF(SIZE(PTABV,1) > SIZE(PTABV,2))THEN
          IF(JLOOPJ <= SIZE(PTABV,2))THEN
             WRITE(INUM,'(I5,2X,E15.8,1X,I4,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ), &
            JLOOPJ,XZWORKZ(1,JLOOPJ)
          ELSE
            WRITE(INUM,'(I5,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ)
          ENDIF
        ELSE IF(SIZE(PTABV,2) > SIZE(PTABV,1))THEN
          IF(JLOOPJ <= SIZE(PTABV,1))THEN
            WRITE(INUM,'(I5,2X,E15.8,1X,I4,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ), &
            JLOOPJ,XZWORKZ(1,JLOOPJ)
          ELSE
            WRITE(INUM,'(23X,I4,2X,E15.8)')JLOOPJ,XZWORKZ(1,JLOOPJ)
          ENDIF
        ELSE
          WRITE(INUM,'(I5,2X,E15.8,1X,I4,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ), &
          JLOOPJ,XZWORKZ(1,JLOOPJ)
        ENDIF
      ENDDO
      WRITE(INUM,'(1X,43(1H*))')

    ELSE

      DO JLOOPI=1,ILOOP
        IF(JLOOPI == 1)THEN
          IDEB=1; IFIN=5
        ELSE
          IDEB=IFIN+1; IFIN=IFIN+5
        ENDIF
        IF(JLOOPI == ILOOP)THEN
          IFIN=SIZE(PTABV,1)
        ENDIF

        WRITE(INUM,'(''TEMPS - ALTITUDES '')')
        WRITE(INUM,'(1X,79(1H*))')
!       WRITE(INUM,'("  K  I->  ",I5,5X,4(5X,I5,5X))')
        ALLOCATE(IE(IFIN-IDEB+1))
        DO III=IDEB,IFIN
        IE(III-IDEB+1)=III
        ENDDO
        WRITE(INUM,'("  K  I->  ",I5,5X,4(5X,I5,5X))')IE
!       WRITE(INUM,'("  K  I->  ",I5,5X,4(5X,I5,5X))')(/(III,III=IDEB,IFIN)/)
        DEALLOCATE(IE)
        WRITE(INUM,'(1X,79(1H.))')
        WRITE(INUM,'("   . TIME->",F7.0,3X,4(4X,F7.0,4X))')(XZZDS(II),II=IDEB,IFIN)
!       WRITE(INUM,'("           ")')
!       WRITE(INUM,'(F7.0,3X,4(4X,F7.0,4X))')(XZZDS(II),II=IDEB,IFIN)
        WRITE(INUM,'(''.'',79(1H*))')
        DO JLOOPJ=SIZE(PTABV,2),1,-1
          WRITE(INUM,'(I4,2X,5(1X,E14.7))')JLOOPJ,(XZWORKZ(II,JLOOPJ),II=IDEB,IFIN)
!         WRITE(INUM,'(I3,2X,5E15.8)')JLOOPJ,(XZWORKZ(II,JLOOPJ),II=IDEB,IFIN)
        ENDDO
        WRITE(INUM,'(1X,79(1H*))')
      ENDDO
    ENDIF

  ENDIF
! ENDIF                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENDIF
!-----------------------------------------------------------------------------
!
!*    2.2       If required, the user provides Max and Min of the field
!*              to be plotted (within section)
! 
ZINT=PINT

IF(NIMNMX == 0 .OR. NIMNMX == 1)THEN

! Modifs for Diachro
!
!CALL GMNMX(ZMIN,ZMAX,ZINT)
  LISOK=.FALSE.
  ZMIN=0.; ZMAX=0.
  CALL READMNMXINT_ISO(NIMNMX,YTEXT(1:LEN_TRIM(YTEXT)),ZMIN,ZMAX,ZINT)

ELSE IF(NIMNMX == 2)THEN
  CALL READXISOLEVP(YTEXT(1:LEN_TRIM(YTEXT)),ILE,ZISOLEVP)
  IF(NVERBIA > 5)THEN
    print *,' IMCOU YTEXT,ILE,ZISOLEVP ',YTEXT(1:LEN_TRIM(YTEXT)),ILE,ZISOLEVP(1:ILE)
  ENDIF

ELSE IF (NIMNMX==3) THEN  ! compute contour values from XISOREF and XDIAINT
  ZISOLEVP(:)=9999.
  ZMN=MINVAL(PTABV,MASK=PTABV/=XSPVAL) 
  ZMX=MAXVAL(PTABV,MASK=PTABV/=XSPVAL)
  CALL READREFINT_ISO(YTEXT(1:LEN_TRIM(YTEXT)),ZMN,ZMX,ZINT,ZISOLEVP)
ENDIF

IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT)THEN
! min + max matrice
if(nverbia >0)then
print *,' ** imcou NLMAX ',NLMAX
endif
ZMN=PTABV(NLMAX/2,SIZE(PTABV,2)/2)
ZMX=PTABV(NLMAX/2,SIZE(PTABV,2)/2)
if(nverbia >0)then
print *,' ** imcou AP ZMN=PTABV(NLMAX/2,SIZE(PTABV,2)/2); ZM...'
endif
ELSE
II2=MAX(1,SIZE(PTABV,1)/2); IJ2=MAX(1,SIZE(PTABV,2)/2)
ZMN=PTABV(II2,IJ2); ZMX=ZMN
!ZMN=999999.; ZMX=-999999.
ENDIF
!-----------------------------------------------------------------------------
IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT)THEN
DO JILOOP=1,NLMAX
  DO JKLOOP=1,IKU
    IF(LPRESY)THEN
! en log(pression)
      IF(XZWORKZ(JILOOP,JKLOOP) < XHMAX)CYCLE
      IF(XZWORKZ(JILOOP,JKLOOP) > XHMIN)CYCLE
    ELSE
      IF(XZWORKZ(JILOOP,JKLOOP) > XHMAX)CYCLE
      IF(XZWORKZ(JILOOP,JKLOOP) < XHMIN)CYCLE
    ENDIF
    IF(PTABV(JILOOP,JKLOOP) == XSPVAL)CYCLE
    IF(PTABV(JILOOP,JKLOOP) < ZMN)ZMN=PTABV(JILOOP,JKLOOP)
    IF(PTABV(JILOOP,JKLOOP) > ZMX)ZMX=PTABV(JILOOP,JKLOOP)
  ENDDO
ENDDO
!-----------------------------------------------------------------------------
ELSE
IF(.NOT.LPXT .AND..NOT.LPYT)THEN
DO JILOOP=1,SIZE(PTABV,1)
  DO JKLOOP=1,SIZE(PTABV,2)
    IF(LPRESY)THEN
! en log(pression)
      IF(XZWORKZ(JILOOP,JKLOOP) < XHMAX)CYCLE
      IF(XZWORKZ(JILOOP,JKLOOP) > XHMIN)CYCLE
    ELSE
      IF(XZWORKZ(JILOOP,JKLOOP) > XHMAX)CYCLE
      IF(XZWORKZ(JILOOP,JKLOOP) < XHMIN)CYCLE
    ENDIF
    IF(PTABV(JILOOP,JKLOOP) == XSPVAL)CYCLE
    IF(PTABV(JILOOP,JKLOOP) < ZMN)ZMN=PTABV(JILOOP,JKLOOP)
    IF(PTABV(JILOOP,JKLOOP) > ZMX)ZMX=PTABV(JILOOP,JKLOOP)
  ENDDO
ENDDO
ELSE
  ZMN=MINVAL(PTABV)
  ZMX=MAXVAL(PTABV)
ENDIF
ENDIF
!-----------------------------------------------------------------------------
YLBL(1:5)='(Min:'
WRITE(YLBL(6:15),'(E10.3)')ZMN
YLBL(16:21)=', Max:'
WRITE(YLBL(22:31),'(E10.3)')ZMX
YLBL(32:32)=')'
!
!*    2.3       Conpack display options 
!
CALL GSLWSC(1.)             ! Line width
!
!
!*    2.4       Contour selection rules
!
!print *,' ** imcou AV SELECT CASE(NIMNMX) '
SELECT CASE(NIMNMX)
  CASE(-1)             ! Automatic contour scanning
    CALL CPSETI('CLS',+16)
    IF((LHACH1 .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))) .OR. &
       (LHACH2 .AND. NSUPER == 2)                                   .OR. &
       (LHACH3 .AND. NSUPER == 3)                                   .OR. &
       (LHACH4 .AND. NSUPER == 4))CALL CPSETI('CLS',+7)

    CALL CPSETR('CIS',-ZINT)
!
  CASE(0)               ! Automatic range with given increment
    CALL CPSETI('CLS',16)
    CALL CPSETR('CIS',ZINT)
    CALL CPSETI('LIS',NULBLL+1)
    CALL CPSETR('CMN',100000000000.)
!   CALL CPSETR('CMN',MAXVAL(PTAB))
    CALL CPSETR('CMX',10000000000.)
!   CALL CPSETR('CMX',MINVAL(PTAB))
!
  CASE(1)               ! Given min, max and increment
    IF(ZMAX == ZMIN)THEN
      ICL=1
      CALL CPSETI('NCL',ICL)
    ELSE
    ICL=NINT((ZMAX-ZMIN)/ZINT)
    IF(ZMIN + ICL*ZINT <= ZMAX)ICL=ICL+1
    CALL CPSETI('NCL',ICL)
!   IF(LCOLAREA .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))CALL CPSETI('NCL',ICL+1)
    ENDIF
    CALL CPSETI('CLS',0)
    ZISO=ZMIN-ZINT
    DO I=1,ICL
    CALL CPSETI('PAI',I)
    CALL CPSETI('AIA',I+1)
    CALL CPSETI('AIB',I)
    ZISO=ZISO+ZINT
    IF(ABS(ZISO)<1.E-20)ZISO=0.
    CALL CPSETR('CLV',ZISO)
    CALL CPSETR('CLU',1.)
    IF(.NOT.LSUPER.OR. (LSUPER .AND. NSUPER == 1))THEN
      IF(LBLUSER1)THEN
        DO JLBL=1,SIZE(XLBLUSER1)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER1(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             if(nverbia > 0)then
             print *,' ISO LABELLE ',ZISO
             endif
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE IF(NSUPER == 2)THEN
      IF(LBLUSER2)THEN
        DO JLBL=1,SIZE(XLBLUSER2)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER2(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE IF(NSUPER == 3)THEN
      IF(LBLUSER3)THEN
        DO JLBL=1,SIZE(XLBLUSER3)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER3(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE IF(NSUPER == 4)THEN
      IF(LBLUSER4)THEN
        DO JLBL=1,SIZE(XLBLUSER4)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER4(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE
      IF(.NOT.LABEL1)THEN
        IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
      ELSE
        IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
      ENDIF
    ENDIF
    ENDDO

  CASE(2,3)                            ! Given contour values     
    ICL=0
    DO I=1,1000
      ICL=ICL+1
! modifs for diachro
      IF(NIMNMX==3 .OR. (NIMNMX==2 .AND.LISOLEVP))THEN
        ZLEV(ICL)=ZISOLEVP(ICL)
        IF(NVERBIA > 5)then
          print *,' ** imcou ICL ZLEV ',ICL,ZLEV(ICL)
        ENDIF
      ELSE  IF (NIMNMX==2 .AND. .NOT.LISOLEVP) THEN 
        IF(I == 1 .AND. XISOLEV(1) == 9999.)THEN
          print *,' NIMNMX=2 . ABSENCE DE VALEURS DANS XISOLEV='
          print *,' RENTREZ LES AU CLAVIER PAR ORDRE CROISSANT ET A RAISON D''1'
          print *,' VALEUR PAR LIGNE. TERMINEZ PAR 9999.'
          print *,' (REMARQUE : elles ne sont pas memorisees et donc valides pour le seul parametre'
          print *,' en cours :',YTEXT(1:LEN_TRIM(YTEXT)),')'
        ENDIF
        IF(XISOLEV(1) == 9999.)THEN
          READ(5,*)ZLEV(ICL)
        ELSE
          ZLEV(ICL)=XISOLEV(ICL)
        ENDIF
      ENDIF
      IF(ZLEV(ICL) == 9999.)EXIT
    ENDDO
    IF(NVERBIA > 5) PRINT*,'ICL= ',ICL
    ICL=ICL-1
    CALL CPSETI('NCL',ICL)
!   IF(LCOLAREA .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))CALL CPSETI('NCL',ICL+1)
    CALL CPSETI('CLS',0)
    DO I=1,ICL
      CALL CPSETI('PAI',I)
      CALL CPSETI('AIA',I+1)
      CALL CPSETI('AIB',I)
      CALL CPSETR('CLV',ZLEV(I))
      CALL CPSETR('CLU',1.)
      IF(.NOT.LSUPER.OR. (LSUPER .AND. NSUPER == 1))THEN
        IF(LBLUSER1)THEN
          DO JLBL=1,SIZE(XLBLUSER1)
           DO JL=-20,20,1
             IF(ZLEV(I) == XLBLUSER1(JLBL)*10.**FLOAT(JL))THEN
               CALL CPSETR('CLU',3.)
               if(nverbia > 0)then
                 print *,' ISO LABELLE ',ZLEV(I)
               endif
               EXIT
             ENDIF
           ENDDO
          ENDDO
        ELSE
          IF(.NOT.LABEL1)THEN
            IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
          ELSE
            IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
          ENDIF
        ENDIF
      ELSE IF(NSUPER == 2)THEN
        IF(LBLUSER2)THEN
          DO JLBL=1,SIZE(XLBLUSER2)
           DO JL=-20,20,1
             IF(ZLEV(I) == XLBLUSER2(JLBL)*10.**FLOAT(JL))THEN
               CALL CPSETR('CLU',3.)
               EXIT
             ENDIF
           ENDDO
          ENDDO
        ELSE
          IF(.NOT.LABEL1)THEN
            IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
          ELSE
            IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
          ENDIF
        ENDIF
      ELSE IF(NSUPER == 3)THEN
        IF(LBLUSER3)THEN
          DO JLBL=1,SIZE(XLBLUSER3)
           DO JL=-20,20,1
             IF(ZLEV(I) == XLBLUSER3(JLBL)*10.**FLOAT(JL))THEN
               CALL CPSETR('CLU',3.)
               EXIT
             ENDIF
           ENDDO
          ENDDO
        ELSE
          IF(.NOT.LABEL1)THEN
            IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
          ELSE
            IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
          ENDIF
        ENDIF
      ELSE IF(NSUPER == 4)THEN
        IF(LBLUSER4)THEN
          DO JLBL=1,SIZE(XLBLUSER4)
           DO JL=-20,20,1
             IF(ZLEV(I) == XLBLUSER4(JLBL)*10.**FLOAT(JL))THEN
               CALL CPSETR('CLU',3.)
               EXIT
             ENDIF
           ENDDO
          ENDDO
        ELSE
          IF(.NOT.LABEL1)THEN
            IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
          ELSE
            IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
          ENDIF
        ENDIF
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(I,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(I-1,NULBLL+1)==0).OR.I==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ENDDO
!
END SELECT
!
!*    2.5    Further Conpack cosmetic parameters
!
SELECT CASE(NIOFFD)
  CASE(0)                 !! No label normalisation, decimal point kept
    III=8                 !
    CALL CPSETI('NEU',III)! 'Numeric exponent use flag'
    CALL CPSETI('NOF',7)! 
    CALL CPSETI('NET',0)  ! Exponent shown as "E"
			  ! III > 0 --> decimal point kept if the number of
                          ! significant digits is  << III; else form requiring
                          ! the fewest character is used
    IF(NSD /= 0)THEN
      CALL CPSETI('NSD',-NSD)  ! Nb de digits significatifs
    ELSE
      CALL CPSETI('NSD',-6)  ! Nb de digits significatifs
    ENDIF
  CASE DEFAULT            !! Label normalization, exponent to the right 
    CALL CPSETI('NEU',-2) ! Exponent notation forced in any case
    CALL CPSETI('NOF',7)! 
    CALL CPSETI('NET',0)  ! Exponent shown as "E"
END SELECT
!
!*   2.6      Special value handling
!
SELECT CASE(NIOFFP)
    
  CASE(0)                     ! No special value used
    CALL CPSETR('SPV',0.)
  CASE DEFAULT                ! XSPVAL used as a special value
    CALL CPSETR('SPV',XSPVAL)

END SELECT
!
!*   2.7     Information label under the plot
!
SELECT CASE(NIOFFM)
    
  CASE(0)                    ! a label is printed under the plot
  CASE DEFAULT               ! no label
    CALL CPSETC('ILT',' ')

END SELECT

ZTEMV=PTABV
CALL CPSETR('SPV',XSPVAL)
!
!*   2.8      Conpack initialization
!
!-----------------------------------------------------------------------------
  IF(LPVT .OR. LPXT .OR. LPYT)THEN
    ILMAX=NLMAX
    NLMAX=SIZE(PTABV,1)
  ENDIF
!-----------------------------------------------------------------------------
IF(NIMNMX <= 0)THEN

  ZTEMV2=ZTEMV
  IF(.NOT.LPXT .AND. .NOT.LPYT)THEN
    IF(LPRESY)THEN
! En log(P)
      WHERE(XZWORKZ(1:NLMAX,1:SIZE(ZTEMV2,2)) < XHMAX+ZDIXEPS)
      ZTEMV2=XSPVAL
      END WHERE
      WHERE(XZWORKZ(1:NLMAX,1:SIZE(ZTEMV2,2)) > XHMIN-ZDIXEPS)
      ZTEMV2=XSPVAL
      END WHERE
    ELSE
      WHERE(XZWORKZ(1:NLMAX,1:SIZE(ZTEMV2,2)) > XHMAX+ZDIXEPS)
      ZTEMV2=XSPVAL
      END WHERE
      WHERE(XZWORKZ(1:NLMAX,1:SIZE(ZTEMV2,2)) < XHMIN-ZDIXEPS)
      ZTEMV2=XSPVAL
      END WHERE
    ENDIF
  ENDIF

!print *,' ZTEMV2'
!print *,ZTEMV2
!print *,' XHMIN  XHMAX ',XHMIN-ZDIXEPS,XHMAX+ZDIXEPS
!print *,XZWORKZ(1,1:IKU)

if(nverbia > 0)then
  print *,' BALISE1 IMCOU'
endif
  CALL CPRECT(ZTEMV2,NLMAX,NLMAX,SIZE(ZTEMV2,2),ZRWRK,JPLRWK,IWRK,JPLIWK)
! CALL CPRECT(ZTEMV2,NLMAX,NLMAX,IKU,ZRWRK,JPLRWK,IWRK,JPLIWK)
  CALL CPPKCL(ZTEMV2,ZRWRK,IWRK)
  CALL CPGETI('NCL',INCL2)
!Janv 2001
! print *,' INCL2 ZTEMV2 ',INCL2
  IF(ALLOCATED(ZZCLV2))THEN
    DEALLOCATE(ZZCLV2)
  ENDIF
  ALLOCATE(ZZCLV2(INCL2))
!Janv 2001
  DO J=1,INCL2
    CALL CPSETI('PAI',J)
    CALL CPGETR('CLV',ZCLV2)
!Janv 2001
!   PRINT *,' ZCLV2 ',ZCLV2
    ZZCLV2(J)=ZCLV2
!Janv 2001
    IF(J == 1)ZCLVD=ZCLV2
    IF(J == INCL2)ZCLVF=ZCLV2
  ENDDO
END IF
!Janv 2001
!print *,' ZCLVD ZCLVF ',ZCLVD,ZCLVF

CALL CPRECT(ZTEMV,NLMAX,NLMAX,SIZE(ZTEMV,2),ZRWRK,JPLRWK,IWRK,JPLIWK)

!CALL CPRECT(ZTEMV,NLMAX,NLMAX,IKU,ZRWRK,JPLRWK,IWRK,JPLIWK)
CALL CPSETR('CWM',XSIZEL/.01)
if(nverbia > 0)then
  print *,' BALISE2 IMCOU NLMAX',NLMAX
endif
!-----------------------------------------------------------------------------
IF(LPVT .OR. LPXT .OR. LPYT)THEN
  NLMAX=ILMAX
ENDIF
if(nverbia > 0)then
  print *,' BALISE3 IMCOU INCL2= ',INCL2
endif
!-----------------------------------------------------------------------------
INCL=0
CALL CPPKCL(ZTEMV,ZRWRK,IWRK)
! Janv 2001
!CALL CPGETI('NCL',INCL)
IF(LCVZOOM)THEN
  IF(NIMNMX <= 0)THEN
    CALL CPSETI('CLS',0)
    IF(INCL2==0)THEN
      CALL CPSETI('NCL',1)
    ELSE
      CALL CPSETI('NCL',INCL2)
    ENDIF
    DO J=1,INCL2
      CALL CPSETI('PAI',J)
      CALL CPSETR('CLV',ZZCLV2(J))
    ENDDO
  ENDIF
! DEALLOCATE(ZZCLV2)
ENDIF
CALL CPGETI('NCL',INCL)
! Janv 2001
if(nverbia > 0)then
  print *,' BALISE3a IMCOU LCVZOOM= ',LCVZOOM
endif
!
!*   2.9      High and low handling
!
SELECT CASE(NHI)
    
  CASE(0)                           ! H + L   are displayed
    IF(INCL /= 0)THEN
      CALL CPLBDR(ZTEMV,ZRWRK,IWRK)
    ENDIF
  CASE DEFAULT                      ! TO BE REVISED*********************
			            ! <0  --> no action (:-1 to be set)
			            ! >0  --> gridpoint value displayed
                                    !         (1: to be set)
END SELECT
!
!print *,' ZTEMV in IMCOU_FORDIACHRO 2.9'    ! Technical message for developper's need
!!print *,ZTEMV
!*   2.10     Line style and color handling 
!
! Janv 2001
IF(NIMNMX <= 0)THEN
!IF(NIMNMX < 0)THEN
  DO J=1,INCL
    CALL CPSETI('PAI',J)
    CALL CPSETR('CLU',1.)
    CALL CPGETR('CLV',ZISO)
    IF(.NOT.LSUPER.OR. (LSUPER .AND. NSUPER == 1))THEN
      IF(LBLUSER1)THEN
        DO JLBL=1,SIZE(XLBLUSER1)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER1(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             if(nverbia > 0)then
               print *,' ISO LABELLE ',ZISO
             endif
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(J,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(J-1,NULBLL+1)==0).OR.J==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE IF(NSUPER == 2)THEN
      IF(LBLUSER2)THEN
        DO JLBL=1,SIZE(XLBLUSER2)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER2(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(J,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(J-1,NULBLL+1)==0).OR.J==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE IF(NSUPER == 3)THEN
      IF(LBLUSER3)THEN
        DO JLBL=1,SIZE(XLBLUSER3)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER3(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(J,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(J-1,NULBLL+1)==0).OR.J==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE IF(NSUPER == 4)THEN
      IF(LBLUSER4)THEN
        DO JLBL=1,SIZE(XLBLUSER4)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER4(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
             EXIT
           ENDIF
         ENDDO
        ENDDO
      ELSE
        IF(.NOT.LABEL1)THEN
          IF((MOD(J,NULBLL+1)==0))CALL CPSETR('CLU',3.)
        ELSE
          IF((MOD(J-1,NULBLL+1)==0).OR.J==1)CALL CPSETR('CLU',3.)
        ENDIF
      ENDIF
    ELSE
      IF(.NOT.LABEL1)THEN
        IF((MOD(J,NULBLL+1)==0))CALL CPSETR('CLU',3.)
      ELSE
        IF((MOD(J-1,NULBLL+1)==0).OR.J==1)CALL CPSETR('CLU',3.)
      ENDIF
    ENDIF
  ENDDO
END IF

if(nverbia > 0)then
  print *,' BALISE3b IMCOU '
endif
SELECT CASE(NDOT)
  
  CASE(0,1,1023,65535)        ! Solid line
      DO J=1,INCL
        CALL CPSETI('PAI',J)
        CALL CPSETI('CLD',65535)
      ENDDO
  CASE (:-1)                  !<0 Dashed negative values, 
                              !   solid positive values
    ICLD=ABS(NDOT)
!     write(0,*)' NDOT',NDOT,' INCL ',INCL
      DO J=1,INCL
        CALL CPSETI('PAI',J)
        CALL CPGETR('CLV',ZCLV)
        IF(ZCLV.GE.0.)CALL CPSETI('CLD',65535)
        IF(ZCLV.LT.0.)CALL CPSETI('CLD',ICLD)
!         write(0,*)' J ZCLV',J,ZCLV
      ENDDO

  CASE DEFAULT                ! NDOT used as a dash pattern
    ICLD=ABS(NDOT)
      DO J=1,INCL
        CALL CPSETI('PAI',J)
        CALL CPSETI('CLD',ICLD)
      ENDDO

END SELECT
!-----------------------------------------------------------------------------
!
! **************************************************************************
! Surfaces en hachures ou/et grises; LHACHx=.TRUE. avec x=1 ou 2 ou 3 ou 4)
! **************************************************************************

IF((LHACH1 .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))) .OR. &
   (LHACH2 .AND. NSUPER == 2)                                   .OR. &
   (LHACH3 .AND. NSUPER == 3)                                   .OR. &
   (LHACH4 .AND. NSUPER == 4))THEN !++++++++++++++++++++++++++++++++++++++++++

  IF(NSUPER > 1)THEN
    IH=IH+1
!   print *,' IHT IH ',IHT,IH
  ENDIF

  WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' VALEURS:'
  DO J=1,INCL
    CALL CPSETI('PAI',J)
    CALL CPSETI('AIB',J)
    CALL CPSETI('AIA',J+1)
    CALL CPGETR('CLV',ZCLV)
    ZLEV(J)=ZCLV
    CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
  ENDDO

  IF(.NOT.LHACHSEL)THEN
    IF(INCL+1 <= 8)THEN
      DO J=1,INCL
        IHACH(J)=INDHACHREF(J)
      ENDDO
      IHACH(INCL+1)=INDHACHREF(8)
    ELSE
      IHACH(1:2)=INDHACHREF(1:2)
      IHACH(3)=INDHACHREF(2)
      IHACH(INCL-1:INCL+1)=INDHACHREF(6:8)

      IF(INCL+1 < 13)THEN
        IHACH(4)=INDHACHREF(3)
      ELSE
        IHACH(4)=INDHACHREF(2)
      ENDIF

      IF(INCL+1 == 9)THEN
        IHACH(5)=INDHACHREF(4)
        IHACH(6)=INDHACHREF(5)
      ELSE
        IHACH(5)=INDHACHREF(3)
        IF(INCL+1 < 13)THEN
          IHACH(6)=INDHACHREF(4)
        ELSE
          IHACH(6)=INDHACHREF(3)
        ENDIF
      ENDIF

      IF(INCL+1 == 10)THEN
        IHACH(7)=INDHACHREF(5)
      ELSE IF(INCL+1 >= 11 .AND. INCL+1 < 14)THEN
        IHACH(7)=INDHACHREF(4)
      ELSE IF(INCL+1 >= 14)THEN
        IHACH(7)=INDHACHREF(3)
      ENDIF

      IF(INCL+1 >= 11 .AND. INCL+1 < 13)THEN
        IHACH(8)=INDHACHREF(5)
      ELSE IF(INCL+1 >= 13)THEN
        IHACH(8)=INDHACHREF(4)
      ENDIF

      IF(INCL+1 >= 12 .AND. INCL+1 < 14)THEN
        IHACH(9)=INDHACHREF(5)
      ELSE IF(INCL+1 >= 14)THEN
        IHACH(9)=INDHACHREF(4)
      ENDIF

      IF(INCL+1 == 13)THEN
        IHACH(10)=INDHACHREF(5)
      ELSE IF(INCL+1 >= 14 .AND. INCL+1 < 15)THEN
        IHACH(10)=INDHACHREF(5)
      ELSE IF(INCL+1 >= 15)THEN
        IHACH(10)=INDHACHREF(4)
      ENDIF

      IF(INCL+1 >= 14)THEN
        IHACH(11)=INDHACHREF(5)
      ENDIF

      IF(INCL+1 >= 15)THEN
        IHACH(12)=INDHACHREF(5)
      ENDIF

      IF(INCL+1 == 16)THEN
        IHACH(13)=INDHACHREF(5)
      ENDIF
    ENDIF

  ELSE

    DO J=1,300
      IHACH(J)=0
    ENDDO
    WRITE(NLUOUT,*)' >>>>>>>SELECTION DES GRISES ET HACHURES PAR L''UTILISATEUR'
    WRITE(NLUOUT,*)' >>>>>>>VOUS DEVEZ FOURNIR ',INCL+1,' INDICES'
    WRITE(NLUOUT,*)' Rentrez sur 1 premiere ligne le nombre d''indices fournis dans la ligne suivante'
    WRITE(NLUOUT,*)' Puis sur la(es) ligne(s) suivante(s) les indices des grises ou hachures' 
    WRITE(NLUOUT,*)' pris dans la table de reference (de grises ou hachures)'
    WRITE(NLUOUT,*)' correspondant aux isocontours ranges par ordre croissant'
    WRITE(NLUOUT,*)' (Entiers separes par 1 blanc)'
    READ(5,*,END=10)INBC
    GO TO 11
    10 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez le nombre d indices '
    READ(5,*)INBC
    11 CONTINUE
    WRITE(YCAR80,*)INBC
    !WRITE(NDIR,'(A80)')YCAR80
    CALL WRITEDIR(NDIR,YCAR80)
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
    READ(5,*,END=12)(IHACH(J),J=1,INBC)
    GO TO 13
    12 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez la valeur des indices '
    READ(5,*)(IHACH(J),J=1,INBC)
    13 CONTINUE
!    WRITE(YCAR160,*)IHACH(1:INBC)
!    YCAR160=ADJUSTL(YCAR160)
!    IF(LEN_TRIM(YCAR160) > 80 .OR. INBC > 20)THEN
     IF(INBC > 20)THEN
!Juillet 99
!      WRITE(YCAR80,'(20I4)')IHACH(1:INBC/2)
!     WRITE(YCAR80,*)IHACH(1:INBC/2)
      !WRITE(NDIR,'(A80)')YCAR80
      CALL WRITEDIR(NDIR,IHACH(1:INBC/2))
!      WRITE(YCAR80,'(20I4)')IHACH(INBC/2+1:INBC)
!     WRITE(YCAR80,*)IHACH(INBC/2+1:INBC)
      !WRITE(NDIR,'(A80)')YCAR80
      CALL WRITEDIR(NDIR,IHACH(INBC/2+1:INBC))
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
    ELSE
 !     WRITE(YCAR80,'(20I4)')IHACH(1:INBC)
!     WRITE(YCAR80,*)IHACH(1:INBC)
      !WRITE(NDIR,'(A80)')YCAR80
      CALL WRITEDIR(NDIR,IHACH(1:INBC))
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
    ENDIF
  ENDIF

  IF(LCOLZERO)THEN
    IHACH(NCOLZERO)=0
  ENDIF
  WRITE(NLUOUT,*)(ZLEV(J),IHACH(J),J=1,INCL)
  WRITE(NLUOUT,*)IHACH(INCL+1)

! Trace des zones hachurees
    CALL GSFAIS(1)
    CALL GSLN(1)
!   CALL GSFACI(1)
    CALL GSPLCI(1)
    CALL ARINAM(IIMAP,JPMAP)
    CALL CPCLAM(ZTEMV,ZRWRK,IWRK,IIMAP)
    CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,SFILLH)
    print *,' Hach: MAP 1 6 5 ',IIMAP(1),IIMAP(6),IIMAP(5)
    CALL GSFAIS(0)
!
! Trace des valeurs

    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
    CALL GSFAIS(1)
    CALL LBSETI('CBL',1)
!   CALL LBSETI('CBL',0)
    DO J=1,INCL
      YLLBS(J)=ADJUSTL(YLLBS(J))
    ENDDO
    IF(NIMNMX <= 0)THEN
      DO J=1,INCL
	IF(ZLEV(J).GT.ZCLVD)EXIT
      ENDDO
      JJD=MAX(1,J-1)
      DO J=INCL,1,-1
	IF(ZLEV(J).LE.ZCLVF)EXIT
      ENDDO
      JJF=MIN(INCL,J)
      INCL2=JJF-JJD+1
    ENDIF
    IF(.NOT.LSUPER .OR. NSUPER == 1)THEN
      IF(ZVR < .8999999)THEN
        print *,' ZVR < .9 ',ZVR
	IF(NIMNMX <= 0)THEN
          CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(ZVR+.2,1.)-ZVR)/10.,MIN(ZVR+.2,1.),&
          ZVB,ZVT,INCL2+1,.15,1.,IHACH(JJD:JJF+1),2,YLLBS(JJD:JJF),INCL2,1)
	ELSE
          CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(ZVR+.2,1.)-ZVR)/10.,MIN(ZVR+.2,1.),ZVB,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,1)
	ENDIF
      ELSE
        IF(INCL <= 8)THEN
          print *,' INCL <= 8 ',INCL
	  IF(NIMNMX <= 0)THEN
            CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB+(ZVT-ZVB)/4.,ZVT,&
            INCL2+1,.15,1.,IHACH(JJD:JJF+1),2,YLLBS(JJD:JJF),INCL2,1)
	  ELSE
            CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB+(ZVT-ZVB)/4.,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,1)
          ENDIF
        ELSE
          print *,' INCL > 8 ',INCL
	  IF(NIMNMX <= 0)THEN
            CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL2+1,.15,1.,IHACH(JJD:JJF+1),2,YLLBS(JJD:JJF),INCL2,1)
	  ELSE
            CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,1)
	  ENDIF
        ENDIF
      ENDIF

    ELSE

!      IF(NSUPERDIA > 2)THEN
!        ZVERA=ZVR-(ZVR-ZVL)/4.
!      ELSE
!        ZVERA=ZVR-(ZVR-ZVL)/3.
!      ENDIF
!      ZINTE=(ZVERA-ZVLDEF)/FLOAT(IHT)
!      IF(IHT == 1)THEN
!	ZD=ZVL; ZF=ZVERA
!      ELSE IF(IHT == 2 .OR. IHT == 3)THEN
!	ZD=ZVLDEF+ZINTE*(IH-1)
!	ZF=ZVLDEF+ZINTE*(IH)-.01
!      ENDIF
      IF(NSUPERDIA > 2)THEN
	ZVLDEF=.05
	ZINTE=.26
      ELSE
	ZVLDEF=.1
	ZINTE=.40
      ENDIF
      ZD=ZVLDEF+ZINTE*(NSUPER-2)
      ZF=ZD+ZINTE-.02
      IF(NIMNMX <= 0)THEN
        IF(INCL2 == 1)THEN
          ZF=ZF-(ZF-ZD)/2.
        ELSE IF(INCL2 <= 4)THEN
          ZF=ZF-(ZF-ZD)/4.
        ENDIF
      ELSE
        IF(INCL == 1)THEN
          ZF=ZF-(ZF-ZD)/2.
        ELSE IF(INCL <= 4)THEN
          ZF=ZF-(ZF-ZD)/4.
        ENDIF
      ENDIF
      IF(NIMNMX <= 0)THEN
        CALL LBLBAR_FORDIACHRO(0,ZD,ZF,ZVT+.01,ZVT+.04,INCL2+1,1.,.33,IHACH(JJD:JJF+1),2,YLLBS(JJD:JJF),INCL2,2)
      ELSE
        CALL LBLBAR_FORDIACHRO(0,ZD,ZF,ZVT+.01,ZVT+.04,INCL+1,1.,.33,IHACH,2,YLLBS,INCL,2)
      ENDIF
    ENDIF

    CALL GSFAIS(0)
!
! Definition de la couleur des isos (0 -> blanc sur papier; 1 -> noir sur papier)
    IF(LISOWHI)CALL GSPLCI(0)
    IF(LISOWHI)CALL GSTXCI(0)

!
!
ELSE IF(LCOLAREA)THEN        !+++++++++++++++++++++++++++++++++++++++++++++++++

! **************************************************************************
! Surfaces couleur (reservees aux dessins avec ou sans superpositions; LCOLAREA=.TRUE.)
! **************************************************************************

  IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN        !00000000000000000000000000000000000000000000

! Selection automatique des couleurs par le programme
! ***************************************************
    IF(.NOT.LCOLAREASEL)THEN     !====================================
       CALL COLOR_FORDIACHRO(INCL+1,1)
       WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' VALEURS:'
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('AIB',J)
         CALL CPSETI('AIA',J+1)
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
         ICOL(J)=J+2
         CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
         if(nverbia >5)then
	   print *,' J ZLEV(J) ICOL(J) A ',J,ZLEV(J),ICOL(J)
	 endif
       ENDDO
       ICOL(INCL+1)=INCL+3
       if(nverbia >0)then
         print *,' ICOL(INCL+1) A ',ICOL(INCL+1),' LCOLBR ',LCOLBR
         print *,' LCOLZERO NCOLZERO ',LCOLZERO,NCOLZERO
       endif
       IF(LCOLBR)THEN
         IF(ZLEV(MAX(1,INCL)) > ZLEV(1) .AND. ICOL(INCL+1) > ICOL(1))THEN
           ALLOCATE(ICOL2(INCL+1))
           if(nverbia >0)then
             print *,' APRES ALLOCATE(ICOL2) '
           endif
           ICOL2(1:INCL+1)=ICOL(INCL+1:1:-1)
           ICOL(1:INCL+1)=ICOL2
!          ICOL(:)=ICOL2
           if(nverbia >0)then
             print *,' AVANT DEALLOCATE(ICOL2) '
           endif
           DEALLOCATE(ICOL2)
         END IF
       END IF
       if(nverbia >0)then
         print *,' LCOLZERO NCOLZERO ',LCOLZERO,NCOLZERO
       endif
       IF(LCOLZERO)THEN
	 ICOL(NCOLZERO)=0
       ENDIF
       if(nverbia >0)then
         print *,' **imcou NLUOUT ',NLUOUT
       endif
       WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
       WRITE(NLUOUT,*)ICOL(INCL+1)
    ELSE                         !====================================

! Selection des couleurs par l'utilisateur
! ****************************************

       IF(LTABCOLDEF)THEN
       ! Choix de la table de couleurs par defaut
         WRITE(NLUOUT,*)' <<< TABCOLDEF >>>'
         CALL TABCOL_FORDIACHRO
       ELSE
       ! Choix d'une table creee par l'utilisateur
         CALL FMLOOK(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
         IF(IRESP == -54)THEN
           YNAMTABCOL(1:32)=' '
           print *,' Entrez le nom de VOTRE TABLE de COULEURS '
! Lecture du nom de la table de couleurs (1 seule fois)
           READ(5,*,END=14)YNAMTABCOL
    GO TO 15
    14 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez le nom de VOTRE TABLE de COULEURS'
    READ(5,*)YNAMTABCOL
    15 CONTINUE
           YNAMTABCOL=ADJUSTL(YNAMTABCOL)
	   !WRITE(NDIR,'(A80)')YNAMTABCOL
           CALL WRITEDIR(NDIR,YNAMTABCOL)
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
! Janv 2001
           CALL FMLOOK(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
           IF(IRESP /= 0)THEN
! Janv 2001
           CALL CREATLINK('DIRCOL',YNAMTABCOL,'CREAT',NVERBIA)
           CALL FMATTR(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
           OPEN(UNIT=ILUCOL,FILE=YNAMTABCOL,FORM='FORMATTED')
! Janv 2001
           ENDIF
! Janv 2001
         END IF

         WRITE(NLUOUT,*)' <<< ',YNAMTABCOL,' >>>'
         REWIND (ILUCOL)
         CALL GQOPS(ISTA)
         CALL GQACWK(1,IER,INB,IWK)
!print *,' COLOR_FORDIACHRO AP GQACWK INB IWK ',INB,IWK
	 CALL GQOPWK(1,IER,INB,IWK)
! Lecture du nb de couleurs de la table, des index de couleur et des
! proportions relatives de rouge, vert, bleu
         READ(ILUCOL,*)INBCT
         DO J=1,INBCT
           READ(ILUCOL,*)IDX,RED,GREEN,BLUE
	   DO JU=1,INB
	   CALL GQOPWK(JU,IER,INBB,IWK)
	   IF(IWK == 9)THEN
	     CYCLE
	   ELSE
             CALL GSCR(IWK,IDX,RED,GREEN,BLUE)
!          CALL GSCR(1,IDX,RED,GREEN,BLUE)
           ENDIF
           ENDDO
         ENDDO
       ENDIF
       WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' VALEURS:'
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('AIB',J)
         CALL CPSETI('AIA',J+1)
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
         if(nverbia >5)then
	   print *,' J ZLEV(J) B ',J,ZLEV(J)
	 endif
         CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
       ENDDO
       DO J=1,300
         ICOL(J)=0
       ENDDO
! Pour 1 dessin donne, lecture du nb d'indices de couleurs et de leur valeur
! sur la ligne suivante
       READ(5,*,END=16)INBC
    GO TO 17
    16 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez le nb d indices de couleur'
    READ(5,*)INBC
    17 CONTINUE
      ! WRITE(YCAR80,*)INBC
       !WRITE(NDIR,'(A80)')YCAR80
       CALL WRITEDIR(NDIR,INBC)
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
       READ(5,*,END=18)(ICOL(J),J=1,INBC)
    GO TO 19
    18 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez la valeur des indices de couleur'
    READ(5,*)(ICOL(J),J=1,INBC)
    19 CONTINUE
!       WRITE(YCAR160,*) ICOL(1:INBC)
!       YCAR160=ADJUSTL(YCAR160)
!       IF(LEN_TRIM(YCAR160) > 80 .OR. INBC > 20)THEN
        IF(INBC > 20)THEN
! Juillet 99
       !  WRITE(YCAR80,'(20I4)')ICOL(1:INBC/2)
!        WRITE(YCAR80,*)ICOL(1:INBC/2)
         !WRITE(NDIR,'(A80)')YCAR80
         CALL WRITEDIR(NDIR,ICOL(1:INBC/2))
        ! WRITE(YCAR80,'(20I4)')ICOL(INBC/2+1:INBC)
!        WRITE(YCAR80,*)ICOL(INBC/2+1:INBC)
         !WRITE(NDIR,'(A80)')YCAR80
         CALL WRITEDIR(NDIR,ICOL(INBC/2+1:INBC))
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
       ELSE
! Juillet 99
       !  WRITE(YCAR80,'(20I4)')ICOL(1:INBC)
!        WRITE(YCAR80,*)ICOL(1:INBC)
         !WRITE(NDIR,'(A80)')YCAR80
         CALL WRITEDIR(NDIR,ICOL(1:INBC))
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
       ENDIF
       WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
       WRITE(NLUOUT,*)ICOL(INCL+1)
! fin de la selection des couleurs par l'utilisateur
    ENDIF                        !====================================
!
! Trace des zones colorees
!*************************
    IF(LMARKER .AND. .NOT. LSPOT)THEN
    ! en etoiles colorees
      !IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT .AND. .NOT.LCVXZ .AND. .NOT.LCVYZ)THEN
      IF(.NOT.LPYT .AND. .NOT.LCVXZ .AND. .NOT.LCVYZ)THEN
      CALL GSMK(3)  ! asterisk is the type of marker
      DO JJ=1,SIZE(ZTEMV,2)
      DO JI=1,SIZE(ZTEMV,1)
	IF(ZTEMV(JI,JJ) /= XSPVAL)THEN
	  IF(ZTEMV(JI,JJ) < ZLEV(1))THEN
	    CALL GSPMCI(ICOL(1))
	  ELSE IF(ZTEMV(JI,JJ) >= ZLEV(INCL))THEN
	    CALL GSPMCI(ICOL(INCL+1))
	  ELSE
	    DO J=1,INCL-1
	      IF(ZTEMV(JI,JJ) >= ZLEV(J) .AND. &
		 ZTEMV(JI,JJ) < ZLEV(J+1))THEN
		CALL GSPMCI(ICOL(J+1))
		EXIT
              ENDIF
	    ENDDO
	  ENDIF
	  ZX=XZZDS(JI)
	  ZY=XZWORKZ(JI,JJ)
	  CALL GPM(1,ZX,ZY)
	ENDIF
      ENDDO
      ENDDO
      ELSE
        print *,'pas de LMARKER teste pour ce type de tracé (PYT, 2D vert //X ou 2D vert //Y)'
        print *,'essayer en modifiant le test IF(.NOT.LPVT... dans imcou_fordiachro'
      ENDIF

    ELSE IF (LSPOT .AND. .NOT. LMARKER) THEN
    ! en paves de couleur
      !IF(.NOT.LPVT .AND. .NOT.LPXT .AND. .NOT.LPYT .AND. .NOT.LCVXZ .AND. .NOT.LCVYZ)THEN
      IF(.NOT.LPYT .AND. .NOT.LCVXZ .AND. .NOT.LCVYZ)THEN
      CALL  GSFAIS(1)  ! solid filling of the polygon
      IND=SIZE(ZTEMV,1)
      ZEPX=(XZZDS(IND/2+1)-XZZDS(IND/2))*0.5
      print *,'LSPOT: contour du pave en noir ?'
      print *,'       (o/O/y/Y recommande pour trace d observations '
      print *,'        epaisseur du contour gere avec XLW1)'
      read(5,*) YREP
      CALL WRITEDIR(NDIR,YREP)
      IF(YREP=='o' .OR. YREP=='O' .OR. YREP=='y' .OR. YREP=='Y') THEN
        ! contour en trait plein noir
        CALL DASHDB(65535)
      END IF
      DO JJ=1,SIZE(ZTEMV,2)-1
      DO JI=1,SIZE(ZTEMV,1)
        IF (JJ==1) THEN
          ZEPYD= XZWORKZ(JI,JJ) - ZWZ(JI+1) ! ZWZ(1:NLMAX+2)
        ELSE
          ZEPYD=XZWORKZ(JI,JJ) - (XZWORKZ(JI,JJ)+XZWORKZ(JI,JJ-1))*0.5
        ENDIF
        IF (JJ==SIZE(ZTEMV,2)-1) THEN
          ZEPYU=0
        ELSE
          ZEPYU=(XZWORKZ(JI,JJ+1)+XZWORKZ(JI,JJ))*0.5 - XZWORKZ(JI,JJ)
        ENDIF
        IF(ZTEMV(JI,JJ) /= XSPVAL)THEN
          IF(ZTEMV(JI,JJ) < ZLEV(1))THEN
            CALL GSFACI(ICOL(1))
	  ELSE IF(ZTEMV(JI,JJ) >= ZLEV(INCL))THEN
	    CALL GSFACI(ICOL(INCL+1))
          ELSE
            DO J=1,INCL-1
              IF(ZTEMV(JI,JJ) >= ZLEV(J) .AND. &
                 ZTEMV(JI,JJ) < ZLEV(J+1))THEN
                CALL GSFACI(ICOL(J+1))
                EXIT
              ENDIF
            ENDDO
          ENDIF
          ZX5(1)=XZZDS(JI)-ZEPX ; ZY5(1)=XZWORKZ(JI,JJ)-ZEPYD
          ZX5(2)=XZZDS(JI)-ZEPX ; ZY5(2)=XZWORKZ(JI,JJ)+ZEPYU
          ZX5(3)=XZZDS(JI)+ZEPX ; ZY5(3)=XZWORKZ(JI,JJ)+ZEPYU
          ZX5(4)=XZZDS(JI)+ZEPX ; ZY5(4)=XZWORKZ(JI,JJ)-ZEPYD
          ZX5(5)=XZZDS(JI)-ZEPX ; ZY5(5)=XZWORKZ(JI,JJ)-ZEPYD
          ! paves
          CALL GFA(5,ZX5,ZY5)
          IF(YREP=='o' .OR. YREP=='O' .OR. YREP=='y' .OR. YREP=='Y') THEN
            ! contour
            CALL GQLWSC(IER,ZWIDTH)
            CALL GSLWSC(XLWIDTH)
            CALL CURVED(ZX5,ZY5,5)
            CALL GSLWSC(ZWIDTH)
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      ELSE
        print *,'pas de LSPOT teste pour ce type de tracé (PYT, 2D vert //X ou 2D vert //Y)'
        print *,'essayer en modifiant le test IF(.NOT.LPVT... dans imcou_fordiachro'
      ENDIF
    ELSE
    ! Trace des surfaces colorees
    CALL GSFAIS(1)
    CALL ARINAM(IIMAP,JPMAP)
    CALL CPCLAM(ZTEMV,ZRWRK,IWRK,IIMAP)
    CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,CCOLR)
    print *,' Col: MAP 1 6 5 ',IIMAP(1),IIMAP(6),IIMAP(5)
    CALL GSPLCI(1)
    CALL GSFAIS(0)
!   CALL GSLN(1)
    ENDIF
    ! Trace de la palette de couleurs (legende)
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
    CALL GSFAIS(1)
    CALL LBSETI('CBL',0)
    DO J=1,INCL
      YLLBS(J)=ADJUSTL(YLLBS(J))
    ENDDO
    IF(NIMNMX <= 0)THEN
      DO J=1,INCL
        IF(ZLEV(J).GT.ZCLVD)EXIT
      ENDDO
      JJD=MAX(1,J-1)
      DO J=INCL,1,-1
        IF(ZLEV(J).LE.ZCLVF)EXIT
      ENDDO
      JJF=MIN(INCL,J)
      INCL2=JJF-JJD+1
!print *,'ZLEV(1:INCL) ',ZLEV(1:INCL)
!print *,' JJD JJF ZLEV(JJD:JJF) ',ZLEV(JJD:JJF)
      CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(1.-ZVR,.2))/10.,1.,ZVB,ZVT,INCL2+1,.15,1.,ICOL(JJD:JJF+1),1,YLLBS(JJD:JJF),INCL2,1)
!     CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL2+1,.15,1.,ICOL(JJD:JJF+1),1,YLLBS(JJD:JJF),INCL2,1)
    ELSE
      CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(1.-ZVR,.2))/10.,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
!     CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
    END IF
    CALL GSFAIS(0)
!
! Definition de la couleur des isos (0 -> blanc sur papier; 1 -> noir sur papier)
    IF(LISOWHI)CALL GSPLCI(0)
    IF(LISOWHI)CALL GSTXCI(0)
!
  ELSE IF(LCOLINE)THEN       !00000000000000000000000000000000000000000000

! Traits couleur dans le cas de superpositions (LCOLAREA=.TRUE. et LCOLINE=.TRUE.)
! **************************************************************************
! Modifs 220396
    IF((LSUPER .AND. NSUPER == 1) .OR. .NOT.LSUPER)CALL TABCOL_FORDIACHRO
    IF(LSUPER)THEN
!Mars 2000
      IF(LCOLISONE)THEN
	IF(NSUPER == 1)CALL GSPLCI(NCOLISONE1)
	IF(NSUPER == 1)CALL GSTXCI(NCOLISONE1)
	IF(NSUPER == 2)CALL GSPLCI(NCOLISONE2)
	IF(NSUPER == 2)CALL GSTXCI(NCOLISONE2)
	IF(NSUPER == 3)CALL GSPLCI(NCOLISONE3)
	IF(NSUPER == 3)CALL GSTXCI(NCOLISONE3)
	IF(NSUPER == 4)CALL GSPLCI(NCOLISONE4)
	IF(NSUPER == 4)CALL GSTXCI(NCOLISONE4)
	IF(NSUPER == 5)CALL GSPLCI(NCOLISONE5)
	IF(NSUPER == 5)CALL GSTXCI(NCOLISONE5)
      ELSE
!Mars 2000
      IF(NSUPER == 1)CALL GSPLCI(2)
      IF(NSUPER == 1)CALL GSTXCI(2)
      IF(NSUPER == 2)CALL GSPLCI(4)
      IF(NSUPER == 2)CALL GSTXCI(4)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==2)CALL GSPLCI(2)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==2)CALL GSTXCI(2)
      IF(NSUPER == 3)CALL GSPLCI(3)
      IF(NSUPER == 3)CALL GSTXCI(3)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==3)CALL GSPLCI(4)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==3)CALL GSTXCI(4)
      IF(NSUPER == 4)CALL GSPLCI(7)
      IF(NSUPER == 4)CALL GSTXCI(7)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==4)CALL GSPLCI(3)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==4)CALL GSTXCI(3)
!!!!!!!! PROVI
!CALL FRSTPT(XDS(1,NMGRID),XHMIN)
!CALL VECTOR(XDS(1,NMGRID),XHMAX)
!CALL VECTOR(XDS(NLMAX,NMGRID),XHMAX)
!CALL VECTOR(XDS(NLMAX,NMGRID),XHMIN)
!CALL VECTOR(XDS(1,NMGRID),XHMIN)
!!!!!!!! PROVI
!Mars 2000
      ENDIF
!Mars 2000
    END IF
  ELSE                       !00000000000000000000000000000000000000000000

! Traits noir et blanc dans le cas de superpositions (LCOLAREA=.TRUE. et LCOLINE=.FALSE.)
! ********************************************************************************
if(nverbia > 0)then
  print *,' BALISE3c IMCOU '
endif

    CALL GSPLCI(1)
    CALL GSLN(1)
    IF(LSUPER)THEN
      IF(NSUPER == 1)CALL GSLN(1)
      IF(NSUPER == 2)CALL GSLN(1)

      IF(LINVPTIR)THEN

        IF(NSUPER == 3)THEN
	  CALL GSLN(2)
	  IF((LCOLAREA.OR.LHACH1) .AND. LHACH2)CALL GSLN(1)
        ENDIF
        IF(NSUPER == 4)CALL GSLN(3)

      ELSE

        IF(NSUPER == 3)THEN
	  CALL GSLN(3)
	  IF((LCOLAREA.OR.LHACH1) .AND. LHACH2)CALL GSLN(1)
        ENDIF
        IF(NSUPER == 4)CALL GSLN(2)

      ENDIF

    END IF

  END IF                     !00000000000000000000000000000000000000000000

ELSE IF( LGREY .AND. .NOT.LCOLAREA )   THEN !++++++++++++++++++++++++++++++
! **************************************************************
! Surfaces en grises ( LGREY=.TRUE.)
!  En cas de superpositions, obligatoirement le 1er dessin
! **************************************************************
  IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN        !000000000000000000
!
! Selection automatique des grises par le programme
! **************************************************
!
  CALL COLOR_FORDIACHRO(INCL+1,2)
  WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' VALEURS:'
  DO J=1,INCL
    CALL CPSETI('PAI',J)
    CALL CPSETI('AIB',J)
    CALL CPSETI('AIA',J+1)
    CALL CPGETR('CLV',ZCLV)
    ZLEV(J)=ZCLV
    ICOL(J)=J+2
    CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
  ENDDO
  ICOL(INCL+1)=INCL+3
       if(nverbia >0)then
         print *,' Grey: ICOL(INCL+1) A ',ICOL(INCL+1),' LCOLBR ',LCOLBR
       endif
  IF(LCOLBR)THEN
    IF(ZLEV(MAX(1,INCL)) > ZLEV(1) .AND. ICOL(INCL+1) > ICOL(1))THEN
      ALLOCATE(ICOL2(INCL+1))
      ICOL2(1:INCL+1)=ICOL(INCL+1:1:-1)
      ICOL(1:INCL+1)=ICOL2
!          ICOL(:)=ICOL2
      DEALLOCATE(ICOL2)
    END IF
  END IF
       if(nverbia >0)then
         print *,' Grey: LCOLZERO NCOLZERO ',LCOLZERO,NCOLZERO
       endif
  IF(LCOLZERO)THEN
    ICOL(NCOLZERO)=0
  ENDIF
  WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
  WRITE(NLUOUT,*)ICOL(INCL+1)
  ! Trace des surfaces grisees
  CALL GSFAIS(1)
  CALL ARINAM(IIMAP,JPMAP)
  CALL CPCLAM(ZTEMV,ZRWRK,IWRK,IIMAP)
  CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,CCOLR)
  print *,' Grey: MAP 1 6 5 ',IIMAP(1),IIMAP(6),IIMAP(5)
  CALL GSPLCI(1)
  CALL GSFAIS(0)
!   CALL GSLN(1)
  ! Trace de la palette de couleurs (legende)
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  CALL GSFAIS(1)
  CALL LBSETI('CBL',0)
  DO J=1,INCL
    YLLBS(J)=ADJUSTL(YLLBS(J))
  ENDDO
  IF(NIMNMX <= 0)THEN
    DO J=1,INCL
      IF(ZLEV(J).GT.ZCLVD)EXIT
    ENDDO
    JJD=MAX(1,J-1)
    DO J=INCL,1,-1
      IF(ZLEV(J).LE.ZCLVF)EXIT
    ENDDO
    JJF=MIN(INCL,J)
    INCL2=JJF-JJD+1
    CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(1.-ZVR,.2))/10.,1.,ZVB,ZVT,INCL2+1,.15,1.,ICOL(JJD:JJF+1),1,YLLBS(JJD:JJF),INCL2,1)
!   CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL2+1,.15,1.,ICOL(JJD:JJF+1),1,YLLBS(JJD:JJF),INCL2,1)
  ELSE
    CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(1.-ZVR,.2))/10.,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
!   CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
  ENDIF
  CALL GSFAIS(0)
!
! Definition de la couleur des isos (0 -> blanc sur papier; 1 -> noir sur papier)
      IF(LISOWHI)CALL GSPLCI(0)
      IF(LISOWHI)CALL GSTXCI(0)
  
  ELSE IF(LCOLINE)THEN       !00000000000000000000000000000000000000000000

! Traits couleur dans le cas de superpositions (LGREY=.TRUE. et LCOLINE=.TRUE.)
! **************************************************************************
    CALL TABCOL_FORDIACHRO
    IF(LSUPER)THEN
!Mars 2000
      IF(LCOLISONE)THEN
	IF(NSUPER == 1)CALL GSPLCI(NCOLISONE1)
	IF(NSUPER == 1)CALL GSTXCI(NCOLISONE1)
	IF(NSUPER == 2)CALL GSPLCI(NCOLISONE2)
	IF(NSUPER == 2)CALL GSTXCI(NCOLISONE2)
	IF(NSUPER == 3)CALL GSPLCI(NCOLISONE3)
	IF(NSUPER == 3)CALL GSTXCI(NCOLISONE3)
	IF(NSUPER == 4)CALL GSPLCI(NCOLISONE4)
	IF(NSUPER == 4)CALL GSTXCI(NCOLISONE4)
	IF(NSUPER == 5)CALL GSPLCI(NCOLISONE5)
	IF(NSUPER == 5)CALL GSTXCI(NCOLISONE5)
      ELSE
!Mars 2000
      IF(NSUPER == 1)CALL GSPLCI(2)
      IF(NSUPER == 1)CALL GSTXCI(2)
      IF(NSUPER == 2)CALL GSPLCI(4)
      IF(NSUPER == 2)CALL GSTXCI(4)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==2)CALL GSPLCI(2)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==2)CALL GSTXCI(2)
      IF(NSUPER == 3)CALL GSPLCI(3)
      IF(NSUPER == 3)CALL GSTXCI(3)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==3)CALL GSPLCI(4)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==3)CALL GSTXCI(4)
      IF(NSUPER == 4)CALL GSPLCI(7)
      IF(NSUPER == 4)CALL GSTXCI(7)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==4)CALL GSPLCI(3)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==4)CALL GSTXCI(3)
!!!!!!!! PROVI
!CALL FRSTPT(XDS(1,NMGRID),XHMIN)
!CALL VECTOR(XDS(1,NMGRID),XHMAX)
!CALL VECTOR(XDS(NLMAX,NMGRID),XHMAX)
!CALL VECTOR(XDS(NLMAX,NMGRID),XHMIN)
!CALL VECTOR(XDS(1,NMGRID),XHMIN)
!!!!!!!! PROVI
!Mars 2000
      ENDIF
!Mars 2000
    END IF

  ELSE                       !00000000000000000000000000000000000000000000

! Traits noir et blanc dans le cas de superpositions (LGREY=.TRUE. et LCOLINE=.FALSE.)
! ********************************************************************************

    CALL GSPLCI(1)
    CALL GSLN(1)
    IF(LSUPER)THEN
      IF(NSUPER == 1)CALL GSLN(1)
      IF(NSUPER == 2)CALL GSLN(1)

      IF(LINVPTIR)THEN

        IF(NSUPER == 3)THEN
	  CALL GSLN(2)
	  IF((LGREY.OR.LHACH1) .AND. LHACH2)CALL GSLN(1)
        ENDIF
        IF(NSUPER == 4)CALL GSLN(3)

      ELSE

        IF(NSUPER == 3)THEN
	  CALL GSLN(3)
	  IF((LGREY.OR.LHACH1) .AND. LHACH2)CALL GSLN(1)
        ENDIF
        IF(NSUPER == 4)CALL GSLN(2)

      ENDIF

    END IF

  END IF                     !00000000000000000000000000000000000000000000
!

ELSE IF(LCOLINE)THEN    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
! **********************************************
! Traits couleur   (LCOLAREA=.FALSE. et LCOLINE=.TRUE.)
! **********************************************

! Cas de superpositions
! *********************
! Modifs 220395=6
  CALL TABCOL_FORDIACHRO
!   IF((LSUPER .AND. NSUPER == 1) .OR. .NOT.LSUPER)CALL TABCOL_FORDIACHRO
! Modifs 270198
! IF(LSUPER)THEN             !............................................
  IF(LSUPER .AND. &          !............................................
    !.NOT.((LHACH1.OR.LHACH2) .AND. NSUPERDIA == 2))THEN
     .NOT.((LHACH1.OR.LHACH2) .AND. NSUPERDIA == 2) .AND. &
     .NOT.( LARROVL .AND. NSUPERDIA == 2          )       )THEN

!Mars 2000
      IF(LCOLISONE)THEN
	IF(NSUPER == 1)CALL GSPLCI(NCOLISONE1)
	IF(NSUPER == 1)CALL GSTXCI(NCOLISONE1)
	IF(NSUPER == 2)CALL GSPLCI(NCOLISONE2)
	IF(NSUPER == 2)CALL GSTXCI(NCOLISONE2)
	IF(NSUPER == 3)CALL GSPLCI(NCOLISONE3)
	IF(NSUPER == 3)CALL GSTXCI(NCOLISONE3)
	IF(NSUPER == 4)CALL GSPLCI(NCOLISONE4)
	IF(NSUPER == 4)CALL GSTXCI(NCOLISONE4)
	IF(NSUPER == 5)CALL GSPLCI(NCOLISONE5)
	IF(NSUPER == 5)CALL GSTXCI(NCOLISONE5)
      ELSE
!Mars 2000

    IF(NSUPER == 1)CALL GSPLCI(2)
    IF(NSUPER == 1)CALL GSTXCI(2)
    IF(NSUPER == 2)CALL GSPLCI(4)
    IF(NSUPER == 2)CALL GSTXCI(4)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==2)CALL GSPLCI(2)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==2)CALL GSTXCI(2)
    IF(NSUPER == 3)CALL GSPLCI(3)
    IF(NSUPER == 3)CALL GSTXCI(3)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==3)CALL GSPLCI(4)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==3)CALL GSTXCI(4)
    IF(NSUPER == 4)CALL GSPLCI(7)
    IF(NSUPER == 4)CALL GSTXCI(7)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==4)CALL GSPLCI(3)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==4)CALL GSTXCI(3)
  
!Mars 2000
      ENDIF
!Mars 2000
  ELSE                       !............................................
! Pas de superpositions
! *********************

! Selection automatique des couleurs par le programme
! ***************************************************

    IF(.NOT.LCOLINESEL)THEN      !::::::::::::::::::::::::::::::::::::

!Mars 2000
       IF(LCOLISONE)THEN
	 ICOL(1:INCL)=NCOLISONE1
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('CLC',ICOL(J))
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
       ENDDO
       WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' COULEUR UNIQUE : ',ICOL(1)
       WRITE(NLUOUT,*)(ZLEV(J),J=1,INCL)
       ELSE
!Mars 2000

       CALL COLOR_FORDIACHRO(INCL,1)
       WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' VALEURS:'
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
         ICOL(J)=J+2
	 if(nverbia > 5)then
	   print *,' J ZLEV(J) ICOL(J) C ',J,ZLEV(J),ICOL(J)
	 endif
         CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
       ENDDO
       IF(LCOLBR)THEN
         IF(ZLEV(MAX(1,INCL)) > ZLEV(1) .AND. ICOL(INCL) > ICOL(1))THEN
           ALLOCATE(ICOL2(INCL))
           ICOL2(1:INCL)=ICOL(INCL:1:-1)
           ICOL(1:INCL)=ICOL2
!          ICOL(:)=ICOL2
           DEALLOCATE(ICOL2)
         END IF
       END IF
       WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('CLC',ICOL(J))
       ENDDO
!Mars 2000
     ENDIF
!Mars 2000

    ELSE                         !::::::::::::::::::::::::::::::::::::

! Selection des couleurs par l'utilisateur
! ****************************************

! Choix de la table de couleurs par defaut
! ****************************************

       IF(LTABCOLDEF)THEN
         WRITE(NLUOUT,*)' <<< TABCOLDEF >>>'
         CALL TABCOL_FORDIACHRO

       ELSE

! Choix d'une table creee par l'utilisateur
! *****************************************

         CALL FMLOOK(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
         IF(IRESP == -54)THEN
           YNAMTABCOL(1:32)=' '
! Lecture du nom de la table de couleurs (1 seule fois)
           print *,' Entrez le nom de VOTRE TABLE de COULEURS '
           READ(5,*,END=20)YNAMTABCOL
    GO TO 21
    20 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez le nom de VOTRE TABLE de COULEURS'
    READ(5,*)YNAMTABCOL
    21 CONTINUE
           YNAMTABCOL=ADJUSTL(YNAMTABCOL)
	   !WRITE(NDIR,'(A80)')YNAMTABCOL
           CALL WRITEDIR(NDIR,YNAMTABCOL)
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
! Janv 2001
           CALL FMLOOK(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
           IF(IRESP /= 0)THEN
! Janv 2001
           CALL CREATLINK('DIRCOL',YNAMTABCOL,'CREAT',NVERBIA)
           CALL FMATTR(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
           OPEN(UNIT=ILUCOL,FILE=YNAMTABCOL,FORM='FORMATTED')
! Janv 2001
           ENDIF
! Janv 2001
         END IF
         WRITE(NLUOUT,*)' <<< ',YNAMTABCOL,' >>>'
         REWIND (ILUCOL)
         CALL GQOPS(ISTA)
         CALL GQACWK(1,IER,INB,IWK)
!print *,' COLOR_FORDIACHRO AP GQACWK INB IWK ',INB,IWK
	 CALL GQOPWK(1,IER,INB,IWK)
! Lecture du nb de couleurs de la table, des index de couleur et des
! proportions relatives de rouge, vert, bleu
         READ(ILUCOL,*)INBCT
         DO J=1,INBCT
           READ(ILUCOL,*)IDX,RED,GREEN,BLUE
	   DO JU=1,INB
	   CALL GQOPWK(JU,IER,INBB,IWK)
	   IF(IWK == 9)THEN
	     CYCLE
	   ELSE
             CALL GSCR(IWK,IDX,RED,GREEN,BLUE)
!          CALL GSCR(1,IDX,RED,GREEN,BLUE)
           ENDIF
           ENDDO
         ENDDO
       END IF
! Pour 1 dessin donne, lecture du nb d'indices de couleurs et de leur valeur
! sur la ligne suivante
         DO J=1,300
           ICOL(J)=1
         ENDDO
         READ(5,*,END=22)INBC
    GO TO 23
    22 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez le nb d indices de couleur'
    READ(5,*)INBC
    23 CONTINUE
         !WRITE(YCAR80,*)INBC
         !WRITE(NDIR,'(A80)')YCAR80
         CALL WRITEDIR(NDIR,INBC)
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
         READ(5,*,END=24)(ICOL(J),J=1,INBC)
    GO TO 25
    24 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",YCAR20)
    YCAR20=ADJUSTL(YCAR20)
    OPEN(5,FILE=YCAR20)
    print *,' INTERACTIF : Entrez la valeur des indices de couleur'
    READ(5,*)(ICOL(J),J=1,INBC)
    25 CONTINUE
!         WRITE(YCAR160,*)ICOL(1:INBC)
!         YCAR160=ADJUSTL(YCAR160)
!         IF(LEN_TRIM(YCAR160) > 80 .OR. INBC > 20)THEN
          IF(INBC > 20)THEN

! Juillet 99
         !  WRITE(YCAR80,'(20I4)')ICOL(1:INBC/2)
!          WRITE(YCAR80,*)ICOL(1:INBC/2)
           !WRITE(NDIR,'(A80)')YCAR80
           CALL WRITEDIR(NDIR,ICOL(1:INBC/2))
           !WRITE(YCAR80,'(20I4)')ICOL(INBC/2+1:INBC)
!          WRITE(YCAR80,*)ICOL(INBC/2+1:INBC)
           !WRITE(NDIR,'(A80)')YCAR80
           CALL WRITEDIR(NDIR,ICOL(INBC/2+1:INBC))
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
         ELSE
          ! WRITE(YCAR80,'(20I4)')ICOL(1:INBC)
!          WRITE(YCAR80,*)ICOL(1:INBC)
           !WRITE(NDIR,'(A80)')YCAR80
           CALL WRITEDIR(NDIR,ICOL(1:INBC))
#ifdef RHODES
    CALL FLUSH(NDIR,ISTAF)
#else
    CALL FLUSH(NDIR)
#endif
         ENDIF
         DO J=1,INCL
           CALL CPSETI('PAI',J)
           CALL CPSETI('CLC',ICOL(J))
           CALL CPGETR('CLV',ZCLV)
           ZLEV(J)=ZCLV
	   if(nverbia > 5)then
	     print *,' J ZLEV(J) ICOL(J) D ',J,ZLEV(J),ICOL(J)
	   endif
           CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
         ENDDO
         WRITE(NLUOUT,*)' >>>>>>>IMCOU_FORDIACHRO VARIABLE : ',HTEXT,' NB ISOC. : ',INCL,' VALEURS:'
         WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)

    END IF                       !::::::::::::::::::::::::::::::::::::

!Mars 2000
       IF(LCOLISONE)THEN
       ELSE
!Mars 2000
       CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
       CALL GSFAIS(0)
       CALL SETUSV('MI',1)
       CALL SET(ZVR,1.,ZVB,ZVT,ZVR,1.,ZVB,ZVT,1)
       IF(NIMNMX <= 0)THEN
         DO J=1,INCL
           IF(ZLEV(J).GE.ZCLVD)EXIT
         ENDDO
         JJD=MAX(1,J)
         DO J=INCL,1,-1
           IF(ZLEV(J).LE.ZCLVF)EXIT
         ENDDO
         JJF=MIN(INCL,J)
         INCL2=JJF-JJD+1
	 IF(INCL2 <= 1)THEN
	   ZINTERV=0.
         ELSE
           ZINTERV=(ZVT-ZVB-.009)/(INCL2-1)
         ENDIF
	 CALL GSCLIP(0)
         DO J=JJD,JJF
           YLLBS(J)=ADJUSTL(YLLBS(J))
           CALL GSPLCI(ICOL(J))
           CALL GSTXCI(ICOL(J))
           if(nverbia > 0)then
             print *,' BALISE3d IMCOU '
           endif
	   IF(ZVR < .9 .AND. INCL < 25)THEN
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.015,0.,-1.)
	   ELSEIF (ZVR < .9 .AND. INCL < 30 .AND. INCL >= 25)THEN
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.012,0.,-1.)
	   ELSEIF (ZVR >= .95 )THEN
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.007,0.,-1.)
	   ELSE
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.009,0.,-1.)
           ENDIF
!          CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.007,0.,-1.)
         ENDDO
	 CALL GSCLIP(1)
       ELSE
	 IF(INCL <= 1)THEN
           ZINTERV=0.
         ELSE
           ZINTERV=(ZVT-ZVB-.009)/(INCL-1)
         ENDIF
	 CALL GSCLIP(0)
           if(nverbia > 0)then
             print *,' BALISE3e IMCOU '
           endif
         DO J=1,INCL
           YLLBS(J)=ADJUSTL(YLLBS(J))
           CALL GSPLCI(ICOL(J))
           CALL GSTXCI(ICOL(J))

	   IF(ZVR < .9 .AND. INCL < 25)THEN
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.015,0.,-1.)
	   ELSEIF (ZVR < .9 .AND. INCL < 30 .AND. INCL >= 25)THEN
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.012,0.,-1.)
	   ELSEIF (ZVR >= .95 )THEN
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.007,0.,-1.)
	   ELSE
             CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.009,0.,-1.)
	   ENDIF
!          CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.007,0.,-1.)
         ENDDO
	 CALL GSCLIP(1)
       END IF
       CALL SETUSV('MI',IMI)
       CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!Mars 2000
      ENDIF
!Mars 2000
       CALL GSTXCI(1)
       CALL GSPLCI(1)
       

  END IF                     !............................................

ELSE                    !+++++++++++++++++++++++++++++++++++++++++++++++++++++
if(nverbia > 0)then
  print *,' BALISE3f IMCOU'
endif

!***************************************************
! Traits noir et blanc (LCOLAREA=.FALSE. et LCOLINE=.FALSE.)
!***************************************************

  CALL GSPLCI(1)

  IF(LSUPER)THEN                   !!!  Overlay case


    IF(NSUPER == 1)THEN            ! If first plot of an overlay: default 
      CALL GSLN(1)                 ! Line is solid

    ELSE                           ! If subsequent plots of an overlay: default
      
      IF(LINVPTIR)THEN

        IF(NSUPER ==2)CALL GSLN(2)    ! line is a special dash type
        IF((LARROVL .OR. (LCOLAREA .OR. LHACH1)) .AND. NSUPER ==2)CALL GSLN(1)
        IF(NSUPER ==3)CALL GSLN(3)
        IF((LARROVL .OR. (LCOLAREA .OR. LHACH1)) .AND. NSUPER ==3)THEN
	  CALL GSLN(1)
	  CALL GSLN(2)
	  IF(LHACH2)CALL GSLN(1)
        ENDIF

      ELSE

        IF(NSUPER ==2)CALL GSLN(3)    ! line is a special dash type
        IF((LARROVL .OR. (LCOLAREA .OR. LHACH1)) .AND. NSUPER ==2)CALL GSLN(1)
        IF(NSUPER ==3)CALL GSLN(2)
        IF((LARROVL .OR. (LCOLAREA .OR. LHACH1)) .AND. NSUPER ==3)THEN
	  CALL GSLN(1)
	  CALL GSLN(3)
	  IF(LHACH2)CALL GSLN(1)
        ENDIF

      ENDIF

    END IF

  END IF                           !!!  Not an overlay case
!
END IF                  !+++++++++++++++++++++++++++++++++++++++++++++++++++++
if(nverbia > 0)then
  print *,' BALISE3g IMCOU'
endif
!
!*      2.11 High and low handling
!
SELECT CASE(NHI)
    
CASE(0)                          ! H + L ara displayed
    IF(INCL /=0)THEN
      CALL CPLBDR(ZTEMV,ZRWRK,IWRK)
    ENDIF
CASE DEFAULT                     ! TO BE REVISED ********************
			         ! <0  --> no action (:-1 to be set)
			         ! >0  --> gridpoint value displayed (1: to be set)
END SELECT
    
!
!*      2.12      Effective contour drawing, perimeter box, grid and labels
!
IF((LCOLAREA .AND. .NOT.LISO .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))&
   .OR.(LHACH1 .AND. .NOT.LISO .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))&
   .OR. (LHACH2 .AND. .NOT.LISO .AND. NSUPER == 2) &
   .OR. (LHACH3 .AND. .NOT.LISO .AND. NSUPER == 3) &
   .OR. (LHACH4 .AND. .NOT.LISO .AND. NSUPER == 4))THEN
if(nverbia > 0)then
  print *,' BALISE3ha IMCOU'
endif

ELSE

if(nverbia > 0)then
  print *,' BALISE3h IMCOU XLWIDTH ',XLWIDTH
endif
  CALL GSLWSC(XLWIDTH)
if(nverbia > 0)then
  print *,' BALISE3ha IMCOU APXLWIDTH '
endif
  IF(NSUPER == 2 .AND. LISOWHI2)THEN
    CALL GSLN(1)
    CALL GSPLCI(0)
    CALL GSTXCI(0)
  ELSE IF(NSUPER == 3 .AND. LISOWHI3)THEN
    CALL GSLN(1)
    CALL GSPLCI(0)
    CALL GSTXCI(0)
  ENDIF
if(nverbia > 0)then
  print *,' BALISE3ha IMCOU AV CPCLDR '
endif
  CALL CPCLDR(ZTEMV,ZRWRK,IWRK)
if(nverbia > 0)then
  print *,' BALISE3hb IMCOU AP CPCLDR '
endif
END IF
IF((NSUPER == 2 .AND. LISOWHI2) .OR. (NSUPER == 3 .AND. LISOWHI3))THEN
! CALL GSPLCI(1)
  CALL GSTXCI(1)
ENDIF
if(nverbia > 0)then
  print *,' BALISE3I IMCOU INCL',INCL
endif
IF(INCL == 0)THEN
  CALL CPLBDR(ZTEMV,ZRWRK,IWRK)
ENDIF
IF(nverbia > 0)THEN
  print *,' **IMCOU AV GETSET'
endif
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL SETUSV('MI',1)
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
IF((LHACH1 .AND. NSUPER == 1) .OR. (LHACH2 .AND. NSUPER == 2) .OR. &
   (LHACH3 .AND. NSUPER == 3) .OR. (LHACH4 .AND. NSUPER == 4))THEN
ELSE
  IF(LSUPER .AND. NSUPER > 1)THEN
  IF((LCOLAREA .AND. NSUPER > 1) .OR. &
     (.NOT.LCOLAREA  .AND. &
      .NOT.((LHACH1.OR.LHACH2) .AND. NSUPERDIA == 2)))THEN
    ILENT=LEN_TRIM(HTEXT)+2
    IF(LPVT)THEN
      IF(NSUPERDIA >= 2 .AND. (LHACH2.OR.LHACH3))THEN
        CALL FRSTPT(.1+(NSUPER-2)*.24,ZVT+.036)
        CALL VECTOR(.1+(NSUPER-2)*.24+.03,ZVT+.036)
      ELSE
        CALL FRSTPT(.1+(NSUPER-2)*.24,ZVT+.016)
        CALL VECTOR(.1+(NSUPER-2)*.24+.03,ZVT+.016)
      ENDIF
    ELSE
      CALL GSLWSC(XLWIDTH)
      IF(NSUPERDIA >= 2 .AND. (LHACH2.OR.LHACH3))THEN
        CALL FRSTPT(.1+(NSUPER-2)*.24+ILENT*.009,ZVT+.05)
        CALL VECTOR(.1+(NSUPER-2)*.24+ILENT*.009+.03,ZVT+.05)
      ELSE
        CALL FRSTPT(.1+(NSUPER-2)*.24+ILENT*.009,ZVT+.03)
        CALL VECTOR(.1+(NSUPER-2)*.24+ILENT*.009+.03,ZVT+.03)
      ENDIF
    ENDIF
  ENDIF
  ENDIF
ENDIF

CALL SETUSV('MI',IMI)
!IF(LPRESY)THEN
! CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,2)
 if(nverbia > 0)then
  print *,' ** imcou vers FIN ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,2,ID ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
 endif
!ELSE
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!ENDIF
CALL GSLWSC(1.)
CALL GSLN(1)
CALL GSPLCI(1)
CALL GSTXCI(1)
!
CALL GSCLIP(0)
CALL GASETI('LTY',1)
! Mai 2000 Abscisses tps en heures si LHEURX=T
IF(LPVT .AND. LHEURX)THEN
! CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL/3600.,ZWR/3600.,ZWB,ZWT,ID)
  FORMAX='          '
  IF(LFMTAXEX)THEN
    FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
  ELSE
    FORMAX='(F8.0)'
  ENDIF
  FORMAY='          '
  IF(LFMTAXEY)THEN
    FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
  ELSE
    FORMAY='(F7.0)'
  ENDIF
  CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
! CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
! CALL LABMOD('(F8.0)','(F7.0)',0,0,10,10,0,0,0)
!!!!!!!Avril 2002
  IF(LMYHEURX)THEN
    ZH=NHEURXGRAD*3600.
  ELSE
!!!!!!!Avril 2002

  IF((ZWR-ZWL)/3600. > 24.)THEN
    ZH=10800.
  ELSE
    ZH=3600.
  ENDIF
!!!!!!!Avril 2002
  ENDIF
!!!!!!!Avril 2002

! CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  DO J=INT(ZWL),INT(ZWR)
    ZJ=J
!!!!!!!Avril 2002
  IF(LMYHEURX)THEN
    ZJJ=ZJ/ZH*NHEURXGRAD
    ZINTT=NHEURXLBL
  ELSE
!!!!!!!Avril 2002

      IF(ZH == 10800.)THEN
	ZJJ=ZJ/ZH*3.
	ZINTT=6.
      ELSE
	ZJJ=ZJ/ZH
	ZINTT=3.
      ENDIF
!!!!!!!Avril 2002
  ENDIF
!!!!!!!Avril 2002
!!!! Mars 2009 pour labels = hhHmm .besoin fournir les extremes sous
!!!! Mars 2009 forme reelle avec OBLIG. 2 decimales pour minutes ex 9.45
!!!! Mars 2009 pour eviter superposition ticks differents
  IF(LHEURX .AND. LAXEXUSER .AND. LNOLABELX)THEN
  ELSE
!!!! Mars 2009 pour eviter superposition ticks differents
    IF(MOD(ZJ,ZH) == 0.)THEN
      CALL FRSTPT(ZJ,ZWB)
      IF(LPRESY)THEN
        CALL VECTOR(ZJ,ZWB+(ZWT-ZWB)/60.)
      ELSE
        IF(MOD(ZJJ,ZINTT) == 0.)THEN
          CALL VECTOR(ZJ,ZWB+(ZWT-ZWB)/90.)
          if(nverbia > 0)then
            print *,' Ap VECTOR A IMCOU'
          endif
        ELSE
          CALL VECTOR(ZJ,ZWB+(ZWT-ZWB)/120.)
          if(nverbia > 0)then
            print *,' Ap VECTOR B IMCOU'
          endif
        ENDIF
      ENDIF
!!!! Mars 2009
  ENDIF
!!!! Mars 2009
 if(nverbia > 0)then
  print *,' ** imcou vers FIN ZJ ZJJ ZINT ',ZJ,ZJJ,ZINTT
 endif

      ZWBBB=ZWB-((ZWT-ZWB)/((ZVT-ZVB)/.02))
      IF(LPRESY)THEN
        ZWBBB=ZWB-((ZWT-ZWB)/((ZVT-ZVB)/.05))
      ENDIF
      IF(.NOT.LNOLABELX)THEN
      IF(MOD(ZJJ,ZINTT) == 0.)THEN
	IF(ZJJ < 10.)THEN
	  YC2='  '
	  WRITE(YC2,'(F2.0)')ZJJ
	  CALL PLCHHQ(ZJ,ZWBBB,YC2,.010,0.,0.)
	ELSEIF(ZJJ < 100.)THEN
	  YC3='   '
	  WRITE(YC3,'(F3.0)')ZJJ
	  CALL PLCHHQ(ZJ,ZWBBB,YC3,.010,0.,0.)
	ELSE
	  YC4='    '
	  WRITE(YC4,'(F4.0)')ZJJ
	  CALL PLCHHQ(ZJ,ZWBBB,YC4,.010,0.,0.)
	ENDIF
      ENDIF
      ENDIF
    ENDIF
  ENDDO
! Mars 2001
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
  IF(LFACTAXEX)THEN
    IF(LFACTAXEY)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL*XFACTAXEX,ZWRR*XFACTAXEX,&
	       ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
    ELSE
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL*XFACTAXEX,ZWRR*XFACTAXEX,&
	       ZWBB,ZWTT,IDD)
    ENDIF
  ELSEIF(LFACTAXEY)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
	       ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
  ELSEIF(LAXEXUSER)THEN
    IF(LAXEYUSER)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,XAXEXUSERD,XAXEXUSERF,&
	       XAXEYUSERD,XAXEYUSERF,IDD)
    ELSE
      CALL SET(ZVL,ZVR,ZVB,ZVT,XAXEXUSERD,XAXEXUSERF,&
	       ZWBB,ZWTT,IDD)
    ENDIF
  ELSEIF(LAXEYUSER)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
	       XAXEYUSERD,XAXEYUSERF,IDD)
  ENDIF
! Mars 2001
  IF(LPRESY)THEN
    CALL AXELOGPRES(XHMIN,XHMAX)
    CALL GRIDAL(0,0,0,0,0,0,5,0.,0.)
!   CALL GRIDAL(0,0,1,9,0,1,5,0.,0.)
 if(nverbia > 0)then
    print *,' **imcou ap GRIDAL(0,0,2,10,0,1,5,0.,0.)'
 endif
  ELSE
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0.)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
      CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0.)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0.)
    ELSE
      CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0.)
    ENDIF
!Avril 2002
  ENDIF
!!!!!!!!Mars 2009 pour ecrire des heures sous forme hhHmm sur axe X
!!!  Besoin de passer les extremes en valeurs réelles
!!! dans XAXEXUSERD et XAXEXUSERF avec 2chiffres decimaux OBLIGATOIREMENT
!!! pour les minutes
!!!  LAXEXUSER=T LHEURX=T LNOLABELX=T Obligatoires . 5 intervalles prevus
!!!  Intervenir sur NCVITVXMJ pour changer ce nb d'intervalles
  IF(LAXEXUSER .AND. LHEURX .AND. LNOLABELX)THEN
!   Conversion extremes en minutes
    ZTA=AINT(XAXEXUSERD)
    ZTB=(XAXEXUSERD-ZTA)*100.
    ZTD=ZTA*60+ZTB
    ZTA=AINT(XAXEXUSERF)
    ZTB=(XAXEXUSERF-ZTA)*100.
    ZTF=ZTA*60+ZTB
    ZTINT=(ZTF-ZTD)/NCVITVXMJ
    ALLOCATE( ZTDX(NCVITVXMJ))
    DO IA=2,NCVITVXMJ
      ZTDX(IA)=ZTD+ZTINT*(IA-1)
    ENDDO
    ZTDX(1)=ZTD
    ZINTV=(XAXEXUSERF-XAXEXUSERD)/NCVITVXMJ
       CALL GRIDAL(NCVITVXMJ,NCVITVXMN,0,0,0,1,5,0.,0.)
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWLC,ZWRC,ZWBC,ZWTC,IDD)
!   
    DO IA=1,NCVITVXMJ+1
      IF(IA == NCVITVXMJ+1)THEN
        ZTD=ZTF
      ELSE
        ZTD=ZTDX(IA)
      ENDIF
      ZTA=AINT(ZTD/60.)
      ZTB=ZTD-(ZTA*60.)
      IF(L24H)THEN
        DO IB=1,10
            if(nverbia >0)print *,' IB ',IB
          IF(ZTA > 24.)THEN
            ZTA=ZTA-24.
            if(nverbia >0)print *,' ZTA A ',ZTA
          ELSE
            IF(ZTA == 24. .AND. ZTB /= 0.)ZTA=ZTA-24.
            if(nverbia >0)print *,' ZTA B ',ZTA
!           CYCLE
          ENDIF
        ENDDO
      ENDIF
      WRITE(YFORMAT,'(I2.2,"H",I2.2)')NINT(ZTA),NINT(ZTB)     
      CALL PLCHHQ(ZWLC+(IA-1)*ZINTV,ZWBC-(ZWTC-ZWBC)/40.,YFORMAT,.01,0.,0.)
    ENDDO
    DEALLOCATE(ZTDX)
  ENDIF
!!!!!!!!Mars 2009 pour ecrire des heures  (Fin)
  IF(LFACTAXEX .OR. LFACTAXEY .OR. LAXEXUSER .OR. LAXEYUSER)THEN
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
  ENDIF

ELSE

  FORMAX='          '
  IF(LFMTAXEX)THEN
    FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
  ELSE
    FORMAX='(F8.0)'
  ENDIF
  FORMAY='          '
  IF(LFMTAXEY)THEN
    FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
  ELSE
    FORMAY='(F8.0)'
  ENDIF

  IF(ABS(ZWB) > 999999. .OR. ABS(ZWT) > 999999.)THEN
    CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!   CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
! CALL LABMOD('(F8.0)','(F8.0)',0,0,10,10,0,0,0)
  ELSE
    FORMAY='          '
    IF(LFMTAXEY)THEN
      FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
    ELSE
      FORMAY='(F8.0)'
!     FORMAY='(F7.0)'
    ENDIF
    CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!   CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
! CALL LABMOD('(F8.0)','(F7.0)',0,0,10,10,0,0,0)
  ENDIF
  IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1) .OR. &
  (NSUPER == 2 .AND. LISOWHI2) .OR. (NSUPER == 3 .AND. LISOWHI3))THEN
! Mars 2001
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
  IF(LFACTAXEX)THEN
    IF(LFACTAXEY)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL*XFACTAXEX,ZWRR*XFACTAXEX,&
	       ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
    ELSE
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL*XFACTAXEX,ZWRR*XFACTAXEX,&
	       ZWBB,ZWTT,IDD)
    ENDIF
  ELSEIF(LFACTAXEY)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
	       ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
  ELSEIF(LAXEXUSER)THEN
    IF(LAXEYUSER)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,XAXEXUSERD,XAXEXUSERF,&
	       XAXEYUSERD,XAXEYUSERF,IDD)
    ELSE
      CALL SET(ZVL,ZVR,ZVB,ZVT,XAXEXUSERD,XAXEXUSERF,&
	       ZWBB,ZWTT,IDD)
    ENDIF
  ELSEIF(LAXEYUSER)THEN
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
	       XAXEYUSERD,XAXEYUSERF,IDD)
  ENDIF
! Mars 2001
  IF(LPRESY)THEN
    CALL AXELOGPRES(XHMIN,XHMAX)
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,0,0,0,0,5,0.,0.)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,0,0,0,0,5,0.,0.)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,0,0,1,0,5,0.,0.)
    ELSE
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,0,0,1,0,5,0.,0.)
    ENDIF
!Avril 2002
  ELSE
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0.)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0.)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,1,0,5,0.,0.)
    ELSE
      CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,1,1,5,0.,0.)
    ENDIF
!Avril 2002
  ENDIF
! CALL GRIDAL(5,0,10,0,1,1,5,0.,0.)
  IF(LFACTAXEX .OR. LFACTAXEY .OR. LAXEXUSER .OR. LAXEYUSER)THEN
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
  ENDIF
  ENDIF
ENDIF
!
    IF(.NOT.LDEFCV2CC)THEN              !%%%%%%%%%%%%%%%%%%%%%%

    IF(NLANGLE.EQ.0.AND.XIDEBCOU.EQ.-999..AND.LXZ)THEN
      CALL GSCLIP(0)
      CALL TRACEXZ
      CALL GSCLIP(1)
    END IF

    ENDIF                               !%%%%%%%%%%%%%%%%%%%%%%
!
!*      2.13      General NCAR parameter reset
!
CALL CPSETI('CLS',16)
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == NSUPERDIA))THEN
  CALL CPRSET
ENDIF
CALL GSLN(1)
!
!*      2.14      Final touch: page information labels
!
IF(nverbia > 0)THEN
  print *,' **IMCOU AV GETSET 2'
endif
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
IF(LANIMT)THEN
  CALL PLCHHQ((0.002-ZVL)*(ZWR-ZWL)/(ZVR-ZVL),(0.050-ZVB)*(ZWT-ZWB)/(ZVT-ZVB), &
	      CTIMEC,.009,0.,-1.)
ENDIF
CALL SETUSV('MI',1)
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!IF(LFACTIMP)THEN
! CALL FACTIMP
!ENDIF
!
!!!!!!!Debut Titres pour NSUPER = 1!!!!!!!!!!!!!!!!!!!
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN
! Mars 2000
IF(LFACTIMP)THEN
  CALL FACTIMP
ENDIF
! Modifs for diachro
!
  IF(LXYO )THEN
! IF(LXYO .AND. XIDEBCOU == -999. .AND. XJDEBCOU == -999.)THEN
    YXYO(1:LEN(YXYO))=' '
    IO=1
    YXYO(IO:IO)='('
!   ZX=XXX(NIDEBCOU,NMGRID)
!   ZY=XXY(NJDEBCOU,NMGRID)
    ZX=XDSX(1,NMGRID)
    ZY=XDSY(1,NMGRID)
    ZXE=XDSX(NLMAX,NMGRID)
    ZYE=XDSY(NLMAX,NMGRID)
    YC8(1:LEN(YC8))=' '
    WRITE(YC8,'(F8.0)')ZX
    YC8=ADJUSTL(YC8)
    IO=IO+1
    YXYO(IO:IO+LEN_TRIM(YC8)-1)=YC8(1:LEN_TRIM(YC8))
    IO=IO+LEN_TRIM(YC8)
    YXYO(IO:IO)=','
    IO=IO+1
    YC8(1:LEN(YC8))=' '
    WRITE(YC8,'(F8.0)')ZY
    YC8=ADJUSTL(YC8)
    IO=IO+1
    YXYO(IO:IO+LEN_TRIM(YC8)-1)=YC8(1:LEN_TRIM(YC8))
    IO=IO+LEN_TRIM(YC8)
    YXYO(IO:IO)=')'
    CALL PLCHHQ(ZVL-.009,ZVB-(ZVB/7.1),YXYO(1:LEN_TRIM(YXYO)),.007,0.,-1.)
    YXYO(1:LEN(YXYO))=' '
    IO=1
    YXYO(IO:IO)='('
    YC8(1:LEN(YC8))=' '
    WRITE(YC8,'(F8.0)')ZXE
    YC8=ADJUSTL(YC8)
    IO=IO+1
    YXYO(IO:IO+LEN_TRIM(YC8)-1)=YC8(1:LEN_TRIM(YC8))
    IO=IO+LEN_TRIM(YC8)
    YXYO(IO:IO)=','
    IO=IO+1
    YC8(1:LEN(YC8))=' '
    WRITE(YC8,'(F8.0)')ZYE
    YC8=ADJUSTL(YC8)
    IO=IO+1
    YXYO(IO:IO+LEN_TRIM(YC8)-1)=YC8(1:LEN_TRIM(YC8))
    IO=IO+LEN_TRIM(YC8)
    YXYO(IO:IO)=')'
    CALL PLCHHQ(ZVR,ZVB-(ZVB/7.1),YXYO(1:LEN_TRIM(YXYO)),.007,0.,+1.)
  ENDIF
! Remodifs le 17/05/96
!
! Titres en X
!
if(nverbia > 0)then
  print *,' BALISE4 IMCOU NLMAX',NLMAX
endif
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXL',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXL',YTEM)
    IF(XSZTITXL /= 0.)THEN
      CALL PLCHHQ(ZVL,ZVB-MIN(ZVB/2.,.05),YTEM,XSZTITXL,0.,-1.)
!     CALL PLCHHQ(ZVL,ZVB/2.,YTEM,XSZTITXL,0.,-1.)
    ELSE
      CALL PLCHHQ(ZVL,ZVB-MIN(ZVB/2.,.05),YTEM,.008,0.,-1.)
!     CALL PLCHHQ(ZVL,ZVB/2.,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXM',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXM',YTEM)
    IF(XSZTITXM /= 0.)THEN
      CALL PLCHHQ((ZVL+ZVR)/2.,ZVB-MIN(ZVB/2.,.05),YTEM(1:LEN_TRIM(YTEM)),XSZTITXM,0.,0.)
!     CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),XSZTITXM,0.,0.)
    ELSE
      CALL PLCHHQ((ZVL+ZVR)/2.,ZVB-MIN(ZVB/2.,.05),YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
!     CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
    ENDIF
  ENDIF
!
! Titres en Y
!
IF(nverbia > 0)THEN
  print *,' **IMCOU AV TITRES Y'
endif
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYM',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYB',ZVL,ZVR,ZVB,ZVT,YTEM)
  ZXPOSTITB1=.002
  ZXYPOSTITB1=.005
  IF(XPOSTITB1 /= 0.)THEN
    ZXPOSTITB1=XPOSTITB1
  ENDIF
  IF(XYPOSTITB1 /= 0.)THEN
    ZXYPOSTITB1=XYPOSTITB1
  ENDIF
  CALL RESOLV_TIT('CTITB1',HLEGEND)
  IF(XSZTITB1 /= 0.)THEN
    CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HLEGEND,XSZTITB1,0.,-1.)
!   CALL PLCHHQ(0.002,0.005,HLEGEND,XSZTITB1,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HLEGEND,.007,0.,-1.)
!   CALL PLCHHQ(0.002,0.005,HLEGEND,.007,0.,-1.)
  ENDIF
!
! Titres TOP 
!
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITT2',YTEM)

  ZXPOSTITT2=.002
  ZXYPOSTITT2=.95
  IF(XPOSTITT2 /= 0.)THEN
    ZXPOSTITT2=XPOSTITT2
  ENDIF
  IF(XYPOSTITT2 /= 0.)THEN
    ZXYPOSTITT2=XYPOSTITT2
  ENDIF
!! Oct 2001
  IF(YTEM == ' ' .OR. YTEM == 'DEFAULT')THEN
!!!Mars 2009 + NPROFILE /= 0
    IF(LPVT .AND. NPROFILE /= 1 .AND. NPROFILE /= 0)THEN          
      YTEM(1:LEN(YTEM))=' '
      WRITE(YTEM,1024)NPROFILE
    ENDIF
  ENDIF

  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    IF(XSZTITT2 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,XSZTITT2,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YTEM,XSZTITT2,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,.008,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  ZXPOSTITT3=.002
  ZXYPOSTITT3=.93
  IF(XPOSTITT3 /= 0.)THEN
    ZXPOSTITT3=XPOSTITT3
  ENDIF
  IF(XYPOSTITT3 /= 0.)THEN
    ZXYPOSTITT3=XYPOSTITT3
  ENDIF
  CALL RESOLV_TIT('CTITT3',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    IF(XSZTITT3 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM(1:LEN_TRIM(YTEM)),XSZTITT3,0.,-1.)
!     CALL PLCHHQ(0.002,0.93,YTEM(1:LEN_TRIM(YTEM)),XSZTITT3,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM(1:LEN_TRIM(YTEM)),.008,0.,-1.)
!     CALL PLCHHQ(0.002,0.93,YTEM(1:LEN_TRIM(YTEM)),.008,0.,-1.)
    ENDIF
  ENDIF
  
  YCARCOU(1:LEN(YCARCOU))=' '
  YCAR(1:LEN(YCAR))=' '

  IF(.NOT.LPVT .AND..NOT.LPXT .AND..NOT.LPYT)THEN
  
    YTEM(1:LEN(YTEM))=' '
    IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
    ENDIF
    IF(.NOT.LANIMT)THEN
  ! YTEM(1:LEN(YTEM))=' '
    ZXPOSTITB3=.002
    ZXYPOSTITB3=.045
    IF(XPOSTITB3 /= 0.)THEN
      ZXPOSTITB3=XPOSTITB3
    ENDIF
    IF(XYPOSTITB3 /= 0.)THEN
      ZXYPOSTITB3=XYPOSTITB3
    ENDIF

      IF(LMINUS .OR. LPLUS)THEN

        IF(.NOT.LTITDEFM .AND. CTITB3MEM /= 'DEFAULT' .AND. &
	   CTITB3MEM /= 'default' .AND. CTITB3MEM /= 'DEFAUT' .AND. &
	   CTITB3MEM /= 'defaut')THEN
	  IF(CTITB3MEM /= ' ' .AND. CTITB3MEM /= 'WHITE' .AND. &
	     CTITB3MEM /= 'white' .AND. CTITB3MEM /= 'BLANC' .AND. &
	     CTITB3MEM /= 'blanc')THEN
	    IF(XSZTITB3 /= 0.)THEN
	      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3MEM(1:LEN_TRIM(CTITB3MEM)),XSZTITB3,0.,-1.)
	    ELSE
	      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3MEM(1:LEN_TRIM(CTITB3MEM)),.009,0.,-1.)
	    ENDIF
	  ENDIF

        ELSE

        CALL RESOLV_TIT('CTITB3',CTITB3)
        IF(CTITB3 /= ' ')THEN
          IF(XSZTITB3 /= 0.)THEN
            CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3(1:LEN_TRIM(CTITB3)),XSZTITB3,0.,-1.)
!           CALL PLCHHQ(0.002,0.050,CTITB3(1:LEN_TRIM(CTITB3)),XSZTITB3,0.,-1.)
          ELSE
            CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3(1:LEN_TRIM(CTITB3)),.009,0.,-1.)
!           CALL PLCHHQ(0.002,0.050,CTITB3(1:LEN_TRIM(CTITB3)),.009,0.,-1.)
	  ENDIF
        ENDIF

        ENDIF

      ELSE
        if(nverbia > 0)then
        print *,' **imcou CTIMEC,YTEM ',CTIMEC,YTEM
        endif
        YTEM(1:LEN(YTEM))=' '
	YTEM=CTIMEC
	YTEM=ADJUSTL(YTEM)
        if(nverbia > 0)then
          print *,' **imcou CTIMEC,YTEM ',CTIMEC,YTEM
        endif
        CALL RESOLV_TIT('CTITB3',YTEM)
!       CALL RESOLV_TIT('CTITB3',CTIMEC)
        IF(YTEM /= ' ')THEN
          IF(XSZTITB3 /= 0.)THEN
            CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM(1:LEN_TRIM(YTEM)),XSZTITB3,0.,-1.)
          ELSE
            CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM(1:LEN_TRIM(YTEM)),.009,0.,-1.)
          ENDIF
        ENDIF
        if(nverbia > 0)then
        print *,' **imcou CTIMEC,YTEM ',CTIMEC,YTEM
        endif
!        IF(CTIMEC /= ' ')THEN
!          IF(XSZTITB3 /= 0.)THEN
!            CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTIMEC,XSZTITB3,0.,-1.)
!          ELSE
!            CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTIMEC,.009,0.,-1.)
!          ENDIF
!        ENDIF
      ENDIF
    ENDIF
    CALL RESOLV_TIT('CTITB2',CLEGEND2)
    CLEGEND2=ADJUSTL(CLEGEND2)
    if(nverbia > 0)then
      print *,' **imcou CLEGEND2 ',CLEGEND2(1:LEN_TRIM(CLEGEND2))
      print *,' **imcou CTITB2 ',CTITB2(1:LEN_TRIM(CTITB2))
    endif
    ZXPOSTITB2=.002
    ZXYPOSTITB2=.025
    IF(XPOSTITB2 /= 0.)THEN
      ZXPOSTITB2=XPOSTITB2
    ENDIF
    IF(XYPOSTITB2 /= 0.)THEN
      ZXYPOSTITB2=XYPOSTITB2
    ENDIF
    IF(CLEGEND2 /= ' ')THEN
      IF(XSZTITB2 /= 0.)THEN
        CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,XSZTITB2,0.,-1.)
!       CALL PLCHHQ(0.002,0.025,CLEGEND2,XSZTITB2,0.,-1.)
      ELSE
        CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,.007,0.,-1.)
!       CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
      ENDIF
    ENDIF
    IF(XIDEBCOU.NE.-999.)THEN

      IF(LDEFCV2CC)THEN           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IF(LDEFCV2IND)THEN
	  WRITE(YCARCOU,1018)NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
	ELSE IF(LDEFCV2LL)THEN
	  WRITE(YCARCOU,1019)XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
	ELSE
	  WRITE(YCARCOU,1020)XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
	ENDIF
      ELSE                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(XIDEBCOU < 99999.)THEN
        IF(XJDEBCOU < 99999.)THEN
          WRITE(YCARCOU,1001)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
        ELSE
          WRITE(YCARCOU,1002)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
        END IF
      ELSE
        IF(XJDEBCOU < 99999.)THEN
          WRITE(YCARCOU,1003)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
        ELSE
          WRITE(YCARCOU,1004)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
        END IF
      END IF
      ENDIF                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELSE
      WRITE(YCARCOU,1000)NIDEBCOU,NJDEBCOU,NLANGLE,NLMAX
    END IF
    if(nverbia > 0)then
    print *,' IMCOU AV RESOLVTIT 1 ',YCARCOU(1:LEN_TRIM(YCARCOU))
    endif
    CALL RESOLV_TIT('CTITT1',YCARCOU)
    if(nverbia > 0)then
    print *,' IMCOU AP RESOLVTIT 1'
    endif
    ZXPOSTITT1=.002
    ZXYPOSTITT1=.98
    IF(XPOSTITT1 /= 0.)THEN
      ZXPOSTITT1=XPOSTITT1
    ENDIF
    IF(XYPOSTITT1 /= 0.)THEN
      ZXYPOSTITT1=XYPOSTITT1
    ENDIF
    IF(YCARCOU /= ' ')THEN
      IF(XSZTITT1 /= 0.)THEN
        CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!       CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
      ELSE
        CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.012,0.,-1.)
!       CALL PLCHHQ(0.002,0.98,YCARCOU,.012,0.,-1.)
      ENDIF
    ENDIF
    YTEM(1:LEN(YTEM))=' '
    CALL RESOLV_TIT('CTITXR',YTEM)
    IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
      CALL RESOLV_TIT('CTITXR',YTEM)
      IF(XSZTITXR /= 0.)THEN
        CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/2.,.05),YTEM,XSZTITXR,0.,-1.)
!       CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,XSZTITXR,0.,-1.)
      ELSE
        CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/2.,.05),YTEM,.008,0.,-1.)
!       CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
      ENDIF
    ENDIF

  ELSE

    IND=INDEX(CLEGEND2(10:LEN_TRIM(CLEGEND2)),'DATE')
    IF(IND == 0)THEN
      CLEGEND2(1:LEN_TRIM(CLEGEND2))=' '
    ELSE
      IND=IND+10-1
      CLEGEND2(IND:LEN_TRIM(CLEGEND2))=' '
    ENDIF
    CALL RESOLV_TIT('CTITB2',CLEGEND2)
    ZXPOSTITB2=.002
    ZXYPOSTITB2=.025
    IF(XPOSTITB2 /= 0.)THEN
      ZXPOSTITB2=XPOSTITB2
    ENDIF
    IF(XYPOSTITB2 /= 0.)THEN
      ZXYPOSTITB2=XYPOSTITB2
    ENDIF
    IF(CLEGEND2 /= ' ')THEN
      IF(XSZTITB2 /= 0.)THEN
        CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,XSZTITB2,0.,-1.)
!       CALL PLCHHQ(0.002,0.025,CLEGEND2,XSZTITB2,0.,-1.)
      ELSE
        CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,.007,0.,-1.)
!       CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
      ENDIF
    ENDIF
    YTEM(1:LEN(YTEM))=' '
    IF(LPVT .OR. LPYT .OR. (LPXT .AND..NOT.LXABSC))THEN
      IF(LPVT .AND. LHEURX)THEN
        YTEM='(H)'
      ELSE
        YTEM='(Sec.)'
      ENDIF
    ELSE IF(LPXT .AND. LXABSC)THEN
      YTEM='(X)'
    ENDIF
    CALL RESOLV_TIT('CTITXR',YTEM)
    IF(YTEM /= ' ')THEN
      IF(XSZTITXR /= 0.)THEN
        CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/2.,.05),YTEM,XSZTITXR,0.,-1.)
!       CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,XSZTITXR,0.,-1.)
      ELSE
        CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/2.,.05),YTEM,.008,0.,-1.)
!       CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
      ENDIF
    ENDIF
    YTEM(1:LEN(YTEM))=' '
    SELECT CASE(CTYPE)
      CASE('CART')
        IF(LPXT .AND..NOT.LXABSC)THEN
          YTEM='(X)'
        ELSE IF(LPXT .AND. LXABSC)THEN
          YTEM='(S)'
        ELSE IF(LPYT)THEN
          YTEM='(Y)'
        ELSE
IF(nverbia > 0)THEN
  print *,' **IMCOU AV Model;Levels;(M)'
endif
          YTEM='Model;Levels;(M)'
          IF(LPRESY)THEN
            YTEM='Pressure;(Mbs)'
          ENDIF
        ENDIF
      CASE DEFAULT
IF(nverbia > 0)THEN
  print *,' **IMCOU AV Levels;(M)'
endif
        YTEM='Levels;(M)'
    END SELECT
    CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
    !CALL PLCHHQ(0.,ZVT-1*.015,'Model',.008,0.,-1.)
    !CALL PLCHHQ(0.,ZVT-2*.015,'Levels',.008,0.,-1.)
    !CALL PLCHHQ(0.,ZVT-3*.015,'(M)',.008,0.,-1.)

    IF(L1DT)THEN
      SELECT CASE(CTYPE)
        CASE('CART')
          IF(LPXT)THEN
            WRITE(YCARCOU,1016)NIINF,NISUP
          ELSE IF(LPYT)THEN
            WRITE(YCARCOU,1017)NJINF,NJSUP
          ELSE
            WRITE(YCARCOU,1012)
          ENDIF
        CASE('SSOL')
  	YCARCOU(1:LEN(YCARCOU))=' '
  	YCARCOU(1:7)='SSOL N.'
  	WRITE(YCARCOU(8:10),'(I3)')NLOOPN
  	YCARCOU(11:13)='  ('
  	WRITE(YCARCOU(14:18),'(F5.0)')XTRAJX(1,1,NLOOPN)
  	YCARCOU(19:19)=','
  	WRITE(YCARCOU(20:24),'(F5.0)')XTRAJY(1,1,NLOOPN)
  	YCARCOU(25:27)=')  '
  	ISUIT=28
  	ISUI=8
  	IF(ALLOCATED(ISTM))THEN
  	  DEALLOCATE(ISTM)
          ENDIF
  	  ALLOCATE(ISTM(NSUPERDIA))
! 20 Nov 2000
  	INDISTM=1
! 20 Nov 2000
  	ISTM(INDISTM)=NLOOPN
        CASE DEFAULT
  	YCARCOU(1:LEN(YCARCOU))=' '
  	YCARCOU(1:4)=CTYPE
  	YCARCOU(5:7)=' N.'
  	WRITE(YCARCOU(8:10),'(I3)')NLOOPN
        if(nverbia > 0)then
        print *,' ** IMCOU YCARCOU AP WRI NLOOPN ',YCARCOU
        endif
  	ISUIT=11
  	IF(ALLOCATED(ISTM))THEN
  	  DEALLOCATE(ISTM)
          ENDIF
  	  ALLOCATE(ISTM(NSUPERDIA))
        if(nverbia > 0)then
        print *,' ** IMCOU NSUPERDIA ISTM ',NSUPERDIA
        endif
  	INDISTM=1
  	ISTM(INDISTM)=NLOOPN
        if(nverbia > 0)then
        print *,' ** IMCOU NSUPERDIA ISTM ',NSUPERDIA,ISTM
        endif
      END SELECT
IF(nverbia > 0)THEN
  print *,' **IMCOU FIN IF(L1DT)'
endif

    ELSE
IF(nverbia > 0)THEN
  print *,' **IMCOU FIN IF(L1DT) et AP ELSE'
endif

      IF(LPXT)THEN
        WRITE(YCARCOU,1016)NIINF,NISUP
      ELSE IF(LPYT)THEN
        WRITE(YCARCOU,1017)NJINF,NJSUP
      ELSE
        IF(XIDEBCOU.NE.-999.)THEN
          IF(LDEFCV2CC)THEN           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            YCAR(1:LEN(YCAR))=' '
	    IF(LDEFCV2IND)THEN
              IF(LPVT .AND. NPROFILE == 1)THEN
	        WRITE(YCARCOU,1023)NIDEBCV,NJDEBCV
              ELSE
                IF(LPVT .AND. NPROFILE /= 1)THEN
	          WRITE(YCARCOU,1018)NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
	          WRITE(YCAR,1024)NPROFILE
                ELSE
	          WRITE(YCARCOU,1018)NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
                ENDIF
              ENDIF
	    ELSE IF(LDEFCV2LL)THEN
              IF(LPVT .AND. NPROFILE == 1)THEN
	        WRITE(YCARCOU,1021)XIDEBCVLL,XJDEBCVLL
              ELSE
                IF(LPVT .AND. NPROFILE /= 1)THEN
	          WRITE(YCARCOU,1019)XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
	          WRITE(YCAR,1024)NPROFILE
                ELSE
	          WRITE(YCARCOU,1019)XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
                ENDIF
              ENDIF
	    ELSE
              IF(LPVT .AND. NPROFILE == 1)THEN
	        WRITE(YCARCOU,1022)XIDEBCV,XJDEBCV
              ELSE
                IF(LPVT .AND. NPROFILE /= 1)THEN
	          WRITE(YCARCOU,1020)XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
	          WRITE(YCAR,1024)NPROFILE
                ELSE
	          WRITE(YCARCOU,1020)XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
                ENDIF
              ENDIF
	    ENDIF
          ELSE                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          IF(XIDEBCOU < 99999.)THEN
            IF(XJDEBCOU < 99999.)THEN
              WRITE(YCARCOU,1011)XIDEBCOU,XJDEBCOU,NLANGLE,NPROFILE
            ELSE
              WRITE(YCARCOU,1013)XIDEBCOU,XJDEBCOU,NLANGLE,NPROFILE
            END IF
          ELSE
            IF(XJDEBCOU < 99999.)THEN
              WRITE(YCARCOU,1014)XIDEBCOU,XJDEBCOU,NLANGLE,NPROFILE
            ELSE
              WRITE(YCARCOU,1015)XIDEBCOU,XJDEBCOU,NLANGLE,NPROFILE
            END IF
          END IF
          ENDIF                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ELSE
          WRITE(YCARCOU,1010)NIDEBCOU,NJDEBCOU,NLANGLE,NPROFILE
        ENDIF
      ENDIF

    END IF

    if(nverbia > 0)then
    print *,' IMCOU AV RESOLVTIT ',YCARCOU(1:LEN_TRIM(YCARCOU))
    endif
    CALL RESOLV_TIT('CTITT1',YCARCOU)
    if(nverbia > 0)then
    print *,' IMCOU AP RESOLVTIT '
    endif
    ZXPOSTITT1=.002
    ZXYPOSTITT1=.98
    IF(XPOSTITT1 /= 0.)THEN
      ZXPOSTITT1=XPOSTITT1
    ENDIF
    IF(XYPOSTITT1 /= 0.)THEN
      ZXYPOSTITT1=XYPOSTITT1
    ENDIF
    IF(YCARCOU /= ' ')THEN
      IF(LSUPER)THEN
        SELECT CASE(CTYPE)
  	CASE ('CART','MASK')
          IF(XSZTITT1 /= 0.)THEN
            CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!           CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
          ELSE
            CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.009,0.,-1.)
!           CALL PLCHHQ(0.002,0.98,YCARCOU,.009,0.,-1.)
          ENDIF
  	CASE DEFAULT
        END SELECT
      ELSE
        IF(XSZTITT1 /= 0.)THEN
          CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!         CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
        ELSE
          CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.012,0.,-1.)
!         CALL PLCHHQ(0.002,0.98,YCARCOU,.012,0.,-1.)
	ENDIF
      ENDIF
    ENDIF

  ENDIF    ! Fin .NOT.LPVT

  IF(.NOT.LPVT .AND..NOT.LPXT .AND..NOT.LPYT)THEN
    IF(LDATFILE)CALL DATFILE_FORDIACHRO
  ENDIF

ENDIF    ! Fin .NOT.SUPER  .OR.  (LSUPER ...
!!!!!!!Fin  Titres pour NSUPER = 1!!!!!!!!!!!!!!!!!!!

IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN

  IF(NLOOPSUPER == 1)THEN
    CALL RESOLV_TIT('CTITVAR1',HTEXT)
  ELSE IF(NLOOPSUPER == 2)THEN
    CALL RESOLV_TIT('CTITVAR2',HTEXT)
  ELSE IF(NLOOPSUPER == 3)THEN
    CALL RESOLV_TIT('CTITVAR3',HTEXT)
  ELSE IF(NLOOPSUPER == 4)THEN
    CALL RESOLV_TIT('CTITVAR4',HTEXT)
  ELSE IF(NLOOPSUPER == 5)THEN
    CALL RESOLV_TIT('CTITVAR5',HTEXT)
  ELSE IF(NLOOPSUPER == 6)THEN
    CALL RESOLV_TIT('CTITVAR6',HTEXT)
  ELSE IF(NLOOPSUPER == 7)THEN
    CALL RESOLV_TIT('CTITVAR7',HTEXT)
  ELSE IF(NLOOPSUPER == 8)THEN
    CALL RESOLV_TIT('CTITVAR8',HTEXT)
  ENDIF

  if(nverbia > 0)then
    print *,' ** IMCOU HTEXT LENTRIM(HTEXT) ',LEN_TRIM(HTEXT),' ',&
    HTEXT(1:LEN_TRIM(HTEXT)),' NLOOPSUPER ',NLOOPSUPER,' NSUPER ',NSUPER,&
    ' NSUPERDIA ',NSUPERDIA
    print *,' XSZTITVAR1 ',XSZTITVAR1
  endif
  IF(HTEXT /= ' ')THEN
      ZSC=.009
      IF(XSZTITVAR1 /= 0. .AND. NLOOPSUPER == 1)ZSC=XSZTITVAR1
      IF(XSZTITVAR2 /= 0. .AND. NLOOPSUPER == 2)ZSC=XSZTITVAR2
      IF(XSZTITVAR3 /= 0. .AND. NLOOPSUPER == 3)ZSC=XSZTITVAR3
      IF(XSZTITVAR4 /= 0. .AND. NLOOPSUPER == 4)ZSC=XSZTITVAR4
      IF(XSZTITVAR5 /= 0. .AND. NLOOPSUPER == 5)ZSC=XSZTITVAR5
      IF(XSZTITVAR6 /= 0. .AND. NLOOPSUPER == 6)ZSC=XSZTITVAR6
      IF(XSZTITVAR7 /= 0. .AND. NLOOPSUPER == 7)ZSC=XSZTITVAR7
      IF(XSZTITVAR8 /= 0. .AND. NLOOPSUPER == 8)ZSC=XSZTITVAR8
  if(nverbia > 0)then
    print *,' ZSC ',ZSC
  endif
    CALL PLCHHQ(MAX(ZVR,.99),0.007,HTEXT,ZSC,0.,+1.)
  ENDIF

  IF(LMINMAX)THEN
    CALL PCSETC('FC','/')
    IF(NSUPERDIA == 1)THEN
        CAll PLCHHQ(ZVR,ZVT+.03,YLBL,.011,0.,+1.)
    ELSE
      CAll PLCHHQ(.98,ZVT+.01+(NSUPER-1)*.02,YLBL,.007,0.,+1.)
    ENDIF
    CALL PCSETC('FC',':')
  ENDIF

ELSE

  SELECT CASE(CTYPE)
    CASE('SSOL','DRST','RSPL','RAPL')
  if(nverbia > 0)then
    print *,' ** IMCOU AP CASE SSOL INDISTM ... ',INDISTM
  endif
      WRITE(YTEX(1:4),'(I4)')NLOOPN
      YTEX(1+5:LEN_TRIM(HTEXT)+5)=HTEXT(1:LEN_TRIM(HTEXT))
      YTEX=ADJUSTL(ADJUSTR(YTEX))
      IF(NSUPER > 1)THEN
	ISTOK=0
	DO JB=1,INDISTM
	  IF(NLOOPN == ISTM(JB))THEN
	    ISTOK=1
	  ENDIF
	ENDDO
	IF(ISTOK == 1)THEN
	ELSE
	  INDISTM=INDISTM+1
	  ISTM(INDISTM)=NLOOPN
	  IF(CTYPE == 'SSOL')THEN
	    IF(ISUIT > 50)THEN
	      WRITE(YCAR(ISUI:ISUI+3),'(I4)')NLOOPN
	      YCAR(ISUI+4:ISUI+6)='  ('
	      WRITE(YCAR(ISUI+7:ISUI+11),'(F5.0)')XTRAJX(1,1,NLOOPN)
	      ISUI=ISUI+12
	      YCAR(ISUI:ISUI)=','
	      ISUI=ISUI+1
	      WRITE(YCAR(ISUI:ISUI+4),'(F5.0)')XTRAJY(1,1,NLOOPN)
	      ISUI=ISUI+5
	      YCAR(ISUI:ISUI+2)=')  '
	      ISUI=ISUI+3
	    ELSE
	      WRITE(YCARCOU(ISUIT:ISUIT+3),'(I4)')NLOOPN
	      YCARCOU(ISUIT+4:ISUIT+6)='  ('
	      WRITE(YCARCOU(ISUIT+7:ISUIT+11),'(F5.0)')XTRAJX(1,1,NLOOPN)
	      ISUIT=ISUIT+12
	      YCARCOU(ISUIT:ISUIT)=','
	      ISUIT=ISUIT+1
	      WRITE(YCARCOU(ISUIT:ISUIT+4),'(F5.0)')XTRAJY(1,1,NLOOPN)
	      ISUIT=ISUIT+5
	      YCARCOU(ISUIT:ISUIT+2)=')  '
	      ISUIT=ISUIT+3
	    ENDIF
	  ELSE
	    WRITE(YCARCOU(ISUIT:ISUIT+4),'(I5)')NLOOPN
	    ISUIT=ISUIT+5
	  ENDIF
	ENDIF
      ENDIF
      IF(NSUPER == NSUPERDIA)THEN
        CALL RESOLV_TIT('CTITT1',YCARCOU)
        ZXPOSTITT1=.002
        ZXYPOSTITT1=.98
        IF(XPOSTITT1 /= 0.)THEN
          ZXPOSTITT1=XPOSTITT1
        ENDIF
        IF(XYPOSTITT1 /= 0.)THEN
          ZXYPOSTITT1=XYPOSTITT1
        ENDIF
        IF(YCARCOU /= ' ')THEN
          IF(XSZTITT1 /= 0.)THEN
            CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!           CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
          ELSE
            CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.009,0.,-1.)
!           CALL PLCHHQ(0.002,0.98,YCARCOU,.009,0.,-1.)
	  ENDIF
	  CALL RESOLV_TIT('CTITT2',YCAR)
           ZXPOSTITT2=.002
           ZXYPOSTITT2=.95
           IF(XPOSTITT2 /= 0.)THEN
             ZXPOSTITT2=XPOSTITT2
           ENDIF
           IF(XYPOSTITT2 /= 0.)THEN
             ZXYPOSTITT2=XYPOSTITT2
           ENDIF
	  IF(YCAR /= ' ')THEN
            IF(XSZTITT2 /= 0.)THEN
              CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,XSZTITT2,0.,-1.)
!             CALL PLCHHQ(0.002,0.95,YCAR,XSZTITT2,0.,-1.)
	    ELSE
              CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,.009,0.,-1.)
!             CALL PLCHHQ(0.002,0.95,YCAR,.009,0.,-1.)
	    ENDIF
	  ENDIF
	ENDIF
	IF(ALLOCATED(ISTM))THEN
	  DEALLOCATE(ISTM)
	ENDIF
      ENDIF
    CASE DEFAULT
      YTEX=ADJUSTL(HTEXT)
  END SELECT

  IF(NLOOPSUPER == 1)THEN
    CALL RESOLV_TIT('CTITVAR1',YTEX)
  ELSE IF(NLOOPSUPER == 2)THEN
    CALL RESOLV_TIT('CTITVAR2',YTEX)
  ELSE IF(NLOOPSUPER == 3)THEN
    CALL RESOLV_TIT('CTITVAR3',YTEX)
  ELSE IF(NLOOPSUPER == 4)THEN
    CALL RESOLV_TIT('CTITVAR4',YTEX)
  ELSE IF(NLOOPSUPER == 5)THEN
    CALL RESOLV_TIT('CTITVAR5',YTEX)
  ELSE IF(NLOOPSUPER == 6)THEN
    CALL RESOLV_TIT('CTITVAR6',YTEX)
  ELSE IF(NLOOPSUPER == 7)THEN
    CALL RESOLV_TIT('CTITVAR7',YTEX)
  ELSE IF(NLOOPSUPER == 8)THEN
    CALL RESOLV_TIT('CTITVAR8',YTEX)
  ENDIF

  IF(YTEX /= ' ')THEN !************************************************

  IF(LEN_TRIM(YTEX) > 25)THEN
    IF(NSUPERDIA >= 2 .AND. (LHACH2.OR.LHACH3))THEN
      CALL PLCHHQ(0.1+(NSUPER-2)*.24,ZVT+0.05,YTEX(1:LEN_TRIM(YTEX)),.005,0.,-1.)
    ELSE
      CALL PLCHHQ(0.1+(NSUPER-2)*.24,ZVT+0.03,YTEX(1:LEN_TRIM(YTEX)),.005,0.,-1.)
    ENDIF
  ELSE
    IF(NSUPERDIA >= 2 .AND. (LHACH2.OR.LHACH3))THEN
      CALL PLCHHQ(0.1+(NSUPER-2)*.24,ZVT+0.05,YTEX(1:LEN_TRIM(YTEX)),.005,0.,-1.)
    ELSE
!   CALL PLCHHQ(0.1+(NSUPER-1)*.24,ZVT+0.03,HTEXT,.007,0.,-1.)
      CALL PLCHHQ(0.1+(NSUPER-2)*.24,ZVT+0.03,YTEX(1:LEN_TRIM(YTEX)),.007,0.,-1.)
    ENDIF
  ENDIF
! CALL PLCHHQ(0.1+(NSUPER-1)*.24,ZVT+0.03,HTEXT,.009,0.,-1.)

  IF(.NOT.LPVT)THEN
    IF(NSUPERDIA >= 2 .AND. (LHACH2.OR.LHACH3))THEN
    ELSE
      CALL PLCHHQ(0.1+(NSUPER-2)*.24,ZVT+0.01,ADJUSTL(CTIMEC(8:15))//'s',.007,0.,-1.)
    ENDIF
  ENDIF

  ENDIF    !**********************************************************
  IF(LMINMAX)THEN
    IF(LPLUS .OR. LMINUS)THEN
      CALL PCSETC('FC','/')
        CAll PLCHHQ(ZVR,ZVT+.03,YLBL(1:LEN_TRIM(YLBL)),.009,0.,+1.)
!       CAll PLCHHQ(0.68,ZVT+.03,YLBL,.009,0.,-1.)
      CALL PCSETC('FC',':')
    ELSE
      CALL PCSETC('FC','/')
!       CAll PLCHHQ(0.1+(NSUPER-1)*.24,ZVT+.01,YLBL,.007,0.,-1.)
      CAll PLCHHQ(.98,ZVT+.01+(NSUPER-1)*.02,YLBL,.007,0.,+1.)
      CALL PCSETC('FC',':')
    ENDIF
  ENDIF

END IF
!
CALL SETUSV('MI',IMI)
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL GSCLIP(1) 
CALL GSLN(1)
CALL GSPLCI(1)
CALL GSTXCI(1)
!
!*       2.14      Heading formats
!
1000 FORMAT('Vertical section IDEB=',I4,' JDEB=',I4,' ANG.=',I3,' NBPTS=',I4)
1001 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1002 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1003 FORMAT('Vertical section XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1004 FORMAT('Vertical section XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1010 FORMAT('Vertical section IDEB=',I4,' JDEB=',I4,' ANG.=',I3,' IPRO=',I4)
1011 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' IPRO=',I4)
1012 FORMAT('Vertical profile (1D)')
1013 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1014 FORMAT('Vertical section XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1015 FORMAT('Vertical section XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1016 FORMAT('Horiz. profile  NIINF=',I5,' NISUP=',I5)
1017 FORMAT('Horiz. profile  NJINF=',I5,' NJSUP=',I5)
1018 FORMAT('Vertical section IND I,J (BEGIN)-(END)=(',I4,',',I4,')-(',I4,',',I4,')')
1019 FORMAT('Vertical section LAT,LON (BEGIN)-(END)=(',F5.1,',',F5.1,')-(',F5.1,',',F5.1,')')
1020 FORMAT('Vertical section CONF. COORD.(BEGIN)-(END)=(',F8.0,',',F8.0,')-(',F8.0,',',F8.0,')')
1021 FORMAT('Vertical profile LAT,LON =(',F5.1,',',F5.1,')')
1022 FORMAT('Vertical profile CONF. COORD.=(',F8.0,',',F8.0,')')
1023 FORMAT('Vertical profile IND I,J =(',I4,',',I4,')')
1024 FORMAT('Profile =',I4)
!
!-----------------------------------------------------------------------------
!
!*       3.        EXIT
!                  ----
!
RETURN
END SUBROUTINE IMCOU_FORDIACHRO
