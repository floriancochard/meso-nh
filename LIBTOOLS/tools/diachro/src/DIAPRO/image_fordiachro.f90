!     ######spl
      SUBROUTINE IMAGE_FORDIACHRO(PTAB,KLREF,PTABINT,KNHI,KNDOT,HTEXTE)
!     #################################################################
!
!!****  *IMAGE_FORDIACHRO* - Isoncontour plots manager for horizontal 
!!                           cross-sections
!!
!!    PURPOSE
!!    -------
!       Calls the NCAR contour routines and defines the display environment
!    for the horizontal cross-section case
!
!!**  METHOD
!!    ------
!!      First, the field is checked for extrema, and the plot geometry is
!!   generated, drawing a cartographic stencil and the continental/state 
!!   outlines when required by the 'LCARTESIAN' parameters. Next, NCAR 
!!   variables are set according to the user requests, and  contours are 
!!   drawn by a call to Conpack utilities (CPRECT/CPCLDR). If a 'Z' section
!!   is requested, the topography outlines are examined to mask the contours 
!!   where map altitude intercepts terrain.
!!
!!     Notice that a TRACE-provided CPMPXY routine is used within the NCAR
!!   Conpack call to map the contoured array matrix onto the stretched model 
!!   cartographic space. The plotted data are NOT interpolated onto a regular 
!!   grid before plotting, instead a coordinate stretching technique is used.
!!   Basically, the contour calculation are made in a "grid index space"
!!   where the meshsize is uniform and equal to 1 between successive model 
!!   points (this corresponds to the x_hat_* and y_hat_* coordinates of the
!!   Meso-NH technical specification book, page 41). In this "grid index space"
!!   contourlines points are located by two floating-point index coordinates 
!!   vaying between 1 and the corresponding array dimension. This "grid index"
!!   coordinates are latter converted back to screen coordinates by CPMPXY to 
!!   obtain a correct display.  
!!
!!    EXTERNAL
!!    --------
!!      GMNMX     : computes min, max and contour increment for current field
!!      BCGRD     : when a cartographic projection applies, defines displayed
!!                  window and draws the continent/state outlines
!!      DEFENETRE : when cartesian geometry applies, defines the display window
!!      TRACEXY   : draws the model gridpont stencil as a dashline overlay
!!
!!      CPSETI !                                          INTEGER   !
!!      CPSETR !  : sets the value of a NCAR parameter,   REEL      !
!!      CPSETC !                                          CHARACTER ! NCAR
!!                                                                  !
!!      CPGETI !                                          INTEGER   !
!!      CPGETR !  : gets the value of a NCAR parameter,   REEL      !
!!      CPGETC !                                          CHARACTER !
!!                                                                  !
!!      CPRECT    : Conpack initialization                          !
!!      CPPKCL    : contour level selection                         !
!!      CPCLDR    : draws contours                                  ! Routines
!!      GSLWSC    : sets line width                                 !
!!                                                                  !
!!      ARINAM    : initialize the contour calculation as a subset  !
!!                  of areas, which may be adressed individually to !
!!                  modify their display characteristics (used for  !
!!                  topography masking here).                       ! 
!!      ARSCAM    : scans the plotting domain and defines the       !
!!                  different areas, then performs the processing   !
!!                  defined in the SFILL routine (here, hatch fill) ! 
!!      CPCLAM    : adds contour in a  previously defined area      ! NCAR
!!                                                                  ! 
!!      SET       : defines the display window in normalized and    !
!!                  user NCAR coordinates                           !
!!      GETSET    : retrieves the normalized and user NCAR          !
!!                  coordinates of a previously used window         ! Routines
!!      PLCHHQ    : prints high-quality character strings           !
!!      GSCLIP    : clips items getting out of the drawing window   !
!!
!!      CPMPXY    : TRACE provided FORTRAN-77 routine directly called
!!                  within CONPACK to map the array space onto the
!!                  cartographic space
!!      SFILL     : TRACE provided FORTAN-77 routine directly called 
!!                  CONPACK to define the hatched area used to locate
!!                  points  where the plot level intercepts topography
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_TITLE  : Declares heading variables for the plots (TRACE)
!!         NCONT  :  Current plot number
!!         CLEGEND:  Current plot heading title
!!
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!       XXX,XXY  : coordinate values for all the MESO-NH grids
!!       XXZS     : topography values for all the MESO_NH grids
!!
!!      Module MODD_CONF   : declares configuration variables of all models 
!!       LCARTESIAN: Logical for cartesian geometry :
!!                   .TRUE.  = cartesian geometry
!!                   .FALSE. = conformal projection
!!
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID    : Current MESO-NH grid indicator
!!
!!      Module MODN_PARA   : defines NAM_DOMAIN_POS namelist
!!         LHORIZ    : must be .TRUE. to perform horizontal cross esctions
!!         LVERTI    : must be .FALSE. to perform horizontal cross sections
!!         Module MODD_DIM1   : Contains dimensions
!!            NIMAX, NJMAX :  x, and y array dimensions
!!            NIINF, NISUP :  Lower and upper array bounds in x direction
!!            NJINF, NJSUP :  Lower bound and upper bound  in y direction
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist
!!                         (former NCAR common)
!!        NIOFFD     : Label normalisation (=0 none, =/=0 active)
!!        NULBLL     : Nb of contours between 2 labelled contours
!!        NIOFFM     : =0    --> message at picture bottom
!!                     =/= 0 --> no message
!!        NIOFFP     : Special point value detection
!!                    (=0 none, =/=0 active)
!!        CTYPHOR    : Horizontal cross-section type
!!                     (='K' --> model level section;
!!                      ='Z' --> constant-altitude section;
!!                      ='P' --> isobar section (planned)
!!                      ='T' --> isentrope section (planned)
!!        XSPVAL     : Special value
!!        XSIZEL     : Label size
!!        LXY        : If =.TRUE., plots  a grid-mesh stencil background
!!
!!      Module MODD_OUT       : Defines a log. unit for printing
!!        NIMAXT : x-size of the displayed section of the model array
!!        NJMAXT : y-size of the displayed section of the model array
!!
!!      Module MODD_SUPER   : defines plot overlay control variables
!!         LSUPER   : =.TRUE. --> plot overlay is active
!!                    =.FALSE. --> plot overlay is not active
!!         NSUPER   : Rank of the current plot in the overlay
!!                    sequence. The initial plot is rank 1.
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
!!      Updated   PM   06/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef NAGf95
USE F90_UNIX  ! for FLUSH and GETENV
#endif

USE MODD_TITLE
USE MODD_MASK3D
USE MODD_COORD
USE MODD_NMGRID
USE MODD_CONF
USE MODN_PARA
USE MODN_NCAR
USE MODD_TIME
USE MODD_TIME1
USE MODD_OUT
USE MODD_SUPER
USE MODD_LUNIT1
USE MODD_RESOLVCAR
USE MODD_HACH
USE MODD_TIT
USE MODD_ALLOC_FORDIACHRO
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODI_READMNMXINT_ISO
USE MODI_READREFINT_ISO
USE MODI_READXISOLEVP
USE MODD_CTL_AXES_AND_STYL
USE MODD_RSISOCOL
!
USE MODI_CREATLINK
USE MODI_WRITEDIR
!
IMPLICIT NONE
!
!        0.0   TRACE interface with the "CPMPXY" routine of the NCAR package
!
! NOTICE:  The CPMPXY and the NCAR graphical utilities are NOT written
! ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE 
!          does not follow the Meso-NH usual rules: it has to be made using 
!          a COMMON stack with  static memory allocation of XZZXX and
!          XZZXY arrays.
!
COMMON/TEMH/XZZXX,XZZXY,NIIMAX,NIJMAX
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
COMMON/COLAREA/ICOL(300)
COMMON/HACHAREA/IHACH(300)
#include "big.h"
!
REAL,DIMENSION(N2DVERTX) :: XZZXX
REAL,DIMENSION(N2DVERTX) ::  XZZXY
INTEGER             :: NIIMAX, NIJMAX
LOGICAL             :: LVERT, LHOR, LPT, LXABS
INTEGER             :: ICOL
INTEGER             :: IHACH
!
!*       0.1   Declarations of dummy arguments and results
!
INTEGER             :: KNHI       ! Extrema processing option
INTEGER             :: KNDOT      ! Line style option 
INTEGER             :: KLREF      ! Cross-section altitude (or Model Level
                                  ! or Pressure depending on user's vertical
                                  ! coordinate choice)

CHARACTER(LEN=*)    :: HTEXTE     ! Plot heading with variable name

REAL                :: PTABINT    ! Contour increments for current plot

REAL,DIMENSION(:,:) :: PTAB       ! Variable array to be plotted

!
!*       0.2   Local variables
!
INTEGER             :: IM, IL, ILE
INTEGER             :: J, JJ, JI, JU, JK
INTEGER             :: JLBL, JL
INTEGER             :: I, ICLD, INCL
INTEGER             :: INBC
INTEGER             :: INBX,INBY
INTEGER,SAVE        :: IDX
INTEGER,SAVE        :: INBCT
INTEGER,SAVE        :: ILUCOL, IRESP
INTEGER,DIMENSION(:),ALLOCATABLE :: ICOL2
INTEGER             :: ILENT
INTEGER             :: ISTA, IER, IWK, INB, INBB
INTEGER             :: INUM, ILOOP, JLOOPI, JLOOPJ, IDEB, IFIN, II
INTEGER,SAVE        :: IH, IHT
INTEGER,DIMENSION(32):: INDHACHREF=(/0,54,52,60,14,59,58,1,57,56,55,54,53,52,51,50, &
                        1,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35/)
#ifdef RHODES
INTEGER          :: ISTAF
#endif
CHARACTER(LEN=80)               :: YCAR80 
CHARACTER(LEN=320)               :: YCAR320
CHARACTER(LEN=70)               :: YPLANH 
CHARACTER(LEN=100)               :: YTEM
CHARACTER(LEN=40)                :: YTEM40
CHARACTER(LEN=8),DIMENSION(300) :: YLLBS  
CHARACTER(LEN=32),SAVE          :: YNAMTABCOL
CHARACTER(LEN=32)               :: YLBL
CHARACTER(LEN=32)               :: YTEXT
CHARACTER(LEN=20)               :: YCAR20
CHARACTER(LEN=4)                :: YC4, YC42
CHARACTER(LEN=1)                :: YREP

LOGICAL :: GISO

REAL,DIMENSION(300) :: ZLEV, ZISOLEVP
REAL                :: ZTABMIN, ZTABMAX, ZTABINT
REAL                :: ZTABMN, ZTABMX
REAL                :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
REAL                :: ZZSPVAL, ZISO
REAL                :: ZLREF,ZWIDTH
REAL                :: ZCLV
REAL                :: RED,GREEN,BLUE
REAL                :: ZINTERV
REAL                :: ZMIN, ZMAX
REAL,SAVE           :: ZSC
REAL,SAVE           :: ZVLDEF, ZVRDEF, ZVBDEF, ZVTDEF
REAL,SAVE           :: ZD, ZF, ZVERA, ZINTE
REAL                :: ZX, ZY
REAL                :: ZXPOSTITT1, ZXYPOSTITT1
REAL                :: ZXPOSTITT2, ZXYPOSTITT2
REAL                :: ZXPOSTITT3, ZXYPOSTITT3
REAL                :: ZXPOSTITB1, ZXYPOSTITB1
REAL                :: ZXPOSTITB2, ZXYPOSTITB2
REAL,SAVE           :: ZXPOSTITB3, ZXYPOSTITB3
REAL                :: ZSZTITVAR1, ZSZTITVAR
REAL,DIMENSION(5)   :: ZX5, ZY5
REAL                :: ZEPX, ZEPY
!
!       0.3    Work arrays for NCAR 
!
! See aforementioned notice. The dimensions of these arrays are 
! subject to possible tuning, but have to be prescribed. Add
! extra size if necessary.
!

INTEGER                 :: ID, ICL, III
INTEGER,PARAMETER       :: JPLRWK=50000, JPLIWK=50000
INTEGER,PARAMETER       :: JPMAP=8000000, JPAREAGRP=300, JPWRK=50000

INTEGER,DIMENSION(JPLIWK)   :: IWRK
INTEGER,DIMENSION(JPMAP)    :: IIMAP
INTEGER,DIMENSION(JPAREAGRP):: IAREA, IGRP

REAL,DIMENSION(JPLRWK)      :: ZRWRK
REAL,DIMENSION(JPWRK)       :: ZXWRK, ZYWRK
!
! SFILL subroutine declared as external provides area control
! in some parts of the contour plot.
!
EXTERNAL SFILL
EXTERNAL SFILLH
EXTERNAL CCOLR
!
!---------------------------------------------------------------------------
!
!*      1.    DISPLAY ENVIRONMENT SETUP
!             -------------------------
!
! Recuperation du nom du processus dans YTEXT
!
NLUOUT=6
YTEXT(1:LEN(YTEXT))=' '
HTEXTE=ADJUSTL(HTEXTE)
DO JJ=1,LEN_TRIM(HTEXTE)
  IF(HTEXTE(JJ:JJ) == ' ')THEN
    YTEXT(1:JJ-1)=HTEXTE(1:JJ-1)
    EXIT
  ENDIF
  IF(JJ == LEN_TRIM(HTEXTE))THEN
    YTEXT=HTEXTE
  ENDIF
ENDDO
YTEXT=ADJUSTL(YTEXT)
!
!*      1.1   Size computations and gridpoint location loading for NCAR  
!
IM=SIZE(PTAB,1)
IL=SIZE(PTAB,2)
ZTABINT=PTABINT
LHORIZ=.TRUE.; LVERTI=.FALSE.
LVERT=LVERTI
LHOR=LHORIZ
LPT=LPXT
! Min and max
ZMIN=PTAB(IM/2,IL/2); ZMAX=PTAB(IM/2,IL/2)
IF(ZMIN == XSPVAL)ZMIN=1.E16
IF(ZMAX == XSPVAL)ZMAX=-1.E16
!ZMIN=999999.; ZMAX=-999999.
if(nverbia > 0)then
  print *,' ** image AV DO JJ=1,IL'
endif
DO JJ=1,IL
  DO JI=1,IM
    IF(PTAB(JI,JJ) /= 888. .AND. PTAB(JI,JJ) /= XSPVAL)THEN
      IF(PTAB(JI,JJ) < ZMIN)ZMIN=PTAB(JI,JJ)
      IF(PTAB(JI,JJ) > ZMAX)ZMAX=PTAB(JI,JJ)
    ENDIF
  ENDDO
ENDDO
YLBL(1:5)='(Min:'
WRITE(YLBL(6:15),'(E10.3)')ZMIN
YLBL(16:21)=', Max:'
WRITE(YLBL(22:31),'(E10.3)')ZMAX
YLBL(32:32)=')'

!

!NIIMAX=NIMAXT
!NIJMAX=NJMAXT
NIIMAX=SIZE(PTAB,1)
NIJMAX=SIZE(PTAB,2)
XZZXX(1:NIIMAX)=XXX(NIINF:NISUP,NMGRID)
XZZXY(1:NIJMAX)=XXY(NJINF:NJSUP,NMGRID)
!

IF(LPRINT)THEN
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=SIZE(PTAB,1)/5
 
!!Octobre 2001 Cas des trajectoires ??
  print *,' ** image, ILOOP,NLOOPT, SIZE(PTAB,1) ',ILOOP,NLOOPT, SIZE(PTAB,1)
!!Octobre 2001
  IF(ILOOP * 5 < SIZE(PTAB,1))ILOOP=ILOOP+1
  WRITE(INUM,'(''CH  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
& CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
  IF(LMINUS .OR. LPLUS)THEN
    WRITE(INUM,'(A55,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITB3(1:55)
  ELSE
    WRITE(INUM,'(A40,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITGAL
  ENDIF
  WRITE(INUM,'(''niinf'',i4,'' njinf'',i4,'' nisup'',i4,'' njsup'',i4,&
&''   '',A1,'' '',i6)')&
  &NIINF,NJINF,NISUP,NJSUP,CTYPHOR,KLREF
  WRITE(INUM,'(''NBVAL en I '',i4,''  NBVAL en J '',i4,''   iter'',i3)') &
  &NISUP-NIINF+1,NJSUP-NJINF+1,ILOOP
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**IMAGE XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
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
! JUin 2001 Ecriture des dates 
  DO JLOOPI=1,ILOOP
    IF(JLOOPI == 1)THEN
      IDEB=1; IFIN=5
      IDEB=IDEB+NIINF-1; IFIN=IFIN+NIINF-1
    ELSE
      IDEB=IFIN+1; IFIN=IFIN+5
    ENDIF
    IF(JLOOPI == ILOOP)THEN
      IFIN=SIZE(PTAB,1)+NIINF-1
    ENDIF
    
    WRITE(INUM,'(1X,78(1H*))')
    WRITE(INUM,'('' J   I-> '',3X,I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
    WRITE(INUM,'(''.'',78(1H*))')
    DO JLOOPJ=SIZE(PTAB,2),1,-1
      WRITE(INUM,'(I4,2X,5(1X,E14.7))')JLOOPJ+NJINF-1,(PTAB(II,JLOOPJ),II=IDEB-NIINF+1,IFIN-NIINF+1)
!     WRITE(INUM,'(I3,2X,5E15.8)')JLOOPJ+NJINF-1,(PTAB(II,JLOOPJ),II=IDEB-NIINF+1,IFIN-NIINF+1)
    ENDDO
    WRITE(INUM,'(1X,78(1H*))')
  ENDDO
ENDIF
IF(LPRINTXY)THEN
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  WRITE(INUM,'(''CH XY  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
& CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
  IF(LMINUS .OR. LPLUS)THEN
    WRITE(INUM,'(A55,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITB3(1:55)
  ELSE
    WRITE(INUM,'(A40,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITGAL
  ENDIF
  WRITE(INUM,'(''niinf'',i4,'' njinf'',i4,'' nisup'',i4,'' njsup'',i4,&
&'' '',A1,'' '',i6)')&
  &NIINF,NJINF,NISUP,NJSUP,CTYPHOR,KLREF
  WRITE(INUM,'(''NBVAL en I '',i4,''  NBVAL en J '',i4)') &
  &NISUP-NIINF+1,NJSUP-NJINF+1

  II=MAX(SIZE(PTAB,1),SIZE(PTAB,2))
  WRITE(INUM,'(1X,73(1H*))')
  WRITE(INUM,'(26X,''X'',38X,''Y'')')
  WRITE(INUM,'(1X,73(1H*))')
  DO JLOOPJ=1,II
    IF(JLOOPJ ==1)THEN
	YC4='    '
	YC42='    '
	WRITE(YC4,'(I4,'')'')')NIINF
	WRITE(YC42,'(I4,'')'')')NJINF
	WRITE(INUM,'(''NIINF('',A4,I4,5X,E15.8,5X,''NJINF('',A4,I4,5X,E15.8)') &
	YC4,JLOOPJ,XZZXX(JLOOPJ),YC42,JLOOPJ,XZZXY(JLOOPJ)
	YC4='    '
	YC42='    '
	WRITE(YC4,'(I4,'')'')')NISUP
	WRITE(YC42,'(I4,'')'')')NJSUP
    ELSE
	IF(SIZE(PTAB,1) > SIZE(PTAB,2))THEN
	  IF(JLOOPJ < SIZE(PTAB,2))THEN
	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZXX(JLOOPJ), &
	    JLOOPJ,XZZXY(JLOOPJ)
	  ELSE IF(JLOOPJ == SIZE(PTAB,1))THEN
	    WRITE(INUM,'(''NISUP('',A4,I4,5X,E15.8)')YC4,JLOOPJ,XZZXX(JLOOPJ)
            WRITE(INUM,'(1X,73(1H*))')
	  ELSE IF(JLOOPJ == SIZE(PTAB,2))THEN
	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,''NJSUP('',A4,I4,5X,E15.8)')&
	    JLOOPJ,XZZXX(JLOOPJ), &
	    YC42,JLOOPJ,XZZXY(JLOOPJ)
	  ELSE IF(JLOOPJ > SIZE(PTAB,2))THEN
	    WRITE(INUM,'(5X,I9,5X,E15.8)')JLOOPJ,XZZXX(JLOOPJ)
	  ENDIF
	ELSE IF(SIZE(PTAB,2) > SIZE(PTAB,1))THEN
	  IF(JLOOPJ < SIZE(PTAB,1))THEN
	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZXX(JLOOPJ), &
	    JLOOPJ,XZZXY(JLOOPJ)
	  ELSE IF(JLOOPJ == SIZE(PTAB,2))THEN
	    WRITE(INUM,'(29X,5X,5X,''NJSUP('',A4,I4,5X,E15.8)') &
	    YC42,JLOOPJ,XZZXY(JLOOPJ)
            WRITE(INUM,'(1X,73(1H*))')
	  ELSE IF(JLOOPJ > SIZE(PTAB,1))THEN
	    WRITE(INUM,'(29X,5X,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZXY(JLOOPJ)
	  ELSE
	    WRITE(INUM,'(''NISUP('',A4,I4,5X,E15.8,5X,5X,I9,5X,E15.8)') &
	    YC4,JLOOPJ,XZZXX(JLOOPJ), &
	    JLOOPJ,XZZXY(JLOOPJ)
	  ENDIF
	ELSE
	  IF(JLOOPJ == SIZE(PTAB,2))THEN
	    WRITE(INUM,'(''NISUP('',A4,I4,5X,E15.8,5X,''NJSUP('',A4,I4,5X,E15.8)') &
	    YC4,JLOOPJ,XZZXX(JLOOPJ), &
	    YC42,JLOOPJ,XZZXY(JLOOPJ)
            WRITE(INUM,'(1X,73(1H*))')
	  ELSE
	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZXX(JLOOPJ), &
	    JLOOPJ,XZZXY(JLOOPJ)
	  ENDIF
	ENDIF
    ENDIF
  ENDDO
ENDIF
!
!*      1.2  Scans for data extrema. Selects display window.
!            If required by LCARTESIAN: selects cartographic projection 
!            and draws coastlines. 
!            If required by LXY: draws a gripoint stencil over the contours.
!

! Modifs for diachro
!
!CALL GMNMX(ZTABMIN,ZTABMAX,ZTABINT)

if(nverbia > 0)then
  print *,' ** image IF(NIMNMX == '
endif
IF(NIMNMX == 0 .OR. NIMNMX == 1)THEN
  LISOK=.FALSE.
  ZTABMIN=0.; ZTABMAX=0.
  CALL READMNMXINT_ISO(NIMNMX,YTEXT(1:LEN_TRIM(YTEXT)),ZTABMIN,ZTABMAX,ZTABINT)

ELSE IF(NIMNMX == 2)THEN
  ZISOLEVP(:)=9999.
  CALL READXISOLEVP(YTEXT(1:LEN_TRIM(YTEXT)),ILE,ZISOLEVP)
  IF(NVERBIA > 5)THEN
    print *,' IMAGE YTEXT,ILE,ZISOLEVP ',YTEXT(1:LEN_TRIM(YTEXT)),ILE,ZISOLEVP(1:ILE)
  ENDIF

ELSE IF (NIMNMX==3) THEN  ! compute contour values from XISOREF and XDIAINT
  ZISOLEVP(:)=9999.
  ZTABMN=MINVAL(PTAB,MASK=PTAB/=XSPVAL) 
  ZTABMX=MAXVAL(PTAB,MASK=PTAB/=XSPVAL)
  CALL READREFINT_ISO(YTEXT(1:LEN_TRIM(YTEXT)),ZTABMN,ZTABMX,ZTABINT,ZISOLEVP)
ENDIF

IF(LCARTESIAN)THEN
  ZVLDEF=.1; ZVRDEF=.9; ZVBDEF=.1; ZVTDEF=.9
ELSE
  ZVLDEF=.05; ZVRDEF=.95; ZVBDEF=.05; ZVTDEF=.95
ENDIF
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
  IF(NSUPER == 1)THEN
    IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(1)   
    IF(LCARTESIAN)CALL DEFENETRE
  END IF
ELSE
  IF(XLW >= 0)THEN
    XLWIDTH=XLW
  ENDIF
  IF(XLW1 >= 0)THEN
    XLWIDTH=XLW1
  ENDIF
  IH=0; IHT=0
  IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(1)   
  IF(LCARTESIAN)CALL DEFENETRE
END IF
!
!IF(LXY)THEN
! CALL GSCLIP(0)
! CALL TRACEXY
!END IF
!
if(nverbia > 0)then
  print *,' ** image AV CALL GSLWSC(1.)'
endif
CALL GSLWSC(1.)
!CALL CPSETI('CFC',1)
!
!
!*      1.3  Selects contour range and increment according to NIMNMX
!
SELECT CASE(NIMNMX)
    
  CASE(-1)                           ! Fully automatic scanning
    CALL CPSETI('CLS',+16)
    IF((LHACH1 .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))) .OR. &
      (LHACH2 .AND. NSUPER == 2)                                   .OR. &
      (LHACH3 .AND. NSUPER == 3)                                   .OR. &
      (LHACH4 .AND. NSUPER == 4))CALL CPSETI('CLS',+7)
    CALL CPSETR('CIS',-ZTABINT)

  CASE(0)                            ! Automatic range and given increment
    CALL CPSETI('CLS',16)
    CALL CPSETR('CIS',ZTABINT)
    CALL CPSETI('LIS',NULBLL+1)
    CALL CPSETR('CMN',10000000000.)
!   CALL CPSETR('CMN',MAXVAL(PTAB))
    CALL CPSETR('CMX',1000000000.)
!   CALL CPSETR('CMX',MINVAL(PTAB))

  CASE(1)                            ! Given range and increment
    IF(ZTABMAX == ZTABMIN)THEN
      ICL=1
      CALL CPSETI('NCL',ICL)
    ELSE
      ICL=NINT((ZTABMAX-ZTABMIN)/ZTABINT)
      IF(NVERBIA >= 5)THEN
      print *,' ztabmin  max, int,ICL ',ZTABMIN,ZTABMAX,ZTABINT,ICL
      ENDIF
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      IF(ZTABMIN + ICL*ZTABINT <= ZTABMAX)ICL=ICL+1
      IF(NVERBIA >= 5)THEN
      print *,' ztabmin  max, int,ICL ',ZTABMIN,ZTABMAX,ZTABINT,ICL
      ENDIF
!     IF(ZTABMIN + ICL*ZTABINT < ZTABMAX)ICL=ICL+1
      CALL CPSETI('NCL',ICL)
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
!     IF((LCOLAREA .OR. LHACH1) .AND. (.NOT.LSUPER.OR. (LSUPER .AND. NSUPER == 1))) CALL CPSETI('NCL',ICL+1)
    ENDIF
    CALL CPSETI('CLS',0)
    ZISO=ZTABMIN-ZTABINT
    DO I=1,ICL
    CALL CPSETI('PAI',I)
    CALL CPSETI('AIA',I+1)
    CALL CPSETI('AIB',I)
    ZISO=ZISO+ZTABINT
    IF(ABS(ZISO)<1.E-20)ZISO=0.
    CALL CPSETR('CLV',ZISO)
    CALL CPSETR('CLU',1.)
    IF(.NOT.LSUPER.OR. (LSUPER .AND. NSUPER == 1))THEN
      IF(LBLUSER1)THEN
        DO JLBL=1,SIZE(XLBLUSER1)
         DO JL=-20,20,1
           IF(ZISO == XLBLUSER1(JLBL)*10.**FLOAT(JL))THEN
             CALL CPSETR('CLU',3.)
!            print *,' ISO LABELLE ',ZISO
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
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
!   IF(ICL /= 1)THEN
!     IF((LCOLAREA .OR. LHACH1) .AND. (.NOT.LSUPER.OR. (LSUPER .AND. NSUPER == 1)))THEN
!       ICL=ICL+1
!       CALL CPSETI('PAI',ICL)
!       CALL CPSETI('AIB',ICL)
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
!       CALL CPSETI('AIA',ICL+1)
!       ZISO=ZISO+ZTABINT
!       IF(ABS(ZISO)<1.E-20)ZISO=0.
!       CALL CPSETR('CLV',ZISO)
!     END IF
!   END IF

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
  CASE(2,3)                   
    ICL=0
    DO I=1,10000
      ICL=ICL+1
      IF(NIMNMX==3 .OR. (NIMNMX==2 .AND.LISOLEVP))THEN
        ZLEV(ICL)=ZISOLEVP(ICL)
        IF(NVERBIA > 5)then
          print *,' ICL ZLEV ',ICL,ZLEV(ICL)
        ENDIF
      ELSE IF (NIMNMX==2 .AND. .NOT.LISOLEVP) THEN ! Given contour values     
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
    !
    ICL=ICL-1
    CALL CPSETI('NCL',ICL)
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
!            print *,' ISO LABELLE ',ZLEV(I)
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
END SELECT
!
!*      1.4  A few cosmetic parameter settings
!
! Label format and normalization
!
if(nverbia > 0)then
  print *,' ** image AV CASE(NIOFFD)',NIOFFD
endif
SELECT CASE(NIOFFD)        
    
CASE(0)                     !! No label normalisation, decimal point kept
    III=9                   ! 'Numeric exponent use flag'
    CALL CPSETI('NEU',III)  ! III > 0 --> decimal point kept if the number of
                            ! significant digits < III; else form requiring the 
                            ! fewest character is used
    CALL CPSETI('NOF',7)
    IF(NSD /= 0)THEN
      CALL CPSETI('NSD',-NSD)
    ELSE
      CALL CPSETI('NSD',-6)
    ENDIF
CASE DEFAULT                !! Label normalisation, scale factor right of the plot
    CALL CPSETI('NEU',-2)   ! Exponential notation forced, in any case 
    CALL CPSETI('NOF',7)
    CALL CPSETI('NET',0)    ! Exponent shown as a "E"

END SELECT
!
! Special value handling
!
SELECT CASE(NIOFFP)
    
CASE(0)                          ! No special value used
    CALL CPSETR('SPV',0.)
CASE DEFAULT                     ! XSPVAL used as a special value
    CALL CPSETR('SPV',XSPVAL)

END SELECT
!
! Information label under the plot
!
SELECT CASE(NIOFFM)
    
CASE(0)               ! A label is printed to the plot bottom
CASE DEFAULT          ! No label
    CALL CPSETC('ILT',' ')

END SELECT
!
!!!!!!!! PROVI
CALL GSCLIP(1)              ! Display clipping activated
!CALL GSCLIP(0)              ! Display clipping activated
!!!!!!!! PROVI
CALL CPSETI('MAP',4)        ! A specific map projection is used, as provided in
                            ! the user-provided "CPMPXY" routine. This important
                            ! parameter informs Conpack of the kind of geographic
                            ! transformation actually made.
CALL CPSETI('SET',0)        ! No "SET" issued by conpack
CALL CPSETR('SPV',XSPVAL)
!
!-------------------------------------------------------------------------------
!
!*      3.   FIELD CONTOURS DRAWING
!            ----------------------       
!
!*      3.1  Conpack initialization
!
if(nverbia > 0)then
  print *,' ** image AV CPRECT(PTAB,IM,IM',IM,IL
endif
CALL CPRECT(PTAB,IM,IM,IL,ZRWRK,JPLRWK,IWRK,JPLIWK)
CALL CPSETR('CWM',XSIZEL/.01)

INCL=0
CALL CPPKCL(PTAB,ZRWRK,IWRK)
CALL CPGETI('NCL',INCL)

!
!*    3.1a     High and low handling
!
SELECT CASE(KNHI)
    
  CASE(0)                           ! H + L   are displayed
! Test rajoute pour eviter la superposition de CONSTANT FIELD ici et ensuite
! avec le 2eme CPLBDR utile en cas de surfaces colorees
    IF(INCL /= 0)THEN
      CALL CPLBDR(PTAB,ZRWRK,IWRK)
    ENDIF
  CASE DEFAULT                      ! TO BE REVISED*********************
			            ! <0  --> no action (:-1 to be set)
			            ! >0  --> gridpoint value displayed
                                    !         (1: to be set)
END SELECT
!
!*     3.2   Line style and color handling
!
!INCL=0
!CALL CPPKCL(PTAB,ZRWRK,IWRK)
!CALL CPGETI('NCL',INCL)
IF(NIMNMX < 0)THEN
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
             print *,' ISO LABELLE ',ZISO
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
SELECT CASE(KNDOT)
  CASE(0,1,1023,65535)           ! Solid lines
    DO J=1,INCL
      CALL CPSETI('PAI',J)
      CALL CPSETI('CLD',65535)
    ENDDO
  CASE (:-1)                     ! <0 Negative value dashed, positive value solid
      ICLD=ABS(KNDOT)
!     write(0,*)' KNDOT',KNDOT,' INCL ',INCL
        DO J=1,INCL
          CALL CPSETI('PAI',J)
          CALL CPGETR('CLV',ZCLV)
          IF(ZCLV.GE.0.)CALL CPSETI('CLD',65535)
          IF(ZCLV.LT.0.)CALL CPSETI('CLD',ICLD)
!         write(0,*)' J ZCLV',I,ZCLV
        ENDDO
  CASE DEFAULT                   ! KNDOT used as a dash pattern
      ICLD=ABS(KNDOT)
        DO J=1,INCL
          CALL CPSETI('PAI',J)
          CALL CPSETI('CLD',ICLD)
        ENDDO
END SELECT

!
! **************************************************************
! Surfaces en hachures ; LHACHx=.TRUE. (avec x=1 ou 2 ou 3 ou 4)
! **************************************************************

IF((LHACH1 .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))) .OR. &
   (LHACH2 .AND. NSUPER == 2)                                   .OR. &
   (LHACH3 .AND. NSUPER == 3)                                   .OR. &
   (LHACH4 .AND. NSUPER == 4))THEN !++++++++++++++++++++++++++++++++++++++++++

  IF(NSUPER > 1)THEN
    IH=IH+1
    if(nverbia >0)then
      print *,' image: HACHures IHT IH ',IHT,IH
    endif
  ENDIF

  WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' VALEURS:'
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
    WRITE(NLUOUT,*)' >>>>>>>SELECTION DES HACHURES PAR L''UTILISATEUR'
    WRITE(NLUOUT,*)' >>>>>>>VOUS DEVEZ FOURNIR ',INCL+1,' INDICES'
    WRITE(NLUOUT,*)' Rentrez sur 1 premiere ligne le nombre d''indices fournis dans la ligne suivante'
    WRITE(NLUOUT,*)' Puis sur la(es) ligne(s) suivante(s) les indices des figures pris dans la table' 
    WRITE(NLUOUT,*)' de reference correspondant aux isocontours ranges par ordre croissant'
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
    !WRITE(YCAR80,*)INBC
    !WRITE(NDIR,'(A80)')YCAR80
    CALL WRITEDIR(NDIR,INBC)
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
   ! WRITE(YCAR320,*)IHACH(1:INBC)
   ! YCAR320=ADJUSTL(YCAR320)
   ! ILENT=LEN_TRIM(YCAR320)
   ILENT=INBC*4
   !! car plantage dans le cas ELSE si ILENT=80 !!
    IF(ILENT == 80 ) THEN
     ! YCAR320=TRIM(YCAR320)//' '
      ILENT=ILENT+1
    END IF
    IF(ILENT > 240 )THEN
      !WRITE(YCAR80,*)IHACH(1:INBC/4)
      CALL WRITEDIR(NDIR,IHACH(1:INBC/4))
      !WRITE(YCAR80,*)IHACH(INBC/4+1:INBC/2)
      CALL WRITEDIR(NDIR,IHACH(INBC/4+1:INBC/2))
      !WRITE(YCAR80,*)IHACH(INBC/2+1:3*INBC/4)
      CALL WRITEDIR(NDIR,IHACH(INBC/2+1:3*INBC/4))
      !WRITE(YCAR80,*)IHACH(3*INBC/4+1:INBC)
      CALL WRITEDIR(NDIR,IHACH(3*INBC/4+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE IF(ILENT > 160 )THEN
     ! WRITE(YCAR80,*)IHACH(1:INBC/3)
      CALL WRITEDIR(NDIR,IHACH(1:INBC/3))
     ! WRITE(YCAR80,*)IHACH(INBC/3+1:2*INBC/3)
      CALL WRITEDIR(NDIR,IHACH(INBC/3+1:2*INBC/3))
     ! WRITE(YCAR80,*)IHACH(2*INBC/3+1:INBC)
      CALL WRITEDIR(NDIR,IHACH(2*INBC/3+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE IF(ILENT > 80 )THEN
    !  WRITE(YCAR80,*)IHACH(1:INBC/2)
      CALL WRITEDIR(NDIR,IHACH(1:INBC/2))
   !   WRITE(YCAR80,*)IHACH(INBC/2+1:INBC)
      CALL WRITEDIR(NDIR,IHACH(INBC/2+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE
      !WRITE(YCAR80,*)IHACH(1:INBC)
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
!   call mapbla(iimap)
    CALL CPCLAM(PTAB,ZRWRK,IWRK,IIMAP)
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
    IF(.NOT.LSUPER .OR. NSUPER == 1 .OR. (NSUPER == 2 .AND. LARROVL .AND. NSUPERDIA == 2))THEN
    IF(ZVR < .8999999)THEN
      print *,' ZVR < .9 ',ZVR
      CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(ZVR+.2,1.)-ZVR)/10.,MIN(ZVR+.2,1.),ZVB,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,1)
    ELSE
        IF(INCL <= 8)THEN
	  if(nverbia >0)then
          print *,' INCL <= 8 ',INCL
	  endif
          CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB+(ZVT-ZVB)/4.,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,1)
        ELSE
	  if(nverbia >0)then
          print *,' INCL > 8 ',INCL
	  endif
          CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,1)
        ENDIF
!       CALL LBLBAR_FORDIACHRO(1,ZVR,1.,ZVB,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,2)
    ENDIF

    ELSE

      ZVERA=ZVR-(ZVR-ZVL)/3.

      IF(IHT == 0)THEN
        IF(NSUPER == 2 .AND. LARROVL .AND. NSUPERDIA > 2)THEN
          ZD=ZVL; ZF=ZVERA
          IF(INCL == 1)THEN
            ZF=ZF-(ZF-ZD)/2.
          ELSE IF(INCL <= 4)THEN
            ZF=ZF-(ZF-ZD)/4.
          ENDIF
          CALL LBLBAR_FORDIACHRO(0,ZD,ZF,ZVT+.01,ZVT+.04,INCL+1,1.,.33,IHACH,2,YLLBS,INCL,2)
        ELSE
        print *,' ** Image IHT=0 -> pas de trace de la table de hachures. Cas imprevu .. A voir.. '
        ENDIF
      ELSE

      ZINTE=(ZVERA-ZVLDEF)/FLOAT(IHT)
      IF(IHT == 1)THEN
	ZD=ZVL; ZF=ZVERA
      ELSE IF(IHT == 2 .OR. IHT == 3)THEN
	ZD=ZVLDEF+ZINTE*(IH-1)
	ZF=ZVLDEF+ZINTE*(IH)-.01
      ENDIF
      IF(INCL == 1)THEN
        ZF=ZF-(ZF-ZD)/2.
      ELSE IF(INCL <= 4)THEN
        ZF=ZF-(ZF-ZD)/4.
      ENDIF
      CALL LBLBAR_FORDIACHRO(0,ZD,ZF,ZVT+.01,ZVT+.04,INCL+1,1.,.33,IHACH,2,YLLBS,INCL,2)

      ENDIF
    ENDIF
    CALL GSFAIS(0)
!
! Definition de la couleur des isos (0 -> blanc sur papier; 1 -> noir sur papier)
    IF(LISOWHI)CALL GSPLCI(0)
    IF(LISOWHI)CALL GSTXCI(0)

!
ELSE IF(LCOLAREA)THEN   !+++++++++++++++++++++++++++++++++++++++++++++++++++++

! **************************************************************************
! Surfaces couleur (reservees aux dessins avec ou sans superpositions;
! LCOLAREA=.TRUE.) . En cas de superpositions, obligatoirement le 1er dessin
! **************************************************************************

  IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN        !00000000000000000000000000000000000000000000

    IF(.NOT.LCOLAREASEL)THEN     !====================================
!
! Selection automatique des couleurs par le programme
! ***************************************************
!
if(nverbia > 0)then
  print *,' ** image AV COLOR_FORDIACHRO(INCL+1) ,INCL',INCL
endif
       CALL COLOR_FORDIACHRO(INCL+1,1)
       WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' VALEURS:'
       IF(INCL /= 0)then
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('AIB',J)
         CALL CPSETI('AIA',J+1)
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
         ICOL(J)=J+2
if(nverbia > 2)then
  print *,' ** image AV GENFORMAT ZCLV ',ZCLV
endif
         CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
       ENDDO
       ENDIF
       ICOL(INCL+1)=INCL+3
if(nverbia > 0)then
  print *,' ** image ICOL(INCL+1) ',ICOL(INCL+1)
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
       IF(LCOLZERO)THEN
	 ICOL(NCOLZERO)=0
       ENDIF
       WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
       WRITE(NLUOUT,*)ICOL(INCL+1)
    ELSE                         !====================================
!
! Selection des couleurs par l'utilisateur
! ****************************************
!
      IF(LTABCOLDEF)THEN
! Choix de la table de couleurs par defaut
        WRITE(NLUOUT,*)' <<< TABCOLDEF >>>'
        CALL TABCOL_FORDIACHRO
      ELSE
! Choix d'une table creee par l'utilisateur
        CALL FMLOOK(YNAMTABCOL,CLUOUT,ILUCOL,IRESP)
        IF(IRESP == -54)THEN
          YNAMTABCOL(1:32)=' '
! Lecture du nom de la table de couleurs (1 seule fois)
          print *,' Entrez le nom de VOTRE TABLE de COULEURS '
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
! Lecture du nb de couleurs de la table, des index de couleur et des
! proportions relatives de rouge, vert, bleu
        CALL GQOPS(ISTA)
        CALL GQACWK(1,IER,INB,IWK)
!print *,' COLOR_FORDIACHRO AP GQACWK INB IWK ',INB,IWK
        CALL GQOPWK(1,IER,INB,IWK)
        READ(ILUCOL,*)INBCT
        DO J=1,INBCT
          READ(ILUCOL,*)IDX,RED,GREEN,BLUE
          DO JU=1,INB
	    CALL GQOPWK(JU,IER,INBB,IWK)
	    IF(IWK == 9)THEN
	      CYCLE
	    ELSE
              CALL GSCR(IWK,IDX,RED,GREEN,BLUE)
!             CALL GSCR(1,IDX,RED,GREEN,BLUE)
            ENDIF
	  ENDDO
        ENDDO
      ENDIF ! fin d'une table creee par l'utilisateur
      WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' VALEURS:'
      DO J=1,INCL
        CALL CPSETI('PAI',J)
        CALL CPSETI('AIB',J)
        CALL CPSETI('AIA',J+1)
        CALL CPGETR('CLV',ZCLV)
        ZLEV(J)=ZCLV
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
      !WRITE(YCAR80,*)INBC
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
   ! WRITE(YCAR320,*)ICOL(1:INBC)
   ! YCAR320=ADJUSTL(YCAR320)
   ! ILENT=LEN_TRIM(YCAR320)
   ! print*,"YCAR320=",YCAR320
   ! print*,"ILENT=",ILENT
   ILENT=INBC*4
    IF(ILENT == 80 ) THEN
     ! YCAR320=TRIM(YCAR320)//' '
      ILENT=ILENT+1
    END IF
    IF(ILENT > 240 )THEN
     ! WRITE(YCAR80,*)ICOL(1:INBC/4)
      CALL WRITEDIR(NDIR,ICOL(1:INBC/4))
     ! WRITE(YCAR80,*)ICOL(INBC/4+1:INBC/2)
      CALL WRITEDIR(NDIR,ICOL(INBC/4+1:INBC/2))
     ! WRITE(YCAR80,*)ICOL(INBC/2+1:3*INBC/4)
      CALL WRITEDIR(NDIR,ICOL(INBC/2+1:3*INBC/4))
     ! WRITE(YCAR80,*)ICOL(3*INBC/4+1:INBC)
      CALL WRITEDIR(NDIR,ICOL(3*INBC/4+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE IF(ILENT > 160 )THEN
     ! WRITE(YCAR80,*)ICOL(1:INBC/3)
      CALL WRITEDIR(NDIR,ICOL(1:INBC/3))
     ! WRITE(YCAR80,*)ICOL(INBC/3+1:2*INBC/3)
      CALL WRITEDIR(NDIR,ICOL(INBC/3+1:2*INBC/3))
     ! WRITE(YCAR80,*)ICOL(2*INBC/3+1:INBC)
      CALL WRITEDIR(NDIR,ICOL(2*INBC/3+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE IF(ILENT > 80 )THEN
      !WRITE(YCAR80,*)ICOL(1:INBC/2)
      CALL WRITEDIR(NDIR,ICOL(1:INBC/2))
     ! WRITE(YCAR80,*)ICOL(INBC/2+1:INBC)
      CALL WRITEDIR(NDIR,ICOL(INBC/2+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE
      !WRITE(YCAR80,*)ICOL(1:INBC)
      CALL WRITEDIR(NDIR,ICOL(1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ENDIF
      print*,(ZLEV(J),ICOL(J),J=1,INCL)
      print*,ICOL(INCL+1)

      WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
      WRITE(NLUOUT,*)ICOL(INCL+1)
! fin de la selection des couleurs par l'utilisateur
    ENDIF                        !====================================
!
! Trace des zones colorees
! ************************
    !IF(LMSKTOP .AND. LMARKER)THEN
    IF(LMARKER .AND. .NOT. LSPOT)THEN
    ! en etoiles colorees
      CALL GSMK(3)  ! asterisk is the type of marker
      DO JJ=1,NIJMAX
      DO JI=1,NIIMAX
	IF(PTAB(JI,JJ) /= XSPVAL)THEN
	  IF(PTAB(JI,JJ) < ZLEV(1))THEN
	    CALL GSPMCI(ICOL(1))
	  ELSE IF(PTAB(JI,JJ) >= ZLEV(INCL))THEN
	    CALL GSPMCI(ICOL(INCL+1))
	  ELSE
	    DO JK=1,INCL-1
	      IF(PTAB(JI,JJ) >= ZLEV(JK) .AND. &
		 PTAB(JI,JJ) < ZLEV(JK+1))THEN
		CALL GSPMCI(ICOL(JK+1))
		EXIT
              ENDIF
	    ENDDO
	  ENDIF
	  ZX=XZZXX(JI)
	  ZY=XZZXY(JJ)
	  CALL GPM(1,ZX,ZY)
	ENDIF
      ENDDO
      ENDDO

    ELSE IF (LSPOT .AND. .NOT. LMARKER) THEN
    ! en paves de couleur
      CALL  GSFAIS(1)  ! solid filling of the polygon
      ZEPX=(XZZXX(NIIMAX/2+1)-XZZXX(NIIMAX/2))*0.5
      ZEPY=(XZZXY(NIJMAX/2+1)-XZZXY(NIJMAX/2))*0.5
      print *,'LSPOT: taille differente de la maille?'
      print *,'       (n/N recommande pour trace de champs modeles)'
      print *,'       (avec contour: o/O/y/Y recommande pour trace d observations '
      print *,'        epaisseur du contour gere avec XLW1)'
      print *,'       (sans contour: a/A recommande pour trace d observations)'
      read(5,*) YREP
      CALL WRITEDIR(NDIR,YREP)
      IF(YREP=='o' .OR. YREP=='O' .OR. YREP=='y' .OR. YREP=='Y' .OR.&
         YREP=='a' .OR. YREP=='A'                               ) THEN
        ! essai de redimensionnement
        print *,'taille du pixel: NIMAX/nx et NJMAX/ny'
        print *,'indiquez nx et ny (2 entiers) ?'
        print *,'      si <=0 le defaut (50) est utilise'
        read(5,*) INBX,INBY
        CALL WRITEDIR(NDIR,INBX)
        CALL WRITEDIR(NDIR,INBY)
        IF(INBX<=0) INBX=50
        IF(INBY<=0) INBY=50
        ZEPX=ZEPX*NIIMAX/INBX ; ZEPY=ZEPY*NIJMAX/INBY
        ! contour en trait plein noir
        CALL DASHDB(65535)
      ENDIF
      DO JJ=1,NIJMAX
      DO JI=1,NIIMAX
        IF(PTAB(JI,JJ) /= XSPVAL)THEN
          IF(PTAB(JI,JJ) < ZLEV(1))THEN
            CALL GSFACI(ICOL(1))
          ELSE IF(PTAB(JI,JJ) >= ZLEV(INCL)) THEN
            CALL GSFACI(ICOL(INCL+1))
          ELSE
            DO JK=1,INCL-1
              IF(PTAB(JI,JJ) >= ZLEV(JK) .AND. &
                 PTAB(JI,JJ) < ZLEV(JK+1))THEN
                CALL GSFACI(ICOL(JK+1))
                EXIT
              ENDIF
            ENDDO
          ENDIF
          ZX5(1)=XZZXX(JI)-ZEPX ; ZY5(1)=XZZXY(JJ)-ZEPY
          ZX5(2)=XZZXX(JI)-ZEPX ; ZY5(2)=XZZXY(JJ)+ZEPY
          ZX5(3)=XZZXX(JI)+ZEPX ; ZY5(3)=XZZXY(JJ)+ZEPY
          ZX5(4)=XZZXX(JI)+ZEPX ; ZY5(4)=XZZXY(JJ)-ZEPY
          ZX5(5)=XZZXX(JI)-ZEPX ; ZY5(5)=XZZXY(JJ)-ZEPY
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
    ! Trace des surfaces couleurs
      CALL GSFAIS(1)
if(nverbia > 0)then
  print *,' ** image AV CALL ARINAM ',JPMAP
endif
      CALL ARINAM(IIMAP,JPMAP)
!     call mapbla(iimap)
      CALL CPCLAM(PTAB,ZRWRK,IWRK,IIMAP)
      CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,CCOLR)
      print *,' Col: MAP 1 6 5 ',IIMAP(1),IIMAP(6),IIMAP(5)
      CALL GSPLCI(1)
      CALL GSFAIS(0)
    ENDIF
!   CALL GSLN(1)
    ! Trace des valeurs (legende)
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
    CALL GSFAIS(1)
    CALL LBSETI('CBL',0)
    DO J=1,INCL
      YLLBS(J)=ADJUSTL(YLLBS(J))
    ENDDO
    IF(ZVR < .9)THEN
      CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(ZVR+.2,1.)-ZVR)/10.,MIN(ZVR+.2,1.),ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
    ELSE
      CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
!     CALL LBLBAR_FORDIACHRO(1,ZVR,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,2)
    ENDIF
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
    CALL TABCOL_FORDIACHRO
!   IF((LSUPER .AND. NSUPER == 1) .OR. .NOT.LSUPER)CALL TABCOL_FORDIACHRO
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
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==2) &
        CALL GSPLCI(2)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==2) &
        CALL GSTXCI(2)
      IF(NSUPER == 3)CALL GSPLCI(3)
      IF(NSUPER == 3)CALL GSTXCI(3)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==3) &
        CALL GSPLCI(4)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==3) &
        CALL GSTXCI(4)
      IF(NSUPER == 4)CALL GSPLCI(7)
      IF(NSUPER == 4)CALL GSTXCI(7)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==4) &
        CALL GSPLCI(3)
      IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==4) &
        CALL GSTXCI(3)
      IF(NSUPER > 4)CALL GSPLCI(NSUPER*2-1)
      IF(NSUPER > 4)CALL GSTXCI(NSUPER*2-1)
!!!!!!!! PROVI
!CALL FRSTPT(XXX(NIINF,NMGRID),XXY(NJINF,NMGRID))
!CALL VECTOR(XXX(NIINF,NMGRID),XXY(NJSUP,NMGRID))
!CALL VECTOR(XXX(NISUP,NMGRID),XXY(NJSUP,NMGRID))
!CALL VECTOR(XXX(NISUP,NMGRID),XXY(NJINF,NMGRID))
!CALL VECTOR(XXX(NIINF,NMGRID),XXY(NJINF,NMGRID))
!!!!!!!! PROVI
       ENDIF

    END IF
  ELSE                       !00000000000000000000000000000000000000000000

! Traits noir et blanc dans le cas de superpositions (LCOLAREA=.TRUE. et LCOLINE=.FALSE.)
! ********************************************************************************

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

ELSE IF( LGREY .AND. .NOT.LCOLAREA ) THEN
! **************************************************************
! Surfaces en grises ( LGREY=.TRUE.)
!  En cas de superpositions, obligatoirement le 1er dessin
! **************************************************************
  IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN        !000000000000000000
!
! Selection automatique des grises par le programme
! **************************************************
!
if(nverbia > 0)then
  print *,' ** image GREY av COLOR_FORDIACHRO(INCL+1,2) ,INCL',INCL
endif
     CALL COLOR_FORDIACHRO(INCL+1,2)
     WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' VALEURS:'
     IF(INCL /= 0)then
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('AIB',J)
         CALL CPSETI('AIA',J+1)
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
         ICOL(J)=J+2
         CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
       ENDDO
     ENDIF
     ICOL(INCL+1)=INCL+3
if(nverbia > 0)then
  print *,' ** image ICOL(INCL+1) ',ICOL(INCL+1)
endif
     IF(LCOLBR)THEN
       IF(ZLEV(MAX(1,INCL)) > ZLEV(1) .AND. ICOL(INCL+1) > ICOL(1))THEN
             print*,zlev(incl),zlev(1),icol(incl+1),icol(1)
         ALLOCATE(ICOL2(INCL+1))
         ICOL2(1:INCL+1)=ICOL(INCL+1:1:-1)
         ICOL(1:INCL+1)=ICOL2
!          ICOL(:)=ICOL2
         DEALLOCATE(ICOL2)
       END IF
     END IF
     IF(LCOLZERO)THEN
       ICOL(NCOLZERO)=0
     ENDIF
     WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)
     WRITE(NLUOUT,*)ICOL(INCL+1)
    ! Trace des zones grisees
    CALL GSFAIS(1)
    CALL ARINAM(IIMAP,JPMAP)
!   call mapbla(iimap)
    CALL CPCLAM(PTAB,ZRWRK,IWRK,IIMAP)
    CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,CCOLR)
    print *,' Grey: MAP 1 6 5 ',IIMAP(1),IIMAP(6),IIMAP(5)
    CALL GSPLCI(1)
    CALL GSFAIS(0)
!   CALL GSLN(1)
    ! Trace des valeurs (legende)
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
    CALL GSFAIS(1)
    CALL LBSETI('CBL',0)
    DO J=1,INCL
      YLLBS(J)=ADJUSTL(YLLBS(J))
    ENDDO
      IF(ZVR < .8999999)THEN
        print *,' ZVR < .9 ',ZVR
        CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(ZVR+.2,1.)-ZVR)/10.,MIN(ZVR+.2,1.),ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
      ELSE
        IF(INCL <= 8)THEN
	  if(nverbia >0)then
          print *,' INCL <= 8 ',INCL
	  endif
          CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB+(ZVT-ZVB)/4.,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
        ELSE
	  if(nverbia >0)then
          print *,' INCL > 8 ',INCL
	  endif
          CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,INCL+1,.15,1.,ICOL,1,YLLBS,INCL,1)
        ENDIF
!       CALL LBLBAR_FORDIACHRO(1,ZVR,1.,ZVB,ZVT,INCL+1,.15,1.,IHACH,2,YLLBS,INCL,2)
      ENDIF
      CALL GSFAIS(0)
!
! Definition de la couleur des isos (0 -> blanc sur papier; 1 -> noir sur papier)
    IF(LISOWHI)CALL GSPLCI(0)
    IF(LISOWHI)CALL GSTXCI(0)

  ELSE IF(LCOLINE)THEN       !00000000000000000000000000000000000000000000

! Traits couleur dans le cas de superpositions (LGREY=.TRUE. et LCOLINE=.TRUE.)
! **************************************************************************

! Modifs 220396
    CALL TABCOL_FORDIACHRO
!   IF((LSUPER .AND. NSUPER == 1) .OR. .NOT.LSUPER)CALL TABCOL_FORDIACHRO
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
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==2) &
        CALL GSPLCI(2)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==2) &
        CALL GSTXCI(2)
      IF(NSUPER == 3)CALL GSPLCI(3)
      IF(NSUPER == 3)CALL GSTXCI(3)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==3) &
        CALL GSPLCI(4)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==3) &
        CALL GSTXCI(4)
      IF(NSUPER == 4)CALL GSPLCI(7)
      IF(NSUPER == 4)CALL GSTXCI(7)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==4) &
        CALL GSPLCI(3)
      IF((LARROVL .OR. LGREY .OR. LHACH1) .AND. NSUPER ==4) &
        CALL GSTXCI(3)
      IF(NSUPER > 4)CALL GSPLCI(NSUPER*2-1)
      IF(NSUPER > 4)CALL GSTXCI(NSUPER*2-1)
!!!!!!!! PROVI
!CALL FRSTPT(XXX(NIINF,NMGRID),XXY(NJINF,NMGRID))
!CALL VECTOR(XXX(NIINF,NMGRID),XXY(NJSUP,NMGRID))
!CALL VECTOR(XXX(NISUP,NMGRID),XXY(NJSUP,NMGRID))
!CALL VECTOR(XXX(NISUP,NMGRID),XXY(NJINF,NMGRID))
!CALL VECTOR(XXX(NIINF,NMGRID),XXY(NJINF,NMGRID))
!!!!!!!! PROVI
       ENDIF

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

ELSE IF(LCOLINE)THEN    !+++++++++++++++++++++++++++++++++++++++++++++++++++++

! **********************************************
! Traits couleur   (LCOLAREA=.FALSE. et LCOLINE=.TRUE.)
! **********************************************
! Cas de superpositions
! *********************

! Modifs 220396
    CALL TABCOL_FORDIACHRO
!   IF((LSUPER .AND. NSUPER == 1) .OR. .NOT.LSUPER)CALL TABCOL_FORDIACHRO
! Modifs 260198
! IF(LSUPER)THEN             !............................................
  IF(LSUPER .AND. &          !............................................
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
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==2) &
      CALL GSPLCI(2)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==2) &
      CALL GSTXCI(2)
    IF(NSUPER == 3)CALL GSPLCI(3)
    IF(NSUPER == 3)CALL GSTXCI(3)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==3) &
      CALL GSPLCI(4)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==3) &
      CALL GSTXCI(4)
    IF(NSUPER == 4)CALL GSPLCI(7)
    IF(NSUPER == 4)CALL GSTXCI(7)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==4) &
      CALL GSPLCI(3)
    IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==4) &
      CALL GSTXCI(3)
    IF(NSUPER > 4)CALL GSPLCI(NSUPER*2-1)
    IF(NSUPER > 4)CALL GSTXCI(NSUPER*2-1)

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
       WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' COULEUR UNIQUE : ',ICOL(1)
       WRITE(NLUOUT,*)(ZLEV(J),J=1,INCL)
       ELSE
!Mars 2000

       CALL COLOR_FORDIACHRO(INCL,1)
       WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' VALEURS:'
       DO J=1,INCL
         CALL CPSETI('PAI',J)
         CALL CPSETI('CLC',J+2)
         CALL CPGETR('CLV',ZCLV)
         ZLEV(J)=ZCLV
         ICOL(J)=J+2
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
        ! WRITE(YCAR320,*)ICOL(1:INBC)
        ! YCAR320=ADJUSTL(YCAR320)
        ! ILENT=LEN_TRIM(YCAR320)
        ILENT=INBC*4
    IF(ILENT == 80 ) THEN
     ! YCAR320=TRIM(YCAR320)//' '
      ILENT=ILENT+1
    END IF
    IF(ILENT > 240 )THEN
     ! WRITE(YCAR80,*)ICOL(1:INBC/4)
      CALL WRITEDIR(NDIR,ICOL(1:INBC/4))
     ! WRITE(YCAR80,*)ICOL(INBC/4+1:INBC/2)
      CALL WRITEDIR(NDIR,ICOL(INBC/4+1:INBC/2))
     ! WRITE(YCAR80,*)ICOL(INBC/2+1:3*INBC/4)
      CALL WRITEDIR(NDIR,ICOL(INBC/2+1:3*INBC/4))
     ! WRITE(YCAR80,*)ICOL(3*INBC/4+1:INBC)
      CALL WRITEDIR(NDIR,ICOL(3*INBC/4+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE IF(ILENT > 160 )THEN
     ! WRITE(YCAR80,*)ICOL(1:INBC/3)
      CALL WRITEDIR(NDIR,ICOL(1:INBC/3))
     ! WRITE(YCAR80,*)ICOL(INBC/3+1:2*INBC/3)
      CALL WRITEDIR(NDIR,ICOL(INBC/3+1:2*INBC/3))
     ! WRITE(YCAR80,*)ICOL(2*INBC/3+1:INBC)
      CALL WRITEDIR(NDIR,ICOL(2*INBC/3+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE IF(ILENT > 80 )THEN
     ! WRITE(YCAR80,*)ICOL(1:INBC/2)
      CALL WRITEDIR(NDIR,ICOL(1:INBC/2))
     ! WRITE(YCAR80,*)ICOL(INBC/2+1:INBC)
      CALL WRITEDIR(NDIR,ICOL(INBC/2+1:INBC))
#ifdef RHODES
      CALL FLUSH(NDIR,ISTAF)
#else
      CALL FLUSH(NDIR)
#endif
    ELSE
     ! WRITE(YCAR80,*)ICOL(1:INBC)
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
           CALL GENFORMAT_FORDIACHRO(ZCLV,YLLBS(J))
         ENDDO
         WRITE(NLUOUT,*)' >>>>>>>IMAGE_FORDIACHRO VARIABLE : ',HTEXTE,' NB ISOC. : ',INCL,' VALEURS:'
         WRITE(NLUOUT,*)(ZLEV(J),ICOL(J),J=1,INCL)

    END IF                       !::::::::::::::::::::::::::::::::::::

!Mai 2009
      IF(LNOLBLBAR)THEN
      ELSE
!Mai 2009
!Mars 2000
       IF(LCOLISONE)THEN
       ELSE
!Mars 2000
       CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
       CALL GSFAIS(0)
       CALL SET(ZVR,1.,ZVB,ZVT,ZVR,1.,ZVB,ZVT,1)
       IF(INCL <= 1)THEN
	 ZINTERV=0.
       ELSE
         ZINTERV=(ZVT-ZVB-.009)/(INCL-1)
       ENDIF
       CALL GSCLIP(0)
       DO J=1,INCL
         YLLBS(J)=ADJUSTL(YLLBS(J))
         CALL GSPLCI(ICOL(J))
         CALL GSTXCI(ICOL(J))
	 IF(ZVR < .9 .AND. INCL < 25)THEN
           CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.015,0.,-1.)
         ELSEIF(ZVR < .9 .AND. INCL < 30 .AND. INCL >= 25)THEN
           CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.012,0.,-1.)
         ELSEIF(ZVR >= .95 )THEN
           CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.007,0.,-1.)
	 ELSE
           CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.009,0.,-1.)
	 ENDIF
!        CALL PLCHHQ(ZVR+(MIN(1.-ZVR,.2))/10.,ZVB+.004+(J-1)*ZINTERV,YLLBS(J),.007,0.,-1.)
       ENDDO
       CALL GSCLIP(1)
       CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!Mars 2000
       ENDIF
!Mars 2000
!Mai 2009
       ENDIF
!Mai 2009
       CALL GSTXCI(1)
       CALL GSPLCI(1)
       

  END IF                     !............................................

ELSE                    !+++++++++++++++++++++++++++++++++++++++++++++++++++++

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
        IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==2)CALL GSLN(1)
        IF(NSUPER ==3)CALL GSLN(3)
        IF((LARROVL .OR. LCOLAREA .OR. LHACH1) .AND. NSUPER ==3)THEN
          CALL GSLN(1)
          CALL GSLN(2)
          IF(LHACH2)CALL GSLN(1)
        ENDIF

      ELSE

        IF(NSUPER ==2)CALL GSLN(3)    ! line is a special dash type
        IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==2) &
          CALL GSLN(1)
        IF(NSUPER ==3)CALL GSLN(2)
        IF((LARROVL .OR. LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER ==3)THEN
          CALL GSLN(1)
          CALL GSLN(3)
          IF(LHACH2)CALL GSLN(1)
        ENDIF

      ENDIF

    END IF

  END IF                           !!!  Not an overlay case
!
END IF                  !+++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!*    3.3     High and low handling
!
IF (nverbia >=5) THEN
  print *,'image KNHI=',KNHI
END IF
SELECT CASE(KNHI)
    
  CASE(0)                           ! H + L   are displayed
! Test rajoute pour eviter la superposition de CONSTANT FIELD ici et ensuite
! avec le 2eme CPLBDR utile en cas de surfaces colorees
    IF(INCL /= 0)THEN
      CALL CPLBDR(PTAB,ZRWRK,IWRK)
    ENDIF
  CASE DEFAULT                      ! TO BE REVISED*********************
			            ! <0  --> no action (:-1 to be set)
			            ! >0  --> gridpoint value displayed
                                    !         (1: to be set)
END SELECT
!
!*     3.4   Effective contour drawing and line width selection
!    
IF(ZMIN == 999999. .AND. ZMAX == -999999.)THEN
  CALL CPSETC('CFT','CONSTANT FIELD - SPECIAL VALUE 999.')
ENDIF
GISO=LISO .AND. .NOT.(LSPOT .OR. LMARKER)
IF((LCOLAREA .AND. .NOT.GISO .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))&
  .OR.(LHACH1 .AND. .NOT.LISO .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))&
  .OR.(LGREY .AND. .NOT.LISO .AND. (.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1)))&
  .OR. (LHACH2 .AND. .NOT.LISO .AND. NSUPER == 2) &
  .OR. (LHACH3 .AND. .NOT.LISO .AND. NSUPER == 3) &
  .OR. (LHACH4 .AND. .NOT.LISO .AND. NSUPER == 4) ) THEN
ELSE
  CALL GSLWSC(XLWIDTH)
  IF(NSUPER == 2 .AND. LISOWHI2)THEN
    CALL GSLN(1)
    CALL GSPLCI(0)
    CALL GSTXCI(0)
  ELSE IF(NSUPER == 3 .AND. LISOWHI3)THEN
    CALL GSLN(1)
    CALL GSPLCI(0)
    CALL GSTXCI(0)
  ENDIF
  IF (nverbia >=5) THEN
    print *,'image av CPCLDR'
  END IF
  CALL CPCLDR(PTAB,ZRWRK,IWRK)
  ! message d erreur pour grd tableau: comment corriger ??
  !CPGIWS   50100 WORDS REQUESTED   50000 WORDS AVAILABLE
  IF (nverbia >=5) THEN
    print *,'image ap CPCLDR'
  END IF
END IF
IF((NSUPER == 2 .AND. LISOWHI2) .OR. (NSUPER == 3 .AND. LISOWHI3))THEN
! CALL GSPLCI(1)
  CALL GSTXCI(1)
ENDIF
IF(INCL == 0)THEN
  CALL CPLBDR(PTAB,ZRWRK,IWRK)
ENDIF

IF (nverbia >=5) THEN
  print *,'image avant CALL GSCLIP '
END IF
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,1)
CALL GSCLIP(0)

YTEM40(1:LEN(YTEM40))=' '
IF(NLOOPSUPER == 1)THEN
  CALL RESOLV_TIT('CTITVAR1',YTEM40)
ELSE IF(NLOOPSUPER == 2)THEN
  CALL RESOLV_TIT('CTITVAR2',YTEM40)
ELSE IF(NLOOPSUPER == 3)THEN
  CALL RESOLV_TIT('CTITVAR3',YTEM40)
ELSE IF(NLOOPSUPER == 4)THEN
  CALL RESOLV_TIT('CTITVAR4',YTEM40)
ELSE IF(NLOOPSUPER == 5)THEN
  CALL RESOLV_TIT('CTITVAR5',YTEM40)
ELSE IF(NLOOPSUPER == 6)THEN
  CALL RESOLV_TIT('CTITVAR6',YTEM40)
ELSE IF(NLOOPSUPER == 7)THEN
  CALL RESOLV_TIT('CTITVAR7',YTEM40)
ELSE IF(NLOOPSUPER == 8)THEN
  CALL RESOLV_TIT('CTITVAR8',YTEM40)
ENDIF
if(nverbia > 0)then
  print *,' image  CTITVAR ',YTEM40(1:LEN_TRIM(YTEM40))
endif

  IF(NSUPER < 4)THEN

    IF((LHACH1 .AND. NSUPER == 1) .OR. (LHACH2 .AND. NSUPER == 2) .OR. &
       (LHACH3 .AND. NSUPER == 3) .OR. (LHACH4 .AND. NSUPER == 4) ) THEN
    ELSE
      IF((LCOLAREA .AND. NSUPER > 1) .OR. &
         (.NOT.LCOLAREA  .AND. &
          .NOT.((LHACH1.OR.LHACH2) .AND. NSUPERDIA == 2)))THEN
        CALL GSLWSC(XLWIDTH)

	IF(YTEM40  /= ' ')THEN
        CALL FRSTPT(.95,.007+(NSUPER-1)*.017)
        CALL VECTOR(.95+.03,.007+(NSUPER-1)*.017)
	ENDIF

      ENDIF
    ENDIF

  ELSE

      IF((LCOLAREA .AND. NSUPER > 1) .OR. &
         (.NOT.LCOLAREA  .AND. &
          .NOT.((LHACH1.OR.LHACH2) .AND. NSUPERDIA == 2)))THEN

	IF(YTEM40  /= ' ')THEN
          CALL PLCHHQ(ZVLDEF+(NSUPER-4)*.25,ZVT+.01,ADJUSTL(CTIMEC(8:15)),.007,0.,-1.)
          CALL FRSTPT(ZVLDEF+(NSUPER-4)*.25+.08,ZVT+.01)
          CALL VECTOR(ZVLDEF+(NSUPER-4)*.25+.08+.03,ZVT+.01)
        ENDIF

      ENDIF

  ENDIF

CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
CALL GSLWSC(1.)
CALL GSLN(1)
CALL GSPLCI(1)
CALL GSTXCI(1)
IF(NSUPER == 1 .OR. .NOT.LSUPER .OR. (NSUPER == 2 .AND. LISOWHI2) .OR.  &
(NSUPER == 3 .AND. LISOWHI3))THEN
  IF(LCARTESIAN)THEN
    CALL DEFENETRE
  ELSE
    CALL BCGRD_FORDIACHRO(2)
  END IF
  IF(LXY)THEN
    CALL GSCLIP(0)
    CALL TRACEXY
  END IF
END IF
!------------------------------------------------------------------------------
!
!*     4.  TOPOGRAPHY MASKING WHEN PLOTTED LEVEL INTERCEPTS TERRAIN
!          --------------------------------------------------------
!
! Initialization of a topographic mask using 
! the NCAR "area" features (see NCAR manual)
!
if(nverbia > 0)then
  print *,' ** image AV CTYPHOR.EQ.Z'
endif
IF(CTYPHOR.EQ.'Z' .AND. (.NOT.LSUPER .OR. NSUPER == 1))THEN
  ZLREF=KLREF
  !                            ! If terrain higher -> a 888. mask value is forced
  DO J=NIINF,NISUP
     DO JJ=NJINF,NJSUP
        IF(ZLREF.LT.XXZS(J,JJ,NMGRID))PTAB(J-NIINF+1,JJ-NJINF+1)=888.
     ENDDO
  ENDDO
  !
  ICL=1                        ! A single contour will be drawn
  CALL CPSETI('CLS',0)         ! User provided contour value
  CALL CPSETI('HCF',1)         ! Area within contour will be hatched
  CALL CPSETC('CFT',' ')       ! No 'CONSTANT FIELD' message issued
  CALL CPSETI('NCL',ICL)       ! A single contour will be drawn
  CALL CPSETI('PAI',ICL)       ! A single contour will be drawn
  CALL CPSETI('AIA',ICL+1)     ! Area number where field values are > 888.
  CALL CPSETI('AIB',ICL)       ! Area number where field values are < 888. 
  CALL CPSETI('CLU',1)         ! Area without contour, if =1 unlabeled contour
  CALL CPSETR('SPV',0.)        ! Resets SPV, erases the special value setting
  CALL CPSETR('CLV',888.)      ! Value of the single contour drawn
  !
  ! As the topography-intercepted area has been set to 888., the rest of the
  ! field array is set to ZZSPVAL to hide it in the subsequent processing
  !
  ZZSPVAL=7777.
    WHERE(PTAB(:,:)/=888.)PTAB(:,:)=ZZSPVAL
    WHERE(PTAB(:,::2)==888.)PTAB(:,::2)=PTAB(:,::2)+1.E-3
  CALL CPSETR('SPV',ZZSPVAL)    ! Special value =  ZZSPVAL
  !
  ! Effective area computation and contour drawing
  !
  CALL ARINAM(IIMAP,JPMAP)                              ! Initialize areas
!   call mapbla(iimap)
if(nverbia > 0)then
  print *,' ** image AV CPRECT'
endif
  CALL CPRECT(PTAB,IM,IM,IL,ZRWRK,JPLRWK,IWRK,JPLIWK)   ! Initialize conpack
  CALL CPCLAM(PTAB,ZRWRK,IWRK,IIMAP)                    ! Contours terrain area
  CALL CPCLDR(PTAB,ZRWRK,IWRK)                          ! Contours outside field
  CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,SFILL)! Hatches 
  !                                                              !terrain area
END IF
!
!-----------------------------------------------------------------------------
!
!*    5.    COMPLETING THE PLOT
!           -------------------
!
!*    5.1   Page information labels
!
CALL GSCLIP(0)
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT

CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,1)
IF(CTYPHOR == 'T')THEN
  IF(.NOT.LTHSTAB)THEN
    CALL PLCHHQ(ZVL+.04,ZVT-.04,'*** UNSTABLE THETA ***',.011,0.,-1.)
  ENDIF
ELSE IF(CTYPHOR == 'E')THEN
  IF(.NOT.LTHSTAB)THEN
      CALL PLCHHQ(ZVL+.04,ZVT-.04,'*** VORTICITE NON MONOTONE ***',.011,0.,-1.)
  ENDIF
ELSE IF(CTYPHOR == 'V')THEN
  IF(.NOT.LTHSTAB)THEN
      CALL PLCHHQ(ZVL+.04,ZVT-.04,'*** FONCTION NON MONOTONE ***',.011,0.,-1.)
  ENDIF

ENDIF
IF(.NOT.LSUPER)THEN

! Modifs du 03/04/96
  IF(LEN_TRIM(HTEXTE) > 25)THEN                      !+++++++++++++
    ZSZTITVAR1=.009
  ELSE
    ZSZTITVAR1=.011
  ENDIF
  IF(XSZTITVAR1 /= 0.)THEN
    ZSZTITVAR1=XSZTITVAR1
  ENDIF
  IF(LCOLAREA .OR. LHACH1 .OR. LGREY)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL RESOLV_TIT('CTITVAR1',HTEXTE)
    IF(HTEXTE /= ' ')THEN
      CALL PLCHHQ(MAX(ZVR,.99),.007,HTEXTE,ZSZTITVAR1,0.,+1.)
!     CALL PLCHHQ(MAX(ZVR,.99),.007,HTEXTE,.011,0.,+1.)
    ENDIF

  ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL RESOLV_TIT('CTITVAR1',HTEXTE)
    IF(HTEXTE /= ' ')THEN
      CALL PLCHHQ(.93,.007,HTEXTE,ZSZTITVAR1,0.,+1.)
!     CALL PLCHHQ(.93,.007,HTEXTE,.011,0.,+1.)
    ENDIF

  ENDIF
  IF(LMINMAX)THEN
    CALL PCSETC('FC','/')
    CAll PLCHHQ(ZVR,ZVT+.03,YLBL,.009,0.,+1.)
    CALL PCSETC('FC',':')
  ENDIF

ELSE

  ZSZTITVAR=0.
  IF(NLOOPSUPER == 1)THEN
    CALL RESOLV_TIT('CTITVAR1',HTEXTE)
    IF(XSZTITVAR1 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR1
    ENDIF
  ELSE IF(NLOOPSUPER == 2)THEN
    CALL RESOLV_TIT('CTITVAR2',HTEXTE)
    IF(XSZTITVAR2 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR2
    ENDIF
  ELSE IF(NLOOPSUPER == 3)THEN
    CALL RESOLV_TIT('CTITVAR3',HTEXTE)
    IF(XSZTITVAR3 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR3
    ENDIF
  ELSE IF(NLOOPSUPER == 4)THEN
    CALL RESOLV_TIT('CTITVAR4',HTEXTE)
    IF(XSZTITVAR4 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR4
    ENDIF
  ELSE IF(NLOOPSUPER == 5)THEN
    CALL RESOLV_TIT('CTITVAR5',HTEXTE)
    IF(XSZTITVAR5 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR5
    ENDIF
  ELSE IF(NLOOPSUPER == 6)THEN
    CALL RESOLV_TIT('CTITVAR6',HTEXTE)
    IF(XSZTITVAR6 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR6
    ENDIF
  ELSE IF(NLOOPSUPER == 7)THEN
    CALL RESOLV_TIT('CTITVAR7',HTEXTE)
    IF(XSZTITVAR7 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR7
    ENDIF
  ELSE IF(NLOOPSUPER == 8)THEN
    CALL RESOLV_TIT('CTITVAR8',HTEXTE)
    IF(XSZTITVAR8 /= 0.)THEN
      ZSZTITVAR=XSZTITVAR8
    ENDIF
  ENDIF
if(nverbia > 0)then
  print *,' image  CTITVAR ',HTEXTE(1:LEN_TRIM(HTEXTE))
endif

! Modifs du 03/04/96 NON NON REFLECHIR EN CAS DE SUPERPOSITIONS
  IF(NSUPER < 4)THEN

    IF(NSUPER == 1)ZSC=999.
    IF(LEN_TRIM(HTEXTE) > 25)THEN                      !+++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF((LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER == 1)THEN !000000000000

        IF(HTEXTE /= ' ')THEN
	  IF(ZSZTITVAR /= 0.)THEN
            CALL PLCHHQ(MAX(ZVR,.99),.007+(NSUPER-1)*.017,HTEXTE,ZSZTITVAR,0.,+1.)
	  ELSE
            CALL PLCHHQ(MAX(ZVR,.99),.007+(NSUPER-1)*.017,HTEXTE,.005,0.,+1.)
	  ENDIF
        ENDIF

      ELSE                                       !00000000000000000000
	IF((LHACH2 .AND. NSUPER == 2) .OR. (LHACH3 .AND. NSUPER == 3) .OR. &
           (LHACH4 .AND. NSUPER == 4) ) THEN

          IF(IHT == 1)THEN
            IF(HTEXTE /= ' ')THEN
	      IF(ZSZTITVAR /= 0.)THEN
                CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,ZSZTITVAR,0.,-1.)
	      ELSE
                CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,.007,0.,-1.)
              ENDIF
            ENDIF
          ELSE
            IF(HTEXTE /= ' ')THEN
	      IF(ZSZTITVAR /= 0.)THEN
                CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,ZSZTITVAR,0.,-1.)
	      ELSE
                CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,.005,0.,-1.)
              ENDIF
            ENDIF
          ENDIF
	ELSE

          IF(HTEXTE /= ' ')THEN
	    IF(ZSZTITVAR /= 0.)THEN
              CALL PLCHHQ(.93,.007+(NSUPER-1)*.017,HTEXTE,ZSZTITVAR,0.,+1.)
	    ELSE
              CALL PLCHHQ(.93,.007+(NSUPER-1)*.017,HTEXTE,.005,0.,+1.)
	    ENDIF
	  ENDIF

        ENDIF
      ENDIF                                      !0000000000000000000

      ZSC=.005

    ELSE                                               !+++++++++++++

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF((LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER == 1)THEN

        IF(HTEXTE /= ' ')THEN
	  IF(ZSZTITVAR /= 0.)THEN
            CALL PLCHHQ(MAX(ZVR,.99),.007+(NSUPER-1)*.017,HTEXTE,ZSZTITVAR,0.,+1.)
	  ELSE
            CALL PLCHHQ(MAX(ZVR,.99),.007+(NSUPER-1)*.017,HTEXTE,.007,0.,+1.)
	  ENDIF
	ENDIF

      ELSE

	IF((LHACH2 .AND. NSUPER == 2) .OR. (LHACH3 .AND. NSUPER == 3) .OR. &
           (LHACH4 .AND. NSUPER == 4))THEN

          IF(HTEXTE /= ' ')THEN
	    IF(ZSZTITVAR /= 0.)THEN
              CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,ZSZTITVAR,0.,-1.)
	    ELSE
              CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,.005,0.,-1.)
	    ENDIF
	  ENDIF

	ELSE

          IF(HTEXTE /= ' ')THEN
	    IF(ZSZTITVAR /= 0.)THEN
              CALL PLCHHQ(.93,.007+(NSUPER-1)*.017,HTEXTE,ZSZTITVAR,0.,+1.)
	    ELSE
              CALL PLCHHQ(.93,.007+(NSUPER-1)*.017,HTEXTE,.007,0.,+1.)
	    ENDIF
	  ENDIF

        ENDIF
      ENDIF

    ENDIF                                              !+++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF((LCOLAREA .OR. LHACH1 .OR. LGREY) .AND. NSUPER == 1)THEN

      IF(HTEXTE /= ' ')THEN
        CALL PLCHHQ(1.-((LEN_TRIM(HTEXTE)+5)*.007),.007+(NSUPER-1)*.017,CTIMEC(8:15),.007,0.,+1.)
      ENDIF

    ELSE

      IF((LHACH2 .AND. NSUPER == 2) .OR. (LHACH3 .AND. NSUPER == 3) .OR. &
         (LHACH4 .AND. NSUPER == 4))THEN
!!!!!!! REFLECHIR
!       CALL PLCHHQ(ZD,ZVT+.04,HTEXTE,.005,0.,-1.)
      ELSE
        IF(HTEXTE /= ' ')THEN
          CALL PLCHHQ(.93-((LEN_TRIM(HTEXTE)+4)*.007),.007+(NSUPER-1)*.017,CTIMEC(8:15),.007,0.,+1.)
        ENDIF
      ENDIF

    ENDIF

    IF(LMINMAX)THEN
      CALL PCSETC('FC','/')
      CAll PLCHHQ(ZVRDEF,ZVT+.01+(NSUPER-1)*.02,YLBL,.007,0.,+1.)
      CALL PCSETC('FC',':')
    ENDIF

  ELSE

    IF(ZSC /= 999.)THEN
      IF(HTEXTE /= ' ')THEN
        CALL PLCHHQ(ZVLDEF+(NSUPER-4)*.25,ZVT+.03,HTEXTE,ZSC,0.,-1.)
      ENDIF
    ELSE
      IF(HTEXTE /= ' ')THEN
        CALL PLCHHQ(ZVLDEF+(NSUPER-4)*.25,ZVT+.03,HTEXTE,.007,0.,-1.)
      ENDIF
    ENDIF

  ENDIF


END IF
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
CALL GSLWSC(1.)
CALL GSLN(1)
CALL GSPLCI(1)
CALL GSTXCI(1)
! Oct 99

!IF(LFACTIMP)THEN
! CALL FACTIMP
!ENDIF
! Oct 99
if(nverbia > 0)then
  print *,' ** image AV NOT LSUPER'
endif
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN
! Mars 2000
IF(LFACTIMP)THEN
  CALL FACTIMP
ENDIF
! Modifs for diachro
! Remodifs le 170596
! Titres en X
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
!     CALL PLCHHQ((ZVL+ZVR)/2.-ZVB/2.,ZVB/2.,YTEM,XSZTITXM,0.,-1.)
    ELSE
      CALL PLCHHQ((ZVL+ZVR)/2.,ZVB-MIN(ZVB/2.,.05),YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
!     CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
!     CALL PLCHHQ((ZVL+ZVR)/2.-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXR',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXR',YTEM)
    IF(XSZTITXR /= 0.)THEN
      CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/2.,.05),YTEM,XSZTITXR,0.,-1.)
!     CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,XSZTITXR,0.,-1.)
    ELSE
      CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/2.,.05),YTEM,.008,0.,-1.)
!     CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
! Titres en Y
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYM',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYB',ZVL,ZVR,ZVB,ZVT,YTEM)
! Titres  TOP
  YTEM(1:LEN(YTEM))=' '
  ZXPOSTITT2=.002
  ZXYPOSTITT2=.95
  IF(XPOSTITT2 /= 0.)THEN
    ZXPOSTITT2=XPOSTITT2
  ENDIF
  IF(XYPOSTITT2 /= 0.)THEN
    ZXYPOSTITT2=XYPOSTITT2
  ENDIF
  CALL RESOLV_TIT('CTITT2',YTEM)
  IF(YTEM /= ' ')THEN
    IF(XSZTITT2 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,XSZTITT2,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YTEM,XSZTITT2,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,.008,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
  ZXPOSTITT3=.002
  ZXYPOSTITT3=.93
  IF(XPOSTITT3 /= 0.)THEN
    ZXPOSTITT3=XPOSTITT3
  ENDIF
  IF(XYPOSTITT3 /= 0.)THEN
    ZXYPOSTITT3=XYPOSTITT3
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITT3',YTEM)
  IF(YTEM /= ' ')THEN
    IF(XSZTITT3 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,XSZTITT3,0.,-1.)
!     CALL PLCHHQ(0.002,0.93,YTEM,XSZTITT3,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,.008,0.,-1.)
!     CALL PLCHHQ(0.002,0.93,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF

! Titre N1 BOTTOM
  ZXPOSTITB1=.002
  ZXYPOSTITB1=.005
  IF(XPOSTITB1 /= 0.)THEN
    ZXPOSTITB1=XPOSTITB1
  ENDIF
  IF(XYPOSTITB1 /= 0.)THEN
    ZXYPOSTITB1=XYPOSTITB1
  ENDIF
  CALL RESOLV_TIT('CTITB1',CLEGEND)
  IF(CLEGEND /= ' ')THEN
    IF(XSZTITB1 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,CLEGEND,XSZTITB1,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,CLEGEND,.007,0.,-1.)
    ENDIF
  ENDIF
! Titre N3 BOTTOM
  ZXPOSTITB3=.002
  ZXYPOSTITB3=.045
  IF(XPOSTITB3 /= 0.)THEN
    ZXPOSTITB3=XPOSTITB3
  ENDIF
  IF(XYPOSTITB3 /= 0.)THEN
    ZXYPOSTITB3=XYPOSTITB3
  ENDIF
  IF(LCNCUM .OR. LCNSUM)THEN
    CALL RESOLV_TIT('CTITB3',CTIMECS)
    IF(CTIMECS /= ' ')THEN
      IF(XSZTITB3 /= 0.)THEN
        CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTIMECS,XSZTITB3,0.,-1.)
      ELSE
        CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTIMECS,.009,0.,-1.)
      ENDIF
    ENDIF
  ELSE
    IF(LMINUS .OR. LPLUS)THEN
      IF(.NOT.LTITDEFM .AND. CTITB3MEM /= 'DEFAULT' .AND. &
      CTITB3MEM /= 'default' .AND. CTITB3MEM /= 'DEFAUT' .AND. &
      CTITB3MEM /= 'defaut')THEN
! Il ne faut pas mettre l'instruction suivante
!       CALL RESOLV_TIT('CTITB3',CTITB3MEM)
	  if(nverbia > 0)then
	  print *,' image  CTITB3MEM ',CTITB3MEM(1:LEN_TRIM(CTITB3MEM))
	  endif
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
! ******************** 200697 ***************
          CALL RESOLV_TIT('CTITB3',CTITB3)
	  if(nverbia > 0)then
	  print *,' image  CTITB3 ',CTITB3(1:LEN_TRIM(CTITB3))
	  endif
          IF(CTITB3 /= ' ')THEN
            IF(XSZTITB3 /= 0.)THEN
              CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3,XSZTITB3,0.,-1.)
            ELSE
              CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3,.009,0.,-1.)
	    ENDIF
          ENDIF
      ENDIF
! ******************** 200697 ***************
    ELSE

      IF(CSTORAGE_TYPE /= 'PG')THEN
! NBPMT=nb de + et -
	IF(NBPMT == 0)THEN
          YTEM(1:LEN(YTEM))=' '
	  YTEM=CTIMEC
	  YTEM=ADJUSTL(YTEM)
          CALL RESOLV_TIT('CTITB3',YTEM)
	  if(nverbia > 0)then
	  print *,' image LEN et CTIMEC ',LEN(CTIMEC),CTIMEC
	  print *,' image LEN et YTEM ',LEN(YTEM),YTEM
	  endif
          IF(YTEM/= ' ')THEN
            IF(XSZTITB3 /= 0.)THEN
              CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,XSZTITB3,0.,-1.)
            ELSE
              CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,.009,0.,-1.)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

    ENDIF
  ENDIF
! Titre N2 BOTTOM
  ZXPOSTITB2=.002
  ZXYPOSTITB2=.025
  IF(XPOSTITB2 /= 0.)THEN
    ZXPOSTITB2=XPOSTITB2
  ENDIF
  IF(XYPOSTITB2 /= 0.)THEN
    ZXYPOSTITB2=XYPOSTITB2
  ENDIF
  CALL RESOLV_TIT('CTITB2',CLEGEND2)
  IF(CLEGEND2 /= ' ')THEN
    IF(XSZTITB2 /= 0.)THEN
      CALL PLCHHQ(0.002,0.025,CLEGEND2,XSZTITB2,0.,-1.)
    ELSE
      CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
    ENDIF
  ENDIF
! Titre N1 TOP
  ZXPOSTITT1=.002
  ZXYPOSTITT1=.98
  IF(XPOSTITT1 /= 0.)THEN
    ZXPOSTITT1=XPOSTITT1
  ENDIF
  IF(XYPOSTITT1 /= 0.)THEN
    ZXYPOSTITT1=XYPOSTITT1
  ENDIF
  WRITE(YPLANH,1001)NIINF,NISUP,NJINF,NJSUP
  CALL RESOLV_TIT('CTITT1',YPLANH)
  IF(YPLANH /= ' ')THEN
    IF(XSZTITT1 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YPLANH,XSZTITT1,0.,-1.)
!     CALL PLCHHQ(0.002,0.98,YPLANH,XSZTITT1,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YPLANH,.012,0.,-1.)
!     CALL PLCHHQ(0.002,0.98,YPLANH,.012,0.,-1.)
    ENDIF
  ENDIF
  IF(LDATFILE)CALL DATFILE_FORDIACHRO
ENDIF
!
1001 FORMAT('HORIZONTAL SECTION NIINF=',I4,' NISUP=',I4,' NJINF=',I4,' NJSUP=',I4)
!
CALL GSCLIP(1)
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!
!*    5.2   NCAR parameter reset
!
CALL CPSETI('CLS',16)
CALL CPRSET
CALL GSLN(1)
!
!--------------------------------------------------------------------------------
!
!*    6.    EXIT
!           ----
!
RETURN
END SUBROUTINE IMAGE_FORDIACHRO
