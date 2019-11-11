!     ######spl
      SUBROUTINE BCGRD_FORDIACHRO(K)
!     ##############################
!
!!****  *BCGRD* - Displays a cartographic background in horizontal mode
!!
!!    PURPOSE
!!    -------
!       Displays a cartographic background for horizontal cross-section
!     contour or arrow maps when the cartographic projection option is 
!     active.  
!       The geographical display window is defined, a grid of latitude-
!     longitude lines, a set of continental/state outlines and, optionaly,
!     a series of landmarks, are plotted on this background. 
!
!!**  METHOD
!!    ------
!!     
!!     The conformal projection routines of MODE_GRIDPROJ are used to
!!    compute the latitude-longitude coordinates of the display box.
!!    Next, the NCAR Ezmap projection parameters are set up to 
!!    correspond to the Meso-NH projection, and a grid of latitude-
!     longitude lines, a set of continental/state outlines and, optionaly,
!     a series of landmarks, are plotted as an overlay on the current map.
!!
!!    EXTERNAL
!!    --------
!!
!!      MAPSTI ! set an NCAR parameter to a valuei, type  INTEGER   !
!!      MAPSTC ! (cartographic projection package)        CHARACTER !
!!      MAPROJ   selects a type cartographic projection             !
!!      MAPDRW   draws a map as specified by the user parameter     !
!!               choice                                             !
!!      MAPIT    draws a polyline on a map, using map coordinates   ! NCAR
!!      MAPIQ    terminates a line drawn by MAPIT                   !
!!      MAPSET   defines the plot window using map coordinates      !
!!      MAPTRN   projects a point onto a geographic map using       !
!!               latitude-longitude to locate the point             !
!!                                                                  !
!!      PWRITX   prints a text                                      !
!!      LABMOD   defines the axes label formats (paired with PERIM) !Routines 
!!      GRIDAL   draws grid lines and labels                        !
!!      PERIM    draws an unlabeled plot perimeter                  !
!!      SET      defines the plot window and viewport using user    !
!!               and normalized NCAR coordinates                    !
!!      GETSET   retrieves the NCAR and user coordinate definitions !
!!      PLCHHQ   high quality printing facility                     !
!!      GSCLIP   clips the plot using the window limits             !
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID      : Current MESO-NH grid indicator
!!
!!      Module MODE_GRIDPROJ:  packages a set of cartographic
!!                             module-procedures
!!         SM_LATLON   : to compute geographic  from conformal (cartographic)
!!                       cartesian coordinates;
!!         SM_XYHAT    : to compute conformal (cartographic) cartesian from
!!                       geographic coordinates;
!!         LATREF2     : to compute the second reference latitude
!!                       in the case of Lambert conformal projection
!!
!!      Module MODD_COORD      : declares gridpoint coordinates (TRACE use)
!!         XXX,XXY     : coordinate values FOR ALL  the MESO-NH grids
!!
!!      Module MODD_GRID1      : declares grid variables (Model module)
!!         XXHAT, XYHAT  : x, y cartographic coordinates of the model grid
!!         XLONOR,XLATOR : longitude and latitude of the (1,1,1) point of
!!                         the model mass grid
!!
!!      Module MODD_GRID    : declaration of grid variables for all models
!!         XLON0,XLAT0 : reference longitude and latitude for the conformal
!!                       projection
!!         XBETA,XRPK  : rotation angle and projection parameter for the 
!!                       conformal projection
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist
!!                         (former NCAR common)
!!         NIFDC   : Coastline data style (0 none, 1 NCAR, 2 IGN)
!!         NLPCAR  : Number of land-mark points to be plotted
!!         XLONCAR :  Longitude of land-mark points
!!         XLATCAR :  Latitude  of land-mark points
!!
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist (former PARA common)
!!       Module MODD_DIM1 : contains dimensions of data arrays
!!         NIINF, NISUP : lower and upper bounds of arrays
!!                        to be plotted in x direction
!!         NJINF, NJSUP : lower and upper bounds of arrays
!!                        to be plotted in y direction
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
!!      Updated   PM   12/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODD_RADAR
USE MODE_GRIDPROJ
USE MODD_COORD
USE MODD_MASK3D
USE MODD_TRAJ3D
USE MODD_RESOLVCAR
USE MODD_GRID1
USE MODD_GRID
USE MODD_CTL_AXES_AND_STYL
USE MODN_NCAR
USE MODN_PARA
USE MODI_CREATLINK
USE MODI_WRITEDIR

IMPLICIT NONE

COMMON/EPAISCONT/ZLWCONT
COMMON/FDC/IFDC

INTEGER :: K
!
!*       0.1   Local variables
!
REAL :: ZLWCONT
REAL :: ZLAT2, ZLAT, ZLON
REAL,SAVE :: ZZLAT, ZZLON
REAL :: ZPL1, ZPL2, ZPL3, ZPL4
REAL :: ZX1, ZX2, ZY1, ZY2, ZXX1, ZXX2, ZYY1, ZYY2
REAL :: ZU, ZV, ZSZ, ZPOS, ZCENT
REAL :: ZI, ZJ, ZX, ZY
INTEGER :: ICONVI, ICONVJ
REAL :: ZXMIN, ZXMAX,ZYMIN, ZYMAX
REAL :: ZWIDTH
CHARACTER(LEN=40),SAVE :: YCAR40=' '
CHARACTER(LEN=80),SAVE :: YCAR80=' '
CHARACTER(LEN=1)       :: YSYMB
CHARACTER(LEN=20)      :: YNOM
CHARACTER(LEN=10) :: FORMAX, FORMAY

INTEGER :: JIP, IT, IDUM, IRPK, JLPCAR,JIJCAR, J,IIT
INTEGER :: IERR, IPOS, ICOLS, ICOLN
INTEGER :: IFDC
!!!!!!!!!!!!!! Modif VD (29/10/2003)
INTEGER :: IDOT,IPT,IDOT0,IPT0,JLOOP
REAL, DIMENSION(200000) :: ZZU,ZZV,ZZU0,ZZV0
!!!!!!!!!!!!!! fin Modif VD
LOGICAL :: GIND,GCONF
!
!-------------------------------------------------------------------------------
!
!*       1.   SETS CARTOGRAPHIC PROJECTION AND DRAWS BACKGROUND MAP
!             -----------------------------------------------------
!
!
!*       1.1  If Lambert case, computes the second reference latitude
!            (required by the NCAR framework for Lambert) 
!
IF(L2CONT)THEN
  IFDC=NIFDC
ELSE
  IF(K == 1)THEN
    IFDC=0
  ELSE
    IFDC=NIFDC
  ENDIF
ENDIF
!!!!IFDC=NIFDC
IF(ABS(XRPK).GT.0..AND.ABS(XRPK).LT.1.)THEN
  IF(NVERBIA >= 5)THEN
    print *,' bcgrd XLAT0,XRPK ',XLAT0,XRPK
  ENDIF
  ZLAT2=LATREF2(XLAT0,XRPK)
  IF(NVERBIA >= 5)THEN
    print *,' bcgrd ZLAT2 ',ZLAT2
  ENDIF
ENDIF
!
!*       1.2  Convert display window diagonal to cartographic coordinates
!
! (The main diagonal of the displayed domain is given by Meso-NH
!  indexes NIINF-NJINF, NISUP-NJSUP)
!
!ZXMIN=100000.
ZXMIN=XXX(NIINF,NMGRID)
ZYMIN=XXY(NJINF,NMGRID)
!ZXMAX=2500000.
ZXMAX=XXX(NISUP,NMGRID)
ZYMAX=XXY(NJSUP,NMGRID)
IF(NVERBIA >= 2)THEN
  print *,' ** bcg NIINF,NJINF,NMGRID,NISUP,NJSUP ',NIINF,NJINF,NMGRID,NISUP,NJSUP
ENDIF
!
CALL SM_LATLON_S(XLATORI,XLONORI,ZXMIN,ZYMIN,ZPL1,ZPL2)
CALL SM_LATLON_S(XLATORI,XLONORI,ZXMAX,ZYMAX,ZPL3,ZPL4)
IF(NVERBIA >= 2)THEN
  print *,' ZXMIN,ZYMIN,ZXMAX,ZYMAX ',ZXMIN,ZYMIN,ZXMAX,ZYMAX
  print *,' XLATORI,XLONORI,ZPL1,ZPL2,ZPL3,ZPL4 ',XLATORI,XLONORI,ZPL1,ZPL2,ZPL3,ZPL4
  print *,' XLATO,XLONO ',XLAT0,XLON0
ENDIF
!
!*       1.3   Selects a standard NCAR continental/state outline mode
!*             and visual details
!
! -> NCAR default : call mapstc('OU','PO')
! ->     None      : call mapstc('OU','NO')
!
!IF (NIFDC.NE.1)THEN
IF (NIFDC.EQ.1 .OR. NIFDC.EQ.3)THEN
  CALL MAPSTC('OU','PO')
ELSE
  CALL MAPSTC('OU','NO')
ENDIF
!
CALL MAPSTI('DO',0)        ! Solid coastlines
!CALL MAPSTI('DO',1)        ! Dotted coastlines
CALL MAPSTI('RE',10000)    ! Plotter resolution
CALL MAPSTI('DL',0)        ! MAPIT draws solid lines
!CALL MAPSTI('DL',1)        ! MAPIT draws dotted lines
!CALL MAPSTI('GR',NIGRNC)   ! Grid spacing in degrees
if(nverbia > 0)then
  print *,' **bcgrd AV CALL MAPSTI(GR,0)'
endif
IF(K == 1)THEN
  CALL MAPSTI('GR',0)   ! Grid spacing in degrees
ELSE IF(K == 2)THEN
  IF(LANIMK )THEN
  ELSE
  CALL MAPSTI('GR',NIGRNC)   ! Grid spacing in degrees
  ENDIF
ENDIF
!
!*       1.4   Selects NCAR cartographic projection
!
IRPK=2
IF(XRPK.EQ.0.)IRPK=0
! Oct 99 Pole Sud Proj. stereog.
IF(ABS(XRPK).EQ.1.)IRPK=1
! Oct 99 Pole Sud Proj. stereog.
!IF(XRPK.EQ.1.)IRPK=1
!write(0,*)' BCGRD IRPK ',IRPK
!
SELECT CASE(IRPK)
  CASE(0)  
    CALL MAPROJ('ME',0.,XLON0,XBETA)               ! Mercator
  CASE(1)
    CALL MAPROJ('ST',90.,XLON0,-XBETA)             ! Polar Stereographic
! Oct 99 Pole Sud Proj. stereog.
!  BESOIN DE VERIFIER si dans ce cas on met XBETA ou -XBETA
    IF(XRPK < 0.)CALL MAPROJ('ST',-90.,XLON0,-XBETA)
! Oct 99 Pole Sud Proj. stereog.
  CASE DEFAULT
    CALL MAPROJ('LC',XLAT0,XLON0+XBETA/XRPK,ZLAT2) ! Lambert
END SELECT
!
!*       1.5   Sets map transformation, map display window
!*             and draws lat-lon grid
!
IF(LVPTUSER)THEN
  CALL MAPPOS(XVPTL,XVPTR,XVPTB,XVPTT)
ELSE
  CALL MAPPOS(.05,.95,.05,.95)
ENDIF
CALL MAPSET('CO',ZPL1,ZPL2,ZPL3,ZPL4)
IF(XLWCONT /= 0.)THEN
  ZLWCONT=XLWCONT
ELSE
  ZLWCONT=5.
ENDIF
!
! Pour V4.1.1 A la place de CALL MAPDRW a mettre en commentaire
! Non c'est fait EN PRINCIPE dans MAPDRW qui est inclus dans le fichier frame
!CALL MAPINT
!CALL MAPGRD
!CALL MAPLBL
!CALL MPLNDR('Earth..1',3)
if(nverbia > 0)then
  print *,' **bcgrd AV CALL MAPDRW'
endif
CALL MAPDRW
!
!*      1.6    Use of non-NCAR coastline data sets if available 
!*             (ex. IGN ones) on fortran unit 1
!
! NOTICE: The use of fortran unit 1 here does not
!         fit Meso-NH file access norm
!
IF((NIFDC.EQ.2 .OR. NIFDC.EQ.3) .AND. K.EQ.2)THEN
  IF(YCAR40(1:LEN(YCAR40)) == ' ')THEN
    print *,'ENTREZ le nom du fichier des contours (geograp. ou polit....) '
    !print *,' avec un PATH ABSOLU  (40 caracteres maximum) et entre quotes'
    print *,' entre quotes (40 caracteres maximum)'
    READ(5,*)YCAR40
    YCAR40=ADJUSTL(YCAR40)
    YCAR80(1:1)="'"
    YCAR80(2:LEN_TRIM(YCAR40)+1)=YCAR40(1:LEN_TRIM(YCAR40))
    YCAR80(LEN_TRIM(YCAR40)+2:LEN_TRIM(YCAR40)+2)="'"
    !WRITE(NDIR,'(A80)')YCAR80
    CALL WRITEDIR(NDIR,YCAR80)
    CALL CREATLINK('DIRFDC',YCAR40(1:LEN_TRIM(YCAR40)),'CREAT',NVERBIA)
! print *,YCAR40
  ENDIF
  OPEN(1,FILE=YCAR40(1:LEN_TRIM(YCAR40)),FORM='FORMATTED',STATUS='OLD')  ! Opens coastline file
!  OPEN(1,FILE='/u/m/mrmh/mrmh005/mesonh/data/cotign')  ! Opens coastline file
  CALL GSCLIP(0)
  CALL GQLWSC(IERR,ZWIDTH)
  IF(XLWCONT /= 0.)THEN
    ZLWCONT=XLWCONT
  ELSE
    ZLWCONT=4.
  ENDIF
  CALL GSLWSC(ZLWCONT)
!
!!!!!!!!! MODIF VD TO introduce dashed lines with NIFDC=2  (29/10/2003)
! Initial coordinate transformation saved
  CALL GETSET(ZX1,ZX2,ZY1,ZY2,ZXX1,ZXX2,ZYY1,ZYY2,IDUM)
! Initial coordinate transformation restored
  CALL SET(ZX1,ZX2,ZY1,ZY2,ZXMIN,ZXMAX,ZYMIN,ZYMAX,IDUM)
   IPT=0
   IPT0=0
   IDOT=838860 ! dashed pattern used for dashed lines (IT=2 or 3)
   IDOT0=65535 ! dashed pattern used for solid lines (IT=0 or 1)
    DO JIP=1,200000
      READ(1,*,END=50)ZLAT,ZLON,IT                ! Reads coastline file
      IF(JIP == 1)print *,' 1er enr. ',ZLAT,ZLON,IT 
!     IF(ABS(ZZLAT-ZLAT) > .2 .OR. ABS(ZZLON-ZLON) > .2)THEN
!       print *,'ZZLAT,ZLAT,ZZLON,ZLON ',ZZLAT,ZLAT,ZZLON,ZLON
!       IT=0
!       CALL MAPIT(ZLAT,ZLON,IT)             ! Draws IGN one coastline point
!     ELSE
!       CALL MAPIT(ZLAT,ZLON,IT)             ! Draws IGN one coastline point
!     ENDIF
      !ZZLAT=ZLAT
      !ZZLON=ZLON
      CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZU,ZV)
!
      IF (IT==2 .OR. IT==3) THEN
        IF (IT==2) THEN
          IF (IPT>0) THEN
           CALL DASHDB(IDOT) 
           CALL CURVED(ZZU,ZZV,IPT)
          ENDIF
          IPT=0
          IF ((ZU>= ZXMIN).AND.(ZU<=ZXMAX).AND.(ZV>=ZYMIN).AND.(ZV<=ZYMAX)) THEN
            IPT=IPT+1
            ZZU(IPT)=ZU
            ZZV(IPT)=ZV 
          ENDIF
        ELSE
          IF ((ZU>= ZXMIN).AND.(ZU<=ZXMAX).AND.(ZV>=ZYMIN).AND.(ZV<=ZYMAX)) THEN
            IPT=IPT+1
            ZZU(IPT)=ZU
            ZZV(IPT)=ZV 
          END IF
        ENDIF
      ELSE
 
        IF (IT==0) THEN  ! begin of the definition of the 
          IF (IPT0>0) THEN
            CALL DASHDB(IDOT0) 
            CALL CURVED(ZZU0,ZZV0,IPT0)
          ENDIF
          IPT0=0
          IF ((ZU>= ZXMIN).AND.(ZU<=ZXMAX).AND.(ZV>=ZYMIN).AND.(ZV<=ZYMAX)) THEN
            IPT0=IPT0+1
            ZZU0(IPT0)=ZU
            ZZV0(IPT0)=ZV 
          ENDIF
        ELSE
          IF ((ZU>= ZXMIN).AND.(ZU<=ZXMAX).AND.(ZV>=ZYMIN).AND.(ZV<=ZYMAX)) THEN
            IPT0=IPT0+1
            ZZU0(IPT0)=ZU
            ZZV0(IPT0)=ZV 
          END IF
        ENDIF
      ENDIF
!
    ENDDO
50 CONTINUE
! finish to draw the last curves :  
      print *,' Dernier enr. ',ZLAT,ZLON,IT 
    !CALL MAPIQ
    IF (IPT>0) THEN
      CALL DASHDB(IDOT)
      CALL CURVED(ZZU,ZZV,IPT)
    ENDIF
    IF (IPT0>0) THEN
      CALL DASHDB(IDOT0)
      CALL CURVED(ZZU0,ZZV0,IPT0)
    ENDIF
!!!!!!!!!!!!!!!!!!! fin modif VD
  CALL GSCLIP(1)                           ! Clipping of extra coastline
  CLOSE(1)
  CALL GSLN(1)                             ! restore solid line
  CALL GSLWSC(ZWIDTH)
ENDIF
!
!*      1.7    Formats and write Map axes with appropriate labels
!*             and axes scale labels 
!
! Initial coordinate transformation saved
CALL GETSET(ZX1,ZX2,ZY1,ZY2,ZXX1,ZXX2,ZYY1,ZYY2,IDUM)
! Sets NCAR user coordinates
GIND=.NOT.LGEOG .OR. &
!!!!!!!!!!!!!!! JOEL!!!!!!!!!!!!
     (.NOT.LGEOG .AND. &
      (LXYZ00 .OR. LMASK3D .OR. LMASK3D_XY .OR. LMASK3D_XZ .OR. LMASK3D_YZ &
!      .OR. LMARKER .OR. LTRAJ3D .OR. LFLUX3D)
       .OR. LMSKTOP .OR. LTRAJ3D .OR. LFLUX3D) .AND. LINDAX )
GCONF= .NOT.LGEOG .AND. &
       (LXYZ00 .OR. LMASK3D .OR. LMASK3D_XY .OR. LMASK3D_XZ .OR. LMASK3D_YZ &
!      .OR. LMARKER .OR. LTRAJ3D .OR. LFLUX3D)
       .OR. LMSKTOP .OR. LTRAJ3D .OR. LFLUX3D) .AND. .NOT.LINDAX 
IF (GCONF) GIND=.FALSE.
!!!!!!!!!!!!!!! JOEL!!!!!!!!!!!!

! limites du domaine en indices de grille
IF(GIND)THEN
   CALL SET(ZX1,ZX2,ZY1,ZY2,FLOAT(NIINF),FLOAT(NISUP),  &
        FLOAT(NJINF),FLOAT(NJSUP),IDUM)
!>>>>>>>>>>>>This section is to be revised***********************

  FORMAX='          '
  IF(LFMTAXEX)THEN
    FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
  ELSE
    FORMAX='(F5.1)'
  ENDIF
  FORMAY='          '
  IF(LFMTAXEY)THEN
    FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
  ELSE
    FORMAY='(F5.1)'
  ENDIF
  CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
! CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
! CALL LABMOD('(F5.1)','(F5.1)',0,0,10,10,0,0,0)
!CALL GASETI('LTY',1)

  IF(NCHPCITVXMJ /= 0 .OR. NCHPCITVYMJ /=0 .OR. NCHPCITVXMN /= 0 .OR. &
     NCHPCITVXMN /= 0)THEN
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,0,0,5,0.,0.)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
      CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,0,1,5,0.,0.)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,1,0,5,0.,0.)
    ELSE
      CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,1,1,5,0.,0.)
    ENDIF
!Avril 2002

  ELSE
    IF(NISUP > 99)THEN
      FORMAX='          '
      IF(LFMTAXEX)THEN
        FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
      ELSE
        FORMAX='(I4)'
      ENDIF
      FORMAY='          '
      IF(LFMTAXEY)THEN
        FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
      ELSE
        FORMAY='(I2)'
      ENDIF
      CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!     CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(I3)','(I2)',0,0,10,10,0,0,0)
!     CALL LABMOD('(I3)','(I2)',3,2,10,10,0,0,0)
      IF(NJSUP > 99)THEN
        FORMAY='          '
        IF(LFMTAXEY)THEN
          FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
        ELSE
          FORMAY='(I4)'
        ENDIF
        CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!       CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
!       CALL LABMOD('(I3)','(I3)',0,0,10,10,0,0,0)
!       CALL LABMOD('(I3)','(I3)',3,3,10,10,0,0,0)
      ENDIF
    ELSE  
      FORMAX='          '
      IF(LFMTAXEX)THEN
        FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
      ELSE
        FORMAX='(I2)'
      ENDIF
      FORMAY='          '
      IF(LFMTAXEY)THEN
        FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
      ELSE
        FORMAY='(I2)'
      ENDIF
      CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!     CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(I2)','(I2)',0,0,10,10,0,0,0)
!     CALL LABMOD('(I2)','(I2)',2,2,10,10,0,0,0)
      IF(NJSUP > 99)THEN
        FORMAY='          '
        IF(LFMTAXEY)THEN
          FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
        ELSE
          FORMAY='(I4)'
        ENDIF
        CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!       CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
!       CALL LABMOD('(I2)','(I3)',0,0,10,10,0,0,0)
!       CALL LABMOD('(I2)','(I3)',2,3,10,10,0,0,0)
      ENDIF
    ENDIF
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(1,NISUP-NIINF,1,NJSUP-NJINF,0,0,5,0.,0.)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
      CALL GRIDAL(1,NISUP-NIINF,1,NJSUP-NJINF,0,1,5,0.,0.)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
      CALL GRIDAL(1,NISUP-NIINF,1,NJSUP-NJINF,1,0,5,0.,0.)
    ELSE
      CALL GRIDAL(1,NISUP-NIINF,1,NJSUP-NJINF,1,1,5,0.,0.)
      !CALL PERIM(NISUP-NIINF,1,NJSUP-NJINF,1)
    ENDIF
!Avril 2002
  ENDIF
ENDIF
!
!!!!!!!!!!!!!!! JOEL!!!!!!!!!!!!
! limites du domaine en coord. conf. (pour lachers de part. LMASK3D)
IF(GCONF) THEN
  CALL SET(ZX1,ZX2,ZY1,ZY2,ZXMIN,ZXMAX,ZYMIN,ZYMAX,1)
  CALL LABMOD('(F8.0)','(F8.0)',0,0,NSZLBX,NSZLBY,12,0,0)
  CALL GRIDAL(1,NISUP-NIINF,1,NJSUP-NJINF,1,1,5,0.,0.)
ENDIF
!!!!!!!!!!!!!!! JOEL!!!!!!!!!!!!
!
! limites du domaine en lat/lon
IF (LGEOG) THEN
  CALL SET(ZX1,ZX2,ZY1,ZY2,ZPL2,ZPL4,ZPL1,ZPL3,IDUM)
  FORMAY='          '
  IF(LFMTAXEY)THEN
    FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
  ELSE
    FORMAY='(F5.1)'
  ENDIF
  IF(ZPL2 < -99. .OR. ZPL4 < -99.)THEN
    FORMAX='          '
    IF(LFMTAXEX)THEN
      FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
    ELSE
      FORMAX='(F6.1)'
    ENDIF
! Ai mis 12 pour rapprocher les labels Y de l'axe; sinon troncature
    CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,12,0,0)
!   CALL LABMOD(FORMAX,FORMAY,0,0,10,10,12,0,0)
!   CALL LABMOD('(F6.1)','(F5.1)',0,0,10,10,12,0,0)
!   CALL LABMOD('(F6.1)','(F5.1)',6,5,10,10,0,0,0)
  ELSE
    FORMAX='          '
    IF(LFMTAXEX)THEN
      FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
    ELSE
      FORMAX='(F6.2)'
    ENDIF
    CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,12,0,0)
!   CALL LABMOD(FORMAX,FORMAY,0,0,10,10,12,0,0)
!   CALL LABMOD('(F6.2)','(F5.1)',0,0,10,10,12,0,0)
!   CALL LABMOD('(F6.2)','(F5.1)',6,5,10,10,0,0,0)
  ENDIF
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
  CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,0,0,5,0.,0.)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
  CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,0,1,5,0.,0.)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
  CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,1,0,5,0.,0.)
    ELSE
  CALL GRIDAL(NCHPCITVXMJ,NCHPCITVXMN,NCHPCITVYMJ,NCHPCITVXMN,1,1,5,0.,0.)
    ENDIF
!Avril 2002
ENDIF
! Initial coordinate transformation restored
CALL SET(ZX1,ZX2,ZY1,ZY2,ZXX1,ZXX2,ZYY1,ZYY2,IDUM)
!
!*     1.8    A series of landmarks is added to the plot when required
!
!!! Enleve le 30/8/99 pour travailler avec les coordonnees conformes ci-apres
!   verifie que idem
!IF(NLPCAR.GE.1)THEN
! DO JLPCAR=1,NLPCAR
!   CALL MAPTRN(XLATCAR(JLPCAR),XLONCAR(JLPCAR),ZU,ZV)
!>>>>>>>May be, this section is to be revised*******************
!   CALL NGWSYM('N',8,ZU,ZV,.012,1,0)
! Obsolete   CALL PWRITX(ZU,ZV,'''KGU''-',6,20,0,0)
! ENDDO
!ENDIF
! Initial coordinate transformation restored
CALL SET(ZX1,ZX2,ZY1,ZY2,ZXMIN,ZXMAX,ZYMIN,ZYMAX,IDUM)
if(nverbia > 0)then
  print *,' **bcgrd AP CALL SET'
endif

IF(K == 2)THEN
IF(NLPCAR.GE.1)THEN
  IF(.NOT.LCOLAREA .AND. .NOT.LCOLINE)THEN
    call tabcol_fordiachro
  ENDIF
  IF(LUMVM .OR. LUTVT .AND. NSUPERDIA == 1)THEN
    call tabcol_fordiachro
  ENDIF
  DO JLPCAR=1,NLPCAR
    ZLAT=XLATCAR(JLPCAR)
    ZLON=XLONCAR(JLPCAR)
    YSYMB=CSYMCAR(JLPCAR)
    ZPOS=XPOSNOM(JLPCAR)
    ICOLS=ICOLSYM(JLPCAR)
    ICOLN=ICOLNOM(JLPCAR)
    IF(XSZSYM(JLPCAR) /= 0.)THEN
      ZSZ=XSZSYM(JLPCAR)
      IF(ZSZ == 9999.)ZSZ=.012
    ELSE
      ZSZ=.012
    ENDIF
    CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZU,ZV)
!   CALL GSTXCI(ICOLS)
    CALL PCSETI('OC',ICOLS)
    IF(YSYMB == '.')THEN
      CALL NGWSYM('N',8,ZU,ZV,ZSZ,ICOLS,0)
!     CALL NGWSYM('N',8,ZU,ZV,ZSZ,1,0)
    ELSE
      CALL PCSETI('OF',2)
      CALL PCSETR('OL',1.5)
      CALL PLCHHQ(ZU,ZV,YSYMB,ZSZ,0.,0.)
      CALL PCSETI('OF',0)
      CALL PCSETR('OL',0.)
    ENDIF
    CALL PCSETI('OC',1)
    IF(XSZNOM(JLPCAR) /= 0.)THEN
      ZSZ=XSZNOM(JLPCAR)
      IF(ZSZ == 9999.)ZSZ=.012
    ELSE
      ZSZ=.012
    ENDIF
    IPOS=ZPOS
!   print *,' ZSZ NOM ',ZSZ
    SELECT CASE(IPOS)
      CASE(0)
	ZCENT=-1.
	ZU=ZU+ZSZ*1.1*(ZXMAX-ZXMIN)
      CASE(45)
	ZCENT=-1.
	ZU=ZU+ZSZ*1.0*(ZXMAX-ZXMIN)
	ZV=ZV+ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
      CASE(90)
	ZCENT=0.
	ZV=ZV+ZSZ*1.5*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!       ZV=ZV+ZSZ*1.5*(ZYMAX-ZYMIN)
      CASE(135)
	ZCENT=1.
	ZU=ZU-ZSZ*1.0*(ZXMAX-ZXMIN)
	ZV=ZV+ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!       ZV=ZV+ZSZ*1.0*(ZYMAX-ZYMIN)
      CASE(180)
	ZCENT=1.
	ZU=ZU-ZSZ*1.1*(ZXMAX-ZXMIN)
      CASE(225)
	ZCENT=1.
	ZU=ZU-ZSZ*1.0*(ZXMAX-ZXMIN)
	ZV=ZV-ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!       ZV=ZV-ZSZ*1.0*(ZYMAX-ZYMIN)
      CASE(270)
	ZCENT=0.
	ZV=ZV-ZSZ*1.5*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!       ZV=ZV-ZSZ*1.5*(ZYMAX-ZYMIN)
      CASE(315)
	ZCENT=-1.
	ZU=ZU+ZSZ*1.0*(ZXMAX-ZXMIN)
	ZV=ZV-ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!       ZV=ZV-ZSZ*1.0*(ZYMAX-ZYMIN)
    END SELECT 
    IF(CNOMCAR(JLPCAR) /= ' ')THEN
      YNOM=CNOMCAR(JLPCAR)
      YNOM=ADJUSTL(YNOM)
      CALL PCSETI('OF',2)
      CALL PCSETI('OC',ICOLN)
      !CALL PCSETR('OL',1.5)
      !MODIF SYLVIE D.: epaisseur des caracteres de CNOMSYM -> XLWNOM
      CALL PCSETR('OL',XLWCONT)
!     CALL GSTXCI(ICOLN)
!     CALL GSPLCI(ICOLN)
      CALL PLCHHQ(ZU,ZV,YNOM(1:LEN_TRIM(YNOM)),ZSZ,0.,ZCENT)
!     CALL PLCHHQ(ZU,ZV+ZSZ*1.5*(ZYMAX-ZYMIN),YNOM(1:LEN_TRIM(YNOM)),ZSZ,0.,ZCENT)
    ENDIF
    CALL PCSETI('OF',0)
    CALL PCSETR('OL',0.)
    CALL PCSETI('OC',1)
    CALL GSTXCI(1)
  ENDDO
ENDIF
IF(NIJCAR.GE.1)THEN
  IF(.NOT.LCOLAREA .AND. .NOT.LCOLINE)THEN
    call tabcol_fordiachro
  ENDIF
  DO JIJCAR=1,NIJCAR
    ZI=XICAR(JIJCAR)
    ZJ=XJCAR(JIJCAR)
    print *,' **bcgrd_fordiachro ZI,ZJ ',ZI,ZJ
    YSYMB=CSYMCAR(JIJCAR)
    ZPOS=XPOSNOM(JIJCAR)
    ICOLS=ICOLSYM(JIJCAR)
    ICOLN=ICOLNOM(JIJCAR)
    IF(XSZSYM(JIJCAR) /= 0.)THEN
      ZSZ=XSZSYM(JIJCAR)
      IF(ZSZ == 9999.)ZSZ=.012
    ELSE
      ZSZ=.012
    ENDIF
    ICONVI=INT(ZI)
    ICONVJ=INT(ZJ)
    if(nverbia > 0)then
    print *,' **bcgrd_fordiachro ICONVI, ICONVJ ',ICONVI,ICONVJ
    endif
    ZX=XXX(ICONVI,NMGRID)+(XXX(MIN(ICONVI+1,SIZE(XXX,1)),NMGRID)-XXX(ICONVI,NMGRID))*(ZI-FLOAT(ICONVI))
    ZY=XXY(ICONVJ,NMGRID)+(XXY(MIN(ICONVJ+1,SIZE(XXY,1)),NMGRID)-XXY(ICONVJ,NMGRID))*(ZJ-FLOAT(ICONVJ))
    if(nverbia > 0)then
    print *,' **bcgrd_fordiachro ZX,ZY ',ZX,ZY
    endif
    CALL PCSETI('OC',ICOLS)
    IF(YSYMB == '.')THEN
      CALL NGWSYM('N',8,ZX,ZY,ZSZ,ICOLS,0)
    ELSE
      CALL PCSETI('OF',2)
      CALL PCSETR('OL',1.5)
      CALL PLCHHQ(ZX,ZY,YSYMB,ZSZ,0.,0.)
      CALL PCSETI('OF',0)
      CALL PCSETR('OL',0.)
    ENDIF
    CALL PCSETI('OC',1)
    IF(XSZNOM(JIJCAR) /= 0.)THEN
      ZSZ=XSZNOM(JIJCAR)
      IF(ZSZ == 9999.)ZSZ=.012
    ELSE
      ZSZ=.012
    ENDIF
    IPOS=ZPOS
    SELECT CASE(IPOS)
      CASE(0)
	ZCENT=-1.
	ZX=ZX+ZSZ*1.1*(ZXMAX-ZXMIN)
      CASE(45)
	ZCENT=-1.
	ZX=ZX+ZSZ*1.0*(ZXMAX-ZXMIN)
	ZY=ZY+ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
      CASE(90)
	ZCENT=0.
	ZY=ZY+ZSZ*1.5*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
      CASE(135)
	ZCENT=1.
	ZX=ZX-ZSZ*1.0*(ZXMAX-ZXMIN)
	ZY=ZY+ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
      CASE(180)
	ZCENT=1.
	ZX=ZX-ZSZ*1.1*(ZXMAX-ZXMIN)
      CASE(225)
	ZCENT=1.
	ZX=ZX-ZSZ*1.0*(ZXMAX-ZXMIN)
	ZY=ZY-ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
      CASE(270)
	ZCENT=0.
	ZY=ZY-ZSZ*1.5*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
      CASE(315)
	ZCENT=-1.
	ZX=ZX+ZSZ*1.0*(ZXMAX-ZXMIN)
	ZY=ZY-ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
    END SELECT 
    IF(CNOMCAR(JIJCAR) /= ' ')THEN
      YNOM=CNOMCAR(JIJCAR)
      YNOM=ADJUSTL(YNOM)
      CALL PCSETI('OF',2)
      CALL PCSETI('OC',ICOLN)
      CALL PCSETR('OL',1.5)
      CALL PLCHHQ(ZX,ZY,YNOM(1:LEN_TRIM(YNOM)),ZSZ,0.,ZCENT)
    ENDIF
    CALL PCSETI('OF',0)
    CALL PCSETR('OL',0.)
    CALL PCSETI('OC',1)
    CALL GSTXCI(1)
  ENDDO
ENDIF
IF(LRADAR)THEN
  CALL GQLWSC(IERR,ZWIDTH)
  ZSZ=.012
  CALL GSLWSC(3.)
  IF(NPORTRAD1 /= 0)THEN
    ZLAT=XLATRAD1
    ZLON=XLONRAD1
    YSYMB=CSYMRAD1
    CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZU,ZV)
    CALL PLCHHQ(ZU,ZV,YSYMB,ZSZ,0.,0.)
    DO J=1,NPORTRAD1
      CALL TRACIRCLE(ZU,ZV,XPORTRAD1(J),XLWRAD1(J))
      CALL SFLUSH
    ENDDO
  ENDIF
  IF(NPORTRAD2 /= 0)THEN
    ZLAT=XLATRAD2
    ZLON=XLONRAD2
    YSYMB=CSYMRAD2
    CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZU,ZV)
    CALL PLCHHQ(ZU,ZV,YSYMB,ZSZ,0.,0.)
    DO J=1,NPORTRAD2
      CALL TRACIRCLE(ZU,ZV,XPORTRAD2(J),XLWRAD2(J))
      CALL SFLUSH
    ENDDO
  ENDIF
  IF(NPORTRAD3 /= 0)THEN
    ZLAT=XLATRAD3
    ZLON=XLONRAD3
    YSYMB=CSYMRAD3
    CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZU,ZV)
    CALL PLCHHQ(ZU,ZV,YSYMB,ZSZ,0.,0.)
    DO J=1,NPORTRAD3
      CALL TRACIRCLE(ZU,ZV,XPORTRAD3(J),XLWRAD3(J))
      CALL SFLUSH
    ENDDO
  ENDIF
  IF(NPORTRAD4 /= 0)THEN
    ZLAT=XLATRAD4
    ZLON=XLONRAD4
    YSYMB=CSYMRAD4
    CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZU,ZV)
    CALL PLCHHQ(ZU,ZV,YSYMB,ZSZ,0.,0.)
    DO J=1,NPORTRAD4
      CALL TRACIRCLE(ZU,ZV,XPORTRAD4(J),XLWRAD4(J))
      CALL SFLUSH
    ENDDO
  ENDIF
  CALL GSLWSC(ZWIDTH)
ENDIF

ENDIF
!
!----------------------------------------------------------------------
!
!*    2.     EXIT
!            ----
!
RETURN
END SUBROUTINE BCGRD_FORDIACHRO
