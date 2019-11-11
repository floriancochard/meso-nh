!     ######spl
      SUBROUTINE VERIFLEN_FORDIACHRO
!     ##############################
!
!!****  *VERIFLEN_FORDIACHRO* - Computes the length of the abscissa axis for the vertical
!!                   cross-sections and checks whether it gets out of the 
!!                   display boundaries
!!
!!    PURPOSE
!!    -------
!       Computes the meshsizes along the abscissa axis of vertical 
!     cross-sections and checks the requested number of points gets
!     out of the display boundaries. The calculation is made for all
!     the possible grids
!
!!**  METHOD
!!    ------
!!      -NA-
!!
!!    EXTERNAL
!!    --------
!!      LENMAILLD : locates the four corners of the x-y gridbox  containing
!!                  the starting point of a vertical cross section. This
!!                  information is a prerequisite to calculate the meshsizes
!!                  along a vertical cross-section abscissa axis.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_COORD   : declares gridpoint coordinates 
!!                           (TRACE use only)
!!         XXDXHAT, XXDYHAT : Mesh size arrays (meters), for all grid locations
!!         XXX, XXY   : XHAT, YHAT values (meters) for all grid locations
!!         XXDS       : Mesh size (meters) along the horizontal axis of an 
!!                      oblique vertical cross-section, for all grid locations
!!         XDS        : Abscissa array along the horizontal axis of an oblique
!!                      vertical cross-section (meters), for all grid locations
!!         XDSX, XDSY : Projections on the MESO-NH cartesian axes of the XDS 
!!                      oblique abscissa (meters), for all grid locations
!!
!!      Module MODD_GRID1      : declares grid variables (Model module)
!!         XXHAT, XYHAT : x, y in the conformal or cartesian plane
!!
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist (former PARA common)
!!         NIDEBCOU,NJDEBCOU :  Origin of a vertical cross-section
!!                              in grid index integer values
!!                              (XIDEBCOU and XJDEBCOU must 
!!                              be = to -999.)
!!         XIDEBCOU,XJDEBCOU :  Origin of a vertical cross-section
!!                              in cartesian (or conformal) real values
!!         NLANGLE           :  Angle between X Meso-NH axis and
!!                              cross-section direction in degrees
!!                              (Integer value anticlockwise)
!!         NLMAX,            :  Number of points horizontally along
!!                              the vertical section
!!         Module MODD_DIM1 : contains dimensions of data arrays
!!                  NIMAX,NKMAX    :  x, and z array dimensions
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!          JPHEXT : Horizontal external points number
!!          JPVEXT : Vertical external points number
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
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   14/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_COORD
USE MODD_CONF
USE MODD_DIM1
USE MODD_TYPE_AND_LH
USE MODD_NMGRID
USE MODD_GRID1
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODN_PARA
USE MODD_PARAMETERS
USE MODD_ALLOC_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_DEFCV
USE MODD_NMGRID
USE MODE_GRIDPROJ

IMPLICIT NONE

REAL :: ZRANGLE, ZCANGLE, ZZSANGLE, ZIX, ZIY, ZZSIC, ZZSIS
REAL,SAVE :: ZCVX1, ZCVX2, ZCVY1, ZCVY2
REAL,SAVE :: ZXREF1, ZXREF2, ZYREF1, ZYREF2
REAL :: ZME, ZMEY, ZMEX
REAL :: ZZLA, ZZLO

INTEGER :: J2LOOP, JILOOP, IIA, IJA, ISIZ
INTEGER :: IIB, IIE, IJB, IJE, IIU, IJU
INTEGER :: IID,  IJD, IIF, IJF
INTEGER :: ICV,  IINF, IJINF, ISUP, IJSUP
INTEGER :: IMODIF

LOGICAL :: GEND, GEND2
LOGICAL :: GOKX, GOKY
!
!-------------------------------------------------------------------------------
IIB=1+JPHEXT
IIE=NIMAX+JPHEXT
IIU=NIMAX+2*JPHEXT
IJB=1+JPHEXT
IJE=NJMAX+JPHEXT
IJU=NJMAX+2*JPHEXT
!
!*       1.   LOCATING THE STARTING GRIDBOX AND CHECKING FOR LOCATION
!*            OUT OF THE DISPLAY BOUNDARIES
!             -------------------------------------------------------
!
!*       1.0  array allocations
!
ISIZ=MAX(SIZE(XXHAT),SIZE(XYHAT))
IF(ALLOCATED(XDS))THEN
  DEALLOCATE(XDS)
END IF
  ALLOCATE(XDS(ISIZ+100,7))
IF(ALLOCATED(XXDS))THEN
  DEALLOCATE(XXDS)
END IF
  ALLOCATE(XXDS(ISIZ+100,7))
IF(ALLOCATED(XDSX))THEN
  DEALLOCATE(XDSX)
END IF
  ALLOCATE(XDSX(ISIZ+100,7))
IF(ALLOCATED(XDSY))THEN
  DEALLOCATE(XDSY)
END IF
  ALLOCATE(XDSY(ISIZ+100,7))
! Avril 2002
IF(LCV .AND. .NOT.LCARTESIAN)THEN
IF(ALLOCATED(XLATCV))THEN
  DEALLOCATE(XLATCV)
ENDIF
ALLOCATE(XLATCV(ISIZ+100))
IF(ALLOCATED(XLONCV))THEN
  DEALLOCATE(XLONCV)
ENDIF
ALLOCATE(XLONCV(ISIZ+100))
ENDIF
! Avril 2002
!
if(nverbia > 0)then
  print *,' ** veriflen LDEFCV2 LDEFCV2LL LDEFCV2IND ',LDEFCV2,LDEFCV2ll,LDEFCV2IND
endif

IF(LDEFCV2)THEN
ZCVX1=XIDEBCV; ZCVX2=XIFINCV; ZCVY1=XJDEBCV; ZCVY2=XJFINCV
LDEFCV2CC=.TRUE.
ELSE IF(LDEFCV2LL)THEN
CALL SM_XYHAT_S(XLATORI,XLONORI,XIDEBCVLL,XJDEBCVLL,ZCVX1,ZCVY1)
CALL SM_XYHAT_S(XLATORI,XLONORI,XIFINCVLL,XJFINCVLL,ZCVX2,ZCVY2)
LDEFCV2CC=.TRUE.
ELSE IF(LDEFCV2IND)THEN
ZCVX1=XXX(NIDEBCV,NMGRID)
ZCVY1=XXY(NJDEBCV,NMGRID)
ZCVX2=XXX(NIFINCV,NMGRID)
ZCVY2=XXY(NJFINCV,NMGRID)
LDEFCV2CC=.TRUE.
ELSE
LDEFCV2CC=.FALSE.
ENDIF
IF(LDEFCV2CC)THEN
  IINF=NIINF; ISUP=NISUP; IJINF=NJINF; IJSUP=NJSUP
  IF(LCV)THEN
    ICV=1
  ELSE
    ICV=0
    LCV=.TRUE.
  ENDIF
  CALL RESOLV_NIJINF_NIJSUP
ENDIF
!
!*       1.1   Checking successive gridbox locations along axis
!
IF(LDEFCV2CC)THEN
  ZRANGLE=ATAN2((ZCVY2-ZCVY1),(ZCVX2-ZCVX1))
! print *,' ** veriflen ZRANGLE,ZCANGLE,ZZSANGLE ',ZRANGLE
  IF(ZCVY2 == ZCVY1 .AND. ABS(ZRANGLE) < 1.E-6)THEN
    ZRANGLE=0.
  ENDIF
  XANGLECV=ZRANGLE
ELSE
  ZRANGLE=FLOAT(NLANGLE)*ACOS(-1.)/180.
ENDIF
ZCANGLE=COS(ZRANGLE)
ZZSANGLE=SIN(ZRANGLE)
if(nverbia > 0)then
  print *,' ** veriflen ZRANGLE,ZCANGLE,ZZSANGLE ',ZRANGLE,ZCANGLE,ZZSANGLE
endif
IF(.NOT.LDEFCV2CC)THEN
  IF(NLANGLE.EQ.0.OR.NLANGLE.EQ.180)ZZSANGLE=0.
  IF(NLANGLE.EQ.90.OR.NLANGLE.EQ.270)ZCANGLE=0.
ELSE
  IF(XANGLECV == 0. .OR. XANGLECV/ACOS(-1.)*180. == 180.)ZZSANGLE=0.
  IF(XANGLECV/ACOS(-1.)*180. == 90. .OR.XANGLECV/ACOS(-1.)*180. == 270.)ZCANGLE=0.
ENDIF
ZZSIC=SIGN(1.,ZCANGLE)
ZZSIS=SIGN(1.,ZZSANGLE)
IF(LDEFCV2CC)THEN
  XIDEBCOU=ZCVX1; XJDEBCOU=ZCVY1
  NLMAX=500
  ZXREF1=MIN(ZCVX1,ZCVX2)
  ZXREF2=MAX(ZCVX1,ZCVX2)
  ZYREF1=MIN(ZCVY1,ZCVY2)
  ZYREF2=MAX(ZCVY1,ZCVY2)
  if(nverbia > 0)then
    print *,' *** veriflen XIDEBCOU XJDEBCOU NLMAX AV calcul ZXREF1,ZXREF2,ZYREF1,ZYREF2'
    print *,' *** veriflen',XIDEBCOU,XJDEBCOU,NLMAX,ZXREF1,ZXREF2,ZYREF1,ZYREF2
  print *,' ** veriflen ZRANGLE,ZCANGLE,ZZSANGLE ',ZRANGLE,ZCANGLE,ZZSANGLE
  endif
ENDIF
!
! Verification origine OK
!
IF(XIDEBCOU.EQ.-999.)THEN
  IF(NIDEBCOU >= NIL .AND. NIDEBCOU <= NIH)THEN
    GOKX=.TRUE.
  ELSE
    print *,' NIDEBCOU EN DEHORS DES LIMITES en X ',NIDEBCOU,' (',NIL,' - ', &
    NIH,')'
    GOKX=.FALSE.
  ENDIF
  IF(NJDEBCOU >= NJL .AND. NJDEBCOU <= NJH)THEN
    GOKY=.TRUE.
  ELSE
    print *,' NJDEBCOU EN DEHORS DES LIMITES en Y ',NJDEBCOU,' (',NJL,' - ', &
    NJH,')'
    GOKY=.FALSE.
  ENDIF
ELSE
  IF(XIDEBCOU >= XXX(NIL,NMGRID) .AND. XIDEBCOU <= XXX(NIH,NMGRID))THEN
    GOKX=.TRUE.
  ELSE
    print *,' XIDEBCOU EN DEHORS DES LIMITES en X ',XIDEBCOU,' (',  &
    XXX(NIL,NMGRID),' - ', &
    XXX(NIH,NMGRID),')'
    GOKX=.FALSE.
  ENDIF
  IF(XJDEBCOU >= XXY(NJL,NMGRID) .AND. XJDEBCOU <= XXY(NJH,NMGRID))THEN
    GOKY=.TRUE.
  ELSE
    print *,' XJDEBCOU EN DEHORS DES LIMITES en Y ',XJDEBCOU,' (',  &
    XXY(NJL,NMGRID),' - ', &
    XXY(NJH,NMGRID),')'
    GOKY=.FALSE.
  ENDIF
ENDIF
IF(.NOT.GOKX .OR. .NOT.GOKY)THEN
  print *,' -> ABORT: REDEFINISSEZ L'' ORIGINE DE LA COUPE '
  LPBREAD=.TRUE.
  !RETURN
  STOP    
ENDIF
!
! Scanning all the existing grids                           
! J2LOOP --> NGRID
!
IMODIF=0
DO J2LOOP=1,7                                           !do 1 (grid loop)
GEND=.FALSE.
GEND2=.FALSE.
!print *,' GRILLE NLMAX ',J2LOOP,' ',NLMAX
  IF(XIDEBCOU.EQ.-999.)THEN    ! Section defined by indexes             
    ZIX=XXDXHAT(NIDEBCOU,J2LOOP)
      IF(ZZSIC.LT.0.)ZIX=XXDXHAT(MAX(NIL,NIDEBCOU-1),J2LOOP)
!     IF(ZZSIC.LT.0.)ZIX=XXDXHAT(MAX(1,NIDEBCOU-1),J2LOOP)
    ZIY=XXDYHAT(NJDEBCOU,J2LOOP)
      IF(ZZSIS.LT.0.)ZIY=XXDYHAT(MAX(NJL,NJDEBCOU-1),J2LOOP)
!     IF(ZZSIS.LT.0.)ZIY=XXDYHAT(MAX(1,NJDEBCOU-1),J2LOOP)
    XDSX(1,J2LOOP)=XXX(NIDEBCOU,J2LOOP)
    XDSY(1,J2LOOP)=XXY(NJDEBCOU,J2LOOP)
  ELSE                         ! Section defined by range
    XDSX(1,J2LOOP)=XIDEBCOU
    XDSY(1,J2LOOP)=XJDEBCOU
    CALL LENMAILLD(XIDEBCOU,XJDEBCOU,IIA,IJA,1,J2LOOP)
    ZIX=XXDXHAT(IIA,J2LOOP)
    ZIY=XXDYHAT(IJA,J2LOOP)
    if(nverbia > 0)then
      print *,' veriflen XIDEBCOU,XJDEBCOU,ZIX,ZIY ',XIDEBCOU,XJDEBCOU,ZIX,ZIY
    endif
  END IF
!
! Scans oblique abscissa from origin to end.
! XDS  ---> X along oblique cross-section
! XXDS ---> X-meshsize along X of oblique cross-section
!
  XDS(1,J2LOOP)=0.
! print *,' TINY ',TINY(1.)
    DO JILOOP=2,NLMAX                                   ! do 2 (abscissa loop)
      XXDS(JILOOP-1,J2LOOP)=ABS(ZIX*ZCANGLE)+ABS(ZIY*ZZSANGLE)
      if(nverbia >8)then
	print *,' **** veriflen boucle DO JILOOP=2,NLMAX, XXDS(JILOOP-1,J2LOOP)',XXDS(JILOOP-1,J2LOOP),JILOOP-1,XXDS(1,J2LOOP)
      endif
      XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+XXDS(JILOOP-1,J2LOOP)
      XDSX(JILOOP,J2LOOP)=XDSX(JILOOP-1,J2LOOP)+XXDS(JILOOP-1,J2LOOP)*ZCANGLE
      XDSY(JILOOP,J2LOOP)=XDSY(JILOOP-1,J2LOOP)+XXDS(JILOOP-1,J2LOOP)*ZZSANGLE
!
! Checks whether the section length fits into the displayed domain?
!
        IF(LDEFCV2CC)THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        IF(ABS((ZCVX1-ZCVX2)/zcvx1) > 1.E-7)THEN
!!*******************************************************
         IF(ZCVX1 /= ZCVX2)THEN
!!*******************************************************
           
          if(nverbia > 0)then
	    print '('' +++veriflen ZCVX1,ZCVX2,ZXREF1,ZXREF2,XDSX(JILOOP,J2LOOP)'',5(2X,F12.3))',ZCVX1,ZCVX2, &
	    ZXREF1,ZXREF2,XDSX(JILOOP,J2LOOP)
	    print '('' +++veriflen ZCVY1,ZCVY2,ZYREF1,ZYREF2,XDSY(JILOOP,J2LOOP)'',5(2X,F12.3))',ZCVY1,ZCVY2, &
	    ZYREF1,ZYREF2,XDSY(JILOOP,J2LOOP)
	  endif

          IF(XDSX(JILOOP,J2LOOP) < ZXREF1) THEN
!         IF(XDSX(JILOOP,J2LOOP) <= ZXREF1) THEN
	    XDSX(JILOOP,J2LOOP) = ZXREF1
	    IF(ZXREF1 == ZCVX1)THEN
	      XDSY(JILOOP,J2LOOP)= ZCVY1
	    ELSE
	      XDSY(JILOOP,J2LOOP)= ZCVY2
	    ENDIF
	    ZMEY=ABS(XDSY(JILOOP,J2LOOP)-XDSY(JILOOP-1,J2LOOP))
	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
	    ZME=ABS(ZMEX*ZCANGLE) + ABS(ZMEY*ZZSANGLE)
            if(NVERBIA > 0)THEN
	      print *,' AP IF(XDSX(JILOOP,J2LOOP) < ZXREF1 Longueur de la derniere maille calculee ',ZME
            endif
	    XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+ZME
	    XXDS(JILOOP-1,J2LOOP)=ZME
	    NLMAX=JILOOP
            if(NVERBIA > 0)THEN
	      print *,' Controles . NLMAX calcule : ',NLMAX,' Grille N.',J2LOOP
	      print *,' Controles . Coord. conformes des extremites de la coupe demandees :'
	      print *,' (',ZCVX1,',',ZCVY1,')   (',ZCVX2,',',ZCVY2,')'
	      print *,' Controles . Coord. conformes des extremites de la coupe calculees :'
	      print *,' (',XDSX(1,J2LOOP),',',XDSY(1,J2LOOP),')   (',XDSX(NLMAX,J2LOOP),',',XDSY(NLMAX,J2LOOP),')'
	      print *,' xds xdsx xdsy ZCANGLE ZSANGLE ', ZCANGLE,ZZSANGLE 
	      print *,' **** XDS'
	      print *,xds(1:nlmax,j2loop)
	      print *,' **** XXDS'
	      print *,xxds(1:nlmax,j2loop)
	      print *,' **** XDSX'
	      print *,xdsx(1:nlmax,j2loop)
	      print *,' **** XDSY'
	      print *,xdsy(1:nlmax,j2loop)
            endif
	    EXIT
          ELSE IF(XDSX(JILOOP,J2LOOP) > ZXREF2)THEN
!         ELSE IF(XDSX(JILOOP,J2LOOP) >= ZXREF2)THEN
	    XDSX(JILOOP,J2LOOP) = ZXREF2
	    IF(ZXREF2 == ZCVX1)THEN
	      XDSY(JILOOP,J2LOOP)= ZCVY1
	    ELSE
	      XDSY(JILOOP,J2LOOP)= ZCVY2
	    ENDIF
	    ZMEY=ABS(XDSY(JILOOP,J2LOOP)-XDSY(JILOOP-1,J2LOOP))
!	    IF(ABS(ZZSANGLE) < 1.E-32)THEN
!	    ZME=ZMEY
!	    ELSE
!	    ZME=ABS(ZMEY/ZZSANGLE)
!	    ENDIF
	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
	    ZME=ABS(ZMEX*ZCANGLE) + ABS(ZMEY*ZZSANGLE)
	    IF(NVERBIA > 0)THEN
	      print *,' AP IF(XDSX(JILOOP,J2LOOP) > ZXREF2  Longueur de la derniere maille calculee ',ZME
	    ENDIF
	    XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+ZME
	    XXDS(JILOOP-1,J2LOOP)=ZME
	    NLMAX=JILOOP
            if(NVERBIA > 0)THEN
	      print *,' Controles . NLMAX calcule : ',NLMAX,' Grille N.',J2LOOP
	      print *,' Controles . Coord. conformes des extremites de la coupe demandees :'
	      print *,' (',ZCVX1,',',ZCVY1,')   (',ZCVX2,',',ZCVY2,')'
	      print *,' Controles . Coord. conformes des extremites de la coupe calculees :'
	      print *,' (',XDSX(1,J2LOOP),',',XDSY(1,J2LOOP),')   (',XDSX(NLMAX,J2LOOP),',',XDSY(NLMAX,J2LOOP),')'
	      print *,' xds xdsx xdsy ZCANGLE ZSANGLE ', ZCANGLE,ZZSANGLE 
	      print *,' **** XDS'
	      print *,xds(1:nlmax,j2loop)
	      print *,' **** XXDS'
	      print *,xxds(1:nlmax,j2loop)
	      print *,' **** XDSX'
	      print *,xdsx(1:nlmax,j2loop)
	      print *,' **** XDSY'
	      print *,xdsy(1:nlmax,j2loop)
            endif
	    EXIT
!!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          ELSE IF(XDSY(JILOOP,J2LOOP) < ZYREF1) THEN
!         ELSEIF(XDSY(JILOOP,J2LOOP) <= ZYREF1) THEN
	    XDSY(JILOOP,J2LOOP) = ZYREF1
	    IF(ZYREF1 == ZYREF2)THEN
	      IF(ABS(XDSX(JILOOP-1,J2LOOP)-ZCVX1) < &
	      ABS(XDSX(JILOOP-1,J2LOOP)-ZCVX2))THEN
		XDSX(JILOOP,J2LOOP)= ZCVX1
	      ELSE
		XDSX(JILOOP,J2LOOP)= ZCVX2
	      ENDIF
	    ELSE
	    IF(ZYREF1 == ZCVY1)THEN
	      XDSX(JILOOP,J2LOOP)= ZCVX1
	    ELSE
	      XDSX(JILOOP,J2LOOP)= ZCVX2
	    ENDIF
	    ENDIF
	    ZMEY=ABS(XDSY(JILOOP,J2LOOP)-XDSY(JILOOP-1,J2LOOP))
	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
	    ZME=ABS(ZMEX*ZCANGLE) + ABS(ZMEY*ZZSANGLE)
	    IF(NVERBIA > 0)THEN
	      print *,' AP IF(XDSY(JILOOP,J2LOOP) <= ZYREF1 Longueur de la derniere maille calculee ',ZME
	    ENDIF
	    XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+ZME
	    XXDS(JILOOP-1,J2LOOP)=ZME
	    NLMAX=JILOOP
            if(NVERBIA > 0)THEN
	      print *,' Controles . NLMAX calcule : ',NLMAX,' Grille N.',J2LOOP
	      print *,' Controles . Coord. conformes des extremites de la coupe demandees :'
	      print *,' (',ZCVX1,',',ZCVY1,')   (',ZCVX2,',',ZCVY2,')'
	      print *,' Controles . Coord. conformes des extremites de la coupe calculees :'
	      print *,' (',XDSX(1,J2LOOP),',',XDSY(1,J2LOOP),')   (',XDSX(NLMAX,J2LOOP),',',XDSY(NLMAX,J2LOOP),')'
	      print *,' xds xdsx xdsy ZCANGLE ZSANGLE ', ZCANGLE,ZZSANGLE 
	      print *,' **** XDS'
	      print *,xds(1:nlmax,j2loop)
	      print *,' **** XXDS'
	      print *,xxds(1:nlmax,j2loop)
	      print *,' **** XDSX'
	      print *,xdsx(1:nlmax,j2loop)
	      print *,' **** XDSY'
	      print *,xdsy(1:nlmax,j2loop)
            endif
	    EXIT
	  ELSE IF(XDSY(JILOOP,J2LOOP) > ZYREF2)THEN
!         ELSE IF(XDSY(JILOOP,J2LOOP) >= ZYREF2)THEN
	    IF(ZYREF1 == ZYREF2)THEN
	      IF(ABS(XDSX(JILOOP-1,J2LOOP)-ZCVX1) < &
	      ABS(XDSX(JILOOP-1,J2LOOP)-ZCVX2))THEN
		XDSX(JILOOP,J2LOOP)= ZCVX1
	      ELSE
		XDSX(JILOOP,J2LOOP)= ZCVX2
	      ENDIF
	    ELSE
	    XDSY(JILOOP,J2LOOP) = ZYREF2
	    IF(ZYREF2 == ZCVY1)THEN
	      XDSX(JILOOP,J2LOOP)= ZCVX1
	    ELSE
	      XDSX(JILOOP,J2LOOP)= ZCVX2
	    ENDIF
	    ENDIF
	    ZMEY=ABS(XDSY(JILOOP,J2LOOP)-XDSY(JILOOP-1,J2LOOP))
	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
	    ZME=ABS(ZMEX*ZCANGLE) + ABS(ZMEY*ZZSANGLE)
            if(NVERBIA > 0)THEN
	      print *,' AP ELSE IF(XDSY(JILOOP,J2LOOP) >= ZYREF2  Longueur de la derniere maille calculee ',ZME
            endif
	    XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+ZME
	    XXDS(JILOOP-1,J2LOOP)=ZME
	    NLMAX=JILOOP
            if(NVERBIA > 0)THEN
	      print *,' Controles . NLMAX calcule : ',NLMAX,' Grille N.',J2LOOP
	      print *,' Controles . Coord. conformes des extremites de la coupe demandees :'
	      print *,' (',ZCVX1,',',ZCVY1,')   (',ZCVX2,',',ZCVY2,')'
	      print *,' Controles . Coord. conformes des extremites de la coupe calculees :'
	      print *,' (',XDSX(1,J2LOOP),',',XDSY(1,J2LOOP),')   (',XDSX(NLMAX,J2LOOP),',',XDSY(NLMAX,J2LOOP),')'
	      print *,' xds xdsx xdsy ZCANGLE ZSANGLE ', ZCANGLE,ZZSANGLE 
	      print *,' **** XDS'
	      print *,xds(1:nlmax,j2loop)
	      print *,' **** XXDS'
	      print *,xxds(1:nlmax,j2loop)
	      print *,' **** XDSX'
	      print *,xdsx(1:nlmax,j2loop)
	      print *,' **** XDSY'
	      print *,xdsy(1:nlmax,j2loop)
            endif
	    EXIT
!!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          
	  ENDIF

!!*******************************************************
	  ELSE
!!*******************************************************

          IF(XDSY(JILOOP,J2LOOP) < ZYREF1) THEN
!         IF(XDSY(JILOOP,J2LOOP) <= ZYREF1) THEN
	    XDSY(JILOOP,J2LOOP) = ZYREF1
	    IF(ZYREF1 == ZCVY1)THEN
	      XDSX(JILOOP,J2LOOP)= ZCVX1
	    ELSE
	      XDSX(JILOOP,J2LOOP)= ZCVX2
	    ENDIF
	    ZMEY=ABS(XDSY(JILOOP,J2LOOP)-XDSY(JILOOP-1,J2LOOP))
	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
	    ZME=ABS(ZMEX*ZCANGLE) + ABS(ZMEY*ZZSANGLE)
	    IF(NVERBIA > 0)THEN
	      print *,' AP IF(XDSY(JILOOP,J2LOOP) < ZYREF1 Longueur de la derniere maille calculee ',ZME
	    ENDIF
	    XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+ZME
	    XXDS(JILOOP-1,J2LOOP)=ZME
	    NLMAX=JILOOP
	    IF(NVERBIA > 0)THEN
	      print *,' Controles . NLMAX calcule : ',NLMAX,' Grille N.',J2LOOP
	      print *,' Controles . Coord. conformes des extremites de la coupe demandees :'
	      print *,' (',ZCVX1,',',ZCVY1,')   (',ZCVX2,',',ZCVY2,')'
	      print *,' Controles . Coord. conformes des extremites de la coupe calculees :'
	      print *,' (',XDSX(1,J2LOOP),',',XDSY(1,J2LOOP),')   (',XDSX(NLMAX,J2LOOP),',',XDSY(NLMAX,J2LOOP),')'
	      print *,' xds xdsx xdsy ZCANGLE ZSANGLE ', ZCANGLE,ZZSANGLE 
	      print *,' **** XDS'
	      print *,xds(1:nlmax,j2loop)
	      print *,' **** XXDS'
	      print *,xxds(1:nlmax,j2loop)
	      print *,' **** XDSX'
	      print *,xdsx(1:nlmax,j2loop)
	      print *,' **** XDSY'
	      print *,xdsy(1:nlmax,j2loop)
            endif
	    EXIT
	  ELSE IF(XDSY(JILOOP,J2LOOP) > ZYREF2)THEN
!         ELSE IF(XDSY(JILOOP,J2LOOP) >= ZYREF2)THEN
	    XDSY(JILOOP,J2LOOP) = ZYREF2
	    IF(ZYREF2 == ZCVY1)THEN
	      XDSX(JILOOP,J2LOOP)= ZCVX1
	    ELSE
	      XDSX(JILOOP,J2LOOP)= ZCVX2
	    ENDIF
	    ZMEY=ABS(XDSY(JILOOP,J2LOOP)-XDSY(JILOOP-1,J2LOOP))
	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
	    ZME=ABS(ZMEX*ZCANGLE) + ABS(ZMEY*ZZSANGLE)
            if(NVERBIA > 0)THEN
	      print *,' AP ELSE IF(XDSY(JILOOP,J2LOOP) > ZYREF2  Longueur de la derniere maille calculee ',ZME
	    ENDIF
	    XDS(JILOOP,J2LOOP)=XDS(JILOOP-1,J2LOOP)+ZME
	    XXDS(JILOOP-1,J2LOOP)=ZME
!	    ZMEX=ABS(XDSX(JILOOP,J2LOOP)-XDSX(JILOOP-1,J2LOOP))
!	    IF(ABS(ZCANGLE) < 1.E-32)THEN
!	    ZME=ZMEX
!	    ELSE
!	    ZME=ABS(ZMEX/ZCANGLE)
!	    ENDIF
!	    IF(NVERBIA > 0)THEN
!	      print *,' Longueur de la derniere maille calculee avec COS pour controle ',ZME
!	    ENDIF
	    NLMAX=JILOOP
            if(NVERBIA > 0)THEN
	      print *,' Controles . NLMAX calcule : ',NLMAX,' Grille N.',J2LOOP
	      print *,' Controles . Coord. conformes des extremites de la coupe demandees :'
	      print *,' (',ZCVX1,',',ZCVY1,')   (',ZCVX2,',',ZCVY2,')'
	      print *,' Controles . Coord. conformes des extremites de la coupe calculees :'
	      print *,' (',XDSX(1,J2LOOP),',',XDSY(1,J2LOOP),')   (',XDSX(NLMAX,J2LOOP),',',XDSY(NLMAX,J2LOOP),')'
	      print *,' xds xdsx xdsy ZCANGLE ZSANGLE ', ZCANGLE,ZZSANGLE 
	      print *,' **** XDS'
	      print *,xds(1:nlmax,j2loop)
	      print *,' **** XXDS'
	      print *,xxds(1:nlmax,j2loop)
	      print *,' **** XDSX'
	      print *,xdsx(1:nlmax,j2loop)
	      print *,' **** XDSY'
	      print *,xdsy(1:nlmax,j2loop)
            endif
	    EXIT
          
	  ENDIF

!!*******************************************************
	  ENDIF
!!*******************************************************

	ELSE              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IF(LPOINTG)THEN
	  IID=1
	  IIF=IIU
	ELSE
	  IID=IIB
	  IIF=IIE
	ENDIF
        IF(XDSX(JILOOP,J2LOOP).LT.XXX(MAX(NIL,IID),J2LOOP).OR.  &
           XDSX(JILOOP,J2LOOP).GT.  &
        XXX(MIN(NIH,IIF),J2LOOP))THEN
          print *,' Vertical section overflows postprocessing window: ', &
          'X boundary reached after ',JILOOP,' points.'
          print *,' Requested number of points: ',NLMAX
          print *,' Computed X : ',XDSX(JILOOP,J2LOOP),' XMIN(NIL): ',  &
          XXX(NIL,J2LOOP),' XMAX(NIH): ',XXX(NIH,J2LOOP),' XMIN(1 or IIB): ', &
          XXX(IID,J2LOOP),' XMAX(IIE or IIU): ',XXX(IIF,J2LOOP)
!         STOP
          GEND=.TRUE.
!         print *,' NLMAX AVANT MODIF GRILLE ',NLMAX,' ',J2LOOP
	  IF(LPOINTG)THEN
            NLMAX=JILOOP
          ELSE
            NLMAX=JILOOP-1
	  ENDIF
          print *,' NLMAX APRES MODIF, NIDEBCOU NIL NIH ',NLMAX,NIDEBCOU, &
          NIL,NIH
          EXIT
        END IF
	IF(LPOINTG)THEN
	  IJD=1
	  IJF=IJU
	ELSE
	  IJD=IJB
	  IJF=IJE
	ENDIF
        IF(XDSY(JILOOP,J2LOOP).LT.XXY(MAX(NJL,IJD),J2LOOP).OR.  &
           XDSY(JILOOP,J2LOOP).GT.  &
        XXY(MIN(NJH,IJF),J2LOOP))THEN
          print *,' Vertical section overflows postprocessing window: ', &
          'Y  boundary reached after ',JILOOP,' points.'
          print *,' Requested number of points : ',NLMAX
          print *,' Computed Y : ',XDSY(JILOOP,J2LOOP),' YMIN(NJL): ',  &
          XXY(NJL,J2LOOP),' YMAX(NJH): ',XXY(NJH,J2LOOP),' YMIN(1 or IJB): ', &
          XXY(IJD,J2LOOP),' YMAX(IJE or IJU): ',XXY(IJF,J2LOOP)
!         STOP
          GEND=.TRUE.
!         print *,' NLMAX AVANT MODIF   GRILLE ',NLMAX,' ',J2LOOP
	  IF(LPOINTG)THEN
	    NLMAX=JILOOP
	  ELSE
            NLMAX=JILOOP-1
	  ENDIF
          print *,' NLMAX APRES MODIF, NJDEBCOU NJL NJH ',NLMAX,NJDEBCOU, &
          NJL,NJH
          EXIT
        END IF

	ENDIF       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gets IIA, IJA indexes to move forward to  next meshbox
!
        CALL LENMAILLD(XDSX(JILOOP,J2LOOP),XDSY(JILOOP,J2LOOP),IIA,IJA,  &
        JILOOP,J2LOOP)

        IF(GEND)THEN
          print *,' NLMAX AVANT MODIF ',NLMAX,' pour grille ',J2LOOP
          NLMAX=JILOOP
          print *,' NLMAX APRES MODIF ',NLMAX
          IF(XIDEBCOU.EQ.-999.)THEN
            print *,' NIDEBCOU,NJDEBCOU,NIL,NIH,NJL,NJH ',NIDEBCOU, &
            NJDEBCOU,NIL,NIH,NJL,NJH
          ENDIF
          EXIT
        ENDIF
        IF(GEND2)THEN
          print *,' NLMAX AVANT MODIF ',NLMAX,' pour grille ',J2LOOP
          NLMAX=JILOOP-1
          print *,' NLMAX APRES MODIF ',NLMAX
          IMODIF=J2LOOP ! car GEND2 remis a f pour grille suivante
          EXIT
        ENDIF
!
        ZIX=XXDXHAT(IIA,J2LOOP)
        ZIY=XXDYHAT(IJA,J2LOOP)

    ENDDO                                                        ! enddo 2
    !
ENDDO                                                            ! enddo 1
! Avril 2002 Calcul lat,lon de la coupe
IF(LCV .AND. .NOT.LCARTESIAN)THEN
  DO J2LOOP=1,1                                           !do 1 (grid loop)
    DO JILOOP=1,NLMAX                                       !do 2
      CALL SM_LATLON_S(XLATORI,XLONORI,XDSX(JILOOP,J2LOOP),&
      XDSY(JILOOP,J2LOOP),ZZLA,ZZLO)
      XLATCV(JILOOP)=ZZLA
      XLONCV(JILOOP)=ZZLO
    ENDDO                                                            ! enddo 2
if(nverbia > 0)then
  print *,' *** LATCV ',XLATCV(1:NLMAX)
  print *,' *** LONCV ',XLONCV(1:NLMAX)
endif
  ENDDO                                                            ! enddo 1
  IF (IMODIF/=0 .AND. LDEFCV2LL) THEN
    ! prise en compte du chgt d extremite 
    XIFINCVLL=XLATCV(NLMAX) ; XJFINCVLL=XLONCV(NLMAX)
  END IF
ENDIF
! Avril 2002
IF(LDEFCV2CC)THEN
  NIINF=IINF; NISUP=ISUP; NJINF=IJINF; NJSUP=IJSUP
  IF(ICV == 0)THEN
    LCV=.FALSE.
  ENDIF
ENDIF
!
CONTAINS
!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
!*    2.     CONTAINED ROUTINE LENMAILLD
!            ---------------------------
!-------------------------------------------------------------------------- 
!     ###################################################
      SUBROUTINE LENMAILLD(PSX,PSY,KIA,KJA,KILOOP,K2LOOP)
!     ###################################################
!
!!****  *LENMAILLD* - Gets the I,J indexes of the gribox containing the current
!!****                point along the abscissa of a vertical cross-section.
!!
!!    PURPOSE
!!    -------
!       Computes the (KIA,KJA) indexes of the gridbox containing the current
!     (PSX,PSY) point along the abscissa of a vertical cross-section, and
!     checks whether the point is within the limits of the postprocessing
!     window. Test is made using grid number K2LOOP to locate the gridpoints.
!
!!**  METHOD
!!    ------
!!     -NA-
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    
!!      None
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
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   14/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!
REAL :: PSX, PSY               ! Given gridpoint location (meters)
INTEGER :: KIA, KJA            ! Return indexes to be used for the relevant
                               ! gridbox containing the given point
INTEGER :: K2LOOP              ! Selector of the grid to be used
INTEGER :: KILOOP              ! Current value of the "oblique index" along
                               ! the oblique vertical cross-section
!
!*       0.2   Local variables
!
INTEGER :: JI, JJ
!
!-------------------------------------------------------------------------------
!
!*        1.    LOCATES CIRCUMSCRIBING GRIDBOX AND CHECKS FOR OVERFLOW
!               ------------------------------------------------------
! X scanning
!
DO JI=NIL,NIH
  IF(PSX.LE.XXX(JI,K2LOOP))GO TO 1
ENDDO
!
print *,' Out of TRACE window X=',PSX,' XMAX=',XXX(NIH,K2LOOP)
!! no more STOP
!STOP
KIA=NIH
GEND2=.TRUE.
!! no more STOP
!
1 CONTINUE
!
!! no more STOP
IF (.NOT. GEND2) THEN
!! no more STOP
IF(ABS(PSX-XXX(JI,K2LOOP)).LE.ABS(PSX-XXX(MAX(NIL,JI-1),K2LOOP)))THEN
  IF(ZZSIC.GT.0.)KIA=JI
  IF(ZZSIC.LT.0.)KIA=MAX(NIL,JI-1)
ELSE
  IF(ZZSIC.GT.0.)KIA=MAX(NIL,JI-1)
  IF(ZZSIC.LT.0.)KIA=MAX(NIL,JI-2)
END IF
!! no more STOP
END IF
!! no more STOP
!
! Y scanning
!
DO JJ=NJL,NJH
  IF(PSY.LE.XXY(JJ,K2LOOP))GO TO 2
ENDDO
!
print *,' Out of TRACE window Y=',PSY,' YMAX=',XXY(NJH,K2LOOP)
!! no more STOP
!STOP
KJA=NJH
GEND2=.TRUE.
!! no more STOP
!
2 CONTINUE
!
!! no more STOP
IF (.NOT. GEND2) THEN
!! no more STOP
IF(ABS(PSY-XXY(JJ,K2LOOP)).LE.ABS(PSY-XXY(MAX(NJL,JJ-1),K2LOOP)))THEN
  IF(ZZSIC.GT.0.)KJA=JJ
  IF(ZZSIC.LT.0.)KJA=MAX(NJL,JJ-1)
ELSE
  IF(ZZSIC.GT.0.)KJA=MAX(NJL,JJ-1)
  IF(ZZSIC.LT.0.)KJA=MAX(NJL,JJ-2)
END IF
!! no more STOP
END IF
!! no more STOP
!
! Index range control
!
IF(KIA.GE.NIH.AND.KILOOP.NE.NLMAX.AND.ZCANGLE.NE.0.)THEN
  print *,' Out of TRACE window, X limit reached',  &
  ' after ',KILOOP,' points.'
  print *,' Requested number of points : ',NLMAX
  print *,' Computed X : ',XDSX(KILOOP,K2LOOP),' XMIN : ',  &
  XXX(NIL,K2LOOP),' XMAX : ',XXX(NIH,K2LOOP)
  GEND=.TRUE.
! EXIT
! STOP
END IF
IF(KJA.GE.NJH.AND.KILOOP.NE.NLMAX.AND.ZZSANGLE.NE.0.)THEN
  print *,' Out of TRACE window, Y limit reached',  &
  '  after',KILOOP,' points.'
  print *,' Requested number of points : ',NLMAX
  print *,' Computed Y : ',XDSY(KILOOP,K2LOOP),' YMIN : ',  &
  XXY(NJL,K2LOOP),' YMAX : ',XXY(NJH,K2LOOP)
  GEND=.TRUE.
! EXIT
! STOP
END IF
!
!------------------------------------------------------------------------------
!
!*     2.    EXIT
!            ----
! 
END SUBROUTINE  LENMAILLD
!------------------------------------------------------------------------------
END SUBROUTINE  VERIFLEN_FORDIACHRO
