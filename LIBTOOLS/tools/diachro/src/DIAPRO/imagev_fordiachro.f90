!     ######spl
      SUBROUTINE IMAGEV_FORDIACHRO(PU,PV,KLREF,HTEXTE)
!     ################################################
!
!!****  *IMAGEV_FORDIACHRO* - Draws a vector arrow plot for an horizontal cross-section
!!
!!    PURPOSE
!!    -------
!       Draws an arrow plot of a UV vector field re-colocated at the
!     mass gridpoint for an horizontal cross-section
!
!!**  METHOD
!!    ------
!!     
!!     Assumption is made that wind components were re-colocated onto the mass
!!   gridpoint location prior to calling IMAGEV_FORDIACHRO. The horizontal coordinates
!!   of the mass gridpoint are first collected, next the gridmap overlay
!!   background and the display window are computed according to user requests,
!!   the visual characteritics of the plot are prescribed, and the wind
!!   arrows are plotted accounting for map projection using the VVECTR NCAR
!!   utility. If IMAGEV_FORDIACHRO works on a constant altitude or pressure level, areas 
!!   where the plotting level intercepts the terrain are hatched and wind 
!!   vector are hidden. Finally, various information labels are printed on
!!   the plot.
!!     
!!     Notice that a TRACE-provided VVUMXY routine is used within the NCAR
!!   vector VVECTR utility to map the wind vectors onto the stretched
!!   MESO-NH model space.  Wind vectors are given in m/s and scaled by VVUMXY
!!   to obtain arrow sizes in "NCAR fractional coordinate" (NCAR User Guide
!!   "Fundamentals", Appendix A, p345 section 1), notice this is different
!!   from what is required for Conpack... The final result is an automatic
!!   arrow scale selection giving a maximum arrow size equal to the meshlength
!!   on the plot. If a different procedure has to be followed VVUMXY should
!!   be updated accordingly. The parameters of the NCAR VVECTR utility can
!!   be printed online by typing "man vectors_params", these feature are not
!!   really documented elsewhere in NCAR user guide.
!!    
!!     Further, notice that the Meso-NH model usually provides the so-called 
!!   covariant wind components in the LFIFM files (multiplied by rho_~_*).
!!   If this assumption is made, the wind modulus of the displayed wind is 
!!   equal to the modulus of the real meteorological wind on the spherical 
!!   earth. 
!!
!!    EXTERNAL
!!    --------
!!      DEFENETRE : when cartesian geometry applies, defines the    !
!!                  display window                                  !
!!      BCGRD     : when a cartographic projection applies, defines !
!!                  displayed                                       !
!!                  window and draws the continent/state outlines   !
!!      GSCLIP    : clips items getting out of the drawing window   ! 
!!      GETSET    : retrieves the normalized and user NCAR          !
!!                  coordinates of a previously used window         ! 
!!      PLCHHQ    : prints high-quality character strings           !
!!                                                                  !
!!      VVSETR !  : gets the value of a NCAR parameter,   REEL      !
!!      VVSETI !                                          INTEGER   !
!!      VVINIT    : initialize a vector plot (arrows)               !
!!      VVECTR    : draws the arrows for a vector plot              !
!!                                                                  !
!!      CPSETI !                                          INTEGER   !
!!      CPSETR !  : sets the value of a NCAR parameter,   REEL      !
!!      CPSETC !                                          CHARACTER ! NCAR
!!                                                                  !
!!      CPGETI !                                          INTEGER   !
!!      CPGETR !  : gets the value of a NCAR parameter,   REEL      !
!!      CPGETC !                                          CHARACTER !
!!                                                                  !
!!      CPRECT    : Conpack initialization (contours)               !
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
!!      CPCLAM    : adds contour in a  previously defined area      ! 
!!      CPRSET    : resets Conpack parameters to default values     !
!!
!!
!!      VVUMXY    : TRACE provided FORTRAN-77 routine directly called
!!                  within the VVECTR NCAR utility to to map the wind
!!                  vectors onto the stretched MESO-NH model space.
!!      CPMPXY    : TRACE provided FORTRAN-77 routine directly called
!!                  within CONPACK to map the array space onto the
!!                  cartographic space
!!      SFILL     : TRACE provided FORTAN-77 routine directly called 
!!                  CONPACK to define the hatched area used to locate
!!                  points  where the plot level intercepts topography
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_TITLE  : Declares heading variables for the plots (TRACE)
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
!!        CTYPHOR    : Horizontal cross-section type
!!                     (='K' --> model level section;
!!                      ='Z' --> constant-altitude section;
!!                      ='P' --> isobar section (planned)
!!                      ='T' --> isentrope section (planned)
!!        XSPVAL     : Special value
!!        NISKIP     : Sampling rate for drawing velocity vectors
!!
!!      Module MODD_OUT       : Defines a log. unit for printing
!!        NIMAXT : x-size of the displayed section of the model array
!!        NJMAXT : y-size of the displayed section of the model array
!!
!!      Module MODD_TIME   ! To be checked, useless..
!!      Module MODD_TIME1  ! To be checked, useless.
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
!!      Updated   PM   13/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_COORD
USE MODD_CONF
USE MODD_GRID
USE MODD_GRID1
USE MODE_GRIDPROJ
USE MODD_TITLE
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_ALLOC_FORDIACHRO
USE MODD_OUT
USE MODN_PARA
USE MODN_NCAR
USE MODD_TIME
USE MODD_TIME1
USE MODD_SUPER
USE MODD_RESOLVCAR
USE MODD_TIT
USE MODD_PVT
USE MODD_MEMCV
USE MODD_CTL_AXES_AND_STYL
USE MODI_RESOLV_TIT
USE MODI_RESOLV_TITY
USE MODI_COMPUTEDIR
!
IMPLICIT NONE

INTERFACE
      SUBROUTINE IMAGE_FORDIACHRO(PTAB,KLREF,PTABINT,KNHI,KNDOT,HTEXTE)
      CHARACTER(LEN=*)   :: HTEXTE
      REAL                :: PTABINT
      REAL,DIMENSION(:,:) :: PTAB
      INTEGER :: KNHI, KNDOT, KLREF
      END SUBROUTINE IMAGE_FORDIACHRO
END INTERFACE
!
!*       0.0   TRACE interface with the "VVUMXY" routine of the NCAR package
!
! NOTICE:  The TRACE provided VVUMXY routine and the NCAR graphical utilities 
! ------   are NOT written in Fortran 90, but in Fortran 77.. This sub-section
!          of TRACE does not follow the Meso-NH usual rules: it has to be made 
!          using a COMMON stack with  static memory allocation of XZZXX and
!          XZZXY arrays.
!
!
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
COMMON/TEMH/XZZX,XZZY,NIIMAX,NIJMAX
#include "big.h"
REAL,DIMENSION(N2DVERTX) :: XZZX
REAL,DIMENSION(N2DVERTX) :: XZZY
INTEGER :: NIIMAX, NIJMAX
LOGICAL :: LVERT, LHOR,LPT, LXABS
!
!*       0.1   NCAR work arrays
!
! See aforementioned notice. The dimensions of these arrays are
! subject to possible tuning, but have to be prescribed. Add
! extra size if necessary.
!
INTEGER,PARAMETER       :: JPLRWK=50000, JPLIWK=50000
INTEGER,PARAMETER       :: JPRSCR=10000, JPISCR=10000
INTEGER,PARAMETER       :: JPMAP=800000, JPAREAGRP=300, JPWRK=50000

REAL,DIMENSION(JPLRWK):: ZRWRK
INTEGER,DIMENSION(JPLIWK):: IWRK
!REAL,DIMENSION(JPRSCR):: ZRSCR
!INTEGER,DIMENSION(JPISCR):: ISCR
INTEGER,DIMENSION(JPMAP):: IIMAP
INTEGER,DIMENSION(JPAREAGRP):: IAREA, IGRP
REAL,DIMENSION(JPWRK)   :: ZXWRK, ZYWRK
!
!*       0.2   Dummy arguments and results
!
INTEGER                 :: KLREF  ! Cross-section altitude (or Model Level
                                  ! or Pressure depending on user's vertical
                                  ! coordinate choice)
CHARACTER(LEN=*) :: HTEXTE       ! Plot heading contataining field name
REAL,DIMENSION(:,:) :: PU,PV      ! Arrays of "wind components" to be plotted
!
!*       0.3   Local variables
!
INTEGER :: JILOOP, JJLOOP, IUB1, IUB2, ID, J, IJ, JA
INTEGER                 :: ICL

INTEGER                 :: IZS

INTEGER                 :: ITER, JTER, ISKIP, IGRNC
INTEGER                 :: II, INUM, IRESP, ILOOP, IDEB, IFIN
INTEGER                 :: JLOOPI, JLOOPJ

CHARACTER(LEN=70) ::   YPLANH, YTEM 
CHARACTER(LEN=40) ::   YTEXTE
CHARACTER(LEN=4)  ::   YTE, YC4, YC42
!
REAL :: ZLREF, ZZSPVAL, ZY, ZINTX, ZINTY
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZZU,ZZV
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZSTRU,ZSTRV
!REAL,DIMENSION(SIZE(PU,1),SIZE(PU,2)) :: ZZU, ZZV
REAL,DIMENSION(:),ALLOCATABLE,SAVE ::  ZZY, ZTEMX,ZTEMY
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZX, ZLAT, ZLON, ZYY
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZDIRU, ZDIRV, ZLA, ZLO
REAL :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
REAL :: ZVINT, ZVY
REAL :: ZXPOSTITT1, ZXYPOSTITT1
REAL :: ZXPOSTITT2, ZXYPOSTITT2
REAL :: ZXPOSTITT3, ZXYPOSTITT3
REAL :: ZXPOSTITB1, ZXYPOSTITB1
REAL :: ZXPOSTITB2, ZXYPOSTITB2
REAL,SAVE :: ZXPOSTITB3, ZXYPOSTITB3
REAL,DIMENSION(18) :: ZCOL

REAL,DIMENSION(:),ALLOCATABLE,SAVE ::  ZSTR1

INTEGER,DIMENSION(18) :: ICOL
INTEGER :: ICOL1,IER
LOGICAL,SAVE       :: GVSUPSCA
!
!*       0.4   External for NCAR use
!
! SFILL subroutine declared as external provides area control
! in some parts of the contour plot.
!
EXTERNAL SFILL
EXTERNAL STUMXY
!
!-------------------------------------------------------------------------------
!
!*       1.    DISPLAY ENVIRONMENT SETUP AND ARROWS PLOTTING
!              ---------------------------------------------
!
!*       1.1   Array sizes calculation and default field value
!
IUB1=UBOUND(PU,1)
IUB2=UBOUND(PU,2)
IF(ALLOCATED(ZZU))THEN
  DEALLOCATE(ZZU)
ENDIF
IF(ALLOCATED(ZZV))THEN
  DEALLOCATE(ZZV)
ENDIF
ALLOCATE(ZZU(SIZE(PU,1),SIZE(PU,2)),ZZV(SIZE(PU,1),SIZE(PU,2)))

!
!DO JJLOOP=1,NJMAXT
DO JJLOOP=1,IUB2
! DO JILOOP=1,NIMAXT
  DO JILOOP=1,IUB1
    ZZU(JILOOP,JJLOOP)=XSPVAL
    ZZV(JILOOP,JJLOOP)=XSPVAL
  ENDDO
ENDDO
!
!*       1.2  Collects XHAT and YHAT values at mass gridpoints (NGRID=1)
!*            where wind components have been relocated in TRACEH
!
DO JILOOP=NIINF,NISUP
  XZZX(JILOOP-NIINF+1)=XXX(JILOOP,1)
!XZZX(JILOOP-NIINF+1)=XXX(JILOOP,NMGRID)
ENDDO
DO JJLOOP=NJINF,NJSUP
  XZZY(JJLOOP-NJINF+1)=XXY(JJLOOP,1)
!XZZY(JJLOOP-NJINF+1)=XXY(JJLOOP,NMGRID)
ENDDO
!
!*       1.3  Collects wind values within the user postprocessing
!*            window with a sampling rate of NISKIP outside values 
!*            are kept to default
!
!DO JJLOOP=1,NJMAXT,NISKIP
DO JJLOOP=1,IUB2,NISKIP
! DO JILOOP=1,NIMAXT,NISKIP
  DO JILOOP=1,IUB1,NISKIP
    ZZU(JILOOP,JJLOOP)=PU(JILOOP,JJLOOP)
    ZZV(JILOOP,JJLOOP)=PV(JILOOP,JJLOOP)
  ENDDO
ENDDO
!!!!!!!!!!!!!!!STREAM
IF(LSTREAM)THEN
  ITER=IUB1/NISKIP+1
  IF(1+(ITER-1)*NISKIP > IUB1)ITER=ITER-1
  JTER=IUB2/NISKIP+1
  IF(1+(JTER-1)*NISKIP > IUB2)JTER=JTER-1
  ALLOCATE(ZDIRU(ITER,JTER),ZDIRV(ITER,JTER))
  ALLOCATE(ZX(ITER,1),ZZY(JTER))
  ALLOCATE(ZTEMX(IUB1),ZTEMY(IUB2))
  ZTEMX(1:IUB1)=XZZX(1:IUB1)
  ZTEMY(1:IUB2)=XZZY(1:IUB2)
  ZX(:,1)=XZZX(1:IUB1:NISKIP)
  ZZY=XZZY(1:IUB2:NISKIP)
  ZDIRU=PU(1:IUB1:NISKIP,1:IUB2:NISKIP)
  ZDIRV=PV(1:IUB1:NISKIP,1:IUB2:NISKIP)
! print *,' **deallocate ZZU ZZV'
   ALLOCATE(ZSTRU(ITER,JTER),ZSTRV(ITER,JTER))
  
  DO JJLOOP=1,JTER
  DO JILOOP=1,ITER
    ZSTRU(JILOOP,JJLOOP)=ZDIRU(JILOOP,JJLOOP)
    ZSTRV(JILOOP,JJLOOP)=ZDIRV(JILOOP,JJLOOP)
  ENDDO
  ENDDO
  XZZX(1:ITER)=ZX(:,1)
  XZZY(1:JTER)=ZZY
! IUB1=ITER
! IUB2=JTER
  DEALLOCATE(ZDIRU,ZDIRV,ZX,ZZY)
!!!!!!!!!!!!!!!STREAM
ALLOCATE(ZSTR1(4*ITER*JTER))
!!!!!!!!!!!!!!!STREAM
ENDIF
!!!!!!!!!!!!!!!STREAM
!
IF(LDIRWIND)THEN
  ISKIP=NISKIP
  NISKIP=1
  IGRNC=NIGRNC
  NIGRNC=5
ENDIF

!000000000000000000000000000000000000000000000000000000000000000
IF(LDIRWIND)THEN
!000000000000000000000000000000000000000000000000000000000000000
  print *,' imagev LDIRWIND ',LDIRWIND
  YTEXTE(1:LEN(YTEXTE))=' '
  YTEXTE='WIND-DIRECTION'
  YTEXTE=ADJUSTL(YTEXTE)
  ITER=IUB1/NISKIP+1
  IF(1+(ITER-1)*NISKIP > IUB1)ITER=ITER-1
  JTER=IUB2/NISKIP+1
  IF(1+(JTER-1)*NISKIP > IUB2)JTER=JTER-1
  ALLOCATE(ZDIRU(ITER,JTER),ZDIRV(ITER,JTER))
  ALLOCATE(ZX(ITER,1),ZZY(JTER))
  ZX(:,1)=XZZX(1:IUB1:NISKIP)
  ZZY=XZZY(1:IUB2:NISKIP)
  ZDIRU=PU(1:IUB1:NISKIP,1:IUB2:NISKIP)
  ZDIRV=PV(1:IUB1:NISKIP,1:IUB2:NISKIP)
   print*,'imagev dd ',minval(ZDIRU),maxval(ZDIRU),minval(ZDIRV), maxval(ZDIRV)
  CALL COMPUTEDIR(ITER,JTER,IUB1,IUB2,NISKIP,ZDIRU,ZDIRV)
   print*,'imagev dd ',minval(ZDIRV), maxval(ZDIRV)
!! Supprime en nov 2001 Appel routine COMPUTEDIR
!! Supprime en nov 2001 Appel routine COMPUTEDIR
  IF(LSUPER)THEN
    NSUPER=NSUPER+1
    print *,' ** imagev DIRWIND NSUPER ',NSUPER
    IF(NSUPER == 1)THEN
      IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(1)   
      IF(LCARTESIAN)CALL DEFENETRE
    END IF
  ELSE
    IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(1)   
    IF(LCARTESIAN)CALL DEFENETRE
  END IF
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! CALL SET(ZVL,ZVR,ZVB,ZVT,ZX(1,1),ZX(ITER,1),ZZY(1),ZZY(JTER),1)

  CALL TABCOL_FORDIACHRO
  IJ=1
  DO J=15,345,30
    IJ=IJ+1
    ZCOL(IJ)=J
  ENDDO
  ZCOL(1)=0.
  IJ=IJ+1
  ZCOL(IJ)=360.
  ICOL(1)=4; ICOL(13)=4; ICOL(2)=88; ICOL(3)=79; ICOL(4)=7
  ICOL(5)=52; ICOL(6)=25; ICOL(7)=2; ICOL(8)=20; ICOL(9)=24
  ICOL(10)=3; ICOL(11)=124; ICOL(12)=5; ICOL(13)=4
  DO JJLOOP=1,JTER
    DO JILOOP=1,ITER
	IF(ZDIRV(JILOOP,JJLOOP) == XSPVAL)THEN
!       print *,J,' CYCLE  ',ZDIRV(JILOOP,JJLOOP),ZCOL(J),ZCOL(J-1)
	  CYCLE
	ENDIF
      DO J=2,IJ
!       print *,J,' ',ZDIRV(JILOOP,JJLOOP),ZCOL(J),ZCOL(J-1)
        
	IF(ZDIRV(JILOOP,JJLOOP) == 0. .OR. ZDIRV(JILOOP,JJLOOP) == 360.)THEN
	  CALL GSPMCI(ICOL(1))
!     print *,' ZDIRV(JILOOP,JJLOOP) J+2 ',ZDIRV(JILOOP,JJLOOP),ICOL(1)
	  EXIT
	ELSE IF(ZDIRV(JILOOP,JJLOOP) < ZCOL(J).AND. &
		ZDIRV(JILOOP,JJLOOP) >= ZCOL(J-1))THEN
	  CALL GSPMCI(ICOL(J-1))
!     print *,' ZDIRV(JILOOP,JJLOOP) J+1 ',ZDIRV(JILOOP,JJLOOP),ICOL(J)
	  EXIT
	ENDIF
      ENDDO
      CALL GSMK(2)
      ZINTX=ZX(JILOOP,1)
      ZINTY=ZZY(JJLOOP)
      CALL GPM(1,ZINTX,ZINTY)
      CALL GSMK(3)
      CALL GPM(1,ZINTX,ZINTY)
      CALL GSMK(5)
      CALL GPM(1,ZINTX,ZINTY)
    ENDDO
  ENDDO
! 
! Legende couleurs
  CALL GSCLIP(0)
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,1)
  ZVINT=(ZVT-ZVB)/12.
  ZVY=ZVB
  YTE='    '
  WRITE(YTE,'(F4.0)')ZCOL(1)
  CALL PLCHHQ(ZVR+.0002,ZVY,YTE,.012,0.,-1.)
! print *,'ZVR,ZVY,YTE ',ZVR,ZVY,YTE
  DO J=1,6
    CALL GSPMCI(ICOL(1))
    ZINTX=ZVR+.005*J
    ZINTY=ZVY+.015
    CALL GSMK(2)
    CALL GPM(1,ZINTX,ZINTY)
    CALL GSMK(3)
    CALL GPM(1,ZINTX,ZINTY)
    CALL GSMK(5)
    CALL GPM(1,ZINTX,ZINTY)
  ENDDO
  ZVY=ZVY+ZVINT/2.
  YTE='    '
  WRITE(YTE,'(F4.0)')ZCOL(2)
  CALL PLCHHQ(ZVR+.0002,ZVY,YTE,.012,0.,-1.)
! print *,'ZVR,ZVY,YTE ',ZVR,ZVY,YTE
  DO J=1,6
    CALL GSPMCI(ICOL(2))
    ZINTX=ZVR+.005*J
    ZINTY=ZVY+.015
    CALL GSMK(2)
    CALL GPM(1,ZINTX,ZINTY)
    CALL GSMK(3)
    CALL GPM(1,ZINTX,ZINTY)
    CALL GSMK(5)
    CALL GPM(1,ZINTX,ZINTY)
  ENDDO
  DO J=3,13
    ZVY=ZVY+ZVINT
    YTE='    '
    WRITE(YTE,'(F4.0)')ZCOL(J)
    CALL PLCHHQ(ZVR+.0002,ZVY,YTE,.012,0.,-1.)
!   print *,'ZVR,ZVY,YTE ',ZVR,ZVY,YTE
    DO JA=1,6
      CALL GSPMCI(ICOL(J))
      ZINTX=ZVR+.005*JA
      ZINTY=ZVY+.015
      CALL GSMK(2)
      CALL GPM(1,ZINTX,ZINTY)
      CALL GSMK(3)
      CALL GPM(1,ZINTX,ZINTY)
      CALL GSMK(5)
      CALL GPM(1,ZINTX,ZINTY)
    ENDDO
  ENDDO
  ZVY=ZVY+ZVINT/2.
  YTE='    '
  WRITE(YTE,'(F4.0)')ZCOL(14)
  CALL PLCHHQ(ZVR+.0002,ZVY,YTE,.012,0.,-1.)
!
! Titre N1 TOP
!! WRITE(YPLANH,1001)NIINF,NISUP,NJINF,NJSUP
!!  ZXPOSTITT1=.002
!!  ZXYPOSTITT1=.98
!!  IF(XPOSTITT1 /= 0.)THEN
!!    ZXPOSTITT1=XPOSTITT1
!!  ENDIF
!!  IF(XYPOSTITT1 /= 0.)THEN
!!    ZXYPOSTITT1=XYPOSTITT1
!!  ENDIF
!!  CALL RESOLV_TIT('CTITT1',YPLANH)
!!  IF(YPLANH /= ' ')THEN
!!    IF(XSZTITT1 /= 0.)THEN
!!      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YPLANH,XSZTITT1,0.,-1.)
!!!     CALL PLCHHQ(0.002,0.98,YPLANH,XSZTITT1,0.,-1.)
!!    ELSE
!!      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YPLANH,.012,0.,-1.)
!!!     CALL PLCHHQ(0.002,0.98,YPLANH,.012,0.,-1.)
!!    ENDIF
!!  ENDIF
!!  IF(LDATFILE)CALL DATFILE_FORDIACHRO

  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(2)
  if(nverbia > 0)then
    print *,'**imagev AP CALL BCGRD_FORDIACHRO(2) 1 '
  endif
  CALL TABCOL_FORDIACHRO
 
  IF(LPRINT)THEN
    CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
    IF(IRESP /= 0)THEN
      CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
      OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
      PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
    ENDIF
    ILOOP=SIZE(ZDIRV,1)/5
    IF(ILOOP * 5 < SIZE(ZDIRV,1))ILOOP=ILOOP+1
    WRITE(INUM,'(''CH  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
  & CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    IF(LMINUS .OR. LPLUS)THEN
      WRITE(INUM,'(A55,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITB3(1:55)
    ELSE
      WRITE(INUM,'(''WIND-DIRECTION'',26X,''(NIINF-NISUP,NJINF-NJSUP)'')')
  !   WRITE(INUM,'(A40,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITGAL
    ENDIF
    WRITE(INUM,'(''niinf'',i4,'' njinf'',i4,'' nisup'',i4,'' njsup'',i4,&
  &''   '',A1,'' '',i6)')&
    &NIINF,NJINF,NISUP,NJSUP,CTYPHOR,KLREF
    WRITE(INUM,'(''NBVAL en I '',i4,''  NBVAL en J '',i4,''   iter'',i3)') &
    &NISUP-NIINF+1,NJSUP-NJINF+1,ILOOP
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**IMAGEV XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
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
        IFIN=SIZE(ZDIRV,1)+NIINF-1
      ENDIF
      
      WRITE(INUM,'(1X,78(1H*))')
      WRITE(INUM,'('' J   I-> '',3X,I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',78(1H*))')
      DO JLOOPJ=SIZE(ZDIRV,2),1,-1
        WRITE(INUM,'(I4,2X,5(1X,E14.7))')JLOOPJ+NJINF-1,(ZDIRV(II,JLOOPJ),II=IDEB-NIINF+1,IFIN-NIINF+1)

  !     WRITE(INUM,'(I3,2X,5E15.8)')JLOOPJ+NJINF-1,(ZDIRV(II,JLOOPJ),II=IDEB-NIINF+1,IFIN-NIINF+1)
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
      WRITE(INUM,'(''WIND-DIRECTION'',26X,''(NIINF-NISUP,NJINF-NJSUP)'')')
  !   WRITE(INUM,'(A40,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITGAL
    ENDIF
    WRITE(INUM,'(''niinf'',i4,'' njinf'',i4,'' nisup'',i4,'' njsup'',i4,&
  &'' '',A1,'' '',i6)')&
    &NIINF,NJINF,NISUP,NJSUP,CTYPHOR,KLREF
    WRITE(INUM,'(''NBVAL en I '',i4,''  NBVAL en J '',i4)') &
    &NISUP-NIINF+1,NJSUP-NJINF+1
  
    II=MAX(SIZE(ZDIRV,1),SIZE(ZDIRV,2))
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
  	YC4,JLOOPJ,XZZX(JLOOPJ),YC42,JLOOPJ,XZZY(JLOOPJ)
  	YC4='    '
  	YC42='    '
  	WRITE(YC4,'(I4,'')'')')NISUP
  	WRITE(YC42,'(I4,'')'')')NJSUP
      ELSE
  	IF(SIZE(ZDIRV,1) > SIZE(ZDIRV,2))THEN
  	  IF(JLOOPJ < SIZE(ZDIRV,2))THEN
  	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZX(JLOOPJ), &
  	    JLOOPJ,XZZY(JLOOPJ)
  	  ELSE IF(JLOOPJ == SIZE(ZDIRV,1))THEN
  	    WRITE(INUM,'(''NISUP('',A4,I4,5X,E15.8)')YC4,JLOOPJ,XZZX(JLOOPJ)
              WRITE(INUM,'(1X,73(1H*))')
  	  ELSE IF(JLOOPJ == SIZE(ZDIRV,2))THEN
  	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,''NJSUP('',A4,I4,5X,E15.8)')&
  	    JLOOPJ,XZZX(JLOOPJ), &
  	    YC42,JLOOPJ,XZZY(JLOOPJ)
  	  ELSE IF(JLOOPJ > SIZE(ZDIRV,2))THEN
  	    WRITE(INUM,'(5X,I9,5X,E15.8)')JLOOPJ,XZZX(JLOOPJ)
  	  ENDIF
  	ELSE IF(SIZE(ZDIRV,2) > SIZE(ZDIRV,1))THEN
  	  IF(JLOOPJ < SIZE(ZDIRV,1))THEN
  	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZX(JLOOPJ), &
  	    JLOOPJ,XZZY(JLOOPJ)
  	  ELSE IF(JLOOPJ == SIZE(ZDIRV,2))THEN
  	    WRITE(INUM,'(29X,5X,5X,''NJSUP('',A4,I4,5X,E15.8)') &
  	    YC42,JLOOPJ,XZZY(JLOOPJ)
              WRITE(INUM,'(1X,73(1H*))')
  	  ELSE IF(JLOOPJ > SIZE(ZDIRV,1))THEN
  	    WRITE(INUM,'(29X,5X,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZY(JLOOPJ)
  	  ELSE
  	    WRITE(INUM,'(''NISUP('',A4,I4,5X,E15.8,5X,5X,I9,5X,E15.8)') &
  	    YC4,JLOOPJ,XZZX(JLOOPJ), &
  	    JLOOPJ,XZZY(JLOOPJ)
  	  ENDIF
  	ELSE
  	  IF(JLOOPJ == SIZE(ZDIRV,2))THEN
  	    WRITE(INUM,'(''NISUP('',A4,I4,5X,E15.8,5X,''NJSUP('',A4,I4,5X,E15.8)') &
  	    YC4,JLOOPJ,XZZX(JLOOPJ), &
  	    YC42,JLOOPJ,XZZY(JLOOPJ)
              WRITE(INUM,'(1X,73(1H*))')
  	  ELSE
  	    WRITE(INUM,'(5X,I9,5X,E15.8,5X,5X,I9,5X,E15.8)')JLOOPJ,XZZX(JLOOPJ), &
  	    JLOOPJ,XZZY(JLOOPJ)
  	  ENDIF
  	ENDIF
      ENDIF
    ENDDO
  ENDIF

  NISKIP=ISKIP
  NIGRNC=IGRNC
  DEALLOCATE(ZX,ZZY,ZDIRU,ZDIRV)
! DEALLOCATE(ZX,ZZY,ZYY,ZLAT,ZLON,ZLA,ZLO,ZDIRU,ZDIRV)
  IF(ALLOCATED(ZYY))DEALLOCATE(ZYY)
  IF(ALLOCATED(ZLAT))DEALLOCATE(ZLAT)
  IF(ALLOCATED(ZLON))DEALLOCATE(ZLON)
  IF(ALLOCATED(ZLA))DEALLOCATE(ZLA)
  IF(ALLOCATED(ZLO))DEALLOCATE(ZLO)
       
!000000000000000000000000000000000000000000000000000000000000000
ELSE
!000000000000000000000000000000000000000000000000000000000000000
!
!*       1.4  Selects display window as requested by LCARTESIAN
!*            Sets Map projection, overlays coastlines and landmarks
!*            if required
!
!
  CALL GSLN(1)
  CALL GSPLCI(1)
  CALL GSTXCI(1)

  IF(LSUPER)THEN
    NSUPER=NSUPER+1
!   print *,' ** imagev NSUPER ',NSUPER

    IF(NSUPER == 1)THEN
      NCOLUVG=NCOLUV1
    ELSEIF(NSUPER == 2)THEN
      NCOLUVG=NCOLUV2
    ELSEIF(NSUPER == 3)THEN
      NCOLUVG=NCOLUV3
    ELSEIF(NSUPER == 4)THEN
      NCOLUVG=NCOLUV4
    ELSEIF(NSUPER == 5)THEN
      NCOLUVG=NCOLUV5
    ELSE
      NCOLUVG=1
    ENDIF
    IF(NSUPER == 1)THEN
      IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(2)   
      IF(LCARTESIAN)CALL DEFENETRE
      if(nverbia > 0)then
        print *,' **imagev AP CALL BCGRD_FORDIACHRO(2) 2 '
      endif
    ENDIF
  ELSE
    IF(.NOT.LCARTESIAN)CALL BCGRD_FORDIACHRO(2)   
    IF(LCARTESIAN)CALL DEFENETRE
    NCOLUVG=NCOLUV1
  ENDIF
!
!*       1.5  Routine VVUMXY of provided by TRACE to locate and scale wind
!*            arrows on the display
!
  LHOR=LHORIZ
  LVERT=LVERTI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STREAM
! GO TO 1000
IF(.NOT.LSTREAM)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STREAM
  CALL VVSETI('MAP',4)
  CALL VVSETI('SET',0)
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  if(nverbia > 0)then
    print *,' **imagev ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT
  endif
  CALL VVSETR('VPL',ZVL)    
  CALL VVSETR('VPR',ZVR)
  CALL VVSETR('VPB',ZVB)
  CALL VVSETR('VPT',ZVT)
  !CALL VVSETR('WDL',100000.)
  CALL VVSETR('WDL',ZWL)
  !CALL VVSETR('WDR',2500000.)
  CALL VVSETR('WDR',ZWR)
  CALL VVSETR('WDB',ZWB)
  CALL VVSETR('WDT',ZWT)
  
! CALL SET(ZVL,ZVR,ZVB,ZVT,100000.,2500000.,ZWB,ZWT,ID)
! Sortie statistiques
  IF(LVST)THEN
    CALL VVSETI('VST',1)
  ELSE
    CALL VVSETI('VST',0)
  ENDIF
  CALL VVSETR('AMX',XAMX)
  CALL VVSETR('VHC',XVHC)
  CALL VVSETR('VRL',XVRL)
  CALL VVSETR('VLC',XVLC)
  IF(XVHC < 0. )THEN
    CALL VVSETC('MXT',' ')
    CALL VVSETC('MXT','Scale')
  END IF
!
!*      1.6   Masks vectors where wind coponents have XSPVAL values
!
  CALL VVSETI('SVF',3)
  CALL VVSETR('USV',XSPVAL)
  CALL VVSETR('VSV',XSPVAL)
!
!*      1.6   Selects look and feel options for the vector display
!             (Text strings, etc..)
!
  if(nverbia > 0)then
    print *,' **imagev AP VVSETR(VSV,XSPVAL)'
  endif
  CALL VVSETI('MNP',-4)
  CALL VVSETR('MNX',(-ZVL+.002)/(ZVR-ZVL))
!
! ZY=-1./5.
! IF(ZVB-(ZVT-ZVB)/5..LT.0.05)ZY=(0.05-ZVB)/(ZVT-ZVB)
! Oct 2000 Essai de repositionnement des fleches min et max
  IF(ZVB <= .1)THEN
    ZY=(-ZVB+0.0395)/(ZVT-ZVB)
  ELSE
    ZY=(-ZVB+0.0545)/(ZVT-ZVB)
  ENDIF
  CALL VVSETR('MNY',ZY)
  CALL VVSETI('MXP',-4)
  CALL VVSETR('MXX',(-ZVL+.14+.002)/(ZVR-ZVL))
  CALL VVSETR('MXY',ZY)
  CALL VVSETR('MXS',.008*.9/(ZVR-ZVL))
! CALL VVSETR('MXS',.008)
  CALL VVSETR('MNS',.008*.9/(ZVR-ZVL))
! CALL VVSETR('MNS',.008)
! Elimination de la legende des fleches si LEGVECT=F
  IF(.NOT.LEGVECT)THEN
    CALL VVSETC('MXT',' ')
    CALL VVSETC('MNT',' ')
  ENDIF
  IF(XVHC >= 0.)THEN
! Janv 2001
    GVSUPSCA=LVSUPSCA
    LVSUPSCA=.FALSE.
  ENDIF
!
!*     1.7    Draws the arrows
!
  IF(XLWV > 0.)THEN
    CALL VVSETR('LWD',XLWV)
  ELSE
    CALL VVSETR('LWD',XLWVDEF)
  ENDIF
  CALL GSCLIP(0)                                     ! Clipping off
  CALL VVSETI('VPO',1)
! CALL GSCLIP(1)                                     ! Clipping off
! if(nverbia > 0)then
! Oct 2000 La ligne suivante est obligatoire sinon plantage avec visu
! dans certains cas -> besoin de revenir sur le pb un jour
  print *,' **imagev AV VVINIT '
!endif
  CALL VVINIT(ZZU,IUB1,ZZV,IUB1,0.,0,IUB1,IUB2,0.,0) ! Initializes VVECTR
  CALL VVECTR(ZZU,ZZV,0.,0,0,0.)                     ! Draws arrows
  CALL GSCLIP(1)                                     ! Clipping back on
  CALL GSLWSC(1.)
  CALL VVRSET
  if(nverbia > 0)then
    print *,' **imagev AP VVRSET '
  endif
! Janv 2001
  IF(XVHC >= 0.)THEN
    LVSUPSCA=GVSUPSCA
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STREAM
!1000 CONTINUE
ELSE
NIIMAX=ITER
!NIIMAX=NIMAXT
NIJMAX=JTER
  CALL STSETI('MAP',4)
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  if(nverbia > 0)then
    print *,' **imagev ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT,NIIMAX,NIJMAX
    print *,' **imagev ap getset ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT,NIIMAX,NIJMAX
  endif
  CALL STSETI('SET',0)
  CALL STSETR('VPL',ZVL)    
  CALL STSETR('VPR',ZVR)
  CALL STSETR('VPB',ZVB)
  CALL STSETR('VPT',ZVT)
  CALL STSETR('WDL',ZWL)
  CALL STSETR('WDR',ZWR)
  CALL STSETR('WDB',ZWB)
  CALL STSETR('WDT',ZWT)
  if(nverbia > 0)then
    print *,' **imagev ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT
  endif
  
  CALL STSETI('AGD',NARSTR)
  CALL STSETI('GBS',0)
  CALL STSETI('CPM',0)
! CALL STSETR('ARL',.009)
  CALL STSETR('ARL',XARLSTR)
  CALL STSETR('DFM',.02)
  CALL STSETR('CDS',1.)
  CALL STSETR('SSP',XSSP)
! CALL STSETR('SSP',.004)
  CALL STSETR('LWD',XLWSTR)
  CALL STSETI('MSK',0)
  CALL STSETI('SVF',3)
  CALL STSETR('USV',XSPVAL)
  CALL STSETR('VSV',XSPVAL)
  CALL GQPLCI(IER,ICOL1)
  CALL GSPLCI(NCOLUVG)
  IZS=4*ITER*JTER
  CALL STINIT(ZSTRU,ITER,ZSTRV,ITER,0.,0,ITER,JTER,ZSTR1,IZS) ! Initializes VVECTR
! CALL STINIT(ZSTRU,ITER,ZSTRV,ITER,ZTEM,ITER,ITER,JTER,ZSTR1,IZS) ! Initializes VVECTR
  CALL STREAM(ZSTRU,ZSTRV,0.,0,STUMXY,ZSTR1)                     ! Draws arrows
  CALL STRSET
  CALL GSPLCI(ICOL1)
  XZZX(1:IUB1)=ZTEMX(1:IUB1)
  XZZY(1:IUB2)=ZTEMY(1:IUB2)
  DEALLOCATE(ZSTR1,ZSTRU,ZSTRV,ZTEMX,ZTEMY)
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STREAM
!
!000000000000000000000000000000000000000000000000000000000000000
ENDIF
!000000000000000000000000000000000000000000000000000000000000000
!------------------------------------------------------------------------------
!
!*     2.  TOPOGRAPHY MASKING WHEN PLOTTED LEVEL INTERCEPTS TERRAIN
!          --------------------------------------------------------
!
!
!*     2.1  Initialization of a topographic mask using
!*          the NCAR "area" features (see NCAR manual)
!
LVERT=LVERTI
LHOR=LHORIZ
if(nverbia >0)then
  print *,' **imagev LVERT, LHOR ',LVERT,LHOR
endif
CALL CPSETI('MAP',4)
CALL CPSETI('SET',0)
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
NIIMAX=IUB1
!NIIMAX=NIMAXT
NIJMAX=IUB2
!NIJMAX=NJMAXT
!print *,' **NIIMAX,NIJMAX ',NIIMAX,NIJMAX
!           
IF(CTYPHOR.EQ.'Z' .AND. (.NOT.LSUPER .OR. NSUPER == 1))THEN
  ZLREF=KLREF
  !
  DO JILOOP=NIINF,NISUP
     DO JJLOOP=NJINF,NJSUP
        !                      If terrain higher than topography  
        !                      a 888. mask value is forced
        !
        IF(ZLREF.LT.XXZS(JILOOP,JJLOOP,1))PU(JILOOP-NIINF+1,JJLOOP-NJINF+1)=888.
     ENDDO
  ENDDO
  !
  ICL=1                        ! A single contour is drawn
  CALL CPSETI('CLS',0)         ! Contour value forced
  CALL CPSETI('HCF',1)         ! All contoured areas will be hatched
  CALL CPSETC('CFT',' ')       ! No 'CONSTANT FIELD' message
  CALL CPSETI('NCL',ICL)       ! A single contour is drawn
  CALL CPSETI('PAI',ICL)       ! A single contour is drawn
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
    WHERE(PU(:,:)/=888.)PU(:,:)=ZZSPVAL
    WHERE(PU(:,::2)==888.)PU(:,::2)=PU(:,::2)+1.E-3
  CALL CPSETR('SPV',ZZSPVAL)    ! Valeur speciale = ZZSPVAL
!
!*      2.2    Effective area computation and contour drawing
!
  CALL ARINAM(IIMAP,JPMAP)                               !Initialize areas
  CALL CPRECT(PU,IUB1,IUB1,IUB2,ZRWRK,JPLRWK,IWRK,JPLIWK)!Initialize conpack
  CALL CPCLAM(PU,ZRWRK,IWRK,IIMAP)                       !Contours terrain area
  CALL CPCLDR(PU,ZRWRK,IWRK)                             !Contours outside field
  CALL ARSCAM(IIMAP,ZXWRK,ZYWRK,JPWRK,IAREA,IGRP,JPAREAGRP,SFILL)!Hatches
  !                                                              !terrain area
ENDIF
!
!-----------------------------------------------------------------------------
!
!*    3.    COMPLETING THE PLOT
!           -------------------
!
!*    3.1   Page information labels
!

CALL GSCLIP(0)
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,1)
if(nverbia > 0)then
  print *,' **imagev 2 ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT
endif

IF(NLOOPSUPER == 1)THEN
  CALL RESOLV_TIT('CTITVAR1',HTEXTE)
ELSE IF(NLOOPSUPER == 2)THEN
  CALL RESOLV_TIT('CTITVAR2',HTEXTE)
ELSE IF(NLOOPSUPER == 3)THEN
  CALL RESOLV_TIT('CTITVAR3',HTEXTE)
ELSE IF(NLOOPSUPER == 4)THEN
  CALL RESOLV_TIT('CTITVAR4',HTEXTE)
ELSE IF(NLOOPSUPER == 5)THEN
  CALL RESOLV_TIT('CTITVAR5',HTEXTE)
ELSE IF(NLOOPSUPER == 6)THEN
  CALL RESOLV_TIT('CTITVAR6',HTEXTE)
ELSE IF(NLOOPSUPER == 7)THEN
  CALL RESOLV_TIT('CTITVAR7',HTEXTE)
ELSE IF(NLOOPSUPER == 8)THEN
  CALL RESOLV_TIT('CTITVAR8',HTEXTE)
ENDIF

IF(.NOT.LSUPER)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(HTEXTE /= ' ')THEN
  CALL PLCHHQ(MAX(ZVR,.99),0.007,HTEXTE(1:LEN_TRIM(HTEXTE)),.011,0.,+1.)
  ENDIF
! CALL PLCHHQ(ZVR-(ZVR-ZVL)/4.,0.007,HTEXTE,.011,0.,-1.)
ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(HTEXTE /= ' ')THEN
  CALL PLCHHQ(MAX(ZVR,.99),0.007+(NSUPER-1)*.017,HTEXTE(1:LEN_TRIM(HTEXTE)),.009,0.,+1.)
  ENDIF
! CALL PLCHHQ(ZVR-(ZVR-ZVL)/4.,0.007+(NSUPER-1)*.017,HTEXTE,.009,0.,-1.)
ENDIF

IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN

  CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)

! Modifs for diachro
! Titres en X
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXL',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXL',YTEM)
    IF(XSZTITXL /= 0.)THEN
      CALL PLCHHQ(ZVL,ZVB/2.,YTEM,XSZTITXL,0.,-1.)
    ELSE
      CALL PLCHHQ(ZVL,ZVB/2.,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXM',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXM',YTEM)
    IF(XSZTITXM /= 0.)THEN
      CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),XSZTITXM,0.,0.)
!     CALL PLCHHQ((ZVL+ZVR)/2.-ZVB/2.,ZVB/2.,YTEM,XSZTITXM,0.,-1.)
    ELSE
      CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
!     CALL PLCHHQ((ZVL+ZVR)/2.-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXR',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXR',YTEM)
    IF(XSZTITXR /= 0.)THEN
      CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,XSZTITXR,0.,-1.)
    ELSE
      CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
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
! Top2
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
  IF(YTEM /= ' ')THEN
    IF(XSZTITT2 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,XSZTITT2,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YTEM,XSZTITT2,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,.008,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF
! Top3
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITT3',YTEM)
  ZXPOSTITT3=.002
  ZXYPOSTITT3=.93
  IF(XPOSTITT3 /= 0.)THEN
    ZXPOSTITT3=XPOSTITT3
  ENDIF
  IF(XYPOSTITT3 /= 0.)THEN
    ZXYPOSTITT3=XYPOSTITT3
  ENDIF
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
  CALL RESOLV_TIT('CTITB1',CLEGEND)
  ZXPOSTITB1=.002
  ZXYPOSTITB1=.005
  IF(XPOSTITB1 /= 0.)THEN
    ZXPOSTITB1=XPOSTITB1
  ENDIF
  IF(XYPOSTITB1 /= 0.)THEN
    ZXYPOSTITB1=XYPOSTITB1
  ENDIF
  IF(CLEGEND /= ' ')THEN
    CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,CLEGEND,.007,0.,-1.)
!   CALL PLCHHQ(0.002,0.005,CLEGEND,.007,0.,-1.)
  ENDIF
! Titre N2 BOTTOM
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
    CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,.007,0.,-1.)
!   CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
  ENDIF
! Titre N3 BOTTOM
  YTEM(1:LEN(YTEM))=' '
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
       if(nverbia > 0)then
         print *,' imagev  CTITB3MEM ',CTITB3MEM(1:LEN_TRIM(CTITB3MEM))
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
!     print *,' **imagev CTITB3 AV RESOLV_TIT ',CTITB3
      CALL RESOLV_TIT('CTITB3',CTITB3)
!     print *,' **imagev CTITB3 AP RESOLV_TIT ',CTITB3
      IF(CTITB3 /= ' ')THEN
        CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3,.009,0.,-1.)
      ENDIF
    ENDIF
  ELSE
    CALL RESOLV_TIT('CTITB3',YTEM)
    IF(YTEM /= ' ')THEN
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,.009,0.,-1.)
!     CALL PLCHHQ(0.002,0.050,YTEM,.009,0.,-1.)
    ENDIF
  ENDIF

! Titre N1 TOP
! Top1
  WRITE(YPLANH,1001)NIINF,NISUP,NJINF,NJSUP
  ZXPOSTITT1=.002
  ZXYPOSTITT1=.98
  IF(XPOSTITT1 /= 0.)THEN
    ZXPOSTITT1=XPOSTITT1
  ENDIF
  IF(XYPOSTITT1 /= 0.)THEN
    ZXYPOSTITT1=XYPOSTITT1
  ENDIF
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

IF(LMINUS .OR. LPLUS)THEN

  ZXPOSTITB3=.002
  ZXYPOSTITB3=.045
  IF(XPOSTITB3 /= 0.)THEN
    ZXPOSTITB3=XPOSTITB3
  ENDIF
  IF(XYPOSTITB3 /= 0.)THEN
    ZXYPOSTITB3=XYPOSTITB3
  ENDIF

  IF(.NOT.LTITDEFM .AND. CTITB3MEM /= 'DEFAULT' .AND. &
     CTITB3MEM /= 'default' .AND. CTITB3MEM /= 'DEFAUT' .AND. &
     CTITB3MEM /= 'defaut')THEN
     if(nverbia > 0)then
       print *,' imagev  CTITB3MEM ',CTITB3MEM(1:LEN_TRIM(CTITB3MEM))
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

!   print *,' **imagev CTITB3 AV RESOLV_TIT ',CTITB3
    CALL RESOLV_TIT('CTITB3',CTITB3)
!   print *,' **imagev CTITB3 AP RESOLV_TIT ',CTITB3
    IF(CTITB3 /= ' ')THEN
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTITB3,.009,0.,-1.)
    ENDIF
  ENDIF

ENDIF

1001 FORMAT('HORIZONTAL SECTION NIINF=',I4,' NISUP=',I4, &
            ' NJINF=',I4,' NJSUP=',I4)
CALL GSCLIP(1)
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
!IF(.NOT.LDIRWIND)THEN
! Conservation de la valeur du logique suivant pour la direction du vent
! pour beneficier des traits pleins en cas de superposition (Mai 99)
IF(LSUPER)THEN
  LARROVL=.TRUE.
ELSE
  LARROVL=.FALSE.
ENDIF
!
IF(LDIRWIND)THEN
! LDIRWIND=.FALSE.
ENDIF
!
!*    3.2   NCAR parameter reset
!
CALL CPSETI('CLS',16)
CALL CPRSET
!
!-------------------------------------------------------------------------
!
!*    4.    EXIT
!           ----
!
if(nverbia > 0)then
print *,' **imagev Sortie'
endif
RETURN
END SUBROUTINE  IMAGEV_FORDIACHRO

