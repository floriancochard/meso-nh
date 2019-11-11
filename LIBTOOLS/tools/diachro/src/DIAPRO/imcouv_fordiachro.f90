!     ######spl
      SUBROUTINE IMCOUV_FORDIACHRO(PU,PW,HLEGEND,HTEXT)
!     #################################################
!
!!****  *IMCOUV_FORDIACHRO* - Draws a vector arrow plot for a vertical cross-section
!!
!!    PURPOSE
!!    -------
!       Draws an arrow plot of a UW vector field re-colocated at the
!     mass gridpoint for a vertical cross-section
!
!!**  METHOD
!!    ------
!!     
!!     Assumption is made that wind components were re-colocated onto the mass
!!   gridpoint location prior to calling IMCOUV.
!!   The wind arrows are plotted using the VVECTR NCAR utility.
!!     
!!     Notice that a TRACE-provided VVUMXY routine is used within the NCAR
!!   vector VVECTR utility to map the wind vectors onto the stretched
!!   MESO-NH model space.  Wind vectors are given in m/s and scaled by VVUMXY
!!   to obtain arrow sizes in "NCAR fractional coordinate" (NCAR User Guide
!!   "Fundamentals", Appendix A, p345 section 1), notice this is different
!!   from what is required for Conpack... The final result is an automatic
!!   arrow scale selection on the plot.
!!   If a different procedure has to be followed VVUMXY should
!!   be updated accordingly. The parameters of the NCAR VVECTR utility can
!!   be printed online by typing "man vectors_params", these feature are not
!!   really documented elsewhere in NCAR user guide.
!!    
!!
!!    EXTERNAL
!!    --------
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
!!      GSLWSC    : sets line width                                 !
!!      VVRSET    : resets VVECTR parameters to default values     !
!!
!!
!!      VVUMXY    : TRACE provided FORTRAN-77 routine directly called
!!                  within the VVECTR NCAR utility to to map the wind
!!                  vectors onto the stretched MESO-NH model space.
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
!!         LHORIZ    : must be .FALSE. to perform vertical cross esctions
!!         LVERTI    : must be .TRUE. to perform vertical cross sections
!!         Module MODD_DIM1   : Contains dimensions
!!            NIMAX, NJMAX :  x, and y array dimensions
!!            NIINF, NISUP :  Lower and upper array bounds in x direction
!!            NJINF, NJSUP :  Lower bound and upper bound  in y direction
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist
!!                         (former NCAR common)
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
!!      Original       19/09/95
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF 
USE MODD_COORD
USE MODD_ALLOC_FORDIACHRO
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_GRID 
USE MODD_GRID1
USE MODD_PARAMETERS
USE MODD_NMGRID
USE MODD_FIELD1_CV2D
USE MODD_SUPER
USE MODD_TITLE
USE MODD_OUT
USE MODN_PARA
USE MODN_NCAR
USE MODD_LUNIT1
USE MODD_CVERT
USE MODD_CTL_AXES_AND_STYL
USE MODD_RESOLVCAR
USE MODD_TIT
USE MODD_PVT
USE MODD_MEMCV
USE MODD_DEFCV
USE MODE_GRIDPROJ
USE MODI_RESOLV_TIT
USE MODI_RESOLV_TITY
!
IMPLICIT NONE

INTERFACE
      SUBROUTINE INTERPOLW(PZZU, PZZW, PSTRU, PSTRW)
      REAL,DIMENSION(:,:) :: PZZU, PZZW, PSTRU, PSTRW
      END SUBROUTINE INTERPOLW
END INTERFACE

!
!*       0.0   TRACE interface with the "VVUMXY" routine of the NCAR package
!
! NOTICE:  The TRACE provided VVUMXY routine and the NCAR graphical utilities 
! ------   are NOT written in Fortran 90, but in Fortran 77.. This sub-section
!          of TRACE does not follow the Meso-NH usual rules: it has to be made 
!          using a COMMON stack with  static memory allocation of XZWORKZ and
!          XZZDS arrays.
!
!
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
COMMON/TEMV/XZWORKZ,XZZDS,NINX,NINY
#include "big.h"
REAL,DIMENSION(N2DVERTX,2500):: XZWORKZ
!REAL,DIMENSION(1000,400):: XZWORKZ
!REAL,DIMENSION(200,200) :: XZWORKZ
REAL,DIMENSION(N2DVERTX):: XZZDS
!REAL,DIMENSION(1000):: XZZDS
!REAL,DIMENSION(200) :: XZZDS
INTEGER :: NINX, NINY
LOGICAL :: LVERT, LHOR, LPT, LXABS
!
!*       0.1   NCAR work arrays
!
! See aforementioned notice. The dimensions of these arrays are
! subject to possible tuning, but have to be prescribed. Add
! extra size if necessary.
!
INTEGER,PARAMETER       :: JPRSCR=50000, JPISCR=50000

REAL,DIMENSION(JPRSCR):: ZRSCR
INTEGER,DIMENSION(JPISCR):: ISCR
!
!*       0.2   Dummy arguments and results
!
REAL,DIMENSION(:,:) :: PU, PW
CHARACTER(LEN=*) :: HTEXT       ! Plot heading containing field name
CHARACTER(LEN=*) :: HLEGEND
!
!*       0.3   Local variables
!
INTEGER :: JILOOP, JKLOOP, ID, IDD
INTEGER :: IKB, IKE, IKU
INTEGER,SAVE :: IKL, ISKIPX, ISKIPY, ISKIPXM
!!! Janvier 2001
INTEGER :: IUB1, IUB2, ITER, JTER, II,JJ, ITERM
INTEGER :: IJ, J, JA, JILOOPD, JILOOPF, I
INTEGER :: JLOOPI, JLOOPJ, III
INTEGER :: ILMAX
INTEGER :: INUM, IRESP, ILOOP, IDEB, IFIN
INTEGER :: IER,ICOL1
INTEGER,DIMENSION(18) :: ICOL
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: IE
!
!! Avec interpol en Z
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZSTRU, ZSTRW
REAL :: ZT,ZRA
INTEGER :: IZS, IU1, JU1,JU2, ISEUIL
!! Avec interpol en Z
!
REAL :: ZMI, ZMA, ZMIG, ZMAG, ZLATB, ZLONB
REAL :: ZRPK, ZBETA, ZLON0, ZVINT, ZVY, ZINTX, ZINTY
REAL,DIMENSION(18) :: ZCOL
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZX,ZLAT,ZLON,ZZY,ZYY
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZDIRU,ZDIRV,ZLA,ZLO
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZSTR1
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZLAB,ZLOB
CHARACTER(LEN=4) :: YTE
!!! Janvier 2001
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE  :: ZZU, ZZW
!!REAL,DIMENSION(NLMAX,SIZE(PU,2)) :: ZZU, ZZW
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZZDS
REAL,DIMENSION(N2DVERTX+20)       :: ZDS, ZWZ
!REAL,DIMENSION(1020)       :: ZDS, ZWZ
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZDS2, ZWZ2
REAL :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
REAL ::ZWLL, ZWRR, ZWBB, ZWTT
REAL :: ZY
REAL :: ZXB,ZYB
REAL :: ZDMX, ZVMX
REAL :: ZRAP
REAL :: ZXPOSTITT1, ZXYPOSTITT1
REAL :: ZXPOSTITT2, ZXYPOSTITT2
REAL :: ZXPOSTITT3, ZXYPOSTITT3
REAL :: ZXPOSTITB1, ZXYPOSTITB1
REAL :: ZXPOSTITB2, ZXYPOSTITB2
REAL :: ZXPOSTITB3, ZXYPOSTITB3
REAL,DIMENSION(1000) :: ZYYY
REAL :: ZU,ZW,ZM,ZUMN,ZWMN,ZMN,ZUMX,ZWMX,ZMX

CHARACTER(LEN=82) :: YCARCOU, YTEM
CHARACTER(LEN=80) :: YCAR
CHARACTER(LEN=40) :: YLBL
CHARACTER(LEN=10) :: YLBLMN,YLBLMX
CHARACTER(LEN=10) :: FORMAX, FORMAY

LOGICAL,SAVE :: GVSUPSCA
!
!*       0.4   External for NCAR use
!
! SFILL subroutine declared as external provides area control
! in some parts of the contour plot.
!
!EXTERNAL SFILL
EXTERNAL STUMXY
!
!-------------------------------------------------------------------------------
!
!*       1.    DISPLAY ENVIRONMENT SETUP AND ARROWS PLOTTING
!              ---------------------------------------------
!
!*       1.1   Array sizes calculation and default field value
!
!
IKU=NKMAX+2*JPVEXT
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
!print *,'size ZZU ZZW ZRSCR ISCR ',SIZE(ZZU,1),SIZE(ZZU,2),SIZE(ZZW,1), &
!SIZE(ZZW,2),SIZE(ZRSCR),SIZE(ISCR)
IF(ALLOCATED(ZZU))THEN
  DEALLOCATE(ZZU)
ENDIF
IF(ALLOCATED(ZZW))THEN
  DEALLOCATE(ZZW)
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(LSTREAM .AND. NISKIP /= 1)THEN
ILMAX=NLMAX/NISKIP
IF(ILMAX*NISKIP < NLMAX)ILMAX=ILMAX+1
ALLOCATE(ZZU(ILMAX,SIZE(PU,2)))
ALLOCATE(ZZW(ILMAX,SIZE(PU,2)))
DO JILOOP=1,ILMAX
  DO JKLOOP=1,IKU
    ZZU(JILOOP,JKLOOP)=XSPVAL
    ZZW(JILOOP,JKLOOP)=XSPVAL
  ENDDO
ENDDO
I=0
DO JILOOP=1,NLMAX,NISKIP
  I=I+1
  XZZDS(I)=XDS(JILOOP,NMGRID)
ENDDO
IF(I == ILMAX)THEN
ELSE
  I=I+1
  XZZDS(I)=XDS(NLMAX,NMGRID)
ENDIF
I=0
DO JILOOP=1,NLMAX,NISKIP
  I=I+1
    XZWORKZ(I,1:IKU)=XWORKZ(JILOOP,1:IKU,NMGRID)
ENDDO
IF(I == ILMAX)THEN
ELSE
  I=I+1
  XZWORKZ(I,1:IKU)=XWORKZ(NLMAX,1:IKU,NMGRID)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE

ALLOCATE(ZZU(NLMAX,SIZE(PU,2)))
ALLOCATE(ZZW(NLMAX,SIZE(PU,2)))
DO JILOOP=1,NLMAX
  DO JKLOOP=1,IKU
    ZZU(JILOOP,JKLOOP)=XSPVAL
    ZZW(JILOOP,JKLOOP)=XSPVAL
  ENDDO
ENDDO
!
!*       1.2  Collects X and Z values 
!
DO JILOOP=1,NLMAX
  XZZDS(JILOOP)=XDS(JILOOP,NMGRID)
ENDDO
DO JILOOP=1,NLMAX
  DO JKLOOP=1,IKU
    XZWORKZ(JILOOP,JKLOOP)=XWORKZ(JILOOP,JKLOOP,NMGRID)
  ENDDO
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*       1.3  Window definition and plot
!
!ZVL=.1
!ZVR=.9
!ZVB=.1
!ZVT=.9
IF(LVPTVUSER)THEN
  ZVL=XVPTVL
  ZVR=XVPTVR
  ZVB=XVPTVB
  ZVT=XVPTVT
ELSE
  ZVL=.1
  ZVR=.9
  ZVB=.1
  ZVT=.9
ENDIF
ZWL=XDS(1,NMGRID)
ZWR=XDS(NLMAX,NMGRID)
! 130101
IF((XHMIN==0..AND.XHMAX==0.).OR.(XHMAX<=XHMIN))THEN
  XHMIN=XWORKZ(1,IKB,NMGRID)
  XHMAX=XWORKZ(1,IKE,NMGRID)
END IF
ZWB=XHMIN
ZWT=XHMAX

LVERTI=.TRUE. ; LHORIZ=.FALSE.
LVERT=LVERTI
LHOR=LHORIZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(LSTREAM .AND. NISKIP /= 1)THEN
NINX=ILMAX
ELSE
NINX=NLMAX
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NINY=IKU
!print *,' **gsclip N1 0 '
CALL GSCLIP(0)

CALL GSLN(1)
CALL GSPLCI(1)
CALL GSTXCI(1)
IF(LSUPER)THEN
  NSUPER=NSUPER+1
! print *,' ***IMCOUV NSUPER*** ',NSUPER
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
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
  ELSE
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  END IF
ELSE
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
    NCOLUVG=NCOLUV1
ENDIF

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
!CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
CALL GASETI('LTY',1)
! Janvier 2001
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
IF(LPV)THEN
  IF(LFACTAXEY)THEN
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
	   ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
  ELSEIF(LAXEYUSER)THEN
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
	     XAXEYUSERD,XAXEYUSERF,IDD)
  ENDIF
!Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0)
  ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0)
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0)
  ELSE
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0)
  ENDIF
!Avril 2002
  CALL FRSTPT((ZWL+ZWR)/2,ZWB)
  CALL VECTOR((ZWL+ZWR)/2,ZWT)
ELSE
! Mars 2001
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
!Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0)
  ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0)
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,1,0,5,0.,0)
  ELSE
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,1,1,5,0.,0)
  ENDIF
!Avril 2002
ENDIF
! Mars 2001
IF(LFACTAXEX .OR. LFACTAXEY .OR. LAXEXUSER .OR. LAXEYUSER)THEN
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
ENDIF
! Mars 2001
!
!*       1.4  Collects wind values within the user postprocessing
!*            window with a sampling rate of NISKIP outside values 
!*            are kept to default
!
! Janvier 2001 On prevoit Vecteurs et direction Vent Horizontal en CV!
! Dans ce cas ZZW contient la composante V passee en argument

! Partie commune de LPRINT
IF(LPRINT .AND..NOT.LULMWM .AND..NOT.LULTWT)THEN
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=SIZE(PU,1)/5
  IF(ILOOP * 5 < SIZE(PU,1))ILOOP=ILOOP+1
  IF(LPV)ILOOP=1

  IF(.NOT.LPVT)THEN
    IF(LPV)THEN
      WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'',''   (NPROFILE,1-IKU)'')')CGROUP,&
&     CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ELSE
      WRITE(INUM,'(''CV  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'',''   (1-NLMAX,1-IKU)'')')CGROUP,&
&     CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ENDIF
  ELSE
    WRITE(INUM,'(''CV  '',''G:'',A16,'' P:'',A25)')CGROUP,&
&   CTITRE(NLOOPP)(1:25)
  ENDIF

  IF(LMINUS .OR. LPLUS)THEN
    WRITE(INUM,'(A70)')CTITB3
  ELSE
    WRITE(INUM,'(A40)')CTITGAL
  ENDIF

  IF(.NOT.LPV)THEN
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2)THEN
        WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,SIZE(PU,2),ILOOP
      ELSE IF(LDEFCV2LL)THEN
        WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,SIZE(PU,2),ILOOP
      ELSE IF(LDEFCV2IND)THEN
        WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,SIZE(PU,2),ILOOP
      ENDIF
    ELSE
      IF(XIDEBCOU /= -999.)THEN
        WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
      ELSE
        WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
      ENDIF
    ENDIF
  ELSE
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2)THEN
        WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,SIZE(PU,2),ILOOP
	WRITE(INUM,'(''nprofile='',I4)')NPROFILE
      ELSE IF(LDEFCV2LL)THEN
        WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,SIZE(PU,2),ILOOP
	WRITE(INUM,'(''nprofile='',I4)')NPROFILE
      ELSE IF(LDEFCV2IND)THEN
        WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
    &'' iku'',i4,'' iter'',i3)')&
       &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,SIZE(PU,2),ILOOP
	WRITE(INUM,'(''nprofile='',I4)')NPROFILE
      ENDIF
    ELSE
      IF(XIDEBCOU /= -999.)THEN
        WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
	WRITE(INUM,'(''nprofile='',I4)')NPROFILE
      ELSE
        WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
    &'' iku'',i4,''    iter'',i3)')&
       &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
	WRITE(INUM,'(''nprofile='',I4)')NPROFILE
      ENDIF
    ENDIF
  ENDIF
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**IMCOUV XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
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
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IF(.NOT.LDIRWIND)THEN
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DO JLOOPI=1,ILOOP
    IF(JLOOPI == 1)THEN
      IDEB=1; IFIN=5
    ELSE
      IDEB=IFIN+1; IFIN=IFIN+5
    ENDIF
    IF(JLOOPI == ILOOP)THEN
      IFIN=SIZE(PU,1)
    ENDIF
    IF(LPV)THEN
      IDEB=NPROFILE;IFIN=NPROFILE
    ENDIF
    
    WRITE(INUM,'(1X,25(1H*),'' U Component '',41(1H*))')
!   WRITE(INUM,'(1X,79(1H*))')
    WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
    WRITE(INUM,'(''.'',79(1H*))')
    DO JLOOPJ=SIZE(PU,2),1,-1
      WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(PU(II,JLOOPJ),II=IDEB,IFIN)
!     WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(PU(II,JLOOPJ),II=IDEB,IFIN)
    ENDDO
    WRITE(INUM,'(1X,79(1H*))')
  ENDDO
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DO JLOOPI=1,ILOOP
    IF(JLOOPI == 1)THEN
      IDEB=1; IFIN=5
    ELSE
      IDEB=IFIN+1; IFIN=IFIN+5
    ENDIF
    IF(JLOOPI == ILOOP)THEN
      IFIN=SIZE(PU,1)
    ENDIF
    IF(LPV)THEN
      IDEB=NPROFILE;IFIN=NPROFILE
    ENDIF
    
    WRITE(INUM,'(1X,25(1H*),'' V Component '',41(1H*))')
!   WRITE(INUM,'(1X,79(1H*))')
    WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
    WRITE(INUM,'(''.'',79(1H*))')
    DO JLOOPJ=SIZE(PW,2),1,-1
      WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(PW(II,JLOOPJ),II=IDEB,IFIN)
!     WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(PW(II,JLOOPJ),II=IDEB,IFIN)
    ENDDO
    WRITE(INUM,'(1X,79(1H*))')
  ENDDO
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(LDIRWIND .OR. LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT)THEN
  ZRAP=1
  ISKIPX=NISKIPVX
  ISKIPY=NISKIPVY
  IF(LPV)THEN
    ISKIPX=1
  ENDIF
ELSE
! Dilatation de la composante W par ZRAP
  IF(LDILW)THEN
    ZRAP=((ZWR-ZWL)/(ZVR-ZVL))/((ZWT-ZWB)/(ZVT-ZVB))
  ELSE
    ZRAP=1
  ENDIF
  !ISKIPX=NISKIPVX
  !ISKIPY=NISKIPVY
  ISKIPX=NISKIP
  ISKIPY=NISKIP
  IF(LSTREAM)ISKIPY=1
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nverbia <0)then
  print *,' ***IMCOUV ZWR,ZWL,ZVR,ZVL,ZWT,ZWB,ZVT,ZVB,ZRAP '
  print *,ZWR,ZWL,ZVR,ZVL,ZWT,ZWB,ZVT,ZVB,ZRAP
endif
! Determination egalement du min et max reels
ZUMN=999.;ZUMX=-999.;ZWMN=999.;ZWMX=-999.;ZMN=999.;ZMX=-999.
!print *,' IMCOUV NLMAX ',NLMAX,NISKIP

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Janvier 2001
DO JKLOOP=IKB,IKE,ISKIPY

  IF(LPV)THEN
    JILOOPD=NPROFILE
    JILOOPF=NPROFILE
  ELSE
    JILOOPD=1
    JILOOPF=NLMAX
  ENDIF
  I=0

  DO JILOOP=JILOOPD,JILOOPF,ISKIPX

    IF(.NOT.LSTREAM)THEN
    IF(XZWORKZ(JILOOP,JKLOOP) > XHMAX .OR. XZWORKZ(JILOOP,JKLOOP) < XHMIN)CYCLE
    ENDIF
    ZU=PU(JILOOP,JKLOOP)
    ZW=PW(JILOOP,JKLOOP)
    IF(ZU /= XSPVAL .AND. ZW /= XSPVAL)THEN   
      ZM=SQRT(ZU*ZU+ZW*ZW)
!
      IF(ZM.LT.ZMN)THEN
        ZMN=ZM;ZUMN=ZU;ZWMN=ZW
      ENDIF
      IF(ZM.GT.ZMX)THEN
        ZMX=ZM;ZUMX=ZU;ZWMX=ZW
      ENDIF
!       
 
      IF(LSTREAM .AND. NISKIP /= 1)THEN
          I=I+1
          if(nverbia <0)then
          print *,' **JILOOP,NISKIP,I,ILMAX ',JILOOP,NISKIP,I,ILMAX
          endif
          ZZU(I,JKLOOP)=PU(JILOOP,JKLOOP)
          ZZW(I,JKLOOP)=PW(JILOOP,JKLOOP)
          IF((JILOOP == JILOOPF .OR. JILOOP > JILOOPF-ISKIPX)&
               .AND. I /= ILMAX)THEN
            I=I+1
          if(nverbia <0)then
          print *,' **JILOOP,JILOOPF,NISKIP,I,ILMAX ',JILOOP,JILOOPF,&
          NISKIP,I,ILMAX
          endif
            ZZU(I,JKLOOP)=PU(JILOOP,JKLOOP)
            ZZW(I,JKLOOP)=PW(JILOOP,JKLOOP)
            EXIT
          ENDIF

      ELSE

        ZZU(JILOOP,JKLOOP)=PU(JILOOP,JKLOOP)
        if(nverbia <0)then
        if(JKloop == IKB)THEN
          print *,' ***IMCOUV PW ',PW(JILOOP,JKLOOP)
        ENDIF
        ENDIF

        IF(.NOT.LSTREAM)THEN
          ZZW(JILOOP,JKLOOP)=PW(JILOOP,JKLOOP)*ZRAP
        ELSE
          ZZW(JILOOP,JKLOOP)=PW(JILOOP,JKLOOP)
        ENDIF

      ENDIF
      if(nverbia <0)then
        if(JKloop == IKB)THEN
          print *,' ***IMCOUV ZW ',ZZW(JILOOP,JKLOOP)
        ENDIF
      ENDIF

    ENDIF

  ENDDO

ENDDO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!print *,' JILOOP,JKLOOP ',JILOOP,JKLOOP
!print *,' ZRAP,ZMN,ZMX,ZUMN,ZWMN,ZUMX,ZWMX ',ZRAP,ZMN,ZMX,ZUMN,ZWMN,ZUMX,ZWMX
!
!        1.41 Topography
!

!print *,' **gsclip N2 1 '
CALL GSCLIP(1)

! Janvier 2001
IF(.NOT. LPV)THEN
! Janvier 2001
IF(.NOT.LSUPER .OR. NSUPER == 1)THEN
  ZDS(1)=XDS(1,NMGRID)
  ZWZ(1)=XHMIN
  DO JILOOP=2,NLMAX+1
    ZDS(JILOOP)=XDS(JILOOP-1,NMGRID)
    ZWZ(JILOOP)=XWZ(JILOOP-1,NMGRID)
  ENDDO
  ZDS(NLMAX+2)=ZDS(NLMAX+1)
  ZWZ(NLMAX+2)=XHMIN
  IF(ALLOCATED(ZDS2))THEN
    DEALLOCATE(ZDS2)
  ENDIF
  IF(ALLOCATED(ZWZ2))THEN
    DEALLOCATE(ZWZ2)
  ENDIF
  ALLOCATE(ZWZ2(NLMAX+2))
  ALLOCATE(ZDS2(NLMAX+2))
  ZDS2=ZDS(1:NLMAX+2)
  ZWZ2=ZWZ(1:NLMAX+2)

  CALL CURVE(ZDS2,ZWZ2,NLMAX+2)
  CALL SFSETR('SP',.008)
  CALL SFSETR('AN',45.)
  CALL SFSETI('DO',0)
  CALL SFWRLD(ZDS2,ZWZ2,NLMAX+2,ZRSCR,JPRSCR,ISCR,JPISCR)
ENDIF
! Janvier 2001
ENDIF
! Janvier 2001

!print *,' **gsclip N3 0 '
CALL GSCLIP(0)

!
!       If required draw a model-level background
!
!
IF(LXZ)THEN
  DO JKLOOP=IKU,1,-1
    IF(ZZU(1,JKLOOP) /= XSPVAL)EXIT
  ENDDO
  IKL=JKLOOP
  CALL GSLN(3)
  DO JKLOOP=1,IKL
    ZYYY(1:NLMAX)=XZWORKZ(1:NLMAX,JKLOOP)
    CALL GPL(NLMAX,XZZDS,ZYYY)
  ENDDO
  CALL GSLN(1)
ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Janvier 2001 + LDIRWIND
IF(LDIRWIND)THEN
  print *,' imcouv LDIRWIND ',LDIRWIND
  IUB1=SIZE(PU,1)
  ITERM=IUB1/ISKIPX+1
  IF(1+(ITERM-1)*ISKIPX > IUB1)ITERM=ITERM-1
  ITER=IUB1
  ISKIPXM=ISKIPX
  ISKIPX=1
  IUB2=SIZE(PU,2)
! 130101
!!! Essai de conservation de 1 a IKU en Y (Pour LPRINT) mais
!!! de 1 a NLMAX par pas de NISKIPX en X
!!!  JTER=(IUB2-IKB)/ISKIPY+1
!!!  IF(IKB+(JTER-1)*ISKIPY > IUB2)JTER=JTER-1
  JTER=IUB2
!!!
  ALLOCATE(ZX(ITER,1),ZZY(ITER,JTER),ZYY(ITER,1),ZLAT(ITER,1),ZLON(ITER,1))
  ALLOCATE(ZLA(ITER,JTER),ZLO(ITER,JTER),ZDIRU(ITER,JTER),ZDIRV(ITER,JTER))
  ALLOCATE(ZZDS(ITER))
! 130101
! print *,' IIIIIMCOUV IUB1, ISKIPX, ITER, IUB2, ISKIPY, JTER,LPV ',IUB1,ISKIPX,ITER,IUB2,ISKIPY,JTER,LPV
!!!
!!!  ZDIRU=PU(1:IUB1:ISKIPX,IKB:IUB2:ISKIPY)
!!!  ZDIRV=PW(1:IUB1:ISKIPX,IKB:IUB2:ISKIPY)
  ZDIRU=XSPVAL
  ZDIRV=XSPVAL
  ZDIRU=PU(1:IUB1:ISKIPX,1:IUB2:1)
  ZDIRV=PW(1:IUB1:ISKIPX,1:IUB2:1)
!!!

  ZZDS=XDS(1:IUB1:ISKIPX,1)

! print *,' IIIIIMCOUV XDSX(1:IUB1) ',XDSX(1:IUB1,1)
! print *,' IIIIIMCOUV ZX(:,1) ',ZX(:,1)

! 130101
  JJ=0

!!!
!!!  DO JKLOOP=IKB,IUB2,ISKIPY
  DO JKLOOP=1,IUB2
!!!
    JJ=JJ+1
    II=0
    DO JILOOP=1,IUB1,ISKIPX
      II=II+1
      ZZY(II,JJ)=XZWORKZ(JILOOP,JKLOOP)
    ENDDO
  ENDDO

! 130101
! print *,' IIIIMCOUV IUB1,ISKIPX,IKB,IUB2,ISKIPY ',IUB1,ISKIPX,IKB,IUB2
! print *,' IIIIMCOUV XZWORKZ(1:NLMAX,IKB) ',XZWORKZ(1:NLMAX,IKB)
! print *,' IIIIMCOUV ZZY(:,1) ',ZZY(:,1)
! print *,' IIIIMCOUV XZWORKZ(1:NLMAX,IKB+1) ',XZWORKZ(1:NLMAX,IKB+1)
! print *,' IIIIMCOUV ZZY(:,2) ',ZZY(:,2)

! 130101
  ZX(:,1)=XDSX(1:IUB1:ISKIPX,1)
  ZYY(:,1)=XDSY(1:IUB1:ISKIPX,1)

  DO JKLOOP=1,JTER
    CALL SM_LATLON_A(XLATORI,XLONORI,ZX,ZYY,ZLAT,ZLON)
    ZLA(:,JKLOOP)=ZLAT(:,1)
    ZLO(:,JKLOOP)=ZLON(:,1)
  ENDDO

  where(zdiru /= xspval .AND. zdirv /= xspval)
    ZDIRU=ATAN2(ZDIRV,ZDIRU)*180./ACOS(-1.)
  endwhere

  if(nverbia > 0)then
    print *,' ZDIRU 1,1 ITER/2,1 1,JTER/2 ITER/2,JTER/2 ITER,JTER '
    print *,ZDIRU(1,1),  ZDIRU(ITER/2,1), ZDIRU(1,JTER/2), ZDIRU(ITER/2,JTER/2), &
    ZDIRU(ITER,JTER)
  endif

  ZRPK=XRPK
  ZBETA=XBETA
  ZLON0=XLON0
  where(zdiru /= xspval .AND. zdirv /= xspval)
    ZDIRU=ZDIRU - (ZRPK*(ZLO-ZLON0)-ZBETA) + 90.
  endwhere
  WHERE(ZDIRU < 0.)ZDIRU=ZDIRU+360.
  WHERE(ZDIRU > 360. .AND. ZDIRU /= XSPVAL)ZDIRU=ZDIRU-360.

  if(nverbia > 0)then
    print *,' ZDIRU 1,1 ITER/2,1 1,JTER/2 ITER/2,JTER/2 ITER,JTER '
    print *,ZDIRU(1,1),  ZDIRU(ITER/2,1), ZDIRU(1,JTER/2), ZDIRU(ITER/2,JTER/2), &
    ZDIRU(ITER,JTER)
  endif

  where(zdiru /= xspval .AND. zdirv /= xspval)
    ZDIRV=360.-ZDIRU
  elsewhere
    ZDIRV=XSPVAL
  endwhere

  if(nverbia > 0)then
    print *,' ZDIRV 1,1 ITER/2,1 1,JTER/2 ITER/2,JTER/2 ITER,JTER '
    print *,ZDIRV(1,1),  ZDIRV(ITER/2,1), ZDIRV(1,JTER/2), ZDIRV(ITER/2,JTER/2), &
    ZDIRV(ITER,JTER)
  endif

  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)

        if(nverbia <0)then
  print *,' ** imcouv ap getset ZWL,ZWR,XDS(1,1),XDS(NLMAX,1),ZX(1,1),ZX(ITER,1) ',ZWL,ZWR,XDS(1,1),XDS(NLMAX,1),ZX(1,1),ZX(ITER,1)
        endif
  IF(ITERM > 6)THEN
    CALL GSCLIP(1)
  ELSE
    CALL GSCLIP(0)
  ENDIF

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

  IF(LPV)THEN
    JILOOPD=NPROFILE
    JILOOPF=NPROFILE
  ELSE
    JILOOPD=1
    JILOOPF=ITER
  ENDIF

!!!
!!!  DO JKLOOP=1,JTER
  DO JKLOOP=IKB,JTER,ISKIPY
!!!

    DO JILOOP=JILOOPD,JILOOPF,ISKIPXM
!   DO JILOOP=JILOOPD,JILOOPF

      IF(ZDIRV(JILOOP,JKLOOP) == XSPVAL)THEN
!       print *,J,' CYCLE  ',ZDIRV(JILOOP,JKLOOP),ZCOL(J),ZCOL(J-1)
	CYCLE
      ENDIF

      DO J=2,IJ
!       print *,J,' ',ZDIRV(JILOOP,JKLOOP),ZCOL(J),ZCOL(J-1)
        
	IF(ZDIRV(JILOOP,JKLOOP) == 0. .OR. ZDIRV(JILOOP,JKLOOP) == 360.)THEN
	  CALL GSPMCI(ICOL(1))
!         print *,' ZDIRV(JILOOP,JKLOOP) J+2 ',ZDIRV(JILOOP,JKLOOP),ICOL(1)
	  EXIT
	ELSE IF(ZDIRV(JILOOP,JKLOOP) < ZCOL(J).AND. &
		ZDIRV(JILOOP,JKLOOP) >= ZCOL(J-1))THEN
	  CALL GSPMCI(ICOL(J-1))
!         print *,' ZDIRV(JILOOP,JKLOOP) J+1 ',ZDIRV(JILOOP,JKLOOP),ICOL(J)
	  EXIT
	ENDIF
      ENDDO

      CALL GSMK(2)

!!! Janvier 2001
      IF(LPV)THEN
        ZINTX=(ZWL+ZWR)/2
      ELSE
        ZINTX=ZZDS(JILOOP)
      ENDIF

      ZINTY=ZZY(JILOOP,JKLOOP)
      IF(ZINTY > XHMAX .OR. ZINTY <XHMIN)THEN
	CYCLE
      ENDIF

      CALL GPM(1,ZINTX,ZINTY)
      CALL GSMK(3)
      CALL GPM(1,ZINTX,ZINTY)
      CALL GSMK(5)
      CALL GPM(1,ZINTX,ZINTY)

    ENDDO
    CALL SFLUSH

  ENDDO

!print *,' **gsclip N4 0 '
  CALL GSCLIP(0)

! Legende couleurs
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
! print *,'ZVR,ZVY,YTE ',ZVR,ZVY,YTE
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

  IF(LPRINT)THEN
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DO JLOOPI=1,ILOOP
      IF(JLOOPI == 1)THEN
        IDEB=1; IFIN=5
      ELSE
        IDEB=IFIN+1; IFIN=IFIN+5
      ENDIF
      IF(JLOOPI == ILOOP)THEN
        IFIN=SIZE(PU,1)
      ENDIF
      IF(LPV)THEN
        IDEB=NPROFILE;IFIN=NPROFILE
      ENDIF
      
      WRITE(INUM,'(1X,79(1H*))')
      WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',79(1H*))')
      DO JLOOPJ=SIZE(ZDIRV,2),1,-1
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(ZDIRV(II,JLOOPJ),II=IDEB,IFIN)
  !     WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(ZDIRV(II,JLOOPJ),II=IDEB,IFIN)
      ENDDO
      WRITE(INUM,'(1X,79(1H*))')
    ENDDO
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF

!print *,' **gsclip N5 0 '
  CALL GSCLIP(0)
  DEALLOCATE(ZX,ZZY,ZYY,ZLAT,ZLON,ZLA,ZLO,ZDIRU,ZDIRV,ZZDS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELSE
!!!! Janvier 2001 + LDIRWIND

IF (.NOT.LSTREAM)THEN
!
!
!*       1.5  Routine VVUMXY of provided by TRACE to locate and scale wind
!*            arrows on the display
!
CALL VVSETI('MAP',4)
CALL VVSETI('SET',0)
CALL VVSETR('VPL',ZVL)    
CALL VVSETR('VPR',ZVR)
CALL VVSETR('VPB',ZVB)
CALL VVSETR('VPT',ZVT)
CALL VVSETR('WDL',ZWL)
CALL VVSETR('WDR',ZWR)
CALL VVSETR('WDB',ZWB)
CALL VVSETR('WDT',ZWT)

! Sortie de statistiques si LVST=T
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

! Elimination de la legende des fleches si LEGVECT=F
IF(.NOT.LEGVECT)THEN
  CALL VVSETC('MXT',' ')
  CALL VVSETC('MNT',' ')
ENDIF

! Janv 2001 Si XVHC <0 (Scale) conservation tout de meme des valeurs > xvhc
! Intervention ds vvectr rajoute a frame (12/1/2001 je n'ai pas fait gd chose
! Besoin peut-etre de reintervenir)
IF(XVHC >= 0.)THEN
  GVSUPSCA=LVSUPSCA
  LVSUPSCA=.FALSE.
ENDIF
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
CALL VVSETI('MNP',-4)
CALL VVSETI('MXP',-4)
CALL VVSETR('MNX',.75)
!CALL VVSETR('MNX',-ZVL)

IF(ZVB <= .15)THEN
  ZY=(-.08)/(ZVT-ZVB)
ELSE
  ZY=(-.13)/(ZVT-ZVB)
ENDIF

CALL VVSETR('MNY',ZY)

IF(ZVR-ZVL >= .78)THEN
  CALL VVSETR('MXX',.75+.16)
ELSE
  CALL VVSETR('MXX',.75+.27)
ENDIF

CALL VVSETR('MXY',ZY)
CALL VVSETR('MXS',.008*.9/(ZVR-ZVL))
CALL VVSETR('MNS',.008*.9/(ZVR-ZVL))
!
!*     1.7    Draws the arrows
!
IF(XLWV > 0.)THEN
  CALL VVSETR('LWD',XLWV)
ELSE
  CALL VVSETR('LWD',XLWVDEF)
ENDIF

!print *,' **gsclip N6 0 '
CALL GSCLIP(0)                                     ! Clipping off
CALL VVSETI('VPO',1)
CALL VVINIT(ZZU,NLMAX,ZZW,NLMAX,0.,0,NLMAX,IKU,0.,0) ! Initializes VVECTR
CALL VVECTR(ZZU,ZZW,0.,0,0,0.)                     ! Draws arrows
!print *,' **gsclip N7 1 '
CALL GSCLIP(1)                                     ! Clipping back on
!
CALL VVRSET

!!!! Janvier 2001 + LDIRWIND
IF(XVHC >= 0.)THEN
  LVSUPSCA=GVSUPSCA
ENDIF

!!!!!!!!!!!!!!!!!!!!STREAM
ELSE

  NSGD=1
  IF(LINTERPOLSTR)THEN
!! Avec interpol en Z
    IF(NISKIP /= 1)THEN
      ALLOCATE(ZSTR1(4*ILMAX*NZSTR))
      ALLOCATE(ZSTRU(ILMAX,NZSTR))
      ALLOCATE(ZSTRW(ILMAX,NZSTR))
    ELSE
      ALLOCATE(ZSTR1(4*NLMAX*NZSTR))
      ALLOCATE(ZSTRU(NLMAX,NZSTR))
      ALLOCATE(ZSTRW(NLMAX,NZSTR))
    ENDIF
    ZSTR1=0.
if(nverbia >0)then
  print *,' Appel interpolw '
  endif

  CALL INTERPOLW(ZZU,ZZW,ZSTRU,ZSTRW)
  if(nverbia >0)then
  print *,' Apres Appel interpolw '
  endif
!! Avec interpol en Z

  ELSE

    IF(NISKIP /= 1)THEN
      ALLOCATE(ZSTR1(4*ILMAX*IKU))
    ELSE
      ALLOCATE(ZSTR1(4*NLMAX*IKU))
    ENDIF
    ZSTR1=0.
!! Recherche d'un seuil pour choisir la frequence de depart 
!! d'1 streamline
   IU1=MAX(JPHEXT+1,2); JU1=MAX(JPVEXT+1,2)
   JU2=(NINY-JPVEXT)  
   ZT=XZWORKZ(IU1,JU2)-XZWORKZ(IU1,JU1)
   DO J=1+JPVEXT+1,NINY-JPVEXT

    ZRA= (XZWORKZ(IU1,J)-XZWORKZ(IU1,J-1))/(ZT)
    IF(ZRA >= .05)THEN
      NSGD=2
      NSEUIL=J
      ISEUIL=NSEUIL
  if(nverbia <0)then
      print *, '** imcouv RAP NSEUIL ',ZRA,NSEUIL
  endif
      EXIT
    ENDIF
   
   ENDDO
   
  ENDIF
  CALL STSETI('MAP',4)
! CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  if(nverbia > 0)then
    print *,' **imcouv ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT,NINX,NINY
    print *,' **imcouv ap getset ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT,NINX,NINY
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
    print *,' **imcouv ZWL,ZWR,ZWB,ZWT ',ZWL,ZWR,ZWB,ZWT
  endif

  CALL STSETR('GBS',0.)
! CALL STSETI('AGD',2)
! pour suppression de la fleche de depart d'1 streamline
! CALL STSETR('AMD',.005)
  CALL STSETI('AGD',NARSTR)
! defaut AGD=4
  CALL STSETI('SGD',2)
  CALL STSETI('CPM',0)
  CALL STSETI('TRT',0)
  CALL STSETI('TRP',0)
  CALL STSETI('CKX',1)
! CALL STSETI('CKP',30)
  CALL STSETR('ARL',XARLSTR)
! defaut ARL=.009
  CALL STSETR('DFM',.02)
  CALL STSETR('CDS',2.)
  CALL STSETR('LWD',XLWSTR)
! defaut LWD=1
  IF(LVERT)THEN
  CALL STSETR('SSP',XSSP)
! defaut SSP=.004
  ELSE
  CALL STSETR('SSP',XSSP)
  ENDIF
  CALL STSETI('MSK',0)
  CALL STSETI('SVF',3)
  CALL STSETR('USV',XSPVAL)
  CALL STSETR('VSV',XSPVAL)

  IF(LINTERPOLSTR)THEN
!! Avec interpol en Z
    ZSTRW=ZSTRW*ZRAP
    IF(NISKIP /= 1)THEN
      IZS=4*NZSTR*ILMAX
      CALL STINIT(ZSTRU,ILMAX,ZSTRW,ILMAX,0.,0,ILMAX,NZSTR,ZSTR1,IZS) !
    ELSE
      IZS=4*NZSTR*NLMAX
      CALL STINIT(ZSTRU,NLMAX,ZSTRW,NLMAX,0.,0,NLMAX,NZSTR,ZSTR1,IZS) !
    ENDIF
    CALL GQPLCI(IER,ICOL1)
    CALL GSPLCI(NCOLUVG)
    CALL STREAM(ZSTRU,ZSTRW,0.,0,0.,ZSTR1)                     ! Draws arrows
!   print *,' **incouv AP STREAM '

  ELSE

!print *,' **gsclip N8 0 '
    CALL GSCLIP(0)                                     ! NO Clipping 
   IF(NISKIP /= 1)THEN
     IZS=4*IKU*ILMAX
   ELSE
     IZS=4*IKU*NLMAX
   ENDIF
   ZZW=ZZW*ZRAP
   IF(NSGD == 2)THEN
     IF(NISKIP /= 1)THEN
       CALL STINIT(ZZU,ILMAX,ZZW,ILMAX,0.,0,ILMAX,IKU,ZSTR1,IZS) ! Initializes VVECTR
     ELSE
       CALL STINIT(ZZU,NLMAX,ZZW,NLMAX,0.,0,NLMAX,IKU,ZSTR1,IZS) ! Initializes VVECTR
     ENDIF
     CALL GQPLCI(IER,ICOL1)
     CALL GSPLCI(NCOLUVG)
     CALL STREAM(ZZU,ZZW,0.,0,0.,ZSTR1)                     ! Draws arrows
! CALL STREAM(ZZU,ZZW,0.,0,STUMXY,ZSTR1)                     ! Draws arrows
   ELSE
     IF(NISKIP /= 1)THEN
       CALL STINIT(ZZU,ILMAX,ZZW,ILMAX,0.,0,ILMAX,IKU,ZSTR1,IZS) ! Initializes VVECTR
     ELSE
       CALL STINIT(ZZU,NLMAX,ZZW,NLMAX,0.,0,NLMAX,IKU,ZSTR1,IZS) ! Initializes VVECTR
     ENDIF
     CALL GQPLCI(IER,ICOL1)
     CALL GSPLCI(NCOLUVG)
     CALL STREAM(ZZU,ZZW,0.,0,0.,ZSTR1)                     ! Draws arrows
! CALL STREAM(ZZU,ZZW,0.,0,STUMXY,ZSTR1)                     ! Draws arrows
  ENDIF
  IF(NSGD == 2)THEN
! ZSTR1=0.
  CALL STSETI('SGD',1)
  NSEUIL=ISEUIL
!   ZSTR1=0.
  IF(NISKIP /= 1)THEN
    CALL STINIT(ZZU,ILMAX,ZZW,ILMAX,0.,0,ILMAX,IKU,ZSTR1,IZS) ! Initializes VVECTR
  ELSE
    CALL STINIT(ZZU,NLMAX,ZZW,NLMAX,0.,0,NLMAX,IKU,ZSTR1,IZS) ! Initializes VVECTR
  ENDIF
  CALL GQPLCI(IER,ICOL1)
  CALL GSPLCI(NCOLUVG)
  CALL STREAM(ZZU,ZZW,0.,0,0.,ZSTR1)                     ! Draws arrows
! CALL STREAM(ZZU,ZZW,0.,0,STUMXY,ZSTR1)                     ! Draws arrows
  ENDIF

  ENDIF

  DEALLOCATE(ZSTR1)
  IF(LINTERPOLSTR)THEN
  DEALLOCATE(ZSTRU)
  DEALLOCATE(ZSTRW)
  ENDIF
!print *,' **gsclip N9 1 '
CALL GSCLIP(1)                                     ! Clipping back on
  CALL STRSET
  CALL GSPLCI(ICOL1)
ENDIF
!!!!!!!!!!!!!!!!!!!!STREAM

ENDIF
!!!! Janvier 2001 + LDIRWIND
IF(LPRINTXY .AND..NOT.LULMWM .AND..NOT.LULTWT)THEN
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=SIZE(PU,1)/5
  IF(ILOOP * 5 < SIZE(PU,1))ILOOP=ILOOP+1
  IF(LPV)ILOOP=1

  IF(.NOT. LPVT)THEN

    IF(.NOT.LPV)THEN
      WRITE(INUM,'(''CV XZ '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'','' (1-NLMAX,1-IKU)'')')CGROUP, &
&     CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ELSE
      WRITE(INUM,'(''PV XZ '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'','' (NPROFILE,1-IKU)'')')CGROUP, &
&     CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
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
      ALLOCATE(ZLAB(NLMAX),ZLOB(NLMAX))
      DO J=1,NLMAX
	ZXB=XDSX(J,NMGRID)
	ZYB=XDSY(J,NMGRID)
	CALL SM_LATLON_S(XLATORI,XLONORI,ZXB,ZYB,ZLATB,ZLONB)
	ZLAB(J)=ZLATB
	ZLOB(J)=ZLONB
      ENDDO
      IF(LDEFCV2LL)THEN
	ZLAB(1)=XIDEBCVLL
	ZLOB(1)=XJDEBCVLL
      ENDIF
      if(nverbia > 0)then
!     print *,' ZLA'
!     print *,ZLA
!     print *,' ZLO'
!     print *,ZLO
      endif
    ENDIF

    IF(.NOT.LPV)THEN
      IF(LDEFCV2CC)THEN
        IF(LDEFCV2)THEN
          WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
      &'' iku'',i4,'' iter'',i3)')&
         &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,SIZE(PU,2),ILOOP
        ELSE IF(LDEFCV2LL)THEN
          WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
      &'' iku'',i4,'' iter'',i3)')&
         &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,SIZE(PU,2),ILOOP
        ELSE IF(LDEFCV2IND)THEN
          WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
      &'' iku'',i4,'' iter'',i3)')&
         &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,SIZE(PU,2),ILOOP
        ENDIF
      ELSE
        IF(XIDEBCOU /= -999.)THEN
          WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
      &'' iku'',i4,''    iter'',i3)')&
         &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
        ELSE
      WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4, &
  &  '' iku'',i4,''    iter'',i3)') &
    & NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
        ENDIF
      ENDIF
    ELSE
      IF(LDEFCV2CC)THEN
        IF(LDEFCV2)THEN
          WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
      &'' iku'',i4,'' iter'',i3)')&
         &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,SIZE(PU,2),ILOOP
	  WRITE(INUM,'(''nprofile='',I4)')NPROFILE
        ELSE IF(LDEFCV2LL)THEN
          WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
      &'' iku'',i4,'' iter'',i3)')&
         &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,SIZE(PU,2),ILOOP
	  WRITE(INUM,'(''nprofile='',I4)')NPROFILE
        ELSE IF(LDEFCV2IND)THEN
          WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
      &'' iku'',i4,'' iter'',i3)')&
         &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,SIZE(PU,2),ILOOP
	  WRITE(INUM,'(''nprofile='',I4)')NPROFILE
        ENDIF
      ELSE
        IF(XIDEBCOU /= -999.)THEN
          WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
      &'' iku'',i4,''    iter'',i3)')&
         &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
	  WRITE(INUM,'(''nprofile='',I4)')NPROFILE
        ELSE
      WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4, &
  &  '' iku'',i4,''    iter'',i3)') &
    & NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,SIZE(PU,2),ILOOP
	  WRITE(INUM,'(''nprofile='',I4)')NPROFILE
        ENDIF
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
          WRITE(INUM,'(''   1 '',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
          XWZ(JLOOPI,NMGRID),ZLAB(JLOOPI),ZLOB(JLOOPI)
        ELSE IF(JLOOPI == NLMAX)THEN
          WRITE(INUM,'(''NLMAX'',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
          XWZ(JLOOPI,NMGRID),ZLAB(JLOOPI),ZLOB(JLOOPI)
        ELSE
          WRITE(INUM,'(''     '',I5,2(1X,E15.8),2(2X,F10.5))')JLOOPI,XDS(JLOOPI,NMGRID), &
          XWZ(JLOOPI,NMGRID),ZLAB(JLOOPI),ZLOB(JLOOPI)
        ENDIF
      ENDDO
      WRITE(INUM,'(1X,66(1H*))')
      DEALLOCATE(ZLAB,ZLOB)
    ENDIF
  
    DO JLOOPI=1,ILOOP
      IF(JLOOPI == 1)THEN
        IDEB=1; IFIN=5
      ELSE
        IDEB=IFIN+1; IFIN=IFIN+5
      ENDIF
      IF(JLOOPI == ILOOP)THEN
        IFIN=SIZE(PU,1)
      ENDIF
      IF(LPV)THEN
	IDEB=NPROFILE; IFIN=NPROFILE
      ENDIF
      
      IF(LPV)THEN
        WRITE(INUM,'(''ALTITUDES   (NPROFILE,1-IKU)'')')
      ELSE
        WRITE(INUM,'(''ALTITUDES   (1-NLMAX,1-IKU)'')')
      ENDIF
      WRITE(INUM,'(1X,79(1H*))')
      WRITE(INUM,'(''  K  X->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',79(1H*))')
      DO JLOOPJ=SIZE(PU,2),1,-1
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(XZWORKZ(II,JLOOPJ),II=IDEB,IFIN)
!       WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(XWORKZ(II,JLOOPJ,NMGRID),II=IDEB,IFIN)
!       WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(XWORKZ(II,JLOOPJ,NMGRID),II=IDEB,IFIN)
      ENDDO
      WRITE(INUM,'(1X,79(1H*))')
    ENDDO

  ENDIF

! ENDIF                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENDIF
!!!! Janvier 2001 + LDIRWIND
!------------------------------------------------------------------------------
!
!*    2.    COMPLETING THE PLOT
!           -------------------
!
!*    2.1   Page information labels
!

!print *,' **gsclip N10 0 '
CALL GSCLIP(0)
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
!print *,' getset ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN

!
! Titres en X
!
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
!
! Titres en Y
!
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYM',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYB',ZVL,ZVR,ZVB,ZVT,YTEM)

CALL RESOLV_TIT('CTITB1',HLEGEND)
ZXPOSTITB1=.002
ZXYPOSTITB1=.005
IF(XPOSTITB1 /= 0.)THEN
  ZXPOSTITB1=XPOSTITB1
ENDIF
IF(XYPOSTITB1 /= 0.)THEN
  ZXYPOSTITB1=XYPOSTITB1
ENDIF

IF(HLEGEND /= ' ')THEN
  IF(XSZTITB1 /= 0.)THEN
    CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HLEGEND,XSZTITB1,0.,-1.)
!   CALL PLCHHQ(0.002,0.005,HLEGEND,XSZTITB1,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HLEGEND,.007,0.,-1.)
!   CALL PLCHHQ(0.002,0.005,HLEGEND,.007,0.,-1.)
  ENDIF
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
    CALL PLCHHQ(0.002,0.025,CLEGEND2,XSZTITB2,0.,-1.)
  ELSE
    CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
  ENDIF
ENDIF
YTEM(1:LEN(YTEM))=' '
CALL RESOLV_TIT('CTITB3',YTEM)
ZXPOSTITB3=.002
ZXYPOSTITB3=.0450
IF(XPOSTITB3 /= 0.)THEN
  ZXPOSTITB3=XPOSTITB3
ENDIF
IF(XYPOSTITB3 /= 0.)THEN
  ZXYPOSTITB3=XYPOSTITB3
ENDIF
IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
  IF(XSZTITB3 /= 0.)THEN
    CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,XSZTITB3,0.,-1.)
!   CALL PLCHHQ(0.002,0.050,YTEM,XSZTITB3,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,.009,0.,-1.)
!   CALL PLCHHQ(0.002,0.050,YTEM,.009,0.,-1.)
  ENDIF
ENDIF
  IF(XIDEBCOU.NE.-999.)THEN
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2IND)THEN
        WRITE(YCARCOU,1018)NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
      ELSE IF(LDEFCV2LL)THEN
        WRITE(YCARCOU,1019)XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
      ELSE
        WRITE(YCARCOU,1020)XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
      ENDIF
    ELSE
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
    END IF
  ELSE
    WRITE(YCARCOU,1000)NIDEBCOU,NJDEBCOU,NLANGLE,NLMAX
  END IF
! Janvier 2001
  IF(LPV)THEN
    YCAR(1:LEN(YCAR))=' '
    WRITE(YCAR,1006)NPROFILE
  ENDIF
! Janvier 2001
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
!   CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.012,0.,-1.)
!   CALL PLCHHQ(0.002,0.98,YCARCOU,.012,0.,-1.)
  ENDIF
ENDIF
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
IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
  IF(XSZTITT2 /= 0.)THEN
    CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,XSZTITT2,0.,-1.)
!   CALL PLCHHQ(0.002,0.95,YTEM,XSZTITT2,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,.008,0.,-1.)
!   CALL PLCHHQ(0.002,0.95,YTEM,.008,0.,-1.)
  ENDIF
! Janvier 2001
ELSE
  IF(LPV)THEN
    CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,.008,0.,-1.)
  ENDIF
ENDIF
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
IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
  IF(XSZTITT3 /= 0.)THEN
    CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,XSZTITT3,0.,-1.)
!   CALL PLCHHQ(0.002,0.93,YTEM,XSZTITT3,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,.008,0.,-1.)
!   CALL PLCHHQ(0.002,0.93,YTEM,.008,0.,-1.)
  ENDIF
ENDIF
IF(LDATFILE)CALL DATFILE_FORDIACHRO
ENDIF

! 15 SEPT 2000 Je mets NSUPER a la place de NSUPER
  IF(NSUPER == 1)THEN
    CALL RESOLV_TIT('CTITVAR1',HTEXT)
  ELSE IF(NSUPER == 2)THEN
    CALL RESOLV_TIT('CTITVAR2',HTEXT)
  ELSE IF(NSUPER == 3)THEN
    CALL RESOLV_TIT('CTITVAR3',HTEXT)
  ELSE IF(NSUPER == 4)THEN
    CALL RESOLV_TIT('CTITVAR4',HTEXT)
  ELSE IF(NSUPER == 5)THEN
    CALL RESOLV_TIT('CTITVAR5',HTEXT)
  ELSE IF(NSUPER == 6)THEN
    CALL RESOLV_TIT('CTITVAR6',HTEXT)
  ELSE IF(NSUPER == 7)THEN
    CALL RESOLV_TIT('CTITVAR7',HTEXT)
  ELSE IF(NSUPER == 8)THEN
    CALL RESOLV_TIT('CTITVAR8',HTEXT)
  ENDIF


IF(HTEXT /= ' ')THEN

IF(.NOT.LSUPER)THEN
  CALL PLCHHQ(0.1,ZVT+0.03,HTEXT,.011,0.,-1.)
ELSE
  CALL PLCHHQ(0.1+(NSUPER-1)*.26,ZVT+0.03,HTEXT,.011,0.,-1.)
ENDIF

ENDIF

!CALL PLCHHQ(0.1,ZVT+0.03,HTEXT,.011,0.,-1.)
IF(LVECTMNMX)THEN
  IF(.NOT.LDIRWIND .AND..NOT.LUMVM .AND..NOT.LUTVT .AND..NOT.LSUMVM &
    .AND..NOT.LSUTVT .AND.LDILW)THEN
  CALL PLCHHQ(.1,ZVT+0.010,'(Vertical component upscaled by domain aspect ratio)',.009,0.,-1.) 
  ENDIF
  IF(.NOT.LDIRWIND)THEN
  YLBLMN='          '
  YLBLMX='          '
  WRITE(YLBLMN,'(E10.3)')ZMN
  WRITE(YLBLMX,'(E10.3)')ZMX
  IF(LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. .NOT.LDILW)THEN
    YLBL(1:4)='Min:'
    YLBL(5:14)=YLBLMN
    YLBL(15:20)=', max:'
    YLBL(21:30)=YLBLMX
    YLBL(31:40)=' '
  ELSE
    YLBL(1:13)='Unscaled min:'
    YLBL(14:23)=YLBLMN
    YLBL(24:29)=', max:'
    YLBL(30:39)=YLBLMX
    YLBL(40:40)=' '
  ENDIF
  CALL PCSETC('FC','/')
  CALL PLCHHQ(.99,ZVT+.010,YLBL,.007,0.,+1.)
! CALL PLCHHQ(.69,.047,YLBL,.007,0.,-1.)
  CALL PCSETC('FC',':')
  ENDIF
ENDIF
IF(LSUPER)THEN
  LARROVL=.TRUE.
ELSE
  LARROVL=.FALSE.
END IF
!
!
!*       2.14      Heading formats
!
1000 FORMAT('Vertical section IDEB=',I4,' JDEB=',I4,' ANG.=',I3,' NBPTS=',I4)
1001 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1002 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1003 FORMAT('Vertical section XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1004 FORMAT('Vertical section XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1006 FORMAT('Vertical profile IPRO=',I4)
1018 FORMAT('Vertical section IND I,J (BEGIN)-(END)=(',I4,',',I4,')-(',I4,',',I4,')')
1019 FORMAT('Vertical section LAT,LON (BEGIN)-(END)=(',F4.1,',',F5.1,')-(',F4.1,',',F5.1,')')
1020 FORMAT('Vertical section CONF. COORD.(BEGIN)-(END)=(',F8.0,',',F8.0,')-(',F8.0,',',F8.0,')')
!
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!print *,'imcouv ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
CALL GSCLIP(1)
!
!-------------------------------------------------------------------------
!
!*    3.    EXIT
!           ----
!
RETURN
END SUBROUTINE  IMCOUV_FORDIACHRO
