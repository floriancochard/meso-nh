!     ######spl
      SUBROUTINE IMCOUPV_FORDIACHRO(PU,PW,HLEGEND,HTEXT)
!     #################################################
!
!!****  *IMCOUPV_FORDIACHRO* - Draws a vector arrow plot for a vertical cross-section
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
!!   gridpoint location prior to calling IMCOUPV.
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
USE MODD_COORD
USE MODD_ALLOC_FORDIACHRO
USE MODD_PARAMETERS
USE MODD_NMGRID
USE MODD_GRID
USE MODD_GRID1 
USE MODD_FIELD1_CV2D
USE MODD_SUPER
USE MODD_TITLE
USE MODD_OUT
USE MODN_PARA
USE MODN_NCAR
USE MODD_LUNIT1
USE MODD_CVERT
USE MODD_PVT
USE MODD_TYPE_AND_LH
USE MODD_CTL_AXES_AND_STYL
USE MODD_RESOLVCAR
USE MODD_TIT
USE MODD_DEFCV
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODE_GRIDPROJ
USE MODI_RESOLV_TIT
USE MODI_RESOLV_TITY
!
IMPLICIT NONE
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
INTERFACE

SUBROUTINE GENFORMAT_FORDIACHRO(PCLV,HLLBS)
REAL                :: PCLV
CHARACTER(LEN=*)    :: HLLBS
END SUBROUTINE
!
END INTERFACE
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
INTEGER :: JLOOPI, JLOOPJ, ILOOP, INUM, IRESP,IDEB,IFIN
INTEGER :: JILOOP, JKLOOP, ID, J
INTEGER :: IKB, IKE, IKU
INTEGER      :: IKL, ILMAX, JLMAX
INTEGER      :: ILENYC, ILENHT
INTEGER      :: INBCOL, IIBID
INTEGER      :: JA, JILOOPD, JILOOPF
INTEGER      :: JJ, IJ, II, IUB1, IUB2, ITER, JTER
INTEGER      :: ISKIPX, ISKIPY, ITERM, ISKIPXM
INTEGER,DIMENSION(:),ALLOCATABLE      :: ICOL
!
REAL,DIMENSION(SIZE(PU,2),SIZE(PU,1)) :: ZZU, ZZV
REAL :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
REAL :: ZY, ZJ, ZH, ZJJ, ZWBB
REAL :: ZDMX, ZVMX
REAL :: ZRAP
REAL :: ZXPOSTITT1, ZXYPOSTITT1
REAL :: ZXPOSTITT2, ZXYPOSTITT2
REAL :: ZXPOSTITT3, ZXYPOSTITT3
REAL :: ZXPOSTITB1, ZXYPOSTITB1
REAL :: ZXPOSTITB2, ZXYPOSTITB2
REAL :: ZXPOSTITB3, ZXYPOSTITB3
REAL,DIMENSION(1000) :: ZYY
REAL :: ZW,ZM,ZUMN,ZWMN,ZMN,ZWMX,ZMX
REAL,DIMENSION(:),ALLOCATABLE      :: ZPARCOLUV
REAL :: ZTEM, ZINT, ZRPK, ZLON0, ZBETA
REAL :: ZVINT, ZVY, ZINTX, ZINTY
REAL,DIMENSION(:,:),ALLOCATABLE    :: ZX, ZLAT, ZLON, ZZY,ZZYY
CHARACTER(LEN=4) :: YTE
REAL,DIMENSION(:,:),ALLOCATABLE    :: ZDIRU, ZDIRV, ZLA, ZLO  
REAL,DIMENSION(:),ALLOCATABLE    :: ZZDS
REAL,DIMENSION(18) :: ZCOL

CHARACTER(LEN=82) :: YCARCOU, YTEM
CHARACTER(LEN=80) :: YCAR
CHARACTER(LEN=40) :: YLBL
CHARACTER(LEN=40) :: YTIT
CHARACTER(LEN=8),DIMENSION(:),ALLOCATABLE :: YLBS
CHARACTER(LEN=8) :: YLBSTEM
CHARACTER(LEN=2) :: YC2
CHARACTER(LEN=3) :: YC3
CHARACTER(LEN=4) :: YC4
CHARACTER(LEN=10) :: YLBLMN,YLBLMX
CHARACTER(LEN=10) :: FORMAX, FORMAY
!
!*       0.4   External for NCAR use
!
! SFILL subroutine declared as external provides area control
! in some parts of the contour plot.
!
!EXTERNAL SFILL
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

!!!! ATTENTION  En entree ICI,PU (U) et PW (V ICI) ont comme 1ere dimension
!!!! Z (1:IKU) et comme 2eme le temps (qui au trace sera en X) -> besoin
!!!! de retablir l'ordre habituel : (Tps,Z) ce qui est fait ds ZZU et ZZV

ILMAX=SIZE(PU,2)
JLMAX=SIZE(PU,1)
if(nverbia > 0)then
print *, ' ENTREE imcoupv ',ILMAX,JLMAX
endif

ZZU=XSPVAL
ZZV=XSPVAL

! Janvier 2001
!IF(.NOT.LUMVMPV)THEN

  DO JKLOOP=1,JLMAX
  DO JILOOP=1,ILMAX
    ZZU(JILOOP,JKLOOP)=PU(JKLOOP,JILOOP)
    ZZV(JILOOP,JKLOOP)=PW(JKLOOP,JILOOP)
  ENDDO
  ENDDO

!ELSE

! Janvier 2001
! DO JKLOOP=1,JLMAX,NISKIPVY
! DO JILOOP=1,ILMAX,NISKIPVX
!   ZZU(JILOOP,JKLOOP)=PU(JKLOOP,JILOOP)
!   ZZV(JILOOP,JKLOOP)=PW(JKLOOP,JILOOP)
! ENDDO
! ENDDO

! Janvier 2001
!ENDIF
! Janvier 2001
!
!
!*       1.2  Collects X and Z values 
!
!*       1.3  Window definition and plot
!

LVERTI=.TRUE. ; LHORIZ=.FALSE.
LVERT=LVERTI
LHOR=LHORIZ

CALL GSCLIP(0)

CALL GSLN(1)
CALL GSPLCI(1)
CALL GSTXCI(1)

!IF(LSUPER)THEN
! NSUPER=NSUPER+1
! IF(NSUPER == 1)THEN
!   CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
! ELSE
!   CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! END IF
!ELSE
! CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
!ENDIF

CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)

!!!!!!!!!!!!!!!
FORMAX='          '
IF(LFMTAXEX)THEN
  FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
ELSE
  FORMAX='(F8.1)'
ENDIF

FORMAY='          '
IF(LFMTAXEY)THEN
  FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
ELSE
  FORMAY='(F7.0)'
ENDIF
!!!!!!!OCt 2001
!IF(ZWL == ZWR)ZWR=ZWL*2
!!!!!!!OCt 2001

IF(LHEURX)THEN
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL/3600.,ZWR/3600.,ZWB,ZWT,ID)
  CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)

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

ELSE

  CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
ENDIF
!!!!!!!!!!!!!!!

! Utilisation de PLCHHQ pour ecriture des labels (sinon 0= WTSTR)
CALL GASETI('LTY',1)

IF(.NOT.LHEURX)THEN
! Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0)
  ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0)
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,1,0,5,0.,0)
  ELSE
    CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NCVITVYMJ,NCVITVYMN,1,1,5,0.,0)
  ENDIF
! Avril 2002
ENDIF

!!!!!!!!!!!!!!!
IF(LHEURX)THEN
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  DO J=INT(ZWL),INT(ZWR)
    ZJ=J

    IF(MOD(ZJ,ZH) == 0.)THEN
      CALL FRSTPT(ZJ,ZWB)
      CALL VECTOR(ZJ,ZWB+(ZWT-ZWB)/90.)

!!!!!!!Avril 2002
  IF(LMYHEURX)THEN
    ZJJ=ZJ/ZH*NHEURXGRAD
    ZINT=NHEURXLBL
  ELSE
!!!!!!!Avril 2002
      IF(ZH == 10800.)THEN
        ZJJ=ZJ/ZH*3.
        ZINT=6.
      ELSE
        ZJJ=ZJ/ZH
        ZINT=3.
      ENDIF

!!!!!!!Avril 2002
  ENDIF
!!!!!!!Avril 2002
      ZWBB=ZWB-((ZWT-ZWB)/((ZVT-ZVB)/.02))

      IF(.NOT. LNOLABELX)THEN
      IF(MOD(ZJJ,ZINT) == 0.)THEN
        IF(ZJJ <10.)THEN
          WRITE(YC2,'(F2.0)')ZJJ
          CALL PLCHHQ(ZJ,ZWBB,YC2,.010,0.,0.)
        ELSEIF(ZJJ <100.)THEN
          WRITE(YC3,'(F3.0)')ZJJ
          CALL PLCHHQ(ZJ,ZWBB,YC3,.010,0.,0.)
        ELSE
          WRITE(YC4,'(F4.0)')ZJJ
          CALL PLCHHQ(ZJ,ZWBB,YC4,.010,0.,0.)
        ENDIF
      ENDIF
      ENDIF

    ENDIF
  ENDDO
! CALL GRIDAL(1,0,NCVITVYMJ,NCVITVYMN,1,1,5,0.,0)
! Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0)
  ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
    CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0)
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,0,5,0.,0)
  ELSE
    CALL GRIDAL(0,0,NCVITVYMJ,NCVITVYMN,0,1,5,0.,0)
  ENDIF
! Avril 2002
ENDIF
!!!!!!!!!!!!!!!

! Janvier 2001
!!! Partie commune de LPRINT
IF(LPRINT)THEN
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=SIZE(ZZU,1)/5
  IF(ILOOP * 5 < SIZE(ZZU,1))ILOOP=ILOOP+1

  IF(.NOT.LPVT)THEN
    WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'',''   (1-1,1-IKU)'')')CGROUP,&
&   CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
  ELSE
    IUB1=SIZE(ZZU,1)
    WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' TD-TF:'',F8.0,''-'',F8.0,''s'')')CGROUP,&
    CTITRE(NLOOPP)(1:25),XZZDS(1),XZZDS(IUB1)
    WRITE(INUM,'('' (1-NBTIME,1-IKU)'')')
  ENDIF

  IF(LMINUS .OR. LPLUS)THEN
    WRITE(INUM,'(A70)')CTITB3
  ELSE
!   WRITE(INUM,'(A40)')CTITGAL
  ENDIF

  IF(LUMVMPV)THEN
    WRITE(INUM,'(''I='',I4,''J='',I4)')&
    NIL,NJL
  ELSE
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2)THEN
        WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
    &'' profile='',I4)')&
       &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,NPROFILE
      ELSE IF(LDEFCV2LL)THEN
        WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
    &'' profile='',I4)')&
       &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,NPROFILE
      ELSE IF(LDEFCV2IND)THEN
        WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
    &'' profile='',i4)')&
       &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,NPROFILE
      ENDIF
    ELSE
      IF(XIDEBCOU /= -999.)THEN
        WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
    &'' profile='',i4)')&
       &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,NPROFILE
      ELSE
        WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
    &'' profile='',i4)')&
       &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,NPROFILE
      ENDIF
    ENDIF
!   WRITE(INUM,'(''nprofile='',I4)')NPROFILE
  ENDIF

    WRITE(INUM,'(''NBVAL en I (TIME): '',i4, &
&  '' NBVAL en K (Z)'',i4,''    iter'',i3)') &
  & SIZE(ZZU,1),SIZE(ZZU,2),ILOOP
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**IMCOUPV XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
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
ENDIF 

!!!! Janvier 2001 + LDIRWIND
IF(LDIRWIND)THEN
  if(nverbia > 0)then
  print *,' imcoupv LDIRWIND ',LDIRWIND
  endif
  ISKIPX=NISKIPVX
  ISKIPY=NISKIPVY
  IUB1=SIZE(ZZU,1)
!!30/01/01
! ITER=IUB1/ISKIPX+1
! IF(1+(ITER-1)*ISKIPX > IUB1)ITER=ITER-1
  ITERM=IUB1/ISKIPX+1
  IF(1+(ITERM-1)*ISKIPX > IUB1)ITERM=ITERM-1
  ITER=IUB1
  ISKIPXM=ISKIPX
  ISKIPX=1
!!30/01/01
  IUB2=SIZE(ZZU,2)
! 130101
!!! Essai de conservation de 1 a IKU en Y (pour LPRINT) mais
!!! de 1 a ITER en X
!!!  JTER=(IUB2-IKB)/ISKIPY+1
!!!  IF(IKB+(JTER-1)*ISKIPY > IUB2)JTER=JTER-1
  JTER=IUB2
!!!
  ALLOCATE(ZX(ITER,1),ZZY(ITER,JTER),ZZYY(ITER,1),ZLAT(ITER,1),ZLON(ITER,1))
  ALLOCATE(ZLA(ITER,JTER),ZLO(ITER,JTER),ZDIRU(ITER,JTER),ZDIRV(ITER,JTER))
  ALLOCATE(ZZDS(ITER))
! 130101
! print *,' IIIIIMCOUPV IUB1, ISKIPX, ITER, IUB2, ISKIPY, JTER,LPV ',IUB1,ISKIPX,ITER,IUB2,ISKIPY,JTER,LPV

!!!
  ZDIRU=XSPVAL
  ZDIRV=XSPVAL
!!!  ZDIRU=ZZU(1:IUB1:ISKIPX,IKB:IUB2:ISKIPY)
!!!  ZDIRV=ZZV(1:IUB1:ISKIPX,IKB:IUB2:ISKIPY)
  ZDIRU=ZZU(1:IUB1:ISKIPX,1:IUB2:1)
  ZDIRV=ZZV(1:IUB1:ISKIPX,1:IUB2:1)
!!!
  if(nverbia > 0)then
    print *,' ZDIRU AP CHARG. ZZU'
    print *,ZDIRU 
    print *,' ZDIRV AP CHARG. ZZV'
    print *,ZDIRV
  endif

! Chargement des temps ICI .
  ZZDS=XTDIRWIND(1:IUB1:ISKIPX)
! print *,' IIIIIMCOUPV XDSX(1:IUB1) ',XDSX(1:IUB1,1)
! print *,' IIIIIMCOUPV ZX(:,1) ',ZX(:,1)
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
! print *,' IIIIMCOUPV IUB1,ISKIPX,IKB,IUB2,ISKIPY ',IUB1,ISKIPX,IKB,IUB2
! print *,' IIIIMCOUPV XZWORKZ(1:NLMAX,IKB) ',XZWORKZ(1:NLMAX,IKB)
! print *,' IIIIMCOUPV ZZY(:,1) ',ZZY(:,1)
! print *,' IIIIMCOUPV XZWORKZ(1:NLMAX,IKB+1) ',XZWORKZ(1:NLMAX,IKB+1)
! print *,' IIIIMCOUPV ZZY(:,2) ',ZZY(:,2)

! 130101
  ZX(:,1)=XDSX(1,1)
  ZZYY(:,1)=XDSY(1,1)

  IF(ALLOCATED(ICOL))THEN
    DEALLOCATE(ICOL)
  ENDIF
  ALLOCATE(ICOL(18))

  DO JKLOOP=1,JTER
    CALL SM_LATLON_A(XLATORI,XLONORI,ZX,ZZYY,ZLAT,ZLON)
    ZLA(:,JKLOOP)=ZLAT(:,1)
    ZLO(:,JKLOOP)=ZLON(:,1)
  ENDDO

  where(zdiru /= xspval .AND. zdirv /= xspval)
    ZDIRU=ATAN2(ZDIRV,ZDIRU)*180./ACOS(-1.)
  endwhere

  if(nverbia > 0)then
    print *,' ZDIRU AP ATAN2 '
    print *,ZDIRU 
    print *,' ZDIRU 1,1 ITER/2,1 1,JTER/2 ITER/2,JTER/2 ITER,JTER  '
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
    print *,' ZDIRU AP WHERE(ZDIRU < 0.'
    print *,ZDIRU 
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
      if(nverbia > 0)then
	print *,' AV LPRINT DIRWIND ZDIRU '
	print *, ZDIRU
	print *,' AV LPRINT DIRWIND ZDIRV '
	print *, ZDIRV
      endif

  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)

if(nverbia > 0)then
  print *,' ** imcoupv ap getset ZWL,ZWR,XDS(1,1),XDS(NLMAX,1),ZX(1,1),ZX(ITER,1) ',ZWL,ZWR,XDS(1,1),XDS(NLMAX,1),ZX(1,1),ZX(ITER,1)
endif
!! 30/01/01
  IF(ITERM > 6)THEN
! IF(ITER > 6)THEN
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
!! 30/01/01
    DO JILOOP=JILOOPD,JILOOPF,ISKIPXM
!   DO JILOOP=JILOOPD,JILOOPF
!! 30/01/01
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
!       ZINTX=(ZWL+ZWR)/2
        ZINTX=ZZDS(JILOOP)
      ELSE
        ZINTX=ZZDS(JILOOP)
!       print *,' **imcoupv ZINTX ',ZINTX
      ENDIF

      ZINTY=ZZY(JILOOP,JKLOOP)
      IF(ZINTY < XHMIN .OR. ZINTY > XHMAX)THEN
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
        IFIN=SIZE(ZZU,1)
      ENDIF
      
      if(nverbia > 0)then
	print *,' ds LPRINT DIRWIND ZDIRU '
	print *, ZDIRU
	print *,' ds LPRINT DIRWIND ZDIRV '
	print *, ZDIRV
      endif
      WRITE(INUM,'(1X,79(1H*))')
      WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',79(1H*))')
      DO JLOOPJ=SIZE(ZZU,2),1,-1
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(ZDIRV(II,JLOOPJ),II=IDEB,IFIN)
  !     WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(ZDIRV(II,JLOOPJ),II=IDEB,IFIN)
      ENDDO
      WRITE(INUM,'(1X,79(1H*))')
    ENDDO
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF

  IF(LPRINTXY)THEN
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
    IF(IRESP /= 0)THEN
      CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
      OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
      PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
    ENDIF

    IF(.NOT.LPVT)THEN
      WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'',''   (1-1,1-IKU)'')')CGROUP,&
  &   CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ELSE
      WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' TD-TF:'',F8.0,''-'',F8.0,''s'')')CGROUP,&
      CTITRE(NLOOPP)(1:25),XZZDS(1),XZZDS(IUB1)
      WRITE(INUM,'('' (1-NBTIME,1-IKU)'')')
    ENDIF
  
    IF(LMINUS .OR. LPLUS)THEN
      WRITE(INUM,'(A70)')CTITB3
    ELSE
  !   WRITE(INUM,'(A40)')CTITGAL
    ENDIF
  
    IF(LUMVMPV)THEN
      WRITE(INUM,'(''I='',I4,''J='',I4)')&
      NIL,NJL
    ELSE
      IF(LDEFCV2CC)THEN
        IF(LDEFCV2)THEN
          WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
      &'' profile='',I4)')&
         &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,NPROFILE
        ELSE IF(LDEFCV2LL)THEN
          WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
      &'' profile='',I4)')&
         &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,NPROFILE
        ELSE IF(LDEFCV2IND)THEN
          WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
      &'' profile='',i4)')&
         &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,NPROFILE
        ENDIF
      ELSE
        IF(XIDEBCOU /= -999.)THEN
          WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
      &'' profile='',i4)')&
         &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,NPROFILE
        ELSE
          WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
      &'' profile='',i4)')&
         &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,NPROFILE
        ENDIF
      ENDIF
  !   WRITE(INUM,'(''nprofile='',I4)')NPROFILE
    ENDIF
  
      WRITE(INUM,'(''NBVAL en I (TIME): '',i4, &
  &  '' NBVAL en K (Z)'',i4,''    iter'',i3)') &
    & SIZE(ZZU,1),SIZE(ZZU,2),ILOOP
  
      II=MAX(SIZE(ZZU,1),SIZE(ZZU,2))
      WRITE(INUM,'(1X,43(1H*))')
      WRITE(INUM,'(2X,''  I'',7X,''TIME'',10X,''K'',9X,''Z'')')
      WRITE(INUM,'(1X,43(1H*))')
      DO JLOOPJ=1,II
        IF(SIZE(ZZU,1) > SIZE(ZZU,2))THEN
          IF(JLOOPJ <= SIZE(ZZU,2))THEN
             WRITE(INUM,'(I5,2X,E15.8,1X,I4,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ), &
            JLOOPJ,XZWORKZ(1,JLOOPJ)
          ELSE
            WRITE(INUM,'(I5,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ)
          ENDIF
        ELSE IF(SIZE(ZZU,2) > SIZE(ZZU,1))THEN
          IF(JLOOPJ <= SIZE(ZZU,1))THEN
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
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF

  CALL GSCLIP(0)
  DEALLOCATE(ZX,ZZY,ZZYY,ZLAT,ZLON,ZLA,ZLO,ZDIRU,ZDIRV,ICOL,ZZDS)

ELSE

!!!! Janvier 2001 + LDIRWIND
  IF(LPRINT)THEN
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DO JLOOPI=1,ILOOP
      IF(JLOOPI == 1)THEN
        IDEB=1; IFIN=5
      ELSE
        IDEB=IFIN+1; IFIN=IFIN+5
      ENDIF
      IF(JLOOPI == ILOOP)THEN
        IFIN=SIZE(ZZU,1)
      ENDIF
      
      IF(INDEX(CGROUP,'UM') /= 0)THEN
        WRITE(INUM,'(1X,20(1H*),'' UM  component '',34(1H*))')
      ELSE
        WRITE(INUM,'(1X,20(1H*),'' UT  component '',34(1H*))')
      ENDIF
      if(nverbia > 0)then
	print *,' ds LPRINT ZZU'
	print *, ZZU
      endif
!     WRITE(INUM,'(1X,79(1H*))')
      WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',79(1H*))')
      DO JLOOPJ=SIZE(ZZU,2),1,-1
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(ZZU(II,JLOOPJ),II=IDEB,IFIN)
  !     WRITE(INUM,'(I4,1X,5E15.8)')JLOOPJ,(ZZU(II,JLOOPJ),II=IDEB,IFIN)
      ENDDO
      WRITE(INUM,'(1X,79(1H*))')

      IF(INDEX(CGROUP,'VM') /= 0)THEN
        WRITE(INUM,'(1X,20(1H*),'' VM  component '',34(1H*))')
      ELSE
        WRITE(INUM,'(1X,20(1H*),'' VT  component '',34(1H*))')
      ENDIF
      WRITE(INUM,'(''  K  I->   '',I4,6X,4(6X,I4,6X))')(/(II,II=IDEB,IFIN)/)
      WRITE(INUM,'(''.'',79(1H*))')
      DO JLOOPJ=SIZE(ZZV,2),1,-1
        WRITE(INUM,'(I4,1X,5(1X,E14.7))')JLOOPJ,(ZZV(II,JLOOPJ),II=IDEB,IFIN)
      ENDDO
      WRITE(INUM,'(1X,79(1H*))')
    ENDDO
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF

  IF(LPRINTXY)THEN
                                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
    IF(IRESP /= 0)THEN
      CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
      OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
      PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
    ENDIF

    IF(.NOT.LPVT)THEN
      WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'',''   (1-1,1-IKU)'')')CGROUP,&
  &   CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
    ELSE
      WRITE(INUM,'(''PV  '',''G:'',A16,'' P:'',A25,'' TD-TF:'',F8.0,''-'',F8.0,''s'')')CGROUP,&
      CTITRE(NLOOPP)(1:25),XZZDS(1),XZZDS(IUB1)
      WRITE(INUM,'('' (1-NBTIME,1-IKU)'')')
    ENDIF
  
    IF(LMINUS .OR. LPLUS)THEN
      WRITE(INUM,'(A70)')CTITB3
    ELSE
  !   WRITE(INUM,'(A40)')CTITGAL
    ENDIF
  
    IF(LUMVMPV)THEN
      WRITE(INUM,'(''I='',I4,''J='',I4)')&
      NIL,NJL
    ELSE
      IF(LDEFCV2CC)THEN
        IF(LDEFCV2)THEN
          WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
      &'' profile='',I4)')&
         &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,NPROFILE
        ELSE IF(LDEFCV2LL)THEN
          WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
      &'' profile='',I4)')&
         &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,NPROFILE
        ELSE IF(LDEFCV2IND)THEN
          WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
      &'' profile='',i4)')&
         &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,NPROFILE
        ENDIF
      ELSE
        IF(XIDEBCOU /= -999.)THEN
          WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
      &'' profile='',i4)')&
         &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,NPROFILE
        ELSE
          WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
      &'' profile='',i4)')&
         &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,NPROFILE
        ENDIF
      ENDIF
  !   WRITE(INUM,'(''nprofile='',I4)')NPROFILE
    ENDIF
  
      WRITE(INUM,'(''NBVAL en I (TIME): '',i4, &
  &  '' NBVAL en K (Z)'',i4,''    iter'',i3)') &
    & SIZE(ZZU,1),SIZE(ZZU,2),ILOOP
  
      II=MAX(SIZE(ZZU,1),SIZE(ZZU,2))
      WRITE(INUM,'(1X,43(1H*))')
      WRITE(INUM,'(2X,''  I'',7X,''TIME'',10X,''K'',9X,''Z'')')
      WRITE(INUM,'(1X,43(1H*))')
      DO JLOOPJ=1,II
        IF(SIZE(ZZU,1) > SIZE(ZZU,2))THEN
          IF(JLOOPJ <= SIZE(ZZU,2))THEN
             WRITE(INUM,'(I5,2X,E15.8,1X,I4,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ), &
            JLOOPJ,XZWORKZ(1,JLOOPJ)
          ELSE
            WRITE(INUM,'(I5,2X,E15.8)')JLOOPJ,XZZDS(JLOOPJ)
          ENDIF
        ELSE IF(SIZE(ZZU,2) > SIZE(ZZU,1))THEN
          IF(JLOOPJ <= SIZE(ZZU,1))THEN
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
                                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ENDIF
! Janvier 2001

  ZZU=XSPVAL
  ZZV=XSPVAL
  IF(.NOT.LUMVMPV)THEN
    DO JKLOOP=IKB,JLMAX,NISKIPVY
    DO JILOOP=1,ILMAX,NISKIPVX
      ZZU(JILOOP,JKLOOP)=PU(JKLOOP,JILOOP)
      ZZV(JILOOP,JKLOOP)=PW(JKLOOP,JILOOP)
    ENDDO
    ENDDO

  ELSE

    DO JKLOOP=1,JLMAX,NISKIPVY
    DO JILOOP=1,ILMAX,NISKIPVX
      ZZU(JILOOP,JKLOOP)=PU(JKLOOP,JILOOP)
      ZZV(JILOOP,JKLOOP)=PW(JKLOOP,JILOOP)
    ENDDO
    ENDDO

  ENDIF
! Janvier 2001

!
!*       1.4  Collects wind values within the user postprocessing
!*            window with a sampling rate of NISKIP outside values 
!*            are kept to default
!

CALL GSCLIP(0)
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
CALL VVSETI('MNP',-4)
CALL VVSETI('MXP',-4)
CALL VVSETR('MNX',.75)
!CALL VVSETR('MNX',-ZVL)
!ZY=-1./5.
!ZY=-MIN(0.12,ZVB+.02)
IF(ZVB <= .15)THEN
  ZY=-ZVB-.020
! ZY=(-.08)/(ZVT-ZVB)
ELSE
!!! Octobre 2001
! ZY=(-.10)/(ZVT-ZVB)
  ZY=(-.13)/(ZVT-ZVB)
!!! Octobre 2001
ENDIF
!IF(ZVB-(ZVT-ZVB)/5..LT.0.05)ZY=(0.05-ZVB)/(ZVT-ZVB)
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

IF(ILMAX > 6)THEN
CALL GSCLIP(1)                                     ! Clipping off
ENDIF
CALL VVSETI('VPO',1)
CALL VVINIT(ZZU,ILMAX,ZZV,ILMAX,0.,0,ILMAX,IKU,0.,0) ! Initializes VVECTR
CALL VVECTR(ZZU,ZZV,0.,0,0,0.)                     ! Draws arrows
CALL GSCLIP(0)                                     ! Clipping back on
!
CALL VVRSET
!------------------------------------------------------------------------------
!
!*    2.    COMPLETING THE PLOT
!           -------------------
!
!*    2.1   Page information labels
!

CALL GSCLIP(0)

CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
!print *,' getset ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT

IF(LCOLPVT)THEN
!print *,' ** imcoupv AP LCOLPVT '

   IF(LCOLUSERUV)THEN
     INBCOL=NBPARCOLUV
     IF(ALLOCATED(ICOL))THEN
       DEALLOCATE(ICOL)
     ENDIF
     ALLOCATE(ICOL(NBCOLUV))
     ALLOCATE(YLBS(NBCOLUV-1))
     ALLOCATE(ZPARCOLUV(NBCOLUV-1))
     ICOL(:)=NINDCOLUV(1:NBCOLUV)
     ZPARCOLUV=XPARCOLUV(1:NBCOLUV-1)
   ELSE
     INBCOL=NBPARCOLUVSTD
     IF(ALLOCATED(ICOL))THEN
       DEALLOCATE(ICOL)
     ENDIF
     ALLOCATE(ICOL(NBCOLUVSTD))
     ALLOCATE(YLBS(NBCOLUVSTD-1))
     ALLOCATE(ZPARCOLUV(NBCOLUVSTD-1))
     ICOL(:)=NCOLUVSTD(1:NBCOLUVSTD)
     ZPARCOLUV=XPARCOLUVSTD(1:NBCOLUVSTD-1)
   ENDIF

   YLBS(:)=' '
!print *,' ** imcoupv AV GENFORMAT '

   DO J=1,INBCOL
     ZTEM=ZPARCOLUV(J)
     CALL GENFORMAT_FORDIACHRO(ZTEM,YLBSTEM)
!    CALL GENFORMAT_FORDIACHRO(ZPARCOLUV(J),YLBS(J))
     YLBS(J)=YLBSTEM
   ENDDO

!print *,' ** imcoupv AP GENFORMAT '
   CALL GSFAIS(1)
   CALL LBLBAR_FORDIACHRO(0,ZVL,ZVR,ZVT+.01,ZVT+.05,INBCOL+1,1.,.15,ICOL,&
   1,YLBS,INBCOL,2)
   CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,ID)
   YTIT(1:LEN(YTIT))=' '
   YTIT=CTITRE(NLOOPP)

   YTIT=ADJUSTR(YTIT)
!  print *,' **imcoupv YTIT NLOOPP ',YTIT,NLOOPP
      CALL PLCHHQ(MIN(ZVR+.1,1.),ZVT+.02,YTIT(1:LEN_TRIM(YTIT)),.007,0.,+1.)
   CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
   DEALLOCATE(ICOL)
   DEALLOCATE(YLBS)
   DEALLOCATE(ZPARCOLUV)
ENDIF

!!! Janvier 2001
ENDIF
!print *,' **imcoupv AV SET(0.,1.,0.,1.,0.,1.,0.,1.,1)'
!!! Janvier 2001

CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN

!
! Titres en X
!
!-------------------------------------------------------------------
  YTEM(1:LEN(YTEM))=' '
  YTEM=ADJUSTL(YTEM)
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
  IF(LHEURX)THEN
    YTEM='(H)'
  ELSE
    YTEM='(sec)'
  ENDIF
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
!-------------------------------------------------------------------
  YTEM(1:LEN(YTEM))=' '
  YTEM='Altitude;(ms)'
  CALL RESOLV_TITY('CTITYM',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYB',ZVL,ZVR,ZVB,ZVT,YTEM)

! Titres Bottom
!-------------------------------------------------------------------
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
! Octobre 2001
YTEM=CTIMEC
YTEM=ADJUSTL(CTIMEC)
! Octobre 2001
CALL RESOLV_TIT('CTITB3',YTEM)
ZXPOSTITB3=.002
ZXYPOSTITB3=.050
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
! Titres Top
!-------------------------------------------------------------------
! Janv 2001
   IF(.NOT.LUMVMPV)THEN
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

   ELSE
     WRITE(YCARCOU,1021)NIl,NJL
   ENDIF

! Janvier 2001
! Conversion METERS/SECONDE en M/S
IIBID=INDEX(HTEXT,'METERS/SECONDE')
ILENHT=LEN_TRIM(HTEXT)
IF(IIBID /= 0)THEN
IF(HTEXT(IIBID:ILENHT) == 'METERS/SECONDE')THEN
  HTEXT(IIBID:ILENHT)=' '
  HTEXT(IIBID:IIBID+2)='M/S '
ENDIF
ENDIF

IF(LUMVMPV)THEN
! Janvier 2001
IF(HTEXT/= ' ')THEN
! print *,' ** imcoupv CUNITGAL ',CUNITGAL
  ILENYC=LEN_TRIM(YCARCOU)
  ILENHT=LEN_TRIM(HTEXT)
  YCARCOU(ILENYC+1:ILENYC+3)=' '
  YCARCOU(ILENYC+4:ILENYC+ILENHT+4-1)=HTEXT(1:ILENHT)
! ILENYC=LEN_TRIM(YCARCOU)
! ILENHT=LEN_TRIM(CUNITGAL)
! YCARCOU(ILENYC+1:ILENYC+1)=' '
ENDIF
! Janvier 2001
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
  IF(.NOT.LUMVMPV)THEN
    YCAR(1:LEN(YCAR))=' '
    WRITE(YCAR,1006)NPROFILE
    CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,.008,0.,-1.)
  ENDIF
! Janvier 2001
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
!-------------------------------------------------------------------
IF(LDATFILE)CALL DATFILE_FORDIACHRO
ENDIF

!-------------------------------------------------------------------
  IF(NSUPERDIA == 1)THEN
    CALL RESOLV_TIT('CTITVAR1',HTEXT)
  ELSE IF(NSUPERDIA == 2)THEN
    CALL RESOLV_TIT('CTITVAR2',HTEXT)
  ELSE IF(NSUPERDIA == 3)THEN
    CALL RESOLV_TIT('CTITVAR3',HTEXT)
  ELSE IF(NSUPERDIA == 4)THEN
    CALL RESOLV_TIT('CTITVAR4',HTEXT)
  ELSE IF(NSUPERDIA == 5)THEN
    CALL RESOLV_TIT('CTITVAR5',HTEXT)
  ELSE IF(NSUPERDIA == 6)THEN
    CALL RESOLV_TIT('CTITVAR6',HTEXT)
  ELSE IF(NSUPERDIA == 7)THEN
    CALL RESOLV_TIT('CTITVAR7',HTEXT)
  ELSE IF(NSUPERDIA == 8)THEN
    CALL RESOLV_TIT('CTITVAR8',HTEXT)
  ENDIF


! Janvier 2001
 IF(.NOT.LUMVMPV)THEN
! Janvier 2001
 IF(HTEXT /= ' ')THEN
 IF(.NOT.LSUPER)THEN
  IF(XSZTITVAR1 /= 0.)THEN
    CALL PLCHHQ(0.1,ZVT+0.03,HTEXT,XSZTITVAR1,0.,-1.)
  ELSE
    CALL PLCHHQ(0.1,ZVT+0.03,HTEXT,.011,0.,-1.)
  ENDIF
 ELSE
  IF(XSZTITVAR1 /= 0. .AND. NSUPER == 1)THEN
    CALL PLCHHQ(0.1+(NSUPER-1)*.26,ZVT+0.03,HTEXT,XSZTITVAR1,0.,-1.)
  ELSE
    CALL PLCHHQ(0.1+(NSUPER-1)*.26,ZVT+0.03,HTEXT,.011,0.,-1.)
  ENDIF
 ENDIF
 ENDIF
! Janvier 2001
 ENDIF
! Janvier 2001
!-------------------------------------------------------------------
IF(LSUPER)THEN
  LARROVL=.TRUE.
ELSE
  LARROVL=.FALSE.
END IF
!
!!!!!!!!!!!!!!!!!!!!!!
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!
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
1021 FORMAT('Vertical profile I=',I4,' J=',I4)
!
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!print *,'imcoupv ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
CALL GSCLIP(1)
!
!-------------------------------------------------------------------------
!
!*    3.    EXIT
!           ----
!
RETURN
END SUBROUTINE  IMCOUPV_FORDIACHRO
