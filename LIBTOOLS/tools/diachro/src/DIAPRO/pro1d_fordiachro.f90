!     ######spl
      SUBROUTINE PRO1D_FORDIACHRO(KPRO,PPRO,PTAB,PTABMIN,PTABMAX,KXDOT, &
      HLEGEND,HTEXT)
!     ####################################################################
!
!!****  *PRO1D_FORDIACHRO* - Draws vertical profiles
!!
!!    PURPOSE
!!    -------
!       Draws vertical profiles.
!
!!**  METHOD
!!    ------
!!      The  NCAR autograph utility is called with appropriate
!!      scaling and headers.
!!
!!    EXTERNAL
!!    --------
!!      SET      : defines the display window in normalized device  !
!!                 coordinate and user coordinate.                  !
!!      LABMOD   : defines the label formats                        ! NCAR
!!      GRIDAL   : draws axes, perimeter, ticks, and labels         !
!!      GSCLIP   : prevents out of window plotting                  !
!!      GSFAIS   : color filling  iusing  GKS                       !
!!      PLCHHQ   : prints high quality texts on graphics            ! routines 
!!      EZXY     : compact utility to draw a Y=f(X) function plot   ! 
!!      AGSETF   : sets an NCAR parameter to a el value (AUTOGRAPH) !
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
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
!!         XHMIN             :  Altitude of the vert. cross-section
!!                              bottom (in meters above sea-level)
!!         XHMAX             :  Altitude of the vert. cross-section
!!                              top (in meters above sea-level)
!!         Module MODD_DIM1 : contains dimensions of data arrays
!!               NKMAX       : z array dimension
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!         JPVEXT : Vertical external points number
!!
!!      Module MODD_ALLVAR  : contains generic variables arrays and structures
!!         L1DT  : Logical identifying the current generic variable as a 1D
!!                    scalar variable  when .TRUE. No signification otherwise.
!!
!!      Module MODD_GRID1      : declares grid variables (Model module)
!!         XZZ      : true gridpoint z altitude
!!
!!      Module MODD_SUPER   : defines plot overlay control variables
!!         LSUPER   : =.TRUE. --> plot overlay is active
!!                    =.FALSE. --> plot overlay is not active
!!         NSUPER   : Rank of the current plot in the overlay
!!                    sequence. The initial plot is rank 1.
!!      Module MODD_TITLE  : Declares heading variables for the plots (TRACE)
!!         CLEGEND2 : Current plot heading title
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
!!      Updated   PM   13/01/95
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!

USE MODN_NCAR
USE MODN_PARA
USE MODD_PARAMETERS
USE MODD_RESOLVCAR 
USE MODD_ALLVAR
USE MODD_GRID1
USE MODD_CONF
USE MODD_DEFCV
USE MODD_SUPER
USE MODD_TITLE
USE MODD_OUT
USE MODD_TYPE_AND_LH
USE MODD_TIT
USE MODD_EXPERIM
USE MODD_ALLOC_FORDIACHRO
USE MODD_CTL_AXES_AND_STYL

IMPLICIT NONE
!
!*       0.1   Dummy arguments and results
!
INTEGER           :: KPRO      ! Profile gridpont index along section x-axis
REAL,DIMENSION(:) :: PPRO      ! Data array to be plotted
REAL,DIMENSION(:) :: PTAB      ! Altitude array for the profile
REAL              :: PTABMIN   ! Minimum altitude of the profile
REAL              :: PTABMAX   ! Maximum altitude of the profile 
INTEGER           :: KXDOT     ! Number of major division along abscissa
CHARACTER(LEN=*)  :: HLEGEND   ! Name of the variable header
CHARACTER(LEN=*)  :: HTEXT     ! General header
!
!*       0.2   Local variables
!
INTEGER      :: INTERVAL
INTEGER      :: IKE,IKB
INTEGER      :: IK
INTEGER,SAVE :: ICOL, ISTYL, IERR
INTEGER,SAVE :: I1D, NSUP
INTEGER      :: ID, IND1, ILENC
INTEGER      :: IKL, IKH, JB
INTEGER,SAVE :: ISUIT, ISUI, INDISTM, ISTOK
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: ISTM
!
REAL         :: Z1, Z2, ZX1, ZX2
REAL,SAVE    :: ZLWSC, ZSC
REAL,SAVE    :: ZSCMIN
REAL,SAVE    :: ZHMIN, ZHMAX
REAL,SAVE    :: ZMNM, ZMXM
REAL         :: ZVL, ZVR, ZVB, ZVT, ZWL, ZWR, ZWB, ZWT
REAL         :: ZDEBX, ZECART1, ZECART2
REAL         :: ZXPOSTITT1, ZXYPOSTITT1
REAL         :: ZXPOSTITT2, ZXYPOSTITT2
REAL         :: ZXPOSTITT3, ZXYPOSTITT3
REAL         :: ZXPOSTITB1, ZXYPOSTITB1
REAL         :: ZXPOSTITB2, ZXYPOSTITB2
REAL         :: ZXPOSTITB3, ZXYPOSTITB3
REAL,DIMENSION(2)  :: ZX(2), ZY(2)
REAL,DIMENSION(2)  :: ZXZERO(2), ZYZERO(2)
!
CHARACTER(LEN=80),SAVE :: YCARCOU
CHARACTER(LEN=80),SAVE :: YCAR
CHARACTER(LEN=40) :: YTEX
CHARACTER(LEN=100) :: YTEM
CHARACTER(LEN=100) :: YTITB3
CHARACTER(LEN=8)  :: YT
CHARACTER(LEN=10) :: FORMAX, FORMAY
!
!----------------------------------------------------------------------
!!!!!!!!!!! 110797
!ZHMIN=XHMIN; ZHMAX=XHMAX
!!!!!!!!!!! 110797
!CALL GQLN  (IERR,ISTYL)
print *,' +++pro1d entree   ISTYL ',ISTYL,' CVARNPV1 ',CVARNPV1(1:LEN(CVARNPV1))
YTEX(1:LEN(YTEX))=' '
!
!*      1.     DISPLAY ENVIRONMENT SETUP AND PROFILE DRAWING
!              ---------------------------------------------
!
!*      1.1    Array size calculation
!
!IK=SIZE(PPRO,1)
SELECT CASE(CTYPE)
  CASE('CART','MASK','SPXY')
    IKB=1+JPVEXT
    IKE=NKMAX+JPVEXT
    IK=(MIN(NKH,IKE)-MAX(NKL,IKB))+1
  CASE('SSOL','DRST','RSPL','RAPL')
    IKB=1
    IKE=SIZE(PPRO)
    IK=IKE
    IKL=NKL
    IKH=NKH
    NKL=IKB
    NKH=IKE
END SELECT
!
!WRITE(YCARCOU,1000)NIDEBCOU,NJDEBCOU,NLANGLE,KPRO
!1000 FORMAT(' Oblique section IDEB=',I2,' JDEB=',I2,' ANG.=',I3,  &
!' IPRO=',I3)
!
!*     1.2    Sets NCAR viewport and window
!
IF(LVPTPVUSER)THEN
  ZX1=XVPTPVL; ZX2=XVPTPVR; Z1=XVPTPVB; Z2=XVPTPVT
ELSE
  Z1=0.1
  Z2=0.9
!Z2=0.1+AMIN1(0.85,(XHMAX-XHMIN)/10000.)
  ZX1=0.13
  ZX2=0.9
ENDIF
!
IF(XHMAX > XHMIN)THEN
ELSE
SELECT CASE(CTYPE)
  CASE('CART','MASK','SPXY')
    XHMIN=0.
    XHMAX=XZZ(1,1,IKE)
  CASE('SSOL','DRST','RSPL')
    XHMIN=PPRO(1)
    XHMAX=PPRO(IK)
  CASE('RAPL')
    IF(PPRO(1) < PPRO(IK))THEN
      XHMIN=PPRO(1)
      XHMAX=PPRO(IK)
    ELSE
      XHMIN=PPRO(IK)
      XHMAX=PPRO(1)
    ENDIF
END SELECT
  
END IF
CALL SET(ZX1,ZX2,Z1,Z2,PTABMIN,PTABMAX,XHMIN,XHMAX,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD
!
!*    1.3     Actually draws the profile
!
CALL AGSETF('SET.',4.)    ! autograph uses the last SET values
CALL AGSETF('BAC.',4.)    ! no perimeter drawn
CALL AGSETF('FRA.',2.)    ! no frame advance
!!!!!Oct 99
CALL GQLN  (IERR,ISTYL)
!print *,' +++pro1d    ISTYL ',ISTYL
ISTYL=ISTYL+NGSLNP
IF(ISTYL == 1)CALL AGSETR('DAS/PA/1.',65535.)
IF(ISTYL == 2)CALL AGSETR('DAS/PA/1.',30583.)
IF(ISTYL == 3)CALL AGSETR('DAS/PA/1.',21845.)
IF(ISTYL == 4)CALL AGSETR('DAS/PA/1.',10023.)
IF(ISTYL == 5)CALL AGSETR('DAS/PA/1.',16191.)
IF(ISTYL == 6)CALL AGSETR('DAS/PA/1.',990.)
IF(ISTYL == 7)CALL AGSETR('DAS/PA/1.',3855.)
IF(ISTYL == 8)CALL AGSETR('DAS/PA/1.',24415.)
IF(ISTYL == 9)CALL AGSETR('DAS/PA/1.',13107.)
IF(ISTYL == 10)CALL AGSETR('DAS/PA/1.',63903.)
call gsln(1)
CALL AGSETR('DAS/SE.',1.)
if(nverbia >0)then
print *,' +++pro1d AV EZXY    ISTYL ',ISTYL
endif
!!!!!Oct 99
CALL EZXY(PTAB(MAX(NKL,IKB):MIN(NKH,IKE)), &
	       PPRO(MAX(NKL,IKB):MIN(NKH,IKE)),IK,0) ! calls autograph
!!!!!!!!!!!!!!JD Mars 2009 Ligne zero sur PV
IF(LINZEROPV)THEN
  IF(NSTYLINZEROPV == 1)CALL AGSETR('DAS/PA/1.',65535.)
  IF(NSTYLINZEROPV == 2)CALL AGSETR('DAS/PA/1.',30583.)
  IF(NSTYLINZEROPV == 3)CALL AGSETR('DAS/PA/1.',21845.)
  IF(NSTYLINZEROPV == 4)CALL AGSETR('DAS/PA/1.',10023.)
  IF(NSTYLINZEROPV == 5)CALL AGSETR('DAS/PA/1.',16191.)
  IF(NSTYLINZEROPV == 6)CALL AGSETR('DAS/PA/1.',990.)
  IF(NSTYLINZEROPV == 7)CALL AGSETR('DAS/PA/1.',3855.)
  IF(NSTYLINZEROPV == 8)CALL AGSETR('DAS/PA/1.',24415.)
  IF(NSTYLINZEROPV == 9)CALL AGSETR('DAS/PA/1.',13107.)
  IF(NSTYLINZEROPV == 10)CALL AGSETR('DAS/PA/1.',63903.)
CALL GSLN(NSTYLINZEROPV)
  ZXZERO(1:2)=0.
  ZYZERO(1)=XHMIN
  ZYZERO(2)=XHMAX
  CALL CURVED(ZXZERO,ZYZERO,2)
ENDIF
!!!!!!!!!!!!!!JD Mars 2009 Ligne zero sur PV
!!!!!Oct 99
CALL GSLN(ISTYL)
!!!!!Oct 99
IF(LSUPER)THEN
  CALL GQPLCI(IERR,ICOL)
! CALL GQLN  (IERR,ISTYL)
END IF
CALL GQLN  (IERR,ISTYL)
CALL GQLWSC(IERR,ZLWSC)
CALL GSLN(1)              ! solid line restored
CALL GSPLCI(1)
CALL GSTXCI(1)
CALL GSLWSC(1.)
SELECT CASE(CTYPE)
  CASE('SSOL','DRST','RSPL','RAPL')
    NKL=IKL
    NKH=IKH
END SELECT
!
!*    1.4     Prints the tick labels
!
IF(NPVITVYMJ /= 0)THEN
  INTERVAL=NPVITVYMJ
ELSE
IF(XHMAX-XHMIN < 2000.)THEN
  INTERVAL=5
ELSE
  INTERVAL=NINT((XHMAX-XHMIN)/1000.)
ENDIF
ENDIF
FORMAX='          '
IF(LFMTAXEX)THEN
  FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
ENDIF
  FORMAY='          '
IF(LFMTAXEY)THEN
  FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
ELSE
  FORMAY='(F6.0)'
ENDIF
FORMAY=ADJUSTL(FORMAY)
!print *,' FORMAX,FORMAY ',FORMAX
!print *,' FORMAX,FORMAY ',FORMAY

IF(.NOT.LSUPER)THEN       !00000000000000000000000000000000

  IF(LFMTAXEX .AND. LFMTAXEY)THEN      !Aout 2000
    CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) ! sets label format ...
!   CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0) ! sets label format ...
  ELSE IF(LFMTAXEX)THEN
    CALL LABMOD(FORMAX,'(F6.0)',0,0,NSZLBX,NSZLBY,0,0,0) ! sets label format ...
!   CALL LABMOD(FORMAX,'(F6.0)',0,0,10,10,0,0,0) ! sets label format ...
  ELSE
  IF(PTABMAX /= 0.)THEN
! IF(LOG10(ABS(PTABMAX)).GE.6. .OR. LOG10(ABS(PTABMAX)).LE. -1.)THEN
  IF(LOG10(ABS(PTABMAX)).GE.6. .OR. LOG10(ABS(PTABMAX)).LT. 0.)THEN
    CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) ! sets label format ...
!   CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,0,0) ! sets label format ...
!   CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0) ! sets label format ...
  ELSE
    IF(PTABMIN /= 0. .AND. (LOG10(ABS(PTABMIN)).GE.5. .OR. LOG10(ABS(PTABMIN)) &
      .LT. 0.))THEN
!     .LE. -1.))THEN
      CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
!     CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0)
    ELSE IF(ABS(PTABMAX-PTABMIN) <= 1.)THEN
!   ELSE IF(ABS(PTABMAX-PTABMIN) < 1.)THEN
      CALL LABMOD('(F8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
      CALL LABMOD('(F8.2)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.2)','(F6.0)',0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.2)','(F6.0)',8,6,10,10,0,0,0)
    ELSE
      IF(PTABMIN <0)THEN
      CALL LABMOD('(F9.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
!     CALL LABMOD('(F9.1)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.0)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.0)','(F6.0)',0,0,10,10,0,0,0)
      ELSE
      CALL LABMOD('(F8.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
!     CALL LABMOD('(F8.1)',FORMAY,0,0,10,10,0,0,0)
!     print *,' format F8.1  ************pro1d'
!     CALL LABMOD('(F8.1)','(F6.0)',0,0,10,10,0,0,0)
      ENDIF
    END IF
  END IF
  ELSE
!  PTABMAX = 0.
    IF(PTABMIN /= 0. .AND. (LOG10(ABS(PTABMIN)).GE.5. .OR. LOG10(ABS(PTABMIN)) &
      .LT. 0.))THEN
!     .LE. -1.))THEN
      CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
!     CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0)
    ELSE IF(ABS(PTABMAX-PTABMIN) <= 1.)THEN
!   ELSE IF(ABS(PTABMAX-PTABMIN) < 1.)THEN
      CALL LABMOD('(F8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
!     CALL LABMOD('(F8.2)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.2)','(F6.0)',0,0,10,10,0,0,0)
    ELSE
      CALL LABMOD('(F9.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) 
!     CALL LABMOD('(F9.1)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.0)',FORMAY,0,0,10,10,0,0,0)
!     CALL LABMOD('(F8.0)','(F6.0)',0,0,10,10,0,0,0)
    END IF
   
  ENDIF

  ENDIF                                !Aout 2000

  CALL GASETI('LTY',1)                       ! High quality perimeter font
! Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,0,0,5,0,0)  ! draws perimeter and labels
  ELSEIF(LNOLABELX .AND. .NOT. LNOLABELY)THEN
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,0,1,5,0,0)  ! draws perimeter and labels
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,1,0,5,0,0)  ! draws perimeter and labels
  ELSE
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,1,1,5,0,0)  ! draws perimeter and labels
  ENDIF
! Avril 2002
! CALL GRIDAL(KXDOT,1,INTERVAL,1,1,1,5,0,0)  ! draws perimeter and labels

ELSE                      !00000000000000000000000000000000

  CALL GASETI('LTY',1)                       ! High quality perimeter font
  SELECT CASE(NSUPER)
  CASE(1)
    NSUP=NSUPER
    ZSCMIN=999.
    ZMNM=PTABMIN
    ZMXM=PTABMAX

    IF(LFMTAXEX)THEN                  ! Aout 2000
      CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!     CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
    ELSE

    IF(PTABMAX /= 0.)THEN
      IF(LOG10(ABS(PTABMAX)).GE.6. .OR. LOG10(ABS(PTABMAX)).LT. 0.)THEN
!     IF(LOG10(ABS(PTABMAX)).GE.6. .OR. LOG10(ABS(PTABMAX)).LE. -1.)THEN
        CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(E8.2)','(F6.0)',8,6,10,10,0,3+(NSUPER-1)*15,0)
      ELSE
        IF(PTABMIN /= 0. .AND. (LOG10(ABS(PTABMIN)).GE.5. .OR. LOG10(ABS(PTABMIN)) &
          .LT. 0.))THEN
!         .LE. -1.))THEN
          CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0)
!         CALL LABMOD('(E8.2)','(F6.0)',8,6,10,10,0,0,0)
        ELSE IF(ABS(PTABMAX-PTABMIN) <= 1.)THEN
!       ELSE IF(ABS(PTABMAX-PTABMIN) < 1.)THEN
          CALL LABMOD('(F8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F8.2)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F8.2)','(F6.0)',0,0,10,10,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F8.2)','(F6.0)',8,6,10,10,0,3+(NSUPER-1)*15,0)
        ELSE
          CALL LABMOD('(F9.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F9.1)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F8.0)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F8.0)','(F6.0)',0,0,10,10,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD('(F8.0)','(F6.0)',8,6,10,10,0,3+(NSUPER-1)*15,0)
        END IF
      END IF
    ELSE
  !  PTABMAX = 0.
      IF(PTABMIN /= 0. .AND. (LOG10(ABS(PTABMIN)).GE.5. .OR. LOG10(ABS(PTABMIN)) &
        .LT. 0.))THEN
!       .LE. -1.))THEN
        CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0)
!       CALL LABMOD('(E8.2)','(F6.0)',8,6,10,10,0,0,0)
      ELSE IF(ABS(PTABMAX-PTABMIN) <= 1.)THEN
!     ELSE IF(ABS(PTABMAX-PTABMIN) < 1.)THEN
        CALL LABMOD('(F8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(F8.2)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(F8.2)','(F6.0)',0,0,10,10,0,0,0)
!       CALL LABMOD('(F8.2)','(F6.0)',8,6,10,10,0,0,0)
      ELSE
        CALL LABMOD('(F9.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(F9.1)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(F8.0)',FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
!       CALL LABMOD('(F8.0)','(F6.0)',0,0,10,10,0,0,0)
!       CALL LABMOD('(F8.0)','(F6.0)',8,6,10,10,0,0,0)
      END IF
     
    ENDIF

    ENDIF                            ! Aout 2000

!   CALL GRIDAL(KXDOT,1,INTERVAL,1,1,1,5,0,0)  ! draws perimeter and labels
! Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,0,0,5,0,0)  ! draws perimeter and labels
  ELSEIF(LNOLABELX .AND. .NOT. LNOLABELY)THEN
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,0,1,5,0,0)  ! draws perimeter and labels
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,1,0,5,0,0)  ! draws perimeter and labels
  ELSE
    CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,1,1,5,0,0)  ! draws perimeter and labels
  ENDIF
! Avril 2002
  CASE(2:)
    IF(PTABMIN == ZMNM .AND. PTABMAX == ZMXM)THEN
    ELSE
      NSUP=NSUP+1
      IF (NSUP > 3)THEN
        WRITE(NLUOUT,*)' ** PRO1D_FORDIACHRO NB DE SUPERPOSITIONS TROP ELEVE. IMPOSSIBILITE D''ECRIRE LES BORNES '
        WRITE(NLUOUT,*)'    DES VARIABLES INSCRITES A DROITE DU DESSIN'
        WRITE(NLUOUT,*)' ** IL SUFFIRAIT PEUT-ETRE DE METTRE EN TETE DES VAR. a SUPERPOSER CELLE DONT '
        WRITE(NLUOUT,*)'    LES BORNES ENGLOBENT LES LIMITES DES AUTRES'
      ELSE

        IF(LFMTAXEX)THEN                  ! Aout 2000
          CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUPER-1)*15,0)
!         CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,3+(NSUPER-1)*15,0)
        ELSE

        IF(PTABMAX /= 0.)THEN
        IF(LOG10(ABS(PTABMAX)).GE.6. .OR. LOG10(ABS(PTABMAX)).LT. 0.)THEN
!       IF(LOG10(ABS(PTABMAX)).GE.6. .OR. LOG10(ABS(PTABMAX)).LE. -1.)THEN
          CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!         CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!         CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,3+(NSUP-1)*15,0)
!         CALL LABMOD('(E8.2)','(F6.0)',8,6,10,10,0,3+(NSUP-1)*15,0)
        ELSE
          IF(PTABMIN /= 0. .AND. (LOG10(ABS(PTABMIN)).GE.5. .OR. LOG10(ABS(PTABMIN)) &
            .LT. 0.))THEN
!           .LE. -1.))THEN
            CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0)
!           CALL LABMOD('(E8.2)','(F6.0)',8,6,10,10,0,0,0)
          ELSE IF(ABS(PTABMAX-PTABMIN) <= 1.)THEN
!         ELSE IF(ABS(PTABMAX-PTABMIN) < 1.)THEN
            CALL LABMOD('(F8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.2)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.2)','(F6.0)',0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.2)','(F6.0)',8,6,10,10,0,3+(NSUP-1)*15,0)
          ELSE
            CALL LABMOD('(F9.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F9.1)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           print *,' format f9.1 **********pro1d'
!           CALL LABMOD('(F8.0)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.0)','(F6.0)',0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.0)','(F6.0)',8,6,10,10,0,3+(NSUP-1)*15,0)
          END IF
        END IF
        ELSE
      !  PTABMAX = 0.
          IF(PTABMIN /= 0. .AND. (LOG10(ABS(PTABMIN)).GE.5. .OR. LOG10(ABS(PTABMIN)) &
            .LT. 0.))THEN
!           .LE. -1.))THEN
            CALL LABMOD('(E8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(E8.2)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(E8.2)','(F6.0)',0,0,10,10,0,0,0)
!           CALL LABMOD('(E8.2)','(F6.0)',8,6,10,10,0,0,0)
          ELSE IF(ABS(PTABMAX-PTABMIN) <= 1.)THEN
!         ELSE IF(ABS(PTABMAX-PTABMIN) < 1.)THEN
            CALL LABMOD('(F8.2)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.2)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.2)','(F6.0)',0,0,10,10,0,0,0)
!           CALL LABMOD('(F8.2)','(F6.0)',8,6,10,10,0,0,0)
          ELSE
            CALL LABMOD('(F9.1)',FORMAY,0,0,NSZLBX,NSZLBY,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F9.1)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           CALL LABMOD('(F8.0)',FORMAY,0,0,10,10,0,3+(NSUP-1)*15,0)
!           print *,' format f9.1 **********pro1d'
!           CALL LABMOD('(F8.0)','(F6.0)',0,0,10,10,0,0,0)
!           CALL LABMOD('(F8.0)','(F6.0)',8,6,10,10,0,0,0)
          END IF
         
        ENDIF

        ENDIF                             ! Aout 2000

!       CALL GRIDAL(KXDOT,1,INTERVAL,1,1,1,5,0,0)  ! draws perimeter and labels
! Avril 2002
        IF(LNOLABELX .AND. LNOLABELY)THEN
          CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,0,0,5,0,0)  ! draws perimeter and labels
        ELSEIF(LNOLABELX .AND. .NOT. LNOLABELY)THEN
          CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,0,1,5,0,0)  ! draws perimeter and labels
        ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
          CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,1,0,5,0,0)  ! draws perimeter and labels
        ELSE
          CALL GRIDAL(KXDOT,NPVITVXMN,INTERVAL,NPVITVYMN,1,1,5,0,0)  ! draws perimeter and labels
        ENDIF
! Avril 2002
        CALL GSPLCI(ICOL)
! Oct 99
!       CALL GSLN(ISTYL)
        CALL GSLN(1)
        CALL AGSETR('DAS/SE.',1.)
! Oct 99
        CALL GSLWSC(ZLWSC)
        ZX(1)=PTABMIN+((MIN(0.06+ZX2,.96)-ZX1)*ABS(PTABMAX-PTABMIN)/(ZX2-ZX1))
!       ZX(1)=PTABMIN+((0.96-ZX1)*ABS(PTABMAX-PTABMIN)/(ZX2-ZX1))
        ZX(2)=PTABMIN+((MIN(0.10+ZX2,1.)-ZX1)*ABS(PTABMAX-PTABMIN)/(ZX2-ZX1))
!       ZX(2)=PTABMIN+((1.00-ZX1)*ABS(PTABMAX-PTABMIN)/(ZX2-ZX1))
        ZY(1)=XHMIN-ABS(((XHMAX-XHMIN)*(3+10+(NSUP-1)*15))/((Z2-Z1)*1024.))
        ZY(2)=ZY(1)
        CALL GSCLIP(0)
! Oct 99
    if(nverbia > 0)then
      print *,' AV CURVED'
    endif
        CALL CURVED(ZX,ZY,2)
!       CALL GPL(2,ZX,ZY)
! Semble inutile
!       CALL GSLN(ISTYL)
! Semble inutile
! Oct 99
!       CALL GSCLIP(1)
        CALL GSLWSC(1.)
      END IF
        CALL GSPLCI(1)
        CALL GSLN(1)
        CALL GSLWSC(1.)
    END IF
  END SELECT
END IF                    !00000000000000000000000000000000
CALL GSCLIP(0)                             ! suppress clipping
!CALL PLCHHQ((PTABMIN-(PTABMAX-PTABMIN)/(ZX2-ZX1)*ZX1),XHMIN+(XHMAX-XHMIN)/2.,  &
!'ALTITUDE',.012,0.,-1.)
!CALL PLCHHQ((PTABMIN-(PTABMAX-PTABMIN)/(ZX2-ZX1)*ZX1),XHMIN+(XHMAX-XHMIN)/2.4, &
!'(M)',.012,0.,-1.)
!
!*    1.5     Headers printing with pretty font,
!*            and possible overlay
!
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)           
!IF(LFACTIMP)THEN
! CALL FACTIMP
!ENDIF
! 
! Page headers
!
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN

IF(LFACTIMP )THEN
  CALL FACTIMP
ENDIF
  ZXPOSTITB1=.002
  ZXYPOSTITB1=.005
  IF(XPOSTITB1 /= 0.)THEN
    ZXPOSTITB1=XPOSTITB1
  ENDIF
  IF(XYPOSTITB1 /= 0.)THEN
    ZXYPOSTITB1=XYPOSTITB1
  ENDIF

  CALL RESOLV_TIT('CTITB1',HLEGEND)
  IF(HLEGEND /= ' ')THEN
    IF(XSZTITB1 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HLEGEND,XSZTITB1,0.,-1.)
!     CALL PLCHHQ(0.002,0.005,HLEGEND,XSZTITB1,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HLEGEND,.007,0.,-1.)
!     CALL PLCHHQ(0.002,0.005,HLEGEND,.007,0.,-1.)
    ENDIF
  ENDIF
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
      CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,XSZTITB2,0.,-1.)
!     CALL PLCHHQ(0.002,0.025,CLEGEND2,XSZTITB2,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,.007,0.,-1.)
!     CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
    ENDIF
  ENDIF
! Titres en X
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXL',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXL',YTEM)
  CALL PLCHHQ(ZVL,ZVB-MIN(ZVB/2.,.05),YTEM,XSZTITXL,0.,-1.)
! CALL PLCHHQ(ZVL,ZVB/2.,YTEM,.008,0.,-1.)
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXM',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXM',YTEM)
  CALL PLCHHQ((ZVL+ZVR)/2.,ZVB-MIN(ZVB/2.,.05),YTEM(1:LEN_TRIM(YTEM)),XSZTITXM,0.,0.)
! CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXR',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXR',YTEM)
    CALL PLCHHQ(ZVR,ZVB-MIN(ZVB/2.,.05),YTEM,XSZTITXR,0.,+1.)
! CALL PLCHHQ(ZVR-ZVB/2.,ZVB/2.,YTEM,.008,0.,-1.)
  ENDIF
    if(nverbia > 0)then
      print *,' ***pro1d 627'
    endif
! Titres en Y
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  YTEM='ALTITUDE;(M)'
  CALL RESOLV_TITY('CTITYM',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYB',ZVL,ZVR,ZVB,ZVT,YTEM)
! Titres  TOP
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
!     CALL PLCHHQ(0.002,0.93,YTEM,XSZTITT3,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,.008,0.,-1.)
!     CALL PLCHHQ(0.002,0.93,YTEM,.008,0.,-1.)
    ENDIF
  ENDIF

ENDIF
!
! Profile location
!
IF(L1DT)THEN
      if(nverbia > 0)then
	print *,' ***pro1d  L1DT=T'
      endif
  SELECT CASE(CTYPE)
    CASE('CART','MASK')
      IF(NIL == 1 .OR. NJL == 1)THEN
      WRITE(YCARCOU,1002)
      ELSE
      WRITE(YCARCOU,1012)NIL,NJL
      ENDIF
        YCAR(1:LEN(YCAR))=' '
    if(nverbia > 0)then
      print *,' ***pro1d 675'
    endif

    CASE('SSOL')
      IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN
        YCARCOU(1:LEN(YCARCOU))=' '
        YCAR(1:LEN(YCAR))=' '
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
        INDISTM=1
        ISTM(INDISTM)=NLOOPN
      ELSE IF(LSUPER .AND. NSUPER > 1)THEN
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
	ENDIF
      ENDIF
    CASE DEFAULT
      if(nverbia > 0)then
	print *,' ***pro1d  CASE DEFAULT'
      endif
      IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 1))THEN
        YCARCOU(1:LEN(YCARCOU))=' '
        YCAR(1:LEN(YCAR))=' '
        YCARCOU(1:4)=CTYPE
        YCARCOU(5:7)=' N.'
        WRITE(YCARCOU(8:10),'(I3)')NLOOPN
      if(nverbia > 0)then
	print *,' ***pro1d  YCARCOU',YCARCOU(1:LEN_TRIM(YCARCOU))
      endif
        ISUIT=11
        IF(ALLOCATED(ISTM))THEN
  	DEALLOCATE(ISTM)
        ENDIF
  	ALLOCATE(ISTM(NSUPERDIA))
        INDISTM=1
        ISTM(INDISTM)=NLOOPN
      if(nverbia > 0)then
	print *,' ***pro1d  YCARCOU',YCARCOU(1:LEN_TRIM(YCARCOU))
      endif
      ELSE IF(LSUPER .AND. NSUPER > 1)THEN
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
	  WRITE(YCARCOU(ISUIT:ISUIT+4),'(I5)')NLOOPN
	  ISUIT=ISUIT+5
	ENDIF
      ENDIF
   END SELECT
ELSE
  YCAR(1:LEN(YCAR))=' '
  IF(XIDEBCOU /= -999.)THEN
      IF(LDEFCV2CC)THEN           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IF(LDEFCV2IND)THEN
	  WRITE(YCARCOU,1018)NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
	  WRITE(YCAR,1006)KPRO,XIPROFV,XJPROFV
	ELSE IF(LDEFCV2LL)THEN
	  WRITE(YCARCOU,1019)XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
	  WRITE(YCAR,1006)KPRO,XIPROFV,XJPROFV
	ELSE
	  WRITE(YCARCOU,1020)XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
	  WRITE(YCAR,1006)KPRO,XIPROFV,XJPROFV
	ENDIF
      ELSE                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF(XIDEBCOU < 99999.)THEN
        IF(XJDEBCOU < 99999.)THEN
          WRITE(YCARCOU,1001)XIDEBCOU,XJDEBCOU,NLANGLE,KPRO
	  IF(.NOT.LCARTESIAN)THEN
	    WRITE(YCAR,1006)KPRO,XIPROFV,XJPROFV
	  ENDIF
        ELSE
          WRITE(YCARCOU,1003)XIDEBCOU,XJDEBCOU,NLANGLE,KPRO
	  IF(.NOT.LCARTESIAN)THEN
	    WRITE(YCAR,1006)KPRO,XIPROFV,XJPROFV
	  ENDIF
        END IF
      ELSE
        IF(XJDEBCOU < 99999.)THEN
          WRITE(YCARCOU,1004)XIDEBCOU,XJDEBCOU,NLANGLE,KPRO
        ELSE
          WRITE(YCARCOU,1005)XIDEBCOU,XJDEBCOU,NLANGLE,KPRO
        END IF
      END IF
      ENDIF                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ELSE
      if(nverbia > 0)then
	print *,' ***pro1d  AV YCARCOU',YCARCOU(1:LEN_TRIM(YCARCOU))
      endif
    WRITE(YCARCOU,1000)NIDEBCOU,NJDEBCOU,NLANGLE,KPRO
  END IF
END IF
    if(nverbia > 0)then
      print *,' ***pro1d 815'
    endif

!IF(L1DT)THEN
IF(L1DT .AND. NIL == 1 .AND. NIH == 1 .AND. NJL == 1 .AND. NJH == 1)THEN
      if(nverbia > 0)then
	print *,' ***pro1d  L1DT AV PCSETI'
      endif
CALL PCSETI('BF',1)               ! Fills a box around characters 
CALL PCSETR('BL',2.)              ! heavy line plotted
CALL PCSETR('BM',.3)              ! sets a box margin
CALL PCSETI('BC(1)',1)            ! sets box color for prints
ENDIF

ZXPOSTITT1=.002
ZXYPOSTITT1=.98
IF(XPOSTITT1 /= 0.)THEN
  ZXPOSTITT1=XPOSTITT1
ENDIF
IF(XYPOSTITT1 /= 0.)THEN
  ZXYPOSTITT1=XYPOSTITT1
ENDIF

ZXPOSTITT2=.002
ZXYPOSTITT2=.95
IF(XPOSTITT2 /= 0.)THEN
  ZXPOSTITT2=XPOSTITT2
ENDIF
IF(XYPOSTITT2 /= 0.)THEN
  ZXYPOSTITT2=XYPOSTITT2
ENDIF


IF(.NOT.LSUPER)THEN
      if(nverbia > 0)then
	print *,' ***pro1d AV RESOLV_TIT(CTITT1,YCARCOU)'
      endif
  CALL RESOLV_TIT('CTITT1',YCARCOU)
  IF(XSZTITT1 /= 0.)THEN
    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!   CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
  ELSE
    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.012,0.,-1.)
!   CALL PLCHHQ(0.002,0.98,YCARCOU,.012,0.,-1.)
  ENDIF
  IF(YCAR /= ' ')THEN
    IF(XSZTITT2 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,XSZTITT2,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YCAR,XSZTITT2,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,.012,0.,-1.)
!     CALL PLCHHQ(0.002,0.95,YCAR,.012,0.,-1.)
    ENDIF
  ENDIF
ELSE
  
  SELECT CASE(CTYPE)

    CASE('CART','MASK')

      IF(NSUPER == 1)THEN
        I1D=2
        IF(L1DT)I1D=1
        CALL RESOLV_TIT('CTITT1',YCARCOU)
        IF(XSZTITT1 /= 0.)THEN
          CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!         CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
	ELSE
          CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.012,0.,-1.)
!         CALL PLCHHQ(0.002,0.98,YCARCOU,.012,0.,-1.)
	ENDIF
    if(nverbia > 0)then
      print *,' ***pro1d 887'
    endif
	IF(.NOT.L1DT)THEN
! Mars 2000
              CALL RESOLV_TIT('CTITT2',YCAR)
          IF(XSZTITT2 /= 0.)THEN
            CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,XSZTITT2,0.,-1.)
!           CALL PLCHHQ(0.002,0.95,YCAR,XSZTITT2,0.,-1.)
	  ELSE
            CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,.012,0.,-1.)
!           CALL PLCHHQ(0.002,0.95,YCAR,.012,0.,-1.)
	  ENDIF
	ENDIF
      ELSE IF(NSUPER >= 2)THEN
        SELECT CASE(I1D)
          CASE(1)
            IF(.NOT.L1DT)THEN
    ! Titres  TOP
              CALL RESOLV_TIT('CTITT2',YCARCOU)
              ZXPOSTITT2=.002
              ZXYPOSTITT2=.92
              IF(XPOSTITT2 /= 0.)THEN
                ZXPOSTITT2=XPOSTITT2
              ENDIF
              IF(XYPOSTITT2 /= 0.)THEN
                ZXYPOSTITT2=XYPOSTITT2
              ENDIF
              IF(XSZTITT2 /= 0.)THEN
                CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCARCOU,XSZTITT2,0.,-1.)
!               CALL PLCHHQ(0.002,0.92,YCARCOU,XSZTITT2,0.,-1.)
	      ELSE
                CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCARCOU,.012,0.,-1.)
!               CALL PLCHHQ(0.002,0.92,YCARCOU,.012,0.,-1.)
	      ENDIF
              I1D=3
            END IF
          CASE(2)
            IF(L1DT)THEN
    ! Titres  TOP
      YTEM(1:LEN(YTEM))=' '
              CALL RESOLV_TIT('CTITT2',YTEM)
              IF(XSZTITT2 /= 0.)THEN
                CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCARCOU,XSZTITT2,0.,-1.)
!               CALL PLCHHQ(0.002,0.95,YCARCOU,XSZTITT2,0.,-1.)
	      ELSE
                CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCARCOU,.012,0.,-1.)
!               CALL PLCHHQ(0.002,0.95,YCARCOU,.012,0.,-1.)
              ENDIF
              I1D=3
            END IF
          CASE(3)
        END SELECT
      END IF

    CASE DEFAULT

	IF(NSUPER == NSUPERDIA)THEN
          CALL RESOLV_TIT('CTITT1',YCARCOU)
	  IF(YCARCOU /= ' ')THEN
            IF(XSZTITT1 /= 0.)THEN
              CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,XSZTITT1,0.,-1.)
!             CALL PLCHHQ(0.002,0.98,YCARCOU,XSZTITT1,0.,-1.)
	    ELSE
              CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU,.010,0.,-1.)
!             CALL PLCHHQ(0.002,0.98,YCARCOU,.010,0.,-1.)
	    ENDIF
	  ENDIF
          CALL RESOLV_TIT('CTITT2',YCAR)
	  IF(YCAR /= ' ')THEN
            IF(XSZTITT2 /= 0.)THEN
              CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,XSZTITT2,0.,-1.)
!             CALL PLCHHQ(0.002,0.95,YCAR,XSZTITT2,0.,-1.)
	    ELSE
              CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YCAR,.010,0.,-1.)
!             CALL PLCHHQ(0.002,0.95,YCAR,.010,0.,-1.)
	    ENDIF
	  ENDIF
	  IF(ALLOCATED(ISTM))THEN
	    DEALLOCATE(ISTM)
	  ENDIF
	ENDIF
  END SELECT
    if(nverbia > 0)then
      print *,' ***pro1d 970'
    endif

END IF
CALL GSFAIS(0)      
CALL PCSETI('BF',0)  ! desactivates text outline option
!
! Variable names
!
ILENC=LEN_TRIM(CTIMEC)
IF(ILENC < LEN(CTIMEC))THEN
  IF(CTIMEC(ILENC:ILENC) == '.')THEN
    CTIMEC(ILENC:ILENC)='s'
  ELSE
    ILENC=ILENC+1
    CTIMEC(ILENC:ILENC)='s'
  ENDIF
ENDIF
YT(1:LEN(YT))=' '
IND1=INDEX(CTIMEC,'=')
YT=CTIMEC(IND1+1:ILENC)
YT=ADJUSTL(YT)
ZXPOSTITB3=.75
ZXYPOSTITB3=.025
IF(XPOSTITB3 /= 0.)THEN
  ZXPOSTITB3=XPOSTITB3
ENDIF
IF(XYPOSTITB3 /= 0.)THEN
  ZXYPOSTITB3=XYPOSTITB3
ENDIF
YTEM(1:LEN(YTEM))=' '
if(nverbia > 0)then
print *,' **pro1d CTITB3 CTITB3MEM ',CTITB3, CTITB3MEM
endif
!!!!!!!!!!!!!!!!=================================================
IF(.NOT.LSUPER)THEN
      if(nverbia > 0)then
         print *,' ***pro1d AV CALL PLCHHQ(0.75,0.007,HTEXT,.011,0.,-1.) '
      endif
  CALL PLCHHQ(0.75,0.007,HTEXT,.011,0.,-1.)
      if(nverbia > 0)then
	print *,' ***pro1d AP CALL PLCHHQ(0.75,0.007,HTEXT,.011,0.,-1.) '
      endif

!! nov 2001
  IF(.NOT.LTITDEFM)THEN
    YTITB3=' '
    YTITB3=CTITB3
    CTITB3=' '
    CTITB3=CTITB3MEM
    CTITB3=ADJUSTL(CTITB3)
if(nverbia > 0)then
print *,' **pro1d CTITB3 CTITB3MEM ',CTITB3, CTITB3MEM
endif
    CALL RESOLV_TIT('CTITB3',YTEM)
!   CTITB3=YTITB3
  ELSE
!! nov 2001
    CALL RESOLV_TIT('CTITB3',YTEM)
  ENDIF

if(nverbia > 0)then
  print *,' YTEM++++++++ ',YTEM,' CTITB3 ',CTITB3
endif

  IF(LTITDEFM)THEN
! ELSE
!! Nov 2001
    IF(XSZTITB3 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTIMEC,XSZTITB3,0.,-1.)
!     CALL PLCHHQ(0.75,0.025,CTIMEC,XSZTITB3,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,CTIMEC,.011,0.,-1.)
!     CALL PLCHHQ(0.75,0.025,CTIMEC,.011,0.,-1.)
    ENDIF

  ELSEIF(YTEM /= ' ')THEN

    IF(XSZTITB3 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM(1:LEN_TRIM(YTEM)),XSZTITB3,0.,-1.)
!     CALL PLCHHQ(0.75,0.025,YTEM(1:LEN_TRIM(YTEM)),XSZTITB3,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM(1:LEN_TRIM(YTEM)),.011,0.,-1.)
!     CALL PLCHHQ(0.75,0.025,YTEM(1:LEN_TRIM(YTEM)),.011,0.,-1.)
    ENDIF
!! Nov 2001
  ENDIF
  IF(.NOT.LTITDEFM)THEN
    CTITB3=YTITB3
  ENDIF

!!!!!!!!!!!!!!!!=================================================
ELSE
!!!!!!!!!!!!!!!!=================================================

!! nov 2001
  IF(.NOT.LTITDEFM)THEN
    YTITB3=' '
    YTITB3=CTITB3
    CTITB3=' '
    CTITB3=CTITB3MEM
    CTITB3=ADJUSTL(CTITB3)
    CALL RESOLV_TIT('CTITB3',YTEM)
!   CTITB3=YTITB3
  ELSE
!   CTITB3=' '
    CALL RESOLV_TIT('CTITB3',YTEM)
  ENDIF
if(nverbia > 0)then
  print *,' YTEM2++++++++ ',YTEM
endif
  IF(YTEM /= 'DEFAULT' .AND. YTEM /= ' ')THEN
    IF(XSZTITB3 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM(1:LEN_TRIM(YTEM)),XSZTITB3,0.,-1.)
!     CALL PLCHHQ(0.75,0.025,YTEM(1:LEN_TRIM(YTEM)),XSZTITB3,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM(1:LEN_TRIM(YTEM)),.011,0.,-1.)
      if(nverbia >0)then
        print *,' ***pro1d CTITB3*******',CTITB3
      endif
!     CALL PLCHHQ(0.75,0.025,YTEM(1:LEN_TRIM(YTEM)),.011,0.,-1.)
    ENDIF
  ENDIF
  IF(.NOT.LTITDEFM)THEN
    CTITB3=YTITB3
  ENDIF
  SELECT CASE(CTYPE)
    CASE('SSOL','DRST','RSPL','RAPL')
      WRITE(YTEX(1:4),'(I4)')NLOOPN
      YTEX(1+5:MIN(LEN(YTEX),LEN_TRIM(HTEXT)+5))=HTEXT(1:MIN(LEN(YTEX),LEN_TRIM(HTEXT)))
      YTEX=ADJUSTL(ADJUSTR(YTEX))
      if(nverbia > 0)then
       print *,' PRO1D**** YTEX LEN_TRIM(HTEXT) ',LEN_TRIM(HTEXT),' ',YTEX
      endif
    CASE DEFAULT
      YTEX(1:MIN(LEN(YTEX),LEN_TRIM(HTEXT)))=HTEXT(1:MIN(LEN(YTEX),LEN_TRIM(HTEXT)))
      YTEX=ADJUSTL(ADJUSTR(YTEX))
  END SELECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD
  IF(NSUPER >4)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD BOTTOM TITRES
!!!!!!!!!!!!!!!
    IF(LVARNPVUSER)THEN
      CALL GSPLCI(ICOL)
      CALL GSTXCI(ICOL)
        IF(XSZVARNPVBOT /=0.)THEN
          ZSCMIN=XSZVARNPVBOT
        ELSE
          ZSCMIN=.008
        ENDIF
      IF(NSUPER == 5)THEN
        IF(CVARNPV5 == 'WHITE' .OR. CVARNPV5 == 'white')THEN
!         CVARNPV5(1:LEN_TRIM(CVARNPV5))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV5 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV5
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 6)THEN
        IF(CVARNPV6 == 'WHITE' .OR. CVARNPV6 == 'white')THEN
!         CVARNPV6(1:LEN_TRIM(CVARNPV6))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV6 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV6
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 7)THEN
        IF(CVARNPV7 == 'WHITE' .OR. CVARNPV7 == 'white')THEN
!         CVARNPV7(1:LEN_TRIM(CVARNPV7))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV7 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV7
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 8)THEN
        IF(CVARNPV8 == 'WHITE' .OR. CVARNPV8 == 'white')THEN
!         CVARNPV8(1:LEN_TRIM(CVARNPV8))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV8 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV8
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 9)THEN
        IF(CVARNPV9 == 'WHITE' .OR. CVARNPV9 == 'white')THEN
!         CVARNPV9(1:LEN_TRIM(CVARNPV9))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV9 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV9
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 10)THEN
        IF(CVARNPV10 == 'WHITE' .OR. CVARNPV10 == 'white')THEN
!         CVARNPV10(1:LEN_TRIM(CVARNPV10))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV10 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV10
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 11)THEN
        IF(CVARNPV11 == 'WHITE' .OR. CVARNPV11 == 'white')THEN
!         CVARNPV11(1:LEN_TRIM(CVARNPV11))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV11 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV11
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 12)THEN
        IF(CVARNPV12 == 'WHITE' .OR. CVARNPV12 == 'white')THEN
!         CVARNPV12(1:LEN_TRIM(CVARNPV12))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV12 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV12
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 13)THEN
        IF(CVARNPV13 == 'WHITE' .OR. CVARNPV13 == 'white')THEN
!         CVARNPV13(1:LEN_TRIM(CVARNPV13))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV13 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV13
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 14)THEN
        IF(CVARNPV14 == 'WHITE' .OR. CVARNPV14 == 'white')THEN
!         CVARNPV14(1:LEN_TRIM(CVARNPV14))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV14 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV14
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 15)THEN
        IF(CVARNPV15 == 'WHITE' .OR. CVARNPV15 == 'white')THEN
!         CVARNPV15(1:LEN_TRIM(CVARNPV15))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV15 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV15
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ENDIF
      print *,' NSUPER YTEX ',NSUPER,YTEX
      IF(XPOSXVARNPV5BOT /= 0.)THEN
        IF(XPOSYVARNPV5BOT == 0.)THEN
          CALL PLCHHQ(XPOSXVARNPV5BOT,.005+(NSUPER-5)*.017,YTEX(1:LEN_TRIM(YTEX)),ZSCMIN,0.,-1.)
        ELSE
          CALL PLCHHQ(XPOSXVARNPV5BOT,XPOSYVARNPV5BOT+(NSUPER-5)*.017,YTEX(1:LEN_TRIM(YTEX)),ZSCMIN,0.,-1.)
        ENDIF
      ELSEIF(XPOSYVARNPV5BOT /= 0.)THEN
          CALL PLCHHQ(.75,XPOSYVARNPV5BOT+(NSUPER-5)*.017,YTEX(1:LEN_TRIM(YTEX)),ZSCMIN,0.,-1.)
      ELSE
        CALL PLCHHQ(.75,.005+(NSUPER-5)*.017,YTEX(1:LEN_TRIM(YTEX)),ZSCMIN,0.,-1.)
      ENDIF
    ELSE
!!!!!!!!!!!!!!!
      CALL GSLN(1)
      CALL GSPLCI(1)
      CALL GSTXCI(1)
      CALL GSLWSC(1.)
      if(nverbia >0)then
        print *,' YTEX BOTTOM ',YTEX(1:LEN_TRIM(YTEX))
        print *,' YT BOTTOM ',YT 
      endif
      IF(ZSCMIN /= 999.)THEN
        CALL PLCHHQ(.75,.005+(NSUPER-5)*.017,YTEX(1:LEN_TRIM(YTEX)),ZSCMIN,0.,-1.)
!       CALL PLCHHQ(.75,.005+(NSUPER-5)*.017,HTEXT,ZSCMIN,0.,-1.)
      ELSE
        CALL PLCHHQ(.75,.005+(NSUPER-5)*.017,YTEX(1:LEN_TRIM(YTEX)),.007,0.,-1.)
!       CALL PLCHHQ(.75,.005+(NSUPER-5)*.017,HTEXT,.007,0.,-1.)
      ENDIF
      CALL PLCHHQ(.62,.005+(NSUPER-5)*.017,YT,.007,0.,-1.)
!     CALL PLCHHQ(.60,.005+(NSUPER-5)*.017,YT,.007,0.,-1.)
      if(nverbia > 0)then
        print *,' ***pro1d 1065'
      endif
      CALL GSPLCI(ICOL)
      CALL GSTXCI(ICOL)
! Oct 99
!   CALL GSLN(ISTYL)
      CALL GSLN(1)
      CALL AGSETR('DAS/SE.',1.)
! Oct 99
      CALL GSLWSC(ZLWSC)
      ZX(1)=.69
!   ZX(1)=.67
      ZX(2)=ZX(1)+.03
      ZY(1)=0.005+(NSUPER-5)*.017
      ZY(2)=ZY(1)
! Oct 99
!     CALL GPL(2,ZX,ZY)
      CALL CURVED(ZX,ZY,2)
!!!!!!!!!!!!!!!
    ENDIF
!!!!!!!!!!!!!!!
! Oct 99
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD
  ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD TOP
    IF(LVARNPVUSER)THEN
        IF(XSZVARNPVTOP /=0.)THEN
          ZSC=XSZVARNPVTOP
        ELSE
          ZSC=.008
        ENDIF
      IF(NSUPER == 1)THEN
        IF(CVARNPV1 == 'WHITE' .OR. CVARNPV1 == 'white')THEN
!         CVARNPV1(1:LEN_TRIM(CVARNPV1))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV1 /= ' ')THEN
          print *,' ***pro1d CVARNPV1 ',CVARNPV1
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=ADJUSTL(CVARNPV1)
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 2)THEN
        IF(CVARNPV2 == 'WHITE' .OR. CVARNPV2 == 'white')THEN
!         CVARNPV2(1:LEN_TRIM(CVARNPV2))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
          print *,' NSUPER=2 YTEX ',YTEX
        ELSEIF(CVARNPV2 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV2
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 3)THEN
        IF(CVARNPV3 == 'WHITE' .OR. CVARNPV3 == 'white')THEN
!         CVARNPV3(1:LEN_TRIM(CVARNPV3))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV3 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV3
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ELSEIF(NSUPER == 4)THEN
        IF(CVARNPV4 == 'WHITE' .OR. CVARNPV4 == 'white')THEN
!         CVARNPV4(1:LEN_TRIM(CVARNPV4))=' '
          YTEX(1:LEN_TRIM(YTEX))=' '
        ELSEIF(CVARNPV4 /= ' ')THEN
          YTEX(1:LEN_TRIM(YTEX))=' '
          YTEX=CVARNPV4
          YTEX=ADJUSTL(YTEX)
        ENDIF
      ENDIF
      CALL GSPLCI(ICOL)
      CALL GSTXCI(ICOL)
      IF(XPOSXVARNPV1TOP /= 0.)THEN
        ZDEBX=XPOSXVARNPV1TOP
      ELSE
        ZDEBX=ZVL
      ENDIF
      IF(XPOSYVARNPV1TOP /= 0.)THEN
        ZECART2=XPOSYVARNPV1TOP-ZVT
      ELSE
        ZECART2=.02
      ENDIF
      print *,' pro1d ZSC ',ZSC,' YTEX ',YTEX(1:LEN_TRIM(YTEX)),' YT ',YT
!     STOP
      CALL PLCHHQ(ZDEBX+(NSUPER-1)*.21,ZVT+ZECART2,YTEX(1:LEN_TRIM(YTEX)),ZSC,0.,-1.)
      IF(YTEX == ' ')THEN
      ELSE
      CALL PLCHHQ(ZDEBX+(NSUPER-1)*.21,ZVT+ZECART2+.02,YT,ZSC,0.,-1.)
      ENDIF
!     CALL PLCHHQ(ZDEBX,.95,YT,ZSC,0.,-1.)
    ELSE
!!!!!!!!!!!!!!
      CALL GSLN(1)
      CALL GSPLCI(1)
      CALL GSTXCI(1)
      CALL GSLWSC(1.)
      ZSC=.007
      IF(LEN_TRIM(HTEXT) >25)THEN
        ZSC=.006
!     ZSC=.005
      ELSE IF(LEN_TRIM(HTEXT) >20)THEN
        ZSC=.007
      ENDIF
      IF(NSUPERDIA > 3)THEN
        ZDEBX=.1
      ELSE
        ZDEBX=ZVL
      ENDIF
      IF(ZVT >= .9)THEN
        ZECART1=.01; ZECART2=.03
      ELSE
        ZECART1=.02; ZECART2=.04
      ENDIF
      ZSCMIN=MIN(ZSCMIN,ZSC)
      if(nverbia > 0)then
        print *,' ***pro1d YTEX TOP ',YTEX(1:LEN_TRIM(YTEX))
      endif
      CALL PLCHHQ(ZDEBX+(NSUPER-1)*.21,ZVT+ZECART2,YTEX(1:LEN_TRIM(YTEX)),ZSC,0.,-1.)
!     CALL PLCHHQ(ZVL+(NSUPER-1)*.21,ZVT+.03,YTEX(1:LEN_TRIM(YTEX)),ZSC,0.,-1.)
!     CALL PLCHHQ(ZVL+(NSUPER-1)*.21,ZVT+.03,HTEXT,.007,0.,-1.)
      CALL PLCHHQ(ZDEBX+(NSUPER-1)*.21,ZVT+ZECART1,YT,.006,0.,-1.)
      if(nverbia > 0)then
        print *,' ***pro1d 1113'
      endif
!     CALL PLCHHQ(ZVL+(NSUPER-1)*.21,ZVT+.01,YT,.006,0.,-1.)
      CALL GSPLCI(ICOL)
      CALL GSTXCI(ICOL)
! Oct 99
!     CALL GSLN(ISTYL)
      CALL GSLN(1)
      CALL AGSETR('DAS/SE.',1.)
! Oct 99
      CALL GSLWSC(ZLWSC)
      ZX(1)= ZDEBX+(NSUPER-1)*.21+(LEN_TRIM(YT)+2)*.010
!     ZX(1)= ZVL+(NSUPER-1)*.21+(LEN_TRIM(YT)+2)*.010
!     ZX(1)= ZVL+(NSUPER-1)*.21+(LEN_TRIM(HTEXT)+2)*.011
      ZX(2)=ZX(1)+.03
      ZY(1)=ZVT+ZECART1
!     ZY(1)=ZVT+.01
!     ZY(1)=ZVT+.02
!     ZY(1)=ZVT+.03
      ZY(2)=ZY(1)
! Oct 99
      CALL CURVED(ZX,ZY,2)
      if(nverbia > 0)then
      print *,' ***pro1d AP CURVED'
      endif
!!!!!!!!!!!!!!!!!!!
    ENDIF
!!!!!!!!!!!!!!!!!!!
!   CALL GPL(2,ZX,ZY)
! Oct 99
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!JDJD
  CALL GSLN  (1)
  CALL GSPLCI(1)
  CALL GSTXCI(1)
  CALL GSLWSC(1.)
END IF
!!!!!!!!!!!!!!!!=================================================
IF(.NOT.LSUPER .OR. NSUPER == 1)THEN
  IF(LDATFILE)CALL DATFILE_FORDIACHRO
ENDIF
1000 FORMAT('Vertical section IDEB=',I4,' JDEB=',I4,' ANG.=',I3,' IPRO=',I4)
1001 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' IPRO=',I4)
1002 FORMAT('Vertical profile (1D)')
1012 FORMAT('Vertical profile (1D) I=',I4,' J=',I4)
1003 FORMAT('Vertical section XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' IPRO=',I4)
1004 FORMAT('Vertical section XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' IPRO=',I4)
1005 FORMAT('Vertical section XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' IPRO=',I4)
1006 FORMAT('Vertical profile IPRO=',I4,'  --> LAT=',F10.5,' ,LON=',F10.5)
1018 FORMAT('Vertical section IND I,J (BEGIN)-(END)=(',I4,',',I4,')-(',I4,',',I4,')')
1019 FORMAT('Vertical section LAT,LON (BEGIN)-(END)=(',F5.1,',',F6.1,')-(',F5.1,',',F6.1,')')
1020 FORMAT('Vert. section CONF. COORD.(BEGIN)-(END)=(',F8.0,',',F8.0,')-(',F8.0,',',F8.0,')')
!
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL GSCLIP(1) ! Restores window clipping
!!!!!!!!!!! 110797
IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == NSUPERDIA))THEN
  XHMIN=ZHMIN; XHMAX=ZHMAX
ENDIF
!!!!!!!!!!! 110797
RETURN
!
!-----------------------------------------------------------------------------
!
!    2.       EXIT
!             ----
!
END SUBROUTINE  PRO1D_FORDIACHRO
