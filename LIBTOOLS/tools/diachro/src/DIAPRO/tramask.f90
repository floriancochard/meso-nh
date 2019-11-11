!     ######spl
      MODULE MODI_TRAMASK
!     ###################
!
INTERFACE
!
SUBROUTINE TRAMASK(PTEM,KLOOP)
INTEGER           :: KLOOP
REAL,DIMENSION(:,:,:)           :: PTEM
END SUBROUTINE TRAMASK
!
END INTERFACE
END MODULE MODI_TRAMASK
!     ######spl
      SUBROUTINE TRAMASK(PTEM,KLOOP)
!     ##############################
!
!!****  *TRAMASK* - 
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!     
!!    EXTERNAL
!!    --------
!!      CLSGKS    : closes NCAR and GKS graphics
!!      COMPCOORD : computes gridpoint locations, meshsizes and topography
!!                  for all the possible grids, and true altitude where
!!                  required.
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID    : Current MESO-NH grid indicator
!!
!!
!!      Module MODN_PARA: Defines NAM_DOMAIN_POS namelist (former PARA common)
!!         NLMAX            :  Number of points horizontally along
!!                             the vertical section
!!         Module MODD_DIM1 : contains dimensions of data arrays
!!              NKMAX      : z array dimension
!!
!!
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!         XXZ      : Gal-Chen z coordinate values for all the MESO-NH grids
!!
!!      Module MODD_GRID1      : declares grid variables (Model module)
!!         XZZ      : true gridpoint z altitude
!!
!!      Module MODD_SUPER   : defines plot overlay control variables
!!         LSUPER   : =.TRUE. --> plot overlay is active
!!                    =.FALSE. --> plot overlay is not active
!!
!-------------------------------------------------------------------------------
!
!*      0.     DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODN_PARA
USE MODN_NCAR
USE MODD_COORD  
USE MODD_GRID1  
USE MODD_SUPER  
USE MODD_RESOLVCAR
USE MODD_TIT
USE MODD_TITLE
USE MODD_CTL_AXES_AND_STYL
!
IMPLICIT NONE
!
!*     0.1    interface declarations
!
INTERFACE
      SUBROUTINE VALMNMX(PMIN,PMAX)
      REAL                :: PMIN, PMAX
      END SUBROUTINE VALMNMX
END INTERFACE
!
!*      0.1    Dummy arguments 
!
INTEGER    :: KLOOP
REAL,DIMENSION(:,:,:)  :: PTEM
!
!*      0.2    local variables 
!
!
INTEGER           :: ID
INTEGER           :: ILENT
INTEGER           :: ISTA, IER, INB, IWK
!
REAL,SAVE         :: ZMIN, ZMAX
REAL,SAVE         :: ZVL, ZVR, ZVB, ZVT, ZWL, ZWR, ZWB,ZWT
!
CHARACTER(LEN=20) :: YNOM
CHARACTER(LEN=40):: YTEXTE
CHARACTER(LEN=60):: YTEM
!
!-------------------------------------------------------------------------------
!
!*      1.     PRELIMINARY CALCULATIONS
!              ------------------------
!
YTEXTE(1:LEN(YTEXTE)) = ' '
ILENT=LEN_TRIM(CTITGAL)
YTEXTE=ADJUSTL(CTITGAL)
YTEXTE=ADJUSTL(YTEXTE)
!
!
!*      1.4  
!
!
YNOM=ADJUSTL(CGROUP)
IF(YNOM.EQ.'QUIT')THEN      
!
!*       1.5    End of job: EXIT
!
  CALL GQOPS(ISTA)
  CALL GQACWK(1,IER,INB,IWK)
  IF(ISTA >1 .AND. INB >1)THEN
    CALL GDAWK(2)
    CALL GCLWK(2)
  ENDIF
! CALL FRAME
  CALL NGPICT(1,1)
  CALL CLSGKS
  STOP
ENDIF
!
!*       1.6   Ooverlay control
!
IF(NSUPERDIA > 1)THEN
  LSUPER=.TRUE.
ELSE
  LSUPER=.FALSE.
ENDIF
IF(KLOOP == 1)NSUPER=0
!print *,' KLOOP NSUPER ',KLOOP,NSUPER
!
!
!*       1.8    Line width and color changes to differentiate the 
!*              successive plots in an overlay sequence 
!
CALL GSCLIP(1)
IF(LSUPER)THEN

  NSUPER=NSUPER+1
  IF(NSUPER == 1)CALL TABCOL_FORDIACHRO
  IF(LCOLINE)THEN
    CALL GSLN(1)
    CALL GSPLCI(NSUPER+1)
    CALL GSTXCI(NSUPER+1)
  ELSE
    CALL GSPLCI(1)
    CALL GSTXCI(1)
    SELECT CASE(NSUPER)
      CASE(:4)
	CALL GSLWSC(1.)
      CASE(5:8)
	CALL GSLWSC(2.)
      CASE(9:12)
	CALL GSLWSC(3.)
      CASE(13:16)
	CALL GSLWSC(4.)
      CASE DEFAULT
	CALL GSLWSC(1.)
    END SELECT
    CALL GSLN(MOD(NSUPER,4))
    IF(MOD(NSUPER,4) == 0)CALL GSLN(4)
  ENDIF

ELSE

  CALL GSLN(1)              ! Solid line if no overlay
  CALL GSPLCI(1)
  CALL GSTXCI(1)

END IF
!
IF(NSUPER <= 1)THEN
  CALL AGSETF('SET.',4.)
  CALL AGSETF('BAC.',4.)
  CALL AGSETF('FRA.',2.)
!print *,' AGSETF '
  ZMIN=MINVAL(PTEM)
  ZMAX=MAXVAL(PTEM)
  CALL VALMNMX(ZMIN,ZMAX)
  IF(ABS(ZMAX-ZMIN) <1.E-4)THEN
    ZMAX=ZMAX+1.
    ZMIN=ZMIN-1.
  ENDIF
! ZMIN=-.5; ZMAX=1.5
  ZWB=ZMIN; ZWT=ZMAX
ENDIF
!print *,' SIZE(PTEM) ',SIZE(PTEM,1),SIZE(PTEM,2),SIZE(PTEM,3)
IF(SIZE(PTEM,1) == 1)THEN
  ZWL=XXY(NJINF,NMGRID); ZWR=XXY(NJSUP,NMGRID)
  CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
  CALL EZXY(XXY(NJINF:NJSUP,NMGRID),PTEM(1,:,1),NJSUP-NJINF+1,0)
ELSE IF(SIZE(PTEM,2) == 1)THEN
  ZWL=XXX(NIINF,NMGRID); ZWR=XXX(NISUP,NMGRID)
! print *,' ZWL ZWR ',ZWL,ZWR
  CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
  CALL EZXY(XXX(NIINF:NISUP,NMGRID),PTEM(:,1,1),NISUP-NIINF+1,0)
ENDIF

CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
!print *,' ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
CALL GSCLIP(0)
IF(LSUPER)THEN
  IF(NSUPER < 4)THEN
    CALL FRSTPT(ZVR-(ZVR-ZVL)/4.+MAX(.18,ILENT*.009),.007+(NSUPER-1)*.017)
    CALL VECTOR(ZVR-(ZVR-ZVL)/4.+MAX(.18,ILENT*.009)+.03,.007+(NSUPER-1)*.017)
  ELSE
    CALL PLCHHQ(ZVL+(NSUPER-4)*.25,ZVT+.01,ADJUSTL(CTIMEC(8:15)),.007,0.,-1.)
    CALL FRSTPT(ZVL+(NSUPER-4)*.25+.08,ZVT+.01)
    CALL VECTOR(ZVL+(NSUPER-4)*.25+.08+.03,ZVT+.01)
  ENDIF
ENDIF

CALL GSPLCI(1)
CALL GSTXCI(1)
CALL GSLN(1)
CALL GSLWSC(1.)
IF(NSUPER <= 1)THEN
! ******************************************************************
  CALL FORMATXY(ZWL,ZWR,ZWB,ZWT)
  CALL GRIDAL(NMASKITVXMJ,NMASKITVXMN,NMASKITVYMJ,NMASKITVYMN,1,1,5,0,0)
! CALL GRIDAL(5,1,5,1,1,1,5,0,0)
ENDIF
CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,1)
IF(.NOT.LSUPER)THEN
  ILENT=ILENT+2
  YTEXTE(ILENT:ILENT+15-8+1)=CTIMEC(8:15)
  CALL PLCHHQ(MAX(ZVR,.99),.007,YTEXTE(1:ILENT+15-8+1),.011,0.,+1.)
ELSE
  IF(NSUPER < 4)THEN
    CALL PLCHHQ(ZVR-(ZVR-ZVL)/4.-.04,.007+(NSUPER-1)*.017,YTEXTE(1:ILENT),  &
    .009,0.,-1.)
    CALL PLCHHQ(ZVR-(ZVR-ZVL)/4.-.12,.007+(NSUPER-1)*.017,CTIMEC(8:15),  &
    .007,0.,-1.)
  ELSE
    CALL PLCHHQ(ZVL+(NSUPER-4)*.25,ZVT+.03,YTEXTE(1:ILENT),  &
    .009,0.,-1.)
  ENDIF
ENDIF

CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
IF(LFACTIMP)THEN
  CALL FACTIMP
ENDIF
! Titres en X
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXL',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXL',YTEM)
  CALL PLCHHQ(ZVL,ZVB-MIN(ZVB/3.,.05),YTEM,.008,0.,-1.)
! CALL PLCHHQ(ZVL,ZVB/3.,YTEM,.008,0.,-1.)
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITXM',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXM',YTEM)
  CALL PLCHHQ((ZVL+ZVR)/2.,ZVB-MIN(ZVB/2.,.05),YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
! CALL PLCHHQ((ZVL+ZVR)/2.,ZVB/2.,YTEM(1:LEN_TRIM(YTEM)),.008,0.,0.)
! CALL PLCHHQ((ZVL+ZVR)/2.-ZVB/3.,ZVB/2.,YTEM,.008,0.,-1.)
  ENDIF
  YTEM(1:LEN(YTEM))=' '
! YTEM='(Sec.)'
  CALL RESOLV_TIT('CTITXR',YTEM)
  IF(YTEM /= ' ' .AND. YTEM /= 'DEFAULT')THEN
    CALL RESOLV_TIT('CTITXR',YTEM)
  CALL PLCHHQ(ZVR-ZVB/2.,ZVB-MIN(ZVB/3.,.05),YTEM,.008,0.,-1.)
! CALL PLCHHQ(ZVR-ZVB/2.,ZVB/3.,YTEM,.008,0.,-1.)
  ENDIF
! Titres en Y
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYT',ZVL,ZVR,ZVB,ZVT,YTEM)
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYM',ZVL,ZVR,ZVB,ZVT,YTEM)
  IF(LCNSUM)THEN
    YTEM='SUM(.TRUE.=1)'
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TITY('CTITYB',ZVL,ZVR,ZVB,ZVT,YTEM)
! Titres  TOP
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITT3',YTEM)
  IF(CTITT3 /= ' ')THEN
    CALL PLCHHQ(0.002,0.93,YTEM,.008,0.,-1.)
  ENDIF
  YTEM(1:LEN(YTEM))=' '
  CALL RESOLV_TIT('CTITT2',YTEM)
  IF(CTITT2 /= ' ')THEN
    CALL PLCHHQ(0.002,0.95,YTEM,.008,0.,-1.)
  ENDIF
    YTEM(1:LEN(YTEM))=' '
    CALL RESOLV_TIT('CTITT1',YTEM)
    IF(CTITT1 /= ' ')THEN
      CALL PLCHHQ(0.002,0.98,YTEM,.012,0.,-1.)
    ENDIF
! Titres  BOTTOM
! YTEM(1:LEN(YTEM))=' '
! CALL RESOLV_TIT('CTITB3',YTEM)
! IF(CTITB3 /= ' ')THEN
!   CALL PLCHHQ(0.002,0.05,YTEM,.008,0.,-1.)
! ENDIF
! YTEM(1:LEN(YTEM))=' '
! CALL RESOLV_TIT('CTITB2',YTEM)
! IF(CTITB2 /= ' ')THEN
!   CALL PLCHHQ(0.002,0.025,YTEM,.007,0.,-1.)
! ENDIF
! YTEM(1:LEN(YTEM))=' '
! CALL RESOLV_TIT('CTITB1',YTEM)
! IF(CTITB1 /= ' ')THEN
!   CALL PLCHHQ(0.002,0.005,YTEM,.007,0.,-1.)
! ENDIF
! Titre N1 BOTTOM
  CALL RESOLV_TIT('CTITB1',CLEGEND)
  CALL PLCHHQ(0.002,0.005,CLEGEND,.007,0.,-1.)
  IF(LCNCUM .OR. LCNSUM)THEN
! Titre N3 BOTTOM
  CALL RESOLV_TIT('CTITB3',CTIMECS)
  CALL PLCHHQ(0.002,0.050,CTIMECS,.009,0.,-1.)
  ELSE
  IF(LMINUS .OR. LPLUS)THEN
    IF(.NOT.LTITDEFM .AND. CTITB3MEM /= 'DEFAULT' .AND. &
    CTITB3MEM /= 'default' .AND. CTITB3MEM /= 'DEFAUT' .AND. &
    CTITB3MEM /= 'defaut')THEN
      IF(CTITB3MEM /= ' ' .AND. CTITB3MEM /= 'WHITE' .AND. &
      CTITB3MEM /= 'white' .AND. CTITB3MEM /= 'BLANC' .AND. &
      CTITB3MEM /= 'blanc')THEN
        CALL PLCHHQ(0.002,0.050,CTITB3MEM(1:LEN_TRIM(CTITB3MEM)),.009,0.,-1.)
      ENDIF
    ELSE
! ******************** 200697 ***************
    CALL RESOLV_TIT('CTITB3',CTITB3)
    IF(CTITB3 /= ' ')THEN
      CALL PLCHHQ(0.002,0.050,CTITB3,.009,0.,-1.)
    ENDIF
    ENDIF
! ******************** 200697 ***************
  ELSE
    YTEM(1:LEN(YTEM))=' '
    YTEM=CTIMEC
    YTEM=ADJUSTL(YTEM)
    CALL RESOLV_TIT('CTITB3',YTEM)
!   CALL RESOLV_TIT('CTITB3',CTIMEC)
    IF(YTEM /= ' ')THEN
      CALL PLCHHQ(0.002,0.050,YTEM(1:LEN_TRIM(YTEM)),.009,0.,-1.)
!     CALL PLCHHQ(0.002,0.050,CTIMEC,.009,0.,-1.)
    ENDIF
  ENDIF
  ENDIF
! Titre N2 BOTTOM
  CALL RESOLV_TIT('CTITB2',CLEGEND2)
  IF(CLEGEND2 /= ' ')THEN
    CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
  ENDIF
IF(LDATFILE)CALL DATFILE_FORDIACHRO
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)
!
!
!----------------------------------------------------------------------------
!
!*       4.     EXIT
!               ----
!
END SUBROUTINE  TRAMASK
