!     #########################
      SUBROUTINE SUBSPXY(KLOOP)
!     #########################
!
!
!!
!!    PURPOSE
!!    -------
!
!     Traite les informations de type SPXY et envoyees sous forme
!     d'un vecteur de coefficients spectraux ou d'un plan
!     Partie retiree de OPER_PROCESS devenue trop volumineuse pour
!     la compilation
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist 
!!                         (former NCAR common)
!!
!!       NIOFFD     : Label normalisation (=0 none, =/=0 active)
!!       NULBLL     : Nb of contours between 2 labelled contours
!!       NIOFFM     : =0    --> message at picture bottom
!!                    =/= 0 --> no message
!!       NIOFFP     : Special point value detection
!!                    (=0 none, =/=0 active)
!!       NHI        : Extrema detection
!!                    (=0 --> H+L, <0 nothing)
!!       NINITA     : For streamlimes
!!       NINITB     : Not yet implemented
!!       NIGRNC     : Not yet implemented
!!       NDOT       : Line style
!!                    (=0|1|1023|65535 --> solid lines;
!!                    <0 --> solid lines for positive values and
!!                    dotted lines(ABS(NDOT))for negative values;
!!                    >0 --> dotted lines(ABS(NDOT)) )
!!       NIFDC      : Coastline data style (0 none, 1 NCAR, 2 IGN)
!!       NLPCAR     : Number of land-mark points to be plotted
!!       NIMNMX     : Contour selection option
!!                    (=-1 Min, max and inc. automatically set;
!!                    =0 Min, max automatically set; inc. given;
!!                    >0 Min, max, inc. given by user)
!!       NISKIP     : Rate for drawing velocity vectors
!!       CTYPHOR    : Horizontal cross-section type
!!                    (='K' --> model level section;
!!                     ='Z' --> constant-altitude section;
!!                     ='P' --> isobar section (planned)
!!                     ='T' --> isentrope section (planned)
!!       XSPVAL     : Special value
!!       XSIZEL     : Label size
!!       XLATCAR, XLONCAR :  Lat. and Long. of land-mark points
!!       LXY        : If =.TRUE., plots  a grid-mesh stencil background
!!       LXZ        : If =.TRUE., plots  a model-level stencil background 
!!
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist 
!!                          (former PARA common)
!!
!!       XIDEBCOU, XJDEBCOU : Origin of a vertical cross-section
!!                            in cartesian (or conformal) real values
!!       XHMIN      : Altitude of the vert. cross-section
!!                    bottom (in meters above sea-level)
!!       XHMAX      : Altitude of the vert. cross-section
!!                    top (in meters above sea-level)
!!
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       30/05/00
!!      Updated   PM   02/12/94
!!                VM   05/04/06 abscisse:2pi/j*OMEGA ET j*OMEGA
!!                             et Module (apres Phase dans les cas PHALO,PHAO)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_FILES_DIACHRO
USE MODD_ALLOC_FORDIACHRO
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODI_TRACEH_FORDIACHRO
USE MODD_TYPE_AND_LH
USE MODD_DIM1
USE MODD_TIT
USE MODD_GRID1
USE MODD_NMGRID
USE MODD_CVERT
USE MODD_MASK3D
USE MODD_TITLE
USE MODD_PARAMETERS
USE MODD_EXPERIM
USE MODN_NCAR    
USE MODN_PARA    
USE MODI_PRECOU_FORDIACHRO
USE MODI_TRACEV_FORDIACHRO
!USE MODI_VARFCT
!USE MODI_PVFCT
!USE MODI_CLOSF
USE MODI_LOADUNITIT
!USE MODI_TRAPRO_FORDIACHRO
USE MODD_COORD
USE MODD_CONF
USE MODD_SUPER
USE MODD_CST
USE MODD_PVT
USE MODD_DEFCV
USE MODE_GRIDPROJ

IMPLICIT NONE

INTERFACE
	      SUBROUTINE IMCOU_FORDIACHRO(PTABV,PINT,HLEGEND,HTEXT)
	      REAL,DIMENSION(:,:) :: PTABV
	      REAL                :: PINT
	      CHARACTER(LEN=*)    :: HTEXT, HLEGEND
	      END SUBROUTINE IMCOU_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE INTERP_FORDIACHRO(KLREF,KD,KF,PTAB,PTABREF)
      REAL,DIMENSION(:,:,:), INTENT(IN)         :: PTAB 
      REAL,DIMENSION(SIZE(PTAB,1),SIZE(PTAB,2)) :: PTABREF
      INTEGER :: KLREF
      END SUBROUTINE INTERP_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE IMAGE_FORDIACHRO(PTAB,KLREF,PTABINT,KNHI,KNDOT,HTEXTE)
      CHARACTER(LEN=*)   :: HTEXTE
      REAL                :: PTABINT
      REAL,DIMENSION(:,:) :: PTAB
      INTEGER :: KNHI, KNDOT, KLREF
      END SUBROUTINE IMAGE_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE TRAXY(PTEMX,PTEMY,KLOOP,HTITX,HTITY,PTIMED,PTIMEF)
      INTEGER    :: KLOOP
      REAL,DIMENSION(:)  :: PTEMX, PTEMY
      REAL               :: PTIMED, PTIMEF
      CHARACTER(LEN=*) :: HTITX, HTITY
      END SUBROUTINE TRAXY
END INTERFACE
COMMON/TEMV/XZWORKZ,XZZDS,NINX,NINY
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
#include "big.h"
REAL,DIMENSION(N2DVERTX,2500) :: XZWORKZ
!REAL,DIMENSION(1000,400) :: XZWORKZ
REAL,DIMENSION(N2DVERTX)     :: XZZDS
!REAL,DIMENSION(1000)     :: XZZDS
INTEGER                 :: NINX, NINY
LOGICAL                 :: LVERT, LHOR, LPT, LXABS
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER           :: KLOOP

!
!*       0.1   Local variables
!              ---------------
!
INTEGER   ::   J, JJ
INTEGER   ::   II, IJ, IK, IKU, IKB, IKE, IIU, IJU
INTEGER   ::   JLOOPP, JLOOPN, JLOOPT, JLOOPK, JLOOPZ
INTEGER   ::   IZ, IOMEGA, IEGAL
INTEGER   ::   ILENT, ILENU
INTEGER   ::   ISUP, IJSUP, IINF, IJINF
INTEGER   ::   IIB, IIE, IJB, IJE, IL
INTEGER   ::   ID
INTEGER,SAVE   ::   INUM, IRESP

REAL      ::   ZWL, ZWR, ZWB, ZWT
REAL      ::   ZVL, ZVR, ZVB, ZVT
REAL      ::   ZOMEGA
REAL      ::   ZMIN, ZMAX, ZZMIN,  ZZMAX 
REAL      ::   ZXPOSTITT1, ZXYPOSTITT1, ZXPOSTITT2, ZXYPOSTITT2
REAL      ::   ZXPOSTITT3, ZXYPOSTITT3
REAL      ::   ZXPOSTITB1, ZXYPOSTITB1, ZXPOSTITB2, ZXYPOSTITB2
REAL      ::   ZXPOSTITB3, ZXYPOSTITB3


REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZWORK3D
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: ZTEMCV
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: ZTEM1, ZTEMV
REAL,DIMENSION(:),ALLOCATABLE,SAVE     :: ZTEM1D, ZWORKZ, ZTEMLO
REAL,DIMENSION(:),ALLOCATABLE,SAVE     :: ZTE, ZWO

CHARACTER(LEN=40)  :: YTEXTE
!CHARACTER(LEN=LEN(CTITGAL))  :: YTITGAL
CHARACTER(LEN=4)   :: YC5S3='-5/3' 
CHARACTER(LEN=16)  :: YTITX, YTITY
CHARACTER(LEN=60)  :: YTEM
!CHARACTER(LEN=60)  :: YTEMP

LOGICAL            :: GOMEGAZOK, GOMEGAXOK, GOMEGAYOK
!------------------------------------------------------------------------------

!*****************************************************************************
!*****************************************************************************
!   CASE('SPXY')

	print *,'SUBSPSY entree NIMAX NJMAX NIL NJL NIH NJH NKL NKH ',NIMAX,NJMAX,NIL,NJL,NIH,NJH,NKL,NKH
	GOMEGAXOK=.FALSE.
	GOMEGAYOK=.FALSE.
	GOMEGAZOK=.FALSE.
	LSPX=.FALSE.
	LSPY=.FALSE.
	LSPZ=.FALSE.
	LSPSECTXY=.FALSE.
	LSPSECTXZ=.FALSE.
	LSPSECTYZ=.FALSE.
	IIB=1+JPHEXT; IIE=NIMAX+JPHEXT
	IJB=1+JPHEXT; IJE=NJMAX+JPHEXT
	IKU=NKMAX+2*JPVEXT
	IKB=1+JPVEXT; IKE=IKU-JPVEXT

	II=SIZE(XVAR,1)
	IJ=SIZE(XVAR,2)
	IK=SIZE(XVAR,3)  
	 
!!!!! UNIDIMENSIONNELS (Eventuellement sur plusieurs niveaux)
        IF(.NOT. LSPSECT)THEN  !iiiiiiiiiiiiiiiiiiiiiiii
!*************************************************************************
! PV // Z
!*************************************************************************
        IF(II == 1 .AND. IJ == 1 .AND. IK /= 1)THEN
          print *,' unidimensionnel1: II,IJ,IK=',II,IJ,IK

	  LSPZ=.TRUE.
	  ALLOCATE(ZTEM1D(SIZE(XVAR,3)),ZWORKZ(SIZE(XVAR,3)))

!+++++++++ Boucle processus +++++++++++++++++++++++++++++++++++

	  DO JLOOPP=1,NBPROCDIA(KLOOP)

	    NLOOPP=NPROCDIA(JLOOPP,KLOOP)
	    CALL LOADUNITIT(JLOOPP,KLOOP)
	    IOMEGA=INDEX(CCOMMENT(NLOOPP),'DOMEGAZ')
	    IF(IOMEGA == 0)THEN
	      IOMEGA=INDEX(CCOMMENT(NLOOPP),'Domegaz')
	      IF(IOMEGA == 0)THEN
		IOMEGA=INDEX(CCOMMENT(NLOOPP),'domegaz')
	      ENDIF
	    ENDIF
	    IF(IOMEGA == 0)THEN
	      PRINT *,' Delta OmegaZ (pulsation) non trouve dans le champ commentaire '
	      PRINT *,' On trace en indices de tableau'          
	      DO J=1,SIZE(ZTEM1D)
		ZTEM1D(J)=J
	      ENDDO
	      GOMEGAZOK=.FALSE.
              ZOMEGA=1.

	    ELSE

	      IEGAL=INDEX(CCOMMENT(NLOOPP)(IOMEGA:LEN_TRIM(CCOMMENT(NLOOPP))),'=')
	      READ(CCOMMENT(NLOOPP)(IOMEGA+IEGAL:LEN_TRIM(CCOMMENT(NLOOPP))),*)XOMEGAZ
    
              IF(XOMEGAZ == 0.)THEN
	        PRINT *,' Delta OmegaZ (pulsation)  =  0'
	        PRINT *,' On trace en indices de tableau'          
	        DO J=1,SIZE(ZTEM1D)
		  ZTEM1D(J)=J
	        ENDDO
	        GOMEGAZOK=.FALSE.
                ZOMEGA=1.

              ELSE

	        DO J=1,SIZE(ZTEM1D)
	  	  ZTEM1D(J)=J*XOMEGAZ
	        ENDDO
	        GOMEGAZOK=.TRUE.
                ZOMEGA=XOMEGAZ
	      ENDIF
	    ENDIF

	    IF(.NOT.LTINCRDIA(KLOOP,1))THEN         !TTTTTTTTTTTTTTTTTTTTTT

!+++++++++ Boucle temps +++++++++++++++++++++++++++++++++++

	      DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
		CALL RESOLV_TIMES(NLOOPT)
		WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NLOOPT,1)
! Partie reelle
	        ZWORKZ(:)=XVAR(1,1,:,NLOOPT,1,NLOOPP)
		ZMIN=MINVAL(ZTEM1D);ZMAX=MAXVAL(ZTEM1D)
		ZZMIN=MINVAL(ZWORKZ);ZZMAX=MAXVAL(ZWORKZ)
		CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		CALL AGSETF('FRA.',2.)
		CALL AGSETC('LAB/NAME.','B')
		CALL AGSETR('LAB/SU.',1.)
		CALL AGSETC('LAB/NAME.','L')
		CALL AGSETR('LAB/SU.',1.)
                CALL PCSETC('FC',':')
		IF(GOMEGAZOK)THEN                    !......................

!------
! _SPO_
!------
		  IF(LSPO)THEN

		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		    CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		    call plchhq(.9,.05,':PGL:X:PRL:Z:',.015,0.,1.)
		    call plchhq(0., .85,':PRU:(R):',.015,0.,-1.)

! Traitement de la partie imaginaire
		    IF(SIZE(XVAR,5) == 2)THEN
		      CALL FRAME
		      ZWORKZ(:)=XVAR(1,1,:,NLOOPT,2,NLOOPP)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PGL:X:PRL:Z:',.015,0.,1.)
		      call plchhq(0., .85,':PRU:(I):',.015,0.,-1.)
		    ENDIF

!---------------
! _OSPLO_ (/log)
!---------------
		  ELSE IF(LOSPLO)THEN

		    ZMIN=LOG10(ZMIN)
		    ZMAX=LOG10(ZMAX)
		    ZZMIN=MINVAL(ZWORKZ(:)*ZOMEGA)
		    ZZMAX=MAXVAL(ZWORKZ(:)*ZOMEGA)
		    CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		    CALL AGSETF('SET.',4.)
                    CALL EZXY(LOG10(ZTEM1D),ZWORKZ*ZOMEGA,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PGL:X:PRU:*(R):',.015,0.,-1.)

!------------------
! _LSPLO_ (Log/log)
!------------------
		  ELSE IF(LSPLO)THEN

		    CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		    CALL AGSETF('SET.',2.)
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PRU:Log(R):',.015,0.,-1.)
		    IF(LM5S3)THEN
		     CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		     CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		     CALL FRSTPT(3.,0.)
		     CALL VECTOR(0.,5.)
		     CALL GSCHH(.2)
		     CALL GTX(0.+.5,5.-.4,YC5S3)
		   ENDIF

!---------------
! _PHALO_ (/log)
!---------------
		  ELSE IF(LPHALO)THEN

                    ZMIN=LOG10(ZMIN)
		    ZMAX=LOG10(ZMAX)
		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
	              ZWORKZ(:)=ATAN2(XVAR(1,1,:,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		    ENDIF

!-------
! _PHAO_
!-------
		  ELSE IF(LPHAO)THEN

		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
	              ZWORKZ(:)=ATAN2(XVAR(1,1,:,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PGL:X:PRL:Z:',.015,0.,1.)
		      call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		    ENDIF
		  ENDIF

		ELSE                                !......................

		  IF(LSPO)THEN
		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		    CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		    call plchhq(.9,.05,':PRU:Ind(:PRL:Z:PRU:):',.015,0.,1.)
		    call plchhq(0., .85,':PRU:(R):',.015,0.,-1.)

! Traitement de la partie imaginaire
		    IF(SIZE(XVAR,5) == 2)THEN
		      CALL FRAME
		      ZWORKZ(:)=XVAR(1,1,:,NLOOPT,2,NLOOPP)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Ind(:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PRU:(I):',.015,0.,-1.)
		    ENDIF
		  ELSE
		  ENDIF

		ENDIF                               !......................
		CALL FRAME

!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!
                IF(LPRINT)THEN

                  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
                  IF(IRESP /= 0)THEN
                    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
                    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
                    PRINT  '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
                  ENDIF

                   WRITE(INUM,'(''SP  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
& CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
                   IF(SIZE(XVAR,5) < 2)THEN
                     IF(GOMEGAZOK)THEN
                       WRITE(INUM,'(''Partie reelle uniquement     DOMEGAZ= '',F6.1)')ZOMEGA
                     ELSE
                       WRITE(INUM,'(''Partie reelle uniquement     DOMEGAZ= '',F6.1,'' -> Trace en indices de grille'')')ZOMEGA
                     ENDIF
                   ELSE
                     IF(GOMEGAZOK)THEN
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAZ= '',F6.1)')ZOMEGA
                     ELSE
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAZ= '',F6.1,'' -> Trace en indices de grille'')')ZOMEGA
                     ENDIF
                   ENDIF
                   WRITE(INUM,'(''NBVAL en K '',i4 )')SIZE(ZTEM1D,1)

                   IF(SIZE(XVAR,5) < 2)THEN
                      WRITE(INUM,'(36(''*''))')
                      WRITE(INUM,'(10X,''X(K)'',9X,''Y(VAL.R)'')')
                      WRITE(INUM,'(36(''*''))')
                   DO J=1,SIZE(ZTEM1D,1)
                       WRITE(INUM,'(I4,2X,F8.1,(5X,E15.8))')J,ZTEM1D(J),ZWORKZ(J) 
                   ENDDO
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                      WRITE(INUM,'(10X,''X(=K)'',8X,''Y(VAL.R)'',11X,''Y(VAL.Im)'')')
                      WRITE(INUM,'(55(''*''))')
                   DO J=1,SIZE(ZTEM1D,1)
                       WRITE(INUM,'(I4,2X,F8.1,(5X,E15.8))')J,ZTEM1D(J),XVAR(1,1,J,NLOOPT,1,NLOOPP),ZWORKZ(J) 
                   ENDDO
                   ENDIF
                   IF(SIZE(XVAR,5) < 2)THEN

                      WRITE(INUM,'(36(''*''))')
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                   ENDIF

                ENDIF
!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!
              ENDDO

            ELSE                                    !TTTTTTTTTTTTTTTTTTTTTT

	      DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		NLOOPT=JLOOPT
		CALL RESOLV_TIMES(NLOOPT)
		WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NLOOPT,1)
! Partie reelle
	        ZWORKZ(:)=XVAR(1,1,:,NLOOPT,1,NLOOPP)
		ZMIN=MINVAL(ZTEM1D);ZMAX=MAXVAL(ZTEM1D)
		ZZMIN=MINVAL(ZWORKZ);ZZMAX=MAXVAL(ZWORKZ)
		CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		CALL AGSETF('FRA.',2.)
		CALL AGSETC('LAB/NAME.','B')
		CALL AGSETR('LAB/SU.',1.)
		CALL AGSETC('LAB/NAME.','L')
		CALL AGSETR('LAB/SU.',1.)
                CALL PCSETC('FC',':')
		IF(GOMEGAZOK)THEN                    !......................

!------
! _SPO_
!------
		  IF(LSPO)THEN

		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		    CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		    call plchhq(.9,.05,':PGL:X:PRL:Z:',.015,0.,1.)
		    call plchhq(0., .85,':PRU:(R):',.015,0.,-1.)

! Traitement de la partie imaginaire
		    IF(SIZE(XVAR,5) == 2)THEN
		      ZWORKZ(:)=XVAR(1,1,:,NLOOPT,2,NLOOPP)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PGL:X:PRL:Z:',.015,0.,1.)
		      call plchhq(0., .85,':PRU:(I):',.015,0.,-1.)
		    ENDIF

!---------------
! _OSPLO_ (/log)
!---------------
		  ELSE IF(LOSPLO)THEN

		    ZMIN=LOG10(ZMIN)
		    ZMAX=LOG10(ZMAX)
		    ZZMIN=MINVAL(ZWORKZ(:)*ZOMEGA)
		    ZZMAX=MAXVAL(ZWORKZ(:)*ZOMEGA)
		    CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		    CALL AGSETF('SET.',4.)
                    CALL EZXY(LOG10(ZTEM1D),ZWORKZ*ZOMEGA,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PGL:X:PRU:*(R):',.015,0.,-1.)

!------------------
! _LSPLO_ (Log/log)
!------------------
		  ELSE IF(LSPLO)THEN

		    CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		    CALL AGSETF('SET.',2.)
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PRU:Log(R):',.015,0.,-1.)
	            IF(LM5S3)THEN
		     CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		     CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
                     CALL FRSTPT(3.,0.)
                     CALL VECTOR(0.,5.)
                     CALL GSCHH(.2)
		     CALL GTX(0.+.5,5.-.4,YC5S3)
	            ENDIF

!---------------
! _PHALO_ (/log)
!---------------
		  ELSE IF(LPHALO)THEN

		    ZMIN=LOG10(ZMIN)
		    ZMAX=LOG10(ZMAX)
		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
	              ZWORKZ(:)=ATAN2(XVAR(1,1,:,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		    ENDIF

!-------
! _PHAO_
!-------
		  ELSE IF(LPHAO)THEN

		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
	              ZWORKZ(:)=ATAN2(XVAR(1,1,:,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PGL:X:PRL:Z:',.015,0.,1.)
		      call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		    ENDIF
		  ENDIF

		ELSE                                !......................

		  IF(LSPO)THEN
		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		    CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		    call plchhq(.9,.05,':PRU:Ind(:PRL:Z:PRU:):',.015,0.,1.)
		    call plchhq(0., .85,':PRU:(R):',.015,0.,-1.)

! Traitement de la partie imaginaire
		    IF(SIZE(XVAR,5) == 2)THEN
		      ZWORKZ(:)=XVAR(1,1,:,NLOOPT,2,NLOOPP)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres
		      CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
		      call plchhq(.9,.05,':PRU:Ind(:PRL:Z:PRU:):',.015,0.,1.)
		      call plchhq(0., .85,':PRU:(I):',.015,0.,-1.)
		    ENDIF
		  ELSE
		  ENDIF

		ENDIF                               !......................
		CALL FRAME

!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!
                IF(LPRINT)THEN

                  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
                  IF(IRESP /= 0)THEN
                    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
                    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
                    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
                  ENDIF

                   WRITE(INUM,'(''SP  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
& CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
                   IF(SIZE(XVAR,5) < 2)THEN
                     IF(GOMEGAZOK)THEN
                       WRITE(INUM,'(''Partie reelle uniquement     DOMEGAZ= '',F6.1)')ZOMEGA
                     ELSE
                       WRITE(INUM,'(''Partie reelle uniquement     DOMEGAZ= '',F6.1,'' -> Trace en indices de grille'')')ZOMEGA
                     ENDIF
                   ELSE
                     IF(GOMEGAZOK)THEN
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAZ= '',F6.1)')ZOMEGA
                     ELSE
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAZ= '',F6.1,'' -> Trace en indices de grille'')')ZOMEGA
                     ENDIF
                   ENDIF
                   WRITE(INUM,'(''NBVAL en K '',i4 )')SIZE(ZTEM1D,1)

                   IF(SIZE(XVAR,5) < 2)THEN
                      WRITE(INUM,'(36(''*''))')
                      WRITE(INUM,'(10X,''X(K)'',9X,''Y(VAL.R)'')')
                      WRITE(INUM,'(36(''*''))')
                   DO J=1,SIZE(ZTEM1D,1)
                       WRITE(INUM,'(I4,2X,F8.1,(5X,E15.8))')J,ZTEM1D(J),ZWORKZ(J) 
                   ENDDO
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                      WRITE(INUM,'(10X,''X(=K)'',8X,''Y(VAL.R)'',11X,''Y(VAL.Im)'')')
                      WRITE(INUM,'(55(''*''))')
                   DO J=1,SIZE(ZTEM1D,1)
                       WRITE(INUM,'(I4,2X,F8.1,(5X,E15.8))')J,ZTEM1D(J),XVAR(1,1,J,NLOOPT,1,NLOOPP),ZWORKZ(J) 
                   ENDDO
                   ENDIF
                   IF(SIZE(XVAR,5) < 2)THEN

                      WRITE(INUM,'(36(''*''))')
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                   ENDIF

                ENDIF
!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!
	      ENDDO

!+++++++++ Boucle temps +++++++++++++++++++++++++++++++++++

	    ENDIF                                   !TTTTTTTTTTTTTTTTTTTTTT
	  ENDDO

!+++++++++ Boucle processus +++++++++++++++++++++++++++++++++++

	  DEALLOCATE(ZWORKZ,ZTEM1D)

!*************************************************************************
! PH // X ou // Y
!*************************************************************************
	ELSE IF((II /= 1 .AND. IJ == 1) .OR. (II == 1 .AND. IJ /= 1))THEN
!       ELSE IF(II /= 1 .AND. IJ == 1 .AND. IK == 1)THEN
          print *,' unidimensionnel2: II,IJ,IK=',II,IJ,IK
          
	  IF(IJ == 1)THEN
! Disposition particuliere pour l'exploitation d'un fichier mal enregistre
! Juin 2001 c.a.d que le vecteur// Y est sur l'indice 1 de XVAR alors que
! NIL=NIH et NJL =/= NJH
            IF(NJL == NJH)THEN
! Cas normal
	      LSPX=.TRUE.
	    ELSE
! Cas anormal
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
	      LSPY=.TRUE.
	    ENDIF
	  ELSEIF(II == 1)THEN
	    LSPY=.TRUE.
	  ENDIF
	  IF(LSPX)THEN
            print*,'cas LPSX=',LSPX
	    ALLOCATE(ZTEM1D(SIZE(XVAR,1)),ZTEMLO(SIZE(XVAR,1)),ZWORKZ(SIZE(XVAR,1)))
	  ELSE
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
              print*,'cas anormal LPSY=',LSPY
	    ALLOCATE(ZTEM1D(SIZE(XVAR,1)),ZTEMLO(SIZE(XVAR,1)),ZWORKZ(SIZE(XVAR,1)))
	    ELSE
!ooooooooooooooooo
! Cas normal
              print*,'cas normal LPSY=',LSPY
            ALLOCATE(ZTEM1D(SIZE(XVAR,2)),ZTEMLO(SIZE(XVAR,2)),ZWORKZ(SIZE(XVAR,2)))
	    ENDIF
	  ENDIF
          if(nverbia > 0)then
           print *,' **subspxy LSPX,LSPY ',LSPX,LSPY
          endif

!+++++++++ Boucle processus +++++++++++++++++++++++++++++++++++

	  DO JLOOPP=1,NBPROCDIA(KLOOP)

	    NLOOPP=NPROCDIA(JLOOPP,KLOOP)
 	    print *,'subspxy NLOOPP',NLOOPP
	    CALL LOADUNITIT(JLOOPP,KLOOP)

!..............
	    IF(LSPX)THEN
	      IOMEGA=INDEX(CCOMMENT(NLOOPP),'DOMEGAX')
	        IF(IOMEGA == 0)THEN
	          IOMEGA=INDEX(CCOMMENT(NLOOPP),'Domegax')
	          IF(IOMEGA == 0)THEN
		    IOMEGA=INDEX(CCOMMENT(NLOOPP),'domegax')
	          ENDIF
	        ENDIF
	      IF(IOMEGA == 0)THEN
	        PRINT *,' Delta OmegaX (pulsation) non trouve dans le champ commentaire '
	        PRINT *,' On trace en indices de tableau'          
	        DO J=1,SIZE(ZTEM1D)
		  ZTEM1D(J)=J
	        ENDDO
	        GOMEGAXOK=.FALSE.
	        ZOMEGA=1.

	      ELSE

	        IEGAL=INDEX(CCOMMENT(NLOOPP)(IOMEGA:LEN_TRIM(CCOMMENT(NLOOPP))),'=')
	        READ(CCOMMENT(NLOOPP)(IOMEGA+IEGAL:LEN_TRIM(CCOMMENT(NLOOPP))),*)XOMEGAX
                print *,' tracé abscisse:j*OMEGAX ou 2pi/j*OMEGAX avec OMEGAX=',XOMEGAX
                IF(XOMEGAX == 0.)THEN
                  PRINT *,' Delta OmegaX (pulsation) = 0'
                  PRINT *,' On trace en indices de tableau'
                  DO J=1,SIZE(ZTEM1D)
                    ZTEM1D(J)=J
                  ENDDO
                  ZTEMLO(:)=ZTEM1D(:)
                  GOMEGAXOK=.FALSE.
                  ZOMEGA=1.
                ELSE
	          DO J=1,SIZE(ZTEM1D)
		    ZTEM1D(J)=J*XOMEGAX            ! lambda pour lin
		    ZTEMLO(J)=2*XPI/(J*XOMEGAX)    ! 2pi/lambda pour log
	          ENDDO
	          ZOMEGA=XOMEGAX
	          GOMEGAXOK=.TRUE.
	        ENDIF
	      ENDIF

!..............
	    ELSE

	    IOMEGA=INDEX(CCOMMENT(NLOOPP),'DOMEGAY')
	    IF(IOMEGA == 0)THEN
	      IOMEGA=INDEX(CCOMMENT(NLOOPP),'Domegay')
	      IF(IOMEGA == 0)THEN
		IOMEGA=INDEX(CCOMMENT(NLOOPP),'domegay')
	      ENDIF
	    ENDIF
          if(nverbia > 0)then
           print *,' **subspxy IOMEGA ',IOMEGA
          endif
	    IF(IOMEGA == 0)THEN
	      PRINT *,' Delta OmegaY (pulsation) non trouve dans le champ commentaire '
	      PRINT *,' On trace en indices de tableau'          
	      DO J=1,SIZE(ZTEM1D)
		ZTEM1D(J)=J
	      ENDDO
	      GOMEGAYOK=.FALSE.
              ZOMEGA=1.
	    ELSE
	      IEGAL=INDEX(CCOMMENT(NLOOPP)(IOMEGA:LEN_TRIM(CCOMMENT(NLOOPP))),'=')
	      READ(CCOMMENT(NLOOPP)(IOMEGA+IEGAL:LEN_TRIM(CCOMMENT(NLOOPP))),*)XOMEGAY
              print *,' tracé abscisse:j*OMEGAY ou 2pi/j*OMEGAY avec OMEGAY=',XOMEGAY
              IF(XOMEGAY == 0.)THEN
                PRINT *,' Delta OmegaY (pulsation) = 0 '
                PRINT *,' On trace en indices de tableau'
                DO J=1,SIZE(ZTEM1D)
                  ZTEM1D(J)=J
                ENDDO
                ZTEMLO(:)=ZTEM1D(:)
                GOMEGAYOK=.FALSE.
                ZOMEGA=1
              ELSE
	        DO J=1,SIZE(ZTEM1D)
		  ZTEM1D(J)=J*XOMEGAY            ! lambda pour lin
		  ZTEMLO(J)=2*XPI/(J*XOMEGAY)    ! 2pi/lambda pour Log
	        ENDDO
	        ZOMEGA=XOMEGAY
	        GOMEGAYOK=.TRUE.
              ENDIF
	    ENDIF
	    ENDIF
!..............
            IF(GOMEGAXOK .OR. GOMEGAYOK) THEN
              IF (LSPO .OR. LPHAO) THEN ! lin
	        ZMIN=MINVAL(ZTEM1D);ZMAX=MAXVAL(ZTEM1D)
              ELSE IF (LSPLO .OR. LOSPLO .OR. LPHALO) THEN ! Log
                ZMAX=MAXVAL(ZTEMLO)
                ! Elimination des valeurs <=0 a cause du Log
                ZMIN=ZMAX
 	        DO J=1,SIZE(ZTEMLO)
	          IF(ZTEMLO(J) > 0.)THEN
	            ZMIN=MIN(ZMIN,ZTEMLO(J))
	          ENDIF
 	        ENDDO
	        where(ZTEMLO <= 0.)ZTEMLO=1.e36
              END IF
	      print *,' ZMIN,ZMAX ',ZMIN,ZMAX
	    ENDIF
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++ Boucle sur N +++++++++++++++++++++++++++++++++++
!++++++++++Parties reelle et imaginaire++++++++++++++++++++

            DO JLOOPN=1,NBNDIA(KLOOP)
              NLOOPN=NNDIA(JLOOPN,KLOOP)
 	print *,'subspxy NLOOPN',NLOOPN

!+++++++++ Boucle sur K +++++++++++++++++++++++++++++++++++
            DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)
	      NLOOPK=NLVLKDIA(JLOOPK,KLOOP,NLOOPN)
!	print *,'subspxy jloopk,NLOOPK ',JLOOPK,NLOOPK

	    IF(.NOT.LTINCRDIA(KLOOP,1))THEN         !TTTTTTTTTTTTTTTTTTTTTT
        print *,'subspxy temps ',LTINCRDIA(KLOOP,1)

!+++++++++ Boucle temps +++++++++++++++++++++++++++++++++++

	      DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
		CALL RESOLV_TIMES(NLOOPT)
		WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NLOOPT,1)

! Partie reelle et imaginaire (suivant valeur de NLOOPN)
		IF(LSPX)THEN
	          ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,NLOOPN,NLOOPP)
		ELSE
!ooooooooooooooooooooooo
! PROVI pour lire vecteurs // Y mal ecrits chez VM
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	          ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,NLOOPN,NLOOPP)
            ELSE
! Cas normal
                  ZWORKZ(:)=XVAR(1,:,NLOOPK,NLOOPT,NLOOPN,NLOOPP)
            ENDIF
		ENDIF
               	ZZMIN=MINVAL(ZWORKZ);ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZZMIN,ZZMAX initiaux ',ZZMIN,ZZMAX
		IF(LVPTUSER)THEN
          	  CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		ELSE
          	  CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		ENDIF
		CALL AGSETF('FRA.',2.)
		CALL AGSETC('LAB/NAME.','B')
		CALL AGSETR('LAB/SU.',1.)
		CALL AGSETC('LAB/NAME.','L')
		CALL AGSETR('LAB/SU.',1.)
                CALL PCSETC('FC',':')

		IF((GOMEGAXOK .AND. LSPX) .OR. (GOMEGAYOK .AND. LSPY))THEN   !......................

!------
! _SPO_
!------
		  IF(LSPO)THEN

		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle (ou imaginaire)
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Ligne 0
		    CALL GSCLIP(1)
		    CALL GSLN(2)
		    CALL FRSTPT(ZMIN,0.)
		    CALL VECTOR(ZMAX,0.)
		    CALL SFLUSH
		    CALL GSLN(1)
! Titres


!---------------
! _OSPLO_ (/log)
!---------------
		  ELSE IF(LOSPLO)THEN

		    ZZMIN=MINVAL(ZWORKZ(:)*ZOMEGA)
		    ZZMAX=MAXVAL(ZWORKZ(:)*ZOMEGA)
		print *,' ZZMIN,ZZMAX *omega ',ZZMIN,ZZMAX
		    IF(LVPTUSER)THEN
          	      CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
		    ELSE
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
                    ENDIF
!Ds AGSETF, 4 signifie que l'on prend en compte les parametres de SET
! 2 que l'on prend en compte le seul dernier parametre et les 4 1ers
!Ds SET: 3 -> X Log + Y lin. 1 -> X+Y lin. 2-> X lin. + Y log 4 -> X log + Y Log
!                   CALL AGSETF('SET.',2.)
 	            CALL AGSETF('SET.',4.)
                    CALL EZXY(ZTEMLO,ZWORKZ*ZOMEGA,SIZE(ZTEMLO),0)
! Ligne 0
		    CALL GSCLIP(1)
		    CALL GSLN(2)
		    CALL FRSTPT(ZMIN,0.)
		    CALL VECTOR(ZMAX,0.)
		    CALL SFLUSH
		    CALL GSLN(1)
! Titres


!------------------
! _LSPLO_ (Log/log)
!------------------
		  ELSE IF(LSPLO)THEN
                    IF (ZZMAX <=0.) THEN
                      IF (NLOOPN==1) PRINT*,' LSPLO partie reelle <=0'
                      IF (NLOOPN==2) PRINT*,' LSPLO partie imaginaire <=0'
                      CYCLE 
                    END IF
                   ! Elimination des valeurs <=0 a cause du Log
                   ZZMIN=ZZMAX
 	           DO J=1,SIZE(ZWORKZ)
	             IF(ZWORKZ(J) > 0.)THEN
	               ZZMIN=MIN(ZZMIN,ZWORKZ(J))
	             ENDIF
 	           ENDDO
                    IF (ZZMIN ==ZZMAX) THEN
                     IF (NLOOPN==1) PRINT*,' LSPLO partie reelle>0 cst ',ZZMIN
                     IF (NLOOPN==2) PRINT*,' LSPLO partie imaginaire>0 cst ',ZZMIN
                     CYCLE 
                    END IF
                   where(ZWORKZ <= 0.)ZWORKZ=1.e36
		print *,' ZZMIN,ZZMAX corrigés ',ZZMIN,ZZMAX

		    IF(LVPTUSER)THEN
          	      CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		    ELSE
 	              CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		    ENDIF
! Verifier qu'avec 4 les limites sont mieux (NON)
                    CALL AGSETF('SET.',2.)
!                   CALL AGSETF('SET.',4.)
                    CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
		IF(LM5S3)THEN
		  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		  CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		  CALL FRSTPT(3.,0.)
		  CALL VECTOR(0.,5.)
		  CALL GSCHH(.2)
		  CALL GTX(0.+.5,5.-.4,YC5S3)
		ENDIF
! Titres

!---------------
! _PHALO_ (/log)
!---------------
		  ELSE IF(LPHALO)THEN
		    !!VM IF(NLOOPN == 2)exit

		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
		     IF(NLOOPN == 1) THEN ! Phase
		      IF(LSPX)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ELSE
!ooooooooooooooooooooooooooo
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
	    ELSE		 
! Cas normal
                        ZWORKZ(:)=ATAN2(-XVAR(1,:,NLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ENDIF
		      ENDIF
!                     print *,' PHALO ZWORKZ ',ZWORKZ
!                     print *,' PHALO ZWORKZ EN DEGRES ',ZWORKZ*45./ATAN(1.)
		      DO J=2,SIZE(ZWORKZ)
			IF(ABS(ZWORKZ(J-1) - ZWORKZ(J)) >= ATAN(1.)*8.)THEN
			  IF(ZWORKZ(J) >  0.)ZWORKZ(J)=ZWORKZ(J)+ATAN(1.)*8.
			  IF(ZWORKZ(J) <  0.)ZWORKZ(J)=ZWORKZ(J)-ATAN(1.)*8.
			ENDIF
		      ENDDO
!                     print *,' PHALO ZWORKZ AP DEROULEMENT PHASE ',ZWORKZ
!                     print *,' PHALO ZWORKZ AP DEROULEMENT PHASE EN DEGRES ',ZWORKZ*45./ATAN(1.)
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      print *,' ZZMIN,ZZMAX de la phase ',ZZMIN,ZZMAX
		      IF(LVPTUSER)THEN
          	        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
		      ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
		      ENDIF
		      CALL AGSETF('SET.',4.)
!                     print *,' PHALO AV EZXY '
                      CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
!                     print *,' PHALO AP EZXY '
                     ELSE IF(NLOOPN == 2) THEN ! Module
		      IF(LSPX)THEN
	                ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)   )
		      ELSE
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)   )
	    ELSE
                        ZWORKZ(:)=XVAR(1,:,NLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(1,:,NLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(1,:,NLOOPK,NLOOPT,2,NLOOPP)   )
            ENDIF
		      ENDIF
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      print *,' ZZMIN,ZZMAX du module ',ZZMIN,ZZMAX
                      ! 4 (log X, log Y) plutot que 3 (log X, linear Y)
		      IF(LVPTUSER)THEN
          	        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		      ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
		IF(LM5S3)THEN
		  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		  CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		  CALL FRSTPT(3.,0.)
		  CALL VECTOR(0.,5.)
		  CALL GSCHH(.2)
		  CALL GTX(0.+.5,5.-.4,YC5S3)
		ENDIF
		     ENDIF ! fin NLOOPN
		    ENDIF ! fin (SIZE(XVAR,5) < 2)

!-------
! _PHAO_
!-------
		  ELSE IF(LPHAO)THEN
		    !!VM IF(NLOOPN == 2)exit

		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
		     IF(NLOOPN == 1) THEN ! Phase
		      IF(LSPX)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ELSE
!ooooooooooooooooooooooooooo
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ELSE
! Cas normal
                        ZWORKZ(:)=ATAN2(-XVAR(1,:,NLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ENDIF
		      ENDIF
!                     print *,' PHAO ZWORKZ ',ZWORKZ
		      DO J=2,SIZE(ZWORKZ)
			IF(ABS(ZWORKZ(J-1) - ZWORKZ(J)) >= ATAN(1.)*8.)THEN
			  IF(ZWORKZ(J) >  0.)ZWORKZ(J)=ZWORKZ(J)+ATAN(1.)*8.
			  IF(ZWORKZ(J) <  0.)ZWORKZ(J)=ZWORKZ(J)-ATAN(1.)*8.
			ENDIF
		      ENDDO
!                     print *,' PHAO ZWORKZ AP DEROULEMENT PHASE ',ZWORKZ
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		      print *,' ZZMIN,ZZMAX de la phase ',ZZMIN,ZZMAX
		      IF(LVPTUSER)THEN
          	        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
                     ELSE IF(NLOOPN == 2) THEN ! Module *k
		      IF(LSPX)THEN
	                ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)   )
		      ELSE
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,NLOOPK,NLOOPT,2,NLOOPP)   )
            ELSE
                        ZWORKZ(:)=XVAR(1,:,NLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(1,:,NLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(1,:,NLOOPK,NLOOPT,2,NLOOPP)   )
            ENDIF
		      ENDIF
                      ! Module * k
                      ZWORKZ(:)=ZWORKZ(:)*ZTEMLO(:)
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZZMIN,ZZMAX du Module *K ',ZZMIN,ZZMAX
                      ! 4 (log X, log Y) 
		      IF(LVPTUSER)THEN
          	        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		      ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
		IF(LM5S3)THEN
		  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		  CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		  CALL FRSTPT(3.,0.)
		  CALL VECTOR(0.,5.)
		  CALL GSCHH(.2)
		  CALL GTX(0.+.5,5.-.4,YC5S3)
		ENDIF
		     ENDIF ! fin boucle NLOOPN
		    ENDIF ! fin (SIZE(XVAR,5) < 2)
		   ENDIF ! fin LSPO,LOSPLO,LSPLO,LPHALO,LPHAO

		ELSE                                !......................

		  IF(LSPO)THEN
		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Titres

                  ELSE
                  ENDIF

		ENDIF                               !......................

		IF(GOMEGAXOK .OR. GOMEGAYOK)THEN   !GGGGGGGGGGGGGGG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Titres 
	        CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
	        call gsclip(0)
! Titres en X et Y
                IF(LSPO)THEN
		  IF(GOMEGAXOK)THEN
		    call plchhq(.9,.05,':PGL:X:PRL:X:',.015,0.,1.)
		  ELSEIF(GOMEGAYOK)THEN
		    call plchhq(.9,.05,':PGL:X:PRL:Y:',.015,0.,1.)
                  ELSE
		    IF(LSPX)THEN
		      call plchhq(.9,.05,':PRU:Ind(:PRL:X:PRU:):',.015,0.,1.)
		    ELSE
		      call plchhq(.9,.05,':PRU:Ind(:PRL:Y:PRU:):',.015,0.,1.)
		    ENDIF
                  ENDIF
                ELSEIF(LOSPLO)THEN
		  IF(LSPX)THEN
		    !!VM call plchhq(.9,.05,':PGL:X:PRL:X:PRU:',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:X:PRU:',.015,0.,1.)
                  ELSE
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Y:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:Y:PRU:',.015,0.,1.)
		  ENDIF
		ELSEIF(LSPLO)THEN
		  IF(LSPX)THEN
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:X:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:X:PRU:',.015,0.,1.)
                  ELSE
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Y:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:Y:PRU:',.015,0.,1.)
		  ENDIF
		ELSEIF(LPHALO)THEN
		  IF(LSPX)THEN
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:X:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:X:PRU:',.015,0.,1.)
                  ELSE
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Y:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:Y:PRU:',.015,0.,1.)
		  ENDIF
		ELSEIF(LPHAO)THEN
		  IF(LSPX)THEN
		    call plchhq(.9,.05,':PGL:X:PRL:X:',.015,0.,1.)
                  ELSE
		    call plchhq(.9,.05,':PGL:X:PRL:Y:',.015,0.,1.)
		  ENDIF
		ENDIF
	        IF(NLOOPN == 2)THEN
                  IF(LSPO)THEN
		    call plchhq(0., .87,':PRU:(I):',.015,0.,-1.)
                  ELSEIF(LOSPLO)THEN
		    call plchhq(0., .87,':PGL:X:PRU:*(I):',.015,0.,-1.)
		  ELSEIF(LSPLO)THEN
		    call plchhq(0., .87,':PRU:(I):',.015,0.,-1.)
		  ELSEIF(LPHALO)THEN
		    call plchhq(0., .85,':PRU:Module:',.015,0.,-1.)
		  ELSEIF(LPHAO)THEN
		    call plchhq(0., .85,':PRU:K*Module:',.015,0.,-1.)
		  ENDIF
	        ELSE
                  IF(LSPO)THEN
		    call plchhq(0., .87,':PRU:(R):',.015,0.,-1.)
                  ELSEIF(LOSPLO)THEN
		    call plchhq(0., .87,':PGL:X:PRU:*(R):',.015,0.,-1.)
		  ELSEIF(LSPLO)THEN
		    call plchhq(0., .87,':PRU:(R):',.015,0.,-1.)
		  ELSEIF(LPHALO)THEN
		    call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		  ELSEIF(LPHAO)THEN
		    call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		  ENDIF
	        ENDIF
! Titres top et bottom
! Top1
		YTEM(1:LEN(YTEM))=' '
		CALL RESOLV_TIT('CTITT1',YTEM)
		ZXPOSTITT1=.002
		ZXYPOSTITT1=.98
		IF(XPOSTITT1 /= 0.)THEN
		  ZXPOSTITT1=XPOSTITT1
	        ENDIF
		IF(XYPOSTITT1 /= 0.)THEN
		  ZXYPOSTITT1=XYPOSTITT1
                ENDIF
		IF(YTEM /= ' ')THEN
		  IF(XSZTITT1 /= 0.)THEN
		    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM,XSZTITT1,0.,-1.)
		  ELSE
		    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM,.008,0.,-1.)
		  ENDIF
		ELSE
		  YTEM=CGROUP(1:LEN_TRIM(CGROUP))
		  YTEM=ADJUSTL(YTEM)
		  IL=LEN_TRIM(YTEM)
		  YTEM(IL+3:IL+5)='K ='
		  IL=IL+6
		  WRITE(YTEM(IL:IL+2),'(I3)')NLOOPK
	          call plchhq(.05,.98,YTEM(1:LEN_TRIM(YTEM)),.015,0.,-1.)
!                 call plchhq(.05,.98,CGROUP(1:LEN_TRIM(CGROUP)),.015,0.,-1.)
                ENDIF
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
		  ELSE
		    CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,.008,0.,-1.)
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
		  ELSE
		    CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,.008,0.,-1.)
		  ENDIF
                ENDIF
! Titres Bottom
! Titre N1 BOTTOM
		YTEM(1:LEN(YTEM))=' '
		YTEM=CTIMEC
		YTEM=ADJUSTL(YTEM)
                CALL RESOLV_TIT('CTITB1',YTEM)
                ZXPOSTITB1=.002
                ZXYPOSTITB1=.005
                IF(XPOSTITB1 /= 0.)THEN
                  ZXPOSTITB1=XPOSTITB1
                ENDIF
                IF(XYPOSTITB1 /= 0.)THEN
                  ZXYPOSTITB1=XYPOSTITB1
                ENDIF
                IF(YTEM /= ' ')THEN
                  CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,YTEM,.009,0.,-1.)
              !   CALL PLCHHQ(0.002,0.005,YTEM,.007,0.,-1.)
                ENDIF
! Titre N2 BOTTOM
		YTEM(1:LEN(YTEM))=' '
                CALL RESOLV_TIT('CTITB2',YTEM)
                ZXPOSTITB2=.002
                ZXYPOSTITB2=.025
                IF(XPOSTITB2 /= 0.)THEN
                  ZXPOSTITB2=XPOSTITB2
                ENDIF
                IF(XYPOSTITB2 /= 0.)THEN
                  ZXYPOSTITB2=XYPOSTITB2
                ENDIF
                IF(YTEM /= ' ')THEN
                  CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,YTEM,.007,0.,-1.)
              !   CALL PLCHHQ(0.002,0.025,YTEM,.007,0.,-1.)
                ENDIF
! Titre N3 BOTTOM
                YTEM(1:LEN(YTEM))=' '
                CALL RESOLV_TIT('CTITB3',YTEM)
                ZXPOSTITB3=.002
                ZXYPOSTITB3=.045
                IF(XPOSTITB3 /= 0.)THEN
                  ZXPOSTITB3=XPOSTITB3
                ENDIF
                IF(XYPOSTITB3 /= 0.)THEN
                  ZXYPOSTITB3=XYPOSTITB3
                ENDIF
                IF(YTEM /= ' ')THEN
                  CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,.007,0.,-1.)
                ENDIF
	        IF(LDATFILE)CALL DATFILE_FORDIACHRO
	        call gsclip(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ENDIF        !GGGGGGGGGGGGGGG
		CALL FRAME
!               print *,'subspxy ap frame '
!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!
                IF(LPRINT)THEN
                  IF(SIZE(XVAR,5) == 2 .AND. NLOOPN == 1)CYCLE
                  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
                  IF(IRESP /= 0)THEN
                    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
                    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
                    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
                  ENDIF

                   WRITE(INUM,'(''SP  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
& CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
!........
                   IF(SIZE(XVAR,5) < 2)THEN
                     IF(GOMEGAXOK)THEN
                       WRITE(INUM,'(''Partie reelle uniquement     DOMEGAX= '',F6.1)')ZOMEGA

                     ELSE
                      
                       IF(GOMEGAYOK)THEN
                         WRITE(INUM,'(''Partie reelle uniquement     DOMEGAY= '',F6.1)')ZOMEGA
                       ELSE

                         IF(LSPX)THEN
                           WRITE(INUM,'(''Partie reelle uniquement     DOMEGAX= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAX
                         ELSE
                           WRITE(INUM,'(''Partie reelle uniquement     DOMEGAY= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAY
                         ENDIF
                       ENDIF
                     ENDIF

                   ELSE
!........
                     IF(GOMEGAXOK)THEN
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAX= '',F6.1)')XOMEGAX
                     ELSE

                       IF(GOMEGAYOK)THEN
                         WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAY= '',F6.1)')XOMEGAY
                       ELSE
                         
                         IF(LSPX)THEN
                           WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAX= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAX
                         ELSE
                           WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAY= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAY
                         ENDIF
                       ENDIF
                     ENDIF
                   ENDIF
!........    
                   IF(LSPX)THEN
                     WRITE(INUM,'(''NBVAL en I '',i4 )')SIZE(ZTEM1D,1)
                   ELSE
                     WRITE(INUM,'(''NBVAL en J '',i4 )')SIZE(ZTEM1D,1)
                   ENDIF

                   IF(SIZE(XVAR,5) < 2)THEN

                      WRITE(INUM,'(36(''*''))')
                      IF(LSPX)THEN
                        WRITE(INUM,'(10X,''X(I)'',9X,''Y(VAL.R)'')')
                      ELSE
                        WRITE(INUM,'(10X,''X(J)'',9X,''Y(VAL.R)'')')
                      ENDIF
                      WRITE(INUM,'(36(''*''))')
                      DO J=1,SIZE(ZTEM1D,1)
                        WRITE(INUM,'(I4,2X,F8.1,(5X,E15.8))')J,ZTEM1D(J),ZWORKZ(J) 
                      ENDDO
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                      IF(LSPX)THEN
                        WRITE(INUM,'(10X,''X(=I)'',8X,''Y(VAL.R)'',11X,''Y(VAL.Im)'')')
                      ELSE
                        WRITE(INUM,'(10X,''X(=J)'',8X,''Y(VAL.R)'',11X,''Y(VAL.Im)'')')
                      ENDIF
                      WRITE(INUM,'(55(''*''))')
                      DO J=1,SIZE(ZTEM1D,1)
                        IF(LSPX)THEN
                          WRITE(INUM,'(I4,2X,F8.1,2(5X,E15.8))')J,ZTEM1D(J),XVAR(J,1,NLOOPK,NLOOPT,1,NLOOPP),ZWORKZ(J) 
                        ELSE
                          WRITE(INUM,'(I4,2X,F8.1,2(5X,E15.8))')J,ZTEM1D(J),XVAR(1,J,NLOOPK,NLOOPT,1,NLOOPP),ZWORKZ(J) 
                        ENDIF
                      ENDDO
                    ENDIF
                   IF(SIZE(XVAR,5) < 2)THEN

                      WRITE(INUM,'(36(''*''))')
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                   ENDIF

                ENDIF
!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!

              ENDDO


            ELSE                                    !TTTTTTTTTTTTTTTTTTTTTT
             print *,'subspxy boucle temps ',LTINCRDIA(KLOOP,1)


	      DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		NLOOPT=JLOOPT
		CALL RESOLV_TIMES(NLOOPT)
		WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NLOOPT,1)
! Partie reelle et imaginaire
		IF(LSPX)THEN
                  ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,NLOOPN,NLOOPP)
		ELSE
!oooooooooooooooooooooo
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
                  ZWORKZ(:)=XVAR(:,1,NLOOPK,NLOOPT,NLOOPN,NLOOPP)
            ELSE
! Cas normal
                  ZWORKZ(:)=XVAR(1,:,NLOOPK,NLOOPT,NLOOPN,NLOOPP)
            ENDIF
		ENDIF
                IF (LSPO .OR. LPHAO) THEN ! lin
		  ZMIN=MINVAL(ZTEM1D);ZMAX=MAXVAL(ZTEM1D)
                ELSE IF (LSPLO .OR. LOSPLO .OR. LPHALO) THEN ! Log
		  ZMAX=MAXVAL(ZTEMLO)
                  ZMIN=ZMAX
		  DO J=1,SIZE(ZTEMLO)
		    IF(ZTEMLO(J) > 0.)THEN
			ZMIN=MIN(ZMIN,ZTEMLO(J))
		    ENDIF
		  ENDDO
		  where(ZTEMLO <= 0.)ZTEMLO=1.e36
                END IF
		ZZMIN=MINVAL(ZWORKZ);ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZMIN,ZMAX,ZZMIN,ZZMAX initiaux ',ZMIN,ZMAX,ZZMIN,ZZMAX
	        IF(LVPTUSER)THEN
                  CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
	        ELSE
		  CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		ENDIF
		CALL AGSETF('FRA.',2.)
		CALL AGSETC('LAB/NAME.','B')
		CALL AGSETR('LAB/SU.',1.)
		CALL AGSETC('LAB/NAME.','L')
		CALL AGSETR('LAB/SU.',1.)
                CALL PCSETC('FC',':')

		IF((GOMEGAXOK .AND. LSPX) .OR. (GOMEGAYOK .AND. LSPY))THEN      !......................

!------
! _SPO_
!------
		  IF(LSPO)THEN

		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle et imaginaire (suivant la valeur de N)
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
! Ligne 0
		    CALL GSCLIP(1)
		    CALL GSLN(2)
		    CALL FRSTPT(ZMIN,0.)
		    CALL VECTOR(ZMAX,0.)
		    CALL SFLUSH
		    CALL GSLN(1)
! Titres


!---------------
! _OSPLO_ (/log)
!---------------
		  ELSE IF(LOSPLO)THEN

		    ZZMIN=MINVAL(ZWORKZ(:)*ZOMEGA)
		    ZZMAX=MAXVAL(ZWORKZ(:)*ZOMEGA)
		print *,' ZZMIN,ZZMAX *omega ',ZZMIN,ZZMAX
	            IF(LVPTUSER)THEN
                      CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
	            ELSE
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
                    ENDIF
		    CALL AGSETF('SET.',4.)
                    CALL EZXY(ZTEMLO,ZWORKZ*ZOMEGA,SIZE(ZTEMLO),0)
! Ligne 0
		    CALL GSCLIP(1)
		    CALL GSLN(2)
		    CALL FRSTPT(ZMIN,0.)
		    CALL VECTOR(ZMAX,0.)
		    CALL SFLUSH
		    CALL GSLN(1)
! Titres

!------------------
! _LSPLO_ (Log/log)
!------------------
		  ELSE IF(LSPLO)THEN
                    IF (ZZMAX <=0.) THEN
                      IF (NLOOPN==1) PRINT*,' LSPLO partie reelle <=0'
                      IF (NLOOPN==2) PRINT*,' LSPLO partie imaginaire <=0'
                      CYCLE 
                    END IF
                   ! Elimination des valeurs <=0 a cause du Log
                   ZZMIN=ZZMAX
 	           DO J=1,SIZE(ZWORKZ)
	             IF(ZWORKZ(J) > 0.)THEN
	               ZZMIN=MIN(ZZMIN,ZWORKZ(J))
	             ENDIF
 	           ENDDO
                    IF (ZZMIN ==ZZMAX) THEN
                     IF (NLOOPN==1) PRINT*,' LSPLO partie reelle>0 cst ',ZZMIN
                     IF (NLOOPN==2) PRINT*,' LSPLO partie imaginaire>0 cst ',ZZMIN
                     CYCLE 
                    END IF
                    where(ZWORKZ <= 0.)ZWORKZ=1.e36
		print *,' ZZMIN,ZZMAX corrigés ',ZZMIN,ZZMAX
                    IF(LVPTUSER)THEN
                      CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
	            ELSE
		      CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		    ENDIF
!                   CALL AGSETF('SET.',4.)
        	    CALL AGSETF('SET.',2.)
                    CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
		IF(LM5S3)THEN
		  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		  CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		  CALL FRSTPT(3.,0.)
		  CALL VECTOR(0.,5.)
		  CALL GSCHH(.2)
		  CALL GTX(0.+.5,5.-.4,YC5S3)
		ENDIF
! Titres

!---------------
! _PHALO_ (/log)
!---------------
		  ELSE IF(LPHALO)THEN

		    !!VM IF(NLOOPN == 2)exit

		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
                     IF(NLOOPN==1) THEN ! Phase
		      IF(LSPX)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ELSE
!ooooooooooooooooooooooooooo
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ELSE
! Cas normal
                        ZWORKZ(:)=ATAN2(-XVAR(1,:,JLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ENDIF
		      ENDIF
		      DO J=2,SIZE(ZWORKZ)
			IF(ABS(ZWORKZ(J-1) - ZWORKZ(J)) >= ATAN(1.)*8.)THEN
			  IF(ZWORKZ(J) >  0.)ZWORKZ(J)=ZWORKZ(J)+ATAN(1.)*8.
			  IF(ZWORKZ(J) <  0.)ZWORKZ(J)=ZWORKZ(J)-ATAN(1.)*8.
			ENDIF
		      ENDDO
!                     print *,' PHALO ZWORKZ AP DEROULEMENT PHASE ',ZWORKZ
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZZMIN,ZZMAX de la phase ',ZZMIN,ZZMAX
	              IF(LVPTUSER)THEN
                        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
	              ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,3)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEM1D),0)
                     ELSE IF(NLOOPN==2) THEN ! Module
		      IF(LSPX)THEN
	                ZWORKZ(:)=XVAR(:,1,JLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)   )
		      ELSE
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=XVAR(:,1,JLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)   )
            ELSE
                        ZWORKZ(:)=XVAR(1,:,JLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(1,:,JLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(1,:,JLOOPK,NLOOPT,2,NLOOPP)   )
            ENDIF
		      ENDIF
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZZMIN,ZZMAX du module ',ZZMIN,ZZMAX
                      ! 4 (log X, log Y) plutot que 3 (log X, linear Y)
	              IF(LVPTUSER)THEN
                        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
	              ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
		IF(LM5S3)THEN
		  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		  CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		  CALL FRSTPT(3.,0.)
		  CALL VECTOR(0.,5.)
		  CALL GSCHH(.2)
		  CALL GTX(0.+.5,5.-.4,YC5S3)
		ENDIF
                     ENDIF ! fin boucle NLOOPN
		    ENDIF

!-------
! _PHAO_
!-------
		  ELSE IF(LPHAO)THEN

		    !!VM IF(NLOOPN == 2)exit

		    IF(SIZE(XVAR,5) < 2)THEN
		      print *,' Absence partie imaginaire. Representation impossible sous cette forme'
		    ELSE
                     IF(NLOOPN==1) THEN ! Phase
		      IF(LSPX)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
		      ELSE
!ooooooooooooooooooooooooooo
! Disposition particuliere pour le traitement des vecteurs // Y  mal enreg.
! Cas anormal
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=ATAN2(-XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ELSE
! Cas normal
                        ZWORKZ(:)=ATAN2(-XVAR(1,:,JLOOPK,NLOOPT,2,NLOOPP),ZWORKZ(:))
            ENDIF
		      ENDIF
!                     print *,' PHALO ZWORKZ ',ZWORKZ
		      DO J=2,SIZE(ZWORKZ)
			IF(ABS(ZWORKZ(J-1) - ZWORKZ(J)) >= ATAN(1.)*8.)THEN
			  IF(ZWORKZ(J) >  0.)ZWORKZ(J)=ZWORKZ(J)+ATAN(1.)*8.
			  IF(ZWORKZ(J) <  0.)ZWORKZ(J)=ZWORKZ(J)-ATAN(1.)*8.
			ENDIF
		      ENDDO
!                     print *,' PHALO ZWORKZ AP DEROULEMENT PHASE ',ZWORKZ
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZZMIN,ZZMAX de la phase ',ZZMIN,ZZMAX
	              IF(LVPTUSER)THEN
                        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
	              ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,1)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
                     ELSE IF(NLOOPN==2) THEN ! Module *k
		      IF(LSPX)THEN
	                ZWORKZ(:)=XVAR(:,1,JLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)   )
		      ELSE
            IF(IJ == 1 .AND. NJL /= NJH)THEN
	                ZWORKZ(:)=XVAR(:,1,JLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(:,1,JLOOPK,NLOOPT,2,NLOOPP)   )
            ELSE
                        ZWORKZ(:)=XVAR(1,:,JLOOPK,NLOOPT,1,NLOOPP)
	                ZWORKZ(:)=SQRT(ZWORKZ(:)*ZWORKZ(:)              + &
                                       XVAR(1,:,JLOOPK,NLOOPT,2,NLOOPP)*&
                                       XVAR(1,:,JLOOPK,NLOOPT,2,NLOOPP)   )
            ENDIF
		      ENDIF
                      ! Module * k
                      ZWORKZ(:)=ZWORKZ(:)*ZTEMLO(:)
		      ZZMIN=MINVAL(ZWORKZ)
		      ZZMAX=MAXVAL(ZWORKZ)
		print *,' ZZMIN,ZZMAX de K*Module ',ZZMIN,ZZMAX
                      ! 4 (log X, log Y) 
	              IF(LVPTUSER)THEN
                        CALL SET(XVPTL,XVPTR,XVPTB,XVPTT,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
	              ELSE
		        CALL SET(.1,.9,.1,.9,ZMIN,ZMAX,ZZMIN,ZZMAX,4)
		      ENDIF
		      CALL AGSETF('SET.',4.)
                      CALL EZXY(ZTEMLO,ZWORKZ,SIZE(ZTEMLO),0)
		IF(LM5S3)THEN
		  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
                      print*,'out GETSET',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID
		  CALL SET(ZVR-.3,ZVR-.05,ZVT-.3,ZVT-.05,0.,5.,0.,5.,1)
		  CALL FRSTPT(3.,0.)
		  CALL VECTOR(0.,5.)
		  CALL GSCHH(.2)
		  CALL GTX(0.+.5,5.-.4,YC5S3)
		ENDIF
                     ENDIF ! fin boucle NLOOPN
		    ENDIF
		  ENDIF

		ELSE                                !......................

		  IF(LSPO)THEN
		    CALL AGSETF('SET.',4.)
! Traitement de la partie reelle
                    CALL EZXY(ZTEM1D,ZWORKZ,SIZE(ZTEM1D),0)
		  ELSE
		  ENDIF

		ENDIF                               !......................

		IF(GOMEGAXOK .OR. GOMEGAYOK)THEN     !GGGGGGGGGGGGGGG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Titres 
	        CALL SET(0.,.9,0.,.9,0.,.9,0.,.9,1)
	        call gsclip(0)
! Titres en X et Y
                IF(LSPO)THEN
		  IF(GOMEGAXOK)THEN
		    call plchhq(.9,.05,':PGL:X:PRL:X:',.015,0.,1.)
		  ELSEIF(GOMEGAYOK)THEN
		    call plchhq(.9,.05,':PGL:X:PRL:Y:',.015,0.,1.)
		  ELSE
		    IF(LSPX)THEN
		      call plchhq(.9,.05,':PRU:Ind(:PRL:X:PRU:):',.015,0.,1.)
		    ELSE
		      call plchhq(.9,.05,':PRU:Ind(:PRL:Y:PRU:):',.015,0.,1.)
		    ENDIF
		  ENDIF
                ELSEIF(LOSPLO)THEN
		  IF(LSPX)THEN
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:X:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:X:PRU:',.015,0.,1.)
		  ELSE
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Y:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:Y:PRU:',.015,0.,1.)
		  ENDIF
		ELSEIF(LSPLO)THEN
		  IF(LSPX)THEN
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:X:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:X:PRU:',.015,0.,1.)
		  ELSE
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Y:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:Y:PRU:',.015,0.,1.)
		  ENDIF
		ELSEIF(LPHALO)THEN
		  IF(LSPX)THEN
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:X:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:X:PRU:',.015,0.,1.)
		  ELSE
		    !!VM call plchhq(.9,.05,':PRU:Log(:PGL:X:PRL:Y:PRU:):',.015,0.,1.)
		    call plchhq(.9,.05,':PRL:K:PRL:Y:PRU:',.015,0.,1.)
		  ENDIF
		ELSEIF(LPHAO)THEN
		  IF(LSPX)THEN
		    call plchhq(.9,.05,':PGL:X:PRL:X:',.015,0.,1.)
		  ELSE
		    call plchhq(.9,.05,':PGL:X:PRL:Y:',.015,0.,1.)
		  ENDIF
		ENDIF
	        IF(NLOOPN == 2)THEN
                  IF(LSPO)THEN
		    call plchhq(0., .87,':PRU:(I):',.015,0.,-1.)
                  ELSEIF(LOSPLO)THEN
		    call plchhq(0., .87,':PGL:X:PRU:*(I):',.015,0.,-1.)
		  ELSEIF(LSPLO)THEN
		    call plchhq(0., .87,':PRU:(I):',.015,0.,-1.)
		  ELSEIF(LPHALO)THEN
		    call plchhq(0., .85,':PRU:Module:',.015,0.,-1.)
		  ELSEIF(LPHAO)THEN
		    call plchhq(0., .85,':PRU:K*Module:',.015,0.,-1.)
		  ENDIF
	        ELSE
                  IF(LSPO)THEN
		    call plchhq(0., .87,':PRU:(R):',.015,0.,-1.)
                  ELSEIF(LOSPLO)THEN
		    call plchhq(0., .87,':PGL:X:PRU:*(R):',.015,0.,-1.)
		  ELSEIF(LSPLO)THEN
		    call plchhq(0., .87,':PRU:(R):',.015,0.,-1.)
		  ELSEIF(LPHALO)THEN
		    call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		  ELSEIF(LPHAO)THEN
		    call plchhq(0., .85,':PRU:Phase:',.015,0.,-1.)
		  ENDIF
	        ENDIF
! Titres top et bottom
! Top1
		YTEM(1:LEN(YTEM))=' '
		CALL RESOLV_TIT('CTITT1',YTEM)
		ZXPOSTITT1=.002
		ZXYPOSTITT1=.98
		IF(XPOSTITT1 /= 0.)THEN
		  ZXPOSTITT1=XPOSTITT1
	        ENDIF
		IF(XYPOSTITT1 /= 0.)THEN
		  ZXYPOSTITT1=XYPOSTITT1
                ENDIF
		IF(YTEM /= ' ')THEN
		  IF(XSZTITT1 /= 0.)THEN
		    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM,XSZTITT1,0.,-1.)
		  ELSE
		    CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM,.008,0.,-1.)
		  ENDIF
		ELSE
		  YTEM=CGROUP(1:LEN_TRIM(CGROUP))
		  YTEM=ADJUSTL(YTEM)
		  IL=LEN_TRIM(YTEM)
		  YTEM(IL+3:IL+5)='K ='
		  IL=IL+6
		  WRITE(YTEM(IL:IL+2),'(I3)')NLOOPK
	          call plchhq(.05,.98,YTEM(1:LEN_TRIM(YTEM)),.015,0.,-1.)
!                 call plchhq(.05,.98,CGROUP(1:LEN_TRIM(CGROUP)),.015,0.,-1.)
                ENDIF
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
		  ELSE
		    CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM,.008,0.,-1.)
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
		  ELSE
		    CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM,.008,0.,-1.)
		  ENDIF
                ENDIF
! Titres Bottom
! Titre N1 BOTTOM
		YTEM(1:LEN(YTEM))=' '
		YTEM=CTIMEC
		YTEM=ADJUSTL(YTEM)
                CALL RESOLV_TIT('CTITB1',YTEM)
                ZXPOSTITB1=.002
                ZXYPOSTITB1=.005
                IF(XPOSTITB1 /= 0.)THEN
                  ZXPOSTITB1=XPOSTITB1
                ENDIF
                IF(XYPOSTITB1 /= 0.)THEN
                  ZXYPOSTITB1=XYPOSTITB1
                ENDIF
                IF(YTEM /= ' ')THEN
                  CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,YTEM,.009,0.,-1.)
              !   CALL PLCHHQ(0.002,0.005,YTEM,.007,0.,-1.)
                ENDIF
! Titre N2 BOTTOM
		YTEM(1:LEN(YTEM))=' '
                CALL RESOLV_TIT('CTITB2',YTEM)
                ZXPOSTITB2=.002
                ZXYPOSTITB2=.025
                IF(XPOSTITB2 /= 0.)THEN
                  ZXPOSTITB2=XPOSTITB2
                ENDIF
                IF(XYPOSTITB2 /= 0.)THEN
                  ZXYPOSTITB2=XYPOSTITB2
                ENDIF
                IF(YTEM /= ' ')THEN
                  CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,YTEM,.007,0.,-1.)
              !   CALL PLCHHQ(0.002,0.025,YTEM,.007,0.,-1.)
                ENDIF
! Titre N3 BOTTOM
                YTEM(1:LEN(YTEM))=' '
                CALL RESOLV_TIT('CTITB3',YTEM)
                ZXPOSTITB3=.002
                ZXYPOSTITB3=.045
                IF(XPOSTITB3 /= 0.)THEN
                  ZXPOSTITB3=XPOSTITB3
                ENDIF
                IF(XYPOSTITB3 /= 0.)THEN
                  ZXYPOSTITB3=XYPOSTITB3
                ENDIF
                IF(YTEM /= ' ')THEN
                  CALL PLCHHQ(ZXPOSTITB3,ZXYPOSTITB3,YTEM,.007,0.,-1.)
                ENDIF
	        IF(LDATFILE)CALL DATFILE_FORDIACHRO
	        call gsclip(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	        print *,'subspxy ap frame '
		ENDIF        !GGGGGGGGGGGGGGG
		CALL FRAME
!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!
                IF(LPRINT)THEN
                  IF(SIZE(XVAR,5) == 2 .AND. NLOOPN == 1)CYCLE

                  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
                  IF(IRESP /= 0)THEN
                    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
                    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
                    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
                  ENDIF

                   WRITE(INUM,'(''SP  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')')CGROUP,&
& CTITRE(NLOOPP)(1:25),XTRAJT(NLOOPT,1)
!........
                   IF(SIZE(XVAR,5) < 2)THEN
                     IF(GOMEGAXOK)THEN
                       WRITE(INUM,'(''Partie reelle uniquement     DOMEGAX= '',F6.1)')ZOMEGA

                     ELSE

                       IF(GOMEGAYOK)THEN
                         WRITE(INUM,'(''Partie reelle uniquement     DOMEGAY= '',F6.1)')ZOMEGA
                       ELSE

                         IF(LSPX)THEN
                           WRITE(INUM,'(''Partie reelle uniquement     DOMEGAX= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAX
                         ELSE
                           WRITE(INUM,'(''Partie reelle uniquement     DOMEGAY= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAY
                         ENDIF
                       ENDIF
                     ENDIF

                   ELSE
!........
                     IF(GOMEGAXOK)THEN
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAX= '',F6.1)')ZOMEGA
                     ELSE

                       IF(GOMEGAYOK)THEN
                       WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAY= '',F6.1)')ZOMEGA
                       ELSE
                         IF(LSPX)THEN
                           WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAX= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAX
                         ELSE
                           WRITE(INUM,'(''Parties reelle + imaginaire  DOMEGAY= '',F6.1,'' -> Trace en indices de grille'')')XOMEGAY
                         ENDIF
                       ENDIF
                     ENDIF
                   ENDIF
!........
                   IF(LSPX)THEN
                     WRITE(INUM,'(''NBVAL en I '',i4 )')SIZE(ZTEM1D,1)
                   ELSE
                     WRITE(INUM,'(''NBVAL en J '',i4 )')SIZE(ZTEM1D,1)
                   ENDIF

                   IF(SIZE(XVAR,5) < 2)THEN

                      WRITE(INUM,'(36(''*''))')
                      IF(LSPX)THEN
                        WRITE(INUM,'(10X,''X(I)'',9X,''Y(VAL.R)'')')
                      ELSE
                        WRITE(INUM,'(10X,''X(J)'',9X,''Y(VAL.R)'')')
                      ENDIF
                      WRITE(INUM,'(36(''*''))')
                      DO J=1,SIZE(ZTEM1D,1)
                        WRITE(INUM,'(I4,2X,F8.1,(5X,E15.8))')J,ZTEM1D(J),ZWORKZ(J) 
                      ENDDO
                    ELSE
                      WRITE(INUM,'(55(''*''))')
                      IF(LSPX)THEN
                        WRITE(INUM,'(10X,''X(=I)'',8X,''Y(VAL.R)'',11X,''Y(VAL.Im)'')')
                      ELSE
                        WRITE(INUM,'(10X,''X(=J)'',8X,''Y(VAL.R)'',11X,''Y(VAL.Im)'')')
                      ENDIF
                      WRITE(INUM,'(55(''*''))')
                      DO J=1,SIZE(ZTEM1D,1)
                        IF(LSPX)THEN
                          WRITE(INUM,'(I4,2X,F8.1,2(5X,E15.8))')J,ZTEM1D(J),XVAR(J,1,NLOOPK,NLOOPT,1,NLOOPP),ZWORKZ(J) 
                        ELSE
                          WRITE(INUM,'(I4,2X,F8.1,2(5X,E15.8))')J,ZTEM1D(J),XVAR(1,J,NLOOPK,NLOOPT,1,NLOOPP),ZWORKZ(J) 
                        ENDIF
                      ENDDO
                   ENDIF
                   IF(SIZE(XVAR,5) < 2)THEN

                      WRITE(INUM,'(36(''*''))')
                   ELSE
                      WRITE(INUM,'(55(''*''))')
                   ENDIF

                ENDIF
!!!!!!!!!!!!!!!!!Mai 2002!!!!!!!!!!!!!!!!!!!!!!!!

	      ENDDO

!+++++++++ Boucle temps +++++++++++++++++++++++++++++++++++

	    ENDIF                                   !TTTTTTTTTTTTTTTTTTTTTT

	  ENDDO
!+++++++++ Boucle sur K +++++++++++++++++++++++++++++++++++
	  ENDDO
!+++++++++ Boucle sur N +++++++++++++++++++++++++++++++++++
	  ENDDO

!+++++++++ Boucle processus +++++++++++++++++++++++++++++++++++

	  DEALLOCATE(ZWORKZ,ZTEM1D,ZTEMLO)

!*************************************************************************
	  ENDIF
!!!!! BIDIMENSIONNELS 
!*************************************************************************
! Plan ( horizontal ou vertical // X ou vertical // Y)
!*************************************************************************
	ELSE                 !iiiiiiiiiiiiiiiiiiiiiiiiii
          print *,' bidimensionnel: II,IJ,IK=',II,IJ,IK


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CH Positionnement NIINF, NJINF, NISUP, NJSUP
! Defaut : NIINF=MAX(IIB,NIL), NJINF=MAX(IJB,NJL), NISUP=MIN(IIE,NIH), 
!          NJSUP=MIN(IJE,NJH)
! Sinon valeurs fournies par l'utilisateur dans les limites (NIL,NJL NIH,
! NJH)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CV Positionnement NIINF, NJINF, NISUP, NJSUP
! CV Positionnement LHORIZ et LVERTI
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		CALL RESOLV_NIJINF_NIJSUP
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CH + CV Allocation matrice 3D de reception des valeurs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            	ALLOCATE (ZWORK3D(1:NISUP-NIINF+1,1:NJSUP-NJINF+1, &
    		                  1:NKH-NKL+1))

!       	print *,' NBPROCDIA(KLOOP) ',NBPROCDIA(KLOOP)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle externe sur les processus
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        	DO JLOOPP=1,NBPROCDIA(KLOOP)
                  NLOOPP=NPROCDIA(JLOOPP,KLOOP)

		  CALL LOADUNITIT(JLOOPP,KLOOP)


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les numeros de R + I
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!       	print *,' NBNDIA(KLOOP) ',NBNDIA(KLOOP)

        	  DO JLOOPN=1,NBNDIA(KLOOP)
		    NLOOPN=NNDIA(JLOOPN,KLOOP)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  temps (Formulation sequentielle)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        	    IF(.NOT.LTINCRDIA(KLOOP,1))THEN

!       	      print *,' NBTIMEDIA(KLOOP,1) ',NBTIMEDIA(KLOOP,1)

        	      DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		        NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)


		        CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))

                        ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
		                       NJINF-NJL+1:NJSUP-NJL+1, &
        	                     :,NTIMEDIA(JLOOPT,KLOOP,1),JLOOPN, &
				     NPROCDIA(JLOOPP,KLOOP))
!                      WRITE(CLEGEND2(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
                       WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
!!!!!!!!!!!!!!!!!!!!!!!!!    CH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        	        IF(LCH)THEN

          		  IF(NBLVLKDIA(KLOOP,1) == 0)THEN

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  altitudes Z (Formulation sequentielle)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			    IF(.NOT.LZINCRDIA(KLOOP))THEN
          		    DO JLOOPZ=1,NBLVLZDIA(KLOOP)

			      IZ=XLVLZDIA(JLOOPZ,KLOOP)
          		      CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)

                              IF(KLOOP == NSUPERDIA)CALL FRAME

          		    ENDDO

			    ELSE

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  altitudes Z (Formulation incrementale)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          		    DO JLOOPZ=INT(XLVLZDIA(1,KLOOP)),INT(XLVLZDIA(2,KLOOP)), &
				      INT(XLVLZDIA(3,KLOOP))
			      IZ=JLOOPZ
          		      CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
  			        IF(KLOOP == NSUPERDIA)CALL FRAME

          		    ENDDO

			    ENDIF
        
          		  ELSE
        
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les niveaux de modele (Formulation sequentielle)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          		    DO JLOOPK=1,NBLVLKDIA(KLOOP,1)
          		      CALL TRACEH_FORDIACHRO(NLVLKDIA(JLOOPK, &
						     KLOOP,1),ZWORK3D,KLOOP)
  			        IF(KLOOP == NSUPERDIA)CALL FRAME
          		    ENDDO

          		  ENDIF
!                         CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)


!!!!!!!!!!!!!!!!!!!!!!!!!    CV    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			ELSE IF(LCV)THEN

			  IF(NLMAX <= 1 .OR. (NLANGLE<0 .OR. NLANGLE>360) .OR. &
			  (NIDEBCOU <=0 .AND. XIDEBCOU == -999.) .OR. &
			  (NJDEBCOU <=0 .AND. XJDEBCOU == -999.))THEN
			    PRINT *,' DEFINISSEZ D''ABORD NIDEBCOU, NJDEBCOU,',&
&                           ' NLMAX, NLANGLE (Pour CV + PV), PROFILE (Pour PV)'
                            PRINT *,'                  ou XIDEBCOU, XJDEBCOU'
			    PRINT *,' PUIS RENTREZ A NOUVEAU VOTRE DIRECTIVE '
!                           print *,' (Pour le 1D, mettre Obligatoirement ',&
!&                           'NLMAX=2 et LPOINTG=T'
			    PRINT *,' VALEURS ACTUELLES: '
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                           I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLMAX,NLANGLE,NPROFILE
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
                          ELSE
			    PRINT *,' VALEURS DES PARAMETRES DE DEFINITION DE LA COUPE',&
&                           ' ou DU PROFIL :'
			    IF(XIDEBCOU == -999. .AND. XJDEBCOU == -999.)THEN
			      PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                             NJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            ELSE
			      PRINT '('' XIDEBCOU:'',F7.1,'' XJDEBCOU:'',F7.1,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',XIDEBCOU, &
&                             XJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    ENDIF
                          ENDIF

			  CALL VERIFLEN_FORDIACHRO
		          ALLOCATE (ZTEMCV(NLMAX,1:IKU))
			  CALL PRECOU_FORDIACHRO(ZWORK3D,ZTEMCV)
			    ILENT=LEN_TRIM(CTITGAL)
			    ILENU=LEN_TRIM(CUNITGAL)
			    YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
			    YTEXTE(ILENT+1:ILENT+1)=' '
			    YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
                            CALL TRACEV_FORDIACHRO(ZTEMCV,KLOOP,YTEXTE(1: &
                            LEN_TRIM(YTEXTE)))
                              IF(KLOOP == NSUPERDIA)CALL FRAME

			  DEALLOCATE(ZTEMCV)
			  DEALLOCATE(XWORKZ,XWZ)

        	        ENDIF
        	      ENDDO
        
        	    ELSE
        
!       	      print *,' NBTIMEDIA(KLOOP,1) ',NBTIMEDIA(KLOOP,1)
!       	      print *,' NTIMEDIA(1 et 2,KLOOP,1) ',NTIMEDIA(1,KLOOP,1), &
!                     NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  temps (Formulation incrementale)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        	      DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1), &
						  NTIMEDIA(3,KLOOP,1)
                        NLOOPT=JLOOPT


		        CALL RESOLV_TIMES(JLOOPT)
        	          ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1,    &
				     NJINF-NJL+1:NJSUP-NJL+1,    &
        	                     :,JLOOPT,JLOOPN,NPROCDIA(JLOOPP,KLOOP))
!                       WRITE(CLEGEND2(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
                        WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)

!!!!!!!!!!!!!!!!!!!!!!!!!    CH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        	        IF(LCH)THEN

          		  IF(NBLVLKDIA(KLOOP,1) == 0)THEN

			    IF(.NOT.LZINCRDIA(KLOOP))THEN
          		    DO JLOOPZ=1,NBLVLZDIA(KLOOP)
                              IZ=XLVLZDIA(JLOOPZ,KLOOP)
          		      CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
  			        IF(KLOOP == NSUPERDIA)CALL FRAME
          		    ENDDO

			    ELSE

			    DO JLOOPZ=INT(XLVLZDIA(1,KLOOP)),INT(XLVLZDIA(2,KLOOP)), &
                                      INT(XLVLZDIA(3,KLOOP))
                              IZ=JLOOPZ
                              CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
  			        IF(KLOOP == NSUPERDIA)CALL FRAME
                            ENDDO
			    ENDIF
        
          		  ELSE
        
          		    DO JLOOPK=1,NBLVLKDIA(KLOOP,1)

          		      CALL TRACEH_FORDIACHRO(NLVLKDIA(JLOOPK,KLOOP,1), &
						     ZWORK3D,KLOOP)
  			        IF(KLOOP == NSUPERDIA)CALL FRAME
          		    ENDDO

          		  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!    CV    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			ELSE IF(LCV)THEN

			  IF(NLMAX <= 1 .OR. (NLANGLE<0 .OR. NLANGLE>360) .OR. &
			  (NIDEBCOU <=0 .AND. XIDEBCOU == -999.) .OR. &
			  (NJDEBCOU <=0 .AND. XJDEBCOU == -999.))THEN
			    PRINT *,' DEFINISSEZ D''ABORD NIDEBCOU, NJDEBCOU,',&
&                           ' NLMAX, NLANGLE (Pour CV + PV), PROFILE (Pour PV)'
                            PRINT *,'                  ou XIDEBCOU, XJDEBCOU'
			    PRINT *,' PUIS RENTREZ A NOUVEAU VOTRE DIRECTIVE '
!                           print *,' (Pour le 1D, mettre Obligatoirement ',&
!&                           'NLMAX=2 et LPOINTG=T'
			    PRINT *,' VALEURS ACTUELLES: '
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                           I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLMAX,NLANGLE,NPROFILE
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
                          ELSE
			    PRINT *,' VALEURS DES PARAMETRES DE DEFINITION DE LA COUPE',&
&                           ' ou DU PROFIL :'
			    IF(XIDEBCOU == -999. .AND. XJDEBCOU == -999.)THEN
			      PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                             NJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
			    ELSE
			      PRINT '('' XIDEBCOU:'',F7.1,'' XJDEBCOU:'',F7.1,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',XIDEBCOU, &
&                             XJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    ENDIF
                          ENDIF

			  CALL VERIFLEN_FORDIACHRO
		          ALLOCATE (ZTEMCV(NLMAX,1:IKU))
			  CALL PRECOU_FORDIACHRO(ZWORK3D,ZTEMCV)
!                         CALL IMCOU_FORDIACHRO(ZTEMCV,XDIAINT,CLEGEND,YTEXTE( &
!                         1:LEN_TRIM(YTEXTE)))
			    ILENT=LEN_TRIM(CTITGAL)
			    ILENU=LEN_TRIM(CUNITGAL)
			    YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
			    YTEXTE(ILENT+1:ILENT+1)=' '
			    YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
                            CALL TRACEV_FORDIACHRO(ZTEMCV,KLOOP,YTEXTE(1: &
                            LEN_TRIM(YTEXTE)))
    			      IF(KLOOP == NSUPERDIA)CALL FRAME
			  DEALLOCATE(ZTEMCV)
			  DEALLOCATE(XWORKZ,XWZ)

        	        ENDIF
        	      ENDDO
        	    ENDIF
        	  ENDDO
        	ENDDO

        ENDIF

!*****************************************************************************
!*****************************************************************************
!------------------------------------------------------------------------------
RETURN
END SUBROUTINE SUBSPXY
