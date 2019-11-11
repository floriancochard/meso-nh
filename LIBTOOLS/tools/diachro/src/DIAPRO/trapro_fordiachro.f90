!     ######spl
      MODULE MODI_TRAPRO_FORDIACHRO
!     #############################
!
INTERFACE
!
SUBROUTINE TRAPRO_FORDIACHRO(PTEM1D,PWORKZ,KLOOP)
INTEGER    :: KLOOP
REAL,DIMENSION(:)  :: PTEM1D, PWORKZ
END SUBROUTINE  TRAPRO_FORDIACHRO
!
END INTERFACE
END MODULE MODI_TRAPRO_FORDIACHRO
!     ######spl
      SUBROUTINE TRAPRO_FORDIACHRO(PTEM1D,PWORKZ,KLOOP)
!     #################################################
!
!!****  *TRAPRO_FORDIACHRO* - Manager of the 1D vertical profile plots
!!
!!    PURPOSE
!!    -------
!!      Controls 1D vertical profiles of scalar or vector variables.
!!     The displayed variables may be either from the Meso-NH basic
!!     set or generic.
!
!!**  METHOD
!!    ------
!!      Arrays are allocated, interactive dialogue is performed, and
!!    a branching is made either on the 'basic Meso-NH' set section
!!    or on the "generic variable" section where calls are made to the
!!    tracing routine PRO1D.
!!     
!!    EXTERNAL
!!    --------
!!      VALNGRID  : retrieves NGRID, the grid indicator, for the current 
!!                  variable name
!!      PRO1D     : tracing routine for the 1D vertical profiles
!!      OPNGKS    : opens NCAR and GKS graphics
!!      CLSGKS    : closes NCAR and GKS graphics
!!      COMPCOORD : computes gridpoint locations, meshsizes and topography
!!                  for all the possible grids, and true altitude where
!!                  required.
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TITLE  : Declares heading variables for the plots (TRACE)
!!         CLEGEND:  Current plot heading title
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
!!      Module MODD_PARAMETERS :  Contains array border depths
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODN_NCAR
!!>>>>> DRAGOON NOTICE: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!>>>>>         Apparently not used
!!>>>>> DRAGOON NOTICE: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!
!!      Module MODD_CVERT:  Declares work arrays for vertical cross-sections
!!         XWORKZ   : working array for true altitude storage (all grids)
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
#ifdef NAGf95
USE F90_UNIX  ! for FLUSH and GETENV
#endif

USE MODD_TITLE
USE MODD_TIT
USE MODD_NMGRID
USE MODN_PARA
USE MODD_PARAMETERS
USE MODN_NCAR
USE MODD_COORD  
USE MODD_GRID1  
USE MODD_GRID  
USE MODD_CONF  
USE MODE_GRIDPROJ
USE MODD_SUPER  
USE MODD_OUT
USE MODD_DEFCV
USE MODD_TYPE_AND_LH
USE MODD_RESOLVCAR
USE MODD_CTL_AXES_AND_STYL
USE MODI_READMNMX_FT_PVKT
!
USE MODI_WRITEDIR
!
IMPLICIT NONE
!
!*     0.1    interface declarations
!
INTERFACE
      SUBROUTINE PRO1D_FORDIACHRO(KPRO,PPRO,PTAB,PTABMIN,PTABMAX,KXDOT,HLEGEND,HTEXT)
      INTEGER :: KPRO, KXDOT
      REAL    :: PTABMIN, PTABMAX
      REAL,DIMENSION(:) :: PTAB, PPRO
      CHARACTER(LEN=*)  :: HTEXT, HLEGEND
      END SUBROUTINE PRO1D_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE VALMNMX(PMIN,PMAX)
      REAL                :: PMIN, PMAX
      END SUBROUTINE VALMNMX
END INTERFACE
!
!*      0.1    Dummy arguments 
!
INTEGER    :: KLOOP
REAL,DIMENSION(:)  :: PTEM1D, PWORKZ
!
!*      0.2    local variables 
!
!
INTEGER           :: IKU, IKE, IKB
INTEGER           :: IMJ
INTEGER           :: I, IM, IMB, J
INTEGER           :: ILENT, ILENU, IENDTXT
INTEGER           :: ISTA, IER, INB, IWK, ICOL
INTEGER           :: IPROFILE, IKL, IKH
INTEGER           :: INUM, IRESP
INTEGER           :: ISTYL   !, ISTY
!
REAL,SAVE         :: ZMIN, ZMAX, ZMN, ZMX
REAL              :: ZBMIN, ZBMAX
REAL,SAVE         :: ZHMIN, ZHMAX
REAL,SAVE         :: ZMNM, ZMXM
REAL,SAVE         :: ZX, ZY, ZLAT, ZLON
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZLA, ZLO
#ifdef RHODES
INTEGER          :: ISTAF
#endif
!
CHARACTER(LEN=20) :: YNOM
CHARACTER(LEN=40):: YTEXTE
CHARACTER(LEN=80):: YCAR80
CHARACTER(LEN=20):: YCAR20
CHARACTER(LEN=3) :: YREP
!
!-------------------------------------------------------------------------------
!
!*      1.     PRELIMINARY CALCULATIONS
!              ------------------------
!
!!!!!!!!!!! 110797
IF(NLOOPSUPER == 1)THEN
ZHMIN=XHMIN; ZHMAX=XHMAX
ENDIF
!!!!!!!!!!! 110797
NGSLNP=0
IKU=NKMAX+2*JPVEXT
SELECT CASE(CTYPE)
  CASE('CART','MASK')
    IKB=1+JPVEXT
    IKE=IKU-JPVEXT
  CASE('SSOL','DRST','RAPL')
    IKB=1
    IKE=SIZE(PTEM1D)
    IKL=NKL
    IKH=NKH
    NKL=1
    NKH=SIZE(PTEM1D)
END SELECT
YTEXTE(1:LEN(YTEXTE)) = ' '
ILENT=LEN_TRIM(CTITGAL)
ILENU=LEN_TRIM(CUNITGAL)
YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
YTEXTE(ILENT+1:ILENT+1)=' '
YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
IENDTXT=ILENT+2+ILENU
!
  IF(CTYPE == 'CART' .AND. .NOT.L1DT)THEN
    IF(.NOT.LCARTESIAN)THEN
      ALLOCATE(ZLA(NLMAX),ZLO(NLMAX))
      DO J=1,NLMAX
        ZX=XDSX(J,NMGRID)
        ZY=XDSY(J,NMGRID)
        CALL SM_LATLON_S(XLATORI,XLONORI,ZX,ZY,ZLAT,ZLON)
        ZLA(J)=ZLAT
        ZLO(J)=ZLON
      ENDDO
     if(nverbia > 0)then
     print *,' ZLA'
     print *,ZLA
     print *,' ZLO'
     print *,ZLO
     endif
      IF(LDEFCV2LL)THEN
	ZLA(1)=XIDEBCVLL
	ZLO(1)=XJDEBCVLL
      ENDIF
!     print *,' ZLA'
!     print *,ZLA
!     print *,' ZLO'
!     print *,ZLO
      XIPROFV=ZLA(NPROFILE); XJPROFV=ZLO(NPROFILE)
      if(nverbia > 0)then
      print *,' NPROFILE ZLA ZLO ',NPROFILE,ZLA(NPROFILE),ZLO(NPROFILE)
      endif
      DEALLOCATE(ZLA,ZLO)
    ENDIF
  ENDIF
IF(LPRINT)THEN
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
    WRITE(INUM,'(''TRAPRO  '',''G:'',A16,'' P:'',A25,'' '',A16,''  (1-IKU)''&
&   )')CGROUP,CTITGAL(1:25),ADJUSTL(CTIMEC)
  ELSE
    WRITE(INUM,'(''TRAPRO  '',''G:'',A16,'' P:'',A25,'' '',A16,'' NBVAL:'',&
&   I8)')CGROUP,CTITGAL(1:25),ADJUSTL(CTIMEC),SIZE(PTEM1D)
  ENDIF
  IF(LPLUS .OR.LMINUS)THEN
    WRITE(INUM,'(A70,A4)')CTITB3,CTYPE
  ELSE
    WRITE(INUM,'(A40,A4)')CTITGAL,CTYPE
  ENDIF
  IF(CTYPE == 'CART' .AND. .NOT.L1DT)THEN
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
!     print *,' ZLA'
!     print *,ZLA
!     print *,' ZLO'
!     print *,ZLO
      XIPROFV=ZLA(NPROFILE); XJPROFV=ZLO(NPROFILE)
      print *,' NPROFILE ZLA ZLO ',NPROFILE,ZLA(NPROFILE),ZLO(NPROFILE)
      DEALLOCATE(ZLA,ZLO)
    ENDIF
    IF(LDEFCV2CC)THEN
      IF(LDEFCV2)THEN
        WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'','' nlmax'',i5,&
&     '' profile'',i4)')XIDEBCV,XJDEBCV,XIFINCV,XJFINCV,NLMAX,NPROFILE
      ELSE IF(LDEFCV2LL)THEN
        WRITE(INUM,'(''ll(deb)-(fin)=('',F8.3,'','',F8.3,'')-('',F8.3,'','',F8.3,'')'','' nlmax'',i5,&
&      '' profile'',i4)')XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL,NLMAX,NPROFILE
      ELSE IF(LDEFCV2IND)THEN
        WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'','' nlmax'',i5,&
&       '' profile'',i4)')NIDEBCV,NJDEBCV,NIFINCV,NJFINCV,NLMAX,NPROFILE
      ENDIF
    ELSE
      IF(XIDEBCOU /= -999.)THEN
        WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4,&
   &   '' profile'',i4)')XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE,NPROFILE
      ELSE
        WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4,&
    &   '' profile'',i4)')NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,NPROFILE
      ENDIF
    ENDIF
  ENDIF
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**TRAPRO XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
    ELSE
      WRITE(INUM,'(1X,75(''*''))')
      WRITE(INUM,'(1X,''    Dates courante   *     modele      *   experience    *      segment'')')
      WRITE(INUM,'(1X,'' J   An  M  J  Sec.  * An  M  J  Sec.  * An  M  J  Sec.  * An  M  J  Sec.'')')
      WRITE(INUM,'(1X,75(''*''))')
      DO J=1,SIZE(XPRDAT,2)
        WRITE(INUM,'(1X,I3,1X,3(I4,I3,I3,I6,'' *''),I4,I3,I3,I6)')J,INT(XPRDAT(:,J))

      ENDDO
    ENDIF
  ENDIF
! JUin 2001 Ecriture des dates 
  WRITE(INUM,'(1X,45(''*''))')
  WRITE(INUM,'(''   K'',12X,''X'',19X,''Z'',"          NBVAL:",I6)')SIZE(PTEM1D)
  WRITE(INUM,'(1X,45(''*''))')
  DO J=SIZE(PTEM1D),1,-1
    WRITE(INUM,'(I5,4X,E15.8,4X,E15.8)')J,PTEM1D(J),PWORKZ(J)
  ENDDO
ENDIF
!
!
!
!*      1.4    Interactive selection of the profile location, and
!*             field name
!
! Profile point selection
!
IF(.NOT. L1DT)THEN

IF(NPROFILE.GT.NLMAX)THEN
  PRINT *,' This point ',NPROFILE,'  lays out of the section limits..'
  PRINT *,' index has to be smaller than ',NLMAX
  PRINT *,' Enter the gridpoint location for the profile: '
  PRINT *,' i.e. the gridpoint index along the oblique vertical section '
  PRINT *,' starting at (NIDEBCOU,NJDEBCOU or XIDEBCOU,XJDEBCOU or .....)?'
  READ(5,*)NPROFILE
ENDIF

ELSE

  IPROFILE=NPROFILE
  NPROFILE=1

ENDIF
!
! Field name selection 
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
! print *,' Profile''s bounds (min and max )? '
! READ(5,*)ZBMIN,ZBMAX
ELSE
  ZBMIN=0.;ZBMAX=0.
END IF
!
!*       1.6   Ooverlay control
!
IF(NSUPERDIA > 1)THEN
  LSUPER=.TRUE.
ELSE
  LSUPER=.FALSE.
ENDIF
IF(KLOOP == 1)NSUPER=0
!
!
!*       1.8    Line width changes to differentiate the 
!*              successive plots in an overlay sequence 
!
CALL TABCOL_FORDIACHRO
CALL GSLWSC(2.)
CALL GSLN(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(LSUPER)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NSUPER=NSUPER+1
!----------------------------------------------------------
  IF(LCOLINE)THEN
!----------------------------------------------------------
    SELECT CASE(NSUPER)
      CASE(1)
	IF(XLWPV1 /= 0.)CALL GSLWSC(XLWPV1)
      CASE(2)
	IF(XLWPV2 /= 0.)CALL GSLWSC(XLWPV2)
      CASE(3)
	IF(XLWPV3 /= 0.)CALL GSLWSC(XLWPV3)
      CASE(4)
	IF(XLWPV4 /= 0.)CALL GSLWSC(XLWPV4)
      CASE(5)
	IF(XLWPV5 /= 0.)CALL GSLWSC(XLWPV5)
      CASE(6)
	IF(XLWPV6 /= 0.)CALL GSLWSC(XLWPV6)
      CASE(7)
	IF(XLWPV7 /= 0.)CALL GSLWSC(XLWPV7)
      CASE(8)
	IF(XLWPV8 /= 0.)CALL GSLWSC(XLWPV8)
!!!!!!!!!!!!!!!!!
      CASE(9)
	IF(XLWPV9 /= 0.)CALL GSLWSC(XLWPV9)
      CASE(10)
	IF(XLWPV10 /= 0.)CALL GSLWSC(XLWPV10)
      CASE(11)
	IF(XLWPV11 /= 0.)CALL GSLWSC(XLWPV11)
      CASE(12)
	IF(XLWPV12 /= 0.)CALL GSLWSC(XLWPV12)
      CASE(13)
	IF(XLWPV13 /= 0.)CALL GSLWSC(XLWPV13)
      CASE(14)
	IF(XLWPV14 /= 0.)CALL GSLWSC(XLWPV14)
      CASE(15)
	IF(XLWPV15 /= 0.)CALL GSLWSC(XLWPV15)
!!!!!!!!!!!!!!!!!
      CASE DEFAULT
	CALL GSLWSC(2.)
    END SELECT
!+++++++++++++++++++++++++++++++++
    IF(LCOLUSER)THEN
!+++++++++++++++++++++++++++++++++
      IF(NSUPER == 1)THEN
        print *,' VOUS VOULEZ VRAIMENT SELECTIONNER LES COULEURS DES PROFILS ? (y/n) '
        YREP='   '
        READ(5,'(A3)',END=10)YREP
	GO TO 20
	10 CONTINUE
        CLOSE(5)
        CALL GETENV("VARTTY",YCAR20)
	YCAR20=ADJUSTL(YCAR20)
        OPEN(5,FILE=YCAR20)
	READ(5,'(A3)')YREP
	20 CONTINUE
        YCAR80=YREP
        !WRITE(NDIR,'(A80)')YCAR80
        CALL WRITEDIR(NDIR,YCAR80)
#ifdef RHODES
	CALL FLUSH(NDIR,ISTAF)
#else
	CALL FLUSH(NDIR)
#endif
          IF(YREP /= 'y' .AND. YREP /= 'yes' .AND. YREP /= 'Y'   &
	  .AND. YREP /= 'YES' &
    	  .AND. YREP /= 'o' .AND. YREP /= 'oui' .AND. YREP /= 'O' &
    	  .AND. YREP /= 'OUI')THEN
            LCOLUSER=.FALSE.
	    print *,' LCOLUSER REMIS A .FALSE.'
            YCAR20(1:LEN(YCAR20))=' '
            INQUIRE(5,NAME=YCAR20)
            print *,' AP INQUIRE YCAR20 ',YCAR20
            YCAR20=ADJUSTL(YCAR20)
            IF(YCAR20(1:8) /= '/dev/tty')BACKSPACE 5
            CALL GSLN(1)
            CALL GSPLCI(NSUPER+1)
          ELSE
    	    print *,' INDICE DE COULEUR POUR ',CTITGAL(1:ILENT),' ? '
    	    READ(5,*,END=11)ICOL
	    GO TO 21
	    11 CONTINUE
            CLOSE(5)
            CALL GETENV("VARTTY",YCAR20)
	    YCAR20=ADJUSTL(YCAR20)
            OPEN(5,FILE=YCAR20)
	    READ(5,*)ICOL
	    21 CONTINUE
	    WRITE(YCAR80,*)ICOL
	    YCAR80=ADJUSTL(YCAR80)
	    !WRITE(NDIR,'(A80)')YCAR80
            CALL WRITEDIR(NDIR,YCAR80)
#ifdef RHODES
	    CALL FLUSH(NDIR,ISTAF)
#else
	    CALL FLUSH(NDIR)
#endif
!           WRITE(NDIR,*)YCAR80
            CALL GSLN(1)
            CALL GSPLCI(ICOL)
          ENDIF
      ELSE
	print *,' INDICE DE COULEUR POUR ',CTITGAL(1:ILENT),' ? '
	READ(5,*,END=12)ICOL
        GO TO 22
        12 CONTINUE
        CLOSE(5)
        CALL GETENV("VARTTY",YCAR20)
        YCAR20=ADJUSTL(YCAR20)
        OPEN(5,FILE=YCAR20)
        READ(5,*)ICOL
        22 CONTINUE
	WRITE(YCAR80,*)ICOL
	YCAR80=ADJUSTL(YCAR80)
	!WRITE(NDIR,'(A80)')YCAR80
        CALL WRITEDIR(NDIR,YCAR80)
#ifdef RHODES
	CALL FLUSH(NDIR,ISTAF)
#else
	CALL FLUSH(NDIR)
#endif
!       WRITE(NDIR,*)YCAR80
        CALL GSLN(1)
        CALL GSPLCI(ICOL)
      ENDIF
!+++++++++++++++++++++++++++++++++
    ELSE
!+++++++++++++++++++++++++++++++++
      CALL GSLN(1)
      CALL GSPLCI(NSUPER+1)
!+++++++++++++++++++++++++++++++++
    ENDIF
!+++++++++++++++++++++++++++++++++
!----------------------------------------------------------
  ELSE
!----------------------------------------------------------
    CALL GSPLCI(1)
!   IF(NSUPER == 1)CALL GSLN(1)  ! Solid line : first in sequence
!   IF(NSUPER == 2)CALL GSLN(3)  ! Dotted line: second in sequence
!   IF(NSUPER == 3)CALL GSLN(2)  ! Dashed line
!   IF(NSUPER == 4)CALL GSLN(4)  ! Dashed line: further on
    CALL GSTXCI(1)
    SELECT CASE(NSUPER)
      CASE(:4)
	CALL GSLWSC(1.)
	IF(NSUPER == 1 .AND. XLWPV1 /= 0.)CALL GSLWSC(XLWPV1)
	IF(NSUPER == 2 .AND. XLWPV2 /= 0.)CALL GSLWSC(XLWPV2)
	IF(NSUPER == 3 .AND. XLWPV3 /= 0.)CALL GSLWSC(XLWPV3)
	IF(NSUPER == 4 .AND. XLWPV4 /= 0.)CALL GSLWSC(XLWPV4)
      CASE(5:8)
	CALL GSLWSC(2.)
	IF(NSUPER == 5 .AND. XLWPV5 /= 0.)CALL GSLWSC(XLWPV5)
	IF(NSUPER == 6 .AND. XLWPV6 /= 0.)CALL GSLWSC(XLWPV6)
	IF(NSUPER == 7 .AND. XLWPV7 /= 0.)CALL GSLWSC(XLWPV7)
	IF(NSUPER == 8 .AND. XLWPV8 /= 0.)CALL GSLWSC(XLWPV8)
      CASE(9:12)
	CALL GSLWSC(3.)
!!!!!!!!!!
	IF(NSUPER == 9 .AND. XLWPV9 /= 0.)CALL GSLWSC(XLWPV9)
	IF(NSUPER == 10 .AND. XLWPV10 /= 0.)CALL GSLWSC(XLWPV10)
	IF(NSUPER == 11 .AND. XLWPV11 /= 0.)CALL GSLWSC(XLWPV11)
	IF(NSUPER == 12 .AND. XLWPV12 /= 0.)CALL GSLWSC(XLWPV12)
!!!!!!!!!!
      CASE(13:16)
	CALL GSLWSC(4.)
!!!!!!!!!!
	IF(NSUPER == 13 .AND. XLWPV13 /= 0.)CALL GSLWSC(XLWPV13)
	IF(NSUPER == 14 .AND. XLWPV14 /= 0.)CALL GSLWSC(XLWPV14)
	IF(NSUPER == 15 .AND. XLWPV15 /= 0.)CALL GSLWSC(XLWPV15)
!!!!!!!!!!
      CASE DEFAULT
	CALL GSLWSC(1.)
    END SELECT
    NGSLNP=0
    IF(NSUPER == 1 .AND. XSTYLPV1 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV1))
      CALL GQLN(IER,ISTYL)
!     print *,' TRAPRO 1 INT(XSTYLPV1) ISTYL  ',INT(XSTYLPV1),ISTYL
      IF(INT(XSTYLPV1) >4)NGSLNP=XSTYLPV1-1
    ELSEIF(NSUPER == 2 .AND. XSTYLPV2 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV2))
      IF(INT(XSTYLPV2) >4)NGSLNP=XSTYLPV2-1
    ELSEIF(NSUPER == 3 .AND. XSTYLPV3 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV3))
      IF(INT(XSTYLPV3) >4)NGSLNP=XSTYLPV3-1
    ELSEIF(NSUPER == 4 .AND. XSTYLPV4 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV4))
      IF(INT(XSTYLPV4) >4)NGSLNP=XSTYLPV4-1
    ELSEIF(NSUPER == 5 .AND. XSTYLPV5 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV5))
      IF(INT(XSTYLPV5) >4)NGSLNP=XSTYLPV5-1
    ELSEIF(NSUPER == 6 .AND. XSTYLPV6 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV6))
      IF(INT(XSTYLPV6) >4)NGSLNP=XSTYLPV6-1
    ELSEIF(NSUPER == 7 .AND. XSTYLPV7 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV7))
      IF(INT(XSTYLPV7) >4)NGSLNP=XSTYLPV7-1
    ELSEIF(NSUPER == 8 .AND. XSTYLPV8 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV8))
      IF(INT(XSTYLPV8) >4)NGSLNP=XSTYLPV8-1
!!!!!!!!!!
    ELSEIF(NSUPER == 9 .AND. XSTYLPV9 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV9))
      IF(INT(XSTYLPV9) >4)NGSLNP=XSTYLPV9-1
    ELSEIF(NSUPER == 10 .AND. XSTYLPV10 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV10))
      IF(INT(XSTYLPV10) >4)NGSLNP=XSTYLPV10-1
    ELSEIF(NSUPER == 11 .AND. XSTYLPV11 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV11))
      IF(INT(XSTYLPV11) >4)NGSLNP=XSTYLPV11-1
    ELSEIF(NSUPER == 12 .AND. XSTYLPV12 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV12))
      IF(INT(XSTYLPV12) >4)NGSLNP=XSTYLPV12-1
    ELSEIF(NSUPER == 13 .AND. XSTYLPV13 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV13))
      IF(INT(XSTYLPV13) >4)NGSLNP=XSTYLPV13-1
    ELSEIF(NSUPER == 14 .AND. XSTYLPV14 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV14))
      IF(INT(XSTYLPV14) >4)NGSLNP=XSTYLPV14-1
    ELSEIF(NSUPER == 15 .AND. XSTYLPV15 /= 0.)THEN
      CALL GSLN(INT(XSTYLPV15))
      IF(INT(XSTYLPV15) >4)NGSLNP=XSTYLPV15-1
!!!!!!!!!!
    ELSE
      CALL GSLN(MOD(NSUPER,4))
      IF(MOD(NSUPER,4) == 0)CALL GSLN(4)
    ENDIF
!----------------------------------------------------------
  END IF
!----------------------------------------------------------
!     print *,' TRAPRO 1  entre ENDIF et ELSE INT(XSTYLPV1) ISTYL  ',INT(XSTYLPV1),ISTYL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CALL GSLWSC(1.)
  CALL GSLWSC(2.)
  IF(XLWPV1 /= 0.)CALL GSLWSC(XLWPV1)
  CALL GSLN(1)              ! Solid line if no overlay
! ISTY=NINT(XSTYLPV1)
! NGSLNP=0
! IF(XSTYLPV1 /= 0.)THEN
!     CALL GSLN(ISTY)
!     IF(INT(XSTYLPV1) >4)NGSLNP=XSTYLPV1-1
!     print *,' TRAPRO 2 INT(XSTYLPV1) ',INT(XSTYLPV1)
! ENDIF
! CALL GQLN(IER,ISTYL)
! print *,' TRAPRO 2 ISTYL ',ISTYL
!+++++++++++++++++++++++++++++++++
  IF(LCOLUSER)THEN
!+++++++++++++++++++++++++++++++++
	print *,' INDICE DE COULEUR ? '
	READ(5,*,END=82)ICOL
        GO TO 92
        82 CONTINUE
        CLOSE(5)
        CALL GETENV("VARTTY",YCAR20)
        YCAR20=ADJUSTL(YCAR20)
        OPEN(5,FILE=YCAR20)
        READ(5,*)ICOL
        92 CONTINUE
	WRITE(YCAR80,*)ICOL
	YCAR80=ADJUSTL(YCAR80)
        !WRITE(NDIR,'(A80)')YCAR80
        CALL WRITEDIR(NDIR,YCAR80)
#ifdef RHODES
	CALL FLUSH(NDIR,ISTAF)
#else
	CALL FLUSH(NDIR)
#endif
!       WRITE(NDIR,*)YCAR80
        CALL GSLN(1)
        CALL GSPLCI(ICOL)
!+++++++++++++++++++++++++++++++++
  ELSE
!+++++++++++++++++++++++++++++++++
    CALL GSPLCI(1)
    CALL GSTXCI(1)
!+++++++++++++++++++++++++++++++++
  ENDIF
!+++++++++++++++++++++++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------------------------
!
!*       2.     PROCESSING THE BASIC SET OF PROGNOSTIC VARIABLES
!               ------------------------------------------------
!
! On a NMGRID CTITGAL CUNITGAL 
!  TESTER NMGRID DANS OPER POUR LE METTRE A 1 SI IL A UNE VALEUR ABERRANTE
IF(XHMAX-XHMIN == 0.)THEN

SELECT CASE(CTYPE)
  CASE('CART','MASK')
    XHMIN=0.
!   XHMIN=PWORKZ(IKB)
    XHMAX=PWORKZ(IKE)
  CASE('SSOL','DRST','RAPL')
    XHMIN=MIN(0.,PWORKZ(IKB))
    XHMAX=MAX(0.,PWORKZ(IKE))
    IF(XHMIN == 0. .AND. XHMAX == 0.)THEN
      XHMIN=-1.
      XHMAX=1.
    ENDIF
END SELECT

ENDIF
!ZPRO(:)=XWORKZ(IPRO,1:IKU,NMGRID) --> PWORKZ
DO I=IKB,IKE
  IF(XHMAX.LT.PWORKZ(I))THEN
  IM=I
  EXIT
  ENDIF
  IM=I
ENDDO
IM=MIN(IM,IKE)
DO I=IKB,IKE
  IF(XHMIN <= PWORKZ(I))THEN
  IMB=MAX(I-1,IKB)
  EXIT
  ENDIF
  IMB=MAX(I-1,IKB)
ENDDO
IF(NPVITVXMJ /= 0)THEN
  IMJ=NPVITVXMJ
ELSE
  IMJ=4
ENDIF
ZMN=0.; ZMX=0.
!IF(XPVMIN /=0. .OR. XPVMAX /=0.)THEN
!  ZMN=XPVMIN
!  ZMX=XPVMAX
!ENDIF
LOK=.FALSE.
if(nverbia > 0)then
  print *,' TRAPRO AP LOK=F LOK,LMNMXUSER ',LOK,LMNMXUSER
endif
IF(LMNMXUSER)THEN
!666666666666666666666666666666666666666666666666666666666666666666
  IF(XPVMAXT-XPVMINT /= 0)THEN
    LOK=.TRUE.
  ELSE
!666666666666666666666666666666666666666666666666666666666666666666
    print *,' TRAPRO ',CTITGAL
    CALL READMNMX_FT_PVKT(CTITGAL(1:LEN_TRIM(CTITGAL)),ZMN,ZMX)
    if(nverbia > 0)THEN
    print *,' TRAPRO ZMN ZMX ',ZMN,ZMX,LOK
    ENDIF
!666666666666666666666666666666666666666666666666666666666666666666
  ENDIF
!666666666666666666666666666666666666666666666666666666666666666666
ENDIF
ZMIN=MINVAL(PTEM1D(MAX(IMB,NKL):MIN(IM,NKH)))
ZMAX=MAXVAL(PTEM1D(MAX(IMB,NKL):MIN(IM,NKH)))

print *,' TRAPRO ZMIN ZMAX TROUVES ',ZMIN,ZMAX

!

  
SELECT CASE(NSUPER)
CASE(:1)
!66666666666666666666666666666666666666666666666666
      IF(LMNMXUSER .AND. LOK)THEN
        IF(XPVMAXT-XPVMINT /= 0)THEN
          print *,' TRAPRO XPVMINT,XPVMAXT UTILISES :',XPVMINT,XPVMAXT
          CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),XPVMINT,XPVMAXT,IMJ,CLEGEND,YTEXTE&
          (1:LEN_TRIM(YTEXTE)))
        ELSE
          CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),ZMN,ZMX,IMJ,CLEGEND,YTEXTE&
          (1:LEN_TRIM(YTEXTE)))
          ZMNM=ZMN; ZMXM=ZMX
	ENDIF
      ELSE
!66666666666666666666666666666666666666666666666666
        CALL VALMNMX(ZMIN,ZMAX)
        IF(ZMAX-ZMIN == 0.)THEN
          ZMIN=ZMIN-1.
          ZMAX=ZMAX+1.
        ENDIF
	  print *,' TRAPRO CALCUL AUTOMATIQUE DES BORNES: ',ZMIN,ZMAX
!         print *,' TRAPRO  av pro1d INT(XSTYLPV1) ',INT(XSTYLPV1)
!         CALL GQLN(IER,ISTYL)
!         print *,' TRAPRO  av pro1d ISTYL ',ISTYL
          CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),ZMIN,ZMAX,IMJ,CLEGEND,YTEXTE &
          (1:LEN_TRIM(YTEXTE)))
          ZMNM=ZMIN; ZMXM=ZMAX
	  if(nverbia > 0)then
          print *,' TRAPRO ap pro1d INT(XSTYLPV1) ',INT(XSTYLPV1)
	  endif
      ENDIF

!66666666666666666666666666666666666666666666666666
!66666666666666666666666666666666666666666666666666
CASE(2:)
!66666666666666666666666666666666666666666666666666
  IF(LMNMXUSER .AND. LOK)THEN
    IF(XPVMAXT-XPVMINT /= 0)THEN
      print *,' TRAPRO XPVMINT,XPVMAXT UTILISES :',XPVMINT,XPVMAXT
      CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),XPVMINT,XPVMAXT,IMJ,CLEGEND,YTEXTE&
      (1:LEN_TRIM(YTEXTE)))
    ELSE
      CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),ZMN,ZMX,IMJ,CLEGEND,YTEXTE &
      (1:LEN_TRIM(YTEXTE)))
    ENDIF
  ELSE
!66666666666666666666666666666666666666666666666666


    IF(ZMIN >=ZMNM .AND. ZMAX <= ZMXM)THEN
      CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),ZMNM,ZMXM,IMJ,CLEGEND,YTEXTE &
      (1:LEN_TRIM(YTEXTE)))
    ELSE
      CALL VALMNMX(ZMIN,ZMAX)
      IF(ZMAX-ZMIN == 0.)THEN
        ZMIN=ZMIN-1.
        ZMAX=ZMAX+1.
      ENDIF
      CALL PRO1D_FORDIACHRO(NPROFILE,PWORKZ(:),PTEM1D(:),ZMIN,ZMAX,IMJ,CLEGEND,YTEXTE &
      (1:LEN_TRIM(YTEXTE)))
    ENDIF
  END IF
!66666666666666666666666666666666666666666666666666
!66666666666666666666666666666666666666666666666666
END SELECT
SELECT CASE(CTYPE)
  CASE('SSOL','DRST','RAPL')
    NKL=IKL
    NKH=IKH
  if(nverbia > 0)then
  print *,' TRAPRO NKL NKH',NKL,NKH
  endif
END SELECT
IF(L1DT)THEN
  NPROFILE=IPROFILE
  if(nverbia > 0)then
  print *,' TRAPRO NPROFILE ',NPROFILE
  endif
ENDIF
!
1000 FORMAT(5X,I4,3X,A12)
!
!----------------------------------------------------------------------------
!
!*       4.     EXIT
!               ----
IF(.NOT.LSUPER  .OR.  (LSUPER .AND. NSUPER == NSUPERDIA))THEN
  XHMIN=ZHMIN; XHMAX=ZHMAX
ENDIF
if(nverbia > 0)then
print *,' TRAPRO SORTIE XSTYLPV1 ',XSTYLPV1
endif
RETURN
!
END SUBROUTINE  TRAPRO_FORDIACHRO
