!     ######spl
      MODULE  MODI_CLOSF
!     ##################
!
INTERFACE
!
SUBROUTINE CLOSF(KLOOPT,KTIMEND,KSEGD,KSEGM,K)
INTEGER,INTENT(IN)          :: KLOOPT,KTIMEND,KSEGD,KSEGM,K
END SUBROUTINE CLOSF
!
END INTERFACE
!
END MODULE MODI_CLOSF
!     ######spl
      SUBROUTINE CLOSF(KLOOPT,KTIMEND,KSEGD,KSEGM,K)
!     ##############################################
!
!!****  *CLOSF* - 
!!
!!    PURPOSE
!!    -------
!      
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
!!      Module
!!
!!      Module
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
!!      Original       24/11/95
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_MEMCV
USE MODD_NMGRID
USE MODD_COORD 
USE MODD_DEFCV
USE MODD_CONF
USE MODD_CTL_AXES_AND_STYL
USE MODN_NCAR
USE MODN_PARA
USE MODD_ALLOC_FORDIACHRO
USE MODD_TIME
USE MODD_TIME1
USE MODD_GRID1
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODE_GRIDPROJ

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER,INTENT(IN)          :: KLOOPT,KTIMEND,KSEGD,KSEGM,K
!
!*       0.1   Local variables
!              ---------------
!
INTEGER :: JJ, IER, INB, IWK, I, IA, ID, J
INTEGER :: IP, IN, IT, IZ, IPV=0
INTEGER :: KLEN, JI, JIM, ICOLI
INTEGER :: INUM, IRESP, ISEGM, ICOLSEGM
INTEGER :: II, IJ
INTEGER,SAVE :: INBTRACECV=0, INBTOT, ICO
LOGICAL,SAVE :: GGEOG, GVPTUSER
REAL    :: ZVPTL, ZVPTR, ZVPTB, ZVPTT
REAL    :: ZZZXD, ZZZXF, ZZZYD, ZZZYF
REAL    :: ZVL, ZVR, ZVB, ZVT, ZWL, ZWR, ZWB, ZWT
REAL    :: PHA, ZAX, ZAY, ZAU, ZAV
REAL    :: ZWIDTH, ZLAT, ZLON
REAL,DIMENSION(100) :: ZX, ZY
CHARACTER(LEN=25) :: CAR1, CAR2, CAR
CHARACTER(LEN=80) :: YTEM
!
!------------------------------------------------------------------------------
!IF(LANIMT)THEN     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL GFLAS2
!  IF(KLOOPT == KTIMEND)THEN
!    DO JJ=KSEGD,KSEGM
!      CALL GFLAS3(JJ)
!    ENDDO
!    CALL GCLWK(9)
!    CALL NGPICT(1,1)
!				  !!!!!!!!!!!!
!    CALL GQACWK(1,IER,INB,IWK)
!    IF(INB > 1)CALL NGPICT(2,3)
! ENDIF
!ELSE IF(LPXT .OR. LPYT)THEN     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(LPXT .OR. LPYT)THEN     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(LANIMK)THEN                !!LANIMK
    CALL GFLAS2
    IF(NBLVLKDIA(K,1) == 0)THEN
! Alt Niv. PR ou TK ...
      IF(.NOT.LZINCRDIA(K))THEN
! Pas incremental
        IF(XLOOPZ == XLVLZDIA(NBLVLZDIA(K),K))THEN
          DO JJ=KSEGD,KSEGM
            CALL GFLAS3(JJ)
          ENDDO
          CALL GCLWK(9)
          CALL NGPICT(1,1)
      				  !!!!!!!!!!!!
          CALL GQACWK(1,IER,INB,IWK)
          IF(INB > 1)CALL NGPICT(2,3)
	ENDIF
      ELSE
! Incremental
        IF(XLOOPZ == XLVLZDIA(2,K))THEN
          DO JJ=KSEGD,KSEGM
            CALL GFLAS3(JJ)
          ENDDO
          CALL GCLWK(9)
          CALL NGPICT(1,1)
      				  !!!!!!!!!!!!
          CALL GQACWK(1,IER,INB,IWK)
          IF(INB > 1)CALL NGPICT(2,3)
	ENDIF
      ENDIF
    ELSE
! Niveaux K
      IF(NLOOPK == NBLVLKDIA(K,1))THEN
          DO JJ=KSEGD,KSEGM
            CALL GFLAS3(JJ)
          ENDDO
          CALL GCLWK(9)
          CALL NGPICT(1,1)
      				  !!!!!!!!!!!!
          CALL GQACWK(1,IER,INB,IWK)
          IF(INB > 1)CALL NGPICT(2,3)
      ENDIF
    ENDIF

  ELSE                          !!LANIMK

  IF(K == NSUPERDIA)THEN        !+++++++++++++++++++++++++++++++++
! Trace du domaine fils eventuellement
    IF(LDOMAIN .AND. (LCH .OR. LCHXY ) .AND. .NOT.LCV)THEN
      ZZZXD=XXX(NDOMAINL,NMGRID)
      ZZZXF=XXX(NDOMAINR,NMGRID)
      ZZZYD=XXY(NDOMAINB,NMGRID)
      ZZZYF=XXY(NDOMAINT,NMGRID)
      CALL GSLWSC(XLWDOMAIN)
      CALL FRSTPT(ZZZXD,ZZZYD)
      CALL VECTOR(ZZZXF,ZZZYD)
      CALL VECTOR(ZZZXF,ZZZYF)
      CALL VECTOR(ZZZXD,ZZZYF)
      CALL VECTOR(ZZZXD,ZZZYD)
    ENDIF
! Trace de segments eventuellement
    IF(LSEGM .AND. (LCH .OR. LCHXY) .AND. .NOT.LCV)THEN
      CALL GQPLCI(IER,ICOLI)
      DO J=1,NCOLSEGM
      !IF(.NOT.LCOLAREA .AND. .NOT.LCOLINE .AND. NCOLSEGMS(J) > 1)THEN
      IF(NCOLSEGMS(J) > 1)THEN
	CALL TABCOL_FORDIACHRO
	print *,' appel a TABCOL_FORDIACHRO pour le trace de polylines couleur'
	EXIT
      ENDIF
      ENDDO
      CALL GSLWSC(XLWSEGM)
      ISEGM=0
      DO J=1,SIZE(NSEGMS,1)
      ! Conversion en coordonnees conformes
        ZLAT=XSEGMS(J,1)
        ZLON=XSEGMS(J,2)
        IF (NSEGMS(J)==1) THEN           ! XSEGMS
          IF (XCONFSEGMS(J,1)==0. .AND. XCONFSEGMS(J,2)==0.) &
            CALL SM_XYHAT_S(XLATORI,XLONORI, &
                            ZLAT,ZLON,                 &
                            XCONFSEGMS(J,1),XCONFSEGMS(J,2))
        ELSE IF (NSEGMS(J)==-1) THEN     ! ISEGMS
          NSEGMS(J)=1
          II=MAX(MIN(INT(ZLAT),NIMAX+2*JPHEXT-1),1)
          IJ=MAX(MIN(INT(ZLON),NJMAX+2*JPHEXT-1),1)
          XCONFSEGMS(J,1)=XXX(II,NMGRID) +  &
                  (ZLAT-FLOAT(II))*(XXX(II+1,NMGRID) - XXX(II,NMGRID) )
          XCONFSEGMS(J,2)=XXY(IJ,NMGRID) + &
                  (ZLON-FLOAT(IJ))*(XXY(IJ+1,NMGRID) - XXY(IJ,NMGRID) )
        END IF
	IF(J == 1 .AND. NSEGMS(J) == 1)THEN       
!       IF((J == 1 .AND. NSEGMS(J) == 1) .OR. &
!          (J >1 .AND. NSEGMS(J) == 1 .AND. &
!          NSEGMS(J-1) == 0))THEN
!         IF(J > 1)CALL SFLUSH
	  ISEGM=ISEGM+1
	  ICOLSEGM=NCOLSEGMS(ISEGM)
          IF((LCOLAREA .OR. LCOLINE) .AND. ICOLSEGM > 1)THEN
	    print *,' Avec LCOLAREA=T ou LCOLINE=T , attention a la superposition des couleurs'
	    print *,' pour les segments preferez NCOLSEGMS= 0 ou 1 '
            !print *,' valeur trouvee: ',NCOLSEGMS(ISEGM),'FORCEE a 1 '
            !ICOLSEGM=1
          ENDIF
          CALL GSPLCI(ICOLSEGM)
          CALL GSTXCI(ICOLSEGM)
          CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
!       ELSE IF(J > 1 .AND. NSEGMS(J) == 1 .AND. &
        ELSE IF(J > 1 .AND. NSEGMS(J) == 1 )THEN   
	  IF(NSEGMS(J-1)== 1)THEN
          CALL VECTOR(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
          ELSE
            CALL SFLUSH
            ISEGM=ISEGM+1
            ICOLSEGM=NCOLSEGMS(ISEGM)
            IF((LCOLAREA .OR. LCOLINE) .AND. ICOLSEGM > 1)THEN
	      print *,' Avec LCOLAREA=T ou LCOLINE=T , attention a la superposition des couleurs'
	      print *,' pour les segments preferez NCOLSEGMS= 0 ou 1 '
            !print *,' valeur trouvee: ',NCOLSEGMS(ISEGM),'FORCEE a 1 '
            !ICOLSEGM=1
            ENDIF
            CALL GSPLCI(ICOLSEGM)
            CALL GSTXCI(ICOLSEGM)
            CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
          ENDIF
	ENDIF
      ENDDO
      CALL SFLUSH
      CALL GSPLCI(ICOLI)
      CALL GSTXCI(1)
    ENDIF
! Trace de la CV dans CH suivante(s) eventuellement
    IF(LTRACECV .AND. (LCH .OR. LCHXY) .AND. .NOT.LCV)THEN
      CALL GQLWSC(IER,ZWIDTH)
      CALL GSLWSC(XLWTRACECV)
      CALL GSMKSC(2.)
      ICOLSEGM=1
      CALL GSPLCI(ICOLSEGM)
      CALL GSTXCI(ICOLSEGM)
      DO J=1,SIZE(NSEGMS,1)-1
	IF(NSEGMS(J) == 2 .AND. NSEGMS(J+1) ==2)THEN       
                print *,'closf J=',J
          CALL GSMK(4)
          CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
          CALL GSMK(5)
          CALL GPM(1,XCONFSEGMS(J+1,1),XCONFSEGMS(J+1,2))
          CALL CURVED(XCONFSEGMS(J:J+1,1),XCONFSEGMS(J:J+1,2),2)
	ENDIF
      ENDDO
      CALL SFLUSH
      CALL GSLWSC(ZWIDTH)
      CALL GSTXCI(1)
    ENDIF
!   Fermeture du dessin ds les cas =/= PH UMVM
    IF(LCH .AND. LCV .AND. (LUMVM .OR. LUTVT))THEN
    IF(nverbia > 0)then
    print *,' ***closf NLMAX ',NLMAX
    print *,' XTEMCVU ',XTEMCVU
    print *,' XTEMCVV ',XTEMCVV
    endif

    ELSE
      IF(LANIMT)THEN
        CALL GFLAS2
        IF(KLOOPT == KTIMEND)THEN
          DO JJ=KSEGD,KSEGM
          CALL GFLAS3(JJ)
          ENDDO
          CALL GCLWK(9)
          CALL NGPICT(1,1)
				  !!!!!!!!!!!!
          CALL GQACWK(1,IER,INB,IWK)
          IF(INB > 1)CALL NGPICT(2,3)
        ENDIF
      ELSE
        CALL NGPICT(1,1)
        CALL GQACWK(1,IER,INB,IWK)
        IF(INB > 1)CALL NGPICT(2,3)
        if(nverbia == -10)then
          print *,' CCCCCLOSF FRAME'
        endif
      ENDIF
    ENDIF
!   Fermeture du dessin ds les cas =/= PH UV
    IF(LTRACECV .AND. LCV .AND..NOT.L1DT)THEN    !.............................
      INBTRACECV=INBTRACECV+1
	IF(LPV)THEN
	  IPV=IPV+1
	  ZX(IPV)=XDSX(NPROFILE,NMGRID)
	  ZY(IPV)=XDSY(NPROFILE,NMGRID)
	ENDIF
	IF(NVERBIA == 10)THEN
	  print *,' closf INBTRACECV ',INBTRACECV
        ENDIF

      IF(INBTRACECV == 1)THEN     !0000000000000000000000000000000
	IP=NBPROCDIA(NLOOPSUPER)
	IN=NBNDIA(NLOOPSUPER)
	IF(.NOT.LTINCRDIA(NLOOPSUPER,1))THEN
	  IT=NBTIMEDIA(NLOOPSUPER,1)
	ELSE
	  IT=(NTIMEDIA(2,NLOOPSUPER,1)-NTIMEDIA(1,NLOOPSUPER,1))&
	  /NTIMEDIA(3,NLOOPSUPER,1)+1
	ENDIF
!       IF(LVLKDIALL(NLOOPSUPER,1))THEN
!         NBLVLKDIA(NLOOPSUPER,1)=0
!         print *,' **closf LTRACECV=T LCV=T LCH=T NBLVLKDIA(NLOOPSUPER,1) remis a 0 pour eliminer LVLKDIALL=T'
!       ENDIF
	IF(.NOT.LCH)THEN
!         print *,' **closf LCH CDIRCUR ',LCH,CDIRCUR(1:LEN_TRIM(CDIRCUR))
	  IZ=1
        ELSE
	IF(NBLVLKDIA(NLOOPSUPER,1) /= 0)THEN
	  IZ=NBLVLKDIA(NLOOPSUPER,1)
	ELSE
	IF(.NOT.LZINCRDIA(NLOOPSUPER))THEN
	  IZ=NBLVLZDIA(NLOOPSUPER)
	ELSE
	  IZ=(XLVLZDIA(2,NLOOPSUPER)-XLVLZDIA(1,NLOOPSUPER))&
	  /XLVLZDIA(3,NLOOPSUPER)+1
	ENDIF
	ENDIF
	ENDIF
	INBTOT=IP*IN*IT*IZ
	IF(NVERBIA == 10)THEN
	  print *,' closf INBTOT,IP,IN,IT,IZ ',INBTOT,IP,IN,IT,IZ
	ENDIF
      ENDIF                       !0000000000000000000000000000000

      IF(INBTRACECV == INBTOT)THEN
	    IF(LVPTUSER)THEN
	      GVPTUSER=.TRUE.
	      ZVPTL=XVPTL; ZVPTR=XVPTR; ZVPTB=XVPTB; ZVPTT=XVPTT
	    ELSE
	      GVPTUSER=.FALSE.
	    LVPTUSER=.TRUE.
	      ZVPTL=XVPTL; ZVPTR=XVPTR; ZVPTB=XVPTB; ZVPTT=XVPTT
	    XVPTL=.10; XVPTR=.90; XVPTB=.10; XVPTT=.90
	    ENDIF
	
          IF(LCARTESIAN)THEN
	    CALL DEFENETRE
          ELSE
	    IF(LGEOG)THEN
	      GGEOG=.TRUE.
	    ELSE
	      GGEOG=.FALSE.
	    ENDIF
!           LGEOG=.TRUE.
!           XVPTL=.12; XVPTR=.88; XVPTB=.12; XVPTT=.88
            CALL BCGRD_FORDIACHRO(1)
            CALL BCGRD_FORDIACHRO(2)
          ENDIF

          CALL GSLWSC(XLWTRACECV)
          CALL GSMKSC(2.)
          DO I =1,NTRACECV
            CALL GSMK(4)
	    CALL GPM(1,XTRACECV(1,I),XYTRACECV(1,I))
            CALL GSMK(5)
	    CALL GPM(1,XTRACECV(2,I),XYTRACECV(2,I))
!           CALL FRSTPT(XTRACECV(1,I),XYTRACECV(1,I))
!           CALL VECTOR(XTRACECV(2,I),XYTRACECV(2,I))
	    CALL CURVED(XTRACECV(1:2,I),XYTRACECV(1:2,I),2)
	    IF(IPV /= 0)THEN
	      DO IA=1,IPV
                CALL GSMKSC(1.)
		CALL GSMK(5)
		CALL GPM(1,ZX(IA),ZY(IA))
	      ENDDO
	      IPV=0
	    ENDIF
          ENDDO
! Janv 2001
          CALL GSMKSC(1.)
! Janv 2001
          CAR(1:LEN(CAR))=' '
          CAR1(1:LEN(CAR1))=' '
          CAR2(1:LEN(CAR2))=' '
          IF(LDEFCV2CC)THEN

            IF(LDEFCV2LL)THEN
              WRITE(CAR,'(''Latitude,Longitude :'')')
              WRITE(CAR1,'(''('',F6.2,'','',F6.2,'')'')')XIDEBCVLL,XJDEBCVLL
              WRITE(CAR2,'(''('',F6.2,'','',F6.2,'')'')')XIFINCVLL,XJFINCVLL
            ELSE IF(LDEFCV2IND)THEN
              WRITE(CAR,'(''Indices de grille I,J : '')')
              WRITE(CAR1,'(''('',I4,'','',I4,'')'')')NIDEBCV,NJDEBCV
              WRITE(CAR2,'(''('',I4,'','',I4,'')'')')NIFINCV,NJFINCV
            ELSE IF(LDEFCV2)THEN
              WRITE(CAR,'(''Coordonnees conformes : '')')
              WRITE(CAR1,'(''('',F10.2,'','',F10.2,'')'')')XTRACECV(1,1),XYTRACECV(1,1)
              WRITE(CAR2,'(''('',F10.2,'','',F10.2,'')'')')XTRACECV(2,1),XYTRACECV(2,1)
            ENDIF
          ELSE
	    IF(XIDEBCOU == -999.)THEN
	      WRITE(CAR,'(''Indices de grille I,J : '')')
	      WRITE(CAR1,'(''('',I4,'','',I4,'')'')')NIDEBCOU,NJDEBCOU
              WRITE(CAR2,'(''(NLMAX='',I4,'',ANG='',I3,'')'')')NLMAX,NLANGLE
	    ELSE
              WRITE(CAR,'(''Coordonnees conformes : '')')
              WRITE(CAR1,'(''('',F10.2,'','',F10.2,'')'')')XIDEBCOU,NJDEBCOU
              WRITE(CAR2,'(''(NLMAX='',I4,'',ANG='',I3,'')'')')NLMAX,NLANGLE
	    ENDIF
          ENDIF
	  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
	  XCURVPTL=ZVL;XCURVPTR=ZVR;XCURVPTB=ZVB;XCURVPTT=ZVT
!
!  Traitement des PH UMVM ou UTVT . Condition
!  LTRACECV=T _CV__K_ (ou Z etc...) UMVM ou UTVT

          IF(LCH .AND. LCV .AND. (LUMVM .OR. LUTVT))THEN
	    CALL ECHELLEPH(KLEN,PHA)
	    CALL GQLWSC(IER,ZWIDTH)
	    IF(XLWV > 0.)THEN
	      CALL GSLWSC(XLWV)
	    ENDIF
	    JIM=0
	    IF(LCOLINE)THEN
	      print *,' PH couleur fleches ? 1=noir 2=rouge 3=vert 4=bleu ... '
	      read(5,*,ERR=10)ICO
	      CALL GSPLCI(ICO)
	      GO TO 20
	      10 CONTINUE
	      BACKSPACE 5
	      print *,' Mai 2000 PH vecteurs vent horizontal. Si LCOLINE=T, possibilite '
	      print *,' de mettre les fleches en couleur en fournissant un indice apres la requete '
	      print *,' En cas d''absence, elles restent en noir '
	      20 CONTINUE
	    ENDIF
	    DO JI=1,SIZE(XTEMCVU,1),NISKIP
	      JIM=JIM+1
	      ZAX=XDSX(JI,1)
	      ZAY=XDSY(JI,1)
	      ZAU=XTEMCVU(JI,1)
	      ZAV=XTEMCVV(JI,1)
	      CALL FLECHE(ZAX,ZAY,ZAU,ZAV,KLEN,PHA)
	    ENDDO
	    CALL SFLUSH
	    CALL GSLWSC(ZWIDTH)
	    CALL GSPLCI(1)
	    CALL GSTXCI(1)
	    if(nverbia > 0)then
	      print *,' ***closf JIM ',JIM,' NISKIP,SIZE(XTEMCVU,1) ',NISKIP,SIZE(XTEMCVU,1)
	    endif
	    IF(LPRINT)THEN
              CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
              IF(IRESP /= 0)THEN
                CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
                OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
                PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
              ENDIF
              WRITE(INUM,'(''CLOSF  '',''G:'',A16,4X,'' NBVAL:'',I5,'' NLMAX:'',I5,'' NISKIP:'',I5)')CGROUP,JIM,NLMAX,NISKIP
              WRITE(INUM,'(A70)')CDIRCUR
              IF(LDEFCV2CC)THEN
                IF(LDEFCV2)THEN
                  WRITE(INUM,'(''cc(deb)-(fin)=('',F8.0,'','',F8.0,'')-('',F8.0,'','',F8.0,'')'')')&
                  &XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
                ELSE IF(LDEFCV2LL)THEN
                  WRITE(INUM,'(''ll(deb)-(fin)=('',F8.4,'','',F8.4,'')-('',F8.4,'','',F8.4,'')'')')&
                  &XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
                ELSE IF(LDEFCV2IND)THEN
                  WRITE(INUM,'(''ij(deb)-(fin)=('',I4,'','',I4,'')-('',I4,'','',I4,'')'')')&
                  &NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
                ENDIF
              ELSE
                IF(XIDEBCOU /= -999.)THEN
                  WRITE(INUM,'(''xidebcou'',F8.0,'' xjdebcou'',F8.0,'' nlmax'',i5,'' nlangle'',i4)')&
                  &XIDEBCOU,XJDEBCOU,NLMAX,NLANGLE
                ELSE
                  WRITE(INUM,'(''nidebcou'',i4,'' njdebcou'',i4,'' nlmax'',i5,'' nlangle'',i4)')&
                  &NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE
                ENDIF
              ENDIF
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**CLOSF XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
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
              WRITE(INUM,'(1X,78(1H*))')
              WRITE(INUM,'(15X,''U'',17X,''V'',17X,''X'',17X,''Y'')')
!             WRITE(INUM,'(16X,''X'',19X,''Y'')')
              WRITE(INUM,'(1X,78(1H*))')
	      JIM=0
              DO JI=1,SIZE(XTEMCVU,1),NISKIP
		JIM=JIM+1
	        ZAX=XDSX(JI,1)
	        ZAY=XDSY(JI,1)
	        ZAU=XTEMCVU(JI,1)
	        ZAV=XTEMCVV(JI,1)
                WRITE(INUM,'(I5,4(2X,E15.8))')JIM,ZAU,ZAV,ZAX,ZAY
              ENDDO
              WRITE(INUM,'(1X,78(1H*))')
	    ENDIF
	    DEALLOCATE(XTEMCVU,XTEMCVV)
	  ENDIF

          CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
          CALL GSLWSC(2.)
          IF(LDATFILE)CALL DATFILE_FORDIACHRO
          YTEM(1:LEN(YTEM))=' '
          CALL GSLWSC(2.)
	  CALL GSTXFP(-13,0)
          CALL RESOLV_TIT('CTITT1',YTEM)
          IF(YTEM /= ' ')THEN
            CALL PLCHHQ(.001,.98,YTEM(1:LEN_TRIM(YTEM)),.012,0.,-1.)
          ELSE
            CALL PLCHHQ(.001,.98,CDIRCUR(1:LEN_TRIM(CDIRCUR)),.012,0.,-1.)
          ENDIF
	  CALL GSTXFP(-13,2)
          CALL GSLWSC(2.)
          YTEM(1:LEN(YTEM))=' '
          CALL RESOLV_TIT('CTITT2',YTEM)
          IF(YTEM /= ' ')THEN
            CALL PLCHHQ(.001,.95,YTEM(1:LEN_TRIM(YTEM)),.009,0.,-1.)
          ENDIF
          YTEM(1:LEN(YTEM))=' '
          CALL RESOLV_TIT('CTITT3',YTEM)
          IF(YTEM /= ' ')THEN
            CALL PLCHHQ(.001,.93,YTEM(1:LEN_TRIM(YTEM)),.009,0.,-1.)
          ENDIF
          YTEM(1:LEN(YTEM))=' '
          CALL RESOLV_TIT('CTITB1',YTEM)
          IF(YTEM /= ' ')THEN
            CALL PLCHHQ(.001,.001,YTEM(1:LEN_TRIM(YTEM)),.009,0.,-1.)
          ENDIF
          CALL PLCHHQ(.001,.04,CAR(1:LEN_TRIM(CAR)),.012,0.,-1.)
          CALL GSMK(4)
	  CALL GPM(1,.35,.04)
          CALL PLCHHQ(.401,.04,CAR1(1:LEN_TRIM(CAR)),.012,0.,-1.)
          CALL GSMK(5)
	  CALL GPM(1,.70,.04)
          IF(LDEFCV2CC)THEN
            CALL PLCHHQ(.751,.04,CAR2(1:LEN_TRIM(CAR)),.012,0.,-1.)
	  ELSE
            CALL PLCHHQ(.721,.04,CAR2(1:LEN_TRIM(CAR)),.012,0.,-1.)
	  ENDIF
          CALL FRAME
          CALL GSLWSC(1.)
          INBTRACECV=0
	  NTRACECV=0
          IF(GVPTUSER)THEN
            LVPTUSER=.TRUE.
            XVPTL=ZVPTL; XVPTR=ZVPTR; XVPTB=ZVPTB; XVPTT=ZVPTT
          ELSE
            LVPTUSER=.FALSE.
            XVPTL=ZVPTL; XVPTR=ZVPTR; XVPTB=ZVPTB; XVPTT=ZVPTT
          ENDIF
          IF(LCARTESIAN)THEN
	  ELSE
          IF(GGEOG)THEN
            LGEOG=.TRUE.
          ELSE
            LGEOG=.FALSE.
          ENDIF
          ENDIF
!    IF(LCARTESIAN)THEN
!    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!    ENDIF
      ENDIF  !000000000000000000000000000000000000000
    ELSE  !..........................................
      INBTRACECV=0
      NTRACECV=0
    ENDIF !..........................................
  ENDIF   !++++++++++++++++++++++++++++++++++++++++++

  ENDIF                         !!LANIMK

ENDIF   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RETURN
END SUBROUTINE CLOSF
