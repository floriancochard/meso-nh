!     ######spl
      PROGRAM  DIAPROG
!     ################
!
!!****  *DIAPROG* - 
!! 
!!
!!    PURPOSE
!!    -------
! 
!
!!**  METHOD
!!    ------
!!      
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHORS
!!    -------
!!    J. Duron      *Lab. Aerologie* 
!!
!!    Copyright 1994,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/11/95 
!!      Updated  PM 23/11/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef NAGf95
USE F90_UNIX  ! for FLUSH and GETENV
#endif

USE MODD_CST
USE MODD_CONF, ONLY : CPROGRAM
USE MODD_MASK3D
USE MODD_COORD
USE MODD_TYPE_AND_LH
USE MODD_GRID
USE MODD_GRID1
USE MODD_ALLOC_FORDIACHRO
USE MODD_ALLOC2_FORDIACHRO
USE MODN_NCAR
USE MODN_PARA
USE MODD_OUT
USE MODD_NMGRID
USE MODD_FILES_DIACHRO
USE MODD_RESOLVCAR 
USE MODI_EXTRACT_AND_OPEN_FILES
USE MODI_READ_DIMGRIDREF
USE MODI_READ_DIACHRO
USE MODI_CARESOLV
USE MODI_CARMEMORY
USE MODI_LOAD_FMTAXES
USE MODI_LOAD_SEGMENTS
USE MODI_CONVLO2UP
USE MODI_OPER_PROCESS
USE MODI_PRINTS
USE MODI_REALLOC_AND_LOAD
USE MODI_ALLOC2_FORDIACHRO
USE MODI_DIFF_OPER
USE MODD_TIT
USE MODD_PVT
USE MODD_MEMCV
USE MODD_EXPR 
USE MODI_RESOLV_TIT
USE MODI_LOAD_TIT
USE MODI_CONVIJ2XY
USE MODD_SEVERAL_RECORDS
USE MODI_VERIF_GROUP
USE MODI_READ_UVW
USE MODI_READ_TYPE
USE MODD_PT_FOR_CH_FORDIACHRO
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
USE MODD_TRAJ3D
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
USE MODI_WRITEDIR

IMPLICIT NONE
!
!*       0.1   Local variables declarations
!
INTEGER           :: JI, JIA, JJ, J, JM, ITOP1, ITOP2
INTEGER           :: JLOOP, INDEXPR, IMULTDIV
INTEGER           :: INDE, ILENC240, ILENC
INTEGER           :: INDPRI, INDTIT, INDPRIL, IZERO
INTEGER           :: IDIR, IDIRESP, ITITDEF
INTEGER           :: ICONVIJ2XY, ICONVXY2IJ, ICONVALLIJ2LL
INTEGER           :: ICNOMCAR, ICSYMCAR
INTEGER           :: ICGROUPSV3, ILENT, IQUOT
#ifdef RHODES
INTEGER          :: ISTAF
#endif
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
INTEGER           :: ICTRAJ_GROUP
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!

REAL,DIMENSION(:),ALLOCATABLE :: ZBID2
REAL,DIMENSION(:,:),ALLOCATABLE :: ZBID1
CHARACTER(LEN=100)  :: CAR100
CHARACTER(LEN=80)   :: CAR80
CHARACTER(LEN=20)   :: CAR20, VARTTY
CHARACTER(LEN=2400) :: CAR240, YCAR240
CHARACTER(LEN=16)   :: YDIRNAM
CHARACTER(LEN=8)    :: YDAT
CHARACTER(LEN=10)   :: YTEM
CHARACTER(LEN=6)    :: YMULTDIV

INTEGER         :: ILINES=25  ! nb de lignes de directives avec &

LOGICAL           :: GMASK3D, GMASK3D_XY, GMASK3D_XZ, GMASK3D_YZ
!-------------------------------------------------------------------------------
!
!*       1.    P
!              ---------------------------------------
!
CPROGRAM='DIAPRO'
!
! Initialisation des parametres de Namelists
!
CALL INIDEF
CTITALL='DEFAULT'; LTITDEF=.TRUE.
CALL RESOLV_TIT('CTITALL',YTEM)
!
! Ouverture du fichier de conservation des directives
! Son nom:  dir.date
!
!CALL DATE(YDAT)
YTEM='          '
CALL DATE_AND_TIME(YDAT,YTEM)
YDIRNAM(1:4)='dir.'
!YTEM(1:2)=YDAT(7:8); YTEM(4:5)=YDAT(4:5); YTEM(7:8)=YDAT(1:2)
!YTEM(1:2)=YDAT(7:8); YTEM(4:5)=YDAT(5:6); YTEM(7:8)=YDAT(3:4)
!YTEM(3:3)=':'; YTEM(6:6)=':'
YDIRNAM(5:6)=YDAT(7:8)
YDIRNAM(7:8)=YDAT(5:6)
YDIRNAM(9:10)=YDAT(3:4)
YDIRNAM(11:11)=':'
YDIRNAM(12:13)=YTEM(1:2)
YDIRNAM(14:14)=':'
YDIRNAM(15:16)=YTEM(3:4)
!YDIRNAM(5:12)=YTEM
CALL FMATTR(YDIRNAM,YDIRNAM,IDIR,IDIRESP)
OPEN(UNIT=IDIR,FILE=YDIRNAM,FORM='FORMATTED')
NDIR=IDIR
!
! Lecture et interpretation des directives
!
DO JJ = 1,100000
  CAR240(1:LEN(CAR240))=' '
  CGROUP(1:LEN(CGROUP))=' '
  CGROUPS(:)(1:LEN(CGROUPS(1)))=' '
    IF(JJ == 1)THEN
      print *,' ENTREZ VOS DIRECTIVES '
    ELSE
      print *,' DIRECTIVE ? '
    ENDIF
  DO JI = 1,ILINES   ! directive sur ILINES lignes
    CAR80(1:LEN(CAR80))=' '
    CAR100(1:LEN(CAR100))=' '
    READ(5,'(A100)',END=10)CAR100
    CAR100=ADJUSTL(CAR100)
    IF( LEN_TRIM(CAR100)>80 .AND. CAR100(1:1) /= '!' ) THEN
       print *,'-- Directive:'
       print *,TRIM(CAR100)
       print *,' depassant 80 car. : ABORT'
       CAR80='QUIT'
       GO TO 99
    ENDIF 
    READ(CAR100,'(A80)')CAR80
    CAR80=ADJUSTL(CAR80)
    !WRITE(IDIR,'(A80)')CAR80
    CALL WRITEDIR(IDIR,CAR80)
    GO TO 20
    10 CONTINUE
    CLOSE(5)
    CALL GETENV("VARTTY",CAR20)
    CAR20=ADJUSTL(CAR20)
    OPEN(5,FILE=CAR20)
    print *,' diaprog INTERACTIF : ENTREZ VOS DIRECTIVES '
    20 CONTINUE
    CAR80=ADJUSTL(CAR80)
    print *,CAR80(1:LEN_TRIM(CAR80))
! Test  FF et commentaires
    IF(CAR80(1:1) == '!')EXIT 
    IF(CAR80(1:4) == 'QUIT' .OR. CAR80(1:4) == 'quit')GO TO 99
    INDE = INDEX(CAR80,'&')
    IF(INDE == 0)THEN
      ! directive sur une ligne
      ILENC=LEN_TRIM(CAR80)
      ILENC240=LEN_TRIM(CAR240)
      IF (ILENC240+ILENC .LE. LEN(CAR240)) THEN
        CAR240(ILENC240+1:ILENC240+ILENC)=CAR80(1:ILENC)
      ELSE
        print *,'Erreur! '//CAR240(1:20)//'...  depasse ',LEN(CAR240),' caracteres'
        CAR240=' '
      ENDIF
      EXIT
    ELSE
      IF (JI==ILINES) THEN
        print *,'-- Pas plus de ',ILINES,' lignes pour une directive : ABORT'
        CAR80='QUIT'
        GO TO 99
      ENDIF 
      DO JIA=INDE-1,1,-1
	IF(CAR80(JIA:JIA) /= ' ')THEN
	  ! suite des directives ligne suivante
          ILENC240=LEN_TRIM(CAR240)
          ILENC=JIA
          IF (ILENC240+ILENC .LE. LEN(CAR240)) THEN
            CAR240(ILENC240+1:ILENC240+ILENC)=CAR80(1:ILENC)
          ELSE
            print *,'Erreur! '//CAR240(1:20)//'...  depasse ',LEN(CAR240),' caracteres'
            CAR240=' '
          ENDIF
          EXIT
        END IF
      ENDDO
    END IF
  ENDDO
#ifdef RHODES
CALL FLUSH(IDIR,ISTAF)
#else
CALL FLUSH(IDIR)
#endif

IF(LEN_TRIM(CAR240) == 0)CYCLE
!
! Conversion des mots cles des instructions en MAJUSCULES
!
CDIRCUR(1:LEN(CDIRCUR))=' '
CALL CONVLO2UP(CAR240(1:LEN_TRIM(CAR240)),YCAR240)
IF(LPBREAD)THEN
  LPBREAD=.FALSE.
  CYCLE
ENDIF
CDIRCUR(1:LEN_TRIM(YCAR240))=YCAR240(1:LEN_TRIM(YCAR240))
!
CAR240(1:LEN(CAR240))=' '
CAR240=ADJUSTL(YCAR240)
print* ,CAR240(1:LEN_TRIM(CAR240))
!
! Juillet 2001 *  ou / par un processus DEB*****************
!
! Desallocation des tableaux si RM*EXPRx (avec x=1 a 9)
INDEXPR=INDEX(CAR240,'RM*EXPR')
IF(INDEXPR /= 0)THEN
  IZERO=0
  CALL LOAD_EXPR(IZERO,CAR240(1:LEN_TRIM(CAR240)))
  CYCLE
ENDIF
! Desallocation des tableaux si RM/EXPRx (avec x=1 a 9)
INDEXPR=INDEX(CAR240,'RM/EXPR')
IF(INDEXPR /= 0)THEN
  IZERO=0
  CALL LOAD_EXPR(IZERO,CAR240(1:LEN_TRIM(CAR240)))
  CYCLE
ENDIF
! Chargement du processus a * ou /
INDEXPR=INDEX(CAR240,'*EXPR')
IF(INDEXPR /= 0)THEN
  IF(CAR240(INDEXPR+6:INDEXPR+6) == '=')THEN
  IZERO=0
  CALL LOAD_EXPR(IZERO,CAR240(1:LEN_TRIM(CAR240)))
  CYCLE
  ENDIF
ENDIF
INDEXPR=INDEX(CAR240,'/EXPR')
IF(INDEXPR /= 0)THEN
  IF(CAR240(INDEXPR+6:INDEXPR+6) == '=')THEN
  IZERO=0
  CALL LOAD_EXPR(IZERO,CAR240(1:LEN_TRIM(CAR240)))
  CYCLE
  ENDIF
ENDIF
!
! Juillet 2001 *  ou / par un processus FIN*****************
!
! Nov 2001 Deplacement des impressions + haut
!
! Traitement des impressions
!
INDPRI=INDEX(CAR240,'PRINT ')
INDPRIL=INDEX(CAR240,'LPRINT ')
IF(INDPRI /= 0 .AND. INDPRIL == 0)THEN
  CALL PRINTS(CAR240(1:LEN_TRIM(CAR240)))
  CYCLE
ENDIF
!
! Lecture eventuelle du groupe SV3 a utiliser comme coord. vert.
!
ICGROUPSV3=INDEX(CAR240,'CGROUPSV3')
IF(ICGROUPSV3 /= 0)THEN
  CGROUPSV3(1:LEN(CGROUPSV3))=' '
  ILENT=LEN_TRIM(CAR240)
  IQUOT=INDEX(CAR240,"'")
  IF(IQUOT == 0)THEN
    IQUOT=INDEX(CAR240,'"')
  ENDIF
  CGROUPSV3=CAR240(IQUOT+1:ILENT-1)
  CGROUPSV3=ADJUSTL(CGROUPSV3)
  print *,' CGROUPSV3 FOURNI ',CGROUPSV3
  CYCLE
ENDIF
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!
! Lecture eventuelle du groupe TRAJ_GROUP a utiliser pour les trajectoires
!
ICTRAJ_GROUP=INDEX(CAR240,'CTRAJ_GROUP')
IF(ICTRAJ_GROUP /= 0)THEN
  CTRAJ_GROUP(1:LEN(CTRAJ_GROUP))=' '
  ILENT=LEN_TRIM(CAR240)
  IQUOT=INDEX(CAR240,"'")
  IF(IQUOT == 0)THEN
    IQUOT=INDEX(CAR240,'"')
  ENDIF
  CTRAJ_GROUP=CAR240(IQUOT+1:ILENT-1)
  CTRAJ_GROUP=ADJUSTL(CTRAJ_GROUP)
  print *,' CTRAJ_GROUP FOURNI ',CTRAJ_GROUP
  CYCLE
ENDIF
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!!!!!!!!
!
! Conversion d'indices de grille I,J en coord. conf et geographiques
!
ICONVIJ2XY=INDEX(CAR240,'CONVIJ2XY')
IF(ICONVIJ2XY /= 0)THEN
  CALL CONVIJ2XY(CAR240)
  IF(LPBREAD)LPBREAD=.FALSE.
  CYCLE
ENDIF
!
ICONVALLIJ2LL=INDEX(CAR240,'CONVALLIJ2LL')
IF(ICONVALLIJ2LL /= 0)THEN
  CALL CONVALLIJ2LL(CAR240)
  IF(LPBREAD)LPBREAD=.FALSE.
  CYCLE
ENDIF
!
! Conversion de coord. conf et geographiques en indices de grille
!
ICONVXY2IJ=INDEX(CAR240,'CONVXY2IJ')
IF(ICONVXY2IJ /= 0)THEN
  CALL CONVXY2IJ(CAR240)
  IF(LPBREAD)LPBREAD=.FALSE.
  CYCLE
ENDIF
!
! Memorisation des textes et symboles associes a un couple lat,lon
!
ICNOMCAR=INDEX(CAR240,'CNOMCAR')
IF(ICNOMCAR /= 0)THEN
  CALL CARESOLV(CAR240)
! CALL CARESOLV(CAR240(1:LEN_TRIM(CAR240)))
  IF(LPBREAD)LPBREAD=.FALSE.
  NBGUIL=0
  CYCLE
ENDIF
ICSYMCAR=INDEX(CAR240,'CSYMCAR')
IF(ICSYMCAR /= 0)THEN
  CALL CARESOLV(CAR240)
! CALL CARESOLV(CAR240(1:LEN_TRIM(CAR240)))
  IF(LPBREAD)LPBREAD=.FALSE.
  if(nverbia >0)then
    print *,' ***DIAPROG ICSYMCAR > 0 AV CYCLE'
  endif
  NBGUIL=0
  CYCLE
ENDIF
!
! Traitement des eventuels segments de dte a superposer sur une CH en PCart.
!
INDTIT=INDEX(CAR240,'XSEGM')
IF(INDTIT /= 0)THEN
  CALL LOAD_SEGMENTS(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  CYCLE
ENDIF
!
! Traitement des eventuels segments de dte a superposer sur une CH
!
INDTIT=INDEX(CAR240,'ISEGM')
IF(INDTIT /= 0)THEN
  CALL LOAD_SEGMENTS(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  CYCLE
ENDIF
!
! Traitement des eventuels formats des labels des axes
!
INDTIT=INDEX(CAR240,'CFMTAXEX')
IF(INDTIT /= 0)THEN
  CALL LOAD_FMTAXES(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  CYCLE
ENDIF
INDTIT=INDEX(CAR240,'CFMTAXEY')
IF(INDTIT /= 0)THEN
  CALL LOAD_FMTAXES(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  CYCLE
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!! pour les retrotrajectoires                                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Traitement eventuel du format des labels des retrotrajectoires
!
INDTIT=INDEX(CAR240,'CFMTRTRAJ')
IF(INDTIT /= 0)THEN
  CALL LOAD_FMTAXES(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  CYCLE
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Traitement des eventuels titres
!
ITITDEF=0
INDTIT=INDEX(CAR240,'LTITDEF')
IF(INDTIT /= 0)THEN
  CALL LOAD_TIT(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  ITITDEF=1
ENDIF
INDTIT=INDEX(CAR240,'CTIT')
IF(INDTIT /= 0)THEN
  CALL LOAD_TIT(CAR240(1:LEN_TRIM(CAR240)),INDTIT)
  IF(INDTIT == 999)THEN
    INDTIT=0
    CYCLE
  ENDIF
  INDTIT=0
ENDIF
IF(ITITDEF == 1)THEN
  ITITDEF=0
  CYCLE
ENDIF
!print *,CAR240
!print *,LEN_TRIM(CAR240)
!READ(*,*)
!
!  Ajout artificiel d'un niveau de modele avec _MSKTOP_ pour des facilites
!  de programmation (si absent)
!
ITOP1=INDEX(CAR240,'_MSKTOP_')
IF(ITOP1 /= 0)THEN
  ITOP2=INDEX(CAR240,'_K_')
  IF(ITOP2 == 0)THEN
    CAR240=ADJUSTL(ADJUSTR(CAR240)//'_K_2')
    print *,' **diaprog . directive generee volontairement :',CAR240(1:LEN_TRIM(CAR240))
  ENDIF
ENDIF
!
!  Traitement SAV et NOSAV (SAVE et NOSAVE)
!
IF(CAR240(1:LEN_TRIM(CAR240)) == 'NOSAV'  .OR. &
   CAR240(1:LEN_TRIM(CAR240)) == 'NOSAVE')THEN
   CALL GDAWK(1)
   CYCLE
ELSE IF(CAR240(1:LEN_TRIM(CAR240)) == 'SAV'  .OR. &
   CAR240(1:LEN_TRIM(CAR240)) == 'SAVE')THEN
   CALL GACWK(1)
   CYCLE
ENDIF
!
! Extraction des noms de fichiers; eventuelle ouverture; mise a jour du numero
! de fichier courant dans la variable NUMFILECUR. Elimination du nom des
! fichiers et des sequences _FILE_ _FILEx_ FILExx_ des instructions d'entree
!
if (nverbia >0)then
  print *,' ****DIAPROG AV EXTRACT_AND_OPEN_FILES '
  print *,CAR240(1:LEN_TRIM(CAR240))
endif

CALL EXTRACT_AND_OPEN_FILES(CAR240(1:LEN_TRIM(CAR240)),YCAR240)
IF(LPBREAD)THEN
  LPBREAD=.FALSE.
  CYCLE
ENDIF

if (nverbia >0)then
  print *,' AP EXTRACT_AND_OPEN_FILES '
  print *,YCAR240(1:LEN_TRIM(YCAR240))
endif
CAR240(1:LEN_TRIM(CAR240))=' '
CAR240=ADJUSTL(YCAR240)
!
! Memorisation de l'instruction d'entree en vue de reutiliser les specifica-
! -pour le groupe suivant avec l'option IDEM
!
!IF(JJ == 1)THEN
!  CALL CARMEMORY(CAR240,1)
!ENDIF
!
!
! Resolution des temps, processus, niveaux, altitudes ...
!
if (nverbia >0)then
  print *,' AV CARESOLV'
endif
CALL CARESOLV(CAR240(1:LEN_TRIM(CAR240)))
if (nverbia >0)then
  print *,' AP CARESOLV'
endif
IF(LPBREAD)THEN
  LPBREAD=.FALSE.
  CYCLE
ENDIF
!
!
DO JLOOP=1,NSUPERDIA

  NLOOPSUPER=JLOOP
  LXYZ=LXYZT(JLOOP)
! Mars 2000
  LUMVMPV=LUMVMPVT(JLOOP)
! Mars 2000
! Memorisation pour - et + 
  IF(JLOOP == 1)THEN
    LTITDEFM=LTITDEF
    CTITB3MEM=CTITB3
    CTITB3MEM=ADJUSTL(CTITB3MEM)
    if(nverbia >0)print *,' **diaprog LTITDEFM, CTITB3MEM ',LTITDEFM,CTITB3MEM
  ENDIF
!!!!!!Oct 2000 Prise en compte de superposition d'un pH issu du 2D Hor. sur
! une CV
  IF(NHISTORY(JLOOP) == 1)THEN
    LCH=.FALSE.
    LCV=.TRUE.
  ELSEIF(NHISTORY(JLOOP) == 3)THEN
    LCH=.TRUE.
    LCV=.TRUE.
  ENDIF
!!!!!!Oct 2000

  IF(NBPM > 1)THEN
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
!   IF(JLOOP >= 2)THEN
    IF(JLOOP >= 1)THEN
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
      IF(NUMPM(JLOOP) == 1)THEN
	LPLUS=.TRUE.
	LMINUS=.FALSE.
      ELSE IF(NUMPM(JLOOP) == 2)THEN
	LMINUS=.TRUE.
	LPLUS=.FALSE.
      ELSE
	LMINUS=.FALSE.
	LPLUS=.FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
        IF(JLOOP < NSUPERDIA)THEN
          IF(NUMPM(JLOOP+1) == 1)THEN
            LPLUS=.TRUE.
          ELSE IF(NUMPM(JLOOP+1) == 2)THEN
            LMINUS=.TRUE.
          ENDIF
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
      ENDIF
    ENDIF
  ENDIF
! print *,' PG PAL LMINUS LPLUS JLOOP ',LMINUS,LPLUS,JLOOP

  CGROUP=ADJUSTL(CGROUPS(JLOOP))
  IF(CGROUP(1:LEN_TRIM(CGROUP)) == ' ')THEN
    EXIT
  ELSE
    NUMFILECUR=NFILESCUR(JLOOP)
    DO J=1,NBFILES
      IF(NUMFILES(J) == NUMFILECUR)THEN
      JM=J
      ENDIF
    ENDDO
!
! Lecture du type d'informations demandees
!
    CALL READ_TYPE(CFILEDIAS(JM),CLUOUTDIAS(JM),CGROUP(1:LEN_TRIM(CGROUP)))
! Sorties sur surfaces isobares ou isentropes ou emagrammes
! Chargement des infos utiles
    IF(LPBREAD)THEN
      LPBREAD=.FALSE.
      EXIT
    ENDIF
	if(nverbia >0)then
	  print *,' **diaprog AP READ_TYPE LTK,LPR,LEV,LSV3 ',LTK,LPR,LEV,LSV3
	endif
!
! Chargement de la temperature pour sorties en surfaces isentropes et
! emagrammes
!
    ! test pour eviter les messages de READ_TH_PR dans ce cas
    IF(CGROUP(1:LEN_TRIM(CGROUP)) == 'ZS' .OR. CGROUP(1:LEN_TRIM(CGROUP)) == &
      'ZSBIS')THEN
    ELSE

    IF((LTK .OR. LEV .OR. LSV3) .OR. LPR .OR. LPRESY .OR. ((LRS .OR. LRS1) .AND. CTYPE=='CART'))THEN
      NMT=2
      IF(CGROUP(LEN_TRIM(CGROUP):LEN_TRIM(CGROUP)) == 'M')NMT=1
!     NMT=1
!     IF(CGROUP(LEN_TRIM(CGROUP):LEN_TRIM(CGROUP)) == 'T')NMT=2
    ENDIF

! LTK = .TRUE. ou LRS = .TRUE. ou LRS1 = .TRUE. ou LEV=.TRUE.

    IF((LTK .OR. LEV .OR. LSV3) .OR. ((LRS .OR. LRS1) .AND. CTYPE=='CART'))THEN
      CALL READ_TH_PR(CFILEDIAS(JM),CLUOUTDIAS(JM),NMT,1)
      IF(LPBREAD)THEN
	LPBREAD=.FALSE.
!       EXIT
	IF(NMT == 1)THEN
	  NMT=2
	ELSE
	  NMT=1
	ENDIF
        CALL READ_TH_PR(CFILEDIAS(JM),CLUOUTDIAS(JM),NMT,1)
        IF(LPBREAD)THEN
	  LPBREAD=.FALSE.
          EXIT
        ENDIF
      ENDIF
      IF(LSV3 .AND. MAXVAL(XZHAT)/MAXVAL(XTH) > 1.E2)THEN
        IF(.NOT.LXYZ00 .OR. CGROUPSV3 == 'Z00')THEN
	if(nverbia >0)then
	  print *,' **diaprog MAXVAL(XZHAT)/MAXVAL(XTH) ',MAXVAL(XZHAT)/MAXVAL(XTH)
	endif
	WHERE(XTH /= XSPVAL)
	  XTH=XTH*1.E3
	ENDWHERE
	if(nverbia >0)then
	  print *,' **diaprog MAXVAL(XTH) ap *1.E3',MAXVAL(XTH)
	  print *,' **diaprog MINVAL(XTH) ap *1.E3',MINVAL(XTH)
	endif
        ENDIF
      ENDIF
    ENDIF

  ENDIF
!
! LPR = .TRUE. ou LRS = .TRUE.  ou LRS1 = .TRUE. Calcul ou lecture de la pression
!
    IF(LPR .OR. LPRESY .OR. ((LRS .OR. LRS1) .AND. CTYPE=='CART'))THEN
      CALL READ_TH_PR(CFILEDIAS(JM),CLUOUTDIAS(JM),NMT,2)
      IF(LPBREAD)THEN
	LPBREAD=.FALSE.
!       EXIT
	IF(NMT == 1)THEN
	  NMT=2
	ELSE
	  NMT=1
	ENDIF
        CALL READ_TH_PR(CFILEDIAS(JM),CLUOUTDIAS(JM),NMT,2)
        IF(LPBREAD)THEN
	  LPBREAD=.FALSE.
	  IF(LPRESY)THEN
	    print *,' Pression absente (PABSM et PABST) -> LPRESY remis a F '
	    LPRESY=.FALSE.
	  ENDIF
          EXIT
        ENDIF
      ENDIF
    ENDIF
!
! Chargement des composantes du vent dans le cas de combinaisons de celles-ci
!
    IF(LUMVM .OR. LMUMVM .OR. LULM .OR. LVTM .OR. LULMWM .OR. &
       LUTVT .OR. LMUTVT .OR. LULT .OR. LVTT .OR. LULTWT .OR. &
       LDIRWM .OR. LDIRWT .OR. &
       LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
      CALL READ_UVW(CFILEDIAS(JM),CLUOUTDIAS(JM),CGROUP(1:LEN_TRIM(CGROUP)))

! Janvier 2001  Vecteurs vent horizontal et direction en CV
      IF(LPBREAD)THEN
!     IF(LPBREAD .OR. (LUMVM .AND. LCV .AND..NOT. LCH) .OR. (LUTVT .AND. LCV &
!     .AND..NOT. LCH))THEN
	LPBREAD=.FALSE.
!       IF(LUMVM .OR. LUTVT)THEN
!         print *,' VECTEURS VENT HORIZONTAL NON PREVUS EN COUPE VERTICALE'
!       ENDIF
	  IF(ALLOCATED(XU))THEN
	    DEALLOCATE(XU)
	  ENDIF
	  IF(ALLOCATED(XV))THEN
	    DEALLOCATE(XV)
	  ENDIF
	  IF(ALLOCATED(XVAR))THEN
	    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	  ENDIF
	IF(JLOOP > 1)CALL FRAME
	EXIT
      ENDIF
    ENDIF
!
! Lecture des informations autres que pour RS + combinaisons composantes vent
!
    IF((.NOT.LRS .AND. .NOT.LRS1 .AND. .NOT.LUMVM .AND. .NOT.LMUMVM .AND. &
	.NOT.LULM .AND. .NOT.LVTM .AND. .NOT.LULMWM .AND. .NOT.LUTVT .AND.&
	.NOT.LMUTVT .AND. .NOT.LULT .AND. .NOT.LVTT .AND. .NOT.LULTWT .AND.  &
	.NOT.LDIRWM .AND. .NOT.LDIRWT .AND. &
	.NOT.LSUMVM .AND. .NOT.LSUTVT .AND. .NOT.LMLSUMVM .AND. .NOT.LMLSUTVT)& 
	.OR. ((LRS .OR. LRS1) .AND. CTYPE /= 'CART'))THEN

      IF(LXYZ .OR. LMSKTOP)THEN
	IF(XXL == 0. .AND. XXH == 0. .AND. XYL == 0. .AND. XYH == 0. &
	  .AND. XZL == 0. .AND. XZH == 0.)THEN
	  print *,' Definissez une fenetre (en metres) dans XXL= XXH= XYL= XYH= XZL= XZH='
	  print *,' Et rentrez a nouveau votre directive '
	  IF(ALLOCATED(XVAR))THEN
	    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	  ENDIF
	  LXYZ=.FALSE. ; LMSKTOP=.FALSE.
	  EXIT
        ELSE
	  GMASK3D=LMASK3D
	  GMASK3D_XY=LMASK3D_XY
	  GMASK3D_XZ=LMASK3D_XZ
	  GMASK3D_YZ=LMASK3D_YZ
	  LMASK3D=.FALSE.; LMASK3D_XY=.FALSE.; LMASK3D_XZ=.FALSE.
	  LMASK3D_YZ=.FALSE.
          CALL TRAMASK3D
	  LMASK3D=GMASK3D
	  LMASK3D_XY=GMASK3D_XY
	  LMASK3D_XZ=GMASK3D_XZ
	  LMASK3D_YZ=GMASK3D_YZ
        ENDIF
	IF(LPBREAD)THEN
	  LPBREAD=.FALSE.
	  EXIT
	ENDIF
      ENDIF
      CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),CGROUP)
      IF(LPBREAD)THEN
        LPBREAD=.FALSE.
        EXIT
      ENDIF
      IF(LGROUP)THEN
      CALL READ_DIACHRO(CFILEDIAS(JM),CLUOUTDIAS(JM),CGROUP)
      ENDIF
!     print *,'SIZE(XVAR,1,2,3,4,5,6) ',SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3), &
!     SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)
!     print *,' XVAR(1,1,1,1,1,1) ',XVAR(1,1,1,1,1,1)
!     print *,' XVAR(1,1,1,2,1,1) ',XVAR(1,1,1,2,1,1)
!     print *,' XVAR(5,5,5,1,1,1) ',XVAR(5,5,5,1,1,1)
!     print *,' XVAR(5,5,5,2,1,1) ',XVAR(5,5,5,2,1,1)
      IF(LPBREAD)THEN
	IF(LFT .OR. LFT1)THEN
	ALLOCATE(ZBID1(1,1),ZBID2(1))
	CALL VARFCT(ZBID1,ZBID2,1)
	IF(JLOOP >1)CALL FRAME
	DEALLOCATE(ZBID1,ZBID2)
	ENDIF
	LPBREAD=.FALSE.
	IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	ENDIF
	EXIT 
      ENDIF

!     print *,' LFIC1 NBSIMULT ',LFIC1,NBSIMULT

      IF(.NOT.LFIC1)THEN

!       print *,' AV REALLOC_AND_LOAD '
        CALL REALLOC_AND_LOAD(CGROUP)
        IF(LPBREAD)THEN
          LPBREAD=.FALSE.
	  IF(ALLOCATED(XVAR))THEN
	    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	  ENDIF
          EXIT
        ENDIF
!     print *,' AP REALLOC_AND_LOAD '

      ELSE

        NBRECOUV=1
        NRECOUV(1)=1
        NRECOUV(2)=SIZE(XTRAJT,1)

      ENDIF
!     print *,' diaprog XSPVAL ',XSPVAL
        IF(LXYZ)THEN
!         IF(ALLOCATED(LMASK3))THEN
!           WHERE(.NOT.LMASK3)XVAR(:,:,:,1,1,1)=XSPVAL
!         ELSE
!           CALL TRAMASK3D
            WHERE(.NOT.LMASK3)XVAR(:,:,:,:,1,1)=XSPVAL
!           WHERE(.NOT.LMASK3)XVAR(:,:,:,1,1,1)=XSPVAL
!         ENDIF
	ENDIF

    ENDIF

! Pour distinguer 1 profil 1D enregistre comme tel (LPV=T et LCV=F) et 1 profil
! extrait d'une matrice 3D (LPV=T et LCV=t)
!
    IF(LPV .OR. LPVT .OR. LPVKT .OR. LPVKT1 .OR. LPXT .AND. (SIZE(XVAR,1)-1 > 0))LCV=.TRUE.
    IF(LPYT)LCV=.TRUE.
    IF(LPXT .OR. LPYT .AND. LCH)LCV=.FALSE.
    IF(NVERBIA > 0)THEN
      print *,' main LPXT LPYT LCV LCH ',LPXT,LPYT,LCV,LCH
    ENDIF

    if(nverbia >0)print *,' ****diaprog AV  KZTNP'
    CALL KZTNP(JLOOP)
    if(nverbia >0)print *,' ****diaprog AP  KZTNP LPBREAD ',LPBREAD
	if(nverbia >0)then
	  print *,' **diaprog AP KZTNP LTK,LPR,LEV,LSV3 ',LTK,LPR,LEV,LSV3
	endif
    IF(LPBREAD)THEN
      LPBREAD=.FALSE.
      IF(ALLOCATED(XVAR))THEN
        CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      ENDIF
      IF(ALLOCATED(XU))THEN
        DEALLOCATE(XU)
      ENDIF
      IF(ALLOCATED(XV))THEN
        DEALLOCATE(XV)
      ENDIF
      EXIT
    ENDIF

    IF((LRS .OR. LRS1) .AND. CTYPE == 'CART')THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    ENDIF

    IF(ALLOCATED(XVAR) .AND. NOPE(JLOOP) /= 0)THEN
      IF(NOPE(JLOOP) == 1)THEN
        !IF(LSV3 .OR. LXYZ)THEN
        WHERE(XVAR(:,:,:,:,:,:) /= XSPVAL)
          XVAR(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)+XCONSTANTE(JLOOP)
        ENDWHERE
        !ELSE
        !  XVAR(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)+XCONSTANTE(JLOOP)
        !ENDIF
! Janvier 2001
	IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. LSUMVM .OR. &
	   LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
	  IF(ALLOCATED(XU))THEN
	    XU(:,:,:,:,:,:)=XU(:,:,:,:,:,:)+XCONSTANTE(JLOOP)
	  ENDIF
	ENDIF
!       print *,' XCONSTANTE(JLOOP) ',XCONSTANTE(JLOOP)
      ELSE IF(NOPE(JLOOP) == 2)THEN
        !IF(LSV3 .OR. LXYZ)THEN
        WHERE(XVAR(:,:,:,:,:,:) /= XSPVAL)
          XVAR(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)*XCONSTANTE(JLOOP)
        ENDWHERE
        !ELSE
        ! XVAR(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)*XCONSTANTE(JLOOP)
        !ENDIF
! Janvier 2001
	IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. LSUMVM .OR. &
	   LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
	  IF(ALLOCATED(XU))THEN
	    XU(:,:,:,:,:,:)=XU(:,:,:,:,:,:)*XCONSTANTE(JLOOP)
	  ENDIF
	ENDIF
      ELSE IF(NOPE(JLOOP) == 3)THEN
	WHERE(XVAR(:,:,:,:,:,:) /= XSPVAL .AND. XVAR >0.)
	  XVAR(:,:,:,:,:,:)=LOG(XVAR(:,:,:,:,:,:))
        ELSEWHERE
	  XVAR(:,:,:,:,:,:)=XSPVAL
	ENDWHERE
      ENDIF
    ENDIF
! Juillet 2001

    IF(ALLOCATED(XVAR) .AND. NMULTDIV(JLOOP) /= 0)THEN
      print *,' ++diaprog JLOOP,NMULTDIV(JLOOP),CMULTDIV(JLOOP) ',JLOOP,NMULTDIV(JLOOP),CMULTDIV(JLOOP)
      IMULTDIV=NMULTDIV(JLOOP)
      YMULTDIV=' '
      YMULTDIV=CMULTDIV(JLOOP)
      YMULTDIV=ADJUSTL(YMULTDIV)
      CALL LOAD_EXPR(IMULTDIV,YMULTDIV(1:LEN_TRIM(YMULTDIV)))
    ENDIF

! Juillet 2001

! Difference entre 2 champs (ou somme) . Presence de la chaine _MINUS_ (_PLUS_)
!
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
!   IF((LMINUS .OR. LPLUS) .AND. JLOOP == 1)THEN
    IF((LMINUS .OR. LPLUS) .AND. (NUMPM(JLOOP) == 0 .OR. NUMPM(JLOOP) == 3))THEN
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
! On memorise le 1er champ
      IF(NBPROCDIA(JLOOP) == 1)THEN
	NGRIDIAM=NGRIDIA(NPROCDIA(NBPROCDIA(JLOOP),JLOOP))
      ELSE
	print *,' ** diaprog Nb de processus > 1 pour une somme ou difference'
      ENDIF
      CALL ALLOC2_FORDIACHRO(1)
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
        IF(ALLOCATED(XU))THEN
          DEALLOCATE(XU)
        ENDIF
        IF(ALLOCATED(XV))THEN
          DEALLOCATE(XV)
        ENDIF
      CYCLE
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
!   IF((LMINUS .OR. LPLUS) .AND. JLOOP >= 2)THEN
    IF((LMINUS .OR. LPLUS) .AND. (NUMPM(JLOOP) == 1 .OR. NUMPM(JLOOP) == 2))THEN
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
      CALL DIFF_OPER(JLOOP)
      CALL ALLOC2_FORDIACHRO(3)
      IF(LPBREAD)THEN
	LPBREAD=.FALSE.
	CTITB3=CTITB3MEM
	LTITDEF=LTITDEFM
        IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
        ENDIF
        IF(ALLOCATED(XU))THEN
          DEALLOCATE(XU)
        ENDIF
        IF(ALLOCATED(XV))THEN
          DEALLOCATE(XV)
        ENDIF
    IF(ALLOCATED(XUMEM))THEN
      CALL ALLOC2_FORDIACHRO(3)
    ENDIF
	if(nverbia > 0)then
	print *,' ** diaprog LPBREAD=T DEALLOCATE '
	endif
        CYCLE
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
!     IF(JLOOP < NBPM)THEN
      IF(JLOOP < NBPM .AND.(NUMPM(JLOOP+1) == 1 .OR. NUMPM(JLOOP+1) == 2))THEN
!!!!!!!!!!!!!!!!!!!!!!!!020398!!!!!!!!!!!!!!!!!!!
	CALL ALLOC2_FORDIACHRO(1)
	CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
        IF(ALLOCATED(XU))THEN
          DEALLOCATE(XU)
        ENDIF
        IF(ALLOCATED(XV))THEN
          DEALLOCATE(XV)
        ENDIF
	CYCLE
      ENDIF
    ENDIF

    if(nverbia >0)print *,' ****diaprog AV OPER_PROCESS'
	if(nverbia >0)then
	  print *,' **diaprog AV OPER_PROC LTK,LPR,LEV,LSV3 ',LTK,LPR,LEV,LSV3
	endif
    CALL OPER_PROCESS(JLOOP,CTYPE)
    if(nverbia >0)print *,' ****diaprog AP OPER_PROCESS'
    IF(LPBREAD)THEN
      LPBREAD=.FALSE.
    ENDIF
    LDIRWIND=.FALSE.
! Oct 2000
    LCHXY=.FALSE.

!   15022000
!! Nov 2001
    IF(NBPMT > 0 .AND. JLOOP == NSUPERDIA)THEN
!   IF(LMINUS .OR. LPLUS .AND. JLOOP == NSUPERDIA)THEN
!! Nov 2001
!   IF(LMINUS .OR. LPLUS)THEN
      CTITB3=CTITB3MEM
      LTITDEF=LTITDEFM
      if(nverbia > 0)then
	print *,' ** diaprog FIN boucle JLOOP == NSUPERDIA CTITB3 LTITDEF NBPMT ',CTITB3,LTITDEF,NBPMT
      endif
    ENDIF

    IF(ALLOCATED(XVAR))THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    ENDIF
    IF(ALLOCATED(XU))THEN
      DEALLOCATE(XU)
    ENDIF
    IF(ALLOCATED(XV))THEN
      DEALLOCATE(XV)
    ENDIF
    IF(ALLOCATED(XUMEM))THEN
      CALL ALLOC2_FORDIACHRO(3)
    ENDIF
    IF(JLOOP == NSUPERDIA)THEN
      XIDEBCOU=-999.;XJDEBCOU=-999.
!6666666666666666666666666666666666666666666666
!     NIINF=0; NJINF=0; NISUP=0; NJSUP=0
!6666666666666666666666666666666666666666666666
    ENDIF

  ENDIF

ENDDO

CDIRPREC=' '
CDIRPREC=CDIRCUR
CDIRPREC=ADJUSTL(CDIRPREC)

ENDDO
99 CONTINUE
CAR240(1:80)=CAR80
CALL CONVLO2UP(CAR240(1:LEN_TRIM(CAR240)),YCAR240)
CAR240=ADJUSTL(YCAR240)
if (nverbia >0)then
  print *,' ****DIAPROG 2 AV EXTRACT_AND_OPEN_FILES '
  print *,CAR240(1:LEN_TRIM(CAR240))
endif
CALL EXTRACT_AND_OPEN_FILES(CAR240(1:LEN_TRIM(CAR240)),YCAR240)
if (nverbia >0)then
  print *,' ****DIAPROG 2 AP EXTRACT_AND_OPEN_FILES '
  print *,YCAR240(1:LEN_TRIM(YCAR240))
endif
CAR240(1:LEN(CAR240))=' '
CAR240=ADJUSTL(YCAR240)
!READ(*,*)
!CALL CARESOLV(CAR240(1:LEN_TRIM(CAR240)))
CLOSE(IDIR)
if (nverbia >0)then
  print *,' ****DIAPROG 3 AP EXTRACT_AND_OPEN_FILES '
endif

STOP 
END PROGRAM DIAPROG
