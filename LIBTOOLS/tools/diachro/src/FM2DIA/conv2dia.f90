!     ######spl
      PROGRAM  FM2DIACHRO
!     ###################
!
!!****  *FM2DIACHRO* -  Conversion des fichiers synchrones LFIFM en
!!                      fichiers de type diachronique (LFIFM egalement)
!! 
!!
!!    PURPOSE
!!    -------
! 
!       Convertit 1 (ou plusieurs fichiers synchrones correspondant a
!       des sorties successives d'un meme run) en 1 fichier diachronique
!
!!**  METHOD
!!    ------
!!      
!       La routine LFILAF (du logiciel LFI) modifiee (--> JDLFILAF) pour
!       l'ouverture d'un fichier FICJD ecrit dans celui-ci le numero,
!       le nom et la longueur totale des enregistrements.
!       Puis un appel a la routine LFILEC permet de lire dans le 2eme mot
!       de chaque enregistrement la longueur du champ commentaire (qui n'est
!       pas necessairement constante) et donc de deduire par soustraction
!       la longueur du champ physique enregistre 
!       de sorte que l'on possede toutes les informations necessaires a la
!       lecture avec FMREAD des enregistrements d'un fichier LFIFM dont on ne 
!       connait pas a priori le contenu. (du moins pour les infos reelles)
!       Dans un premier temps, on ecrit dans le fichier diachonique avec
!       la routine WRITE_LFIFM1_FORDIACHRO_CV l'entete des fichiers d'entree
!       en particulier les parametres de grille, l'etat de reference ...
!       Puis en bouclant sur le nombre de fichiers a traiter et le nombre
!       d'enregistrements de chacun, on lit chaque champ et on regroupe
!       progressivement dans un enregistrement du fichier diachronique unique
!       pour un meme parametre les differentes echeances trouvees.
!       ACTUELLEMENT (Avril 97) SONT PRIS EN COMPTE LES CHAMPS DE LONGUEUR
!       IIU*IJU*IKU  , IIU*IJU  et  1
!
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
!!      Original    30/01/96 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF          
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT          
USE MODD_DIM1, ONLY: NIMAX,NJMAX,NKMAX          
USE MODD_GRID  ! XLON0,XLAT0, XBETA,XRPK
USE MODD_GRID1 ! XLONOR,XLATOR
USE MODD_TIME1 ! TDTCUR
!
USE MODD_DIACHRO
USE MODD_OUT_DIA
USE MODD_REA_LFI    
USE MODD_DIMGRID_FORDIACHRO
!USE MODI_READ_DESFM
USE MODI_READ_DIMGRIDREF_FM2DIA
USE MODI_WRITE_DIMGRIDREF
USE MODI_WRITE_OTHERSFIELDS
USE MODI_MENU_DIACHRO
USE MODI_INI_CST

IMPLICIT NONE
!
!*       0.1   Local variables declarations
!
INTEGER           :: ILUDES    ! Logical unit number for the DES file
INTEGER           :: INUMER

INTEGER,DIMENSION(100) :: IFICJD

INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: IWORK 

INTEGER           :: INUM, ISIZ, INBM

INTEGER           :: IRESP, IVAR, IEL_OR_SEL
INTEGER           :: INEWSIZE, ITYPCOD

INTEGER           :: JJ, J, JA, J1, J2, I2
INTEGER           :: INB, IID, JI, JIP1, ICODEL, IL, IDA
INTEGER           :: I4, IOK, IKEEP
INTEGER           :: IX,IY,IZ,ixyz     ! resolution degradee
INTEGER           :: IIMAXIN,IJMAXIN,IKMAXIN,IIMAXOUT,IJMAXOUT,IKMAXOUT  !  "
REAL,DIMENSION(:),   ALLOCATABLE  :: ZXHAT,ZYHAT !,ZZHAT  !        "
REAL,DIMENSION(:,:), ALLOCATABLE  :: ZZS                  !        "
INTEGER,DIMENSION(:), ALLOCATABLE  :: IIMAX, IJMAX, IKMAX
REAL,DIMENSION(:), ALLOCATABLE  :: ZTIMECUR,ZLON0,ZRPK,ZLONOR,ZLATOR,ZLAT0,ZBETA
LOGICAL,DIMENSION(:), ALLOCATABLE :: OCARTESIAN
LOGICAL           :: GDTOUT, GOK

CHARACTER(LEN=16) :: YRECFM, YRECFM2
CHARACTER(LEN=3)  :: YREPON
CHARACTER(LEN=16) :: YREF
CHARACTER(LEN=16) :: YCOMMENT
CHARACTER(LEN=80) :: YCAR80  
CHARACTER(LEN=16),DIMENSION(100)      :: YFICJD, YFICJDOUT
CHARACTER(LEN=16),DIMENSION(:), ALLOCATABLE,SAVE :: YRECT, YRECID, YKEEP
CHARACTER(LEN=16),DIMENSION(5)                   :: YPRI
!-------------------------------------------------------------------------------
!
!*       1.    Definition du type de traitement et init du fichier de constantes
!              -----------------------------------------------------------------
!
CPROGRAM='FM2DIA'
!
CCONF='POSTP'
CALL INI_CST
OPEN(80,FILE='dirconv',FORM='FORMATTED')
!
!
!*	 2.    Lecture du nombre de fichiers a regrouper et de leur nom 
!              --------------------------------------------------------
!              Doivent etre dissocies en *.des et *.lfi et
!              rentres en ordre chronologique (1 / 1 ligne)
!
PRINT *,' ENTER NUMBER OF INPUT FM FILES'
READ(5,*)NNBF
YCAR80(1:LEN(YCAR80))=' '
WRITE(YCAR80,*)NNBF
YCAR80=ADJUSTL(YCAR80)
WRITE(80,'(A80)')YCAR80

DO J=1,NNBF
  PRINT *,' ENTER FM FILE NAME'
  READ(5,'(A28)')CNAMFILED(J)   
  YCAR80(1:LEN(YCAR80))=' '
  YCAR80=CNAMFILED(J)
  YCAR80=ADJUSTL(YCAR80)
  WRITE(80,'(A80)')YCAR80
ENDDO
!
!
!*	 3.    Lecture du nom du fichier diachronique a creer
!              ----------------------------------------------
!

PRINT *,' ENTER DIACHRONIC FILE NAME'
READ(5,'(A28)')CFILEDIA
YCAR80(1:LEN(YCAR80))=' '
YCAR80=CFILEDIA     
YCAR80=ADJUSTL(YCAR80)
WRITE(80,'(A80)')YCAR80
!
!*       4.    Ouverture du fichier correspondant au listing
!              ---------------------------------------------
!
CLUOUTD='LISTING_DIA'
CALL FMATTR(CLUOUTD,CLUOUTD,NLUOUTD,NRESP)
OPEN(UNIT=NLUOUTD,FILE=CLUOUTD,FORM='FORMATTED')
!
!*       5.    Boucle sur les fichiers a lire 
!              ------------------------------
!
DO J=1,NNBF

  CLFIFMD(J)=ADJUSTL(ADJUSTR(CNAMFILED(J))//'.lfi')
  CDESFMD(J)=ADJUSTL(ADJUSTR(CNAMFILED(J))//'.des')

!
!*       5.1   Ouverture des fichiers LFIFM et DESFM
!
  CSTATU='OLD'
  NVERB=5
! Modif demandee par Nicole Asencio. 28/9/98
  NFTYPE=2
! NFTYPE=0
  CALL FMOPEN(CNAMFILED(J),CSTATU,CLUOUTD,NNPRAR,NFTYPE,NVERB,NNINAR,NRESP)
  IF(NRESP.NE.0)THEN
    WRITE(0,*)'BUG OPENING LFIFM FILE ',CLFIFMD(J),'  RETURN CODE= ',NRESP
  END IF
!
!*       5.2   Fermeture du fichier DESFM  (ACTUELLEMENT NON INTEGRE DANS LE
!                                           FICHIER DIACHRONIQUE)
! en 5.6 avec LFIFM par FMCLOS
!
!*	 5.3   Lecture du numero, nom et longueur des enregistrements
!              Memorisation dans les tableaux NNUMT,CRECFM2T,NSIZT
!
!
  GDTOUT=.TRUE.
  CALL FMLOOK(CLFIFMD(J),CLUOUTD,INUMER,NRESP)
  CALL JDLFILAF(NRESP,INUMER,GDTOUT)
!
  YFICJD(J)='FICJD'
  YFICJDOUT(J)='FICJDOUT'
  CALL FMATTR(YFICJD(J),YFICJDOUT(J),IFICJD(J),NRESP)
  OPEN(UNIT=IFICJD(J),FILE=YFICJD(J),FORM='FORMATTED',STATUS='OLD')
!
  NNB=0
  DO JJ=1,10000
    READ(IFICJD(J),*,END=99)INUM,YRECFM2,ISIZ
    NNB=NNB+1
  ENDDO
99 CONTINUE

  IF(J == 1)THEN
    INBM=NNB
  ENDIF

  WRITE(NLUOUTD,*)' ******** FICHIER N: ',J,CNAMFILED(J)(1:LEN_TRIM(CNAMFILED(J))), &
  ' NB ENR. ',NNB
  WRITE(NLUOUTD,*)' ******** '
!
  REWIND(IFICJD(J))
!
  IF(J == 1)THEN
    ALLOCATE(NNUMT(NNB+100,100),NSIZT(NNB+100,100),NLENC(NNB+100,100))
    ALLOCATE(CRECFM2T(NNB+100,100))
  ENDIF
!
  DO JJ=1,NNB
    READ(IFICJD(J),*)NNUMT(JJ,J),CRECFM2T(JJ,J),NSIZT(JJ,J)
    ALLOCATE(IWORK(NSIZT(JJ,J)))
    CALL LFILEC(NRESP,INUMER,CRECFM2T(JJ,J),IWORK,NSIZT(JJ,J))
    NLENC(JJ,J)=IWORK(2)     ! longueur de la zone commentaire
! Determination de la longueur de la zone de donnees
! 2 = 1er mot : numero de grille et 2eme mot : longueur de la zone commentaire
    NSIZT(JJ,J)=NSIZT(JJ,J)-2-NLENC(JJ,J)
    CALL GET_COMPHEADER(IWORK(3+NLENC(JJ,J)),NSIZT(JJ,J),INEWSIZE,ITYPCOD)
    IF (INEWSIZE >= 0) THEN ! compressed field found
      WRITE (NLUOUTD,*) TRIM(CRECFM2T(JJ,J)),' is compressed (old/new SIZE):',NSIZT(JJ,J),INEWSIZE
      NSIZT(JJ,J)=INEWSIZE
    END IF
    DEALLOCATE(IWORK)
  ENDDO
!
  CLOSE (IFICJD(J))
  CALL FMFREE(YFICJD(J),YFICJDOUT(J),NRESP)

! Verification de l'egalite du nombre d'enregistrements dans les differents
! fichiers

  IF(J > 1)THEN
    IF(INBM /= NNB)THEN
      WRITE(NLUOUTD,*)' ******************************************'
      WRITE(NLUOUTD,*)' Nb enregistrents different (/ 1er fichier)'
      WRITE(NLUOUTD,*)' ******************************************'
      WRITE(NLUOUTD,*)' ( - = absence par rapport au 1er fichier, + = ajout)'
      WRITE(NLUOUTD,*)' ( + ne sont pas integres dans le fichier diachronique)'
    ENDIF
  ENDIF

! Verification de l'identite des enregistrements dans les differents fichiers

  IF(J > 1)THEN
    IF(INBM /= NNB)THEN
      IF (INBM > NNB)THEN
        DO JJ=1,INBM
    	  GOK=.FALSE.
    	  DO JA=1,NNB
    	    IF(CRECFM2T(JJ,1) == CRECFM2T(JA,J))THEN
    	      GOK=.TRUE.
    	      EXIT
    	    ELSE
    	      CYCLE
    	    ENDIF
    	  ENDDO
  	  IF(.NOT.GOK)THEN
  	    NNUMT(JJ,1)=0
  	    WRITE(NLUOUTD,*)' -  ',CRECFM2T(JJ,1)
  	  ENDIF
        ENDDO

      ELSE

        DO JJ=1,NNB
      	  GOK=.FALSE.
      	  DO JA=1,INBM
      	    IF(CRECFM2T(JJ,J) == CRECFM2T(JA,1))THEN
      	      GOK=.TRUE.
      	      EXIT
      	    ELSE
      	      CYCLE
      	    ENDIF
      	  ENDDO
	  IF(.NOT.GOK)THEN
	    WRITE(NLUOUTD,*)' +  ',CRECFM2T(JJ,J)
	  ENDIF
  	ENDDO
      ENDIF
    ENDIF
  ENDIF
  !
!
!*       5.4   Lecture et ecriture des parametres "intouchables"
!        5.40  lecture
!
  CALL READ_DIMGRIDREF_FM2DIA(J,CNAMFILED(J),CLUOUTD)
!
!        5.41  Writing or checking  DIM., GRID., REF. VARIABLES
!
  IF(J == 1)THEN  ! premier fichier
! resolution degradee
    IX=1 ; IY=1 ; IZ=1
    IF (NIMAX>1) THEN
      print *,'- DO YOU WANT COARSER RESOLUTION along X ? (y/n)'
      READ(5,*)YREPON
      YCAR80(1:LEN(YCAR80))=' '
      YCAR80=YREPON
      YCAR80=ADJUSTL(YCAR80)
      WRITE(80,'(A80)')YCAR80
      IF(YREPON == 'y' .OR. YREPON == 'yes' .OR. YREPON == 'o' .OR. &
         YREPON == 'oui')THEN
        print *,' Enter the ratio IX (1 point on IX points kept) '
        READ(5,*) IX
        YCAR80(1:LEN(YCAR80))=' '
        WRITE(YCAR80,*)IX
        YCAR80=ADJUSTL(YCAR80)
        WRITE(80,'(A80)')YCAR80
      ENDIF
    ENDIF
    IF (IX>1) THEN
      IIMAXOUT=((NIMAX+2*JPHEXT-1)/IX +1) -2*JPHEXT
      IF (IIMAXOUT<=0) THEN
        print *,'TOO COARSER RESOLUTION along X for ',NIMAX,'points'
        STOP
      ENDIF
      print*,'old X physical domain: ',NIMAX,'pts  - new one: ',IIMAXOUT
    ENDIF
    !
    IF (NJMAX>1) THEN
      print *,'- DO YOU WANT COARSER RESOLUTION along Y ? (y/n)'
      READ(5,*)YREPON
      YCAR80(1:LEN(YCAR80))=' '
      YCAR80=YREPON
      YCAR80=ADJUSTL(YCAR80)
      WRITE(80,'(A80)')YCAR80
      IF(YREPON == 'y' .OR. YREPON == 'yes' .OR. YREPON == 'o' .OR. &
         YREPON == 'oui')THEN
        print *,' Enter the ratio IY (1 point on IY points kept) '
        READ(5,*) IY
        YCAR80(1:LEN(YCAR80))=' '
        WRITE(YCAR80,*)IY
        YCAR80=ADJUSTL(YCAR80)
        WRITE(80,'(A80)')YCAR80
      ENDIF
    ENDIF
    IF (IY>1) THEN
      IJMAXOUT=((NJMAX+2*JPHEXT-1)/IY +1) -2*JPHEXT
      IF (IJMAXOUT<=0) THEN
        print *,'TOO COARSER RESOLUTION along Y for ',NJMAX,'points'
        STOP
      ENDIF
      print*,'old Y physical domain: ',NJMAX,'pts  - new one: ',IJMAXOUT
    ENDIF
    !
    !print *,'- DO YOU WANT COARSER RESOLUTION along Z ? (y/n)'
    !READ(5,*)YREPON
    !IF(YREPON == 'y' .OR. YREPON == 'yes' .OR. YREPON == 'o' .OR. &
    !   YREPON == 'oui')THEN
    !  print *,' Enter the ratio IZ (1 point on IZ points kept) '
    !  READ(5,*) IZ
    !ENDIF
    !IF (IZ>1) THEN
    !  IKMAXOUT=((NKMAX+2*JPVEXT-1)/IZ +1) -2*JPVEXT
    !  IF (IKMAXOUT<=0) THEN
    !    print *,'TOO COARSER RESOLUTION along Y for ',NKMAX,'points'
    !    STOP
    !  ENDIF
    !  print*,'old Z physical domain: ',NKMAX,'pts  - new one: ',IKMAXOUT
    !ENDIF
    !
    IF (IX>1) THEN
      ALLOCATE(ZXHAT(SIZE(XXHAT)))
      ZXHAT(:)=XXHAT(:)
      DEALLOCATE(XXHAT)
      ALLOCATE(XXHAT(IIMAXOUT+2*JPHEXT))
      XXHAT(:)=ZXHAT(1:NIMAX+2*JPHEXT:IX)
    ENDIF
    IF (IY>1) THEN
      ALLOCATE(ZYHAT(SIZE(XYHAT)))
      ZYHAT(:)=XYHAT(:)
      DEALLOCATE(XYHAT)
      ALLOCATE(XYHAT(IJMAXOUT+2*JPHEXT))
      XYHAT(:)=ZYHAT(1:NJMAX+2*JPHEXT:IY)
    ENDIF
    ixyz=0
    IF (IX>1) ixyz=1
    IF (IY>1) ixyz=ixyz+10
    IF (ixyz>0) THEN
      ALLOCATE(ZZS(SIZE(XZS,1),SIZE(XZS,2)))
      ZZS(:,:)=XZS(:,:)
      DEALLOCATE(XZS)
    ENDIF
    SELECT CASE(ixyz)
      CASE (1)   !X
        ALLOCATE(XZS(IIMAXOUT+2*JPHEXT,SIZE(ZZS,2)))
        DO J2=1,SIZE(ZZS,2)
          XZS(:,J2)=ZZS(1:NIMAX+2*JPHEXT:IX,J2)
        END DO
        IIMAXIN=NIMAX
        NIMAX  =IIMAXOUT
      CASE (10)  !Y
        ALLOCATE(XZS(SIZE(ZZS,1),IJMAXOUT+2*JPHEXT))
        DO J1=1,SIZE(ZZS,1)
          XZS(J1,:)=ZZS(J1,1:NJMAX+2*JPHEXT:IY)
        END DO
        IJMAXIN=NJMAX
        NJMAX  =IJMAXOUT
      CASE (11)  !X et Y
        ALLOCATE(XZS(IIMAXOUT+2*JPHEXT,IJMAXOUT+2*JPHEXT))
        I2=0
        DO J2=1,SIZE(ZZS,2),IY
          I2=I2+1
          XZS(:,I2)=ZZS(1:NIMAX+2*JPHEXT:IX,J2)
        END DO
        IIMAXIN=NIMAX
        NIMAX  =IIMAXOUT
        IJMAXIN=NJMAX
        NJMAX  =IJMAXOUT
    END SELECT
    !IF (IZ>1) THEN
    !  ALLOCATE(ZZHAT(SIZE(XZHAT)))
    !  ZZHAT(:)=XZHAT(:)
    !  DEALLOCATE(XZHAT)
    !  ALLOCATE(XZHAT(IKMAXOUT+2*JPVEXT))
    !  XZHAT(:)=ZZHAT(1:NKMAX+2*JPVEXT:IZ)
    !  IKMAXIN=NKMAX
    !  NKMAX  =IKMAXOUT
    !ENDIF
    !
    CALL WRITE_DIMGRIDREF
    !  
    IF (IX>1) THEN
      NIMAX=IIMAXIN
      DEALLOCATE(XXHAT)
      ALLOCATE(XXHAT(SIZE(ZXHAT)))
      XXHAT(:)=ZXHAT(:)
      DEALLOCATE(ZXHAT)
    ENDIF
    IF (IY>1) THEN
      NJMAX=IJMAXIN
      DEALLOCATE(XYHAT)
      ALLOCATE(XYHAT(SIZE(ZYHAT)))
      XYHAT(:)=ZYHAT(:)
      DEALLOCATE(ZYHAT)
    ENDIF
    !IF (IZ>1) THEN
    !  NKMAX=IKMAXIN
    !  DEALLOCATE(XZHAT)
    !  ALLOCATE(XZHAT(SIZE(ZZHAT)))
    !  XZHAT(:)=ZZHAT(:)
    !  DEALLOCATE(ZZHAT)
    !ENDIF
    IF (ixyz>0) THEN
      DEALLOCATE(XZS)
      ALLOCATE(XZS(SIZE(ZZS,1),SIZE(ZZS,2)))
      XZS(:,:)=ZZS(:,:)
      DEALLOCATE(ZZS)
    ENDIF
    !
    ALLOCATE(IIMAX(NNBF),IJMAX(NNBF),IKMAX(NNBF),ZTIMECUR(NNBF))
    ALLOCATE(ZLON0(NNBF),ZLAT0(NNBF),ZLONOR(NNBF),ZLATOR(NNBF), &
                                        ZRPK(NNBF),ZBETA(NNBF)  )
    ALLOCATE(OCARTESIAN(NNBF))
  ENDIF
!
  IIMAX(J)=NIMAX ; IJMAX(J)=NJMAX ; IKMAX(J)=NKMAX
  ZTIMECUR(J)=TDTCUR%TIME
  ZLON0(J)=XLON0   ; ZLAT0(J)=XLAT0
  ZLONOR(J)=XLONOR ; ZLATOR(J)=XLATOR
  ZRPK(J)=XRPK     ; ZBETA(J)=XBETA
  OCARTESIAN(J)=LCARTESIAN
!
  IF(J > 1)THEN   ! fichiers suivants
  !
    IF(IIMAX(J) /= IIMAX(1))THEN
      PRINT *,' J IIMAX(J) IIMAX(1) ',J,IIMAX(J),IIMAX(1)
    ENDIF
    IF(IJMAX(J) /= IJMAX(1))THEN
      PRINT *,' J IJMAX(J) IJMAX(1) ',J,IJMAX(J),IJMAX(1)
    ENDIF
    IF(IKMAX(J) /= IKMAX(1))THEN
      PRINT *,' J IKMAX(J) IKMAX(1) ',J,IKMAX(J),IKMAX(1)
    ENDIF
    IF(ZTIMECUR(J) /= ZTIMECUR(1))THEN
      PRINT *,' J ZTIMECUR(J) ZTIMECUR(1) ',J,ZTIMECUR(J),ZTIMECUR(1)
    ENDIF
    IF(ZLON0(J) /= ZLON0(1))THEN
      PRINT *,' J ZLON0(J) ZLON0(1) ',J,ZLON0(J),ZLON0(1)
    ENDIF
    IF(ZRPK(J) /= ZRPK(1))THEN
      PRINT *,' J ZRPK(J) ZRPK(1) ',J,ZRPK(J),ZRPK(1)
    ENDIF
    IF(ZLONOR(J) /= ZLONOR(1))THEN
      PRINT *,' J ZLONOR(J) ZLONOR(1) ',J,ZLONOR(J),ZLONOR(1)
    ENDIF
    IF(ZLATOR(J) /= ZLATOR(1))THEN
      PRINT *,' J ZLATOR(J) ZLATOR(1) ',J,ZLATOR(J),ZLATOR(1)
    ENDIF
    IF(ZLAT0(J) /= ZLAT0(1))THEN
      PRINT *,' J ZLAT0(J) ZLAT0(1) ',J,ZLAT0(J),ZLAT0(1)
    ENDIF
    IF(ZBETA(J) /= ZBETA(1))THEN
      PRINT *,' J ZBETA(J) ZBETA(1) ',J,ZBETA(J),ZBETA(1)
    ENDIF
    IF((OCARTESIAN(J) .AND..NOT. OCARTESIAN(1)) .OR. &
       (.NOT. OCARTESIAN(J) .AND. OCARTESIAN(1)))THEN
      PRINT *,' J OCARTESIAN(J) OCARTESIAN(1) ',J,OCARTESIAN(J),OCARTESIAN(1)
    ENDIF
    !
  ENDIF
!
  IF(J == NNBF)THEN  ! dernier fichier
    DEALLOCATE(IIMAX,IJMAX,IKMAX,ZTIMECUR)
    DEALLOCATE(ZLON0,ZRPK,ZLONOR,ZLATOR,ZLAT0,ZBETA)
    DEALLOCATE(OCARTESIAN)
  END IF
!
!        5.42  Eventuelle eliminination de certains parametres ds le fic. diach.
!
  IF(J == 1)THEN

    print *,'- DELETION OF PARAMETERS AT TIME t-dt ? (enter 1) '
    print *,'- DELETION OF PARAMETERS AT TIME t    ? (enter 2) '
    print *,'- NO DELETION                         ? (enter 0) '
    READ(5,*)ICODEL
    YCAR80(1:LEN(YCAR80))=' '
    WRITE(YCAR80,*)ICODEL
    YCAR80=ADJUSTL(YCAR80)
    WRITE(80,'(A80)')YCAR80
    
    IF(ICODEL == 0)THEN
        IEL_OR_SEL=0                                ! conv2dia.elim
    ELSE
      print *,'- Do you want to ELIM or to SELECT parameters ? (E/S)'
      READ(5,*)YREPON
      YCAR80(1:LEN(YCAR80))=' '
      YCAR80=YREPON
      YCAR80=ADJUSTL(YCAR80)
      WRITE(80,'(A80)')YCAR80
      IF(YREPON == 'E' .OR. YREPON == 'e')THEN
        IEL_OR_SEL=0                                ! conv2dia.elim
        print*,'as conv2dia.elim'
      ELSE IF(YREPON == 'S' .OR. YREPON == 's')THEN
        IEL_OR_SEL=1                                ! conv2dia.select
        print*,'as conv2dia.select'
      ELSE 
        STOP 'Bad answer'
      ENDIF
    ENDIF
    !
    ALLOCATE(YRECT(SIZE(CRECFM2T,1)))
    YRECT(:)(1:LEN(YRECT))=' '
    INB=0
    DO JI=1,NNB
    IF(NNUMT(JI,J) /= 0)THEN
      INB=INB+1
      YRECT(INB)=CRECFM2T(JI,J)
      YRECT(INB)=ADJUSTL(YRECT(INB))
!     print *,' INB, YRECT ',INB,YRECT(INB)
    ENDIF    
    ENDDO
    ALLOCATE(YRECID(INB))
    YRECID(:)(1:LEN(YRECID))=' '
    IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
      ALLOCATE(YKEEP(INB))
      YKEEP(:)(1:LEN(YKEEP))=' '
      IKEEP=1
      YKEEP(IKEEP)='ZS'
      YKEEP(IKEEP)=ADJUSTL(YKEEP(IKEEP))
      IKEEP=IKEEP+1
      YKEEP(IKEEP)='ZSMT'
      YKEEP(IKEEP)=ADJUSTL(YKEEP(IKEEP))
    ENDIF
    !
    IID=0
    DO JI = 1,INB-1
      YREF(1:LEN(YREF))=' '
      IL=LEN_TRIM(YRECT(JI))-1
      YREF(1:IL)=YRECT(JI)(1:IL)
!     YREF=ADJUSTL(YREF)
      IF(YRECT(JI)(IL+1:IL+1) == 'M')THEN
      DO JIP1=JI+1,INB
	IL=LEN_TRIM(YRECT(JIP1))-1
	IF(YRECT(JIP1)(1:IL) == YREF .AND. YRECT(JIP1)(IL+1:IL+1)=='T' )THEN
          IF ( IEL_OR_SEL==0 .OR.            &  ! conv2dia.elim
              (IEL_OR_SEL==1                 &  ! conv2dia.select
               .AND.(YREF(1:IL)=='PABS'      &  ! et PABS
                     .OR.YREF(1:IL)=='POVO'  &  ! ou POVO
                     .OR.YREF(1:IL)=='TH'   ))) THEN  ! ou TH
            IID=IID+1
            YRECID(IID)=' '
            YRECID(IID)=YREF
            YRECID(IID)=ADJUSTL(YRECID(IID))
            !
            IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
              IKEEP=IKEEP+1
              YKEEP(IKEEP)=YREF
              YKEEP(IKEEP)=ADJUSTL(YKEEP(IKEEP))
            ENDIF
            EXIT
          ENDIF
        ENDIF
      ENDDO
      ENDIF
    ENDDO
    IF(ICODEL == 0)THEN
    ELSE IF(ICODEL == 1)THEN
      DO JI=1,IID
        YRECID(JI)=ADJUSTL(ADJUSTR(YRECID(JI))//'M')
      ENDDO
      IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
        DO JI=3,IKEEP
          YKEEP(JI)=ADJUSTL(ADJUSTR(YKEEP(JI))//'T')
        ENDDO
      ENDIF
    ELSE IF(ICODEL == 2)THEN
      DO JI=1,IID
        YRECID(JI)=ADJUSTL(ADJUSTR(YRECID(JI))//'T')
      ENDDO
      IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
        DO JI=3,IKEEP
          YKEEP(JI)=ADJUSTL(ADJUSTR(YKEEP(JI))//'M')
        ENDDO
      ENDIF
    ENDIF
    !
    I4=0
    YPRI=' '
    IF(ICODEL /= 0)THEN

    print *,' PARAMETRES RESTANTS'
    DO JI = 1,NNB
      DO JIP1 = 1,IID
        IF(CRECFM2T(JI,J) == YRECID(JIP1))THEN
  	  NNUMT(JI,J)=0
  	  EXIT
        ENDIF
      ENDDO
      IF(NNUMT(JI,J) /= 0)THEN
	I4=I4+1
	YPRI(I4)=CRECFM2T(JI,J)
	IF(I4 == 5 .OR. JI == NNB)THEN
          print 10,YPRI
	  I4=0
          YPRI=' '
	ENDIF
      ENDIF     
    ENDDO

    IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
      print *,' '
      print *,' Some parameters are automatically recorded (for vert. interpolations):'
      print *,' --> ',(YKEEP(JI)(1:LEN_TRIM(YKEEP(JI))+1),JI=1,IKEEP)
      print *,' '
    ENDIF

    YREPON(1:LEN(YREPON))=' '
    IF (IEL_OR_SEL==0) THEN                    ! conv2dia.elim
      print *,'- Do you want to SUPPRESS others parameters ? (y/n) '
    ELSE IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
      print *,'- Do you want to KEEP others parameters ? (y/n) '
    ENDIF
    READ(5,*)YREPON
    YCAR80(1:LEN(YCAR80))=' '
    YCAR80=YREPON
    YCAR80=ADJUSTL(YCAR80)
    WRITE(80,'(A80)')YCAR80
    IF(YREPON == 'y' .OR. YREPON == 'yes' .OR. YREPON == 'o' .OR. &
       YREPON == 'oui' .OR. YREPON =='Y' .OR. YREPON =='YES'      )THEN
      print *,'- Enter their names in UPPERCASE  (1/1 line) '
      print *,'End by END '
      IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
        IF(ICODEL == 1)THEN
          print *,'  NOTA: if you want to plot RS ,don''t forget : RVT,UT,VT'
        ELSE IF(ICODEL == 2)THEN
          print *,'  NOTA: if you want to plot RS ,don''t forget : RVM,UM,VM'
        ENDIF
        print *,' '
      ENDIF
      DO JI=1,1000
        YREF=' '
        READ(5,*)YREF
        YCAR80(1:LEN(YCAR80))=' '
        IF (IEL_OR_SEL==0) THEN                    ! conv2dia.elim
          IID=IID+1
          YRECID(IID)=' '
          YRECID(IID)=ADJUSTL(YREF)
          YCAR80=YRECID(IID)  
        ELSE IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
          IKEEP=IKEEP+1
          YKEEP(IKEEP)=' '
          YKEEP(IKEEP)=ADJUSTL(YREF)
          YCAR80=YKEEP(IKEEP)  
        ENDIF
        YCAR80=ADJUSTL(YCAR80)
        WRITE(80,'(A80)')YCAR80
        IF(YREF == 'END')THEN
          CLOSE(80)
          EXIT
        ENDIF
      ENDDO
    ENDIF
    !
    10 FORMAT(1X,5A15)
    I4=0
!   YPRI(:)=' '
    DO JI = 1,NNB
      IF (IEL_OR_SEL==0) THEN                    ! conv2dia.elim
        DO JIP1 = 1,IID
          IF(CRECFM2T(JI,J) == YRECID(JIP1))THEN
            NNUMT(JI,J)=0
            EXIT
          ENDIF
        ENDDO
      ELSE IF (IEL_OR_SEL==1) THEN               ! conv2dia.select
        IF(NNUMT(JI,J) /=0) THEN
          IOK=0
          DO JIP1 = 1,IKEEP
            IF(CRECFM2T(JI,J) == YKEEP(JIP1))THEN
              IOK=1
              EXIT
            ENDIF
          ENDDO
          IF(IOK==0)THEN
            NNUMT(JI,J)=0
          ENDIF
        ENDIF
      ENDIF
      IF(NNUMT(JI,J) /= 0)THEN
        IF(I4 == 5)THEN
          print 10,YPRI
          I4=0
!         YPRI(1:5)='    '
        ENDIF
        I4=I4+1
        YPRI(I4)=CRECFM2T(JI,J)
      ENDIF     
      IF(JI == NNB)THEN
        print 10,YPRI(1:I4)
      ENDIF     
    ENDDO

    ENDIF     ! (ICODEL/=0)
!              
!              
!        5.43  Elimination des dates
!              
    DO JI=1,NNB
      IDA=INDEX(CRECFM2T(JI,J),'%TDA')
      IF(IDA /= 0)THEN
        NNUMT(JI,J)=0
      ENDIF
      IDA=INDEX(CRECFM2T(JI,J),'%TIM')
      IF(IDA /= 0)THEN
        NNUMT(JI,J)=0
      ENDIF
!        5.44  Elimination des champs dont le nom depasse 13 caracteres
!        (13 = 16 (=max.LEN(RECFM)=JPNCPN) -3 (=LEN('.TYpe','.DIm','.TItre',
!                              '.UNite','.COmment','.PRoc1','.TRajt','.DAtim'))
      IF (LEN_TRIM(CRECFM2T(JI,J))>13 .AND. NNUMT(JI,J)/=0) THEN
        NNUMT(JI,J)=0
        print*,'Variable ',CRECFM2T(JI,J), ' not written (name too long)'
        WRITE(NLUOUTD,*)'Variable ',CRECFM2T(JI,J), ' not written (name too long)'
      END IF
    ENDDO
  ENDIF    !(J==1)
!
!
!*       5.5   Lecture et ecriture des autres champs
!
  CALL WRITE_OTHERSFIELDS(J,CFILEDIA,CLUOUTDIA,IX,IY,IZ)
!
!*       5.6   Fermeture du Fichier d'entree traite et liberation des unites
!              logiques correspondantes (DES et LFI)
!
  CALL FMCLOS(CNAMFILED(J),'KEEP',CLUOUTD,NRESP)
!
ENDDO
!
!*       6.    Terminaison du fichier diachronique et impression du nom des
!              groupes enregistres
!              -------------------------------------------------------------
!
CALL MENU_DIACHRO(CFILEDIA,CLUOUTDIA,'END')
CALL MENU_DIACHRO(CFILEDIA,CLUOUTDIA,'READ')

CLOSE(NLUOUTD)
CALL FMFREE(CLUOUTD,CLUOUTD,NRESP)
!
!*       7.    Fermeture du fichier diachronique 
!              ---------------------------------
!
CALL FMCLOS(CFILEDIA,'KEEP',CLUOUTDIA,NRESP)
!
!------------------------------------------------------------------------------
!
!*      4.    EPILOGUE
!             --------

STOP

END PROGRAM FM2DIACHRO
