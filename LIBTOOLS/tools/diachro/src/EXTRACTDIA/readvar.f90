!     ######
      SUBROUTINE  READVAR(HLABELCHAMP,HFILENAME,HFLAGFILE,&
       KVERBIA,KRETCODE)
!     ################
!
!!****  *READVAR* - 
!! 
!!
!!    PURPOSE
!!    -------
!     Extraction d un champ du fichier diachronique et initialisation
!     des differents parametres utiles (grille, relief...)
! 
!
!!**  METHOD
!!    ------
!     utilisation des routines de diaprog : le tableau de stockage
!     XVAR est alloué par les routines de lecture.
!
!     au maximum 44 fichiers simultanement ouverts 
!       44 =limite FMOPEN= (JPNXFM-10)/2 avec JPNXFM=99
!
!     HFLAGFILE='OPE' lors de la premiere utilisation du fichier
!     HFLAGFILE='NOP' lors des utilisations suivantes
!     HFLAGFILE='CLO' fermeture du fichier traite ( decremente
!      le nombre de fichiers ouverts comptabilises par FMOPEN)
!
!     KVERBIA= 0 impressions reduites au minimum (entree et sortie de la
!      routine)
!     KVERBIA >0 impressions pour signaler chaque etape de READVAR
!
!     KRETCODE = 0 execution de READVAR correcte
!     KRETCODE = 1 erreur lors de l ouverture du fichier
!     KRETCODE = 2 champ inconnu dans le fichier
!     KRETCODE = 3 Nombre de fichiers ouverts simultanement > limite
!
!!
!!    EXTERNAL
!!    --------
!!          CREATLINK : à l'ouverture du fichier, HFLAGFILE='OPE',
!!                      création d'un lien dans le directory local
!!                      si le fichier existe sous $DIRLFI
!!          TO_COMPUTING_UNITS: passage unites vers unites plus pertinentes 
!!                              pour effectuer des calculs       
!!
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHORS
!!    -------
!!    I. Mallet et N. Asencio * CNRM*
!!
!!    Copyright 2003,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/03/2003
!!      N. Asencio  01/2005    call To_Computing_units
!!      G. TANGUY  03/2010     problème pour les champs sur point de flux 
!                              on remplace les 999 sur les mailles à côtés des bords du domaine 
!                              par la valeur la plus proche dans le domaine zoomé
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! modules MesoNH
USE MODD_PARAMETERS, ONLY: XUNDEF,JPHEXT
USE MODD_DIM1, ONLY: NIMAX,NJMAX,NKMAX
USE MODD_GRID1, ONLY: XZZ
! modules DIACHRO
!                    grille : XXDXHAT(:,1:7) et XXX(:,1:7), XXZS(:,:,1:7)
USE MODD_COORD
USE MODD_TYPE_AND_LH, ONLY: NIL,NIH,NJL,NJH,NKL,NKH,CTYPE,LICP,LJCP,LKCP
!                    XVAR(i,j,k,,,), XMASK,XTRAJ ,XDATIME(16,t) ,CUNITE(p)
USE MODD_ALLOC_FORDIACHRO
!                    nom de fichiers NLUOUT,CLFIFM, CDESFM
USE MODD_OUT
USE MODD_FILES_DIACHRO, ONLY: NBFILES,CFILEDIAS,CLUOUTDIAS,NRESPDIAS, &
                              NLUOUTDIAS, NNPRARDIAS, NFTYPEDIAS,     &
                              NNINARDIAS, NVERBDIAS
!
USE MODD_DIACHRO, ONLY:CFILEDIA       
!
USE MODI_FMREAD
USE MODI_READ_DIACHRO
USE MODI_VERIF_GROUP
USE MODI_ALLOC_FORDIACHRO
!
! modules TOOL
USE MODI_CREATLINK
! modules EXTRACTDIA
USE MODI_TO_COMPUTING_UNITS
USE MODD_READLH
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------
!
CHARACTER(LEN=*), INTENT(IN) :: HLABELCHAMP, HFILENAME ! nom du champ et du fichier
CHARACTER(LEN=3), INTENT(INOUT) :: HFLAGFILE           ! ouverture/ deja ouvert
INTEGER, INTENT(IN)          :: KVERBIA                ! prints de controle
!
INTEGER, INTENT(OUT)         :: KRETCODE   ! Code de retour de la routine 
!
!*       0.2   Local variables
!              ---------------
!
CHARACTER(LEN=13) :: YGP ! limite a 13 (ou 9 si plusieurs procs) 
                         !car read_diachro lit YRECFM(1:16)=YGP//'.PROCnn'
CHARACTER(LEN=32) :: YDESFM
INTEGER           :: JLOOP,JLOOPFIN,JI                              
INTEGER           :: IRESP,ILUDES
INTEGER           :: ILENG, ILENCH, IGRID, ILENDIM, IGROUP
INTEGER           :: idim3
INTEGER,DIMENSION(:),ALLOCATABLE :: ITABCHAR
CHARACTER(LEN=16) :: YRECFM
CHARACTER(LEN=20) :: YCOMMENT
CHARACTER(LEN=16),DIMENSION(:),ALLOCATABLE    :: YGROUP 
! pour traiter les champs budget deja zoomes
REAL , allocatable, dimension(:,:,:,:,:,:):: ZVARSAVE !
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
!      
print *,'---------'
print *,'Beginning of READVAR ',TRIM(HFILENAME),' ',HFLAGFILE,' ',TRIM(HLABELCHAMP)
!
! Code de retour de la routine : 0 = OK
!                                1 = erreur lors de l ouverture du fichier
!                                2 = champ inconnu
!                                3 = erreur sur le nombre de fichier
IF ( HFLAGFILE /= 'OPE' .AND. HFLAGFILE /= 'NOP' .AND. HFLAGFILE /= 'CLO' ) THEN
  KRETCODE=1
  print * ,'erreur d initialisation de HFLAGFILE =', HFLAGFILE
  print * ,'HFLAGFILE peut prendre les valeurs: OPE,NOP,CLO'
  print *,'---------'
  RETURN
ENDIF

KRETCODE=0
! code de retour d erreur des routines diaprog
LPBREAD=.FALSE.
!
IF(ALLOCATED(XVAR))THEN
! desallocation des tableaux alloues dans READ_DIACHRO (via ALLOC_FOR_DIACHRO)
  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
  if (KVERBIA >0)then
    print *,'*after ALLOC_FORDIACHRO(1,1,1,1,1,1,3)'
  endif
ENDIF
!-------------------------------------------------------------------------------
!
!*       2.    CLOSE THE FILE
!              --------------
!      
IF ( HFLAGFILE(1:3) == 'CLO' ) THEN
   CALL FMCLOS(HFILENAME,'KEEP',CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
   !! FMFREE ne relache pas l unite logique pour .lfi car .des deja relache
   DO JLOOP=1,NBFILES
     ! reperage de l indice de CFILEDIAS pour le fichier HFILENAME
     IF (CFILEDIAS(JLOOP) == HFILENAME ) THEN
      ! decalage du tableau CFILEDIAS pour supprimer cet element
        DO JLOOPFIN= JLOOP,NBFILES-1
          CFILEDIAS(JLOOPFIN)=CFILEDIAS(JLOOPFIN+1)
          CLUOUTDIAS(JLOOPFIN)=CLUOUTDIAS(JLOOPFIN+1)
          NLUOUTDIAS(JLOOPFIN)=NLUOUTDIAS(JLOOPFIN+1)                
          NNPRARDIAS(JLOOPFIN)=NNPRARDIAS(JLOOPFIN+1)
          NFTYPEDIAS(JLOOPFIN)=NFTYPEDIAS(JLOOPFIN+1)
          NVERBDIAS(JLOOPFIN)=NVERBDIAS(JLOOPFIN+1)
        ENDDO
        ! suppression du lien
        CALL CREATLINK('DIRLFI',CFILEDIAS(JLOOP),'CLEAN',KVERBIA)
        EXIT
     ENDIF
   ENDDO
   NBFILES=NBFILES-1
   print *,'End of READVAR: close of file ',TRIM(HFILENAME)
   print *,'---------'
   RETURN
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    OPEN THE FILE (first call)
!              --------------------------
!      
IF ( HFLAGFILE(1:3) == 'OPE' ) THEN
!
  if (KVERBIA >0)then
    print'(A23,I2,A17)','*before OPENning file, ',NBFILES,' currently opened'
  endif
!     utilisation de tableaux et de NBFILES pour calquer la methode
!     diaprog et permettre le traitement de plusieurs fichiers simultanement
  NBFILES=NBFILES+1
  !IF (NBFILES > 44 ) THEN
    ! 44 =limite FMOPEN= (JPNXFM-10)/2 avec JPNXFM=99
  !!limite >44 car fmfree de file.des
    !KRETCODE=3
    !print *,' ****READVAR: pour FMOPEN erreur nb de fichiers ouverts >44 ',&
    !          ' nbfiles= ',NBFILES
    !RETURN
  !ENDIF
  IF (NBFILES > size(CFILEDIAS) ) THEN
    KRETCODE=3
    print'(A58,I3,A10,I3)',' ****READVAR: pour diachro erreur nb de fichiers ouverts > ',&
                  size(CFILEDIAS), ' nbfiles= ',NBFILES
    print *,'---------'
    RETURN
  ENDIF
  CFILEDIAS(NBFILES)=HFILENAME
  CLUOUTDIAS(NBFILES)=CLUOUTDIAS(1)
  NNPRARDIAS(NBFILES)=0
  NFTYPEDIAS(NBFILES)= NFTYPEDIAS(1)
  NVERBDIAS(NBFILES)=KVERBIA
  ! listing OUT_DIA
  CALL FMLOOK(CLUOUTDIAS(NBFILES),CLUOUTDIAS(NBFILES),&
              NLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
  IF (NRESPDIAS(NBFILES)/=0) THEN
    ! ouverture du listing
    CALL FMATTR(CLUOUTDIAS(NBFILES),CLUOUTDIAS(NBFILES),&
                NLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
    OPEN(UNIT=NLUOUTDIAS(NBFILES),FILE=CLUOUTDIAS(NBFILES),&
       FORM='FORMATTED')
  END IF                  
  ! fichier diachronique
  CALL CREATLINK('DIRLFI',CFILEDIAS(NBFILES),'CREAT',KVERBIA)
  CALL FMOPEN(CFILEDIAS(NBFILES),'OLD',CLUOUTDIAS(NBFILES),&
              NNPRARDIAS(NBFILES),NFTYPEDIAS(NBFILES),NVERBDIAS(NBFILES),&
              NNINARDIAS(NBFILES),NRESPDIAS(NBFILES))
!      apres cet appel , variables initialisees:
!      NINARDIAS(NBFILES)= nb d articles dans le fichier
!      NRESPDIAS(NBFILES)= code de retour
!      une unite logique pour HFILENAME.des et HFILENAME.lfi
!    
  if (KVERBIA >0)then
    print'(A,A,A,5(I5,X))','*after OPENning files ',&
                    TRIM(CFILEDIAS(NBFILES)),&
                    TRIM(CLUOUTDIAS(NBFILES)),NNPRARDIAS(NBFILES), &
                    NFTYPEDIAS(NBFILES),NVERBDIAS(NBFILES),&
                    NNINARDIAS(NBFILES),NRESPDIAS(NBFILES)
  endif
  !
  IF (NRESPDIAS(NBFILES).NE.0)THEN
    KRETCODE=1
    print'(A52,A20,A6,I3)',' ****READVAR: erreur lors de l ouverture du fichier ',&
            CFILEDIAS (NBFILES), 'code= ',NRESPDIAS(NBFILES)
    print *,'---------'
    RETURN
  ENDIF
  !  
  ! partie DES du fichier: fermeture et unite logique relachee
  !YDESFM(1:LEN(YDESFM))=' '
  !YDESFM=ADJUSTL(ADJUSTR(CFILEDIAS(NBFILES))//'.des')
  !CALL FMLOOK(YDESFM,YDESFM,ILUDES,IRESP)
  !CLOSE(ILUDES)
  !CALL FMFREE(YDESFM,CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
!! ne pas relacher unite logique car compute_r00_pc doit fermer (avec FMCLOS)
!!le fic.  d entree qui a ete amende des var. Lag.
!
! READ JPHEXT
    CALL FMREAD(CFILEDIAS(NBFILES),'JPHEXT',CLUOUTDIAS(NBFILES),1,JPHEXT,IGRID,ILENCH,YCOMMENT,NRESPDIAS(NBFILES))

!*       3.1   Reading head of file
!              --------------------
!      
  CALL READ_FILEHEAD(1,HFILENAME,CLUOUTDIAS(NBFILES))
  if (KVERBIA >0)then
    print'(A41,3(I4,X))','*after READ_FILEHEAD, NIMAX,NJMAX,NKMAX= ',&
    NIMAX,NJMAX,NKMAX
  endif
  !
  ! lecture de MENU_BUDGET.DIM, MENU_BUDGET
  ! appel a INI_CST
  ! appel a READ_DIMGRIDREF: appel a SET_DIM pour lecture de IMAX, J,K-MAX 
  !                                         et calcul de I,J,K-INF,SUP 
  !                         lecture de CARTESIAN,THINSHELL,STORAGE_TYPE     
  !                         appel a SET_GRID
  ! appel a COMPCOORD_FORDIACHRO(0): pour les 7 grilles, 
  !          calcul de X,Y,Z-HAT(m) dans XXX,XXY,XXZ(:,1:7)        ! (MODD_COORD)
  !                 de topography altitude values(m):XXZS(:,:,1:7) ! (MODD_COORD)
  !                 de meshsize values XXDXHAT,XXDYHAT(:,1:7)      ! (MODD_COORD)
        
  !    apres cette lecture, les variables suivantes sont disponibles:
  !    NIMAX,NJMAX,NKMAX , apres SETDIM, LCARTESIAN, LTHINSHELL,CSTORAGE_TYPE,
  !    NGRID
  !    XXHAT(IIU)   pour la grille de U
  !    XYHAT(IJU)   pour la grille de V
  !    XZHAT(IIU)
  !    XMAP(IIU,IJU)
  !    XLAT(IIU,IJU)   pour la grille de masse
  !    XLON(IIU,IJU)   pour la grille de masse
  !    XDXHAT(IIU),XDYHAT(IJU)
  !    XZS(IIU,IJU)
  !    XZZ(IIU,IJU,IKU)  pour la grille W
  !    TDTMOD,TDTCUR,TDTEXP,TDTSEG,
  !    NSTOP,NOUT_TIMES,NOUT_NUMB, XTSTEP,XSEGLEN,
  ! 
  CALL COMPCOORD_FORDIACHRO(4)  ! NGRID set to 4 then XZZ is the true height
                                !of w-point as in the model
  if (KVERBIA >0)then
    print *,'*after COMPCOORD_FORDIACHRO(4)'
  endif
  !
  ! indiquera  au prochain appel de READVAR que le fichier courant 
  !est deja ouvert  (lecture du champ sans init des modules)
  HFLAGFILE(1:3)='NOP'
ENDIF
!-------------------------------------------------------------------------------
!
!*       4.    LIST OF GROUPS
!              --------------
!
IF(HLABELCHAMP(1:5)=='GROUP')THEN
  print *,'*following groups are present in the file ',TRIM(HFILENAME)
  ILENDIM=1
  YRECFM='MENU_BUDGET.DIM'
  CALL FMREAD(HFILENAME,YRECFM,CLUOUTDIAS(NBFILES),ILENDIM,ILENG,&
  IGRID,ILENCH,YCOMMENT,NRESPDIAS(NBFILES))
  IF(NRESPDIAS(NBFILES) == -47)THEN
    print *,' No record MENU_BUDGET '
    RETURN
  ENDIF
  ALLOCATE(ITABCHAR(ILENG))
  YRECFM='MENU_BUDGET'
  CALL FMREAD(HFILENAME,YRECFM,CLUOUTDIAS(NBFILES),ILENG,ITABCHAR, &
  IGRID,ILENCH,YCOMMENT,NRESPDIAS(NBFILES))
  IGROUP=ILENG/16
  ALLOCATE(YGROUP(IGROUP))
  DO JLOOP=1,IGROUP
    DO JI= 1,16
      YGROUP(JLOOP)(JI:JI)=CHAR(ITABCHAR(16*(JLOOP-1)+JI))
    ENDDO
  ENDDO
  print *,'****************************** GROUPS *****************************'
  print 100,(YGROUP(JLOOP),JLOOP=1,IGROUP)
100 FORMAT(1X,5A15)
  DEALLOCATE(ITABCHAR,YGROUP)
!
ELSE
!-------------------------------------------------------------------------------
!
!*       5.    TEST IF GROUP EXISTS 
!              --------------------
!
YGP=HLABELCHAMP
CALL VERIF_GROUP(HFILENAME,CLUOUTDIAS(NBFILES),YGP)
IF(LPBREAD)THEN
  print *,' ****READVAR: Groupe ',TRIM(YGP),' inconnu dans le fichier ', &
          TRIM(HFILENAME)
  KRETCODE=2
  LPBREAD=.FALSE.
  print *,'---------'
  RETURN
ENDIF
CFILEDIA=HFILENAME
!
!-------------------------------------------------------------------------------
!
!*       6.   READ GROUP
!             ----------
!
if (KVERBIA >0)then
  print *,'*before READ_DIACHRO'
endif
!
CALL READ_DIACHRO(HFILENAME,CLUOUTDIAS(NBFILES),YGP)
if (KVERBIA >0)then
  print *,'*after READ_DIACHRO'
endif
!
! lecture d'un enregistrement de nom CGROUP (en fait plusieurs enregistrements 
!lus dans les variables suivantes:
!CGROUP//'.TYPE' => CTYPE('CART','MASK','SPXY','SSOL','RSPL','DRST','RAPL')
                                                     ! MODD_TYPE_AND_LH
!CGROUP//'.DIM'  si CTYPE='CART','MASK','SPXY'
!             NIL,NJL,NKL,NIH,NJH,NKH,LICP,LJCP,LKCP ! MODD_TYPE_AND_LH
! = zoom inside the complete x-y-zgrid
!                appel de ALLOC_FORDIACHRO pour allouer les var. suivantes
!CGROUP//'.TITRE'  =>CTITRE(p)                       ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.UNITE'  =>CUNITE(p)                       ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.COMMENT' =>COMMENT(p)                     ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.PROCp' =>XVAR(i,j,k,t,n,p),NGRIDIA(p)     ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.TRAJT' =>XTRAJT(t,n)                      ! MODD_ALLOC_FORDIACHRO
! 
!CGROUP//'.TRAJX' =>XTRAJX(k,t,n)  optional          ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.TRAJY' =>XTRAJY(k,t,n)    "               ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.TRAJZ' =>XTRAJZ(k,t,n)    "               ! MODD_ALLOC_FORDIACHRO
!CGROUP//'.MASK'  =>XMASK(i,j,1,t,n,1)  " (si CTYPE='MASK')! MODD_ALLOC_FORDIACHRO
!CGROUP//'.DATIM' =>XDATIME(16,t)                    ! MODD_ALLOC_FORDIACHRO
! EXP.YEAR=XDATIME(1,t); EXP.MONTH=XDATIME(2,t) 
! EXP.DAY=XDATIME(3,t) ; EXP.TIME=XDATIME(4,t)
! SEG.YEAR=XDATIME(5,t); SEG.MONTH=XDATIME(6,t)
! SEG.DAY=XDATIME(7,t);  SEG.TIME=XDATIME(8,t)
! MOD.YEAR=XDATIME(9,t); MOD.MONTH=XDATIME(10,t) 
! MOD.DAY=XDATIME(11,t) ; MOD.TIME=XDATIME(12,t)
! CUR.YEAR=XDATIME(13,t); CUR.MONTH=XDATIME(14,t)
! CUR.DAY=XDATIME(15,t);  CUR.TIME=XDATIME(16,t)
!

! Passage a des unites plus pertinentes pour calculs si necessaire
CALL TO_COMPUTING_UNITS(YGP,CUNITE(1))
!
! Traitement d un champ eventuellement zoome
!
IF (CTYPE == 'CART' .AND. .NOT. LICP .AND. .NOT. LJCP ) THEN
  IF( SIZE(XVAR,1) /= SIZE(XZZ,1) .OR. SIZE(XVAR,2) /= SIZE(XZZ,2) )THEN
        ! replace le zoom dans le domaine total avant tout autre traitement
        !pour avoir les memes indices pour XLON,XLAT et ZHAT et XVAR
        if (KVERBIA > 0 ) then
          print *,' Replace un champ zoome dans le domaine total:'
          print'(A19,3(I4,X))','NIMAX,NJMAX,NKMAX= ',NIMAX,NJMAX,NKMAX
          print'(A25,6(I4,X))','nil,nih,njl,njh,nkl,nkh= ',nil,nih,njl,njh,nkl,nkh
        endif
        ! sauve XVAR
        ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
                          size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
        ZVARSAVE=XVAR
        if (KVERBIA > 0 ) then
          print *,'dimensions 4 5 6 :'
          print'(3(I5,x))',size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)
        endif
        DEALLOCATE(XVAR)
        idim3=SIZE(XZZ,3)
        IF (SIZE(ZVARSAVE,3) /= SIZE(XZZ,3)) THEN
          IF (SIZE(ZVARSAVE,3)/=1 )THEN
            !champ 3D zoome selon k
            idim3=SIZE(XZZ,3)
          ELSE
            !champ 2D
            idim3=SIZE(ZVARSAVE,3)
          ENDIF
        ENDIF
        ! nouveau XVAR= domaine total
        ALLOCATE(XVAR(SIZE(XZZ,1),SIZE(XZZ,2),idim3,&
                      SIZE(ZVARSAVE,4),SIZE(ZVARSAVE,5),SIZE(ZVARSAVE,6)))
        ! init seulement du zoom lu stocke dans ZVARSAVE
        XVAR=XUNDEF
        XVAR(nil:nih,njl:njh,nkl:nkh,:,:,:)=ZVARSAVE(:,:,:,:,:,:)
        DEALLOCATE (ZVARSAVE)

        !! GAELLE mars 2010
        IF (nil /= 1) THEN
           XVAR(nil-1,:,:,:,:,:)=XVAR(nil,:,:,:,:,:)
        ENDIF
        IF (nih /= SIZE(XZZ,1) ) THEN 
            XVAR(nih+1,:,:,:,:,:)= XVAR(nih,:,:,:,:,:)
        ENDIF
        IF (njl /= 1) THEN
           XVAR(:,njl-1,:,:,:,:)=XVAR(:,njl,:,:,:,:)
        ENDIF
        IF(njh /= SIZE(XZZ,2) ) THEN
           XVAR(:,njh+1,:,:,:,:)=XVAR(:,njh,:,:,:,:)
        ENDIF
        IF (nkl /= 1) THEN
           XVAR(:,:,nkl-1,:,:,:)=XVAR(:,:,nkl,:,:,:)
        ENDIF
        IF(nkh /= idim3) THEN
           XVAR(:,:,nkh+1,:,:,:)=XVAR(:,:,nkh,:,:,:)
        ENDIF
        !! GAELLE mars 2010

!     ENDIF
  ENDIF
ENDIF
!
! Traitement d un champ partiellement ecrit
!
IF (CTYPE == 'CART' .AND. .NOT. LKCP) THEN
  IF( SIZE(XVAR,3) /= SIZE(XZZ,3) )THEN
        if (KVERBIA > 0 ) then
          print *,' Replace un champ partiellement ecrit dans le domaine total:'
          print'(A7,I3)','NKMAX= ',NKMAX
          print'(A9,2(I3,X))','nkl,nkh= ',nkl,nkh
        endif
    ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
                      size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
    ZVARSAVE=XVAR
    IF (SIZE(ZVARSAVE,3)/=1 )THEN
      !champ 3D zoome selon k
      idim3=SIZE(XZZ,3)
    ELSE
      !champ 2D
      idim3=SIZE(ZVARSAVE,3)
    ENDIF
    print*,idim3
    DEALLOCATE(XVAR)
    ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),idim3,&
                  SIZE(ZVARSAVE,4),SIZE(ZVARSAVE,5),SIZE(ZVARSAVE,6)))
    XVAR=XUNDEF
    XVAR(:,:,nkl:nkh,:,:,:)=ZVARSAVE(:,:,:,:,:,:)
    !! GAELLE mars 2010
    IF (nkl /= 1) THEN
    XVAR(:,:,nkl-1,:,:,:)=XVAR(:,:,nkl,:,:,:)
    ENDIF
    print*,nkh,idim3
    IF(nkh /= idim3) THEN
    XVAR(:,:,nkh+1,:,:,:)=XVAR(:,:,nkh,:,:,:)
    ENDIF
    !! GAELLE mars 2010

    DEALLOCATE (ZVARSAVE)
  ENDIF
ENDIF
!
NREADIL=1 ; NREADIH=SIZE(XVAR,1)
NREADJL=1 ; NREADJH=SIZE(XVAR,2)
NREADKL=1 ; NREADKH=SIZE(XVAR,3)
IF (CTYPE=='CART' .OR. CTYPE=='MASK' .OR. CTYPE=='SPXY') THEN 
  IF (.NOT. LICP) THEN
    NREADIL=NIL ; NREADIH=NIH
  END IF
  IF (.NOT. LJCP) THEN
    NREADJL=NJL ; NREADJH=NJH
  END IF
  IF (.NOT. LKCP) THEN
    NREADKL=NKL ; NREADKH=NKH
  END IF
ENDIF
if (KVERBIA >= 0) then
  print*,'End of READVAR: the group ',&
          TRIM(YGP),' of file ',TRIM(HFILENAME),&
          ' is available in the XVAR array with sizes'
  print'(A4,I4,5(A5,I4))','  1:',SIZE(XVAR,1),',1:',SIZE(XVAR,2),',1:',SIZE(XVAR,3),&
           ',1:',SIZE(XVAR,4),',1:',SIZE(XVAR,5),',1:',SIZE(XVAR,6)
  IF (CTYPE=='CART' .OR. CTYPE=='MASK' .OR. CTYPE=='SPXY') THEN 
    print'(A90,6(I4,A))',&
         '(initialized in the zoom (NREADIL:NREADIH,NREADJL:NREADJH,NREADKL:NREADKH)= ',&
         NREADIL,':',NREADIH,',',NREADJL,':',NREADJH,',',NREADKL,':',NREADKH,')'
  END IF
endif
!
ENDIF  ! HLABELCHAMP(1:5)/='GROUP'
print *,'---------'
!
END SUBROUTINE READVAR
