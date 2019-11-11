!-----------------------------------------------------------------
!     ####################
      SUBROUTINE TRATRAJ3D
!     ####################
!
!!****  *TRATRAJ3D* - (Demande Joel Stein,Nicole Asencio, Francois Gheusi)
!!                    (Avril 00)
!!
!!    PURPOSE
!!    -------
!       Materialisation du positionnement de particules a divers instants
!       issues d'une position initiale connue ,
!       par transport de leurs coordonnees initiales dans les tableaux
!       scalaires SVx1, SVx2, SVx3
!
!       Conjointement :
!         ecriture a chaque point de la trajectoire d'un champ donne 
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron  et J. Stein  * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       12/04/00
!!      21/11/03  J. Stein Modification of the test for the field
!!                    computation along the backward trajectories
!!      10/03/04  JD  Ajout titres standard et possibilite de modification de
!!                    ceux-ci
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TRAJ3D
USE MODD_TITLE
USE MODD_TIT
USE MODI_INTERPXYZ
USE MODD_MASK3D
USE MODD_RESOLVCAR
USE MODD_CONF
USE MODD_COORD
USE MODD_GRID1
USE MODD_NMGRID
USE MODD_DIM1
USE MODD_PARAMETERS
USE MODD_SEVERAL_RECORDS
USE MODD_FILES_DIACHRO
USE MODD_ALLOC_FORDIACHRO
USE MODI_REALLOC_AND_LOAD
USE MODN_NCAR
USE MODD_CTL_AXES_AND_STYL
USE MODN_PARA
USE MODI_TIT_TRA3D
USE MODI_WRITEDIR
!
IMPLICIT NONE
!
COMMON/COLAREA/ICOL(300)
!
!*       0.1   Local variables
!
INTEGER           :: JKLOOP,JILOOP , JJLOOP, J, ID, IGRID, JTLOOP, JI
INTEGER           :: IIB, IIE, IJB, IJE, IKB, IKE
INTEGER           :: ICL, ICOL, ILOOP, IDEB, IFIN, INUM, IRESP
!
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZSVM1, ZSVM2, ZSVM3, ZCHAMP
REAL :: ZVL, ZVR, ZVB, ZVT, ZWL, ZWR, ZWB, ZWT
REAL :: ZVLL, ZVRL, ZVBL, ZVTL
REAL :: ZMINZ, ZMAXZ, ZINTZ, ZISO
REAL,DIMENSION(300) :: ZLEV
!CHARACTER(LEN=8),DIMENSION(300) :: YLLBS
CHARACTER(LEN=16) :: YGROUP
CHARACTER(LEN=75) :: YCAR
CHARACTER(LEN=12) :: YCHAMP
CHARACTER(LEN=100),SAVE  :: YTEM2
CHARACTER(LEN=110),SAVE  :: YTEM1
INTEGER  :: JPART,ICOLOR
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZXPOS,ZYPOS,ZZPOS, ZCHAMP_POS  ! positions aux
!   instants correspondants aux differents fichiers
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: GPART_IN   ! particule in the
                                                   ! computational domain?
!
!-------------------------------------------------------------------------------
IGRID=NMGRID
NMGRID=1

!
! boucle generale sur les fichiers
!
DO JTLOOP=1,NBFILES
! on lit les champs X0,Y0 et Z0 de la trajectoire pour tous les fichiers
!
! partie selon X
  YGROUP='LGXM'
  CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
  IF(LPBREAD)THEN
    YGROUP='LGXT'
    LPBREAD=.FALSE.
    CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
    IF(LPBREAD)THEN
      YGROUP='SVM001'
      LPBREAD=.FALSE.
      CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
      IF(LPBREAD)THEN
        YGROUP='SVT001'
        LPBREAD=.FALSE.
        CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
        IF(LPBREAD)THEN
          YGROUP='SVM1'
          LPBREAD=.FALSE.
          CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
          IF(LPBREAD)THEN
            YGROUP='SVT1'
            LPBREAD=.FALSE. 
            CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
            !
            IF(LPBREAD)THEN
          print *,' Absence de variable LGXM, SVM001, LGXT ou SVT001 .. Operation impossible'
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  !
  IF (LGROUP) THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    CALL READ_DIACHRO(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
  ENDIF
  !
  IF (.NOT. ALLOCATED(ZSVM1)) THEN
    ALLOCATE(ZSVM1(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
    ZSVM1=11111.
  ENDIF
  IF(MAXVAL(XXHAT)/MAXVAL(XVAR) > 1.E2)THEN
    print *,' ** Tratraj3D MAXVAL(XXHAT),MAXVAL(XVAR),*1000(KM->M) ',MAXVAL(XXHAT),MAXVAL(XVAR)
    WHERE(XVAR(:,:,:,1,1,1) /= XSPVAL)
      ZSVM1(:,:,:)=XVAR(:,:,:,1,1,1)*1000.
    ELSEWHERE
      ZSVM1(:,:,:)=XVAR(:,:,:,1,1,1)
    ENDWHERE
  ELSE
    ZSVM1(:,:,:)=XVAR(:,:,:,1,1,1)
  ENDIF
  !
  ! Chargement clegend clegend2
  CALL RESOLV_TIMES(1)
  YTEM2=' '
  YTEM1=' '
  YTEM2=CLEGEND2
  ! Elimination volontaire de 104 a 108 charge ds resolv_times pour RS
  YTEM1=CLEGEND(1:103)
  !
  IF(.NOT.LFIC1)THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
!    CALL REALLOC_AND_LOAD(YGROUP)
    IF(LPBREAD)THEN
      print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
      ' L''UN DES FICHIERS '
      IF(ALLOCATED(XVAR))THEN
        CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      ENDIF
      RETURN
    ENDIF
  ENDIF
!
! partie selon Y
  YGROUP='LGYM'
  CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
  IF(LPBREAD)THEN
    YGROUP='LGYT'
    LPBREAD=.FALSE.
    CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
    IF(LPBREAD)THEN
      YGROUP='SVM002'
      LPBREAD=.FALSE.
      CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
      IF(LPBREAD)THEN
        YGROUP='SVT002'
        LPBREAD=.FALSE.
        CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
        IF(LPBREAD)THEN
          YGROUP='SVM2'
          LPBREAD=.FALSE.
          CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
          IF(LPBREAD)THEN
            YGROUP='SVT2'
            LPBREAD=.FALSE. 
            CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
            !
            IF(LPBREAD)THEN
          print *,' Absence de variable LGYM, SVM002, LGYT ou SVT002 .. Operation impossible'
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  !
  IF (LGROUP) THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    CALL READ_DIACHRO(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
  ENDIF
  !
  IF (.NOT. ALLOCATED(ZSVM2)) THEN
    ALLOCATE(ZSVM2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
    ZSVM2=11111.
  ENDIF
  IF(MAXVAL(XYHAT)/MAXVAL(XVAR) > 1.E2)THEN
    print *,' ** Tratraj3D MAXVAL(XYHAT),MAXVAL(XVAR),*1000(KM->M) ',MAXVAL(XYHAT),MAXVAL(XVAR)
    WHERE(XVAR(:,:,:,1,1,1) /= XSPVAL)
      ZSVM2(:,:,:)=XVAR(:,:,:,1,1,1)*1000.
    ELSEWHERE
      ZSVM2(:,:,:)=XVAR(:,:,:,1,1,1)
    ENDWHERE
  ELSE
    ZSVM2(:,:,:)=XVAR(:,:,:,1,1,1)
  ENDIF
  !
  IF(.NOT.LFIC1)THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    IF(LPBREAD)THEN
      print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
      ' L''UN DES FICHIERS '
      IF(ALLOCATED(XVAR))THEN
        CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      ENDIF
      RETURN
    ENDIF
  ENDIF
!
! partie selon Z
  YGROUP='LGZM'
  CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
  IF(LPBREAD)THEN
    YGROUP='LGZT'
    LPBREAD=.FALSE.
    CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
    IF(LPBREAD)THEN
      YGROUP='SVM003'
      LPBREAD=.FALSE.
      CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
      IF(LPBREAD)THEN
        YGROUP='SVT003'
        LPBREAD=.FALSE.
        CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
        IF(LPBREAD)THEN
          YGROUP='SVM3'
          LPBREAD=.FALSE. 
          CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
          IF(LPBREAD)THEN
            YGROUP='SVT3'
            LPBREAD=.FALSE. 
            CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
            !
            IF(LPBREAD)THEN
          print *,' Absence de variable LGZM, SVM003, LGZT ou SVT003 .. Operation impossible'
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  !
  IF (LGROUP) THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    CALL READ_DIACHRO(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),YGROUP)
  ENDIF
  !
  IF (.NOT. ALLOCATED(ZSVM3)) THEN
    ALLOCATE(ZSVM3(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
    ZSVM3=11111.
  ENDIF
  IF(MAXVAL(XZHAT)/MAXVAL(XVAR) > 1.E2)THEN
    print *,' ** Tratraj3D MAXVAL(XZHAT),MAXVAL(XVAR),*1000(KM->M) ',MAXVAL(XZHAT),MAXVAL(XVAR)
    WHERE(XVAR(:,:,:,1,1,1) /= XSPVAL)
      ZSVM3(:,:,:)=XVAR(:,:,:,1,1,1)*1000.
    ELSEWHERE
      ZSVM3(:,:,:)=XVAR(:,:,:,1,1,1)
    ENDWHERE
  ELSE
    ZSVM3(:,:,:)=XVAR(:,:,:,1,1,1)
  ENDIF
  !
  IF(.NOT.LFIC1)THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    IF(LPBREAD)THEN
      print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
      ' L''UN DES FICHIERS '
      IF(ALLOCATED(XVAR))THEN
        CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      ENDIF
      RETURN
    ENDIF
  ENDIF
  IF (JTLOOP==1) THEN
    ! on calcule ici les grilles verticales a cause du cas du champ ALT
    ! qui pose un probleme car il est situe sur un niveau de w
    IIB=1+JPHEXT; IIE=SIZE(ZSVM1,1)-JPHEXT
    IJB=1+JPHEXT; IJE=SIZE(ZSVM1,2)-JPHEXT
    IKB=1+JPVEXT; IKE=SIZE(ZSVM1,3)-JPVEXT
    !
    ! Calcul des altitudes pour la grille 1 dans XZZ
    !
    CALL COMPCOORD_FORDIACHRO(1)
    !
  ENDIF
!
! on lit un champ supplementaire pour le tracer sur la trajectoire
!
  IF (LTRAJ_GROUP) THEN
   IF ( CTRAJ_GROUP=='ALT') THEN
     IF (.NOT. ALLOCATED(ZCHAMP)) THEN
      ALLOCATE(ZCHAMP(SIZE(ZSVM3,1),SIZE(ZSVM3,2),SIZE(ZSVM3,3)))
      ZCHAMP=11111.
     ENDIF
     IF (JTLOOP==1) ZCHAMP(:,:,:)=XZZ(:,:,:)
   ELSE
    CALL VERIF_GROUP(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),CTRAJ_GROUP)
    IF(LPBREAD)THEN
      print *,' Absence de variable CTRAJ_GROUP .. Operation impossible'
      RETURN
    ENDIF
    !
    IF (LGROUP) THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      CALL READ_DIACHRO(CFILEDIAS(JTLOOP),CLUOUTDIAS(JTLOOP),CTRAJ_GROUP)
    ENDIF
    !
    IF (.NOT. ALLOCATED(ZCHAMP)) THEN
      ALLOCATE(ZCHAMP(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
      ZCHAMP=11111.
    ENDIF
    !
    ZCHAMP(:,:,:)=XVAR(:,:,:,1,1,1)
    !
    IF(.NOT.LFIC1)THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      IF(LPBREAD)THEN
        print *,' REQUETE IMPOSSIBLE .',CTRAJ_GROUP,' N''EXISTE PAS DANS', &
        ' L''UN DES FICHIERS '
        IF(ALLOCATED(XVAR))THEN
          CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
        ENDIF
        RETURN
      ENDIF
    ENDIF
   END IF
  ELSE
    IF (.NOT. ALLOCATED(ZCHAMP)) ALLOCATE(ZCHAMP(0,0,0))
  ENDIF
!
! on recherche la valeur R0 d'origine pour le point courant R
!
  IF (JTLOOP==1) THEN
   ALLOCATE(ZXPOS(NPART,NBFILES+1))
   ALLOCATE(ZYPOS(NPART,NBFILES+1))
   ALLOCATE(ZZPOS(NPART,NBFILES+1))
   ALLOCATE(GPART_IN(NPART,NBFILES+1))
   GPART_IN=.TRUE.
   IF (LTRAJ_GROUP) THEN
     ALLOCATE(ZCHAMP_POS(NPART,NBFILES+1))
   ELSE
!!!Octobre 2001
     ALLOCATE(ZCHAMP_POS(NPART,NBFILES+1))
!    ALLOCATE(ZCHAMP_POS(1,2))
!!!Octobre 2001
!    ALLOCATE(ZCHAMP_POS(0,0))
   END IF
   !
   ZXPOS(:,1)=XXPART(1:NPART)
   ZYPOS(:,1)=XYPART(1:NPART)
   ZZPOS(:,1)=XZPART(1:NPART)
   !
   DO JPART=1,NPART
     IF (ZXPOS(JPART,1).LT.XXX(IIB,1) .OR. ZXPOS(JPART,1).GT.XXX(IIE,1) .OR.   &
         ZYPOS(JPART,1).LT.XXY(IJB,1) .OR. ZYPOS(JPART,1).GT.XXY(IJE,1)        &
        ) THEN
       ZXPOS(JPART,1)=MIN(XXX(IIE,1),MAX(XXX(IIB,1),ZXPOS(JPART,1)))
       ZYPOS(JPART,1)=MIN(XXY(IJE,1),MAX(XXY(IJB,1),ZYPOS(JPART,1)))
       print *,' la particule ',JPART,' est sortie du domaine'
       print *,'nouvelles valeurs de XXPART et XYPART:'
       print *,'XXPART=',ZXPOS(JPART,1),'XYPART=',ZYPOS(JPART,1)
     END IF
   END DO
  ENDIF
!
!
  DO JPART=1,NPART
    !
    IF(GPART_IN(JPART,JTLOOP)) THEN
         ! the particule is in the simulation domain
      CALL INTERPXYZ(ZSVM1(:,:,:),               &
                     ZSVM2(:,:,:),               &
                     ZSVM3(:,:,:),               &
                     ZCHAMP(:,:,:),              &
                     ZXPOS(JPART,JTLOOP),        &
                     ZYPOS(JPART,JTLOOP),        &
                     ZZPOS(JPART,JTLOOP),        &
                     XXX(2,1),XXY(2,1),          & 
                     XXDXHAT(3,1),XXDYHAT(3,1),  &
                     XZZ,LTRAJ_GROUP,            &
                     ZXPOS(JPART,JTLOOP+1),      &
                     ZYPOS(JPART,JTLOOP+1),      &
                     ZZPOS(JPART,JTLOOP+1),      &
                     ZCHAMP_POS(JPART,JTLOOP)  )
         !
      IF (ZXPOS(JPART,JTLOOP+1).LT.XXX(IIB,1) .OR.   &
          ZXPOS(JPART,JTLOOP+1).GT.XXX(IIE,1) .OR.   &
          ZYPOS(JPART,JTLOOP+1).LT.XXY(IJB,1) .OR.   &
          ZYPOS(JPART,JTLOOP+1).GT.XXY(IJE,1)        &
         )  THEN
         ! it is the first time the particule has been gone out
        GPART_IN(JPART,JTLOOP+1)=.FALSE.
        ZXPOS(JPART,JTLOOP+1)=ZXPOS(JPART,JTLOOP)
        ZYPOS(JPART,JTLOOP+1)=ZYPOS(JPART,JTLOOP)
        ZZPOS(JPART,JTLOOP+1)=ZZPOS(JPART,JTLOOP)
        print *,'la particule ',JPART,' est sortie du domaine apres ',JTLOOP+1,' avances' 
      ENDIF
    ELSE
         ! the particule is out of the simulation domain
        GPART_IN(JPART,JTLOOP+1)=.FALSE.
        ZXPOS(JPART,JTLOOP+1)=ZXPOS(JPART,JTLOOP)
        ZYPOS(JPART,JTLOOP+1)=ZYPOS(JPART,JTLOOP)
        ZZPOS(JPART,JTLOOP+1)=ZZPOS(JPART,JTLOOP)
        ZCHAMP_POS(JPART,JTLOOP)=ZCHAMP_POS(JPART,JTLOOP-1)         
    END IF
    ! fin de la boucle sur les particules
  ENDDO
!
! fin de la boucle sur les instants de la trajectoire
!
ENDDO
!
DEALLOCATE(ZSVM1,ZSVM2,ZSVM3,ZCHAMP,GPART_IN)   ! dealloc des champs
!
! sortie des trajectoires
IF(LPRINT)THEN
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  ILOOP=NPART/5
  IF(ILOOP * 5 < NPART)ILOOP=ILOOP+1
ENDIF
DO JTLOOP=1,NBFILES+1
  print*,'*****************'
  print*,'JTLOOP= ', JTLOOP
  print*,'*****************'
  print*,'XPOS= ',ZXPOS(1:NPART,JTLOOP)
  print*,'YPOS= ',ZYPOS(1:NPART,JTLOOP)
  print*,'ZPOS= ',ZZPOS(1:NPART,JTLOOP)
  IF (LTRAJ_GROUP) print*,'CHAMPPOS= ',ZCHAMP_POS(1:NPART,JTLOOP)
  IF(LPRINT)THEN
    WRITE(INUM,'(A,I3)') 'LOOP= ',JTLOOP
    DO JI=1,ILOOP
      IF (JI==1) THEN
        IDEB=1 ; IFIN=4
      ELSE  
        IDEB=IFIN+1 ; IFIN=IFIN+5
      ENDIF
      IF (JI==ILOOP) THEN
        IFIN=NPART
      ENDIF
      IF (JI==1) THEN
        WRITE(INUM,'(A12,4(3X,E12.6))')' XPOS=',ZXPOS(IDEB:IFIN,JTLOOP)
      ELSE
        WRITE(INUM,'(4(E12.6,3X),E12.6)') ZXPOS(IDEB:IFIN,JTLOOP)
      ENDIF
    END DO
    DO JI=1,ILOOP
      IF (JI==1) THEN
        IDEB=1 ; IFIN=4
      ELSE  
        IDEB=IFIN+1 ; IFIN=IFIN+5
      ENDIF
      IF (JI==1) THEN
        WRITE(INUM,'(A12,4(3X,E12.6))')' YPOS=',ZYPOS(IDEB:IFIN,JTLOOP)
      ELSE
        WRITE(INUM,'(4(E12.6,3X),E12.6)') ZYPOS(IDEB:IFIN,JTLOOP)
      ENDIF
    END DO
    DO JI=1,ILOOP
      IF (JI==1) THEN
        IDEB=1 ; IFIN=4
      ELSE  
        IDEB=IFIN+1 ; IFIN=IFIN+5
      ENDIF
      IF (JI==1) THEN
        WRITE(INUM,'(A12,4(3X,E12.6))')' ZPOS=',ZZPOS(IDEB:IFIN,JTLOOP)
      ELSE
        WRITE(INUM,'(4(E12.6,3X),E12.6)') ZZPOS(IDEB:IFIN,JTLOOP)
      ENDIF
      IF (JI==ILOOP) WRITE(INUM,*)
    END DO
  ENDIF
END DO
!
!-------------------------------------------------------------------------------
!
!!!!!!!!!!!!JOEL!!!!!!!!!!
!!!!!!!!!!!!JOEL!!!!!!!!!!
! Visualisation des trajectoires sur XY, XZ, YZ
!!!!!!!!!!!!JOEL!!!!!!!!!!
!!!!!!!!!!!!JOEL!!!!!!!!!!
!
! Recuperation de la fenetre d'affichage courante pour restauration en fin de
! routine
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! 
! Determination de NIINF NJINF NISUP NJSUP si non initialises par l'utilisateur
IF(NIINF == 0 .AND. NISUP == 0 .AND. NJINF == 0 .AND. NJSUP == 0)THEN
  CALL RESOLV_NIJINF_NIJSUP
ENDIF

!
!!!!!! XY 
!
YCAR(1:LEN_TRIM(YCAR))=' '
WRITE(YCAR,'(''TRAJ **XY** '')')
IF( LTRAJ_GROUP) THEN
  ! car TIT_TRA3D ne trace rien sur la 1e image dans le cas LTRAJ_GROUP ...!
  CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
  CALL PCSETC('FC','/')
  CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
  CALL PCSETC('FC',':')
ELSE
  CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)
ENDIF

IF(LDATFILE)CALL DATFILE_FORDIACHRO

IF(LCARTESIAN)THEN
  CALL DEFENETRE
ELSE
  ! trace de la grille lat-lon
  CALL GSLWSC(1.)
  CALL GSTXCI(1)
  CALL GSPLCI(1)
  CALL BCGRD_FORDIACHRO(2)
  !CALL BCGRD_FORDIACHRO(1)
ENDIF
!
! couleur en fct de l alt ZZPOS (15 intervalles)
ICL=15
CALL COLOR_FORDIACHRO(ICL+2,1)
CALL TABCOL_FORDIACHRO
ZMAXZ=MAXVAL(ZZPOS) ; ZMINZ=MINVAL(ZZPOS)
ZINTZ=NINT(ZMAXZ-ZMINZ)/15
IF(ZMINZ + ICL*ZINTZ <= ZMAXZ)ICL=ICL+1
CALL CPSETI('NCL',ICL)
CALL CPSETI('CLS',0)
ZISO=ZMINZ-ZINTZ
DO JI=1,ICL
  CALL CPSETI('PAI',JI)
  CALL CPSETI('AIA',JI+1)
  CALL CPSETI('AIB',JI)
  ZISO=ZISO+ZINTZ
  IF(ABS(ZISO)<1.E-20)ZISO=0.
  CALL CPSETR('CLV',ZISO)
  CALL CPSETR('CLU',1.)
  ZLEV(JI)=ZISO
  !CALL GENFORMAT_FORDIACHRO(ZISO,YLLBS(JI))
  ICOL(JI)=JI
END DO
!
IF (.NOT.LCOLINE) THEN
  print *,' LCOLINE=F: Retro-trajectoires et marqueurs noirs dans le plan XY'
ENDIF
!
CALL GSLWSC(3.)
DO JPART=1,NPART
  CALL GSMK(4)  
  IF (.NOT.LCOLINE) THEN
    ICOLOR=1
    CALL GSPMCI(1)
  ELSE  
    ICOLOR= 1+ MOD((JPART-1),16)   ! boucle sur les 16 premieres couleurs 
    ! couleur du marker en fct de l alt ZZPOS
    IF(ZZPOS(JPART,1) <ZLEV(1))THEN
      CALL GSPMCI(1)
    ELSEIF(ZZPOS(JPART,1) >=ZLEV(ICL))THEN
      CALL GSPMCI(ICL+1)
    ELSE
      DO JI=1,ICL-1
        IF(ZZPOS(JPART,1) >= ZLEV(JI) .AND. &
          ZZPOS(JPART,1) < ZLEV(JI+1))THEN
          CALL GSPMCI(JI+1)
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  CALL GSTXCI(ICOLOR)
  CALL GSPLCI(ICOLOR)
  CALL GPM(1,ZXPOS(JPART,1),ZYPOS(JPART,1))
  CALL FRSTD(ZXPOS(JPART,1),ZYPOS(JPART,1))  
  CALL GSMK(3)  
  DO JTLOOP=2,NBFILES+1
    IF (LCOLINE) THEN ! couleur du marker en fct de l alt ZZPOS
      IF(ZZPOS(JPART,JTLOOP) <ZLEV(1))THEN
        CALL GSPMCI(1)
      ELSEIF(ZZPOS(JPART,JTLOOP) >=ZLEV(ICL))THEN
        CALL GSPMCI(ICL+1)
      ELSE
        DO JI=1,ICL-1
          IF(ZZPOS(JPART,JTLOOP) >= ZLEV(JI) .AND. &
             ZZPOS(JPART,JTLOOP) < ZLEV(JI+1))THEN
            CALL GSPMCI(JI+1)
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    CALL VECTD(ZXPOS(JPART,JTLOOP),ZYPOS(JPART,JTLOOP))
    CALL GPM(1,ZXPOS(JPART,JTLOOP),ZYPOS(JPART,JTLOOP))
  ENDDO
  CALL LASTD
ENDDO
! Trace des valeurs de ZZPOS en legende: A revoir...
!CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!CALL GSFAIS(1)
!CALL LBSETI('CBL',0)
!DO JI=1,ICL
!  YLLBS(JI)=ADJUSTL(YLLBS(JI))
!ENDDO
!IF(ZVR < .9)THEN
!  CALL LBLBAR_FORDIACHRO(1,ZVR+(MIN(ZVR+.2,1.)-ZVR)/10.,MIN(ZVR+.2,1.),ZVB,ZVT,ICL+1,.15,1.,ICOL,1,YLLBS,ICL,1)
!ELSE
!  CALL LBLBAR_FORDIACHRO(1,ZVR+(1.-ZVR)/10.,1.,ZVB,ZVT,ICL+1,.15,1.,ICOL,1,YLLBS,ICL,1)
!ENDIF
!
CALL FRAME
!
!
IF( LTRAJ_GROUP) THEN
  CALL GSLWSC(1.)
  CALL GSTXCI(1)
  CALL GSPLCI(1)
  CALL GSTXCI(1)
  YCAR(1:LEN_TRIM(YCAR))=' '
  WRITE(YCAR,'(''TRAJ **XY**   '',A16)') CTRAJ_GROUP
  CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
  !CALL PCSETC('FC','/')
  !CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
  !CALL PCSETC('FC',':')
  CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

  IF(LDATFILE)CALL DATFILE_FORDIACHRO

  IF(LCARTESIAN)THEN
    CALL DEFENETRE
  ELSE
    CALL BCGRD_FORDIACHRO(1)
  ENDIF

  CALL GSLWSC(3.)
  DO JPART=1,NPART
    CALL GSMK(4)  
    ICOLOR= 1+ MOD((JPART-1),16)   ! boucle sur les 16 premieres couleurs 
    CALL GSTXCI(ICOLOR)
    CALL GSPLCI(ICOLOR)
    CALL GSPMCI(ICOLOR)
    CALL GPM(1,ZXPOS(JPART,1),ZYPOS(JPART,1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   WRITE(YCHAMP,'(E10.5)') ZCHAMP_POS(JPART,1)
  ! CALL PLCHHQ(ZXPOS(JPART,1),ZYPOS(JPART,1),YCHAMP,10.,0.,-1.)
  
   WRITE(YCHAMP,CFMTRTRAJ) ZCHAMP_POS(JPART,1)
   CALL PLCHHQ(ZXPOS(JPART,1),ZYPOS(JPART,1),YCHAMP,NSZRTRAJ,0.,-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL FRSTD(ZXPOS(JPART,1),ZYPOS(JPART,1))  
    CALL GSMK(3)  
    DO JTLOOP=2,NBFILES+1
      CALL VECTD(ZXPOS(JPART,JTLOOP),ZYPOS(JPART,JTLOOP))
      CALL GPM(1,ZXPOS(JPART,JTLOOP),ZYPOS(JPART,JTLOOP))
      IF (JTLOOP<NBFILES+1) THEN
      ! le dernier point pour CHAMP se rapporte a l'echeance precedente
      ! donc il ne peut pas etre calcule et trace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !     WRITE(YCHAMP,'(E10.5)') ZCHAMP_POS(JPART,JTLOOP)
   !     CALL PLCHHQ(ZXPOS(JPART,JTLOOP),ZYPOS(JPART,JTLOOP),YCHAMP,10.,0.,-1.)

         WRITE(YCHAMP,CFMTRTRAJ) ZCHAMP_POS(JPART,JTLOOP)
         CALL PLCHHQ(ZXPOS(JPART,JTLOOP),ZYPOS(JPART,JTLOOP),YCHAMP,NSZRTRAJ,0.,-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF
    ENDDO
    CALL LASTD
  ENDDO
  !
  ! trace de la grille lat-lon
  CALL GSLWSC(1.)
  CALL GSTXCI(1)
  CALL GSPLCI(1)
  CALL BCGRD_FORDIACHRO(2)
  CALL FRAME
ENDIF
!
!!!!!! XZ 
!
CALL GSLWSC(1.)
CALL GSTXCI(1)
CALL GSPLCI(1)
CALL GSTXCI(1)
WRITE(YCAR,'(''TRAJ **XZ** '')')
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!CALL PCSETC('FC','/')
!CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!CALL PCSETC('FC',':')
CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

IF(LDATFILE)CALL DATFILE_FORDIACHRO

CALL SET(.1,.9,.1,.9,XXX(NIINF,1),XXX(NISUP,1), &
XHMIN,XHMAX,1)
CALL LABMOD('(F8.0)','(F6.0)',9,6,10,10,0,0,0)
CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NINT((XHMAX-XHMIN)/1000.),1,1,1,5,0.,0.)
!
CALL GSLWSC(3.)
DO JPART=1,NPART
  CALL GSMK(4)  
  ICOLOR= 1+ MOD((JPART-1),16)   ! boucle sur les 16 premieres couleurs 
  CALL GSPLCI(ICOLOR)
  CALL GSTXCI(ICOLOR)
  CALL GSPMCI(ICOLOR)
  CALL GPM(1,ZXPOS(JPART,1),ZZPOS(JPART,1))
  CALL FRSTD(ZXPOS(JPART,1),ZZPOS(JPART,1))  
  CALL GSMK(3)  
  DO JTLOOP=2,NBFILES+1
    CALL VECTD(ZXPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
    CALL GPM(1,ZXPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
  ENDDO
  CALL LASTD
ENDDO
!
CALL FRAME
!
!
IF (LTRAJ_GROUP) THEN
  CALL GSLWSC(1.)
  CALL GSTXCI(1)
  CALL GSPLCI(1)
  CALL GSTXCI(1)
  WRITE(YCAR,'(''TRAJ **XZ**     '',A16)') CTRAJ_GROUP
  CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!  CALL PCSETC('FC','/')
!  CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!  CALL PCSETC('FC',':')
  CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

  IF(LDATFILE)CALL DATFILE_FORDIACHRO

  CALL SET(.1,.9,.1,.9,XXX(NIINF,1),XXX(NISUP,1), &
  XHMIN,XHMAX,1)
  CALL LABMOD('(F8.0)','(F6.0)',9,6,10,10,0,0,0)
  CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NINT((XHMAX-XHMIN)/1000.),1,1,1,5,0.,0.)
  !
  CALL GSLWSC(3.)
  DO JPART=1,NPART
    CALL GSMK(4)  
    ICOLOR= 1+ MOD((JPART-1),16)   ! boucle sur les 16 premieres couleurs 
    CALL GSPLCI(ICOLOR)
    CALL GSTXCI(ICOLOR)
    CALL GSPMCI(ICOLOR)
    CALL GPM(1,ZXPOS(JPART,1),ZZPOS(JPART,1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  WRITE(YCHAMP,'(E10.5)') ZCHAMP_POS(JPART,1)
  !  CALL PLCHHQ(ZXPOS(JPART,1),ZZPOS(JPART,1),YCHAMP,10.,0.,-1.)
  
   WRITE(YCHAMP,CFMTRTRAJ) ZCHAMP_POS(JPART,1)
   CALL PLCHHQ(ZXPOS(JPART,1),ZZPOS(JPART,1),YCHAMP,NSZRTRAJ,0.,-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL FRSTD(ZXPOS(JPART,1),ZZPOS(JPART,1))  
    CALL GSMK(3)  
    DO JTLOOP=2,NBFILES+1
      CALL VECTD(ZXPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
      CALL GPM(1,ZXPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
      IF (JTLOOP<NBFILES+1) THEN
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    WRITE(YCHAMP,'(E10.5)') ZCHAMP_POS(JPART,JTLOOP)
  !    CALL PLCHHQ(ZXPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP),YCHAMP,10.,0.,-1.)
  
   WRITE(YCHAMP,CFMTRTRAJ) ZCHAMP_POS(JPART,JTLOOP)
   CALL PLCHHQ(ZXPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP),YCHAMP,NSZRTRAJ,0.,-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF
    ENDDO
    CALL LASTD
  ENDDO
  !
  CALL FRAME
END IF
!
!!!!!! YZ 
!
CALL GSLWSC(1.)
CALL GSPLCI(1)
CALL GSTXCI(1)
WRITE(YCAR,'(''TRAJ **YZ** '')')
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!CALL PCSETC('FC','/')
!CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!CALL PCSETC('FC',':')
CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)
  
IF(LDATFILE)CALL DATFILE_FORDIACHRO

CALL SET(.1,.9,.1,.9,XXY(NJINF,1),XXY(NJSUP,1), &
XHMIN,XHMAX,1)
CALL LABMOD('(F8.0)','(F6.0)',9,6,10,10,0,0,0)
CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NINT((XHMAX-XHMIN)/1000.),1,1,1,5,0.,0.)
!
CALL GSLWSC(3.)
DO JPART=1,NPART
  CALL GSMK(4)  
  ICOLOR= 1+ MOD((JPART-1),16)   ! boucle sur les 16 premieres couleurs 
  CALL GSPLCI(ICOLOR)
  CALL GSTXCI(ICOLOR)
  CALL GSPMCI(ICOLOR)
  CALL GPM(1,ZYPOS(JPART,1),ZZPOS(JPART,1))
  CALL FRSTD(ZYPOS(JPART,1),ZZPOS(JPART,1))  
  CALL GSMK(3)  
  DO JTLOOP=2,NBFILES+1
    CALL VECTD(ZYPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
    CALL GPM(1,ZYPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
  ENDDO
  CALL LASTD
ENDDO
!
CALL FRAME
!
IF (LTRAJ_GROUP) THEN
  CALL GSLWSC(1.)
  CALL GSPLCI(1)
  CALL GSTXCI(1)
    WRITE(YCAR,'(''TRAJ **YZ**     '',A16)') CTRAJ_GROUP
  CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!  CALL PCSETC('FC','/')
!  CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!  CALL PCSETC('FC',':')
  CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

  IF(LDATFILE)CALL DATFILE_FORDIACHRO

  CALL SET(.1,.9,.1,.9,XXY(NJINF,1),XXY(NJSUP,1), &
  XHMIN,XHMAX,1)
  CALL LABMOD('(F8.0)','(F6.0)',9,6,10,10,0,0,0)
  CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NINT((XHMAX-XHMIN)/1000.),1,1,1,5,0.,0.)
  !
  CALL GSLWSC(3.)
  DO JPART=1,NPART
    CALL GSMK(4)  
    ICOLOR= 1+ MOD((JPART-1),16)   ! boucle sur les 16 premieres couleurs 
    CALL GSPLCI(ICOLOR)
    CALL GSTXCI(ICOLOR)
    CALL GSPMCI(ICOLOR)
    CALL GPM(1,ZYPOS(JPART,1),ZZPOS(JPART,1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  WRITE(YCHAMP,'(E10.5)') ZCHAMP_POS(JPART,1)
  !CALL PLCHHQ(ZYPOS(JPART,1),ZZPOS(JPART,1),YCHAMP,10.,0.,-1.)

   WRITE(YCHAMP,CFMTRTRAJ) ZCHAMP_POS(JPART,1)
   CALL PLCHHQ(ZYPOS(JPART,1),ZZPOS(JPART,1),YCHAMP,NSZRTRAJ,0.,-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL FRSTD(ZYPOS(JPART,1),ZZPOS(JPART,1))  
    CALL GSMK(3)  
    DO JTLOOP=2,NBFILES+1
      CALL VECTD(ZYPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
      CALL GPM(1,ZYPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP))
      IF (JTLOOP<NBFILES+1) THEN
             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! WRITE(YCHAMP,'(E10.5)') ZCHAMP_POS(JPART,JTLOOP)
   ! CALL PLCHHQ(ZYPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP),YCHAMP,10.,0.,-1.)
   
   WRITE(YCHAMP,CFMTRTRAJ) ZCHAMP_POS(JPART,JTLOOP)
   CALL PLCHHQ(ZYPOS(JPART,JTLOOP),ZZPOS(JPART,JTLOOP),YCHAMP,NSZRTRAJ,0.,-1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ENDIF
    ENDDO
    CALL LASTD
  ENDDO
  !
  CALL FRAME
END IF
!
!
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)  
!
!
CALL GSTXCI(1)
CALL GSPLCI(1)
CALL GSLWSC(1.)
CALL GSLN(1)
DEALLOCATE(ZXPOS,ZYPOS,ZZPOS,ZCHAMP_POS)   ! dealloc des champs
NMGRID=IGRID
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
!
RETURN
!
END SUBROUTINE TRATRAJ3D 
