!-----------------------------------------------------------------
!     ####################
      SUBROUTINE TRAMASK3D
!     ####################
!
!!****  *TRAMASK3D* - (Demande Joel Stein,Nicole Asencio, Francois Gheusi)
!!                    (Mai 99)
!!
!!    PURPOSE
!!    -------
!       Materialisation du positionnement de particules a un instant donne
!       issues d'une position initiale connue ,
!       par transport de leurs coordonnees initiales dans les tableaux
!       scalaires SVx1, SVx2, SVx3
!       L'utilisateur definit une fenetre spatiale dans les limites
!       XXL= XXH= XYL= XYH= XZL= XZH= (en metres) correspondant a une
!       position initiale et recherche dans les =/= enr. de ces tableaux
!       (<-> a des termes d'evolution temporelle) les points correspondant
!       a cette fenetre -> occurences vraies d'un masque.
!
!       Si LMASK3D=T , visualisation de la projection de ces occurences
!       sur XY, XZ, YZ.
!
!       Conjointement :
!       thetae_msktop_ (Valeurs <-> surface des occurences.T. du masque)
!                       (a partir du sommet)
!       thetae_xyz__z_7000 (Extraction des valeurs de thetae corresp. aux
!       occurences .T. du masque en affectant aux autres points la valeur
!       XSPVAL puis trace comme habituellement d'une coupe horizontale, ici
!       d'altitude donnee)
!       thetae_sv3_5000,4000 (Trace d'une coupe horizontale d'altitudes
!       donnees SVx3. Le masque n'intervient pas dans ce cas)
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
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       29/04/99
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
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
!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!
!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!
USE MODD_CTL_AXES_AND_STYL
USE MODN_PARA
USE MODD_TRAJ3D
!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!
!!!!!!!!!!!!!!JOEL!!!!!!!!!!!!!!
USE MODD_TITLE
USE MODI_TIT_TRA3D
      USE MODD_ALLOC_FORDIACHRO

!
IMPLICIT NONE
!
!*       0.1   Local variables
!
INTEGER           :: JKLOOP,JILOOP , JJLOOP, J, JM, ID, IGRID, JTLOOP
INTEGER           :: IIB, IIE, IJB, IJE, IKB, IKE
INTEGER           :: IDBID
INTEGER           :: INUM,IRESP,ILOOP
!
REAL,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE :: ZSVM1, ZSVM2, ZSVM3
REAL :: ZVL, ZVR, ZVB, ZVT, ZWL, ZWR, ZWB, ZWT, ZX, ZY
REAL :: ZWLBID, ZWRBID, ZWBBID, ZWTBID
CHARACTER(LEN=16) :: YGROUP
CHARACTER(LEN=75) :: YCAR
CHARACTER(LEN=10) :: YFORMAX, YFORMAY
CHARACTER(LEN=2)  :: YNUMBER     ! number of the start for the lag. var
INTEGER           :: ILENTRIMSV3 ! length of the CGROUPSV3 var.
REAL,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE :: ZFIELD_LAG
LOGICAL  :: GLAG
CHARACTER(LEN=16) :: YSTO_CGROUPSV3 ! storage of CGROUPSV3
CHARACTER(LEN=100),SAVE  :: YTEM2
CHARACTER(LEN=110),SAVE  :: YTEM1
!
!-------------------------------------------------------------------------------
IGRID=NMGRID
NMGRID=1
CALL TABCOL_FORDIACHRO
DO J=1,NBFILES
  IF(NUMFILES(J) == NUMFILECUR)THEN
    JM=J
  ENDIF
ENDDO
!
IF(LXYZ00)THEN
  YSTO_CGROUPSV3=CGROUPSV3    
  ILENTRIMSV3=LEN(TRIM(CGROUPSV3))
  YNUMBER=CGROUPSV3(ILENTRIMSV3-1:ILENTRIMSV3)
  ! on verifie que CGROUPSV3 contient une variable lagrangienne
  ! pertinente sinon on remet Z000 pour cette routine puis on remet
  ! CROUPSV3 a ce qu il etait avant de rentrer dans cette routine
  IF (ICHAR(YNUMBER(1:1))<48 .OR. ICHAR(YNUMBER(1:1))>57 .OR. &
      ICHAR(YNUMBER(2:2))<48 .OR. ICHAR(YNUMBER(2:2))>57       ) THEN
    CGROUPSV3='Z000'
    PRINT*,'**TRAMASK3D: CGROUPSV3 force a Z000'
    PRINT*,'son ancienne valeur ',YSTO_CGROUPSV3, &
           ' sera remise a la sortie de tramask3d'
    ILENTRIMSV3=LEN(TRIM(CGROUPSV3))
    YNUMBER=CGROUPSV3(ILENTRIMSV3-1:ILENTRIMSV3)
  ENDIF
ENDIF
!
!
! Lecture des X0 -> chargement dans ZSVM1
!
IF(LXYZ00)THEN
  YGROUP='X0'//YNUMBER
ELSE
  YGROUP='LGXM'
ENDIF
CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
IF(LPBREAD)THEN
  IF(LXYZ00)THEN
    print *,' Absence de variable X00 .. Operation impossible'
    RETURN
  ELSE
    YGROUP='LGXT'
    LPBREAD=.FALSE.
    CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
    IF(LPBREAD)THEN
      YGROUP='SVM001'
      LPBREAD=.FALSE.
      CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
      IF(LPBREAD)THEN
        YGROUP='SVT001'
        LPBREAD=.FALSE.
        CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
        IF(LPBREAD)THEN
          YGROUP='SVM1'
          LPBREAD=.FALSE.
          CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
          IF(LPBREAD)THEN
            YGROUP='SVT1'
            LPBREAD=.FALSE.
            CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
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
ENDIF
IF(LGROUP)THEN
  IF(LMASK3D)THEN
  print *,' **TRAMASK3D utilisation de ',YGROUP
  ENDIF
  CALL READ_DIACHRO(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
ENDIF
! Chargement clegend clegend2
CALL RESOLV_TIMES(1)
YTEM2=' '
YTEM2=CLEGEND2
YTEM1=' '
! Elimination volontaire de 104 a 108 charge ds resolv_times pour RS
YTEM1=CLEGEND(1:103)
!
IF(.NOT.LFIC1)THEN
  CALL REALLOC_AND_LOAD(YGROUP)
  IF(LPBREAD)THEN
    print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
    ' L''UN DES FICHIERS '
    IF(ALLOCATED(XVAR))THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    ENDIF
    RETURN
  ENDIF
ENDIF
ALLOCATE(ZSVM1(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4)))
!IF(YGROUP == 'SVM1')THEN
IF(MAXVAL(XXHAT)/MAXVAL(XVAR) > 1.E2)THEN
  print *,' ** Tramask3D MAXVAL(XXHAT),MAXVAL(XVAR),*1000(KM->M) ',MAXVAL(XXHAT),MAXVAL(XVAR)
  WHERE(XVAR(:,:,:,:,1,1) /= XSPVAL)
    ZSVM1(:,:,:,:)=XVAR(:,:,:,:,1,1)*1000.
  ELSEWHERE
    ZSVM1(:,:,:,:)=XVAR(:,:,:,:,1,1)
  ENDWHERE
ELSE
  ZSVM1(:,:,:,:)=XVAR(:,:,:,:,1,1)
ENDIF
CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
!
! Lecture des Y0 -> chargement dans ZSVM2
!
IF(LXYZ00)THEN
  YGROUP='Y0'//YNUMBER
ELSE
  YGROUP='LGYM'
ENDIF
CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
IF(LPBREAD)THEN
  IF(LXYZ00)THEN
    print *,' Absence de variable Y00 .. Operation impossible'
    RETURN
  ELSE
    YGROUP='LGYT'
    LPBREAD=.FALSE.
    CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
    IF(LPBREAD)THEN
      YGROUP='SVM002'
      LPBREAD=.FALSE.
      CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
      IF(LPBREAD)THEN
        YGROUP='SVT002'
        LPBREAD=.FALSE.
        CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
        IF(LPBREAD)THEN
          YGROUP='SVM2'
          LPBREAD=.FALSE.
          CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
          IF(LPBREAD)THEN
            YGROUP='SVT2'
            LPBREAD=.FALSE.
            CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
            !
            IF(LPBREAD)THEN
        print *,' Absence de variable LGYM ou SVM002 ou LGYT ou SVT002 .. Operation impossible'
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDIF
IF(LGROUP)THEN
  IF(LMASK3D)THEN
  print *,' **TRAMASK3D utilisation de ',YGROUP
  ENDIF
  CALL READ_DIACHRO(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
ENDIF
IF(.NOT.LFIC1)THEN
  CALL REALLOC_AND_LOAD(YGROUP)
  print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
  ' L''UN DES FICHIERS '
  IF(ALLOCATED(XVAR))THEN
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
  ENDIF
  RETURN
ENDIF
ALLOCATE(ZSVM2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4)))
!IF(YGROUP == 'SVM2')THEN
IF(MAXVAL(XYHAT)/MAXVAL(XVAR) > 1.E2)THEN
  print *,' ** Tramask3D MAXVAL(XYHAT),MAXVAL(XVAR),*1000(KM->M) ',MAXVAL(XYHAT),MAXVAL(XVAR)
  WHERE(XVAR(:,:,:,:,1,1) /= XSPVAL)
    ZSVM2(:,:,:,:)=XVAR(:,:,:,:,1,1)*1000.
  ELSEWHERE
    ZSVM2(:,:,:,:)=XVAR(:,:,:,:,1,1)
  ENDWHERE
ELSE
  ZSVM2(:,:,:,:)=XVAR(:,:,:,:,1,1)
ENDIF
CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
!
! Lecture des Z0 -> chargement dans ZSVM3
!
IF(LXYZ00)THEN
  YGROUP='Z0'//YNUMBER
ELSE
  YGROUP='LGZM'
ENDIF
CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
IF(LPBREAD)THEN
  IF(LXYZ00)THEN
    print *,' Absence de variable Z00 .. Operation impossible'
    RETURN
  ELSE
    YGROUP='LGZT'
    LPBREAD=.FALSE.
    CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
    IF(LPBREAD)THEN
      YGROUP='SVM003'
      LPBREAD=.FALSE.
      CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
      IF(LPBREAD)THEN
        YGROUP='SVT003'
        LPBREAD=.FALSE.
        CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
        IF(LPBREAD)THEN
          YGROUP='SVM3'
          LPBREAD=.FALSE.
          CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
          IF(LPBREAD)THEN
            YGROUP='SVT3'
            LPBREAD=.FALSE.
            CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
            !
            IF(LPBREAD)THEN
        print *,' Absence de variable LGZM ou SVM003 ou LGZT ou SVT003 .. Operation impossible'
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDIF
IF(LGROUP)THEN
  IF(LMASK3D)THEN
  print *,' **TRAMASK3D utilisation de ',YGROUP
  ENDIF
  CALL READ_DIACHRO(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
ENDIF
IF(.NOT.LFIC1)THEN
  CALL REALLOC_AND_LOAD(YGROUP)
  IF(LPBREAD)THEN
    print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
    ' L''UN DES FICHIERS '
    IF(ALLOCATED(XVAR))THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    ENDIF
    RETURN
  ENDIF
ENDIF
ALLOCATE(ZSVM3(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4)))
!IF(YGROUP == 'SVM3')THEN
IF(MAXVAL(XZHAT)/MAXVAL(XVAR) > 1.E2)THEN
  print *,' ** Tramask3D MAXVAL(XZHAT),MAXVAL(XVAR),*1000(KM->M) ',MAXVAL(XZHAT),MAXVAL(XVAR)
  WHERE(XVAR(:,:,:,:,1,1) /= XSPVAL)
    ZSVM3(:,:,:,:)=XVAR(:,:,:,:,1,1)*1000.
  ELSEWHERE
    ZSVM3(:,:,:,:)=XVAR(:,:,:,:,1,1)
  ENDWHERE
ELSE
  ZSVM3(:,:,:,:)=XVAR(:,:,:,:,1,1)
ENDIF
CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
!
! Lecture du champ lagrangien suppl -> chargement dans ZFIELD_LAG
!
GLAG=LXYZ00 .AND.                                         &
    (CGROUPSV3(1:2).NE.'SV' .AND. CGROUPSV3(1:2).NE.'LG'  &
                            .AND. CGROUPSV3(1:2).NE.'Z0'  &
    )
!
IF( GLAG )THEN
  YGROUP=CGROUPSV3
  CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
  IF(LPBREAD)THEN
    print *,' Absence de variable ',CGROUPSV3,' .. Operation impossible'
    RETURN
  ENDIF
  IF(LGROUP)THEN
    IF(LMASK3D)THEN
    print *,' **TRAMASK3D utilisation suppl. de ',YGROUP
    ENDIF
    CALL READ_DIACHRO(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
  ENDIF
  IF(.NOT.LFIC1)THEN
    CALL REALLOC_AND_LOAD(YGROUP)
    IF(LPBREAD)THEN
      print *,' REQUETE IMPOSSIBLE .',YGROUP,' N''EXISTE PAS DANS', &
      ' L''UN DES FICHIERS '
      IF(ALLOCATED(XVAR))THEN
        CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
      ENDIF
      RETURN
    ENDIF
  ENDIF
  ALLOCATE(ZFIELD_LAG(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4)))
  ZFIELD_LAG(:,:,:,:)=XVAR(:,:,:,:,1,1)
  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
ENDIF
!
! Determination du masque en fonction de la fenetre XXL,XXH,XYL,XYH,XZL,XZH
!
IF(ALLOCATED(LMASK3))THEN
  DEALLOCATE(LMASK3)
ENDIF
ALLOCATE(LMASK3(SIZE(ZSVM1,1),SIZE(ZSVM1,2),SIZE(ZSVM1,3),SIZE(ZSVM1,4)))
LMASK3=.FALSE.
!
IF (GLAG) THEN
 LMASK3=(XXL < ZSVM1 .AND. XXH >ZSVM1) .AND. (XYL < ZSVM2 .AND. XYH > ZSVM2) &
 .AND. (XZL < ZFIELD_LAG .AND. XZH > ZFIELD_LAG)
ELSE
 LMASK3=(XXL < ZSVM1 .AND. XXH >ZSVM1) .AND. (XYL < ZSVM2 .AND. XYH > ZSVM2) &
 .AND. (XZL < ZSVM3 .AND. XZH > ZSVM3)
ENDIF      
!
! Calcul des altitudes pour la grille 1 dans XZZ
!
CALL COMPCOORD_FORDIACHRO(1)
!
!-------------------------------------------------------------------------------
!
! Visualisation du masque sur XY, XZ, YZ
!
IF(LMASK3D .OR. LMASK3D_XY .OR. LMASK3D_XZ .OR. LMASK3D_YZ)THEN
!IF(LMASK3D)THEN
!
  IF(LPRINT)THEN
    CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
    IF(IRESP /= 0)THEN
      CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
      OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
      PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
    ENDIF
  ENDIF

! Recuperation de la fenetre d'affichage courante pour restauration en fin de
! routine
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! 
! Determination de NIINF NJINF NISUP NJSUP si non initialises par l'utilisateur
IF(NIINF == 0 .AND. NISUP == 0 .AND. NJINF == 0 .AND. NJSUP == 0)THEN
  CALL RESOLV_NIJINF_NIJSUP
ENDIF

IIB=1+JPHEXT; IIE=SIZE(ZSVM1,1)-JPHEXT
IJB=1+JPHEXT; IJE=SIZE(ZSVM1,2)-JPHEXT
IKB=1+JPVEXT; IKE=SIZE(ZSVM1,3)-JPVEXT

DO JTLOOP=1,SIZE(ZSVM1,4)
if(nverbia >0)then
print *,' ** TRAMASK3D JTLOOP ',JTLOOP
endif
!
!!!!!! XY 
!
IF(LMASK3D_XY)THEN

IF(NJMAX /= 1)THEN
  IF(LPRINT)THEN
   ILOOP=SIZE(ZSVM1,1)/5
   IF(ILOOP * 5 < SIZE(ZSVM1,1)) ILOOP=ILOOP+1
   WRITE(INUM,'(''CH  '',''G:'',A16,'' P:'',A25,'' T:'',F8.0,''s'')') &
        CGROUPSV3,CTITRE(1)(1:25),XTRAJT(JTLOOP,1)
   WRITE(INUM,'(A40,''(NIINF-NISUP,NJINF-NJSUP)'')')CTITGAL
   WRITE(INUM,'(''niinf'',i4,'' njinf'',i4,'' nisup'',i4,'' njsup'',i4,&
         &''   '',A1,'' '',i6)')&
         &NIINF,NJINF,NISUP,NJSUP,CTYPHOR,IKE
   WRITE(INUM,'(''NBVAL en I '',i4,''  NBVAL en J '',i4,''   iter'',i3)') &
        &NISUP-NIINF+1,NJSUP-NJINF+1,ILOOP
  ENDIF

YCAR(1:LEN_TRIM(YCAR))=' '
WRITE(YCAR,'(''MASK **XY-  ** window:('',F8.0,'':'',F8.0,'','',F8.0,'':'',F8.0,'','',F6.0,'':'',F6.0,'')'')') &
            XXL,XXH,XYL,XYH,XZL,XZH
YCAR(11:12)=YGROUP(3:4)
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!CALL PCSETC('FC','/')
!CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!CALL PCSETC('FC',':')
CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

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
! trace du masque (etoiles colorees)
CALL GSMK(3)
DO JKLOOP=IKE,IKB,-1
DO JILOOP=IIB,IIE
 DO JJLOOP=IJB,IJE
   IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP))THEN
     ZX=XXX(JILOOP,1)
     ZY=XXY(JJLOOP,1)
     CALL GPM(1,ZX,ZY)
   ENDIF
 ENDDO
ENDDO
ENDDO
! trace de la boite de lacher
CALL GSPLCI(4)
CALL GSLWSC(3.)
CALL FRSTPT(XXL,XYL)
CALL VECTOR(XXH,XYL)
CALL VECTOR(XXH,XYH)
CALL VECTOR(XXL,XYH)
CALL VECTOR(XXL,XYL)
CALL FRAME
ENDIF

ENDIF
CALL GSLWSC(1.)
CALL GSPLCI(1)
!
!!!!!! XZ 
!
IF(LMASK3D_XZ)THEN
WRITE(YCAR,'(''MASK **XZ-  ** window:('',F8.0,'':'',F8.0,'','',F8.0,'':'',F8.0,'','',F6.0,'':'',F6.0,'')'')') &
            XXL,XXH,XYL,XYH,XZL,XZH
!
IF(GLAG) THEN
  YCAR(11:12)=YNUMBER
ELSE
  YCAR(11:12)=YGROUP(3:4)
ENDIF
!
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!CALL PCSETC('FC','/')
!CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!CALL PCSETC('FC',':')
CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

IF(LDATFILE)CALL DATFILE_FORDIACHRO

CALL SET(.1,.9,.1,.9,XXX(NIINF,1),XXX(NISUP,1), &
         XHMIN,XHMAX,1)
YFORMAX='          '
IF(LFMTAXEX)THEN
  YFORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
ELSE
  YFORMAX='(F8.0)'
ENDIF
YFORMAY='          '
IF(LFMTAXEY)THEN
  YFORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
ELSE
  YFORMAY='(F6.0)'
ENDIF
!
CALL LABMOD(YFORMAX,YFORMAY,0,0,10,10,0,0,0)
!CALL LABMOD('(F8.0)','(F6.0)',9,6,10,10,0,0,0)
CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NINT((XHMAX-XHMIN)/1000.),1,1,1,5,0.,0.)
!
IF (GLAG) THEN
  ! trace du masque (etoiles colorees)
  CALL GSMK(3)
  CALL GSPMCI(1)
  DO JILOOP=IIB,IIE
  DO JJLOOP=IJB,IJE
  DO JKLOOP=IKB,IKE
    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP) )THEN
      ZX=XXX(JILOOP,1)
      ZY=XZZ(JILOOP,JJLOOP,JKLOOP)
      CALL GPM(1,ZX,ZY)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  !
  ! trace de la zone de lacher (cercles)
  CALL GSMK(4)
  CALL GSPMCI(3)
  DO JILOOP=IIB,IIE
  DO JJLOOP=IJB,IJE
  DO JKLOOP=IKB,IKE
    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP) )THEN
      ZX=ZSVM1(JILOOP,JJLOOP,JKLOOP,JTLOOP)
      ZY=ZSVM3(JILOOP,JJLOOP,JKLOOP,JTLOOP)
      CALL GPM(1,ZX,ZY)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! 
ELSE
  ! trace du masque (etoiles colorees)
  DO JILOOP=IIB,IIE
  DO JJLOOP=IJB,IJE
  DO JKLOOP=IKB,IKE
    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP) )THEN
      ZX=XXX(JILOOP,1)
      ZY=XZZ(JILOOP,JJLOOP,JKLOOP)
      CALL GPM(1,ZX,ZY)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! trace de la boite de lacher
  CALL GSPLCI(3)
  CALL GSLWSC(3.)
  CALL FRSTPT(XXL,XZL)
  CALL VECTOR(XXH,XZL)
  CALL VECTOR(XXH,XZH)
  CALL VECTOR(XXL,XZH)
  CALL VECTOR(XXL,XZL)
ENDIF
!
CALL FRAME
ENDIF
CALL GSLWSC(1.)
CALL GSPLCI(1)
!
!!!!!! YZ 
!
IF(LMASK3D_YZ)THEN

IF(NJMAX /= 1)THEN
WRITE(YCAR,'(''MASK **YZ-  ** window:('',F8.0,'':'',F8.0,'','',F8.0,'':'',F8.0,'','',F6.0,'':'',F6.0,'')'')') &
      XXL,XXH,XYL,XYH,XZL,XZH
IF(GLAG) THEN
  YCAR(11:12)=YNUMBER
ELSE
  YCAR(11:12)=YGROUP(3:4)
ENDIF
!
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!CALL PCSETC('FC','/')
!CALL PLCHHQ(.002,.98,YCAR(1:LEN_TRIM(YCAR)),.012,0.,-1.)
!CALL PCSETC('FC',':')
CALL TIT_TRA3D(YCAR,YTEM1,YTEM2,ZVR)

IF(LDATFILE)CALL DATFILE_FORDIACHRO

CALL SET(.1,.9,.1,.9,XXY(NJINF,1),XXY(NJSUP,1), &
         XHMIN,XHMAX,1)
YFORMAX='          '
IF(LFMTAXEX)THEN
  YFORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
ELSE
  YFORMAX='(F8.0)'
ENDIF
YFORMAY='          '
IF(LFMTAXEY)THEN
  YFORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
ELSE
  YFORMAY='(F6.0)'
ENDIF

CALL LABMOD(YFORMAX,YFORMAY,0,0,10,10,0,0,0)
!CALL LABMOD('(F8.0)','(F6.0)',9,6,10,10,0,0,0)
CALL GRIDAL(NCVITVXMJ,NCVITVXMN,NINT((XHMAX-XHMIN)/1000.),1,1,1,5,0.,0.)
IF (GLAG) THEN
  ! trace du masque (etoiles colorees)
  CALL GSMK(3)
  CALL GSPMCI(1)
  DO JILOOP=IIB,IIE
  DO JJLOOP=IJB,IJE
  DO JKLOOP=IKB,IKE
    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP) )THEN
      ZX=XXY(JJLOOP,1)
      ZY=XZZ(JILOOP,JJLOOP,JKLOOP)
      CALL GPM(1,ZX,ZY)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  !
  ! trace de la zone de lacher (cercles)
  CALL GSMK(4)
  CALL GSPMCI(2)
  DO JILOOP=IIB,IIE
  DO JJLOOP=IJB,IJE
  DO JKLOOP=IKB,IKE
    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP) )THEN
      ZX=ZSVM2(JILOOP,JJLOOP,JKLOOP,JTLOOP)
      ZY=ZSVM3(JILOOP,JJLOOP,JKLOOP,JTLOOP)
      CALL GPM(1,ZX,ZY)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  !
ELSE
  ! trace du masque (etoiles colorees)
  DO JILOOP=IIB,IIE
  DO JJLOOP=IJB,IJE
  DO JKLOOP=IKB,IKE
    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,JTLOOP) )THEN
      ZX=XXY(JJLOOP,1)
      ZY=XZZ(JILOOP,JJLOOP,JKLOOP)
      CALL GPM(1,ZX,ZY)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! trace de la boite de lacher
  CALL GSPLCI(2)
  CALL GSLWSC(3.)
  CALL FRSTPT(XYL,XZL)
  CALL VECTOR(XYH,XZL)
  CALL VECTOR(XYH,XZH)
  CALL VECTOR(XYL,XZH)
  CALL VECTOR(XYL,XZL)
ENDIF
!
CALL FRAME
ENDIF
ENDIF

ENDDO
! Recuperation du viewport courant pour son eventuelle impression
CALL GETSET(XCURVPTL,XCURVPTR,XCURVPTB,XCURVPTT,ZWLBID,ZWRBID,ZWBBID,ZWTBID,IDBID)
! Restauration de la fenetre a l'entree de la routine
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
ENDIF

CALL GSPLCI(1)
CALL GSPMCI(1)
CALL GSLWSC(1.)
CALL GSLN(1)
DEALLOCATE(ZSVM1,ZSVM2,ZSVM3)
IF (GLAG) DEALLOCATE(ZFIELD_LAG)
IF (LXYZ00) CGROUPSV3=YSTO_CGROUPSV3
NMGRID=IGRID
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
!
!
RETURN
END SUBROUTINE TRAMASK3D
