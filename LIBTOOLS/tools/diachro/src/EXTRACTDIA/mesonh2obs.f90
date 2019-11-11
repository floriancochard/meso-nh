      PROGRAM  MESONH2OBS
!     ###################
!
!!****  *MESONH2OBS* -  Interpolation d un champ Mesonh sur les points
!                       d'observation donnees en entree
!                       en sortie un fichier ascii lon lat alt valeur_interpolee_modele
!! 
!!
!!    PURPOSE
!!    -------
! 
!     Lecture en entree:
!       d'un fichier ascii contenant les localisations (lon,lat ou lat,lon)
!       d'un fichier diachronique a traiter (boucle sur les fichiers)
!       du champ modèle a interpoler (boucle sur les champs)
!      
!     Ecriture en sortie:
!       d'un fichier au format 
!         lon lat alt new_val_modele   avec alt=altitude d un niveau de modele
!      ou lat lon alt new_val_modele    ou  alt=Z constante
!                                       ou  alt=P constante
!
!!**  METHOD
!!    ------
!      3.2.a Creation de la grille des obs en X, Y et Z
!           lecture du fichier de localisations
!         (LLHV:lon-lat, LLZV:lon-lat-alt(metres), LLPV:lon-lat-pres(hPa))
!         (llhv:lat-lon, llzv:lat-lon-alt(metres), llpv:lat-lon-pres(hPa)
!           calcul des X et Y correspondants  
!      3.2.b Interpolation verticale du champ (3D) MesoNH (LLZV ou LLPV)
!      3.2.c Interpolation horizontale sur la grille des obs
!          c1     "            "       du champ MesoNH
!          c2     "            "       du tableau de grille vert. 
!                                   (champ 2D ou champ 3D en LLHV)
!      3.4. Ecriture par writellhv
!      
!     
!!      
!!    EXTERNAL
!!    --------
!!          CREATLINK : a l'ouverture du fichier, HFLAGFILE='OPE',
!!                      creation d'un lien dans le directory local
!!                      si le fichier existe sous $DIROBS
!!          DD et FF  : calcul de dd et ff a partir des composantes U et V
!!          READVAR   : lecture d unchamp du fichier diachronique
!!          WRITELLHV : ecriture format lon lat alt val
!!          SYSTEM    : renommer le fichier de sortie avec un nom > 28 carateres
!!          zinter    : interpolation verticale en Z=cst
!!          pinter    : interpolation verticale en P=cst
!!          hor_interp_4pts : interpolation horizontale
!!          SM_XYHAT  : creation de la grille des Obs
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHORS
!!    -------
!!    N. Asencio
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/09/2003
!     09/10/2003 use XXX(:,NGRID) et XXY(:,NGRID) for hor_interp4pts 
!                and SM_XYHAT calls
!     04/05/2005 add a control for the min and max of the field before
!                and after interpolation(s)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
#ifdef NAGf95
USE F90_UNIX  ! for FLUSH
USE F90_UNIX_PROC  ! for SYSTEM
#endif      
! modules MesoNH
USE MODD_CST
USE MODD_PARAMETERS, ONLY:JPHEXT,JPVEXT,XUNDEF
USE MODD_DIM1, ONLY:NIMAX,NJMAX,NKMAX, NIINF, NISUP ,NJINF,NJSUP
USE MODD_GRID,  ONLY: XLATORI,XLONORI
USE MODD_GRID1, ONLY: XXHAT,XYHAT,XZZ
USE MODE_GRIDPROJ  ! subroutine SM_XYHAT 
! modules DIACHRO
USE MODN_NCAR,  ONLY: XSPVAL      
USE MODD_COORD
!                    XVAR(i,j,k,,,), XMASK,XTRAJ ,XDATIME(16,t)
!                     et NGRIDIA , NGRIDIAM ( appel interp_grids)
USE MODD_ALLOC_FORDIACHRO
!                      nverbia, CGROUP
USE MODD_RESOLVCAR 
!
! modules tools
USE MODI_HOR_INTERP_4PTS 
USE MODI_ZINTER      
USE MODI_PINTER      
USE MODI_WRITELLHV
USE MODI_DD
USE MODI_FF
USE MODI_UV_TO_ZONAL_AND_MERID
USE MODI_CREATLINK
USE MODI_LOW2UP
USE MODI_WRITEDIR
! modules extractdia
USE MODD_READLH                 ! domaine initialise par READVAR: 
                                !NREADIL,NREADIH, NREADJL,NREADJH,
                                !NREADKL,NREADKH
!
IMPLICIT NONE                       
!    
!*       0.1   Local variables declarations
!
! indices de boucle
INTEGER           :: JILOOP,JLOOPFILE,JGR,JNobsLOOPsite,JNobsLOOPz,JNobsLOOPtriplet
! zoom  suivant les 6 dimensions des champs diachro
INTEGER           :: ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin
INTEGER           :: ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup
REAL , allocatable, dimension(:,:,:) :: ZVAR3D
REAL , allocatable, dimension(:,:)   :: ZDATIME
!  pour l interpolation horizontale  
REAL , allocatable, dimension(:,:) :: ZOBSLAT,ZOBSLON,ZOBSX,ZOBSY
REAL , allocatable, dimension(:) :: ZOBSALT
REAL , allocatable, dimension(:) :: ZOBSLATlu,ZOBSLONlu,ZOBSALTlu
REAL , allocatable, dimension(:,:,:) :: ZVARNEWH
! pour l interpolation verticale: zinter
REAL , allocatable, dimension(:,:,:) :: ZVARZCST
INTEGER :: ikdebzint ! premier niveau a traiter      
! pour l appel a WRITELLHV
REAL , allocatable, dimension(:,:,:) :: ZALT
! calcul dd ff
REAL , allocatable, dimension(:,:,:) :: zwork3d,zwork3d2,zffvent,zdirvent
! pour pinter : interpolation a P=cst
REAL , allocatable, dimension(:,:,:) :: zpabs      
INTEGER :: iskip,IGRID,IGRIDOUT,ILUDIR
REAL :: zmini,zmaxi
!                                   
INTEGER :: iret,iret2,ilocverbia,inbvalxy,idimlonlat,inbvalz,inbvalxyz, &
           inbvalz3d,inbvalxyz3d
!! **** la taille des variables caracteres contenant les noms
!!      de fichiers est obligatoirement de 28 ****
!!      pour toutes les routines diachro
CHARACTER(LEN=28) :: YFILEIN,YFILEIN2,YDUMMYFILE
CHARACTER(LEN=100):: YFILEOBS, YSAVEFILEOBS
! **** la longueur du nom ne doit pas depasser 13 car. si le fichier
! contient des groupes a un seul PROCessus, ou 9 si plusieurs PROCessus ****
CHARACTER(LEN=13) :: YGROUP
CHARACTER(LEN=4)  :: YTYPEOUT
CHARACTER(LEN=5)  :: YFLAGREADVAR,YFLAGWRITE
CHARACTER(LEN=9)  :: ygrillevert     ! type de grille verticale selon
                                     ! champ2D/3D et YTYPEOUT
CHARACTER(LEN=36) :: YFILEOUT
CHARACTER(LEN=3)  :: YREP
CHARACTER(LEN=100):: ycommand                                     
CHARACTER(LEN=11) :: YLUDIR      !  Name of the dir file
!-------------------------------------------------------------------------------
!
!*       1.     Init
!               ----
!
YFILEIN2=''
! active(1) ou desactive(0) les prints de controle dans les routines
! READVAR et WRITE
ilocverbia=3
! 
! dans mesonh Xundef est utilise 
! dans les routines diachro XSPVAL est utilise
 XSPVAL=XUNDEF                                    
!
! ouverture d un fichier dir ou vont s ecrire les entrees clavier
YLUDIR='dirmnh2obs'
CALL FMATTR(YLUDIR,YLUDIR,ILUDIR,iret)
OPEN(UNIT=ILUDIR,FILE=YLUDIR,FORM='FORMATTED')
!
!                                       
PRINT*, '- Format of the output file: '
PRINT*, ' (and of the input observation file with positions)'
PRINT*, ' Lon-Lat-Height(MNH)-Value= LLHV'
PRINT*, ' lat-lon-height(MNH)-value= llhv'
PRINT*, ' Lon-Lat-Z(m)-Value       = LLZV'
PRINT*, ' lat-lon-Z(m)-value       = llzv'
PRINT*, ' Lon-Lat-P(hPa)-Value     = LLPV'
PRINT*, ' lat-lon-P(hPa)-value     = llpv'
PRINT*, '?'
READ(5,'(A)')YTYPEOUT
CALL WRITEDIR(ILUDIR,YTYPEOUT)
!
SELECT CASE (YTYPEOUT(1:2))  ! type de coordonnées lon,lat ou lat,lon
  CASE('LL')
    PRINT*,'-> positions in the observation file are given in  lon lat'
  CASE('ll')
    PRINT*,'-> positions in the observation file are given in  lat lon'
END SELECT
    inbvalz=1
SELECT CASE (YTYPEOUT(3:3))
  CASE('Z','z','P','p')                                   
    PRINT*,'- Are the vertical levels included in the input observation file ?'
    PRINT*,'   Y= format of the obs file=coord1 coord2 level'
    PRINT*,'   N= format of the obs file=coord1 coord2 '
    PRINT*,'      and levels provided interactively '
    READ(5,'(A)') YREP
    CALL WRITEDIR(ILUDIR,YREP)
    YREP=ADJUSTL(YREP)
    SELECT CASE (YREP(1:1))
      CASE('O','o','Y','y') 
        inbvalz=0
      CASE DEFAULT 
        PRINT*, '- Number of vertical levels for the interpolation ', YTYPEOUT(3:3),' ?'
        READ(5,*) inbvalz
        CALL WRITEDIR(ILUDIR,inbvalz)
        PRINT*, '- List of these levels (in meters or in hPa): exemple 500 1500 ?'
        allocate (ZOBSALTlu(inbvalz))
        READ(5,*) ZOBSALTlu
        DO JILOOP=1,inbvalz
          CALL WRITEDIR(ILUDIR,ZOBSALTlu(JILOOP))
        END DO
        PRINT*, ' interpolation for the following levels ',YTYPEOUT(3:3),'='
        PRINT*, ZOBSALTlu
    END SELECT
  CASE('H','h')
    PRINT*,'-> the vertical levels will be the same as in the model'
END SELECT
!
PRINT*, '- Name of the file which contains the localisation of the obs ?'
READ(5,'(A)',END=99) YFILEOBS
CALL WRITEDIR(ILUDIR,YFILEOBS)
!
PRINT*, '- Prints : 0= mini 1=mode debug in mesonh2obs'
PRINT*, '                   3=debug mode in dichro routines'
PRINT*, '?'
READ(5,*)ilocverbia
CALL WRITEDIR(ILUDIR,ilocverbia)
PRINT*, ' output prints= ',ilocverbia 
IF (ilocverbia >2) nverbia=ilocverbia
!
!-------------------------------------------------------------------------------
!
!*       2.     Boucle sur les fichiers a traiter
!               ---------------------
DO JLOOPFILE=1,100000
  !
  !*       2.1    Lecture Nom de fichier et type de sortie
  !              ----------------------
  PRINT*, '- Name of the diachro file (without .lfi) (END to stop) ?'
  IF (LEN_TRIM(YFILEIN2)/=0) PRINT*, ' other than ',TRIM(YFILEIN2)
  READ(5,'(A28)',END=99) YFILEIN
  CALL WRITEDIR(ILUDIR,YFILEIN)
  IF ( YFILEIN(1:3) == 'END' .OR. YFILEIN(1:3) == 'end' ) GO TO 99
  !
  !  indique que le fichier d entree lu doit etre ouvert dans READVAR
  YFLAGREADVAR='OPE'
  !  indique que le fichier de sortie doit etre ouvert dans WRITELLHV
  !  et que l entete sera ecrite uniquement lors de la premiere ecriture
  YFLAGWRITE='NEW1H'
  !
  IF (YTYPEOUT(1:4)=='LLPV' .OR. YTYPEOUT(1:4)=='llpv') THEN
    CALL READVAR('PABST',YFILEIN,YFLAGREADVAR,ilocverbia,iret)
    IF ( iret /= 0 ) then
      print *, '- PABST not found, name of the pressure variable ?'
      read *,YGROUP
      CALL WRITEDIR(ILUDIR,YGROUP)
      CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
      IF ( iret /= 0 ) then
        print *,' interpolation at P=cst not possible because PABST and ',TRIM(YGROUP), ' are not available'
        STOP
      ENDIF
    ENDIF
    ! stockage de ZPABS utilise par pinter
    ALLOCATE ( ZPABS(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3))) 
    ZPABS(:,:,:)=XVAR(:,:,:,1,1,1)
  END IF
  PRINT*, 'Input file: ', TRIM(YFILEIN), ', type of output is: ',YTYPEOUT
  !
  NIINF=0
  NISUP=0
  NJINF=0
  NJSUP=0
  !
  !*       3.     Boucle sur les champs a traiter dans le fichier
  !               ----------------------
  !
  DO JGR=1,10000
    !
    PRINT*, '-  Name of the group in upper case (13 characters max.)'
    PRINT*, ' (ex: THT ou DD ou FF ou DD10 ou FF10 )'
    PRINT*, '(GROUP for the list of groups, END to stop)?'
    READ(5,'(A13)',END=88) CGROUP
    CALL WRITEDIR(ILUDIR,CGROUP)
    CGROUP=ADJUSTL(CGROUP)
    CALL LOW2UP(CGROUP)
    IF (CGROUP=='END') GO TO 88
    IF (CGROUP(1:5)/='GROUP') &
      PRINT*,'you asked for the following record: ',TRIM(CGROUP)
    !
    IGRIDOUT=-1
    !
    !*     3.1  Lecture et initialisation de XVAR (MODD_ALLOC_FORDIACHRO)
    !
    SELECT CASE (CGROUP(1:2))                                   
      !
      CASE('DD','FF','UM','VM','UT','VT')
        !
        ! Lecture du champ UM et VM apres traitement de UM (voir en 3.2)
        IF (LEN(TRIM(CGROUP)) ==2) THEN
          YGROUP='UT'
        ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
          YGROUP='UT'//CGROUP(3:4)
        ELSE
          print*,'** problem with the name of group: ',CGROUP
          CYCLE
        ENDIF
        CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
        IF ( iret /= 0 ) then
          print *,TRIM(CGROUP),': ',TRIM(YGROUP),' not available'
          IF (LEN(TRIM(CGROUP)) ==2) THEN
            YGROUP='UM'
          ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
            YGROUP='UM'//CGROUP(3:4)
          ENDIF
          CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret2)
          IF ( iret2 /= 0 ) then
            print *,'** no processing for ',TRIM(CGROUP), &
                    ' because UT and ',TRIM(YGROUP),' not available'
            CYCLE
          ENDIF
          iret=iret2
        ENDIF
        IGRIDOUT=1  ! le champ DD,FF,UZON ou VMED sera en grille masse
        !
        !  3.1.1  traitement sup  du tableau XVAR si DD ou FF ou UM ou VM
        !
        ! Allocation des tableaux de stockage de la premiere composante
        ALLOCATE(zwork3d(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
        zwork3d(:,:,:)=XVAR(:,:,:,1,1,1)
        IF (LEN(TRIM(CGROUP)) ==2) THEN
          YGROUP='VT'
        ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
          YGROUP='VT'//CGROUP(3:4)
        ENDIF
        CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
        if ( iret /= 0 ) then
          print *,TRIM(CGROUP),': ',TRIM(YGROUP),' not available'
          IF (LEN(TRIM(CGROUP)) ==2) THEN
            YGROUP='VM'
          ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
            YGROUP='VM'//CGROUP(3:4)
          ENDIF
          CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret2)
          IF ( iret2 /= 0 ) then
            print *,'** traitement of ',TRIM(CGROUP), &
                    ' not possible because VT and ',TRIM(YGROUP), &
                    ' are not available'
            CYCLE
          ENDIF
          iret=iret2
!          CYCLE
        endif
        ! Allocation des tableaux de calcul
        ALLOCATE(zffvent(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
        ALLOCATE(zdirvent(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
        zffvent=XSPVAL
        zdirvent=XSPVAL
        !
        !   Calcul de dd ff
        !
        IF (CGROUP(1:2) == 'FF' .OR. CGROUP(1:2) == 'DD') THEN
          ! call ff (zwork3d,zwork3d2,zffvent,kvext,khext,kgrid)
          !!CALL FF(zwork3d(:,:,:),XVAR(:,:,:,1,1,1),zffvent,0,0,3)
          IGRID=NGRIDIA(SIZE(XVAR,6))
          print *,'avant ff dd:JPVEXT,JPHEXT,IGRID', JPVEXT,JPHEXT,IGRID
          CALL FF(zwork3d(:,:,:),XVAR(:,:,:,1,1,1),zffvent,JPVEXT,JPHEXT,IGRID)
          ! tous les points de grille: iskip=1
          iskip=1
          ! call dd(zwork3d,zwork3d2,zdirvent,iskip,kgrid,PLON=ZOBSLON)
          CALL DD(zwork3d(:,:,:),XVAR(:,:,:,1,1,1),zdirvent,iskip,3)
          print *,' End of computation of dd and ff'
          !
          ! Stockage dans le tableau XVAR qui est le tableau ecrit
          !
          IF (CGROUP(1:2) == 'FF' ) THEN
            XVAR(:,:,:,1,1,1)=zffvent
            IF (LEN(TRIM(CGROUP)) ==2) THEN
              YGROUP='VENTFF'
            ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
              YGROUP='VENT'//CGROUP(3:4)//'FF'
            ENDIF
          ELSE IF (CGROUP(1:2) == 'DD') THEN
            XVAR(:,:,:,1,1,1)=zdirvent
            IF (LEN(TRIM(CGROUP)) ==2) THEN
              YGROUP='VENTDD'
            ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
              YGROUP='VENT'//CGROUP(3:4)//'DD'
            ENDIF
               CUNITE(1)='degrees'
          ENDIF
          CGROUP=YGROUP
          CTITRE(1)=YGROUP
          NGRIDIA(1)=IGRIDOUT ! dd et ff en grille de masse
        !
        !   Calcul des composantes zonale et meridienne
        !
        ELSE IF (CGROUP(1:2) == 'UM' .OR. CGROUP(1:2) == 'VM' .OR. &
                 CGROUP(1:2) == 'UT' .OR. CGROUP(1:2) == 'VT'      ) THEN
          CALL UV_TO_ZONAL_AND_MERID(zwork3d(:,:,:),  &
                                     XVAR(:,:,:,1,1,1), &
                                     23,PZC=zffvent,PMC=zdirvent)
          IF (CGROUP(1:1) == 'U' ) THEN
            XVAR(:,:,:,1,1,1)=zffvent(:,:,:)
            IF (LEN(TRIM(CGROUP)) ==2) THEN
              YGROUP='UZON'
            ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
              YGROUP='U'//CGROUP(3:4)//'ZON'
            ENDIF
            CTITRE(1)='U zonal wind component'
          ELSE IF (CGROUP(1:1) == 'V' ) THEN
            XVAR(:,:,:,1,1,1)=zdirvent(:,:,:)
            IF (LEN(TRIM(CGROUP)) ==2) THEN
              YGROUP='VMED'
            ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
              YGROUP='V'//CGROUP(3:4)//'MED'
            END IF
            CTITRE(1)='V meridian wind component'
          ENDIF
          CGROUP=YGROUP
          NGRIDIA(1)=IGRIDOUT ! UZON et VMED en grille de masse
        ENDIF
        DEALLOCATE(zwork3d)
        DEALLOCATE(zffvent,zdirvent)
      !
      CASE default
        !
        ! Lecture du  champ CGROUP et stockage dans XVAR
        ! + Initialisation (si YFLAGREADVAR='OPE') des variables
        ! des modules (cf USE en debut de programme)
        !
        CALL READVAR(CGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
        IF (CGROUP(1:5)=='GROUP') CYCLE
    END SELECT    
    !
    IF ( iret == 0 ) THEN
      IF(SIZE(NGRIDIA,1)/=1) THEN
        print *,'** no processing for ',TRIM(CGROUP), &
                ' because several processus'
        CYCLE
      ENDIF
      IGRID=NGRIDIA(1)
      IF(IGRIDOUT==-1) IGRIDOUT=NGRIDIA(1)
      !
      zmini=MINVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
      zmaxi=MAXVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
      print * ,' After reading, min,max of the field ',TRIM(CGROUP),'=', zmini,zmaxi
      if (ilocverbia >= 0 ) then
        print *,' Size of array read= ',SIZE(XVAR,1),SIZE(XVAR,2),&
        SIZE(XVAR,3),SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)
        PRINT*, 'NIINF,NISUP,NJINF,NJSUP', NIINF,NISUP,NJINF,NJSUP
      endif
      ! 
      IF (NIMAX==1 .AND. NJMAX==1) THEN
        PRINT *,'** 1D group: use rather extractdia since the observation file will not be taken into account' 
        CYCLE
      ENDIF
      !
      !*       3.2  Traitement du tableau XVAR sur les points lat lon
      !             --------------------------------------------------
      ! 
      !        3.2.a Creation de la grille  des obs en X et Y et Z
      !
      IF ( .NOT. ALLOCATED ( ZOBSX) ) THEN
        !  creation de la grille  des obs en X et Y
        ! realisee apres la premiere lecture d un champ modele
        ! pour avoir les tableaux XXX XXY initialises 
        print *,' Creation of the X Y grid of the obs'
        ! 
        ! Lecture des sites lon lat
        YSAVEFILEOBS=YFILEOBS
        CALL CREATLINK('DIROBS',YFILEOBS,'CREAT',ilocverbia)
        OPEN(UNIT=8,FILE=TRIM(ADJUSTL(YFILEOBS)),STATUS='OLD',&
             FORM='FORMATTED')
        idimlonlat=SIZE(XVAR,1)*SIZE(XVAR,2)
        ALLOCATE ( ZOBSLATlu(idimlonlat) )
        ALLOCATE ( ZOBSLONlu(idimlonlat) )
        IF (.NOT.ALLOCATED (ZOBSALTlu)) THEN
           ALLOCATE ( ZOBSALTlu(idimlonlat) )
           ZOBSALTlu=0.
        ENDIF
        inbvalxy=0
        DO JILOOP=1,idimlonlat
          SELECT CASE ( YTYPEOUT(1:4))
            CASE ('LLZV','LLPV')
              IF ( inbvalz == 0) THEN
                READ (8,*,END=888) &
                    ZOBSLONlu(JILOOP),ZOBSLATlu(JILOOP),ZOBSALTlu(JILOOP)
              ELSE
                ! niveaux vert. deja lus en interactif
                READ (8,*,END=888) &
                    ZOBSLONlu(JILOOP),ZOBSLATlu(JILOOP)
              ENDIF
            CASE ('llzv','llpv')
              IF ( inbvalz == 0) THEN
                READ (8,*,END=888) &
                    ZOBSLATlu(JILOOP),ZOBSLONlu(JILOOP),ZOBSALTlu(JILOOP)
              ELSE
                ! niveaux vert. deja lus en interactif
                READ (8,*,END=888) &
                    ZOBSLATlu(JILOOP),ZOBSLONlu(JILOOP)
              ENDIF
            CASE ('LLHV')
              READ (8,*,END=888) &
                    ZOBSLONlu(JILOOP),ZOBSLATlu(JILOOP)
              inbvalz=SIZE(XVAR,3)
            CASE ('llhv')
              READ (8,*,END=888) &
                    ZOBSLATlu(JILOOP),ZOBSLONlu(JILOOP)
              inbvalz=SIZE(XVAR,3)
          END SELECT
          inbvalxy=inbvalxy+1
        ENDDO
        print *,' The program can take into account ', idimlonlat,' positions at the maximum'
        print *,' next values of the file ',TRIM(YFILEOBS) ,' will not be read'
  888      CONTINUE
        print *, ' End of reading of the observations localisation file ',TRIM(YFILEOBS)
        CALL CREATLINK('DIROBS',YFILEOBS,'CLEAN',ilocverbia)
        !
        IF ( inbvalz == 0) THEN
            ! niveaux vert. lus avec les coordonnees
            inbvalz=inbvalxy
            ! nombre de triplets= nombre de valeurs lues
            inbvalxyz=inbvalxy
        ELSE
            ! nombre de triplets= coordonnées lues * niveaux vert. interactifs
            inbvalxyz=inbvalxy*inbvalz
        ENDIF
        print *, ' Number of positions = ', inbvalxy
        print *, ' Number of vertical levels = ', inbvalz
        if (ilocverbia >= 4 ) then
            print *, 'lon, lat read :'
            print *,ZOBSLONlu,ZOBSLATlu
        endif
        !
        ! preparation des arguments pour SM_XYHAT : tableaux 2D
        ALLOCATE ( ZOBSLAT(inbvalxy,1), ZOBSLON(inbvalxy,1), ZOBSALT(inbvalz) )
        ZOBSLAT(1:inbvalxy,1)=ZOBSLATlu(1:inbvalxy)
        ZOBSLON(1:inbvalxy,1)=ZOBSLONlu(1:inbvalxy)
        ZOBSALT(1:inbvalz)  =ZOBSALTlu(1:inbvalz)
        DEALLOCATE (ZOBSLATlu,ZOBSLONlu,ZOBSALTlu)
        ALLOCATE ( ZOBSX(size(ZOBSLAT,1),size(ZOBSLAT,2)) )
        ALLOCATE ( ZOBSY(size(ZOBSLAT,1),size(ZOBSLAT,2)) )
        ! les 2 premiers arg. doivent etre XXHAT et XYHAT (pas XXX et XXY)
        !! peu importe en masdev4_6 car plus utilises.. 
        !CALL SM_XYHAT(XXHAT,XYHAT,XLATORI,XLONORI,&
        !! XXHAT,XYHAT supprimes en masdev4_7
        CALL SM_XYHAT(XLATORI,XLONORI,&
                      ZOBSLAT,ZOBSLON,ZOBSX,ZOBSY)
        if (ilocverbia >= 4 ) then
          ! XXX= XXHAT et XXY=XYHAT pour les 7 grilles
          print *, ' after SM_XYHAT, old limits X ',XXX(1,IGRID), XXX(SIZE(XVAR,1),IGRID)
          print *, 'new limits X ',ZOBSX(1,1),ZOBSX(inbvalxy,1)
          print *, 'old limits Y ',XXY(1,IGRID), XXY(SIZE(XVAR,2),IGRID)
          print *, 'new limits Y ',ZOBSY(1,1),ZOBSY(inbvalxy,1)
          DO JILOOP= 1,SIZE(XVAR,1) 
            print *, 'XXHAT ZOBSX ',XXX(JILOOP,IGRID),ZOBSX(JILOOP,1) 
          ENDDO
          DO JILOOP= 1,SIZE(XVAR,2) 
            print *, 'XYHAT ZOBSY ',XXY(JILOOP,IGRID),ZOBSY(1,JILOOP) 
          ENDDO
        endif
      ENDIF ! fin grille ZOBSX deja allouee
      !
      ! 777 = debut du traitement du tableau XVAR : utilise si DD ou FF
      ! pour reprise du traitement sur la deuxieme composante
777  CONTINUE
      ! 
      !        3.2.b interpolation  selon la verticale du champ Mesonh
      !              --------------------------------------------------
      ! cette interpolation verticale est realisee avant tout
      ! changement de la grille horizontale par l interpolation horizontale
      IF ( SIZE(XVAR,3)>1 .AND. SIZE(XVAR,2)>1 .AND. SIZE(XVAR,1)>1 ) THEN
      ! champ 3D 
        IF ( IGRID /=4 ) THEN
          print * , ' init of the model altitudes XZZ for NGRID=',IGRID
          ! init de XZZ pour cette grille
          ! car la routine readvar initialise XZZ pour NGRID=4 
          CALL COMPCOORD_FORDIACHRO(IGRID)
        ENDIF
        SELECT CASE ( YTYPEOUT(1:4))
          CASE ('LLZV','llzv','LLPV','llpv')
            ! interpolation  selon la verticale 
            print*,' Interpolation on ',YTYPEOUT(3:3),'=cst ',inbvalz,' levels'
            if (ilocverbia >= 1 ) then
              print *, 'levels= ',ZOBSALT 
            endif
            ALLOCATE (ZVARZCST(SIZE(XVAR,1),SIZE(XVAR,2),inbvalz))
            ALLOCATE (ZVAR3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
            ikdebzint=2
            ZVAR3D(:,:,:)=XVAR(:,:,:,1,1,1)
            IF ( YTYPEOUT(1:4)=='LLZV' .OR. YTYPEOUT(1:4)=='llzv' ) THEN
              CALL ZINTER(ZVAR3D,XZZ,ZVARZCST,ZOBSALT,ikdebzint,XSPVAL)
            ELSE IF ( YTYPEOUT(1:4)=='LLPV' .OR. YTYPEOUT(1:4)=='llpv' ) THEN
              CALL PINTER(ZVAR3D,IGRID,XSPVAL,ZOBSALT,ZVARZCST,ZPABS)
            ENDIF
            DEALLOCATE(XVAR)
            ALLOCATE(XVAR(SIZE(ZVARZCST,1),SIZE(ZVARZCST,2),SIZE(ZVARZCST,3),1,1,1))
            XVAR(:,:,:,1,1,1)=ZVARZCST
            zmini=MINVAL(XVAR(:,:,:,1,1,1),MASK=XVAR(:,:,:,1,1,1)/=XSPVAL)
            zmaxi=MAXVAL(XVAR(:,:,:,1,1,1),MASK=XVAR(:,:,:,1,1,1)/=XSPVAL)
            print * ,' After vertical interpolation, min,max of the field ',TRIM(CGROUP),'=', zmini,zmaxi
            DEALLOCATE(ZVARZCST,ZVAR3D)
            !
            ! ZOBSALT = grille verticale, tableau 1D passe en argument a zinter
            ! mise a jour du tableau 3D ZALT  passe en argument de WRITELLHV
            if ( ALLOCATED(ZALT) ) DEALLOCATE (ZALT)
            ALLOCATE ( ZALT(1,1,size(ZOBSALT,1)) )
            ZALT(1,1,:)=ZOBSALT
            ygrillevert='listevert'
            !
          CASE ('LLHV','llhv')  
            ! pas d interpolation verticale (h=grille modele)
            ygrillevert='XZZ'
            inbvalz=SIZE(XVAR,3)
            inbvalxyz=inbvalxy*inbvalz
            ! l interpolation horizontale sera faite apres l init de
            ! la nouvelle grille horizontale
        END SELECT
      !
      ELSE
      ! champ 2D  : pas d interpolation verticale
      ! la grille verticale utilisee est XXZS (i,j,NGRID)
        ygrillevert='XXZS'
        inbvalz3d=inbvalz  ! sauvegarde du nombre de niveaux verticaux
        inbvalxyz3d=inbvalxyz  ! sauvegarde du nombre de triplets
        inbvalz=1
        inbvalxyz=inbvalxy*inbvalz
        ! l interpolation horizontale sera faite apres l init de
        ! la nouvelle grille horizontale
      ENDIF
      !
      !        3.2.c interp. horizontale sur la nouvelle grille XY des obs 
      !              ----------------------------------------------------
      !
      !        3.2.c.1 interpolation horizontale du champ Mesonh
      !
      print *,' Interpolation to the new lat-lon grid of the field ',TRIM(CGROUP)
      ALLOCATE ( ZVAR3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
      ZVAR3D(:,:,:)=XVAR(:,:,:,1,1,1)
      ALLOCATE ( ZVARNEWH(size(ZOBSX,1),size(ZOBSX,2),SIZE(XVAR,3)) )
      if (ilocverbia >= 1 ) then
        print *, ' before HOR_INTERP_4PTS'
      endif
      CALL HOR_INTERP_4PTS (XXX(:,IGRID),XXY(:,IGRID),ZVAR3D,ZOBSX,ZOBSY,ZVARNEWH)
      DEALLOCATE(XVAR)
      ALLOCATE( XVAR(size(ZVARNEWH,1),size(ZVARNEWH,2),SIZE(ZVARNEWH,3),1,1,1) )
      XVAR(:,:,:,1,1,1)=ZVARNEWH(:,:,:)
      DEALLOCATE(ZVARNEWH,ZVAR3D)
      if (ilocverbia >= 1 ) then
        print *, ' after HOR_INTERP_4PTS'
      endif
      zmini=MINVAL(XVAR(:,:,:,1,1,1),MASK=XVAR(:,:,:,1,1,1)/=XSPVAL)
      zmaxi=MAXVAL(XVAR(:,:,:,1,1,1),MASK=XVAR(:,:,:,1,1,1)/=XSPVAL)
      print * ,'after horizontal interpolation, min,max of the field ',TRIM(CGROUP),'=', zmini,zmaxi
      !
      !        3.2.c.2 interpolation horizontale du tableau 3D de grille verticale
      SELECT CASE (ygrillevert (1:3))
        CASE ('XXZ')   ! champs 2D
          print *, ' Horizontal interpolation of XXZS for NGRID=',IGRID
          ALLOCATE ( ZVAR3D(SIZE(XXZS,1),SIZE(XXZS,2),1 ))
          ZVAR3D(:,:,1)=XXZS(:,:,IGRID)
          zmini=MINVAL(ZVAR3D(:,:,1),MASK=ZVAR3D(:,:,1)/=XSPVAL)
          zmaxi=MAXVAL(ZVAR3D(:,:,1),MASK=ZVAR3D(:,:,1)/=XSPVAL)
          print * ,'min,max of the vertical grid XXZS=', zmini,zmaxi
          ALLOCATE ( ZVARNEWH(size(ZOBSX,1),size(ZOBSX,2),1) )
          if (ilocverbia >= 1 ) then
            print *, ' before HOR_INTERP_4PTS'
          endif
          CALL HOR_INTERP_4PTS (XXX(:,IGRID),XXY(:,IGRID),ZVAR3D,ZOBSX,ZOBSY,ZVARNEWH)
          if ( ALLOCATED(ZALT) ) DEALLOCATE (ZALT)
          ALLOCATE( ZALT(size(ZVARNEWH,1),size(ZVARNEWH,2),IGRIDOUT) )
          ZALT(:,:,IGRIDOUT)=ZVARNEWH(:,:,1)
          DEALLOCATE(ZVARNEWH,ZVAR3D)
          zmini=MINVAL(ZALT(:,:,IGRIDOUT),MASK=ZALT(:,:,IGRIDOUT)/=XSPVAL)
          zmaxi=MAXVAL(ZALT(:,:,IGRIDOUT),MASK=ZALT(:,:,IGRIDOUT)/=XSPVAL)
          print * ,'after horizontal interpolation, min,max of the vertical grid =', zmini,zmaxi
        CASE ('XZZ')   ! champs 3D (LLHV)
          print *, ' Horizontal interpolation of XZZ'
          ALLOCATE ( ZVAR3D(SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)))
          ZVAR3D(:,:,:)=XZZ(:,:,:)
          zmini=MINVAL(ZVAR3D(:,:,:),MASK=ZVAR3D(:,:,:)/=XSPVAL)
          zmaxi=MAXVAL(ZVAR3D(:,:,:),MASK=ZVAR3D(:,:,:)/=XSPVAL)
          print * ,'min,max of the vertical grid XZZ=', zmini,zmaxi
          ALLOCATE ( ZVARNEWH(size(ZOBSX,1),size(ZOBSX,2),SIZE(XZZ,3)) )
          if (ilocverbia >= 1 ) then
            print *,' before HOR_INTERP_4PTS'
          endif
          CALL HOR_INTERP_4PTS (XXX(:,IGRID),XXY(:,IGRID),ZVAR3D,ZOBSX,ZOBSY,ZVARNEWH)
          IF ( ALLOCATED(ZALT) ) DEALLOCATE (ZALT)
          ALLOCATE( ZALT(size(ZVARNEWH,1),size(ZVARNEWH,2),SIZE(ZVARNEWH,3)) )
          ZALT(:,:,:)=ZVARNEWH(:,:,:)
          DEALLOCATE(ZVARNEWH,ZVAR3D)
          zmini=MINVAL(ZALT(:,:,:),MASK=ZALT(:,:,:)/=XSPVAL)
          zmaxi=MAXVAL(ZALT(:,:,:),MASK=ZALT(:,:,:)/=XSPVAL)
          print * ,'after horizontal interpolation, min,max of the vertical grid =', zmini,zmaxi
        CASE ('lis')   ! champs 3D (LLZV,LLPV)
          ! Pas d interpolation horizontale du tableau 
          !contenant la liste des niveaux verticaux
        CASE DEFAULT
          print *,'** type of vertical grid= ',TRIM(ygrillevert),' not correct'
          STOP
      END SELECT
      !
      !*      3.4  traitement sup si pluies cumulees
      !            -------------------------
      !
      IF (INDEX(CGROUP,'AC') /=0 ) THEN
        IF (.NOT.ALLOCATED(zwork3d))  THEN
          PRINT*, '- ACcumulated rain, do you want to make difference with a previous instant (o/O/y/Y/n/N) ?'
          READ(5,'(A1)')YREP
          CALL WRITEDIR(ILUDIR,YREP)
          CALL LOW2UP(YREP)
          IF (YREP=='Y' .OR. YREP=='O')THEN
            PRINT*, '- Name of diachro file (without .lfi) ?'
            READ(5,'(A28)',END=99) YFILEIN2
            CALL WRITEDIR(ILUDIR,YFILEIN2)
            ALLOCATE(zwork3d(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
            zwork3d(:,:,:)=XVAR(:,:,:,1,1,1)
            ALLOCATE(ZDATIME(16,SIZE(XDATIME,2)))
            ZDATIME(:,:)=XDATIME(:,:)
            YFLAGREADVAR='OPE'
            CALL READVAR(CGROUP,YFILEIN2,YFLAGREADVAR,ilocverbia,iret)
            if ( iret /= 0 ) then
              print *,TRIM(CGROUP),' not available'
              YFLAGREADVAR='CLO'
              CALL READVAR(CGROUP,YFILEIN2,YFLAGREADVAR,ilocverbia,iret2)
              YFLAGREADVAR='NOP'
              CYCLE
            endif
            ! pour traiter le deuxieme champ
            GO TO 777 
          ENDIF
        ENDIF
        IF (ALLOCATED(zwork3d) .AND. .NOT.ALLOCATED(zwork3d2))  THEN
          ALLOCATE(zwork3d2(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
          zwork3d2=XVAR(:,:,:,1,1,1)
          ! Stockage dans le tableau XVAR qui est le tableau ecrit
          XVAR(:,:,:,1,1,1)=XUNDEF
          WHERE( zwork3d(:,:,:) /= XUNDEF .AND. zwork3d2(:,:,:) /= XUNDEF) &
                 XVAR(:,:,:,1,1,1)=zwork3d(:,:,:)-zwork3d2(:,:,:)
          ! sauvegarde de la valeur de CGROUP
          YGROUP=CGROUP
          CGROUP=ADJUSTL( ADJUSTR(CGROUP)//'diff')
          ! pour avoir le temps du 1er fichier
          XDATIME(:,:)=ZDATIME(:,:)
          DEALLOCATE(zwork3d,ZDATIME)
        ENDIF
      ENDIF
      !
      !*      3.5  Ecriture  du tableau XVAR (module MODD_ALLOC_FORDIACHRO)
      !            -------------------------
      !
      print *,' Format of writing= ', YTYPEOUT(1:4)
      print *,'size of XVAR ',SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)
      ivarideb=NREADIL
      ivarifin=NREADIH
      ivarjdeb=NREADJL
      ivarjfin=NREADJH
      ivarkdeb=NREADKL
      ivarkfin=NREADKH
      ivartinf=1
      ivartsup=1
      ivartrajinf=1
      ivartrajsup=1
      ivarprocinf=1
      ivarprocsup=1
      if (ilocverbia >= 1 ) then
        print *,'size of ZALT ',SIZE(ZALT,1),SIZE(ZALT,2),SIZE(ZALT,3)
        IF(SIZE(ZALT,3)<=10)THEN
          print *,ZALT(1:SIZE(ZALT,1),1:SIZE(ZALT,2),1:SIZE(ZALT,3))
        ELSE
          print *,ZALT(:,:,1:10)
        ENDIF
      endif
      !
      ! Ecriture triplet par triplet lat,lon,alt  traites
      !
      print *,' number of triplets taken into account =',inbvalxyz
      IF ( inbvalxyz == inbvalxy*inbvalz) THEN
      ! cas fichier d obs contient seulement les coordonnees
      ! les niveaux sont passes en interactif: double boucle sites
      ! puis niveaux
        DO JNobsLOOPsite=1,inbvalxy
          DO JNobsLOOPz=1,inbvalz
           if (ilocverbia >= 0 ) then
            print *,' obs ',JNobsLOOPsite,' lat lon alt',ZOBSLAT(JNobsLOOPsite:JNobsLOOPsite,1),&
            ZOBSLON(JNobsLOOPsite:JNobsLOOPsite,1),ZALT(1,1,JNobsLOOPz:JNobsLOOPz)
              if (ilocverbia >= 1 ) then
                print *,' size XVAR', SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)
              endif
           endif
            CALL WRITELLHV(JNobsLOOPsite,JNobsLOOPsite,1,1,JNobsLOOPz,JNobsLOOPz,&
                       ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,ilocverbia,iret,&
                       HFILENAME_SUP='obs',&
                       PLON=ZOBSLON,PLAT=ZOBSLAT,PALT=ZALT)
             print *, ' WRITELLHV, return value=',iret
            !  indiquera  a WRITELLHV que le fichier courant en ecriture est 
            !  deja ouvert et de ne pas ecrire l entete
            YFLAGWRITE='OLDNH'
          ENDDO
          IF ( inbvalz == 1) THEN
                 ! une seule valeur par site donc pas d entete entre 2 sites
            YFLAGWRITE='OLDNH'
          ELSE
                 ! nouvelle entete pour le site suivant
            YFLAGWRITE='OLD1H'
          ENDIF
        ENDDO
      !
      ELSE
      ! cas fichier d obs contient les coordonnees et les altitudes
      ! simple boucle sur le nombre de triplets
        DO JNobsLOOPtriplet=1,inbvalxy
          if (ilocverbia >= 0 ) then
              print *,' obs ',JNobsLOOPtriplet,'lat lon alt',ZOBSLAT(JNobsLOOPtriplet:JNobsLOOPtriplet,1),&
              ZOBSLON(JNobsLOOPtriplet:JNobsLOOPtriplet,1),ZALT(1,1,JNobsLOOPtriplet:JNobsLOOPtriplet)
            if (ilocverbia >= 1 ) then
              print *,' size XVAR', SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)
            endif
          endif
          CALL WRITELLHV(JNobsLOOPtriplet,JNobsLOOPtriplet,1,1,JNobsLOOPtriplet,JNobsLOOPtriplet,&
                     ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                     CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,ilocverbia,iret,&
                     HFILENAME_SUP='obs',&
                     PLON=ZOBSLON,PLAT=ZOBSLAT,PALT=ZALT)
          print *, ' WRITELLHV, return value=',iret
            !  indiquera  a WRITELLHV que le fichier courant en ecriture est 
            !  deja ouvert et de ne pas ecrire l entete
          YFLAGWRITE='OLDNH'
        ENDDO
      ENDIF
      ! restore le nombre de niveaux verticaux et de triplets
      IF (ygrillevert=='XXZS') THEN
        inbvalz=inbvalz3d
        inbvalxyz=inbvalxyz3d
      END IF
      ! fermeture du 2e fichier ouvert pour diff de pluie cumulee
      IF (INDEX(CGROUP,'AC') /=0 .AND.  ALLOCATED(zwork3d2))  THEN
        DEALLOCATE(zwork3d2)
         !CALL READVAR(YGROUP,YFILEIN2,'CLO',ilocverbia,iret2)
         ! il faut close avec un champ forcement present mais
         ! pas AC...diff donc avec YGROUP qui memorise le nom de groupe
         ! existant sans diff
         ! ce close fait planter le prog: erreur non trouvée apres
         ! 1/j de recherche
         !  de toute facon la fermeture se fera avec la fin de programme
      ENDIF
    ELSE  ! iret /=0
      print *, ' READVAR, return value=',iret
    ENDIF 
    ! Pour indiquer l ecriture de l entete pour la variable suivante
  YFLAGWRITE='OLD1H'
  ENDDO ! boucle champ a traiter
  ! pour clore le traitement meme si la liste des champs est
  ! incomplete ( non terminee par END)
  88  CONTINUE
  CGROUP='END'
!
!---------------------------------------------------------------------------
!
!*       4.    Fermeture fichiers
!              ------------------
!
  IF ( CGROUP(1:3) == 'END' .AND. YFLAGWRITE(1:3)/='NEW') THEN
    PRINT*, 'END -> Close the output file'
    YFLAGWRITE='CLOSE'
    ! dans cet appel seul l argument YFLAGWRITE est pris en compte, tous
    ! les autres arguments sont ignores
    SELECT CASE(YTYPEOUT(1:4))
      CASE('LLHV','llhv','LLPV','llpv','LLZV','llzv')
        CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                      ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                      CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,ilocverbia,iret,&
                      HFILENAME_SUP='obs')
      CASE DEFAULT
        PRINT*, 'Closure of output type ',YTYPEOUT ,' not coded'
    END SELECT
    ! renomme le fichier de sortie en ajoutant le nom du fichier d obs
    ! effectue en fin de traitement car pour les routines FM, les noms de fichiers 
    !sont limites a 28 caracteres
    YFILEOUT=ADJUSTR(YFILEIN)//ADJUSTL(YTYPEOUT(1:4))
    YFILEOUT=ADJUSTL( ADJUSTR(YFILEOUT)//'_obs')
    ycommand='mv '//TRIM(YFILEOUT)//' '//TRIM(TRIM(YSAVEFILEOBS)//'_'//TRIM(YFILEOUT))
    print *,'command= ',ycommand
    call SYSTEM ( TRIM(ycommand) )    
  ENDIF
!
ENDDO ! fin boucle des fichiers a traiter
!-------------------------------------------------------------------------------
!
!*       5.    Fin de boucle sur les fichiers
!              ------------------
!
99 CONTINUE
!
!   Suppression de tous les liens eventuellemnet crees
YDUMMYFILE=''
CALL CREATLINK(' ',YDUMMYFILE,'CLEAN',ilocverbia)
PRINT*, 'The file ',TRIM(YLUDIR),' stores all the input directives'
PRINT*, ' you must give a new name to use it again'
CLOSE(ILUDIR)
!
IF ( YFLAGWRITE(1:3)/='NEW') THEN
  PRINT*, 'Output files ',TRIM(YSAVEFILEOBS),'*obs are available'  
END IF
!
END PROGRAM MESONH2OBS
