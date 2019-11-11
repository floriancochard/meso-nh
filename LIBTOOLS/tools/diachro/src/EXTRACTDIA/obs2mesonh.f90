!     ###################
      PROGRAM  OBS2MESONH
!     ###################
!
!!****  *OBS2MESONH* -  Interpolation d une liste de valeurs observées
!                       sur la grille Mesonh
!                       en entrée un fichier ASCII au format [ ] pour optionnel
!                                 [YYYMMJJHHMMSS]
!                                 lon lat [alt] valeur_obs  
!                              ou lat lon [alt] valeur_obs 
!                       en sortie un fichier diachronique
!! 
!!
!!    PURPOSE
!!    -------
! 
!
!!**  METHOD
!!    ------
!!      
!     Lecture en entree:
!  d un fichier ascii contenant les localisations (lon,lat alt,valeur) a traiter
!  d un fichier modele pour recuperer la grille XYZ
!     Pour chaque obs lue, 1) recherche du point de grille xy la contenant,
!                          2) recherche sur la verticale du niveau K la contenant
!     Si plusieurs obs sont contenues dans un meme point de grille, calcul de la moyenne des ces obs
!     Pour certaines variables (unite dBz, nom de champ_WVBT ou _IRBT), 
!   passage a des unites plus pertinentes pour effectuer les moyennes
!   et retour aux unites d origine avant ecriture (voir les 2 routines
! symetriques To_computing_units et From_computing_units)
!     Mise a XSPVAL des points de grille ne contenant pas d'obs
!      
!     Ecriture en sortie:
!  d un fichier  diachronique ( utiliser LSPOT=T dans diaprog pour visualiser 
!       toutes les points grilles non XSPVAL)
!
!!
!!    EXTERNAL
!!    --------
!!          CREATLINK : à l'ouverture du fichier, HYFLAGFILE='OPE',
!!                      création d'un lien dans le directory local
!!                      si le fichier existe sous $DIROBS
!!          READVAR   : lecture d unchamp du fichier diachronique
!!          WRITEVAR : ecriture format lon lat alt val
!!          SM_XYHAT  : création de la grille des Obs
!!          TO_COMPUTING_UNITS: passage unites vers unites plus pertinentes 
!!                              pour effectuer des calculs
!!          FROM_COMPUTING_UNITS: passage inverse avant ecriture
!!                               (appele par writevar)
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHORS
!!    -------
!!    N. Asencio and J. Stein
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/11/2003
!!      Fev 2005: ajout de champs diachroniques ALT_champ N_champ
!!                changement de grille pour le vent (zonal,meridien->
!!                grille Mesonh)
!!      04/05/2005 add a control for the min and max of the field before
!                and after interpolation(s)
!                  observations outside the mesonh domain are rejected
!       19/09/2005 G.Jaubert CNRM  
!                  1) Nom du fichier .lfi en output demande 
!                  2) l'enregistrement peut ne pas contenir alt si champ 2D
!                  3) si le fichier de donnees commence par une date,
!                     reinitialisation de la date dans le lfi de sortie
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
USE MODD_GRID, ONLY: XLON0,XLAT0,XLONORI,XLATORI
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
USE MODD_READLH !NREADIL,IH,...
!
!
USE MODI_WRITEVAR
USE MODI_CREATLINK
USE MODI_LOW2UP      
USE MODI_WRITEDIR      
USE MODI_UV_TO_ZONAL_AND_MERID
USE MODI_TO_COMPUTING_UNITS
!
IMPLICIT NONE                       
!
!*       0.1   Local variables declaration
!
!! indices de boucle
INTEGER     :: JILOOP,JOBS
! zoom  suivant les 6 dimensions des champs diachro
INTEGER     :: ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin
INTEGER     :: ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup
INTEGER     :: IL   ! indice de positionnement dans yligne
INTEGER     :: IIMNH,IJMNH,IKMNH,IKMAX,IGRID
!
INTEGER     :: IAN, IMOIS, IJOUR, IHEUR, IMINU, ISEC ! date observation
REAL        :: XSEC ! heure observation en secondes 
!                                   
! stockage
REAL, allocatable, dimension(:,:,:)    ::  ZOBSinMNH,ZALTOBS
INTEGER , allocatable, dimension(:,:,:) ::  ICPTOBSinMNH    
REAL              :: ZXOBS,ZYOBS,ZOBSLONLU,ZOBSLATLU,ZOBSALTLU,ZVALOBS
! pour le passage composantes meridienne,zonale a la grille Mesonh
REAL, allocatable, dimension(:,:,:,:,:,:):: ZVENTSAVE
REAL, allocatable, dimension(:,:,:)      :: ZWORK3D,ZWORK3D2
!
REAL :: zmini,zmaxi
REAL :: ZVAR1, ZVAR2, ZVAR3
LOGICAL :: galtobs
INTEGER           :: ILUDIR,iret,ilocverbia,inbval,inbvalrej,IKM1
!! **** la taille des variables caracteres contenant les noms
!!      de fichiers est obligatoirement de 28 ****
!!      pour toutes les routines diachro
CHARACTER(LEN=28) :: YFILEGRID,YDUMMYFILE , YFILEOUTNAME
CHARACTER(LEN=100):: YFILEOBS
CHARACTER(LEN=3)  :: YFLAGREADVAR ,YFLAGWRITE
CHARACTER (LEN=10):: YUNITE,YUNITEMAJ
CHARACTER (LEN=3) :: YSTOCK,YFILEOUTSUFFIX
CHARACTER (LEN=4) :: YLL
CHARACTER(LEN=11) :: YLUDIR      !  Name of the dir file
CHARACTER(LEN=14) :: YDATEOBS    ! observation  date (YYYYMMDDHHMISS) 
CHARACTER(LEN=100) :: YLIGNE    ! 
!-------------------------------------------------------------------------------
!
!*       1.    Init
!              ----
!
! active(1) ou desactive(0) les prints de controle dans les routines
! READVAR et WRITEVAR
ilocverbia=0
! 
! dans mesonh Xundef est utilise 
! dans les routines diachro XSPVAL est utilisé
 XSPVAL=XUNDEF                                    
!
! ouverture d un fichier dir ou vont s ecrire les entrees clavier
YLUDIR='dirobs2mnh'
CALL FMATTR(YLUDIR,YLUDIR,ILUDIR,iret)
OPEN(UNIT=ILUDIR,FILE=YLUDIR,FORM='FORMATTED')
!
NIINF=0
NISUP=0
NJINF=0
NJSUP=0
iret=0
!
!*      2.  Lecture et initialisation des modules Mesonh
!          ----------------------------

PRINT*, '- Name of the diachro file to read the grid '&
      ,'(without .lfi) ?'
READ(5,'(A)',END=99) YFILEGRID
YFILEGRID=ADJUSTL(YFILEGRID)      
CALL WRITEDIR(ILUDIR,YFILEGRID)
!
PRINT*, '- Prints : 0= mini 1=debug mode in obs2mesonh'
PRINT*, '                   2= print of input values'
PRINT*, '                   3=debug mode in diachro routines'
PRINT*, '?'
READ(5,*)ilocverbia
CALL WRITEDIR(ILUDIR,ilocverbia)
PRINT*, ' output prints= ',ilocverbia 
IF (ilocverbia >2) nverbia=ilocverbia
!
! Lecture du  champ ZSBIS  pour  obtenir
! l Initialisation des variables
! des modules (cf USE en debut de programme)
!
!  indique que le fichier lu doit etre ouvert dans READVAR
YFLAGREADVAR='OPE'
CALL READVAR('ZSBIS',YFILEGRID,YFLAGREADVAR,ilocverbia,iret)
print *, 'READVAR(zsbis), return value= ',iret
IF ( iret /= 0 ) THEN 
  print *,'** Error when reading the grid in the FM diachro file: ',TRIM(YFILEGRID)
  STOP
ENDIF
!
if (ilocverbia >= 0 ) then
  print *,' Size of input array(zs)= ',SIZE(XVAR,1),SIZE(XVAR,2),&
            SIZE(XVAR,3),SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)
  PRINT*, 'NIINF,NISUP,NJINF,NJSUP', NIINF,NISUP,NJINF,NJSUP
  PRINT*, 'LATORI,LONORI ', XLATORI,XLONORI
  zmini= 0.5 * (XXHAT(1)+XXHAT(2))
  zmaxi= 0.5 * (XYHAT(1)+XYHAT(2))
  CALL SM_LATLON(XLATORI,XLONORI,zmini,zmaxi,ZVAR1,ZVAR2)
  PRINT*, 'LATOR,LONOR ', ZVAR1,ZVAR2
endif
!
!
PRINT*,'- Name of the output file ? (the default [CR or empty line], "obs" will be added to the input FM file'
READ(5,'(A28)',END=88) YFILEOUTNAME
YFILEOUTNAME=ADJUSTL(YFILEOUTNAME)
CALL WRITEDIR(ILUDIR,YFILEOUTNAME)
IF (YFILEOUTNAME(1:2) == 'll' .OR. YFILEOUTNAME(1:2) == 'LL' ) THEN
      print * ,'** OBS2MESONH: the format of the input file was modified Oct2005'
      print * ,' the 3rd line is for the name of the output file'
      print * ,'instead of the format of the input file: modify your directives'
      STOP
ELSE
    ! input phase avec obs2mesonh post Oct2005
  IF (LEN_TRIM(YFILEOUTNAME) == 0) THEN
    YFILEOUTSUFFIX='obs'
    YFILEOUTNAME=YFILEGRID
  ELSE
    YFILEOUTSUFFIX='NEN'
    PRINT*,'Output file name=', TRIM(YFILEOUTNAME)
  ENDIF
ENDIF
!
! fichier a creer dans WRITEvar
YFLAGWRITE='NEW'
! ecriture de ZSBIS dans FILEOUT
CALL WRITEVAR(NREADIL,NREADIH,NREADJL,NREADJH,NREADKL,NREADKL,&
              1,SIZE(XVAR,4),1,SIZE(XVAR,5),1,SIZE(XVAR,6), &
              'ZSBIS',YFILEOUTNAME,YFLAGWRITE,YFILEOUTSUFFIX,ilocverbia,iret)
!  indiquera  a WRITEVAR que le fichier courant en ecriture est deja ouvert
YFLAGWRITE='OLD'
!-------------------------------------------------------------------------------
!
!*       3.    LOOP on obs files
!              -----------------
!    
PRINT*,' Loop on the observation files:'
DO JOBS=1,10000
  777 CONTINUE  ! Point de reprise pour le traitement de la 2e composante vent
  !
  PRINT*,'- Format of the input observation file:' 
  PRINT*,'  LL=   n lines Lon,Lat,val '
  PRINT*,'  ll=   n lines lat,lon,val '
  PRINT*,'  DLL=  date (YYYYMMDDHHMISS) then n lines Lon,Lat,val '
  PRINT*,'  Dll=  date (YYYYMMDDHHMISS) then n lines lat,lon,val '
  PRINT*,'  LLa=  n lines Lon,Lat,alt(m),val'
  PRINT*,'  lla=  n lines lat,lon,alt(m),val'
  PRINT*,'  DLLa= date (YYYYMMDDHHMISS) then n lines Lon,Lat,alt(m),val'
  PRINT*,'  Dlla= date (YYYYMMDDHHMISS) then n lines lat,lon,alt(m),val'
  PRINT*,'(END to stop)?'
  READ(5,'(A)',END=88) YLL
  YLL=ADJUSTL(YLL)
  CALL WRITEDIR(ILUDIR,YLL)
  IF (YLL(1:3)=='end' .OR. YLL(1:3)=='END') GO TO 88
  IF (YLL(1:3)=='LLa') THEN
    print*, 'format Lon,Lat,alt(m),value'
    galtobs=.true.
  ELSE IF (YLL(1:3)=='lla') THEN
    print*, 'format lat,lon,alt(m),value'
    galtobs=.true.
  ELSE IF (YLL(1:4)=='DLLa') THEN
    print*, 'format Date then Lon,Lat,alt(m),valeur'
    galtobs=.true.
  ELSE IF (YLL(1:4)=='Dlla') THEN
    print*, 'format Date then lat,lon,alt(m),valeur'
    galtobs=.true.
  ELSE IF (YLL(1:2)=='LL') THEN
    print*, 'format Lon,Lat,valeur'
    galtobs=.false.
  ELSE IF (YLL(1:2)=='ll') THEN
    print*, 'format lat,lon,valeur'
    galtobs=.false.
  ELSE IF (YLL(1:3)=='DLL') THEN
    print*, 'format Date then Lon,Lat,value'
    galtobs=.false.
  ELSE IF (YLL(1:3)=='Dll') THEN
    print*, 'format Date then lat,lon,value'    
    galtobs=.false.
  ELSE
    print*, '** incorrect format ',YLL
    CYCLE
  ENDIF
  PRINT*,'- Name of the input observation file ?'
  READ(5,'(A)',END=88) YFILEOBS
  YFILEOBS=ADJUSTL(YFILEOBS)
  CALL WRITEDIR(ILUDIR,YFILEOBS)
  !
  !*       3.1   Lecture du fichier d obs a traiter
  !              ----------------------
  PRINT*, '- Name of the new field to be created:'
  PRINT*, '(if it is a wind field you have to name the field : '
  PRINT*, ' WTxx: the field is localised at vertical flux points, ',&
          'otherwise at mass points (example : WT10) '
  PRINT*, ' UTxx: the field (U-component for zonal) will be converted to ',&
          'MesoNH wind components'
  PRINT*, 'the V-component must be provided immediately after with VTxx'
  PRINT*, '?'
  READ(5,'(A9)',END=88) CGROUP
  CGROUP=ADJUSTL(CGROUP)
  CALL WRITEDIR(ILUDIR,CGROUP)
  CALL LOW2UP(CGROUP)      
  PRINT*, '- Unit of the new field ?'
  READ(5,'(A)') YUNITE
  YUNITE=ADJUSTL(YUNITE)
  CALL WRITEDIR(ILUDIR,YUNITE)
  PRINT*, '- Profil of the new field :'
  PRINT*, '   3D=XYZ '
  PRINT*, '   2D=XY  (obs altitudes not taken into account)'
  PRINT*, '   1D=Z   (vertical profil (_PV_ for diaprog) localised at ',&
           'lat-lon of the 1st obs'
  PRINT*, ' 1D/2D/3D ?'        
  READ(5,'(A)') YSTOCK      
  YSTOCK=ADJUSTL(YSTOCK)
  CALL WRITEDIR(ILUDIR,YSTOCK)
  IF ( (YSTOCK == '3D' .OR. YSTOCK == '1D') .AND. .NOT.(galtobs) ) THEN
      print * ,'** It is not possible to store ',TRIM(YSTOCK),' profil ',&
               'because no altitude was provided in the input obs file'
      print *, ' change your inputs:', TRIM(YLL)
      STOP
  ENDIF
  ! tableau de stockage des valeurs des obs et compteur de ces valeurs stockees
  IF(ALLOCATED(ZOBSinMNH)) DEALLOCATE(ZOBSinMNH)
  IF(ALLOCATED(ZALTOBS)) DEALLOCATE(ZALTOBS)
  IF(ALLOCATED(ICPTOBSinMNH)) DEALLOCATE(ICPTOBSinMNH)
  ! XVAR = futur tableau a ecrire via writevar
  IF(ALLOCATED(XVAR)) DEALLOCATE(XVAR)
  IF ( YSTOCK == '3D' .OR. YSTOCK == '1D' ) THEN
    ALLOCATE(XVAR( SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3),1,1,1))
    ALLOCATE(ZOBSinMNH( SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)))
    ALLOCATE(ZALTOBS( SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)))
    ALLOCATE(ICPTOBSinMNH( SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)))
  ELSE
    ALLOCATE(XVAR( SIZE(XZZ,1),SIZE(XZZ,2),1,1,1,1))
    ALLOCATE(ZOBSinMNH( SIZE(XZZ,1),SIZE(XZZ,2),1))
    ALLOCATE(ZALTOBS( SIZE(XZZ,1),SIZE(XZZ,2),1))
    ALLOCATE(ICPTOBSinMNH( SIZE(XZZ,1),SIZE(XZZ,2),1))
  END IF
  !
  PRINT*, 'Mesonh field to be created: ', TRIM(CGROUP),' ',TRIM(YUNITE),' ',TRIM(YSTOCK)
  ! init de la grille verticale Mesonh suivant le nom de variable
  SELECT CASE (CGROUP(1:1))
    CASE ('W') 
      IGRID=4        ! champ d obs sur la grille W
      ! init du tableau des altitudes XZZ pour la grille masse
      CALL COMPCOORD_FORDIACHRO(1)            
      if (ilocverbia > 0 ) then
        print *,' after COMPCOORD_FORDIACHRO (mass grid for W field)'
      endif
      !                             -----------XZZ(k)   grille masse
      !                                 x W(k)
      !                             -----------XZZ(k-1) grille masse
    CASE default 
      IGRID=1        ! champ d obs sur la grille de masse 
      ! init du tableau des altitudes XZZ pour la grille W 
      CALL COMPCOORD_FORDIACHRO(4)            
      if (ilocverbia > 0 ) then
        print *,' after COMPCOORD_FORDIACHRO (W grid for mass field)'
      endif
      !                             -----------XZZ(k+1) grille W
      !                                 x T(k)
      !                             -----------XZZ(k)   grille W
  END SELECT
  !
  !* 3.2.  pour chaque obs lue, recherche de la maille mesonh
  !        contenant cette obs : cumul dans cette maille
  !                              mis a jour du compteur d obs par maille
  !              ----------------------
  ! 
  print *,'YFILEOBS=', TRIM(YFILEOBS)
  CALL CREATLINK('DIROBS',YFILEOBS,'CREAT',ilocverbia)
  OPEN (UNIT=8,FILE=TRIM(ADJUSTL(YFILEOBS)),STATUS='OLD',FORM='FORMATTED')
  !
  ZOBSINMNH(:,:,:)=0.
  ICPTOBSinMNH(:,:,:)=0
  IKMAX=1
  !
  inbval=0
  inbvalrej=0
  if (ilocverbia >= 2 ) print *, 'before reading of input obs file'
  IF (YLL(1:1)=='D' ) THEN
    ! lecture de la date
    READ (8,'(A14)',ERR=886) YDATEOBS
    YDATEOBS=ADJUSTL(YDATEOBS)
    ! verification YDATEOBS est une date
    IF (YDATEOBS(1:4)<='1900' .OR. YDATEOBS(1:4)>='2020' &
      .OR. YDATEOBS(5:6)<'01' .OR. YDATEOBS(5:6)>'12' &
      .OR. YDATEOBS(7:8)<'01' .OR. YDATEOBS(7:8)>'31' &
      .OR. YDATEOBS(9:10)<'00' .OR. YDATEOBS(9:10)>'23' &
      .OR. YDATEOBS(11:12)<'00' .OR. YDATEOBS(11:12)>'59' &
      .OR. YDATEOBS(13:14)<'00' .OR. YDATEOBS(13:14)>'59' ) GO TO 887
    ! Pourquoi cette init ? voir GJ iret=1
    ! reinitialisation de XDATIME
    READ(YDATEOBS,'(i4,5I2)') IAN, IMOIS, IJOUR, IHEUR, IMINU, ISEC
    XSEC=IHEUR*3600.+IMINU*60.+ISEC
    DO JILOOP=1,13,4
      XDATIME(JILOOP,1)=IAN
      XDATIME(JILOOP+1,1)=IMOIS
      XDATIME(JILOOP+2,1)=IJOUR
      XDATIME(JILOOP+3,1)=XSEC
    END DO
    if (ilocverbia >= 2 ) print *,'XDATIME initialised to:',XDATIME(1:4,1)
  ENDIF
  ! lecture de la position et valeur des observations
  ! boucle infinie : arret sur fin de fichier
  ! init de la position à -999
  IIMNH=-999
  IJMNH=-999
! modif GJ deduction format  JILOOP=0
  DO
! modif GJ deduction format    IF ( JILOOP == 0 .AND. YSTOCK =='2D'  ) THEN
! modif GJ deduction format       ! le format de donnees peut ne pas contenir alt 
! modif GJ deduction format       ! lecture du premier enregistrement en caracteres puis decodage
! modif GJ deduction format      if (ilocverbia >= 2 ) print *,'Recherche du nombre de variables dans un enregistrement'
! modif GJ deduction format      READ (8,'(A100)') YLIGNE
! modif GJ deduction format      YLIGNE=ADJUSTL(YLIGNE)
! modif GJ deduction format      il=index(YLIGNE,' ')-1
! modif GJ deduction format      READ(Yligne(1:il),*) ZVAR1
! modif GJ deduction format      YLIGNE=ADJUSTL(YLIGNE(il+1:))
! modif GJ deduction format      il=index(YLIGNE,' ')-1
! modif GJ deduction format      READ(Yligne(1:il),*) ZVAR2
! modif GJ deduction format      IF (TRIM(YLL)=='LL' .or. TRIM(YLL)=='DLL') THEN
! modif GJ deduction format        ZOBSLONlu=ZVAR1
! modif GJ deduction format        ZOBSLATlu=ZVAR2
! modif GJ deduction format      ELSE IF (YLL(1:2)=='ll' .or. YLL(1:3)=='Dll') THEN
! modif GJ deduction format        ZOBSLATlu=ZVAR1
! modif GJ deduction format        ZOBSLONlu=ZVAR2
! modif GJ deduction format      ENDIF
! modif GJ deduction format      YLIGNE=ADJUSTL(YLIGNE(il+1:))
! modif GJ deduction format      il=index(YLIGNE,' ')-1
! modif GJ deduction format      READ(Yligne(1:il),*) ZVAR3
! modif GJ deduction format      IF ( LEN_TRIM(YLIGNE(il+1:)) == 0 ) THEN
! modif GJ deduction format        print *,' Champ 2D avec enregistrement sans altitude'
! modif GJ deduction format        ZOBSALTlu=-999
! modif GJ deduction format        ZVALOBS=ZVAR3
! modif GJ deduction format        YLL(LEN_TRIM(YLL):LEN_TRIM(YLL))='A'
! modif GJ deduction format      ELSE
! modif GJ deduction format        ZOBSALTlu=ZVAR3
! modif GJ deduction format        YLIGNE=ADJUSTL(YLIGNE(il+1:))
! modif GJ deduction format        READ(Yligne,*) ZVALOBS
! modif GJ deduction format      ENDIF
! modif GJ deduction format      if (ilocverbia >= 2 ) THEN
! modif GJ deduction format        print *,'1er enregistrement: lat=',ZOBSLATlu,' lon=',ZOBSLONlu
! modif GJ deduction format        print *,'                    alt=',ZOBSALTlu,' var=',ZVALOBS
! modif GJ deduction format      endif
! modif GJ deduction format    ELSE
! modif GJ deduction format      ! lecture d'un enregistrement
! modif GJ deduction format      IF (TRIM(YLL)=='LL' .or. TRIM(YLL)=='DLL') THEN
! modif GJ deduction format        READ (8,*,ERR=887,END=888) ZOBSLONlu,ZOBSLATlu,ZOBSALTlu,ZVALOBS
! modif GJ deduction format      ELSE IF (YLL(1:2)=='ll' .or. YLL(1:3)=='Dll') THEN
! modif GJ deduction format        READ (8,*,ERR=887,END=888) ZOBSLATlu,ZOBSLONlu,ZOBSALTlu,ZVALOBS
! modif GJ deduction format      ELSE IF (TRIM(YLL)=='LA' .or. TRIM(YLL)=='DLA') THEN
! modif GJ deduction format        READ (8,*,ERR=887,END=888) ZOBSLONlu,ZOBSLATlu,ZVALOBS
! modif GJ deduction format        ZOBSALTlu=-999
! modif GJ deduction format      ELSE IF (YLL(1:2)=='lA' .or. YLL(1:3)=='DlA') THEN
! modif GJ deduction format        READ (8,*,ERR=887,END=888) ZOBSLATlu,ZOBSLONlu,ZVALOBS
! modif GJ deduction format        ZOBSALTlu=-999
! modif GJ deduction format      ELSE
! modif GJ deduction format        print * ,' Format des obs =',YLL(1:4),' valeur incorrecte'
! modif GJ deduction format        print *, 'valeurs possibles: ll ou LL ou Dll ou DLL ou llh ou LLh ou Dllh ou DLLh'
! modif GJ deduction format        STOP
! modif GJ deduction format      ENDIF
! modif GJ deduction format    ENDIF
! modif GJ deduction format    JILOOP=1
    IF (YLL(1:3)=='LLa' .OR. YLL(1:4)=='DLLa') THEN
      READ (8,*,END=888) ZOBSLONlu,ZOBSLATlu,ZOBSALTlu,ZVALOBS
    ELSE IF (YLL(1:3)=='lla'.OR. YLL(1:4)=='Dlla' ) THEN
      READ (8,*,END=888) ZOBSLATlu,ZOBSLONlu,ZOBSALTlu,ZVALOBS
    ELSE IF (YLL(1:2)=='LL'.OR.YLL(1:3)=='DLL' ) THEN
      READ (8,*,END=888) ZOBSLONlu,ZOBSLATlu,ZVALOBS
      ZOBSALTlu= XSPVAL
    ELSE IF (YLL(1:2)=='ll' .OR. YLL(1:3)=='Dll' ) THEN
      READ (8,*,END=888) ZOBSLATlu,ZOBSLONlu,ZVALOBS
      ZOBSALTlu= XSPVAL
    ELSE
      print * ,'** Obs format=',YLL(1:4),' is an incorrect value'
      print *, 'correct values are: ll or LL or Dll or DLL or lla or LLa or Dlla or DLLa'
      STOP
    ENDIF  

    ! recupere les coordonnées de l obs sur le plan conforme
    IF (YSTOCK == '3D' .OR. YSTOCK == '2D' .OR. &
                         (YSTOCK == '1D' .AND. IIMNH == -999) ) THEN
      ! recupere pour chaque obs si 2D ou 3D , pour la premiere obs si 1D
      !(les 2 premiers arg. doivent etre XXHAT et XYHAT (pas XXX et XXY))
      !! peu importe en masdev4_6 car plus utilises..
      !CALL SM_XYHAT(XXHAT,XYHAT,XLATORI,XLONORI, &
      !! XXHAT,XYHAT supprimes en masdev4_7
      CALL SM_XYHAT(XLATORI,XLONORI, &
                     ZOBSLATlu,ZOBSLONlu,ZXOBS,ZYOBS)
      ! quelle est la maille horizontale mesonh qui contient cette obs ?
      ! XXHAT(I),XXHAT(I+1) = limites X de la maille I
      ! XYHAT(J),XYHAT(J+1) = limites Y de la maille J
      IF ( ZXOBS >= XXHAT(2) .AND. ZXOBS <= XXHAT(NIMAX+2-1) .AND.&
           ZYOBS >= XYHAT(2) .AND. ZYOBS <= XYHAT(NJMAX+2-1) ) THEN
        IIMNH=MAX(MIN(COUNT(XXHAT(:)<ZXOBS),NIMAX+2-1),2)
        IJMNH=MAX(MIN(COUNT(XYHAT(:)<ZYOBS),NJMAX+2-1),2)
      ELSE
        print * ,'*** The observation at lat,lon ',ZOBSLATlu,ZOBSLONlu,&
                 'is out of the Mesonh domain, not treated ***'
        inbvalrej=inbvalrej+1
        CYCLE
      ENDIF
    ELSE
     if (ilocverbia >= 2 ) then
       print *, ' Profil ', YSTOCK, ': following obs are at the same localisation',&
                ' as the first one ', ZOBSLATlu,ZOBSLONlu
       print *, 'i,j=', IIMNH,IJMNH
     endif
    ENDIF

    if (ilocverbia >= 3 ) then
       print *, ZXOBS,IIMNH &
               ,XXX(IIMNH,IGRID),XXX(IIMNH-1,IGRID),XXHAT(IIMNH),&
                XXHAT(IIMNH-1),XXHAT(IIMNH+1)
       print *, ZYOBS,IJMNH &
               ,XXY(IJMNH,IGRID), XXY(IJMNH-1,IGRID),XYHAT(IJMNH),&
                XYHAT(IJMNH-1),XYHAT(IJMNH+1)
    endif
    IF ( YSTOCK == '3D' .OR. YSTOCK == '1D' ) THEN
      ! quelle est la maille verticale mesonh qui contient cette obs ?
      ! cas des obs a localiser sur la grille de masse
      ! XZZ_W (K) , XZZ_W(K+1) = limites Z de la maille_masse K
      ! cas des obs à localiser sur la grille de W
      ! XZZ_masse (K-1) , XZZ_masse(K) = limites Z de la maille_W K
      IKMNH=MIN(COUNT(XZZ(IIMNH,IJMNH,:)< ZOBSALTlu),NKMAX+2-1)
      IF ( IKMNH == 0 .AND. ZVALOBS /= XSPVAL ) THEN
        print *,'obs under the first model level, stored at',&
             ' k=1', ZOBSLONlu,ZOBSLATlu,ZOBSALTlu,ZVALOBS
        IKMNH=1
      ENDIF
      IF ( IGRID == 4 ) THEN  ! champ d obs sur la grille W
        IKMNH=IKMNH+1
      ENDIF
      ! stocke le niveau max pour minimiser la taille du tableau a ecrire
      IKMAX=MAX(IKMAX,IKMNH)
    ELSE
      IKMNH=1
    ENDIF
    ! stockage
    if (ilocverbia >= 2 ) then
      IKM1=MAX(IKMNH-1,1)
      print *, ZOBSLONlu,ZOBSLATlu,ZOBSALTlu,ZVALOBS,IIMNH,IJMNH,IKMNH &
             , XZZ(IIMNH,IJMNH,IKMNH),XZZ(IIMNH,IJMNH,IKM1)
    endif
    IF (ZVALOBS /= XSPVAL ) THEN
      if (ilocverbia >= 2 ) then
        print *,'before TO_COMPUTING_UNITS', ZVALOBS
      endif
      CALL TO_COMPUTING_UNITS(CGROUP,YUNITE,ZVALOBS) 
      if (ilocverbia >= 2 ) then
        print *,'after TO_COMPUTING_UNITS', ZVALOBS
      endif
      ! Voir une amelioration en moyenne ponderee avec la distance :
      ! serait utile pour des mailles tres grandes
      if (ilocverbia >=3 ) then
        print *, 'Storage indexes i,j,k=',IIMNH,IJMNH,IKMNH
      endif
      ZALTOBS(IIMNH,IJMNH,IKMNH)=ZOBSALTlu
      ZOBSinMNH(IIMNH,IJMNH,IKMNH)=ZOBSinMNH(IIMNH,IJMNH,IKMNH)+ZVALOBS
      ICPTOBSinMNH(IIMNH,IJMNH,IKMNH)=ICPTOBSinMNH(IIMNH,IJMNH,IKMNH)+1
    ENDIF
    !
    inbval=inbval+1
  ENDDO   ! fin de boucle de lecture du fichier d obs
GO TO 888
886 CONTINUE
  print *,' *** WARNING: in reading the date in the obs file ***'
  print *,' not enough rows (4)'
  GO TO 888
887 CONTINUE
  print *,' *** WARNING: in reading the obs file ***'
  print *,'             every record must contains 4 values'
  print *,'             or 3 in 2D'
888   CONTINUE
  !
  print *, 'End of reading the input obs file'
  CLOSE (UNIT=8)
  ! suppression du lien
  CALL CREATLINK('DIROBS',YFILEOBS,'CLEAN',ilocverbia)
  print *, 'number of obs taken into account in the model grid= ', inbval
  print *, 'number of obs out of domain not taken into account= ', inbvalrej
  !
  ! mise a indef des mailles MNH non concernées par les obs
  !
  WHERE ( ICPTOBSinMNH(:,:,:) == 0)  
     ZOBSinMNH(:,:,:)= XSPVAL
     ZALTOBS(:,:,:)= XSPVAL
  END WHERE
  print *, 'number of meshes set to indef= ', COUNT(ICPTOBSinMNH(:,:,:) ==0) 
  print *, 'number of meshes initialised= ', COUNT(ICPTOBSinMNH(:,:,:) > 0) 

  IF ( (COUNT (ICPTOBSinMNH(:,:,:) > 0) ) == 0 ) THEN
     print *, '**** no observation is localised into the model grid'
     print *, ' the field is not written in the output diachronic file'
     CYCLE
  ENDIF
  !
  ! calcul eventuel de la moyenne des obs incluses dans les mailles mesonh
  WHERE ( ICPTOBSinMNH(:,:,:) > 0) &
    ZOBSinMNH(:,:,:)=ZOBSinMNH(:,:,:)/ICPTOBSinMNH(:,:,:)
  print *, 'end of computation of the average on ',&
    COUNT(ICPTOBSinMNH(:,:,:) >0) , ' meshes'
  !
  ! traitement particulier des composantes du vent
  SELECT CASE (CGROUP(1:1))
    CASE ('U','V')
      IF ( .NOT. ALLOCATED (ZVENTSAVE) ) THEN
        ALLOCATE(ZVENTSAVE( SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3), &
                            SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)   ))
        ALLOCATE(ZWORK3D ( SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3) )) 
        ALLOCATE(ZWORK3D2( SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3) )) 
        print *, 'Treatment for the wind: storage of the zonal component 1'
        ZWORK3D(:,:,:)=ZOBSinMNH(:,:,:)
        print *, '  and treatement of the Obs file for the 2d component'
        GO TO 777
      ELSE
        print *, 'Treatment for the wind: storage of the meridional component 2'
        ZWORK3D2(:,:,:)=ZOBSinMNH(:,:,:)
        CALL UV_TO_ZONAL_AND_MERID(ZWORK3D,ZWORK3D2,0,      &
                                PZC=ZVENTSAVE(:,:,:,1,1,1), &
                                PMC=XVAR(:,:,:,1,1,1)       )
        
        print *,' after UV_TO_ZONAL_AND_MERID'
        DEALLOCATE( ZWORK3D,ZWORK3D2)
      ENDIF
      ! Fin traitement particulier des composantes du vent
    CASE DEFAULT
      ! init du champ  passe par module a writevar
       XVAR(:,:,:,1,1,1)=ZOBSinMNH(:,:,:)
  ENDSELECT
  ! 
  ! init des variables passees par module a writevar
  NGRIDIA(1)=IGRID
  CTITRE(1)=CGROUP
  CCOMMENT(1)='from '//ADJUSTL(YFILEOBS)
  CUNITE(1)=YUNITE
  !
  !*      3.3 Ecriture  du tableau XVAR (module MODD_ALLOC_FORDIACHRO)
  !           --------------------------------------------------
  zmini=MINVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
  zmaxi=MAXVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
  print * ,' After treatment, min,max of the field ',TRIM(CGROUP),'=', zmini,zmaxi
  print *,' Writing in diachronic format'
  if (ilocverbia >= 1 ) then
    print *,'dimensions of XVAR ', SIZE(XVAR,1) , SIZE(XVAR,2), SIZE(XVAR,3)
  endif
  !
  ivarideb=1
  ivarifin=SIZE(XVAR,1)
  ivarjdeb=1
  ivarjfin=SIZE(XVAR,2)
  ivarkdeb=1
  if ( IKMAX <= 2 ) THEN
    ivarkdeb= IKMAX
  endif
  ivarkfin=IKMAX
  ivartinf=1
  ivartsup=1
  ivartrajinf=1
  ivartrajsup=1
  ivarprocinf=1
  ivarprocsup=1
  IF ( YSTOCK == '1D' ) THEN
    ! tableaux 1D stockés pour permettre un trace diaprog en profil vertical
    ivarideb=IIMNH
    ivarifin=IIMNH
    ivarjdeb=IJMNH
    ivarjfin=IJMNH
    print * ,' Storage of 1D profil, position i,j in the grid=',ivarideb,ivarjdeb
  ENDIF
  if (ilocverbia >= 2 ) then
    print *,'before WRITEVAR',' input arguments ',&
             ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
             ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
             TRIM(CGROUP),' ',TRIM(YFILEOUTNAME),' ',TRIM(YFLAGWRITE),' ',&
             TRIM(YFILEOUTSUFFIX),&
             ilocverbia,iret
  endif
  CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
                CGROUP,YFILEOUTNAME,YFLAGWRITE,YFILEOUTSUFFIX,ilocverbia,iret)

  print *, ' WRITEVAR, return value for (',TRIM(CGROUP),')= ',iret
  IF ( iret /= 0 ) THEN 
    print *,'** Error when writing in the file: ',TRIM(YFILEOUTNAME)
    STOP
  ENDIF              
  !
  ! traitement eventuel de la 2e composante du vent
  IF ( ALLOCATED (ZVENTSAVE) ) THEN
    XVAR(:,:,:,:,:,:)= ZVENTSAVE(:,:,:,:,:,:)
    CGROUP='U'//CGROUP(2:)
    CTITRE(1)='U'//CGROUP(2:)
    if (ilocverbia >= 2 ) then
     print *,'before WRITEVAR',' input arguments ',&
             ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
             ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
             TRIM(CGROUP),' ',TRIM(YFILEOUTNAME),' ',TRIM(YFLAGWRITE),' ',&
             TRIM(YFILEOUTSUFFIX),&
             ilocverbia,iret
    endif
    CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
                CGROUP,YFILEOUTNAME,YFLAGWRITE,YFILEOUTSUFFIX,ilocverbia,iret)

     print *, ' WRITEVAR, return value for (',TRIM(CGROUP),')= ',iret
     IF ( iret /= 0 ) THEN 
       print *,'** Error when writing in the file: ',TRIM(YFILEOUTNAME)
       STOP
     ENDIF              
     DEALLOCATE (ZVENTSAVE)
  ENDIF
  !
  IF (YSTOCK == '2D' ) THEN
   !IF (COUNT(ZALTOBS(:,:,:) /= XSPVAL) /= 0) THEN
   IF (galtobs) THEN
    ! stockage egalement de l altitude des obs comme champ diachronique
    XVAR(:,:,:,1,1,1)=ZALTOBS(:,:,:)          
    NGRIDIA(1)=1
    CTITRE(1)='ALT_'//ADJUSTL(CGROUP)
    CCOMMENT(1)='from '//ADJUSTL(YFILEOBS)
    CUNITE(1)='m'          
    if (ilocverbia >= 2 ) then
      print *,'before WRITEVAR',' input arguments ',&
           ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
           ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
           TRIM(CGROUP),' ',TRIM(YFILEOUTNAME),' ',TRIM(YFLAGWRITE),' ',&
           TRIM(YFILEOUTSUFFIX),ilocverbia,iret
    endif
    CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
              ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
              CTITRE(1),YFILEOUTNAME,YFLAGWRITE,YFILEOUTSUFFIX,ilocverbia,iret)

    print *, 'WRITEVAR, return value (',TRIM(CTITRE(1)),')= ',iret
    IF ( iret /= 0 ) THEN 
      print *,'** Error when writing in the file: ',TRIM(YFILEOUTNAME)
      STOP
    ENDIF              
   ELSE
     print * , ' No altitudes in the Obs file: no field ALT_'//ADJUSTL(CGROUP)
   ENDIF
  ENDIF
  !
  ! + stockage du nombre d obs par point de grille comme champ diachronique
  XVAR(:,:,:,1,1,1)=ICPTOBSinMNH(:,:,:)          
  NGRIDIA(1)=1
  CTITRE(1)='N_'//ADJUSTL(CGROUP)
  CCOMMENT(1)='from '//ADJUSTL(YFILEOBS)
  CUNITE(1)='count'          
  if (ilocverbia >= 2 ) then
    print *,'before WRITEVAR',' input arguments ',&
         ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
         ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
         TRIM(CTITRE(1)),' ',TRIM(YFILEOUTNAME),' ',TRIM(YFLAGWRITE),' ',&
         TRIM(YFILEOUTSUFFIX),&
         ilocverbia,iret
  endif
  CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
            ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
            CTITRE(1),YFILEOUTNAME,YFLAGWRITE,YFILEOUTSUFFIX,ilocverbia,iret)

  print *, 'WRITEVAR return value (',TRIM(CTITRE(1)),')= ',iret              
  IF ( iret /= 0 ) THEN 
      print *,'** Error when writing in the file: ',TRIM(YFILEOUTNAME)
      STOP
  ENDIF            
  !
ENDDO ! boucle fichier obs a traiter
!
! Fin de boucle sur les fichiers d obs
! pour clore le traitement meme si la liste des champs est
! incomplete ( non terminee par END)
88  CONTINUE
YFILEOBS='END'
!
!---------------------------------------------------------------------------
!
!*       4.    Fermeture fichiers
!              ------------------
!
IF ( YFILEOBS(1:3) == 'END' ) THEN
  PRINT*, 'END -> Close the output file'
  YFLAGWRITE='CLO'
  ! dans cet appel seul l argument YFLAGWRITE est pris en compte, tous
  ! les autres arguments sont ignorés
  CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
               ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
               CGROUP,YFILEOUTNAME,YFLAGWRITE,YFILEOUTSUFFIX,ilocverbia,iret)
  print *, 'WRITEVAR, return value=',iret
  IF ( iret > 0 ) THEN 
    print *,'** Error when closing the file: ',TRIM(YFILEOUTNAME)
  ENDIF            
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.    Fin de programme
!              ------------------
!
99  CONTINUE
!
!   Suppression de tous les liens eventuellement crees
YDUMMYFILE=''
CALL CREATLINK(' ',YDUMMYFILE,'CLEAN',ilocverbia)
PRINT*, 'The file ',TRIM(YLUDIR),' stores all the input directives'
PRINT*, ' you must give a new name to use it again'
CLOSE(ILUDIR)
!
IF (iret==0) THEN
  print *,'================'
  IF  (YFILEOUTSUFFIX /= 'NEN' ) THEN
    PRINT*, 'Output files *',TRIM(YFILEOUTSUFFIX), '.lfi are available'  
  ELSE
    PRINT*, 'Output file ', TRIM(YFILEOUTNAME), '.lfi is available'  
  ENDIF
  PRINT*, ' Use LCOLAREA=T and LSPOT=T in diaprog to plot the fields'
ENDIF
!
END PROGRAM OBS2MESONH
