      PROGRAM  EXTRACTDIA
!     ###################
!
!!****  *EXTRACTDIA* -  lecture d'enregistrements dans fichier diachronique,
!                         traitement,
!                         ecriture (11 types de format de fichier possibles)
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
!       d'une liste de fichiers diachroniques 
!       du format de sortie
!       d'une liste de champs a traiter pour chaque fichier diachronique
!       d'un zoom selon toutes les directions inclu dans le champ a traiter
!         ( seul le zoom selon i,j,k est possible pour le format DIAC)
!
!     Ecriture en sortie:
!       d'un fichier  au format fonction de TYPEOUT c.a.d 
!         DIAC= type diachro (un seul fichier contenant toutes
!                                       les variables selectionnées)
!         LLHV= lon lat alt val (un seul fichier contenant toutes
!                                       les variables selectionnées) 
!         llhv= lat lon alt val (un seul fichier contenant toutes
!                                       les variables selectionnées) 
!         ll ou LL zv lon lat  niveau Z val
!
!         ll ou LL pv lon lat  niveau P val
!
!----------------------
! AJOUT NOVEMBRE 2009:
!----------------------
!         IJHV= i j alt val (un seul fichier contenant toutes
!                                       les variables selectionnées) 
!         jihv= j i alt val (un seul fichier contenant toutes
!                                       les variables selectionnées) 
!         IJ ou ji zv lon lat  niveau Z val
!
!         IJ ou ji pv lon lat  niveau P val
!----------------------
!         FREE= format libre a choisir par l utilisateur (un fichier par variable)
!         KCDL ou ZCDL ou PCDL= format CDL (à convertir en netcdf via "tonetcdf")
!                               (un seul fichier contenant toutes
!                                       les variables selectionnées)
!           KCDL si les niveaux verticaux sont les niveaux du modele
!           ZCDL si les niveaux verticaux sont des niveaux Z=constante donnes au programme
!           PCDL si les niveaux verticaux sont des niveaux P=constante donnes au programme
!
!  pour les formats *CDL,*Z*,*P*, 2 types de grille horizontale sont possibles:
!    'CONF' grille reguliere sur le plan de projection (conforme ou cartesien)
!    'LALO' grille reguliere en lat-lon
!             dans ce cas les composantes du vent sont transformees
!             en composantes zonales et méridiennes.
! sauf pour IJPV, IJZV, jipv, jizv :  CONF obligatoire
!!
!!    EXTERNAL
!!    --------
!!          FROM_COMPUTING_UNITS: retour aux unites initiales  avant ecriture
!!                               = passage inverse a celui realise par
!!                                 TO_COMPUTING_UNITS      
!!          appele par writevar,writecdl,writellhv 
!!              et par extractdia avant l ecriture au format FREE
!!    REFERENCE
!!    ---------
!!
!!    AUTHORS
!!    -------
!!    I. Mallet , N. Asencio, J. Stein
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/03/2003
!       call to dd and ff routines
!       call to writeLLHV if LLHV
!       clean writevar to delete choice LLHV inside this routine
!       add PCDL,LLZV,llzv,LLPV,llpv cases
!       allow a zoom 0,0,jdeb,jfin or ideb,ifin,0,0 or 0,0,0,0  05/2005
!        add ALT 3Dfield if KCDL, add the LAT and LON 3Dfields if CONF and *CDL
!       04/11/2009 (G. Tanguy) : add case IJHV,IJZV, IJPV , JIHV, JIZV, JIPV
!       29/03/2011 (G. TANGUY) : add case ZGRB PGRB
!       11/07/2014 (G. TANGUY) : correction pour les donnees LES de type SSOl
!                                (vlev et field ne correspondaient pas suite à
!                                mauvais zoom)
!       16/12/2014 (G.DELAUTIER) : ajout cas LLAV llav : altitude au dessus du
!       sol
!       18/02/2015 (G.DELAUTIER) : ajout cas AGRB : altitude au dessus du
!       sol
!       Avril 2015 (G.DELAUTIER) : ajout CFIXRESOL pour car GRIB +correction
!       pour FF10MAX
! -----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! modules MesoNH
USE MODD_CONF, ONLY: NVERB,LCARTESIAN
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT,XUNDEF,NUNDEF
USE MODD_DIM1, ONLY: NIMAX,NJMAX,NKMAX
USE MODD_GRID, ONLY: XLATORI,XLONORI
USE MODD_GRID1, ONLY: XZS,XZZ,XLAT,XLON,XXHAT,XYHAT
USE MODE_GRIDPROJ  ! subroutines SM_XYHAT et SM_LATLON 
USE MODI_UV_TO_ZONAL_AND_MERID
USE MODI_HOR_INTERP_4PTS
USE MODI_ZINTER
USE MODI_PINTER                          
! modules DIACHRO
USE MODD_FILES_DIACHRO
USE MODN_NCAR,  ONLY: XSPVAL  
USE MODD_ALLOC_FORDIACHRO, ONLY: XVAR, &         ! XVAR(i,j,k,t,n,p)
                                 XTRAJZ, &       ! XTRAJZ(k,t,n)
                                 XDATIME, &      ! XDATIME(16,t)
                                 CTITRE, CUNITE,&! CTITRE(p),CUNITE(p)
                                 NGRIDIA, & ! NGRIDIA(p)
                                 NGRID
USE MODD_COORD, ONLY: XXX,XXY,XXZS, & !  XXX(:,1:7), XXY(:,1:7), XXZS(:,:,1:7)
                      XXDXHAT,XXDYHAT ! XXDXHAT(:,1:7), XXDYHAT(:,1:7)
USE MODD_RESOLVCAR, ONLY: CGROUP, NVERBIA, &
                          NNDIA, NPROCDIA, NBPROCDIA !pour appel a interp_grids
USE MODD_TYPE_AND_LH, ONLY: NIL,NIH,NJL,NJH,NKL,NKH,CTYPE,LICP,LJCP
! modules tools
USE MODI_CHANGE_A_GRID 
USE MODI_LOW2UP 
USE MODI_CREATLINK 
USE MODI_DD
USE MODI_FF
USE MODI_WRITEDIR                                 
USE MODI_WRITELLHV
USE MODI_WRITEGRIB
USE MODI_WRITECDL                                 
USE MODI_WRITEVAR                                 
USE MODI_FROM_COMPUTING_UNITS
USE MODD_READLH
USE MODI_INI2LALO
USE MODI_INT2LALO
!                                 
IMPLICIT NONE                       
!
!*       0.1   Local variables declarations
!
INTEGER           :: I
INTEGER           :: ILUDIR,IRESP
INTEGER           :: JLOOP,JI,JJ,JK,J5,J6,J4,JA,JGR,ii
! zoom lu pour les 6 dimensions possibles
INTEGER           :: iideb,iifin,ijdeb,ijfin,ikdeb,ikfin
REAL              :: zideb,zifin,zjdeb,zjfin
INTEGER, dimension(2) :: iloc
INTEGER           :: itinf,itsup,itrajinf,itrajsup,iprocinf,iprocsup
! zoom recalcule en fonction des dimensions du champ traite
INTEGER           :: ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin
INTEGER           :: ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup
INTEGER           :: ivarzmin,ivarzmax
INTEGER           :: inbvertz,IND_VERT,IND_LL,IND_IJ
REAL , allocatable, dimension(:,:,:):: ZWORK3D,ZWORK3D2,zffvent,zdirvent
REAL , allocatable, dimension(:,:)  :: zwork2d,zwork2d2
REAL , allocatable, dimension(:,:)  :: ZLAT,ZLON
REAL , allocatable, dimension(:,:)  :: ZDIFFLON,ZDIFFLAT
! pour traiter les champs budget deja zoomes
REAL , allocatable, dimension(:,:,:,:,:,:):: ZVARSAVE                                 
! pour l interpolation verticale a P=cst : pinter
REAL , allocatable, dimension(:,:,:) :: ZPABS                          
! pour les interpolations verticales a P ou Z=cst
REAL , allocatable, dimension(:,:,:) :: ZVARZCST
REAL , allocatable, dimension(:) :: zlistevert
INTEGER :: ikdebzint ! premier niveau a traiter
!  pour l interpolation sur grille reguliere lat lon
REAL , allocatable, dimension(:,:) :: ZNEWLAT,ZNEWLON,ZNEWX,ZNEWY
REAL              :: ZDELTALAT,ZDELTALON
REAL :: zmini,zmaxi
INTEGER           :: inetadd ! compteur de champs supp dans le fichier Netcdf
INTEGER           :: IFLAGzcst,IGRID
INTEGER           :: IDIM1,IDIM2,I1,I2,IZOOMIDEB,IZOOMIFIN,IZOOMJDEB,IZOOMJFIN
INTEGER           :: IAN,IMOIS,IJOUR,IHEURE,IMINUTE,ISECONDE
! 
INTEGER           :: ilocverbia,iret,iret2,iskip,ISAVENGRIDIA,iarg,INDX,IK
CHARACTER(LEN=3)  :: YK
!                    flag pour initialiser/ne pas initialiser le zoom d
!                      d ecriture : 
!                      ne pas initialiser quand ajout par le programme
!                      des champs ALT LAT LON qui doivent conserver le
!                      zoom de l utilisateur
INTEGER           :: ino_init_zoom
! **** la taille des variables caracteres contenant les noms
!      de fichiers est obligatoirement de 28 ****
CHARACTER(LEN=28) :: YFILEIN,YFILEOUT
! **** la longueur du nom ne doit pas depasser 13 car. si le fichier
! contient des groupes a un seul PROCessus, ou 9 si plusieurs PROCessus ****
CHARACTER(LEN=13) :: YGROUP,YGROUP_OLD
CHARACTER(LEN=20) :: YGROUP_SAVE
CHARACTER(LEN=4)  :: YTYPEOUT
CHARACTER(LEN=1)  :: YTYPEOUT3
CHARACTER(LEN=3)  :: YSUFFIX_file
CHARACTER(LEN=250):: YFMTFREE   ! format ecriture des champs si YTYPEOUT='FREE'
CHARACTER(LEN=45) :: YFILEOUTFREE ! nom du fichier de sortie si YTYPEOUT='FREE'
CHARACTER(LEN=5)  :: YFLAGREADVAR ,YFLAGWRITE
CHARACTER(LEN=4)  :: YOUTGRID  ! grille en sortie:
                  !CONF pour rester dans le plan conforme,
                  ! (le logiciel graphique devra réaliser la projection)
                  !LALO pour passer à lat,lon réguliers
CHARACTER(LEN=28) :: YDUMMYFILE
CHARACTER(LEN=11) :: YLUDIR      !  Name of the dir file
REAL   , DIMENSION(:,:)  ,ALLOCATABLE        :: ZX,ZY             
! GRIB
INTEGER :: IND_GRB
INTEGER :: ICODCOD ! Parameter grib code
INTEGER :: ICODLEV ! grib code for the Type of Level 
INTEGER :: ICODOLL ! bottom level if layer 
INTEGER :: ICODOLH ! level or top of level if layer
CHARACTER(LEN=256):: YINPLINE    ! input agregation line read from Namelist
LOGICAL :: LVAR2D
INTEGER :: ILEVEL2D ! en option : altitude du champ 2D à coder dans le fichier GRIB
LOGICAL :: LLEVEL2D 
REAL,DIMENSION(4) :: ZLATLON
INTEGER,DIMENSION(4) :: ILATLON
INTEGER :: INX,INY
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZALT
REAL,DIMENSION(:,:,:),ALLOCATABLE :: zlistevert3D
INTEGER :: IZLIST
CHARACTER(LEN=1) :: CFIXRESOL
REAL :: ZDX_GRB,ZDY_GRB,ZCONTROL
!-------------------------------------------------------------------------------
!
!*       1.     INIT
!               ----
!
!
inetadd=0  !compteur de champs supp dans le fichier Netcdf
!
!Prints : 0=mini 1=debug mode in extractdia, readvar and writevar , writecdl, writellhv
!                3=debug mode in routines diachro'
! nverbia= controle des prints dans les routines diachro
ilocverbia=0
! 
! dans mesonh Xundef est utilise =999.
! dans les routines diachro XSPVAL est utilisé                  
XSPVAL=XUNDEF                                    
!
! ouverture d un fichier dir ou vont s ecrire les entrees clavier
YLUDIR='dirextract'
CALL FMATTR(YLUDIR,YLUDIR,ILUDIR,IRESP)
OPEN(UNIT=ILUDIR,FILE=YLUDIR,FORM='FORMATTED')
!
! Possibilite de definir un zoom d ecriture 
!  definition locale du zoom pour extractdia et writevar, writecdl, writellhv
iideb=0
iifin=0
ijdeb=0
ijfin=0
ikdeb=0
ikfin=0
itinf=0
itsup=0
itrajinf=0
itrajsup=0
iprocinf=0
iprocsup=0
!
!-------------------------------------------------------------------------------
!
!*       2.     INPUT FILE AND FORMAT
!               ---------------------
!
!*       2.1   name of file and output format
!              ------------------------------
!
PRINT*, '- Name of the diachro file (without .lfi) ?'
READ(5,'(A28)') YFILEIN
CALL WRITEDIR(ILUDIR,YFILEIN)
!
PRINT*, '- type of the output file ?'
PRINT*, '(DIAC/llhv/llzv/llpv/llav/LLHV/LLZV/LLPV/LLAV/IJHV/IJZV/IJPV/jihv/jizv/jipv/FREE/KCDL/ZCDL/PCDL/ZGRB/PGRB/AGRB)'
READ(5,'(A4)')YTYPEOUT
CALL WRITEDIR(ILUDIR,YTYPEOUT)
PRINT*,'the file ',TRIM(YFILEIN),' will be converted in type ',YTYPEOUT
!
PRINT*, '- Prints : 0=mini 1=debug mode in extractdia'
PRINT*, '                  3=debug mode in routines diachro'
PRINT*, '?'
READ(5,*)ilocverbia
CALL WRITEDIR(ILUDIR,ilocverbia)
PRINT*, ' output prints= ',ilocverbia
if ( ilocverbia > 2) nverbia=ilocverbia   ! verbosity of diachro routines
NVERB=ilocverbia                          ! verbosity of mesonh routines
!
!*       2.2   other parameters
!              ----------------
!
SELECT CASE (YTYPEOUT)                                   
  CASE('LLHV','llhv','DIAC','FREE','KCDL','ZCDL','PCDL','llzv','LLZV',&
          &'llpv','LLPV','IJHV','IJZV','IJPV','jihv','jizv','jipv','ZGRB','PGRB','AGRB',&
          &'llav','LLAV') ! lecture des choix de l utilisateur
    IF ( YTYPEOUT == 'FREE' ) THEN
      PRINT*, '- format of writing for fields ? '
      PRINT*, '    (fortran syntaxe of FMT in WRITE)'
      PRINT*,'exemple: (10F9.3) or (8F0.3)'
      PRINT*, '?'
      READ(5,'(A)') YFMTFREE
      CALL WRITEDIR(ILUDIR,YFMTFREE)
      PRINT*, ' format=', TRIM(YFMTFREE)
    ENDIF
    ! lecture du zoom
    IND_VERT= INDEX(YTYPEOUT(1:4),'Z') + INDEX(YTYPEOUT(1:4),'P') + &
              INDEX(YTYPEOUT(1:4),'z') + INDEX(YTYPEOUT(1:4),'p') + &
              INDEX(YTYPEOUT(1:4),'a') + INDEX(YTYPEOUT(1:4),'A')
    IND_LL= INDEX(YTYPEOUT(1:2),'L') + INDEX(YTYPEOUT(1:2),'l') 
    IND_IJ= INDEX(YTYPEOUT(1:2),'IJ') + INDEX(YTYPEOUT(1:2),'ji') 
    IND_GRB=INDEX(YTYPEOUT(1:4),'GRB')
print*,YTYPEOUT,IND_IJ
    IF (IND_LL==0 .AND. IND_GRB==0) THEN
      IF (IND_VERT/=0) THEN
        ! cas 'ZCDL','PCDL','jizv','jipv','IJZV','IJPV'
        PRINT*, '- zoom on the 2 first dimensions: '
        PRINT*, '              ideb,ifin,jdeb,jfin'
        PRINT*, '0,0,0,0 for the whole physical domain'
        PRINT*, '-1,-1,-1,-1 for the whole domain'
        PRINT*, '?'
        READ(5,*) iideb,iifin,ijdeb,ijfin
        CALL WRITEDIR(ILUDIR,iideb)
        CALL WRITEDIR(ILUDIR,iifin)
        CALL WRITEDIR(ILUDIR,ijdeb)
        CALL WRITEDIR(ILUDIR,ijfin)
      ELSE 
        ! cas 'DIAC','FREE','KCDL','IJHV','jihv'
        PRINT*, '- zoom on the 3 first dimensions: '
        PRINT*, '              ideb,ifin,jdeb,jfin,kdeb,kfin'
        PRINT*, '0,0,0,0,0,0 for the whole physical domain'
        PRINT*, '-1,-1,-1,-1,-1,-1 for the whole domain'
        PRINT*, '?'
        READ(5,*) iideb,iifin,ijdeb,ijfin,ikdeb,ikfin
        CALL WRITEDIR(ILUDIR,iideb)
        CALL WRITEDIR(ILUDIR,iifin)
        CALL WRITEDIR(ILUDIR,ijdeb)
        CALL WRITEDIR(ILUDIR,ijfin)
        CALL WRITEDIR(ILUDIR,ikdeb)
        CALL WRITEDIR(ILUDIR,ikfin)
      END IF
    ELSE
      ! cas 'llzv','LLZV','llpv','LLPV','llhv','LLHV','llav' 'LLAV'
      !      'ZGRB','PGRB','AGRB'
      PRINT*, '- zoom on the 2 first directions: '
      PRINT*, '              lonmin,lonmax,latmin,latmax'
      PRINT*, '0.,0.,0.,0. for the whole physical domain'
      PRINT*, '-1.,-1.,-1.,-1. for the whole domain'
      PRINT*, '?'
      READ(5,*) zideb,zifin,zjdeb,zjfin
      CALL WRITEDIR(ILUDIR,zideb)
      CALL WRITEDIR(ILUDIR,zifin)
      CALL WRITEDIR(ILUDIR,zjdeb)
      CALL WRITEDIR(ILUDIR,zjfin)
      if(zideb==0. .AND. zifin==0.) then
        iideb=0 ; iifin=0
      else if(zideb==-1. .AND. zifin==-1.) then
        iideb=-1 ; iifin=-1
      else
        iideb=-2 ; iifin=-2
      endif
      if(zjdeb==0. .AND. zjfin==0.) then
        ijdeb=0 ; ijfin=0
      else if(zjdeb==-1. .AND. zjfin==-1.) then
        ijdeb=-1 ; ijfin=-1
      else
        ijdeb=-2 ; ijfin=-2
      endif
      IF (IND_GRB/=0) THEN
        PRINT*,'Do you want to fix resolution in x and y ? (y/n)'
        PRINT*,'(only available with LALO)'
        READ(5,*) CFIXRESOL
        CALL WRITEDIR(ILUDIR,CFIXRESOL)
        IF (CFIXRESOL=='y') THEN
          PRINT*,'Enter x resolution (in millidegrees)'
          READ(5,*) ZDX_GRB
          PRINT*,'Enter y resolution (in millidegrees)'
          READ(5,*) ZDY_GRB
          CALL WRITEDIR(ILUDIR,ZDX_GRB)
          CALL WRITEDIR(ILUDIR,ZDY_GRB)
        ENDIF

      ENDIF
      IF (IND_VERT==0) THEN
        ! cas 'llhv','LLHV'
        PRINT*, '- zoom on the 3rd dimension: '
        PRINT*, '                 kdeb,kfin'
        PRINT*, '0,0 for the whole physical domain'
        PRINT*, '-1,-1 for the whole domain'
        PRINT*, '?'
        READ(5,*) ikdeb,ikfin
        CALL WRITEDIR(ILUDIR,ikdeb)
        CALL WRITEDIR(ILUDIR,ikfin)
      END IF
    END IF
    PRINT*, '- zoom on the 3 last dimensions : '
    PRINT*, '   itinf,itsup,itrajinf,itrajsup,iprocinf,iprocsup'
    PRINT*, '0,0,0,0,0,0 for the whole last dimensions'
    PRINT*, '?'
    READ(5,*) itinf,itsup,itrajinf,itrajsup,iprocinf,iprocsup
    CALL WRITEDIR(ILUDIR,itinf)
    CALL WRITEDIR(ILUDIR,itsup)
    CALL WRITEDIR(ILUDIR,itrajinf)
    CALL WRITEDIR(ILUDIR,itrajsup)
    CALL WRITEDIR(ILUDIR,iprocinf)
    CALL WRITEDIR(ILUDIR,iprocsup)
    IF ((iideb==-2) .AND. (ijdeb==-2)) THEN
      PRINT'(A6,4(E10.4,X),2(I4,X),2(I5,X),4(I4,X))', ' zoom= ',zideb,zifin,zjdeb,zjfin,ikdeb,ikfin&
                      ,itinf,itsup,itrajinf,itrajsup,iprocinf,iprocsup
    ELSE
      PRINT'(A6,6(I4,X),2(I5,X),4(I4,X))', ' zoom= ',iideb,iifin,ijdeb,ijfin,ikdeb,ikfin&
                      ,itinf,itsup,itrajinf,itrajsup,iprocinf,iprocsup
    END IF
    IF (IND_VERT/=0) THEN
      PRINT*, '- Number of vertical levels for ',YTYPEOUT(IND_VERT:IND_VERT),' interpolation ?'
      READ(5,*) inbvertz
      CALL WRITEDIR(ILUDIR,inbvertz)
      PRINT*, '- Ordered list of these levels (in meters or in hPa): exemple 500 1500 ?'
      allocate (zlistevert(inbvertz))
      READ(5,*) zlistevert
      DO JI=1,inbvertz
        CALL WRITEDIR(ILUDIR,zlistevert(JI))
      END DO
      PRINT*, ' interpolation for the following ',YTYPEOUT(IND_VERT:IND_VERT),' levels='
      PRINT*, zlistevert 
    ENDIF
    YOUTGRID='CONF'
    IF (YTYPEOUT/='DIAC' .AND. YTYPEOUT/='llhv' .AND. YTYPEOUT/='LLHV' .AND.&
   & IND_IJ==0) THEN
      PRINT *,'- Fields in regular LAt/LOn grid'
      PRINT *,'  or    in regular grid on CONFormal plan (native MesoNH grid) ?'
      PRINT *,'LALO/CONF ?'
      READ(5,*) YOUTGRID
      CALL WRITEDIR(ILUDIR,YOUTGRID)
      PRINT*, ' Output grid= ', YOUTGRID
      PRINT*, ''
      YSUFFIX_file=YTYPEOUT(1:2)//YTYPEOUT(4:4)
      IF ( YTYPEOUT(2:4) == 'CDL') THEN
        PRINT*, '!!!!!!!! Warning !!!!!!!!'
        PRINT*, 'For the CDL type, the dimensions are initialised'
        PRINT*, ' with those of the first field:'
        PRINT*, 'the values of the 6 dimensions must be the maximum that'
        PRINT*, ' will be treated '
        PRINT*, '!!!!!!!! Warning !!!!!!!!'
        PRINT*, 'For the CDL type, the coordinates must be the same'
        PRINT*, ' for all fields'
        PRINT*, '(stored in the output file with LAT/LON/VLEV groups)'
        PRINT*, '!!!!!!!!'
      ENDIF
    ELSE IF (IND_IJ/=0) THEN ! dans le cas des points de grille on prend les
                             !  coordonnees conformes
      YOUTGRID='CONF'
    ENDIF
  CASE DEFAULT
    PRINT*, 'Incorrect value for the output type:',YTYPEOUT
    PRINT*, 'the following ones are currently available :'
    PRINT*, 'DIAC,LLHV,llhv,FREE,KCDL,ZCDL,PCDL,llzv,LLZV,llpv,LLPV,llav,LLAV'
    PRINT*, 'IJHV,IJZV,IJPV,jihv,jizv,jipv,ZGRB,PGRB,AGRB'
    STOP
END SELECT
! 
!*       2.3   init for input file and output file
!              -----------------------------------
! in READVAR, input file must be opened before reading
YFLAGREADVAR='OPE'
! in WRITE routine, output file is new
YFLAGWRITE='NEW'
! 
!*       2.4   lecture de la pression pour interpolation
!              -----------------------------------------
IF (INDEX(YTYPEOUT(1:4),'p')/=0 .OR. INDEX(YTYPEOUT(1:4),'P')/=0 )THEN
  CALL READVAR('PABST',YFILEIN,YFLAGREADVAR,ilocverbia,iret)
  IF ( iret /= 0 ) then
    print *, '- PABST not found, name of the pressure variable ? '
    read *,YGROUP
    CALL WRITEDIR(ILUDIR,YGROUP)
    CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
    IF ( iret /= 0 ) then
      print *,' interpolation at P=cst not possible because PABST and ',TRIM(YGROUP),' are not available'
      STOP
    ENDIF
  ENDIF
  ! stockage de ZPABS utilise par pinter
  ALLOCATE ( ZPABS(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3))) 
  ZPABS(:,:,:)=XVAR(:,:,:,1,1,1) 
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    LOOP ON GROUPS IN THE FILE
!              --------------------------
!
DO JGR=1,10000
  !    
  !*      3.0  preparation pour la lecture du champ suivant
  !
  ino_init_zoom=0
  IF (IND_GRB==0) THEN
    PRINT*,'- Name of the group in upper case (13 characters max.)'
    PRINT*,' (ex: THT or DD or FF or DD10 or FF10 or LAT or LON or VLEV)'
    PRINT*,'(GROUP for the list of groups, END to stop)?'
    READ(5,'(A13)',END=88) CGROUP
    CALL WRITEDIR(ILUDIR,CGROUP)
    CGROUP=ADJUSTL(CGROUP)
    CALL LOW2UP(CGROUP)
  ELSE ! CASE ZGRB or PGRB
    LLEVEL2D=.FALSE.      
    LVAR2D=.FALSE. 
    PRINT*,'- Name of the group in upper case (13 characters max.)'
    PRINT*,' MesoNH field name, grib parameter indicator'
    PRINT*,' (ex: UT 131, VT 132, GROUP for the list of groups, END to stop)'
    PRINT*,' optional : you can add FOR 2D FIELDS ONLY the altitude (in meters)'
    PRINT*,' of the field after  the grib parameter indicator exple : UT10 131 10'
    READ(5,'(A)') YINPLINE
    YINPLINE= TRIM(ADJUSTL(YINPLINE))
    IF (LEN_TRIM(YINPLINE) == 0) CYCLE ! skip blank line
    CALL WRITEDIR(ILUDIR,YINPLINE)
    CALL TAB2SPACE(YINPLINE)
    ! extract field name
    INDX= INDEX(YINPLINE,' ')
    CGROUP= YINPLINE(1:INDX-1)
    IF (CGROUP=='END') GO TO 88
    ! 
    IF (CGROUP /='GROUP') THEN
      ICODLEV=NUNDEF
      ICODOLH=NUNDEF
      ICODOLL=NUNDEF
      YINPLINE= ADJUSTL(YINPLINE(INDX:))
      INDX= INDEX(YINPLINE,' ')
      IF (INDX == 1 ) THEN
        PRINT*, ' Parameter indicator is missing. ',CGROUP,' not treated.'
        CYCLE
      END IF
      READ(YINPLINE(1:INDX-1),*) ICODCOD
      IF (NVERB>=5) print*, ' Parameter indicator: ',ICODCOD
      YINPLINE= ADJUSTL(YINPLINE(INDX:))
      INDX= INDEX(YINPLINE,' ')
      IF (INDX /= 1 ) THEN
        READ(YINPLINE(1:INDX-1),*) ILEVEL2D     
        PRINT*, 'Level found : ',ILEVEL2D
        PRINT*, 'it will be only used if the field ',CGROUP,' is 2D'
        LLEVEL2D=.TRUE.
      END IF

    ENDIF
  ENDIF
  IF (CGROUP=='END') GO TO 88
  ! point de reprise pour forcer l ecriture des champs VLEV,LAT,LON 
  ! dans les fichiers netcdf
77 CONTINUE
  YGROUP_SAVE=CGROUP(1:13)
  YK=''
  INDX=INDEX(CGROUP,'_K_')
  IF (INDX/=0) THEN
    CGROUP=YGROUP_SAVE(1:INDX-1)
    YK(1:3)=YGROUP_SAVE(INDX+3:INDX+5)
    READ(YK,'(I3)') IK
  END IF
  IF (CGROUP(1:5)/='GROUP') &
    PRINT*,'you asked for the following record: ',TRIM(CGROUP)
  !
  !*      3.1  Lecture et initialisation du tableau XVAR
  !            passé en module MODD_ALLOC_FORDIACHRO
  !
  !
  !      3.1.1 Cas particulier pour le vent
  !
  IF ( CGROUP(1:2) == 'UT' .OR. &
       CGROUP(1:2) == 'VT' .OR. &
       CGROUP(1:2) == 'DD' .OR. &
       CGROUP(1:2) == 'FF' .AND. CGROUP(1:7) /= 'FF10MAX'     )  THEN
    !
    IF ( (CGROUP(1:2)=='UT'.OR.CGROUP(1:2)=='VT') .AND. &
          YOUTGRID(1:4) /= 'LALO'                       ) THEN
      ! Lecture du champ U ou V sans calcul 
      ! les composantes du vent restent dans le plan conforme
      CALL READVAR(CGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
    ELSE
      ! Lecture des 2 composantes du vent  : commence par UM...
      !(stockees dans les tableaux ZWORK3D et ZWORK3D2)
      ! max 13 car.
      YGROUP='UT'//CGROUP(3:13) 
      CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
      IF ( iret /= 0 ) then
        print *,TRIM(CGROUP),': ',TRIM(YGROUP),' not available'
        ! echec , on tente UM....
        YGROUP='UM'//CGROUP(3:13)
        CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret2)
        IF ( iret2 /= 0 ) then
          print *,'** no processing for ',TRIM(CGROUP), &
                  ' because UT and ',TRIM(YGROUP),' are not available'
          CYCLE
        ENDIF
      ENDIF
      ! allocation du tableau de stockage de la 1e composante du vent
      ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
                        size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
      ZVARSAVE=XVAR
      !
      ! deuxieme composante VT....
      YGROUP='VT'//CGROUP(3:13)
      CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
      IF ( iret /= 0 ) then
        print *,TRIM(CGROUP),': ',TRIM(YGROUP),' not available'
        ! echec , on tente VM....
        YGROUP='VM'//CGROUP(3:13)
   CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret2)
        IF ( iret2 /= 0 ) then
          print *,'** no processing for ',TRIM(CGROUP), &
                  ' because VT and ',TRIM(YGROUP),' are not available'
          CYCLE
        ENDIF
        iret=iret2
      ENDIF
      !
      ! Calcul de ff
      IF (CGROUP(1:2) == 'FF' ) THEN
        IF (LEN(TRIM(CGROUP)) ==2) THEN
          YGROUP='VENTFF'
        ELSE IF (LEN(TRIM(CGROUP)) ==3) THEN
          YGROUP='VENT'//CGROUP(3:3)//'FF'
        ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
          YGROUP='VENT'//CGROUP(3:4)//'FF'
        ELSE
          ! 13 car max
          YGROUP='VENTFF'//CGROUP(3:9)
        ENDIF
        ! allocation du tableau de calcul
        IF (allocated(ZWORK3D)) DEALLOCATE(ZWORK3D)
        ALLOCATE(ZWORK3D(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
        ZWORK3D(:,:,:)=XSPVAL
        DO J6=1,SIZE(XVAR,6)
          IGRID=NGRIDIA(J6)
          DO J5=1,SIZE(XVAR,5)
          DO J4=1,SIZE(XVAR,4)
            CALL FF (ZVARSAVE(:,:,:,J4,J5,J6),XVAR(:,:,:,J4,J5,J6),ZWORK3D, &
                     JPVEXT,JPHEXT,IGRID)
            XVAR(:,:,:,J4,J5,J6)=ZWORK3D(:,:,:)
          END DO
          END DO
          ! initialisation des variables necessaires a l ecriture
          CGROUP=YGROUP
          CTITRE(J6)=YGROUP
          NGRIDIA(J6)=1
        END DO
        DEALLOCATE(ZWORK3D)
        ! Calcul de dd par rapport au Nord geographique
      ELSE IF (CGROUP(1:2) == 'DD') THEN
        IF (CTYPE=='CART' .OR. CTYPE=='MASK' .OR. CTYPE=='SPXY') THEN 
          IF (LEN(TRIM(CGROUP)) ==2) THEN
            YGROUP='VENTDD'
          ELSE IF (LEN(TRIM(CGROUP)) ==3) THEN
            YGROUP='VENT'//CGROUP(3:3)//'DD'
          ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
            YGROUP='VENT'//CGROUP(3:4)//'DD'
          ELSE
          ! 13 car max
            YGROUP='VENTDD'//CGROUP(3:9) 
          ENDIF
          ! allocation du tableau de calcul
          IF (allocated(ZWORK3D)) DEALLOCATE(ZWORK3D)
          ALLOCATE(ZWORK3D(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
          DO J6=1,SIZE(XVAR,6)
            IGRID=NGRIDIA(J6)
            DO J5=1,SIZE(XVAR,5)
            DO J4=1,SIZE(XVAR,4)
              iskip=1 ! tous les points de grille
              CALL DD(ZVARSAVE(:,:,:,J4,J5,J6),XVAR(:,:,:,J4,J5,J6),ZWORK3D, &
                      iskip,IGRID,PLON=XLON(NIL:NIH,NJL:NJH))
              XVAR(:,:,:,J4,J5,J6)=ZWORK3D(:,:,:)
            END DO
            END DO
            ! initialisation des variables necessaires a l ecriture
            CGROUP=YGROUP
            CTITRE(J6)=YGROUP
            CUNITE(J6)='degrees'
            NGRIDIA(J6)=1
          END DO
          DEALLOCATE(ZWORK3D)
        ELSE
          print *,'** processing of ',TRIM(CGROUP),' is not performed for CTYPE= ',CTYPE
          CYCLE
        ENDIF
      ELSE IF (CGROUP(1:2) == 'UT' .OR. CGROUP(1:2) == 'VT') THEN
        IF (CTYPE=='CART' .OR. CTYPE=='MASK' .OR. CTYPE=='SPXY') THEN 
        ! Calcul des composantes zonale et meridienne
        !(YOUTGRID(1:4) == 'LALO') avec la routine UV_TO_ZONAL_AND_MERID
          print*,' Translate to meridional and zonal wind components'
          ALLOCATE(ZWORK3D(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
          ALLOCATE(ZWORK3D2(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
          IF (ilocverbia >= 3 ) then
            print *,'before UV_TO_ZONAL_AND_MERID KGRID=23'
            print'(A31,3(I5,X))',' dimensions of the input arrays',size(ZVARSAVE,1),&
                                      size(ZVARSAVE,2),size(ZVARSAVE,3)
            print'(3(I5,X))',size(XVAR,1),size(XVAR,2),size(XVAR,3)
            print'(A32,3(I5,X))',' dimensions of the output arrays',size(ZWORK3D,1),&
                                       size(ZWORK3D,2),size(ZWORK3D,3)
            print'(3(I5,X))',size(ZWORK3D2,1),size(ZWORK3D2,2),size(ZWORK3D2,3)
          ENDIF
          DO J6=1,SIZE(XVAR,6)
            DO J5=1,SIZE(XVAR,5)
            DO J4=1,SIZE(XVAR,4)
              CALL UV_TO_ZONAL_AND_MERID(ZVARSAVE(:,:,:,J4,J5,J6), &
                                         XVAR(:,:,:,J4,J5,J6),     &
                                         23,PZC=ZWORK3D,PMC=ZWORK3D2)
              IF (CGROUP(1:1) == 'U' ) THEN
                XVAR(:,:,:,J4,J5,J6)=ZWORK3D(:,:,:)
              ENDIF
              IF (CGROUP(1:1) == 'V' ) THEN
                XVAR(:,:,:,J4,J5,J6)=ZWORK3D2(:,:,:)
              ENDIF
            END DO
            END DO
          END DO
          IF (ilocverbia >= 3 ) then
            print *,'after UV_TO_ZONAL_AND_MERID KGRID=23'
          END IF
          ! Stockage dans le tableau XVAR qui est le tableau ecrit
          ! de la composante souhaitée
          IF (CGROUP(1:1) == 'U' ) THEN
            print *, ' U zonal wind component'
            IF (LEN(TRIM(CGROUP)) ==2) THEN
              YGROUP='UZON'
            ELSE IF (LEN(TRIM(CGROUP)) ==3) THEN
              YGROUP='U'//CGROUP(3:3)//'ZON'
            ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
              YGROUP='U'//CGROUP(3:4)//'ZON'
            ELSE
              ! 13 car max
              YGROUP='UZON'//CGROUP(3:9)
            ENDIF
            CTITRE(:)='U zonal wind component'
          ELSE IF (CGROUP(1:1) == 'V' ) THEN
            print *, ' V meridian wind component'
            IF (LEN(TRIM(CGROUP)) ==2) THEN
              YGROUP='VMED'
            ELSE IF (LEN(TRIM(CGROUP)) ==3) THEN
              YGROUP='V'//CGROUP(3:3)//'MED'  
            ELSE IF (LEN(TRIM(CGROUP)) ==4) THEN
              YGROUP='V'//CGROUP(3:4)//'MED'
            ELSE
              ! 13 car max
              YGROUP='VZON'//CGROUP(3:9)
            END IF
            CTITRE(:)='V meridian wind component'
          ENDIF
          CGROUP=YGROUP
          NGRIDIA(:)=1  ! UZON et VMED en grille de masse
          DEALLOCATE(ZWORK3D,ZWORK3D2)
        ELSE
          print *,' No processing of UZON and VMED for CTYPE= ',CTYPE
          CYCLE
        ENDIF
      ENDIF
      DEALLOCATE(ZVARSAVE)
    ENDIF
  !
  !      3.1.2 LATitude ou LONgitude de chaque point de la grille conforme
  !
  ELSE IF (CGROUP(1:3)=='LAT' .OR. CGROUP(1:3)=='LON') THEN
    print *, 'LAT/LON asked and YFLAGREADVAR=', YFLAGREADVAR
   IF ( YFLAGREADVAR /= 'NOP') THEN
    ! Lecture d un champ 2D quelconque pour initialiser XLAT et XLON
    CALL READVAR('ZSBIS',YFILEIN,YFLAGREADVAR,ilocverbia,iret)
    IF ( iret /= 0 ) then
     ! cas de fichier diachronique sans ZSBIS
      print *, '- Name of one group in upper case '
      read *,YGROUP
      CALL WRITEDIR(ILUDIR,YGROUP)
      CALL LOW2UP(YGROUP)
      CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
      IF ( iret /= 0 ) then
         print * ,'**group ', TRIM(YGROUP) , 'not found'
         stop
      ENDIF
    ENDIF
   ENDIF
    ! init du tableau XVAR au champ souhaite
    DEALLOCATE(XVAR)
    ALLOCATE(XVAR(size(XLAT,1),size(XLAT,2),1,1,1,1) )
    IF (CGROUP(1:3)=='LAT') THEN
      XVAR(:,:,1,1,1,1)=XLAT(:,:)
      CTITRE(1)='latitudes'
      CUNITE(1)='degrees_north'
    ELSE IF (CGROUP(1:3)=='LON') THEN
      XVAR(:,:,1,1,1,1)=XLON(:,:)
      CTITRE(1)='longitudes'
      CUNITE(1)='degrees_east'
    ENDIF
  !
  !      3.1.3 Altitude de chaque point de la grille conforme
  !
  ELSE IF (CGROUP(1:4)=='VLEV') THEN
    print *, 'VLEV asked and YFLAGREADVAR=', YFLAGREADVAR
    IF(CTYPE=='SSOL'.OR.CTYPE=='DRST'.OR.CTYPE=='RAPL'.OR.CTYPE=='RSPL') THEN
      IF ( YFLAGREADVAR == 'NOP') THEN
      ! altitude des niveaux du groupe precedent dans XTRAJZ
        print *,'warning, for CTYPE=',CTYPE,' Vertical LEVels of previous group (',TRIM(YGROUP_OLD),')'
        DEALLOCATE(XVAR)
        ALLOCATE(XVAR(1,1,size(XTRAJZ,1),1,1,1))
        XVAR(1,1,:,1,1,1)=XTRAJZ(:,1,1)
      ELSE
        print*,'** no processing with VLEV at the first group'
        GOTO 99
      ENDIF
    ELSE
      IF ( YFLAGREADVAR /= 'NOP') THEN
        ! Lecture d un champ 2D quelconque pour initialiser les tableaux XZZ
        CALL READVAR('ZSBIS',YFILEIN,YFLAGREADVAR,ilocverbia,iret)
        IF ( iret /= 0 ) then
          ! cas de fichier diachronique sans ZSBIS
          print *, '- Name of one group in upper case '
          read *,YGROUP
          CALL WRITEDIR(ILUDIR,YGROUP)
          CALL LOW2UP(YGROUP)
          CALL READVAR(YGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
          IF ( iret /= 0 ) then
            print * ,'** group ', TRIM(YGROUP) , 'not found'
            stop
          ENDIF
        ENDIF
      ENDIF
      ! init de XZZ a la grille de masse ( par defaut readvar 
      ! l initialise a la grille 4 des  vitesse verticales W)
      CALL COMPCOORD_FORDIACHRO(1)
      ! init du tableau XVAR au champ souhaite
      DEALLOCATE(XVAR)
      ALLOCATE(XVAR(size(XZZ,1),size(XZZ,2),size(XZZ,3),1,1,1))
      XVAR(:,:,:,1,1,1)=XZZ(:,:,:)
      ! retour au XZZ grille 4
      CALL COMPCOORD_FORDIACHRO(4)
    ENDIF
    CTITRE(1)='model levels altitudes ASL'
    CUNITE(1)='meters'
  !
  !      3.1.4 Default case
  !
  ELSE
    !
    ! Lecture du  champ CGROUP et stockage dans XVAR
    ! + Initialisation (si YFLAGREADVAR='OPE') des variables
    ! des modules (cf USE en debut de programme)
    ! Appel a menu_diachro pour la liste des groupes si CGROUP(1:5)=='GROUP'
    CALL READVAR(CGROUP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
    IF (CGROUP(1:5)=='GROUP') CYCLE
    !
  ENDIF
  !
  IF ( iret == 0 ) THEN
    zmini=MINVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
    zmaxi=MAXVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
    print * ,' After read, min,max of the variable ',TRIM(CGROUP),'=', zmini,zmaxi
    !       
    !*      3.2  Init des bornes min max du zoom en fonction des
    !            dimensions du tableau XVAR traite
    !
   IF ( ino_init_zoom == 0) THEN
    IF (iideb == 0 .AND. iifin == 0 ) THEN
      ivarideb=NREADIL ; ivarifin=NREADIH
      IF (ivarideb/=ivarifin) THEN  ! domI/=1
        ivarideb=MAX(1+JPHEXT,NREADIL) 
        ivarifin=MIN(SIZE(XVAR,1)-JPHEXT,NREADIH)
      ENDIF
    ELSE IF (iideb == -1 .AND. iifin == -1 ) THEN
      ivarideb=MAX(1,NREADIL) 
      ivarifin=MIN(SIZE(XVAR,1),NREADIH)
    ELSE IF (iideb == -2 .AND. iifin == -2 ) THEN
      ivarideb=-2
      iideb=1+JPHEXT
      IF (zideb >= minval(XLON)) THEN
        DO JJ=1,SIZE(XLON,2)
          ivarideb=MAX(MIN(COUNT(XLON(:,JJ)<zideb),SIZE(XLON,1)),iideb)
          iideb=ivarideb
        END DO
      ENDIF
      ivarifin=-2
      iifin=1+JPHEXT
      IF (zifin <= maxval(XLON)) THEN
        DO JJ=1,SIZE(XLON,2)
          ivarifin=MAX(MIN(COUNT(XLON(:,JJ)<zifin),SIZE(XLON,1)),iifin)
          iifin=ivarifin
        END DO
      ENDIF
    ELSE
      ivarideb=max(iideb,NREADIL)
      ivarifin=min(iifin,NREADIH)
      ivarideb=min(ivarideb,ivarifin)
    ENDIF
    IF(ijdeb == 0 .AND. ijfin == 0) THEN
      ivarjdeb=NREADJL ; ivarjfin=NREADJH
      IF (ivarjdeb/=ivarjfin) THEN  ! domJ/=1
        ivarjdeb=MAX(1+JPHEXT,NREADJL)
        ivarjfin=MIN(SIZE(XVAR,2)-JPHEXT,NREADJH)
      ENDIF
    ELSE IF (ijdeb == -1 .AND. ijfin == -1 ) THEN
      ivarjdeb=MAX(1,NREADJL)
      ivarjfin=MIN(SIZE(XVAR,2),NREADJH)
    ELSE IF (ijdeb == -2 .AND. ijfin == -2 ) THEN
      ivarjdeb=-2
      ijdeb=1+JPHEXT
      IF (zjdeb >= minval(XLAT)) THEN
        DO JI=1,SIZE(XLAT,1)
          ivarjdeb=MAX(MIN(COUNT(XLAT(JI,:)<zjdeb),SIZE(XLAT,2)),ijdeb)
          ijdeb=ivarjdeb
        END DO
      ENDIF
      ivarjfin=-2
      ijfin=1+JPHEXT
      IF (zjfin <= maxval(XLAT)) THEN
        DO JI=1,SIZE(XLAT,1)
          ivarjfin=MAX(MIN(COUNT(XLAT(JI,:)<zjfin),SIZE(XLAT,2)),ijfin)
          ijfin=ivarjfin
        END DO
      ENDIF
    ELSE
      ivarjdeb=max(ijdeb,NREADJL)
      ivarjfin=min(ijfin,NREADJH)
      ivarjdeb=min(ivarjdeb,ivarjfin)
    ENDIF
    IF(ivarideb==-2 .OR. ivarifin==-2 .OR. ivarjdeb==-2 .OR.  ivarjfin==-2) THEN
      print *,'****zoom provided is not included in the FM-file grid'
      print *,'LON (zoom: ',zideb,zifin,') (file: ',minval(XLON),maxval(XLON)
      print *,'LAT (zoom: ',zjdeb,zjfin,') (file: ',minval(XLAT),maxval(XLAT)
      GOTO 99
    ENDIF
    IF (IND_VERT/=0) THEN
      ivarzmin=1   ; ivarzmax=inbvertz
    ELSE
      ivarzmin=MAX(1,NREADKL)  ; ivarzmax=MIN(SIZE(XVAR,3),NREADKH)
      inbvertz=ivarzmax-ivarzmin+1
    ENDIF
    IF (ikdeb == 0 .AND. ikfin == 0 ) THEN
      ivarkdeb=NREADKL ; ivarkfin=NREADKH
      IF (ivarkdeb/=ivarkfin .AND. CTYPE/='SSOL') THEN  ! domK/=1
        ivarkdeb=MAX(1+JPVEXT,NREADKL)
        ivarkfin=min(ivarzmax,SIZE(XVAR,3)-JPVEXT)
      ENDIF
    ELSEIF (ikdeb == -1 .AND. ikfin ==-1 ) THEN
      ivarkdeb=ivarzmin
      ivarkfin=ivarzmax
    ELSE
      ivarkdeb=max(ikdeb,ivarzmin)
      ivarkfin=min(ikfin,ivarzmax)
      ivarkdeb=min(ivarkdeb,ivarkfin)
    ENDIF   
    IF (INDX/=0) THEN
      ivarkdeb=IK ; ivarkfin=IK
    END IF
   ENDIF

    IF (itinf == 0 .AND. itsup == 0 ) THEN
      ivartinf=1 ; ivartsup=SIZE(XVAR,4)
    ELSE
      ivartinf=max(itinf,1)
      ivartsup=min(itsup,SIZE(XVAR,4))
      ivartinf=min(ivartinf,ivartsup)
    ENDIF
    IF (itrajinf == 0 .AND. itrajsup == 0 ) THEN
      ivartrajinf=1 ; ivartrajsup=SIZE(XVAR,5)
    ELSE
      ivartrajinf=max(itrajinf,1)
      ivartrajsup=min(itrajsup,SIZE(XVAR,5))
      ivartrajinf=min(ivartrajinf,ivartrajsup)
    ENDIF
    IF (iprocinf == 0 .AND. iprocsup == 0 ) THEN
      ivarprocinf=1 ; ivarprocsup=SIZE(XVAR,6)
    ELSE
      ivarprocinf=max(iprocinf,1)
      ivarprocsup=min(iprocsup,SIZE(XVAR,6))
      ivarprocinf=min(ivarprocinf,ivarprocsup)
    ENDIF
    if (ilocverbia > 0 ) then
      PRINT*,' Zoom limits initialized with:'
      PRINT'(A53,6(I4,X))','ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin',&
            ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin 
      PRINT'(A53,6(I4,X))','ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocfin',&
            ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup 
    endif
    !
    !*      3.3  Ecriture  du tableau XVAR (module MODD_ALLOC_FORDIACHRO) 
    !
    print *,' Write with the format ', YTYPEOUT(1:4)
    SELECT CASE(YTYPEOUT(1:4))
      !
      CASE('DIAC')
        CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                      ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,  &
                      CGROUP,YFILEIN,YFLAGWRITE,'2  ',ilocverbia,iret)
        if (ilocverbia > 0 ) then
          print'(A17,I2))','WRITEVAR return= ',iret
        end if
      !
      CASE('FREE')
        if (ilocverbia >= 0 ) then
          print*,' format ',YTYPEOUT
          print'(A53,X,A50,6(I4,X),2(I6,X),4(I4,X))',&
                 ' domaine for writting : ideb,ifin,jdeb,jfin,kdeb,kfin', &
                 ',itinf,itsup,itrajinf,itrajsup,iprocinf,iprocsup= ', &
              ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
              ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup 
        endif
        !  Retour aux unites initiales si necessaire
        CALL FROM_COMPUTING_UNITS(CGROUP,CUNITE(1)) 
        !
        YFILEOUTFREE=ADJUSTL(ADJUSTR(YFILEIN)//'.'//ADJUSTL(ADJUSTR(CGROUP)))
        OPEN (UNIT=7,STATUS='NEW',FORM='FORMATTED',FILE=YFILEOUTFREE)
        ! a. Ecriture de l entete
        !temps courant
        IAN=XDATIME(13,1)
        IMOIS=XDATIME(14,1)
        IJOUR=XDATIME(15,1)
        IHEURE=XDATIME(16,1)/3600
        IMINUTE=(XDATIME(16,1)-(IHEURE*3600))/60
        ISECONDE=ISECONDE-(IHEURE*3600)-(IMINUTE*60)
        WRITE(7,FMT='(6(I4,X),2(I6,X),4(I4,X),4(I4,X),A42,A33)') ivarideb,&
                   ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                   ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                   IAN,IMOIS,IJOUR,IHEURE,IMINUTE ,&
                   'format ligne1=  12 Indices (.deb .fin) du ',&
                   'tableau  an mois jour hUTC minute'
        ! b. ecriture des données au fmt choisi par l utilisateur
        WRITE(7,FMT=YFMTFREE) &
         XVAR(ivarideb:ivarifin,ivarjdeb:ivarjfin,ivarkdeb:ivarkfin,&
              ivartinf:ivartsup,ivartrajinf:ivartrajsup,ivarprocinf:ivarprocsup)
        PRINT*,'File ',TRIM(YFILEOUTFREE),' available'
        CLOSE(7)
      !
      CASE('LLHV','llhv','IJHV','jihv')
        IF (CTYPE == 'SSOL') THEN
          ALLOCATE(ZALT(1,1,SIZE(XTRAJZ,1)))
          ZALT(1,1,:)=XTRAJZ(:,1,1)
          CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,&
                       ilocverbia,iret,PALT=ZALT)
        ELSE
          CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,&
                       ilocverbia,iret)     
        ENDIF  
        if (ilocverbia > 0 ) then
          print*,' WRITELLHV return= ',iret
        end if
      !
      CASE('KCDL','ZCDL','PCDL','LLZV','LLPV','llpv','llzv',&
             & 'IJZV','jizv','IJPV','jipv','llav','LLAV')
        ! replace field at mass points
        IF ( CGROUP /= 'VLEV' ) THEN
          If (ALLOCATED(ZWORK3D))DEALLOCATE(ZWORK3D)
          If (ALLOCATED(ZWORK3D2))DEALLOCATE(ZWORK3D2)
          ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
          ALLOCATE(ZWORK3D2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
          DO J6=ivarprocinf,ivarprocsup
            IGRID=NGRIDIA(J6)
            IF(SIZE(XVAR,3)/=1 .OR. IGRID/=4) THEN 
              ! pas d interpolation verticale pour champ 2D
              DO J5=ivartrajinf,ivartrajsup
                DO J4=ivartinf,ivartsup
                  ZWORK3D(:,:,:)=XVAR(:,:,:,J4,J5,J6)
                  print'(A29,3(X,I4))',' mass point grid for J4,J5,J6=',J4,J5,J6
                  CALL CHANGE_A_GRID(ZWORK3D,IGRID,ZWORK3D2)
                  NGRIDIA(J6)=IGRID
                  ! IGRID=1 en sortie de change_a_grid
                  XVAR(:,:,:,J4,J5,J6)=ZWORK3D2(:,:,:)               
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(ZWORK3D,ZWORK3D2)
        ENDIF
        !
        ! a. reinit avant ecriture de la grille verticale correspondant a la
        !grille de masse sur laquelle le champ a ete interpole
        IFLAGzcst=0
        IF (IND_VERT/=0) THEN
          IF ( CGROUP == 'VLEV' ) THEN
            ! ecriture de la liste des niveaux verticaux 
            IFLAGzcst=1
            DEALLOCATE(XVAR)
            allocate(XVAR(1,1,inbvertz,1,1,1))
            XVAR(1,1,:,1,1,1)=zlistevert
            ivarideb=1 ; ivarifin=1
            ivarjdeb=1 ; ivarjfin=1
            ivarkdeb=1 ; ivarkfin=inbvertz
            CTITRE(1)='vertical_levels'
            CUNITE(1)='user choice'
            IF ( YTYPEOUT(IND_VERT:IND_VERT) == 'z' .OR.  YTYPEOUT(IND_VERT:IND_VERT) == 'Z' ) THEN
                CUNITE(1)='km'
                XVAR=XVAR*0.001
            ENDIF
            IF ( YTYPEOUT(IND_VERT:IND_VERT) == 'p' .OR.  YTYPEOUT(IND_VERT:IND_VERT) == 'P' ) THEN
              CUNITE(1)='hPa'
            ENDIF
          ENDIF
        ! b. interpolation eventuelle selon la verticale
          IF( SIZE(XVAR,3)>1 .AND. CGROUP /= 'VLEV' ) THEN
            ! VLEV, LON, LAT et chps 2D ne passent pas cette partie 
            if (ilocverbia >= 0 ) then
              print*,' Interpolations on ',inbvertz,' ', &
                     YTYPEOUT(IND_VERT:IND_VERT),'-levels'
            endif
            if (ilocverbia >= 1 .AND. IND_VERT/=0) THEN
              print*,'levels= ',zlistevert 
            endif
            ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
                              size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
            ZVARSAVE=XVAR
            ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
            ALLOCATE(ZVARZCST(SIZE(XVAR,1),SIZE(XVAR,2),inbvertz))
            DEALLOCATE(XVAR)
            ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),SIZE(ZVARZCST,3),&
                          size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)))
            DO J6=ivarprocinf,ivarprocsup
              IGRID=NGRIDIA(J6)
              ! init du tableau des altitudes  XZZ pour la grille= IGRID
              CALL COMPCOORD_FORDIACHRO(IGRID)
              DO J5=ivartrajinf,ivartrajsup
                DO J4=ivartinf,ivartsup
                  ZWORK3D(:,:,:)=ZVARSAVE(:,:,:,J4,J5,J6)
                  ikdebzint=2
                  IF (INDEX(YTYPEOUT(1:4),'Z')/=0 .OR. INDEX(YTYPEOUT(1:4),'z')/=0) THEN
                    CALL ZINTER(ZWORK3D,XZZ,ZVARZCST,zlistevert,ikdebzint,XSPVAL)
                  ELSE IF (INDEX(YTYPEOUT(1:4),'A')/=0 .OR. INDEX(YTYPEOUT(1:4),'a')/=0) THEN
                    IF (.NOT. ALLOCATED(zlistevert3D)) THEN
                      ALLOCATE(zlistevert3D(SIZE(XZS,1),SIZE(XZS,2),SIZE(zlistevert)))
                      DO IZLIST=1,SIZE(zlistevert)
                        zlistevert3D(:,:,IZLIST)=XZS(:,:)+zlistevert(IZLIST)
                      ENDDO
                    ENDIF     
                    CALL SINTER(ZWORK3D,XZZ,ZVARZCST,zlistevert3D,ikdebzint,XSPVAL)
                  ELSE IF (INDEX(YTYPEOUT(1:4),'P')/=0 .OR. INDEX(YTYPEOUT(1:4),'p')/=0) THEN
                    CALL PINTER(ZWORK3D,0,XSPVAL,zlistevert,ZVARZCST,ZPABS)
                  ELSE IF (INDEX(YTYPEOUT(1:4),'H')/=0 .OR. INDEX(YTYPEOUT(1:4),'h')/=0) THEN
                    ZVARZCST(:,:,:)=ZWORK3D(:,:,:)
                  ELSE
                    print*,'** ERROR in vertical interpolations with ',YTYPEOUT
                  ENDIF
                  XVAR(:,:,:,J4,J5,J6)=ZVARZCST
                END DO
              END DO
            END DO
            DEALLOCATE(ZVARSAVE,ZVARZCST,ZWORK3D)
            zmini=MINVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
            zmaxi=MAXVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
            print * ,' After vertical interpolation, min,max of the variable ',TRIM(CGROUP),'=', zmini,zmaxi
            ivarkdeb=1
            ivarkfin=inbvertz
            IF (ilocverbia >= 5 ) then
              print*,'ivarkdeb,ivarkfin= ',ivarkdeb,ivarkfin 
            ENDIF
          ENDIF
        ENDIF
        ! c. interpolation eventuelle sur l horizontale
        IF ( YOUTGRID(1:4) == 'LALO' ) THEN
          if (ilocverbia >= 0 ) then
            print *,'Translate to a regular lat lon grid '
          end if
          IF ( .NOT. ALLOCATED (ZNEWX) ) THEN
            IF ( IFLAGzcst == 1 ) THEN
              print*,'** no processing with VLEV at the first group'
              GOTO 99
            ELSE
            ! c.1. creation de la grille réguliere en lat lon
              if (ilocverbia >= 2 ) then
                print *,'grid creation, size of XLON: ',SIZE(XLON,1),SIZE(XLON,2) 
              end if
              ! calcul des coord X Y des points de la grille lat-lon reguliere
              ! determine le maximum d espacement en lat et lon sur le zoom
              ALLOCATE(ZDIFFLON(SIZE(XLON,1)-1,SIZE(XLON,2)))
              ALLOCATE(ZDIFFLAT(SIZE(XLAT,1),SIZE(XLAT,2)-1))
              
              DO ii=1,SIZE(XLON,1)-1
                DO jj=1,SIZE(XLON,2)
                     ZDIFFLON(ii,jj)=XLON(ii+1,jj)-XLON(ii,jj)
                END DO
              END DO

              DO ii=1,SIZE(XLAT,1)
                DO jj=1,SIZE(XLAT,2)-1
                     ZDIFFLAT(ii,jj)=XLAT(ii,jj+1)-XLAT(ii,jj)
                END DO
              END DO
!              ZDELTALON=NINT(maxval(ZDIFFLON)*1000.)
!              ZDELTALAT=NINT(maxval(ZDIFFLAT)*1000.)
                           ZDELTALON=maxval(ZDIFFLON)
              ZDELTALAT=maxval(ZDIFFLAT)
              DEALLOCATE(ZDIFFLON)
              DEALLOCATE(ZDIFFLAT)
              if (ZDELTALON == 0 .OR. ZDELTALAT == 0 ) THEN
                print *,' error during ZDELTALON,ZDELTALAT computation=', ZDELTALON,ZDELTALAT
                print *,'XLON(ivarideb+1,ivarjdeb)-XLON(ivarideb,ivarjdeb)'&
                        ,'XLON(ivarifin,ivarjfin)-XLON(ivarifin-1,ivarjfin)'&
                        ,'XLAT(ivarideb,ivarjdeb+1)-XLAT(ivarideb,ivarjdeb)'&
                        ,'XLAT(ivarifin,ivarjfin)-XLAT(ivarifin,ivarjfin-1)'
                print *,XLON(ivarideb+1,ivarjdeb)-XLON(ivarideb,ivarjdeb)&
                       ,XLON(ivarifin,ivarjfin)-XLON(ivarifin-1,ivarjfin)&
                       ,XLAT(ivarideb,ivarjdeb+1)-XLAT(ivarideb,ivarjdeb)&
                       ,XLAT(ivarifin,ivarjfin)-XLAT(ivarifin,ivarjfin-1)           
                print *, 'ivarideb+1,ivarjdeb,ivarifin-1,ivarjfin',ivarideb+1,ivarjdeb,ivarifin-1,ivarjfin
                print *,'Verify the fields LAT LON of the FM file'
                ALLOCATE(ZX(SIZE(XLAT,1),SIZE(XLAT,2)),ZY(SIZE(XLAT,1),SIZE(XLAT,2)))
                ZX(1:SIZE(XZZ,1),1) = XXX(1:SIZE(XZZ,1),IGRID)
                ZX(:,2:SIZE(XZZ,2)) = SPREAD(ZX(:,1),2,SIZE(XZZ,2)-1)
                ZY(1,1:SIZE(XZZ,2)) = XXY(1:SIZE(XZZ,2),IGRID)
                ZY(2:SIZE(XZZ,1),:) = SPREAD(ZY(1,:),1,SIZE(XZZ,1)-1)
                CALL SM_LATLON(XLATORI,XLONORI,ZX,ZY,XLAT,XLON)
                ZDELTALON=max(XLON(ivarideb+1,ivarjdeb)-XLON(ivarideb,ivarjdeb)&
                           ,XLON(ivarifin,ivarjfin)-XLON(ivarifin-1,ivarjfin))
                ZDELTALAT=max(XLAT(ivarideb,ivarjdeb+1)-XLAT(ivarideb,ivarjdeb)&
                           ,XLAT(ivarifin,ivarjfin)-XLAT(ivarifin,ivarjfin-1))
                print *,' After Model Grid computation: ZDELTALON,ZDELTALAT=', ZDELTALON,ZDELTALAT
              endif
              IDIM1=(maxval(XLON)-minval(XLON))/ZDELTALON
              IDIM2=(maxval(XLAT)-minval(XLAT))/ZDELTALAT
              ALLOCATE (ZNEWLAT(IDIM1,IDIM2),ZNEWLON(IDIM1,IDIM2) )
              if (ilocverbia >= 1 ) then
                print*,' ZDELTALON,ZDELTALAT= ',ZDELTALON,ZDELTALAT
              endif
              if (ilocverbia >= 2 ) then
                print*,' IDIM1,IDIM2= ',IDIM1,IDIM2
              endif
              ! depart de la nouvelle grille : coin Sud Ouest
              DO JI=1,IDIM1
                ZNEWLON(JI,:)=minval(XLON) + (JI-1) *ZDELTALON
              ENDDO
              DO JJ=1,IDIM2
                ZNEWLAT(:,JJ)=minval(XLAT) + (JJ-1) *ZDELTALAT
              ENDDO
              if (ilocverbia >= 4 ) then
                print*, 'new lat lon grid=',ZNEWLAT(1,:)
                print*, ZNEWLON(:,1)
              endif

              ALLOCATE (ZNEWX(IDIM1,IDIM2))
              ALLOCATE (ZNEWY(IDIM1,IDIM2))
              CALL SM_XYHAT(XLATORI,XLONORI,ZNEWLAT,ZNEWLON,ZNEWX,ZNEWY)
              if (ilocverbia >= 4 ) then
                ! XXX= XXHAT et XXY=XYHAT pour les 7 grilles
                print*,' After SM_XYHAT old limits X: ', &
                       XXX(1,IGRID),XXX(SIZE(XVAR,1),IGRID)
                print*,'                new limits X: ', &
                       ZNEWX(1,1),ZNEWX(IDIM1,IDIM2)
                print*,'                old limits Y: ', &
                       XXY(1,IGRID),XXY(SIZE(XVAR,2),IGRID)
                print*,'                new limits Y: ', &
                       ZNEWY(1,1),ZNEWY(IDIM1,IDIM2)
              endif
              if (ilocverbia >= 5 ) then
                DO JI= 1,SIZE(XVAR,1) 
                  print*,'XXHAT ZNEWX',XXX(JI,IGRID),ZNEWX(JI,1),ZNEWX(JI,IDIM2)
                ENDDO
                DO JJ= 1,SIZE(XVAR,2) 
                  print*,'XYHAT ZNEWY',XXY(JJ,IGRID),ZNEWY(1,JJ),ZNEWX(IDIM1,JJ)
                ENDDO
              endif
              ! calcul de la section de tableau correspondant au zoom

!===================================================================================================================
              I1=(maxval(XLON(ivarideb:ivarifin,ivarjdeb:ivarjfin)) &
                 -minval(XLON(ivarideb:ivarifin,ivarjdeb:ivarjfin)) )/ZDELTALON
              I2=(maxval(XLAT(ivarideb:ivarifin,ivarjdeb:ivarjfin)) &
                 -minval(XLAT(ivarideb:ivarifin,ivarjdeb:ivarjfin)) )/ZDELTALAT
              IZOOMIDEB=MAX(MIN(COUNT(ZNEWLON(:,1)<minval(XLON(ivarideb:ivarifin,ivarjdeb:ivarjfin))),IDIM1),1)
              IZOOMJDEB=MAX(MIN(COUNT(ZNEWLAT(1,:)<minval(XLAT(ivarideb:ivarifin,ivarjdeb:ivarjfin))),IDIM2),1)
              IZOOMIFIN=MAX(MIN(COUNT(ZNEWLON(:,1)<maxval(XLON(ivarideb:ivarifin,ivarjdeb:ivarjfin))),IDIM1),1)
              IZOOMJFIN=MAX(MIN(COUNT(ZNEWLAT(1,:)<maxval(XLAT(ivarideb:ivarifin,ivarjdeb:ivarjfin))),IDIM2),1)
!=====================================================================================================================


              if (ilocverbia >= 2 ) then
                print*,' ZOOM along i in the LON-LAT grid: ', &
                       IZOOMIDEB,IZOOMIFIN,I1
                print*,'            j                    : ', &
                       IZOOMJDEB,IZOOMJFIN,I2
              endif
            ENDIF
          ENDIF ! fin grille ZNEWX deja allouee
          ! c.2. interpolation sur la nouvelle grille
          IF( IFLAGzcst/= 1 .AND. (NREADIH-NREADIL)>0 .AND. (NREADJH-NREADJL)>0 )THEN
            ! interpolation vers la nouvelle grille réguliere en lat lon 
            !sauf la grille verticale definie en niveaux Z et champs 1D
            if (ilocverbia >= 1 ) then
              print*,' interpolation for the variable  ',trim(CGROUP)
            end if
            allocate(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
            allocate(ZWORK3D2(IDIM1,IDIM2,SIZE(XVAR,3)))
            ! stockage des champs interpoles dans la nouvelle grille
            if (allocated (ZVARSAVE)) DEALLOCATE(ZVARSAVE)
            allocate(ZVARSAVE(IDIM1,IDIM2,SIZE(XVAR,3),&
                     SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)))
            ! boucle sur les dimensions 4 5 6
            DO J6=ivarprocinf,ivarprocsup
              DO J5=ivartrajinf,ivartrajsup
                DO J4=ivartinf,ivartsup
                  ZWORK3D(:,:,:)=XVAR(:,:,:,J4,J5,J6)
                  if (ilocverbia >= 2 ) then
                    print *,'before HOR_INTERP_4PTS J4,J5,J6=', J4,J5,J6
                  end if
                  CALL HOR_INTERP_4PTS(XXX(:,IGRID),XXY(:,IGRID),ZWORK3D, &
                                       ZNEWX,ZNEWY,ZWORK3D2)
                  ZVARSAVE(:,:,:,J4,J5,J6)=ZWORK3D2(:,:,:)
                END DO
              END DO
            END DO
            ! resultat dans XVAR passe en module
            DEALLOCATE (XVAR)
            ALLOCATE(XVAR(IDIM1,IDIM2,SIZE(ZVARSAVE,3),&
                    SIZE(ZVARSAVE,4),SIZE(ZVARSAVE,5),SIZE(ZVARSAVE,6)))
            XVAR=XSPVAL
            XVAR(:,:,:,ivartinf:ivartsup,ivartrajinf:ivartrajsup,ivarprocinf:ivarprocsup)= &
            ZVARSAVE(:,:,:,ivartinf:ivartsup,ivartrajinf:ivartrajsup,ivarprocinf:ivarprocsup)
            DEALLOCATE (ZVARSAVE)
            zmini=MINVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
            zmaxi=MAXVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
            print * ,' After horizontal interpolation, min,max of the variable ',TRIM(CGROUP),'=', zmini,zmaxi
            if (ilocverbia >= 2 ) then
              print*, 'After HOR_INTERP_4PTS all the dim 4,5,6'
            endif
            IF (allocated(ZWORK3D)) DEALLOCATE(ZWORK3D)
            IF (allocated(ZWORK3D2)) DEALLOCATE(ZWORK3D2)
          ENDIF
        ENDIF
        ! d. ecriture des donnees au format cdl ou llz/llp
        IF ( YOUTGRID(1:4) == 'LALO' ) THEN
          IF ( IFLAGzcst /= 1 ) THEN
            ivarideb=IZOOMIDEB
            ivarifin=IZOOMIFIN
            ivarjdeb=IZOOMJDEB
            ivarjfin=IZOOMJFIN
          ENDIF
          SELECT CASE(YTYPEOUT(1:4))
          CASE('LLZV','llzv','LLPV','llpv','LLAV','llav')
            IF (allocated(ZWORK3D)) DEALLOCATE(ZWORK3D)
            ALLOCATE(ZWORK3D(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
            IF (SIZE(XVAR,3)==inbvertz) THEN
              ZWORK3D(1,1,:)=zlistevert
            ELSE
              ZWORK3D(1,1,:)=XSPVAL
            ENDIF
            CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,&
                       ilocverbia,iret,PLON=ZNEWLON,PLAT=ZNEWLAT,&
                       PALT=ZWORK3D)       
            if (ilocverbia > 0 ) then
              print*,'WRITELLHV LALO return= ', YTYPEOUT,'= ',iret
            end if
            DEALLOCATE(ZWORK3D)
        !
          CASE('KCDL','ZCDL','PCDL')
           YGROUP=ADJUSTL(ADJUSTR(CGROUP)//ADJUSTL(YK))
           CALL WRITECDL(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                        ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                        ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
                        YGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YSUFFIX_file, &
                        ilocverbia,iret,PGRIDX=ZNEWLON(:,1),PGRIDY=ZNEWLAT(1,:))
           IF (ilocverbia >= 1 ) print *,' counter of added fields=',inetadd
           if ( inetadd == 0) then
            IF( SIZE(XZZ,3)<=1 .OR. SIZE(XZZ,2)<=1 .OR. SIZE(XZZ,1)<=1 ) THEN
            ! VLEV, LON, LAT et chps 2D ne passent pas cette partie 
             print *,' *****The program could not add the VLEV 3Dfield to the netcdf file****'
            ELSE
             print *,' The program adds the VLEV 3Dfield to the netcdf file'
             YGROUP_OLD=CGROUP(1:13)
             CGROUP='VLEV'
             inetadd=inetadd+1
             YFLAGWRITE='OLD'
             ino_init_zoom=1
             GO TO 77
            ENDIF
           endif
           if ( inetadd == 1 .AND. YOUTGRID(1:4) == 'CONF' )THEN
             print *,' The program adds the LAT 3Dfield to the netcdf file'      
             CGROUP='LAT'
             inetadd=inetadd+1
             YFLAGWRITE='OLD'
             ino_init_zoom=1
             GO TO 77
           endif
           if ( inetadd == 2 .AND. YOUTGRID(1:4) == 'CONF' )THEN
             print *,' The program adds the LON 3Dfield to the netcdf file'      
             CGROUP='LON'
             inetadd=inetadd+1
             YFLAGWRITE='OLD'
             ino_init_zoom=1
             GO TO 77
           endif
           
          END SELECT
        ELSE ! pas d interpolation horizontale
          SELECT CASE(YTYPEOUT(1:4))
          CASE('LLZV','llzv','LLPV','llpv','IJZV','jizv','IJPV','jipv','LLAV','llav')
            IF (SIZE(XVAR,3)==inbvertz) THEN  ! champ 3D
              IF (allocated(ZWORK3D)) DEALLOCATE(ZWORK3D)
              ALLOCATE(ZWORK3D(size(XVAR,1),size(XVAR,2),size(XVAR,3)))
              ZWORK3D(1,1,:)=zlistevert
              CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,&
                       ilocverbia,iret,&
                       PALT=ZWORK3D)       
            ELSE                              ! champ 2D
              IF((YTYPEOUT(3:3)=='z').OR.(YTYPEOUT(3:3)=='p')) YTYPEOUT3='h'
              IF((YTYPEOUT(3:3)=='Z').OR.(YTYPEOUT(3:3)=='P')) YTYPEOUT3='H'
              IF((YTYPEOUT(3:3)=='a').OR.(YTYPEOUT(3:3)=='A')) YTYPEOUT3='H'
              CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
                       CGROUP,YFILEIN,YFLAGWRITE, &
                       YTYPEOUT(1:2)//YTYPEOUT3//YTYPEOUT(4:4), &
                       ilocverbia,iret)
            ENDIF
            if (ilocverbia > 0 ) then
              print*,' WRITELLHV for ', YTYPEOUT,', return value= ',iret
            end if
            IF (allocated(ZWORK3D)) DEALLOCATE(ZWORK3D)
        !
          CASE('KCDL','ZCDL','PCDL')
            YGROUP=ADJUSTL(ADJUSTR(CGROUP)//ADJUSTL(YK))
            CALL WRITECDL(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                        ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                        ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
                        YGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YSUFFIX_file, &
                        ilocverbia,iret,PGRIDX=XXX(:,IGRID),PGRIDY=XXY(:,IGRID))
            IF (ilocverbia >= 1 ) print *,' counter of added fields=',inetadd
           if ( inetadd == 0) then
             if (ivarkdeb == ivarkfin .AND. ivarkdeb == 1 ) THEN
               print *, 'No VLEV field for only one vertical position'
             else
             print *,' The program adds the VLEV 3Dfield to the netcdf file'
             YGROUP_OLD=CGROUP(1:13)
             CGROUP='VLEV'
             inetadd=inetadd+1
             YFLAGWRITE='OLD'
             ino_init_zoom=1
             GO TO 77
             endif
           endif
           if ( inetadd == 1 .AND. YOUTGRID(1:4) == 'CONF' )THEN
             if (ivarideb /= ivarifin .AND. ivarjdeb /= ivarjfin ) THEN

             print *,' The program adds the LAT 3Dfield to the netcdf file'      
             CGROUP='LAT'
             inetadd=inetadd+1
             ino_init_zoom=1
             GO TO 77
             else
               print *, ' No LAT field for only one location', ivarideb,ivarifin,ivarjdeb,ivarjfin
             endif
           endif
           if ( inetadd == 2 .AND. YOUTGRID(1:4) == 'CONF' )THEN
             if (ivarideb /= ivarifin .AND. ivarjdeb /= ivarjfin ) THEN      
             print *,' The program adds the LON 3Dfield to the netcdf file'      
             CGROUP='LON'
             inetadd=inetadd+1
             ino_init_zoom=1
             GO TO 77
             else
               print *, ' No LON field for only one location', ivarideb,ivarifin,ivarjdeb,ivarjfin
             endif
           endif
          END SELECT              
        ENDIF
        ! retour a XZZ pour NGRID a 4 (cf readvar)
        CALL COMPCOORD_FORDIACHRO(4)
!============================================   
   CASE('ZGRB','PGRB','AGRB')
        IF(SIZE(XVAR,3)==1) THEN
           LVAR2D=.TRUE.
        ENDIF
        ! replace field at mass points
        IF ( CGROUP /= 'VLEV' ) THEN
          If (ALLOCATED(ZWORK3D))DEALLOCATE(ZWORK3D)
          If (ALLOCATED(ZWORK3D2))DEALLOCATE(ZWORK3D2)
          ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
          ALLOCATE(ZWORK3D2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
          DO J6=ivarprocinf,ivarprocsup
            IGRID=NGRIDIA(J6)
            IF(SIZE(XVAR,3)/=1 .OR. IGRID/=4) THEN 
              ! pas d interpolation verticale pour champ 2D
              DO J5=ivartrajinf,ivartrajsup
                DO J4=ivartinf,ivartsup
                  ZWORK3D(:,:,:)=XVAR(:,:,:,J4,J5,J6)
                  print'(A29,3(X,I4))',' mass point grid for J4,J5,J6=',J4,J5,J6
                  CALL CHANGE_A_GRID(ZWORK3D,IGRID,ZWORK3D2)
                  ! IGRID=1 en sortie de change_a_grid
                  NGRIDIA(J6)=IGRID
                  XVAR(:,:,:,J4,J5,J6)=ZWORK3D2(:,:,:)               
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(ZWORK3D,ZWORK3D2)
        ENDIF
        !
        ! a. reinit avant ecriture de la grille verticale correspondant a la
        !grille de masse sur laquelle le champ a ete interpole
        IFLAGzcst=0
        IF (IND_VERT/=0) THEN
          ! b. interpolation eventuelle selon la verticale 
          IF( SIZE(XVAR,3)>1 .AND. SIZE(XVAR,2)>1 .AND. SIZE(XVAR,1)>1 ) THEN
            ! VLEV, LON, LAT et chps 2D ne passent pas cette partie 
            if (ilocverbia >= 0 ) then
              print*,' Interpolations on ',inbvertz,' ', &
                     YTYPEOUT(IND_VERT:IND_VERT),'-levels'
            endif
            if (ilocverbia >= 1 .AND. IND_VERT/=0) THEN
              print*,'levels= ',zlistevert 
            endif
            ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
                              size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
            ZVARSAVE=XVAR
            ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
            ALLOCATE(ZVARZCST(SIZE(XVAR,1),SIZE(XVAR,2),inbvertz))
            DEALLOCATE(XVAR)
            ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),SIZE(ZVARZCST,3),&
                          size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)))
            DO J6=ivarprocinf,ivarprocsup
              IGRID=NGRIDIA(J6)
              ! init du tableau des altitudes  XZZ pour la grille= IGRID
              CALL COMPCOORD_FORDIACHRO(IGRID)
              DO J5=ivartrajinf,ivartrajsup
                DO J4=ivartinf,ivartsup
                  ZWORK3D(:,:,:)=ZVARSAVE(:,:,:,J4,J5,J6)
                  ikdebzint=2
                  IF (INDEX(YTYPEOUT(1:4),'Z')/=0 .OR. INDEX(YTYPEOUT(1:4),'z')/=0) THEN
                    CALL ZINTER(ZWORK3D,XZZ,ZVARZCST,zlistevert,ikdebzint,XSPVAL)
                  ELSE IF (INDEX(YTYPEOUT(1:4),'A')/=0 .OR. INDEX(YTYPEOUT(1:4),'a')/=0) THEN
                    IF (.NOT. ALLOCATED(zlistevert3D)) THEN
                      ALLOCATE(zlistevert3D(SIZE(XZS,1),SIZE(XZS,2),SIZE(zlistevert)))
                      DO IZLIST=1,SIZE(zlistevert)
                        zlistevert3D(:,:,IZLIST)=XZS(:,:)+zlistevert(IZLIST)
                      ENDDO
                    ENDIF     
                    CALL SINTER(ZWORK3D,XZZ,ZVARZCST,zlistevert3D,ikdebzint,XSPVAL)
                  ELSE IF (INDEX(YTYPEOUT(1:4),'P')/=0 .OR. INDEX(YTYPEOUT(1:4),'p')/=0) THEN
                    CALL PINTER(ZWORK3D,0,XSPVAL,zlistevert,ZVARZCST,ZPABS)
                  ELSE IF (INDEX(YTYPEOUT(1:4),'H')/=0 .OR. INDEX(YTYPEOUT(1:4),'h')/=0) THEN
                    ZVARZCST(:,:,:)=ZWORK3D(:,:,:)
                  ELSE
                    print*,'** ERROR in vertical interpolations with ',YTYPEOUT
                  ENDIF
                  XVAR(:,:,:,J4,J5,J6)=ZVARZCST
                END DO
              END DO
            END DO
            DEALLOCATE(ZVARSAVE,ZVARZCST,ZWORK3D)
            zmini=MINVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
            zmaxi=MAXVAL(XVAR(:,:,:,:,:,:),MASK=XVAR(:,:,:,:,:,:)/=XSPVAL)
            print * ,' After vertical interpolation, min,max of the variable ',TRIM(CGROUP),'=', zmini,zmaxi
            ivarkdeb=1
            ivarkfin=inbvertz
            IF (ilocverbia >= 5 ) then
              print*,'ivarkdeb,ivarkfin= ',ivarkdeb,ivarkfin 
            ENDIF
          ENDIF
        ENDIF
        ! c. interpolation eventuelle sur l horizontale
        IF ( YOUTGRID(1:4) == 'LALO' ) THEN
          ZLATLON(1)=MAXVAL(XLAT)*1000.
          ZLATLON(2)=MINVAL(XLAT)*1000.
          ZLATLON(3)=MINVAL(XLON)*1000.
          ZLATLON(4)=MAXVAL(XLON)*1000.
      
          IF (ZJFIN /=0 .AND. ZJFIN/=-1) ZLATLON(1)=zjfin*1000.
          IF (ZJDEB /=0 .AND. ZJDEB/=-1) ZLATLON(2)=zjdeb*1000.
          IF (ZIDEB /=0 .AND. ZIDEB/=-1) ZLATLON(3)=zideb*1000.
          IF (ZIFIN /=0 .AND. ZIFIN/=-1) ZLATLON(4)=zifin*1000.
         
          ILATLON(:)=ZLATLON(:)
          
          IF (ILATLON(1)>  ZLATLON(1)) ILATLON(1)=ILATLON(1)-1
          IF (ILATLON(2)<  ZLATLON(2)) ILATLON(2)=ILATLON(2)+1
          IF (ILATLON(3)<  ZLATLON(3)) ILATLON(3)=ILATLON(3)+1
          IF (ILATLON(4)>  ZLATLON(4)) ILATLON(4)=ILATLON(4)-1
 
          ZLATLON=ILATLON
          IF (CFIXRESOL=="y") THEN
            INX=(ZLATLON(4)-ZLATLON(3))/ZDX_GRB +1 
            INY=(ZLATLON(1)-ZLATLON(2))/ZDY_GRB +1 
            ZCONTROL=ZLATLON(3)+(INX-1)*ZDX_GRB
            IF (ZCONTROL/=ZLATLON(4)) THEN
              print*,"warning : need to change E longitude"
              ZLATLON(4)=ZCONTROL
              print*,"lon min" ,ZLATLON(3)
              print*,"lon max" ,ZLATLON(4)
              print*,"dx",ZDX_GRB
            ENDIF
            ZCONTROL=ZLATLON(2)+(INY-1)*ZDY_GRB
            IF (ZCONTROL/=ZLATLON(1)) THEN
              print*,"warning : need to change N latitude"
              ZLATLON(1)=ZCONTROL
              print*,"lat min" ,ZLATLON(2)
              print*,"lat max" ,ZLATLON(1)
              print*,"dy",ZDY_GRB
            ENDIF
            print*,INX,INY,ZLATLON,ZDX_GRB,ZDY_GRB
          ELSE
            CALL INI2LALO(ZLATLON,INX,INY)
          ENDIF
          ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
                              size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
          ZVARSAVE=XVAR
          ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
          ALLOCATE(ZVARZCST(INX,INY,size(XVAR,3)))
          DEALLOCATE(XVAR)
          ALLOCATE(XVAR(INX,INY,SIZE(ZVARSAVE,3),&
                          size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)))

          DO J6=ivarprocinf,ivarprocsup
            DO J5=ivartrajinf,ivartrajsup
              DO J4=ivartinf,ivartsup
                ZWORK3D(:,:,:)=ZVARSAVE(:,:,:,J4,J5,J6)
                CALL INT2LALO('BILI',ZWORK3D,ZLATLON,XSPVAL,ZVARZCST)
                XVAR(:,:,:,J4,J5,J6)=ZVARZCST
              END DO
            END DO
          END DO
          DEALLOCATE(ZVARSAVE,ZVARZCST,ZWORK3D)
        ENDIF
       print*," ZLATLON apres INT2lalo",ZLATLON
                 ! d. ecriture des donnees au format  GRIB
        IF ( YOUTGRID(1:4) == 'LALO' ) THEN
          IF ( IFLAGzcst /= 1 ) THEN
            ivarideb=1
            ivarifin=SIZE(XVAR,1)
            ivarjdeb=1
            ivarjfin=SIZE(XVAR,2)
          ENDIF
          IF (LLEVEL2D) THEN
                 CALL WRITEGRIB(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YTYPEOUT,&
                       ilocverbia,iret,ICODCOD,&
                       zlistevert,LVAR2D,KLEVEL2D=ILEVEL2D,PLATLON=ZLATLON)        
          ELSE
                 CALL WRITEGRIB(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YTYPEOUT,&
                       ilocverbia,iret,ICODCOD,&
                       zlistevert,LVAR2D,PLATLON=ZLATLON)        
          ENDIF 
          if (ilocverbia > 0 ) then
             print*,'WRITEGRIB LALO return= ', YTYPEOUT,'= ',iret
          end if
        ELSE ! pas d interpolation horizontale (CONF)
           IF (LCARTESIAN) THEN
                 PRINT*,"===================================="
                 PRINT*,"WARNING : WITH LCARTESIAN=TRUE PLEASE ASK LALO"
                 PRINT*,"===================================="
                 STOP      
           ENDIF
           IF (LLEVEL2D) THEN
                 CALL WRITEGRIB(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YTYPEOUT,&
                       ilocverbia,iret,ICODCOD,zlistevert,LVAR2D,KLEVEL2D=ILEVEL2D)
           ELSE
                 CALL WRITEGRIB(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                       ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                       ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                       CGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YTYPEOUT,&
                       ilocverbia,iret,ICODCOD,zlistevert,LVAR2D)      
           ENDIF         
           if (ilocverbia > 0 ) then
                 print*,'WRITEGRIB CONF return= ', YTYPEOUT,'= ',iret
           end if
        ENDIF
        ! retour a XZZ pour NGRID a 4 (cf readvar)
        CALL COMPCOORD_FORDIACHRO(4)
     END SELECT
     ! indiquera aux routines d ecriture que le fichier courant est deja ouvert
     YFLAGWRITE='OLD'
  ! 
  ELSE   ! iret /=0
    print *, ' READVAR return= ',iret
  ENDIF  
END DO ! boucle champ a traiter
!
!
!---------------------------------------------------------------------------
!
!*       4.    CLOSURE OF OUTPUT FILE
!              ----------------------
!
!pour clore le traitement meme si la liste des champs est non terminee par END
88 CONTINUE
!
IF (ALLOCATED(ZNEWX))   DEALLOCATE(ZNEWX,ZNEWY)
IF (ALLOCATED(ZNEWLAT)) DEALLOCATE(ZNEWLAT,ZNEWLON)
IF (ALLOCATED(ZWORK2D)) DEALLOCATE(ZWORK2D,ZWORK2D2)
!
PRINT*, 'END ->  Close the output file'
YFLAGWRITE='CLO'
SELECT CASE(YTYPEOUT(1:4))
  CASE('DIAC')
    CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                  ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,  &
                 CGROUP,YFILEIN,YFLAGWRITE,'2  ',ilocverbia,iret)
  CASE('LLHV','llhv','LLZV','llzv','LLPV','llpv','LLAV','llav',&
          'IJHV','IJZV','IJPV','jihv','jizv','jipv')             
    CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                 ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                 CGROUP,YFILEIN,YFLAGWRITE,YTYPEOUT,ilocverbia,iret)      
  CASE('KCDL','ZCDL','PCDL')
    CALL WRITECDL(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                  ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                  CGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YSUFFIX_file,      &
                  ilocverbia,iret,PGRIDX=XXX(:,IGRID),PGRIDY=XXY(:,IGRID))
  CASE('ZGRB','PGRB','AGRB')
    CALL WRITEGRIB(ivarideb,ivarifin,ivarjdeb,ivarjfin, &
                   ivarkdeb,ivarkfin,ivartinf,ivartsup, &
                   ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                   CGROUP,YFILEIN,YFLAGWRITE,YOUTGRID,YTYPEOUT,&
                   ilocverbia,iret,ICODCOD,zlistevert,LVAR2D)   
  CASE DEFAULT
    PRINT*, 'Closure of output type ',YTYPEOUT ,' not coded'
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       5.    END
!              ---
!
99 CONTINUE
PRINT*, 'Delete the links if necessary'
YDUMMYFILE=''
CALL CREATLINK(' ',YDUMMYFILE,'CLEAN',ILOCVERBIA)
PRINT*, 'The file ',TRIM(YLUDIR),' stores all the input directives '
PRINT*, ' you must give a new name to use it again'
CLOSE(ILUDIR)
!
!-------------------------------------------------------------------------------
!
CONTAINS 
!
!------------------------------------------------------------------------------
!
SUBROUTINE TAB2SPACE(HTEXT)
IMPLICIT NONE
CHARACTER(len=*),INTENT(INOUT) :: HTEXT

CHARACTER, PARAMETER :: YPTAB = CHAR(9) ! TAB character is ASCII : 9
CHARACTER, PARAMETER :: YPCOM = CHAR(44)! COMma character is ASCII : 44
INTEGER              :: JI

DO JI=1,LEN_TRIM(HTEXT)
  IF (HTEXT(JI:JI)==YPTAB .OR. HTEXT(JI:JI)==YPCOM) HTEXT(JI:JI) = ' '
END DO
END SUBROUTINE TAB2SPACE
!------------------------------------------------------------------------------

END PROGRAM EXTRACTDIA
!
