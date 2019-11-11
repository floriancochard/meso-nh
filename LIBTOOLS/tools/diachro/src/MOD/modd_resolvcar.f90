!     ######spl
      MODULE  MODD_RESOLVCAR
!     ######################
!
!!****  *MODD_RESOLVCAR* - 
!!
!!    PURPOSE
!!    -------
!
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        24/11/95
!!     updated   PM    21/11/94
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE

!
! Unite logique fichier de directives
!
INTEGER           :: NDIR

!
! NLOOPN = JLOOPN OPER memorise pour SSOL + DRST + RAPL + RSPL
!
INTEGER           :: NLOOPN

!
! NLOOPP = JLOOPP OPER memorise pour CART
!
INTEGER           :: NLOOPP

!
! NLOOPK = JLOOPK OPER memorise pour CART et LANIMK
! XLOOPZ = JLOOPZ               "                pour les niveaux =/= niv.modele
!
INTEGER           :: NLOOPK
REAL,SAVE         :: XLOOPZ

!
! Logiques de gestion des differents types d'operations
!
LOGICAL           :: LCH, LCV, LPV, LPH, LPVT, LCN, LFT, LFT1, LPVKT, L1K
LOGICAL           :: LTK, LPR, LPVKT1, LZT, LXT, LYT, LXYDIA, LZTPVKT1
LOGICAL           :: LXYWINCUR=.FALSE.
LOGICAL           :: LPXT, LPYT
LOGICAL           :: LEV       ! Potential vorticity
LOGICAL           :: LMINUS, LPLUS
LOGICAL           :: LANIMK, LANIMT
LOGICAL           :: LCNCUM, LCNSUM, LCHXY, LCVXZ, LCVYZ, LRS, LRS1
LOGICAL           :: L1DT
LOGICAL           :: LPRINT=.FALSE., LPRINTXY=.FALSE.
! Ecriture des dates dans le fichier FICVAL . Associe a XPRDAT
LOGICAL           :: LPRDAT=.FALSE.
LOGICAL           :: LPOINTG=.FALSE.
LOGICAL           :: L2DBX=.FALSE., L2DBY=.FALSE.
LOGICAL           :: LXYO=.FALSE.
LOGICAL,SAVE      :: LMNMXLOC=.FALSE.
LOGICAL           :: LXABSC=.FALSE. ! Cas V(x,t) --> X en abscisse ou non
LOGICAL           :: LXMINTOP=.FALSE. ! Cas V(x,t) --> Min X en ordonnee en haut
                                  ! ou non
! Streamlines
LOGICAL,SAVE      :: LSTREAM=.FALSE.
LOGICAL,SAVE      :: LINTERPOLSTR=.FALSE.
INTEGER,SAVE      :: NZSTR=80, NARSTR=4
INTEGER,SAVE      :: NSGD, NSEUIL
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XZSTR
REAl,SAVE         :: XLWSTR=1., XARLSTR=.009, XSSP=.004

! Logique presence ou non des labels sur les axes .T. -> absence .F.presence
LOGICAL           :: LNOLABELX=.FALSE.
LOGICAL           :: LNOLABELY=.FALSE.
!
! logiques de gestion des combinaisons des composantes du vent
!
LOGICAL           :: LUMVM, LUTVT, LMUMVM, LMUTVT
LOGICAL           :: LULM, LULT, LVTM, LVTT, LULMWM, LULTWT
LOGICAL           :: LSUMVM, LSUTVT, LMLSUMVM, LMLSUTVT
! Representation anterieure de ULM et VTM
LOGICAL,SAVE      :: LULMVTMOLD=.FALSE.
! CH orientation pour calcul ULM et VTM
REAl,SAVE         :: XANGULVT
! CH direction du vent
LOGICAL,SAVE      :: LDIRWIND, LDIRWM, LDIRWT
! Cas LDIRWIND _PVT_  Memorisation des temps
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XTDIRWIND
! Logique de gestion des statistiques des vecteurs vents
LOGICAL,SAVE      :: LVST=.FALSE.
! Logique pour supprimer la dilatation de la la composante W ds 
! representation ULMWM ou ULTWT
LOGICAL,SAVE      :: LDILW=.TRUE.
! Logique pour conserver les valeurs > XVHC ds le cas ou XVHC est <0
! c.a.d ou = scale
LOGICAL,SAVE      :: LVSUPSCA=.FALSE.
!
! logique de gestion des couleurs en cas de zoom en CV
! =F -> iso +couleurs zoom ident. au graphique integral
! =T ->        "       "   fction du min et du max du zoom
!
LOGICAL,SAVE      :: LCVZOOM=.FALSE.
!
! logique de gestion des couleurs a 0 (background)
!
LOGICAL           :: LCOLZERO=.FALSE.
!
! rang de l'index de couleur a mettre a 0
!
INTEGER           :: NCOLZERO
!
! RS
!
INTEGER           :: NIRS, NJRS
REAL              :: XIRS, XJRS, XIRSCC, XJRSCC
LOGICAL           :: LNOUVRS=.FALSE.
!
! Ajout (ou - ou * ) constante entre ()
!
CHARACTER(LEN=20),DIMENSION(100) :: CFACT=' '  !(Nom gpe + fact * ou + ou -)
REAL,DIMENSION(100)    :: XCONSTANTE=0.
INTEGER,DIMENSION(100)    :: NOPE(100)
INTEGER :: NPARG=0, NPARD=0              ! Parentheses Gauche et dte
INTEGER :: NOPEL=0                         ! Compteur fact * ou + ou -
LOGICAL :: LFACTIMP=.TRUE.               ! Impression fact * ou + ou -
!
! Superpositions
!
INTEGER :: NLOOPSUPER     ! Indice de boucle des superpositions dans pg pal
LOGICAL :: LSUPERDIA
INTEGER :: NSUPERDIA
! Fev 2001
CHARACTER(LEN=600),DIMENSION(100) :: CARSUP
!CHARACTER(LEN=240),DIMENSION(50) :: CARSUP
INTEGER,DIMENSION(100) :: NFILESCUR
!
! Cas superpositions CV 3D + PH 2D Hor. (Oct 2000)
! D'une maniere generale superp. coupes =/=
! Conventions actuelles . 1 pour 1 CV 1+2=3 pour CV+K=PH
!
INTEGER,DIMENSION(100),SAVE :: NHISTORY=0
REAL, SAVE   :: XLWPH1=2,XLWPH2=2,XLWPH3=2,XLWPH4=2,XLWPH5=2,XLWPH6=2
REAL, SAVE   :: XLWPH7=2,XLWPH8=2
!
! Temps
!
! Les 2 derniers indices = superpositions + n.traj ou station...
! 23/04/03 dim. n augmentee de 20 a 45
! 17/01/05 dim. n augmentee de 45 a 100
LOGICAL,DIMENSION(100,100)           :: LTIMEDIALL, LTINCRDIA     

INTEGER,DIMENSION(100,100)           :: NBTIMEDIA

INTEGER, DIMENSION(120,100,100)  :: NTIMEDIA

REAL, DIMENSION(120,100,100)     :: XTIMEDIA
!
! Processus
!
! Dernier indice = superpositions
LOGICAL,DIMENSION(100)           :: LPROCDIALL, LPINCRDIA

INTEGER,DIMENSION(100)           :: NBPROCDIA

INTEGER, DIMENSION(120,100)  :: NPROCDIA
!
! Niveaux K
!
! Les 2 derniers indices = superpositions + n.traj ou station...
! 10/10/07 dim. 1 augmentee de 120 a 160 (nb de niveaux K)
LOGICAL,DIMENSION(100,100)           :: LVLKDIALL, LKINCRDIA

INTEGER,DIMENSION(100,100)           :: NBLVLKDIA

INTEGER, DIMENSION(160,100,100)  :: NLVLKDIA
!
! Niveaux Z
!
LOGICAL,DIMENSION(100)           :: LZINCRDIA

INTEGER,DIMENSION(100)           :: NBLVLZDIA

INTEGER, DIMENSION(120,100)  :: NLVLZDIA

REAL, DIMENSION(120,100)  :: XLVLZDIA
!
! Numeros masques ou trajectoires
!
LOGICAL,DIMENSION(100)           :: LNDIALL, LNINCRDIA

INTEGER,DIMENSION(100)           :: NBNDIA

INTEGER, DIMENSION(120,100)  :: NNDIA
!
! Nom du groupe
!
CHARACTER(LEN=100) :: CTIMECS
CHARACTER(LEN=16) :: CGROUP, CTIMEC
CHARACTER(LEN=22) :: CUNITGAL
CHARACTER(LEN=40) :: CTITGAL
CHARACTER(LEN=16),DIMENSION(100) :: CGROUPS
!
! Intervalle des isocontours , extremes ou valeurs
!
REAL      :: XDIAINT
REAL      :: XISOMIN, XISOMAX
REAL,DIMENSION(300) :: XISOLEV
REAL      :: XISOREF
!
! Nb chiffres signicatifs pour les High and Low isocontours + cste(champ cst)
!
INTEGER,SAVE :: NSD=0
CHARACTER(LEN=10),SAVE :: CFMTMNMX=' '
!
! Formats axe X et axe Y et possibilite de * par un facteur ou donner bornes
!
LOGICAL,SAVE :: LFMTAXEX=.FALSE., LFMTAXEY=.FALSE.
CHARACTER(LEN=10),SAVE :: CFMTAXEX=' ', CFMTAXEY=' '
! Taille des labels= NSZLBX/1024 et NSZLBY/1024
INTEGER,SAVE :: NSZLBX=10, NSZLBY=10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!! pour les retrotrajectoires                                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Formats labels pour retrotrajectoire (ctraj3d_group)
!
LOGICAL,SAVE :: LFMTRTRAJ=.FALSE.
CHARACTER(LEN=10),SAVE :: CFMTRTRAJ='(E10.5)'
REAL,SAVE :: NSZRTRAJ=10.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Mars 2001
LOGICAL,SAVE :: LFACTAXEX=.FALSE., LFACTAXEY=.FALSE.
LOGICAL,SAVE :: LAXEXUSER=.FALSE., LAXEYUSER=.FALSE.
REAL,   SAVE :: XFACTAXEX=1., XFACTAXEY=1.
REAL,   SAVE :: XAXEXUSERD=1., XAXEXUSERF=1.
REAL,   SAVE :: XAXEYUSERD=1., XAXEYUSERF=1.
!
! Profil vertical: indice en X dans la CV (NLMAX,IKU)
!
INTEGER   :: NPROFILE
!
! PV Bornes en X
REAL      :: XPVMINTRUE, XPVMAXTRUE  ! Fournies par l'utilisateur
REAL      :: XPVMINT, XPVMAXT        ! Fournies par l'utilisateur (Bornes ident.
!
! PV Epaisseur traits des =/= profils et figures
REAL,save :: XLWPV1=0., XLWPV2=0., XLWPV3=0., XLWPV4=0.
REAL,save :: XLWPV5=0., XLWPV6=0., XLWPV7=0., XLWPV8=0.
!*JD*Mars2009 Pour les budgets
REAL,save :: XLWPV9=0., XLWPV10=0., XLWPV11=0., XLWPV12=0.
REAL,save :: XLWPV13=0., XLWPV14=0., XLWPV15=0.
!*JD*Mars2009 Pour les budgets
REAL,save :: XSTYLPV1=0., XSTYLPV2=0., XSTYLPV3=0., XSTYLPV4=0.
REAL,save :: XSTYLPV5=0., XSTYLPV6=0., XSTYLPV7=0., XSTYLPV8=0.
!*JD*Mars2009 Pour les budgets
REAL,save :: XSTYLPV9=0., XSTYLPV10=0., XSTYLPV11=0., XSTYLPV12=0.
REAL,save :: XSTYLPV13=0., XSTYLPV14=0., XSTYLPV15=0.
!*JD*Mars2009 Pour les budgets
! PV =/= entre parametre de GSLN pour 1 PV si > 4 (valeur max autorisee) et 1
! pour assurer le passage de cette valeur entre trapro et pro1d
INTEGER,save :: NGSLNP=0
!
!*JD*Mars2009 Gestion nom variables + position + taille
!*JD*Mars2009
LOGICAL,SAVE :: LVARNPVUSER=.FALSE.
CHARACTER(LEN=22),SAVE :: CVARNPV1=' ',CVARNPV2=' ',CVARNPV3=' ',CVARNPV4=' '
CHARACTER(LEN=22),SAVE :: CVARNPV5=' ',CVARNPV6=' ',CVARNPV7=' ',CVARNPV8=' '
CHARACTER(LEN=22),SAVE :: CVARNPV9=' ',CVARNPV10=' ',CVARNPV11=' ',CVARNPV12=' '
CHARACTER(LEN=22),SAVE :: CVARNPV13=' ',CVARNPV14=' ',CVARNPV15=' '
REAL,save :: XSZVARNPVTOP=0.,XSZVARNPVBOT=0.
REAL,save :: XPOSXVARNPV1TOP=0.,XPOSXVARNPV5BOT=0.
REAL,save :: XPOSYVARNPV1TOP=0.,XPOSYVARNPV5BOT=0.
!*JD*Mars2009  Ligne Zero sur PV
LOGICAL,SAVE :: LINZEROPV=.FALSE.
INTEGER,SAVE    :: NSTYLINZEROPV=1
!*JD*Mars2009 
LOGICAL,SAVE :: LCONVG2MASS
!*JD*Mars2009 
LOGICAL,SAVE :: LVARNPHUSER=.FALSE.
CHARACTER(LEN=22),SAVE :: CVARNPH1=' ',CVARNPH2=' ',CVARNPH3=' ',CVARNPH4=' '
CHARACTER(LEN=22),SAVE :: CVARNPH5=' ',CVARNPH6=' ',CVARNPH7=' ',CVARNPH8=' '
! pour plusieurs variables
!
! PVKT + FT  Bornes en Y des processus
REAL      :: XPVMIN, XPVMAX          ! Calculees par le programme
! 
! Zones de recouvrement
INTEGER,SAVE   :: NBRECOUV
INTEGER,DIMENSION(20),SAVE  :: NRECOUV
!Mai 2000
! FT + PVKT + FT1 Epaisseur des traits de l'ensemble des traces
! Presence de valeurs manquantes
!
REAL,SAVE      :: XLWFTALL=2.
REAL,SAVE      :: XSPVALT
LOGICAL,SAVE   :: LSPVALT=.FALSE.
!
! FT + PVKT + FT1 Bornes des variables en X (Temps)
LOGICAL,SAVE        :: LTIMEUSER=.FALSE.
LOGICAL,SAVE        :: LFTCLIP=.TRUE.  ! pour desactiver le clipping 
! cas LTIMEUSER=F et LMNMXUSER=F pour eviter la disparition de traits 
!aux bornes lors conversion en PS
REAL,SAVE      :: XTIMEMIN, XTIMEMAX
!
! FT + PVKT + FT1 Bornes des variables en Y, Noms des variables, Nb de variables
REAL,SAVE      :: XFTMIN, XFTMAX, XPVKTMIN, XPVKTMAX ! Fournies par l'utilisateur
REAL,SAVE      :: XFT1MIN, XFT1MAX ! Fournies par l'utilisateur cas +sieurs var.
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XFTMN, XFTMX
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: NCOLI
CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: CFTMN, CFTMX, CCOLI
INTEGER,SAVE   :: NBFTMN=0, NBFTMX=0, NBCOLI=0
INTEGER,SAVE   :: NCOLIVAL
LOGICAL,SAVE        :: LMNMXUSER=.FALSE.
LOGICAL,SAVE        :: LCOLUSER=.FALSE.
LOGICAL,SAVE   :: LOK=.FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FT1 Format courbes par USER (JD 240209) et nom variables representees
! le tout controle par LFT1LUSER=T (et LCOLINE=T pour les entiers et reels)
!
LOGICAL,SAVE        :: LFT1LUSER=.FALSE.
INTEGER,SAVE   :: NFT1COL1=0, NFT1COL2=0, NFT1COL3=0, NFT1COL4=0, NFT1COL5=0
INTEGER,SAVE   :: NFT1COL6=0, NFT1COL7=0, NFT1COL8=0, NFT1COL9=0, NFT1COL10=0
INTEGER,SAVE   :: NFT1COL11=0, NFT1COL12=0, NFT1COL13=0, NFT1COL14=0, NFT1COL15=0
!
REAL,SAVE      :: XFT1LW1=2.,XFT1LW2=2.,XFT1LW3=2.,XFT1LW4=2.,XFT1LW5=2.
REAL,SAVE      :: XFT1LW6=2.,XFT1LW7=2.,XFT1LW8=2.,XFT1LW9=2.,XFT1LW10=2.
REAL,SAVE      :: XFT1LW11=2.,XFT1LW12=2.,XFT1LW13=2.,XFT1LW14=2.,XFT1LW15=2.
!
INTEGER,SAVE      :: NFT1STY1=1,NFT1STY2=1,NFT1STY3=1,NFT1STY4=1,NFT1STY5=1
INTEGER,SAVE      :: NFT1STY6=1,NFT1STY7=1,NFT1STY8=1,NFT1STY9=1,NFT1STY10=1
INTEGER,SAVE      :: NFT1STY11=1,NFT1STY12=1,NFT1STY13=1,NFT1STY14=1,NFT1STY15=1
!
CHARACTER(LEN=10) :: CFT1TIT1='          ',CFT1TIT2='          ', &
CFT1TIT3='          ',CFT1TIT4='          ',CFT1TIT5='          '
CHARACTER(LEN=10) :: CFT1TIT6='          ',CFT1TIT7='          ', &
CFT1TIT8='          ',CFT1TIT9='          ',CFT1TIT10='          '
CHARACTER(LEN=10) :: CFT1TIT11='          ',CFT1TIT12='          ', &
CFT1TIT13='          ',CFT1TIT14='          ',CFT1TIT15='          '
!
! FT1 gestion de la fenetre par l utilisateur
!
LOGICAL,SAVE   :: LVPTFT1USER=.FALSE.
REAL,SAVE      :: XVPTFT1L,XVPTFT1R,XVPTFT1B,XVPTFT1T
!
! Pour supprimer les labels a droite du dessin cas FT1 (je ne sais pas pour FT)
!
LOGICAL,SAVE   :: LBLFT1SUP=.FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! JD Avril 2009
! Trace filaire 2D Suppression des noms de var. et leur figure en Top
LOGICAL,SAVE   :: LXYNVARTOP=.TRUE.
LOGICAL,SAVE   :: LXYSTYLTOP=.TRUE.
LOGICAL,SAVE   :: LPHCOLUSER=.FALSE.
LOGICAL,SAVE   :: LPHSTYUSER=.FALSE.
INTEGER,SAVE      :: NPHSTY1=1,NPHSTY2=1,NPHSTY3=1,NPHSTY4=1,NPHSTY5=1
INTEGER,SAVE      :: NPHSTY6=1,NPHSTY7=1,NPHSTY8=1
INTEGER,SAVE      :: NPHCOL1=1,NPHCOL2=1,NPHCOL3=1,NPHCOL4=1,NPHCOL5=1
INTEGER,SAVE      :: NPHCOL6=1,NPHCOL7=1,NPHCOL8=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Isocontours : cas NIMNMX=1 
!
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XISOMN, XISOMX, XISOINT
CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: CISOMN, CISOMX, CISOINT
INTEGER,SAVE   :: NBISOMN, NBISOMX, NBISOINT
LOGICAL        :: LISOK=.FALSE.
!
! Isocontours : cas NIMNMX=2 
!
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XISOLEVP
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: NLENP
CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: CISOLEVP
INTEGER,SAVE   :: NBISOLEVP
LOGICAL        :: LISOLEVP=.FALSE.
!
! Isocontours : cas NIMNMX=3 
!
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XISOREFP
CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: CISOREF
INTEGER,SAVE   :: NBISOREF
LOGICAL        :: LISOREF=.FALSE.
!
! Isocontours : epaisseurs des traits en cas de superpositions ou non
!
REAL,save      :: XLW1, XLW2, XLW3, XLW4
! Epaisseur traits continents
REAL,save      :: XLWCONT=0.
!
! LXY=.TRUE.    Bornes en Y
REAL,SAVE      :: XVARMIN=0., XVARMAX=0.
!
! LZT=.TRUE.    Bornes en Y
REAL,SAVE      :: XZTMIN=0., XZTMAX=0.
!
! Sommes et differences
!
INTEGER        :: NBPM          ! Nb de sommes et differences et superpositions 
INTEGER        :: NBPMT         ! Nb de sommes et differences
INTEGER,DIMENSION(99)   :: NUMPM  ! 1 --> +   2 --> -   0 --> rien
!
! Vents /= vent normal : suffixe
!
CHARACTER(LEN=2),SAVE :: CSUFWIND='  '
INTEGER,SAVE          :: NSUFWIND=0
!
LOGICAL :: LSYMBT
INTEGER,SAVE          :: NCOLSYMB, NTYPSYMB
!
! Spectres
!
LOGICAL,SAVE  :: LINDSP=.FALSE.
LOGICAL,SAVE  :: LOGNEP=.TRUE.
LOGICAL,SAVE  :: LM5S3=.FALSE.
LOGICAL,SAVE  :: LSPMNMXUSER=.FALSE., LSPMNMXALLT=.FALSE.
LOGICAL,SAVE  :: LSPLO=.FALSE., LSPO=.FALSE., LOSPLO=.FALSE., LPHALO=.FALSE., LPHAO=.FALSE.
LOGICAL,SAVE  :: LSPSECT=.FALSE., LSPSECTXY=.FALSE., LSPSECTXZ=.FALSE., LSPSECTYZ=.FALSE.
LOGICAL,SAVE  :: LSPX=.FALSE., LSPY=.FALSE., LSPZ=.FALSE.
REAL          :: XOMEGAX, XOMEGAY, XOMEGAZ
REAL,SAVE     :: XSPMIN=0., XSPMAX=0.

LOGICAL,SAVE  :: LBID
!
! Table de couleurs N2
!
LOGICAL,SAVE  :: LTABCOLDEF2=.FALSE.
!
! Trajectoires
!
LOGICAL,SAVE  :: LCONV2XY=.FALSE., LCONT=.FALSE., LRELIEF=.FALSE.
! L2CONT pour trace comme anterieurement des continents 2fois (Septembre 2000)
LOGICAL,SAVE  :: L2CONT=.FALSE.
INTEGER,SAVE  :: NLATLON
!
INTEGER,SAVE  :: NVERBIA=0, NSSPG=0
! LINVWB=.FALSE. (1,0,0.,0.,0.), (1,1,1.,1.,1.)
! LINVWB=.TRUE. (1,1,0.,0.,0.), (1,0,1.,1.,1.)
! Definition Noir et Blanc
LOGICAL,SAVE  :: LINVWB=.TRUE.
!
LOGICAL,SAVE  :: LISOWHI2=.FALSE., LISOWHI3=.FALSE.
!
! Vecteurs vent (Les autres parametres de meme type sont ds MODN_NCAR)
! = valeur en deca de laquelle les vecteurs ne st pas representes
REAL,SAVE :: XVLC=0.
!
! Textes + symboles a ajouter a des plans horizontaux localises
! a XLATCAR,XLONCAR (definis dans MODN_NCAR)
!
CHARACTER(LEN=20),DIMENSION(400) :: CNOMCAR=' '
CHARACTER(LEN=1),DIMENSION(400) :: CSYMCAR='.'
REAL,DIMENSION(400)              :: XPOSNOM=90.
REAL,DIMENSION(400)              :: XSZNOM=.012
REAL,DIMENSION(400)              :: XSZSYM=.012
INTEGER,DIMENSION(400)              :: ICOLSYM=1
INTEGER,DIMENSION(400)              :: ICOLNOM=1
INTEGER                         :: NOMCAR, NSYMCAR, NPOSNOM
INTEGER                         :: NSZNOM, NSZSYM, NCOLSYM, NCOLNOM
!
! Tableau utilise ds imcoupv pour charger les composantes u et v (UMVM_PVT_)
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XTEM2D, XTEM2D2
LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: LUMVMPVT
LOGICAL,SAVE                          :: LUMVMPV=.FALSE.
! 
! Logique de gestion d'interpolation a partir du haut ou du bas
!
LOGICAL,SAVE                          :: LINTERPTOP=.TRUE.
! Logique precisant que l'on demande les niv. des CH en reel
LOGICAL,SAVE                          :: LCHREEL=.FALSE.
!
! Profils horizontaux 07042000 pour trace UTVT ou UMVM  pour recuperer
! les coordonnees du debut des fleches
!
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XTEMCVU, XTEMCVV
! PH UTVT et UMVM pour expression des X en heures
! Ajoute ds imcou pour _PVT_ et LHEURX=T en Mai 2000
LOGICAL,SAVE                          :: LHEURX=.TRUE.
LOGICAL,SAVE                          :: LMYHEURX=.FALSE.
INTEGER,SAVE                          :: NHEURXLBL=2
INTEGER,SAVE                          :: NHEURXGRAD=1
! Avril 2009 ds imcou pour _PVT_
! Possibilite mettre des temps exprimes sous forme hhHmm dont les bornes
! sont fournies par utilisateur dans XAXUSERD= et XAXUSERF=(reels a 2 decimales)
! avec LAXEXUSER=T LHEURX=T et LNOLABELX=T
! Borne de fin toujours > borne debut mais si on veut une expression 
! des heures entre 0 et 24H , on met L24H=T
LOGICAL,SAVE                          :: L24H=.FALSE.
! Avril 2009 ds imcou pour _PVT_
! Mai 2009 ds image + imcou
LOGICAL,SAVE                          :: LNOLBLBAR=.FALSE.
! Mai 2009 ds image + imcou
! Avril 2002 lat,lon CV et PH
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XLATCV, XLONCV
! 
! Limites du domaine fils sur le domaine pere. 
! Fournies en indices de grille du domaine pere
!
LOGICAL,SAVE         :: LDOMAIN=.FALSE.
INTEGER,SAVE         :: NDOMAINL=1, NDOMAINR=1, NDOMAINB=1, NDOMAINT=1
REAL,SAVE            :: XLWDOMAIN=2.
!
! Trace segment de dte sur une coupe horizontale en proj. cart
!
LOGICAL,SAVE         :: LSEGM=.FALSE.
! Elements du tableau entiers mis a 1 si XSEGMx =/= 0
INTEGER,DIMENSION(100),SAVE  :: NSEGMS=0
! nb couleurs lues
INTEGER,SAVE                :: NCOLSEGM=1
! Numeros des couleurs
INTEGER,DIMENSION(30),SAVE                :: NCOLSEGMS=1
! Couples lat,long extremites de segments de dte .
! Si =0,0 discontinuite ds les segments (plume levee!!)
!REAL,DIMENSION(2),SAVE  :: XSEGM1=0.,XSEGM2=0.,XSEGM3=0.,XSEGM4=0.,XSEGM5=0.
!REAL,DIMENSION(2),SAVE  :: XSEGM6=0.,XSEGM7=0.,XSEGM8=0.,XSEGM9=0.,XSEGM10=0.
!REAL,DIMENSION(2),SAVE  :: XSEGM11=0., XSEGM12=0., XSEGM13=0., XSEGM14=0.
REAL,DIMENSION(100,2),SAVE  :: XCONFSEGMS=0., XSEGMS=0.
REAL,SAVE            :: XLWSEGM=2.
!
! Logique d'inversion des pointilles et tiretes pour les isocontours N/B
!
LOGICAL,SAVE         :: LINVPTIR=.FALSE.
! 15052000 
! Pour impression de la fenetre papier courante
!
REAL,SAVE            :: XCURVPTL, XCURVPTR, XCURVPTB, XCURVPTT
! Pour ecriture des dates dans le fichier FICVAL
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE  ::  XPRDAT
!
! Ajout constante de temps pour ch courbe FT PVKT PVKT1
REAL,SAVE            :: XFT_ADTIM1=0
REAL,SAVE            :: XFT_ADTIM2=0
REAL,SAVE            :: XFT_ADTIM3=0
REAL,SAVE            :: XFT_ADTIM4=0
REAL,SAVE            :: XFT_ADTIM5=0
REAL,SAVE            :: XFT_ADTIM6=0
REAL,SAVE            :: XFT_ADTIM7=0
REAL,SAVE            :: XFT_ADTIM8=0
!
! Ajout constante de temps pour ch courbe FT1
REAL,SAVE            :: XFT1_ADTIM1=0
REAL,SAVE            :: XFT1_ADTIM2=0
REAL,SAVE            :: XFT1_ADTIM3=0
REAL,SAVE            :: XFT1_ADTIM4=0
REAL,SAVE            :: XFT1_ADTIM5=0
REAL,SAVE            :: XFT1_ADTIM6=0
REAL,SAVE            :: XFT1_ADTIM7=0
REAL,SAVE            :: XFT1_ADTIM8=0
!
! FT PVKT 3 ou 4 courbes / 1 diagramme (meme parametre avec bornes fixees)
LOGICAL,SAVE         :: LFT3C=.FALSE., LFT4C=.FALSE.
! 
! FT PVKT FT1 PVKT1 bornes calculees avec min et max effectifs
!(pour evol. temp. ds varfct)
LOGICAL,SAVE         :: LFTBAUTO=.FALSE.
LOGICAL,SAVE         :: LFT1BAUTO=.FALSE.
!
!NOVEMBRE 2009 : ajout de l apossibilité de tourner les titres en Y
LOGICAL,SAVE            ::L90TITYT=.FALSE.
LOGICAL,SAVE            ::L90TITYM=.FALSE.
LOGICAL,SAVE            ::L90TITYB=.FALSE.
LOGICAL,SAVE            ::LPATCH=.FALSE.

!
END MODULE MODD_RESOLVCAR
