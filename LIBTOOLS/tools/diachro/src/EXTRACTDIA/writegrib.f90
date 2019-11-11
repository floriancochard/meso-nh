!     #################################
      MODULE MODI_WRITEGRIB
!     #################################
INTERFACE WRITEGRIB
      SUBROUTINE  WRITEGRIB(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
                   kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin,&
                   HLABELCHAMP,HFILENAME,HFLAGFILE,HOUTGRID,HTYPEOUT,  &
                   KVERBIA,KRETCODE,KCODCOD,PLEV,OVAR2D,KLEVEL2D,PLATLON)
!
CHARACTER(LEN=*), intent(inout)  :: HLABELCHAMP         ! nom du champ
                       ! inout pour modifier le nom VLEV en altitude
CHARACTER(LEN=*), intent(in)  ::  HFILENAME             ! nom du fichier
CHARACTER(LEN=*), intent(in)  :: HFLAGFILE              ! NEW=creation 
                                                        ! OLD=ajout 
                                                        ! CLOSE=fermeture
CHARACTER(LEN=*), intent(in) :: HOUTGRID               ! format grille reguliere plan conforme
                                                        !ou lat lon CONF/LALO                                                        
CHARACTER(LEN=*), INTENT(in) :: HTYPEOUT              ! type de fichier sortie                                                        
INTEGER , intent(in)         :: KVERBIA                 ! prints de controle
                                      ! desactive (0) / active (1) les prints
                                      ! limites sur les 6 dimensions
INTEGER , intent(in)         :: kideb,kifin,kjdeb,kjfin,kkdeb,kkfin   
INTEGER , intent(in)         :: kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin
! 
INTEGER , intent(out)        :: KRETCODE   ! Code de retour de la routine 
INTEGER, INTENT(IN)               :: KCODCOD   ! parameter code
REAL, DIMENSION(:), INTENT(IN) :: PLEV !niveaux verticaux
LOGICAL,INTENT(IN)           :: OVAR2D ! champ 2D (surface) si TRUE sinon 3D
INTEGER,OPTIONAL,INTENT(IN) :: KLEVEL2D
REAL,DIMENSION(:),OPTIONAL,INTENT(IN) :: PLATLON
END SUBROUTINE
END INTERFACE
END MODULE MODI_WRITEGRIB
!
!     ################
      SUBROUTINE  WRITEGRIB(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
                   kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin,&
                   HLABELCHAMP,HFILENAME,HFLAGFILE,HOUTGRID,HTYPEOUT,  &
                   KVERBIA,KRETCODE,KCODCOD,PLEV,OVAR2D,KLEVEL2D,PLATLON)
!     ################
! pour getenv et system
#ifdef NAGf95
USE F90_UNIX
USE F90_UNIX_PROC       
#endif
!
USE MODN_NCAR,  ONLY: XSPVAL       
!
!                    grille : XXDXHAT(:,1:7) et XXX(:,1:7), XXZS(:,:,1:7)
USE MODD_COORD
!                     min max des indices selon x et y
USE MODD_TYPE_AND_LH
!                    XVAR(i,j,k,,,), XMASK,XTRAJ ,XDATIME(16,t)   
USE MODD_ALLOC_FORDIACHRO
USE MODD_FILES_DIACHRO, ONLY: NBFILES, CLUOUTDIAS, NRESPDIAS
!          
USE MODI_TEMPORAL_DIST ! interface modules
USE MODI_FROM_COMPUTING_UNITS
USE MODD_CONF
USE MODD_TIME, ONLY: TDTEXP,TDTSEG
USE MODD_TIME1, ONLY: TDTCUR
USE MODD_GRID
USE MODD_GRID1
!
USE MODN_OUTFILE
USE MODE_GRIDPROJ
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!              -----------------
!
CHARACTER(LEN=*), intent(inout)  :: HLABELCHAMP         ! nom du champ
                       ! inout pour modifier le nom VLEV en altitude
CHARACTER(LEN=*), intent(in)  :: HFILENAME              ! nom du fichier
CHARACTER(LEN=*), intent(in)  :: HFLAGFILE              !NEW=creation 
                                                        !OLD=ajout 
                                                        !CLOSE=fermeture
CHARACTER(LEN=*), intent(in) :: HOUTGRID               ! format grille reguliere plan conforme
                                                        !ou lat lon CONF/LALO                                                        
CHARACTER(LEN=*), INTENT(in) :: HTYPEOUT              ! type de fichier sortie
                                                        
INTEGER , intent(in)         :: KVERBIA                 ! prints de controle
                                      !desactive (0) / active (1) les prints
                                      ! limites sur les 6 dimensions
INTEGER , intent(in)         :: kideb,kifin,kjdeb,kjfin,kkdeb,kkfin   
INTEGER , intent(in)         :: kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin
! 
INTEGER , intent(out)        :: KRETCODE   ! Code de retour de la routine 
INTEGER, INTENT(IN)               :: KCODCOD   ! parameter code
REAL, DIMENSION(:), INTENT(IN) :: PLEV !niveaux verticaux
LOGICAL,INTENT(IN)           :: OVAR2D ! champ 2D (surface) si TRUE sinon 3D
INTEGER,OPTIONAL,INTENT(IN) :: KLEVEL2D
REAL,DIMENSION(:),OPTIONAL,INTENT(IN) :: PLATLON
!
!
!
!
INTEGER           :: ILOOP,JLOOP,KLOOP,KLOOP4,KLOOP5,KLOOP6, iret
INTEGER           :: IAN,IMOIS,IJOUR,ISECONDE,ibasetime
INTEGER           :: IAN2,IMOIS2,IJOUR2,ISECONDE2,IANREF
INTEGER, dimension(:), ALLOCATABLE :: ioffset_time
INTEGER  :: zbasetime
!DOUBLE PRECISION  :: zbasetime

!
REAL              :: zmini ,zmaxi
!
! taille=100  et 28 cf diaprog 
CHARACTER (LEN=100) :: YSAVETITRE, YSAVECOMMENT, YSAVEUNITE 
CHARACTER (LEN=40)  :: YFILEOUT   ! Fichier de sortie
CHARACTER (LEN=100) :: ycommand, ytextdim
CHARACTER (LEN=13), save :: YLIBELLEDIM1,YLIBELLEDIM2
CHARACTER (LEN=5)   :: YNUM
CHARACTER (LEN=28)  :: YLABELCHAMPnew
INTEGER :: ikdeb,ikfin,iitdeb,iitfin,iitrdeb,iitrfin,JK
CHARACTER (LEN=15)  :: YNETCDFCHAMP
CHARACTER  (LEN=8) :: YDATE
CHARACTER  (LEN=10) :: YTIME
CHARACTER  (LEN=5) :: YZONE
INTEGER,DIMENSION(8) :: IVALUES
REAL,DIMENSION(:,:,:,:,:), ALLOCATABLE:: XVAR2
INTEGER :: II,IJ,IK,IT,IM
!
!
INTEGER :: IGRIBFILE  ! logical unit for grib file
INTEGER :: IRESP
CHARACTER (LEN=22)  :: YFIELDGRIB
CHARACTER (LEN=6)  :: YLEV
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZFIELD
INTEGER :: JI,JJ,JJI,JJJ,IIX,IJY
CHARACTER (LEN=22)  ::YSUFFIX
!
! POUR GRIBEX
INTEGER, DIMENSION(2)             :: ISEC0    ! see gribex documentation
INTEGER, DIMENSION(1024)          :: ISEC1
INTEGER, DIMENSION(1024)          :: ISEC2
INTEGER, DIMENSION(2)             :: ISEC3
INTEGER, DIMENSION(512)           :: ISEC4
!
REAL,    DIMENSION(512)         :: ZSEC2
REAL,    DIMENSION(2)           :: ZSEC3
!
REAL,    DIMENSION(:),ALLOCATABLE :: ZSEC4
INTEGER                         :: IPUNP ! length of data array ZSEC4
INTEGER                         :: INBITS ! number of bits for coding
INTEGER, DIMENSION(:),ALLOCATABLE :: INBUFF ! grib buffer
INTEGER                         :: IPACK ! length of grib buffer INBUFF
CHARACTER(LEN=1)                :: YOPER ! requested function
INTEGER                         :: IWORD ! number of words of INBUFF occupied by coded data
INTEGER                         :: IERR ! return gribex code
!
REAL    :: ZLENGTH    ! length of forecast in seconds
REAL    :: ZLATREF2   ! second reference latitude in Lambert projection
REAL    :: ZMAP60     ! map factor at 60³ parallel nearest of the pole
INTEGER :: ITIME
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
! 
print *,' --------- '
print *,'Entree WRITEGRIB ',TRIM(HFILENAME),' ',TRIM(HLABELCHAMP),' ', &
                              TRIM(HFLAGFILE),' ',KVERBIA

KRETCODE=0
CALL FROM_COMPUTING_UNITS(HLABELCHAMP,CUNITE(1)) 
LPBREAD=.FALSE.                                                        
print'(A41,6(I4,X))','WRITEGRIB: ideb,ifin,jdeb,jfin,kdeb,kfin= ',&
          kideb,kifin,kjdeb,kjfin,kkdeb,kkfin
print'(A42,2(I10,X),4(I4,X))','          tdeb,tfin,trdeb,trfin,pdeb,pfin= ',&
          kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin 
print'(A26,6(I4,X))','          nil,nih,njl,njh,nkl,nkh=',nil,nih,njl,njh,nkl,nkh

YFILEOUT=TRIM(HFILENAME)//'.'//HTYPEOUT(2:4)
print*,'fichier de sortie   YFILEOUT= ',YFILEOUT
!   
!-------------------------------------------------------------------------------
!
!*       2.1   OUVERTURE DES FICHIERS DE SORTIE
!            -------------------
!
IF ( HFLAGFILE(1:3) == 'NEW' ) THEN
! Open the MULTIGRIB file if necessary     
  print*,'The output GRIB file is named: ', YFILEOUT
  CALL PBOPEN(IGRIBFILE,YFILEOUT,"W",IRESP)
  IF (IRESP /= 0) print*, 'ERROR when opening file, IRESP=',IRESP
END IF
!-------------------------------------------------------------------------------
!
!*       3.    ECRITURE du champ 
!              --------
!
IF ( HFLAGFILE(1:3) /= 'CLO' ) THEN
! Pour l'instant on ne triate qu'un seul temps, trajectoire/mask ou processus
! si on en veut un en particulier il faut le faire avec le zoom
! probleme s'il y a plusieurs temps dans un fichier ils auront tous le meme 
! dans le fichier grib donc il vaut mieux ne pas concatener les fichiers
  IF (kitdeb/=kitfin) THEN
   PRINT*,"=== WARNING ==="
   PRINT*," you are asking for several times : (",kitdeb,":",kitfin,")"
   PRINT*," only the first one (",kitdeb,") will be take into account"
  ENDIF
  IF (kitrdeb/=kitrfin) THEN
   PRINT*,"=== WARNING ==="
   PRINT*," you are asking for several trajectories : (",kitrdeb,":",kitrfin,")"
   PRINT*," only the first one (",kitrdeb,") will be take into account"
  ENDIF
  IF (kipdeb/=kipfin) THEN
   PRINT*,"=== WARNING ==="
   PRINT*," you are asking for several processus : (",kipdeb,":",kipfin,")"
   PRINT*," only the first one (",kipdeb,") will be take into account"
  ENDIF
!
  !=========================================
  ! ecriture de la section 1 du GRIB  
  !=========================================
  ISEC1(:)=0
  ISEC1(1)=1
  ISEC1(2)=85      !  Idendification of center : French Weather Service
  ISEC1(3)=96      ! Generating process identification number : MESONH identifier
  ISEC1(4)=255     ! Grid definition : non-standard grid definition
  ISEC1(5)=192     ! section 2 included, section 3 included (missing value)
  ISEC1(6)=KCODCOD ! parameter indicator
  ISEC1(10)=TDTEXP%TDATE%YEAR-100*(TDTEXP%TDATE%YEAR/100) ! year of century
  ISEC1(11)=TDTEXP%TDATE%MONTH ! month of reference date (start of experiment)
  ISEC1(12)=TDTEXP%TDATE%DAY ! day of reference date (start of experiment)
  ISEC1(13)=NINT(TDTEXP%TIME)/3600 ! hour of reference date (start of experiment)
  ISEC1(14)=NINT(TDTEXP%TIME)/60 - 60*ISEC1(13) ! minutes of reference date (start of exper0,0,0,0,0,0iment)
  ISEC1(15)=1                    ! time unit: hour
  CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH,TDTCUR%TDATE%DAY,TDTCUR%TIME, &
     TDTEXP%TDATE%YEAR,TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME, &
     ZLENGTH) 
  ISEC1(16)=NINT(ZLENGTH/3600)    ! length of forecast (period of time since the start of experiment
  !IF (NVERB>=5) print*, 'ZLENGTH=',ZLENGTH,'ISEC1=',ISEC1(16)
  print*, 'ZLENGTH=',ZLENGTH,'ISEC1=',ISEC1(16)
  ISEC1(21)=TDTEXP%TDATE%YEAR/100+1 ! century of data
  !=========================================
   
  ! zoom sur les dimensions

  ALLOCATE(ZFIELD(kifin-kideb+1,kjfin-kjdeb+1,kkfin-kkdeb+1))
  IIX=kifin-kideb+1
  IJY=kjfin-kjdeb+1
  IK=kkfin-kkdeb+1
!  print*,IIX,IJY,IK
!  print*,SHAPE(XVAR)
  ZFIELD(:,:,:)=XVAR(kideb:kifin,kjdeb:kjfin,kkdeb:kkfin,kitdeb,kitrdeb,kipdeb)      
!===========================================================================
!===========================================================================
!                       GRILLE LAT/LON REGULIERE
!===========================================================================
!===========================================================================
  IF (HOUTGRID=='LALO' )THEN 
    ISEC2(:)=0
    ISEC2(1)=0                          ! lat/lon regular grid
    ISEC2(2)=IIX                 
    ISEC2(3)=IJY
    ISEC2(4)=PLATLON(1)
    ISEC2(5)=PLATLON(3)
    ISEC2(6)= 128
    ISEC2(7)=PLATLON(2) 
    ISEC2(8)=PLATLON(4)
!   print*,"ISEC2(2),ISEC2(3)",ISEC2(2),ISEC2(3)
!   print*,"ISEC2(4),ISEC2(5), ISEC2(7),ISEC2(8)",ISEC2(4),ISEC2(5), ISEC2(7),ISEC2(8)
    ISEC2(9)= (ISEC2(8)-ISEC2(5))/(IIX-1)
    IF (ISEC2(9)<0) ISEC2(9)=(ISEC2(8)-ISEC2(5)+360000.)/(IIX-1)
    ISEC2(10)=(ISEC2(4)-ISEC2(7))/(IJY-1)
!   print*,"ISEC2(9),ISEC2(10)",ISEC2(9),ISEC2(10)
    !
    ! quelques verif de coherence 
    !
    IF (ISEC2(7)/= (ISEC2(4)-ISEC2(10)*(IJY-1))) THEN
            print*,"ERREUR : ISEC2(7)/= (ISEC2(4)-ISEC2(10)*(IJY-1)))"
            print*,"ISEC2(7)=",ISEC2(7)
            print*,"ISEC2(4)=",ISEC2(4)
            print*,"ISEC2(10)=",ISEC2(10)
            print*,"IJY=",IJY
            STOP
    ENDIF
    IF (ISEC2(8)/= (ISEC2(5)+ISEC2(9)*(IIX-1))) THEN
            print*,"ERREUR : ISEC2(8)/= (ISEC2(5)+ISEC2(9)*(IIX-1)))"
            print*,"ISEC2(8)=",ISEC2(8)
            print*,"ISEC2(5)=",ISEC2(5)
            print*,"ISEC2(9)=",ISEC2(9)
            print*,"IIX=",IIX
            STOP
    ENDIF

    ISEC2(11)=0                   ! scanning: i+, j-
    ISEC2(12)=IK              ! number of vertical levels
!===========================================================================
!===========================================================================
!            ON RESTE SUR LA PROJECTION CONFORME 
!===========================================================================
!===========================================================================
  ELSE IF (HOUTGRID=='CONF') THEN
     !     print*,"XRPK=",XRPK
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXX   Mercator               XXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          
    IF (ABS(XRPK)<1.E-10) THEN   
      ISEC2(:)=0
      ISEC2(1)=1
      ISEC2(2)=IIX                 
      ISEC2(3)=IJY
      ISEC2(4)=1000.*XLAT(kideb,kjdeb)   ! latitude of first point
      ISEC2(5)=1000.*(XLON(kideb,kjdeb) &! longitude of first point
           -360.*NINT(XLON(kideb,kjdeb)/360.))
      ISEC2(6)=0
      ISEC2(7)=1000.*XLAT(kifin,kjfin)   ! latitude of last point
      ISEC2(8)=1000.*(XLON(kifin,kjfin) &! longitude of last point
           -360.*NINT(XLON(kifin,kjfin)/360.))
      ISEC2(9)=1000.*XLAT0
      ISEC2(11)=64
      ISEC2(12)=IK
      ISEC2(13)=XXHAT(kideb+1)-XXHAT(kideb)  ! DX at XLAT0
      ISEC2(14)=XYHAT(kjdeb+1)-XYHAT(kjdeb)  ! DY at XLAT0
      ISEC2(19)=8
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXX   Polar Stereographic    XXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ELSE IF ( ABS(XRPK)==1. ) THEN 
      ISEC2(:)=0
      ISEC2(1)=5               ! polar stereographic projection
      ISEC2(2)=IIX       ! number of points along x
      ISEC2(3)=IJY       ! number of points along y
      ISEC2(4)=1000.*XLAT(kideb,kjdeb)     ! latitude of first point
      ISEC2(5)=1000.*(XLON(kideb,kjdeb) &  ! longitude of first point
           -360.*NINT(XLON(kideb,kjdeb)/360.))
      ISEC2(7)=1000.*(XLON0-360.*NINT(XLON0/360.))! longitude of the reference meridian
      IF (XRPK>0.) THEN
        ZMAP60=( COS(XLAT0*XPI/180.)   / COS(XPI/3.)   )**(1.-ABS(XRPK)) &
         *((1+SIN(XLAT0*XPI/180.))/(1+SIN(XPI/3.)))**(ABS(XRPK))
      ELSE IF (XRPK<0.) THEN
        ZMAP60=( COS(-XLAT0*XPI/180.)   / COS(-XPI/3.)   )**(1.-ABS(XRPK)) &
           *((1+SIN(-XLAT0*XPI/180.))/(1+SIN(-XPI/3.)))**(ABS(XRPK))
      END IF
      ISEC2(9)=(XXHAT(kideb+1)-XXHAT(kideb)) /ZMAP60  ! DX at 60³
      ISEC2(10)=(XYHAT(kjdeb+1)-XYHAT(kjdeb))/ZMAP60  ! DY at 60³
      ISEC2(11)=64                           ! scanning I+, J+, (I,J) (01000000)
      ISEC2(12)=IK
      IF (XRPK>1.-1.E-10) THEN
        ISEC2(13)=0                          ! North pole in the domain
      ELSE
        ISEC2(13)=1                          ! South pole in the domain
        ! bizarre normalement c'est 128 d'apres la doc mais ca ne marche pas
        ! par contre ok avec 1
      END IF
      ISEC2(19)=8
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXX    Conformal Lambert   XXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ELSE IF ( ABS(XRPK)<1. .AND. ABS(XRPK)>1.E-10 ) THEN 
      ISEC2(:)=0
      ISEC2(1)=3
      ISEC2(2)=IIX                 
      ISEC2(3)=IJY
      ISEC2(4)=1000.*XLAT(kideb,kjdeb)   ! latitude of first point
      ISEC2(5)=1000.*(XLON(kideb,kjdeb) &! longitude of first point
           -360.*NINT(XLON(kideb,kjdeb)/360.))
      print*,"ISEC2(4),ISEC2(5) :",ISEC2(4),ISEC2(5)
      ISEC2(6)=128
      ISEC2(7)=1000.*(XLON0-360.*NINT(XLON0/360.)) ! reference longitude
      ISEC2(8)=0
      ISEC2(9)=(XXHAT(kideb+1)-XXHAT(kideb)) 
      ISEC2(10)=(XYHAT(kjdeb+1)-XYHAT(kjdeb))
      print*,"ISEC2(9),ISEC2(10) :",ISEC2(9),ISEC2(10)
      ISEC2(11)=64                   
      ISEC2(12)=IK   
      IF (XRPK>0.) THEN
        ISEC2(13)=0                          ! North pole in the projection plane
        ZLATREF2=LATREF2(XLAT0,XRPK)
        ISEC2(14)=1000*MAX(XLAT0,ZLATREF2)
        ISEC2(15)=1000*MIN(XLAT0,ZLATREF2)
      ELSE
        ISEC2(13)=128                        ! South pole in the projection plane
        ZLATREF2=LATREF2(XLAT0,XRPK)
        ISEC2(14)=1000*MIN(XLAT0,ZLATREF2)
        ISEC2(15)=1000*MAX(XLAT0,ZLATREF2)
      END IF
      ISEC2(19)=8                            ! U and V along x and y axes
      ISEC2(20)=-90000                       ! latitude of south pole
    ELSE
      print*,"ERREUR : seules les projection stereographique,"
      print*," mercator ou lambert son reconnues"
      STOP
    END IF    
  ELSE 
    print*,"HOUTGRID=",HOUTGRID," non reconnu"
  ENDIF
!
  ZSEC2(:)=0
  ZSEC2(1)=XBETA*XPI/180.                ! angle of rotation (unit ????, supposed radian)
  ZSEC2(2)=1.
  DO JK=1,IK
    ZSEC2(JK+10)=PLEV(JK)
  END DO
!
  ISEC3(:)=0
  ISEC3(1)=0       ! missing data is considered
  ISEC3(2)=-1
  ZSEC3(:)=0
  ZSEC3(1)=0       ! not used 
  ZSEC3(2)=999. ! value for missing data
  !
  ALLOCATE(INBUFF((SIZE(ZFIELD,1)*SIZE(ZFIELD,2)*4)+4202))
  IPUNP=IIX*IJY
  INBITS=24
  IPACK=((IIX*IJY*INBITS/8)+(2101*2))
  !
  ISEC4(1)=IIX*IJY                   ! number of data to be packed
  ISEC4(2)=INBITS                      ! number of bits used for each value
  ISEC4(3)=0                         ! 0 since section 2 is present
  ISEC4(4)=0                         ! simple packing
  ISEC4(5)=0                         ! floating point data representation
  ISEC4(6:42)=0 
  ! 
  ALLOCATE(ZSEC4(IIX*IJY))
  ZSEC4(:)=0.
  !
  DO JK=kkdeb,kkfin
    IF (OVAR2D) THEN
       IF (PRESENT(KLEVEL2D)) THEN
         ISEC1(7)=105 ! type of level  : altitude
         ISEC1(8)= KLEVEL2D! value of level
         ISEC1(9)=0 ! bottom level if layer
         ISEC2(12)=1
         ZSEC2(11)=KLEVEL2D       
       ELSE        
         ISEC1(7)=105 ! type of level  : altitude
         ISEC1(8)=XZHAT(2) ! value of level
         ISEC1(9)=0 ! bottom level if layer
         ISEC2(12)=1
         ZSEC2(11)=XZHAT(2) 
       ENDIF       
    ELSE
      IF (HTYPEOUT(1:1) == 'P') THEN
        ISEC1(7)=100 ! type of level  : isobaric surfac
        ISEC1(8)=NINT(PLEV(JK)) ! value of level
        ISEC1(9)=0 ! bottom level if layer 
      ELSEIF (HTYPEOUT(1:1) == 'A') THEN
        ISEC1(7)=105 ! type of level  : isobaric surfac
        ISEC1(8)=NINT(PLEV(JK)) ! value of level
        ISEC1(9)=0 ! bottom level if layer 
      ELSE     ! code as height levels
        ISEC1(7)=103 ! type of level  : altitude
        ISEC1(8)=NINT(PLEV(JK)) ! value of level
        ISEC1(9)=0 ! bottom level if layer     
      ENDIF
    ENDIF

        ! 
    ZSEC4(1:IIX*IJY)=RESHAPE(ZFIELD(:,:,JK),(/IIX*IJY/))
    IF (NVERB>=10)  CALL GRSDBG(1)   ! switch ON(1)/OFF(0) debug printing
    CALL GRSDBG(0)                  ! pas de redirection possible...
    YOPER = 'C'      !  for coding 
    IERR = 1
    INBUFF(:)=0.
    CALL GRIBEX (ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4, &
                 ZSEC4,IPUNP,INBUFF,IPACK,IWORD,YOPER,IERR)
    print'(A,I3,A,I5,A,I7,A,I7)', 'FIELD= ',KCODCOD,' LEVEL= ',NINT(PLEV(JK)), &
                                 ' IPACK= ',IPACK,' IWORD= ',IWORD
       
    CALL PBWRITE(IGRIBFILE,INBUFF,ISEC0(1),IERR)
    print*, 'in unit IGRIBFILE=',IGRIBFILE, &
           ' number of bytes: ',IERR
    IF (IERR < 0) THEN
      print*, 'ERROR when writing in GRIB file:  IERR=',IERR
      STOP
    ENDIF
  END DO
  DEALLOCATE(ZFIELD)
  DEALLOCATE(ZSEC4)
  DEALLOCATE(INBUFF)
ENDIF
!-------------------------------------------------------------------------------
!
!*       4.    FERMETURE des fichiers de sortie
!              --------------------------------
!
IF ( HFLAGFILE(1:3) == 'CLO' ) THEN
   print*,'WRITEGRIB: avant fermeture fichier de sortie ',YFILEOUT
   CALL PBCLOSE(IGRIBFILE,IRESP)
   print*, 'After close of ',YFILEOUT,' IRESP=',IRESP
ENDIF


END SUBROUTINE WRITEGRIB                                      
