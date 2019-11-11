!     #################################
      MODULE MODI_WRITECDL
!     #################################
INTERFACE WRITECDL
      SUBROUTINE  WRITECDL(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
                   kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin,&
                   HLABELCHAMP,HFILENAME,HFLAGFILE,HTYPEGRID,  &
                   HFILENAME_SUP,KVERBIA,KRETCODE,             &
                   PGRIDX,PGRIDY                               )
!
CHARACTER(LEN=*), intent(inout)  :: HLABELCHAMP         ! nom du champ
                       ! inout pour modifier le nom VLEV en altitude
CHARACTER(LEN=*), intent(in)  ::  HFILENAME             ! nom du fichier
CHARACTER(LEN=*), intent(in)  :: HFLAGFILE              ! NEW=creation 
                                                        ! OLD=ajout 
                                                        ! CLOSE=fermeture
CHARACTER(LEN=3)              :: HFILENAME_SUP          ! chaine de caracteres
                                                        ! a rajouter a
                                                        ! HFILENAME
CHARACTER(LEN=*), intent(in) :: HTYPEGRID               !  format grille reguliere plan conforme
                                                        ! ou lat lon CONF/LALO
INTEGER , intent(in)         :: KVERBIA                 ! prints de controle
                                      ! desactive (0) / active (1) les prints
                                      ! limites sur les 6 dimensions
INTEGER , intent(in)         :: kideb,kifin,kjdeb,kjfin,kkdeb,kkfin   
INTEGER , intent(in)         :: kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin
! 
INTEGER , intent(out)        :: KRETCODE   ! Code de retour de la routine 
REAL, DIMENSION(:), INTENT(IN) :: PGRIDX, PGRIDY
END SUBROUTINE
END INTERFACE
END MODULE MODI_WRITECDL
!
!     ################
      SUBROUTINE  WRITECDL(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
                   kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin,&
                   HLABELCHAMP,HFILENAME,HFLAGFILE,HTYPEGRID,  &
                   HFILENAME_SUP,KVERBIA,KRETCODE,             &
                   PGRIDX,PGRIDY                               )
!     ################
!
!!****  *writedcdl* - 
!! 
!!
!!    PURPOSE
!!    -------
!     Ecriture d'un fichier  de type CDL pour etre transformé en netcdf
!     via ncgen -b file.cdl
! 
!
!!**  METHOD
!!    ------
!   Ecriture ascii de 2 fichiers en parallele:
! un fichier pour l entete
! un fichier pour les données
!   Chaque appel de la routine writecdl complete le fichier d entete
! et le fichier de données.
!   Ces 2 fichiers seront concatenes avant d'appeler ncgen ( outil netcdf
! qui cree un fichier netcdf a partir d un fichier ascii de format CDL).
!   Voir le script tonetcdf  ci-dessous:
!# concatenation de l entete et des données
!# 
!cat ${FILE}hcl ${FILE}dcl > ${FILE}cdl
!#
!# outil netcdf : ncgen 
!#
!ncgen -b ${FILE}cdl      
!
!     XVAR est alloué avant l appel a writecdl
!
!     HFLAGFILE='NEW' lors de la premiere utilisation du fichier
!     HFLAGFILE='OLD' lors des utilisations suivantes
!     HFLAGFILE='CLO' pour la fermeture du fichier de sortie
!      ( fin de mise a jour du menu )
!
!     KVERBIA= 0 impressions reduites au minimum (entree et sortie de la
!      routine)
!     KVERBIA >0 impressions pour signaler chaque etape de READVAR
!
!     KRETCODE = 0 execution de writecdl correcte
!     KRETCODE = 1 erreur lors de l ouverture du fichier
!     KRETCODE = 2 erreur lors de la fermeture du fichier
!
!     kideb,kifin,kjdeb,kjfin,kkdeb,kkfin = limites en indices i,j,k du
!                                                   domaine à traiter dans XVAR
!     kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin = limites en indices
!                                                  des dimensions 4,5,6 de XVAR
!!      
!!
!!    EXTERNAL
!!    --------
!!
!!          FROM_COMPUTING_UNITS: retour aux unites initiales  avant ecriture
!!                               = passage inverse a celui realise par
!!                                 TO_COMPUTING_UNITS
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
!!     N. Asencio * CNRM*
!!
!!    Copyright 2003,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!!    23/06/2009 G. TANGUY * CNRM*
!! ajout du champ _Fillvalue pour les valeurs indéfinies
!! modification de l'ecriture de "time" : type int et la référence est prise au
!! premier janvier deux ans auparavant
!! ecriture de la dimension de vertical_levels quand il n'y a qu'un seul niveau
!! demandé
!! ajout de la variable YNETCDFCHAMP pour remplacer HLABELCHAMP dans ce
!! programme ce qui évite de tronquer vertical_levels
!! ajout du champ global attributes pour préciser la simulation dans l'entête
!! 18/02/2010 : time doit etre ecrit en premier puisqu'il est UNLIMITED
!!              changement de l'ordre avec le mask
!! Nov 2010 : ajout des paramètres de cartes (LON0,LAT0,LONOR,LATOR,RPK,BETA)
!!            pour les projections conformes (utile sous NCL pour retracer la carte)
!!            Passage des coordonnées en metres au lieu de km (coord conformes
!!            et niveaux verticaux)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
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
USE MODI_TEMPORAL_DIST_FOR_EXT
USE MODI_FROM_COMPUTING_UNITS
USE MODD_CONF, ONLY: CEXP
USE MODD_TIME, ONLY: TDTEXP,TDTSEG
USE MODD_TIME1, ONLY: TDTCUR
USE MODD_GRID

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
CHARACTER(LEN=3)              :: HFILENAME_SUP          ! chaine de caracteres
                                                        !a rajouter a HFILENAME
CHARACTER(LEN=*), intent(in) :: HTYPEGRID               ! format grille reguliere plan conforme
                                                        !ou lat lon CONF/LALO
INTEGER , intent(in)         :: KVERBIA                 ! prints de controle
                                      !desactive (0) / active (1) les prints
                                      ! limites sur les 6 dimensions
INTEGER , intent(in)         :: kideb,kifin,kjdeb,kjfin,kkdeb,kkfin   
INTEGER , intent(in)         :: kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin
! 
INTEGER , intent(out)        :: KRETCODE   ! Code de retour de la routine 
REAL, DIMENSION(:), INTENT(IN) :: PGRIDX, PGRIDY
!
!*       0.2   Variables locales
!              -----------------
!
INTEGER           :: ILOOP,JLOOP,KLOOP,KLOOP4,KLOOP5,KLOOP6, iret
INTEGER,save      :: ILUOUT1HEAD,ILUOUT2DATA  ! unites logiques de sortie 
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
CHARACTER (LEN=28)  :: YFILEOUT,YFILEOUT1,YFILEOUT2   ! Fichiers de sortie
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
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
! 
  IAN=XDATIME(13,1)
  IMOIS=XDATIME(14,1)
  IJOUR=XDATIME(15,1)
  ISECONDE=XDATIME(16,1)
  IANREF=IAN-2
!
YNETCDFCHAMP=HLABELCHAMP
if (KVERBIA >= 0) then
   print *,' --------- '
   print *,'Entree WRITECDL ',TRIM(HFILENAME),' ',TRIM(YNETCDFCHAMP),' ', &
                              TRIM(HFLAGFILE),' ',TRIM(HTYPEGRID),' ', &
                              TRIM(HFILENAME_SUP),' ',KVERBIA
endif
!
! Code de retour de la routine : 0 = OK
!                                1 = erreur lors de l ouverture du fichier
!                                2 = erreur lors de la fermeture du fichier
KRETCODE=0
!
!  Retour aux unites initiales si necessaire
CALL FROM_COMPUTING_UNITS(YNETCDFCHAMP,CUNITE(1)) 
!
!
! code de retour d erreur des routines diaprog
LPBREAD=.FALSE.                                                        
!
if (KVERBIA > 0) then
  print'(A41,6(I4,X))','WRITECDL: ideb,ifin,jdeb,jfin,kdeb,kfin= ',&
          kideb,kifin,kjdeb,kjfin,kkdeb,kkfin
  print'(A42,2(I10,X),4(I4,X))','          tdeb,tfin,trdeb,trfin,pdeb,pfin= ',&
          kitdeb,kitfin,kitrdeb,kitrfin,kipdeb,kipfin 
  print'(A26,6(I4,X))','  nil,nih,njl,njh,nkl,nkh=',nil,nih,njl,njh,nkl,nkh
endif
!
!*       1.1    nom des fichiers de sortie (ajout d un suffixe hkcl/dkcl
!                                                           ou hzcl/dzcl)
!
YFILEOUT=ADJUSTL(ADJUSTR(HFILENAME(1:LEN(HFILENAME)-1))//HFILENAME_SUP)
YFILEOUT1=ADJUSTL(ADJUSTR(HFILENAME(1:LEN(HFILENAME)-1))//'h'//&
          ADJUSTL(HFILENAME_SUP))
YFILEOUT2=ADJUSTL(ADJUSTR(HFILENAME(1:LEN(HFILENAME)-1))//'d'//&
          ADJUSTL(HFILENAME_SUP))
if (KVERBIA > 0) then
    print*,'fichier d entete   YFILEOUT1= ',YFILEOUT1
    print*,'fichier de donnees YFILEOUT2= ',YFILEOUT2
endif
!   
!-------------------------------------------------------------------------------
!
!*       2.1   OUVERTURE DES FICHIERS DE SORTIE
!              -------------------
!
IF ( HFLAGFILE(1:3) == 'NEW' ) THEN
  !
  ! recupere l unite logique et ouverture des fichiers
  !
  !*    2.1.1 Fichier entete : partie commune a toutes les variables
  !           --------------
  CALL FMATTR(YFILEOUT1,CLUOUTDIAS(NBFILES),ILUOUT1HEAD,NRESPDIAS(NBFILES))
  IF (NRESPDIAS(NBFILES).NE.0)THEN
    KRETCODE=1
    print *,' ****WRITECDL: erreur lors de l ouverture du fichier ',&
            TRIM(YFILEOUT1),' code= ',NRESPDIAS(NBFILES)
    RETURN
  ENDIF
  OPEN(UNIT=ILUOUT1HEAD,FILE=YFILEOUT1,STATUS='NEW',FORM='FORMATTED')
  ! creation du debut de l entete
  !nom du fichier
  write(ILUOUT1HEAD,*) 'netcdf ',YFILEOUT,' { '
  !dimensions
  write(ILUOUT1HEAD,*) 'dimensions: '
  SELECT CASE (HTYPEGRID(1:4) ) 
  CASE ('CONF')
     YLIBELLEDIM1='W_E_direction'
     YLIBELLEDIM2='S_N_direction'
  CASE ('LALO')
     YLIBELLEDIM1='longitude'
     YLIBELLEDIM2='latitude'
  CASE DEFAULT
     print*, ' type de grille incorrect: LALO/CONF possibles et non ', HTYPEGRID
  END SELECT
  !
  write(ILUOUT1HEAD,*) '   ',TRIM(YLIBELLEDIM1),'  = ', kifin-kideb +1, ';'
  write(ILUOUT1HEAD,*) '   ',TRIM(YLIBELLEDIM2),'  = ', kjfin-kjdeb +1, ';'
  write(ILUOUT1HEAD,*) '   vertical_levels   = ', kkfin-kkdeb +1, ';'
!  write(ILUOUT1HEAD,*) '       time   =  ',kitfin-kitdeb +1, ';'
  write(ILUOUT1HEAD,*) '   time   =  UNLIMITED ; // (',kitfin-kitdeb +1,' currently) ;'
   write(ILUOUT1HEAD,*) '   mask = ', kitrfin-kitrdeb +1, ';'
  write(ILUOUT1HEAD,*) 'variables: '
  
!  write (ILUOUT1HEAD,*) '        double time(time);'
  write (ILUOUT1HEAD,*) '        int time(time);'
  write(ILUOUT1HEAD,'(A,I4,A)') ' time:units = "seconds since ',IANREF,'-1-1 00:00:00" ;'
  write(ILUOUT1HEAD,'(A,I4,A)') ' time:time_origin = "',IANREF,'-1-1 00:00:00" ;'

  !reference temporelle
!  write (ILUOUT1HEAD,*) '        int base_time ;'
!  write (ILUOUT1HEAD,*)' base_time:units = "seconds since 1970-01-01'&
!                            ,'00:00:00 UTC" ;'
!  write (ILUOUT1HEAD,*) ' base_time:long_name = ',&
!                       '"base time for the file" ;'
  !evolution temporelle / reference
!  write (ILUOUT1HEAD,*) '        int time_offset(time) ;'
!  write (ILUOUT1HEAD,*)' time_offset:units = "seconds" ;'
!  write (ILUOUT1HEAD,*) ' time_offset:long_name = "time offset from'&
!                       ,' base time" ;'
  SELECT CASE (HTYPEGRID(1:4) ) 
  CASE ('CONF')
    !grille réguliere selon x dans le plan conforme
    write (ILUOUT1HEAD,*) '        float W_E_direction(W_E_direction);'
    write (ILUOUT1HEAD,*) '   W_E_direction:units = "km" ;'
    write (ILUOUT1HEAD,*) '   W_E_direction:long_name = "model grid in the conformal projection" ;'
    !grille réguliere selon y dans le plan conforme
    write (ILUOUT1HEAD,*) '        float S_N_direction(S_N_direction);'
    write (ILUOUT1HEAD,*) '   S_N_direction:units = "km" ;'
    write (ILUOUT1HEAD,*) '   S_N_direction:long_name = "model grid in the conformal projection" ;'
    write (ILUOUT1HEAD,*) '        float LON0 ;'
    write (ILUOUT1HEAD,*) '   LON0:units = "degrees_east" ;'
    write (ILUOUT1HEAD,*) '   LON0:long_name = "reference longitude for the conformal projection" ;'
    write (ILUOUT1HEAD,*) '        float LAT0 ;'
    write (ILUOUT1HEAD,*) '   LAT0:units = "degrees_north" ;'
    write (ILUOUT1HEAD,*) '   LAT0:long_name = "reference latitude for the conformal projection" ;'
    write (ILUOUT1HEAD,*) '        float LONOR ;'
    write (ILUOUT1HEAD,*) '   LONOR:units = "degrees_east" ;'
    write (ILUOUT1HEAD,*) '   LONOR:long_name = "longitude of point x=0,y=0 in the conformal projection" ;'
    write (ILUOUT1HEAD,*) '        float LATOR ;'
    write (ILUOUT1HEAD,*) '   LATOR:units = "degrees_north" ;'
    write (ILUOUT1HEAD,*) '   LATOR:long_name = "latitude of point x=0,y=0 in  the conformal projection" ;'
    write (ILUOUT1HEAD,*) '        float BETA ;'
    write (ILUOUT1HEAD,*) '   BETA:units = "degrees" ;'
    write (ILUOUT1HEAD,*) '   BETA:long_name = "Rotation angle for the conformal projection" ;'
    write (ILUOUT1HEAD,*) '        float RPK ;'
    write (ILUOUT1HEAD,*) '   RPK:units = " " ;'
    write (ILUOUT1HEAD,*) '   RPK:long_name = "projection parameter for the conformal projection" ;'

  CASE('LALO')
    !grille réguliere selon x en longitude
    write (ILUOUT1HEAD,*) '        float longitude(longitude);'
    write (ILUOUT1HEAD,*) '   longitude:units = "degrees_east" ;'
    write (ILUOUT1HEAD,*) '   longitude:long_name = "longitudes" ;'
    !grille réguliere selon y en latitude
    write (ILUOUT1HEAD,*) '        float latitude(latitude);'
    write (ILUOUT1HEAD,*) '   latitude:units = "degrees_north" ;'
    write (ILUOUT1HEAD,*) '   latitude:long_name = "latitudes" ;'
  END SELECT
  !
  !*    2.1.2 Fichier contenant les donnees: variables contenant la grille
  !           ------------------------------ 
  CALL FMATTR(YFILEOUT2,CLUOUTDIAS(NBFILES),ILUOUT2DATA,NRESPDIAS(NBFILES))
  IF (NRESPDIAS(NBFILES).NE.0)THEN
    KRETCODE=1
    print *,' ****WRITECDL: erreur lors de l ouverture du fichier ',&
            TRIM(YFILEOUT2),' code= ',NRESPDIAS(NBFILES)
    RETURN
  ENDIF
  OPEN(UNIT=ILUOUT2DATA,FILE=YFILEOUT2,STATUS='NEW',FORM='FORMATTED')
  !
  !calcul  et ecriture du nombre de secondes depuis le 01/01 2 ans auparavant
   zbasetime=0.
  if (KVERBIA > 0) then
    print *,' calcul ibasetime: IAN,IMOIS,IJOUR,ISECONDE,zbasetime'
    print *,IAN,IMOIS,IJOUR,ISECONDE,zbasetime
  endif
   CALL TEMPORAL_DIST_FOR_EXT(IAN,IMOIS,IJOUR,ISECONDE,IANREF,01,01,0,zbasetime)
  if (KVERBIA > 0) then
     print *, IAN,IMOIS,IJOUR,ISECONDE,zbasetime
  endif
  !
  ibasetime=zbasetime
  write(ILUOUT2DATA,*) 'data: '
  write(ILUOUT2DATA,*) 'time = '!,zbasetime !, ' ;'

 ! write(ILUOUT2DATA,*) 'base_time = ',ibasetime, ' ;'
  !ecriture de l instant du fichier= 0 seconde / reference
  !write(ILUOUT2DATA,*) 'time_offset = 0 ;'
!  ytextdim='time_offset = '
!  write(ILUOUT2DATA,*) ytextdim
  ALLOCATE(ioffset_time(kitfin-kitdeb+1)) ; ioffset_time(:)=0
  DO JK=kitdeb,kitfin
    !ibasetime=XTRAJT(JK,1)-XTRAJT(kitdeb,1)  !
    ! cas ou TEXP et TSEG sont faux
    IAN=XDATIME(13,kitdeb)
    IMOIS=XDATIME(14,kitdeb)
    IJOUR=XDATIME(15,kitdeb)
    ISECONDE=XDATIME(16,kitdeb)
    IAN2=XDATIME(13,JK)
    IMOIS2=XDATIME(14,JK)
    IJOUR2=XDATIME(15,JK)
    ISECONDE2=XDATIME(16,JK)
    CALL TEMPORAL_DIST_FOR_EXT(IAN2,IMOIS2,IJOUR2,ISECONDE2,IAN,IMOIS,IJOUR,ISECONDE,zbasetime)
    ioffset_time(jk-kitdeb+1)=ibasetime+zbasetime
  ENDDO
  write(ILUOUT2DATA,1010,advance='no') ioffset_time(1:kitfin-kitdeb+1)
  DEALLOCATE(ioffset_time)
  WRITE(ILUOUT2DATA,'(";")')
  write(ILUOUT2DATA,*) ' '
!------------------------------------------------------------------
  SELECT CASE (HTYPEGRID(1:4) ) 
  CASE ('CONF')
    ! grille régulière selon X en km
    write(ILUOUT2DATA,*) ' W_E_direction ='
    write(ILUOUT2DATA,1000,advance='no') PGRIDX(kideb:kifin)*0.001
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    ! grille régulière selon Y en km
    write(ILUOUT2DATA,*) ' S_N_direction ='
    write(ILUOUT2DATA,1000,advance='no') PGRIDY(kjdeb:kjfin)*0.001
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    !parametre de la grille
    write(ILUOUT2DATA,*) ' LON0 ='
    write(ILUOUT2DATA,1000,advance='no') XLON0
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    write(ILUOUT2DATA,*) ' LAT0 ='
    write(ILUOUT2DATA,1000,advance='no') XLAT0
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    write(ILUOUT2DATA,*) ' LONOR ='
    write(ILUOUT2DATA,1000,advance='no') XLONORI
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    write(ILUOUT2DATA,*) ' LATOR ='
    write(ILUOUT2DATA,1000,advance='no') XLATORI
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    write(ILUOUT2DATA,*) ' BETA ='
    write(ILUOUT2DATA,1000,advance='no') XBETA
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    write(ILUOUT2DATA,*) ' RPK ='
    write(ILUOUT2DATA,1000,advance='no') XRPK
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '


  CASE('LALO')
    write(ILUOUT2DATA,*) 'longitude ='
    write(ILUOUT2DATA,1000,advance='no') PGRIDX(kideb:kifin)
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
    write(ILUOUT2DATA,*) 'latitude ='
    write(ILUOUT2DATA,1000,advance='no') PGRIDY(kjdeb:kjfin)
    WRITE(ILUOUT2DATA,'(";")')
    write(ILUOUT2DATA,*) ' '
  END SELECT
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    ECRITURE du champ dans YFILEOUT2 et de l entete dans YFILEOUT1
!              --------
!
IF ( HFLAGFILE(1:3) /= 'CLO' ) THEN
  !
  if (KVERBIA > 0) then
    print*,'WRITECDL: format CDL ecriture en cours '
  endif
  ! 
  ! Ecriture du champ + lat,lon ,altitude du niveau
  ! 
  !      3.1  liste des dimensions tel que "Last dim varies fastest"
  ! 
  ytextdim=''
  !Process: ecriture d une variable netcdf par processus donc lignes commentees
  !IF ( kipfin-kipdeb > 0) THEN
  !   ytextdim='process '
  !ENDIF
  ! ATTENTION le TEMPS DOIT ETRE LA PREMIERE VARIABLE CAR UNLIMITED  
  !Time
  SELECT CASE (YNETCDFCHAMP)
  CASE ('VLEV')
    if (KVERBIA >= 2) then
      print*,' No temporal dimension for ', YNETCDFCHAMP
    endif
    IF ( SIZE(XVAR,2) > 1 ) THEN
      ! cas du champ 3D pour les altitudes
      ! passage en km pour utilisation Zebra
      YNETCDFCHAMP='VLEV'
      CUNITE(1)='km'
      XVAR=XVAR*0.001
    ELSE
      ! cas d une liste de niveaux verticaux choisis par l utilisateur
      ! on garde l unité donnée par extractdia metres ou hPa
      YNETCDFCHAMP='vertical_levels'
     ENDIF
  CASE ('LAT','LON')
    if (KVERBIA >= 2) then
      print*,' No temporal dimension for ', YNETCDFCHAMP
    endif
  CASE DEFAULT
    ! Les variables doivent avoir la dimension time meme si
    ! cette dimension est egale a 1
    !IF ( kitfin-kitdeb > 0 ) THEN
    IF (ytextdim /= '') ytextdim=ADJUSTL(ADJUSTR(ytextdim)//',')
    ytextdim=ADJUSTL(ADJUSTR(ytextdim)//'time ')
    !ENDIF
  END SELECT

  !Mask
  SELECT CASE (YNETCDFCHAMP)
  CASE ('VLEV','LAT','LON')
  CASE DEFAULT
    IF ( kitrfin-kitrdeb > 0) THEN
      IF (ytextdim /= '') ytextdim=ADJUSTL(ADJUSTR(ytextdim)//',')
      ytextdim=ADJUSTL(ADJUSTR(ytextdim)//'mask ')
    ENDIF      
  END SELECT

  !Z
  SELECT CASE (YNETCDFCHAMP)
  CASE ('LAT','LON')
    if (KVERBIA >= 2) then
      print*,' No vertical dimension for ', YNETCDFCHAMP
    endif
  CASE ('vertical_levels')
    IF (ytextdim /= '') ytextdim=ADJUSTL(ADJUSTR(ytextdim)//',')
    ytextdim=ADJUSTL(ADJUSTR(ytextdim)//'vertical_levels ')
  CASE DEFAULT
    IF ( kkfin-kkdeb > 0) THEN
      IF (ytextdim /= '') ytextdim=ADJUSTL(ADJUSTR(ytextdim)//',')
      ytextdim=ADJUSTL(ADJUSTR(ytextdim)//'vertical_levels ')
    ENDIF
  END SELECT
  !Y
  IF ( kjfin-kjdeb > 0) THEN
    IF (ytextdim /= '') ytextdim=ADJUSTL(ADJUSTR(ytextdim)//',')
    ytextdim=ADJUSTL(ADJUSTR(ytextdim)//ADJUSTL(YLIBELLEDIM2))
  ENDIF
  !X
  IF ( kifin-kideb > 0) THEN
    IF (ytextdim /= '') ytextdim=ADJUSTL(ADJUSTR(ytextdim)//',')
    ytextdim=ADJUSTL(ADJUSTR(ytextdim)//ADJUSTL(YLIBELLEDIM1))
  ENDIF
  !
  if (KVERBIA >= 2) then
    print *,' dimensions du tableau= ', TRIM(ytextdim)
  end if
  !
 ! Ecriture d une variable netcdf par processus
 ! nommée nom_var+pnum_process
 DO  KLOOP6=kipdeb,kipfin
  YLABELCHAMPnew=ADJUSTL(YNETCDFCHAMP)
  IF ( SIZE(XVAR,6)  > 1  ) THEN
    ! ajout du numéro de processus
    WRITE (YNUM,'(I5)') KLOOP6
    YLABELCHAMPnew=ADJUSTL(ADJUSTR(YNETCDFCHAMP)//'p'//ADJUSTL(YNUM))
  ENDIF
  write (ILUOUT1HEAD,*) '        float ',TRIM(YLABELCHAMPnew),'(',TRIM(ytextdim),') ;'
  write (ILUOUT1HEAD,*) TRIM(YLABELCHAMPnew), ':long_name = "',TRIM(CTITRE(kloop6)),'" ;'
  write (ILUOUT1HEAD,*) TRIM(YLABELCHAMPnew), ':units = "',TRIM(CUNITE(kloop6)),'" ;'
  SELECT CASE (YNETCDFCHAMP)
  CASE ('LAT','LON')
    ikdeb=1 ; ikfin=1 ; iitdeb=1 ; iitfin=1 ; iitrdeb=1 ; iitrfin=1
  CASE DEFAULT
    ikdeb=kkdeb ; ikfin=kkfin ; iitdeb=kitdeb ; iitfin=kitfin ; iitrdeb=kitrdeb ; iitrfin=kitrfin
  END SELECT
  IF (ANY(XVAR(kideb:kifin,kjdeb:kjfin,ikdeb:ikfin,iitdeb:iitfin,iitrdeb:iitrfin,kloop6)/=XSPVAL)) THEN
    zmini=MINVAL(XVAR(kideb:kifin,kjdeb:kjfin,ikdeb:ikfin, &
                      iitdeb:iitfin,iitrdeb:iitrfin,kloop6), &
                 MASK=XVAR(kideb:kifin,kjdeb:kjfin,ikdeb:ikfin, &
                           iitdeb:iitfin,iitrdeb:iitrfin,kloop6)/=XSPVAL     )
    zmaxi=MAXVAL(XVAR(kideb:kifin,kjdeb:kjfin,ikdeb:ikfin, &
                      iitdeb:iitfin,iitrdeb:iitrfin,kloop6), &
                 MASK=XVAR(kideb:kifin,kjdeb:kjfin,ikdeb:ikfin, &
                           iitdeb:iitfin,iitrdeb:iitrfin,kloop6)/=XSPVAL     )
  ELSE
    zmini=XSPVAL ; zmaxi=XSPVAL
  ENDIF
  IF (ABS (zmini) > 1.E-05 .AND. ABS(zmaxi) > 1.E-05 ) THEN
    write (ILUOUT1HEAD,FMT=101) TRIM(YLABELCHAMPnew),zmini,zmaxi
  ELSE
    write (ILUOUT1HEAD,FMT=103) TRIM(YLABELCHAMPnew),zmini,zmaxi
  ENDIF
  IF (YNETCDFCHAMP /= 'vertical_levels') THEN
    write (ILUOUT1HEAD,FMT=102) TRIM(YLABELCHAMPnew),XSPVAL
    write (ILUOUT1HEAD,FMT=104) TRIM(YLABELCHAMPnew),XSPVAL
  ENDIF
  ! 
  !      3.2 ecriture des valeurs: Last dim varies fastest
  ! 
! on intervertit la place du temps et la place du mask avant l'ecriture

ALLOCATE(XVAR2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,5),SIZE(XVAR,4)))

DO II=kideb,kifin
  DO IJ=kjdeb,kjfin
    DO IK=ikdeb,ikfin
      DO IT=iitdeb,iitfin 
        DO IM=iitrdeb,iitrfin
          XVAR2(II,IJ,IK,IM,IT)=XVAR(II,IJ,IK,IT,IM,kloop6)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO



  write(ILUOUT2DATA,*) TRIM(YLABELCHAMPnew),' = '
  IF (ABS (zmini) > 1.E-04 .AND. ABS(zmaxi) > 1.E-04 ) THEN
    WRITE(ILUOUT2DATA,FMT=1000,advance='no') XVAR2(kideb:kifin,kjdeb:kjfin,&
                   ikdeb:ikfin,iitrdeb:iitrfin,iitdeb:iitfin)
  ELSE
    WRITE(ILUOUT2DATA,FMT=1001,advance='no') XVAR2(kideb:kifin,kjdeb:kjfin,&
                   ikdeb:ikfin,iitrdeb:iitrfin,iitdeb:iitfin)        
  ENDIF
DEALLOCATE(XVAR2)
  WRITE(ILUOUT2DATA,'(";")')
  write(ILUOUT2DATA,*) ' '
 END DO

  !
101   FORMAT (1H ,A,16H :actual_range = ,F0.5,3Hf ,,F0.5,3Hf ;) 
103   FORMAT (1H ,A,16H :actual_range = ,E11.5,3Hf ,,E11.5,3Hf ;) 
102   FORMAT (1H ,A,18H :missing_value = ,F0.5,3Hf ;)
104   FORMAT (1H ,A,15H :_FillValue = ,F0.5,3Hf ;)
!105   FORMAT (8H time = ,E17.11,3Hf ;)

  ! le ":" est le descripteur de fin d'exploitation d'un format. 
  ! sous f95 et pgf90. D. Gazen
1000  FORMAT (7(F0.5,:,", "))      
1001  FORMAT (7(E11.5,:,", "))      
1010  FORMAT (7(I10,:,", "))      
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.    FERMETURE des fichiers de sortie
!              --------------------------------
!
IF ( HFLAGFILE(1:3) == 'CLO' ) THEN
  ! fin de fichier de données
  WRITE(ILUOUT2DATA,*) '}'  
  if (KVERBIA > 0) then
    print*,'WRITECDL: avant fermeture fichier de sortie ',YFILEOUT
  endif
  ! force les buffers a etre vides pour permettre a l appel
  ! systeme de traiter les fichiers complets
  !CALL FLUSH (ILUOUT1HEAD)
  !CALL FLUSH (ILUOUT2DATA)
  !
  ! fermeture
    write (ILUOUT1HEAD,*)  "// global attributes:"
    write (ILUOUT1HEAD,*)  '  :title = "Meso-NH simulation" ;'
    write (ILUOUT1HEAD,*)  '  :grid_resolution_in_meters = "',  XXDXHAT(1,1),' x ',XXDYHAT(1,1),'" ;'
    write (ILUOUT1HEAD,*)  '  :description = "Data are from the file ', HFILENAME, '" ;'
    write (ILUOUT1HEAD,'(A46,3(I4,X),F12.4,A25,3(I4,X),F12.4,A3)')&
    '  :comments = " Meso-NH  experience starts at ',TDTEXP,' and segment starts at ', TDTSEG,' ";'
    CALL DATE_AND_TIME(YDATE, YTIME, YZONE, IVALUES)
    write (ILUOUT1HEAD,FMT=201) IVALUES(3),IVALUES(2),IVALUES(1),IVALUES(5),IVALUES(6),IVALUES(7)
201   FORMAT ('   :history = "created on  ',I2,'/',I2,'/',I4, ' at ',I2,':',I2,':',I2,'" ;') 


  CLOSE(ILUOUT1HEAD)
  CALL FMFREE(YFILEOUT1,CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
  IF (NRESPDIAS(NBFILES).NE.0)THEN
    KRETCODE=2
    print *,' ****WRITECDL: erreur lors de la fermeture du fichier ',&
            TRIM(YFILEOUT1),' code= ',NRESPDIAS(NBFILES)
  ENDIF
  CLOSE(ILUOUT2DATA)
  CALL FMFREE(YFILEOUT2,CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
  IF (NRESPDIAS(NBFILES).NE.0)THEN
    KRETCODE=2
    print *,' ****WRITECDL: erreur lors de la fermeture du fichier ',&
            TRIM(YFILEOUT2),' code= ',NRESPDIAS(NBFILES)
  ENDIF
  !
  if (KVERBIA > 0) then
    print *,'WRITECDL: before calling tonetcdf'
  end if
  ycommand='tonetcdf '//ADJUSTL(ADJUSTR(HFILENAME))
  call SYSTEM ( TRIM(ycommand) )
  !
  if (KVERBIA >= 0) then
    print*,'Sortie WRITECDL: Fichier ',TRIM(YFILEOUT),' disponible au format cdl'
    print*,' --------- '
  endif
  !
ENDIF
!
!
HLABELCHAMP=YNETCDFCHAMP

END SUBROUTINE WRITECDL
