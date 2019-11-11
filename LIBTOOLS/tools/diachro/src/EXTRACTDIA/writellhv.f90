!     #################################
      MODULE MODI_WRITELLHV
!     #################################
INTERFACE WRITELLHV     
      SUBROUTINE WRITELLHV(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,     &
                           KTDEB,KTFIN,KTRDEB,KTRFIN,KPDEB,KPFIN,   &
                           HLABELCHAMP,HFILENAME,HFLAGFILE,HTYPEOUT,&
                           KVERBIA,KRETCODE,HFILENAME_SUP,PLON,PLAT,PALT )
!
CHARACTER(LEN=*), INTENT(in) :: HLABELCHAMP,HFILENAME ! nom du champ et du fichier
CHARACTER(LEN=*), INTENT(in) :: HFLAGFILE             ! NEW=creation 
                                                      ! OLD=ajout 
                                                      ! CLOSE=fermeture
                                                      ! NEW1H=creation entete speciale
                                                      ! OLDNH= ajout sans entete
CHARACTER(LEN=*), INTENT(in) :: HTYPEOUT              ! type de fichier sortie
                                                      ! LL?V= lon lat alt val
                                                      ! ll?v= lat lon alt val
                                                      !?=H,h alt du niveau k
                                                      !  Z,z alt apres
                                                      !  P,p interpol. verticale
                                                      ! en Z=cst Presssion=cst
INTEGER , INTENT(in)         :: KVERBIA               ! prints de controle
                                      ! desactive (0) / active (1) les prints
                                      ! limites sur les 6 dimensions
INTEGER , INTENT(in)         :: KIDEB,KIFIN,KJDEB,KJFIN,KKDEB,KKFIN
INTEGER , INTENT(in)         :: KTDEB,KTFIN,KTRDEB,KTRFIN,KPDEB,KPFIN
INTEGER , INTENT(out)        :: KRETCODE   ! Code de retour de la routine       
CHARACTER(LEN=3) ,OPTIONAL   :: HFILENAME_SUP    ! chaine de caracteres
                                                 !a rajouter a HFILENAME
REAL, DIMENSION(:,:), INTENT(IN), OPTIONAL   :: PLON,PLAT ! tableaux des lat et
                                                          ! lon si LLZV ou LLPV
REAL, DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: PALT ! tableau des altitudes
END SUBROUTINE       
END INTERFACE
END MODULE MODI_WRITELLHV       
!     ######
      SUBROUTINE WRITELLHV(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,     &
                           KTDEB,KTFIN,KTRDEB,KTRFIN,KPDEB,KPFIN,   &
                           HLABELCHAMP,HFILENAME,HFLAGFILE,HTYPEOUT,&
                           KVERBIA,KRETCODE,HFILENAME_SUP,PLON,PLAT,PALT )
!     ################
!
!!****  *WRITELLHV* - 
!! 
!!
!!    PURPOSE
!!    -------
!     Ecriture d'un fichier de type lon,lat,alt,val (LL) ou lat,lon,alt,val (ll)
!         lon,lat= type LLHV,llhv: position dans la grille modele
!                  type LLZV,llzv/LLPV,llpv: apres interpolation horizontale
!                                                                (PLAT,PLON)
!         alt= type LLHV,llhv: position verticale de la grille du modèle (XZZ)
!                        ou apres interpolation verticale a Z ou P=cst (PALT)
!              type LLZVllzv,/LLPV,llpv: apres interpolation verticale 
!                                         a Z ou P=cst (PALT)
! NB: ces interpolations ont ete realisees avant l'appel de WRITELLHV
! 
!
!!**  METHOD
!!    ------
!     utilisation des routines de diaprog : le tableau de stockage
!     XVAR est alloué avant l appel a WRITELLHV
!
!     HFLAGFILE='NEW' lors de la premiere utilisation du fichier
!     HFLAGFILE='OLD' lors des utilisations suivantes avec nouvelle entete
!     HFLAGFILE='NEW1H' lors de la premiere utilisation du fichier et gestion
!                d une entete speciale (cas mesonh2obs)
!     HFLAGFILE='OLDNH' lors des utilisations suivantes sans nouvelle entete
!                      (cas mesonh2obs)
!     HFLAGFILE='OLD1H' lors des utilisations suivantes du fichier et gestion
!                d une entete speciale (cas mesonh2obs)
!     HFLAGFILE='CLO' pour la fermeture du fichier de sortie
!      ( fin de mise a jour du menu )
!
!     KVERBIA= 0 impressions reduites au minimum (entree et sortie de la
!      routine)
!     KVERBIA >0 impressions pour signaler chaque etape de READVAR
!
!     KRETCODE = 0 execution de WRITELLHV correcte
!     KRETCODE = 1 erreur lors de l ouverture du fichier
!     KRETCODE = 2 erreur lors de la fermeture du fichier 
!
!     kideb,kifin,kjdeb,kjfin,kkdeb,kkfin = limites en indices i,j,k du
!       domaine à traiter dans XVAR       
!     KTDEB,KTFIN,KTRDEB,KTRFIN,KPDEB,KPFIN = limites en indices
!       des dimensions 4,5,6 de XVAR       
!!      
!!    EXTERNAL
!!    --------
!!
!!          FROM_COMPUTING_UNITS: retour aux unites initiales  avant ecriture
!!                               = passage inverse a celui realise par
!!                                 TO_COMPUTING_UNITS
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHORS
!!    -------
!!     N. Asencio  * CNRM*
!!
!!    Copyright 2003,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!     04/11/2009 (G. Tanguy) : add case IJHV,IJZV, IJPV , JIHV, JIZV, JIPV
!     11/07/2014 (G. Tanguy) : correctoin bug IJHV au lieu de JIHV
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! modules MESONH
USE MODD_CST
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT
USE MODE_GRIDPROJ
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_GRID1, ONLY: XZZ,XXHAT,XYHAT
!      
! modules DIACHRO
USE MODN_NCAR,  ONLY: XSPVAL       
!                    XVAR(i,j,k,,,),XMASK,XTRAJT,X,Y,Z,XDATIME(16,t),CUNITE(p)
USE MODD_ALLOC_FORDIACHRO
USE MODD_COORD, ONLY: XXX,XXY,XXZS, & !  XXX(:,1:7), XXY(:,1:7), XXZS(:,:,1:7)
                      XXDXHAT,XXDYHAT ! XXDXHAT(:,1:7), XXDYHAT(:,1:7)
!                    nom de fichiers NLUOUT,CLFIFM, CDESFM
USE MODD_OUT
USE MODD_FILES_DIACHRO, ONLY: NBFILES, CLUOUTDIAS, NRESPDIAS
!                    pour appel a FMATTR et FMCLOS
!USE MODD_DIACHRO, ONLY:CFILEDIA,CLUOUTDIA, &
!                       NLUOUTDIA,NRESPDIA,NNPRARDIA,NFTYPEDIA,NVERBDIA,NNINARDIA
!
!
USE MODI_FROM_COMPUTING_UNITS
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!
CHARACTER(LEN=*), INTENT(IN):: HLABELCHAMP,HFILENAME ! nom du champ et du fichier
CHARACTER(LEN=*), INTENT(IN):: HFLAGFILE             ! NEW=creation 
                                                     ! OLD=ajout 
                                                     ! CLOSE=fermeture
                                                     ! NEW1H=creation entete speciale
                                                     ! OLDNH=ajout  sans entete
CHARACTER(LEN=*), INTENT(IN):: HTYPEOUT              ! type de fichier sortie
                                                     ! LL?V= lon lat alt val
                                                     ! ll?v= lat lon alt val
                                                     !?=H,h alt du niveau k
                                                     !  Z,z alt apres
                                                     !  P,p interpol. verticale
                                                     ! en Z=cst Presssion=cst
INTEGER, INTENT(IN)         :: KVERBIA               ! prints de controle
                                      ! desactive (0) / active (1) les prints
                                                ! limites sur les 6 dimensions
INTEGER, INTENT(IN)         :: KIDEB,KIFIN,KJDEB,KJFIN,KKDEB,KKFIN
INTEGER, INTENT(IN)         :: KTDEB,KTFIN,KTRDEB,KTRFIN,KPDEB,KPFIN
! 
INTEGER , INTENT(OUT)       :: KRETCODE   ! Code de retour de la routine 
CHARACTER(LEN=3) ,OPTIONAL   :: HFILENAME_SUP    ! chaine de caracteres
                                                 !a rajouter a HFILENAME
REAL, DIMENSION(:,:), INTENT(IN), OPTIONAL   :: PLON,PLAT ! tableaux des lat et
                                                          ! lon si LLZV ou LLPV
REAL, DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: PALT ! tableau des altitudes
!
!*       0.2   Declarations des variables locales
!
INTEGER      ::   JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP ! indices de boucle
INTEGER,save ::   ILUOUTLL                        ! unite logique de sortie 
INTEGER      ::   IAN,IMOIS,IJOUR,IHEURE,IMINUTE,ISEC,INBVAL,IGRID
INTEGER      ::   IIU,IJU
! taille= 28 cf routines FM 
CHARACTER (LEN=28)  :: YFILEOUT                        ! Fichier de sortie
REAL   , DIMENSION(:,:)  ,ALLOCATABLE        :: ZLAT,ZLON ! lat et lon
REAL   , DIMENSION(:,:)  ,ALLOCATABLE        :: ZX,ZY
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATION 
!              --------------
!      
if (KVERBIA >= 0) then
  print *,'--------- '
  print *,'Entree WRITELLHV: ',TRIM(HFILENAME),' ',TRIM(HLABELCHAMP),' ', &
                               TRIM(HFLAGFILE),' ',TRIM(HTYPEOUT),' ',KVERBIA
endif
!
! Code de retour de la routine : 0 = OK
!                                1 = erreur lors de l ouverture du fichier
!                                2 = erreur lors de la fermeture du fichier 
KRETCODE=0
!
!  Retour aux unites initiales si necessaire
IF (HFLAGFILE(1:3) /= 'CLO' ) THEN
  IF (HLABELCHAMP/='END') CALL From_Computing_Units(HLABELCHAMP,CUNITE(1)) 
END IF
!
!
! init du zoom
if (KVERBIA > 0 .AND.  HFLAGFILE(1:3) /= 'CLO' ) THEN
  print*,'WRITELLHV: zoom '
  print'(A,6(I4,X))','  ideb,ifin,jdeb,jfin,kdeb,kfin=',&
           kideb,kifin,kjdeb,kjfin,kkdeb,kkfin
  print'(A,2(I8,X),4(I4,X))','  tdeb,tfin,trdeb,trfin,pdeb,pfin=',&
           KTDEB,KTFIN,KTRDEB,KTRFIN,KPDEB,KPFIN
endif
!
!*       1.1   nom du fichier de sortie (ajout d un suffixe LLHV/LLZV/LLPV)
!
SELECT CASE ( HTYPEOUT(1:4) )
 CASE ('LLHV','llhv','LLZV','llzv','LLPV','llpv','jihv','IJHV',&
         'IJZV','jizv','IJPV','jipv','llav','LLAV') 
   YFILEOUT=ADJUSTL(ADJUSTR(HFILENAME(1:LEN(HFILENAME)-1))//HTYPEOUT(1:4))
 CASE DEFAULT
   PRINT*,' ****WRITELLHV: type ', TRIM(HTYPEOUT),' non prevu'
   PRINT*,'types possibles: LLHV/llhv, LLZV/llzv, LLPV/llpv, IJHV/jihv'
   PRINT*,'IJZV/jizv, IJPV/jipv,LLAV/llav'
   KRETCODE=1
   RETURN
END SELECT
IF ( PRESENT(HFILENAME_SUP)) THEN
    IF(HFILENAME_SUP(1:3) /= '  ') THEN
      YFILEOUT=ADJUSTL( ADJUSTR(YFILEOUT)//'_'//ADJUSTL(HFILENAME_SUP) )     
    ENDIF
ENDIF
!
!*       1.2   ouverture du fichier de sortie et allocations
!
IF ( HFLAGFILE(1:3) == 'NEW' ) THEN
  ! recupere l unite logique et ouverture du fichier
  CALL FMATTR(YFILEOUT,CLUOUTDIAS(NBFILES),ILUOUTLL,NRESPDIAS(NBFILES))
  IF (NRESPDIAS(NBFILES)==0 ) THEN
    OPEN(UNIT=ILUOUTLL,FILE=YFILEOUT,STATUS='NEW',FORM='FORMATTED')
  ELSE
    PRINT*,' ****WRITELLHV: error when openning ', TRIM(YFILEOUT), &
           'code= ',NRESPDIAS(NBFILES)
    KRETCODE=1
    RETURN
  ENDIF
ENDIF
!
!*       1.3   test sur les arguments optionnels
!
IF ( HFLAGFILE(1:3) /= 'CLO' ) THEN
IIU=SIZE(XZZ,1) ; IJU=SIZE(XZZ,2)
!
IF (.NOT.PRESENT(PLAT) .AND. .NOT.PRESENT(PLON)) THEN
! utilisation des lat. et lon. de la grille modele
  ALLOCATE(ZX(IIU,IJU),ZY(IIU,IJU))
  ALLOCATE(ZLAT(IIU,IJU),ZLON(IIU,IJU))
  if (KVERBIA>0) print*,'WRITELLHV: LAT et LON de la grille modele '
ELSE ! ( present(PLAT) .or. present(PLON) )
  IF ( (PRESENT(PLAT) .AND. .NOT.PRESENT(PLON)) .OR. &
       (.NOT.PRESENT(PLAT) .AND. PRESENT(PLON)) .OR. &
       .NOT.PRESENT(PALT)                            ) THEN
    PRINT*,' ****WRITELLHV: latitudes ET longitudes doivent etre presentes '
    PRINT*,'               ET altitudes '
    KRETCODE=1
    RETURN
  ENDIF
  ! Cas de passage par argument de PLAT et PLON différents de 
  !ceux de la grille du modele
  IF (PRESENT (PLON)) THEN
    ALLOCATE(ZLON(SIZE(PLON,1),SIZE(PLON,2)))
    ZLON=PLON
  ENDIF
  IF (PRESENT (PLAT)) THEN
    ALLOCATE(ZLAT(SIZE(PLON,1),SIZE(PLON,2)))
    ZLAT=PLAT
  ENDIF
ENDIF
ENDIF
!
!------------------------------------------------------------------------------
!
!*       2.    ECRITURE DU CHAMP DANS LE FICHIER DE SORTIE
!              -------------------------------------------
!
IF ( HFLAGFILE(1:3) /= 'CLO' ) THEN
  if (KVERBIA > 0) then
    print'(A,I4)','WRITELLHV: unite sortie ILUOUTLL= ', ILUOUTLL
  endif
  ! ecriture de la ligne d entete de champ
  !(temps courant)
  IAN=XDATIME(13,1)
  IMOIS=XDATIME(14,1)
  IJOUR=XDATIME(15,1)
  IHEURE=XDATIME(16,1)/3600
  IMINUTE=(XDATIME(16,1)-(IHEURE*3600))/60
  IF ( HFLAGFILE(4:5) /= 'NH') THEN
    ! first line
    write(ILUOUTLL,FMT='(I4,4(I2,X),A,A,A,A)') IAN,IMOIS,IJOUR,IHEURE,IMINUTE,TRIM(HLABELCHAMP),' ',TRIM(CUNITE(1)),&
                    ' first_line_format=Year Month Day UTCHour Minute VARIABLE_NAME UNIT'
    ! second line
    IF ( HFLAGFILE(4:5)== '1H') THEN
    ! entete unique donnant le nombre de valeurs totales ecrites lors de
    ! plusieurs appels avec OLDNH
      write(ILUOUTLL,*) 'second_line_format=values written in the same chronological order than the OBS file' 
    ELSE
    ! entete donnant exactement le nombre de valeurs ecrites lors de cet appel
      write(ILUOUTLL,FMT='(6(I4,X),A)') kkdeb,kkfin,kjdeb,kjfin,kideb,kifin ,&
                'second_line_format=values written from (k=kbeg,kend (j=jbeg,jend (i=ibeg,iend)))'
    ENDIF
  ENDIF
  !
  if (KVERBIA > 0) then
    print'(A,6(I4,X))',' kideb,kifin,kjdeb,kjfin,kkdeb,kkfin= ',kideb,kifin,kjdeb,kjfin,kkdeb,kkfin
    print'(A,2(I6,X),4(I4,X))',' ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin= ',&
    ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin
    print'(A,6(I4,X))',' dimensions de XVAR ',SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),&
                                  SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)

  endif
  ! ecriture du champ + lat,lon ,altitude du niveau
  INBVAL= (kkfin-kkdeb+1) * (kjfin-kjdeb+1) * (kifin-kideb+1)
  DO JPLOOP= KPDEB,KPFIN
    IGRID=NGRIDIA(JPLOOP)
    IF (.NOT.PRESENT(PLAT) .AND. .NOT.PRESENT(PLON)) THEN
      ZX(1:IIU,1) = XXX(1:IIU,IGRID)
      ZX(:,2:IJU) = SPREAD(ZX(:,1),2,IJU-1)
      ZY(1,1:IJU) = XXY(1:IJU,IGRID)
      ZY(2:IIU,:) = SPREAD(ZY(1,:),1,IIU-1)
      ! les 2 premiers arg. doivent etre XXHAT et XYHAT (pas XXX et XXY)
      !! peu importe en masdev4_6 car plus utilises.. 
      !CALL SM_LATLON(XXHAT,XYHAT,XLATORI,XLONORI, &
      !! supprimes en masdev4_7
      CALL SM_LATLON(XLATORI,XLONORI,             &
                     ZX,ZY,ZLAT,ZLON              )
    ENDIF
    ! init de XZZ a la grille du champ (par defaut readvar
    !l initialise a la grille 4 des vitesses verticales W)
    CALL COMPCOORD_FORDIACHRO(IGRID)
    if (KVERBIA > 0) then
      print'(A,I2)','*after COMPCOORD_FORDIACHRO ',IGRID
    endif
    DO JTRLOOP= KTRDEB,KTRFIN
      DO JTLOOP= KTDEB,KTFIN
        IAN=XDATIME(13,JTLOOP)
        IMOIS=XDATIME(14,JTLOOP)
        IJOUR=XDATIME(15,JTLOOP)
        IHEURE=XDATIME(16,JTLOOP)/3600
        IMINUTE=(XDATIME(16,JTLOOP)-(IHEURE*3600))/60
        ISEC=XDATIME(16,JTLOOP)-IHEURE*3600-IMINUTE*60
        IF ( HFLAGFILE(4:5) /= 'NH') THEN
          IF ( HFLAGFILE(4:5) == '1H') THEN       
          ! plusieurs futurs appels avec OLDNH : le nombre de lignes ne peut 
          ! etre connu a cet instant
            write(ILUOUTLL,FMT='(F10.5,X,I6,A,3(I2,X),A,2(I2,X),A,A)') XSPVAL,&
                               JTLOOP,'(',            &
                              IHEURE,IMINUTE,ISEC,')',  &
                              JTRLOOP,JPLOOP, & 
                    ' undef_value for these timenumber,',&
                    ' (UTCHour Min. Sec.), trajectorynumber, processnumber'
          ELSE
            write(ILUOUTLL,FMT='(I7,X,F10.5,X,I6,A,3(I2,X),A,2(I2,X),A,A)') INBVAL,&
                              XSPVAL,JTLOOP,'(',            &
                              IHEURE,IMINUTE,ISEC,')',  &
                              JTRLOOP,JPLOOP, & 
                    'number_of_next_lines, undef_value for these timenumber,',&
                    ' (UTCHour Min. Sec.), trajectorynumber, processnumber'
          ENDIF
        ENDIF
        DO JKLOOP= kkdeb,kkfin
          SELECT CASE ( HTYPEOUT(1:4) )
          CASE ('LLHV','llhv') 
            IF (kkdeb == 1 .AND. kkfin == 1) THEN
              ! champ 2D: altitude donnee par PALT(:,:,IGRID) ou XXZS(:,:,IGRID)
              DO JJLOOP= kjdeb,kjfin
              DO JILOOP= kideb,kifin
                IF (PRESENT (PALT) ) THEN
                  if (KVERBIA > 0) then
                    print'(A,I2,X,F10.5)', 'LLHV 2D igrid PALT(:,:)= ',IGRID, &
                                                       PALT(JILOOP,JJLOOP,IGRID)
                  endif
                  IF (HTYPEOUT(1:4)=='LLHV') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLON(JILOOP,JJLOOP), &
                                            ZLAT(JILOOP,JJLOOP), &
                                            PALT(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='llhv') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLAT(JILOOP,JJLOOP), &
                                            ZLON(JILOOP,JJLOOP), &
                                            PALT(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ELSE
                  IF (HTYPEOUT(1:4)=='LLHV') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLON(JILOOP,JJLOOP), &
                                            ZLAT(JILOOP,JJLOOP), &
                                            XXZS(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='llhv') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLAT(JILOOP,JJLOOP), &
                                            ZLON(JILOOP,JJLOOP), &
                                            XXZS(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ENDIF
              END DO
              END DO
            ELSE
              ! champ 3D
              !altitude des niveaux donnee par XZZ ou PALT
              DO JJLOOP= kjdeb,kjfin
              DO JILOOP= kideb,kifin
                IF (PRESENT (PALT) ) THEN
                  if (KVERBIA > 0 .AND. JILOOP==1 .AND. JJLOOP==1) then
                    print '(A,I4,X,F10.5)', 'LLHV 3D K,PALT(1,1,K)= ',JKLOOP, &
                                                      PALT(JILOOP,JJLOOP,JKLOOP)
                  endif
                  IF (HTYPEOUT(1:4)=='LLHV') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLON(JILOOP,JJLOOP),       &
                                            ZLAT(JILOOP,JJLOOP),       &
                                            PALT(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='llhv') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLAT(JILOOP,JJLOOP),       &
                                            ZLON(JILOOP,JJLOOP),       &
                                            PALT(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ELSE
                  IF (HTYPEOUT(1:4)=='LLHV') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLON(JILOOP,JJLOOP),       &
                                            ZLAT(JILOOP,JJLOOP),       &
                                            XZZ(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='llhv') THEN
                    WRITE(ILUOUTLL,FMT=1000)ZLAT(JILOOP,JJLOOP),       &
                                            ZLON(JILOOP,JJLOOP),       &
                                            XZZ(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ENDIF
              END DO
              END DO
            ENDIF
          CASE ('IJHV','jihv') 
            IF (kkdeb == 1 .AND. kkfin == 1) THEN
              ! champ 2D: altitude donnee par PALT(:,:,IGRID) ou XXZS(:,:,IGRID)
              DO JJLOOP= kjdeb,kjfin
              DO JILOOP= kideb,kifin
                IF (PRESENT (PALT) ) THEN
                  if (KVERBIA > 0) then
                    print '(A,I2,X,F10.5)', 'IJHV 2D igrid PALT(:,:)= ',IGRID, &
                                                       PALT(JILOOP,JJLOOP,IGRID)
                  endif
                  IF (HTYPEOUT(1:4)=='IJHV') THEN
                    WRITE(ILUOUTLL,FMT=1001) JILOOP, &
                                             JJLOOP, &
                                            PALT(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='jihv') THEN
                    WRITE(ILUOUTLL,FMT=1001)JJLOOP, &
                                            JILOOP, &
                                            PALT(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ELSE
                  IF (HTYPEOUT(1:4)=='IJHV') THEN
                    WRITE(ILUOUTLL,FMT=1001)JILOOP, &
                                            JJLOOP, &
                                            XXZS(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='jihv') THEN
                    WRITE(ILUOUTLL,FMT=1001)JJLOOP, &
                                            JILOOP, &
                                            XXZS(JILOOP,JJLOOP,IGRID), & 
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ENDIF
              END DO
              END DO
            ELSE
              ! champ 3D
              !altitude des niveaux donnee par XZZ ou PALT
              DO JJLOOP= kjdeb,kjfin
              DO JILOOP= kideb,kifin
                IF (PRESENT (PALT) ) THEN
                  if (KVERBIA > 0 .AND. JILOOP==1 .AND. JJLOOP==1) then
                    print '(A,I4,X,F10.5)', 'IJHV 3D K,PALT(1,1,K)= ',JKLOOP, &
                                                      PALT(JILOOP,JJLOOP,JKLOOP)
                  endif
                  IF (HTYPEOUT(1:4)=='IJHV') THEN
                    WRITE(ILUOUTLL,FMT=1001)JILOOP,       &
                                            JJLOOP,       &
                                            PALT(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='jihv') THEN
                    WRITE(ILUOUTLL,FMT=1001)JILOOP,       &
                                            JJLOOP,       &
                                            PALT(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ELSE
                  IF (HTYPEOUT(1:4)=='IJHV') THEN
                    WRITE(ILUOUTLL,FMT=1001)JILOOP,       &
                                            JJLOOP,       &
                                            XZZ(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ELSE IF (HTYPEOUT(1:4)=='jihv') THEN
                    WRITE(ILUOUTLL,FMT=1001)JILOOP,       &
                                            JJLOOP,       &
                                            XZZ(JILOOP,JJLOOP,JKLOOP), &
                                            XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
                  ENDIF
                ENDIF
              END DO
              END DO
            ENDIF
  
          CASE ('LLZV','llzv','LLPV','llpv','LLAV','llav')
            IF (PRESENT (PALT) ) THEN
            !altitude des niveaux donnee par PALT
              if (KVERBIA > 0) then
                print'(A,A,I4,X,F10.5)', HTYPEOUT(1:4),' K,PALT(1,1,K)= ',JKLOOP,PALT(1,1,JKLOOP)
              endif
            ELSE
              PRINT*,'** WRITELLHV: les altitudes doivent etre passees par argument'
              PRINT*,'          pour HTYPEOUT= ',HTYPEOUT(1:4)
              KRETCODE=1
              RETURN
            ENDIF
            DO JJLOOP= kjdeb,kjfin
            DO JILOOP= kideb,kifin
              IF (HTYPEOUT(1:2)=='LL') THEN
                WRITE(ILUOUTLL,FMT=1000)ZLON(JILOOP,JJLOOP),       &
                                        ZLAT(JILOOP,JJLOOP),       &
                                        PALT(1,1,JKLOOP),          &
                                        XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
              ELSE IF (HTYPEOUT(1:2)=='ll') THEN
                WRITE(ILUOUTLL,FMT=1000)ZLAT(JILOOP,JJLOOP),       &
                                        ZLON(JILOOP,JJLOOP),       &
                                        PALT(1,1,JKLOOP),          &
                                        XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
              ENDIF
            END DO
            END DO
          CASE ('IJZV','jizv','IJPV','jipv') 
            IF (PRESENT (PALT) ) THEN
            !altitude des niveaux donnee par PALT
              if (KVERBIA > 0) then
                print'(A,A,I4,X,F10.5)', HTYPEOUT(1:4),' K,PALT(1,1,K)= ',JKLOOP,PALT(1,1,JKLOOP)
              endif
            ELSE
              PRINT*,'** WRITELLHV: les altitudes doivent etre passees par argument'
              PRINT*,'          pour HTYPEOUT= ',HTYPEOUT(1:4)
              KRETCODE=1
              RETURN
            ENDIF
            DO JJLOOP= kjdeb,kjfin
            DO JILOOP= kideb,kifin
              IF (HTYPEOUT(1:2)=='IJ') THEN
                WRITE(ILUOUTLL,FMT=1001)JILOOP,       &
                                        JJLOOP,       &
                                        PALT(1,1,JKLOOP),          &
                                        XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
              ELSE IF (HTYPEOUT(1:2)=='ji') THEN
                WRITE(ILUOUTLL,FMT=1001)JJLOOP,       &
                                        JILOOP,       &
                                        PALT(1,1,JKLOOP),          &
                                        XVAR(JILOOP,JJLOOP,JKLOOP,JTLOOP,JTRLOOP,JPLOOP)
              ENDIF
            END DO
            END DO
  
          END SELECT
        END DO
      END DO
    END DO
  END DO
!    
1000  FORMAT ( 2(F11.6,1x),F8.2,1x,E15.9)
1001  FORMAT ( 2(I4,1x),F8.2,1x,E15.9)

  if (KVERBIA >= 0) then
    print*,'WRITELLHV: ecriture de ',TRIM(HLABELCHAMP)
    print*,'--------- '
  endif
ENDIF
!-------------------------------------------------------------------------------
!
!*       3.    FERMETURE DU FICHIER DE SORTIE
!              ------------------------------
!
IF ( HFLAGFILE(1:3) == 'CLO' ) THEN
  if (KVERBIA > 0) then
    print*,'WRITELLHV: before closing file ',TRIM(YFILEOUT),' unit ',iluoutll
  endif
  !
  ! fichier de sortie
  CLOSE(UNIT=ILUOUTLL)
  CALL FMFREE(YFILEOUT,CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
  IF( NRESPDIAS(NBFILES)==0 ) THEN
    if (KVERBIA >= 0) then
            print*,'End of WRITELLHV: File ',TRIM(YFILEOUT),' available with format ',HTYPEOUT 
      print*,'--------- '
    endif
  ELSE
    PRINT*,' ****WRITELLHV: error when closing ', TRIM(YFILEOUT), &
           ' code= ',NRESPDIAS(NBFILES)
    KRETCODE=2
    RETURN
  ENDIF
  !
ENDIF
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITELLHV
