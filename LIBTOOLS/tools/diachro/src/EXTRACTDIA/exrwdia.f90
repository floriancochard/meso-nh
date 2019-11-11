      PROGRAM  EXRWDIA
!     ###################
!
!!****  *EXRWDIA* -  lecture d'enregistrements dans des fichiers diachroniques,
!                    traitement,
!                    ecriture (plusieurs types de format de fichier possibles)
!     DIAC= fichier diachronique utilisable via diaprog (appel à WRITEVAR)
!     LLHV= fichier ascii lon,lat,altitude,valeur (appel à WRITELLHV)
!        les sorties LLZV LLPV llzv llpv sont codées dans extractdia
!        conseil: sortir en format 'DIAC' puis utiliser extractdia
!                 pour des sorties LLZV LLPV llzv llpv
!     KCDL= fichier netcdf (appel à WRITECDL)
!        les sorties ZCDL ou PCDL sont codées dans extractdia
!        conseil: sortir en format 'DIAC' puis utiliser extractdia
!                 pour des sorties ZCDL ou PCDL
!     FREE= fichier ascii , l ecriture est à coder par l utilisateur
! 
! CONTRAINTES:
! Au maximum 50 fichiers ouverts simultanement
! 50 = limite du module MODD_FILES_DIACHRO
! Au max 44 fichiers simultanement ouverts par FMOPEN (c.a.d via
! READVAR et WRITEVAR )
!! 
!!
!!    PURPOSE
!!    -------
! 
!!**  METHOD
!!    ------
!      Exemple de programme simple a adapter aux besoins
!      pour info supplementaires, voir le programme interactif extractdia.f90
!
!  Rappel1: un fichier LFI diachronique (000 ou issu de conv2dia) contient
!          des champs stockes dans un tableau a 6 dimensions (XVAR passe
!          par module a toutes les routines de traitement)
!  le logiciel graphique diaprog interprete les dimensions ainsi:
!          XVAR( dimension1=i=grille horizontale selon x,
!                dimension2=j=grille horizontale selon y,
!                dimension3=k=grille verticale selon z,
!                dimension4=t=echeances temporelles,
!                dimension5=traj=masques des budgets ou trajectoires,
!                dimension6=p/proc=processus )
!
!  Rappel2: les variables sont stockees sur 7 grilles differentes dans
!          les fichiers LFI ( 1=grille masse , 3=grille W...)
!  Voir le book3 de Mesonh pour traiter correctement ces differentes
! localisations que peuvent representer 1/2deltax sur l horizontale et
! 1/2 niveau selon la verticale
!  XVAR( i,j,k,:,:,:) pour U n est pas localise au meme lieu que
!  XVAR( i,j,k,:,:,:) pour V et XVAR( i,j,k,:,:,:) pour Theta
!
! Rappel3: les composantes U et V sont dans le plan de projection Mesonh 
!         (cartesien ou conforme) et ne correspondent pas a Uzonal et Vmeridien.
! Utiliser les routines DD, FF et UV_TO_ZONAL_AND_MERID pour changer de repere.
!
!!
! READVAR : lit un champ et alimente un tableau XVAR + grille
!          et tous les parametres necessaires aux traitements futurs
!           transforme certaines unites pour des traitements plus corrects:
!            les dBz sont passees en Ze , les temp. de brillance en W
! WRITEVAR, WRITECDL, WRITELLHV pour changer de format
!           les routines writevar,writecdl,writellh effectuent la 
!          transformation inverse  sur les unites avant ecriture      
! Voir les routines TO_COMPUTING_UNITS et FROM_COMPUTING_UNITS pour le
! detail des variables traitees.
!!      
!!
!!    REFERENCE
!!    ---------
!!    'CREATION et EXPLOITATION de FICHIERS DIACHRONIQUES' J.Duron oct.2001
!!
!!    AUTHORS
!!    -------
!!    I. Mallet et N. Asencio
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/03/2003
!!      Modifications 01/2005 : Nicole Asencio
!!        ajout de modules et des commentaires pour une utilisation sur 
!!        des fichiers diachroniques 000
!!      17/06/2005 : ajout de commentaires sur l utilisation de XZZ
!!                   et de la routine MOYZ
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! modules MesoNH
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT,XUNDEF
!                    NIMAX,NJMAX,NKMAX, NIINF, NISUP
USE MODD_DIM1
!                    grille : XXDXHAT(:,1:7) et XXX(:,1:7), XXZS(:,:,1:7)
USE MODD_COORD
!                    ref grille: XLON0,XLAT0,XBETA,XRPK
USE MODD_GRID
!                    descriptif grille: XXHAT(:) ,XLAT(:,:),XDXHAT(:),XMAP(:,:)
!                    ,XZS(:,:),XZZ(:,:,:) ,XCOSSLOPE(:,:),XDIRCOSXW(:,:)
USE MODD_GRID1
!
! modules DIACHRO
USE MODN_NCAR,  ONLY: XSPVAL      
USE MODD_ALLOC_FORDIACHRO, ONLY: XVAR, &         ! XVAR(i,j,k,t,n,p)
                                 XDATIME, &      ! XDATIME(16,t)
                                 CCOMMENT,&      ! CCOMMENT(p)
                                 CTITRE, CUNITE,&! CTITRE(p),CUNITE(p)
                                 NGRIDIA,&       ! NGRIDIA(p)
                                 XTRAJT          ! XTRAJT(t,n)
USE MODD_RESOLVCAR, ONLY: CGROUP, NVERBIA, &
                          NNDIA, NPROCDIA, NBPROCDIA !pour appel a interp_grids                                 
USE MODD_COORD, ONLY: XXX,XXY,XXZS, & !  XXX(:,1:7), XXY(:,1:7), XXZS(:,:,1:7)
                      XXDXHAT,XXDYHAT ! XXDXHAT(:,1:7), XXDYHAT(:,1:7)
USE MODD_PVT, ONLY: LPRESYT                                 
USE MODD_TYPE_AND_LH, ONLY: CTYPE,LICP,LJCP,LKCP
! 
! modules tools
USE MODI_CHANGE_A_GRID          ! changement de grille dans les grilles mesonh
USE MODI_ZINTER                 ! interpolation a Z=cst
USE MODI_PINTER                 ! interpolation a P=cst
USE MODI_ZMOY                   ! moyenne sur une couche verticale
USE MODI_DD                     ! calcul dd ,ff a partir de U,V grille mesonh
USE MODI_FF
USE MODI_WRITELLHV              ! routines
USE MODI_WRITECDL               !d
USE MODI_WRITEVAR               !ecriture
USE MODI_FROM_COMPUTING_UNITS   ! voir routine symetrique TO_COMPUTING_UNITS
                                !pour la liste des variables traitees
USE MODI_HOR_INTERP_4PTS        ! interpolation horizontale 4 points
USE MODI_UV_TO_ZONAL_AND_MERID  ! passage composantes vent Mesonh
                                !a Zonal+ Meridien
USE MODI_LOW2UP                 ! conversion en Majuscules
!
! modules extractdia
USE MODD_READLH                 ! domaine initialise par READVAR: 
                                !NREADIL,NREADIH, NREADJL,NREADJH,
                                !NREADKL,NREADKH
!
IMPLICIT NONE
!
!*       0.1   Local variables
!
!
INTEGER           :: JI,JJ,JK,J4,J5,J6,ILECTTRAITE,NBLECTTRAITE
INTEGER           :: ilocverbia,iret,inbvertz,ikdebzint,IGRID,ISKIP
! zoom recalculé en fonction des dimensions du champ traite
INTEGER           :: ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin
INTEGER           :: ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup      
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZWORK3D,ZWORK3D2,ZVARZCST
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZWORK2D
REAL, allocatable, dimension(:,:,:,:,:,:):: ZVARSAVE 
REAL :: ZMIN,ZMAX
! **** la taille des variables caracteres contenant les noms
!      de fichiers est obligatoirement de 28 ****
CHARACTER(LEN=28) :: YFILEIN
! **** la longueur du nom ne doit pas depasser 13 car. si le fichier
! contient des groupes a un seul PROCessus, ou 9 si plusieurs PROCessus ****
CHARACTER(LEN=13) :: YCHAMP
!
CHARACTER(LEN=4)  :: YTYPEOUT
CHARACTER(LEN=3)  :: YFLAGREADVAR ,YFLAGWRITE
CHARACTER(LEN=3)  :: YSUFFIX=''
REAL , allocatable, dimension(:) :: zlistevert      
!
CHARACTER(LEN=13) :: YCHAMP2
CHARACTER(LEN=2), DIMENSION(15) :: LIST
!-------------------------------------------------------------------------------
!
!*       1.    INIT
!              ----
!
! active(1) ou desactive(0) les prints de controle dans les routines
! READVAR et WRITEVAR
ilocverbia=1 
! active(1) ou desactive(0) les prints de controle dans les routines diachro
NVERBIA=0  
! 
XSPVAL=XUNDEF      ! dans mesonh Xundef est utilise 
                   ! dans les routines diachro XSPVAL est utilise
!
!
!*       1.2   Init de parametres pour la lecture
!                                       
! nom du fichier diachronique en supprimant .lfi
YFILEIN='fichier diachronique en supprimant .lfi' 
! indique que le fichier lu doit etre ouvert dans READVAR
!(initialisation des variables des modules documentés en debut de programme)
! rq: si d autres fichiers traites dans ce programme, remettre 'OPE'
!avant le 1er appel a READVAR pour chaque fichier
YFLAGREADVAR='OPE'        
! type du format de sortie (DIAC/LLHV/FREE/KCDL)
YTYPEOUT='DIAC'    
! ouverture du fichier et ecriture de l entete dans les routines WRITExxx
YFLAGWRITE='NEW' 
! nom du champ a lire
YCHAMP='THM' 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!Boucle lecture + traitement----------------------------------------------------
!NBLECTTRAITE=3
!DO ILECTTRAITE=1,NBLECTTRAITE
!
!-------------------------------------------------------------------------------
!
!*       3.     Lecture du champ YCHAMP et stockage dans XVAR
!              ----------------------
!                                       
CALL READVAR(YCHAMP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
!
!!  apres cette lecture, les tableaux suivants sont disponibles:
! a. si YFLAGREADVAR='OPE' lecture de l entete du fichier:
!  X,Y,Z-HAT(m) dans XXX,XXY,XXZ(:,1:7)        ! (MODD_COORD)
!  topography altitude values(m):XXZS(:,:,1:7) ! (MODD_COORD)
!  meshsize values XXDXHAT,XXDYHAT(:,1:7)      ! (MODD_COORD)
! rq: derniere dimension (1-7) fait reference aux 7 grilles de MesoNH
!    (1: masse, 2: flux selon x, 3: flux selon y, 4: flux suivant z ,
!     5:Vertical vorticity , 6:y-component vorticity, 7:x-component vorticity )
!  NIMAX,NJMAX,NKMAX,  LCARTESIAN, LTHINSHELL,CSTORAGE_TYPE,
!  XXHAT(IIU),XYHAT(IJU),XZHAT(IKU))
!  XMAP(IIU,IJU)  XLAT(IIU,IJU) XLON(IIU,IJU)
!  XDXHAT(IIU),XDYHAT(IJU)
!  XZS(IIU,IJU) ,  XZZ(IIU,IJU,IKU) avec XZZ=grille pour numéro de
!                  grille=4 quelle que soit la grille de variable lue
!                  CALL COMPCOORD_FORDIACHRO(IGRID_var) pour initialiser XZZ
!                  avec la grille de la variable traitée
!  TDTMOD,TDTCUR,TDTEXP,TDTSEG,
!  NSTOP,NOUT_TIMES,NOUT_NUMB, XTSTEP,XSEGLEN,                       
! b. variables relatives a l enregistrement demande:
!  XVAR(i,j,k,t,n,p)= champ
!  NGRIDIA(p)= indice de grille des p processus 
!  CTYPE= CART/MASK/SPXY/SSOL/DRST/RSPL/RAPL
!  CTITRE(p)= titre des p processus UTILISE DANS XMIN_nom de diaprog
!  CUNITE(p)= unite des p processus
!  CCOMMENT(p)= commentaire des p processus
!  XDATIME(16,t)= dates relatives au champ
!     XDATIME(1,t)=TDTEXP%TDATE%YEAR; XDATIME(2,t)=TDTEXP%TDATE%MONTH
!     XDATIME(3,t)=TDTEXP%TDATE%DAY;  XDATIME(4,t)=TDTEXP%TIME
!     XDATIME(5,t)=TDTSEG%TDATE%YEAR; XDATIME(6,t)=TDTSEG%TDATE%MONTH
!     XDATIME(7,t)=TDTSEG%TDATE%DAY;  XDATIME(8,t)=TDTSEG%TIME
!     XDATIME(9,t)=TDTMOD%TDATE%YEAR; XDATIME(10,t)=TDTMOD%TDATE%MONTH
!     XDATIME(11,t)=TDTMOD%TDATE%DAY; XDATIME(12,t)=TDTMOD%TIME
!     XDATIME(13,t)=TDTCUR%TDATE%YEAR;XDATIME(14,t)=TDTCUR%TDATE%MONTH
!     XDATIME(15,t)=TDTCUR%TDATE%DAY; XDATIME(16,t)=TDTCUR%TIME      
!     XTRAJT(t,n)= nombre de secondes depuis le debut de la simulation
! optionnels suivant la valeur de CTYPE:
!XTRAJX-Y-Z(k,t,n) XMASK(i,j,t,n)
! rq: p=1 (nb de processus) si fichier pseudo-diachro sorti de conv2dia
! rq: pour plus d infos sur la nature d un enregistrement dans un
! fichier diachronique, voir 'CREATION et EXPLOITATION de FICHIERS
! DIACHRONIQUES' (J. Duron, octobre 2001)
!      
!!
! la routine READVAR a modifie YFLAGREADVAR a 'NOP'
!pour indiquer que le fichier courant est deja ouvert 
!(le prochain champ sera lu sans initialisation des modules relatifs a l entete)
!
! la routine READVAR a transforme certaines unites pour des traitements
!plus corrects:
! les dBz sont passees en Ze , les temp. de brillance en W
!
! les routines writevar,writecdl,writellh effectuent la transformation
!inverse avant ecriture
!
!   Definir le zoom a traiter dans les calculs: 
! valeurs par defaut:
    ivarideb=NREADIL
    ivarifin=NREADIH
    ivarjdeb=NREADJL
    ivarjfin=NREADJH
    ivarkdeb=NREADKL
    ivarkfin=NREADKH
    ivartinf=1
    ivartsup=size(XVAR,4)
    ivartrajinf=1
    ivartrajsup=size(XVAR,5)
    ivarprocinf=1
    ivarprocsup=size(XVAR,6)                   
!
!-------------------------------------------------------------------------------
!
!*       4.    EXEMPLES DE CALCUL  activer les lignes de code en
!            !optionx en début de ligne 
!              ------------------
!
!*      4.1    Interpolation sur la grille de masse Mesonh
!                                       
! replace field at mass points
!option1      ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
!option1      ALLOCATE(ZWORK3D2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
!option1      DO J6=ivarprocinf,ivarprocsup
!option1        IGRID=NGRIDIA(J6)
!option1        IF(SIZE(XVAR,3)/=1 .OR. IGRID/=4) THEN 
!option1          ! pas d interpolation verticale pour champ 2D
!option1          DO J5=ivartrajinf,ivartrajsup
!option1            DO J4=ivartinf,ivartsup
!option1              ZWORK3D(:,:,:)=XVAR(:,:,:,J4,J5,J6)
!option1              CALL CHANGE_A_GRID(ZWORK3D,IGRID,ZWORK3D2)
!option1              ! IGRID=1 en sortie de change_a_grid
!option1              XVAR(:,:,:,J4,J5,J6)=ZWORK3D2(:,:,:)
!option1              NGRIDIA(J6)=1
!option1            ENDDO
!option1          ENDDO
!option1        ENDIF
!option1      ENDDO
!option1      DEALLOCATE(ZWORK3D,ZWORK3D2)
!
!*      4.2   Maximum du champ sur la verticale
!
!option2      ALLOCATE(ZWORK2D(SIZE(XVAR,1),SIZE(XVAR,2)))
!option2      ZWORK2D=0.
!option2      ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
!option2                        size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
!option2      ZVARSAVE=XVAR
!option2      DEALLOCATE(XVAR)
!option2      ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),1,               &
!option2                    size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)))
!option2      DO J6=ivarprocinf,ivarprocsup
!option2      DO J5=ivartrajinf,ivartrajsup
!option2      DO J4=ivartinf,ivartsup
!option2        DO JK=ivarkdeb,ivarkfin
!option2        DO JJ=ivarjdeb,ivarjfin
!option2        DO JI=ivarideb,ivarifin
!option2          IF (ZVARSAVE(JI,JJ,JK,J4,J5,J6) .GT.  ZWORK2D(JI,JJ) ) THEN
!option2            ZWORK2D(JI,JJ)= ZVARSAVE(JI,JJ,JK,J4,J5,J6) 
!option2          ENDIF
!option2        END DO
!option2        END DO
!option2        END DO
!option2        XVAR(:,:,1,J4,J5,J6)=ZWORK2D(:,:)
!option2     END DO
!option2     END DO
!option2     CTITRE(J6)='MAX'//ADJUSTL(TRIM(CTITRE(J6)))
!option2     END DO
!option2     DEALLOCATE(ZVARSAVE,ZWORK2D)
!option2     YCHAMP='MAX'//ADJUSTL(TRIM(YCHAMP))
!option2     CCOMMENT(ivarprocinf:ivarprocsup)=nouveau_comment
!option2     CUNITE(ivarprocinf:ivarprocsup)=nouvelle_unite
!
!*      4.3    Interpolation sur des niveaux Z=cst ou P=cst
!
!option3!     inbvertz=nombre de niveaux verticaux souhaite
!option3!     allocate ( zlistevert(inbvertz))
!option3!     zlistevert= tableau contenant les differentes valeurs de Z en metres
!option3!                                                              P en hPa
!option3!     print * , ' interpolations sur ',inbvertz,' niveaux'
!option3!     print *, 'niveaux=',zlistevert 
!option3     YSUFFIX='zcl'
!option3       IF (YSUFFIX  == 'zcl'  .AND. SIZE(XVAR,3) > 1 .AND. &
!option3           SIZE(XVAR,2) > 1 .AND. SIZE(XVAR,1) > 1                ) THEN
!option3         ! ALT ne passe pas cette partie car ses dimensions 1 et 2 =1
!option3         if (ilocverbia >= 0 ) then
!option3           print*,' interpolations sur Z=cst',inbvertz,' niveaux'
!option3         endif
!option3         if (ilocverbia >= 1 ) then
!option3           print*,'niveaux= ',zlistevert 
!option3         endif
!option3         ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
!option3                           size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
!option3         ZVARSAVE=XVAR
!option3         ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
!option3         ALLOCATE(ZVARZCST(SIZE(XVAR,1),SIZE(XVAR,2),inbvertz))
!option3         DEALLOCATE(XVAR)
!option3         ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),SIZE(ZVARZCST,3),&
!option3                       size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)))
!option3         DO J6=ivarprocinf,ivarprocsup
!option3           IGRID=NGRIDIA(J6)
!option3           ! init du tableau des altitudes  XZZ pour la grille= IGRID
!option3           CALL COMPCOORD_FORDIACHRO(IGRID)
!option3           DO J5=ivartrajinf,ivartrajsup
!option3             DO J4=ivartinf,ivartsup
!option3               ZWORK3D(:,:,:)=ZVARSAVE(:,:,:,J4,J5,J6)
!option3               ikdebzint=2
!option3               CALL ZINTER(ZWORK3D,XZZ,ZVARZCST,zlistevert,ikdebzint,XSPVAL)
!option3               ou bien       a P=cst precede d un READVAR de ZPABS
!option3               CALL PINTER(ZWORK3D,IGRID,XSPVAL,zlistevert,ZVARZCST,ZPABS)
!option3               XVAR(:,:,:,J4,J5,J6)=ZVARZCST
!option3             END DO
!option3           END DO
!option3         END DO
!option3         DEALLOCATE(ZVARSAVE,ZVARZCST,ZWORK3D)
!option3         ivarkdeb=1
!option3         ivarkfin=inbvertz
!option3         IF (ilocverbia >= 5 ) then
!option3           print*,'ivarkdeb,ivarkfin= ',ivarkdeb,ivarkfin 
!option3         ENDIF
!option3       ENDIF
!
!*      4.4    Moyenne verticale entre deux niveaux zmin et zmax
!              pour des variables lues sans prise en compte du volume
!              de chaque maille (*RHO les variables avant l appel si
!              nécessaire)
!
!option4!     zmin=base  
!option4!     zmax=sommet
!option4      ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
!option4                        size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
!option4      ZVARSAVE=XVAR
!option4      ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
!option4      ALLOCATE(ZWORK2D(SIZE(XVAR,1),SIZE(XVAR,2)))
!option4      DEALLOCATE(XVAR)
!option4      ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),1,&
!option4                    size(ZVARSAVE,4),size(ZVARSAVE,5),size(ZVARSAVE,6)))
!option4      DO J6=ivarprocinf,ivarprocsup
!option4        IGRID=NGRIDIA(J6)
!option4        DO J5=ivartrajinf,ivartrajsup
!option4          DO J4=ivartinf,ivartsup
!option4            ZWORK3D(:,:,:)=ZVARSAVE(:,:,:,J4,J5,J6)
!option4            ! JPVEXT, JPHEXT: points a exclure verticalement et horizontalement
!option4            CALL ZMOY(ZWORK3D,IGRID,ZWORK2D,zmin,zmax,XSPVAL,JPVEXT,JPHEXT)
!option4            XVAR(:,:,1,J4,J5,J6)=ZWORK2D(:,:)
!option4          END DO
!option4        END DO
!option4        CTITRE(J6)='MEANZ'//ADJUSTL(TRIM(CTITRE(J6)))
!option4      END DO
!option4      DEALLOCATE(ZVARSAVE,ZWORK2D,ZWORK3D)
!option4      YCHAMP='MEANZ'//ADJUSTL(TRIM(YCHAMP))
!option4      !CCOMMENT(ivarprocinf:ivarprocsup)=nouveau_comment
!option4      !CUNITE(ivarprocinf:ivarprocsup)=nouvelle_unite
!
!        *      4.5    Calcul de la direction ou de la force du vent
!
!option5! lecture de la 1e composante du vent (UT ou UM ou LSUM ou UMxx ou UTxx)
!option5!stockee dans ZVARSAVE
!option5      ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
!option5                        size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
!option5      ZVARSAVE=XVAR
!option5      ! lecture de la 2e composante du vent (VT ou VM ou LSVM ou VMxx ou VTxx)
!option5      !stockee dans XVAR
!option5      !      YCHAMP='2e composante du vent'
!option5      !CALL READVAR(YCHAMP,YFILEIN,YFLAGREADVAR,ilocverbia,iret)
!option5      ! 
!option5      ALLOCATE(ZWORK3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
!option5      ZWORK3D=XSPVAL
!option5      DO J6=ivarprocinf,ivarprocsup
!option5        IGRID=NGRIDIA(J6)
!option5        DO J5=ivartrajinf,ivartrajsup
!option5          DO J4=ivartinf,ivartsup
!option5            iskip=1 ! tous les points de grille
!option5            CALL DD(ZVARSAVE(:,:,:,J4,J5,J6),XVAR(:,:,:,J4,J5,J6),ZWORK3D, &
!option5                    iskip,IGRID)
!option5            !CALL FF(ZVARSAVE(:,:,:,J4,J5,J6),XVAR(:,:,:,J4,J5,J6),ZWORK3D, &
!option5            !        JPVEXT,JPHEXT,IGRID)
!option5            XVAR(:,:,:,J4,J5,J6)=ZWORK3D(:,:,:)
!option5          END DO
!option5        END DO
!option5        NGRIDIA(J6)=1    ! resultat sur la grille de masse
!option5        CUNITE(J6)='degrees' ! pour dd
!option5      END DO
!option5      DEALLOCATE(ZVARSAVE,ZWORK3D)
!    
!
!        *      4.6    Compression de budgets 3D (equivalent au type
!                      CART de MAINPROG=Model)
!
!option6 ! Sauvagerde la variable à traiter precedemment lue
!option6      ALLOCATE(ZVARSAVE(size(XVAR,1),size(XVAR,2),size(XVAR,3),   &
!option6                        size(XVAR,4),size(XVAR,5),size(XVAR,6))   )
!option6      ZVARSAVE=XVAR      
!option6 !  récupère le champ représentant la densité
!option6 !       RJS pour scalaires, RJX pour U, RJY pour V RJZ pour W 
!option6 CALL READVAR("RJS_0001", YFILEIN, YFLAGREADVAR,ilocverbia,iret)
!option6 ALLOCATE(RHODJS(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3)))
!option6 RHODJS(:,:,:)=XVAR(:,:,:,1,1,1)      
!option6 ! donc pour la partie calcul,la variable devient : 
!option6 WHERE(ZVARSAVE(:,:,:) /= XUNDEF) ZVARSAVE(:,:,:)= ZVARSAVE(:,:,:)* RHODJS(:,:,:)
!option6 
!option6 ! exemple pour compression suivant Z 
!option6 ALLOCATE(ZWORK2D(SIZE(XVAR,1),SIZE(XVAR,2)))
!option6 ALLOCATE(ZWORK2DRho(SIZE(XVAR,1),SIZE(XVAR,2)))
!option6 ZWORK2D=0.      
!option6       
!option6 do J6=ivarprocinf, ivarprocsup
!option6 do J5=ivartrajinf, ivartrajsup
!option6 do J4=ivartinf, ivartsup
!option6       ZWORK2D(:,:) = 0.0 
!option6       ZWORK2DRHo(:,:) = 0.0    
!option6       do JK=ivarkdeb, ivarkfin
!option6             ZWORK2DRHo(:,:)=ZWORK2DRHo(:,:)+RHODJS(:,:,JK)
!option6             ZWORK2D(:,:)=ZWORK2D(:,:) + ZVARSAVE(:,:,JK,J4,J5,J6))   
!option6       enddo  !Loop on K
!option6       ! stockage dans XVAR qui sera utilisé pour l ecriture
!option6       ! pour le tracé ecriture de la variable/densité
!option6       ! remq: pour 2 boites 1 et 2 :var*rho1+var*rho2/(rho1+rho2)
!option6       WHERE(abs(ZWORK2DRho(:,:)).gt.0.0)    
!option6         XVAR(:,:,1,J4,J5,J6)=ZWORK2D(:,:)/ZWORK2DRHo(:,:)
!option6       ELSEWHERE
!option6         XVAR(:,:,1,J4,J5,J6)=XUNDEF
!option6       ENDWHERE        
!option6 ENDDO
!option6 ENDDO
!option6 ENDDO
!option6       ! nouvelles limites pour l ecriture      
!option6         ! si ce champ 2D = moyenne ou compression (trace sans relief)
!option6         ! LJCP=.TRUE.  si  compresse selon j 
!option6         ! LICP=.TRUE.  si  compresse selon i
!option6         ! LKCP=.TRUE.  si  compresse selon k
!option6   ivarkdeb=1 
!option6   ivarkfin=1
!option6   LKCP=.TRUE.
!
!-------------------------------------------------------------------------------
!
!*      5.      traitement perso
!              ----------------------
!              ----------------------
!  Si vous effectuez des calculs sur des variables "Rapport de melange"
!  ne pas oublier de *RHODREF si necessaire (hydrometeores, flux, ...)
!      
! Préférer l utilisation d'un tableau de travail qui conserve les
! dimensions de Xvar pour les dimensions 1,2 et 3 et définir ensuite
! le zoom pour écrire un sous tableau via les variables ivar.deb et
! ivar.fin
!
!      
! ..... code utilisateur .....
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!  
!*      6.     Preparation de l ecriture
!              ----------------------                                   
!
!   XVAR contient le champ a ecrire (module MODD_ALLOC_FORDIACHRO):
! vous pouvez modifier les variables suivantes si vous ne conservez pas les 
! valeurs lues par READVAR
!     YCHAMP=nouveau_nom EN MAJUSCULES
! **** la longueur de YCHAMP ne doit pas depasser 13 car. si le fichier
! contient des groupes à un seul PROCessus, ou 9 si plusieurs PROCessus ****
!     NGRIDIA(p)=nouvelle_grille
!     CTITRE(p)=nouveau_nom_p  utilise dans les directives DIAPROG XISOLEV_
!     CCOMMENT(p)=nouveau_comment
!     CUNITE(p)=nouvelle_unite
!     XDATIME(1:16,t)= dates du fichier
!     XTRAJT(t,n)= nombre de secondes depuis le début de la simulation

! VOUS DEVEZ MODIFIER LA TAILLE DE CES TABLEAUX si vous modifiez la
!  dimension 4 ou 5 ou 6 du tableau XVAR a ecrire 
!
! fichier a creer par  les routines d ecriture WRITEVAR, WRITELLHV, WRITECDL 
!(au premier appel de la routine)
!   nom du fichier de sortie= YFILEIN+suffixe (suffixe=2 par défaut) si writevar
!                                    +'LLHV'+suffixe                si writellhv
!                                    +'d'+suffixe                   si writecdl 
!                           et       +'h'+suffixe                   si writecdl 
!
YFLAGWRITE='NEW'      
!
YCHAMP=                                    !nouveau_nom EN MAJUSCULES
NGRIDIA(1:SIZE(XVAR,6))=                   !nouvelle_grille 
CTITRE(1:SIZE(XVAR,6))=                    !nouveau_nom_p pour XMIN XISOLEV diaprog
CCOMMENT(1:SIZE(XVAR,6))=                  !nouveau_comment
CUNITE(1:SIZE(XVAR,6))=                    !nouvelle_unite
XDATIME(1:16,1:SIZE(XVAR,4))=              ! nouvelles dates du fichier
XTRAJT(1:SIZE(XVAR,4),1:SIZE(XVAR,5))=     ! nouveaux timing des champs
XVAR=tableau à ecrire                      ! passé par module aux routines write*
!
!    Redefinir le zoom d ecriture si different du zoom de lecture
!  dans writevar, controle que ce nouveau zoom est inclus dans le zoom de lecture)
!      ivarideb=1
!      ivarifin=1
!      ivarjdeb=1
!      ivarjfin=1
!      ivarkdeb=1
!      ivarkfin=size(XVAR,3)
!      ivartinf=1
!      ivartsup=1
!      ivartrajinf=1
!      ivartrajsup=1
!      ivarprocinf=1
!      ivarprocsup=size(XVAR,6)                   
!
!-------------------------------------------------------------------------------
!
!*      7.     ECRITURE
!              --------
!
SELECT CASE(YTYPEOUT(1:4))
!
  CASE('DIAC')      
      YSUFFIX ='2' ! fichier de sortie= fichier d entree + ysuffix
      !
      ! Traitement par diaprog des champs 2D de type X,Z et Y,Z: 
      ! si le champ 2D XZ correspond à la position j=jpos de la grille 
        !utiliser ivarjdeb=ivarjfin=jpos pour positionner ce champ
      ! idem avec un champ 2D YZ et ipos
      ! si ce champ 2D = moyenne ou compression (trace sans relief)
        ! LJCP=.TRUE.  si  compresse selon j 
        ! LICP=.TRUE.  si  compresse selon i
        ! LKCP=.TRUE.  si  compresse selon k
      IF ( ivarprocinf /= ivarprocsup ) THEN
        ! **** la longueur du nom ne doit pas depasser 9 caracteres si 
        ! plusieurs PROCessus car .PROCn sera ajouté ultérieurement****
        YCHAMP(:)=YCHAMP(1:9) 
      ENDIF
      IF ( SIZE(XVAR,6) /= SIZE(NGRIDIA,1))THEN
        print * ,' *** erreur possible: la dimension6 de XVAR=',SIZE(XVAR,6) ,&
        'est differente de la dimension des tableaux NGRIDIA,CUNIT...'
      ENDIF
      print *,'LICP,LJCP,LKCP,YCHAMP:' , LICP,LJCP,LKCP,YCHAMP,'-'
      print *,'YFLAGWRITE,YFILEIN,YSUFFIX',YFLAGWRITE,YFILEIN,YSUFFIX
  CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup, &
                YCHAMP,YFILEIN,YFLAGWRITE,YSUFFIX,ilocverbia,iret)
  print *, 'retour WRITEVAR=',iret
!
!*     7.2  Ecriture via writellhv
!          ---------------------
!
  CASE('LLHV')
! si YTYPEOUT='LLHV'
! fichier de sortie= fichier d entree + LLHV
   CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                  ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                  YCHAMP,YFILEIN,YFLAGWRITE,YTYPEOUT,&
                  ilocverbia,iret)       
   print *, 'retour WRITELLHV=',iret
!
!*     7.3  Ecriture via writecdl
!          ---------------------
!
  CASE('KCDL')
YSUFFIX='kcl' !fichier de sortie= fichier d entree + ysuffix 
        ! les sorties ZCDL ou PCDL sont codees dans extractdia
        ! conseil: sortir en format 'DIAC' puis utiliser extractdia
        !         pour des sorties ZCDL ou PCDL
   CALL WRITECDL(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                 ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                 YCHAMP,YFILEIN,YFLAGWRITE,'CONF',YSUFFIX,ilocverbia,iret,&
                 PGRIDX=XXX(:,IGRID),PGRIDY=XXY(:,IGRID)                  )
    print *, 'retour WRITEDCL=',iret
!
!*     7.4  Ecriture  format libre
!          ---------------------
!
  CASE('FREE')                 
      ! retour aux unites initiales pour XVAR
      CALL FROM_COMPUTING_UNITS(YCHAMP,CUNITE(1))
      ! coder ici son write Fortran
   print *, 'retour WRITE FREE a coder ='
END SELECT
!
!-------------------------------------------------------------------------------
!
!*     8. boucle possible
! ..... reprise possible des etapes 2 a 8    
     !  pour changer de fichier en lecture : YFLAGREADVAR='OPE' 
     !                                       YFILEIN=' deuxieme fichier'
     !  pour garder le meme fichier en lecture : YFLAGREADVAR='NOP'
     !                                           YCHAMP='autre variable'
YFLAGREADVAR='NOP'
     !  pour changer de fichier en ecriture : YFLAGWRITE='NEW'
     !                                        YSUFFIX='nouveau suffixe'
     !  pour garder le meme fichier en ecriture : YFLAGWRITE='OLD'
YFLAGWRITE='OLD'
!
     ! Pour liberer les unites et ne pas dépasser la limite de 44 fichiers
     ! ouverts simultanement, executer ces 2 lignes des qu un fichier
     ! n est plus utilise
!YFLAGREADVAR='CLO'
!CALL READVAR('',YFILEIN,YFLAGREADVAR,ilocverbia,iret)           
!
!END DO  ! fin boucle lecture+traitement
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       9.    Fermeture fichiers : obligatoire
!              ------------------
!
     ! Au max 44 fichiers simultanement ouverts par FMOPEN (c.a.d via
     ! READVAR et WRITEVAR )
PRINT*, 'Fermeture du fichier d entree'
YFLAGREADVAR='CLO'
CALL READVAR('',YFILEIN,YFLAGREADVAR,ilocverbia,iret)
!
PRINT*, 'Fermeture du fichier de sortie'
YFLAGWRITE='CLO'
SELECT CASE(YTYPEOUT(1:4))
  CASE('DIAC')
    CALL WRITEVAR(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,  &
                YCHAMP,YFILEIN,YFLAGWRITE,YSUFFIX,ilocverbia,iret)
  CASE('LLHV')
    CALL WRITELLHV(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
                   ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
                   YCHAMP,YFILEIN,YFLAGWRITE,YTYPEOUT,&
                   ilocverbia,iret)               
  CASE('KCDL')            
    CALL WRITECDL(ivarideb,ivarifin,ivarjdeb,ivarjfin,ivarkdeb,ivarkfin,&
             ivartinf,ivartsup,ivartrajinf,ivartrajsup,ivarprocinf,ivarprocsup,&
             YCHAMP,YFILEIN,YFLAGWRITE,'CONF',YSUFFIX,ilocverbia,iret, &
             PGRIDX=XXX(:,IGRID),PGRIDY=XXY(:,IGRID)                   )
    ! Remarque: le script "tonetcdf" est lance par writecdl pour obtenir
    !le fichier au format "netcdf" et non au format intermediaire "cdl"
    ! Verifiez que votre PATH donne acces a cette commande 
END SELECT
!
!-------------------------------------------------------------------------------
END PROGRAM EXRWDIA
!
