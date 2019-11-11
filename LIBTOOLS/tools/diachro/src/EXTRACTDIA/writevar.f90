!     #################################
      MODULE MODI_WRITEVAR
!     #################################
INTERFACE WRITEVAR
      SUBROUTINE  WRITEVAR(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
       ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin,  &
       HLABELCHAMP,HFILENAME,HFLAGFILE,HFILENAME_SUP,KVERBIA,KRETCODE)
!
CHARACTER(LEN=*), INTENT(IN) :: HLABELCHAMP, HFILENAME ! nom du champ et du fichier
CHARACTER(LEN=3), INTENT(IN) :: HFLAGFILE              ! NEW=creation 
                                                       ! OLD=ajout 
                                                       ! CLO=fermeture
CHARACTER(LEN=3)             :: HFILENAME_SUP          ! chaine de caracteres
                                                       ! a rajouter a
                                                       ! HFILENAME
                                                       ! si ='NEN' alors HFILENAME
                                                       ! contient le nom complet
INTEGER , INTENT(IN)         :: KVERBIA                ! prints de controle
!
INTEGER , intent(in)         :: kideb,kifin,kjdeb,kjfin,kkdeb,kkfin   
INTEGER , intent(in)         :: ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin
!
INTEGER  , INTENT(OUT)       :: KRETCODE  ! Code de retour de la routine 
!
END SUBROUTINE
END INTERFACE
END MODULE MODI_WRITEVAR
!     ######
      SUBROUTINE  WRITEVAR(kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
       ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin, &
       HLABELCHAMP,HFILENAME,HFLAGFILE,HFILENAME_SUP,KVERBIA,KRETCODE)
!     ################
!
!!****  *WRITEVAR* - 
!! 
!!
!!    PURPOSE
!!    -------
!     Ecriture d'un fichier  de type:
!       diachronique en vue d'un traitement via diaprog
! 
!
!!**  METHOD
!!    ------
!     utilisation des routines de diaprog : le tableau de stockage
!     XVAR est alloué avant l appel a WRITEVAR
!
!     HFLAGFILE='NEW' lors de la premiere utilisation du fichier
!     HFLAGFILE='OLD' lors des utilisations suivantes
!     HFLAGFILE='CLO' pour la fermeture du fichier de sortie
!      ( fin de mise a jour du menu )
!
!     KVERBIA= 0 impressions reduites au minimum (entree et sortie de la
!      routine)
!     KVERBIA >0 impressions pour signaler chaque etape de WRITEVAR
!
!     KRETCODE = 0 execution de WRITEVAR correcte
!     KRETCODE = 1 erreur lors de l ouverture du fichier
!     KRETCODE = 2 erreur lors de l ecriture du champ 
!     KRETCODE = 3 erreur lors de la fermeture du champ 
!     KRETCODE = -1 pas de fermeture car pas d ouverture
!
!     kideb,kifin,kjdeb,kjfin,kkdeb,kkfin = limites en indices i,j,k du
!       domaine à traiter dans XVAR       
!     ktdeb,ktfin,ktrdeb,ktrfin = limites en indices
!       des dimensions 4,5 de XVAR  
!
!!    EXTERNAL
!!    --------
!!          FROM_COMPUTING_UNITS: retour aux unites initiales  avant ecriture
!!                               = passage inverse a celui realise par
!!                                 TO_COMPUTING_UNITS
!!                              
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHORS
!!    -------
!!    I. Mallet , N. Asencio , J. Stein * CNRM*
!!
!!    Copyright 2003,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/03/2003
!       N. Asencio  01/2005 : take in account 2D fields XZ, YZ and
!                             zoomed fields inside the complete x-y-z-grid
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! modules MESONH
USE MODD_CST
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT
!                    NIMAX,NJMAX,NKMAX,NIINF, NISUP
USE MODD_DIM1
USE MODD_GRID, ONLY: XLAT0,XLON0,XRPK,XBETA
!                    descriptif grille: XXHAT(:) ,XLAT(:,:),XDXHAT(:),XMAP(:,:)
!                    ,XZS(:,:),XZZ(:,:,:) ,XCOSSLOPE(:,:),XDIRCOSXW(:,:)
USE MODD_GRID1
!      
! modules DIACHRO
USE MODN_NCAR,  ONLY: XSPVAL    
USE MODD_COORD ! grille : XXDXHAT(:,1:7) et XXX(:,1:7), XXZS(:,:,1:7)
USE MODD_TYPE_AND_LH ! zoom selon x et y et z :  NIL,NIH,NJL,NJH,NKL,NKH,CTYPE
USE MODD_ALLOC_FORDIACHRO ! XVAR(i,j,k,,,), XMASK,XTRAJ ,XDATIME(16,t)   
USE MODD_OUT ! nom de fichiers NLUOUT,CLFIFM, CDESFM
USE MODD_FILES_DIACHRO ! NBFILES + nom des fichiers CFILEDIAS, CLUOUTDIAS
!                    pour l appel a WRITE_DIMGRIDREF, FMATTR et FMCLOS
USE MODD_DIACHRO, ONLY:CFILEDIA,CLUOUTDIA, &
                       NLUOUTDIA,NRESPDIA,NNPRARDIA,NFTYPEDIA,NVERBDIA,NNINARDIA
USE MODD_READLH
!
USE MODI_WRITE_DIMGRIDREF      
USE MODI_WRITE_DIACHRO      
USE MODI_MENU_DIACHRO
USE MODI_FROM_COMPUTING_UNITS
! 
!
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!              ----------------
!
CHARACTER(LEN=*), INTENT(IN) :: HLABELCHAMP, HFILENAME ! nom du champ et du fichier
CHARACTER(LEN=3), INTENT(IN) :: HFLAGFILE              ! NEW=creation 
                                                       ! OLD=ajout 
                                                       ! CLO=fermeture
CHARACTER(LEN=3)             :: HFILENAME_SUP          ! chaine de caracteres
                                                       !a rajouter a HFILENAME
                                                       ! si ='NEN' alors HFILENAME
                                                       ! contient le nom complet
INTEGER , INTENT(IN)         :: KVERBIA                ! prints de controle
!
INTEGER , intent(in)         :: kideb,kifin,kjdeb,kjfin,kkdeb,kkfin   
INTEGER , intent(in)         :: ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin
!
INTEGER , INTENT(OUT)        :: KRETCODE  ! Code de retour de la routine 
!
!*       0.2   Declarations des variables locales
!              -----------------------------------
!
INTEGER           ::  ISAVENGRIDIA,iret
!                 repositionne le zoom/grille si zoom d un champ deja zoome
INTEGER          :: InewIL,InewIH,InewJL,InewJH,InewKL,InewKH
!
REAL ,DIMENSION(:,:,:,:,:,:) , ALLOCATABLE :: ZVARZS,& ! stockage dans
                                                       ! un tableau 6d de ZS 
                                                       ! avant son ecriture
                                              ZVARSAVE ! sauvegarde de XVAR
!
! taille=100  et 28 cf diachro 
CHARACTER (LEN=100) :: YSAVETITRE, YSAVECOMMENT, YSAVEUNITE 
CHARACTER (LEN=28), SAVE  :: YFILEOUT='zadefinir'        ! Fichier de sortie
CHARACTER (LEN=28)  :: YSAVEFILEDIA             ! sauve le contenu de CFILEDIA 
CHARACTER (LEN=3)   :: YFLAGZS 
CHARACTER (LEN=3)   :: YFLAGFILE 
!
INTEGER,SAVE   ::   IGROUP=0  ! pour compter le nb de champs ecrits
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
!      
! Code de retour de la routine : 0 = OK
!                                1 = erreur lors de l ouverture du fichier
KRETCODE=0
!
YFLAGFILE=HFLAGFILE
!
if (KVERBIA >= 0) then
  print *,'--------- '
  print *,'Beginning of WRITEVAR ',TRIM(HFILENAME),' ',TRIM(HLABELCHAMP),' ',&
                             TRIM(YFLAGFILE)  ,' ',&
                             TRIM(HFILENAME_SUP),' ',KVERBIA
endif
!
! code de retour d erreur des routines diaprog
LPBREAD=.FALSE.                                                        
!
!*       1.1    Determine le nom du fichier de sortie au premier passage
!              -------------------
!
IF (YFILEOUT=='zadefinir') THEN
  ! alignement à droite pour que le test LEN(YFILEOUT)-1:LEN(YFILEOUT)) == '.Z' fonctionne
  YFILEOUT=(ADJUSTR(HFILENAME))
  IF (HFILENAME_SUP(1:3) /= 'NEN' ) THEN
  ! cas d un appel obs2mesonh 
  !avec redefinition totale du nom de fichier de sortie (on prend HFILENAME tel quel)
    IF (HFILENAME_SUP(1:3)=='SAM') THEN
    ! cas d un appel dans compute_r00pc
     ! pas d ajout de suffixe (on complete un fichier existant ouvert en 'OLD')
     ! m.a.j. de la liste des enregistrements diachroniques
      CALL MENU_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),'READ')
      IF(YFLAGFILE(1:3)/='CLO') YFLAGFILE='OLD'
    ELSE
      ! ajout d un suffixe 2
      IF (LEN_TRIM(HFILENAME_SUP) == 0) HFILENAME_SUP='2  '
      !  
      IF ( YFILEOUT(LEN(YFILEOUT)-1:LEN(YFILEOUT)) == '.Z' ) THEN
        ! ajout du suffixe devant le .Z
        ! et suppression de .Z car le fichier cree sera non compresse
        YFILEOUT=ADJUSTL(YFILEOUT(1:LEN(YFILEOUT)-2)//HFILENAME_SUP)
      ELSE
        ! ajout en fin de nom 
        YFILEOUT=ADJUSTL(YFILEOUT(1:LEN(YFILEOUT))//HFILENAME_SUP)
      END IF
    END IF
  ENDIF
  YFILEOUT=ADJUSTL(YFILEOUT)
END IF
!    
if (KVERBIA > 0) then
  PRINT*,'WRITEVAR: output diachronic file ',YFILEOUT
endif
! 
!*       1.2    Appel avec fichier courant different du fichier a ecrire
!              -------------------
!          cas possibles dans compute_r00pc, exrwdia et obs2mesonh avec HFILENAME_SUP(1:3) /= 'NEN'
!
IF ( YFLAGFILE(1:3) /= 'CLO'  ) THEN      
!   reinit eventuelle de l entete si fichier courant different du fichier a ecrire
  YSAVEFILEDIA=CFILEDIA
  IF ( YSAVEFILEDIA /= HFILENAME  .AND. HFILENAME_SUP(1:3) /= 'NEN' ) THEN
   ! seul le cas compute_r00pc est concerné
   ! dans le cas  obs2mesonh avec HFILENAME_SUP(1:3) /= 'NEN', la reinit de 
   ! l entete  (date et heure) a été  faite dans obs2mesonh
    if (KVERBIA > 0) then
      print *,'WRITEVAR: fichier courant dans READVAR ',YSAVEFILEDIA
      print *,' different du fichier a ecrire ', HFILENAME
      print *,' seul XVAR est sauve. La grille spatiale est supposée identique.'
    endif
    ISAVENGRIDIA=NGRIDIA(1)
    YSAVETITRE=CTITRE(1)
    YSAVECOMMENT=CCOMMENT(1)
    YSAVEUNITE=CUNITE(1)
    ! lecture d un champ de HFILENAME pour reinitialiser les modules diachro
    !pour creer l en tete du fichier de sortie YFILEOUT(HFILENAME)
    ALLOCATE(ZVARSAVE(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),&
                      SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)) )
    ZVARSAVE=XVAR
    YFLAGZS='NOP'
    CALL READVAR ('ZSBIS',HFILENAME,YFLAGZS,KVERBIA,iret)
    if (KVERBIA > 0) then
      print *,'WRITEVAR: apres reinit des modules pour le fichier ',HFILENAME
    endif
    DEALLOCATE(XVAR)
    ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),SIZE(ZVARSAVE,3),&
                  SIZE(ZVARSAVE,4),SIZE(ZVARSAVE,5),SIZE(ZVARSAVE,6)) )
    XVAR=ZVARSAVE
    NGRIDIA(1)=ISAVENGRIDIA
    CTITRE(1)=YSAVETITRE            
    CCOMMENT(1)=YSAVECOMMENT
    CUNITE(1)=YSAVEUNITE            
  ENDIF
  CFILEDIA=ADJUSTL(YFILEOUT)
ENDIF
!      
!-------------------------------------------------------------------------------
!
!*       2.    Ouverture du fichier de sortie
!              -------------------
!      
IF ( YFLAGFILE(1:3) /= 'CLO' ) THEN
! Repositionne eventuellement le zoom  en I et J , pour K (2 cas)
InewIL=max(NREADIL,kideb)
InewJL=max(NREADJL,kjdeb)
InewKL=max(NREADKL,kkdeb)
InewIH=min(NREADIH,kifin)
InewJH=min(NREADJH,kjfin)
InewKH=min(NREADKH,kkfin)
IF ( NREADKL == NREADKH .AND. SIZE(XVAR,3) > 1 )THEN
   ! en lecture le tableau contient un seul niveau vertical
   ! en ecriture le tableau (autre variable) contient plusieurs niveaux: 
   ! ecriture du zoom utilisateur
   InewKL=kkdeb
   InewKH=kkfin
   print *, '* warning: desaccord sur le zoom selon la verticale'
   print *, ' le zoom lu=',NREADKL,NREADKH ,'et le zoom ecrit=',kkdeb,kkfin
   ! Pour des traces diaprog sur ce nouveau zoom
   NREADKL=kkdeb
   NREADKH=kkfin
ENDIF
  if (KVERBIA > 1) then
    print*,'ancienne localisation du champ/grille :',NREADIL,NREADIH,NREADJL,NREADJH,NREADKL,NREADKH
    print*,' zoom demande: ', kideb,kifin,kjdeb,kjfin,kkdeb,kkfin
    print*,'nouvelle localisation du champ/grille :',&
             InewIL,InewIH,InewJL,InewJH,InewKL,InewKH
  endif
ENDIF
!      
IF ( YFLAGFILE(1:3) == 'NEW' ) THEN
  !
  CLUOUTDIA=CLUOUTDIAS(NBFILES)
  NLUOUTDIA=NLUOUTDIAS(NBFILES)
  if (KVERBIA >0)then
      print *,'WRITEVAR: avant OPEN_FILES ',TRIM(YFILEOUT),' ',TRIM(CFILEDIA), &
                                            ' ',TRIM(CLUOUTDIA)
  endif
  !
  if (KVERBIA > 1) then
    print *,'WRITEVAR: lat0,lon0 ',XLAT0,XLON0
  endif
  !      Ouverture et ecriture de l entete
  CALL WRITE_DIMGRIDREF
  IF (NRESPDIA.NE.0)THEN
    KRETCODE=1
    print *,' ****WRITEVAR: erreur lors de l ouverture du fichier ',&
            YFILEOUT, 'code= ',NRESPDIA
    RETURN
  ENDIF 
  !
  IF (TRIM(HLABELCHAMP)/='ZSBIS') THEN
  ! Ecriture de ZS avec le nom ZSBIS necessaire pour tracer
  !  le champ "ZS" dans diaprog
    ALLOCATE(ZVARZS(SIZE(XZS,1),SIZE(XZS,2),1,1,1,1))
    ZVARZS(:,:,1,1,1,1)=XZS
    ISAVENGRIDIA=NGRIDIA(1)
    YSAVETITRE=CTITRE(1)
    YSAVECOMMENT=CCOMMENT(1)
    YSAVEUNITE=CUNITE(1)
    NGRIDIA(1)=4
    CTITRE(1)='ZSBIS'
    CUNITE(1)='m'
    CCOMMENT(1)='X_Y_ZS (m)' 
    CALL WRITE_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),'ZSBIS','CART',NGRIDIA,&
                       XDATIME,ZVARZS(kideb:kifin,kjdeb:kjfin,:,:,:,:),&
                       XTRAJT,CTITRE,CUNITE,CCOMMENT, &
                       .FALSE.,.FALSE.,.FALSE.,InewIL,InewIH,InewJL,InewJH,1,1)
    if (KVERBIA > 0) then
      print *,'WRITEVAR(zs) size= 1:',size(ZVARZS,1),',1:',size(ZVARZS,2)
      print *,'  InewIL,InewIH,InewJL,InewJH,1,1=', InewIL,InewIH,InewJL,InewJH
    end if
    DEALLOCATE(ZVARZS)
    NGRIDIA(1)=ISAVENGRIDIA
    CTITRE(1)=YSAVETITRE
    CCOMMENT(1)=YSAVECOMMENT
    CUNITE(1)=YSAVEUNITE                   
    if (KVERBIA > 1) then
      print *,'WRITEVAR: apres write_diachro ZSBIS'
    endif
  !
    IF (NRESPDIA.NE.0)THEN
      KRETCODE=2
      print *,' ****WRITEVAR: erreur lors de l ecriture de ZS  dans ',&
            YFILEOUT, ' code= ',NRESPDIA
      RETURN
    ELSE 
      IGROUP=IGROUP+1
    ENDIF 
  !
  ENDIF 
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4     Ecriture du champ dans YFILEOUT
!              -------------------
!
IF ( YFLAGFILE(1:3) /= 'CLO' ) THEN
  !
  if (KVERBIA >= 0) then
    print*,'WRITEVAR: ecriture en cours de ',HLABELCHAMP
  endif
  !  Retour aux unites initiales si necessaire
  CALL FROM_COMPUTING_UNITS(HLABELCHAMP,CUNITE(1)) 
  !
  if (KVERBIA > 1) then
    print*,'WRITEVAR: NGRID,NGRIDIA(:) =',NGRID,NGRIDIA
  endif
  !
  IF ( SIZE(XVAR,6) /= SIZE(NGRIDIA,1))THEN
    print * ,' *** erreur possible: la dimension6 de XVAR=',SIZE(XVAR,6) ,&
             'est differente de la dimension des tableaux NGRIDIA,CUNIT...'
  ENDIF
  IF ( SIZE(XVAR,4) /= SIZE(XDATIME,2))THEN
    print * ,' *** erreur possible: la dimension4 de XVAR=',SIZE(XVAR,4) ,&
             'est differente de la dimension des tableaux XDATIME,XTRAJT...'
  ENDIF
  !
  IF (ALLOCATED(XMASK)) THEN
    ! CTYPE='MASK'
    IF ( SIZE(XVAR,5) /= SIZE(XMASK,5))THEN
      print * ,' *** erreur possible: la dimension5 de XVAR=',SIZE(XVAR,5) ,&
               'est differente de la dimension5 du tableau XMASK'
    ENDIF
    CALL WRITE_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),HLABELCHAMP,CTYPE,   &
                       NGRIDIA(kpdeb:kpfin),XDATIME(:,ktdeb:ktfin),     &
                       XVAR(kideb:kifin,kjdeb:kjfin,kkdeb:kkfin,& 
                            ktdeb:ktfin,ktrdeb:ktrfin,kpdeb:kpfin),&
                       XTRAJT(ktdeb:ktfin,:),CTITRE(kpdeb:kpfin),&
                       CUNITE(kpdeb:kpfin),CCOMMENT(kpdeb:kpfin), &
                       LICP,LJCP,LKCP,InewIL,InewIH,InewJL,InewJH,InewKL,InewKH,&
   !                   LICP,LJCP,LKCP,kideb,kifin,kjdeb,kjfin,kkdeb,kkfin,&
                       PMASK=XMASK)
  ELSE IF (ALLOCATED(XTRAJX).AND.ALLOCATED(XTRAJY).AND.ALLOCATED(XTRAJZ))THEN
    IF ( CTYPE=='SSOL' ) THEN
      CALL WRITE_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),HLABELCHAMP,CTYPE,   &
                       NGRIDIA(kpdeb:kpfin),XDATIME(:,ktdeb:ktfin),  &
                       XVAR(kideb:kifin,kjdeb:kjfin,kkdeb:kkfin,& 
                            ktdeb:ktfin,ktrdeb:ktrfin,kpdeb:kpfin),&
                       XTRAJT(ktdeb:ktfin,:),CTITRE(kpdeb:kpfin),&
                       CUNITE(kpdeb:kpfin),CCOMMENT(kpdeb:kpfin), &
                       PTRAJX=XTRAJX,PTRAJY=XTRAJY,               &
                       PTRAJZ=XTRAJZ(kkdeb:kkfin,1:1,ktrdeb:ktrfin))
    ELSE
    ! CTYPE='DRST' or CTYPE='RSPL' or CTYPE='RAPL'
      CALL WRITE_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),HLABELCHAMP,CTYPE,   &
                       NGRIDIA(kpdeb:kpfin),XDATIME(:,ktdeb:ktfin),  &
                       XVAR(kideb:kifin,kjdeb:kjfin,kkdeb:kkfin,& 
                            ktdeb:ktfin,ktrdeb:ktrfin,kpdeb:kpfin),&
                       XTRAJT(ktdeb:ktfin,:),CTITRE(kpdeb:kpfin),&
                       CUNITE(kpdeb:kpfin),CCOMMENT(kpdeb:kpfin), &
                       PTRAJX=XTRAJX,PTRAJY=XTRAJY,               &
                       PTRAJZ=XTRAJZ(kkdeb:kkfin,ktdeb:ktfin,ktrdeb:ktrfin))
    ENDIF
  ELSE IF (.NOT.ALLOCATED(XTRAJX) .AND. .NOT.ALLOCATED(XTRAJY) .AND. .NOT.ALLOCATED(XTRAJZ))THEN
    ! CTYPE='CART' or CTYPE='SPXY'
    CALL WRITE_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),HLABELCHAMP,CTYPE,   &
                       NGRIDIA(kpdeb:kpfin),XDATIME(:,ktdeb:ktfin),      &
                       XVAR(kideb:kifin,kjdeb:kjfin,kkdeb:kkfin,& 
                            ktdeb:ktfin,ktrdeb:ktrfin,kpdeb:kpfin),&
                       XTRAJT(ktdeb:ktfin,:),CTITRE(kpdeb:kpfin),&
                       CUNITE(kpdeb:kpfin),CCOMMENT(kpdeb:kpfin), &
                       LICP,LJCP,LKCP,InewIL,InewIH,InewJL,InewJH,InewKL,InewKH)
    !                  LICP,LJCP,LKCP,kideb,kifin,kjdeb,kjfin,kkdeb,kkfin)
  ELSE
    KRETCODE=2
    print *,' ****WRITEVAR: cas d ecriture non prevu pour ',HLABELCHAMP
    RETURN
  ENDIF 
  if (KVERBIA > 0) then
    print *,'WRITEVAR(champ)'
    print *,'  ideb,ifin,jdeb,jfin,kdeb,kfin=', &
        kideb,kifin,kjdeb,kjfin,kkdeb,kkfin
    print *,'  tdeb,tfin,trdeb,trfin,pdeb,pfin=',&
        ktdeb,ktfin,ktrdeb,ktrfin,kpdeb,kpfin
  end if
  if (KVERBIA > 1) then
    print*,'WRITEVAR: apres write_diachro, CTYPE=',CTYPE,' xdatime(16,ktdeb:ktfin)'
    do iret=ktdeb,ktfin
     print*, iret,' ',XDATIME(1:4,iret)
     print*, XDATIME(5:8,iret)
     print*, XDATIME(9:12,iret)
     print*, XDATIME(13:16,iret)
    end do
  endif
  IF (NRESPDIA.NE.0)THEN
    KRETCODE=2
    print *,' ****WRITEVAR: erreur lors de l ecriture de ',HLABELCHAMP,&
            ' dans ',YFILEOUT, ' code= ',NRESPDIA
    RETURN
  ELSE 
    IGROUP=IGROUP+1
  ENDIF 
  !
  CFILEDIA=YSAVEFILEDIA
  IF ( YSAVEFILEDIA /= HFILENAME .AND. HFILENAME_SUP(1:3) /= 'NEN') THEN
    ! retablit les infos du fichier courant
    if (KVERBIA > 0) then
      print *,'WRITEVAR: avant retour aux infos des modules pour ',&
              ' le fichier courant ', YSAVEFILEDIA
    endif
    !      
    YFLAGZS='NOP'
    CALL READVAR ('ZSBIS',YSAVEFILEDIA,YFLAGZS,KVERBIA,iret)
    DEALLOCATE(XVAR)
    ALLOCATE(XVAR(SIZE(ZVARSAVE,1),SIZE(ZVARSAVE,2),SIZE(ZVARSAVE,3),&
             SIZE(ZVARSAVE,4),SIZE(ZVARSAVE,5),SIZE(ZVARSAVE,6)) )
    XVAR=ZVARSAVE
    DEALLOCATE(ZVARSAVE)
  ENDIF
  if (KVERBIA >= 0) then
    print *,'--------- '
  endif
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4     FERMETURE  des fichiers de sortie
!              ---------------------------------
!
IF ( YFLAGFILE(1:3) == 'CLO' ) THEN
  if (KVERBIA > 0 .AND. IGROUP>0) then
    print *,'WRITEVAR: before closing the output file ',TRIM(YFILEOUT)
    print *,' List of the ',IGROUP,' variables :'
  endif
  !
  ! fichier de sortie
  CALL MENU_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),'END')
  if (KVERBIA > 0 .AND. IGROUP>0) then
    CALL MENU_DIACHRO(YFILEOUT,CLUOUTDIAS(NBFILES),'READ')
  endif
  IF (IGROUP>0) THEN
    CALL FMCLOS(YFILEOUT,'KEEP',CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
    !           
    if (NRESPDIAS(NBFILES)==0) then
      print*,'End of WRITEVAR: file ',TRIM(YFILEOUT),' available '
      print *,'--------- '
    else
      print *,' ****WRITEVAR: error when closing the file ',&
              TRIM(YFILEOUT), ' code= ',NRESPDIAS(NBFILES)
      KRETCODE=3
    endif
  ELSE
    print *,' ****WRITEVAR: file not opened, so no closing'
    KRETCODE=-1
  END IF
  ! pour determination du nom du fichier de sortie au prochain appel
  YFILEOUT='zadefinir'
  IGROUP=0
  !
ENDIF
!
END SUBROUTINE WRITEVAR
