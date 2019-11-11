!     #################################
      MODULE MODI_CREATLINK
!     #################################
INTERFACE CREATLINK
      SUBROUTINE  CREATLINK (HVARDIR,HFILENAME,HFLAGCREAT,KVERB)
!
CHARACTER(LEN=*)   , INTENT(in)    :: HVARDIR
CHARACTER(LEN=*) , INTENT(inout) :: HFILENAME ! FILENAME (1:28) sera reinit
CHARACTER(LEN=*), INTENT(in)    :: HFLAGCREAT
INTEGER,          INTENT(in)    :: KVERB
!
END SUBROUTINE
END INTERFACE
END MODULE MODI_CREATLINK
!     ################
      SUBROUTINE  CREATLINK (HVARDIR,HFILENAME,HFLAGCREAT,KVERB)
!     ################
!
!!****  *CREATLINK* - 
!! 
!!
!!    PURPOSE
!!    -------
!  Si necessaire, cree un lien symbolique entre le fichier
!  VARDIR/FILENAME et le directory courant ./FILENAME
!  necessaire pour diachro qui ne traite que les fichiers presents
!  dans le directory courant
!
!!**  METHOD
!! 
!    GETENV pour recuperer la valeur de la variable VARDIR qui
!    contient le nom du directory
!    fabrique les commandes UNIX "ln -s dir/file fileloc" avec fileloc=file(1:28)
!                                 rmlink fileloc dir/file"
!    execution de la premiere commande par CALL SYSTEM
!    execution de la seconde commande si HFLAGCREAT=CLEAN
!!    AUTHORS
!!    -------
!!     N. Asencio * CNRM*
!!
!!    Copyright 2003,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!      N. Asencio  sept. 2003  tronque le nom du fichier local à 28 car.
!                             (limite max des routines FM)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
#ifdef NAGf95
USE F90_UNIX  ! for FLUSH
USE F90_UNIX_PROC  ! for SYSTEM
#endif
!
IMPLICIT NONE
!
!*       0.1   Arguments
! 
CHARACTER(LEN=*)   , INTENT(in)    :: HVARDIR
CHARACTER(LEN=*) , INTENT(inout) :: HFILENAME ! FILENAME (1:28) sera reinit
CHARACTER(LEN=*), INTENT(in)    :: HFLAGCREAT
INTEGER,          INTENT(in)    :: KVERB
!
!*       0.2   Local variables
!
INTEGER :: II
CHARACTER (LEN=28) :: yficloc
CHARACTER(LEN=200) :: ydirloc
CHARACTER(LEN=350) :: ycommandlfi
! longueur commande = 'ln -s ' + dirloc +'/'+FILENAME + '.lfi .'
CHARACTER(LEN=350) :: ycommand
! stckage des commandes rm pour l appel avec clean
CHARACTER(LEN=350), dimension (200) , SAVE :: ycleancommand=''
INTEGER,                              SAVE :: icomptclean=0
!
!-------------------------------------------------------------------------------
!
! nom sous le directory local au plus de 28 caracteres
! voir les limites des routines FM de Mesonh
yficloc=ADJUSTL(HFILENAME)
!
!
!*       1.  CLEAN THE LINK
!            --------------
!
IF ( HFLAGCREAT(1:5) == 'CLEAN') THEN
  !
  IF ( HVARDIR == '' .AND. HFILENAME == '' ) THEN
  ! supprime tous les liens
    DO II=1,icomptclean
      IF ( ycleancommand(II) /= '') then
        print *,' execution de ',TRIM(ycleancommand(II))
        CALL SYSTEM (ycleancommand(II))
      END IF
    END DO
  ELSE
    print *,' creatlink option supprime un seul lien ', &
            TRIM(HVARDIR),' ',TRIM(HFILENAME)
    ! supprime un seul lien
    DO II=1,icomptclean
    ! recherche du lien a supprimer, execution de la commande et reinit
      ycommand=ycleancommand(II)
      if ( ycommand(1:29) == 'rmlink ./'//TRIM(yficloc) ) then
        print *,' execution de ',TRIM(ycleancommand(II))
        CALL SYSTEM (ycleancommand(II))
        ycleancommand(II)=''
      else
        IF (KVERB >= 5) THEN
        print *,'ycommand(1:29)= ', ycommand(1:29)
        print *,'rmlink ./'//TRIM(yficloc)
        ENDIF
      endif
    END DO
  !
  ENDIF
!
!*       2.   CREATE THE LINK
!             ---------------
!
ELSE
!
  icomptclean=icomptclean+1
  !  
  ! recupere la valeur de la variable d environnement $VARDIR
  ydirloc= ' '
  CALL GETENV(HVARDIR,ydirloc)
  print *,TRIM(HVARDIR),'=',TRIM(ydirloc)
  !
  IF (ydirloc(1:1) /= ' ' .AND. ydirloc(1:1) /= '.' ) THEN
    ! fichier sous un directory different du directory courant
    IF (HVARDIR == 'DIRLFI') THEN
      ! ajoute .lfi au nom de fichier ( dans ce cas le nom verifie la
      !                                contrainte de 28 car.    )
      ! prepare la creation 
      ycommandlfi=ADJUSTR(HFILENAME)//'.lfi'
      ycommand=ADJUSTR(ydirloc)//'/'//ADJUSTL(ycommandlfi)
      ycommand=TRIM(ycommand)//' .'
      ycommand='ln -s '//ADJUSTL(ycommand)
      ! prepare le nettoyage
      ycleancommand(icomptclean)='rmlink ./'//ADJUSTL(ycommandlfi)
      ycleancommand(icomptclean)=TRIM(ycleancommand(icomptclean))//' '//ADJUSTL(ADJUSTR(ydirloc))
      ycleancommand(icomptclean)=TRIM(ycleancommand(icomptclean))//'/'//ADJUSTL(ycommandlfi)
      IF (KVERB >= 5) THEN
        print *,'cleancommand=' ,TRIM(ycleancommand(icomptclean)) 
      ENDIF
    ELSE
      ! prepare la creation en tronquant a 28 car. le nom local
      ycommand=ADJUSTR(HFILENAME)//' '//ADJUSTL(yficloc)
      ycommand=ADJUSTR(ydirloc)//'/'//ADJUSTL(ycommand)
      ycommand='ln -s '//ADJUSTL(ycommand)
      ! prepare le nettoyage
      !ycleancommand(icomptclean)='rmlink ./'//ADJUSTL(ADJUSTR(yficloc))//&
      !                           ' '//ADJUSTL( TRIM(ydirloc)//'/'//ADJUSTL(ADJUSTR(HFILENAME)) )
      ycleancommand(icomptclean)='rmlink ./'//ADJUSTL(ADJUSTR(yficloc))
      ycommandlfi=TRIM(ydirloc)//'/'//ADJUSTL(ADJUSTR(HFILENAME))
      ycleancommand(icomptclean)=TRIM(ycleancommand(icomptclean))//' '&
                               //ADJUSTL(ycommandlfi)
      print *,'cleancommand=' ,TRIM(ycleancommand(icomptclean)) 
    ENDIF
    print *,' creation du lien :',TRIM(ycommand)
    CALL SYSTEM(ycommand)
 ELSE
   ! fichier deja sous le directory courant: 
   !si longueur du nom est >28 car. creation du lien avec un nom tronque
   !
   IF ( LEN_TRIM(HFILENAME) > 28) THEN
     ! prepare la creation en tronquant a 28 car. le nom local
     ydirloc='.'
     ycommand=ADJUSTR(HFILENAME)//' '//ADJUSTL(yficloc)
     ycommand=ADJUSTR(ydirloc)//'/'//ADJUSTL(ycommand)
     ycommand='ln -s '//ADJUSTL(ycommand)
     ! prepare le nettoyage
     ycleancommand(icomptclean)=TRIM(ydirloc)//'/'//ADJUSTL(ADJUSTR(HFILENAME)) 
     ycleancommand(icomptclean)=ADJUSTR(yficloc)//' '//ADJUSTL(ycleancommand(icomptclean))
     ycleancommand(icomptclean)='rmlink ./'//ADJUSTL(ycleancommand(icomptclean))
     print *,' creation du lien :',TRIM(ycommand)
     CALL SYSTEM(ycommand)                          
   ELSE
     print *,' pas de creation de lien pour ' ,TRIM(HFILENAME)
   ENDIF

 ENDIF
 IF ( LEN_TRIM(HFILENAME) > 28) THEN
     ! reinitialisation du nom passe en argument
     HFILENAME=' '
     HFILENAME(1:28)=yficloc
     print *,' creatlink: reinit du nom du fichier: ', TRIM(HFILENAME)
 ENDIF
!
ENDIF
!
END SUBROUTINE CREATLINK
