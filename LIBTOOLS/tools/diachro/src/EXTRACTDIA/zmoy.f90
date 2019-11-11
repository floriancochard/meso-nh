!     ##################
      MODULE MODI_ZMOY
!     ##################
!
INTERFACE
      SUBROUTINE ZMOY(pvar,KGRID,pmoyz,pvalmin,pvalmax,pundef,KJPVEXT,KJPHEXT)
!
REAL ,  intent(in), dimension (:,:,:) :: pvar     ! champ3D  a traiter
INTEGER , intent(in) :: KGRID                     ! numero de grille du champ
INTEGER , intent(in) :: KJPvext,KJPhext           ! points a exclure
REAL    , intent(in) :: pvalmin,pvalmax           ! definition de la couche
                                                  ! altitude en mètres
REAL    , intent(in) :: pundef                    ! valeur indefinie
REAL ,   intent(out), dimension (:,:) :: pmoyz    ! champ2D moyenné sur la couche      
END SUBROUTINE ZMOY
END INTERFACE
END MODULE MODI_ZMOY
!
!------------------------------------------------------------------------------
!
!     ####################################################
      SUBROUTINE ZMOY(pvar,KGRID,pmoyz,pvalmin,pvalmax,pundef,KJPVEXT,KJPHEXT)
!     ################
!
!!****  *zmoy* - 
!! 
!!
!!    PURPOSE
!!    -------
!   moyenne sur la couche pvalmin,pvalmax
!   pvar peut etre partiellement indefini ( = pundef)
!
!!**  METHOD
!! 
!!    AUTHORS
!!    -------
!!     N. Asencio * CNRM* d apres  evoltempo.f90 J. Stein
!!
!!    Copyright 2003,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!    appel de zinter avec le parametre optionel KNIVMOD
!    toutes les grilles sont traitées par appel à COMPCOORD_FORDIACHRO 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!                    descriptIF grille: XXHAT ,XLAT,XDXHAT,XMAP,XZS,XZZ 
USE MODD_GRID1, ONLY:XZZ
!
USE MODI_ZINTER
IMPLICIT NONE
!*       0.1   Arguments d'appel
REAL ,  intent(in), dimension (:,:,:) :: pvar     ! champ3D  a traiter
INTEGER , intent(in) :: KGRID                     ! numero de grille du champ
INTEGER , intent(in) :: KJPvext,KJPhext           ! points a exclure
REAL    , intent(in) :: pvalmin,pvalmax           ! definition de la couche
                                                  ! altitude en mètres
REAL    , intent(in) :: pundef                    ! valeur indefinie
REAL ,   intent(out), dimension (:,:) :: pmoyz    ! champ2D moyenné sur la couche      

!*       0.2 variables locales
INTEGER :: ji,jj,jk   ! boucles
INTEGER :: ikmin,ikmax   
REAL , allocatable ,dimension (:,:,:)    :: zinterpomin,zinterpomax 
!!                           champs interpolés aux bornes de la couche
INTEGER , allocatable ,dimension (:,:) :: iknivmin,iknivmax
!!                           stockage des premiers niveaux modele
!!                           situés au dessus de chaque borne de la couche 
REAL  :: zhmin
! specIFique a l appel de zinter
REAL , allocatable ,dimension (:)   :: pnivz ! liste des niveaux verticaux
!  ici un seul niveau utilisé mais zinter s'attEND à un tableau 1D
INTEGER :: ikdebmod ! premier niveau modele au dessus du sol
!
!-------------------------------------------------------------------------------
!
!*       1. interpolation sur z=pvalmin et pvalmax et recuperation
!           des niveaux K correspondant
!
IF (.NOT. ALLOCATEd(zinterpomin))  &
  ALLOCATE(zinterpomin(size(pvar,1),size(pvar,2),1))
IF (.NOT. ALLOCATEd(zinterpomax))  &
  ALLOCATE(zinterpomax(size(pvar,1),size(pvar,2),1))
IF (.NOT. ALLOCATEd(iknivmin))  ALLOCATE(iknivmin(size(pvar,1),size(pvar,2)))
IF (.NOT. ALLOCATEd(iknivmax))  ALLOCATE(iknivmax(size(pvar,1),size(pvar,2)))
IF (.NOT. ALLOCATEd(pnivz))  ALLOCATE(pnivz(1))
!
! init du tableau  XZZ pour la grille= KGRID
CALL COMPCOORD_FORDIACHRO(KGRID)
!
ikdebmod=2
pnivz(1)=pvalmin
CALL ZINTER(pvar,XZZ,zinterpomin,pnivz,ikdebmod,pundef,KNIVMOD=iknivmin)
pnivz(1)=pvalmax
CALL ZINTER(pvar,XZZ,zinterpomax,pnivz,ikdebmod,pundef,KNIVMOD=iknivmax)
!
!    en retour de zinter, knivmax= premiers niveaux modele > pvalmax
! pour obtenir les derniers niveaux inclus dans la couche:
WHERE ( iknivmax /= 1+KJPVEXT ) iknivmax=iknivmax-1
!
!-------------------------------------------------------------------------------
!
!*       2. moyenne verticale sur la couche
!
pmoyz=0.
!
! Cumul
!
DO jj=1+KJPHEXT,SIZE(pvar,2)-KJPHEXT
  DO ji=1+KJPHEXT,SIZE(pvar,1)-KJPHEXT
    ikmin=max(iknivmin(ji,jj),1+KJPVEXT)
    ikmax=iknivmax(ji,jj)
    !
    !  borne inferieure de la couche
    !
   IF ( zinterpomin(ji,jj,1) /= pundef .AND. pvar(ji,jj,ikmin) /= pundef ) then
     pmoyz(ji,jj) = &
      0.5*( zinterpomin(ji,jj,1)+pvar(ji,jj,ikmin) )*(XZZ(ji,jj,ikmin)-pvalmin) 
   ENDIF
   !
   !  borne superieure de la couche
   !
   IF ( zinterpomax(ji,jj,1) /= pundef .AND. pvar(ji,jj,ikmax) /= pundef ) then
     pmoyz(ji,jj) = pmoyz(ji,jj) + &
      0.5*( zinterpomax(ji,jj,1)+pvar(ji,jj,ikmax) )*(pvalmax-XZZ(ji,jj,ikmax)) 
   ENDIF
   !
   ! tous les niveaux modele inclus dans la couche
   !
   DO jk=ikmin,ikmax-1
     IF ( pvar(ji,jj,jk) /= pundef .AND. pvar(ji,jj,jk+1) /= pundef ) then
       pmoyz(ji,jj) = pmoyz(ji,jj) + &
        0.5*( pvar(ji,jj,jk) + pvar(ji,jj,jk+1))*(XZZ(ji,jj,jk+1)-XZZ(ji,jj,jk))
     ENDIF
   END DO
   !
   ! calcul de la hauteur utile de la couche
   zhmin=max(pvalmin,XZZ(ji,jj,ikdebmod))
   IF ( pmoyz(ji,jj) /= 0.) pmoyz(ji,jj)=pmoyz(ji,jj)/ (pvalmax-zhmin)
   !
  END DO
END DO
!
! passage a indef des zones ou la moyenne est restee a l init 0.
WHERE ( pmoyz == 0. ) pmoyz=pundef
!
! nettoyage
IF ( ALLOCATED(zinterpomin))  DEALLOCATE(zinterpomin)
IF ( ALLOCATED(zinterpomax))  DEALLOCATE(zinterpomax)
IF ( ALLOCATED(iknivmin))  DEALLOCATE(iknivmin)
IF ( ALLOCATED(iknivmax))  DEALLOCATE(iknivmax)
IF ( ALLOCATED(pnivz))  DEALLOCATE(pnivz)

END SUBROUTINE ZMOY
