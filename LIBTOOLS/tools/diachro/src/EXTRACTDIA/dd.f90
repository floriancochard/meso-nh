!     ############################################################
      MODULE MODI_DD
!     ############################################################
!
INTERFACE
      SUBROUTINE DD(pu,pv,pddvent,kiskip,KGRID,PLON)
!
REAL    , intent(in), dimension (:,:,:) :: pu,pv ! composantes u et v
INTEGER , intent(in) :: kiskip                   ! nb points a sauter
INTEGER , intent(in) :: KGRID                    ! grille des champs u et v
REAL    , intent(inout), dimension (:,:,:) :: pddvent ! direction vent
REAL    ,intent(in), dimension (:,:),OPTIONAL   :: PLON ! tableau des lon
!
END SUBROUTINE DD
END INTERFACE
END MODULE MODI_DD
!
!------------------------------------------------------------------------------
!
!     ################
      SUBROUTINE DD(pu,pv,pddvent,kiskip,KGRID,PLON)
!     ################
!
!!****  *DD* - 
!! 
!!
!!    PURPOSE
!!    -------
!  calcul de la direction du vent par rapport au Nord geographique
!  0=360 pour un vent venant du Nord
!
!!**  METHOD
!  Appel de computedir niveau vertical par niveau vertical
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
!!      call to change_a_grid  15/04/2004  (I.Mallet) 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_CHANGE_A_GRID
USE MODI_COMPUTEDIR
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!
REAL    , intent(in), dimension (:,:,:) :: pu,pv ! composantes u et v
INTEGER , intent(in) :: kiskip                   ! nb points a sauter
INTEGER , intent(in) :: KGRID                    ! grille des champs u et v
REAL,  intent(inout), dimension (:,:,:) :: pddvent ! direction vent
REAL    ,intent(in), dimension (:,:) , OPTIONAL  :: PLON ! tableau des lon
!
!*       0.2 variables locales
!
INTEGER :: JK,IGRID
REAL, allocatable , dimension (:,:)   :: zwork2d
REAL, allocatable , dimension (:,:,:) :: zwork3du
!
!-------------------------------------------------------------------------------
!
print *,'entree dd ',kiskip,SIZE(pu,3),SIZE(pu,1),SIZE(pu,2)
!
ALLOCATE(zwork3du(size(pu,1),size(pu,2),size(pu,3)))
IF (KGRID /= 1 ) THEN
  ! les 2 composantes sont dans les grilles U(2) et V(3) Mesonh
  IGRID=2
  CALL CHANGE_A_GRID(PU,IGRID,zwork3du)
  IGRID=3
  CALL CHANGE_A_GRID(PV,IGRID,pddvent)
ELSE
  zwork3du(:,:,:)=PU(:,:,:)
  pddvent (:,:,:)=PV(:,:,:)
ENDIF
! 
! Tableau de travail : 2D pour computedir
ALLOCATE(zwork2d(size(pu,1),size(pu,2)))
! 
! Calcul niveau par niveau et stockage dans le tableau 3D
!
IF (PRESENT(PLON)) THEN
  ! grille lon (passee en arg.) differente de celle de Mesonh
  print *,' dd: grille lon utilisateur'
  do JK=1,SIZE(pu,3)
    zwork2d(:,:)=pddvent(:,:,JK)
    CALL COMPUTEDIR (size(PU,1),size(PU,2),size(PV,1),size(PV,2),   &
                     kiskip,zwork3du(:,:,JK),zwork2d(:,:), PLO=PLON )
    pddvent(:,:,JK)=zwork2d(:,:)
  end do
ELSE
  print *,' dd: grille lat lon mesonh'
  ! computedir recalculera PLO en fonction de XXHAT et XYHAT
  do JK=1,SIZE(pu,3)
    zwork2d(:,:)=pddvent(:,:,JK)
    CALL COMPUTEDIR (size(PU,1),size(PU,2),size(PV,1),size(PV,2), &
                     kiskip,zwork3du(:,:,JK),zwork2d(:,:)         )
    pddvent(:,:,JK)=zwork2d(:,:)
  end do       
ENDIF
DEALLOCATE(zwork2d,zwork3du)
!
END SUBROUTINE DD
