!     ############################################################
      MODULE MODI_FF
!     ############################################################
!
INTERFACE
      SUBROUTINE FF(PU,PV,PFFVENT,KVEXT,KHEXT,KGRID)
!
REAL, INTENT(IN), DIMENSION (:,:,:) :: PU,PV    ! composantes u et V
INTEGER, INTENT(IN) :: KVEXT,KHEXT              ! points a exclure
REAL, INTENT(INOUT), DIMENSION (:,:,:) :: PFFVENT ! module vent
INTEGER, INTENT(IN) :: KGRID                    ! grille des champs PU,PV
!
END SUBROUTINE FF
END INTERFACE
END MODULE MODI_FF
!
!------------------------------------------------------------------------------
!

!     ################
      SUBROUTINE FF(PU,PV,PFFVENT,KVEXT,KHEXT,KGRID)
!     ################
!
!!****  *FF* - 
!! 
!!
!!    PURPOSE
!!    -------
!  calcul du module du vent
!
!!**  METHOD
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_NCAR,  ONLY: XSPVAL
!
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!
REAL, INTENT(IN), DIMENSION (:,:,:) :: PU,PV    ! composantes u et V
INTEGER, INTENT(IN) :: KVEXT,KHEXT              ! points a exclure
REAL, INTENT(INOUT), DIMENSION (:,:,:) :: PFFVENT ! module vent
INTEGER, INTENT(IN) :: KGRID                    ! grille des champs PU,PV
!
!*       0.2 variables locales
!
INTEGER :: JI,JJ,JK  ! loop indexes
INTEGER :: JK1,JK2
!
!-------------------------------------------------------------------------------
!
IF (SIZE(PU,3) == 1) THEN
    JK1=1 
    JK2=1
ELSE
    JK1=1+KVEXT
    JK2=SIZE(PU,3)-KVEXT
ENDIF
IF (KGRID == 1 ) THEN
  ! les 2 composantes sont au point de masse UM10,VM10 ou colocalisées
  ! apres interpolation horizontale
  DO JK=JK1,JK2
    DO JJ=1+KHEXT,SIZE(PU,2)-KHEXT
      DO JI=1+KHEXT,SIZE(PU,1)-KHEXT
      ! calcul de la force du vent
        IF ( PU(JI,JJ,JK) /= XSPVAL .AND. PV(JI,JJ,JK) /= XSPVAL) then
          PFFVENT(JI,JJ,JK)=sqrt( PU(JI,JJ,JK)**2+ PV(JI,JJ,JK)**2 )
        ELSE
          PFFVENT(JI,JJ,JK)=XSPVAL
        ENDIF
      end DO
    end DO
  end DO
ELSE
  ! les 2 composantes sont dans les grilles U et V Mesonh
  DO JK=JK1,JK2
    DO JJ=1+KHEXT,SIZE(PU,2)-KHEXT
      DO JI=1+KHEXT,SIZE(PU,1)-KHEXT
      ! calcul de la force du vent
        IF (PU(JI,JJ,JK) /= XSPVAL .AND. PU(JI+1,JJ,JK) /= XSPVAL &
           .AND. PV(JI,JJ,JK) /= XSPVAL .AND. PV(JI,JJ+1,JK) /= XSPVAL) then
          PFFVENT(JI,JJ,JK) = sqrt(0.25*(PU(JI,JJ,JK)+PU(JI+1,JJ,JK))**2+   &
                                   0.25*(PV(JI,JJ,JK)+PV(JI,JJ+1,JK))**2    )
        ELSE
          PFFVENT(JI,JJ,JK) = XSPVAL
        ENDIF
      end DO
    end DO
  end DO
ENDIF
!
END SUBROUTINE FF
