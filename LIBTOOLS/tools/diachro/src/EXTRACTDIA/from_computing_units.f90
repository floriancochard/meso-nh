!     ############################################################
      MODULE MODI_FROM_COMPUTING_UNITS
!     ############################################################
!
INTERFACE
      SUBROUTINE FROM_COMPUTING_UNITS(HCHAMP,HUNITS)
!
CHARACTER(LEN=*) , intent(in)    :: HCHAMP     ! Nom du champ 
CHARACTER(LEN=*) , intent(inout) :: HUNITS     ! Unite
!
END SUBROUTINE FROM_COMPUTING_UNITS
END INTERFACE
END MODULE MODI_FROM_COMPUTING_UNITS
!
!------------------------------------------------------------------------------
!
!     ################
      SUBROUTINE FROM_COMPUTING_UNITS(HCHAMP,HUNITS)
!     ################
!
!!****  *FROM_COMPUTING_UNITS* - 
!! 
!!
!!    PURPOSE
!!    -------
!  Retour vers l'unite initiale apres un passage a une unite adaptee au calcul
!  dans la routine To_Computing_Units(YCHAMP,CUNIT) 
!
!!**  METHOD
!   mettre a jour suivant les variables Mesonh qui necessitent ce passage
! AU 01/2005 : les reflectivités radarexprimees en dBz
!              les temperatures de brillance
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
!!      Original  25/01/2005  (N. Asencio)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY:XUNDEF
USE MODD_ALLOC_FORDIACHRO, ONLY: XVAR
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!
CHARACTER(LEN=*) , intent(in)    :: HCHAMP     ! Nom du champ 
CHARACTER(LEN=*) , intent(inout) :: HUNITS     ! Unite
!
!*       0.2 variables locales
!
!
!-------------------------------------------------------------------------------
!
!print *,'entree FROM_COMPUTING_UNITS ',TRIM(HCHAMP),' ',TRIM(HUNITS)
!
!
! Critère= Unite modifiee dans To_Computing_Units
!
SELECT CASE (HUNITS)
  CASE ('Ze_to_DBZ','Ze_listOBS')
    ! Reflectivités radar
    WHERE ( XVAR <= 0. ) XVAR=XUNDEF
    WHERE ( XVAR /= XUNDEF ) XVAR=10*alog10(XVAR) 
    ! Retour a l unite initiale
    HUNITS='dBZ'
    print *,'**** FROM_COMPUTING_UNITS:Passage Ze a DBZ avant ecriture ****'
  CASE ('W_to_C')
         ! finir les modd voulus et utiliser cet appel
         ! Mesonh
         ! passage rad -> temp brillance pour le satellite KGEO
         ! call MAKE_RADSAT(KYEARF, KMONTHF, KDAYF, PSECF, &
         !                  KGEO, KLON, PRADB, PRADF)
         ! Viviane
         !ZOBS est en radiance, je la transforme en tempe de brillance
         ! IF (ZRADMOY > 0. .AND. (ALOG(ZRADMOY)-PCOEFA) /=  0. ) THEN
         !  ZOBS(JILOOP,JJLOOP)=PCOEFB/(ALOG(ZRADMOY)-PCOEFA)


         !WHERE ( XVAR /= XUNDEF .AND. XVAR > 0. .AND. (ALOG(XVAR)-PCOEFA) /=  0.) &
         !XVAR=PCOEFB/(ALOG(XVAR)-PCOEFA)
                 XVAR=XVAR
         ! Retour à l unité initiale
    HUNITS='C'
    print *,'****FROM_COMPUTING_UNITS:Passage Radiance vers Temperature de Brillance avant ecriture ****'
    print *, ' Ce passage est inactif pour l instant'

END SELECT
!
END SUBROUTINE FROM_COMPUTING_UNITS
