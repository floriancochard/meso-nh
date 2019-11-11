!     ############################################################
      MODULE MODI_TO_COMPUTING_UNITS
!     ############################################################
!
INTERFACE
      SUBROUTINE TO_COMPUTING_UNITS(HCHAMP,HUNITS,PVALOBS)
!
CHARACTER(LEN=*) , intent(in) :: HCHAMP                    ! Nom du champ 
CHARACTER(LEN=*) , intent(inout) :: HUNITS                 ! Unité
REAL , intent(inout) ,OPTIONAL:: PVALOBS           ! cas de traitement 1 valeur
!
END SUBROUTINE TO_COMPUTING_UNITS
END INTERFACE
END MODULE MODI_TO_COMPUTING_UNITS
!
!------------------------------------------------------------------------------
!
!     ################
      SUBROUTINE TO_COMPUTING_UNITS(HCHAMP,HUNITS,PVALOBS)
!     ################
!
!!****  *TO_COMPUTING_UNITS* - 
!! 
!!
!!    PURPOSE
!!    -------
!  Passage vers une unité adaptee au calcul
!  appel a From_Computing_Units(YCHAMP,CUNIT) pour revenir a l unite initiale
!
!!**  METHOD
!  par defaut traite le tableau XVAR passe en module
!ou PVALOBS passe en argument
!
!  Changement du nom d unite pour diagnostiquer le traitement inverse
! dans From_Computing_Units (routine symetrique)
!   mettre a jour suivant les variables Mesonh qui necessitent ce passage
! AU 01/2005 : les reflectivités radar exprimees en dBz
!              les températures de brillance
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

USE MODI_LOW2UP
IMPLICIT NONE
!
!*       0.1   Arguments d'appel
!
CHARACTER(LEN=*) , intent(in) :: HCHAMP                    ! Nom du champ 
CHARACTER(LEN=*) , intent(inout) :: HUNITS                 ! Unité
REAL , intent(inout), OPTIONAL:: PVALOBS           ! cas de traitement 1 valeur
!
!*       0.2 variables locales
!
INTEGER :: ILOOP
CHARACTER (LEN=13) :: YNAME
CHARACTER (LEN=10) :: YUNIT
! provisoire pour passer la compil
REAL :: PCOEFA,PCOEFB
!
!-------------------------------------------------------------------------------
!
!print *,'entree TO_COMPUTING_UNITS ',TRIM(HCHAMP),' ',TRIM(HUNITS)
!
! passage en majuscules
YNAME=HCHAMP
CALL LOW2UP(YNAME)
YUNIT=HUNITS
CALL LOW2UP(YUNIT)
!
! Critere= nom de variable
IF (INDEX(HCHAMP,'_IRBT')/=0 .OR. INDEX(HCHAMP,'_WVBT')/=0) THEN
         ! Prevoir la routine inverse a MAKE_RADSAT
         ! Mesonh
         ! passage rad -> temp brillance pour le satellite KGEO
         ! call MAKE_RADSAT(KYEARF, KMONTHF, KDAYF, PSECF, &
         !                  KGEO, KLON, PRADB, PRADF)
         ! Viviane
         !ZOBS est en radiance, je la transforme en tempe de brillance
         ! IF (ZRADMOY > 0. .AND. (ALOG(ZRADMOY)-PCOEFA) /=  0. ) THEN
         !  ZOBS(JILOOP,JJLOOP)=PCOEFB/(ALOG(ZRADMOY)-PCOEFA)
         ! Viviane
         ! transformation des tempe de brillance en radiance
         ! IF ( TAB_OBS(IOBS)%PTROBS%XVALOBS /= ZUNDEF .AND. &
         ! TAB_OBS(IOBS)%PTROBS%XVALOBS /= 0.) THEN
         ! TAB_OBS(IOBS)%PTROBS%XVALOBS = EXP(PCOEFA+PCOEFB/&
         !                    TAB_OBS(IOBS)%PTROBS%XVALOBS)
    IF (PRESENT (PVALOBS)) THEN
    !    IF (PVALOBS/= XUNDEF .AND. PVALOBS /= 0.) PVALOBS=EXP(PCOEFA+PCOEFB/PVALOBS) 
         PVALOBS=PVALOBS
    ELSE
    !    WHERE (XVAR /= XUNDEF .AND. XVAR /= 0.) XVAR=EXP(PCOEFA+PCOEFB/XVAR)
         XVAR=XVAR
    ENDIF
         !
         ! Pour indiquer le travail inverse dans From_Computing_Units
         HUNITS='W_to_C'
         print *,'****TO_COMPUTING_UNITS: Passage Temperature de Brillance vers Radiance avant calcul ****'
         print *,' Ce passage est inactif pour l instant'
ENDIF
!
! Critere = unite
SELECT CASE (YUNIT)
 CASE ('DBZ','dBz','dBZ','ZE_LISTOBS')
         ! Reflectivites radar
   IF (PRESENT (PVALOBS)) THEN     
     IF (PVALOBS /= XUNDEF ) PVALOBS=10.0**(PVALOBS/10.0)
     ! Pour indiquer le travail inverse dans From_Computing_Units
     HUNITS='Ze_listOBS'
   ELSE        
     WHERE (XVAR /= XUNDEF ) XVAR=10.0**(XVAR/10.0)
     ! Pour indiquer le travail inverse dans From_Computing_Units
     HUNITS='Ze_to_DBZ'
  ENDIF
  IF ( YUNIT /= 'ZE_LISTOBS' ) THEN
    ! print pour la premiere Obs traitee seulement
    print *,'****TO_COMPUTING_UNITS: Passage DBZ a Ze avant calcul ****'
  ENDIF
END SELECT
!
END SUBROUTINE TO_COMPUTING_UNITS
