!     ######spl
      MODULE  MODI_ALLOC_FORDIACHRO
!     #############################
!
INTERFACE
!
SUBROUTINE ALLOC_FORDIACHRO(KI,KJ,KK,KT,KN,KP,KOP,KNTRAJT,KKTRAJX,  &
 KKTRAJY,KKTRAJZ,KTTRAJX,KTTRAJY,KTTRAJZ,KNTRAJX,KNTRAJY,KNTRAJZ,KIMASK, &
       KJMASK,KKMASK,KTMASK,KNMASK,KPMASK)
INTEGER :: KI,KJ,KK,KT,KN,KP,KOP
INTEGER,OPTIONAL :: KNTRAJT,KKTRAJX,KKTRAJY,KKTRAJZ,KTTRAJX, &
		    KTTRAJY,KTTRAJZ,KNTRAJX,KNTRAJY,KNTRAJZ,KIMASK, &
                    KJMASK,KKMASK,KTMASK,KNMASK,KPMASK
END SUBROUTINE ALLOC_FORDIACHRO
!
END INTERFACE
!
END MODULE MODI_ALLOC_FORDIACHRO
!     ######spl
      SUBROUTINE ALLOC_FORDIACHRO(KI,KJ,KK,KT,KN,KP,KOP,KNTRAJT,KKTRAJX,  &
      KKTRAJY,KKTRAJZ,KTTRAJX,KTTRAJY,KTTRAJZ,KNTRAJX,KNTRAJY,KNTRAJZ,KIMASK, &
      KJMASK,KKMASK,KTMASK,KNMASK,KPMASK)
!     #########################################################################
!
!!****  *ALLOC_FORDIACHRO* - Allocation de tableaux dont les dimensions
!                            sont fournies en arguments  de la routine
!       (VALABLE UNIQUEMENT DANS LE CADRE DU TRAITEMENT D'1 FICHIER
!        DIACHRONIQUE : lecture ou/et ecriture)
!!
!!    PURPOSE
!!    -------
!       En fonction d'un code operation transmis dans l'argument KOP
!       alloue ou desalloue des tableaux
!
!!**  METHOD
!!    ------
!!     
!!     KOP=1
!      Alloue des tableaux (en utilisant les 6 1ers arguments qui sont des
!      dimensions fournies par l'utilisateur) destines a etre charges par
!      l'utilisateur et ecrits dans 1 fichier diachronique.
!      Le nombre, le nom et le profil de ces tableaux est dependant du
!      type d'informations a ecrire (CTYPE du MODULE : MODD_TYPE_AND_LH)
!     
!      KOP=2
!      Alloue des tableaux (dont les dimensions ont ete lues dans un
!      enregistrement d'1 fichier diachronique et transmis en arguments)
!      destines a lire les valeurs du champ correspondant au groupe
!      demande dans le meme fichier diachronique
!
!      KOP=3
!      Desalloue les tableaux alloues avec KOP=1 ou 2
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/01/96
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ALLOC_FORDIACHRO
USE MODD_TYPE_AND_LH

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER :: KI, KJ, KK, KT, KN,KP, KOP
INTEGER,OPTIONAL :: KNTRAJT,KKTRAJX,KKTRAJY,KKTRAJZ,KTTRAJX, &
		    KTTRAJY,KTTRAJZ,KNTRAJX,KNTRAJY,KNTRAJZ,KIMASK, &
		    KJMASK,KKMASK,KTMASK,KNMASK,KPMASK
!
!*       0.1   Local variables
!              ---------------

!
!------------------------------------------------------------------------------
!
IF (KOP == 1)THEN

  ALLOCATE(XDATIME(16,KT))
  ALLOCATE(NGRIDIA(KP))
  SELECT CASE (CTYPE)
    CASE ('CART','SPXY')
      ALLOCATE(XVAR(KI,KJ,KK,KT,KN,KP))
      ALLOCATE(XTRAJT(KT,KN))
    CASE ('MASK','SSOL')
      ALLOCATE(XVAR(1,1,KK,KT,KN,KP))
      ALLOCATE(XTRAJT(KT,1))
    CASE ('DRST','RAPL')
      ALLOCATE(XVAR(1,1,KK,KT,KN,KP))
      ALLOCATE(XTRAJT(KT,KN))
    CASE ('RSPL')
      ALLOCATE(XVAR(1,1,1,KT,KN,KP))
      ALLOCATE(XTRAJT(KT,KN))
  END SELECT

  ALLOCATE(CTITRE(KP),CUNITE(KP),CCOMMENT(KP))

  IF (CTYPE == 'SSOL')THEN
    ALLOCATE(XTRAJX(1,1,KN),XTRAJY(1,1,KN),XTRAJZ(KK,1,KN))
  ENDIF
  IF (CTYPE == 'DRST')THEN
    ALLOCATE(XTRAJX(1,KT,KN),XTRAJY(1,KT,KN),XTRAJZ(KK,KT,KN))
  ENDIF
  IF (CTYPE == 'RSPL')THEN
    ALLOCATE(XTRAJX(1,KT,KN),XTRAJY(1,KT,KN),XTRAJZ(1,KT,KN))
  ENDIF
  IF (CTYPE == 'RAPL')THEN
    ALLOCATE(XTRAJX(KK,KT,KN),XTRAJY(KK,KT,KN),XTRAJZ(KK,KT,KN))
  ENDIF

  IF (CTYPE == 'MASK')THEN
    ALLOCATE(XMASK(KI,KJ,1,KT,KN,1))
  ENDIF

ELSE IF(KOP == 2)THEN

  ALLOCATE(XDATIME(16,KT))
  ALLOCATE(XVAR(KI,KJ,KK,KT,KN,KP))
  ALLOCATE(XTRAJT(KT,KNTRAJT))
  ALLOCATE(CTITRE(KP),CUNITE(KP),CCOMMENT(KP))
  CTITRE(:)(1:LEN(CTITRE))=' '
  CUNITE(:)(1:LEN(CUNITE))=' '
  CCOMMENT(:)(1:LEN(CCOMMENT))=' '
  ALLOCATE(NGRIDIA(KP))
  IF(KKTRAJX /= 0 .AND. KTTRAJX /= 0 .AND. KNTRAJX /=0 )THEN
    ALLOCATE(XTRAJX(KKTRAJX,KTTRAJX,KNTRAJX))
  ENDIF
  IF(KKTRAJY /= 0 .AND. KTTRAJY /= 0 .AND. KNTRAJY /=0 )THEN
    ALLOCATE(XTRAJY(KKTRAJY,KTTRAJY,KNTRAJY))
  ENDIF
  IF(KKTRAJZ /= 0 .AND. KTTRAJZ /= 0 .AND. KNTRAJZ /=0 )THEN
    ALLOCATE(XTRAJZ(KKTRAJZ,KTTRAJZ,KNTRAJZ))
  ENDIF
  IF(KIMASK /= 0 .AND. KJMASK /= 0 .AND. KKMASK/= 0 .AND. &
     KTMASK /= 0 .AND. KNMASK /= 0 .AND. KPMASK/= 0 )THEN
     ALLOCATE(XMASK(KIMASK,KJMASK,KKMASK,KTMASK,KNMASK,KPMASK))
  ENDIF

ELSE

  IF(ALLOCATED(XDATIME))DEALLOCATE(XDATIME)
  IF(ALLOCATED(XVAR))DEALLOCATE(XVAR)
  IF(ALLOCATED(XTRAJT))DEALLOCATE(XTRAJT)
  IF(ALLOCATED(CTITRE))DEALLOCATE(CTITRE)
  IF(ALLOCATED(CUNITE))DEALLOCATE(CUNITE)
  IF(ALLOCATED(CCOMMENT))DEALLOCATE(CCOMMENT)
! DEALLOCATE(XVAR,XTRAJT,CTITRE,CUNITE,CCOMMENT)

  IF (CTYPE == 'SSOL' .OR. &
      CTYPE == 'DRST' .OR. &
      CTYPE == 'RSPL' .OR. &
      CTYPE == 'RAPL')THEN
    IF(ALLOCATED(XTRAJX)) DEALLOCATE(XTRAJX)
    IF(ALLOCATED(XTRAJY)) DEALLOCATE(XTRAJY)
    IF(ALLOCATED(XTRAJZ)) DEALLOCATE(XTRAJZ)
  ENDIF

  IF (CTYPE == 'MASK')THEN
    IF(ALLOCATED(XMASK)) DEALLOCATE(XMASK)
  ENDIF

  IF(ALLOCATED(NGRIDIA))THEN
    DEALLOCATE(NGRIDIA)
  ENDIF

ENDIF

!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE ALLOC_FORDIACHRO
