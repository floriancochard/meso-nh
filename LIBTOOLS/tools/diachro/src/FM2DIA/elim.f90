!     ######spl
      SUBROUTINE ELIM(HRECFM)
!     #######################
!
!!****  *ELIM* - Mise a 0 des numeros d'enregistrements lus correspondant
!                aux parametres "intouchables"
!!
!!    PURPOSE
!!    -------
!       On met arbitrairement a 0 les numeros d'enregistrements lus correspon-
!       -dant aux parametres "intouchables" pour les eliminer du traitement
!       realise dans la routine WRITE_OTHERFIELDS 
!
!!**  METHOD
!!    ------
!       On met a une valeur nulle les numeros d'enregistrements correspon-
!       -dant aux parametres "intouchables" ecrits dans le fichier
!      diachronique une fois pour toutes avec la routine 
!      WRITE_LFIFM1_FORDIACHRO_CV pour ne pas les prendre en compte dans
!      la routine WRITE_OTHERFIELDS
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHORS
!!    -------
!!    J. Duron      *Lab. Aerologie* 
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/01/96 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE  MODD_DIMGRID_FORDIACHRO        
!
!*       0.1   Dummy arguments
!
CHARACTER(LEN=*)  :: HRECFM 
!
!*       0.2   Local variables declarations
!
INTEGER           :: J
!
!----------------------------------------------------------------------------
!
DO J=1,NNB
  IF(HRECFM == CRECFM2T(J,1))NNUMT(J,1)=0
ENDDO
!
!----------------------------------------------------------------------------

RETURN

END SUBROUTINE ELIM
