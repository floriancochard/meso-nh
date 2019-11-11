!     ######spl
      SUBROUTINE LOAD_XPRDAT(KIND,KLOOPT)
!     ################################
!
!!****  *LOAD_FMTAXES* - 
!!
!!    PURPOSE
!!    -------
!       Charger dans XPRDAT les dates modele, exp., segment et courante
!       pour ecriture dans le fichier FICVAL (G.Jaubert JUIN 2001)
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
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
!!      Original       20/06/01
!!      Updated   PM   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO

IMPLICIT NONE
!
!*       0.1   Dummy arguments

INTEGER          :: KIND,KLOOPT
!
!------------------------------------------------------------------------------
IF(.NOT.ALLOCATED(XPRDAT))THEN
  RETURN
ENDIF
! Chargement des dates courante , modele, experience et segment
XPRDAT(1:4,KIND)=XDATIME(13:16,KLOOPT)
XPRDAT(5:8,KIND)=XDATIME(9:12,KLOOPT)
XPRDAT(9:12,KIND)=XDATIME(1:4,KLOOPT)
XPRDAT(13:16,KIND)=XDATIME(5:8,KLOOPT)
RETURN
END SUBROUTINE LOAD_XPRDAT 

