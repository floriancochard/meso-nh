!     ######spl
      MODULE  MODD_MEMGRIUV
!     #############################
!
!!****  *MODD_MEMGRIUV* - 
!!
!!    PURPOSE
!!    -------
!       Memorisation du numero de grille de U et V ds read_uvw
!       pour test dans precou pour faire ou non l'interpolation
!       sur la grille de masse
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        28/11/01
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE

INTEGER,SAVE :: NGRIU=2 , NGRIV=3

END MODULE MODD_MEMGRIUV
