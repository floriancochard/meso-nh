!     ######spl
      MODULE  MODD_DIACHRO
!     ####################
!
!!****  *MODD_DIACHRO* - 
!!
!!    PURPOSE
!!    -------
!
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
!!     original  JD    08/01/96
!!     updated   PM   
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


! Nom fichier diachronique
!
CHARACTER(LEN=28),SAVE :: CFILEDIA                
!
! Listing associe au traitement diachronique
!
CHARACTER(LEN=16),SAVE :: CLUOUTDIA='OUT_DIA'               
!
! Numero logique du listing et parametres d'ouverture du fichier
!
INTEGER,SAVE           :: NLUOUTDIA, NNPRARDIA, NFTYPEDIA=2, NVERBDIA,  &
                          NNINARDIA, NRESPDIA

CHARACTER(LEN=28),SAVE :: CMY_NAME_DIA            
CHARACTER(LEN=28),SAVE :: CDAD_NAME_DIA

!
END MODULE MODD_DIACHRO
