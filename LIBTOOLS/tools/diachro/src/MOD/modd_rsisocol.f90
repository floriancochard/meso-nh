!     ######spl
      MODULE  MODD_RSISOCOL
!     #############################
!
!!****  *MODD_RSISOCOL* - Contient les parametres de gestion de la couleur
!!      des RS et isocontours dans le cas ou on veut une seule couleur
!!      et en trait plein
!!
!!    PURPOSE
!!    -------
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
!!     original        01/02/96
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE

! Cas _RS_ et _RS1_
! 
LOGICAL, SAVE   :: LCOLRSONE=.FALSE.
INTEGER, SAVE   :: NCOLRSONE=0
LOGICAL, SAVE   :: LCOLRS1ONE=.FALSE.
INTEGER, SAVE   :: NCOLRS1ONE1=0
INTEGER, SAVE   :: NCOLRS1ONE2=0
INTEGER, SAVE   :: NCOLRS1ONE3=0
INTEGER, SAVE   :: NCOLRS1ONE4=0
INTEGER, SAVE   :: NCOLRS1ONE5=0
!
! Pour recuperer les altitudes  des RS sur la grille 1 pour les noter
! sur les profils de vent
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XALTRS
!
! Isocontours
LOGICAL,SAVE :: LCOLISONE=.FALSE.
!
INTEGER,SAVE    :: NCOLISONE1=0           
INTEGER,SAVE    :: NCOLISONE2=0           
INTEGER,SAVE    :: NCOLISONE3=0           
INTEGER,SAVE    :: NCOLISONE4=0           
INTEGER,SAVE    :: NCOLISONE5=0           
!
END MODULE MODD_RSISOCOL
