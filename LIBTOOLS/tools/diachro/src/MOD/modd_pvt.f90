!     ######spl
      MODULE  MODD_PVT
!     #############################
!
!!****  *MODD_PVT* - Contient les parametres de gestion de la couleur
!!      des fleches du vent horizontal  UV (couleur induite par 1 3eme
!!      parametre) dans le seul cas a ce jour (22/3/2000) d'un PV
!!      enregistre dans un fic. diachronique  de type 'CART'
!!      Le gpe contient U, V et d'autres parametres
!!      U et V doivent avoir ete enr. sur la grille 1
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

! Cas _UMVM_ et _PVT_ (L1DT=T)
! Tableau avec les indices de couleur mis a jour ds OPER et utilise ds
! VVUMXY
INTEGER,DIMENSION(:,:), ALLOCATABLE,SAVE  :: NCOL2DUV  
!
! Logique mis a T par pg si _UMVM_ et _PVT_ et 3 processus fournis
LOGICAL,SAVE :: LCOLPVT=.FALSE.

INTEGER,DIMENSION(7),SAVE :: NCOLUVSTD=(/15,4,5,7,3,2,10/)
INTEGER,SAVE                          :: NBCOLUVSTD=7, NBCOLUV
INTEGER,SAVE                          :: NBPARCOLUVSTD=6, NBPARCOLUV

REAL,DIMENSION(6), SAVE  :: XPARCOLUVSTD

! User
INTEGER,DIMENSION(50),SAVE :: NINDCOLUV
REAL,DIMENSION(50), SAVE  :: XPARCOLUV

LOGICAL,SAVE :: LCOLUSERUV=.FALSE.
INTEGER,SAVE :: NISKIPVX=1, NISKIPVY=1

! Septembre 2000  (Pour commodite)
! Ajout pour les couleurs de fleches de imagev et imcouv
INTEGER,SAVE :: NCOLUVG=1, NCOLUV1=1, NCOLUV2=1, NCOLUV3=1, NCOLUV4=1,NCOLUV5=1
!
! Octobre 2000 (Pour Jerome -> coordonnee verticale=pression pour _PVT_)
LOGICAL,SAVE :: LPRESY=.FALSE., LPRESYT=.FALSE.
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XPRESM
REAL,SAVE :: XPMIN=0.,XPMAX=0.,XPINT=0.
!
END MODULE MODD_PVT
