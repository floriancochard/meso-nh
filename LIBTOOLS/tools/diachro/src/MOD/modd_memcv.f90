!     ######spl
      MODULE  MODD_MEMCV
!     ####################
!
!!****  *MODD_MEMCV* - 
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
!!     original        10/11/96
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


! 
! Info. pour superposer en INTERACTIF (LSTI=T)
! 1 symbole (LSYMB=T)
! 1 texte   (LTEXTG=T) sur le graphique (LTEXTIT=T) hors du graphique
! LSYMBTEXTG=T 1 symbole + 1 texte sur le graphique
LOGICAL,SAVE :: LSYMB=.FALSE., LTEXTG=.FALSE., LTEXTIT=.FALSE., &
		LSYMBTEXTG=.FALSE., LSTI=.FALSE.
!
! Info. pour tracer la trace d'une CV (et PH) dans un plan horizontal
LOGICAL,SAVE :: LTRACECV=.FALSE.
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE    :: XTRACECV, XYTRACECV
REAL,SAVE    :: XLWTRACECV=3.
INTEGER,SAVE :: NTRACECV
!
! Memorisation de la directive courante
CHARACTER(LEN=2400) :: CDIRCUR, CDIRPREC   ! LEN=LEN(CAR240) de diaprog.f90
!
! Longueur en fraction axe X de la fleche de l'echelle (<-> 20m/s ou a XVHCPH
! s'il est =/= de 20 = PHA/4 = IPHAS4 dans echelle.f90 .Peut etre module dans
! echelleph.f90)
! dans le cas d'un PH vecteurs (LCH+LCV+LUMVM+LTRACECV)
!
REAL,SAVE  :: XVRLPH=-1.
REAL,SAVE  :: XVHCPH=20.
!
! Logique d'eventuelle elimination de la legende des fleches en CH+CV
! (Defaut : T)
LOGICAL,SAVE :: LEGVECT=.TRUE.
END MODULE MODD_MEMCV
