!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_allvar.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_ALLVAR
!     ###################
!
!!****  *MODD_ALLVAR* - Declaration des tableaux de travail pour les 
!                       variables autres que prognostiques
!                       et des types de variables permettant la memorisation
!                       du nom de ces variables, du parametre NGRID, des unites
!!
!!    PURPOSE
!!    -------
!       Declare des tableaux de travail pour des variables 3D, 2D, 1D,
!     scalaires ou vectorielles ne figurant pas parmi les champs de base
!     du modele
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!       
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/06/94                      
!!      Updated   PM     /11/94  
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_TYPE_ALLVAR
!
IMPLICIT NONE
!
INTEGER,SAVE  :: NVAR3D, NVAR2D

LOGICAL :: LSCAL1D, LSCAL2D

REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XWORK3D

REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XWORKX3D
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XWORKY3D
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XWORKZ3D

REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XWORK2D

REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XWORKX2D
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XWORKY2D
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XWORKZ2D

REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XWORK1D

TYPE (X_Y_Z_)     :: XT1
TYPE (X_Y_)       :: XT2
TYPE (VX_VY_VZ_)  :: XT3
TYPE (VX_VY_)     :: XT4
TYPE (Z_)         :: XT5
END MODULE MODD_ALLVAR
