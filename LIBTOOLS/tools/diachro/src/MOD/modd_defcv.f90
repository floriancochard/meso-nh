!     ######spl
      MODULE  MODD_DEFCV
!     ####################
!
!!****  *MODD_DEFCV* - 
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


! Definition  des limites CV en coord. conformes
LOGICAL,SAVE :: LDEFCV2
REAL,SAVE    :: XIDEBCV, XIFINCV, XJDEBCV, XJFINCV
! Definition  des limites CV en Lat/lon
LOGICAL,SAVE :: LDEFCV2LL
REAL,SAVE    :: XIDEBCVLL, XIFINCVLL, XJDEBCVLL, XJFINCVLL
! Definition  des limites CV en indices de points de grille
LOGICAL,SAVE :: LDEFCV2IND
INTEGER,SAVE :: NIDEBCV, NIFINCV, NJDEBCV, NJFINCV
!
! Logique general pour moi
LOGICAL,SAVE :: LDEFCV2CC
!
! Angle de la coupe en valeur reelle (/axe des X)
REAL,SAVE    :: XANGLECV
!
! PV : localisation en indices de grille, LL et CC (Transmission entre trapro
! et pro1d)
INTEGER, SAVE :: NIPROFV, NJPROFV
REAL,SAVE     :: XIPROFV, XJPROFV
END MODULE MODD_DEFCV
