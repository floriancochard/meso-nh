!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:48
!-----------------------------------------------------------------
MODULE YOERRTAM

#include "tsmbkind.h"

USE OPARRRTM


IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOERRTAM* - RRTM definition of atmospheric profiles
!     includes all commons with arrays depending on vertical coordinate
!     ------------------------------------------------------------------

REAL_B   , ALLOCATABLE, DIMENSION(:,:,:) :: TAU

REAL_B   , ALLOCATABLE, DIMENSION(:,:) :: TAUAERL

REAL_B   , ALLOCATABLE, DIMENSION(:) :: FAC00 , FAC01 , FAC10 , FAC11
INTEGER_M, ALLOCATABLE, DIMENSION(:) :: JP    , JT    , JT1
REAL_B    :: ONEMINUS
REAL_B   , ALLOCATABLE, DIMENSION(:) :: COLH2O, COLCO2, COLO3 &
                                    &, COLN2O, COLCH4, COLO2 &
                                    &, CO2MULT 
INTEGER_M :: LAYTROP, LAYSWTCH, LAYLOW
REAL_B   , ALLOCATABLE, DIMENSION(:) :: PAVEL , TAVEL
REAL_B   , ALLOCATABLE, DIMENSION(:) :: PZ    , TZ
REAL_B    :: TBOUND
INTEGER_M :: NLAYERS
REAL_B   , ALLOCATABLE, DIMENSION(:) :: SELFFAC, SELFFRAC
INTEGER_M, ALLOCATABLE, DIMENSION(:) :: INDSELF
REAL_B   , ALLOCATABLE, DIMENSION(:,:,:) :: PFRAC
REAL_B   , ALLOCATABLE, DIMENSION(:) :: SEMISS
REAL_B    :: SEMISLW
INTEGER_M :: IREFLECT
INTEGER_M :: NUMANGS, IOUT, ISTART, IEND
REAL_B   , ALLOCATABLE, DIMENSION(:) :: COLDRY , WBRODL
REAL_B   , ALLOCATABLE, DIMENSION(:) :: CLDFRAC
REAL_B   , ALLOCATABLE, DIMENSION(:,:) :: TAUCLDU, TAUCLDD
INTEGER_M :: NMOL, IXSECT, NXMOL
REAL_B   , ALLOCATABLE, DIMENSION(:,:) :: WKL
INTEGER_M, ALLOCATABLE, DIMENSION(:) :: IXINDX
REAL_B   , ALLOCATABLE, DIMENSION(:,:) :: WX
REAL_B   , ALLOCATABLE, DIMENSION(:) :: FORFAC
INTEGER_M, ALLOCATABLE, DIMENSION(:) :: INDLAY
INTEGER_M, ALLOCATABLE, DIMENSION(:) :: INDLEV

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
!     -----------------------------------------------------------------
END MODULE YOERRTAM
