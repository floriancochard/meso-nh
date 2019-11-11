!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:47
!-----------------------------------------------------------------
MODULE YOEAER


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOEAER* - RADIATIVE CHARACTERISTICS OF THE AEROSOLS
!     -----------------------------------------------------------------


!*       1.    SHORTWAVE COEFFICIENTS
!              ----------------------

REAL_B , DIMENSION(4,5)   :: RTAUA=RESHAPE((/&
 &.730719_JPRB,.730719_JPRB,.730719_JPRB,.730719_JPRB,&
 &.912819_JPRB,.912819_JPRB,.912819_JPRB,.912819_JPRB,&
 &.725059_JPRB,.725059_JPRB,.725059_JPRB,.725059_JPRB,&
 &.745405_JPRB,.745405_JPRB,.745405_JPRB,.745405_JPRB,&
 &.682188_JPRB,.682188_JPRB,.682188_JPRB,.682188_JPRB /)&
 &,SHAPE=(/4,5/))
REAL_B , DIMENSION(4,5)   :: RPIZA=RESHAPE((/&
 &.872212_JPRB,.872212_JPRB,.872212_JPRB,.872212_JPRB,&
 &.982545_JPRB,.982545_JPRB,.982545_JPRB,.982545_JPRB,&
 &.623143_JPRB,.623143_JPRB,.623143_JPRB,.623143_JPRB,&
 &.944887_JPRB,.944887_JPRB,.944887_JPRB,.944887_JPRB,&
 &.997975_JPRB,.997975_JPRB,.997975_JPRB,.997975_JPRB /)&
 &,SHAPE=(/4,5/))
REAL_B , DIMENSION(4,5)   :: RCGA=RESHAPE((/&
 &.647596_JPRB,.647596_JPRB,.647596_JPRB,.647596_JPRB,&
 &.739002_JPRB,.739002_JPRB,.739002_JPRB,.739002_JPRB,&
 &.580845_JPRB,.580845_JPRB,.580845_JPRB,.580845_JPRB,&
 &.662657_JPRB,.662657_JPRB,.662657_JPRB,.662657_JPRB,&
 &.624246_JPRB,.624246_JPRB,.624246_JPRB,.624246_JPRB /)&
 &,SHAPE=(/4,5/))

!*       2.    LONGWAVE COEFFICIENTS
!              ---------------------


REAL_B , DIMENSION(5,5)   :: RAER=RESHAPE((/&
 &.038520_JPRB, .037196_JPRB, .040532_JPRB, .054934_JPRB, .038520_JPRB ,&
 &.12613_JPRB , .18313_JPRB , .10357_JPRB , .064106_JPRB, .126130_JPRB ,&
 &.012579_JPRB, .013649_JPRB, .018652_JPRB, .025181_JPRB, .012579_JPRB ,&
 &.011890_JPRB, .016142_JPRB, .021105_JPRB, .028908_JPRB, .011890_JPRB ,&
 &.013792_JPRB, .026810_JPRB, .052203_JPRB, .066338_JPRB, .013792_JPRB /&
       &),SHAPE=(/5,5/))


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : -------
!  RTAUA :  REAL     S.W. NORMALIZED OPTICAL THICKNESS AT 0.55 MICRON
!  RPIZA :  REAL     S.W. SINGLE SCATTERING ALBEDO
!  RCGA  :  REAL     S.W. ASSYMETRY FACTOR
!  RAER  :  REAL     L.W. ABSORPTION COEFFICIENTS
!     -----------------------------------------------------------------
END MODULE YOEAER
