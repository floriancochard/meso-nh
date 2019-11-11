!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARKIND1
!
!     *** Define usual kinds for strong typing ***
!     J.Escobar : 9/06/2015, for I*8 compilation force JPIM to default size
!     J.-P. Chaboureau: 14/10/2016, adding logical kind JPLM for RTTOV
!     J.Escobar 30/03/2017  : Management of compilation of ECMWF_RAD in REAL*8 with MNH_REAL=R4
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIT = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: JPIS = SELECTED_INT_KIND(4)
INTEGER :: JINT_DEF
INTEGER, PARAMETER :: JPIM = KIND(JINT_DEF) ! SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)
!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JPRT = SELECTED_REAL_KIND(2,1)
INTEGER, PARAMETER :: JPRS = SELECTED_REAL_KIND(4,2)
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)
REAL               :: REAL_DEF_JPRB
INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300) !  KIND(REAL_DEF_JPRB) 
INTEGER, PARAMETER :: JPRB_DEF = KIND(REAL_DEF_JPRB) 
!
!     Logical Kinds
!     -------------
INTEGER, PARAMETER :: JPLM = KIND(.TRUE.)               !Standard logical type
!
END MODULE PARKIND1
