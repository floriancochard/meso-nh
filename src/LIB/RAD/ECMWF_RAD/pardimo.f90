!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARDIMO


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!****-------------------------------------------------------------------
!****  CD PARDIMO : PARAMETERS RELATED TO OBSERVATIONS
!****  ----------
!****  Auteurs    : J.PAILLEUX, D.VASILJEVIC, P.CAILLE    90/05-91/09
!****-------------------------------------------------------------------
!*  OBSERVATIONS ARRAYS DIMENSIONS
!*  ------------------------------
!*    JPNOTP  : MAXIMUM NUMBER OF OBSERVATIONS TYPES
!*    JPXTIM  : MAXIMUM NUMBER OF TIME SLOTS
!*
!*  OBSERVATIONS PROCESSING
!*  -----------------------
!*    JPNBTM  : MAXIMUM NUMBER OF MESSAGES TYPES
!*    JPNUMV  : NUMBER OF VARIABLES IN THE ARRAY NVNUMB
!*    JPXSOBT : MAXIMUM NUMBER OF SUB-OBSERVATIONS TYPES
!*    JPXVAR  : MAXIMUM NUMBER OF VARIABLES
!*    JPXAREA : MAXIMUM NUMBER OF AREAS
!*    JPXDEP  : MAXIMUM NUMBER OF ADDITIONAL OBS. DEP.
!*    JPXUPD  : MAXIMUM NUMBER OF UPDATES ALLOWED FOR CMA
!*    JPXDEL  : LEN. OF THE OBS. SEP. LIST FOR SAVING
!*
!****-------------------------------------------------------------------
INTEGER_M, PARAMETER :: JPNOTP=10
INTEGER_M, PARAMETER :: JPXTIM=25

INTEGER_M, PARAMETER :: JPNBTM=36
INTEGER_M, PARAMETER :: JPNUMV=65
INTEGER_M, PARAMETER :: JPXSOBT=3
INTEGER_M, PARAMETER :: JPXVAR=25
INTEGER_M, PARAMETER :: JPXAREA=20
INTEGER_M, PARAMETER :: JPXDEP=100
INTEGER_M, PARAMETER :: JPXDEL=639
INTEGER_M, PARAMETER :: JPXUPD=100

!****-------------------------------------------------------------------
END MODULE PARDIMO
