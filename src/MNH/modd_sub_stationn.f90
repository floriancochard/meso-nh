!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 12:28:56
!-----------------------------------------------------------------
!     ############################
      MODULE MODD_SUB_STATION_n
!     ############################
!
!!****  *MODD_STATION* - declaration of stations
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the different stations types.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!       
!!    AUTHOR
!!    ------
!!	P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/01/02
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE SUB_STATION_t
!
  LOGICAL :: GSTATFIRSTCALL  = .TRUE.     ! 
!
  INTEGER,DIMENSION(:), POINTER :: II=>NULL()       ! mass station position (x index)
  INTEGER,DIMENSION(:), POINTER :: IJ=>NULL()       ! mass station position (y index)
  INTEGER,DIMENSION(:), POINTER :: IU=>NULL()       ! U flux point station position (x index)
  INTEGER,DIMENSION(:), POINTER :: IV=>NULL()       ! V flux point station position (y index)
!
  REAL, DIMENSION(:), POINTER  :: ZTHIS_PROCS=>NULL()     ! 
!
  REAL,DIMENSION(:), POINTER :: ZXCOEF=>NULL()      ! X direction interpolation coefficient
  REAL,DIMENSION(:), POINTER :: ZUCOEF=>NULL()      ! X direction interpolation coefficient (for U)
  REAL,DIMENSION(:), POINTER :: ZYCOEF=>NULL()      ! Y direction interpolation coefficient
  REAL,DIMENSION(:), POINTER :: ZVCOEF=>NULL()      ! Y direction interpolation coefficient (for V)
!

END TYPE SUB_STATION_t

TYPE(SUB_STATION_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: SUB_STATION_MODEL

LOGICAL, POINTER :: GSTATFIRSTCALL=>NULL()
INTEGER,DIMENSION(:), POINTER :: II=>NULL()
INTEGER,DIMENSION(:), POINTER :: IJ=>NULL()
INTEGER,DIMENSION(:), POINTER :: IU=>NULL()
INTEGER,DIMENSION(:), POINTER :: IV=>NULL()
REAL, DIMENSION(:), POINTER  :: ZTHIS_PROCS=>NULL()
REAL,DIMENSION(:), POINTER :: ZXCOEF=>NULL()
REAL,DIMENSION(:), POINTER :: ZUCOEF=>NULL()
REAL,DIMENSION(:), POINTER :: ZYCOEF=>NULL()
REAL,DIMENSION(:), POINTER :: ZVCOEF=>NULL()

CONTAINS

SUBROUTINE SUB_STATION_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
SUB_STATION_MODEL(KFROM)%II=>II
SUB_STATION_MODEL(KFROM)%IJ=>IJ
SUB_STATION_MODEL(KFROM)%IU=>IU
SUB_STATION_MODEL(KFROM)%IV=>IV
SUB_STATION_MODEL(KFROM)%ZTHIS_PROCS=>ZTHIS_PROCS
SUB_STATION_MODEL(KFROM)%ZXCOEF=>ZXCOEF
SUB_STATION_MODEL(KFROM)%ZUCOEF=>ZUCOEF
SUB_STATION_MODEL(KFROM)%ZYCOEF=>ZYCOEF
SUB_STATION_MODEL(KFROM)%ZVCOEF=>ZVCOEF
!
! Current model is set to model KTO
GSTATFIRSTCALL=>SUB_STATION_MODEL(KTO)%GSTATFIRSTCALL
II=>SUB_STATION_MODEL(KTO)%II
IJ=>SUB_STATION_MODEL(KTO)%IJ
IU=>SUB_STATION_MODEL(KTO)%IU
IV=>SUB_STATION_MODEL(KTO)%IV
ZTHIS_PROCS=>SUB_STATION_MODEL(KTO)%ZTHIS_PROCS
ZXCOEF=>SUB_STATION_MODEL(KTO)%ZXCOEF
ZUCOEF=>SUB_STATION_MODEL(KTO)%ZUCOEF
ZYCOEF=>SUB_STATION_MODEL(KTO)%ZYCOEF
ZVCOEF=>SUB_STATION_MODEL(KTO)%ZVCOEF

END SUBROUTINE SUB_STATION_GOTO_MODEL

END MODULE MODD_SUB_STATION_n
