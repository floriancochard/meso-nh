!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_IO_SURF_MNH
!     ##################
!
!!****  *MODD_IO_SURF_MNH - 
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	S.Malardel   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!	M.Faivre 2014
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!
!*       0.   DECLARATIONS
!
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_PARAMETERS, ONLY: JPMODELMAX

IMPLICIT NONE

INTEGER                              :: NHALO = 0

TYPE IO_SURF_MNH_t
  TYPE(TFILEDATA),POINTER        :: TPINFILE => NULL() ! Input FM-file
  CHARACTER(LEN=28)              :: COUTFILE    ! Name of the output FM-file
  TYPE(TFILEDATA),POINTER        :: TOUT => NULL() ! Output_listing file
  CHARACTER(LEN=6)               :: CMASK
  INTEGER, DIMENSION(:), POINTER :: NMASK=>NULL()     ! 1D mask to read only interesting surface
  !                                           ! points on current processor
  INTEGER, DIMENSION(:), POINTER :: NMASK_ALL=>NULL() ! 1D mask to read all surface points all processors
  !
  CHARACTER(LEN=5)               :: CACTION = '     '! action being done ('READ ','WRITE')
  !
  ! number of points in each direction on current processor
  INTEGER                              :: NIU,NJU
  ! indices of physical points in each direction on current processor
  INTEGER                              :: NIB,NJB,NIE,NJE
  ! number of points in each direction on all processors
  INTEGER                              :: NIU_ALL,NJU_ALL
  ! indices of physical points in each direction on all processors
  INTEGER                              :: NIB_ALL,NJB_ALL,NIE_ALL,NJE_ALL
END type IO_SURF_MNH_t
!
TYPE(IO_SURF_MNH_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: IO_SURF_MNH_MODEL
!
!!!!!!!!!!!!!!!!!!!! LOCAL VARIABLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE(TFILEDATA)       ,POINTER :: TPINFILE => NULL()  ! Input FM-file
CHARACTER(LEN=28)     ,POINTER :: COUTFILE =>NULL()   ! Name of the output FM-file
TYPE(TFILEDATA)       ,POINTER :: TOUT => NULL()      ! Output_listing file
CHARACTER(LEN=6)      ,POINTER :: CMASK =>NULL()
INTEGER, DIMENSION(:), POINTER :: NMASK=>NULL()     ! 1D mask to read only interesting surface
!                                           ! points on current processor
INTEGER, DIMENSION(:), POINTER :: NMASK_ALL=>NULL() ! 1D mask to read all surface points all processors
!
CHARACTER(LEN=5)      ,POINTER :: CACTION => NULL() ! action being done ('READ ','WRITE')
!
! number of points in each direction on current processor
INTEGER             , POINTER  :: NIU=>NULL(),NJU=>NULL()
! indices of physical points in each direction on current processor
INTEGER             , POINTER  :: NIB=>NULL(),NJB=>NULL(),NIE=>NULL(),NJE=>NULL()
! number of points in each direction on all processors
INTEGER             , POINTER  :: NIU_ALL=>NULL(),NJU_ALL=>NULL()
! indices of physical points in each direction on all processors
INTEGER             , POINTER  :: NIB_ALL=>NULL(),NJB_ALL=>NULL(),NIE_ALL=>NULL(),NJE_ALL=>NULL()
!
CONTAINS

SUBROUTINE IO_SURF_MNH_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
! save curretnt state for allocated arrays
IO_SURF_MNH_MODEL(KFROM)%NMASK=>NMASK
IO_SURF_MNH_MODEL(KFROM)%NMASK_ALL=>NMASK_ALL

! current model is set for model KTO 
TPINFILE=>IO_SURF_MNH_MODEL(KTO)%TPINFILE
COUTFILE=>IO_SURF_MNH_MODEL(KTO)%COUTFILE
TOUT=>IO_SURF_MNH_MODEL(KTO)%TOUT
CMASK=>IO_SURF_MNH_MODEL(KTO)%CMASK
NMASK=>IO_SURF_MNH_MODEL(KTO)%NMASK
NMASK_ALL=>IO_SURF_MNH_MODEL(KTO)%NMASK_ALL
CACTION=>IO_SURF_MNH_MODEL(KTO)%CACTION
NIU=>IO_SURF_MNH_MODEL(KTO)%NIU
NJU=>IO_SURF_MNH_MODEL(KTO)%NJU
NIB=>IO_SURF_MNH_MODEL(KTO)%NIB
NJB=>IO_SURF_MNH_MODEL(KTO)%NJB
NIE=>IO_SURF_MNH_MODEL(KTO)%NIE
NJE=>IO_SURF_MNH_MODEL(KTO)%NJE
NIU_ALL=>IO_SURF_MNH_MODEL(KTO)%NIU_ALL
NJU_ALL=>IO_SURF_MNH_MODEL(KTO)%NJU_ALL
NIB_ALL=>IO_SURF_MNH_MODEL(KTO)%NIB_ALL
NJB_ALL=>IO_SURF_MNH_MODEL(KTO)%NJB_ALL
NIE_ALL=>IO_SURF_MNH_MODEL(KTO)%NIE_ALL
NJE_ALL=>IO_SURF_MNH_MODEL(KTO)%NJE_ALL
END SUBROUTINE IO_SURF_MNH_GOTO_MODEL

END MODULE MODD_IO_SURF_MNH
