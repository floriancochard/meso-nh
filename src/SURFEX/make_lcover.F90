!SFX_LIC Copyright 1997-2018 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.

!     #########
      SUBROUTINE MAKE_LCOVER(OCOVER)
!     ##############################################################
!
!!**** *PGD_COVER* monitor for averaging and interpolations of cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NPROC, NCOMM
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef SFX_MNH
USE MODD_IO_ll, ONLY : ISIOP, ISP, ISNPROC
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
#endif
!
IMPLICIT NONE
!
#if defined(SFX_MPI) || defined(SFX_MNH)
INCLUDE "mpif.h"
#endif
!
!*    0.1    Declaration of arguments
!            ------------------------
!
LOGICAL, DIMENSION(:), INTENT(INOUT) :: OCOVER
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: INFOMPI, JPROC, JCOVER
!
INTEGER :: IRANK_SAVE, IPROC_SAVE, IPIO_SAVE, ICOMM_SAVE
!
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GCOVER_ALL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('MAKE_LCOVER',0,ZHOOK_HANDLE)
!
#ifdef SFX_MNH
IRANK_SAVE = NRANK
IPROC_SAVE = NPROC
IPIO_SAVE = NPIO
ICOMM_SAVE = NCOMM
!
! on met les infos de mésonh
NRANK = ISP-1
NPROC = ISNPROC
NPIO = ISIOP-1
NCOMM = NMNH_COMM_WORLD
#endif
!
ALLOCATE(GCOVER_ALL(SIZE(OCOVER),0:NPROC-1))
!
!
IF (NPROC>1) THEN
#if defined(SFX_MPI) || defined(SFX_MNH)
  CALL MPI_ALLGATHER(OCOVER,SIZE(OCOVER),MPI_LOGICAL,GCOVER_ALL,SIZE(OCOVER),&
                  MPI_LOGICAL,NCOMM,INFOMPI)
#endif
ELSE
  GCOVER_ALL(:,0) = OCOVER(:)
ENDIF
!
!
OCOVER(:) = .FALSE.
DO JPROC = 0,NPROC-1
  DO JCOVER=1,SIZE(OCOVER)
    IF (GCOVER_ALL(JCOVER,JPROC)) OCOVER(JCOVER) = .TRUE.
  ENDDO
ENDDO
!
DEALLOCATE(GCOVER_ALL)
!
!
IF (NPROC>1) THEN
#if defined(SFX_MPI) || defined(SFX_MNH)
  CALL MPI_BCAST(OCOVER,SIZE(OCOVER),MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
#endif
ENDIF
!
#ifdef SFX_MNH
NRANK = IRANK_SAVE
NPROC = IPROC_SAVE
NPIO = IPIO_SAVE
NCOMM = ICOMM_SAVE
#endif
!
IF (LHOOK) CALL DR_HOOK('MAKE_LCOVER',1,ZHOOK_HANDLE)
!
END SUBROUTINE MAKE_LCOVER
