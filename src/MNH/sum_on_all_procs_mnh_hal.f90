!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #########
      SUBROUTINE SUM_ON_ALL_PROCS_MNH_HAL(KSIZE,KIN,KOUT)
!     #######################################################
!
!
!!****  *SUM_ON_ALL_PROCS* - sums the values of the integers provided on each processor
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson    *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2011 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_IO_SURF_MNH, ONLY : NHALO
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                   INTENT(IN) :: KSIZE  ! sim of integer array
INTEGER, DIMENSION(KSIZE), INTENT(IN) :: KIN    ! array of integer to sum
INTEGER,                   INTENT(OUT):: KOUT   ! sum on all processors 
!                                               ! (excluding halos)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL    :: ZIN
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: NIMAX, NJMAX
INTEGER :: JI, JJ
INTEGER :: IINDEX
!
INTEGER :: IRESP ! return code
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
NIMAX=IIE-IIB+1
NJMAX=IJE-IJB+1
!
ZIN = 0.
DO JJ=1,NJMAX+2*NHALO
  DO JI=1,NIMAX+2*NHALO
    IINDEX = JI + (JJ-1) * (NIMAX+2*NHALO)
    ZIN = ZIN + FLOAT(KIN(IINDEX))
  END DO
END DO
!
CALL REDUCESUM_ll(ZIN,IRESP)
KOUT = NINT(ZIN)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SUM_ON_ALL_PROCS_MNH_HAL
