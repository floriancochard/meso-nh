!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 surfex 2006/05/23 16:02:23
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_MNHGET_SIZE_FULL_n
!     #########################
INTERFACE
      SUBROUTINE MNHGET_SIZE_FULL_n(HPROGRAM,KDIM_FULL,KSIZE_FULL)
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER         ,  INTENT(IN)  :: KDIM_FULL  ! total number of points
INTEGER         ,  INTENT(OUT) :: KSIZE_FULL ! total number of points on this proc
!
END SUBROUTINE MNHGET_SIZE_FULL_n
!
END INTERFACE
END MODULE MODI_MNHGET_SIZE_FULL_n
!
!     #######################################################
      SUBROUTINE MNHGET_SIZE_FULL_n(HPROGRAM,KDIM_FULL,KSIZE_FULL)
!     #######################################################
!
!!****  *MNHGET_SIZE_FULL_n* - routine to get number of points on this processor
!!                             (MESONH universe)
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
!!      Original    09/2003 
!!                  02/2015 (M.Moge) case('PGD') to compute KSIZE_FULL
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODE_ll
!
USE MODD_DIM_n,      ONLY : NIMAX, NJMAX
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODD_CONF,       ONLY : CPROGRAM
!
USE MODD_IO_SURF_MNH, ONLY : NHALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER         ,  INTENT(IN)  :: KDIM_FULL  ! total number of points
INTEGER         ,  INTENT(OUT) :: KSIZE_FULL ! total number of points on this proc
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IIB, IIE, IJB, IJE
!-------------------------------------------------------------------------------
!
SELECT CASE(CPROGRAM)
  CASE ('NESPGD')
    IIB = JPHEXT + 1
    IIE = JPHEXT + NIMAX
    IJB = JPHEXT + 1
    IJE = JPHEXT + NJMAX
  CASE DEFAULT
    CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
END SELECT
! 
SELECT CASE(CPROGRAM)
  CASE ('PGD')
    KSIZE_FULL = (IIE-IIB+1)*(IJE-IJB+1)
  CASE DEFAULT
    KSIZE_FULL = (IIE-IIB+1+2*NHALO)*(IJE-IJB+1+2*NHALO)
END SELECT
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHGET_SIZE_FULL_n
