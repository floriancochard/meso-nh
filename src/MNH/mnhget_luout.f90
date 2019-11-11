!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_MNHGET_LUOUT
!     #############################
INTERFACE
      SUBROUTINE MNHGET_LUOUT(HPROGRAM,KLUOUT)
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling GROUND
INTEGER,           INTENT(OUT) :: KLUOUT   ! Logical unit of output listing
!
END SUBROUTINE MNHGET_LUOUT
!
END INTERFACE
END MODULE MODI_MNHGET_LUOUT
!
!     #######################################################
      SUBROUTINE MNHGET_LUOUT(HPROGRAM,KLUOUT)
!     #######################################################
!
!!****  *MNHGET_LUOUT* - get output listing logical unit
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,    ONLY: CPROGRAM
USE MODE_ll
USE MODD_LUNIT,   ONLY: TLUOUT0
USE MODD_LUNIT_n, ONLY: LUNIT_MODEL,TLUOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling GROUND
INTEGER,           INTENT(OUT) :: KLUOUT   ! Logical unit of output listing
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IMI ! model index
!
!-------------------------------------------------------------------------------
!
SELECT CASE (CPROGRAM)
  CASE ('REAL  ','PGD   ','NESPGD')
    KLUOUT = TLUOUT0%NLU
  CASE ('IDEAL ')
    KLUOUT = TLUOUT%NLU
  CASE ('MESONH','DIAG  ','SPAWN ')
    CALL GET_MODEL_NUMBER_ll(IMI)
    KLUOUT = LUNIT_MODEL(IMI)%TLUOUT%NLU
  CASE DEFAULT
    KLUOUT = TLUOUT0%NLU
END SELECT
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHGET_LUOUT
