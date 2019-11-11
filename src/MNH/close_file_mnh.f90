!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_CLOSE_FILE_MNH
!     #############################
INTERFACE
      SUBROUTINE CLOSE_FILE_MNH(HPROGRAM,KUNIT)
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KUNIT    ! logical unit of file
!
END SUBROUTINE CLOSE_FILE_MNH
!
END INTERFACE
END MODULE MODI_CLOSE_FILE_MNH
!
!     #######################################################
      SUBROUTINE CLOSE_FILE_MNH(HPROGRAM,KUNIT)
!     #######################################################
!
!!****  *CLOSE_FILE_MNH* - closes file read by surface in MESOHN
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
!!      Original    01/2003 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,             ONLY: CPROGRAM
USE MODD_IO_ll,            ONLY: TFILEDATA
USE MODD_IO_NAM,           ONLY: TFILE
USE MODD_LUNIT,            ONLY: TLUOUT0
USE MODD_LUNIT_n,          ONLY: TLUOUT
!
USE MODE_FM,               ONLY: IO_FILE_CLOSE_ll
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KUNIT    ! logical unit of file
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILUOUT         ! output listing logical unit
TYPE(TFILEDATA),POINTER :: TZFILE
!-------------------------------------------------------------------------------
!
SELECT CASE(CPROGRAM)
  CASE('REAL  ','IDEAL ','DIAG  ','PGD   ')
    TZFILE => TLUOUT0
    ILUOUT = TLUOUT0%NLU
  CASE('MESONH','SPAWN ')
    TZFILE => TLUOUT
    ILUOUT = TLUOUT%NLU
  CASE DEFAULT
    TZFILE => NULL()
    ILUOUT = -1
END SELECT
!
!-------------------------------------------------------------------------------
!
!* special case: closing of the output listing file
!  ------------------------------------------------
!
IF (ILUOUT==KUNIT) THEN
  CALL PRINT_MSG(NVERB_DEBUG,'IO','CLOSE_FILE_MNH','called for '//TRIM(TZFILE%CNAME))
  CALL IO_FILE_CLOSE_ll(TZFILE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!* closes the namelist
!  -------------------
!
IF (.NOT.ASSOCIATED(TFILE)) CALL PRINT_MSG(NVERB_FATAL,'IO','CLOSE_FILE_MNH','TFILE not associated')
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','CLOSE_FILE_MNH','called for '//TRIM(TFILE%CNAME))
!
IF (TFILE%NLU==KUNIT) THEN
  CALL IO_FILE_CLOSE_ll(TFILE)
  TFILE => NULL()
ELSE
  WRITE(ILUOUT,*) 'Error for closing a file: '
  WRITE(ILUOUT,*) 'logical unit ',KUNIT,' does not correspond to file', TFILE%CNAME
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'IO','CLOSE_FILE_MNH','')
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLOSE_FILE_MNH
