!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_MNHCLOSE_NAMELIST
!     #############################
INTERFACE
      SUBROUTINE MNHCLOSE_NAMELIST(HPROGRAM,KLUNAM)
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KLUNAM   ! logical unit of namelist
!
END SUBROUTINE MNHCLOSE_NAMELIST
!
END INTERFACE
END MODULE MODI_MNHCLOSE_NAMELIST
!
!     #######################################################
      SUBROUTINE MNHCLOSE_NAMELIST(HPROGRAM,KLUNAM)
!     #######################################################
!
!!****  *MNHCLOSE_NAMELIST* - closes namelists read by surface in MESOHN
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
USE MODD_IO_NAM,           ONLY: TNAM
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
INTEGER,           INTENT(IN)  :: KLUNAM   ! logical unit of namelist
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears 
                                    ! at the open of the file in LFI  routines 
!
INTEGER           :: IMI            ! model index
INTEGER           :: ILUOUT         ! output listing logical unit
CHARACTER(LEN=16) :: YLUOUT         ! output listing file name
!-------------------------------------------------------------------------------
!
IF (.NOT.ASSOCIATED(TNAM)) CALL PRINT_MSG(NVERB_FATAL,'IO','CLOSE_FILE_MNH','TNAM not associated')
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','MNHCLOSE_NAMELIST','called for '//TRIM(TNAM%CNAME))
!
!-------------------------------------------------------------------------------
!
!* closes the namelist
!  -------------------
!
IF (TNAM%NLU==KLUNAM) THEN
  CALL IO_FILE_CLOSE_ll(TNAM)
  TNAM => NULL()
ELSE
  SELECT CASE(CPROGRAM)
    CASE('REAL  ','IDEAL ','DIAG  ')
      ILUOUT = TLUOUT0%NLU
    CASE('MESONH','SPAWN ')
      ILUOUT = TLUOUT%NLU
    CASE DEFAULT
      CALL PRINT_MSG(NVERB_FATAL,'GEN','MNHCLOSE_NAMELIST','CPROGRAM '//TRIM(CPROGRAM)//' not allowed')
  END SELECT
  !
  WRITE(ILUOUT,*) 'Error for closing a namelist file: '
  WRITE(ILUOUT,*) 'logical unit ',KLUNAM,' does not correspond to namelist file', TNAM%CNAME
  !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','MNHCLOSE_NAMELIST','')
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHCLOSE_NAMELIST
