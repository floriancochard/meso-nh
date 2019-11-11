!MNH_LIC Copyright 2003-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_OPEN_FILE_MNH
!     #########################
INTERFACE
      SUBROUTINE OPEN_FILE_MNH(KUNIT,HFILE,HFORM,HACTION,HACCESS,KRECL)
!
INTEGER,           INTENT(OUT):: KUNIT    ! logical unit
CHARACTER(LEN=28), INTENT(IN) :: HFILE    ! file to open
CHARACTER(LEN=11), INTENT(IN) :: HFORM    ! type of file
CHARACTER(LEN=9),  INTENT(IN) :: HACTION  ! action
CHARACTER(LEN=6),  INTENT(IN) :: HACCESS  ! access type
INTEGER,           INTENT(IN) :: KRECL    ! record length
!
END SUBROUTINE OPEN_FILE_MNH
!
END INTERFACE
END MODULE MODI_OPEN_FILE_MNH
!
!     #######################################################
      SUBROUTINE OPEN_FILE_MNH(KUNIT,HFILE,HFORM,HACTION,HACCESS,KRECL)
!     #######################################################
!
!!****  *OPEN_FILE_MNH* - routine to open a file
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
USE MODD_IO_NAM,           ONLY: TFILE
!
USE MODE_FM,               ONLY: IO_FILE_OPEN_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,           INTENT(OUT):: KUNIT    ! logical unit
CHARACTER(LEN=28), INTENT(IN) :: HFILE    ! file to open
CHARACTER(LEN=11), INTENT(IN) :: HFORM    ! type of file
CHARACTER(LEN=9),  INTENT(IN) :: HACTION  ! action
CHARACTER(LEN=6),  INTENT(IN) :: HACCESS  ! access type
INTEGER,           INTENT(IN) :: KRECL    ! record length
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IRESP
!
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','OPEN_FILE_MNH','called for '//TRIM(HFILE))
!
CALL IO_FILE_ADD2LIST(TFILE,TRIM(HFILE),'SURFACE_DATA',HACTION, &
                      HFORM=HFORM,HACCESS=HACCESS,KRECL=KRECL,  &
                      OOLD=.TRUE.) !OOLD=T because the file may already be in list
CALL IO_FILE_OPEN_ll(TFILE)
!
KUNIT = TFILE%NLU
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE OPEN_FILE_MNH
