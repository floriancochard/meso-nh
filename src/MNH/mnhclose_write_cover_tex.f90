!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################################
      SUBROUTINE MNHCLOSE_WRITE_COVER_TEX
!     ##################################
!
!!****  *MNHCLOSE_WRITE_COVER_TEX* - closes cover listing file
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
!
USE MODE_FM,               ONLY: IO_FILE_CLOSE_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
!
USE MODI_TRANSFER_FILE
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=*),PARAMETER :: YTEX = 'class_cover_data.tex' ! name of tex file
!
INTEGER                 :: IRESP
TYPE(TFILEDATA),POINTER :: TZFILE
!-------------------------------------------------------------------------------
!
!*       5.     Prints of cover parameters in a tex file
!               ----------------------------------------
!
TZFILE => NULL()
!
IF (TRIM(CPROGRAM)=='PGD') THEN
  CALL IO_FILE_FIND_BYNAME(YTEX,TZFILE,IRESP)
  CALL IO_FILE_CLOSE_ll(TZFILE)
  CALL TRANSFER_FILE('fujitransfer.x','NIL',YTEX)
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHCLOSE_WRITE_COVER_TEX
