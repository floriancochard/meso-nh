!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################################
      SUBROUTINE MNHOPEN_WRITE_COVER_TEX(KTEX)
!     ##################################
!
!!****  *MNHOPEN_WRITE_COVER_TEX* - opens cover listing file (in MESONH universe)
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
USE MODE_FM,               ONLY: IO_FILE_OPEN_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, INTENT(OUT) :: KTEX ! logical unit of Tex file
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=*),PARAMETER :: YTEX = 'class_cover_data.tex' ! name of tex file
!
TYPE(TFILEDATA),POINTER :: TZFILE
!-------------------------------------------------------------------------------
!
!*       5.     Prints of cover parameters in a tex file
!               ----------------------------------------
!
TZFILE => NULL()
!
IF (TRIM(CPROGRAM)=='PGD') THEN
  CALL IO_FILE_ADD2LIST(TZFILE,YTEX,'TXT','WRITE')
  CALL IO_FILE_OPEN_ll(TZFILE,HPOSITION='REWIND')
  KTEX = TZFILE%NLU
ELSE
  KTEX=0
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHOPEN_WRITE_COVER_TEX
