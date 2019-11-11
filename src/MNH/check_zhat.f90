!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#####################
MODULE MODI_CHECK_ZHAT
!#####################
INTERFACE
      SUBROUTINE CHECK_ZHAT(HFMFILE,HDAD_NAME)
!
CHARACTER(LEN=*),    INTENT(IN)    :: HFMFILE   ! name of the Mesonh input file
CHARACTER(LEN=*),    INTENT(INOUT) :: HDAD_NAME ! true name of the Mesonh input file
!
END SUBROUTINE CHECK_ZHAT
END INTERFACE
END MODULE MODI_CHECK_ZHAT
!     ########################################
      SUBROUTINE CHECK_ZHAT(HFMFILE,HDAD_NAME)
!     ########################################
!
!!****  *CHECK_ZHAT* - checks coherence between large scale and fine
!!                     vertical grids for nesting purposes
!! 
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!  1 Resolution ratios during previous spawning are read in FM file.
!!  2 The 2 orographies are averaged on the grid with coarser resolution.
!!  3 If the 2 orographies on coarse grids are identical, then DAD_NAME
!!    is kept; if not, it is initialized to ' ', and nesting wont be
!!    allowed between the output file and its father.
!!
!!    EXTERNAL
!!    --------
!!
!!    function FMREAD
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     : contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1 
!!         XZHAT
!!      Module MODD_DIM1 
!!         NKMAX
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/09/96
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_DIM_n
USE MODD_GRID_n
USE MODD_IO_ll,            ONLY: TFILEDATA
USE MODD_LUNIT,            ONLY: TLUOUT0
USE MODD_PARAMETERS
!
USE MODE_FMREAD
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
CHARACTER(LEN=*),    INTENT(IN)    :: HFMFILE   ! name of the Mesonh input file
CHARACTER(LEN=*),    INTENT(INOUT) :: HDAD_NAME ! true name of the Mesonh input file
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER             :: IKMAX                ! vertical dimension in input file
REAL, DIMENSION(:), ALLOCATABLE :: ZZHAT    ! vertical grid in input file
LOGICAL             :: GSLEVE               ! flag for sleve coordinate
REAL                :: ZLEN1                ! Decay scale for smooth topography
REAL                :: ZLEN2                ! Decay scale for small-scale topography deviation
!
INTEGER             :: IRESP                ! return-code if problems occured
INTEGER             :: ILUOUT0              ! logical unit for file CLUOUT0
LOGICAL             :: GTHINSHELL
TYPE(TFILEDATA),POINTER :: TZFMFILE
!
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!-------------------------------------------------------------------------------
!
!*            1. Reading grid and dimension
!                --------------------------
!
CALL IO_FILE_FIND_BYNAME(TRIM(HFMFILE),TZFMFILE,IRESP)
!
CALL IO_READ_FIELD(TZFMFILE,'KMAX',IKMAX)
ALLOCATE(ZZHAT(IKMAX+2*JPVEXT))
CALL IO_READ_FIELD(TZFMFILE,'ZHAT',ZZHAT)
CALL IO_READ_FIELD(TZFMFILE,'THINSHELL',GTHINSHELL)
IF ( TZFMFILE%NMNHVERSION(1)<4 .OR. (TZFMFILE%NMNHVERSION(1)==4 .AND. TZFMFILE%NMNHVERSION(2)<=6) ) THEN
  GSLEVE = .FALSE.
ELSE
  CALL IO_READ_FIELD(TZFMFILE,'SLEVE',GSLEVE)
ENDIF
!
!*            2. Check dimensions
!                ----------------
!
IF ( IKMAX /= NKMAX ) THEN
  HDAD_NAME=' '
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'Vertical grid has a new number of levels; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
  RETURN
END IF
!
!*            3. Check the vertical grid
!                -----------------------
!
IF ( ANY(ABS(XZHAT(:)-ZZHAT(:))>1.E-10) ) THEN
  HDAD_NAME=' '
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'Vertical grid has been changed; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
END IF
!
!*            4. Check the thinshell approximation
!                ---------------------------------
!
IF ( (GTHINSHELL .OR. LTHINSHELL) .AND. (.NOT. GTHINSHELL .OR. .NOT. LTHINSHELL) ) THEN
  HDAD_NAME=' '
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'thinshell approximation changed; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
END IF
!
!*            5. Check the type of vertical grid
!                -------------------------------
!
IF ( (GSLEVE .OR. LSLEVE) .AND. (.NOT. GSLEVE .OR. .NOT. LSLEVE) ) THEN
  HDAD_NAME=' '
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'type of grid (SLEVE or GAL-CHEN) changed; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
END IF
!
!*            6. Check the SLEVE coordinate parameters
!                -------------------------------------
!
IF ( GSLEVE .AND. LSLEVE ) THEN
  CALL IO_READ_FIELD(TZFMFILE,'LEN1',ZLEN1)
  CALL IO_READ_FIELD(TZFMFILE,'LEN2',ZLEN2)
  IF (ZLEN1 /= XLEN1 .OR. ZLEN2 /= XLEN2) THEN
    HDAD_NAME=' '
    WRITE (ILUOUT0,*) '********************************************************'
    WRITE (ILUOUT0,*) 'Decay scales for SLEVE coordinate changed; no nesting allowed'
    WRITE (ILUOUT0,*) '********************************************************'
  END IF
END IF
!
DEALLOCATE(ZZHAT)
!-------------------------------------------------------------------------------
!
WRITE (ILUOUT0,*) 'Routine CHECK_ZHAT completed'
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CHECK_ZHAT
