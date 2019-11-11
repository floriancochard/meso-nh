!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_MNHOPEN_NAMELIST
!     #############################
INTERFACE
      SUBROUTINE MNHOPEN_NAMELIST(HPROGRAM,KLUNAM,HFILE)
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(OUT) :: KLUNAM   ! logical unit of namelist
CHARACTER(LEN=28), INTENT(IN)  :: HFILE ! ASCII file to open
!
END SUBROUTINE MNHOPEN_NAMELIST
!
END INTERFACE
END MODULE MODI_MNHOPEN_NAMELIST
!
!     #######################################################
      SUBROUTINE MNHOPEN_NAMELIST(HPROGRAM,KLUNAM,HFILE)
!     #######################################################
!
!!****  *MNHOPEN_NAMELIST* - opens namelists files for surface (MESONH universe)
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
USE MODD_CONF,             ONLY: CPROGRAM
USE MODD_IO_NAM,           ONLY: TNAM
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
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(OUT) :: KLUNAM   ! logical unit of namelist
CHARACTER(LEN=28), INTENT(IN)  :: HFILE ! ASCII file to open
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears 
                                    ! at the open of the file in LFI  routines 
INTEGER           :: IMI            ! model index
!
CHARACTER(LEN=28) :: YNAM ! name of namelist
!-------------------------------------------------------------------------------
!
!* reading of namelist
!  -------------------
!
IF (LEN_TRIM(HFILE)>0) THEN
  YNAM = HFILE
ELSE
 SELECT CASE(CPROGRAM)
  CASE('PGD   ')
    YNAM='PRE_PGD1.nam '
  CASE('REAL  ')
    YNAM='PRE_REAL1.nam'
  CASE('IDEAL ')
    YNAM='PRE_IDEA1.nam'
  CASE('MESONH')
    CALL GET_MODEL_NUMBER_ll(IMI)
    WRITE(YNAM,FMT='(A5,I1,A22)') 'EXSEG',IMI,'.nam                  '
  CASE('DIAG  ')
    YNAM='DIAG1.nam    '
  CASE('SPAWN ')
    YNAM='SPAWN1.nam   '
  CASE('NESPGD')
    YNAM='PRE_NEST_PGD1.nam'
  CASE('ZOOMPG')
    YNAM='PRE_ZOOM1.nam'
  CASE('SPEC ')
    YNAM='SPEC1.nam'
  CASE DEFAULT
    CALL PRINT_MSG(NVERB_FATAL,'IO','MNHOPEN_NAMELIST','CPROGRAM '//TRIM(CPROGRAM)//' not allowed')
 END SELECT
END IF
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','MNHOPEN_NAMELIST','called for '//TRIM(YNAM))
!
CALL IO_FILE_ADD2LIST(TNAM,TRIM(YNAM),'NML','READ',OOLD=.TRUE.) !OOLD=T because the file may already be in list
CALL IO_FILE_OPEN_ll(TNAM)
!
KLUNAM = TNAM%NLU
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHOPEN_NAMELIST
