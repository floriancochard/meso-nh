!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_MNHGET_DESFM_n
!     #########################
INTERFACE
      SUBROUTINE MNHGET_DESFM_n(HACTION,KLUDES)
!
CHARACTER(LEN=5), INTENT(IN)  :: HACTION ! 'READ ', 'WRITE'
INTEGER,          INTENT(OUT) :: KLUDES  ! logical unit of .des file
!
END SUBROUTINE MNHGET_DESFM_n
!
END INTERFACE
END MODULE MODI_MNHGET_DESFM_n
!
!     #######################################################
      SUBROUTINE MNHGET_DESFM_n(HACTION,KLUDES)
!     #######################################################
!
!!****  *MNHGET_DESFM* - routine to open .des file
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
!!      S.Malardel   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2003
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CONF,           ONLY : CPROGRAM
USE MODD_LUNIT_n,        ONLY : TINIFILE
USE MODD_LUNIT,          ONLY : TPGDFILE,TOUTDATAFILE
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=5), INTENT(IN)  :: HACTION ! 'READ ', 'WRITE'
INTEGER,          INTENT(OUT) :: KLUDES  ! logical unit of .des file
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!NONE
!-------------------------------------------------------------------------------
!
!*       1.    Return logical unit of .des files
!
KLUDES=0
!
IF (HACTION=='READ ') THEN
  SELECT CASE(CPROGRAM)
    CASE('MESONH','DIAG  ')
      KLUDES = TINIFILE%TDESFILE%NLU
    CASE('REAL  ')
      KLUDES = TPGDFILE%TDESFILE%NLU
    CASE('IDEAL ')
      KLUDES = 0
  END SELECT
ELSE IF (HACTION=='WRITE') THEN
  IF (CPROGRAM == 'PGD   ' .OR. CPROGRAM =='NESPGD' .OR. &
      CPROGRAM == 'ZOOMPG' .OR. CPROGRAM =='DIAG  '      ) THEN
    KLUDES = 0
  ELSE
    KLUDES = TOUTDATAFILE%TDESFILE%NLU
  END IF
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHGET_DESFM_n
