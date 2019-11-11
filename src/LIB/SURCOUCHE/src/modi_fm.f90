!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------

MODULE MODI_FM_ll
!
INTERFACE 
!
SUBROUTINE SET_FMPACK_ll(O1D,O2D,OPACK)
LOGICAL, INTENT(IN) :: O1D,O2D,OPACK
END SUBROUTINE SET_FMPACK_ll
!
SUBROUTINE IO_FILE_OPEN_ll(TPFILE,KRESP,OPARALLELIO,HPOSITION,HSTATUS,HPROGRAM_ORIG)
USE MODD_IO_ll, ONLY: TFILEDATA
TYPE(TFILEDATA),POINTER,INTENT(INOUT)         :: TPFILE ! File structure
INTEGER,                INTENT(OUT), OPTIONAL :: KRESP  ! Return code
LOGICAL,                INTENT(IN),  OPTIONAL :: OPARALLELIO
CHARACTER(LEN=*),       INTENT(IN),  OPTIONAL :: HPOSITION
CHARACTER(LEN=*),       INTENT(IN),  OPTIONAL :: HSTATUS
CHARACTER(LEN=*),       INTENT(IN),  OPTIONAL :: HPROGRAM_ORIG !To emulate a file coming from this program
END SUBROUTINE IO_FILE_OPEN_ll
!
SUBROUTINE IO_FILE_CLOSE_ll(TPFILE,KRESP,OPARALLELIO,HPROGRAM_ORIG)
USE MODD_IO_ll, ONLY: TFILEDATA
TYPE(TFILEDATA),  INTENT(INOUT)         :: TPFILE ! File structure
INTEGER,          INTENT(OUT), OPTIONAL :: KRESP  ! Return code
LOGICAL,          INTENT(IN),  OPTIONAL :: OPARALLELIO
CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: HPROGRAM_ORIG !To emulate a file coming from this program
END SUBROUTINE IO_FILE_CLOSE_ll
!
END INTERFACE
END MODULE MODI_FM_ll
