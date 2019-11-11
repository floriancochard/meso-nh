!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
MODULE MODE_MODELN_HANDLER
IMPLICIT NONE 

INTEGER, SAVE, PRIVATE     :: ICURRENT_MODEL = -1

CONTAINS 

FUNCTION GET_CURRENT_MODEL_INDEX()
INTEGER :: GET_CURRENT_MODEL_INDEX
!!
GET_CURRENT_MODEL_INDEX = ICURRENT_MODEL
!!
END FUNCTION GET_CURRENT_MODEL_INDEX

SUBROUTINE GOTO_MODEL(KMI, ONOFIELDLIST)
!JUAN
USE MODI_GOTO_MODEL_WRAPPER
!JUAN
INTEGER,           INTENT(IN) :: KMI
LOGICAL, OPTIONAL, INTENT(IN) :: ONOFIELDLIST
!!
IF (ICURRENT_MODEL == -1) THEN
  ICURRENT_MODEL = 1 ! Default model index
  CALL GOTO_MODEL_WRAPPER(ICURRENT_MODEL, KMI, ONOFIELDLIST)
  ICURRENT_MODEL = KMI
ELSE
  IF (ICURRENT_MODEL /= KMI) THEN
!   Switch to model KMI, only if necessary
    CALL GOTO_MODEL_WRAPPER(ICURRENT_MODEL, KMI, ONOFIELDLIST)
    ICURRENT_MODEL = KMI 
  END IF
END IF
!!
END SUBROUTINE GOTO_MODEL

END MODULE MODE_MODELN_HANDLER
