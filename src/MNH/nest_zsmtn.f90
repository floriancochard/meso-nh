!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_NEST_ZSMT_n
!     #######################
!
INTERFACE
!
      SUBROUTINE NEST_ZSMT_n(YFIELD)
!
CHARACTER(LEN=6),   INTENT(IN)  :: YFIELD ! name of the field to nest
!
END SUBROUTINE NEST_ZSMT_n
!
END INTERFACE
!
END MODULE MODI_NEST_ZSMT_n
!
!
!
!     #############################
      SUBROUTINE NEST_ZSMT_n(YFIELD)
!     #############################
!
!!****  *NEST_ZSMT_n* - make smooth orography field coherent between model $n and its sons
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        12/01/05
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CONF, ONLY: NMODEL, CPROGRAM
USE MODD_NESTING, ONLY: NDAD
USE MODD_GRID_n, ONLY: XZSMT, XZS
!
USE MODI_FILL_ZSMTn
USE MODE_MODELN_HANDLER
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: YFIELD ! name of the field to nest
!
!
!*       0.2   declarations of local variables
!
INTEGER :: JMI ! model index loop
INTEGER :: IMI ! current model index
REAL, DIMENSION(:,:),  POINTER  :: DPTR_XZSMT ! dummy pointer to correct ifort bug
!
!-------------------------------------------------------------------------------
!
!*       1.    Modify son field
!              ----------------
!
IMI = GET_CURRENT_MODEL_INDEX()
DO JMI=1,NMODEL
  IF (IMI /= NDAD(JMI)) CYCLE
  DPTR_XZSMT=>XZSMT
  CALL FILL_ZSMT_n(YFIELD,DPTR_XZSMT,JMI)
END DO
CALL MPPDB_CHECK2D(XZS,"nest_zsmt_n:XZS",PRECISION)
CALL MPPDB_CHECK2D(XZSMT,"nest_zsmt_n:XZSMT",PRECISION)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE NEST_ZSMT_n
