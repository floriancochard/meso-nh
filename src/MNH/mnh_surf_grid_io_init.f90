!MNH_LIC Copyright 2015-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#######################
MODULE MODI_MNH_SURF_GRID_IO_INIT
  !#######################
  !
  INTERFACE
    !     ###############################
          SUBROUTINE MNH_SURF_GRID_IO_INIT(KIMAX,KJMAX)
    !     ###############################
    !!
    !!    PURPOSE
    !!    -------
    !!
    !!    Initializes parallel data structures for I/O in surfex (for grid reading)
    !!
    !!    METHOD
    !!    ------
    !!
    !!    EXTERNAL
    !!    --------
    !!
    !!
    !!    IMPLICIT ARGUMENTS
    !!    ------------------
    !!
    !!
    !!    REFERENCE
    !!    ---------
    !!
    !!    AUTHOR
    !!    ------
    !!
    !!    M.Moge                   CNRS - LA
    !!
    !!    MODIFICATION
    !!    ------------
    !!
    !!    Original      19/03/2015
    !!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
    !----------------------------------------------------------------------------
    !
    !*    0.     DECLARATION
    !            -----------
    !
    USE MODE_ll
    USE MODE_FM
    USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT, JPMODELMAX
    USE MODD_CONF,       ONLY : CPROGRAM, L1D, L2D, LPACK
    !
    USE MODE_SPLITTINGZ_ll
    !
    USE MODI_GET_SURF_GRID_DIM_N
    USE MODI_GET_LUOUT
    !
    IMPLICIT NONE
    !
    !*    0.1    Declaration of dummy arguments
    !            ------------------------------
    !
    INTEGER,               INTENT(IN)    :: KIMAX ! number of points in X direction
    INTEGER,               INTENT(IN)    :: KJMAX ! number of points in Y direction
          END SUBROUTINE MNH_SURF_GRID_IO_INIT
  !
  END INTERFACE
END MODULE MODI_MNH_SURF_GRID_IO_INIT
!     ###############################
      SUBROUTINE MNH_SURF_GRID_IO_INIT(KIMAX,KJMAX)
!     ###############################
!!
!!    PURPOSE
!!    -------
!!
!!    Initializes parallel data structures for I/O in surfex (for grid reading)
!!
!!    METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    M.Moge                   CNRS - LA
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original      19/03/2015
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODE_ll
USE MODE_FM
USE MODE_IO_ll
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT, JPMODELMAX
USE MODD_CONF,       ONLY : CPROGRAM, L1D, L2D, LPACK
!
!JUANZ
USE MODE_SPLITTINGZ_ll
!JUANZ
!
USE MODI_GET_SURF_GRID_DIM_N
USE MODI_GET_LUOUT
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
INTEGER,               INTENT(IN)    :: KIMAX ! number of points in X direction
INTEGER,               INTENT(IN)    :: KJMAX ! number of points in Y direction
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IINFO_ll ! return code of // routines
INTEGER :: ILUOUT   ! output listing logical unit
!
!------------------------------------------------------------------------------
!
IF (CPROGRAM=='IDEAL ' .OR. CPROGRAM=='SPAWN ' .OR. CPROGRAM=='REAL  ') RETURN
!
L1D=(KIMAX==1).AND.(KJMAX==1)
L2D=(KIMAX/=1).AND.(KJMAX==1)
LPACK=L1D.OR.L2D
CALL SET_FMPACK_ll(L1D,L2D,LPACK)
CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
CALL SET_DIM_ll(KIMAX, KJMAX, 1)
CALL SET_LBX_ll('OPEN',1)
CALL SET_LBY_ll('OPEN', 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(KIMAX+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(KJMAX+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
CALL INI_PARAZ_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNH_SURF_GRID_IO_INIT
