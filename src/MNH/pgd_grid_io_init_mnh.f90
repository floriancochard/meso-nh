!MNH_LIC Copyright 2004-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#######################
MODULE MODI_PGD_GRID_IO_INIT_MNH
  !#######################
  !
  INTERFACE
    !     ###############################
#ifdef MNH_PARALLEL
          SUBROUTINE PGD_GRID_IO_INIT_MNH(UG,KGRID_PAR,PGRID_PAR,HGRID,ORECT,KIMAX,KJMAX,KDXRATIO,KDYRATIO)
#else
      SUBROUTINE PGD_GRID_IO_INIT_MNH(UG)
#endif
    !     ###############################
    !!
    !!    PURPOSE
    !!    -------
    !!
    !!    Initializes parallel routines for further I/O
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
    !!    V. Masson                   Meteo-France
    !!
    !!    MODIFICATION
    !!    ------------
    !!
    !!    Original      01/2004
    !!    10/10/2011  J.Escobar call INI_PARAZ_ll
    !!    2014        M.Faivre
    !!    07/2015     M.Moge when initializing a child model from a father model (with PREP_PGD), 
    !!                we need to initialize the parallel data structures using a modified version
    !!                of INI_PARAZ_ll/INI_CHILD : INI_PARAZ_CHILD_ll
    !!                In this case, when entering PGD_GRID_IO_INIT_MNH we have only one model : the father
    !!                When exiting, we have only one model : the child
    !!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
    !----------------------------------------------------------------------------
    !
    !*    0.     DECLARATION
    !            -----------
    !
    USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
    !
    USE MODE_ll
    USE MODE_FM
    USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT, JPMODELMAX
    USE MODD_CONF,       ONLY : CPROGRAM, L1D, L2D, LPACK
    !
    USE MODE_IO_ll
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
#ifdef MNH_PARALLEL
    TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
    INTEGER,                         INTENT(IN)    :: KGRID_PAR ! size of PGRID_PAR
    REAL,    DIMENSION(KGRID_PAR),   INTENT(IN)    :: PGRID_PAR ! grid parameters
    CHARACTER(LEN=10),     INTENT(IN), OPTIONAL    :: HGRID
    LOGICAL,               INTENT(IN), OPTIONAL    :: ORECT
    ! if KIMAX,KJMAX,KDXRATIO,KDYRATIO present, this means we are in PREP_PGD, and we only initialise the child model, 
    ! using a father model read from a file and previously initialized with INI_PARAZ_ll
    INTEGER,               INTENT(IN), OPTIONAL    :: KIMAX
    INTEGER,               INTENT(IN), OPTIONAL    :: KJMAX
    INTEGER,               INTENT(IN), OPTIONAL    :: KDXRATIO ! ratio in X direction
    INTEGER,               INTENT(IN), OPTIONAL    :: KDYRATIO ! ratio in Y direction
#endif
          END SUBROUTINE PGD_GRID_IO_INIT_MNH
  !
  END INTERFACE
END MODULE MODI_PGD_GRID_IO_INIT_MNH
!     ###############################
#ifdef MNH_PARALLEL
      SUBROUTINE PGD_GRID_IO_INIT_MNH(UG,KGRID_PAR,PGRID_PAR,HGRID,ORECT,KIMAX,KJMAX,KDXRATIO,KDYRATIO)
#else
      SUBROUTINE PGD_GRID_IO_INIT_MNH(UG)
#endif
!     ###############################
!!
!!    PURPOSE
!!    -------
!!
!!    Initializes parallel routines for further I/O
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
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original      01/2004
!!    10/10/2011  J.Escobar call INI_PARAZ_ll
!!    2014        M.Faivre
!!  06/2016     (G.Delautier) phasage surfex 8
!!  01/2018      (G.Delautier) SURFEX 8.1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
!
USE MODE_ll
USE MODE_FM
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT, JPMODELMAX
USE MODD_CONF,       ONLY : CPROGRAM, L1D, L2D, LPACK
USE MODD_DIM_n,      ONLY : NIMAX_ll, NJMAX_ll, NKMAX
!
!JUANZ
USE MODE_SPLITTINGZ_ll
!JUANZ
!
USE MODI_GET_SURF_GRID_DIM_N
USE MODI_GET_LUOUT
!
USE MODD_MNH_SURFEX_n
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
#ifdef MNH_PARALLEL
INTEGER,                         INTENT(IN)    :: KGRID_PAR ! size of PGRID_PAR
REAL,    DIMENSION(KGRID_PAR),   INTENT(IN)    :: PGRID_PAR ! grid parameters
CHARACTER(LEN=10),     INTENT(IN), OPTIONAL    :: HGRID
LOGICAL,               INTENT(IN), OPTIONAL    :: ORECT
! if KIMAX,KJMAX,KDXRATIO,KDYRATIO present, this means we are in PREP_PGD, and we only initialise the child model, 
! using a father model read from a file and previously initialized with INI_PARAZ_ll
INTEGER,               INTENT(IN), OPTIONAL    :: KIMAX
INTEGER,               INTENT(IN), OPTIONAL    :: KJMAX
INTEGER,               INTENT(IN), OPTIONAL    :: KDXRATIO ! ratio in X direction
INTEGER,               INTENT(IN), OPTIONAL    :: KDYRATIO ! ratio in Y direction
#endif
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IINFO_ll ! return code of // routines
INTEGER :: IIMAX    ! number of points in X direction
INTEGER :: IJMAX    ! number of points in Y direction
INTEGER :: IDXRATIO ! ratio in X direction
INTEGER :: IDYRATIO ! ratio in Y direction
INTEGER :: ILUOUT   ! output listing logical unit
!
LOGICAL :: GRECT           ! true when grid is rectangular
CHARACTER(LEN=10) :: YGRID ! grid type
!
!------------------------------------------------------------------------------
!
IF (CPROGRAM=='IDEAL ' .OR. CPROGRAM=='SPAWN ') RETURN
!
!
#ifdef MNH_PARALLEL
IF ( PRESENT(KIMAX) .AND. PRESENT(KJMAX) .AND. PRESENT(HGRID) .AND. PRESENT(ORECT) &
  .AND. PRESENT(KDXRATIO) .AND. PRESENT(KDYRATIO) ) THEN
  YGRID = HGRID
  GRECT = ORECT
  IIMAX = KIMAX
  IJMAX = KJMAX
  IDXRATIO = KDXRATIO
  IDYRATIO = KDYRATIO
ELSE
  CALL GET_SURF_GRID_DIM_n(UG,YGRID,GRECT,IIMAX,IJMAX,KGRID_PAR,PGRID_PAR)
  IDXRATIO = 1
  IDYRATIO = 1
ENDIF
#else
  CALL GET_SURF_GRID_DIM_n(UG,YGRID,GRECT,IIMAX,IJMAX)
  IDXRATIO = 1
  IDYRATIO = 1
#endif
!
!
IF (YGRID/='CONF PROJ ' .AND. YGRID/='CARTESIAN') THEN
  CALL GET_LUOUT('MESONH',ILUOUT)
  WRITE(ILUOUT,*) "Error, grid type GRID=",YGRID, &
                  " is not supported by MESONH"
END IF
!------------------------------------------------------------------------------
!
!
L1D=(IIMAX==1).AND.(IJMAX==1)
L2D=(IIMAX/=1).AND.(IJMAX==1)
LPACK=L1D.OR.L2D
CALL SET_FMPACK_ll(L1D,L2D,LPACK)
CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
NIMAX_ll = IIMAX
NJMAX_ll = IJMAX
NKMAX    = 1
CALL SET_DIM_ll(IIMAX, IJMAX, 1)
CALL SET_LBX_ll('OPEN',1)
CALL SET_LBY_ll('OPEN', 1)
CALL SET_XRATIO_ll(IDXRATIO, 1)
CALL SET_YRATIO_ll(IDYRATIO, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(IIMAX+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(IJMAX+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
!JUANZ CALL INI_PARA_ll(IINFO_ll)
! for PREP_PGD, when constructing a son grid from a father grid,
! we DON'T want to call INI_PARAZ_ll for the child domain if it has already been called on the father domain :
! INI_PARAZ_ll would split the global son grid without taking into account the RATIO, so it will SPLIT in the middle
! of the cells of the father.
! To avoid this, we call a modified INI_PARAZ_CHILD_ll, that will split the father domain and the use the ratio to 
! get the son splitting.

#ifdef MNH_PARALLEL
IF ( PRESENT(KIMAX) .AND. PRESENT(KJMAX) .AND. PRESENT(HGRID) .AND. PRESENT(ORECT) &
  .AND. PRESENT(KDXRATIO) .AND. PRESENT(KDYRATIO) ) THEN
  CALL INI_PARAZ_CHILD_ll(IINFO_ll)
  CALL SET_XRATIO_ll(1, 1)  ! il faut faire ça dans le cas PREP_PGD sur le modele fils car dans ce cas on ne 
  CALL SET_YRATIO_ll(1, 1)  ! voit en fait plus qu'un seul modele, le modele pere n'existe plus vraiment dans la suite
                            ! donc le ratio n'a plus de sens, et doit etre a 1
ELSE
  CALL INI_PARAZ_ll(IINFO_ll)
ENDIF
#else
CALL INI_PARAZ_ll(IINFO_ll)
#endif
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_GRID_IO_INIT_MNH
