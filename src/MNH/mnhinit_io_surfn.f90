!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_MNHINIT_IO_SURF_n
!     #########################
INTERFACE
      SUBROUTINE MNHINIT_IO_SURF_n(HPROGRAM,HMASK,HACTION)
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    ! type of mask
CHARACTER(LEN=5),  INTENT(IN)  :: HACTION  ! action performed ('READ ','WRITE')
!
END SUBROUTINE MNHINIT_IO_SURF_n
!
END INTERFACE
END MODULE MODI_MNHINIT_IO_SURF_n
!
!     #######################################################
      SUBROUTINE MNHINIT_IO_SURF_n(HPROGRAM,HMASK,HACTION)
!     #######################################################
!
!!****  *MNHINIT_IO_SURF_n* - routine to open surface IO files (MESONH universe)
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------surfex/mnhinit_io_surfn.f90
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
!!	S.Malardel   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2003 
!!  06/2016     (G.Delautier) phasage surfex 8
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,             ONLY: CPROGRAM
USE MODD_DIM_n,            ONLY: NIMAX_ll, NJMAX_ll
USE MODD_IO_SURF_MNH,      ONLY: TOUT, TPINFILE, COUTFILE,                             &
                                 NMASK, CMASK, NIU, NJU, NIB, NJB, NIE, NJE, CACTION,  &
                                 NMASK_ALL, NIU_ALL, NJU_ALL, NIB_ALL, NJB_ALL,        &
                                 NIE_ALL, NJE_ALL, NHALO
USE MODD_LUNIT,            ONLY: CLUOUT0, TPGDFILE, TLUOUT0, TOUTDATAFILE
USE MODD_LUNIT_n,          ONLY: CMASK_SURFEX, TINIFILE, TINIFILEPGD, TLUOUT
USE MODD_MNH_SURFEX_n,     ONLY: YSURF_CUR
USE MODD_PARAMETERS,       ONLY: JPHEXT
!
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODI_GET_1D_MASK
USE MODI_GET_SURF_MASK_N
USE MODI_GET_TYPE_DIM_N
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
CHARACTER(LEN=6),  INTENT(IN)  :: HMASK
CHARACTER(LEN=5),  INTENT(IN)  :: HACTION  ! action performed ('READ ','WRITE')
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILU            ! 1D physical dimension of entire surface on all processors
INTEGER           :: ILU_ALL        ! 1D physical dimension of entire surface on all processors
INTEGER           :: ILM            ! 1D physical dimension of any surface type
REAL, DIMENSION(:),   ALLOCATABLE :: ZFULL  ! total cover
!-------------------------------------------------------------------------------
!
CACTION = HACTION
!
!*       1.    initialization of output listing name
!
SELECT CASE(CPROGRAM)
  CASE('MESONH','SPAWN ')
    TOUT => TLUOUT
  CASE DEFAULT
    TOUT => TLUOUT0
END SELECT
!
!
!*       2.    initialization of surface file
!
IF (HACTION=='READ ') THEN
  SELECT CASE(CPROGRAM)
    CASE('MESONH','DIAG  ')
      IF(CMASK_SURFEX=="PGD") THEN
        TPINFILE=>TINIFILEPGD
      ELSE
        TPINFILE=>TINIFILE
      ENDIF
    CASE('REAL  ','IDEAL ','NESPGD','SPAWN ','ZOOMPG')
      TPINFILE => TPGDFILE
  END SELECT
ELSE IF (HACTION=='WRITE') THEN
  IF (ASSOCIATED(TOUTDATAFILE)) THEN
    COUTFILE = TOUTDATAFILE%CNAME
  ELSE
    COUTFILE = ''
  END IF
END IF
!
!*       3.    initialisation of 2D arrays
! 
CALL GET_DIM_EXT_ll('B',NIU,NJU)
CALL GET_INDICE_ll (NIB,NJB,NIE,NJE)
CALL GET_GLOBALDIMS_ll (NIMAX_ll,NJMAX_ll)
NIB_ALL = 1 + JPHEXT
NIE_ALL = NIMAX_ll + JPHEXT
NJB_ALL = 1 + JPHEXT
NJE_ALL = NJMAX_ll + JPHEXT
NIU_ALL = NIMAX_ll + 2* JPHEXT
NJU_ALL = NJMAX_ll + 2* JPHEXT
!
!
!*       4.    initialisation 1D physical dimension and mask
! 
ILU = (NIE-NIB+1+2*NHALO)*(NJE-NJB+1+2*NHALO)
ILU_ALL = (NIE_ALL-NIB_ALL+1)*(NJE_ALL-NJB_ALL+1)
!
CMASK=HMASK
ILM = ILU
!
!
!* dimension and mask for entire field
ALLOCATE(NMASK_ALL(ILU_ALL))
ALLOCATE(ZFULL(ILU_ALL))
ZFULL=1.
CALL GET_1D_MASK(ILU_ALL,ILU_ALL,ZFULL,NMASK_ALL)
DEALLOCATE(ZFULL)
!
!* dimension and mask for distributed field on processors
CALL GET_TYPE_DIM_n(YSURF_CUR%DTCO,YSURF_CUR%U,HMASK,ILM)
ALLOCATE(NMASK(ILM))
CALL GET_SURF_MASK_n(YSURF_CUR%DTCO,YSURF_CUR%U,HMASK,ILM,NMASK,ILU,TOUT%NLU)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHINIT_IO_SURF_n
