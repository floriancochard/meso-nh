!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_READ_HGRID
!     ######################
INTERFACE
      SUBROUTINE READ_HGRID(KMI,TPFMFILE,HMY_NAME,HDAD_NAME,HSTORAGE_TYPE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
INTEGER,           INTENT(IN)  :: KMI          ! model index
TYPE(TFILEDATA),   INTENT(IN)  :: TPFMFILE     ! file n
CHARACTER(LEN=28), INTENT(OUT) :: HMY_NAME     ! True Name of FM-file
CHARACTER(LEN=28), INTENT(OUT) :: HDAD_NAME    ! Name of father
CHARACTER(LEN=2) , INTENT(OUT) :: HSTORAGE_TYPE
!
END SUBROUTINE READ_HGRID
END INTERFACE
END MODULE MODI_READ_HGRID
!
!     ####################################################################
      SUBROUTINE READ_HGRID(KMI,TPFMFILE,HMY_NAME,HDAD_NAME,HSTORAGE_TYPE)
!     ####################################################################
!
!!****  *READ_HGRID* - to read grid information in FM file into PGD modules
!!                     (KMI==0) or for model KMI
!!
!!    PURPOSE
!!    -------
!!
!!    CAUTION : if used to fill the MODD_PGD... modules, the projection
!!              definition module MODD_GRID will be eraised without test:
!!              In this case (KMI==0), it is used to define the
!!              projection.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      FMREAD   : to read data in LFIFM file
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_GRID : contains projection definition
!!        XLAT0
!!        XLON0
!!        XRPK
!!        XBETA
!!        XLATORI
!!        XLONORI
!!      Module MODD_GRID$n : contains domain definition
!!        XXHAT
!!        XYHAT
!!      Module MODD_DIM$n : contains domain size
!!        NIMAX
!!        NJMAX
!!      Module MODD_PARAMETERS :
!!        JPHEXT
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
!!      Original        26/09/96
!!            M.Faivre      2014
!!            G.Delautier 2017 BUG for MNH2LPDM
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CONF, ONLY : CPROGRAM
USE MODD_GRID
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_PGDDIM
USE MODD_PGDGRID
!
USE MODE_FIELD, ONLY : TFIELDDATA,TFIELDLIST,FIND_FIELD_ID_FROM_MNHNAME
USE MODE_FM,    ONLY : SET_FMPACK_ll
USE MODE_FMREAD
USE MODE_GRIDPROJ
USE MODE_IO_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
USE MODI_READ_HGRID_n
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,           INTENT(IN)  :: KMI          ! model index
TYPE(TFILEDATA),   INTENT(IN)  :: TPFMFILE     ! file n
CHARACTER(LEN=28), INTENT(OUT) :: HMY_NAME     ! True Name of FM-file
CHARACTER(LEN=28), INTENT(OUT) :: HDAD_NAME    ! Name of father
CHARACTER(LEN=2) , INTENT(OUT) :: HSTORAGE_TYPE
!
!
!*       0.2   declarations of local variables
!
INTEGER                :: IRESP
INTEGER                :: IID, IMI
LOGICAL                :: G1D,G2D,GPACK
INTEGER                :: IINFO_ll
REAL                   :: ZLATOR, ZLONOR, ZXHATM, ZYHATM
TYPE(TFIELDDATA)       :: TZFIELD
!-------------------------------------------------------------------------------
!
!*       1.     TEST ON MODEL INDEX
!               -------------------
!
! KMI may be 0 !
IF (KMI<0 .OR. KMI>JPMODELMAX) CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_HGRID','KMI<0 .OR. KMI>JPMODELMAX')
IF (KMI/=0) THEN
  IF (CPROGRAM/='M2LPDM') THEN
    IMI = GET_CURRENT_MODEL_INDEX()
    CALL GOTO_MODEL(KMI)
    CALL GO_TOMODEL_ll(KMI, IINFO_ll)
    CALL READ_HGRID_n(TPFMFILE,HMY_NAME,HDAD_NAME,HSTORAGE_TYPE)
    CALL GO_TOMODEL_ll(IMI, IINFO_ll)
    CALL GOTO_MODEL(IMI)
    RETURN
  ELSE
    IMI = GET_CURRENT_MODEL_INDEX()
    CALL GOTO_MODEL(KMI)
    CALL READ_HGRID_n(TPFMFILE,HMY_NAME,HDAD_NAME,HSTORAGE_TYPE)
    CALL GOTO_MODEL(IMI)
    RETURN
    RETURN    
  END IF
END IF
!
!*       2.     READING IN MODD_PGD...
!               ----------------------
!
!*       2.1    General information :
!               -------------------
!
CALL IO_READ_FIELD(TPFMFILE,'MY_NAME',     HMY_NAME)
CALL IO_READ_FIELD(TPFMFILE,'DAD_NAME',    HDAD_NAME)
CALL IO_READ_FIELD(TPFMFILE,'STORAGE_TYPE',HSTORAGE_TYPE)
!
!*       2.2    Grid information :
!               ----------------
!
!20131010 recompute properly NPGDIMAX NPGDJMAX
!GET_DIM_PHYS_ll impact => 1st one no visible impact
CALL GET_DIM_PHYS_ll ( 'B',NPGDIMAX,NPGDJMAX)
!
CALL IO_READ_FIELD(TPFMFILE,'LAT0',  XLAT0)
CALL IO_READ_FIELD(TPFMFILE,'LON0',  XLON0)
CALL IO_READ_FIELD(TPFMFILE,'RPK',   XRPK)
CALL IO_READ_FIELD(TPFMFILE,'BETA',  XBETA)
CALL IO_READ_FIELD(TPFMFILE,'LATORI',XPGDLATOR)
CALL IO_READ_FIELD(TPFMFILE,'LONORI',XPGDLONOR)
CALL IO_READ_FIELD(TPFMFILE,'IMAX',  NPGDIMAX)
CALL IO_READ_FIELD(TPFMFILE,'JMAX',  NPGDJMAX)
!
!20131010 recompute properly NPGDIMAX NPGDJMAX
!GET_DIM_PHYS_ll impact 2nd one => prevent run failures
CALL GET_DIM_PHYS_ll ( 'B',NPGDIMAX,NPGDJMAX)
!
IF (.NOT.(ALLOCATED(XPGDXHAT))) ALLOCATE(XPGDXHAT(NPGDIMAX+2*JPHEXT))
IF (.NOT.(ALLOCATED(XPGDYHAT))) ALLOCATE(XPGDYHAT(NPGDJMAX+2*JPHEXT))
!20131023 change FMREAD option '--' -> 'XX' ou 'YY' for // reading
CALL IO_READ_FIELD(TPFMFILE,'XHAT',XPGDXHAT)
CALL IO_READ_FIELD(TPFMFILE,'YHAT',XPGDYHAT)
!
!*       3.   Read the configuration (MODD_CONF)
!
CALL IO_READ_FIELD(TPFMFILE,'L1D',G1D,IRESP)
IF (IRESP/=0) THEN
  G1D=.FALSE.
  IF( (NPGDIMAX == 1).AND.(NPGDJMAX == 1) ) G1D=.TRUE.
ENDIF
!
CALL IO_READ_FIELD(TPFMFILE,'L2D',G2D,IRESP)
IF (IRESP/=0) THEN
  G2D=.FALSE.
  IF( (NPGDIMAX /= 1).AND.(NPGDJMAX == 1) ) G2D=.TRUE.
ENDIF
!
CALL IO_READ_FIELD(TPFMFILE,'PACK',GPACK,IRESP)
IF (IRESP/=0) GPACK=.TRUE.
!
CALL SET_FMPACK_ll(G1D,G2D,GPACK)
!-------------------------------------------------------------------------------
IF (TPFMFILE%NMNHVERSION(1)<4 .OR. (TPFMFILE%NMNHVERSION(1)==4 .AND. TPFMFILE%NMNHVERSION(2)<=5)) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME('LONORI',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CMNHNAME = 'LONOR'
  CALL IO_READ_FIELD(TPFMFILE,TZFIELD,XPGDLONOR)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('LATORI',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CMNHNAME = 'LATOR'
  CALL IO_READ_FIELD(TPFMFILE,TZFIELD,XPGDLATOR)
  !
  ZXHATM = - 0.5 * (XPGDXHAT(1)+XPGDXHAT(2))
  ZYHATM = - 0.5 * (XPGDYHAT(1)+XPGDYHAT(2))
  CALL SM_LATLON(XPGDLATOR,XPGDLONOR,ZXHATM,ZYHATM,ZLATOR,ZLONOR)
  XPGDLATOR = ZLATOR
  XPGDLONOR = ZLONOR
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_HGRID
