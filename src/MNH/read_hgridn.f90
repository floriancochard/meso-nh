!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_READ_HGRID_n
!     #######################
!
INTERFACE
      SUBROUTINE READ_HGRID_n(TPFMFILE,HMY_NAME,HDAD_NAME,HSTORAGE_TYPE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),   INTENT(IN)  :: TPFMFILE     ! file n
CHARACTER(LEN=28), INTENT(OUT) :: HMY_NAME     ! True Name of FM-file
CHARACTER(LEN=28), INTENT(OUT) :: HDAD_NAME    ! Name of father
CHARACTER(LEN=2) , INTENT(OUT) :: HSTORAGE_TYPE
!
END SUBROUTINE READ_HGRID_n
!
END INTERFACE
END MODULE MODI_READ_HGRID_n
!
!     #################################################################
      SUBROUTINE READ_HGRID_n(TPFMFILE,HMY_NAME,HDAD_NAME,HSTORAGE_TYPE)
!     #################################################################
!
!!****  *READ_HGRID_n* - to read grid information in FM file of model $n
!!
!!    PURPOSE
!!    -------
!
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
!!      Module MODD_GRID_n : contains domain definition
!!        XXHAT
!!        XYHAT
!!      Module MODD_DIM_n : contains domain size
!!        NIMAX
!!        NJMAX
!!      Module MODD_PARAMETERS :
!!        JPHEXT
!!      Module MODD_LUNIT :
!!        CLUOUT
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
!!         M.Faivre     2014
!!         M.Moge       06/2015 case ( CPROGRAM .EQ. "NESPGD"  .OR. CPROGRAM .EQ. "SPAWN ")
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CONF
USE MODD_DIM_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LUNIT_n
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT, JPMODELMAX
!
USE MODE_FIELD, ONLY : TFIELDDATA,TFIELDLIST,FIND_FIELD_ID_FROM_MNHNAME
USE MODE_FMREAD
USE MODE_GRIDPROJ
USE MODE_IO_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
TYPE(TFILEDATA),   INTENT(IN)  :: TPFMFILE     ! file n
CHARACTER(LEN=28), INTENT(OUT) :: HMY_NAME     ! True Name of FM-file
CHARACTER(LEN=28), INTENT(OUT) :: HDAD_NAME    ! Name of father
CHARACTER(LEN=2) , INTENT(OUT) :: HSTORAGE_TYPE
!
!*       0.2   declarations of local variables
!
INTEGER             :: ILUOUT
INTEGER             :: IRESP
REAL                :: ZLAT0,ZLON0,ZRPK,ZBETA
REAL                :: ZEPS = 1.E-10
INTEGER             :: IID, IMI
!
!-------------------------------------------------------------------------------
REAL :: ZLATOR, ZLONOR, ZXHATM, ZYHATM
!-------------------------------------------------------------------------------
!JUAN REALZ
INTEGER             :: IIU,IJU
!JUAN REALZ
INTEGER             :: IXOR, IYOR, IXEND, IYEND
INTEGER             :: IJPHEXT
TYPE(TFIELDDATA)    :: TZFIELD
!
ILUOUT = TLUOUT%NLU
!
!*       1.     General information :
!               -------------------
!
CALL IO_READ_FIELD(TPFMFILE,'MY_NAME',     HMY_NAME)
CALL IO_READ_FIELD(TPFMFILE,'DAD_NAME',    HDAD_NAME)
CALL IO_READ_FIELD(TPFMFILE,'STORAGE_TYPE',HSTORAGE_TYPE)
!
!*       2.     Grid information :
!               ----------------
!
IF( (TPFMFILE%NMNHVERSION(1)<4 .OR. (TPFMFILE%NMNHVERSION(1)==4 .AND. TPFMFILE%NMNHVERSION(2)<=5)) &
   .AND. HSTORAGE_TYPE == 'PG') THEN
  LCARTESIAN=.FALSE.
ELSE
  CALL IO_READ_FIELD(TPFMFILE,'CARTESIAN',LCARTESIAN)
ENDIF
CALL IO_READ_FIELD(TPFMFILE,'LAT0',ZLAT0)
CALL IO_READ_FIELD(TPFMFILE,'LON0',ZLON0)
CALL IO_READ_FIELD(TPFMFILE,'BETA',ZBETA,IRESP)
IF(IRESP/=0) ZBETA=0.
IF (.NOT.LCARTESIAN ) THEN
  CALL IO_READ_FIELD(TPFMFILE,'RPK',   ZRPK)
  CALL IO_READ_FIELD(TPFMFILE,'LATORI',XLATORI)
  CALL IO_READ_FIELD(TPFMFILE,'LONORI',XLONORI)
ENDIF
!
IMI = GET_CURRENT_MODEL_INDEX()
IF (IMI == 1) THEN
  XLAT0=ZLAT0
  XLON0=ZLON0
  XBETA=ZBETA
  IF (.NOT.LCARTESIAN) XRPK=ZRPK
ELSE
  IF (     ABS(XLAT0-ZLAT0)> ZEPS .OR. ABS(XLON0-ZLON0)> ZEPS  &
                                  .OR. ABS(XBETA-ZBETA)> ZEPS  ) THEN
    WRITE(ILUOUT,*) 'projections are different in the two input files:'
    WRITE(ILUOUT,*) 'model ',IMI,' : XLAT0= ',ZLAT0,' XLON0= ',ZLON0, &
                                               ' XBETA= ',ZBETA
    WRITE(ILUOUT,*) 'model 1 : XLAT0= ',XLAT0,' XLON0= ',XLON0, &
                                              ' XBETA= ',XBETA
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_HGRID_n','')
  END IF
  IF (.NOT.LCARTESIAN ) THEN
    IF ( ABS(XRPK-ZRPK)> ZEPS ) THEN
      WRITE(ILUOUT,*) 'projections are different in the two input files:'
      WRITE(ILUOUT,*) 'model ',IMI,' : XRPK= ',ZRPK
      WRITE(ILUOUT,*) 'model 1 : XRPK= ',XRPK
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_HGRID_n','')
    END IF
  END IF
END IF
!
IF (CPROGRAM/='IDEAL ') THEN
  !* WARNING : the following initialization of dimensions is ONLY valid for 
  !            monoprocessor runs, or if :
  !            a) NIMAX_ll, NJMAX_ll, and corresponding NIMAX_ll, NJMAX_ll are
  !               correctly initialized in later routines (e.g. spawn_model2.f90)
  !            b) and arrays XXHAT, XYHAT, XZS, XZSMT are deallocated after this 
  !               routine (as in ini_size_spawn.f90)
  !$20140506 try 'XX','YY' it is FMREADN0_LL scalar reading so leave '--'
  CALL IO_READ_FIELD(TPFMFILE,'IMAX',  NIMAX)
  CALL IO_READ_FIELD(TPFMFILE,'JMAX',  NJMAX)
  CALL IO_READ_FIELD(TPFMFILE,'JPHEXT',IJPHEXT)
  IF ( IJPHEXT .NE. JPHEXT ) THEN
     IF (CPROGRAM == 'REAL' ) THEN
        WRITE(ILUOUT,FMT=*) ' READ_HGRID_N : JPHEXT in PRE_REAL1.nam/NAM_REAL_CONF ( or default value )&
           & JPHEXT=',JPHEXT
     ELSE
        WRITE(ILUOUT,FMT=*) ' READ_HGRID_N : JPHEXT in PRE_NEST_PGD1.nam/NAM_CONF_NEST ( or default value )&
           & JPHEXT=',JPHEXT
     END IF

     WRITE(ILUOUT,FMT=*) ' different from PGD files=',TPFMFILE%CNAME ,' value JPHEXT=',IJPHEXT
     WRITE(ILUOUT,FMT=*) '-> JOB ABORTED'
     CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_HGRID_n','')
  END IF
END IF
!
!*       2.1  Read the configuration (MODD_CONF)
!
IF (IMI == 1) THEN   
  CALL IO_READ_FIELD(TPFMFILE,'L1D',L1D,IRESP)
  IF (IRESP/=0) THEN
    L1D=.FALSE.
    IF( (NIMAX == 1).AND.(NJMAX == 1) ) L1D=.TRUE.
  ENDIF
!
  CALL IO_READ_FIELD(TPFMFILE,'L2D',L2D,IRESP)
  IF (IRESP/=0) THEN
    L2D=.FALSE.
    IF( (NIMAX /= 1).AND.(NJMAX == 1) ) L2D=.TRUE.
  ENDIF
!
  CALL IO_READ_FIELD(TPFMFILE,'PACK',LPACK,IRESP)
  IF (IRESP/=0) LPACK=.TRUE.
!  CALL SET_FMPACK_ll(L1D,L2D,LPACK)
END IF
!
!*       2.2    Grid information :
!               ----------------
!JUAN REALZ
IF ( CPROGRAM .EQ. "REAL  " ) THEN
  CALL GET_DIM_EXT_ll('B',IIU,IJU)
  CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
  IF (.NOT. (ASSOCIATED(XXHAT))) ALLOCATE(XXHAT(IIU))
  IF (.NOT. (ASSOCIATED(XYHAT))) ALLOCATE(XYHAT(IJU))
ELSE IF ( CPROGRAM .EQ. "NESPGD"  .OR. CPROGRAM .EQ. "SPAWN ") THEN
  NIMAX_ll = NIMAX
  NJMAX_ll = NJMAX
  CALL GET_INDICE_ll( IXOR, IYOR, IXEND, IYEND )
  NIMAX = IXEND - IXOR + 1
  NJMAX = IYEND - IYOR + 1
  IIU = NIMAX+2*JPHEXT
  IJU = NJMAX+2*JPHEXT
  IF (.NOT. (ASSOCIATED(XXHAT))) ALLOCATE(XXHAT(IIU))
  IF (.NOT. (ASSOCIATED(XYHAT))) ALLOCATE(XYHAT(IJU))
ELSE
  IF (.NOT. (ASSOCIATED(XXHAT))) ALLOCATE(XXHAT(NIMAX+2*JPHEXT))
  IF (.NOT. (ASSOCIATED(XYHAT))) ALLOCATE(XYHAT(NJMAX+2*JPHEXT))
ENDIF
!JUAN REALZ

CALL IO_READ_FIELD(TPFMFILE,'XHAT',XXHAT)
CALL IO_READ_FIELD(TPFMFILE,'YHAT',XYHAT)
!
!JUAN REALZ
IF ( CPROGRAM .EQ. "REAL  " ) THEN
IF (.NOT. (ASSOCIATED(XZS))) ALLOCATE(XZS(IIU,IJU))
ELSE
IF (.NOT. (ASSOCIATED(XZS))) ALLOCATE(XZS(NIMAX+2*JPHEXT,NJMAX+2*JPHEXT))
ENDIF
!JUAN REALZ

CALL IO_READ_FIELD(TPFMFILE,'ZS',XZS)
!
!JUAN REALZ
IF ( CPROGRAM .EQ. "REAL  " ) THEN
IF (.NOT. (ASSOCIATED(XZSMT))) ALLOCATE(XZSMT(IIU,IJU))
ELSE
IF (.NOT. (ASSOCIATED(XZSMT))) ALLOCATE(XZSMT(NIMAX+2*JPHEXT,NJMAX+2*JPHEXT))
ENDIF
!JUAN REALZ

IF (TPFMFILE%NMNHVERSION(1)<4 .OR. (TPFMFILE%NMNHVERSION(1)==4 .AND. TPFMFILE%NMNHVERSION(2)<=6)) THEN
  XZSMT = XZS
ELSE
  CALL IO_READ_FIELD(TPFMFILE,'ZSMT',XZSMT)
!
END IF
!
!-------------------------------------------------------------------------------
IF (TPFMFILE%NMNHVERSION(1)<4 .OR. (TPFMFILE%NMNHVERSION(1)==4 .AND. TPFMFILE%NMNHVERSION(2)<=5)) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME('LONORI',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CMNHNAME = 'LONOR'
  CALL IO_READ_FIELD(TPFMFILE,TZFIELD,XLONORI)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('LATORI',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CMNHNAME = 'LATOR'
  CALL IO_READ_FIELD(TPFMFILE,TZFIELD,XLATORI)
  !
  ZXHATM = - 0.5 * (XXHAT(1)+XXHAT(2))
  ZYHATM = - 0.5 * (XYHAT(1)+XYHAT(2))
  CALL SM_LATLON(XLATORI,XLONORI,ZXHATM,ZYHATM,ZLATOR,ZLONOR)
  XLATORI = ZLATOR
  XLONORI = ZLONOR
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_HGRID_n
