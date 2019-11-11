!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!###########################
MODULE MODI_INIT_FOR_CONVLFI
!###########################
!
!
INTERFACE
      SUBROUTINE INIT_FOR_CONVLFI(TPINIFILE)
!
USE MODD_IO_ll,ONLY: TFILEDATA
!
TYPE(TFILEDATA),        INTENT(IN)    :: TPINIFILE   ! file being read
!
END SUBROUTINE INIT_FOR_CONVLFI
END INTERFACE
END MODULE MODI_INIT_FOR_CONVLFI
!
!     ############################################
      SUBROUTINE INIT_FOR_CONVLFI(TPINIFILE)
!     ############################################
!
!!****  *INIT_FOR_CONVLFI * - light monitor to initialize the variables 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize some variables   
!     necessary in the conversion program.
!
!!**  METHOD
!!    ------
!!      This initialization takes some parts of the whole initialization modules
!!    of monitor INIT: 
!!        geometry and dimensions from ini_sizen
!!        grids, metric coefficients, dates and times from set_grid
!!        reading of the pressure field
!!             
!!
!!    EXTERNAL
!!    --------
!!      INI_CST    : to initialize physical constants
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	I. Mallet       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                 20/02/01 
!!      J.-P. Pinty and D. Gazen 31/03/04 Add the 2D capability for V5D plots
!!    10/10/2011  J.Escobar call INI_PARAZ_ll
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_IO_ll,ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_CST
USE MODD_DIM_n
USE MODD_FIELD_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_TIME
USE MODD_TIME_n
USE MODD_VAR_ll, ONLY : NPROC
!
USE MODE_FIELD, ONLY : TFIELDDATA,TFIELDLIST,FIND_FIELD_ID_FROM_MNHNAME
USE MODE_TIME
USE MODE_GRIDPROJ
USE MODE_GRIDCART
!
USE MODE_FM
USE MODE_FMREAD
USE MODE_GATHER_ll
USE MODE_IO_ll
USE MODE_ll
!
USE MODI_INI_CST
!JUANZ
USE MODE_SPLITTINGZ_ll
!JUANZ
!
IMPLICIT NONE
!
!*       0.1   Arguments variables
!
TYPE(TFILEDATA),        INTENT(IN)    :: TPINIFILE   ! file being read
!
!*       0.2   Local variables
!
INTEGER  :: IRESP
CHARACTER (LEN=40)     :: YTITLE               ! Title for date print
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZJ      ! Jacobian
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
REAL                         :: ZXHATM,ZYHATM    ! coordinates of mass point 
REAL                         :: ZLATORI, ZLONORI ! lat and lon of left-bottom point
!
INTEGER             :: IIU,IJU       ! Upper dimension in x,y direction (local)
INTEGER             :: IKU           ! Upper dimension in z direction
INTEGER             :: IINFO_ll      ! return code of // routines
INTEGER             :: IID
TYPE(TFIELDDATA)    :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZE EACH MODEL SIZES AND DEPENDENCY (ini_sizen)
!              ------------------------------------------
!
!*       1.1   Read the geometry kind in the LFIFM file (Cartesian or spherical)
!
CALL IO_READ_FIELD(TPINIFILE,'CARTESIAN',LCARTESIAN)
!
!*       1.2  Read configuration and dimensions in initial file and initialize
!             subdomain dimensions and parallel variables
!
CALL IO_READ_FIELD(TPINIFILE,'IMAX',NIMAX_ll)
CALL IO_READ_FIELD(TPINIFILE,'JMAX',NJMAX_ll)
!
CALL IO_READ_FIELD(TPINIFILE,'L1D',L1D,IRESP)
IF (IRESP/=0) THEN
  L1D=.FALSE.
  IF( (NIMAX_ll == 1).AND.(NJMAX_ll == 1) ) L1D=.TRUE.
ENDIF  
!
CALL IO_READ_FIELD(TPINIFILE,'L2D',L2D,IRESP)
IF (IRESP/=0) THEN
  L2D=.FALSE.
  IF( (NIMAX_ll /= 1).AND.(NJMAX_ll == 1) ) L2D=.TRUE.
ENDIF  
!
CALL IO_READ_FIELD(TPINIFILE,'PACK',LPACK,IRESP)
IF (IRESP/=0) LPACK=.TRUE.
!
CALL SET_FMPACK_ll(L1D,L2D,LPACK)
!
CALL IO_READ_FIELD(TPINIFILE,'KMAX',NKMAX)
!
CSPLIT ='BSPLITTING' ; NHALO = 1
CALL SET_SPLITTING_ll(CSPLIT)
CALL SET_JP_ll(1,JPHEXT,JPVEXT, NHALO)
CALL SET_DAD0_ll()
CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
CALL SET_FMPACK_ll(L1D,L2D,LPACK)
CALL SET_LBX_ll('OPEN', 1)
CALL SET_LBY_ll('OPEN', 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(NIMAX_ll+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(NJMAX_ll+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
!JUANZ CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!
!*       1.4  Compute sizes of arrays of the extended sub-domain (ini_modeln)
!
IKU=NKMAX + 2*JPVEXT
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
!
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZE GRIDS AND METRIC COEFFICIENTS (set_grid)
!              ---------------------
!
!        2.1  reading
!
CALL IO_READ_FIELD(TPINIFILE,'LAT0',XLAT0)
CALL IO_READ_FIELD(TPINIFILE,'LON0',XLON0)
CALL IO_READ_FIELD(TPINIFILE,'BETA',XBETA)
CALL IO_READ_FIELD(TPINIFILE,'XHAT',XXHAT)
CALL IO_READ_FIELD(TPINIFILE,'YHAT',XYHAT)
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_READ_FIELD(TPINIFILE,'RPK',XRPK)
  CALL IO_READ_FIELD(TPINIFILE,'LONORI',XLONORI)
  CALL IO_READ_FIELD(TPINIFILE,'LATORI',XLATORI)
  !
  IF (TPINIFILE%NMNHVERSION(1)<4 .OR. (TPINIFILE%NMNHVERSION(1)==4 .AND. TPINIFILE%NMNHVERSION(2)<=5)) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('LONORI',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CMNHNAME = 'LONOR'
    CALL IO_READ_FIELD(TPINIFILE,TZFIELD,XLONORI)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('LATORI',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CMNHNAME = 'LATOR'
    CALL IO_READ_FIELD(TPINIFILE,TZFIELD,XLATORI)
    !
    ALLOCATE(ZXHAT_ll(NIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(NJMAX_ll+2 * JPHEXT))
    CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP) !//
    CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP) !//
    ZXHATM = - 0.5 * (ZXHAT_ll(1)+ZXHAT_ll(2))
    ZYHATM = - 0.5 * (ZYHAT_ll(1)+ZYHAT_ll(2))
    CALL SM_LATLON(XLATORI,XLONORI,ZXHATM,ZYHATM,ZLATORI,ZLONORI)
    DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
    XLATORI = ZLATORI
    XLONORI = ZLONORI
  END IF
END IF
!
ALLOCATE(XZS(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'ZS',XZS,IRESP)
IF (IRESP/=0) XZS(:,:)=0.
!
ALLOCATE(XZSMT(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'ZSMT',XZSMT,IRESP)
IF (IRESP/=0) XZSMT(:,:)=XZS(:,:)
!
ALLOCATE(XZHAT(IKU))
CALL IO_READ_FIELD(TPINIFILE,'ZHAT',XZHAT)
CALL IO_READ_FIELD(TPINIFILE,'ZTOP',XZTOP)
!
CALL IO_READ_FIELD(TPINIFILE,'SLEVE',LSLEVE,IRESP)
IF (IRESP/=0) LSLEVE = .FALSE.
!
IF (LSLEVE) THEN
  CALL IO_READ_FIELD(TPINIFILE,'LEN1',XLEN1)
  CALL IO_READ_FIELD(TPINIFILE,'LEN2',XLEN2)
END IF
!
CALL IO_READ_FIELD(TPINIFILE,'DTEXP',TDTEXP)
CALL IO_READ_FIELD(TPINIFILE,'DTMOD',TDTMOD)
CALL IO_READ_FIELD(TPINIFILE,'DTSEG',TDTSEG)
CALL IO_READ_FIELD(TPINIFILE,'DTCUR',TDTCUR)
!
YTITLE='CURRENT DATE AND TIME'
CALL SM_PRINT_TIME(TDTCUR,TLUOUT,YTITLE)
!
!*       2.2    Spatial grid
! 
ALLOCATE(XDXHAT(IIU))
ALLOCATE(XDYHAT(IJU))
ALLOCATE(XZZ(IIU,IJU,IKU))
ALLOCATE(ZJ(IIU,IJU,IKU))
!
CALL INI_CST
!
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XDXHAT,XDYHAT,XZZ,ZJ) 
ELSE
  ALLOCATE(XLON(IIU,IJU))
  ALLOCATE(XLAT(IIU,IJU))
  ALLOCATE(XMAP(IIU,IJU))
  CALL SM_GRIDPROJ(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XLATORI,XLONORI, &
                   XMAP,XLAT,XLON,XDXHAT,XDYHAT,XZZ,ZJ)  
END IF    
!
!-------------------------------------------------------------------------------
!
!*       3.    INITIALIZE THE PROGNOSTIC AND SURFACE FIELDS (read_field)
!              --------------------------------------------
ALLOCATE(XPABST(IIU,IJU,IKU))
CALL IO_READ_FIELD(TPINIFILE,'PABST',XPABST)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_FOR_CONVLFI
