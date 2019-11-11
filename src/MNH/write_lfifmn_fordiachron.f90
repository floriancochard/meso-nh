!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################################
      SUBROUTINE WRITE_LFIFMN_FORDIACHRO_n(TPFILE)
!     #############################################
!
!!****  *WRITE_LFIFM_FORDIACHRO_n* - routine to write a LFIFM file for model _n
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write an initial LFIFM File
!     of name YFMFILE//'.lfi' with the FM routines.
!
!!**  METHOD
!!    ------
!!      The data are written in the LFIFM file :
!!        - dimensions
!!        - grid variables
!!        - configuration variables
!!        - prognostic variables at time t and t-dt
!!        - 1D anelastic reference state
!!
!!      The localization on the model grid is also indicated :
!!
!!        IGRID = 1 for mass grid point
!!        IGRID = 2 for U grid point
!!        IGRID = 3 for V grid point
!!        IGRID = 4 for w grid point
!!        IGRID = 0 for meaningless case
!!
!!
!!    EXTERNAL
!!    --------
!!      FMWRIT : FM-routine to write a record
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DIM_n   : contains dimensions
!!      Module MODD_TIME_n   : contains time variables and uses MODD_TIME
!!      Module MODD_GRID    : contains spatial grid variables for all models
!!      Module MODD_GRID_n : contains spatial grid variables
!!      Module MODD_REF     : contains reference state variables
!!      Module MODD_LUNIT_n: contains logical unit variables.
!!      Module MODD_CONF    : contains configuration variables for all models
!!      Module MODD_CONF_n  : contains configuration variables
!!      Module MODD_PARAM_n    : contains parameterization options
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/05/94
!!       V. Ducrocq    27/06/94
!!       J.Stein       20/10/94 (name of the FMFILE)
!!       J.Stein       06/12/94 add the LS fields
!!       J.P. Lafore   09/01/95 add the DRYMASST
!!       J.Stein       20/01/95 add TKE and change the ycomment for the water
!!                              variables
!!       J.Stein       23/01/95 add a TKE switch and MODD_PARAM1
!!       J.Stein       16/03/95 remove R from the historical variables
!!       J.Stein       20/03/95 add the EPS var.
!!       J.Stein       30/06/95 add the variables related to the subgrid condens
!!       S. Belair     01/09/95 add surface variables and ground parameters
!!       J.-P. Pinty   15/09/95 add the radiation parameters
!!       J. DURON      24/06/99 add GPACK varaible to disable pack option
!!       October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!       J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!       P.Wautelet: 11/07/2016 : removed MNH_NCWRIT define
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_TIME_n
USE MODD_PARAM_n
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_REF
USE MODD_REF_n, ONLY: XTHVREF, XRHODREF
USE MODD_LUNIT_n
USE MODD_TIME
USE MODD_TYPE_DATE
USE MODD_NESTING
USE MODD_IO_ll, ONLY: TFILEDATA
!
USE MODE_FMWRIT
USE MODE_GATHER_ll
USE MODE_GRIDPROJ
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA),INTENT(IN) :: TPFILE
!
!*       0.2   Declarations of local variables
!
INTEGER                           :: IRESP          ! return-code!
!
LOGICAL                           :: GPACK
!
REAL                              :: ZLATOR, ZLONOR ! geographical coordinates of 1st mass point
REAL                              :: ZXHATM, ZYHATM ! conformal    coordinates of 1st mass point
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll       !  Position x in the conformal
                                                    ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll       !   Position y in the conformal
                                                    ! plane (array on the complete domain)
!
!-------------------------------------------------------------------------------
!
GPACK=LPACK
LPACK=.FALSE.
!
!*       1.0    Version :
!
CALL IO_WRITE_FIELD(TPFILE,'L1D', L1D)
CALL IO_WRITE_FIELD(TPFILE,'L2D', L2D)
CALL IO_WRITE_FIELD(TPFILE,'PACK',LPACK)
CALL IO_WRITE_FIELD(TPFILE,'SURF',CSURF)
!
!*       1.1    Dimensions :
!
CALL IO_WRITE_FIELD(TPFILE,'IMAX',NIMAX_ll)
CALL IO_WRITE_FIELD(TPFILE,'JMAX',NJMAX_ll)
CALL IO_WRITE_FIELD(TPFILE,'KMAX',NKMAX)
!
CALL IO_WRITE_FIELD(TPFILE,'JPHEXT',JPHEXT)
!
!*       1.2    Grid variables :
!
IF (.NOT.LCARTESIAN) THEN
!
  CALL IO_WRITE_FIELD(TPFILE,'RPK',   XRPK)
  CALL IO_WRITE_FIELD(TPFILE,'LONORI',XLONORI)
  CALL IO_WRITE_FIELD(TPFILE,'LATORI',XLATORI)
  !
  !* diagnostic of 1st mass point
  !
  ALLOCATE(ZXHAT_ll(NIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(NJMAX_ll+2 * JPHEXT))
  CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP) !//
  CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP) !//
  ZXHATM = 0.5 * (ZXHAT_ll(1)+ZXHAT_ll(2))
  ZYHATM = 0.5 * (ZYHAT_ll(1)+ZYHAT_ll(2))
  CALL SM_LATLON(XLATORI,XLONORI,ZXHATM,ZYHATM,ZLATOR,ZLONOR)
  DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
!
  CALL IO_WRITE_FIELD(TPFILE,'LONOR',ZLONOR)
  CALL IO_WRITE_FIELD(TPFILE,'LATOR',ZLATOR)
END IF
!
CALL IO_WRITE_FIELD(TPFILE,'THINSHELL',LTHINSHELL)
CALL IO_WRITE_FIELD(TPFILE,'LAT0',XLAT0)
CALL IO_WRITE_FIELD(TPFILE,'LON0',XLON0)
CALL IO_WRITE_FIELD(TPFILE,'BETA',XBETA)
!
CALL IO_WRITE_FIELD(TPFILE,'XHAT',XXHAT)
CALL IO_WRITE_FIELD(TPFILE,'YHAT',XYHAT)
CALL IO_WRITE_FIELD(TPFILE,'ZHAT',XZHAT)
CALL IO_WRITE_FIELD(TPFILE,'ZTOP',XZTOP)
!
CALL IO_WRITE_FIELD(TPFILE,'ZS',   XZS)
CALL IO_WRITE_FIELD(TPFILE,'ZSMT', XZSMT)
CALL IO_WRITE_FIELD(TPFILE,'SLEVE',LSLEVE)
!
IF (LSLEVE) THEN
  CALL IO_WRITE_FIELD(TPFILE,'LEN1',XLEN1)
  CALL IO_WRITE_FIELD(TPFILE,'LEN2',XLEN2)
END IF
!
CALL IO_WRITE_FIELD(TPFILE,'DTMOD',TDTMOD)
CALL IO_WRITE_FIELD(TPFILE,'DTCUR',TDTCUR)
CALL IO_WRITE_FIELD(TPFILE,'DTEXP',TDTEXP)
CALL IO_WRITE_FIELD(TPFILE,'DTSEG',TDTSEG)
!
!*       1.3    Configuration  variables :
!
CALL IO_WRITE_FIELD(TPFILE,'CARTESIAN',LCARTESIAN)
CALL IO_WRITE_FIELD(TPFILE,'LBOUSS',   LBOUSS)
!
!*       1.6    Reference state variables :
!
CALL IO_WRITE_FIELD(TPFILE,'RHOREFZ',XRHODREFZ)
CALL IO_WRITE_FIELD(TPFILE,'THVREFZ',XTHVREFZ)
CALL IO_WRITE_FIELD(TPFILE,'EXNTOP', XEXNTOP)
!
CALL IO_WRITE_FIELD(TPFILE,'RHODREF',XRHODREF)
CALL IO_WRITE_FIELD(TPFILE,'THVREF', XTHVREF)
!
LPACK=GPACK
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_LFIFMN_FORDIACHRO_n
