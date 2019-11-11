! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     #############################################
      SUBROUTINE WRITE_LFIFMN_FORDIACHRO_n(HFMFILE)
!     #############################################
!
!!****  *WRITE_LFIFM_FORDIACHRO_n* - routine to write a LFIFM file for model _n
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write an initial LFIFM File
!     of name HFMFILE//'.lfi' with the FM routines.
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
USE MODE_FMWRIT
USE MODE_GRIDPROJ
USE MODE_ll
USE MODI_GATHER_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
CHARACTER(LEN=28), INTENT(IN) :: HFMFILE      ! Name of FM-file to write
!
!*       0.2   Declarations of local variables
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
                                    !  at the open of the file                                                                      !  LFI  routines
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be written
CHARACTER(LEN=100) :: YCOMMENT       ! Comment string
!
LOGICAL                :: GPACK
!
REAL                              :: ZLATOR, ZLONOR ! geographical coordinates of 1st mass point
REAL                              :: ZXHATM, ZYHATM ! conformal    coordinates of 1st mass point
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
!
!-------------------------------------------------------------------------------
!
!JUAN
RETURN
!JUAN
GPACK=LPACK
LPACK=.FALSE.
!
!*       1.0    Version :
!
YRECFM='MASDEV'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',NMASDEV,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='BUGFIX'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',NBUGFIX,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='BIBUSER'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',CBIBUSER,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='PROGRAM'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',CPROGRAM,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='L1D'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',L1D,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='L2D'
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',L2D,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='PACK'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',LPACK,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='SURF'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',CSURF,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='STORAGE_TYPE'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',CSTORAGE_TYPE,IGRID,ILENCH,YCOMMENT,IRESP)
!
!*       1.1    Dimensions :
!
YRECFM='IMAX'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',NIMAX_ll,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='JMAX'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',NJMAX_ll,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='KMAX'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',NKMAX,IGRID,ILENCH,YCOMMENT,IRESP)
!
!*       1.2    Grid variables :
!
IF (.NOT.LCARTESIAN) THEN
!
  YRECFM='RPK'
  !CALL ELIM(YRECFM)
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XRPK,IGRID,ILENCH,YCOMMENT,IRESP)
!
  YRECFM='LONORI'
  !CALL ELIM(YRECFM)
  YCOMMENT='DEGREES'
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XLONORI,IGRID,ILENCH,YCOMMENT,IRESP)
!
  YRECFM='LATORI'
  !CALL ELIM(YRECFM)
  YCOMMENT='DEGREES'
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XLATORI,IGRID,ILENCH,YCOMMENT,IRESP)
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
  YRECFM='LONOR'
  YCOMMENT='DEGREES'
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',ZLONOR,IGRID,ILENCH,YCOMMENT,IRESP)
!
  YRECFM='LATOR'
  YCOMMENT='DEGREES'
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',ZLATOR,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
YRECFM='THINSHELL'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',LTHINSHELL,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='LAT0'
!CALL ELIM(YRECFM)
YCOMMENT='reference latitude for conformal projection (DEGREES)'
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XLAT0,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='LON0'
!CALL ELIM(YRECFM)
YCOMMENT='reference longitude for conformal projection (DEGREES)'
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XLON0,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='BETA'
!CALL ELIM(YRECFM)
YCOMMENT='rotation angle (DEGREES)'
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XBETA,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='XHAT'
!CALL ELIM(YRECFM)
YCOMMENT='Position x in the conformal or cartesian plane (METERS)'
IGRID=2
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'XX',XXHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='YHAT'
!CALL ELIM(YRECFM)
YCOMMENT='Position y in the conformal or cartesian plane (METERS)'
IGRID=3
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'YY',XYHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='ZHAT'
!CALL ELIM(YRECFM)
YCOMMENT='height level without orography (METERS)'
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XZHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='ZS'
!CALL ELIM(YRECFM)
YCOMMENT='orography (METERS)'
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'XY',XZS,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='ZSMT'
!CALL ELIM(YRECFM)
YCOMMENT='smooth orography (METERS)'
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'XY',XZSMT,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='SLEVE'
!CALL ELIM(YRECFM)
YCOMMENT=' '
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',LSLEVE,IGRID,ILENCH,YCOMMENT,IRESP)
!
IF (LSLEVE) THEN
  YRECFM='LEN1'
  !CALL ELIM(YRECFM)
  YCOMMENT='METERS'
  IGRID=4
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XLEN1,IGRID,ILENCH,YCOMMENT,IRESP)
  YRECFM='LEN2'
  !CALL ELIM(YRECFM)
  YCOMMENT='METERS'
  IGRID=4
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XLEN2,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
YRECFM='DTCUR'
!CALL ELIM(YRECFM)
IGRID=0
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',TDTCUR,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='DTEXP'
!CALL ELIM(YRECFM)
IGRID=0
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',TDTEXP,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='DTMOD'
!CALL ELIM(YRECFM)
IGRID=0
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',TDTMOD,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='DTSEG'
!CALL ELIM(YRECFM)
IGRID=0
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',TDTSEG,IGRID,ILENCH,YCOMMENT,IRESP)
!
!*       1.3    Configuration  variables :
!
YRECFM='CARTESIAN'
!CALL ELIM(YRECFM)
YCOMMENT='Logical for cartesian geometry '
IGRID=0
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',LCARTESIAN,IGRID,ILENCH,YCOMMENT,IRESP)
!
!*       1.6    Reference state variables :
!
YRECFM='RHOREFZ'
!CALL ELIM(YRECFM)
YCOMMENT='rhodz for reference state without orography (kg/m3)'
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XRHODREFZ,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='THVREFZ'
!CALL ELIM(YRECFM)
YCOMMENT='thetavz for reference state without orography (K)'
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XTHVREFZ,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='EXNTOP'
!CALL ELIM(YRECFM)
YCOMMENT='Exner function at model top'
IGRID=4
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'--',XEXNTOP,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='RHODREF'
!CALL ELIM(YRECFM)
YCOMMENT='Dry density for reference state with orography (kg/m3)'
IGRID=1
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'XY',XRHODREF,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='THVREF'
YCOMMENT='Thetav for reference state with orography (K)'
YCOMMENT='  '
IGRID=1
ILENCH=LEN(YCOMMENT)
CALL FMWRIT(HFMFILE,YRECFM,CLUOUT,'XY',XTHVREF,IGRID,ILENCH,YCOMMENT,IRESP)
!
LPACK=GPACK
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_LFIFMN_FORDIACHRO_n
