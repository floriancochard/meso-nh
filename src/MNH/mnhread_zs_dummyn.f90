!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_MNHREAD_ZS_DUMMY_n
!     ##########################
INTERFACE
      SUBROUTINE MNHREAD_ZS_DUMMY_n(TPINIFILE)
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),    INTENT(IN)   :: TPINIFILE    !Initial file
!
END SUBROUTINE MNHREAD_ZS_DUMMY_n
!
END INTERFACE
!
END MODULE MODI_MNHREAD_ZS_DUMMY_n
!
!     ##########################################################################
      SUBROUTINE MNHREAD_ZS_DUMMY_n(TPINIFILE)
!     ##########################################################################
!
!!****  *MNHREAD_ZS_DUMMY_n* - reads zs and dummy surface fields
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to read the LFIFM part of 
!!       physiographic data file with the FM routines.
!!
!!
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
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_GRID_n,     ONLY : XZS
USE MODD_GR_FIELD_n, ONLY : XSSO_STDEV, XSSO_ANISOTROPY, XSSO_DIRECTION, XSSO_SLOPE, &
                            XAVG_ZS, XSIL_ZS, XMIN_ZS, XMAX_ZS
USE MODD_IO_ll,      ONLY : TFILEDATA
USE MODD_PARAM_n,    ONLY : CSURF
!
USE MODI_READ_DUMMY_GR_FIELD_n
!
USE MODE_ll
USE MODE_FMREAD
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(TFILEDATA),    INTENT(IN)   :: TPINIFILE    !Initial file
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IIU            ! X array size
INTEGER           :: IJU            ! Y array size
!
!-------------------------------------------------------------------------------
!
!*       1.     READS IN THE LFI FILE
!	        ---------------------
!
!*       1.0    General information :
!               -------------------
!
!*       1.1    Dimensions :
!               ----------
!
!* x and y dimensions in the file
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!
!
!*       1.3    Orography :
!               ---------
IF (.NOT.(ASSOCIATED(XZS))) THEN
  ALLOCATE(XZS(IIU,IJU))
  CALL IO_READ_FIELD(TPINIFILE,'ZS',XZS)
END IF
!
IF (CSURF /='EXTE') RETURN
!
!*       2.     Physiographic data fields:
!               -------------------------
!
!*       2.1    Orographic characteristics :
!               --------------------------
!
IF (.NOT.(ASSOCIATED(XSSO_ANISOTROPY))) ALLOCATE(XSSO_ANISOTROPY(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'SSO_ANIS',XSSO_ANISOTROPY(:,:))
!
IF (.NOT.(ASSOCIATED(XSSO_SLOPE))) ALLOCATE(XSSO_SLOPE(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'SSO_SLOPE',XSSO_SLOPE(:,:))
!
IF (.NOT.(ASSOCIATED(XSSO_DIRECTION))) ALLOCATE(XSSO_DIRECTION(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'SSO_DIR',XSSO_DIRECTION(:,:))
!
IF (.NOT.(ASSOCIATED(XAVG_ZS))) ALLOCATE(XAVG_ZS(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'AVG_ZS',XAVG_ZS(:,:))
!
IF (.NOT.(ASSOCIATED(XSIL_ZS))) ALLOCATE(XSIL_ZS(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'SIL_ZS',XSIL_ZS(:,:))
!
IF (.NOT.(ASSOCIATED(XMAX_ZS))) ALLOCATE(XMAX_ZS(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'MAX_ZS',XMAX_ZS(:,:))
!
IF (.NOT.(ASSOCIATED(XMIN_ZS))) ALLOCATE(XMIN_ZS(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'MIN_ZS',XMIN_ZS(:,:))
!
IF (.NOT.(ASSOCIATED(XSSO_STDEV))) ALLOCATE(XSSO_STDEV(IIU,IJU))
CALL IO_READ_FIELD(TPINIFILE,'SSO_STDEV',XSSO_STDEV(:,:))
!
!-------------------------------------------------------------------------------
!
!*      3.     Dummy fields
!              ------------
!
CALL READ_DUMMY_GR_FIELD_n(TPINIFILE,1,IIU,1,IJU,.TRUE.)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHREAD_ZS_DUMMY_n
