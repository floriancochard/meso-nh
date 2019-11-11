!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_MNHWRITE_ZS_DUMMY_n
!     ##########################
INTERFACE
      SUBROUTINE MNHWRITE_ZS_DUMMY_n(TPFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! File characteristics
!
END SUBROUTINE MNHWRITE_ZS_DUMMY_n
!
END INTERFACE
END MODULE MODI_MNHWRITE_ZS_DUMMY_n
!
!     ###################################################
      SUBROUTINE MNHWRITE_ZS_DUMMY_n(TPFILE)
!     ###################################################
!
!!****  *MNHWRITE_ZS_WH_DUMMY_n* - writes zs and dummy surface fields.
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to write the LFIFM part of 
!!       physiographic data file with the FM routines.
!!
!!
!!**  METHOD
!!    ------
!!          
!!    EXTERNAL
!!    --------
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   01/2004
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_GR_FIELD_n, ONLY : XSSO_STDEV, XSSO_ANISOTROPY, XSSO_DIRECTION, XSSO_SLOPE, &
                            XAVG_ZS, XSIL_ZS, XMIN_ZS, XMAX_ZS
!
USE MODD_PARAM_n,    ONLY : CSURF
USE MODD_IO_ll,      ONLY : TFILEDATA
!
USE MODI_WRITE_DUMMY_GR_FIELD_n
!
USE MODE_FMWRIT
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! File characteristics
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!*       1.     Orography :
!               ---------
!
!***********************************!
! Already written (in write_lfifmn) !
!***********************************!
!
IF (CSURF /='EXTE') RETURN
!-------------------------------------------------------------------------------
!
!*       2.     Orographic characteristics :
!               --------------------------
!
CALL IO_WRITE_FIELD(TPFILE,'SSO_ANIS', XSSO_ANISOTROPY)
CALL IO_WRITE_FIELD(TPFILE,'SSO_SLOPE',XSSO_SLOPE)
CALL IO_WRITE_FIELD(TPFILE,'SSO_DIR',  XSSO_DIRECTION)
CALL IO_WRITE_FIELD(TPFILE,'AVG_ZS',   XAVG_ZS)
CALL IO_WRITE_FIELD(TPFILE,'SIL_ZS',   XSIL_ZS)
CALL IO_WRITE_FIELD(TPFILE,'MAX_ZS',   XMAX_ZS)
CALL IO_WRITE_FIELD(TPFILE,'MIN_ZS',   XMIN_ZS)
CALL IO_WRITE_FIELD(TPFILE,'SSO_STDEV',XSSO_STDEV)
!
!-------------------------------------------------------------------------------
!
!*      3.     Dummy fields
!              ------------
!
CALL WRITE_DUMMY_GR_FIELD_n(TPFILE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHWRITE_ZS_DUMMY_n
