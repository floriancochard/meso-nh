!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
SUBROUTINE SYSTEM_MNH(HCOMMAND)
!!
!!    PURPOSE
!!    -------
!!    This subroutine writes the 1 command line HCOMMAND
!!    in the file pipe_name and flushes the buffer.
!!    Modifications:
!!      Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!!
  USE MODE_IO_ll
!!
!*       0.    DECLARATIONS
!              ------------
!
  IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
  CHARACTER(LEN=*)    :: HCOMMAND
!
!*       0.2   Declaration of local variables
!              ------------------------------
#if defined(MNH_LINUX) || defined(MNH_SP4)
  CHARACTER(LEN=*),PARAMETER :: CFILE="file_for_xtransfer"
#else
#if !defined(MNH_SX5)
  CHARACTER(LEN=*),PARAMETER :: CFILE="file_for_fujitransfer"
#else
  CHARACTER(LEN=*),PARAMETER :: CFILE="file_for_nectransfer"
#endif
#endif
  INTEGER                    :: IUNIT
!
!
!
  IUNIT = -1
  OPEN(NEWUNIT=IUNIT,FILE=CFILE,ACCESS="sequential",FORM="formatted",POSITION="append")
  WRITE(IUNIT,*) HCOMMAND
  CLOSE(IUNIT)

END SUBROUTINE SYSTEM_MNH
