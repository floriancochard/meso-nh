!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!#########################
MODULE MODI_TRANSFER_FILE
!#########################
!
INTERFACE
      SUBROUTINE TRANSFER_FILE(HTRANS,HCPIO,HFILENAME)
!
CHARACTER(LEN=*), INTENT(IN) :: HTRANS    ! unix command for transfer
CHARACTER(LEN=*), INTENT(IN) :: HCPIO     ! CPIO option
CHARACTER(LEN=*), INTENT(IN) :: HFILENAME ! name of the file to transfer
!
END SUBROUTINE TRANSFER_FILE
END INTERFACE
END MODULE MODI_TRANSFER_FILE
!
!     ################################################
      SUBROUTINE TRANSFER_FILE(HTRANS,HCPIO,HFILENAME)
!     ################################################
!
!!****  *TRANSFER_FILE* - writes transfer.x command for a file in the pipe_name
!! 
!!    PURPOSE
!!    -------
!!    This subroutine writes the unix command line HTRANS HCPIO HFILENAME
!!    in the file pipe_name and flushes the buffer.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      Routine FLUSH : to flush the buffer
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/04/95
!!      modified by E.pesin 03/98 (for FUJITSU machine)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
CHARACTER(LEN=*), INTENT(IN) :: HTRANS    ! unix command for transfer
CHARACTER(LEN=*), INTENT(IN) :: HCPIO     ! CPIO option
CHARACTER(LEN=*), INTENT(IN) :: HFILENAME ! name of the file to transfer
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
CHARACTER(LEN=100)            :: YCOMMAND  ! command writen in pipe_name
!
!-------------------------------------------------------------------------------
!
WRITE(YCOMMAND,'(A," ",A," ",A," >> OUTPUT_TRANSFER 2>&1 &")') TRIM(HTRANS),TRIM(HCPIO),TRIM(HFILENAME)
PRINT *,'YCOMMAND =',YCOMMAND
!
print*, 'WARNING: routine TRANSFER_FILE DOES NOT WORK'
!!!!!CALL SYSTEM(YCOMMAND)
print*, 'WARNING: "CALL SYSTEM(YCOMMAND)" is not called in TRANSFER_FILE routine'
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRANSFER_FILE
