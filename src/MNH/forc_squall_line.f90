!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      %Z% Lib:%F%, Version:%I%, Date:%D%, Last modified:%E%
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_FORC_SQUALL_LINE
!     ############################
!
INTERFACE
!
SUBROUTINE FORC_SQUALL_LINE(PRTHS, PRHODJ, PDXHAT, PZHAT)
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRTHS  ! Source of TH
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ ! 
REAL, DIMENSION(:),     INTENT(IN)    :: PDXHAT ! Stretching in direction "x"
REAL, DIMENSION(:),     INTENT(IN)    :: PZHAT  ! Cartesian positions of "z"
!
END SUBROUTINE FORC_SQUALL_LINE    
!
END INTERFACE
!
END MODULE MODI_FORC_SQUALL_LINE    
!
!     #########################################################
      SUBROUTINE FORC_SQUALL_LINE(PRTHS, PRHODJ, PDXHAT, PZHAT)
!     #########################################################
!
!!
!!    PURPOSE
!!    -------
!!      This routine provides the intial disturbance of Theta: a cooling rate 
!!    to generate the cold pool of the squall line
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_BLANK
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!      J-P Pinty, Lab. Aerologie, 25/01/08 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_BLANK, ONLY : XDUMMY1,    & ! cooling rate (K/s)
                       XDUMMY2,    & ! vertical size of the disturbance
                       XDUMMY3,    & ! horizontal size of the disturbance
                       XDUMMY4,    & ! left border of the disturbance
                       XDUMMY5       ! duration (s) of the disturbance
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRTHS  ! Source of TH
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ ! 
REAL, DIMENSION(:),     INTENT(IN)    :: PDXHAT ! Stretching in direction "x"
REAL, DIMENSION(:),     INTENT(IN)    :: PZHAT  ! Cartesian positions of "z"
!
!*       0.2   Declarations of local variables :
!
INTEGER          :: JI,JK          ! Loop indexes along the x,y,z directions
INTEGER          :: JIBEG,JIEND    ! Loop indexes for the cooling area
!
!------------------------------------------------------------------------------
!
!
! SIZE OF THE COLD POOL
!
JIBEG = INT(XDUMMY4*FLOAT(SIZE(PDXHAT)))
JIEND = JIBEG + NINT(XDUMMY3/PDXHAT(JIBEG))
!
DO JK = 1+JPVEXT,SIZE(PZHAT)-JPVEXT
  IF (PZHAT(JK)<=XDUMMY2) THEN
    DO JI = JIBEG,JIEND
      PRTHS(JI,2,JK) = PRTHS(JI,2,JK) - PRHODJ(JI,2,JK)*XDUMMY1
    END DO
  END IF
END DO
!
RETURN
END SUBROUTINE FORC_SQUALL_LINE
