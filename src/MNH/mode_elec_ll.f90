!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###################
      MODULE MODE_ELEC_ll
!     ###################
!
!!    Purpose
!!    -------
!!
!!    Implicit arguments
!!    ------------------
!!
!!    Authors
!!    -------
!!    J. Escobar   * LA *
!!    C. Barthe   * LACy *
!!
!!    Modifications
!!    -------------
!!    Original 08/02/2010
!!
!------------------------------------------------------------------------
!
USE MODD_MPIF
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
!
IMPLICIT NONE
!
!include "mpif.h"
!
!
INTEGER, PARAMETER :: IFIRST_PROC = 0   ! 0/1 to increase numerotation of proc number 
INTEGER, PARAMETER :: MPI_PRECISION = MPI_DOUBLE_PRECISION
!
!
INTERFACE SUM_ELEC_ll
  MODULE PROCEDURE RSUM_ELEC_ll, ISUM_ELEC_ll, ISUM0_ELEC_ll, RSUM0_ELEC_ll
END INTERFACE SUM_ELEC_ll
!
INTERFACE MAX_ELEC_ll
  MODULE PROCEDURE RMAX0_ELEC_ll, IMAX0_ELEC_ll
END INTERFACE MAX_ELEC_ll
!
INTERFACE MIN_ELEC_ll
  MODULE PROCEDURE RMIN0_ELEC_ll, IMIN0_ELEC_ll
END INTERFACE MIN_ELEC_ll
!
INTERFACE EXTREMA_ELEC_ll
  MODULE PROCEDURE EXTREMA_ELEC_IJ_ll, EXTREMA_ELEC_IJK_ll
END INTERFACE EXTREMA_ELEC_ll
!
!------------------------------------------------------------------------
!
CONTAINS
!
!------------------------------------------------------------------------
!
!    ###############################
     SUBROUTINE MYPROC_ELEC_ll(KPROC)
!    ###############################
!
!!
!! Purpose : identify the proc rank
!!
!
IMPLICIT NONE
!
!*     0.     DECLARATIONS
!             ------------
!
!*     0.1    declaration of dummy arguments
!
INTEGER :: KPROC
!
!*     0.2    declaration of local variables
!
INTEGER :: INFO
!
!
CALL MPI_COMM_RANK (NMNH_COMM_WORLD, KPROC, INFO)
!
KPROC = KPROC + IFIRST_PROC

END SUBROUTINE MYPROC_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ###################################
     SUBROUTINE RSUM_ELEC_ll(PSUM_INOUT)
!    ###################################
!
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSUM_INOUT
!
!*     0.2    declaration of local variables
!
REAL, DIMENSION(SIZE(PSUM_INOUT)) :: ZTAB
INTEGER :: IDIM, INFO
!
!
IDIM = SIZE(PSUM_INOUT)
ZTAB = PSUM_INOUT
!
INFO = -1
!
! Sum(Proc)
CALL MPI_ALLREDUCE(ZTAB, PSUM_INOUT, IDIM, MPI_PRECISION, &
                   MPI_SUM, NMNH_COMM_WORLD, INFO)
!
END SUBROUTINE RSUM_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ######################################€####
     SUBROUTINE RMIN0_ELEC_ll(PMIN_INOUT, KPROC)
!    ###########################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
REAL,    INTENT(INOUT) :: PMIN_INOUT
INTEGER, INTENT(INOUT) :: KPROC
!
!*     0.2    declaration of local variables
!
REAL    :: ZTAB
INTEGER :: IDIM, INFO
INTEGER :: KPROC_LOCAL
!
!
IDIM = 1
ZTAB = PMIN_INOUT
!
INFO = -1
!
!*     1.1    max(Proc)
!
CALL MPI_ALLREDUCE(ZTAB, PMIN_INOUT, IDIM, MPI_PRECISION, &
                   MPI_MIN, NMNH_COMM_WORLD, INFO)
!
!*     1.2    find the proc number of the maximum
!
IF (PMIN_INOUT .EQ. ZTAB) THEN
  CALL MYPROC_ELEC_ll(KPROC_LOCAL)
ELSE
  KPROC_LOCAL = -1
ENDIF
!
!*     1.3    broadcast to all proc 
!
CALL MPI_ALLREDUCE(KPROC_LOCAL, KPROC,IDIM , MPI_INTEGER, &
                   MPI_MAX, NMNH_COMM_WORLD, INFO) 
!
END SUBROUTINE RMIN0_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ######################################€####
     SUBROUTINE RMAX0_ELEC_ll(PMAX_INOUT, KPROC)
!    ###########################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
REAL,    INTENT(INOUT) :: PMAX_INOUT
INTEGER, INTENT(INOUT) :: KPROC
!
!*     0.2    declaration of local variables
!
REAL    :: ZTAB
INTEGER :: IDIM, INFO
INTEGER :: KPROC_LOCAL
!
!
IDIM = 1
ZTAB = PMAX_INOUT
!
INFO = -1
!
!*     1.1    max(Proc)
!
CALL MPI_ALLREDUCE(ZTAB, PMAX_INOUT, IDIM, MPI_PRECISION, &
                   MPI_MAX, NMNH_COMM_WORLD, INFO)
!
!*     1.2    find the proc number of the maximum
!
IF (PMAX_INOUT .EQ. ZTAB) THEN
  CALL MYPROC_ELEC_ll(KPROC_LOCAL)
ELSE
  KPROC_LOCAL = -1
ENDIF
!
!*     1.3    broadcast to all proc 
!
CALL MPI_ALLREDUCE(KPROC_LOCAL, KPROC,IDIM , MPI_INTEGER, &
                   MPI_MAX, NMNH_COMM_WORLD, INFO) 
!
END SUBROUTINE RMAX0_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ###########################################
     SUBROUTINE IMIN0_ELEC_ll(KMIN_INOUT, KPROC)
!    ###########################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
INTEGER, INTENT(INOUT) :: KMIN_INOUT
INTEGER, INTENT(INOUT) :: KPROC
!
!*     0.2    declaration of local variables
!
INTEGER :: ITAB
INTEGER :: IDIM, INFO
INTEGER :: IPROC_LOCAL
!
!
IDIM = 1
ITAB = KMIN_INOUT
!
INFO = -1
!
!*     1.1    min(Proc)
!
CALL MPI_ALLREDUCE(ITAB, KMIN_INOUT, IDIM, MPI_INTEGER, &
                   MPI_MIN, NMNH_COMM_WORLD, INFO)
!
!*     1.2    find the proc number of the maximum
!
IF (KMIN_INOUT .EQ. ITAB) THEN
  CALL MYPROC_ELEC_ll(IPROC_LOCAL)
ELSE
  IPROC_LOCAL = -1
ENDIF
!
!*     1.3    broadcast to all proc 
!
CALL MPI_ALLREDUCE(IPROC_LOCAL, KPROC,IDIM, MPI_INTEGER, &
                   MPI_MAX, NMNH_COMM_WORLD, INFO) 
!
END SUBROUTINE IMIN0_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ###########################################
     SUBROUTINE IMAX0_ELEC_ll(KMAX_INOUT, KPROC)
!    ###########################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
INTEGER, INTENT(INOUT) :: KMAX_INOUT
INTEGER, INTENT(INOUT) :: KPROC
!
!*     0.2    declaration of local variables
!
INTEGER :: ITAB
INTEGER :: IDIM, INFO
INTEGER :: IPROC_LOCAL
!
!
IDIM = 1
ITAB = KMAX_INOUT
!
INFO = -1
!
!*     1.1    min(Proc)
!
CALL MPI_ALLREDUCE(ITAB, KMAX_INOUT, IDIM, MPI_INTEGER, &
                   MPI_MAX, NMNH_COMM_WORLD, INFO)
!
!*     1.2    find the proc number of the maximum
!
IF (KMAX_INOUT .EQ. ITAB) THEN
  CALL MYPROC_ELEC_ll(IPROC_LOCAL)
ELSE
  IPROC_LOCAL = -1
ENDIF
!
!*     1.3    brodcast to all proc 
!
CALL MPI_ALLREDUCE(IPROC_LOCAL, KPROC, IDIM, MPI_INTEGER, &
                   MPI_MAX, NMNH_COMM_WORLD, INFO) 
!
END SUBROUTINE IMAX0_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ###################################
     SUBROUTINE ISUM_ELEC_ll(KSUM_INOUT)
!    ###################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
INTEGER, DIMENSION(:), INTENT(INOUT) :: KSUM_INOUT
!
!*     0.2    declaration of local variables
!
INTEGER, DIMENSION(SIZE(KSUM_INOUT)) :: ITAB
INTEGER                              :: IDIM, INFO
!
!
IDIM = SIZE(KSUM_INOUT)
ITAB = KSUM_INOUT
!
INFO = -1
!
!*     1.1    sum(Proc)
!
CALL MPI_ALLREDUCE(ITAB, KSUM_INOUT, IDIM, MPI_INTEGER, &
                   MPI_SUM, NMNH_COMM_WORLD, INFO)
!
END SUBROUTINE ISUM_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ####################################
     SUBROUTINE ISUM0_ELEC_ll(KSUM_INOUT)
!    ####################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
INTEGER, INTENT(INOUT) :: KSUM_INOUT
!
!*     0.2    declaration of local variables
!
INTEGER :: ITAB
INTEGER :: IDIM, INFO
!
!
IDIM = 1
ITAB = KSUM_INOUT
!
INFO = -1
!
!*     1.1    sum(Proc)
!
CALL MPI_ALLREDUCE(ITAB, KSUM_INOUT, IDIM, MPI_INTEGER, &
                   MPI_SUM, NMNH_COMM_WORLD, INFO)
!
END SUBROUTINE ISUM0_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ####################################
     SUBROUTINE RSUM0_ELEC_ll(PSUM_INOUT)
!    ####################################
!
!! Purpose:
!
!*     0.     DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
REAL, INTENT(INOUT) :: PSUM_INOUT
!
!*     0.2    declaration of local variables
!
REAL    :: ZTAB
INTEGER :: IDIM, INFO
!
!
IDIM = 1
ZTAB = PSUM_INOUT
!
INFO = -1
!
!*     1.1    sum(Proc)
!
CALL MPI_ALLREDUCE(ZTAB, PSUM_INOUT, IDIM, MPI_PRECISION, &
                   MPI_SUM, NMNH_COMM_WORLD, INFO)
!
END SUBROUTINE RSUM0_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ###################################################################
     SUBROUTINE CORNER_BOX_ELEC_ll (KPROC_COORD,          &
                                    KI_CENTER, KJ_CENTER, &
                                    KI_RADIUS, KJ_RADIUS, &
                                    KIS, KIE, KJS, KJE,   &
                                    KI_WEST, KI_EST, KJ_SOUTH, KJ_NORTH)
!    ###################################################################
!
!! Purpose:
!
!*     0.     DECLARATION
!             -----------
!
USE MODD_PARAMETERS, ONLY : JPHEXT
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
INTEGER, INTENT(IN)  :: KPROC_COORD
INTEGER, INTENT(IN)  :: KI_CENTER, KJ_CENTER
INTEGER, INTENT(IN)  :: KI_RADIUS, KJ_RADIUS
INTEGER, INTENT(IN)  :: KIS, KIE, KJS, KJE
INTEGER, INTENT(OUT) :: KI_WEST, KI_EST, KJ_SOUTH, KJ_NORTH
!
!*      0.2    declaration of local variables
!
INTEGER :: ICENT_LOC,  JCENT_LOC
INTEGER :: ICENT_GLOB, JCENT_GLOB
!
INTEGER :: INW_LOC, INE_LOC, ISE_LOC, ISW_LOC
INTEGER :: JNW_LOC, JNE_LOC, JSE_LOC, JSW_LOC
!
INTEGER :: INW_GLOB, INE_GLOB, ISE_GLOB, ISW_GLOB
INTEGER :: JNW_GLOB, JNE_GLOB, JSE_GLOB, JSW_GLOB
!
INTEGER :: IS_GLOB, IE_GLOB, JS_GLOB, JE_GLOB
!
INTEGER :: IXOR, IYOR, IERR
!
!
ICENT_LOC = KI_CENTER + JPHEXT
JCENT_LOC = KJ_CENTER + JPHEXT
CALL GET_OR_ll('B', IXOR, IYOR)
ICENT_GLOB = IXOR + ICENT_LOC - 1
JCENT_GLOB = IYOR + JCENT_LOC - 1
!
! The proc with the center of the cell broadcast the global coord of the cell
!
CALL MPI_BCAST(ICENT_GLOB, 1, MPI_INTEGER, KPROC_COORD, NMNH_COMM_WORLD, IERR)
CALL MPI_BCAST(JCENT_GLOB, 1, MPI_INTEGER, KPROC_COORD, NMNH_COMM_WORLD, IERR)
!
IS_GLOB = KIS + IXOR -1
IE_GLOB = KIE + IXOR -1
JS_GLOB = KJS + IYOR -1
JE_GLOB = KJE + IYOR -1
!
INW_GLOB = ICENT_GLOB - KI_RADIUS 
JNW_GLOB = JCENT_GLOB + KJ_RADIUS 
!
INE_GLOB = ICENT_GLOB + KI_RADIUS 
JNE_GLOB = JCENT_GLOB + KJ_RADIUS 
!
ISE_GLOB = ICENT_GLOB + KI_RADIUS 
JSE_GLOB = JCENT_GLOB - KJ_RADIUS 
!
ISW_GLOB = ICENT_GLOB - KI_RADIUS 
JSW_GLOB = JCENT_GLOB - KJ_RADIUS 
!
INW_GLOB = MAX (INW_GLOB, IS_GLOB)  
JNW_GLOB = MIN (JNW_GLOB, JE_GLOB)
!
INE_GLOB = MIN (INE_GLOB, IE_GLOB)
JNE_GLOB = MIN (JNE_GLOB, JE_GLOB)
!
ISE_GLOB = MIN (ISE_GLOB, IE_GLOB)
JSE_GLOB = MAX (JSE_GLOB, JS_GLOB)
!
ISW_GLOB = MAX (ISW_GLOB, IS_GLOB)
JSW_GLOB = MAX (JSW_GLOB, JS_GLOB)
!
KI_WEST  = INW_GLOB - IXOR +1
KJ_NORTH = JNW_GLOB - IYOR +1
KI_EST   = INE_GLOB - IXOR +1
JNE_LOC  = JNE_GLOB - IYOR +1
ISE_LOC  = ISE_GLOB - IXOR +1
KJ_SOUTH = JSE_GLOB - IYOR +1
ISW_LOC  = ISW_GLOB - IXOR +1
JSW_LOC  = JSW_GLOB - IYOR +1
!
END SUBROUTINE CORNER_BOX_ELEC_ll
!
!------------------------------------------------------------------------
!
!    ########################################################
     SUBROUTINE EXTREMA_ELEC_IJ_ll(OMASK, KI_WEST, KI_EAST, &
                                       KJ_SOUTH, KJ_NORTH)
!    ########################################################
!
!
!! Purpose: find the global coordinates of the boundaries 
!! where OMASK = T
!
!
!*     0.     DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY : JPHEXT
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
LOGICAL, DIMENSION(:,:), INTENT(IN  )  :: OMASK
INTEGER,                 INTENT(INOUT) :: KI_WEST, KI_EAST, KJ_SOUTH, KJ_NORTH
!
!*      0.2    declaration of local variables
!
INTEGER :: II, IJ, IPROC_COORD
INTEGER :: I_GLOB, J_GLOB
INTEGER :: IXOR, IYOR, IERR
!
!
! origin's coordinates of each extended subdomain
CALL GET_OR_ll('B', IXOR, IYOR)
KI_WEST = 100000000
KI_EAST = -1
KJ_SOUTH = 100000000
KJ_NORTH = -1 
!
DO II = 1, SIZE(OMASK,1) 
  DO IJ = 1, SIZE(OMASK,2) 
    I_GLOB = II + IXOR - 1 !+ JPHEXT
    J_GLOB = IJ + IYOR - 1 !+ JPHEXT
    IF (OMASK(II,IJ)) THEN
      KI_WEST  = MIN(KI_WEST,  I_GLOB)
      KI_EAST  = MAX(KI_EAST,  I_GLOB)
      KJ_SOUTH = MIN(KJ_SOUTH, J_GLOB)
      KJ_NORTH = MAX(KJ_NORTH, J_GLOB)
    END IF
  END DO
END DO
!
CALL MIN_ELEC_ll(KI_WEST,  IPROC_COORD)
CALL MAX_ELEC_ll(KI_EAST,  IPROC_COORD)
CALL MIN_ELEC_ll(KJ_SOUTH, IPROC_COORD)
CALL MAX_ELEC_ll(KJ_NORTH, IPROC_COORD)
!
END SUBROUTINE EXTREMA_ELEC_IJ_ll
!
!------------------------------------------------------------------------
!
!    ###########################################################
     SUBROUTINE EXTREMA_ELEC_IJK_ll(OMASK, KI_WEST, KI_EAST,   &
                                           KJ_SOUTH, KJ_NORTH, &
                                           KK_BOTTOM, KK_TOP   )
!    ###########################################################
!
!
!! Purpose: find the global coordinates of the boundaries
!! where OMASK = T
!
!
!*     0.     DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY : JPHEXT
!
IMPLICIT NONE
!
!*     0.1    declaration of dummy arguments
!
LOGICAL, DIMENSION(:,:,:), INTENT(IN  )  :: OMASK
INTEGER,                   INTENT(INOUT) :: KI_WEST, KI_EAST, &
                                            KJ_SOUTH, KJ_NORTH, &
                                            KK_BOTTOM, KK_TOP
!
!*      0.2    declaration of local variables
!
INTEGER :: II, IJ, IK, IPROC_COORD
INTEGER :: I_GLOB, J_GLOB
INTEGER :: IXOR, IYOR, IERR
!
!
! origin's coordinates of each extended subdomain
CALL GET_OR_ll('B', IXOR, IYOR)
KI_WEST = 100000000
KI_EAST = -1
KJ_SOUTH = 100000000
KJ_NORTH = -1
KK_BOTTOM = 100000000
KK_TOP    = -1
!
DO II = 1, SIZE(OMASK,1)
  DO IJ = 1, SIZE(OMASK,2)
    DO IK = 1, SIZE(OMASK,3)
      I_GLOB = II + IXOR - 1 !+ JPHEXT
      J_GLOB = IJ + IYOR - 1 !+ JPHEXT
      IF (OMASK(II,IJ,IK)) THEN
        KI_WEST   = MIN(KI_WEST,   I_GLOB)
        KI_EAST   = MAX(KI_EAST,   I_GLOB)
        KJ_SOUTH  = MIN(KJ_SOUTH,  J_GLOB)
        KJ_NORTH  = MAX(KJ_NORTH,  J_GLOB)
        KK_BOTTOM = MIN(KK_BOTTOM, IK)
        KK_TOP    = MAX(KK_TOP,    IK)
      END IF
    END DO
  END DO
END DO
!
CALL MIN_ELEC_ll(KI_WEST,   IPROC_COORD)
CALL MAX_ELEC_ll(KI_EAST,   IPROC_COORD)
CALL MIN_ELEC_ll(KJ_SOUTH,  IPROC_COORD)
CALL MAX_ELEC_ll(KJ_NORTH,  IPROC_COORD)
CALL MIN_ELEC_ll(KK_BOTTOM, IPROC_COORD)
CALL MAX_ELEC_ll(KK_TOP,    IPROC_COORD)
!
END SUBROUTINE EXTREMA_ELEC_IJK_ll
!
!------------------------------------------------------------------------
!
END MODULE MODE_ELEC_ll
