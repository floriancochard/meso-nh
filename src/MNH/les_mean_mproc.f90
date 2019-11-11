!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 les 2006/05/18 13:07:25
!-----------------------------------------------------------------
!      ###################
MODULE MODI_LES_MEAN_MPROC
!      ###################
!
INTERFACE LES_MEAN_MPROC
!
      SUBROUTINE LES_MEAN_MPROC_2D(PA_MEAN,    KAVG_PTS,   KUND_PTS,   &
                                               KAVG_PTS_ll,KUND_PTS_ll )
!
REAL,                    INTENT(INOUT):: PA_MEAN
INTEGER,                 INTENT(IN)   :: KAVG_PTS
INTEGER,                 INTENT(IN)   :: KUND_PTS
!
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_MEAN_MPROC_2D
!
      SUBROUTINE LES_MEAN_MPROC_3D(PA_MEAN,    KAVG_PTS,   KUND_PTS,   &
                                               KAVG_PTS_ll,KUND_PTS_ll )
!
REAL,    DIMENSION(:),   INTENT(INOUT)  :: PA_MEAN
INTEGER, DIMENSION(:),   INTENT(IN)     :: KAVG_PTS
INTEGER, DIMENSION(:),   INTENT(IN)     :: KUND_PTS
!
INTEGER, DIMENSION(:), OPTIONAL,   INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,   INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_MEAN_MPROC_3D
!
END INTERFACE
!
END MODULE MODI_LES_MEAN_MPROC
!
!     ##################################################################
      SUBROUTINE LES_MEAN_MPROC_2D(PA_MEAN,    KAVG_PTS,   KUND_PTS,   &
                                               KAVG_PTS_ll,KUND_PTS_ll )
!     ##################################################################
!
!
!!****  *LES_MEAN_MPROC* computes the average of one field from one processor
!!                       to all processors
!!
!!    PURPOSE
!!    -------
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
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,                    INTENT(INOUT)  :: PA_MEAN
INTEGER,                 INTENT(IN)     :: KAVG_PTS
INTEGER,                 INTENT(IN)     :: KUND_PTS
!
INTEGER, OPTIONAL,       INTENT(OUT)    :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT)    :: KUND_PTS_ll
!
!       0.2  declaration of local variables
!
INTEGER :: INFO_ll
!
REAL :: ZA_SUM   ! sum of the field 1) on current proc. 2) on all proc.
REAL :: ZAVG_PTS ! total aver. pts. 1) on current proc. 2) on all proc.
REAL :: ZUND_PTS ! total unde. pts. 1) on current proc. 2) on all proc.
!
INTEGER :: IAVG_PTS_ll
INTEGER :: IUND_PTS_ll
!-------------------------------------------------------------------------------
!
ZA_SUM   = KAVG_PTS * PA_MEAN
ZAVG_PTS = KAVG_PTS
ZUND_PTS = KUND_PTS
!
CALL REDUCESUM_ll(ZA_SUM,INFO_ll)
CALL REDUCESUM_ll(ZAVG_PTS,INFO_ll)
CALL REDUCESUM_ll(ZUND_PTS,INFO_ll)
!
IAVG_PTS_ll = NINT(ZAVG_PTS)
IUND_PTS_ll = NINT(ZUND_PTS)
!
IF ( IAVG_PTS_ll > 0 ) THEN
  PA_MEAN = ZA_SUM / IAVG_PTS_ll
ELSE
  PA_MEAN = XUNDEF
END IF
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_MPROC_2D
!
!     ##################################################################
      SUBROUTINE LES_MEAN_MPROC_3D(PA_MEAN,    KAVG_PTS,   KUND_PTS,   &
                                               KAVG_PTS_ll,KUND_PTS_ll )
!     ##################################################################
!
!
!!****  *LES_MEAN_MPROC* computes the average of one field from one processor
!!                       to all processors
!!
!!    PURPOSE
!!    -------
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
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:),   INTENT(INOUT)  :: PA_MEAN
INTEGER, DIMENSION(:),   INTENT(IN)     :: KAVG_PTS
INTEGER, DIMENSION(:),   INTENT(IN)     :: KUND_PTS
!
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
INTEGER :: INFO_ll
!
REAL, DIMENSION(SIZE(PA_MEAN,1)) :: ZA_SUM   ! sum of the field 1) on current proc. 2) on all proc.
REAL, DIMENSION(SIZE(PA_MEAN,1)) :: ZAVG_PTS ! total aver. pts. 1) on current proc. 2) on all proc.
REAL, DIMENSION(SIZE(PA_MEAN,1)) :: ZUND_PTS ! total unde. pts. 1) on current proc. 2) on all proc.
INTEGER, DIMENSION(SIZE(PA_MEAN,1)) :: IAVG_PTS_ll ! total aver. pts. 1) on current proc. 2) on all proc.
INTEGER, DIMENSION(SIZE(PA_MEAN,1)) :: IUND_PTS_ll ! total unde. pts. 1) on current proc. 2) on all proc.
!-------------------------------------------------------------------------------
!
ZA_SUM  (:) = KAVG_PTS(:) * PA_MEAN(:)
ZAVG_PTS(:) = KAVG_PTS(:)
ZUND_PTS(:) = KUND_PTS(:)
!
CALL REDUCESUM_ll(ZA_SUM(:),INFO_ll)
CALL REDUCESUM_ll(ZAVG_PTS(:),INFO_ll)
CALL REDUCESUM_ll(ZUND_PTS(:),INFO_ll)
!
IAVG_PTS_ll(:) = NINT(ZAVG_PTS(:))
IUND_PTS_ll(:) = NINT(ZUND_PTS(:))
!
WHERE ( IAVG_PTS_ll(:) > 0 ) 
  PA_MEAN(:) = ZA_SUM(:) / IAVG_PTS_ll(:)
ELSEWHERE
  PA_MEAN(:) = XUNDEF
END WHERE
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_MPROC_3D
