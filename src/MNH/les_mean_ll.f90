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
MODULE MODI_LES_MEAN_ll
!      ###################
!
INTERFACE LES_MEAN_ll
!
      SUBROUTINE LES_MEAN_ll_2D(PA, OMASK,                            &
                                PA_MEAN_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PA_MEAN_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_MEAN_ll_2D
!
      SUBROUTINE LES_MEAN_ll_3D(PA, OMASK,                            &
                                PA_MEAN_ll, KAVG_PTS_ll, KUND_PTS_ll  )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),           INTENT(OUT) :: PA_MEAN_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_MEAN_ll_3D
!
      SUBROUTINE LES_MEAN_ll_3DM(PA, OMASK,                            &
                                PA_MEAN_ll, KAVG_PTS_ll, KUND_PTS_ll  )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),           INTENT(OUT) :: PA_MEAN_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_MEAN_ll_3DM
!
END INTERFACE
!
END MODULE MODI_LES_MEAN_ll
!
!     #################################################################
      SUBROUTINE LES_MEAN_ll_2D(PA, OMASK,                            &
                                PA_MEAN_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     #################################################################
!
!
!!****  *LES_MEAN_ll* computes the average of one field on all processors
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
USE MODI_LES_MEAN_1PROC
USE MODI_LES_MEAN_MPROC
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PA_MEAN_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
!       0.2  declaration of local variables
!
INTEGER :: IAVG_PTS
INTEGER :: IUND_PTS
INTEGER :: IAVG_PTS_ll
INTEGER :: IUND_PTS_ll
!-------------------------------------------------------------------------------
!
CALL LES_MEAN_1PROC(PA, OMASK,                   &
                    PA_MEAN_ll, IAVG_PTS, IUND_PTS  )
!
CALL LES_MEAN_MPROC(PA_MEAN_ll,    IAVG_PTS,    IUND_PTS,   &
                                   IAVG_PTS_ll, IUND_PTS_ll )
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_ll_2D
!
!     #################################################################
      SUBROUTINE LES_MEAN_ll_3D(PA, OMASK,                            &
                                PA_MEAN_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     #################################################################
!
!
!!****  *LES_MEAN_ll* computes the average of one field on all processors
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
USE MODI_LES_MEAN_1PROC
USE MODI_LES_MEAN_MPROC
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),           INTENT(OUT) :: PA_MEAN_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IAVG_PTS
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IUND_PTS
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IAVG_PTS_ll
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IUND_PTS_ll
!-------------------------------------------------------------------------------
!
CALL LES_MEAN_1PROC(PA, OMASK,                      &
                    PA_MEAN_ll, IAVG_PTS, IUND_PTS  )
!
CALL LES_MEAN_MPROC(PA_MEAN_ll, IAVG_PTS,    IUND_PTS,   &
                                IAVG_PTS_ll, IUND_PTS_ll )
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_ll_3D
!
!     #################################################################
      SUBROUTINE LES_MEAN_ll_3DM(PA, OMASK,                           &
                                PA_MEAN_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     #################################################################
!
!
!!****  *LES_MEAN_ll* computes the average of one field on all processors
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
!!      Original         06/11/02
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODI_LES_MEAN_1PROC
USE MODI_LES_MEAN_MPROC
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),           INTENT(OUT) :: PA_MEAN_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IAVG_PTS
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IUND_PTS
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IAVG_PTS_ll
INTEGER, DIMENSION(SIZE(PA_MEAN_ll)) :: IUND_PTS_ll
!-------------------------------------------------------------------------------
!
CALL LES_MEAN_1PROC(PA, OMASK,                      &
                    PA_MEAN_ll, IAVG_PTS, IUND_PTS  )
!
CALL LES_MEAN_MPROC(PA_MEAN_ll, IAVG_PTS,    IUND_PTS,   &
                                IAVG_PTS_ll, IUND_PTS_ll )
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_ll_3DM

