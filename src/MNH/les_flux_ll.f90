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
MODULE MODI_LES_FLUX_ll
!      ###################
!
INTERFACE LES_FLUX_ll
!
      SUBROUTINE LES_FLUX_ll_2D(PA, PB, OMASK,                         &
                                PAB_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PAB_FLUX_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_FLUX_ll_2D
!
      SUBROUTINE LES_FLUX_ll_3D(PA, PB, OMASK,                         &
                                PAB_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PB
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),            INTENT(OUT) :: PAB_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_FLUX_ll_3D
!
      SUBROUTINE LES_FLUX_ll_3DM(PA, PB, OMASK,                        &
                                PAB_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PB
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),            INTENT(OUT) :: PAB_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_FLUX_ll_3DM
!
END INTERFACE
!
END MODULE MODI_LES_FLUX_ll
!
!     ##################################################################
      SUBROUTINE LES_FLUX_ll_2D(PA, PB, OMASK,                         &
                                PAB_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     ##################################################################
!
!
!!****  *LES_FLUX_ll* computes the average of the flux between the
!!                    perturbations of two fields A and B, on all processors
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
!!      V. Masson        06/11/02 direct use of anomalies
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODI_LES_MEAN_ll
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PAB_FLUX_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: ZAB
INTEGER                                :: IAVG_PTS_ll
INTEGER                                :: IUND_PTS_ll
!-------------------------------------------------------------------------------
!
ZAB(:,:) = PA(:,:)*PB(:,:)
CALL LES_MEAN_ll (ZAB(:,:),OMASK(:,:),PAB_FLUX_ll,IAVG_PTS_ll,IUND_PTS_ll)
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_FLUX_ll_2D
!
!     ##################################################################
      SUBROUTINE LES_FLUX_ll_3D(PA, PB, OMASK,                         &
                                PAB_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     ##################################################################
!
!

!!****  *LES_FLUX_ll* computes the average of the flux between the
!!                    perturbations of two fields A and B, on all processors
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
USE MODI_LES_MEAN_ll
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PB
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),              INTENT(OUT) :: PAB_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,    INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,    INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: ZAB
INTEGER, DIMENSION(SIZE(PA,3))                    :: IAVG_PTS_ll
INTEGER, DIMENSION(SIZE(PA,3))                    :: IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
ZAB(:,:,:) = PA(:,:,:)*PB(:,:,:)
!
CALL LES_MEAN_ll (ZAB(:,:,:),OMASK(:,:,:),PAB_FLUX_ll(:),IAVG_PTS_ll(:),IUND_PTS_ll(:))
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_FLUX_ll_3D
!
!     ##################################################################
      SUBROUTINE LES_FLUX_ll_3DM(PA, PB, OMASK,                        &
                                PAB_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     ##################################################################
!
!

!!****  *LES_FLUX_ll* computes the average of the flux between the
!!                    perturbations of two fields A and B, on all processors
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
USE MODI_LES_MEAN_ll
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PB
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),              INTENT(OUT) :: PAB_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,    INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,    INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: ZAB
INTEGER, DIMENSION(SIZE(PA,3))                    :: IAVG_PTS_ll
INTEGER, DIMENSION(SIZE(PA,3))                    :: IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
ZAB(:,:,:) = PA(:,:,:)*PB(:,:,:)
!
CALL LES_MEAN_ll (ZAB(:,:,:),OMASK(:,:),PAB_FLUX_ll(:),IAVG_PTS_ll(:),IUND_PTS_ll(:))
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_FLUX_ll_3DM

