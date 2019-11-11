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
MODULE MODI_LES_3RD_MOMENT_ll
!      ###################
!
INTERFACE LES_3RD_MOMENT_ll
!
      SUBROUTINE LES_3RD_MOMENT_ll_2D(PA, PB, PC,                            &
                                     OMASK,                                  &
                                     PABC_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB
REAL,    DIMENSION(:,:), INTENT(IN)  :: PC
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PABC_FLUX_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_3RD_MOMENT_ll_2D
!
      SUBROUTINE LES_3RD_MOMENT_ll_3D(PA, PB, PC,                           &
                                     OMASK,                                 &
                                     PABC_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PB
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PC
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),            INTENT(OUT) :: PABC_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_3RD_MOMENT_ll_3D
!
      SUBROUTINE LES_3RD_MOMENT_ll_3DM(PA, PB, PC,                          &
                                     OMASK,                                 &
                                     PABC_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PB
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PC
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),            INTENT(OUT) :: PABC_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KUND_PTS_ll
!
END SUBROUTINE LES_3RD_MOMENT_ll_3DM
!
END INTERFACE
!
END MODULE MODI_LES_3RD_MOMENT_ll
!
!     ##################################################################
      SUBROUTINE LES_3RD_MOMENT_ll_2D(PA, PB, PC,                            &
                                     OMASK,                                  &
                                     PABC_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll  )
!     ##################################################################
!
!
!!****  *LES_3RD_MOMENT_ll* computes the average of the flux between the
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
REAL,    DIMENSION(:,:), INTENT(IN)  :: PC
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PABC_FLUX_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KAVG_PTS_ll
INTEGER, OPTIONAL,       INTENT(OUT) :: KUND_PTS_ll
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: ZABC
INTEGER                                :: IAVG_PTS_ll
INTEGER                                :: IUND_PTS_ll
!-------------------------------------------------------------------------------
!
ZABC(:,:) = PA(:,:)*PB(:,:)*PC(:,:)
CALL LES_MEAN_ll (ZABC(:,:),OMASK(:,:),PABC_FLUX_ll,IAVG_PTS_ll,IUND_PTS_ll)
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_3RD_MOMENT_ll_2D
!
!     ##################################################################
      SUBROUTINE LES_3RD_MOMENT_ll_3D(PA, PB, PC,                           &
                                     OMASK,                                 &
                                     PABC_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll )
!     ##################################################################
!
!

!!****  *LES_3RD_MOMENT_ll* computes the average of the flux between the
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
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PC
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),            INTENT(OUT) :: PABC_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: ZABC
INTEGER, DIMENSION(SIZE(PA,3))                    :: IAVG_PTS_ll
INTEGER, DIMENSION(SIZE(PA,3))                    :: IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
ZABC(:,:,:) = PA(:,:,:)*PB(:,:,:)*PC(:,:,:)
!
CALL LES_MEAN_ll (ZABC(:,:,:),OMASK(:,:,:),PABC_FLUX_ll(:),IAVG_PTS_ll(:),IUND_PTS_ll(:))
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_3RD_MOMENT_ll_3D
!
!     ##################################################################
      SUBROUTINE LES_3RD_MOMENT_ll_3DM(PA, PB, PC,                          &
                                     OMASK,                                 &
                                     PABC_FLUX_ll, KAVG_PTS_ll, KUND_PTS_ll )
!     ##################################################################
!
!

!!****  *LES_3RD_MOMENT_ll* computes the average of the flux between the
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
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PC
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),            INTENT(OUT) :: PABC_FLUX_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KAVG_PTS_ll
INTEGER, DIMENSION(:), OPTIONAL,  INTENT(OUT) :: KUND_PTS_ll
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: ZABC
INTEGER, DIMENSION(SIZE(PA,3))                    :: IAVG_PTS_ll
INTEGER, DIMENSION(SIZE(PA,3))                    :: IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
ZABC(:,:,:) = PA(:,:,:)*PB(:,:,:)*PC(:,:,:)
!
CALL LES_MEAN_ll (ZABC(:,:,:),OMASK(:,:),PABC_FLUX_ll(:),IAVG_PTS_ll(:),IUND_PTS_ll(:))
!
IF (PRESENT(KAVG_PTS_ll)) KAVG_PTS_ll = IAVG_PTS_ll
IF (PRESENT(KUND_PTS_ll)) KUND_PTS_ll = IUND_PTS_ll
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_3RD_MOMENT_ll_3DM



