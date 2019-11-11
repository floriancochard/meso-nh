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
!      ######################
MODULE MODI_LES_ANOMALY_FIELD
!      ######################
!
INTERFACE
!
SUBROUTINE LES_ANOMALY_FIELD(PF,PF_ANOM)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PF
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PF_ANOM
!
END SUBROUTINE LES_ANOMALY_FIELD
!
END INTERFACE
!
END MODULE MODI_LES_ANOMALY_FIELD
!
!          #############################
SUBROUTINE LES_ANOMALY_FIELD(PF,PF_ANOM)
!          #############################
!
!
!!****  *LES_ANOMALY_FIELD* - computes anomaly to mean field
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    March 2005
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODI_LES_VER_INT
USE MODI_LES_MEAN_ll
!
USE MODD_LES
!
!* 0.1    declaration of dummy arguments
!         ------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PF
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PF_ANOM
!
!* 0.2    declaration of local variables
!         ------------------------------
!
REAL, DIMENSION(SIZE(PF_ANOM,3)) :: ZMEAN
INTEGER :: JI, JJ

CALL LES_VER_INT(PF, PF_ANOM)
CALL LES_MEAN_ll(PF_ANOM, LLES_CURRENT_CART_MASK, ZMEAN  )
DO JJ=1,SIZE(PF_ANOM,2)
  DO JI=1,SIZE(PF_ANOM,1)
    PF_ANOM(JI,JJ,:) = PF_ANOM(JI,JJ,:) - ZMEAN(:)
  END DO
END DO

END SUBROUTINE LES_ANOMALY_FIELD
