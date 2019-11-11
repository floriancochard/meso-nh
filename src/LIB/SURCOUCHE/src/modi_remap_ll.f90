!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!     ####################
      MODULE MODI_REMAP_ll 
!     ####################
!
INTERFACE
!
!!     ########################################################
       SUBROUTINE REMAP_2WAY_X_ll( PFIELDIN, PFIELDOUT, KINFO )
!!     ########################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
      END SUBROUTINE REMAP_2WAY_X_ll
!
!     ########################################################
      SUBROUTINE REMAP_X_2WAY_ll( PFIELDIN, PFIELDOUT, KINFO )
!     ########################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
      END SUBROUTINE REMAP_X_2WAY_ll
!
!     #####################################################
      SUBROUTINE REMAP_X_Y_ll( PFIELDIN, PFIELDOUT, KINFO )
!     #####################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
      END SUBROUTINE REMAP_X_Y_ll
!
!     #####################################################
      SUBROUTINE REMAP_Y_X_ll( PFIELDIN, PFIELDOUT, KINFO )
!     #####################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
      END SUBROUTINE REMAP_Y_X_ll
!
END INTERFACE
!
END MODULE MODI_REMAP_ll
