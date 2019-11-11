!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #################################
       MODULE MODI_LIMA_DROPLETS_SELF_COLLECTION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION (LDCOMPUTE,                      &
                                             PRHODREF,                       &
                                             PCCT, PLBDC3,                   &
                                             P_CC_SELF,                      &
                                             PA_CC                           )
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC3  ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_SELF
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
!
END SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION
END INTERFACE
END MODULE MODI_LIMA_DROPLETS_SELF_COLLECTION
!
!     ######################################################################
      SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION (LDCOMPUTE,                      &
                                                PRHODREF,                       &
                                                PCCT, PLBDC3,                   &
                                                P_CC_SELF,                      &
                                                PA_CC                           )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the self-collection of cloud droplets rate
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XCTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XSELFC
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT     ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC3   ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_SELF
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PCCT)) :: ZW ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Cloud droplets self collection
!	        ------------------------------
!
!
P_CC_SELF(:)=0.
!
WHERE( PCCT(:)>XCTMIN(2) .AND. LDCOMPUTE(:) )
   ZW(:) = XSELFC*(PCCT(:)/PLBDC3(:))**2 * PRHODREF(:) ! analytical integration
   P_CC_SELF(:) = - ZW(:)
   PA_CC(:) = PA_CC(:) + P_CC_SELF(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION
