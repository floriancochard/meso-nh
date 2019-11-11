!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #################################
       MODULE MODI_LIMA_GRAUPEL_DEPOSITION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_GRAUPEL_DEPOSITION (LDCOMPUTE,                            &
                                       PRGT, PSSI, PLBDG, PAI, PCJ, PLSFACT, &
                                       P_TH_DEPG, P_RG_DEPG,                 &
                                       PA_TH, PA_RV, PA_RG                   )
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DEPG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
!!
END SUBROUTINE LIMA_GRAUPEL_DEPOSITION
END INTERFACE
END MODULE MODI_LIMA_GRAUPEL_DEPOSITION
!
!     ###########################################################################
      SUBROUTINE LIMA_GRAUPEL_DEPOSITION (LDCOMPUTE,                            &
                                          PRGT, PSSI, PLBDG, PAI, PCJ, PLSFACT, &
                                          P_TH_DEPG, P_RG_DEPG,                 &
                                          PA_TH, PA_RV, PA_RG                   )
!     ###########################################################################
!
!!    PURPOSE
!!    -------
!!      Deposition of water vapour on graupel
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
USE MODD_PARAM_LIMA,       ONLY : XRTMIN
USE MODD_PARAM_LIMA_MIXED, ONLY : X0DEPG, XEX0DEPG, X1DEPG, XEX1DEPG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DEPG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
!
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Deposition of vapour on graupel
!	        -------------------------------
!
P_TH_DEPG(:) = 0.0
P_RG_DEPG(:) = 0.0
WHERE ( (PRGT(:)>XRTMIN(6)) .AND. LDCOMPUTE(:) )
   P_RG_DEPG(:) = ( PSSI(:)/(PAI(:)) ) *                               &
        ( X0DEPG*PLBDG(:)**XEX0DEPG + X1DEPG*PCJ(:)*PLBDG(:)**XEX1DEPG )
   P_TH_DEPG(:) = P_RG_DEPG(:)*PLSFACT(:)
END WHERE
!
PA_RV(:) = PA_RV(:) - P_RG_DEPG(:)
PA_RG(:) = PA_RG(:) + P_RG_DEPG(:)
PA_TH(:) = PA_TH(:) + P_TH_DEPG(:)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_GRAUPEL_DEPOSITION
