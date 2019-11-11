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
!#####################
MODULE MODI_LES_RES_TR
!#####################
!
INTERFACE
!

      SUBROUTINE  LES_RES_TR(OUSERV,                                &
                             PMEAN_DUDZ, PMEAN_DVDZ, PMEAN_DWDZ,    &
                             PMEAN_DTHLDZ, PMEAN_DRTDZ, PMEAN_DSVDZ )

LOGICAL,              INTENT(IN) :: OUSERV ! flag to use water vapor
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DUDZ  ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DVDZ  ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DWDZ  ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DThlDZ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DRtDZ ! mean of du/dz
REAL, DIMENSION(:,:), INTENT(IN) :: PMEAN_DSvDZ ! mean of du/dz

END SUBROUTINE LES_RES_TR
!
END INTERFACE
!
END MODULE MODI_LES_RES_TR
!
!     ###############################################################
      SUBROUTINE  LES_RES_TR(OUSERV,                                &
                             PMEAN_DUDZ, PMEAN_DVDZ, PMEAN_DWDZ,    &
                             PMEAN_DTHLDZ, PMEAN_DRTDZ, PMEAN_DSVDZ )
!     ###############################################################
!
!
!!****  *LES_RES_TR* computes the resolved transport terms
!!
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
!!      Original         27/09/02
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
!
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_CONF
USE MODD_NSV, ONLY : NSV
!
USE MODI_LES_MEAN_ll
!
USE MODE_ll
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
LOGICAL,              INTENT(IN) :: OUSERV      ! flag to use water vapor
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DUDZ  ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DVDZ  ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DWDZ  ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DThlDZ! mean of du/dz
REAL, DIMENSION(:),   INTENT(IN) :: PMEAN_DRtDZ ! mean of du/dz
REAL, DIMENSION(:,:), INTENT(IN) :: PMEAN_DSvDZ ! mean of du/dz
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE ::  ZWORK_LES
!
INTEGER :: JSV      ! scalar variables counter
INTEGER :: IIU, IJU ! array sizes
INTEGER :: JI, JJ   ! loop counters
!
!-------------------------------------------------------------------------------
!
!*      1.   Initializations
!            ---------------
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!
ALLOCATE(ZWORK_LES(IIU,IJU,NLES_K))
!-------------------------------------------------------------------------------
!
!*      2.   Dynamical production by mean flow
!            ---------------------------------
!
!
!* Ke
!
DO JJ=1,IJU
  DO JI=1,IIU
    ZWORK_LES(JI,JJ,:) = - XU_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:) * PMEAN_DUDZ(:) &
                         - XV_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:) * PMEAN_DVDZ(:) &
                         - XW_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:) * PMEAN_DWDZ(:)
  END DO
END DO

CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_KE(:,NLES_DP))
!
!* Thl2
!
DO JJ=1,IJU
  DO JI=1,IIU
    ZWORK_LES(JI,JJ,:) = - 2. * XThl_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:) * PMEAN_DTHLDZ(:)
  END DO
END DO

CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_Thl2(:,NLES_DP))
!
!* WThl
!
DO JJ=1,IJU
  DO JI=1,IIU
    ZWORK_LES(JI,JJ,:) = - XW_ANOM(JI,JJ,:) * XThl_ANOM(JI,JJ,:) * PMEAN_DWDZ(:)   &
                         - XW_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:)   * PMEAN_DTHLDZ(:)
  END DO
END DO

CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_WThl(:,NLES_DP))
!
!* Rt2
!
IF (OUSERV) THEN
  DO JJ=1,IJU
    DO JI=1,IIU
      ZWORK_LES(JI,JJ,:) = - 2. * XRt_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:) * PMEAN_DRtDZ(:)
    END DO
  END DO

  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_Rt2(:,NLES_DP))
END IF
!
!* WRt
!
IF (OUSERV) THEN
  DO JJ=1,IJU
    DO JI=1,IIU
      ZWORK_LES(JI,JJ,:) = - XW_ANOM(JI,JJ,:) * XRt_ANOM(JI,JJ,:) * PMEAN_DWDZ(:)   &
                           - XW_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:)  * PMEAN_DRtDZ(:)
    END DO
  END DO

  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_WRt(:,NLES_DP))
END IF
!
!* ThlRt
!
IF (OUSERV) THEN
  DO JJ=1,IJU
    DO JI=1,IIU
      ZWORK_LES(JI,JJ,:) = - XW_ANOM(JI,JJ,:) * XRt_ANOM(JI,JJ,:)  * PMEAN_DTHLDZ(:)   &
                           - XW_ANOM(JI,JJ,:) * XThl_ANOM(JI,JJ,:) * PMEAN_DRtDZ(:)
    END DO
  END DO

  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_ThlRt(:,NLES_DP))
END IF
!
!* Sv2
!
DO JSV=1,NSV
  DO JJ=1,IJU
    DO JI=1,IIU
      ZWORK_LES(JI,JJ,:) = - 2. * XSv_ANOM(JI,JJ,:,JSV) * XW_ANOM(JI,JJ,:) * PMEAN_DSvDZ(:,JSV)
    END DO
  END DO
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_Sv2(:,NLES_DP,JSV))
END DO
!
!* WSv
!
DO JSV=1,NSV
  DO JJ=1,IJU
    DO JI=1,IIU
      ZWORK_LES(JI,JJ,:) = - XW_ANOM(JI,JJ,:) * XSv_ANOM(JI,JJ,:,JSV) * PMEAN_DWDZ(:)     &
                           - XW_ANOM(JI,JJ,:) * XW_ANOM(JI,JJ,:)      * PMEAN_DSvDZ(:,JSV)
    END DO
  END DO
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_WSv(:,NLES_DP,JSV))
END DO
!
!-------------------------------------------------------------------------------
!
!*      3.   Advection by mean flow
!            ----------------------
!
ZWORK_LES = 0.
!
!
CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_KE(:,NLES_ADVM))
CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_Thl2(:,NLES_ADVM))
CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_WThl(:,NLES_ADVM))
IF (OUSERV) THEN
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_Rt2(:,NLES_ADVM))
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_WRt(:,NLES_ADVM))
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_ThlRt(:,NLES_ADVM))
END IF
DO JSV=1,NSV
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_Sv2(:,NLES_ADVM,JSV))
  CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_RES_WSv(:,NLES_ADVM,JSV))
END DO
!
CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, X_LES_BU_SBG_Tke(:,NLES_ADVM))
!
!-------------------------------------------------------------------------------
!
!*      4.   Advection by resolved flow
!            --------------------------
!
X_LES_BU_SBG_Tke(:,NLES_ADVR) = X_LES_BU_SBG_Tke(:,NLES_TOTADV) &
                              - X_LES_BU_SBG_Tke(:,NLES_ADVM)
!
!-------------------------------------------------------------------------------
!
!*      5.   Turbulent transport
!            -------------------
!
X_LES_BU_RES_Ke(:,NLES_TR  ) = X_LES_BU_RES_Ke(:,NLES_TOTADV) &
                             - X_LES_BU_RES_Ke(:,NLES_ADVM)    &
                             - X_LES_BU_RES_Ke(:,NLES_DP)
!
X_LES_BU_RES_Thl2(:,NLES_TR  ) = X_LES_BU_RES_Thl2(:,NLES_TOTADV) &
                               - X_LES_BU_RES_Thl2(:,NLES_ADVM)    &
                               - X_LES_BU_RES_Thl2(:,NLES_DP)
!
X_LES_BU_RES_WThl(:,NLES_TR  ) = X_LES_BU_RES_WThl(:,NLES_TOTADV) &
                               - X_LES_BU_RES_WThl(:,NLES_ADVM)    &
                               - X_LES_BU_RES_WThl(:,NLES_DP)
!
IF (OUSERV) THEN
  X_LES_BU_RES_Rt2(:,NLES_TR  ) = X_LES_BU_RES_Rt2(:,NLES_TOTADV) &
                                - X_LES_BU_RES_Rt2(:,NLES_ADVM)    &
                                - X_LES_BU_RES_Rt2(:,NLES_DP)
  !
  X_LES_BU_RES_WRt(:,NLES_TR  ) = X_LES_BU_RES_WRt(:,NLES_TOTADV) &
                                - X_LES_BU_RES_WRt(:,NLES_ADVM)    &
                                - X_LES_BU_RES_WRt(:,NLES_DP)
  !
  X_LES_BU_RES_ThlRt(:,NLES_TR) = X_LES_BU_RES_ThlRt(:,NLES_TOTADV) &
                                - X_LES_BU_RES_ThlRt(:,NLES_ADVM)    &
                                - X_LES_BU_RES_ThlRt(:,NLES_DP)
END IF
!
DO JSV=1,NSV
  X_LES_BU_RES_Sv2(:,NLES_TR  ,JSV) = X_LES_BU_RES_Sv2(:,NLES_TOTADV,JSV) &
                                    - X_LES_BU_RES_Sv2(:,NLES_ADVM,JSV)    &
                                    - X_LES_BU_RES_Sv2(:,NLES_DP,JSV)
  !
  X_LES_BU_RES_WSv(:,NLES_TR  ,JSV) = X_LES_BU_RES_WSv(:,NLES_TOTADV,JSV) &
                                    - X_LES_BU_RES_WSv(:,NLES_ADVM,JSV)    &
                                    - X_LES_BU_RES_WSv(:,NLES_DP,JSV)
END DO
!
!--------------------------------------------------------------------------------
DEALLOCATE(ZWORK_LES)
!--------------------------------------------------------------------------------
!
END SUBROUTINE  LES_RES_TR
