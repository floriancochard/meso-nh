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
!      #######################
MODULE MODI_LES_BUDGET_TEND_n
!      #######################
!
INTERFACE LES_BUDGET_TEND_n
!
      SUBROUTINE LES_BUDGET_TEND_n
!
!
END SUBROUTINE LES_BUDGET_TEND_n
!
END INTERFACE
!
END MODULE MODI_LES_BUDGET_TEND_n

!     ##############################
      SUBROUTINE  LES_BUDGET_TEND_n
!     ##############################
!
!
!!****  *LES_BUDGET_TEND_n* computes tendencies of LES fluxes, variances
!!                          and covariances
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
!!      Original         07/02/00
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_NSV
USE MODD_LES
USE MODD_LES_n
USE MODD_FIELD_n
USE MODD_REF_n
USE MODD_DYN_n
USE MODD_CONF_n
USE MODD_LES_BUDGET
!
USE MODE_ll
!
USE MODI_LES_VER_INT
USE MODI_THL_RT_FROM_TH_R
USE MODI_LES_MEAN_ll
USE MODI_SHUMAN
USE MODI_LES_ANOMALY_FIELD
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!       0.2  declaration of local variables
!
!
INTEGER :: IIU, IJU                            ! hor. indices
INTEGER :: IKU                                 ! ver. index
!
INTEGER :: JSV                                 ! scalar variables counter
INTEGER :: JRR                                 ! moist variables counter
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZU_ANOM
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZV_ANOM
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZW_ANOM
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL_ANOM
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZRT_ANOM
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZSV_ANOM
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZRT
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZWORK_LES
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZTEND   ! tendency of subgrid Tke
!
INTEGER :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll  ! list of fields to exchange
!-------------------------------------------------------------------------------
!
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=SIZE(XRTHS,3)
!
!* get conservative variables
!  --------------------------
!
ALLOCATE(ZTHL        (IIU,IJU,IKU))
ALLOCATE(ZRT         (IIU,IJU,IKU))
CALL THL_RT_FROM_TH_R(LUSERV, LUSERC, LUSERR,             &
                      LUSERI, LUSERS, LUSERG, LUSERH,     &
                      XCURRENT_L_O_EXN_CP,                &
                      XTHT, XRT,                          &
                      ZTHL, ZRT                           )

!
!* anomalies at time-step 
!  ---------------------------------------------
!
ALLOCATE(ZU_ANOM    (IIU,IJU,NLES_K))
ALLOCATE(ZV_ANOM    (IIU,IJU,NLES_K))
ALLOCATE(ZW_ANOM    (IIU,IJU,NLES_K))
ALLOCATE(ZTHL_ANOM  (IIU,IJU,NLES_K))
ALLOCATE(ZRT_ANOM   (IIU,IJU,NLES_K))
ALLOCATE(ZSV_ANOM   (IIU,IJU,NLES_K,NSV))

CALL LES_ANOMALY_FIELD(MXF(XUT),ZU_ANOM)
CALL LES_ANOMALY_FIELD(MYF(XVT),ZV_ANOM)
CALL LES_ANOMALY_FIELD(MZF(1,IKU,1,XWT),ZW_ANOM)
CALL LES_ANOMALY_FIELD(ZTHL,ZTHL_ANOM)
CALL LES_ANOMALY_FIELD(ZRT,ZRT_ANOM)
DO JSV=1,NSV
  CALL LES_ANOMALY_FIELD(XSVT(:,:,:,JSV),ZSV_ANOM(:,:,:,JSV))
END DO
DEALLOCATE(ZTHL)
DEALLOCATE(ZRT)
!
ALLOCATE(ZWORK_LES(IIU,IJU,NLES_K))
!
!* KE budget
!
ZWORK_LES = (XU_ANOM * XU_ANOM - ZU_ANOM * ZU_ANOM) / XCURRENT_TSTEP &
           +(XV_ANOM * XV_ANOM - ZV_ANOM * ZV_ANOM) / XCURRENT_TSTEP &
           +(XW_ANOM * XW_ANOM - ZW_ANOM * ZW_ANOM) / XCURRENT_TSTEP
CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_KE(:,NLES_TEND))

!* WThl budget

ZWORK_LES = (XW_ANOM * XTHL_ANOM - ZW_ANOM * ZTHL_ANOM) / XCURRENT_TSTEP
CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_WThl(:,NLES_TEND))

!* Thl2 budget

ZWORK_LES = (XTHL_ANOM * XTHL_ANOM - ZTHL_ANOM * ZTHL_ANOM) / XCURRENT_TSTEP
CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_Thl2(:,NLES_TEND))


IF (LUSERV) THEN
!* ThlRt budget
  ZWORK_LES = (XRT_ANOM * XTHL_ANOM - ZRT_ANOM * ZTHL_ANOM) / XCURRENT_TSTEP
  CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_ThlRt(:,NLES_TEND))

!* Rt2 budget
  ZWORK_LES = (XRT_ANOM * XRT_ANOM - ZRT_ANOM * ZRT_ANOM) / XCURRENT_TSTEP
  CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_Rt2(:,NLES_TEND))

!* WRt budget
  ZWORK_LES = (XRT_ANOM * XW_ANOM - ZRT_ANOM * ZW_ANOM) / XCURRENT_TSTEP
  CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_WRt(:,NLES_TEND))

END IF

DO JSV=1,NSV
!* WSv budget
  ZWORK_LES = (XW_ANOM * XSV_ANOM(:,:,:,JSV) - ZW_ANOM * ZSV_ANOM(:,:,:,JSV)) / XCURRENT_TSTEP
  CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_WSv(:,NLES_TEND,JSV))

!* Sv2 budget
  ZWORK_LES = (XSV_ANOM(:,:,:,JSV)**2 - ZSV_ANOM(:,:,:,JSV)**2) / XCURRENT_TSTEP
  CALL LES_MEAN_ll(-ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_RES_Sv2(:,NLES_TEND,JSV))
END DO

!* Tke budget
ALLOCATE(ZTEND(IIU,IJU,IKU))
ZTEND(:,:,:) = -(XRTKES*XCURRENT_TSTEP/XRHODJ-XTKET) / XCURRENT_TSTEP
CALL LES_VER_INT( ZTEND, ZWORK_LES )
CALL LES_MEAN_ll(ZWORK_LES,LLES_CURRENT_CART_MASK,X_LES_BU_SBG_Tke(:,NLES_TEND))
DEALLOCATE(ZTEND)
!
!
DEALLOCATE(ZWORK_LES    )
!
DEALLOCATE(ZU_ANOM)
DEALLOCATE(ZV_ANOM)
DEALLOCATE(ZW_ANOM)
DEALLOCATE(ZTHL_ANOM)
DEALLOCATE(ZRT_ANOM)
DEALLOCATE(ZSV_ANOM)
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_BUDGET_TEND_n   
