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
!    ######################
     MODULE MODI_LES_BUDGET
!    ######################
!
INTERFACE
!
      SUBROUTINE LES_BUDGET(PVARS,KBUDN,HBUVAR)

REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS    ! Source
INTEGER               , INTENT(IN) :: KBUDN    ! variable number
CHARACTER (LEN=*)    , INTENT(IN) :: HBUVAR   ! Identifier of the Budget of the
                                               ! variable that is considered

END SUBROUTINE LES_BUDGET

END INTERFACE

END MODULE MODI_LES_BUDGET
!
!     ####################################
      SUBROUTINE LES_BUDGET(PVARS,KBUDN,HBUVAR)
!     ####################################
!
!!****  *LES_BUDGET* - stores
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
!!      Original    September 19, 2002
!!      25/11/2016  Q.Rodier correction bug variance u'^2  v'^2  w'^2
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_NSV
!
USE MODI_SHUMAN
USE MODI_THL_RT_FROM_TH_R
USE MODI_LES_VER_INT
USE MODI_LES_MEAN_ll
!
USE MODE_ll
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!* 0.2    declaration of local variables
!         ------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS    ! Source
INTEGER               , INTENT(IN) :: KBUDN    ! variable number
CHARACTER (LEN=*)    , INTENT(IN) :: HBUVAR   ! Identifier of the Budget of the
                                               ! variable that is considered

!* 0.2    declaration of local variables
!         ------------------------------
!
INTEGER :: ILES_BU         ! number of process in LES budgets
!
INTEGER :: IIU, IJU, IKU   ! array dimensions
INTEGER :: JSV             ! scalar loop counter
REAL    :: ZTIME1, ZTIME2
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZTEND
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZWORK_LES  ! work array
!
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZANOM     ! field anomaly after process occured
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL_ANOM ! THL anomaly after process occured

REAL, DIMENSION(NLES_K) :: ZLES_PROF

INTEGER :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll  ! list of fields to exchange

INTEGER :: JK
!-------------------------------------------------------------------------------
!
CALL SECOND_MNH(ZTIME1)
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=SIZE(PVARS,3)
!
!
ALLOCATE(ZWORK_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZANOM       (IIU,IJU,NLES_K))
!
!*           LES budget term depending of the current physical process
!            ---------------------------------------------------------
!
CALL LES_BU_PROCESS(HBUVAR,ILES_BU)
!
!*           test on variable type
!            ---------------------
!
!
SELECT CASE (KBUDN)
!
!* u
!
  CASE(1)
    CALL LES_BUDGET_ANOMALY(PVARS,'X',ZANOM)
    !
    !* action in KE budget
    ZWORK_LES = 0.5*( ZANOM ** 2 - XU_ANOM ** 2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_KE(:,ILES_BU) = X_LES_BU_RES_KE(:,ILES_BU) + ZLES_PROF(:)
    !
    !* update fields
    XCURRENT_RUS = PVARS
    XU_ANOM = ZANOM
!
!* v
!
  CASE(2)
    CALL LES_BUDGET_ANOMALY(PVARS,'Y',ZANOM)
    !
    !* action in KE budget
    ZWORK_LES = 0.5*( ZANOM ** 2 - XV_ANOM ** 2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_KE(:,ILES_BU) = X_LES_BU_RES_KE(:,ILES_BU) + ZLES_PROF(:)
    !
    !* update fields
    XCURRENT_RVS = PVARS
    XV_ANOM = ZANOM
!
!* w
!
  CASE(3)
    CALL LES_BUDGET_ANOMALY(PVARS,'Z',ZANOM)
    !
    !* action in KE budget
    ZWORK_LES = 0.5*( ZANOM ** 2 - XW_ANOM ** 2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_KE(:,ILES_BU) = X_LES_BU_RES_KE(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in WTHL budget
    ZWORK_LES = ( ZANOM * XTHL_ANOM - XW_ANOM * XTHL_ANOM ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_WTHL(:,ILES_BU) = X_LES_BU_RES_WTHL(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in WRT budget
    IF (LCURRENT_USERV) THEN
      ZWORK_LES = ( ZANOM * XRT_ANOM - XW_ANOM * XRT_ANOM ) / XCURRENT_TSTEP
      CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
      X_LES_BU_RES_WRT(:,ILES_BU) = X_LES_BU_RES_WRT(:,ILES_BU) + ZLES_PROF(:)
    END IF
    !
    !* action in WSV budget
    DO JSV=1,NSV
      ZWORK_LES = ( ZANOM * XSV_ANOM(:,:,:,JSV) - XW_ANOM * XSV_ANOM(:,:,:,JSV)) / XCURRENT_TSTEP
      CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
      X_LES_BU_RES_WSV(:,ILES_BU,JSV) = X_LES_BU_RES_WSV(:,ILES_BU,JSV) + ZLES_PROF(:)
    END DO
    !
    !* update fields
    XCURRENT_RWS = PVARS
    XW_ANOM = ZANOM
!
!* Th
!
  CASE(4)
    XCURRENT_RTHLS = XCURRENT_RTHLS + PVARS - XCURRENT_RTHS
    CALL LES_BUDGET_ANOMALY(XCURRENT_RTHLS,'-',ZANOM)
    !
    !* action in WTHL budget
    ZWORK_LES = ( ZANOM * XW_ANOM - XW_ANOM * XTHL_ANOM ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_WTHL(:,ILES_BU) = X_LES_BU_RES_WTHL(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in THL2 budget
    ZWORK_LES = ( ZANOM ** 2 - XTHL_ANOM**2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_THL2(:,ILES_BU) = X_LES_BU_RES_THL2(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in THLRT budget
    IF (LCURRENT_USERV) THEN
      ZWORK_LES = ( ZANOM * XRT_ANOM - XRT_ANOM * XTHL_ANOM ) / &
                    XCURRENT_TSTEP
      CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
      X_LES_BU_RES_THLRT(:,ILES_BU) = X_LES_BU_RES_THLRT(:,ILES_BU) + ZLES_PROF(:)
    END IF
    !
    !* update fields
    XCURRENT_RTHS = PVARS
    XTHL_ANOM = ZANOM
!
!* Tke
!
  CASE(5)
    ALLOCATE(ZTEND(IIU,IJU,IKU))
    ZTEND(:,:,:) = (PVARS(:,:,:)-XCURRENT_RTKES(:,:,:)) / XCURRENT_RHODJ
    XCURRENT_RTKES = PVARS
    CALL LES_VER_INT( ZTEND, ZWORK_LES )
    DEALLOCATE(ZTEND)
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_SBG_Tke(:,ILES_BU) = X_LES_BU_SBG_Tke(:,ILES_BU) + ZLES_PROF(:)
!
!* Rv, Rr, Ri, Rs, Rg, Rh
!
  CASE(6,8,9,10,11,12)
    !* transformation into conservative variables: RT
    XCURRENT_RRTS = XCURRENT_RRTS + PVARS(:,:,:) - XCURRENT_RRS(:,:,:,KBUDN-5)
    CALL LES_BUDGET_ANOMALY(XCURRENT_RRTS,'-',ZANOM)
    !
    !* action in WRT budget
    ZWORK_LES = ( ZANOM * XW_ANOM - XW_ANOM * XRT_ANOM ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_WRT(:,ILES_BU) = X_LES_BU_RES_WRT(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in RT2 budget
    ZWORK_LES = ( ZANOM **2 - XRT_ANOM **2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_RT2(:,ILES_BU) = X_LES_BU_RES_RT2(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in THLRT budget
    ZWORK_LES = ( ZANOM * XTHL_ANOM - XTHL_ANOM * XRT_ANOM ) / &
                  XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_THLRT(:,ILES_BU) = X_LES_BU_RES_THLRT(:,ILES_BU) + ZLES_PROF(:)
    !
    !* update fields
    XCURRENT_RRS(:,:,:,KBUDN-5) = PVARS
    XRT_ANOM = ZANOM
!
!* Rc
!
  CASE(7)
    !* transformation into conservative variables: theta_l; RT
    XCURRENT_RRTS  = XCURRENT_RRTS  + PVARS(:,:,:) - XCURRENT_RRS(:,:,:,KBUDN-5)
    XCURRENT_RTHLS = XCURRENT_RTHLS - XCURRENT_L_O_EXN_CP &
                                    * (PVARS(:,:,:) - XCURRENT_RRS(:,:,:,KBUDN-5))

    !* anomaly of THL
    ALLOCATE(ZTHL_ANOM(IIU,IJU,NLES_K))
    CALL LES_BUDGET_ANOMALY(XCURRENT_RTHLS,'-',ZTHL_ANOM)
    !* anomaly of RT
    CALL LES_BUDGET_ANOMALY(XCURRENT_RRTS,'-',ZANOM)
    !
    !* action in WTHL budget
    ZWORK_LES = ( ZTHL_ANOM * XW_ANOM - XTHL_ANOM * XW_ANOM ) / &
                  XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_WTHL(:,ILES_BU) = X_LES_BU_RES_WTHL(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in THL2 budget
    ZWORK_LES = ( ZTHL_ANOM **2 - XTHL_ANOM **2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_THL2(:,ILES_BU) = X_LES_BU_RES_THL2(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in THLRT budget
    ZWORK_LES = ( ZANOM * ZTHL_ANOM - XRT_ANOM * XTHL_ANOM ) / &
                  XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_THLRT(:,ILES_BU) = X_LES_BU_RES_THLRT(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in WRT budget
    ZWORK_LES = ( ZANOM * XW_ANOM - XRT_ANOM * XW_ANOM ) /  &
                  XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_WRT(:,ILES_BU) = X_LES_BU_RES_WRT(:,ILES_BU) + ZLES_PROF(:)
    !
    !* action in RT2 budget
    ZWORK_LES = ( ZANOM **2 - XRT_ANOM **2 ) / XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_RT2(:,ILES_BU) = X_LES_BU_RES_RT2(:,ILES_BU) + ZLES_PROF(:)
    !
    !
    !* update fields
    XCURRENT_RRS(:,:,:,KBUDN-5) = PVARS
    XRT_ANOM = ZANOM
    XTHL_ANOM = ZTHL_ANOM
    DEALLOCATE(ZTHL_ANOM)
!
!* SV
!
  CASE(13:)
    CALL LES_BUDGET_ANOMALY(PVARS,'-',ZANOM)
    !
    !* action in WSV budget
    ZWORK_LES = ( ZANOM * XW_ANOM - XSV_ANOM(:,:,:,KBUDN-12) * XW_ANOM ) / &
                  XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_WSV(:,ILES_BU,KBUDN-12) = X_LES_BU_RES_WSV(:,ILES_BU,KBUDN-12) + ZLES_PROF(:)
    !
    !* action in SV2 budget
    ZWORK_LES = ( ZANOM **2 - XSV_ANOM(:,:,:,KBUDN-12) **2 ) / &
                  XCURRENT_TSTEP
    CALL LES_MEAN_ll( ZWORK_LES, LLES_CURRENT_CART_MASK, ZLES_PROF)
    X_LES_BU_RES_SV2(:,ILES_BU,KBUDN-12) = X_LES_BU_RES_SV2(:,ILES_BU,KBUDN-12) + ZLES_PROF(:)
    !
    !* update fields
    XCURRENT_RSVS(:,:,:,KBUDN-12) = PVARS
    XSV_ANOM(:,:,:,KBUDN-12) = ZANOM

END SELECT
!
!
!* deallocations
!
DEALLOCATE(ZWORK_LES)
DEALLOCATE(ZANOM    )
!
CALL SECOND_MNH(ZTIME2)
!
XTIME_LES_BU_PROCESS = XTIME_LES_BU_PROCESS + ZTIME2 - ZTIME1
XTIME_LES_BU = XTIME_LES_BU + ZTIME2 - ZTIME1
!
!---------------------------------------------------
CONTAINS
!---------------------------------------------------
SUBROUTINE LES_BU_PROCESS(HBU,KLES_BU)
CHARACTER (LEN=*), INTENT(IN)  :: HBU     ! Identifier of the Budget of the
                                          ! variable that is considered
INTEGER,           INTENT(OUT) :: KLES_BU ! LES budget identifier
!
IF (HBU(1:3)=='ADV') THEN
  KLES_BU = NLES_TOTADV
ELSE IF (HBU(1:3)=='REL') THEN
  KLES_BU = NLES_RELA
ELSE IF (HBU(1:5)=='VTURB') THEN
  KLES_BU = NLES_VTURB
ELSE IF (HBU(1:5)=='HTURB') THEN
  KLES_BU = NLES_HTURB
ELSE IF (HBU(1:4)=='GRAV') THEN
  KLES_BU = NLES_GRAV
ELSE IF (HBU(1:4)=='PRES') THEN
  KLES_BU = NLES_PRES
ELSE IF (HBU(1:4)=='PREF') THEN
  KLES_BU = NLES_PREF
ELSE IF (HBU(1:4)=='CURV') THEN
  KLES_BU = NLES_CURV
ELSE IF (HBU(1:3)=='COR') THEN
  KLES_BU = NLES_COR
ELSE IF (HBU(1:2)=='DP') THEN
  KLES_BU = NLES_DP
ELSE IF (HBU(1:2)=='TP') THEN
  KLES_BU = NLES_TP
ELSE IF (HBU(1:2)=='TR') THEN
  KLES_BU = NLES_TR
ELSE IF (HBU(1:4)=='DISS') THEN
  KLES_BU = NLES_DISS
ELSE IF (HBU(1:3)=='DIF') THEN
  KLES_BU = NLES_DIFF
ELSE IF (HBU(1:3)=='RAD') THEN
  KLES_BU = NLES_RAD
ELSE IF (HBU(1:4)=='NEST') THEN
  KLES_BU = NLES_NEST
ELSE
  KLES_BU = NLES_MISC
END IF
!
END SUBROUTINE LES_BU_PROCESS
!
!--------------------------------------------------------------------------------
!
SUBROUTINE LES_BUDGET_ANOMALY(PVARS,HGRID,PANOM)
!
USE MODI_LES_ANOMALY_FIELD
!
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PVARS
CHARACTER(LEN=1),       INTENT(IN)   :: HGRID
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PANOM
!
REAL, DIMENSION(SIZE(PVARS,1),SIZE(PVARS,2),SIZE(PVARS,3)) :: ZS
REAL, DIMENSION(SIZE(PVARS,1),SIZE(PVARS,2),SIZE(PVARS,3)) :: ZRHODJ
INTEGER :: IINFO_ll

    SELECT CASE (HGRID)
      CASE ('X')
        ZRHODJ(:,:,:) = MXM(XCURRENT_RHODJ)
        ZS(:,:,:) =  PVARS(:,:,:) / ZRHODJ * XCURRENT_TSTEP
      CASE ('Y')
        ZRHODJ(:,:,:) = MYM(XCURRENT_RHODJ)
        ZS(:,:,:) =  PVARS(:,:,:) / ZRHODJ * XCURRENT_TSTEP
      CASE ('Z')
        ZRHODJ(:,:,:) = MZM(1,IKU,1,XCURRENT_RHODJ)
        ZS(:,:,:) =  PVARS(:,:,:) / ZRHODJ * XCURRENT_TSTEP
      CASE DEFAULT
        ZRHODJ(:,:,:) =     XCURRENT_RHODJ
        ZS(:,:,:) =  PVARS(:,:,:) / ZRHODJ * XCURRENT_TSTEP
    END SELECT

    NULLIFY(TZFIELDS_ll)
    CALL ADD3DFIELD_ll(TZFIELDS_ll, ZS)
    CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
    CALL CLEANLIST_ll(TZFIELDS_ll)
    !
    SELECT CASE (HGRID)
      CASE ('X')
        ZS(:,:,:) = MXF(ZS)
      CASE ('Y')
        ZS(:,:,:) = MYF(ZS)
      CASE ('Z')
        ZS(:,:,:) = MZF(1,IKU,1,ZS)
    END SELECT

    CALL LES_ANOMALY_FIELD(ZS,PANOM)

END SUBROUTINE LES_BUDGET_ANOMALY
!--------------------------------------------------------------------------------
!
END SUBROUTINE LES_BUDGET
