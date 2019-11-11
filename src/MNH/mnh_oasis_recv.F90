!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########
MODULE MODI_MNH_OASIS_RECV
!     ##########
!
INTERFACE 
!
      SUBROUTINE MNH_OASIS_RECV(HPROGRAM,KI,KSW,PTIMEC,PTSTEP_SURF,   &
                             PZENITH,PSW_BANDS,          &
                             PTSRAD,PDIR_ALB,PSCA_ALB,PEMIS,PTSURF)
!
      CHARACTER(LEN=6),      INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
!
      INTEGER,               INTENT(IN)  :: KI          ! number of points on this proc
      INTEGER,               INTENT(IN)  :: KSW         ! number of short-wave spectral bands
      REAL,                  INTENT(IN)  :: PTIMEC      ! Cumulated run time step (s)
      REAL,                  INTENT(IN)  :: PTSTEP_SURF ! Surfex time step
!
      REAL, DIMENSION(:),    INTENT(IN)  :: PZENITH     ! zenithal angle       (radian from the vertical)
      REAL, DIMENSION(:),    INTENT(IN)  :: PSW_BANDS   ! mean wavelength of each shortwave band (m)
!
      REAL, DIMENSION(:),    INTENT(OUT) :: PTSRAD      ! radiative temperature                 (K)
      REAL, DIMENSION(:,:),  INTENT(OUT) :: PDIR_ALB    ! direct albedo for each spectral band  (-)
      REAL, DIMENSION(:,:),  INTENT(OUT) :: PSCA_ALB    ! diffuse albedo for each spectral band (-)
      REAL, DIMENSION(:),    INTENT(OUT) :: PEMIS       ! emissivity                            (-)
      REAL, DIMENSION(:),    INTENT(OUT) :: PTSURF      ! surface effective temperature         (K)
!
      END SUBROUTINE MNH_OASIS_RECV
!
END INTERFACE
!
END MODULE MODI_MNH_OASIS_RECV
!
!     ####################################################################
SUBROUTINE MNH_OASIS_RECV   (HPROGRAM,KI,KSW,PTIMEC,PTSTEP_SURF,   &
                             PZENITH,PSW_BANDS,          &
                             PTSRAD,PDIR_ALB,PSCA_ALB,PEMIS,PTSURF )
!#############################################
!
!!****  *MNH_OASIS_RECV*
!!
!!    PURPOSE
!!    -------
!!    Meso-NH driver that receive coupling fields from oasis
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	J. Pianezze   *LPO*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2014
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_SFX_OASIS,  ONLY : XTSTEP_CPL_LAND, &
                            XTSTEP_CPL_SEA,  &
                            XTSTEP_CPL_WAVE,  &
                            LWATER
!
USE MODD_SFX_OASIS,  ONLY : LCPL_LAND,         &
                            LCPL_GW,LCPL_FLOOD,&
                            LCPL_SEA,          &
                            LCPL_SEAICE,       &
                            LCPL_WAVE
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_MNH_SURFEX_n
!
USE MODI_GET_LUOUT
USE MODI_SFX_OASIS_RECV
USE MODI_PUT_SFX_LAND
USE MODI_PUT_SFX_SEA
USE MODI_PUT_SFX_WAVE
USE MODI_UPDATE_ESM_SURF_ATM_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),       INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
!
INTEGER,                INTENT(IN)  :: KI          ! number of points on this proc
INTEGER,                INTENT(IN)  :: KSW         ! number of short-wave spectral bands
REAL,                   INTENT(IN)  :: PTIMEC      ! Cumulated run time step (s)
REAL,                   INTENT(IN)  :: PTSTEP_SURF ! Surfex time step
!
REAL, DIMENSION(:),    INTENT(IN)  :: PZENITH   ! zenithal angle       (radian from the vertical)
REAL, DIMENSION(:),    INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
!
REAL, DIMENSION(:),    INTENT(OUT) :: PTSRAD    ! radiative temperature                 (K)
REAL, DIMENSION(:,:),  INTENT(OUT) :: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),  INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:),    INTENT(OUT) :: PEMIS     ! emissivity                            (-)
REAL, DIMENSION(:),    INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI) :: ZLAND_WTD     ! Land water table depth (m)
REAL, DIMENSION(KI) :: ZLAND_FWTD    ! Land grid-cell fraction of water table rise (-)
REAL, DIMENSION(KI) :: ZLAND_FFLOOD  ! Land Floodplains fraction (-)
REAL, DIMENSION(KI) :: ZLAND_PIFLOOD ! Land Potential flood infiltration(kg/m2/s)
REAL, DIMENSION(KI) :: ZSEA_SST      ! Sea surface temperature (K)
REAL, DIMENSION(KI) :: ZSEA_UCU      ! Sea u-current stress (Pa)
REAL, DIMENSION(KI) :: ZSEA_VCU      ! Sea v-current stress (Pa)
REAL, DIMENSION(KI) :: ZSEAICE_SIT   ! Sea-ice Temperature (K)
REAL, DIMENSION(KI) :: ZSEAICE_CVR   ! Sea-ice cover (-)
REAL, DIMENSION(KI) :: ZSEAICE_ALB   ! Sea-ice albedo (-)
REAL, DIMENSION(KI) :: ZWAVE_CHA     ! Charnock coefficient (-)
REAL, DIMENSION(KI) :: ZWAVE_UCU     ! u-current velocity (m/s)
REAL, DIMENSION(KI) :: ZWAVE_VCU     ! v-current velocity (m/s)
REAL, DIMENSION(KI) :: ZWAVE_HS      ! Significant wave height (m)
REAL, DIMENSION(KI) :: ZWAVE_TP      ! Peak period (s)
!
INTEGER             :: ILUOUT
REAL                :: ZTIME_CPL
!
LOGICAL             :: GRECV_LAND
LOGICAL             :: GRECV_FLOOD
LOGICAL             :: GRECV_SEA
LOGICAL             :: GRECV_WAVE
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       1.     init coupling fields:
!               ----------------------------------
!
ZTIME_CPL = PTIMEC-PTSTEP_SURF
!
GRECV_LAND =(LCPL_LAND .AND. MOD(ZTIME_CPL,XTSTEP_CPL_LAND)==0.0)
GRECV_SEA  =(LCPL_SEA  .AND. MOD(ZTIME_CPL,XTSTEP_CPL_SEA )==0.0)
GRECV_WAVE =(LCPL_WAVE .AND. MOD(ZTIME_CPL,XTSTEP_CPL_WAVE)==0.0)
!
IF(GRECV_LAND)THEN
  ZLAND_WTD    (:) = XUNDEF
  ZLAND_FWTD   (:) = XUNDEF
  ZLAND_FFLOOD (:) = XUNDEF
  ZLAND_PIFLOOD(:) = XUNDEF
ENDIF
!
IF(GRECV_SEA)THEN
  ZSEA_SST   (:) = XUNDEF
  ZSEA_UCU   (:) = XUNDEF
  ZSEA_VCU   (:) = XUNDEF
  ZSEAICE_SIT(:) = XUNDEF
  ZSEAICE_CVR(:) = XUNDEF
  ZSEAICE_ALB(:) = XUNDEF
ENDIF
!
IF(GRECV_WAVE)THEN
  ZWAVE_CHA(:) = XUNDEF
  ZWAVE_UCU(:)  = XUNDEF
  ZWAVE_VCU(:)  = XUNDEF
  ZWAVE_HS(:)  = XUNDEF
  ZWAVE_TP(:)  = XUNDEF
ENDIF
!
!
!*       2.     Receive fields to other models proc by proc:
!               --------------------------------------------
!
CALL SFX_OASIS_RECV(HPROGRAM,KI,KSW,ZTIME_CPL,        &
                    GRECV_LAND, GRECV_SEA, GRECV_WAVE, &
                    ZLAND_WTD    (:),ZLAND_FWTD   (:), &
                    ZLAND_FFLOOD (:),ZLAND_PIFLOOD(:), &
                    ZSEA_SST     (:),ZSEA_UCU     (:), &
                    ZSEA_VCU     (:),ZSEAICE_SIT  (:), &
                    ZSEAICE_CVR  (:),ZSEAICE_ALB  (:), &
                    ZWAVE_CHA    (:),ZWAVE_UCU    (:), &
                    ZWAVE_VCU    (:),ZWAVE_HS     (:), &
                    ZWAVE_TP     (:)                    )
!
!
!*       3.     Put definitions for exchange of coupling fields :
!               -------------------------------------------------
!
!-------------------------------------------------------------------------------
! Put variable over land tile
!-------------------------------------------------------------------------------
!
IF(GRECV_LAND)THEN
  CALL PUT_SFX_LAND(YSURF_CUR%IM%O, YSURF_CUR%IM%S, YSURF_CUR%IM%K, &
                    YSURF_CUR%IM%NK, YSURF_CUR%IM%NP, YSURF_CUR%U,      &
                    ILUOUT,LCPL_GW,LCPL_FLOOD,        &
                    ZLAND_WTD   (:),ZLAND_FWTD   (:), &
                    ZLAND_FFLOOD(:),ZLAND_PIFLOOD(:)  )        
ENDIF
!
!-------------------------------------------------------------------------------
! Put variable over sea and/or water tile
!-------------------------------------------------------------------------------
!
IF(GRECV_SEA)THEN
  CALL PUT_SFX_SEA(YSURF_CUR%SM%S, YSURF_CUR%U, YSURF_CUR%WM%W, &
                   ILUOUT,LCPL_SEAICE,LWATER,                   &
                   ZSEA_SST   (:),ZSEA_UCU   (:),               &
                   ZSEA_VCU   (:),ZSEAICE_SIT(:),               &
                   ZSEAICE_CVR(:),ZSEAICE_ALB(:)                )
ENDIF
!
!-------------------------------------------------------------------------------
! Put variable over sea and/or water tile for waves
!-------------------------------------------------------------------------------
!
IF(GRECV_WAVE)THEN
  CALL PUT_SFX_WAVE(YSURF_CUR%SM%S, YSURF_CUR%U,               &
                    ILUOUT,ZWAVE_CHA(:),ZWAVE_UCU(:),          &
                    ZWAVE_VCU(:),ZWAVE_HS(:),ZWAVE_TP(:)       )
ENDIF
!
!-------------------------------------------------------------------------------
! Update radiative properties at time t+1 for radiative scheme
!-------------------------------------------------------------------------------
!
GRECV_FLOOD=(GRECV_LAND.AND.LCPL_FLOOD)
!
IF(GRECV_SEA.OR.GRECV_FLOOD)THEN     
  CALL UPDATE_ESM_SURF_ATM_n(YSURF_CUR%FM%F, YSURF_CUR%IM, YSURF_CUR%SM%S, &
                             YSURF_CUR%U, YSURF_CUR%WM%W, &
                             HPROGRAM, KI, KSW, PZENITH(:), PSW_BANDS, &
                             PTSRAD(:), PDIR_ALB(:,:),     &
                             PSCA_ALB(:,:), PEMIS(:),      &
                             PTSURF(:)                                       ) 
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNH_OASIS_RECV
