!MNH_LIC Copyright 2002-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################
MODULE MODI_WRITE_LES_RT_BUDGET_n
!######################
!
INTERFACE
!
      SUBROUTINE  WRITE_LES_RT_BUDGET_n(TPDIAFILE,HLES_AVG)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE ! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG  ! flag to perform the averages
!                                         ! or normalizations
END SUBROUTINE WRITE_LES_RT_BUDGET_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LES_RT_BUDGET_n
!
!     ######################
      SUBROUTINE  WRITE_LES_RT_BUDGET_n(TPDIAFILE,HLES_AVG)
!     ######################
!
!
!!****  *WRITE_LES_n* writes the LES final diagnostics for model _n 
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
!!      Original         06/11/02
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_LES_n
USE MODD_LES_BUDGET
!
USE MODE_ll
!
USE MODE_LES_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE ! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG  ! flag to perform the averages
!                                         ! or normalizations
!
!
!*      0.2  declaration of local variables
!
INTEGER :: ILES
INTEGER :: ILES_STA
INTEGER :: JLES
INTEGER :: ILES_P1, ILES_P2
!
INTEGER :: JK ! vertical loop counter
INTEGER :: JT ! temporal loop counter
!
CHARACTER(len=9), DIMENSION(:), ALLOCATABLE :: YSUBTITLE
CHARACTER(len=8)                            :: YGROUP
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLES_BUDGET
!
!-------------------------------------------------------------------------------
!
!*          Initializations
!            ---------------
!
ALLOCATE(ZLES_BUDGET(NLES_K,NLES_TIMES,40))
ALLOCATE(YSUBTITLE(40))
!
ZLES_BUDGET(:,:,:) = XUNDEF
!-------------------------------------------------------------------------------
!
!*      2.  total water variance budget
!           ---------------------------
!
!
YGROUP= 'BU_RT2  '
ILES=0
ILES_STA=ILES
!
!
!* 2.1 production by mean gradients
!      ----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
ILES_P1=ILES
!
ZLES_BUDGET(:,:,ILES)= - 2. * XLES_SUBGRID_WRt(:,:,1) * XLES_MEAN_dRtdz(:,:,1)
!
!
!* 2.3 production by resolved gradients
!      --------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
!
ZLES_BUDGET(:,:,ILES)= - 2. * XLES_RES_ddxa_Rt_SBG_UaRt(:,:,1)  &
                          - ZLES_BUDGET(:,:,ILES_P1)
!
!
!* 2.3 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_WRt2(:,:,1)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES) = - ( XLES_SUBGRID_WRt2 (JK+1,:,1)      &
                              -XLES_SUBGRID_WRt2 (JK-1,:,1)    ) &
                           / ( XLES_Z            (JK+1)          &
                              -XLES_Z            (JK-1)        )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
END IF
!
!* 2.4 dissipation
!      -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS    '
!
ZLES_BUDGET(:,:,ILES) =  XLES_SUBGRID_DISS_Rt2(:,:,1)
!
!
!* 2.5 residual of subgrid budget
!      --------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
ILES_STA=ILES
!
!* 2.6 tendency
!      --------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TEND'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_TEND)
!
!
!* 2.7 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_Rt2(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES ADV '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_ADVM)
END IF
!
!* 2.8 forcing
!      -------
!
IF ( ANY(XLES_BU_RES_Rt2(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES FORC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_FORC)
END IF
!
!
!* 2.9 production
!      ----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_DP)

!
!* 2.10 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_TR)
!
!
!* 2.11 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_VTURB) + XLES_BU_RES_Rt2(:,:,NLES_HTURB)
!
!
!* 2.11 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_Rt2(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_DIFF)
END IF
!
!* 2.11 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_Rt2(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_RELA)
END IF
!
!* 2.11 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_Rt2(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_NEST)
END IF
!
!* 2.11 other effects
!       -------------
!
IF ( ANY( XLES_BU_RES_Rt2(:,:,NLES_MISC) &
         +XLES_BU_RES_Rt2(:,:,NLES_MICR)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Rt2(:,:,NLES_MISC) &
                      + XLES_BU_RES_Rt2(:,:,NLES_MICR)
END IF
!
!
!* 2.12 residual of resolved budget
!       ---------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
!* 2.13 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG TEND'
!
ZLES_BUDGET(:,:,ILES) = 0.
IF (NLES_TIMES>2) THEN
  DO JK=1,NLES_K
    DO JT=2,NLES_TIMES-1
      ZLES_BUDGET(JK,JT,ILES) =- ( XLES_SUBGRID_Rt2 (JK  ,JT+1,1) &
                                 - XLES_SUBGRID_Rt2 (JK  ,JT-1,1))&
                              / (2.* XLES_TEMP_SAMPLING)
    END DO
    ZLES_BUDGET(JK,1         ,ILES) = ZLES_BUDGET(JK,2           ,ILES)
    ZLES_BUDGET(JK,NLES_TIMES,ILES) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES)
  END DO
END IF
!
!* 2.14 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)=  -XLES_MEAN_W(JK,:,1)                &
                         * ( XLES_SUBGRID_Rt2(JK+1,:,1)        &
                            -XLES_SUBGRID_Rt2(JK-1,:,1)      ) &
                         / ( XLES_Z          (JK+1)            &
                            -XLES_Z          (JK-1)          )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 2.15 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - ( XLES_RES_W_SBG_Rt2   (JK+1,:,1)    &
                             -XLES_RES_W_SBG_Rt2   (JK-1,:,1)  ) &
                          / ( XLES_Z           (JK+1)            &
                             -XLES_Z           (JK-1)          )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!
!* 2.16 writing
!       -------
!
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),"Rt variance budget"//YSUBTITLE(:ILES),"kg2 kg-2 s-1", &
                       ZLES_BUDGET(:,:,:ILES),HLES_AVG)
!
!
!-------------------------------------------------------------------------------
!
!*      3.  total water flux budget
!           -----------------------
!
!
YGROUP= 'BU_WRT  '
ILES=0
ILES_STA=ILES
!
!* 3.1 production by mean gradients
!     -----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
ILES_P1=ILES
!
ZLES_BUDGET(:,:,ILES) =  - XLES_SUBGRID_W2(:,:,1) * XLES_MEAN_dRtdz(:,:,1)
!
!
!* 3.2 production by gradient of resolved motions
!     -------------------------------------------
!
ILES=ILES+1
YSUBTITLE(2) = ' SBG DP R'
ILES_P2=ILES
!
ZLES_BUDGET(:,:,ILES)=- XLES_RES_ddz_Rt_SBG_W2(:,:,1) &
                      - ZLES_BUDGET(:,:,ILES_P1)
!
!
!
!* 3.3 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_W2Rt(:,:,1)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES) = - ( XLES_SUBGRID_W2Rt (JK+1,:,1)       &
                              -XLES_SUBGRID_W2Rt (JK-1,:,1)     ) &
                           / ( XLES_Z           (JK+1)            &
                              -XLES_Z           (JK-1)          )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
END IF
!
!
!* 3.4 presso-correlations
!      -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG PRES'
!
ZLES_BUDGET(:,:,ILES) =  XLES_SUBGRID_RtPz(:,:,1)
!
!
!* 3.5 thermal production
!      ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TP  '
!
ZLES_BUDGET(:,:,ILES) =  XG * XLES_SUBGRID_RtThv(:,:,1)   &
                            / XLES_MEAN_Thv     (:,:,1)
!
!
!* 3.6 dissipation
!      -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS    '
!
ZLES_BUDGET(:,:,ILES) = 0.
!
!
!* 3.7 residual of subgrid budget
!      --------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
ILES_STA=ILES
!
!* 3.8 tendency
!      --------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TEND'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_TEND)
!
!* 3.9 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_WRt(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES ADV '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_ADVM)
END IF
!
!* 3.10 forcing
!       -------
!
IF ( ANY(XLES_BU_RES_WRt(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES FORC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_FORC)
END IF
!
!* 3.11 production by temperature gradient (and vertical wind gradient)
!       ----------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_DP)
!
!* 3.12 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_TR)
!
!
!* 3.13 presso-correlations
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES PRES'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_PRES)
!
!
!* 3.14 thermal production
!       ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_GRAV)
!
!
!* 3.15 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_VTURB) + XLES_BU_RES_WRt(:,:,NLES_HTURB)
!
!* 3.15 effect of Coriolis
!       ------------------
!
IF ( ANY(XLES_BU_RES_WRt(:,:,NLES_COR)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES CORI'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_COR)
END IF
!
!* 3.15 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_WRt(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_DIFF)
END IF
!
!* 3.15 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_WRt(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_RELA)
END IF
!
!* 3.15 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_WRt(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_NEST)
END IF
!
!* 3.15 other effects
!       -------------
!
IF ( ANY( XLES_BU_RES_WRt(:,:,NLES_MISC) &
         +XLES_BU_RES_WRt(:,:,NLES_MICR) /= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WRt(:,:,NLES_MISC) &
                      + XLES_BU_RES_WRt(:,:,NLES_MICR)
END IF
!
!* 3.16 residual of resolved Wthl budget
!       --------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
!* 3.17 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TEND'
!
ZLES_BUDGET(:,:,ILES) = 0.
IF (NLES_TIMES>2) THEN
  DO JK=1,NLES_K
    DO JT=2,NLES_TIMES-1
      ZLES_BUDGET(JK,JT,ILES) =- ( XLES_SUBGRID_WRt (JK  ,JT+1,1) &
                                 - XLES_SUBGRID_WRt (JK  ,JT-1,1))&
                                 / (2.* XLES_TEMP_SAMPLING)
    END DO
    ZLES_BUDGET(JK,1         ,ILES) = ZLES_BUDGET(JK,2           ,ILES)
    ZLES_BUDGET(JK,NLES_TIMES,ILES) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES)
  END DO
END IF
!
!
!
!* 3.18 neglected terms : advection for subgrid quantity
!       ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - XLES_MEAN_W(JK,:,1)                              &
                              * ( XLES_SUBGRID_WRt(JK+1,:,1)                 &
                                 -XLES_SUBGRID_WRt(JK-1,:,1)               ) &
                             /  ( XLES_Z          (JK+1)                     &
                                 -XLES_Z          (JK-1)                   )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 3.19 neglected terms : advection for subgrid quantity
!       ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - ( XLES_RES_W_SBG_WRt(JK+1,:,1)   &
                             -XLES_RES_W_SBG_WRt(JK-1,:,1) ) &
                          / ( XLES_Z          (JK+1)         &
                             -XLES_Z          (JK-1)       )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 3.20 neglected terms : production by gradient of vertical velocity for subgrid quantity
!       ----------------------------------------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG DPGW'
!
ZLES_BUDGET(:,:,ILES)=- XLES_RES_ddxa_W_SBG_UaRt(:,:,1)
!
!
!* 3.21 neglected terms : production by hor. gradient of Thl for subgrid quantity
!       -------------------------------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG DPGT'
!
ZLES_BUDGET(:,:,ILES)=-XLES_RES_ddxa_Rt_SBG_UaW(:,:,1)       &
                      -ZLES_BUDGET(:,:,ILES_P1) -ZLES_BUDGET(:,:,ILES_P2)
!
!
!* 3.22 writing
!       -------
!
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),"Rt flux budget"//YSUBTITLE(:ILES),"m kg kg-1 s-2", &
                       ZLES_BUDGET(:,:,:ILES),HLES_AVG)
!
!
!-------------------------------------------------------------------------------
!
!*      4.  liquid potential temperature - total water covariance budget
!           ------------------------------------------------------------
!
!
YGROUP= 'BU_THLR '
ILES=0
ILES_STA=ILES
!
!
!* 2.1 production by mean gradients
!      ----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
ILES_P1=ILES
!
ZLES_BUDGET(:,:,ILES)=-XLES_SUBGRID_WRt (:,:,1) * XLES_MEAN_dThldz(:,:,1) &
                      -XLES_SUBGRID_WThl(:,:,1) * XLES_MEAN_dRtdz (:,:,1)
!
!
!* 2.3 production by resolved gradients
!      --------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
!
ZLES_BUDGET(:,:,ILES)= - XLES_RES_ddxa_Rt_SBG_UaThl(:,:,1)  &
                       - XLES_RES_ddxa_Thl_SBG_UaRt(:,:,1)  &
                       - ZLES_BUDGET(:,:,ILES_P1)
!
!
!* 2.3 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_WThlRt(:,:,1)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES) = - ( XLES_SUBGRID_WThlRt (JK+1,:,1)      &
                              -XLES_SUBGRID_WThlRt (JK-1,:,1)    ) &
                           / ( XLES_Z              (JK+1)          &
                              -XLES_Z              (JK-1)        )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
END IF
!
!* 2.4 dissipation
!      -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS    '
!
ZLES_BUDGET(:,:,ILES) =  XLES_SUBGRID_DISS_ThlRt(:,:,1)
!
!
!* 2.5 residual of subgrid budget
!      --------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
ILES_STA=ILES
!
!* 2.6 tendency
!      --------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TEND'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_TEND)
!
!
!* 2.7 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_ThlRt(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(7) = ' RES ADV '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_ADVM)
END IF
!
!* 2.8 forcing
!      -------
!
IF ( ANY(XLES_BU_RES_ThlRt(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES FORC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_FORC)
END IF
!
!* 2.9 production
!      ----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_DP)
!
!* 2.10 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_TR)
!
!
!* 2.11 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_VTURB) + XLES_BU_RES_ThlRt(:,:,NLES_HTURB)
!
!
!* 2.11 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_ThlRt(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_DIFF)
END IF
!
!* 2.11 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_ThlRt(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_RELA)
END IF
!
!* 2.11 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_ThlRt(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_NEST)
END IF
!
!* 2.11 miscellaneous effects
!       ---------------------
!
IF ( ANY( XLES_BU_RES_ThlRt(:,:,NLES_MISC) &
         +XLES_BU_RES_ThlRt(:,:,NLES_PREF) &
         +XLES_BU_RES_ThlRt(:,:,NLES_RAD ) &
         +XLES_BU_RES_ThlRt(:,:,NLES_MICR) /= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_ThlRt(:,:,NLES_MISC) &
                      + XLES_BU_RES_ThlRt(:,:,NLES_PREF) &
                      + XLES_BU_RES_ThlRt(:,:,NLES_RAD ) &
                      + XLES_BU_RES_ThlRt(:,:,NLES_MICR)
END IF
!
!
!* 2.12 residual of resolved budget
!       ---------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
!* 2.13 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG TEND'
!
ZLES_BUDGET(:,:,ILES) = 0.
IF (NLES_TIMES>2) THEN
  DO JK=1,NLES_K
    DO JT=2,NLES_TIMES-1
      ZLES_BUDGET(JK,JT,ILES) =- ( XLES_SUBGRID_ThlRt (JK  ,JT+1,1) &
                                 - XLES_SUBGRID_ThlRt (JK  ,JT-1,1))&
                                / (2.* XLES_TEMP_SAMPLING)
    END DO
    ZLES_BUDGET(JK,1         ,ILES) = ZLES_BUDGET(JK,2           ,ILES)
    ZLES_BUDGET(JK,NLES_TIMES,ILES) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES)
  END DO
END IF
!
!* 2.14 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)=  -XLES_MEAN_W(JK,:,1)                &
                         * ( XLES_SUBGRID_ThlRt(JK+1,:,1)      &
                            -XLES_SUBGRID_ThlRt(JK-1,:,1)    ) &
                         / ( XLES_Z          (JK+1)            &
                            -XLES_Z          (JK-1)          )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 2.15 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - ( XLES_RES_W_SBG_ThlRt (JK+1,:,1)    &
                             -XLES_RES_W_SBG_ThlRt (JK-1,:,1)  ) &
                          / ( XLES_Z           (JK+1)            &
                             -XLES_Z           (JK-1)          )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!
!* 2.16 writing
!       -------
!
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),"Thl-Rt covariance budget"//YSUBTITLE(:ILES), &
                                              "K kg kg-1 s-1",ZLES_BUDGET(:,:,:ILES),HLES_AVG)
!
!-------------------------------------------------------------------------------
!
DEALLOCATE(ZLES_BUDGET)
DEALLOCATE(YSUBTITLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_LES_RT_BUDGET_n 
