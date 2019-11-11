!MNH_LIC Copyright 2000-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################
MODULE MODI_WRITE_LES_BUDGET_n
!######################
!
INTERFACE
!
      SUBROUTINE  WRITE_LES_BUDGET_n(TPDIAFILE,HLES_AVG)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG ! flag to perform the averages
!                                        ! or normalizations
END SUBROUTINE WRITE_LES_BUDGET_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LES_BUDGET_n

!     ######################
      SUBROUTINE  WRITE_LES_BUDGET_n(TPDIAFILE,HLES_AVG)
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
!!      Original   07/02/00
!!                 06/11/02 (V. Masson) new LES budgets
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
USE MODD_CONF_n
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
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG ! flag to perform the averages
!                                        ! or normalizations
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
ZLES_BUDGET=XUNDEF
YSUBTITLE(:)=' '
!-------------------------------------------------------------------------------
!
!*      1.  total (resolved+subgrid) kinetic energy budget
!            ------------------------------------
!
YGROUP= 'BU_KE   '
ILES=0
ILES_STA=ILES
!
!* 1.1 tendency
!     --------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TEND'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_TEND)
!
!
!* 1.2 production by mean wind gradient
!     ---------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
!
ZLES_BUDGET(:,:,ILES)= - XLES_SUBGRID_WU (:,:,1) * XLES_MEAN_DUDZ(:,:,1)  &
                       - XLES_SUBGRID_WV (:,:,1) * XLES_MEAN_DVDZ(:,:,1)  &
                       - XLES_SUBGRID_W2 (:,:,1) * XLES_MEAN_DWDZ(:,:,1)
!
!* 1.3 production by wind gradient of resolved motions
!     ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_DP) - ZLES_BUDGET(:,:,2)
!
!
!
!* 1.4 advection
!      ---------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG ADVM'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_ADVM)
END IF
!
!
!* 1.5 forcing
!      -------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG FORC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_FORC)
END IF
!
!
!* 1.6 turbulent transport
!      -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_TR)
!
!
!* 1.7 advection by resolved motions
!      -----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG ADVR'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_ADVR)
!
!
!* 1.8 presso-correlations
!     ----------------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_PRES)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG PRES'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_PRES)
END IF
!
!
!* 1.9 thermal production
!      ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TP'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_TP)
!
!
!* 1.10 dissipation
!       -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_DISS)
!
!
!* 1.11 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_DIFF)
END IF
!
!* 1.12 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_RELA)
END IF
!
!* 1.13 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_NEST)
END IF
!
!
!* 1.14 other effects
!       -------------
!
IF ( ANY(XLES_BU_SBG_TKE(:,:,NLES_MISC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_SBG_TKE(:,:,NLES_MISC)
END IF
!
!
!* 1.15 residual of subgrid budget
!       --------------------------
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
!* 1.16 tendency
!       --------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TEND'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_TEND)
!
!
!* 1.17 advection
!       ---------
!
IF ( ANY(XLES_BU_RES_KE(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES ADV '
!
ZLES_BUDGET(:,:,ILES) =  XLES_BU_RES_Ke(:,:,NLES_ADVM)
END IF
!
!
!* 1.18 forcing
!       -------
!
IF ( ANY(XLES_BU_RES_KE(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES FORC'
!
ZLES_BUDGET(:,:,ILES) =  XLES_BU_RES_Ke(:,:,NLES_FORC)
END IF
!
!
!* 1.19 production by mean wind gradient
!       --------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Ke(:,:,NLES_DP)
!
!
!* 1.20 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_TR)
!
!
!
!* 1.21 presso-correlations
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES PRES'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_PRES)
!
!
!* 1.22 thermal production
!       ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_GRAV)
!
!
!* 1.23 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_VTURB) + XLES_BU_RES_KE(:,:,NLES_HTURB)
!
!* 1.24 effect of Coriolis (must be zero)
!       ------------------
!
IF ( ANY(XLES_BU_RES_KE(:,:,NLES_COR)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES CORI'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_COR)
END IF
!
!* 1.25 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_KE(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_DIFF)
END IF
!
!* 1.26 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_KE(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_RELA)
END IF
!
!* 1.27 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_KE(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_NEST)
END IF
!
!
!* 1.28 other effects
!       -------------
!
IF ( ANY( XLES_BU_RES_KE(:,:,NLES_MISC) &
         +XLES_BU_RES_KE(:,:,NLES_CURV) /= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_KE(:,:,NLES_MISC)  &
                      + XLES_BU_RES_KE(:,:,NLES_CURV)
END IF
!
!* 1.29 residual of resolved Ke budget
!       ------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RESI'
!
ZLES_BUDGET(:,:,ILES) = 0.
DO JLES=ILES_STA+1,ILES-1
  ZLES_BUDGET(:,:,ILES) = ZLES_BUDGET(:,:,ILES) - ZLES_BUDGET(:,:,JLES)
END DO
!
!* 1.30 writing
!       -------
!
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),"resolved KE budget"//YSUBTITLE(:ILES),"m2 s-3", &
                       ZLES_BUDGET(:,:,:ILES),HLES_AVG)
!
!-------------------------------------------------------------------------------
!
!
!*      2.  temperature variance budget
!           ---------------------------
!
YGROUP= 'BU_THL2 '
ILES=0
!
ILES_STA=ILES
!
!* 2.1 production by mean gradients
!      ----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
ILES_P1=ILES
!
ZLES_BUDGET(:,:,ILES)= - 2. * XLES_SUBGRID_WThl(:,:,1) * XLES_MEAN_dThldz(:,:,1)
!
!
!* 2.2 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_WThl2(:,:,1)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES) = - ( XLES_SUBGRID_WThl2 (JK+1,:,1)      &
                              -XLES_SUBGRID_WThl2 (JK-1,:,1)    ) &
                           / ( XLES_Z            (JK+1)           &
                              -XLES_Z            (JK-1)         )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
END IF
!
!* 2.3 production by resolved gradients
!      --------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
ILES_P2=ILES
!
ZLES_BUDGET(:,:,ILES)= - 2. * XLES_RES_ddxa_Thl_SBG_UaThl(:,:,1)  &
                       - ZLES_BUDGET(:,:,ILES_P1)
!
!
!* 2.4 dissipation
!      -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS    '
!
ZLES_BUDGET(:,:,ILES) =  XLES_SUBGRID_DISS_Thl2(:,:,1)
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
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_TEND)
!
!
!* 2.7 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_Thl2(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES ADV '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_ADVM)
END IF
!
!
!* 2.8 forcing
!      -------
!
IF ( ANY(XLES_BU_RES_Thl2(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES FORC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_FORC)
END IF
!
!
!* 2.9 production
!      ----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_DP)

!
!* 2.10 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_TR)
!
!
!* 2.11 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_VTURB) + XLES_BU_RES_Thl2(:,:,NLES_HTURB)
!
!* 2.12 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_Thl2(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_DIFF)
END IF
!
!* 2.13 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_Thl2(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_RELA)
END IF
!
!* 2.14 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_Thl2(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_NEST)
END IF
!
!* 2.15 other effects
!       -------------
!
IF ( ANY( XLES_BU_RES_Thl2(:,:,NLES_MISC) &
         +XLES_BU_RES_Thl2(:,:,NLES_RAD ) &
         +XLES_BU_RES_Thl2(:,:,NLES_MICR) &
         + XLES_BU_RES_Thl2(:,:,NLES_PREF) /= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_Thl2(:,:,NLES_MISC) &
                      + XLES_BU_RES_Thl2(:,:,NLES_RAD ) &
                      + XLES_BU_RES_Thl2(:,:,NLES_MICR) &
                      + XLES_BU_RES_Thl2(:,:,NLES_PREF)
END IF
!
!
!* 2.16 residual of resolved budget
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
!* 2.17 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG TEND'
!
ZLES_BUDGET(:,:,ILES) = 0.
IF (NLES_TIMES>2) THEN
  DO JK=1,NLES_K
    DO JT=2,NLES_TIMES-1
      ZLES_BUDGET(JK,JT,ILES) =- ( XLES_SUBGRID_Thl2 (JK  ,JT+1,1) &
                                 - XLES_SUBGRID_Thl2 (JK  ,JT-1,1))&
                                / (2.* XLES_TEMP_SAMPLING)
    END DO
    ZLES_BUDGET(JK,1         ,ILES) = ZLES_BUDGET(JK,2           ,ILES)
    ZLES_BUDGET(JK,NLES_TIMES,ILES) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES)
  END DO
END IF
!
!* 2.18 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)=  -XLES_MEAN_W(JK,:,1)                 &
                         * ( XLES_SUBGRID_Thl2(JK+1,:,1)        &
                            -XLES_SUBGRID_Thl2(JK-1,:,1)      ) &
                         / ( XLES_Z           (JK+1)            &
                            -XLES_Z           (JK-1)          )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 2.19 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - ( XLES_RES_W_SBG_Thl2   (JK+1,:,1)    &
                             -XLES_RES_W_SBG_Thl2   (JK-1,:,1)  ) &
                          / ( XLES_Z            (JK+1)            &
                             -XLES_Z            (JK-1)          )
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
CALL LES_DIACHRO_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),"thetal variance budget"//YSUBTITLE(:ILES),"K2 s-1", &
                       ZLES_BUDGET(:,:,:ILES),HLES_AVG)
!
!-------------------------------------------------------------------------------
!
!*      3.  temperature flux budget
!            ---------------------
!
YGROUP= 'BU_WTHL '
ILES=0
!
ILES_STA=ILES
!
!* 3.1 production by mean gradients
!     -----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
ILES_P1=ILES
!
ZLES_BUDGET(:,:,ILES) =  - XLES_SUBGRID_W2(:,:,1) * XLES_MEAN_dThldz(:,:,1)
!
!
!* 3.2 production by gradient of resolved motions
!     -------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
!
ZLES_BUDGET(:,:,ILES)=- XLES_RES_ddz_Thl_SBG_W2(:,:,1) &
                      - ZLES_BUDGET(:,:,ILES_P1)
!
!
!* 3.3 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_W2Thl(:,:,1)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES) = - ( XLES_SUBGRID_W2Thl (JK+1,:,1)       &
                              -XLES_SUBGRID_W2Thl (JK-1,:,1)     ) &
                           / ( XLES_Z            (JK+1)            &
                              -XLES_Z            (JK-1)          )
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
IF ( ANY(XLES_SUBGRID_ThlPz(:,:,1)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG PRES'
!
ZLES_BUDGET(:,:,ILES) =  XLES_SUBGRID_ThlPz(:,:,1)
END IF
!
!
!* 3.5 thermal production
!      ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TP  '
!
IF (LUSERV) THEN
  ZLES_BUDGET(:,:,ILEs) =  XG * XLES_SUBGRID_ThlThv(:,:,1)   &
                              / XLES_MEAN_Thv      (:,:,1)
ELSE
  ZLES_BUDGET(:,:,ILES) =  XG * XLES_SUBGRID_ThlThv(:,:,1)   &
                              / XLES_MEAN_Th       (:,:,1)
END IF
!
!* 3.6 residual of subgrid budget
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
!* 3.7 tendency
!      --------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TEND'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_TEND)
!
!* 3.8 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_WThl(:,:,NLES_ADVM)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES ADV '
!
ZLES_BUDGET(:,:,ILES) =  XLES_BU_RES_WThl(:,:,NLES_ADVM)
END IF
!
!* 3.9  forcing
!       -------
!
IF ( ANY(XLES_BU_RES_WThl(:,:,NLES_FORC)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES FORC'
!
ZLES_BUDGET(:,:,ILES) =  XLES_BU_RES_WThl(:,:,NLES_FORC)
END IF
!
!
!* 3.10 production by temperature gradient (and vertical wind gradient)
!       ----------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
ZLES_BUDGET(:,:,ILES) =  XLES_BU_RES_WThl(:,:,NLES_DP)
!
!* 3.11 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
ZLES_BUDGET(:,:,ILES) =  XLES_BU_RES_WThl(:,:,NLES_TR)
!
!
!* 3.12 presso-correlations
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES PRES'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_PRES)
!
!
!* 3.13 thermal production
!       ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TP  '
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_GRAV)
!
!
!* 3.14 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_VTURB) + XLES_BU_RES_WThl(:,:,NLES_HTURB)
!
!* 3.15 effect of Coriolis
!       ------------------
!
IF ( ANY(XLES_BU_RES_WThl(:,:,NLES_COR)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES CORI'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_COR)
END IF
!
!* 3.16 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_WThl(:,:,NLES_DIFF)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NUMD'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_DIFF)
END IF
!
!* 3.17 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_WThl(:,:,NLES_RELA)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RELA'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_RELA)
END IF
!
!* 3.18 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_WThl(:,:,NLES_NEST)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES NEST'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_NEST)
END IF
!
!* 3.19 other effects
!       -------------
!
IF ( ANY( XLES_BU_RES_WThl(:,:,NLES_MISC) &
         +XLES_BU_RES_WThl(:,:,NLES_RAD ) &
         +XLES_BU_RES_WThl(:,:,NLES_MICR) &
         +XLES_BU_RES_WThl(:,:,NLES_PREF) &
         +XLES_BU_RES_WThl(:,:,NLES_CURV) /= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' RES MISC'
!
ZLES_BUDGET(:,:,ILES) = XLES_BU_RES_WThl(:,:,NLES_MISC) &
                      + XLES_BU_RES_WThl(:,:,NLES_RAD ) &
                      + XLES_BU_RES_WThl(:,:,NLES_MICR) &
                      + XLES_BU_RES_WThl(:,:,NLES_PREF) &
                      + XLES_BU_RES_WThl(:,:,NLES_CURV)
END IF
!
!
!* 3.20 residual of resolved Wthl budget
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
!* 3.21 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG TEND'
!
ZLES_BUDGET(:,:,ILES) = 0.
IF (NLES_TIMES>2) THEN
  DO JK=1,NLES_K
    DO JT=2,NLES_TIMES-1
      ZLES_BUDGET(JK,JT,ILES) =- ( XLES_SUBGRID_WThl (JK  ,JT+1,1) &
                               - XLES_SUBGRID_WThl (JK  ,JT-1,1))&
                               / (2.* XLES_TEMP_SAMPLING)
    END DO
    ZLES_BUDGET(JK,1         ,ILES) = ZLES_BUDGET(JK,2           ,ILES)
    ZLES_BUDGET(JK,NLES_TIMES,ILES) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES)
  END DO
END IF
!
!
!
!* 3.22 neglected terms : advection for subgrid quantity
!       ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - XLES_MEAN_W(JK,:,1)                               &
                              * ( XLES_SUBGRID_WThl(JK+1,:,1)                 &
                                 -XLES_SUBGRID_WThl(JK-1,:,1)               ) &
                             /  ( XLES_Z          (JK+1)                      &
                                 -XLES_Z          (JK-1)                    )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 3.23 neglected terms : advection for subgrid quantity
!       ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JK=2,NLES_K-1
  ZLES_BUDGET(JK,:,ILES)= - ( XLES_RES_W_SBG_WThl(JK+1,:,1)   &
                             -XLES_RES_W_SBG_WThl(JK-1,:,1) ) &
                          / ( XLES_Z          (JK+1)          &
                             -XLES_Z          (JK-1)        )
END DO
!
ZLES_BUDGET(1     ,:,ILES) = ZLES_BUDGET(2       ,:,ILES)
ZLES_BUDGET(NLES_K,:,ILES) = ZLES_BUDGET(NLES_K-1,:,ILES)
!
!* 3.24 neglected terms : production by gradient of vertical velocity for subgrid quantity
!       ----------------------------------------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG DPGW'
!
ZLES_BUDGET(:,:,ILES)=- XLES_RES_ddxa_W_SBG_UaThl(:,:,1)
!
!
!* 3.25 neglected terms : production by hor. gradient of Thl for subgrid quantity
!       -------------------------------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG DPGT'
!
ZLES_BUDGET(:,:,ILES)=-XLES_RES_ddxa_Thl_SBG_UaW(:,:,1)       &
                      -ZLES_BUDGET(:,:,ILES_P1) -ZLES_BUDGET(:,:,ILES_P2)

!
!
!* 3.22 writing
!       -------
!
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),"thetal flux budget"//YSUBTITLE(:ILES),"m K s-2", &
                       ZLES_BUDGET(:,:,:ILES),HLES_AVG)
!
!
!-------------------------------------------------------------------------------
!
DEALLOCATE(ZLES_BUDGET)
DEALLOCATE(YSUBTITLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_LES_BUDGET_n 
