!MNH_LIC Copyright 2002-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################
MODULE MODI_WRITE_LES_SV_BUDGET_n
!######################
!
INTERFACE
!
      SUBROUTINE  WRITE_LES_SV_BUDGET_n(TPDIAFILE,HLES_AVG)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE ! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG  ! flag to perform the averages
!                                         ! or normalizations
END SUBROUTINE WRITE_LES_SV_BUDGET_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LES_SV_BUDGET_n

!     ######################
      SUBROUTINE  WRITE_LES_SV_BUDGET_n(TPDIAFILE,HLES_AVG)
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
USE MODD_CONF_n
USE MODD_LES_BUDGET
USE MODD_NSV
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
!*      0.2  declaration of local variables
!
INTEGER :: ILES
INTEGER :: ILES_STA
INTEGER :: JLES
INTEGER :: ILES_P1, ILES_P2
!
INTEGER :: JK ! vertical loop counter
INTEGER :: JT ! temporal loop counter
INTEGER :: JSV! scalar loop counter
INTEGER :: JP ! process loop counter
!
CHARACTER(len=9), DIMENSION(:), ALLOCATABLE :: YSUBTITLE
CHARACTER(len=8)                            :: YGROUP
CHARACTER(len=20)                           :: YTITLE
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZLES_BUDGET
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSV_BUDGET
!
!-------------------------------------------------------------------------------
!
!*          Initializations
!            ---------------
!
ALLOCATE(ZLES_BUDGET(NLES_K,NLES_TIMES,50,NSV))
ALLOCATE(YSUBTITLE(50))
!
ZLES_BUDGET(:,:,:,:) = XUNDEF
!-------------------------------------------------------------------------------
!
!*      2.  total scalar variance budget
!           ----------------------------
!
!
YGROUP='BU_SV2  '
!
ILES=0
ILES_STA=ILES
!
!* 2.1 production by mean gradients
!      ----------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP M'
ILES_P1=ILES
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV)= - 2. * XLES_SUBGRID_WSv(:,:,1,JSV) * XLES_MEAN_dSvdz(:,:,1,JSV)
END DO
!
!
!* 2.2 production by resolved gradients
!      --------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV)= - 2. * XLES_RES_ddxa_Sv_SBG_UaSv(:,:,1,JSV)  &
                             - ZLES_BUDGET(:,:,ILES_P1,JSV)
END DO
!
!
!* 2.3 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_WSv2(:,:,1,:)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JSV=1,NSV
  DO JK=2,NLES_K-1
    ZLES_BUDGET(JK,:,ILES,JSV) = - ( XLES_SUBGRID_WSv2 (JK+1,:,1,JSV)  &
                                    -XLES_SUBGRID_WSv2 (JK-1,:,1,JSV)) &
                            / ( XLES_Z            (JK+1)          &
                               -XLES_Z            (JK-1)        )
  END DO
  ZLES_BUDGET(1     ,:,ILES,JSV) = ZLES_BUDGET(2       ,:,ILES,JSV)
  ZLES_BUDGET(NLES_K,:,ILES,JSV) = ZLES_BUDGET(NLES_K-1,:,ILES,JSV)
END DO
END IF
!
!
!* 2.4 dissipation
!      -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS    '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) =  XLES_SUBGRID_DISS_Sv2(:,:,1,JSV)
END DO
!
!
!* 2.5 residual of subgrid budget
!      --------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG RESI'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
  DO JLES=ILES_STA+1,ILES-1
    ZLES_BUDGET(:,:,ILES,JSV) = ZLES_BUDGET(:,:,ILES,JSV) - ZLES_BUDGET(:,:,JLES,JSV)
  END DO
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
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_TEND,JSV)
END DO
!
!* 2.7 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_Sv2(:,:,NLES_ADVM,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES ADV '
!
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_ADVM,JSV)
  END DO
END IF
!
!* 2.8 forcing
!      -------
!
IF ( ANY(XLES_BU_RES_Sv2(:,:,NLES_FORC,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES FORC'
!
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_FORC,JSV)
  END DO
END IF
!
!
!* 2.9 production
!      ----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_DP,JSV)
END DO
!
!* 2.10 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_TR,JSV)
END DO
!
!
!* 2.11 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_VTURB,JSV) &
                            + XLES_BU_RES_Sv2(:,:,NLES_HTURB,JSV)
END DO
!
!* 2.11 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_Sv2(:,:,NLES_DIFF,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES NUMD'
!
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_DIFF,JSV)
  END DO
END IF
!
!* 2.11 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_Sv2(:,:,NLES_RELA,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES RELA'
!
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_RELA,JSV)
  END DO
END IF
!
!* 2.11 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_Sv2(:,:,NLES_NEST,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES NEST'
!
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_NEST,JSV)
  END DO
END IF
!
!* 2.11 other effects
!       -------------
!
IF ( ANY(XLES_BU_RES_Sv2(:,:,NLES_MISC,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES MISC'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_Sv2(:,:,NLES_MISC,JSV)
  END DO
END IF
!
!* 2.12 residual of resolved budget
!       ---------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RESI'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
  DO JLES=ILES_STA+1,ILES-1
    ZLES_BUDGET(:,:,ILES,JSV) = ZLES_BUDGET(:,:,ILES,JSV) - ZLES_BUDGET(:,:,JLES,JSV)
  END DO
END DO
!
!* 2.13 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG TEND'
!
IF (NLES_TIMES>2) THEN
  DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
    DO JK=1,NLES_K
      DO JT=2,NLES_TIMES-1
        ZLES_BUDGET(JK,JT,ILES,JSV) =- ( XLES_SUBGRID_Sv2 (JK  ,JT+1,1,JSV) &
                                       - XLES_SUBGRID_Sv2 (JK  ,JT-1,1,JSV))&
                                      / (2.* XLES_TEMP_SAMPLING)
      END DO
      ZLES_BUDGET(JK,1         ,ILES,JSV) = ZLES_BUDGET(JK,2           ,ILES,JSV)
      ZLES_BUDGET(JK,NLES_TIMES,ILES,JSV) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES,JSV)
    END DO
  END DO
END IF
!
!* 2.14 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JSV=1,NSV
  DO JK=2,NLES_K-1
    ZLES_BUDGET(JK,:,ILES,JSV)=  -XLES_MEAN_W(JK,:,1)                &
                               * ( XLES_SUBGRID_Sv2(JK+1,:,1,JSV)    &
                                  -XLES_SUBGRID_Sv2(JK-1,:,1,JSV)  ) &
                               / ( XLES_Z          (JK+1)            &
                                  -XLES_Z          (JK-1)          )
  END DO
  ZLES_BUDGET(1     ,:,ILES,JSV) = ZLES_BUDGET(2       ,:,ILES,JSV)
  ZLES_BUDGET(NLES_K,:,ILES,JSV) = ZLES_BUDGET(NLES_K-1,:,ILES,JSV)
END DO
!
!* 2.15 neglected term: advection for subgrid quantity
!       ----------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JSV=1,NSV
  DO JK=2,NLES_K-1
    ZLES_BUDGET(JK,:,ILES,JSV)=-( XLES_RES_W_SBG_Sv2   (JK+1,:,1,JSV)    &
                                 -XLES_RES_W_SBG_Sv2   (JK-1,:,1,JSV)  ) &
                              / ( XLES_Z           (JK+1)                &
                                 -XLES_Z           (JK-1)              )
  END DO
  ZLES_BUDGET(1     ,:,ILES,JSV) = ZLES_BUDGET(2       ,:,ILES,JSV)
  ZLES_BUDGET(NLES_K,:,ILES,JSV) = ZLES_BUDGET(NLES_K-1,:,ILES,JSV)
END DO
!
!
!* 2.16 writing
!       -------
!
ALLOCATE(ZSV_BUDGET(NLES_K,NLES_TIMES,ILES,NSV))
DO JSV=1,NSV
  DO JP=1,ILES
    ZSV_BUDGET(:,:,JP,JSV) = ZLES_BUDGET(:,:,JP,JSV)
  END DO
END DO

YTITLE = "Sv variance budget  "
CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),YTITLE//YSUBTITLE(:ILES),"kg2/kg2/s",ZSV_BUDGET,HLES_AVG)
!
DEALLOCATE(ZSV_BUDGET)
!-------------------------------------------------------------------------------
!
!*      3.  total water flux budget
!           -----------------------
!
!
YGROUP = 'BU_WSV  '
!
!
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
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) =  - XLES_SUBGRID_W2 (:,:,1)     * XLES_MEAN_DSvDZ(:,:,1,JSV) &
                               - XLES_SUBGRID_WSv(:,:,1,JSV) * XLES_MEAN_DWDZ (:,:,1)
END DO
!
!
!* 3.2 production by gradient of resolved motions
!     -------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG DP R'
ILES_P2=ILES
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV)=- XLES_RES_ddz_Sv_SBG_W2(:,:,1,JSV) &
                            - ZLES_BUDGET(:,:,ILES_P1,JSV)
END DO
!
!
!
!* 3.3 turbulent transport
!      -------------------
!
IF ( ANY(XLES_SUBGRID_W2Sv(:,:,1,:)/= 0.) ) THEN
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TR  '
!
DO JSV=1,NSV
  DO JK=2,NLES_K-1
    ZLES_BUDGET(JK,:,ILES,JSV) = - ( XLES_SUBGRID_W2Sv (JK+1,:,1,JSV)  &
                                    -XLES_SUBGRID_W2Sv (JK-1,:,1,JSV)) &
                                 / ( XLES_Z            (JK+1)          &
                                    -XLES_Z            (JK-1)        )
  END DO
  ZLES_BUDGET(1     ,:,ILES,JSV) = ZLES_BUDGET(2       ,:,ILES,JSV)
  ZLES_BUDGET(NLES_K,:,ILES,JSV) = ZLES_BUDGET(NLES_K-1,:,ILES,JSV)
END DO
END IF
!
!
!* 3.4 presso-correlations
!      -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG PRES'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) =  XLES_SUBGRID_SvPz(:,:,1,JSV)
END DO
!
!
!* 3.5 thermal production
!      ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG TP  '
!
IF (LUSERV) THEN
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) =  XG * XLES_SUBGRID_SvThv(:,:,1,JSV)   &
                                    / XLES_MEAN_Thv     (:,:,1)
  END DO
ELSE
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) =  XG * XLES_SUBGRID_SvThv(:,:,1,JSV)   &
                                    / XLES_MEAN_Th      (:,:,1)
  END DO
END IF
!
!
!* 3.6 dissipation
!      -----------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' DISS    '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
END DO
!
!
!* 3.7 residual of subgrid budget
!      --------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' SBG RESI'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
  DO JLES=ILES_STA+1,ILES-1
    ZLES_BUDGET(:,:,ILES,JSV) = ZLES_BUDGET(:,:,ILES,JSV) - ZLES_BUDGET(:,:,JLES,JSV)
  END DO
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
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_TEND,JSV)
END DO
!
!* 3.9 advection
!      ---------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_ADVM,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES ADV '
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_ADVM,JSV)
  END DO
END IF
!
!* 3.10 forcing
!       -------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_FORC,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES FORC'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_FORC,JSV)
  END DO
END IF
!
!* 3.11 production by temperature gradient (and vertical wind gradient)
!       ----------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES DP  '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_DP,JSV)
END DO
!
!* 3.12 turbulent transport
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TR  '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_TR,JSV)
END DO
!
!
!* 3.13 presso-correlations
!       -------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES PRES'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_PRES,JSV)
END DO
!
!
!* 3.14 thermal production
!       ------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES TP  '
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_GRAV,JSV)
END DO
!
!
!* 3.15 effect of subgrid scale motions on the resolved flow
!       ----------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES SBGT'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_VTURB,JSV) + XLES_BU_RES_WSv(:,:,NLES_HTURB,JSV)
END DO
!
!* 3.15 effect of Coriolis
!       ------------------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_COR,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES CORI'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_COR,JSV)
  END DO
END IF
!
!* 3.15 effect of diffusion
!       -------------------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_DIFF,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES NUMD'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_DIFF,JSV)
  END DO
END IF
!
!* 3.15 effect of relaxation
!       --------------------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_RELA,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES RELA'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_RELA,JSV)
  END DO
END IF
!
!* 3.15 effect of nesting
!       -----------------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_NEST,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES NEST'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_NEST,JSV)
  END DO
END IF
!
!* 3.15 other effects
!       -------------
!
IF ( ANY(XLES_BU_RES_WSv(:,:,NLES_MISC,:)/= 0.) ) THEN
  ILES=ILES+1
  YSUBTITLE(ILES) = ' RES MISC'
  !
  DO JSV=1,NSV
    ZLES_BUDGET(:,:,ILES,JSV) = XLES_BU_RES_WSv(:,:,NLES_MISC,JSV)
  END DO
END IF
!
!* 3.16 residual of resolved WSv budget
!       -------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' RES RESI'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
  DO JLES=ILES_STA+1,ILES-1
    ZLES_BUDGET(:,:,ILES,JSV) = ZLES_BUDGET(:,:,ILES,JSV) - ZLES_BUDGET(:,:,JLES,JSV)
  END DO
END DO
!
!* 3.17 neglected term: tendency
!       ------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG TEND'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV) = 0.
  IF (NLES_TIMES>2) THEN
    DO JK=1,NLES_K
      DO JT=2,NLES_TIMES-1
        ZLES_BUDGET(JK,JT,ILES,JSV) =- ( XLES_SUBGRID_WSv (JK  ,JT+1,1,JSV) &
                                       - XLES_SUBGRID_WSv (JK  ,JT-1,1,JSV))&
                                       / (2.* XLES_TEMP_SAMPLING)
      END DO
      ZLES_BUDGET(JK,1         ,ILES,JSV) = ZLES_BUDGET(JK,2           ,ILES,JSV)
      ZLES_BUDGET(JK,NLES_TIMES,ILES,JSV) = ZLES_BUDGET(JK,NLES_TIMES-1,ILES,JSV)
    END DO
  END IF
END DO
!
!
!
!* 3.18 neglected terms : advection for subgrid quantity
!       ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVM'
!
DO JSV=1,NSV
  DO JK=2,NLES_K-1
    ZLES_BUDGET(JK,:,ILES,JSV)= - XLES_MEAN_W(JK,:,1)                       &
                                  * ( XLES_SUBGRID_WSv(JK+1,:,1,JSV)        &
                                     -XLES_SUBGRID_WSv(JK-1,:,1,JSV)      ) &
                                 /  ( XLES_Z          (JK+1)                &
                                     -XLES_Z          (JK-1)              )
  END DO
  ZLES_BUDGET(1     ,:,ILES,JSV) = ZLES_BUDGET(2       ,:,ILES,JSV)
  ZLES_BUDGET(NLES_K,:,ILES,JSV) = ZLES_BUDGET(NLES_K-1,:,ILES,JSV)
END DO
!
!* 3.19 neglected terms : advection for subgrid quantity
!       ------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG ADVR'
!
DO JSV=1,NSV
  DO JK=2,NLES_K-1
    ZLES_BUDGET(JK,:,ILES,JSV)=-( XLES_RES_W_SBG_WSv(JK+1,:,1,JSV)   &
                                 -XLES_RES_W_SBG_WSv(JK-1,:,1,JSV) ) &
                              / ( XLES_Z          (JK+1)             &
                                 -XLES_Z          (JK-1)           )
  END DO
  ZLES_BUDGET(1     ,:,ILES,JSV) = ZLES_BUDGET(2       ,:,ILES,JSV)
  ZLES_BUDGET(NLES_K,:,ILES,JSV) = ZLES_BUDGET(NLES_K-1,:,ILES,JSV)
END DO
!
!* 3.20 neglected terms : production by gradient of vertical velocity for subgrid quantity
!       ----------------------------------------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG DPGW'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV)=- XLES_RES_ddxa_W_SBG_UaSv(:,:,1,JSV)
END DO
!
!
!* 3.21 neglected terms : production by hor. gradient of Thl for subgrid quantity
!       -------------------------------------------------------------------------
!
ILES=ILES+1
YSUBTITLE(ILES) = ' NSG DPGT'
!
DO JSV=1,NSV
  ZLES_BUDGET(:,:,ILES,JSV)=-XLES_RES_ddxa_Sv_SBG_UaW(:,:,1,JSV)       &
                            -ZLES_BUDGET(:,:,ILES_P1,JSV) -ZLES_BUDGET(:,:,ILES_P2,JSV)
END DO
!
!
!* 3.22 writing
!       -------
!
ALLOCATE(ZSV_BUDGET(NLES_K,NLES_TIMES,ILES,NSV))
DO JSV=1,NSV
  DO JP=1,ILES
    ZSV_BUDGET(:,:,JP,JSV) = ZLES_BUDGET(:,:,JP,JSV)
  END DO
END DO

YTITLE = "Sv flux budget      "
CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,YGROUP,YSUBTITLE(:ILES),YTITLE//YSUBTITLE(:ILES),"mkg/kg/s2",ZSV_BUDGET,HLES_AVG)
!
DEALLOCATE(ZSV_BUDGET)
!-------------------------------------------------------------------------------
!
DEALLOCATE(ZLES_BUDGET)
DEALLOCATE(YSUBTITLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_LES_SV_BUDGET_n 
