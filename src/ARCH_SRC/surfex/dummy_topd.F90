MODULE MODD_BUDGET_COUPL_ROUT
END MODULE MODD_BUDGET_COUPL_ROUT
!     ######################
!     ######################
!     ######################      
MODULE MODD_COUPLING_TOPD
END MODULE MODD_COUPLING_TOPD 
!     ######################
!     ######################
!     ######################        
MODULE MODD_DUMMY_EXP_PROFILE
!
END MODULE MODD_DUMMY_EXP_PROFILE
!     ######################
!     ######################
!     ######################
MODULE MODD_TOPD_PAR
!
REAL,    PARAMETER :: XSTEPK = 0.05 
INTEGER, PARAMETER :: NDIM = 20 
INTEGER, PARAMETER :: JPCAT = 10    
INTEGER :: NUNIT = 19
REAL, DIMENSION(JPCAT) :: XF_PARAM_BV
REAL, DIMENSION(JPCAT) :: XC_DEPTH_RATIO_BV
END MODULE MODD_TOPD_PAR

!     ######################
!     ######################
!     ######################
MODULE MODD_TOPODYN
USE MODD_TOPD_PAR, ONLY : JPCAT
CHARACTER(LEN=15), DIMENSION(JPCAT) :: CCAT     ! base name for topographic files
INTEGER                             :: NNCAT    ! catchments number
!
INTEGER                             :: NNB_TOPD_STEP   ! number of TOPODYN time steps
REAL                                :: XTOPD_STEP      ! TOPODYN time step
!
INTEGER                             :: NMESHT   ! maximal number of catchments meshes

REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDMAXT   ! maximal deficit on TOPODYN grid (m)
REAL, ALLOCATABLE, DIMENSION(:)     :: XDXT     ! catchment grid mesh size (m)
REAL, ALLOCATABLE, DIMENSION(:)     :: XMPARA   ! M parameter on TOPODYN grid (m)

INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNMC     ! catchments pixels number
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XCONN    ! pixels reference number and 
                                                ! connections between
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NLINE    ! second index of the pixel in the array 
                                                ! XCONN
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTANB    ! pixels topographic slope (Tan(Beta))
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XSLOP    ! pixels topographic slope/length flow

!Variables à priori inutiles
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDAREA   ! drainage area (aire drainee)

! Variables defining the catchments

INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNXC     ! number of topographic grid points on 
                                                ! abscissa axis
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNYC     ! number of topographic grid points on ordinate 
                                                ! axis
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNPT     ! number of pixels in the topographic 
                                                ! domain
INTEGER                             :: NPMAX    ! maximal number of pixels in the 
                                                ! topographic grid

REAL, ALLOCATABLE, DIMENSION(:)     :: XX0,XY0  ! coordinates bottom-left pixel of each 
                                                ! topographic domain

REAL, ALLOCATABLE, DIMENSION(:)     :: XNUL     ! undefined value in topographic files

REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTOPD    ! topographic values in topographic files
REAL, DIMENSION(JPCAT)              :: XRTOP_D2 ! depth used by topodyn for lateral transfers
                                                ! (expressed in ratio of isba d2)
                                                !
! Variables used in routing module 
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNISO    ! number of time step for the isochrones
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XCISO    ! isochrones routing constants 

REAL, DIMENSION(JPCAT)              :: XQINIT   ! Initial discharge at the outlet of the catchments
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XQTOT    ! Total discharge at the outlet of the catchments

REAL, DIMENSION(JPCAT)              :: XSPEEDR,XSPEEDH ! River and hillslope speed
REAL, DIMENSION(JPCAT)              :: XSPEEDG         ! Ground speed
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDRIV, XDHIL    ! River and hillslope distances
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDGRD           ! Ground distance
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTIME_TOPD      ! Time to go to the outlet
                                                       ! at the soil surface
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTIME_TOPD_DRAIN! Time to go to the outlet in the ground

INTEGER, ALLOCATABLE, DIMENSION(:)  :: NX_STEP_ROUT   ! number of maximal time step to join the outlet of 
                                                ! any catchment

! Variables used in exfiltration module 
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XLAMBDA  ! pure topographic index
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XCSTOPT  ! hydraulic conductivity at saturation on 
                                                ! TOP-LAT grid
                                                !ludo
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XQB_DR
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XQB_RUN
! for topodyn alone
REAL, ALLOCATABLE, DIMENSION(:)   :: XRI,XRI_PREV! recharge on ISBA grid
REAL, ALLOCATABLE, DIMENSION(:)   :: XSRFULL! reservoir of interception for
!TOPODYN only
REAL, ALLOCATABLE, DIMENSION(:,:) :: XDEFT! pixel deficit

END MODULE MODD_TOPODYN
!     ######################
!     ######################
!     ######################
SUBROUTINE INIT_SURF_TOPD(DEC, IO, S, K, NP, NPE, UG, U, HPROGRAM, KI)
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_NP_t, ISBA_NPE_t
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEC
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=*),  INTENT(IN)     :: HPROGRAM      !
INTEGER,           INTENT(IN)     :: KI
!
END SUBROUTINE INIT_SURF_TOPD

      SUBROUTINE ISBA_TO_TOPD(PVARI,PVART)
REAL, DIMENSION(:), INTENT(IN)      :: PVARI   ! variable from ISBA grid
REAL, DIMENSION(:,:), INTENT(OUT)   :: PVART   ! variable for TOPODYN grid
END SUBROUTINE ISBA_TO_TOPD
!     ######################
!     ######################
!     ######################
SUBROUTINE COUPLING_SURF_TOPD (DE, DEC, DC, DMI, G, IO, S, K, NK, NP, NPE, UG, U,HPROGRAM,KI)
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_SFX_GRID_n, ONLY : GRID_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_NK_t, ISBA_NP_t, ISBA_NPE_t
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
TYPE(DIAG_t), INTENT(INOUT) :: DC
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DE
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEC
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMI
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_NK_t), INTENT(INOUT) :: NK
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENt(INOUT) :: NPE
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
CHARACTER(LEN=6), INTENT(IN)         :: HPROGRAM ! program calling surf. schemes
INTEGER,          INTENT(IN)         :: KI       ! Surfex grid dimension
END SUBROUTINE COUPLING_SURF_TOPD
!     ######################
!     ######################
!     ######################
SUBROUTINE PGD_TOPD(HISBA, HGRID, PGRID_PAR, KDIM_FULL, PSSO_SLOPE, HPROGRAM)
! 
 CHARACTER(LEN=*), INTENT(IN) :: HISBA
 CHARACTER(LEN=*), INTENT(IN) :: HGRID
REAL, DIMENSION(:), INTENT(IN) :: PGRID_PAR
 INTEGER, INTENT(IN) :: KDIM_FULL
 REAL, DIMENSION(:), INTENT(INOUT) :: PSSO_SLOPE
!
CHARACTER(LEN=*),  INTENT(IN)     :: HPROGRAM    !
 
END SUBROUTINE PGD_TOPD
!     ######################
!     ######################
!     ######################
SUBROUTINE READ_NAMELISTS_TOPD(HPROGRAM)
 CHARACTER(LEN=*),  INTENT(IN)     :: HPROGRAM    
END SUBROUTINE READ_NAMELISTS_TOPD


