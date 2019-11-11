!     #########
      SUBROUTINE DEFAULT_BLOWSNW
!     ########################################################################
!
!!****  *DEFAULT_BLOWSNW* - routine to set default values for the configuration for BLOWSNW
!!
!!    PURPOSE
!!    -------
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
!!      Vincent Vionnet CNRM
!!
!!    MODIFICATIONS
!!    -------------
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BLOWSNW_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
! Set initial values of variables. These are modified by namelist


REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEFAULT_BLOWSNW',0,ZHOOK_HANDLE)

!     Grain size distribution is a two parameter gamma distribution
XEMIRADIUS_SNW= 75e-6
XEMIALPHA_SNW = 3.

!     Default value of the ratio between diffusion coefficient for momentum variables and blowing snow variables
XRSNOW_SBL = 4

LSNOW_SALT = .FALSE. ! Flag to fix snow concentration in the saltation layer      
! Snow concentration in the saltation layer (kg_{snow}/m3_{air}) if LSNOW_SALT = T
XCONC_SALT = 0.
!
LSNOW_PRECIP = .FALSE. ! Flag to impose uniform and constant precipitation rate
! Fixed snow precipitation rate SWE (kg/ (m2 s)) if LSNOW_PRECIP = T
XSNOW_PRECIP = 0.

LBLOWSNW_CANOSUBL = .FALSE. ! Flag to activate sublimation of Canopy variables

LBLOWSNW_CANODIAG = .FALSE. ! Flag to get additional diagnostic at Canopy levels : mean radius 
                             ! and number and mass variables in _/m3

LBLOWSNW_ADV    = .TRUE.    ! Flag to account for advection effects on 
                            ! total mass and number in Canopy variables and 
                            ! to compute divergence of
                            ! saltation flux

LSNOW_WIND = .FALSE. ! Flag to activate effects of snow particle on wind profile
! Increase in surface rougnhess due to saltation following Dover (1993) z0_s = z0_ns*(u*/uth*)²
XSNOW_ROUGHNESS = 0.        ! = 0 not activated; =1 activated   
! Buoyancy effect in presence of suspended particles of blowing snow.
XSNOW_BUOYANCY  = 0.        ! = 0 not activated; =1 activated  

CSNOW_SALT = 'SORE' ! Paramaterization to compute particle flux in saltation
!             'POME': Pomeroy and Gray, 1990
!             'SORE': Sorensen (1991) : used at SLF before Doorshot's model
!             'MANN': Concentration in the saltation layer is computed according to Mann
!                     parameterization for particle number density just above the ground (10 cm) 

CSNOW_SEDIM = 'TABC' ! Paramaterization to compute settling velocity
!              'CARR': follow Carrier's drag coefficient
!              'TABC': based on tabuleted values of Carrier's formulation for alpha=3
!              'MITC': based on Mitchell's formulation for settling velocity of ice spheres
!              'NONE': sedimentation is desactivated (for test only!)


IF (LHOOK) CALL DR_HOOK('DEFAULT_BLOWSNW',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_BLOWSNW
