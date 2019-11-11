!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_DYN_n
!     #################
!
!!****  *MODD_DYN$n* - declaration of dynamic control variables 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the dynamic
!     control variables.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DYNn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/05/94                      
!!      Modifications 16/11/94   (Lafore+Pinty)  For NUM_DIFF
!!      Modifications 06/01/95   (Lafore)        For LSTEADY_DMASS
!!      Modifications 28/07/96   (Masson)        Supress LSTEADY_DMASS
!!      Modifications 15/03/98   (Stein)         Add LHO_RELAX for each variables 
!!      Modifications 22/01/01   (Gazen)         Add LHORELAX_SVC2R2, _SVCHEM, _SVLG
!!      Modifications 29/11/02   (Pinty)         Add  LHORELAX_SVC1R3, _SVELEC
!!      Modifications    07/05   (P.Tulet)       Add  relaxation for dust and aerosol
!!      Modifications    05/07   (C.Lac)         Separation of num diffusion
!!      Modifications    07/10   (M.Leriche)     Add relaxation for ice phase chemical
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Modification    07/2017  (V. Vionnet)    Add blowing snow variable
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX, JPSVMAX
IMPLICIT NONE

TYPE DYN_t
!
  INTEGER :: NSTOP          ! Number of time step
  REAL    :: XTSTEP         ! Time step
!
!++++++++++++++++++++++++++++++++++
!PART USED BY THE PRESSURE SOLVER
!++++++++++++++++++++++++++++++++++
!
  REAL, DIMENSION(:,:,:), POINTER :: XBFY=>NULL()     ! Vectors giving the non
  REAL, DIMENSION(:,:,:), POINTER :: XBFB=>NULL()     ! Vectors giving the non
  REAL, DIMENSION(:,:,:), POINTER :: XBF_SXP2_YP1_Z=>NULL()     ! Vectors giving the non
  REAL, DIMENSION(:,:,:), POINTER :: XBF=>NULL()      ! vanishing elements of the 
  REAL, DIMENSION(:),     POINTER :: XAF=>NULL(),XCF=>NULL() ! tri-diag matrix in the pressure equation
!
                                                   ! Arrays of sinus or cosinus
                                                   ! values for the FFT 
  REAL, DIMENSION(:),   POINTER :: XTRIGSX=>NULL()  ! in x-direction   
  REAL, DIMENSION(:),   POINTER :: XTRIGSY=>NULL()  ! in y-direction   
  INTEGER, DIMENSION(:),POINTER :: NIFAXX =>NULL()  ! Decomposition in prime numbers 
  INTEGER, DIMENSION(:),POINTER :: NIFAXY =>NULL()  ! for the FFT in x and y directions 
  CHARACTER(LEN=5)       :: CPRESOPT     ! Choice of the pressure solver
  INTEGER                :: NITR         ! Number of iterations for the 
                                            ! pressure solver
  LOGICAL                :: LITRADJ      ! Choice to adjust the number of 
                                            !solver iterations during
                                            !the simulation
  LOGICAL                :: LRES         ! Choice of a different residual
                                         ! divergence limit
  REAL                   :: XRES         ! Value of residual divergence limit
  REAL                   :: XRELAX       ! relaxation coefficient for the
                                         ! Richardson's method
!
  REAL :: XDXHATM                        ! mean grid increment in the 
  REAL :: XDYHATM                        ! x and  y directions

  REAL, DIMENSION (:), POINTER  :: XRHOM=>NULL()  !  mean of XRHODJ on the plane x y 
                                                 !  localized at a mass level
!
!++++++++++++++++++++++++++++++++++
!PART USED BY THE ABSORBING LAYERS 
!++++++++++++++++++++++++++++++++++
!
  INTEGER      :: NALBOT     ! Vertical index corresponding to the     
                                ! absorbing layer base
!  
  INTEGER      :: NALBAS     ! Vertical index corresponding to the     
                                ! absorbing layer base
!                                
  REAL, DIMENSION(:),   POINTER :: XALK=>NULL()  ! Function of the absorbing
                                ! layer damping coefficient defined for
                                !   u,v,and theta
  REAL, DIMENSION(:),   POINTER :: XALKW=>NULL() ! Idem but defined for w 
!
  REAL, DIMENSION(:),   POINTER :: XALKBAS=>NULL()  ! Function of the absorbing
                                ! layer damping coefficient defined for
                                !   u,v,and theta
  REAL, DIMENSION(:),   POINTER :: XALKWBAS=>NULL() ! Idem but defined for w 
!
  LOGICAL ::  LVE_RELAX ! switch to activate the VErtical RELAXation
  LOGICAL ::  LVE_RELAX_GRD ! switch to activate the VErtical RELAXation
!
!                            switch to activate the HOrizontal RELAXation
!  LOGICAL :: LHORELAX_UVWTH
!
  LOGICAL :: LHORELAX_RV, LHORELAX_RC, LHORELAX_RR, LHORELAX_RI
  LOGICAL :: LHORELAX_RS, LHORELAX_RG, LHORELAX_RH
!
!  LOGICAL :: LHORELAX_TKE
!
  LOGICAL :: LHORELAX_SVC2R2
  LOGICAL :: LHORELAX_SVC1R3
  LOGICAL :: LHORELAX_SVLIMA
  LOGICAL :: LHORELAX_SVELEC
  LOGICAL :: LHORELAX_SVCHEM
  LOGICAL :: LHORELAX_SVCHIC
  LOGICAL :: LHORELAX_SVLG
  LOGICAL :: LHORELAX_SVDST
  LOGICAL :: LHORELAX_SVSLT
  LOGICAL :: LHORELAX_SVAER
  LOGICAL :: LHORELAX_SVPP 
#ifdef MNH_FOREFIRE
  LOGICAL :: LHORELAX_SVFF
#endif
  LOGICAL :: LHORELAX_SVCS 
  LOGICAL :: LHORELAX_SVSNW  
  LOGICAL, DIMENSION(:),POINTER :: LHORELAX_SV =>NULL()
!
  REAL    :: XRIMKMAX   ! Max. value of the horiz. relaxation coeff.
!  INTEGER :: NRIMX,NRIMY! Number of points in the lateral absorbing
                           ! layer in the x and y directions 
! sizes of the West-east total LB area
  INTEGER :: NSIZELBX_ll,NSIZELBXU_ll      ! for T,V,W and u 
  INTEGER :: NSIZELBXTKE_ll                ! for TKE       
  INTEGER :: NSIZELBXR_ll,NSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
  INTEGER :: NSIZELBY_ll,NSIZELBYV_ll      ! for T,U,W  and v
  INTEGER :: NSIZELBYTKE_ll                ! for TKE       
  INTEGER :: NSIZELBYR_ll,NSIZELBYSV_ll    ! for Rx and SV    
  LOGICAL, DIMENSION(:,:),  POINTER   :: LMASK_RELAX=>NULL()  ! Mask for lateral
                           ! relaxation: True where it has to be performed
  REAL, DIMENSION(:,:),     POINTER   :: XKURELAX=>NULL() ! Horizontal relaxation
  REAL, DIMENSION(:,:),     POINTER   :: XKVRELAX=>NULL() ! coefficients for the
  REAL, DIMENSION(:,:),     POINTER   :: XKWRELAX=>NULL() ! u, v and mass locations
!
!++++++++++++++++++++++++++++++++++++
!PART USED BY THE NUMERICAL DIFFUSION
!++++++++++++++++++++++++++++++++++++
!
  REAL              :: XT4DIFU ! Damping time scale for 2*dx wavelength
                               ! specified for the 4nd order num. diffusion
                               ! for momentum
  REAL              :: XT4DIFTH! for theta and mixing ratios                  
  REAL              :: XT4DIFSV! for scalar variables                         
  REAL              :: XDK2U   ! 2nd order num. diffusion coef. /dx2
                               ! for momentum
  REAL              :: XDK4U   ! 4nd order num. diffusion coef. /dx4
                               ! for momentum
  REAL              :: XDK2TH  ! for theta and mixing ratios                                        
  REAL              :: XDK4TH  ! for theta and mixing ratios           
  REAL              :: XDK2SV  ! for scalar variables                                                               
  REAL              :: XDK4SV  ! for scalar variables                  
!
END TYPE DYN_t

TYPE(DYN_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: DYN_MODEL
LOGICAL    , DIMENSION(JPMODELMAX),         SAVE :: DYN_FIRST_CALL = .TRUE.

INTEGER, POINTER :: NSTOP=>NULL()
REAL, POINTER :: XTSTEP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XBFY=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XBFB=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XBF_SXP2_YP1_Z=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XBF=>NULL()
REAL, DIMENSION(:),     POINTER :: XAF=>NULL(),XCF=>NULL()
REAL, DIMENSION(:),   POINTER :: XTRIGSX=>NULL()
REAL, DIMENSION(:),   POINTER :: XTRIGSY=>NULL()
INTEGER, DIMENSION(:), POINTER :: NIFAXX=>NULL()
INTEGER, DIMENSION(:), POINTER :: NIFAXY=>NULL()
CHARACTER(LEN=5), POINTER :: CPRESOPT=>NULL()
INTEGER, POINTER :: NITR=>NULL()
LOGICAL, POINTER :: LITRADJ=>NULL()
LOGICAL, POINTER :: LRES=>NULL()
REAL, POINTER :: XRES=>NULL()
REAL, POINTER :: XRELAX=>NULL()
REAL, POINTER :: XDXHATM=>NULL()
REAL, POINTER :: XDYHATM=>NULL()
REAL, DIMENSION (:), POINTER  :: XRHOM=>NULL()
INTEGER, POINTER :: NALBOT=>NULL()
REAL, DIMENSION(:),   POINTER :: XALK=>NULL()
REAL, DIMENSION(:),   POINTER :: XALKW=>NULL()
INTEGER, POINTER :: NALBAS=>NULL()
REAL, DIMENSION(:),   POINTER :: XALKBAS=>NULL()
REAL, DIMENSION(:),   POINTER :: XALKWBAS=>NULL()
LOGICAL, POINTER :: LVE_RELAX=>NULL()
LOGICAL, POINTER :: LVE_RELAX_GRD=>NULL()
LOGICAL, POINTER :: LHORELAX_UVWTH=>NULL()
LOGICAL, POINTER :: LHORELAX_RV=>NULL(), LHORELAX_RC=>NULL(), LHORELAX_RR=>NULL(), LHORELAX_RI=>NULL()
LOGICAL, POINTER :: LHORELAX_RS=>NULL(), LHORELAX_RG=>NULL(), LHORELAX_RH=>NULL()
LOGICAL, POINTER :: LHORELAX_TKE=>NULL()
LOGICAL, POINTER :: LHORELAX_SVC2R2=>NULL()
LOGICAL, POINTER :: LHORELAX_SVC1R3=>NULL()
LOGICAL, POINTER :: LHORELAX_SVLIMA=>NULL()
LOGICAL, POINTER :: LHORELAX_SVELEC=>NULL()
LOGICAL, POINTER :: LHORELAX_SVCHEM=>NULL()
LOGICAL, POINTER :: LHORELAX_SVCHIC=>NULL()
LOGICAL, POINTER :: LHORELAX_SVLG=>NULL()
LOGICAL, POINTER :: LHORELAX_SVDST=>NULL()
LOGICAL, POINTER :: LHORELAX_SVSLT=>NULL()
LOGICAL, POINTER :: LHORELAX_SVAER=>NULL()
LOGICAL, POINTER :: LHORELAX_SVPP=>NULL()
#ifdef MNH_FOREFIRE
LOGICAL, POINTER :: LHORELAX_SVFF=>NULL()
#endif
LOGICAL, POINTER :: LHORELAX_SVCS=>NULL()
LOGICAL, POINTER :: LHORELAX_SVSNW=>NULL()
LOGICAL, DIMENSION(:), POINTER :: LHORELAX_SV=>NULL()
REAL, POINTER :: XRIMKMAX=>NULL()
INTEGER, POINTER :: NRIMX=>NULL(),NRIMY=>NULL()
INTEGER, POINTER :: NSIZELBX_ll=>NULL(),NSIZELBXU_ll=>NULL()
INTEGER, POINTER :: NSIZELBXTKE_ll=>NULL()
INTEGER, POINTER :: NSIZELBXR_ll=>NULL(),NSIZELBXSV_ll=>NULL()
INTEGER, POINTER :: NSIZELBY_ll=>NULL(),NSIZELBYV_ll=>NULL()
INTEGER, POINTER :: NSIZELBYTKE_ll=>NULL()
INTEGER, POINTER :: NSIZELBYR_ll=>NULL(),NSIZELBYSV_ll=>NULL()
LOGICAL, DIMENSION(:,:),  POINTER   :: LMASK_RELAX=>NULL()
REAL, DIMENSION(:,:),     POINTER   :: XKURELAX=>NULL()
REAL, DIMENSION(:,:),     POINTER   :: XKVRELAX=>NULL()
REAL, DIMENSION(:,:),     POINTER   :: XKWRELAX=>NULL()
REAL, POINTER :: XT4DIFU=>NULL()
REAL, POINTER :: XDK2U=>NULL()
REAL, POINTER :: XDK4U=>NULL()
REAL, POINTER :: XT4DIFTH=>NULL()
REAL, POINTER :: XDK2TH=>NULL()
REAL, POINTER :: XDK4TH=>NULL()
REAL, POINTER :: XT4DIFSV=>NULL()
REAL, POINTER :: XDK2SV=>NULL()
REAL, POINTER :: XDK4SV=>NULL()

CONTAINS

SUBROUTINE DYN_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
IF (DYN_FIRST_CALL(KTO)) THEN
ALLOCATE (DYN_MODEL(KTO)%NIFAXX(19))
ALLOCATE (DYN_MODEL(KTO)%NIFAXY(19))
ALLOCATE (DYN_MODEL(KTO)%LHORELAX_SV(JPSVMAX))
DYN_FIRST_CALL(KTO) = .FALSE.
ENDIF
! Save current state for allocated arrays
DYN_MODEL(KFROM)%XBFY=>XBFY
DYN_MODEL(KFROM)%XBFB=>XBFB
DYN_MODEL(KFROM)%XBF_SXP2_YP1_Z=>XBF_SXP2_YP1_Z
DYN_MODEL(KFROM)%XBF=>XBF
DYN_MODEL(KFROM)%XAF=>XAF
DYN_MODEL(KFROM)%XCF=>XCF
DYN_MODEL(KFROM)%XTRIGSX=>XTRIGSX
DYN_MODEL(KFROM)%XTRIGSY=>XTRIGSY
DYN_MODEL(KFROM)%XRHOM=>XRHOM
DYN_MODEL(KFROM)%XALK=>XALK
DYN_MODEL(KFROM)%XALKW=>XALKW
DYN_MODEL(KFROM)%XALKBAS=>XALKBAS
DYN_MODEL(KFROM)%XALKWBAS=>XALKWBAS
DYN_MODEL(KFROM)%LMASK_RELAX=>LMASK_RELAX
DYN_MODEL(KFROM)%XKURELAX=>XKURELAX
DYN_MODEL(KFROM)%XKVRELAX=>XKVRELAX
DYN_MODEL(KFROM)%XKWRELAX=>XKWRELAX
!
! Current model is set to model KTO
NSTOP=>DYN_MODEL(KTO)%NSTOP
XTSTEP=>DYN_MODEL(KTO)%XTSTEP
XBFY=>DYN_MODEL(KTO)%XBFY
XBFB=>DYN_MODEL(KTO)%XBFB
XBF_SXP2_YP1_Z=>DYN_MODEL(KTO)%XBF_SXP2_YP1_Z
XBF=>DYN_MODEL(KTO)%XBF
XAF=>DYN_MODEL(KTO)%XAF
XCF=>DYN_MODEL(KTO)%XCF
XTRIGSX=>DYN_MODEL(KTO)%XTRIGSX
XTRIGSY=>DYN_MODEL(KTO)%XTRIGSY
NIFAXX=>DYN_MODEL(KTO)%NIFAXX
NIFAXY=>DYN_MODEL(KTO)%NIFAXY
CPRESOPT=>DYN_MODEL(KTO)%CPRESOPT
NITR=>DYN_MODEL(KTO)%NITR
LITRADJ=>DYN_MODEL(KTO)%LITRADJ
LRES=>DYN_MODEL(KTO)%LRES          
XRES=>DYN_MODEL(KTO)%XRES
XRELAX=>DYN_MODEL(KTO)%XRELAX
XDXHATM=>DYN_MODEL(KTO)%XDXHATM
XDYHATM=>DYN_MODEL(KTO)%XDYHATM
XRHOM=>DYN_MODEL(KTO)%XRHOM
NALBOT=>DYN_MODEL(KTO)%NALBOT
XALK=>DYN_MODEL(KTO)%XALK
XALKW=>DYN_MODEL(KTO)%XALKW
NALBAS=>DYN_MODEL(KTO)%NALBAS
XALKBAS=>DYN_MODEL(KTO)%XALKBAS
XALKWBAS=>DYN_MODEL(KTO)%XALKWBAS
LVE_RELAX=>DYN_MODEL(KTO)%LVE_RELAX
LVE_RELAX_GRD=>DYN_MODEL(KTO)%LVE_RELAX_GRD
!LHORELAX_UVWTH=>DYN_MODEL(KTO)%LHORELAX_UVWTH !Done in FIELDLIST_GOTO_MODEL
LHORELAX_RV=>DYN_MODEL(KTO)%LHORELAX_RV
LHORELAX_RC=>DYN_MODEL(KTO)%LHORELAX_RC
LHORELAX_RR=>DYN_MODEL(KTO)%LHORELAX_RR
LHORELAX_RI=>DYN_MODEL(KTO)%LHORELAX_RI
LHORELAX_RS=>DYN_MODEL(KTO)%LHORELAX_RS
LHORELAX_RG=>DYN_MODEL(KTO)%LHORELAX_RG
LHORELAX_RH=>DYN_MODEL(KTO)%LHORELAX_RH
!LHORELAX_TKE=>DYN_MODEL(KTO)%LHORELAX_TKE !Done in FIELDLIST_GOTO_MODEL
LHORELAX_SVC2R2=>DYN_MODEL(KTO)%LHORELAX_SVC2R2
LHORELAX_SVC1R3=>DYN_MODEL(KTO)%LHORELAX_SVC1R3
LHORELAX_SVLIMA=>DYN_MODEL(KTO)%LHORELAX_SVLIMA
LHORELAX_SVELEC=>DYN_MODEL(KTO)%LHORELAX_SVELEC
LHORELAX_SVCHEM=>DYN_MODEL(KTO)%LHORELAX_SVCHEM
LHORELAX_SVCHIC=>DYN_MODEL(KTO)%LHORELAX_SVCHIC
LHORELAX_SVLG=>DYN_MODEL(KTO)%LHORELAX_SVLG
LHORELAX_SVDST=>DYN_MODEL(KTO)%LHORELAX_SVDST
LHORELAX_SVSLT=>DYN_MODEL(KTO)%LHORELAX_SVSLT
LHORELAX_SVAER=>DYN_MODEL(KTO)%LHORELAX_SVAER
LHORELAX_SVPP=>DYN_MODEL(KTO)%LHORELAX_SVPP
#ifdef MNH_FOREFIRE
LHORELAX_SVFF=>DYN_MODEL(KTO)%LHORELAX_SVFF
#endif
LHORELAX_SVCS=>DYN_MODEL(KTO)%LHORELAX_SVCS
LHORELAX_SVSNW=>DYN_MODEL(KTO)%LHORELAX_SVSNW
LHORELAX_SV=>DYN_MODEL(KTO)%LHORELAX_SV
XRIMKMAX=>DYN_MODEL(KTO)%XRIMKMAX
!NRIMX=>DYN_MODEL(KTO)%NRIMX !Done in FIELDLIST_GOTO_MODEL
!NRIMY=>DYN_MODEL(KTO)%NRIMY !Done in FIELDLIST_GOTO_MODEL
NSIZELBX_ll=>DYN_MODEL(KTO)%NSIZELBX_ll
NSIZELBXU_ll=>DYN_MODEL(KTO)%NSIZELBXU_ll
NSIZELBXTKE_ll=>DYN_MODEL(KTO)%NSIZELBXTKE_ll
NSIZELBXR_ll=>DYN_MODEL(KTO)%NSIZELBXR_ll
NSIZELBXSV_ll=>DYN_MODEL(KTO)%NSIZELBXSV_ll
NSIZELBY_ll=>DYN_MODEL(KTO)%NSIZELBY_ll
NSIZELBYV_ll=>DYN_MODEL(KTO)%NSIZELBYV_ll
NSIZELBYTKE_ll=>DYN_MODEL(KTO)%NSIZELBYTKE_ll
NSIZELBYR_ll=>DYN_MODEL(KTO)%NSIZELBYR_ll
NSIZELBYSV_ll=>DYN_MODEL(KTO)%NSIZELBYSV_ll
LMASK_RELAX=>DYN_MODEL(KTO)%LMASK_RELAX
XKURELAX=>DYN_MODEL(KTO)%XKURELAX
XKVRELAX=>DYN_MODEL(KTO)%XKVRELAX
XKWRELAX=>DYN_MODEL(KTO)%XKWRELAX
XT4DIFU=>DYN_MODEL(KTO)%XT4DIFU
XDK2U=>DYN_MODEL(KTO)%XDK2U
XDK4U=>DYN_MODEL(KTO)%XDK4U
XT4DIFTH=>DYN_MODEL(KTO)%XT4DIFTH
XDK2TH=>DYN_MODEL(KTO)%XDK2TH
XDK4TH=>DYN_MODEL(KTO)%XDK4TH
XT4DIFSV=>DYN_MODEL(KTO)%XT4DIFSV
XDK2SV=>DYN_MODEL(KTO)%XDK2SV
XDK4SV=>DYN_MODEL(KTO)%XDK4SV

END SUBROUTINE DYN_GOTO_MODEL

END MODULE MODD_DYN_n
