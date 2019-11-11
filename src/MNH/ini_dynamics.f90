!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_INI_DYNAMICS
!     ########################
INTERFACE
SUBROUTINE INI_DYNAMICS(PLON,PLAT,PRHODJ,PTHVREF,PMAP,PZZ,                   &
               PDXHAT,PDYHAT,PZHAT,HLBCX,HLBCY,PTSTEP,                       &
               OVE_RELAX,OVE_RELAX_GRD,OHORELAX_UVWTH,OHORELAX_RV,           &
               OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,  &
               OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                         &
               OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,&
               OHORELAX_SVCHEM,OHORELAX_SVAER,OHORELAX_SVDST,OHORELAX_SVSLT, &
               OHORELAX_SVPP,OHORELAX_SVCS,  OHORELAX_SVCHIC,OHORELAX_SVSNW, &
#ifdef MNH_FOREFIRE
               OHORELAX_SVFF,                                                &
#endif
               PRIMKMAX,KRIMX,KRIMY,PALKTOP,PALKGRD,PALZBOT,PALZBAS,         &
               PT4DIFU,PT4DIFTH,PT4DIFSV,                                    &
               PCORIOX,PCORIOY,PCORIOZ,PCURVX,PCURVY,                        &
               PDXHATM,PDYHATM,PRHOM,PAF,PBFY,PCF,                           &
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                                &
               PALK,PALKW,KALBOT,PALKBAS,PALKWBAS,KALBAS,                    &
               OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX,                     &
               PDK2U,PDK4U,PDK2TH,PDK4TH,PDK2SV,PDK4SV,OZDIFFU,PZDIFFU_HALO2,& 
               PBFB,&
               PBF_SXP2_YP1_Z) !JUAN Z_SPLITING
!  intent in arguments
!
USE MODE_TYPE_ZDIFFU
IMPLICIT NONE
!
REAL, DIMENSION(:,:),   INTENT(IN)        :: PLON,PLAT !Longitude and latitude
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PRHODJ    ! rho J
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PTHVREF   ! virtual potential 
                                      ! temperature of the reference state
REAL, DIMENSION(:,:),   INTENT(IN)        :: PMAP      ! Map factor
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ       ! height 
REAL, DIMENSION(:),     INTENT(IN)        :: PDXHAT     ! Stretching in x direction
REAL, DIMENSION(:),     INTENT(IN)        :: PDYHAT     ! Stretching in y direction
REAL, DIMENSION(:),     INTENT(IN)        :: PZHAT     ! Gal-Chen Height   
CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HLBCX    ! x-direction LBC type
CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HLBCY    ! y-direction LBC type
LOGICAL,                INTENT(IN)        :: OVE_RELAX ! logical
             ! switch to activate the VErtical  RELAXation
LOGICAL,                INTENT(IN)        :: OVE_RELAX_GRD ! logical
             ! switch to activate the VErtical  RELAXation (ground layer)
LOGICAL,            INTENT(IN) :: OHORELAX_UVWTH  ! switch for the 
                       ! horizontal relaxation for U,V,W,TH
LOGICAL,            INTENT(IN) :: OHORELAX_RV     ! switch for the 
                       ! horizontal relaxation for Rv
LOGICAL,            INTENT(IN) :: OHORELAX_RC     ! switch for the 
                       ! horizontal relaxation for Rc
LOGICAL,            INTENT(IN) :: OHORELAX_RR     ! switch for the 
                       ! horizontal relaxation for Rr
LOGICAL,            INTENT(IN) :: OHORELAX_RI     ! switch for the 
                       ! horizontal relaxation for Ri
LOGICAL,            INTENT(IN) :: OHORELAX_RS     ! switch for the 
                       ! horizontal relaxation for Rs
LOGICAL,            INTENT(IN) :: OHORELAX_RG     ! switch for the 
                       ! horizontal relaxation for Rg
LOGICAL,            INTENT(IN) :: OHORELAX_RH     ! switch for the 
                       ! horizontal relaxation for Rh
LOGICAL,            INTENT(IN) :: OHORELAX_TKE    ! switch for the 
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the 
                       ! horizontal relaxation for sv variables
LOGICAL,             INTENT(IN):: OHORELAX_SVC2R2 ! switch for the 
                       ! horizontal relaxation for c2r2 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVC1R3 ! switch for the 
                       ! horizontal relaxation for c1r3 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVELEC ! switch for the 
                       ! horizontal relaxation for elec variables
LOGICAL,             INTENT(IN):: OHORELAX_SVLG   ! switch for the 
                       ! horizontal relaxation for lg variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHEM ! switch for the 
                       ! horizontal relaxation for chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHIC ! switch for the 
                       ! horizontal relaxation for ice chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVAER  ! switch for the 
                       ! horizontal relaxation for aer variables
LOGICAL,             INTENT(IN):: OHORELAX_SVDST  ! switch for the 
                       ! horizontal relaxation for dst variables
LOGICAL,             INTENT(IN):: OHORELAX_SVSLT  ! switch for the 
                       ! horizontal relaxation for slt variables
LOGICAL,             INTENT(IN):: OHORELAX_SVPP   ! switch for the 
                       ! horizontal relaxation for passive pollutants
LOGICAL,             INTENT(IN):: OHORELAX_SVSNW   ! switch for the
                       ! horizontal relaxation for blowing snow variables    
#ifdef MNH_FOREFIRE
LOGICAL,             INTENT(IN):: OHORELAX_SVFF   ! switch for the 
                       ! horizontal relaxation for ForeFire variables
#endif
LOGICAL,             INTENT(IN):: OHORELAX_SVCS   ! switch for the 
                       ! horizontal relaxation for conditional sampling
REAL,                    INTENT(IN)    :: PRIMKMAX !Max. value of the horiz.
                                     ! relaxation coefficients
INTEGER,                INTENT(IN)     :: KRIMX,KRIMY ! Number of points in 
                                     ! the rim zone in the x and y directions 
REAL,     INTENT(IN)   :: PALKTOP    ! Damping coef. at the top of the absorbing
                                     ! layer
REAL,     INTENT(IN)   :: PALKGRD    ! Damping coef. at the top of the absorbing
                                     ! layer                                    
REAL,     INTENT(IN)   :: PALZBOT    ! Height of the absorbing layer base
REAL,     INTENT(IN)   :: PALZBAS    ! Height of the absorbing layer base
REAL,     INTENT(IN)   :: PT4DIFU    ! Damping time scale for 2*dx wavelength
                                     ! specified for the 4nd order num. diffusion
                                     ! for momentum
REAL,     INTENT(IN)   :: PT4DIFTH   ! for meteorological scalar variables
REAL,     INTENT(IN)   :: PT4DIFSV   ! for tracer scalar variables

REAL,     INTENT(IN)   :: PTSTEP     ! Time step    
!
!  intent out arguments
!
REAL, INTENT(OUT) :: PDXHATM                     ! mean grid increment in the x
                                                 ! direction
REAL, INTENT(OUT) :: PDYHATM                     ! mean grid increment in the y
                                                 ! direction
!
REAL, DIMENSION (:), INTENT(OUT) :: PRHOM        !  mean of XRHODJ on the plane x y 
                                                 !  localized at a mass level
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PCORIOX,PCORIOY ! Hor. Coriolis parameters
REAL, DIMENSION(:,:), INTENT(OUT) :: PCORIOZ         ! Vert. Coriolis parameter 
REAL, DIMENSION(:,:), INTENT(OUT) :: PCURVX,PCURVY   ! Curvature coefficients
! 
REAL, DIMENSION(:),     INTENT(OUT) :: PAF ! vectors giving the non-vanishing        
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFY ! elements of the tri-diag matrix        
                                            ! on an y-slice of global physical domain
REAL, DIMENSION(:),     INTENT(OUT) :: PCF ! in the pressure equation
REAL, DIMENSION(:),     INTENT(OUT) :: PTRIGSX ! Arrays for sinus or cosinus 
REAL, DIMENSION(:),     INTENT(OUT) :: PTRIGSY ! values for the FFT in x and
                                               ! y directions
INTEGER, DIMENSION(:),  INTENT(OUT) :: KIFAXX  ! Decomposition in prime numbers
INTEGER, DIMENSION(:),  INTENT(OUT) :: KIFAXY  ! for the FFT in x and y
                                               ! direction      
INTEGER          ,  INTENT(OUT)  :: KALBOT     ! Vertical index corresponding
                                               ! to the absorbing layer base
!  
REAL, DIMENSION(:),   INTENT(OUT) :: PALK  ! Function of the absorbing
                                           ! layer damping coefficient
                                           ! defined for  u,v,and theta
REAL, DIMENSION(:),   INTENT(OUT) :: PALKW ! Idem but defined for w    
INTEGER          ,  INTENT(OUT)  :: KALBAS     ! Vertical index corresponding
                                               ! to the absorbing layer base
!  
REAL, DIMENSION(:),   INTENT(OUT) :: PALKBAS  ! Function of the absorbing
                                           ! layer damping coefficient
                                           ! defined for  u,v,and theta
REAL, DIMENSION(:),   INTENT(OUT) :: PALKWBAS ! Idem but defined for w
LOGICAL, DIMENSION(:,:),  INTENT(OUT)   :: OMASK_RELAX  ! True where the 
                                           ! lateral relax. has to be performed
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKURELAX  !  Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKVRELAX  !  coefficients for the
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKWRELAX  ! u, v and mass locations
REAL,                 INTENT(OUT) :: PDK2U ! 2nd order num. diffusion coef. /dx2
REAL,                 INTENT(OUT) :: PDK4U ! 4nd order num. diffusion coef. /dx4
                                           ! for momentum
REAL,                 INTENT(OUT) :: PDK2TH! for meteorological scalar variables
REAL,                 INTENT(OUT) :: PDK4TH! 
REAL,                 INTENT(OUT) :: PDK2SV! for tracer scalar variables
REAL,                 INTENT(OUT) :: PDK4SV! 
!
LOGICAL, INTENT(IN) :: OZDIFFU
TYPE(TYPE_ZDIFFU_HALO2)                       :: PZDIFFU_HALO2
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFB ! elements of the tri-diag matrix
                                            ! on an b-slice of global physical domain
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBF_SXP2_YP1_Z ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
END SUBROUTINE INI_DYNAMICS
!
END INTERFACE
!
END MODULE MODI_INI_DYNAMICS
!     ######################################################################
SUBROUTINE INI_DYNAMICS(PLON,PLAT,PRHODJ,PTHVREF,PMAP,PZZ,                   &
               PDXHAT,PDYHAT,PZHAT,HLBCX,HLBCY,PTSTEP,                       &
               OVE_RELAX,OVE_RELAX_GRD,OHORELAX_UVWTH,OHORELAX_RV,           &
               OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,  &
               OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                         &
               OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,&
               OHORELAX_SVCHEM,OHORELAX_SVAER,OHORELAX_SVDST,OHORELAX_SVSLT, &
               OHORELAX_SVPP,OHORELAX_SVCS,  OHORELAX_SVCHIC,OHORELAX_SVSNW, &
#ifdef MNH_FOREFIRE
               OHORELAX_SVFF,                                                &
#endif
               PRIMKMAX,KRIMX,KRIMY,PALKTOP,PALKGRD,PALZBOT,PALZBAS,         &
               PT4DIFU,PT4DIFTH,PT4DIFSV,                                    &
               PCORIOX,PCORIOY,PCORIOZ,PCURVX,PCURVY,                        &
               PDXHATM,PDYHATM,PRHOM,PAF,PBFY,PCF,                           &
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                                &
               PALK,PALKW,KALBOT,PALKBAS,PALKWBAS,KALBAS,                    &
               OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX,                     &
               PDK2U,PDK4U,PDK2TH,PDK4TH,PDK2SV,PDK4SV,OZDIFFU,PZDIFFU_HALO2,& 
               PBFB,&
               PBF_SXP2_YP1_Z)  !JUAN Z_SPLITING  
!     ######################################################################
!
!!****  *INI_DYNAMICS* - routine to initialize the parameters for the dynamics
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set or compute the parameters used
!     by the MESONH dynamics : 
!              * Coriolis parameters
!              * Curvature coefficients
!              * Pressure solver coefficients
!              * Absorbing layer coefficients
!              * Numerical difussion coefficients
!
!!**  METHOD
!!    ------
!!      - Coriolis parameters and curvature terms :
!!           Horizontal Coriolis parameters are not initialized if thinshell
!!    approximation is made (LTHINSHELL=.TRUE.).
!!           Curvature coefficients are not initialized if Cartesian geometry
!!    (LCARTESIAN=.TRUE.)  
!!      - Coefficients and variables for pressure solver :
!!          This is done by TRID
!!      - Coefficients and variables for the absorbing layer 
!!          ( upper and lateral) : This is done by RELAXDEF
!!      - Coefficients for the numerical diffusion
!!      
!!    EXTERNAL
!!    --------   
!!      TRID    : to initialize pressure solver
!!      RELAXDEF: to compute the relaxation coefficients
!!      GET_DIM_EXT_ll : get extended sub-domain sizes
!!
!!     Module MODI_TRID   : interface for routine TRID
!!     Module MODI_RELAXDEF : interface for routine RELAXDEF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF  : contains declaration of configuration variables
!!        
!!         LTHINSHELL : Logical for  THINSHELL approximation                
!!                     .TRUE.  = thinshell approximation
!!         LCARTESIAN : Logical for cartesian geometry :
!!                     .TRUE.  = cartesian geometry 
!!         L1D        : Logical for 1D configuration :
!!                     .TRUE.  = 1D model 
!!
!!      Module MODD_CST   : contains physical constants 
!!
!!         XPI        : Pi
!!         XOMEGA     : Earth rotation
!!
!!      Module MODD_GRID  : contains grid variables
!!
!!         XLON0      :  Reference longitude for the conformal projection
!!         XLAT0      :  Reference latitude for the conformal projection  
!!         XBETA      :  Rotation angle for the conformal projection
!!         XRPK       :  Projection parameter for the conformal projection
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_DYNAMICS)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        01/07/94 
!!      Modification    18/10/94 (J. Stein) to add the abs. layer
!!      Modification    16/11/94 (Lafore+Pinty) to add the num. diffusion
!!      Modification    06/12/94  (J.Stein) add the switch LABSLAYER
!!      Modification    12/12/94  (J.Stein) add the lateral relaxation 
!!      Modification    16/01/95  (J.Stein) conditional CALL to trid for 1D case
!!      Modification    13/08/98  (N.Asencio) add parallel code
!!      Modification    20/05/06  Remove KEPS
!!      Modification    07/2013   (Bosseur & Filippi) Adds Forefire
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!     Vionnet 07/2017 : blow snow
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_CONF
USE MODD_CST
USE MODD_GRID
USE MODD_LUNIT_n, ONLY: TLUOUT
!
USE MODI_RELAXDEF
! USE MODI_TRID
USE MODI_TRIDZ
USE MODI_ZDIFFUSETUP
!
USE MODE_ll
USE MODE_TYPE_ZDIFFU
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!  intent in arguments
!
REAL, DIMENSION(:,:),   INTENT(IN)        :: PLON,PLAT !Longitude and latitude
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PRHODJ    ! rho J
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PTHVREF   ! virtual potential 
                                      ! temperature of the reference state
REAL, DIMENSION(:,:),   INTENT(IN)        :: PMAP      ! Map factor
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ       ! height 
REAL, DIMENSION(:),     INTENT(IN)        :: PDXHAT     ! Stretching in x direction
REAL, DIMENSION(:),     INTENT(IN)        :: PDYHAT     ! Stretching in y direction
REAL, DIMENSION(:),     INTENT(IN)        :: PZHAT     ! Gal-Chen Height   
CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HLBCX    ! x-direction LBC type
CHARACTER(LEN=4), DIMENSION(:), INTENT(IN) :: HLBCY    ! y-direction LBC type
LOGICAL,                INTENT(IN)        :: OVE_RELAX ! logical
             ! switch to activate the VErtical  RELAXation 
LOGICAL,                INTENT(IN)        :: OVE_RELAX_GRD ! logical
             ! switch to activate the VErtical  RELAXation (ground layer)       
LOGICAL,            INTENT(IN) :: OHORELAX_UVWTH  ! switch for the 
                       ! horizontal relaxation for U,V,W,TH
LOGICAL,            INTENT(IN) :: OHORELAX_RV     ! switch for the 
                       ! horizontal relaxation for Rv
LOGICAL,            INTENT(IN) :: OHORELAX_RC     ! switch for the 
                       ! horizontal relaxation for Rc
LOGICAL,            INTENT(IN) :: OHORELAX_RR     ! switch for the 
                       ! horizontal relaxation for Rr
LOGICAL,            INTENT(IN) :: OHORELAX_RI     ! switch for the 
                       ! horizontal relaxation for Ri
LOGICAL,            INTENT(IN) :: OHORELAX_RS     ! switch for the 
                       ! horizontal relaxation for Rs
LOGICAL,            INTENT(IN) :: OHORELAX_RG     ! switch for the 
                       ! horizontal relaxation for Rg
LOGICAL,            INTENT(IN) :: OHORELAX_RH     ! switch for the 
                       ! horizontal relaxation for Rh
LOGICAL,            INTENT(IN) :: OHORELAX_TKE    ! switch for the 
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the 
                       ! horizontal relaxation for sv variables
LOGICAL,             INTENT(IN):: OHORELAX_SVC2R2 ! switch for the 
                       ! horizontal relaxation for c2r2 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVC1R3 ! switch for the 
                       ! horizontal relaxation for c1r3 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVELEC ! switch for the 
                       ! horizontal relaxation for elec variables
LOGICAL,             INTENT(IN):: OHORELAX_SVLG   ! switch for the 
                       ! horizontal relaxation for lg variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHEM ! switch for the 
                       ! horizontal relaxation for chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHIC ! switch for the
                       ! horizontal relaxation for ice chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVAER  ! switch for the 
                       ! horizontal relaxation for aer variables
LOGICAL,             INTENT(IN):: OHORELAX_SVDST  ! switch for the 
                       ! horizontal relaxation for dst variables
LOGICAL,             INTENT(IN):: OHORELAX_SVSLT  ! switch for the 
                       ! horizontal relaxation for slt variables
LOGICAL,             INTENT(IN):: OHORELAX_SVPP   ! switch for the 
                       ! horizontal relaxation for passive pollutants
LOGICAL,             INTENT(IN):: OHORELAX_SVSNW   ! switch for the
                       ! horizontal relaxation for blowing snow variables
#ifdef MNH_FOREFIRE
LOGICAL,             INTENT(IN):: OHORELAX_SVFF   ! switch for the 
                       ! horizontal relaxation for ForeFire variables
#endif
LOGICAL,             INTENT(IN):: OHORELAX_SVCS   ! switch for the 
                       ! horizontal relaxation for conditional sampling
REAL,                    INTENT(IN)    :: PRIMKMAX !Max. value of the horiz.
                                     ! relaxation coefficients
INTEGER,                INTENT(IN)     :: KRIMX,KRIMY ! Number of points in 
                                     ! the rim zone in the x and y directions 
REAL,     INTENT(IN)   :: PALKTOP    ! Damping coef. at the top of the absorbing
                                     ! layer
REAL,     INTENT(IN)   :: PALZBOT    ! Height of the absorbing layer base
REAL,     INTENT(IN)   :: PALKGRD    ! Damping coef. at the top of the absorbing
                                     ! layer
REAL,     INTENT(IN)   :: PALZBAS    ! Height of the absorbing layer base
REAL,     INTENT(IN)   :: PT4DIFU    ! Damping time scale for 2*dx wavelength
                                     ! specified for the 4nd order num. diffusion
                                     ! for momentum
REAL,     INTENT(IN)   :: PT4DIFTH   ! for meteorological scalar variables
REAL,     INTENT(IN)   :: PT4DIFSV   ! for tracer scalar variables

REAL,     INTENT(IN)   :: PTSTEP     ! Time step    
!
!  intent out arguments
!
REAL, INTENT(OUT) :: PDXHATM                     ! mean grid increment in the x
                                                 ! direction
REAL, INTENT(OUT) :: PDYHATM                     ! mean grid increment in the y
                                                 ! direction
!
REAL, DIMENSION (:), INTENT(OUT) :: PRHOM        !  mean of XRHODJ on the plane x y 
                                                 !  localized at a mass level
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PCORIOX,PCORIOY ! Hor. Coriolis parameters
REAL, DIMENSION(:,:), INTENT(OUT) :: PCORIOZ         ! Vert. Coriolis parameter 
REAL, DIMENSION(:,:), INTENT(OUT) :: PCURVX,PCURVY   ! Curvature coefficients
! 
REAL, DIMENSION(:),     INTENT(OUT) :: PAF ! vectors giving the non-vanishing        
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFY ! elements of the tri-diag matrix
                                            ! on an y-slice of global physical domain
REAL, DIMENSION(:),     INTENT(OUT) :: PCF ! in the pressure equation
REAL, DIMENSION(:),     INTENT(OUT) :: PTRIGSX ! Arrays for sinus or cosinus 
REAL, DIMENSION(:),     INTENT(OUT) :: PTRIGSY ! values for the FFT in x and
                                               ! y directions
INTEGER, DIMENSION(:),  INTENT(OUT) :: KIFAXX  ! Decomposition in prime numbers
INTEGER, DIMENSION(:),  INTENT(OUT) :: KIFAXY  ! for the FFT in x and y
                                               ! direction      
INTEGER          ,  INTENT(OUT)  :: KALBOT     ! Vertical index corresponding
                                               ! to the absorbing layer base
!  
REAL, DIMENSION(:),   INTENT(OUT) :: PALK  ! Function of the absorbing
                                           ! layer damping coefficient
                                           ! defined for  u,v,and theta
REAL, DIMENSION(:),   INTENT(OUT) :: PALKW ! Idem but defined for w     
INTEGER          ,  INTENT(OUT)  :: KALBAS     ! Vertical index corresponding
                                               ! to the absorbing layer base
!  
REAL, DIMENSION(:),   INTENT(OUT) :: PALKBAS  ! Function of the absorbing
                                           ! layer damping coefficient
                                           ! defined for  u,v,and theta
REAL, DIMENSION(:),   INTENT(OUT) :: PALKWBAS ! Idem but defined for w
LOGICAL, DIMENSION(:,:),  INTENT(OUT)   :: OMASK_RELAX  ! True where the 
                                           ! lateral relax. has to be performed
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKURELAX  !  Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKVRELAX  !  coefficients for the
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKWRELAX  ! u, v and mass locations
REAL,                 INTENT(OUT) :: PDK2U ! 2nd order num. diffusion coef. /dx2
REAL,                 INTENT(OUT) :: PDK4U ! 4nd order num. diffusion coef. /dx4
                                           ! for momentum
REAL,                 INTENT(OUT) :: PDK2TH! for meteorological scalar variables
REAL,                 INTENT(OUT) :: PDK4TH! 
REAL,                 INTENT(OUT) :: PDK2SV! for tracer scalar variables
REAL,                 INTENT(OUT) :: PDK4SV! 
LOGICAL, INTENT(IN) :: OZDIFFU
TYPE(TYPE_ZDIFFU_HALO2)                       :: PZDIFFU_HALO2
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFB ! elements of the tri-diag matrix
                                            ! on an b-slice of global physical domain
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBF_SXP2_YP1_Z ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
!
!*       0.2   declarations of local variables
!
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2)) :: ZGAMMA ! Gamma =K(lambda-lambda0) - beta
REAL                                       :: ZMBETA ! -beta
REAL                                       :: ZCDR   ! to convert degrees in
                                                     ! radians 
INTEGER                                    :: ILUOUT ! Logical unit number for output_listing file
INTEGER                                    :: IIU,IJU !  Upper bounds in x,y directions
LOGICAL                                    :: GHORELAX
LOGICAL, DIMENSION(7) :: GHORELAXR ! local array of logical
#ifdef MNH_FOREFIRE
LOGICAL, DIMENSION(13):: GHORELAXSV! local array of logical
#else
LOGICAL, DIMENSION(12):: GHORELAXSV! local array of logical
#endif
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTES CORIOLIS PARAMETERS AND CURVATURE COEFFICIENTS
!              -------------------------------------------------------
! 
ZCDR = XPI/180.
IF (.NOT.LCARTESIAN) THEN
  ZGAMMA(:,:) = XRPK * ((PLON(:,:) - XLON0)*ZCDR) - (XBETA*ZCDR)
  IF (.NOT.LTHINSHELL) THEN    
    PCORIOX(:,:) = - 2. * XOMEGA * COS(PLAT(:,:)*ZCDR) * SIN(ZGAMMA(:,:))
    PCORIOY(:,:) =   2. * XOMEGA * COS(PLAT(:,:)*ZCDR) * COS(ZGAMMA(:,:))
  END IF
  PCORIOZ(:,:)   =   2. * XOMEGA * SIN(PLAT(:,:)*ZCDR) 
  PCURVX (:,:)   = COS(ZGAMMA(:,:)) * (SIN(PLAT(:,:)*ZCDR) -XRPK)         &
                 / COS(PLAT(:,:)*ZCDR)
  PCURVY (:,:)   = SIN(ZGAMMA(:,:)) * (SIN(PLAT(:,:)*ZCDR) -XRPK)         &
                 / COS(PLAT(:,:)*ZCDR)
ELSE
  ZMBETA = - (XBETA*ZCDR)
  PCORIOX(:,:) = - 2. * XOMEGA * COS(XLAT0*ZCDR) * SIN(ZMBETA)
  PCORIOY(:,:) =   2. * XOMEGA * COS(XLAT0*ZCDR) * COS(ZMBETA)
  PCORIOZ(:,:) =   2. * XOMEGA * SIN(XLAT0*ZCDR) 
END IF           
!
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZATION OF PRESSURE SOLVER
!              ---------------------------------
!
IF (.NOT.L1D) THEN
! CALL TRID(HLBCX,HLBCY,                                                  &
!           PMAP,PDXHAT,PDYHAT,PDXHATM,PDYHATM,PRHOM,PAF,                 &
!           PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                            &
!           PRHODJ,PTHVREF,PZZ,PBFY)
  CALL TRIDZ(HLBCX,HLBCY,                                                 &
            PMAP,PDXHAT,PDYHAT,PDXHATM,PDYHATM,PRHOM,PAF,                 &
            PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                            &
            PRHODJ,PTHVREF,PZZ,PBFY,PBFB,                                 &
            PBF_SXP2_YP1_Z)
END IF
!
!
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE ABSORBING LAYER COEFFICIENTS
!              ----------------------------------------
!
GHORELAXR(1) = OHORELAX_RV
GHORELAXR(2) = OHORELAX_RC
GHORELAXR(3) = OHORELAX_RR
GHORELAXR(4) = OHORELAX_RI
GHORELAXR(5) = OHORELAX_RS
GHORELAXR(6) = OHORELAX_RG
GHORELAXR(7) = OHORELAX_RH
!
GHORELAXSV(1) = OHORELAX_SVC2R2
GHORELAXSV(2) = OHORELAX_SVC1R3
GHORELAXSV(3) = OHORELAX_SVELEC
GHORELAXSV(4) = OHORELAX_SVLG
GHORELAXSV(5) = OHORELAX_SVCHEM
GHORELAXSV(6) = OHORELAX_SVAER
GHORELAXSV(7) = OHORELAX_SVDST
GHORELAXSV(8) = OHORELAX_SVSLT
GHORELAXSV(9) = OHORELAX_SVPP
GHORELAXSV(10)= OHORELAX_SVCS
GHORELAXSV(11) = OHORELAX_SVCHIC
GHORELAXSV(12) = OHORELAX_SVSNW
#ifdef MNH_FOREFIRE
GHORELAXSV(13) = OHORELAX_SVFF
#endif
!
GHORELAX=ANY(GHORELAXR) .OR. ANY(GHORELAXSV) .OR. ANY(OHORELAX_SV) &
                        .OR. OHORELAX_UVWTH  .OR. OHORELAX_TKE 
!
IF (GHORELAX .OR. OVE_RELAX.OR.OVE_RELAX_GRD) THEN
  CALL RELAXDEF( OVE_RELAX,OVE_RELAX_GRD,OHORELAX_UVWTH,OHORELAX_RV,  &
     OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,     &
     OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                            &
     OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,   &
     OHORELAX_SVCHEM, OHORELAX_SVAER, OHORELAX_SVDST, OHORELAX_SVSLT, &
     OHORELAX_SVPP, OHORELAX_SVCS, OHORELAX_SVCHIC,OHORELAX_SVSNW,    &
     PALKTOP,PALKGRD, PALZBOT,PALZBAS,                                &
     PZZ, PZHAT, PTSTEP,                                              &
     PRIMKMAX,KRIMX,KRIMY,                                            &
     PALK, PALKW, KALBOT,                                             &
     PALKBAS, PALKWBAS, KALBAS,                                       &
     OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX                         )
END IF
!
!
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTE THE NUMERICAL DIFFUSION COEFFICIENTS
!              --------------------------------------------
!
PDK4U = 1.0/(16.0*PT4DIFU) ! The damping rate for the 2*dx wavelength is the same
PDK2U = 2.0*PDK4U          ! for the 2nd and the 4th order diffusion schemes
                           ! for momentum
PDK4TH= 1.0/(16.0*PT4DIFTH) ! for meteorological scalar variables                       
PDK2TH= 2.0*PDK4TH         
PDK4SV= 1.0/(16.0*PT4DIFSV) ! for tracer scalar variables                       
PDK2SV= 2.0*PDK4SV         
!
! Call ZDIFFUSETUP if OZDIFFU is true (parameters for truly horizontal diffusion)
!
IF (OZDIFFU) THEN
  CALL ZDIFFUSETUP (PZZ,&
                    PZDIFFU_HALO2)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.    PRINT ON OUTPUT_LISTING
!              -----------------------
!
IF (NVERB >= 10) THEN
  CALL GET_DIM_EXT_ll ('B',IIU,IJU)
  ILUOUT = TLUOUT%NLU
!
  WRITE(ILUOUT,*) 'INI_DYNAMICS : Some PCORIOZ values'
  WRITE(ILUOUT,*) '(1,1)          (IIU/2,IJU/2)         (IIU,IJU)  '
  WRITE(ILUOUT,*) PCORIOZ(1,1),PCORIOZ(IIU/2,IJU/2),PCORIOZ(IIU,IJU) 
!
  IF (.NOT.LTHINSHELL) THEN
    WRITE(ILUOUT,*) 'INI_DYNAMICS : Some PCORIOX values'
    WRITE(ILUOUT,*) '(1,1)          (IIU/2,IJU/2)         (IIU,IJU)  '
    WRITE(ILUOUT,*) PCORIOX(1,1),PCORIOX(IIU/2,IJU/2),PCORIOX(IIU,IJU) 
!
    WRITE(ILUOUT,*) 'INI_DYNAMICS : Some PCORIOY values'
    WRITE(ILUOUT,*) '(1,1)          (IIU/2,IJU/2)         (IIU,IJU)  '
    WRITE(ILUOUT,*) PCORIOY(1,1),PCORIOY(IIU/2,IJU/2),PCORIOY(IIU,IJU) 
  END IF
!
  IF ( .NOT. LCARTESIAN ) THEN
    WRITE(ILUOUT,*) 'INI_DYNAMICS : Some PCURVX values'
    WRITE(ILUOUT,*) '(1,1)          (IIU/2,IJU/2)         (IIU,IJU)  '
    WRITE(ILUOUT,*) PCURVX(1,1),PCURVX(IIU/2,IJU/2),PCURVX(IIU,IJU)
! 
    WRITE(ILUOUT,*) 'INI_DYNAMICS : Some PCURVY values'
    WRITE(ILUOUT,*) '(1,1)          (IIU/2,IJU/2)         (IIU,IJU)  '
    WRITE(ILUOUT,*) PCURVY(1,1),PCURVY(IIU/2,IJU/2),PCURVY(IIU,IJU) 
  END IF
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_DYNAMICS
