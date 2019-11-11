!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_INI_DYNAMICS
!     ########################
INTERFACE
SUBROUTINE INI_DYNAMICS(HLUOUT,PLON,PLAT,PRHODJ,PTHVREF,PMAP,PZZ,           &
               PDXHAT,PDYHAT,PZHAT,HLBCX,HLBCY,PTSTEP,                      &
               OVE_RELAX,OHORELAX_UVWTH,OHORELAX_RV,                        &
               OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG, &
               OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                        &
               OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,&
               OHORELAX_SVCHEM,OHORELAX_SVAER,OHORELAX_SVDST,OHORELAX_SVSLT,&
               PRIMKMAX,KRIMX,KRIMY,PALKTOP,PALZBOT,                        &
               PT4DIFF,                                                     &
               PCORIOX,PCORIOY,PCORIOZ,PCURVX,PCURVY,                       &
               PDXHATM,PDYHATM,PRHOM,PAF,PBFY,PCF,                          &
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                               &
               PALK,PALKW,KALBOT,OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX,  &
               PDK2,PDK4,OZDIFFU,PZDIFFU_HALO2,                             &
               PBFB,&
               PBF_SXP2_YP1_Z) !JUAN Z_SPLITING
!  intent in arguments
!
USE MODE_TYPE_ZDIFFU
IMPLICIT NONE
!
CHARACTER(LEN=*),       INTENT(IN)        :: HLUOUT    ! name for output-listing
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
LOGICAL,             INTENT(IN):: OHORELAX_SVAER  ! switch for the 
                       ! horizontal relaxation for aer variables
LOGICAL,             INTENT(IN):: OHORELAX_SVDST  ! switch for the 
                       ! horizontal relaxation for dst variables
LOGICAL,             INTENT(IN):: OHORELAX_SVSLT  ! switch for the 
                       ! horizontal relaxation for slt variables
REAL,                    INTENT(IN)    :: PRIMKMAX !Max. value of the horiz.
                                     ! relaxation coefficients
INTEGER,                INTENT(IN)     :: KRIMX,KRIMY ! Number of points in 
                                     ! the rim zone in the x and y directions 
REAL,     INTENT(IN)   :: PALKTOP    ! Damping coef. at the top of the absorbing
                                     ! layer
REAL,     INTENT(IN)   :: PALZBOT    ! Height of the absorbing layer base
REAL,     INTENT(IN)   :: PT4DIFF    ! Damping time scale for 2*dx wavelength
                                     ! specified for the 4nd order num. diffusion

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
LOGICAL, DIMENSION(:,:),  INTENT(OUT)   :: OMASK_RELAX  ! True where the 
                                           ! lateral relax. has to be performed
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKURELAX  !  Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKVRELAX  !  coefficients for the
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKWRELAX  ! u, v and mass locations
REAL,                 INTENT(OUT) :: PDK2  ! 2nd order num. diffusion coef. /dx2
REAL,                 INTENT(OUT) :: PDK4  ! 4nd order num. diffusion coef. /dx4
!
LOGICAL, INTENT(IN) :: OZDIFFU
TYPE(TYPE_ZDIFFU_HALO2)                       :: PZDIFFU_HALO2
!
!JUAN Z_SPLITING
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFB ! elements of the tri-diag matrix
                                            ! on an b-slice of global physical domain
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBF_SXP2_YP1_Z ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
!JUAN Z_SPLITING
!
END SUBROUTINE INI_DYNAMICS
!
END INTERFACE
!
END MODULE MODI_INI_DYNAMICS
!
!
!
!     ######################################################################
SUBROUTINE INI_DYNAMICS(HLUOUT,PLON,PLAT,PRHODJ,PTHVREF,PMAP,PZZ,           &
               PDXHAT,PDYHAT,PZHAT,HLBCX,HLBCY,PTSTEP,                      &
               OVE_RELAX,OHORELAX_UVWTH,OHORELAX_RV,                        &
               OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG, &
               OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                        &
               OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,&
               OHORELAX_SVCHEM,OHORELAX_SVAER,OHORELAX_SVDST,OHORELAX_SVSLT,&
               PRIMKMAX,KRIMX,KRIMY,PALKTOP,PALZBOT,                        &
               PT4DIFF,                                                     &
               PCORIOX,PCORIOY,PCORIOZ,PCURVX,PCURVY,                       &
               PDXHATM,PDYHATM,PRHOM,PAF,PBFY,PCF,                          &
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                               &
               PALK,PALKW,KALBOT,OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX,  &
               PDK2,PDK4,OZDIFFU,PZDIFFU_HALO2,                             &
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
!!      FMLOOK  : to retrieve logical unit number linked to a file
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_CONF
USE MODD_CST
USE MODD_GRID
!
USE MODI_TRID
!JUAN Z_SPLITING
USE MODI_TRIDZ
!JUAN Z_SPLITING
USE MODI_RELAXDEF
USE MODI_ZDIFFUSETUP
USE MODE_ll
USE MODE_FM
!
USE MODE_TYPE_ZDIFFU
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!  intent in arguments
!
CHARACTER(LEN=*),       INTENT(IN)        :: HLUOUT    ! name for output-listing
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
LOGICAL,             INTENT(IN):: OHORELAX_SVAER  ! switch for the 
                       ! horizontal relaxation for aer variables
LOGICAL,             INTENT(IN):: OHORELAX_SVDST  ! switch for the 
                       ! horizontal relaxation for dst variables
LOGICAL,             INTENT(IN):: OHORELAX_SVSLT  ! switch for the 
                       ! horizontal relaxation for slt variables
REAL,                    INTENT(IN)    :: PRIMKMAX !Max. value of the horiz.
                                     ! relaxation coefficients
INTEGER,                INTENT(IN)     :: KRIMX,KRIMY ! Number of points in 
                                     ! the rim zone in the x and y directions 
REAL,     INTENT(IN)   :: PALKTOP    ! Damping coef. at the top of the absorbing
                                     ! layer
REAL,     INTENT(IN)   :: PALZBOT    ! Height of the absorbing layer base
REAL,     INTENT(IN)   :: PT4DIFF    ! Damping time scale for 2*dx wavelength
                                     ! specified for the 4nd order num. diffusion

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
LOGICAL, DIMENSION(:,:),  INTENT(OUT)   :: OMASK_RELAX  ! True where the 
                                           ! lateral relax. has to be performed
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKURELAX  !  Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKVRELAX  !  coefficients for the
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKWRELAX  ! u, v and mass locations
REAL,                 INTENT(OUT) :: PDK2  ! 2nd order num. diffusion coef. /dx2
REAL,                 INTENT(OUT) :: PDK4  ! 4nd order num. diffusion coef. /dx4
LOGICAL, INTENT(IN) :: OZDIFFU
TYPE(TYPE_ZDIFFU_HALO2)                       :: PZDIFFU_HALO2
!JUAN Z_SPLITING
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFB ! elements of the tri-diag matrix
                                            ! on an b-slice of global physical domain
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBF_SXP2_YP1_Z ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
!JUAN Z_SPLITING
!
!*       0.2   declarations of local variables
!
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2)) :: ZGAMMA ! Gamma =K(lambda-lambda0) - beta
REAL                                       :: ZMBETA ! -beta
REAL                                       :: ZCDR   ! to convert degrees in
                                                     ! radians 
INTEGER                                    :: ILUOUT,IRESP ! Logical unit number
                                      !  for output_listing file, and return code
INTEGER                                    :: IIU,IJU !  Upper bounds in x,y directions
LOGICAL                                    :: GHORELAX
LOGICAL, DIMENSION(7) :: GHORELAXR ! local array of logical
LOGICAL, DIMENSION(8) :: GHORELAXSV! local array of logical
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
!!$  CALL TRID(HLUOUT,HLBCX,HLBCY,                                           &
!!$            PMAP,PDXHAT,PDYHAT,PDXHATM,PDYHATM,PRHOM,PAF,                 &
!!$            PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                            &
!!$            PRHODJ,PTHVREF,PZZ,PBFY)
!JUAN Z_SPLITING
  CALL TRIDZ(HLUOUT,HLBCX,HLBCY,                                           &
            PMAP,PDXHAT,PDYHAT,PDXHATM,PDYHATM,PRHOM,PAF,                 &
            PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                            &
            PRHODJ,PTHVREF,PZZ,PBFY,PBFB,&
            PBF_SXP2_YP1_Z)
!JUAN Z_SPLITING
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
!
GHORELAX=ANY(GHORELAXR) .OR. ANY(GHORELAXSV) .OR. ANY(OHORELAX_SV) &
                        .OR. OHORELAX_UVWTH  .OR. OHORELAX_TKE 
!
IF (GHORELAX .OR. OVE_RELAX) THEN
  CALL RELAXDEF( OVE_RELAX,OHORELAX_UVWTH,OHORELAX_RV,                &
     OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,     &
     OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                            &
     OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,   &
     OHORELAX_SVCHEM, OHORELAX_SVAER, OHORELAX_SVDST, OHORELAX_SVSLT, &
     PALKTOP, PALZBOT,                                                &
     PZZ, PZHAT, PTSTEP,                                              &
     PRIMKMAX,KRIMX,KRIMY,                                            &
     PALK, PALKW, KALBOT,                                             &
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
PDK4 = 1.0/(16.0*PT4DIFF) ! The damping rate for the 2*dx wavelength is the same
PDK2 = 2.0*PDK4           ! for the 2nd and the 4th order diffusion schemes
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
  CALL FMLOOK_ll(HLUOUT,HLUOUT,ILUOUT,IRESP)
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
