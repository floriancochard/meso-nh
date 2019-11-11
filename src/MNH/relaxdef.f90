!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ####################
      MODULE MODI_RELAXDEF
!     ####################
!
INTERFACE
!
      SUBROUTINE RELAXDEF ( OVE_RELAX,OVE_RELAX_GRD,                         &
            OHORELAX_UVWTH,OHORELAX_RV,                                      &
            OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,     &
            OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                            &
            OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,   &
            OHORELAX_SVCHEM, OHORELAX_SVAER, OHORELAX_SVDST, OHORELAX_SVSLT, &
            OHORELAX_SVPP,OHORELAX_SVCS, OHORELAX_SVCHIC,OHORELAX_SVSNW,     &
            PALKTOP,PALKGRD, PALZBOT,PALZBAS,                                &
            PZZ, PZHAT, PTSTEP,                                              &
            PRIMKMAX,KRIMX, KRIMY,                                           &
            PALK, PALKW, KALBOT,PALKBAS,PALKWBAS,KALBAS,                     &  
            OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX    )
!
LOGICAL,                INTENT(IN)        :: OVE_RELAX ! logical
             ! switch to activate the VErtical  RELAXation 
LOGICAL,                INTENT(IN)        :: OVE_RELAX_GRD ! logical
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
LOGICAL,             INTENT(IN):: OHORELAX_SVCS   ! switch for the 
                       ! horizontal relaxation for conditional sampling 
LOGICAL,             INTENT(IN):: OHORELAX_SVSNW  ! switch for the
                       ! horizontal relaxation for blowing snow variables                        
REAL,                   INTENT(IN)        :: PALKTOP   ! Damping coef. at the 
                                                       ! top of the abs. layer
REAL,                   INTENT(IN)        :: PALKGRD   ! Damping coef. at the 
                                                       ! top of the abs. layer
REAL,                   INTENT(IN)        :: PALZBOT   ! Height of the abs.
                                                       ! layer base
REAL,                   INTENT(IN)        :: PALZBAS   ! Height of the abs.
                                                       ! layer base            
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ       ! Height 
REAL, DIMENSION(:),     INTENT(IN)        :: PZHAT     ! Gal-Chen Height 
REAL,                   INTENT(IN)        :: PTSTEP    ! Time step         
REAL,                     INTENT(IN)    :: PRIMKMAX !Max. value of the horiz.
                                          ! relaxation coefficients
INTEGER,                INTENT(IN)        :: KRIMX,KRIMY ! Number of points in 
                                       ! the rim zone in the x and y directions 
!
REAL, DIMENSION(:),     INTENT(OUT)       :: PALK      ! Function of the damping
                                                       ! coef. defined for u, v
                                                       ! and theta
REAL, DIMENSION(:),     INTENT(OUT)       :: PALKW     ! Idem but defined for w 
INTEGER,                INTENT(OUT)       :: KALBOT    ! Vertical index corres- 
                                                       ! -ponding to the abs. 
                                                       ! layer base
!
REAL, DIMENSION(:),     INTENT(OUT)       :: PALKBAS      ! Function of the damping
                                                       ! coef. defined for u, v
                                                       ! and theta
REAL, DIMENSION(:),     INTENT(OUT)       :: PALKWBAS     ! Idem but defined for w 
INTEGER,                INTENT(OUT)       :: KALBAS    ! Vertical index corres- 
                                                       ! -ponding to the abs. 
                                                       ! layer base
!                                                       
LOGICAL, DIMENSION(:,:),  INTENT(OUT)   :: OMASK_RELAX  ! True where the 
                                         !  lateral relax. has to be performed
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKURELAX !  Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKVRELAX !  coefficients for the
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKWRELAX ! u, v and mass locations
!
END SUBROUTINE RELAXDEF
!
END INTERFACE
!
END MODULE MODI_RELAXDEF
!     ######spl
      SUBROUTINE RELAXDEF ( OVE_RELAX,OVE_RELAX_GRD,                         &
            OHORELAX_UVWTH,OHORELAX_RV,                                      &
            OHORELAX_RC,OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,     &
            OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,                            &
            OHORELAX_SVC2R2,OHORELAX_SVC1R3,OHORELAX_SVELEC,OHORELAX_SVLG,   &
            OHORELAX_SVCHEM, OHORELAX_SVAER, OHORELAX_SVDST, OHORELAX_SVSLT, &
            OHORELAX_SVPP,OHORELAX_SVCS, OHORELAX_SVCHIC,OHORELAX_SVSNW,     &
            PALKTOP,PALKGRD, PALZBOT,PALZBAS,                                &
            PZZ, PZHAT, PTSTEP,                                              &
            PRIMKMAX,KRIMX, KRIMY,                                           &
            PALK, PALKW, KALBOT,PALKBAS,PALKWBAS,KALBAS,                     &
            OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX    )
!     #########################################################################
!
!!****  *RELAXDEF* - routine to set up the damping parameters for both
!!                   the top absorbing layer and the lateral rim zone
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the vertical level at which
!     begins the absorbing layer and the damping coefficient used within the
!     absorbing layer. The routine delivers also the arrays containing the
!     lateral relaxation coefficients for the staggered fields.
!
!!**  METHOD
!!    ------
!       The routine returns a vertical index corresponding to the base
!     of the absorbing layer (KALBOT) and two vectors PALK and PALKW 
!     which are function of the absorbing layer damping coefficient
!     Kal(zhat). PALK and PALKW are idendtical but defined on different
!     vertical grids.
!                                     Kal(zhat)
!        PALK, PALKW --------->  -------------------
!                                 1 + 2 dt Kal(zhat)
!        with
!                                 pi  zhat - zhatalb
!        Kal(zhat) = Kal (H) (sin --- --------------)**2   for zhat > zhatalb
!                                  2   H - zhatalb
!
!       The routine returns also the relaxation coefficients defined only
!     on a band around the boundary.
!                                           Krelax(i,j)
!        PKURELAX, PKVRELAX, PKWRELAX = --------------------
!                                       1 + 2 dt Krelax(i,j)
!        with
!                                    pi  d(i,j)-bound
!        Krelax(i,j) = rimKmax (cos --- --------------)**2
!                                    2   band width 
!
!        with "d(i,j)-bound" being the distance of the current location to
!     the closest boundary.
!     In the corners, d(i,j) represent the elliptical distance to the center of
!     an ellipse located at the intersection point between the innermost lines 
!     of the relaxation bands in the I and J directions.
!
!        The band width is given by the number of normal velocity points in the 
!     involved direction ( u for x direction and v for y direction ).
!     For example, for krimx = 4 :
!
!       u (Ku=rimKmax)  u               u               u               u(Ku=0)
!         
!       x               x               x               x               x              
!
!      (i=iib)                                                    (i=iib+krimx)
!                                   
!!      
!!    EXTERNAL
!!    --------   
!!    FUNCTION RELAX : compute the decrease of the lateral relaxation
!!    in function of the distance to the closest boundary
!!    GET_INDICE_ll       : get physical sub-domain bounds
!!    GET_GLOBALDIMS_ll   : get physical global domain sizes
!!    GET_OR_ll           : get origine coordonates of the physical sub-domain 
!!                          in global indices
!!    GET_INTERSECTION_ll : get the indices of the intersection inside the
!!                          extended global domain
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      Module MODD_CST   : contains physical constants 
!!         XPI        : pi
!!
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT, JPVEXT: define the number of marginal points out of the
!!        physical domain along horizontal and vertical directions respectively
!!       
!!      Module MODD_CONF
!!        L2D           : logical switch for 2 Dimensional verion
!!       
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine RELAXDEF)
!!      
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/12/94 
!!      J. Stein    12/12/94  add switch for horizontal and vertical computation
!!                            + bug in the loop bounds + add PKVRELAX and 
!!                            PKWRELAX
!!      J.Stein     10/01/95  change ZSBAR in ZSHAT
!!    J.P.Pinty     10/02/95  add a 2D switch 
!!      I. Mallet   18/03/96  introduce a function to compute lateral coeffients
!!      N. Asencio  10/08/98  add parallel code
!!      J. Escobar  09/03/2008 bug in extended  dim --> 2*JPHEXT
!!      V. Masson, C.Lac 09/2010 reproducibility : replacement of SUM3D_ll to SUMALL_ll
!!                             and of PZZ(IIB,IJB,IKE+1) to PZHAT(IKE+1)
!!      J.Escobar   30/09/2010  introduction of CPP MACRO(REAL16) for reproductibility test 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_CONF
USE MODI_RELAX
USE MODE_ll
!JUAN
USE MODE_REPRO_SUM
!JUAN
!
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!
LOGICAL,                INTENT(IN)        :: OVE_RELAX ! logical
             ! switch to activate the VErtical  RELAXation 
LOGICAL,                INTENT(IN)        :: OVE_RELAX_GRD ! logical
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
                       ! horizontal relaxation for sv
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
LOGICAL,             INTENT(IN):: OHORELAX_SVCS   ! switch for the 
                       ! horizontal relaxation for conditional sampling 
LOGICAL,             INTENT(IN):: OHORELAX_SVSNW  ! switch for the
                       ! horizontal relaxation for blowing snow variables                       
REAL,                   INTENT(IN)        :: PALKTOP   ! Damping coef. at the 
                                                       ! top of the abs. layer
REAL,                   INTENT(IN)        :: PALKGRD   ! Damping coef. at the 
                                                       ! top of the abs. layer                                                       
REAL,                   INTENT(IN)        :: PALZBOT   ! Height of the abs.
                                                       ! layer base 
REAL,                   INTENT(IN)        :: PALZBAS   ! Height of the abs.
                                                       ! layer base                                                        
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ       ! Height 
REAL, DIMENSION(:),     INTENT(IN)        :: PZHAT     ! Gal-Chen Height 
REAL,                   INTENT(IN)        :: PTSTEP    ! Time step         
REAL,                     INTENT(IN)    :: PRIMKMAX !Max. value of the horiz.
                                          ! relaxation coefficients
INTEGER,                INTENT(IN)        :: KRIMX,KRIMY ! Number of points in 
                                       ! the rim zone in the x and y directions 
!
REAL, DIMENSION(:),     INTENT(OUT)       :: PALK      ! Function of the damping
                                                       ! coef. defined for u, v
                                                       ! and theta
REAL, DIMENSION(:),     INTENT(OUT)       :: PALKW     ! Idem but defined for w 
INTEGER,                INTENT(OUT)       :: KALBOT    ! Vertical index corres- 
                                                       ! -ponding to the abs. 
                                                       ! layer base
!
REAL, DIMENSION(:),     INTENT(OUT)       :: PALKBAS      ! Function of the damping
                                                       ! coef. defined for u, v
                                                       ! and theta
REAL, DIMENSION(:),     INTENT(OUT)       :: PALKWBAS     ! Idem but defined for w 
INTEGER,                INTENT(OUT)       :: KALBAS    ! Vertical index corres- 
!
LOGICAL, DIMENSION(:,:),  INTENT(OUT)   :: OMASK_RELAX  ! True where the 
                                         !  lateral relax. has to be performed
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKURELAX !  Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKVRELAX !  coefficients for the
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PKWRELAX ! u, v and mass locations
!
!
!*       0.2   declarations of local variables
!
INTEGER, DIMENSION(1)      :: IWORK       ! Work array                 
REAL                       :: ZWORK       ! Work variable    
REAL                       :: ZZSHAT        
REAL                       :: ZZTOP       ! Height of the model top         
REAL                       :: ZALZHAT     ! Gal-Chen height of the abs. layer 
                                          ! base
REAL                       :: ZZHATK      ! Gal-Chen height of u(k), v(k), and 
                                          ! theta(k)
!
INTEGER                    :: IKRIMAX     ! Maximum width of the rim zone
                                          !  (number of points)
                                            ! Complete domain:
INTEGER                    :: IIB_ll,IJB_ll ! Lower bounds of the physical
                                            !  global domain in x and y directions
INTEGER                    :: IIE_ll,IJE_ll ! Upper bounds of the physical
                                            ! gloval domain in x and y directions
                                            ! sub-domain:
INTEGER                    :: IIB,IJB       ! Lower bounds of the physical
                                            ! sub-domain in x and y directions
INTEGER                    :: IIE,IJE,IKE   ! Upper bounds of the physical
                                            ! sub-domain in x,y and z directions
INTEGER                    :: IORX,IORY     ! origines of the extended sub-domain 
                                            ! in global landmarks
INTEGER                    :: JI,JJ,JK   ! Loop indexes 
                                            ! Sponge zone in complete domain:
INTEGER                    :: IIWEND_ll     ! Sponge zone internal limit for
INTEGER                    :: IIEEND_ll     ! the West, East, South and
INTEGER                    :: IJSEND_ll     ! North boundaries
INTEGER                    :: IJNEND_ll     ! in global domain landmarks
                                            ! Intersection aera between Sponge zone
                                            ! and sub-domain:
INTEGER                    :: IIWENDINT     !  internal limit for
INTEGER                    :: IIEENDINT     ! the West, East, South and
INTEGER                    :: IJSENDINT     ! North boundaries
INTEGER                    :: IJNENDINT     ! in sub_domain landmarks
INTEGER                    :: IIBINT,IJBINT ! Lower bounds of the intersection 
                                            ! in x and y directions
INTEGER                    :: IIEINT,IJEINT ! Upper bounds of the intersection 
                                            ! in x and y directions
!
REAL                       :: ZXDEPTH     ! Inverse squared length of the rim 
REAL                       :: ZYDEPTH     ! zone in the x and y directions
REAL                       :: ZXUPOS      ! Running squared distance to the
REAL                       :: ZYUPOS      ! lateral boundary for "u" locations
REAL                       :: ZXVPOS      ! Running squared distance to the
REAL                       :: ZYVPOS      ! lateral boundary for "v" locations
REAL                       :: ZXWPOS      ! Running squared distance to the
REAL                       :: ZYWPOS      ! lateral boundary for "w" locations
REAL                       :: ZPOS        ! Normalized distance
                                          ! to the inner boundary
LOGICAL, DIMENSION(7)      :: GHORELAXR   ! local array of logical
LOGICAL, DIMENSION(12)     :: GHORELAXSV  ! local array of logical
INTEGER                    :: IINFO_ll    ! return code of parallel routine
INTEGER                    :: IIMAX_ll,IJMAX_ll ! Number of points of
                                                ! Global physical domain
                                                ! in the x and y directions
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
!              ---------------------------------------
! 
CALL GET_INDICE_ll( IIB,IJB,IIE,IJE)
IKE = SIZE(PZZ,3) - JPVEXT
! Global physical dimensions
CALL GET_GLOBALDIMS_ll ( IIMAX_ll,IJMAX_ll)
!
!-------------------------------------------------------------------------------
!
IF(OVE_RELAX) THEN
!
!*       2.    COMPUTE THE VERTICAL INDEX OF THE ABS. LAYER BASE
!              -------------------------------------------------
!
!*       2.1   compute the mean orography
!
!   global sum on global extended domain at k=2
!JUAN16
  ZZSHAT =  SUM_DD_R2_ll(PZZ(IIB:IIE,IJB:IJE,2))
!JUAN16
!   complete with average
  ZZSHAT = ZZSHAT / (IIMAX_ll*IJMAX_ll)
!
!*       2.2   compute the mean  ZHAT height corresponding to PALZBOT
!
  ZZTOP   = PZHAT(IKE+1)

  ZALZHAT = ZZTOP * (PALZBOT  - ZZSHAT) / (ZZTOP -ZZSHAT)
!
!*       2.3   select for KALBOT the closest vertical level to ZALZHAT
!
  IWORK(:) = MINLOC( ABS (PZHAT(:)-ZALZHAT))
  KALBOT = IWORK(1)
! KALBOT = GMINLOC_ll( ABS (PZHAT(:)-ZALZHAT))
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTES THE DAMPING COEF. OF THE ABSORBING LAYER
!              -------------------------------------------------
!
  ZWORK = 0.5 * XPI / (ZZTOP -PZHAT(KALBOT))
  PALKW(:)=0.
  PALK(:)=0.
!
  DO JK = KALBOT, IKE+1
    PALKW(JK) = PALKTOP * SIN ( ZWORK * (PZHAT(JK) - PZHAT(KALBOT))) **2
    PALKW(JK) = PALKW(JK) / (1. + 2. * PTSTEP * PALKW(JK))
  END DO
!
  DO JK = KALBOT, IKE
    ZZHATK = 0.5 * (PZHAT(JK) +PZHAT(JK+1))
    PALK(JK) = PALKTOP * SIN ( ZWORK * (ZZHATK - PZHAT(KALBOT))) **2
    PALK(JK) = PALK(JK) / (1. + 2. * PTSTEP * PALK(JK))
  END DO
!
END IF
!
IF (OVE_RELAX_GRD) THEN
!
  IWORK(:) = MINLOC( ABS (PZHAT(:)-PALZBAS))
  KALBAS = IWORK(1)
!
!
!
  ZWORK = 0.5 * XPI / PZHAT(KALBAS)
  PALKWBAS(:)=0.
  PALKBAS(:)=0.
!
  DO JK = 1,KALBAS+1
    PALKWBAS(JK) = PALKGRD * SIN ( ZWORK * (-PZHAT(JK) + PZHAT(KALBAS))) **2
    PALKWBAS(JK) = PALKWBAS(JK) / (1. + 2. * PTSTEP * PALKWBAS(JK))
  END DO
!
  DO JK = 1,KALBAS
    ZZHATK = 0.5 * (PZHAT(JK) +PZHAT(JK+1))
    PALKBAS(JK) = PALKGRD * SIN ( ZWORK * (-ZZHATK + PZHAT(KALBAS))) **2
    PALKBAS(JK) = PALKBAS(JK) / (1. + 2. * PTSTEP * PALKBAS(JK))
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTES THE HORIZONTAL RELAXATION COEFFICIENTS
!              -----------------------------------------------
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
!
IF ( ANY(GHORELAXR) .OR. ANY(GHORELAXSV) .OR. ANY(OHORELAX_SV) &
                    .OR. OHORELAX_UVWTH  .OR. OHORELAX_TKE     ) THEN 
!
!*       4.1   some settings
!
  IF ( KRIMX /= 0 ) THEN 
    ZXDEPTH  = (1.0/FLOAT(KRIMX))**2
  ELSE
    ZXDEPTH  = 0.
  END IF
  IF ( KRIMY /= 0 ) THEN 
    ZYDEPTH  = (1.0/FLOAT(KRIMY))**2
  ELSE
    ZYDEPTH  = 0.
  END IF
!
  PKURELAX(:,:) = 0.0
  PKVRELAX(:,:) = 0.0
  PKWRELAX(:,:) = 0.0
!
!             Physical global domain bounds
!
  IIB_ll = 1 + JPHEXT
  IIE_ll = IIMAX_ll + JPHEXT
  IJB_ll = 1 + JPHEXT
  IJE_ll = IJMAX_ll + JPHEXT
!
!             Global limits of the relaxation area
!
  IF (L2D) THEN
    IKRIMAX= KRIMX
    IJSEND_ll = IJB_ll      
    IJNEND_ll = IJE_ll + 1 
    IIWEND_ll = IIB_ll     + KRIMX
    IIEEND_ll = IIE_ll + 1 - KRIMX
  ELSE
    IKRIMAX=MAX(KRIMX,KRIMY) 
    IJSEND_ll = IJB_ll     + KRIMY
    IJNEND_ll = IJE_ll + 1 - KRIMY
    IIWEND_ll = IIB_ll     + KRIMX
    IIEEND_ll = IIE_ll + 1 - KRIMX
  END IF
!
!
!  Four areas are processed  :
!     south : IIB_ll    , IIE_ll       ; IJB_ll    , IJSEND_ll-1
!     north : IIB_ll    , IIE_ll       ; IJNEND_ll , IJE_ll
!     west  : IIB_ll    , IIWEND_ll-1  ; IJSEND_ll , IJNEND_ll-1
!     east  : IIEEND_ll , IIE_ll       ; IJSEND_ll , IJNEND_ll-1
!  initializations are made for each area if the intersection with
!  the physical sub-domain isn't empty
   
!   Get coordonates origines of the extended sub-domain
    CALL GET_OR_ll ('B',IORX,IORY)
!
  !
  ! computations only made in the 3D version
  !
!
!*       4.2   scans the South boundary
!
  CALL GET_INTERSECTION_ll ( IIB_ll,IJB_ll,IIE_ll,IJSEND_ll-1,      &
                             IIBINT,IJBINT,IIEINT,IJSENDINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll /= 1 ) THEN
    ! The sub-domain  contains the south relaxation area:
    ! intersection limits : IIBINT,IIEINT along x and IJBINT,IJSENDINT along y
    ! along x: three zones are processed IIB_ll      to IIWEND_ll-1
    !                                    IIWEND_ll   to IIEEND_ll-1
    !                                    IIEEND_ll-1 to IIE_ll
    ! initializes IIWENDINT et IIEENDINT in local indices
    ! to process one , two or three zones along x sub-domain
    IIWENDINT = MIN( IIEINT , MAX ( IIBINT , (IIWEND_ll - IORX +1) ))
    IIEENDINT = MAX( IIBINT , MIN ( (IIEINT +1) , (IIEEND_ll - IORX +1) ))
!
    DO JJ=IJBINT,IJSENDINT
      !  in global landmarks
      ZYUPOS = (FLOAT(JJ+IORY-1 - IJSEND_ll) + 0.5)**2
      ZYVPOS = (FLOAT(JJ+IORY-1 - IJSEND_ll)      )**2
      ZYWPOS = (FLOAT(JJ+IORY-1 - IJSEND_ll) + 0.5)**2
      !
      DO JI=IIBINT,IIWENDINT-1
        ZXUPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)     )**2
        ZXVPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)+ 0.5)**2
        ZXWPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)+ 0.5)**2
        !
        ZPOS = MIN(1.,SQRT(ZXUPOS*ZXDEPTH+ZYUPOS*ZYDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZXVPOS*ZXDEPTH+ZYVPOS*ZYDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZXWPOS*ZXDEPTH+ZYWPOS*ZYDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      END DO
      !
      DO JI=IIWENDINT,IIEENDINT-1
        ZPOS = MIN(1.,SQRT(ZYUPOS*ZYDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZYVPOS*ZYDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZYWPOS*ZYDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      END DO
      !
      DO JI=IIEENDINT,IIEINT
        ZXUPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)     )**2
        ZXVPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)+ 0.5)**2
        ZXWPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)+ 0.5)**2
        !
        ZPOS = MIN(1.,SQRT(ZXUPOS*ZXDEPTH+ZYUPOS*ZYDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZXVPOS*ZXDEPTH+ZYVPOS*ZYDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZXWPOS*ZXDEPTH+ZYWPOS*ZYDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      !    
      END DO
      !
    END DO
!
  END IF  ! end of if intersection isn t empty
!
!*       4.3   scans the North boundary
!
  CALL GET_INTERSECTION_ll ( IIB_ll,IJNEND_ll,IIE_ll,IJE_ll,   &
                             IIBINT,IJNENDINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll /= 1 ) THEN
    ! The sub-domain  contains the north relaxation area:
    ! intersection limits : IIBINT,IIEINT along x and IJNENDINT,IJEINT along y
    ! along x: three zones are processed IIB_ll      to IIWEND_ll-1
    !                                    IIWEND_ll   to IIEEND_ll-1
    !                                    IIEEND_ll-1 to IIE_ll
    ! initializes IIWENDINT et IIWENDINT in local indices
    ! to process one , two or three zones along x sub-domain
    IIWENDINT = MIN( IIEINT , MAX ( IIBINT , (IIWEND_ll - IORX +1) ))
    IIEENDINT = MAX( IIBINT , MIN ( (IIEINT +1) , (IIEEND_ll - IORX +1) ))
!
    DO JJ=IJNENDINT,IJEINT
      ! in global landmarks
      ZYUPOS = (FLOAT(JJ+IORY-1 - IJNEND_ll) + 0.5)**2
      ZYVPOS = (FLOAT(JJ+IORY-1 - IJNEND_ll)      )**2
      ZYWPOS = (FLOAT(JJ+IORY-1 - IJNEND_ll) + 0.5)**2
      !
      DO JI=IIBINT,IIWENDINT-1
        ZXUPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)     )**2
        ZXVPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)+ 0.5)**2
        ZXWPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)+ 0.5)**2
        !
        ZPOS = MIN(1.,SQRT(ZXUPOS*ZXDEPTH+ZYUPOS*ZYDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
        !
        ZPOS = MIN(1.,SQRT(ZXVPOS*ZXDEPTH+ZYVPOS*ZYDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
       !
        ZPOS = MIN(1.,SQRT(ZXWPOS*ZXDEPTH+ZYWPOS*ZYDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      END DO
      !
      DO JI=IIWENDINT,IIEENDINT-1
       !
        ZPOS = MIN(1.,SQRT(ZYUPOS*ZYDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
       !
        ZPOS = MIN(1.,SQRT(ZYVPOS*ZYDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
       !
        ZPOS = MIN(1.,SQRT(ZYWPOS*ZYDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      END DO
      !
      DO JI=IIEENDINT,IIEINT
        ZXUPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)     )**2
        ZXVPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)+ 0.5)**2
        ZXWPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)+ 0.5)**2
        !
        ZPOS = MIN(1.,SQRT(ZXUPOS*ZXDEPTH+ZYUPOS*ZYDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      !
        ZPOS = MIN(1.,SQRT(ZXVPOS*ZXDEPTH+ZYVPOS*ZYDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
       !
        ZPOS = MIN(1.,SQRT(ZXWPOS*ZXDEPTH+ZYWPOS*ZYDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      END DO
      !
    END DO
    !
  END IF  ! end of if intersection isn t empty
!
!*       4.4   scans the remaining band in the West boundary
!
  CALL GET_INTERSECTION_ll ( IIB_ll,IJSEND_ll,IIWEND_ll-1,IJNEND_ll-1 ,&
                             IIBINT,IJSENDINT,IIWENDINT,IJNENDINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll /= 1 ) THEN
    ! The sub-domain  contains the west relaxation area:
    ! intersection limits : IIBINT,IIWENDINT along x and IJSENDINT,IJNENDINT along y
!   
    DO JI=IIBINT,IIWENDINT
      ZXUPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)     )**2
      ZXVPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)+ 0.5)**2
      ZXWPOS = (FLOAT(JI+IORX-1 - IIWEND_ll)+ 0.5)**2
      DO JJ=IJSENDINT,IJNENDINT
       !
        ZPOS = MIN(1.,SQRT(ZXUPOS*ZXDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      !
        ZPOS = MIN(1.,SQRT(ZXVPOS*ZXDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
       !
        ZPOS = MIN(1.,SQRT(ZXWPOS*ZXDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
     END DO
    END DO
  END IF  ! end of if intersection isn t empty
!
!*       4.5   scans the remaining band in the East boundary
!
  CALL GET_INTERSECTION_ll ( IIEEND_ll,IJSEND_ll,IIE_ll,IJNEND_ll-1, &
                             IIEENDINT,IJSENDINT,IIEINT,IJNENDINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll /= 1 ) THEN
    ! The sub-domain  contains the east relaxation area:
    ! intersection limits : IIENDINT,IIEINT along x and IJSENDINT,IJNENDINT along y
!
    DO JI=IIEENDINT,IIEINT
      ZXUPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)     )**2
      ZXVPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)+ 0.5)**2
      ZXWPOS = (FLOAT(JI+IORX-1 - IIEEND_ll)+ 0.5)**2
      DO JJ=IJSENDINT,IJNENDINT
       !
        ZPOS = MIN(1.,SQRT(ZXUPOS*ZXDEPTH))
        PKURELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      !
        ZPOS = MIN(1.,SQRT(ZXVPOS*ZXDEPTH))
        PKVRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
       !
        ZPOS = MIN(1.,SQRT(ZXWPOS*ZXDEPTH))
        PKWRELAX(JI,JJ) = PRIMKMAX*RELAX(ZPOS,IKRIMAX)
      END DO
    END DO
    !
  END IF  ! end of if intersection isn t empty
  !
END IF
!
!
!*       4.6   prepares the implicit scheme
!
PKURELAX(:,:) = PKURELAX(:,:) / (1.0 + 2.0 * PTSTEP * PKURELAX(:,:))
PKVRELAX(:,:) = PKVRELAX(:,:) / (1.0 + 2.0 * PTSTEP * PKVRELAX(:,:)) 
PKWRELAX(:,:) = PKWRELAX(:,:) / (1.0 + 2.0 * PTSTEP * PKWRELAX(:,:))
!
!*       4.7   sets up the logical mask
!
OMASK_RELAX(:,:) = (PKURELAX(:,:) > 1.E-16) .OR. &
                   (PKVRELAX(:,:) > 1.E-16) .OR. &
                   (PKWRELAX(:,:) > 1.E-16)
!   
!-------------------------------------------------------------------------------
!
END SUBROUTINE RELAXDEF
