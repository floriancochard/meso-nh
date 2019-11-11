!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_ELEC_TRIDZ
!     ######################
!
INTERFACE
!
      SUBROUTINE ELEC_TRIDZ(HLBCX,HLBCY,                                &
                      PMAP,PDXHAT,PDYHAT,PDXHATM,PDYHATM,PRHOM,         &
                      PAF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,            &
                      PRHODJ,PTHVREF,PZZ,PBFY,PEPOTFW_TOP,              &
                      PBFB,PBF_SXP2_YP1_Z)!++cb - Z_SPLITTING
!
IMPLICIT NONE
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! density of reference * J
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF    ! Virtual Potential
                                        ! Temperature of the reference state
!
REAL, DIMENSION(:,:), INTENT(IN) :: PMAP         ! scale factor
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! height z
!
REAL, DIMENSION(:), INTENT(IN) :: PDXHAT          ! Stretching in x direction
REAL, DIMENSION(:), INTENT(IN) :: PDYHAT          ! Stretching in y direction
!
REAL, INTENT(OUT) :: PDXHATM                     ! mean grid increment in the x
                                                 ! direction
REAL, INTENT(OUT) :: PDYHATM                     ! mean grid increment in the y
                                                 ! direction
!
REAL, DIMENSION (:), INTENT(OUT) :: PRHOM        !  mean of XRHODJ on the plane
                                                 !  x y localized at a mass 
                                                 !  level
!
REAL, DIMENSION(:), INTENT(OUT)     :: PAF,PCF  ! vectors giving the nonvanishing
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFY     ! elements (yslice) of the tri-diag.
                                                ! matrix in the pressure eq. 
!++cb - Z_SPLITTING
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFB     ! elements (bsplit slide) of the tri-diag.
                                                ! matrix in the pressure eq. 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBF_SXP2_YP1_Z ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq
!++cb
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PEPOTFW_TOP ! top boundary condition of 
                                                   ! the tri-diag. system; it 
                                                   ! has the dimension of an
                                                   ! electrical potential
!
                                                 ! arrays of sin or cos values
                                                 ! for the FFT :
REAL, DIMENSION(:), INTENT(OUT) :: PTRIGSX       ! - along x
REAL, DIMENSION(:), INTENT(OUT) :: PTRIGSY       ! - along y
!
                                                 ! decomposition in prime
                                                 ! numbers for the FFT:
INTEGER, DIMENSION(19), INTENT(OUT) :: KIFAXX      ! - along x
INTEGER, DIMENSION(19), INTENT(OUT) :: KIFAXY      ! - along y

!
END SUBROUTINE ELEC_TRIDZ
!
END INTERFACE
!
END MODULE MODI_ELEC_TRIDZ
!
!     ###################################################################
      SUBROUTINE ELEC_TRIDZ(HLBCX,HLBCY,                                &
                      PMAP,PDXHAT,PDYHAT,PDXHATM,PDYHATM,PRHOM,         &
                      PAF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,            &
                      PRHODJ,PTHVREF,PZZ,PBFY,PEPOTFW_TOP,              &
                      PBFB,PBF_SXP2_YP1_Z)
!    ####################################################################
!
!!****  *ELEC_TRIDZ * - Compute coefficients for the flat operator to get the
!!                     electric potential from the Gauss equation
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the vertical time independent 
!     coefficients a(k), b(k), c(k) required for the calculation of the "flat"  
!     (i.e. neglecting the orography) operator Laplacian.  RHOJ is averaged on 
!     the whole horizontal domain. The form of the eigenvalues  of the flat 
!     operator depends on the lateral boundary conditions. Furthermore, this
!     routine initializes TRIGS and IFAX arrays required for the FFT transform 
!     used in the routine PRECOND.
!     ELEC_TRIDZ (to invert the Gauss equation) differs from TRID (to solve the
!     pressure equation) by the bottom boundary condition, here Dirichlet 
!     instead of Neumann, because the earth surface is conductive so it is a
!     surface with an electrical equipotential referenced to zero.
!
!!**  METHOD
!!    ------
!!     The forms of the eigenvalues of the horizontal Laplacian are given by:
!!     Cyclic conditions:
!!     -----------------
!!                    <rhoj>        2 ( pi     )       <rhoj>      2 (  pi    )
!!     b(m,n) = -4 -----------   sin  (----- m ) -4 ----------- sin  (----- n )
!!                 <dxx> <dxx>        ( imax   )    <dyy> <dyy>      ( jmax   )
!!
!!     Open conditions:
!!     -----------------
!!                    <rhoj>        2 ( pi     )       <rhoj>      2 (  pi    )
!!     b(m,n) = -4 -----------   sin  (----- m ) -4 ----------- sin  (----- n )
!!                 <dxx> <dxx>        ( 2imax  )    <dyy> <dyy>      ( 2jmax  )
!!      
!!     Cyclic condition along x and open condition along y:
!!     ------------------------------------------------------
!!                    <rhoj>        2 ( pi     )       <rhoj>      2 (  pi    )
!!     b(m,n) = -4 -----------   sin  (----- m ) -4 ----------- sin  (----- n )
!!                 <dxx> <dxx>        ( imax   )    <dyy> <dyy>      ( 2jmax  )
!!
!!     Open condition along x and cyclic condition along y:
!!     ------------------------------------------------------
!!                    <rhoj>        2 ( pi     )       <rhoj>      2 (  pi    )
!!     b(m,n) = -4 -----------   sin  (----- m ) -4 ----------- sin  (----- n )
!!                 <dxx> <dxx>        ( 2imax  )    <dyy> <dyy>      ( jmax   )
!!
!!     where m = 0,1,2....imax-1, n = 0,1,2....jmax-1
!!     Note that rhoj contains the Jacobian J = Deltax*Deltay*Deltaz = volume of
!!     an elementary mesh.

!!
!!    EXTERNAL
!!    --------
!!      Function FFTFAX: initialization of TRIGSX,IFAXX,TRIGSY,IFAXY for  
!!      the FFT transform
!!      GET_DIM_EXT_ll    : get extended sub-domain sizes
!!      GET_INDICE_ll     : get physical sub-domain bounds 
!!      GET_DIM_PHYS_ll   : get physical sub-domain sizes
!!      GET_GLOBALDIMS_ll : get physical global domain sizes
!!      GET_OR_ll         : get origine coordonates of the physical sub-domain 
!!                          in global indices
!!      REDUCESUM_ll      : sum into a scalar variable
!!      GET_SLICE_ll      : get a slice of the global domain
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CST : define constants
!!         XPI : pi
!!         XCPD
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT, JPVEXT: define the number of marginal points out of the 
!!        physical domain along horizontal and vertical directions respectively
!!      Module MODD_CONF: model configurations 
!!         LCARTESIAN: logical for CARTESIAN geometry
!!                     .TRUE. = Cartesian geometry used
!!         L2D: logical for 2D model version
!!  
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine ELEC_TRIDZ)
!!      
!!    AUTHOR
!!    ------
!!	P. HÅreil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/07/94 
!!                  14/04/95  (J. Stein) bug in the ZDZM computation
!!                                        ( stretched case)  
!!                   8/07/96  (P. Jabouille) change the FFT initialization
!!                            which now works for odd number.
!!                  14/01/97 Durran anelastic equation (Stein,Lafore)
!!                  15/06/98  (D.Lugato, R.Guivarch) Parallelisation
!!                  10/08/98 (N. Asencio) add parallel code
!!                                        use PDXHAT, PDYHAT and not PXHAT,PYHAT
!!                                        PBFY is initialized
!!                  20/08/00 (J. Stein, J. Escobar) optimisation of the solver
!!                                        PBFY transposition
!!                  14/03/02 (P. Jabouille) set values for meaningless spectral coefficients
!!                                       (to avoid problem in bouissinesq configuration)
!!                  01/07/12 (J-P. Pinty)  add a non-homogeneous fair-weather 
!!                                         top boundary condition (Neuman)
!!                  16/10/14 (C. Bovalo)   add the z-splitting
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CST
USE MODD_CONF
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAMETERS
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_MSG
!++cb - Z_SPLITTING
USE MODE_SPLITTINGZ_ll , ONLY : GET_DIM_EXTZ_ll,GET_ORZ_ll,LWESTZ_ll,LSOUTHZ_ll
!
USE MODE_REPRO_SUM
!++cb -
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! density of reference * J
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF    ! Virtual Potential
                                        ! Temperature of the reference state
!
REAL, DIMENSION(:,:), INTENT(IN) :: PMAP         ! scale factor
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! height z
!
REAL, DIMENSION(:), INTENT(IN) :: PDXHAT          ! Stretching in x direction
REAL, DIMENSION(:), INTENT(IN) :: PDYHAT          ! Stretching in y direction
!
REAL, INTENT(OUT) :: PDXHATM                     ! mean grid increment in the x
                                                 ! direction
REAL, INTENT(OUT) :: PDYHATM                     ! mean grid increment in the y
                                                 ! direction
!
REAL, DIMENSION (:), INTENT(OUT) :: PRHOM        !  mean of XRHODJ on the plane
                                                 !  x y localized at a mass 
                                                 !  level
!
REAL, DIMENSION(:), INTENT(OUT)     :: PAF,PCF  ! vectors giving the nonvanishing
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFY     ! elements (yslice) of the tri-diag.
! matrix in the pressure eq. which is transposed. PBFY is a y-slices structure
!++cb - Z_SPLITTING
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBFB     ! elements (bsplit slide) of the tri-diag.
                                                ! matrix in the pressure eq. 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PBF_SXP2_YP1_Z ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
!++cb
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PEPOTFW_TOP ! top boundary condition of 
                                                   ! the tri-diag. system; it 
                                                   ! has the dimension of an
                                                   ! electrical potential
!                                 
                                                 ! arrays of sin or cos values
                                                 ! for the FFT :
REAL, DIMENSION(:), INTENT(OUT) :: PTRIGSX       ! - along x
REAL, DIMENSION(:), INTENT(OUT) :: PTRIGSY       ! - along y
!
                                                 ! decomposition in prime
                                                 ! numbers for the FFT:
INTEGER, DIMENSION(19), INTENT(OUT) :: KIFAXX      ! - along x
INTEGER, DIMENSION(19), INTENT(OUT) :: KIFAXY      ! - along y

!
!*       0.2   declarations of local variables
!
INTEGER                         ::  IRESP    ! FM return code
INTEGER                         :: ILUOUT    ! Logical unit number for
                                             ! output-listing   
INTEGER :: IIB,IIE,IJB,IJE,IKB,IKE     ! indice values of the physical subdomain
INTEGER :: IKU                         ! size of the arrays along z
INTEGER :: IIB_ll,IIE_ll,IJB_ll,IJE_ll !  indice values of the physical global domain
INTEGER :: IIMAX,IJMAX       ! Number of points of the physical subdomain
INTEGER :: IIMAX_ll,IJMAX_ll ! Number of points of Global physical domain
!
INTEGER :: JI,JJ,JK                                    ! loop indexes
!
INTEGER :: INN        ! temporary result for the computation of array TRIGS
!
REAL, DIMENSION (:,:), ALLOCATABLE  :: ZEIGEN_ll  ! eigenvalues b(m,n) in global representation
REAL, DIMENSION (:),   ALLOCATABLE  :: ZEIGENX_ll ! used for the computation of ZEIGEN_ll 
!
REAL, DIMENSION( SIZE(PDXHAT))  :: ZWORKX  ! work array to compute PDXHATM
REAL, DIMENSION( SIZE(PDYHAT))  :: ZWORKY  ! work array to compute PDYHATM
!
REAL :: ZGWNX,ZGWNY           ! greater wave numbers allowed by the model 
                              ! configuration in x and y directions respectively
!
REAL, DIMENSION (SIZE(PZZ,3)) :: ZDZM      !  mean of deltaz on the plane x y
                                           !  localized at a w-level
!
REAL :: ZANGLE,ZDEL ! needed for the initialization of the arrays used by the FFT
!
REAL :: ZINVMEAN ! inverse of inner points number in an horizontal grid
!
INTEGER ::  IINFO_ll         ! return code of parallel routine
REAL, DIMENSION (SIZE(PMAP,1)) :: ZXMAP ! extraction of PMAP array along x
REAL, DIMENSION (SIZE(PMAP,2)) :: ZYMAP ! extraction of PMAP array along y
INTEGER :: IORXY_ll,IORYY_ll     ! origin's coordinates of the y-slices subdomain
INTEGER :: IIUY_ll,IJUY_ll       ! dimensions of the y-slices subdomain
INTEGER :: IXMODE_ll,IYMODE_ll   ! number of modes in the x and y direction for global point of view
INTEGER :: IXMODEY_ll,IYMODEY_ll ! number of modes in the x and y direction for y_slice point of view
!++cb - Z_SPLITTING
INTEGER :: IORXB_ll,IORYB_ll     ! origin's coordinates of the b-slices subdomain
INTEGER :: IIUB_ll,IJUB_ll       ! dimensions of the b-slices subdomain
INTEGER :: IXMODEB_ll,IYMODEB_ll ! number of modes in the x and y direction for b_slice point of view
!
INTEGER :: IORX_SXP2_YP1_Z_ll,IORY_SXP2_YP1_Z_ll     ! origin's coordinates of the b-slices subdomain
INTEGER :: IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll       ! dimensions of the b-slices subdomain
INTEGER :: IXMODE_SXP2_YP1_Z_ll,IYMODE_SXP2_YP1_Z_ll ! number of modes in the x and y direction for b_slice point of view
!++cb
!
!++JUAN-16
!TYPE(DOUBLE_DOUBLE)  , DIMENSION (SIZE(PZZ,3)) :: ZRHOM_ll   , ZDZM_ll
REAL, ALLOCATABLE, DIMENSION(:,:)              :: ZRHOM_2D   , ZDZM_2D
!++JUAN-16
!
!
!
!
!
!------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION
!              --------------
!
!*       1.1  retrieve a logical unit number
!             ------------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       1.2  compute loop bounds
!             -------------------
!
!  extended sub-domain
CALL GET_DIM_EXT_ll ('Y',IIUY_ll,IJUY_ll)
!++cb - Z_SPLITTING
CALL GET_DIM_EXT_ll ('B',IIUB_ll,IJUB_ll)
! P1/P2 splitting
CALL GET_DIM_EXTZ_ll('SXP2_YP1_Z',IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll)
!++cb
IKU=SIZE(PRHODJ,3)                        
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1   +JPVEXT                          
IKE=IKU -JPVEXT
!  physical sub-domain
CALL GET_DIM_PHYS_ll ( 'B',IIMAX,IJMAX)
!
!  global physical domain limits
CALL GET_GLOBALDIMS_ll ( IIMAX_ll, IJMAX_ll)
IIB_ll = 1 + JPHEXT
IIE_ll = IIMAX_ll  + JPHEXT
IJB_ll = 1 + JPHEXT
IJE_ll = IJMAX_ll + JPHEXT
!
!     the use of local array ZEIGENX and ZEIGEN would require some technical modifications
!
ALLOCATE (ZEIGENX_ll(IIMAX_ll + 2*JPHEXT))
ALLOCATE (ZEIGEN_ll(IIMAX_ll  + 2*JPHEXT, IJMAX_ll + 2*JPHEXT))
ZEIGEN_ll = 0.0
!  Get the origin coordinates of the extended sub-domain in global landmarks
CALL GET_OR_ll('Y',IORXY_ll,IORYY_ll)
!++cb - Z_SPLITTING
CALL GET_OR_ll('B',IORXB_ll,IORYB_ll)
! P1/P2 Splitting
CALL GET_ORZ_ll('SXP2_YP1_Z',IORX_SXP2_YP1_Z_ll,IORY_SXP2_YP1_Z_ll)
!++cb
!
!*       1.3 allocate x-slice array

!
!*       1.4 variables for the eigenvalues computation
!
ZGWNX = XPI/REAL(IIMAX_ll)
ZGWNY = XPI/REAL(IJMAX_ll)
!
!------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE AVERAGE OF RHOJ*CPD*THETAVREF ALONG XY
!              --------------------------------------------------
!
ZINVMEAN = 1./REAL(IIMAX_ll*IJMAX_ll)
!JUAN-16
ALLOCATE(ZRHOM_2D(IIB:IIE, IJB:IJE))
!
DO JK = 1,SIZE(PZZ,3)
   IF ( CEQNSYS == 'DUR' .OR. CEQNSYS == 'MAE' ) THEN
      DO JJ = IJB,IJE
         DO JI = IIB,IIE
            ZRHOM_2D(JI,JJ) = PRHODJ(JI,JJ,JK)*XCPD*PTHVREF(JI,JJ,JK)*ZINVMEAN       
         END DO
      END DO
   ELSEIF ( CEQNSYS == 'LHE' ) THEN
      DO JJ = IJB,IJE
         DO JI = IIB,IIE
            ZRHOM_2D(JI,JJ) = PRHODJ(JI,JJ,JK)*ZINVMEAN          
         END DO
      END DO
   END IF
   !  global sum
   PRHOM(JK) =  SUM_DD_R2_ll(ZRHOM_2D)
END DO

!
!  global sum
!CALL REDUCESUM_ll(ZRHOM_ll,IINFO_ll)
!PRHOM = ZRHOM_ll
!JUAN-16
!
!------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE MEAN INCREMENT BETWEEN Z LEVELS
!              -------------------------------------------
!
!JUAN-16
!ZDZM_ll = 0.
ALLOCATE(ZDZM_2D(IIB:IIE, IJB:IJE))
!
DO JK = IKB-1,IKE
   DO JJ = IJB,IJE
      DO JI = IIB,IIE
         ZDZM_2D(JI,JJ) = (PZZ(JI,JJ,JK+1)-PZZ(JI,JJ,JK))*ZINVMEAN
      END DO
   END DO
   ZDZM(JK)  =  SUM_DD_R2_ll(ZDZM_2D)
END DO
ZDZM(IKE+1) = ZDZM(IKE)
!
!  global sum
!CALL REDUCESUM_ll(ZDZM_ll,IINFO_ll)
!ZDZM = ZDZM_ll
!JUAN-16
!
!
! vertical average to arrive at a w-level 
DO JK = IKE+1,IKB,-1
  ZDZM(JK) = (ZDZM(JK) + ZDZM(JK-1))*0.5
END DO
!
ZDZM(IKB-1) = -999.
!
!------------------------------------------------------------------------------
!
!*       4.    COMPUTE THE MEAN INCREMENT BETWEEN X LEVELS
!              -------------------------------------------
!
PDXHATM =0.
!     . local sum
IF (LCARTESIAN) THEN
    PDXHATM = SUM_1DFIELD_ll ( PDXHAT,'X',IIB_ll,IIE_ll,IINFO_ll)
ELSE
  ! Extraction of x-slice PMAP at j=(IJB_ll+IJE_ll)/2
  CALL GET_SLICE_ll (PMAP,'X',(IJB_ll+IJE_ll)/2,ZXMAP(IIB:IIE) &
                       ,IIB,IIE,IINFO_ll)
  ! initialize the work array = PDXHAT/ZXMAP
  ZWORKX(IIB:IIE) = PDXHAT(IIB:IIE)/ ZXMAP (IIB:IIE)
  PDXHATM = SUM_1DFIELD_ll ( ZWORKX,'X',IIB_ll,IIE_ll,IINFO_ll)
END IF
!    . division to complete sum
PDXHATM = PDXHATM /  REAL(IIMAX_ll)
!
!------------------------------------------------------------------------------
!
!*       5.    COMPUTE THE MEAN INCREMENT BETWEEN Y LEVELS
!              -------------------------------------------
!
PDYHATM = 0.
IF (LCARTESIAN) THEN
  PDYHATM = SUM_1DFIELD_ll ( PDYHAT,'Y',IJB_ll,IJE_ll,IINFO_ll)
ELSE
  ! Extraction of y-slice PMAP at i=IIB_ll+IIE_ll/2 
  CALL GET_SLICE_ll (PMAP,'Y',(IIB_ll+IIE_ll)/2,ZYMAP(IJB:IJE) &
                       ,IJB,IJE,IINFO_ll)
  ! initialize the work array = PDYHAT / ZYMAP
  ZWORKY(IJB:IJE) = PDYHAT(IJB:IJE) / ZYMAP (IJB:IJE)
  PDYHATM = SUM_1DFIELD_ll ( ZWORKY,'Y',IJB_ll,IJE_ll,IINFO_ll)
END IF
!    . division to complete sum
PDYHATM= PDYHATM / REAL(IJMAX_ll)
!
!------------------------------------------------------------------------------
!
!*      6.    COMPUTE THE OUT-DIAGONAL ELEMENTS A AND C OF THE MATRIX
!             -------------------------------------------------------
!
DO JK = IKB,IKE
  PAF(JK) = 0.5 * ( PRHOM(JK-1) + PRHOM(JK)   ) / ZDZM(JK)   **2 
  PCF(JK) = 0.5 * ( PRHOM(JK)   + PRHOM(JK+1) ) / ZDZM(JK+1) **2 
END DO
!
! at the upper level PAF and PCF are computed using the Neumann 
! condition applying on the vertical component of the gradient
! while at the lower level, the Dirichlet condition is used
!            
!            
! Neumann boundary condition (top of atmosphere)
!            
PAF(IKE+1) = -0.5 * ( PRHOM(IKE)   + PRHOM(IKE+1) ) / ZDZM(IKE+1) **2 
!            
! Dirichlet boundary condition (earth surface)
!            
PCF(IKB-1) = 0.0
!
PAF(IKB-1) = 999.
PCF(IKE+1) = 999.
!
!------------------------------------------------------------------------------
!*       7.    COMPUTE THE DIAGONAL ELEMENTS B OF THE MATRIX
!              ---------------------------------------------
!
!*       7.1   compute the horizontal eigenvalues 
!
!
!*       7.1.1   compute the eigenvalues along the x direction 
!
SELECT CASE (HLBCX(1))
! in the cyclic case, the eigenvalues are the same for two following JM values:
! it corresponds to the real and complex parts of the FFT
  CASE('CYCL')                     ! cyclic case
    IXMODE_ll = IIMAX_ll+2
    IXMODEY_ll = IIUY_ll
    IXMODEB_ll = IIUB_ll !++cb - Z_SPLITTING
!
    DO JI = 1,IXMODE_ll
      ZEIGENX_ll(JI) = - (   2. * SIN ( (JI-1)/2*ZGWNX ) / PDXHATM )**2
    END DO
  CASE DEFAULT                !    other cases
    IXMODE_ll = IIMAX_ll
!
!
    IF (LEAST_ll(HSPLITTING='Y')) THEN
      IXMODEY_ll = IIUY_ll - 2
    ELSE
      IXMODEY_ll = IIUY_ll
    END IF
!++cb - Z_SPLITTING
    IF (LEAST_ll(HSPLITTING='B')) THEN
      IXMODEB_ll = IIUB_ll - 2
    ELSE
      IXMODEB_ll = IIUB_ll
    END IF
!++cb
!
!
    DO JI = 1,IXMODE_ll
      ZEIGENX_ll(JI) = - (   2. *SIN (0.5*REAL(JI-1)*ZGWNX ) / PDXHATM  )**2
    END DO
END SELECT
!
!*       7.1.2   compute the eventual eigenvalues along the y direction
!
IF (.NOT. L2D) THEN
!
! y lateral boundary conditions for three-dimensional cases
!
  SELECT CASE (HLBCY(1))
! in the cyclic case, the eigenvalues are the same for two following JN values:
! it corresponds to the real and complex parts of the FFT result
!
    CASE('CYCL')                      ! 3D cyclic case
      IYMODE_ll = IJMAX_ll+2
      IYMODEY_ll = IJUY_ll
      IYMODEB_ll = IJUB_ll !++cb - Z_SPLITTING
!
      DO JJ = 1,IYMODE_ll
        DO JI = 1,IXMODE_ll
          ZEIGEN_ll(JI,JJ) = ZEIGENX_ll(JI) -    &
                         (  2.* SIN ( (JJ-1)/2*ZGWNY ) / PDYHATM  )**2
        END DO
      END DO
!
    CASE DEFAULT                      ! 3D non-cyclic cases
      IYMODE_ll = IJMAX_ll
      IYMODEY_ll = IJUY_ll - 2
      IYMODEB_ll = IJUB_ll - 2 !++cb - Z_SPLITTING
!
      DO JJ = 1,IYMODE_ll
        DO JI = 1,IXMODE_ll
          ZEIGEN_ll(JI,JJ) = ZEIGENX_ll(JI) - (  2.* SIN (0.5*REAL(JJ-1)*ZGWNY ) / &
                                                 PDYHATM   )**2
        END DO
      END DO
!
  END SELECT
ELSE
!
!  copy the x eigenvalue array in a 2D array for a 2D case
!
  IYMODE_ll = 1
  IYMODEY_ll = 1
  ZEIGEN_ll(1:IXMODE_ll,1)=ZEIGENX_ll(1:IXMODE_ll)
!
END IF
!
DEALLOCATE(ZEIGENX_ll)
!
!*       7.2   compute the matrix diagonal elements
!
!
PBFY = 1.
PBFB = 1.  ! ++cb - Z_SLIDE
PBF_SXP2_YP1_Z = 1.  ! ++cb - Z_SLIDE
!++cb
IF (L2D) THEN
  DO JK= IKB,IKE
    DO JJ= 1, IYMODEY_ll
      DO JI= 1, IXMODEY_ll
        PBFY(JI,JJ,JK) = PRHOM(JK)* ZEIGEN_ll(JI+IORXY_ll-1,JJ+IORYY_ll-1) - 0.5 * &
        ( ( PRHOM(JK-1) + PRHOM(JK)   ) / ZDZM(JK)  **2                          &
         +( PRHOM(JK)   + PRHOM(JK+1) ) / ZDZM(JK+1)**2 )
      END DO
    END DO
  END DO
! at the upper level PBFY is computed using the Neumann condition
! while at the lower level, the Dirichlet condition is used
!
! lower level: Dirichlet condition
  PBFY(1:IXMODEY_ll,1:IYMODEY_ll,IKB) = PBFY(1:IXMODEY_ll,1:IYMODEY_ll,IKB) - &
                                        PAF(IKB)
  PAF(IKB)   = 0.0
  PBFY(1:IXMODEY_ll,1:IYMODEY_ll,IKB-1) = 1.0
  !
! upper level: Neumann condition
  PBFY(1:IXMODEY_ll,1:IYMODEY_ll,IKE+1) = 0.5 * ( PRHOM(IKE)   + PRHOM(IKE+1) ) /  &
                                          ZDZM(IKE+1) **2
  !
ELSE
  DO JK= IKB,IKE
    DO JJ= 1, IYMODEY_ll
      DO JI= 1, IXMODEY_ll
        PBFY(JJ,JI,JK) = PRHOM(JK)* ZEIGEN_ll(JI+IORXY_ll-1,JJ+IORYY_ll-1) - 0.5 * &
        ( ( PRHOM(JK-1) + PRHOM(JK)   ) / ZDZM(JK)  **2                          &
         +( PRHOM(JK)   + PRHOM(JK+1) ) / ZDZM(JK+1)**2 )
      END DO
    END DO
  END DO
! at the upper level PBFY is computed using the Neumann condition
! while at the lower level, the Dirichlet condition is used
!
! lower level: Dirichlet
  PBFY(1:IYMODEY_ll,1:IXMODEY_ll,IKB) = PBFY(1:IYMODEY_ll,1:IXMODEY_ll,IKB) - &
                                        PAF(IKB)
!  PAF(IKB)   = 0.0  ! done later
  PBFY(1:IYMODEY_ll,1:IXMODEY_ll,IKB-1) = 1.0
!
! upper level: Neumann
  PBFY(1:IYMODEY_ll,1:IXMODEY_ll,IKE+1) = 0.5 * ( PRHOM(IKE)   + PRHOM(IKE+1) ) /  &
                                          ZDZM(IKE+1) **2
!
!++cb - Z_SPLITTING
!++cb - for Z splitting we need to do the vertical substitution
!++cb - in Bsplitting slides so need PBFB 
  DO JK= IKB,IKE
    DO JJ= IJB,IJE
      DO JI= IIB, IIE
        PBFB(JI,JJ,JK) = PRHOM(JK)* ZEIGEN_ll(JI+IORXB_ll-IIB_ll,JJ+IORYB_ll-IJB_ll) - 0.5 * &
        ( ( PRHOM(JK-1) + PRHOM(JK)   ) / ZDZM(JK)  **2                          &
         +( PRHOM(JK)   + PRHOM(JK+1) ) / ZDZM(JK+1)**2 )
      END DO
    END DO
  END DO
! at the upper level PBFB is computed using the Neumann condition
! while at the lower level, the Dirichlet is used
!
! lower level: Dirichlet
  PBFB(IIB:IIE,IJB:IJE,IKB) = PBFB(IIB:IIE,IJB:IJE,IKB) - &
                              PAF(IKB)
!
!  PAF(IKB)   = 0.0
  PBFB(IIB:IIE,IJB:IJE,IKB-1) = 1.0
  !
! upper level: Neumann
  PBFB(IIB:IIE,IJB:IJE,IKE+1) = 0.5 * ( PRHOM(IKE)   + PRHOM(IKE+1) ) /  &
                                ZDZM(IKE+1) **2
!
IF (HLBCX(1) == 'CYCL' .AND. .NOT.(L2D) ) THEN
   !++cb
   ! fill unused 2 coef with NI+1 coef (lost in Z transposition elsewhere )
   JI = IXMODE_ll -1
   ZEIGEN_ll(2,:)  = ZEIGEN_ll(JI,:)
END IF
IF (HLBCY(1) == 'CYCL' .AND. .NOT.(L2D) ) THEN
   !++cb
   ! fill unused (:,2,:) coef with NJ+1 coef (lost in Z transposition elsewhere )
   JJ = IYMODE_ll -1
   ZEIGEN_ll(:,2)  = ZEIGEN_ll(:,JJ)
END IF
  !
!++cb - Z_SPLITTING
!++cb - Z_SPLITTING
!++cb - for Z splitting we need to do the vertical substitution
!++cb - in _SXP2_YP1_Zsplitting slides so need PBF_SXP2_YP1_Z
  DO JK=IKB,IKE
    DO JJ= 1, IJU_SXP2_YP1_Z_ll
      DO JI= 1, IIU_SXP2_YP1_Z_ll
        PBF_SXP2_YP1_Z(JI,JJ,JK) = PRHOM(JK)* ZEIGEN_ll(JI+IORX_SXP2_YP1_Z_ll-IIB_ll,JJ+IORY_SXP2_YP1_Z_ll-IJB_ll) - 0.5 * &
        ( ( PRHOM(JK-1) + PRHOM(JK)   ) / ZDZM(JK)  **2                          &
         +( PRHOM(JK)   + PRHOM(JK+1) ) / ZDZM(JK+1)**2 )
      END DO
    END DO
  END DO
! at the upper level PBF_SXP2_YP1_Z is computed using the Neumann condition
! while at the lower level, the Dirichlet condition is used
!
! lower level: Dirichlet
  PBF_SXP2_YP1_Z(1:IIU_SXP2_YP1_Z_ll,1:IJU_SXP2_YP1_Z_ll,IKB) = PBF_SXP2_YP1_Z(1:IIU_SXP2_YP1_Z_ll,1:IJU_SXP2_YP1_Z_ll,IKB) - &
                                                                PAF(IKB)
  PAF(IKB) = 0.0
  PBF_SXP2_YP1_Z(1:IIU_SXP2_YP1_Z_ll,1:IJU_SXP2_YP1_Z_ll,IKB-1) = 1.0
  !
! upper level: Neumann
  PBF_SXP2_YP1_Z(1:IIU_SXP2_YP1_Z_ll,1:IJU_SXP2_YP1_Z_ll,IKE+1) =  0.5 * ( PRHOM(IKE)   + PRHOM(IKE+1) ) /  &
                               ZDZM(IKE+1) **2
!++cb
END IF
!
! second coefficent is meaningless in cyclic case
IF (HLBCX(1) == 'CYCL' .AND. L2D .AND. SIZE(PBFY,1) .GE. 2 ) PBFY(2,:,:)=1.
IF (HLBCX(1) == 'CYCL' .AND. .NOT.(L2D) .AND. LWEST_ll(HSPLITTING='Y') .AND. SIZE(PBFY,2) .GE. 2 ) &
    PBFY(:,2,:)=1.
IF (HLBCY(1) == 'CYCL' .AND. .NOT.(L2D) .AND. SIZE(PBFY,1) .GE. 2  ) PBFY(2,:,:)=1.
!++cb - Z_SPLITTING
! second coefficent is meaningless in cyclic case
IF (HLBCX(1) == 'CYCL' .AND. L2D .AND. SIZE(PBFB,1) .GE. 2  ) PBFB(2,:,:)=1.
IF (HLBCX(1) == 'CYCL' .AND. .NOT.(L2D) .AND. LWEST_ll(HSPLITTING='B') .AND. SIZE(PBFB,2) .GE. 2 ) &
    PBFB(:,2,:)=1.
IF (HLBCY(1) == 'CYCL' .AND. .NOT.(L2D)  .AND. SIZE(PBFB,1) .GE. 2 ) PBFB(2,:,:)=1.
!++cb
!
DEALLOCATE(ZEIGEN_ll)
!
!
!------------------------------------------------------------------------------
!*       8.    INITIALIZATION OF THE TRIGS AND IFAX ARRAYS FOR THE FFT 
!              -------------------------------------------------------
!
!        8.1 x lateral boundary conditions
!                  
CALL SET99(PTRIGSX,KIFAXX,IIMAX_ll)
!
! test on the value of KIMAX: KIMAX must be factorizable as a product
! of powers of 2,3 and 5. KIFAXX(10) is equal to IIMAX if the decomposition 
! is correct, then KIFAXX(1) contains the number of decomposition factors
! of KIMAX.
!
IF (KIFAXX(10) /=  IIMAX_ll) THEN
      WRITE(UNIT=ILUOUT,FMT="('   ERROR',/,                               &
      &'     : THE FORM OF THE FFT USED FOR THE INVERSION OF THE FLAT ',/,&
      &'    OPERATOR REQUIRES THAT KIMAX MUST BE FACTORIZABLE'         ,/,&
      & '         AS A PRODUCT OF POWERS OF 2, 3 AND 5.')")
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','ELEC_TRIDZ','')
END IF 
!
IF (HLBCX(1) /= 'CYCL') THEN
!
!     extra trigs for shifted (co) sine transform (FFT55)
!
  INN=2*(IIMAX_ll)
  ZDEL=ASIN(1.0)/REAL(IIMAX_ll)
  DO JI=1,IIMAX_ll
    ZANGLE=REAL(JI)*ZDEL
    PTRIGSX(INN+JI)=SIN(ZANGLE)
  END DO
!
ENDIF
!
!       8.2 y lateral boundary conditions
!
IF (.NOT. L2D) THEN 
  CALL SET99(PTRIGSY,KIFAXY,IJMAX_ll)
  !
  ! test on the value of KJMAX: KJMAX must be factorizable as a product
  ! of powers of 2,3 and 5. KIFAXY(10) is equal to IJMAX_ll if the decomposition 
  ! is correct, then KIFAXX(1) contains the number of decomposition factors
  ! of IIMAX_ll.
  !
  IF (KIFAXY(10) /= IJMAX_ll) THEN
      WRITE(UNIT=ILUOUT,FMT="('   ERROR',/,                               &
      &'     : THE FORM OF THE FFT USED FOR THE INVERSION OF THE FLAT ',/,&
      &'    OPERATOR REQUIRES THAT KJMAX MUST BE FACTORIZABLE'         ,/,&
      & '         AS A PRODUCT OF POWERS OF 2, 3 AND 5.')")
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','ELEC_TRIDZ','')
 END IF 
 !
 !       8.3 top boundary conditions
 !
 PEPOTFW_TOP(:,:) = PEPOTFW_TOP(:,:) / ZDZM(IKE+1)

 !
 !     other cases 
 !
 IF (HLBCY(1) /= 'CYCL') THEN                
 !
 !     extra trigs for shifted (co) sine transform
 !
   INN=2*(IJMAX_ll)
   ZDEL=ASIN(1.0)/REAL(IJMAX_ll)
   DO JJ=1,IJMAX_ll
     ZANGLE=REAL(JJ)*ZDEL
     PTRIGSY(INN+JJ)=SIN(ZANGLE)
   END DO
 !
 ENDIF
 !
ENDIF
!
!------------------------------------------------------------------------------
!
END SUBROUTINE ELEC_TRIDZ
