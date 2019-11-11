!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 solver 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_FLAT_INV
!     ####################
!
INTERFACE
!
      SUBROUTINE FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF, &
      PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,PY,PF_1_Y)
!
!  
IMPLICIT NONE
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type 
!
REAL, INTENT(IN) :: PDXHATM                     ! mean grid increment in the x
                                                ! direction
REAL, INTENT(IN) :: PDYHATM                     ! mean grid increment in the y
                                                ! direction
!
REAL, DIMENSION (:), INTENT(IN) :: PRHOM         !  mean of XRHODJ on the plane x y 
                                                 !  localized at a mass level
!
REAL, DIMENSION(:), INTENT(IN)     :: PAF,PCF    ! vectors giving the nonvanishing
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF        ! elements of the tri-diag. 
                                                 ! matrix in the pressure eq.
!
                                                 ! arrays of sin or cos values
                                                 ! for the FFT :
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSX        ! - along x
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSY        ! - along y
!
                                                 ! decomposition in prime 
                                                 ! numbers for the FFT:
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXX      ! - along x
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXY      ! - along y
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PY         ! RHS of the equation
!                                                         
REAL, DIMENSION(:,:,:), INTENT(OUT):: PF_1_Y     ! solution of the equation
!
END SUBROUTINE FLAT_INV
!
END INTERFACE
!
END MODULE MODI_FLAT_INV
!     ######################################################################
      SUBROUTINE FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,   &
      PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,PY,PF_1_Y)
!     ######################################################################
!
!!****  *FLAT_INV * - Invert the flat quasi-laplacian operator
!!
!!    PURPOSE
!!    -------
!       This routine solves the following equation:
!                    F ( F_1_Y ) =  Y 
!     where F represents the quasi-laplacian without orography. The solution is 
!     F_1_Y.
!
!!**  METHOD
!!    ------
!!      The horizontal part of F is inverted with a FFT transform. For each 
!!    horizontal direction, the FFT form depends on the lateral boundary 
!!    conditions :
!!    - CRAY intrinsic function RFFTMLT in the cyclic case 
!!    - fast cosine transform called FFT55 for all other boundary condtions.
!!    Then, in the wavenumber space, we invert for each 
!!    horizontal mode i,j a tridiagonal matrix by a classical double sweep 
!!    method.  The singular mean mode (i,j)=(0,0) corresponds to the 
!!    undetermination of the pressure to within a constant and is treated apart.
!!    To fix this degree of freedom, we set the horizontal mean value of the 
!!    pressure perturbation to 0 at the upper level of the model. 
!!       
!!    EXTERNAL
!!    --------
!!      Subroutine FFT55 : aplly multiple fast real staggered (shifted) 
!!    cosine transform
!!      Subroutine RFFTMLT : apply real-to-complex or complex-to-real Fast 
!!    Fourier Transform (FFT) on multiple input vectors.
!!      Subroutine FFT991  : equivalent to RFFTMLT
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT, JPVEXT: define the number of marginal points out of the 
!!        physical domain along horizontal and vertical directions respectively
!!      Module MODD_CONF: model configurations
!!        L2D: logical for 2D model version
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (subroutine FLAT_INV)     
!!
!!    AUTHOR
!!    ------
!!	P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/07/94 
!!      Revision  Jabouille (juillet 96) replace the CRAY intrinsic function
!!                             RFFTMLT by the arpege routine FFT991
!!                  17/07/97  ( J. Stein and V. Masson) initialize the corner
!!                             verticals
!!                  17/07/97  ( J. Stein and V. Masson) initialize the corner
!!                             verticals
!!      Revision  Jabouille (septembre 97) suppress the particular case for
!!                                         tridiagonal inversion
!!                Stein     ( January 98 ) faster computation for the unused 
!!                                points under the ground and out of the domain  
!!      Modification   Lugato, Guivarch (June 1998) Parallelisation
!!                     Escobar, Stein   (July 2000) optimisation
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CONF
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODI_FFT55
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type 
!
REAL, INTENT(IN) :: PDXHATM                     ! mean grid increment in the x
                                                ! direction
REAL, INTENT(IN) :: PDYHATM                     ! mean grid increment in the y
                                                ! direction
!
REAL, DIMENSION (:), INTENT(IN) :: PRHOM         !  mean of XRHODJ on the plane x y 
                                                 !  localized at a mass level
!
REAL, DIMENSION(:), INTENT(IN)     :: PAF,PCF    ! vectors giving the nonvanishing
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF        ! elements of the tri-diag. 
                                                 ! matrix in the pressure eq.
!
                                                 ! arrays of sin or cos values
                                                 ! for the FFT :
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSX        ! - along x
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSY        ! - along y
!
                                                 ! decomposition in prime 
                                                 ! numbers for the FFT:
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXX      ! - along x
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXY      ! - along y
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PY         ! RHS of the equation
!                                                         
REAL, DIMENSION(:,:,:), INTENT(OUT):: PF_1_Y     ! solution of the equation
!
!*       0.2   declaration of local variables
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: ZY ! work array to store 
                                                    ! the RHS of the equation
!
!REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: ZWORK ! work array used by 
! the FFT. It has been enlarged in order to be sufficient for 2D and 3D cases
!
REAL, DIMENSION(SIZE(PBF,1),SIZE(PBF,2),SIZE(PBF,3)) :: ZAF ! work array to
!                                                        ! expand PAF
INTEGER :: IIB          ! indice I for the first inner mass point along x
INTEGER :: IIE          ! indice I for the last inner mass point along x
INTEGER :: IIMAX        ! number of inner mass points along the x direction
INTEGER :: IJB          ! indice J for the first inner mass point along y
INTEGER :: IJE          ! indice J for the last inner mass point along y
INTEGER :: IJMAX        ! number of inner mass points along the y direction
INTEGER :: IKB          ! indice K for the first inner mass point along z
INTEGER :: IKE          ! indice K for the last inner mass point along z
INTEGER :: IKU          ! size of the arrays along z
INTEGER :: IKMAX        ! number of inner mass points along the z direction
!
REAL :: ZDXM2,ZDYM2     ! respectively equal to PDXHATM*PDXHATM 
                        ! and PDYHATM*PDYHATM                         
INTEGER :: JI,JJ,JK     ! loop indexes along x, y, z respectively
!
!
INTEGER :: IIE_INT,IJE_INT   ! highest indice I  and J values for the x y modes.
                             ! They  depend on the l.b.c. !
!
INTEGER :: ILOTX,ILOTY ! number of data vectors along x, y resp. computed 
                       ! in parallel during the FFT process 
!
INTEGER :: INC1X,INC1Y ! increment within each data vector for the FFT along
                       ! x, y resp.
!
INTEGER :: INC2X,INC2Y ! increment between the start of one data vector and 
                       ! the next for the FFT along x,y resp.
!
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORKX ! work array used by
! the FFT. It has been enlarged in order to be sufficient for 2D and 3D cases
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORKY ! work array used by
! the FFT. It has been enlarged in order to be sufficient for 2D and 3D cases
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZGAM
                                                        ! intermediate arrays
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZBETX            ! for the tridiag.
                                                        ! matrix inversion
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_X  ! array in X slices distribution
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_Y  ! array in Y slices distribution
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_YR ! array in Y slices distribution
!
INTEGER :: IINFO_ll           ! return code of parallel routine
!
INTEGER :: IIX,IJX,IIY,IJY ! dimensions of the extended x or y slices subdomain
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_YT  ! array in Y slices distribution transpose
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_YRT ! array in Y slices distribution transpose
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE LOOP BOUNDS
!              -------------------
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll('X',IIX,IJX)
CALL GET_DIM_EXT_ll('Y',IIY,IJY)
IIMAX = IIX-2*JPHEXT
IJMAX = IJY-2*JPHEXT
!
IKU=SIZE(PY,3)
IKB=1+JPVEXT
IKE=IKU - JPVEXT
IKMAX=IKE-IKB+1
!
!!
ALLOCATE(ZBAND_X(IIX,IJX,IKU))
ALLOCATE(ZBAND_Y(IIY,IJY,IKU))
ALLOCATE(ZBAND_YR(IIY,IJY,IKU))
ALLOCATE(ZWORKX(IIX,IJX,IKU))
ALLOCATE(ZWORKY(IIY,IJY,IKU))
ALLOCATE(ZBETX(IIY,IJY))
ALLOCATE(ZGAM(IIY,IJY,IKU))
IF (.NOT. L2D) THEN
  ALLOCATE(ZBAND_YT(IJY,IIY,IKU))
  ALLOCATE(ZBAND_YRT(IJY,IIY,IKU))
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE ARRAY INCREMENTS FOR THE FFT
!              ----------------------------------------
!
IF(.NOT. L2D) THEN
!
  ILOTX = IJX*IKU
  INC1X = 1
  INC2X = IIX
!
  ILOTY = IIY*IKU
  INC1Y = 1
  INC2Y = IJY
!
ELSE
!
  ILOTX = IKU
  INC1X = 1
  INC2X = IIX*IJX
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      3.    FORM HOMOGENEOUS BOUNDARY CONDITIONS FOR A NONCYCLIC CASE
!             ---------------------------------------------------------
!
!
!*      3.1 copy the RHS in a local array REMAP functions will shift the indices for the FFT
!
PF_1_Y = 0.
ZY = PY
!
!*      3.2 form homogeneous boundary condition used by the FFT for non-periodic
!           cases
!
!  modify the RHS in the x direction
!
IF (HLBCX(1) /= 'CYCL') THEN
!
  IF (LWEST_ll(HSPLITTING='B')) THEN
    DO JK=IKB,IKE                    
      DO JJ = IJB, IJE
        ZY(IIB,JJ,JK)     = ZY(IIB,JJ,JK)     + PY(IIB-1,JJ,JK)
      END DO
    END DO
  END IF
!
  IF (LEAST_ll(HSPLITTING='B')) THEN
    DO JK=IKB,IKE
      DO JJ = IJB, IJE
        ZY(IIE,JJ,JK)     = ZY(IIE,JJ,JK)     - PY(IIE+1,JJ,JK)
      END DO
    END DO
  END IF
END IF
!
! modify the RHS in the same way along y
!
IF (HLBCY(1) /= 'CYCL'.AND. (.NOT. L2D)) THEN
  IF (LSOUTH_ll(HSPLITTING='B')) THEN
    DO JK=IKB,IKE
      DO JI = IIB, IIE
        ZY(JI,IJB,JK)     = ZY(JI,IJB,JK)     + PY(JI,IJB-1,JK)
      END DO
    END DO
  END IF
!
  IF (LNORTH_ll(HSPLITTING='B')) THEN
    DO JK=IKB,IKE
      DO JI = IIB, IIE
        ZY(JI,IJE,JK)     = ZY(JI,IJE,JK)     - PY(JI,IJE+1,JK)
      END DO
    END DO
  END IF
END IF
!
!
!*      3.3 2way structure -> xslice structure, + data shift
!
ZBAND_X=0.
CALL REMAP_2WAY_X_ll(ZY,ZBAND_X,IINFO_ll)
!
!  
!-------------------------------------------------------------------------------
!
!*      4.    APPLY A REAL TO COMPLEX FFT
!             ---------------------------
!
!
IF (HLBCX(1) == 'CYCL') THEN
  CALL FFT991(ZBAND_X(1,1,IKB-1),ZWORKX,PTRIGSX,KIFAXX,INC1X,INC2X,  &
              IIMAX,ILOTX,-1 )
ELSE
  CALL FFT55(ZBAND_X(1,1,IKB-1),ZWORKX,PTRIGSX,KIFAXX,INC1X,INC2X,  &
             IIMAX,ILOTX,-1 )
END IF  
!
!
ZBAND_Y=0.
CALL REMAP_X_Y_ll(ZBAND_X,ZBAND_Y,IINFO_ll)
!
IF (.NOT. L2D) THEN                                 
!
! array transposition I --> J
!
  CALL FAST_TRANSPOSE(ZBAND_Y,ZBAND_YT,IIY,IJY,IKU)
!
  IF (HLBCY(1) == 'CYCL') THEN
      CALL FFT991(ZBAND_YT(1,1,IKB-1),ZWORKY,PTRIGSY,KIFAXY,INC1Y,INC2Y,  &
                  IJMAX,ILOTY,-1 )
  ELSE
      CALL FFT55(ZBAND_YT(1,1,IKB-1),ZWORKY,PTRIGSY,KIFAXY,INC1Y,INC2Y,  &
                  IJMAX,ILOTY,-1 )
  END IF 
!
END IF                                        
!
!       singular matrix case : the last term is computed by setting the 
!       average of the pressure field equal to zero.
IF (LWEST_ll(HSPLITTING='Y')) THEN
  IF (L2D) THEN
    ZBAND_Y(1,1,IKE+1)=0
  ELSE
    ZBAND_YT(1,1,IKE+1)=0.
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*      5.    MATRIX INVERSION FOR THE FLAT OPERATOR
!             --------------------------------------
!
CALL FAST_SPREAD(PAF,ZAF,IIY,IJY,IKU)
!
IF (LWEST_ll(HSPLITTING='Y')) THEN
   ZAF(1,1,IKE+1)=0. !singular matrix corresponding to the horizontal average
END IF
!
IF (L2D) THEN
  CALL FAST_SUBSTITUTION_2D(ZBAND_YR,ZBETX,PBF,ZGAM,PCF,ZAF &
                           ,ZBAND_Y,IIY,IJY,IKU)
ELSE
  CALL FAST_SUBSTITUTION_3D(ZBAND_YRT,ZBETX,PBF,ZGAM,PCF,ZAF &
                      ,ZBAND_YT,IIY,IJY,IKU)
END IF
!  
!
!-------------------------------------------------------------------------------
!
!*      6.    APPLY A COMPLEX TO REAL FFT
!             ---------------------------
!
!  
IF (.NOT. L2D) THEN
  IF (HLBCY(1) == 'CYCL') THEN
    CALL FFT991( ZBAND_YRT(1,1,IKB-1),ZWORKY,PTRIGSY,KIFAXY,INC1Y,INC2Y,  &
                  IJMAX,ILOTY,+1 )
  ELSE
      CALL FFT55( ZBAND_YRT(1,1,IKB-1),ZWORKY,PTRIGSY,KIFAXY,INC1Y,INC2Y,  &
                  IJMAX,ILOTY,+1 )
  END IF 
  ! array transposition J --> I
  CALL FAST_TRANSPOSE(ZBAND_YRT,ZBAND_YR,IJY,IIY,IKU)
ENDIF
!
! Transposition Y-> X
!
ZBAND_X=0.
CALL REMAP_Y_X_ll(ZBAND_YR,ZBAND_X,IINFO_ll)
!
!
IF (HLBCX(1) == 'CYCL') THEN
    CALL FFT991( ZBAND_X(1,1,IKB-1),ZWORKX,PTRIGSX,KIFAXX,INC1X,INC2X,  &
                IIMAX,ILOTX,+1 )
ELSE
    CALL FFT55( ZBAND_X(1,1,IKB-1),ZWORKX,PTRIGSX,KIFAXX,INC1X,INC2X,  &
                IIMAX,ILOTX,+1 )
END IF  
!
!-------------------------------------------------------------------------------
!
!*      7.    RETURN TO A NON HOMOGENEOUS NEUMAN CONDITION FOR NON-CYCLIC CASES
!             -----------------------------------------------------------------
!
!*      7.1   Transposition + shift X -> 2way
!
CALL REMAP_X_2WAY_ll(ZBAND_X,PF_1_Y,IINFO_ll)
!
!*      7.2 complete the lateral boundaries
!
IF (HLBCX(1) /= 'CYCL') THEN
!
!*      7.2.1 return to a non-homogeneous case in the x direction
! 
  ZDXM2 = PDXHATM*PDXHATM
!
  IF (LWEST_ll(HSPLITTING='B')) THEN
    DO JK=IKB,IKE
      DO JJ = IJB,IJE
        PF_1_Y(IIB-1,JJ,JK) = PF_1_Y(IIB,JJ,JK) - PY(IIB-1,JJ,JK)*ZDXM2/PRHOM(JK)
      END DO
    END DO
  END IF
!
  IF (LEAST_ll(HSPLITTING='B')) THEN
    DO JK=IKB,IKE
      DO JJ = IJB,IJE
        PF_1_Y(IIE+1,JJ,JK) = PF_1_Y(IIE,JJ,JK) + PY(IIE+1,JJ,JK)*ZDXM2/PRHOM(JK)
      END DO
    END DO
  END IF
!
!    we set the solution at the corner point by the condition:
!    dxm ( P ) = 0
  IF (LWEST_ll(HSPLITTING='B')) THEN
    DO JJ = IJB,IJE
      PF_1_Y(IIB-1,JJ,IKB-1) = PF_1_Y(IIB,JJ,IKB-1) 
      PF_1_Y(IIB-1,JJ,IKE+1) = PF_1_Y(IIB,JJ,IKE+1)
    END DO
  END IF
  IF (LEAST_ll(HSPLITTING='B')) THEN
    DO JJ = IJB,IJE
      PF_1_Y(IIE+1,JJ,IKB-1) = PF_1_Y(IIE,JJ,IKB-1)
      PF_1_Y(IIE+1,JJ,IKE+1) = PF_1_Y(IIE,JJ,IKE+1)
    END DO
  END IF
!
ELSE
!
!*      7.2.2 periodize the pressure function field along the x direction
!
! in fact this part is useless because it is done in the routine
! REMAP_X_2WAY.
!
END IF
!
IF (.NOT.L2D) THEN
  IF (HLBCY(1) /= 'CYCL') THEN
!
!*      7.2.3 return to a non-homogeneous case in the y direction
!
    ZDYM2 = PDYHATM*PDYHATM 
!
    IF (LSOUTH_ll(HSPLITTING='B')) THEN
      DO JK=IKB,IKE
        DO JI = IIB,IIE
          PF_1_Y(JI,IJB-1,JK) = PF_1_Y(JI,IJB,JK) - PY(JI,IJB-1,JK)*ZDYM2/PRHOM(JK)
        END DO
      END DO
    END IF
!
    IF (LNORTH_ll(HSPLITTING='B')) THEN
      DO JK=IKB,IKE
        DO JI = IIB,IIE
          PF_1_Y(JI,IJE+1,JK) = PF_1_Y(JI,IJE,JK) + PY(JI,IJE+1,JK)*ZDYM2/PRHOM(JK)
        END DO
      END DO
    END IF
!    we set the solution at the corner point by the condition:
!    dym ( P )  = 0
!
    IF (LSOUTH_ll(HSPLITTING='B')) THEN
      DO JI = IIB,IIE
        PF_1_Y(JI,IJB-1,IKB-1) = PF_1_Y(JI,IJB,IKB-1)
        PF_1_Y(JI,IJB-1,IKE+1) = PF_1_Y(JI,IJB,IKE+1)
      END DO
    END IF
!
    IF (LNORTH_ll(HSPLITTING='B')) THEN
      DO JI = IIB,IIE
        PF_1_Y(JI,IJE+1,IKB-1) = PF_1_Y(JI,IJE,IKB-1)
        PF_1_Y(JI,IJE+1,IKE+1) = PF_1_Y(JI,IJE,IKE+1)
      END DO
    END IF
  ELSE
!
!*      7.2.4 periodize the pressure function field along the y direction
!
!
! in fact this part is useless because it is done in the routine
! REMAP_X_2WAY.
!
  END IF
!
END IF
!
IF (.NOT. L2D .AND. HLBCX(1)/='CYCL' .AND. HLBCY(1)/='CYCL') THEN
! the following verticals are not used
  IF ( (LWEST_ll(HSPLITTING='B')).AND.(LSOUTH_ll(HSPLITTING='B')) ) THEN
    PF_1_Y(IIB-1,IJB-1,:)=0.
  END IF
!
  IF ( (LWEST_ll(HSPLITTING='B')).AND.(LNORTH_ll(HSPLITTING='B')) ) THEN
    PF_1_Y(IIB-1,IJE+1,:)=0.
  END IF
!
  IF ( (LEAST_ll(HSPLITTING='B')).AND.(LSOUTH_ll(HSPLITTING='B')) ) THEN
    PF_1_Y(IIE+1,IJB-1,:)=0.
  END IF
!
  IF ( (LEAST_ll(HSPLITTING='B')).AND.(LNORTH_ll(HSPLITTING='B')) ) THEN
    PF_1_Y(IIE+1,IJE+1,:)=0.
  END IF
END IF
!
DEALLOCATE(ZBAND_X)
DEALLOCATE(ZBAND_Y)
IF (.NOT. L2D) THEN
  DEALLOCATE(ZBAND_YT)
  DEALLOCATE(ZBAND_YRT)
END IF
DEALLOCATE(ZBAND_YR)
DEALLOCATE(ZWORKX)
DEALLOCATE(ZWORKY)
DEALLOCATE(ZBETX)
DEALLOCATE(ZGAM)
!
!-------------------------------------------------------------------------------
!
CONTAINS
  SUBROUTINE FAST_TRANSPOSE(PX,PXT,KNI,KNJ,KNK)
    INTEGER                      :: KNI,KNJ,KNK ! 3D dimension of X and XT
    REAL, DIMENSION(KNI*KNJ,KNK) :: PX
    REAL, DIMENSION(KNJ*KNI,KNK) :: PXT
    !
    INTEGER                      :: IJI,II,IJ,IIJ ! index in array X and XT
    INTEGER                      :: JK
!
    DO JK=1,KNK
       ! PERMUTATION(PX,PXT)
       !CDIR NODEP
       !OCL NOVREC
       DO IJI = 1, KNJ*KNI
          ! I,J Indice in XT array from linearised index IJI
          II   = 1 +    (IJI-1)/KNJ
          IJ   = IJI - (II-1)*KNJ 
          ! linearised index in X
          IIJ = II + (IJ-1)*KNI 
          ! transposition
          PXT(IJI,JK) = PX(IIJ,JK)
       
       END DO
    END DO
!    
END SUBROUTINE FAST_TRANSPOSE

SUBROUTINE FAST_SUBSTITUTION_3D(PBAND_YR,PBETX,PPBF,PGAM,PPCF,PAF &
                            ,PBAND_Y,KIY,KJY,KKU)
INTEGER                        :: KIY,KJY,KKU
REAL, DIMENSION (KIY*KJY,KKU)  :: PBAND_YR,PBAND_Y,PPBF,PGAM,PAF
REAL, DIMENSION (KIY*KJY)      :: PBETX
REAL, DIMENSION (KKU)          :: PPCF
INTEGER                        :: JK
!
!
!       initialization 
!
!
PBAND_YR = 0.0
PBETX(:) = PPBF(:,IKB-1)
PBAND_YR(:,IKB-1) = PBAND_Y(:,IKB-1)  &
                                      / PBETX(:)
!
!        decomposition and forward substitution
!
DO JK = IKB,IKE+1
  PGAM(:,JK) = PPCF(JK-1) / PBETX(:)
!
  PBETX(:) = PPBF(:,JK) -              &      
                 PAF(:,JK)*PGAM(:,JK)
!
  PBAND_YR(:,JK) = ( PBAND_Y(:,JK) -  &
       PAF(:,JK)*PBAND_YR(:,JK- 1) )  &
                                /PBETX(:)
!
END DO
!
!       backsubstitution
!
DO JK = IKE,IKB-1,-1
  PBAND_YR(:,JK) = PBAND_YR(:,JK) -    &
          PGAM(:,JK+1)*PBAND_YR(:,JK+1)
END DO
!  
!
END SUBROUTINE FAST_SUBSTITUTION_3D
!
SUBROUTINE FAST_SUBSTITUTION_2D(PBAND_YR,PBETX,PPBF,PGAM,PPCF,PAF &
                            ,PBAND_Y,KIY,KJY,KKU)
INTEGER                        :: KIY,KJY,KKU
REAL, DIMENSION (KIY,KJY,KKU)  :: PBAND_YR,PBAND_Y,PPBF,PGAM,PAF
REAL, DIMENSION (KIY,KJY)      :: PBETX
REAL, DIMENSION (KKU)          :: PPCF
INTEGER                        :: JK
!
!
!       initialization 
!
!
PBAND_YR = 0.0
PBETX(:,1) = PPBF(:,1,IKB-1)
PBAND_YR(:,1,IKB-1) = PBAND_Y(:,1,IKB-1)  &
                                     / PBETX(:,1)
!
!        decomposition and forward substitution
!
DO JK = IKB,IKE+1
  PGAM(:,1,JK) = PPCF(JK-1) / PBETX(:,1)
!
  PBETX(:,1) = PPBF(:,1,JK) -              &      
                 PAF(:,1,JK)*PGAM(:,1,JK)
!
  PBAND_YR(:,1,JK) = ( PBAND_Y(:,1,JK) -  &
       PAF(:,1,JK)*PBAND_YR(:,1,JK- 1) )  &
                                /PBETX(:,1)
!
END DO
!
!       backsubstitution
!
DO JK = IKE,IKB-1,-1
  PBAND_YR(:,1,JK) = PBAND_YR(:,1,JK) -    &
          PGAM(:,1,JK+1)*PBAND_YR(:,1,JK+1)
END DO
!  
!
END SUBROUTINE FAST_SUBSTITUTION_2D

SUBROUTINE FAST_SPREAD(PTAB1D,PTAB3D,KIY,KJY,KKU)
INTEGER                        :: KIY,KJY,KKU
REAL, DIMENSION (KKU)          :: PTAB1D
REAL, DIMENSION (KIY*KJY,KKU)  :: PTAB3D

INTEGER                        :: JIJ,JK
!
DO JK=1,KKU
   DO JIJ=1,KIY*KJY
      PTAB3D(JIJ,JK) = PTAB1D(JK)
   ENDDO
ENDDO
!
END SUBROUTINE FAST_SPREAD
!
!------------------------------------------------------------------------------  
END SUBROUTINE FLAT_INV
