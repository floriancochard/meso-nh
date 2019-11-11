!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_CONRESOLZ
!     ####################
!
INTERFACE
!
      SUBROUTINE CONRESOLZ(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV, &
      PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,         &
      KITR,PY,PPHI,                                                            &
      PBFB,&
      PBF_SXP2_YP1_Z) !JUAN Z_SPLITING
!  
IMPLICIT NONE
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type 
!
                                                 ! Metric coefficients:
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! d*xx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! d*yy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX      ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY      ! d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ      ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ    ! density of reference * J
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHETAV   ! virtual pot. temp. at time t
!
REAL, INTENT(IN) :: PDXHATM                     ! mean grid increment in the x
                                                ! direction
REAL, INTENT(IN) :: PDYHATM                     ! mean grid increment in the y
                                                ! direction
!
REAL, DIMENSION (:), INTENT(IN) :: PRHOM         !  XRHODJ mean on the X Y plane
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
INTEGER, INTENT(IN) :: KITR                      ! number of iterations for the
                                                 ! pressure solver
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PY         ! RHS of the equation
!          
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PPHI    ! solution of the equation 
!
!JUAN Z_SPLITING
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBFB       ! elements of the tri-diag. b-slide 
                                                 ! matrix in the pressure eq.
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF_SXP2_YP1_Z       ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
!JUAN Z_SPLITING
END SUBROUTINE CONRESOLZ
!
END INTERFACE
!
END  MODULE MODI_CONRESOLZ
!
!
!
!     #########################################################################
      SUBROUTINE CONRESOLZ(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV, &
      PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,         &
      KITR,PY,PPHI,                                                            &
      PBFB,&
      PBF_SXP2_YP1_Z) !JUAN Z_SPLITING
!     #########################################################################
!
!!****  *CONRESOLZ * - solve an elliptic equation by the conjugate residual
!!       method
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to solve an elliptic equation using 
!     the preconditioned conjugate residual (CR) method. This is a version
!     of the scheme proposed by Skamarock, Smolarkiewicz and Klemp (MWR, 1997).
!     
!!**  METHOD
!!    ------
!!      The equation to be solved reads:
!!
!!             Q (PHI) = Y
!!
!!    where Q is the quasi-Laplacian ( subroutine QLAP ) and PHI the pressure
!!    function.
!!    We precondition the problem by the operator F :
!!             -1               -1 
!!            F   * Q (PHI) = F    (Y)
!!    F represents the flat Laplacian ie. without orography. Its inversion is 
!!    realized in the routine FLAT_INVZ. This equation is solved with a Conjugate
!!    Residual method.
!!    The initial guess is given by the pressure at the previous time step.
!!    The resolution stops after ITR iterations of the solver.
!!
!!    EXTERNAL
!!    --------
!!      Subroutine GDIV: compute J times the divergence of 1/J times a vector
!!      Function QLAP: compute the complete quasi-Laplacian Q
!!      Subroutine FLAT_INVZ : invert the flat quasi-laplacien F
!!      Function DOTPROD: compute the dot product of 2 vectors
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODI_GDIV: interface for the subroutine GDIV
!!      Module MODI_QLAP: interface for the function QLAP
!!      Module MODI_FLAT_INVZ: interface for the subroutine FLAT_INVZ
!!      Module MODI_DOTPROD: interface for the function DOTPROD
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine CONRESOL)
!!      Skamarock, Smolarkiewicz and Klemp (1997) MWR
!!      
!!    AUTHOR
!!    ------
!!	J.-P. Pinty *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/08/99 
!!      J.-P. Pinty & P. Jabouille
!!                  11/07/00 bug in ZALPHA
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_GDIV
USE MODI_QLAP
USE MODI_FLAT_INVZ
USE MODI_DOTPROD
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type 
!
                                                 ! Metric coefficients:
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! d*xx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! d*yy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX      ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY      ! d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ      ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ    ! density of reference * J
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHETAV   ! virtual pot. temp. at time t
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
INTEGER, INTENT(IN) :: KITR                      ! number of iterations for the
                                                 ! pressure solver
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PY         ! RHS of the equation
!          
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PPHI    ! solution of the equation 
!JUAN Z_SPLITING
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBFB       ! elements of the tri-diag. b-slide 
                                                 ! matrix in the pressure eq.
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF_SXP2_YP1_Z       ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
!JUAN Z_SPLITING
!
!*      0.2    declarations of local variables
!
INTEGER :: JM                                    ! loop index   
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZDELTA, ZKSI  
     ! array containing the auxilary fields DELTA and KSI of the CR method
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZP, ZQ  
     ! array containing the auxilary fields P and Q of the CR method
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZRESIDUE
     ! array containing the error field at each iteration Q(PHI) - Y
!
REAL :: ZALPHA, ZLAMBDA      ! amplitude of the descent in the Conjugate
                             ! directions
REAL :: ZDOT_DELTA           ! dot product of ZDELTA by itself
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
!                             
!*       1.1    compute the vector: r^(0) =  Q(PHI) - Y
!    
ZRESIDUE = QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV,PPHI) - PY
!
!*       1.2    compute the vector: p^(0) = F^(-1)*( Q(PHI) - Y )
!    
CALL FLAT_INVZ(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,  &
                     PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZRESIDUE,ZP,&
                     PBFB,&
                     PBF_SXP2_YP1_Z) !JUAN Z_SPLITING 
!JUAN      print*, "size ZP=",SIZE(ZP)
!JUAN      print*, "size ZRESIDUE=",SIZE(ZRESIDUE)
!
!*       1.3    compute the vector: delta^(0) = Q ( p^(0) )
!
ZDELTA = QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV,ZP)
!
!-------------------------------------------------------------------------------
!
!*       2.    ITERATIVE LOOP
!              --------------
!
DO JM = 1,KITR                                                         
!
!*       2.1    compute the step LAMBDA
!
  ZDOT_DELTA = DOTPROD(ZDELTA,  ZDELTA,HLBCX,HLBCY)            ! norm of DELTA
  ZLAMBDA  = - DOTPROD(ZRESIDUE,ZDELTA,HLBCX,HLBCY) / ZDOT_DELTA
!
!*       2.2    update the pressure function PHI
!
  PPHI = PPHI + ZLAMBDA * ZP
!
!
  IF( JM == KITR ) EXIT
!
!
!*       2.3    update the residual error: r
!
  ZRESIDUE = ZRESIDUE + ZLAMBDA * ZDELTA
!                        
!*       2.4    compute the vector: q = F^(-1)*( Q(PHI) - Y )
!    
CALL FLAT_INVZ(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,  &
                 PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZRESIDUE,ZQ,    &
                 PBFB,&
                 PBF_SXP2_YP1_Z) !JUAN Z_SPLITTING
!   
!*       2.5     compute the auxiliary field: ksi = Q ( q )
!
  ZKSI= QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV,ZQ)
!                                         -1
!*       2.6     compute the step ALPHA
!  
  ZALPHA = - DOTPROD(ZKSI,ZDELTA,HLBCX,HLBCY) / ZDOT_DELTA   ! lambda
!
!*       2.7     update p and DELTA
!
  ZP     = ZQ   + ZALPHA * ZP
  ZDELTA = ZKSI + ZALPHA * ZDELTA
!
END DO              ! end of the loop for the iterative solver
!
!  
!-------------------------------------------------------------------------------
!
END SUBROUTINE CONRESOLZ
