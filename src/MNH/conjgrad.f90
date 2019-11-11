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
      MODULE MODI_CONJGRAD
!     ####################
!
INTERFACE
!
      SUBROUTINE CONJGRAD(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV, &
      PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,         &
      KITR,PY,PPHI)
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHETAV   ! virtual potential temp. at time t
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
!
END SUBROUTINE CONJGRAD
!
END INTERFACE
!
END  MODULE MODI_CONJGRAD
!
!
!
!     #########################################################################
      SUBROUTINE CONJGRAD(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV, &
      PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,         &
      KITR,PY,PPHI)
!     #########################################################################
!
!!****  *CONJGRAD * - solve an elliptic equation by the conjugate gradient 
!!       method
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to solve an elliptic equation using 
!     the preditioned conjugate gradient (CG) method. This is a version of the 
!     CG called ORTHOMIN (Young and Jea 1980).
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
!!    realized in the routine FLAT_INV. This equation is solved with a Conjugate 
!!    Gradient method.
!!    The initial guess is given by the pressure at the previous time step.
!!    The resolution stops after ITR iterations of the solver.
!!
!!    EXTERNAL
!!    --------
!!      Subroutine GDIV: compute J times the divergence of 1/J times a vector
!!      Function QLAP: compute the complete quasi-Laplacian Q
!!      Subroutine FLAT_INV : invert the flat quasi-laplacien F
!!      Function DOTPROD: compute the dot product of 2 vectors
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODI_GDIV: interface for the subroutine GDIV
!!      Module MODI_QLAP: interface for the function QLAP
!!      Module MODI_FLAT_INV: interface for the subroutine FLAT_INV
!!      Module MODI_DOTPROD: interface for the function DOTPROD
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine CONJGRAD)
!!      Kapitza and Eppel (1992) Beit. Physik ...
!!      Young and Jea (1980) ....
!!      
!!    AUTHOR
!!    ------
!!	P. HÅreil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/07/94 
!!
!!                  14/01/97 Durran anelastic equation (Stein,Lafore)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_GDIV
USE MODI_QLAP
USE MODI_FLAT_INV
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHETAV   ! virtual potential temp. at time t
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
!
!*      0.2    declarations of local variables
!
INTEGER :: JM                                    ! loop index   
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZDELTA  
     ! array containing the auxilary field DELTA of the CG method
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZP  
     ! array containing the auxilary field P of the CG method
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZWORK ! work 
     ! array containing the source term to be multiplied by the F inverse
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZWORKD ! work 
     ! array containing the result of the F inversion * Q (DELTA)
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZWORKP ! work 
     ! array containing the result of the F inversion * Q (P)
!
REAL :: ZALPHA, ZLAMBDA      ! amplitude of the descent in the Conjugate
                             ! directions
REAL :: ZDOTPP               ! dot product of ZWORKP by itself
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
ZLAMBDA = 0.
ZP = 0.
!
!-------------------------------------------------------------------------------
!
!*       2.    ITERATIVE LOOP
!              --------------
!
DO JM = 1,KITR                                                         
!
!*       2.1    compute the new pressure function
!
  PPHI = PPHI + ZLAMBDA * ZP     !   the case JM =0 is special because 
                                 !   PPHI is not changed
!
!*       2.2    compute the auxiliary field DELTA
!
!                                -1
!   compute the vector DELTA = F  * ( Y - Q ( PHI ) )
!    
  ZWORK = QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV,PPHI) 
                                                               ! Q (PHI)
!                       
  ZWORK = PY - ZWORK                                           ! Y - Q (PHI)
!
  CALL FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF, &!  -1
                  PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZWORK,ZDELTA)   ! F  (Y - Q (PHI)))
!   
!*      2.3     compute the auxiliary field P
!
!                                         -1
!   compute the vector P = DELTA + alpha F  * Q ( DELTA ) 
!  
  IF (JM == 1) THEN
    ZP = ZDELTA          ! P = DELTA at the first solver iteration
  ELSE
    ZWORK  = QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,  &
                PDZZ,PRHODJ,PTHETAV,ZDELTA)                         ! Q ( DELTA )
    CALL FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,  & !  -1
            PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZWORK,ZWORKD)             ! F  * Q ( DELTA ) 
!
    ZALPHA = - DOTPROD(ZWORKD,ZWORKP,HLBCX,HLBCY)/ZDOTPP   ! ZWORKP,ZDOTPP come 
                     ! from the previous solver iteration  (section 2.4)
    ZP = ZDELTA + ZALPHA * ZP                              ! new vector P
!
  END IF
!
!*      2.4     compute LAMBDA
!
!
  ZWORK = QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,&
  PDZZ,PRHODJ,PTHETAV,ZP)                                       ! Q ( P )
  CALL FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,& !  -1
  PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZWORK,ZWORKP)                   ! F  * Q ( P )
!
!    store the scalar product to compute lambda and next P
  ZDOTPP = DOTPROD(ZWORKP,ZWORKP,HLBCX,HLBCY)   
!
  ZLAMBDA = DOTPROD(ZDELTA,ZWORKP,HLBCX,HLBCY) / ZDOTPP   ! lambda
!
!
END DO              ! end of the loop for the iterative solver
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE FINAL PRESSURE FUNCTION
!              -----------------------------------
!
PPHI = PPHI + ZLAMBDA * ZP
!  
!-------------------------------------------------------------------------------
!
END SUBROUTINE CONJGRAD
