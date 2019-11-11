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
!     ######spl
      MODULE MODI_RICHARDSON
!     ######################
!
INTERFACE
!
      SUBROUTINE RICHARDSON(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,          &
      PRHODJ,PTHETAV,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,PTRIGSX,PTRIGSY,    &
      KIFAXX,KIFAXY,KITR,KTCOUNT,PRELAX,PY,PPHI)
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHETAV   ! virtual potential temp. at 
                                                 ! time t
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
INTEGER, INTENT(IN) :: KTCOUNT                   ! present iteration of the 
                                                 ! model temporal counter
!
REAL, INTENT(IN) :: PRELAX                       ! relaxation coefficient 
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PY         ! RHS of the equation       
!          
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PPHI    ! solution of the equation 
!
END SUBROUTINE RICHARDSON
!
END INTERFACE
!
END MODULE MODI_RICHARDSON
!     ######spl
      SUBROUTINE RICHARDSON(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,          &
      PRHODJ,PTHETAV,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF,PTRIGSX,PTRIGSY,    &
      KIFAXX,KIFAXY,KITR,KTCOUNT,PRELAX,PY,PPHI)
!     ########################################################################
!
!!****  *RICHARDSON * - solve the elliptic pressure equation by the Richardson's 
!!      method
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to solve the elliptic pressure equation
!     by the Richardson's method preconditioned by the flat Laplacian.
!     
!!**  METHOD
!!    ------
!!      
!!      The equation to be solved reads:
!!
!!             Q (PHI) = Y
!!
!!    where Q is the quasi-Laplacian ( subroutine QLAP ) and PHI the pressure
!!    function. 
!!    We precondition the problem by the operator F :
!!             -1               -1 
!!            F   * Q (PHI) = F    (Y)
!!    F represents the flat quasi-laplacian ie. without orography. Its 
!!    inversion is realized in the routine FLAT_INV.
!!    This equation is solved with the Richardson's method:
!!
!!           (n+1)      (n)         (  -1                -1    (n)   )
!!        PHI      = PHI    + RELAX ( F   ( Y ) - F   * Q ( PHI    ) )
!!                                  (                                )
!!
!!    EXTERNAL
!!    --------
!!      Subroutine GDIV: compute J times the divergence of 1/J times a vector
!!      Function QLAP: compute the complete quasi-Laplacian Q
!!      Subroutine FLAT_INV : invert the flat quasi-laplacien F
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF: model configuration
!!        CCONF: option for the first time step 
!!               of the segment ( start or restart)
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine RICHARDSON)
!!
!!    AUTHOR
!!    ------
!!	P. HÅreil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/07/94 
!!                  14/01/97 Durran anelastic equation (Stein,Lafore)
!!                  15/06/98  (D.Lugato, R.Guivarch) Parallelisation
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CONF
USE MODI_GDIV
USE MODI_QLAP
USE MODI_FLAT_INV
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHETAV   ! virtual potential temp. at 
                                                 ! time t
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
INTEGER, INTENT(IN) :: KTCOUNT                   ! present iteration of the 
                                                 ! model temporal counter
!
REAL, INTENT(IN) :: PRELAX                       ! relaxation coefficient 
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PY         ! RHS of the equation       
!          
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PPHI    ! solution of the equation 
!
!*       0.2   declarations of local variables
!
INTEGER :: JM                                    ! loop index 
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZF_1_Y 
                !   RHS of the preconditioned problem
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZCORREC
                !   iterative correction to the solution
!
REAL, DIMENSION(SIZE(PPHI,1),SIZE(PPHI,2),SIZE(PPHI,3)) :: ZWORK
                !   quasi-laplacien Q of the iterative solution
!
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE RHS OF THE PRECONDITIONED PROBLEM
!              ---------------------------------------------
!
CALL FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF, &  !  -1
             PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,PY,ZF_1_Y)           ! F  ( Y )
!
!-------------------------------------------------------------------------------
!
!*       2.    ITERATIVE LOOP
!              --------------
!
IF ( KTCOUNT < 1 .OR. ( KTCOUNT == 1 .AND. CCONF == 'START') ) THEN
  PPHI = ZF_1_Y   ! when no first guess is available, we take the solution
                  ! for the flat problem 
END IF
!
DO JM = 1,KITR        
!
  ZWORK = QLAP(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHETAV,PPHI)
                                                                 ! Q (PHI)
!
  ZCORREC = 0.
!
  CALL FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOM,PAF,PBF,PCF, & !  -1
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZWORK,ZCORREC )     ! F  * Q (PHI)
!
!                                                       -1       -1
!    update the iterative solution PHI = PHI + relax* (F  (Y) - F  * Q (PHI))
  PPHI = PPHI + PRELAX * (ZF_1_Y - ZCORREC) 
!        
END DO
!-------------------------------------------------------------------------------
!
END SUBROUTINE RICHARDSON
