!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!    ################## 
      MODULE MODI_CH_SIS
!!    ################## 
!!
!
INTERFACE
!!
SUBROUTINE CH_SIS(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI)
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT  
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN)  :: KMI      ! model number
END SUBROUTINE CH_SIS
!!
END INTERFACE
!!
END MODULE MODI_CH_SIS
!!    ######################################################################
      SUBROUTINE CH_SIS(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI)
!!    ######################################################################
!!
!!*** *CH_SIS*
!!
!!    PURPOSE
!!    -------
!     solve one time-step of the chemical differential euation d/dt C = f(C)
!!
!!**  METHOD
!!    ------
!!    Semi-Implicit Symmetric (SIS):
!!    y^n+1 = y^n + dt*(1-dt/2 J^n)^-1 f^n
!!
!!    REFERENCE
!!    ---------
!!    Ramaroson, R. A., M. Pirre, and D.Cariolle,
!!    A box model for on-line computations of diurnal variations
!!    in a 1-D model: potential for application in multidimensional cases,
!!    Ann. Geophysicae, 10, 416--428, 1992.
!!
!!    AUTHOR
!!    ------
!!    K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/04/95
!!    31/07/96 (K. Suhre) restructured
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_FCN
USE MODI_CH_JAC
USE MODI_CH_GAUSS
!%
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model number
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KVECNPT,KEQ)     :: ZDCDTVECT
REAL, DIMENSION(KEQ)     :: ZDCDT
REAL, DIMENSION(KVECNPT,KEQ,KEQ) :: ZWORKVECT
REAL, DIMENSION(KEQ,KEQ) :: ZWORK,ZINV
INTEGER :: IFAIL, JI,JJ
!
!------------------------------------------------------------------------------
!
!*       1.   CALCULATE FIRST DERIVATIVE AND JACOBIAN
!        --------------------------------------------
!
CALL CH_FCN(PTSIMUL,PCONC,ZDCDTVECT,KMI,KVECNPT,KEQ)
!
!%kloos%
!Matrix ZJAC replaced by ZWORK
!
CALL CH_JAC(PTSIMUL,PCONC,ZWORKVECT,KMI,KVECNPT,KEQ)
!
!*       2.   CALCULATE (1-dt/2 J^n)^-1
!        ------------------------------
!
DO JJ=1,KVECNPT

ZWORK(:,:) = - (0.5*PDTACT) * ZWORKVECT(JJ,:,:)
ZDCDT(:) = ZDCDTVECT(JJ,:)
!%
DO JI = 1, KEQ
  ZWORK(JI,JI) = ZWORK(JI,JI) + 1.0
ENDDO
!
!%kloos%
! Inversion of matrix modified
!
!*       3.   COMPUTE LU FACTORIZATION OF MATRIX (1-dt/2) J^n with LAPACK LIBRARY
!        ------------------------------------------------------------------------
!
IFAIL = 1
CALL CH_GAUSS(ZWORK,ZINV,KEQ,IFAIL)
IF (IFAIL.NE.0) THEN
!callabortstop
CALL ABORT
  STOP 'CH_SIS ERROR: matrix cannot be inverted by CH_GAUSS'
ENDIF
!
!*       4.   CALCULATE (1-dt/2 J^n)^-1 f^n
!        -----------------------------------
!
ZDCDT(:)=MATMUL(ZINV,ZDCDT)
!
!*       5.   CALCULATE y^n+1 = y^n + dt*(1-dt/2 J^n)^-1 f^n
!        ---------------------------------------------------
!
PNEWCONC(JJ,:) = PCONC(JJ,:) + PDTACT * ZDCDT(:)

ENDDO
!%
!
END SUBROUTINE CH_SIS
