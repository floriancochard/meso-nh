!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!!    #####################
      MODULE MODI_CH_SOLVER_n
!!    #####################
!!
!
INTERFACE
SUBROUTINE CH_SOLVER_n(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI)
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
INTEGER, INTENT(IN) :: KMI      ! model number
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
END SUBROUTINE CH_SOLVER_n
END INTERFACE
!!
END MODULE MODI_CH_SOLVER_n
!!    ######################################################################### 
      SUBROUTINE CH_SOLVER_n(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI)
!!    ######################################################################### 
!!
!!*** *CH_SOLVER_n*
!!
!!    PURPOSE
!!    -------
!!      solution of one timestep of the chemical differential equation
!!
!!**  METHOD
!!    ------
!!      Calls the individual solvers and passes the corresponding parameters
!!    on. A fixed integration step is used, however, some solvers internally
!!    use variable time steps.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/04/95
!!    10/11/95 KS: integrate SVODE solver (AFF)
!!    10/04/96 KS: integrate QSSA and EXQSSA solver (AFF)
!!    31/07/96 KS: add TPK to parameterlist (some solver desactivated !@)
!!    01/08/01 (C. Mari)  change routine to _n
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    01/06/07 (P. Tulet)  model number in argument (for AROME)
!!    01/06/07 (JP Pinty & M Leriche) add Rosenbrock solvers
!!
!!    EXTERNAL
!!    --------
!!    calls the different solvers like SIS, LinSSA, NAG, etc 
USE MODI_CH_SIS
USE MODI_CH_LINSSA
USE MODI_CH_CRANCK
!USE MODI_CH_SVODE
USE MODI_CH_QSSA
USE MODI_CH_EXQSSA
!@USE MODI_CH_D02EAF
!@USE MODI_CH_D02EBF
!@USE MODI_CH_D02NBF
USE MODE_RBK90_Integrator
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_SOLVER_n
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
INTEGER, INTENT(IN) :: KMI      ! model number
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
!
!*       0.1  declaration of local variables
!
REAL ::  ZSTART_TIME, ZFINAL_TIME
REAL, DIMENSION(KVECNPT,KEQ) ::  ZCONC_DP
REAL, DIMENSION(KEQ) :: ZRTOL
REAL, DIMENSION(KEQ) :: ZATOL
INTEGER :: ICNTRL_U(20) = 0
!
!------------------------------------------------------------------------------
!
!
!*       1.   CALL THE SOLVERS
!        ---------------------
!
SELECT CASE (CSOLVER)
!
CASE ('SIS')
!
  ! call SIS
  CALL CH_SIS(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT,KMI)
  !
CASE ('LINSSA', 'LinSSA')
!
  ! call LinSSA
  CALL CH_LINSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                 NSSA, NSSAINDEX)
!
CASE ('CRANCK')
!
  ! call Cranck-Nicholson method
  CALL CH_CRANCK(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, & 
                 XALPHA)
!
CASE ('D02EAF')
!
  ! call NAG's stiff-solver D02EAF
  !@CALL CH_D02EAF(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KMI)
!callabortstop
CALL ABORT
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER)'
!
CASE ('D02EBF')
!
  ! call NAG's stiff-solver D02EBF
  !@CALL CH_D02EBF(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KMI)
!callabortstop
CALL ABORT
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER)'
!
CASE ('D02NBF')
!
  ! call NAG's stiff-solver D02NBF
  !@CALL CH_D02NBF(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KMI)
!callabortstop
CALL ABORT
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER)'
!
CASE ('SVODE')
!
  ! call SVODE solver

 ! CALL CH_SVODE(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
 !               XRTOL, XATOL, NPED)
!callabortstop
CALL ABORT
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER) until Masdev47'
!
CASE ('QSSA')
!
  ! call QSSA
  CALL CH_QSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
               XSLOW, XFAST, NQSSAITER)
  !
CASE ('EXQSSA')
!
  ! call EXQSSA
  CALL CH_EXQSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                 XSLOW, XFAST, XDTMAX)
  !
CASE ('ROS1','ROS2','ROS3','ROS4','RODAS3','RODAS4','ROSENBROCK')
! default method ROS1 with ROSENBROCK icntrl_u(3)=0
  ! call ROSENBROCK methods
  !
  ZSTART_TIME = PTSIMUL
  ZFINAL_TIME = PTSIMUL + PDTACT
!  ZRTOL(:) = 1.0d-4
  ZRTOL(:) = 1.0d-3
!  ZATOL(:) = 1.0d-3
  ZATOL(:) = 1.0d-2
  ZCONC_DP(:,:) = PCONC(:,:)
  IF (CSOLVER=='ROS2') ICNTRL_U(3) = 1
  IF (CSOLVER=='ROS3') ICNTRL_U(3) = 2
  IF (CSOLVER=='ROS4') ICNTRL_U(3) = 3
  IF (CSOLVER=='RODAS3') ICNTRL_U(3) = 4
  IF (CSOLVER=='RODAS4') ICNTRL_U(3) = 5
  CALL CH_ROSENBROCK(ZSTART_TIME, ZFINAL_TIME, ZCONC_DP,          &
                     KEQ, KMI, KVECNPT, ZATOL, ZRTOL, ICNTRL_U    )
  PNEWCONC(:,:) = ZCONC_DP(:,:)
!
CASE ('NONE')
!
  ! no integration at all (for debugging purposes)
  PNEWCONC(:,:) = PCONC(:,:)
!
CASE DEFAULT
!callabortstop
CALL ABORT
  STOP 'CH_SOLVER_n ERROR: requested solver not supported (CSOLVER)'
END SELECT
!
END SUBROUTINE CH_SOLVER_n
