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
!**FILE:     ch_svode.f90
!**AUTHOR:   Karsten Suhre
!**DATE:     Fri Nov 10 09:17:45 GMT 1995
!**PURPOSE:  interface for solver SVODE
!**ORIGINAL: K. SUHRE 10/11/95 
!**MODIFIED: 18/02/99 (K. Suhre) vectorize SVODE
!
!!    ############################# 
      MODULE MODI_CH_SVODE
!!    ############################# 

!
INTERFACE

SUBROUTINE CH_SVODE(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                    PRTOL, PATOL, KPED )
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    
				! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC 
				! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model index
REAL,    INTENT(IN) :: PRTOL, PATOL
INTEGER, INTENT(IN) :: KPED
END SUBROUTINE CH_SVODE

END INTERFACE

END MODULE MODI_CH_SVODE

!!    ########################################################################
SUBROUTINE CH_SVODE(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                    PRTOL, PATOL, KPED )
!!    ######################################################################## 
!!
!!****  *CH_SVODE*

!!    PURPOSE
!!    -------
!!    solve one time-step of the chemical differential equation d/dt C = f(C)

!!    METHOD
!!    ------
!!    the stiff-solver SVODE (written in FORTRAN77)
!!    will be used to solve the system
!!    NOTE: this subroutine is not 100% doctorized since we want
!!          to conserve quasi-standard variable names used by SVODE

!!    REFERENCE
!!    ---------
!!    MesoNH book 2 

!!    AUTHOR
!!    ------
!!    K. Suhre

!!    MODIFICATIONS
!!    -------------
!!    Original 10/11/95
!!    01/08/01 (C. Mari)  add arguments
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
!!    SVODE: stiff solver from NETLIB
!!    CH_SVODE_FCN, CH_SVODE_JAC: argument list converters for CH_FCN, CH_JAC

!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    
				! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC 
				! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model index
REAL,    INTENT(IN) :: PRTOL, PATOL
INTEGER, INTENT(IN) :: KPED

!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
REAL, DIMENSION(KEQ) :: ZCONC  ! concentration vector
REAL :: ZTBEGIN, ZTEND         ! begin and end of integration step

INTEGER :: ITOL    = 1      ! type of error, don't play with this parameter
INTEGER :: ITASK   = 1      ! integrate from PTSIMUL to PTSIMUL + PDTACT
INTEGER :: ISTATE           ! integer flag (input and output)
INTEGER :: IOPT    = 0      ! 0 to indicate no optional input used
INTEGER :: IMF              ! 21=use stiff full analytical Jacobian
                            ! 22=use stiff full numerical Jacobian
REAL, DIMENSION(KEQ) :: ZRTOL          ! relative tolerance (XRTOL)
REAL, DIMENSION(KEQ) :: ZATOL          ! absolute tolerance (XATOL)

! workspace declaration
INTEGER :: IRW ! will be set to (22 + 9*KEQ + 2*KEQ*KEQ) in the code
INTEGER :: IIW ! will be set to (30 + KEQ) in the code
REAL, DIMENSION(22 + 9*KEQ + 2*KEQ*KEQ) :: ZWORK ! workspace
INTEGER, DIMENSION(30 + KEQ)            :: IWORK ! workspace

! dummy parameters
INTEGER, DIMENSION(1) :: IPAR ! dummy parameter
REAL, DIMENSION(1)    :: ZPAR ! dummy parameter

INTEGER :: JI ! loop counter

EXTERNAL CH_SVODE_FCN, CH_SVODE_JAC

!*    EXECUTABLE STATEMENTS
!     ---------------------

ZTBEGIN = PTSIMUL
ZTEND = PTSIMUL + PDTACT

ZRTOL(:) = PRTOL
ZATOL(:) = PATOL

! set array dimensions
IRW = (22 + 9*KEQ + 2*KEQ*KEQ) 
IIW = (30 + KEQ) 

! choose calculation of Jacobian
IF (KPED .EQ. 0) THEN
  IMF = 22        ! numerical calculation of Jacobian
ELSE
  IMF = 21        ! analytical calculation of Jacobian using CH_JAC
ENDIF

! at each call to SVODE start a new iteration cycle
ISTATE = 1

! this solver is not vectorized, we loop over all elements
DO JI = 1, KVECNPT

  ZCONC(:) = PCONC(JI,:)

  ! call SVODE solver
  CALL SVODE (CH_SVODE_FCN, KEQ, ZCONC, ZTBEGIN, ZTEND, &
	      ITOL, ZRTOL, ZATOL, ITASK, &
              ISTATE, IOPT, ZWORK, IRW, IWORK, IIW, CH_SVODE_JAC, IMF, &
              ZPAR, IPAR, KMI, JI)

  IF (ISTATE.LT.0) THEN
    PRINT *, "Problems !!! ISTATE = ", ISTATE
    PRINT *, "at vector element ", JI, " out of ", KVECNPT
    STOP "CH_SVODE: program stopped due to SVODE error!"
  ENDIF

  PNEWCONC(JI,:) = ZCONC(:)

END DO

END SUBROUTINE CH_SVODE

SUBROUTINE CH_SVODE_FCN(KEQ, PTIME, PCONC, PFCN, PPAR, KPAR, KMI, KINDEX)
USE MODI_CH_FCN
IMPLICIT NONE
!*  routine to change parameter interface for CH_FCN
INTEGER :: KEQ  ! that is NEQ (dimension of the problem)
REAL    :: PTIME ! integration time
REAL, DIMENSION(KEQ) :: PCONC ! concentration
REAL, DIMENSION(KEQ) :: PFCN  ! first derivative of the system
INTEGER :: KPAR  ! dummy parameter
REAL    :: PPAR  ! dummy parameter
INTEGER, INTENT(IN) :: KMI      ! model number
INTEGER, INTENT(IN) :: KINDEX   ! index of CCSTYPE array to treat 
REAL, DIMENSION(1,KEQ) :: ZCONC ! concentration
REAL, DIMENSION(1,KEQ) :: ZFCN  ! first derivative of the system

! call CH_FCN
ZCONC(1,:) = PCONC(:)
CALL CH_FCN(PTIME,ZCONC,ZFCN,KMI,1,KEQ)
PFCN(:) = ZFCN(1,:)
RETURN
END SUBROUTINE CH_SVODE_FCN

SUBROUTINE CH_SVODE_JAC(KEQ, PTIME, PCONC, KL, KU, PJAC, KROWPD, PPAR, KPAR, KMI, KINDEX)
USE MODI_CH_JAC
IMPLICIT NONE
!*  routine to change parameter interface for CH_JAC
INTEGER :: KEQ  ! that is NEQ (dimension of the problem)
REAL    :: PTIME ! integration time
REAL, DIMENSION(KEQ) :: PCONC     ! concentration
REAL, DIMENSION(KEQ,KEQ) :: PJAC ! Jacobian matrix
REAL :: PPAR                      ! not used, dummy variables
INTEGER :: KL, KU, KROWPD, KPAR   ! not used, dummy variables
INTEGER, INTENT(IN) :: KMI      ! model number
INTEGER, INTENT(IN) :: KINDEX   ! index of CCSTYPE array to treat 
REAL, DIMENSION(1,KEQ) :: ZCONC ! concentration
REAL, DIMENSION(1,KEQ,KEQ) :: ZJAC  ! first derivative of the system
! call CH_JAC
ZCONC(1,:) = PCONC(:)
CALL CH_JAC(PTIME,ZCONC,ZJAC,KMI,1,KEQ)
PJAC(:,:) = ZJAC(1,:,:)
RETURN
END SUBROUTINE CH_SVODE_JAC


