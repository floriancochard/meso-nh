!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!    ######################### 
      MODULE MODN_CH_SOLVER_n
!!    #########################
!!
!!*** *MODN_CH_SOLVER$n*
!!
!!    PURPOSE
!!    -------
!     namelist for parameters of the stiff solvers
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!    01/08/01 (C. Mari)  change routine to $n 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_SOLVER_n, ONLY: &
         CSOLVER_n => CSOLVER, &
         NSSA_n => NSSA, &
         NSSAINDEX_n => NSSAINDEX, &
         XRTOL_n => XRTOL, &
         XATOL_n => XATOL, &
         NRELAB_n => NRELAB, &
         NPED_n => NPED, &
         NMAXORD_n => NMAXORD, &
         LPETZLD_n => LPETZLD, &
         CMETHOD_n => CMETHOD, &
         CNORM_n => CNORM, &
         NTRACE_n => NTRACE, &
         XALPHA_n => XALPHA, &
         XSLOW_n => XSLOW, &
         XFAST_n => XFAST, &
         NQSSAITER_n => NQSSAITER, &
         XDTMIN_n => XDTMIN, &
         XDTMAX_n => XDTMAX, &
         XDTFIRST_n => XDTFIRST
!
IMPLICIT NONE
!
CHARACTER*32  :: CSOLVER
INTEGER  :: NSSA
INTEGER, DIMENSION(1000)  :: NSSAINDEX
REAL  :: XRTOL
REAL  :: XATOL
INTEGER  :: NRELAB
INTEGER  :: NPED
INTEGER  :: NMAXORD
LOGICAL  :: LPETZLD
CHARACTER*1  :: CMETHOD
CHARACTER*1  :: CNORM
INTEGER  :: NTRACE
REAL  :: XALPHA
REAL  :: XSLOW
REAL  :: XFAST
INTEGER  :: NQSSAITER
REAL  :: XDTMIN
REAL  :: XDTMAX
REAL  :: XDTFIRST
!
NAMELIST/NAM_CH_SOLVERn/CSOLVER,NSSA,NSSAINDEX,XRTOL,XATOL,NRELAB,NPED,&
                        NMAXORD,LPETZLD,CMETHOD,CNORM,NTRACE,XALPHA,XSLOW,&
                        XFAST,NQSSAITER,XDTMIN,XDTMAX,XDTFIRST
!
CONTAINS
!
SUBROUTINE INIT_NAM_CH_SOLVERn
  CSOLVER = CSOLVER_n
  NSSA = NSSA_n
  NSSAINDEX = NSSAINDEX_n
  XRTOL = XRTOL_n
  XATOL = XATOL_n
  NRELAB = NRELAB_n
  NPED = NPED_n
  NMAXORD = NMAXORD_n
  LPETZLD = LPETZLD_n
  CMETHOD = CMETHOD_n
  CNORM = CNORM_n
  NTRACE = NTRACE_n
  XALPHA = XALPHA_n
  XSLOW = XSLOW_n
  XFAST = XFAST_n
  NQSSAITER = NQSSAITER_n
  XDTMIN = XDTMIN_n
  XDTMAX = XDTMAX_n
  XDTFIRST = XDTFIRST_n
END SUBROUTINE INIT_NAM_CH_SOLVERn

SUBROUTINE UPDATE_NAM_CH_SOLVERn
  CSOLVER_n = CSOLVER
  NSSA_n = NSSA
  NSSAINDEX_n = NSSAINDEX
  XRTOL_n = XRTOL
  XATOL_n = XATOL
  NRELAB_n = NRELAB
  NPED_n = NPED
  NMAXORD_n = NMAXORD
  LPETZLD_n = LPETZLD
  CMETHOD_n = CMETHOD
  CNORM_n = CNORM
  NTRACE_n = NTRACE
  XALPHA_n = XALPHA
  XSLOW_n = XSLOW
  XFAST_n = XFAST
  NQSSAITER_n = NQSSAITER
  XDTMIN_n = XDTMIN
  XDTMAX_n = XDTMAX
  XDTFIRST_n = XDTFIRST
END SUBROUTINE UPDATE_NAM_CH_SOLVERn

END MODULE MODN_CH_SOLVER_n
