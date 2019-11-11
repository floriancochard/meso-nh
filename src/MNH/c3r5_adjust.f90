!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 microph 2006/06/06 18:25:10
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_C3R5_ADJUST
!     #######################
!
INTERFACE
!
      SUBROUTINE C3R5_ADJUST( KRR, KMI, HRAD,                                  &
                             HTURBDIM, OSUBG_COND, PTSTEP,                     &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, PPABST,         &
                             PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT,         &
                             PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,         &
                             PCCT, PCIT, PCNUCS, PCCS, PINUCS, PCIS,           &
                             PTHS, PSRCS, PCLDFR                               )
         !
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KMI      ! Model index 
CHARACTER*4,              INTENT(IN)    :: HTURBDIM ! Dimensionality of the
                                                    ! turbulence scheme
CHARACTER*4,              INTENT(IN)    :: HRAD     ! Radiation scheme name
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
REAL,                     INTENT(IN)    :: PTSTEP   ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF! Dry density of the 
                                                   ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t     
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRRT ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRIT    ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRST ! Aggregate  m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRGT ! Graupel    m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRHT ! Hail       m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRRS ! Rain water m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRSS ! Aggregate  m.r. at t+1
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRGS ! Graupel    m.r. at t+1
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRHS ! Hail       m.r. at t+1
!
REAL, DIMENSION(:,:,:),   INTENT(IN)       :: PCCT    ! Cloud water conc. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)       :: PCIT    ! Cloud ice   conc. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PCNUCS  ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PCCS    ! Cloud water conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PINUCS  ! Ice Nucl.   conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PCIS    ! Cloud ice   conc. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
!
END SUBROUTINE C3R5_ADJUST
!
END INTERFACE
!
END MODULE MODI_C3R5_ADJUST
!
!     ##########################################################################
      SUBROUTINE C3R5_ADJUST( KRR, KMI, HRAD,                                  &
                             HTURBDIM, OSUBG_COND, PTSTEP,                     &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, PPABST,         &
                             PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT,         &
                             PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,         &
                             PCCT, PCIT, PCNUCS, PCCS, PINUCS, PCIS,           &
                             PTHS, PSRCS, PCLDFR                               )
!     ##########################################################################
!
!!****  *C3R5_ADJUST* -  compute the fast  microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the fast microphysical sources
!!      through an explict scheme and a saturation ajustement procedure.
!!
!!
!!**  METHOD
!!    ------
!!      Reisin et al.,    1996 for the explicit scheme when ice is present
!!      Langlois, Tellus, 1973 for the implict adjustment for the cloud water
!!      (refer also to book 1 of the documentation).
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!         XP00               ! Reference pressure
!!         XMD,XMV            ! Molar mass of dry air and molar mass of vapor
!!         XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
!!         XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
!!         XCL                ! Cl (liquid)
!!         XTT                ! Triple point temperature
!!         XLVTT              ! Vaporization heat constant
!!         XALPW,XBETAW,XGAMW ! Constants for saturation vapor 
!!                            !  pressure  function 
!!      Module  MODD_CONF 
!!         CCONF
!!      Module MODD_BUDGET:
!!         NBUMOD 
!!         CBUTYPE
!!         NBUPROCCTR 
!!         LBU_RTH    
!!         LBU_RRV  
!!         LBU_RRC  
!!      Module MODD_LES : NCTR_LES,LTURB_LES,NMODNBR_LES
!!                        XNA declaration (cloud fraction as global var)
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 1 and Book2 of documentation ( routine FAST_TERMS )
!!      Langlois, Tellus, 1973
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!   
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/12/94 
!!      Modifications: March 1, 1995 (J.M. Carriere) 
!!                                  Introduction of cloud water with order 1
!!                                  formulation
!!      Modifications: June 8, 1995 ( J.Stein )
!!                                  Cleaning 
!!      Modifications: August 30, 1995 ( J.Stein )
!!                                  add Lambda3 for the subgrid condensation
!!   
!!                     October 16, 1995 (J. Stein)     change the budget calls 
!!                     March   16, 1996 (J. Stein)     store the cloud fraction
!!                     April   03, 1996 (J. Stein)     displace the nebulosity
!!                                      computation in the all and nothing case
!!                     April   15, 1996 (J. Stein)     displace the lambda 3 
!!                         multiplication and change the nebulosity threshold
!!                     September 16, 1996  (J. Stein)  bug in the SG cond for
!!                                                     the M computation
!!                     October 10, 1996 (J. Stein)     reformulate the Subgrid
!!                                                     condensation scheme
!!                     October 8,  1996 (Cuxart,Sanchez) Cloud frac. LES diag (XNA)
!!                     December 6, 1996 (J.-P. Pinty)  correction of Delta_2
!!                     November 5, 1996 (J. Stein) remove Rnp<0 values
!!                     November 13 1996 (V. Masson) add prints in test above
!!                     March 11, 1997 (J.-M. Cohard)  C2R2 option
!!                     April  6, 2001 (J.-P. Pinty)   C3R5 option 
!!
!-------------------------------------------------------------------------------
!
PRINT *,'C3R5_ADJUST IS NOT YET DEVELOPPED'
!callabortstop
CALL ABORT
STOP
!
END SUBROUTINE C3R5_ADJUST
