!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 microph 2006/08/10 17:06:04
!-----------------------------------------------------------------
!      ######################
       MODULE MODI_ICE_C1R3
!      ######################
!
INTERFACE
      SUBROUTINE ICE_C1R3 (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,            &
                           PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PW_NU,  &
                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,       &
                           PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                           PCCT, PCRT, PCIT, PCNS, PCCS, PCRS, PINS, PCIS, &
                           PINPRS, PINPRG                                  )
!
!
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
						   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                      ! integration for  ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT    ! Rain water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel/hail m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCNS    ! Cloud  C. nuclei C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINS    ! Ice nuclei C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
!
END SUBROUTINE ICE_C1R3
END INTERFACE
END MODULE MODI_ICE_C1R3
!     ######################################################################
      SUBROUTINE ICE_C1R3 (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,            &
                           PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PW_NU,  &
                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,       &
                           PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                           PCCT, PCRT, PCIT, PCNS, PCCS, PCRS, PINS, PCIS, &
                           PINPRS, PINPRG                                  )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources of cloud water and
!!           rain water concentrations and mixing ratios
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources:
!!    nucleation, sedimentation, autoconversion, accretion, self-collection 
!!    and vaporisation which are parameterized according to Cohard and Pinty 
!!    QJRMS, 2000
!!
!!
!!**  METHOD
!!    ------
!!      The activation of CCN is checked for quasi-saturated air parcels 
!!    to update the cloud droplet number concentration. Then assuming a 
!!    generalized gamma distribution law for the cloud droplets and the 
!!    raindrops, the zeroth and third order moments tendencies are evaluated
!!    for all the coalescence terms by integrating the Stochastic Collection 
!!    Equation. As autoconversion is a process that cannot be resolved 
!!    analytically, the Berry-Reinhardt parameterisation is employed with
!!    modifications to initiate the raindrop spectrum mode. The integration
!!    of the raindrop evaporation of the raindrops below clouds is 
!!    straightformward.
!!
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!
!!      Module MODD_CST     
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure
!!                               function over liquid water
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH        
!!                        .FALSE. = no budget of RTH
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV 
!!                        .FALSE. = no budget of RRV 
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC 
!!                        .FALSE. = no budget of RRC 
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR 
!!                        .FALSE. = no budget of RRR 
!!
!!    REFERENCE
!!    ---------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             31/12/96 
!!      Jean-Pierre PINTY     7/ 4/01  Code cleaning
!!      Jean-Pierre PINTY     7/ 5/01  Bug correction in BERFI
!!      Jean-Pierre PINTY    17/ 5/01  Reset PINS=0 in case of IMLT
!!      Jean-Pierre PINTY    17/ 5/01  Bug in RRCFRIG and RICFRRG
!!      Jean-Pierre PINTY    29/ 5/01  Bug in RCHONI and graupel shedding
!!      Jean-Pierre PINTY    29/ 6/01  Bug in RCHONI and RVHNCI
!!      Jean-Pierre PINTY    29/ 6/01  Add RHHONI process (freezing haze part.)
!!      Jean-Pierre PINTY    13/ 9/01  Recode the RCHONI and RVHNCI processes
!!      Jean-Pierre PINTY    23/ 9/01  Recode the HM processes according to
!!                                     Beheng(1986) and Ovtchin. et al. (2000)
!!                                     and add the S to I conversion rate
!!      Jean-Pierre PINTY     1/10/01  Bug in RVHNCI
!!      Jean-Pierre PINTY     8/10/01  Revise limits in sedim. and review S->I
!!      Jean-Pierre PINTY    18/10/01  Revise Snow to Ice conversion
!!      Jean-Pierre PINTY    18/12/01  Revise Graupel wet growth (limitation)
!!
!-------------------------------------------------------------------------------
!
PRINT *,'ICE_C1R3 IS NOT YET DEVELOPPED'
!callabortstop
CALL ABORT
STOP
!
END SUBROUTINE ICE_C1R3
