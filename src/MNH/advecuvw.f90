!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_ADVECUVW 
!     ####################
!
INTERFACE
!
      SUBROUTINE ADVECUVW ( PUT,  PVT,  PWT,                    &
                            PRUCT, PRVCT, PRWCT,                &
                            PRUS,  PRVS,  PRWS                  )
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT, PWT ! Wind at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS ! Sources of Momentum
!
END SUBROUTINE ADVECUVW
!
END INTERFACE
!
END MODULE MODI_ADVECUVW 
!
!
!
!     ###########################################################
      SUBROUTINE ADVECUVW ( PUT,  PVT,  PWT,                    &
                            PRUCT, PRVCT, PRWCT,                &
                            PRUS,  PRVS,  PRWS                  )
!     ###########################################################
!
!!****  *ADVECUVW * - routine to compute the advection terms of momentum
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the three advection terms
!!    of each component of the momentum, written in flux form.
!!      The advection velocity is taken as the contravariant form of 
!!    the momentum for extension to non-cartesian geometry and 
!!    conformal projection cases. The different sources terms are stored for
!!    the budget computations.
!!     
!!
!!**  METHOD
!!    ------
!!      The left and right lateral EXTernal zones, have been previously
!!    prepared in routine LBC_S, to avoid particular cases close to the
!!    Lateral Boundaries in this routine.
!!      The Shuman functions are used to write the mean and finite 
!!    differences operators.
!!
!!    EXTERNAL
!!    --------
!!      DXM,DYM,DZM : Shuman functions (finite differences operators)
!!      BUDGET      : Stores the different budget components
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPVEXT: define the number of marginal points out of the 
!!        physical domain along the vertical direction.
!!    
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RU       : logical for budget of RU (wind component along x)
!!                        .TRUE. = budget of RU         
!!                        .FALSE. = no budget of RU 
!!         LBU_RV       : logical for budget of RV (wind component along y)
!!                        .TRUE. = budget of RV         
!!                        .FALSE. = no budget of RV 
!!         LBU_RW        : logical for budget of RW (wind component along z)
!!                        .TRUE. = budget of RW         
!!                        .FALSE. = no budget of RW 
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine ADVECUVW )
!!
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty      * Laboratoire d'Aerologie*
!!  	J.-P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/07/94 
!!      Corrections 06/09/94 (J.-P. Lafore)
!!                  02/11/94 (J.Stein)   extrapolation under the ground
!!                  16/03/95 (J.Stein)   remove R from the historical variables
!!                  01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                  16/10/95 (J. Stein)     change the budget calls 
!!                  19/12/96 (J.-P. Pinty)  update the budget calls 
!!                  07/11/02 (V. Masson)    update the budget calls 
!! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_BUDGET
USE MODD_GRID_n
!
USE MODI_SHUMAN
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT, PWT ! Wind at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS ! Sources of Momentum
!
INTEGER :: IKU
!
!  
!-------------------------------------------------------------------------------
!
IKU=SIZE(XZHAT)
!
!*       1.     COMPUTES THE ADVECTIVE TENDANCIES
!	        ---------------------------------
!
PRUS(:,:,:) = PRUS(:,:,:)                              &
             -DXM( MXF(PRUCT(:,:,:))*MXF(PUT(:,:,:)) ) 
IF (LBUDGET_U)  CALL BUDGET (PRUS,1,'ADVX_BU_RU')
!
PRUS(:,:,:) = PRUS(:,:,:)                              &
             -DYF( MXM(PRVCT(:,:,:))*MYM(PUT(:,:,:)) ) 
IF (LBUDGET_U)  CALL BUDGET (PRUS,1,'ADVY_BU_RU')
!
PRUS(:,:,:) = PRUS(:,:,:)                              &
             -DZF(1,IKU,1, MXM(PRWCT(:,:,:))*MZM(1,IKU,1,PUT(:,:,:)) )
IF (LBUDGET_U)  CALL BUDGET (PRUS,1,'ADVZ_BU_RU')
!
!
PRVS(:,:,:) = PRVS(:,:,:)                              &
             -DXF( MYM(PRUCT(:,:,:))*MXM(PVT(:,:,:)) ) 
IF (LBUDGET_V)  CALL BUDGET (PRVS,2,'ADVX_BU_RV')
!
PRVS(:,:,:) = PRVS(:,:,:)                              &
             -DYM( MYF(PRVCT(:,:,:))*MYF(PVT(:,:,:)) )  
IF (LBUDGET_V)  CALL BUDGET (PRVS,2,'ADVY_BU_RV')
!
PRVS(:,:,:) = PRVS(:,:,:)                              &
             -DZF(1,IKU,1, MYM(PRWCT(:,:,:))*MZM(1,IKU,1,PVT(:,:,:)) )
IF (LBUDGET_V)  CALL BUDGET (PRVS,2,'ADVZ_BU_RV')
!
!
PRWS(:,:,:) = PRWS(:,:,:)                              &
             -DXF( MZM(1,IKU,1,PRUCT(:,:,:))*MXM(PWT(:,:,:)) ) 
IF (LBUDGET_W)  CALL BUDGET (PRWS,3,'ADVX_BU_RW')
!
PRWS(:,:,:) = PRWS(:,:,:)                              &
             -DYF( MZM(1,IKU,1,PRVCT(:,:,:))*MYM(PWT(:,:,:)) ) 
IF (LBUDGET_W)  CALL BUDGET (PRWS,3,'ADVY_BU_RW')
!
PRWS(:,:,:) = PRWS(:,:,:)                              &
             -DZM(1,IKU,1, MZF(1,IKU,1,PRWCT(:,:,:))*MZF(1,IKU,1,PWT(:,:,:)) )
IF (LBUDGET_W)  CALL BUDGET (PRWS,3,'ADVZ_BU_RW')
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECUVW
