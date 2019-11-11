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
!     #######################
      MODULE MODI_ADVECSCALAR 
!     #######################
INTERFACE
      SUBROUTINE ADVECSCALAR ( KSV, PSVT, PRUCT, PRVCT, PRWCT, PRSVS  )                     
!
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT 
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum 
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                  ! Sources terms 
END SUBROUTINE ADVECSCALAR
!
END INTERFACE
!
END MODULE MODI_ADVECSCALAR 
!
!
!
!     ######################################################################
      SUBROUTINE ADVECSCALAR ( KSV, PSVT, PRUCT, PRVCT, PRWCT, PRSVS  )                     
!     ######################################################################
!
!!****  *ADVECSCALAR * - routine to compute the advection tendancies of the
!!                       tracer scalar fields.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the total advection 
!!    tendancies of all the scalar fields, written in flux form.
!!      The advection velocity is taken as the contravariant form of 
!!    the momentum for extension to non-cartesian geometry and 
!!    conformal projection cases.
!!
!!
!!**  METHOD
!!    ------
!!      The left and right lateral EXTernal zones, have been previously
!!    prepared in routine LBC_S, to avoid particular cases close to the
!!    Lateral Boundaries in this routine.
!!      The Shuman functions are used to write the mean and finite 
!!    differences operators.
!!      The different sources terms are stored for the budget
!!    computations.
!!
!!    EXTERNAL
!!    --------
!!      MXM,MYM,MZM : Shuman functions (mean operators)
!!      DXM,DYM,DZM : Shuman functions (finite differences operators)
!!      BUDGET      : Stores the different budget components
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RSV      : logical for budget of RSVx (scalar variables)
!!                        .TRUE. = budget of RSV
!!                        .FALSE. = no budget of RSV
!!        
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine ADVECSCALAR )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!  	J.-P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/07/94 
!!      Corrections 06/09/94 (J.-P. Lafore)
!!                  16/03/95 (J. Stein)     remove R from the historical var.
!!                  01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                  16/10/95 (J. Stein)     change the budget calls 
!!                  19/12/96 (J.-P. Pinty)  update the budget calls 
!!                  07/11/02 (V. Masson)    update the budget calls
!!                  24/04/06 (C.Lac)        Split scalar and passive
!!                                          tracer routines
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
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
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum 
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                  ! Sources terms 
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JSV           ! Loop index for Scalar Variables
INTEGER :: IKU
!
!  
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTES THE ADVECTIVE TENDENCIES
!     	        ---------------------------------
!
IKU=SIZE(XZHAT)
! 
                                        ! Case with KSV Scalar Variables
DO JSV=1,KSV
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)                            &
                    -DXF( PRUCT(:,:,:) * MXM (PSVT(:,:,:,JSV)) ) 
END DO
IF (LBUDGET_SV) THEN
  DO JSV=1,KSV
    CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVX_BU_RSV')
  END DO
END IF
!
DO JSV=1,KSV
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)                            &
                    -DYF( PRVCT(:,:,:) * MYM (PSVT(:,:,:,JSV)) ) 
END DO
IF (LBUDGET_SV) THEN
  DO JSV=1,KSV
    CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVY_BU_RSV')
  END DO
END IF
!
DO JSV=1,KSV
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)                            &
                    -DZF(1,IKU,1, PRWCT(:,:,:) * MZM (1,IKU,1, PSVT(:,:,:,JSV)) ) 
END DO
IF (LBUDGET_SV) THEN
  DO JSV=1,KSV
    CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVZ_BU_RSV')
  END DO
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECSCALAR
