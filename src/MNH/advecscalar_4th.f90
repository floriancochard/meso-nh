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
!     ###############################
      MODULE MODI_ADVECSCALAR_4TH
!     ###############################
!
INTERFACE
!
      SUBROUTINE ADVECSCALAR_4TH (HLBCX,HLBCY, KSV, PRUCT, PRVCT, PRWCT,  &
                                  PSVT, PRSVS, TPHALO2LIST                )
!
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                  INTENT(IN)    :: KSV   ! Number of Scalar Variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT ! of momentum
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVT           !  
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS         !
!
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list for diffusion
!
END SUBROUTINE ADVECSCALAR_4TH
!
END INTERFACE
!
END MODULE MODI_ADVECSCALAR_4TH
!
!     ######################################################################
      SUBROUTINE ADVECSCALAR_4TH (HLBCX,HLBCY, KSV, PRUCT, PRVCT, PRWCT,  &
                                  PSVT, PRSVS, TPHALO2LIST                )
!     ######################################################################
!
!!****  *ADVEC_4TH_ORDER * - routine to compute the 4th order centered
!!                           advection tendency of scalar variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to call the ADVEC_4TH_ORDER_ALGO
!!    routine for the horizontal advection and the MZM4 and MZF4 functions for
!!    the vertical advection of each prognostic variable. The code is 
!!    parallelized and works for various boundary conditions.
!!
!!**  METHOD
!!    ------
!!      For each prognostic variable the ADVEC_4TH_ORDER routine calls
!!    the ADVEC_4TH_ORDER_ALGO routine which computes the numerical advection
!!    of any 3D field.
!!      The following variables are passed as argument to ADVEC_4TH_ORDER_ALGO :
!!
!!    -- The variable at t
!!    -- The second layer of the halo of the field at t
!!    -- The horizontal advection fluxes
!!    -- The localisation on the model grid :
!!
!!        IGRID = 1 for mass grid point
!!        IGRID = 2 for U grid point
!!        IGRID = 3 for V grid point
!!        IGRID = 4 for w grid point
!!
!!    EXTERNAL
!!    --------
!!      BUDGET      : Stores the different budget components
!!                    (not used in current version)
!!      ADVEC_4TH_ORDER_ALGO : computes the horizontal advection fluxes
!!      MZF4 and MZM4 : computes the vertical advection fluxes
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    MODULE MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         Switches for budgets activations:
!!
!!         LBU_RSV      : logical for budget of RSVx (scalar variable)
!!                        .TRUE. = budget of RSVx
!!                        .FALSE. = no budget of RSVx
!!
!!    MODULE MODD_ARGSLIST
!!         HALO2LIST_ll : type for a list of "HALO2_lls"
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine ADVEC_4TH_ORDER )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   25/10/05
!!
!! Correction :	
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
USE MODD_GRID_n
USE MODD_CONF
USE MODD_BUDGET
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
USE MODI_SHUMAN
USE MODI_BUDGET
!
! incorporate ADVEC_4TH_ORDER_ALG, MZF4 and MZM4
USE MODI_ADVEC_4TH_ORDER_AUX
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                  INTENT(IN)    :: KSV   ! Number of Scalar Variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT  ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT  !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT  ! of momentum
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PSVT           !    
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS    
!
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list for diffusion
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: JSV           ! Loop index for Scalar Variables
INTEGER:: IIB,IJB        ! Begining useful area  in x,y,z directions
INTEGER:: IIE,IJE        ! End useful area in x,y,z directions
!
TYPE(HALO2LIST_ll), POINTER :: TZHALO2LIST
!
INTEGER :: IGRID ! localisation on the model grid
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZMEANX, ZMEANY ! fluxes
INTEGER :: IKU
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTES THE DOMAIN DIMENSIONS
!               ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKU=SIZE(XZHAT)
!
!-------------------------------------------------------------------------------
!
!*       2.     CALL THE ADVEC_4TH_ORDER_ALGO ROUTINE FOR EACH FIELD
!               ----------------------------------------------------
!
IGRID = 1
!
!
! Case with KSV tracers
!
DO JSV=1,KSV
!
!!$  IF(NHALO == 1) THEN
     TZHALO2LIST => TPHALO2LIST     
     CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PSVT(:,:,:,JSV), IGRID, &
                               ZMEANX, ZMEANY,TZHALO2LIST%HALO2 )
!!$  ELSE
!!$     CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PSVT(:,:,:,JSV), IGRID, ZMEANX, ZMEANY)
!!$  ENDIF
!
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)                            &
                    -DXF( PRUCT(:,:,:) * ZMEANX(:,:,:) )
  IF (LBUDGET_SV) CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVX_BU_RSV')
!
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)                            &
                    -DYF( PRVCT(:,:,:) * ZMEANY(:,:,:) ) 
  IF (LBUDGET_SV) CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVY_BU_RSV')
!
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)                           &
                    -DZF(1,IKU,1,  PRWCT(:,:,:) * MZM4(PSVT(:,:,:,JSV)) )
  IF (LBUDGET_SV) CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVZ_BU_RSV')
ENDDO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECSCALAR_4TH


