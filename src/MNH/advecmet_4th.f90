!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/05/22 19:00:38
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_ADVECMET_4TH
!     ###############################
!
INTERFACE
!
      SUBROUTINE ADVECMET_4TH (HLBCX,HLBCY, KRR,               &
                     PRUCT, PRVCT, PRWCT,                      &
                     PTHT, PTKET, PRT,                         &
                     PRTHS, PRTKES, PRRS, TPHALO2LIST  ) 
!
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                  INTENT(IN)    :: KRR   ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT ! of momentum
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT, PTKET   ! Vars at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES! Source terms
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
!
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list for diffusion
!
END SUBROUTINE ADVECMET_4TH
!
END INTERFACE
!
END MODULE MODI_ADVECMET_4TH
!
!     ######################################################################
      SUBROUTINE ADVECMET_4TH (HLBCX,HLBCY, KRR,               &
                     PRUCT, PRVCT, PRWCT,                      &
                     PTHT, PTKET,  PRT,                        &
                     PRTHS, PRTKES,  PRRS, TPHALO2LIST         )                          
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
!!         LBU_RU       : logical for budget of RU (wind component along x)
!!
!!         LBU_RU       : logical for budget of RU (wind component along x)
!!                        .TRUE. = budget of RU
!!                        .FALSE. = no budget of RU
!!         LBU_RV       : logical for budget of RV (wind component along y)
!!                        .TRUE. = budget of RV
!!                        .FALSE. = no budget of RV
!!         LBU_RW        : logical for budget of RW (wind component along z)
!!                        .TRUE. = budget of RW
!!                        .FALSE. = no budget of RW
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH
!!                        .FALSE. = no budget of RTH
!!         LBU_RTKE     : logical for budget of RTKE (turbulent kinetic energy)
!!                        .TRUE. = budget of RTKE
!!                        .FALSE. = no budget of RTKE
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV
!!                        .FALSE. = no budget of RRV
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC
!!                        .FALSE. = no budget of RRC
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR
!!                        .FALSE. = no budget of RRR
!!         LBU_RRI      : logical for budget of RRI (ice)
!!                        .TRUE. = budget of RRI
!!                        .FALSE. = no budget of RRI
!!         LBU_RRS      : logical for budget of RRS (snow)
!!                        .TRUE. = budget of RRS
!!                        .FALSE. = no budget of RRS
!!         LBU_RRG      : logical for budget of RRG (graupel)
!!                        .TRUE. = budget of RRG
!!                        .FALSE. = no budget of RRG
!!         LBU_RRH      : logical for budget of RRH (hail)
!!                        .TRUE. = budget of RRH
!!                        .FALSE. = no budget of RRH
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
USE MODD_CONF
USE MODD_GRID_n
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
INTEGER,                  INTENT(IN)    :: KRR   ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT  ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT  !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT  ! of momentum
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT, PTKET   ! Vars at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES! Source terms
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
!
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list for diffusion
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: JRR           ! Loop index for  moist variables
INTEGER:: IIB,IJB        ! Begining useful area  in x,y,z directions
INTEGER:: IIE,IJE        ! End useful area in x,y,z directions
INTEGER :: IKU
!
LOGICAL     :: GTKEALLOC                 ! true if TKE arrays are not zero-sized
!
TYPE(HALO2LIST_ll), POINTER :: TZHALO2LIST
!
INTEGER :: IGRID ! localisation on the model grid
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZMEANX, ZMEANY ! fluxes
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTES THE DOMAIN DIMENSIONS
!               ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
GTKEALLOC = SIZE(PTKET,1) /= 0
IKU=SIZE(XZHAT)
!
!-------------------------------------------------------------------------------
!
!*       2.     CALL THE ADVEC_4TH_ORDER_ALGO ROUTINE FOR EACH FIELD
!               ----------------------------------------------------
!
IGRID = 1
!
!!$IF (NHALO == 1) THEN
   TZHALO2LIST => TPHALO2LIST
   CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PTHT, IGRID, ZMEANX, ZMEANY, &
                             TZHALO2LIST%HALO2 )
!!$ELSE
!!$   CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PTHT, IGRID, ZMEANX, ZMEANY)
!!$ENDIF
!
! Thermodynamical variable
!
PRTHS(:,:,:) = PRTHS(:,:,:)                      &
              -DXF( PRUCT(:,:,:) * ZMEANX(:,:,:) ) 
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVX_BU_RTH')
!
PRTHS(:,:,:) = PRTHS(:,:,:)                      &
              -DYF( PRVCT(:,:,:) * ZMEANY(:,:,:) ) 
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVY_BU_RTH')
!
PRTHS(:,:,:) = PRTHS(:,:,:)                           &
              -DZF(1,IKU,1, PRWCT(:,:,:) * MZM4(PTHT(:,:,:)) )
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'ADVZ_BU_RTH')
!
! Turbulence variables
!
IF ( GTKEALLOC ) THEN
!!$  IF(NHALO == 1) THEN
    TZHALO2LIST => TZHALO2LIST%NEXT
    CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PTKET, IGRID, &
                            ZMEANX, ZMEANY, TPHALO2=TZHALO2LIST%HALO2)
!!$  ELSE
!!$    CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PTKET, IGRID, ZMEANX, ZMEANY)
!!$  ENDIF
!
  PRTKES(:,:,:) = PRTKES(:,:,:) 	            &
                 -DXF( PRUCT(:,:,:) * ZMEANX(:,:,:) ) 
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVX_BU_RTKE')
!
  PRTKES(:,:,:) = PRTKES(:,:,:) 	            &
                 -DYF( PRVCT(:,:,:) * ZMEANY(:,:,:) ) 
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVY_BU_RTKE')
!
   PRTKES(:,:,:) = PRTKES(:,:,:) 	                  &
                 -DZF(1,IKU,1, PRWCT(:,:,:) * MZM4(PTKET(:,:,:)) )
  IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'ADVZ_BU_RTKE')
ENDIF
!
!
! Case with KRR moist variables
!
DO JRR=1, KRR
!!$  IF(NHALO == 1) THEN 
    TZHALO2LIST => TZHALO2LIST%NEXT
    CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PRT(:,:,:,JRR), IGRID, &
                              ZMEANX, ZMEANY,TPHALO2=TZHALO2LIST%HALO2 )
!!$  ELSE
!!$    CALL ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PRT(:,:,:,JRR), IGRID, ZMEANX, ZMEANY)
!!$  ENDIF
!
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)                            &
                     -DXF( PRUCT(:,:,:) * ZMEANX(:,:,:) ) 
  IF (JRR==1 .AND. LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVX_BU_RRV')
  IF (JRR==2 .AND. LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVX_BU_RRC')
  IF (JRR==3 .AND. LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVX_BU_RRR')
  IF (JRR==4 .AND. LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVX_BU_RRI')
  IF (JRR==5 .AND. LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'ADVX_BU_RRS')
  IF (JRR==6 .AND. LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'ADVX_BU_RRG')
  IF (JRR==7 .AND. LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'ADVX_BU_RRH')
!
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)                            &
                     -DYF( PRVCT(:,:,:) * ZMEANY(:,:,:) )
  IF (JRR==1 .AND. LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVY_BU_RRV')
  IF (JRR==2 .AND. LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVY_BU_RRC')
  IF (JRR==3 .AND. LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVY_BU_RRR')
  IF (JRR==4 .AND. LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVY_BU_RRI')
  IF (JRR==5 .AND. LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'ADVY_BU_RRS')
  IF (JRR==6 .AND. LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'ADVY_BU_RRG')
  IF (JRR==7 .AND. LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'ADVY_BU_RRH')
!
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)                            &
                    -DZF(1,IKU,1, PRWCT(:,:,:) * MZM4(PRT(:,:,:,JRR)) )
  IF (JRR==1 .AND. LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6 ,'ADVZ_BU_RRV')
  IF (JRR==2 .AND. LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7 ,'ADVZ_BU_RRC')
  IF (JRR==3 .AND. LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8 ,'ADVZ_BU_RRR')
  IF (JRR==4 .AND. LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9 ,'ADVZ_BU_RRI')
  IF (JRR==5 .AND. LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'ADVZ_BU_RRS')
  IF (JRR==6 .AND. LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'ADVZ_BU_RRG')
  IF (JRR==7 .AND. LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'ADVZ_BU_RRH')
ENDDO
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECMET_4TH


