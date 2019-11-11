!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_NUM_DIFF
!     ####################
!
INTERFACE
!
      SUBROUTINE NUM_DIFF ( HLBCX, HLBCY, KRR, KSV, PDK2U, PDK4U,           &
                     PDK2TH, PDK4TH, PDK2SV, PDK4SV, KMI,                   &
                     PUM, PVM, PWM, PTHM, PTKEM, PRM, PSVM,                 &
                     PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PRHODJ,                &
                     PRUS, PRVS, PRWS, PRTHS, PRTKES, PRRS, PRSVS,          &
                     OZDIFFU,ONUMDIFU, ONUMDIFTH, ONUMDIFSV,                &
                     TPHALO2LIST, TPHALO2LSLIST,PZDIFFU_HALO2               )
!
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODE_TYPE_ZDIFFU
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                  INTENT(IN)    :: KRR   ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV   ! Number of Scalar Variables
INTEGER,                  INTENT(IN)    :: KMI   ! Model index
!
REAL,                     INTENT(IN)    :: PDK2U  ! 2nd order dif. coef. /dx2
REAL,                     INTENT(IN)    :: PDK4U  ! 4th order dif. coef. /dx4
                                                  ! for momentum
REAL,                     INTENT(IN)    :: PDK2TH ! for theta and mixing ratios
REAL,                     INTENT(IN)    :: PDK4TH ! and TKE
REAL,                     INTENT(IN)    :: PDK2SV ! for theta and mixing ratios
REAL,                     INTENT(IN)    :: PDK4SV ! and TKE

!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSUM,PLSVM      ! Larger
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSWM,PLSTHM     ! Scales
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSRVM           ! Variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ            ! dry density of
!                                 ! anelastic reference state * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUM, PVM, PWM        ! Variables
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHM, PTKEM          !    at
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRM , PSVM           !   t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS     ! Source terms
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS , PRSVS         !
!
LOGICAL, INTENT(IN) :: OZDIFFU
LOGICAL, INTENT(IN) :: ONUMDIFU
LOGICAL, INTENT(IN) :: ONUMDIFTH
LOGICAL, INTENT(IN) :: ONUMDIFSV
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST, TPHALO2LSLIST ! list for diffusion
TYPE(TYPE_ZDIFFU_HALO2) :: PZDIFFU_HALO2
!
END SUBROUTINE NUM_DIFF
!
END INTERFACE
!
END MODULE MODI_NUM_DIFF
!
!
!
!     ######################################################################
      SUBROUTINE NUM_DIFF ( HLBCX, HLBCY, KRR, KSV, PDK2U, PDK4U,           &
                     PDK2TH, PDK4TH, PDK2SV, PDK4SV, KMI,                   &
                     PUM, PVM, PWM, PTHM, PTKEM, PRM, PSVM,                 &
                     PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PRHODJ,                &
                     PRUS, PRVS, PRWS, PRTHS, PRTKES, PRRS, PRSVS,          &
                     OZDIFFU,ONUMDIFU, ONUMDIFTH, ONUMDIFSV,                &
                     TPHALO2LIST, TPHALO2LSLIST,PZDIFFU_HALO2               )
!     ######################################################################
!
!!****  *NUM_DIFF * - routine to call the NUM_DIFF_ALGO routine
!!                    for each prognostic variable
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to call the NUM_DIFF_ALGO
!!    routine for each prognostic variable. The NUM_DIFF_ALGO routine
!!    applies an horizontal diffusion to the prognostic variables.
!!
!!
!!**  METHOD
!!    ------
!!    For each prognostic variable the NUM_DIFF routine calls
!!    the NUM_DIFF_ALGO routine which computes the numerical diffusion
!!    for a 3D field.
!!      The following variables are passed as argument to NUM_DIFF_ALGO :
!!
!!    -- The source term and the variable at t-dt
!!
!!    -- The density PRHODJ
!!          or the mean operator in x or y or z direction applied to PRHODJ
!!          (as appropriate)
!!
!!    -- The second layer of the halo of the field at t-dt
!!
!!    -- The second layer of the halo of the LS field at t-dt, if necessary
!!
!!    -- The LS field at t-dt, if necessary
!!
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
!!
!!      NUM_DIFF_ALGO : applies an horizontal diffusion to the prognostic
!!                      variables
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
!!         LBU_RSV      : logical for budget of RSVx (scalar variable)
!!                        .TRUE. = budget of RSVx
!!                        .FALSE. = no budget of RSVx
!!
!!    MODULE MODD_PARALLEL2
!!
!!         HALO2_ll     : type for a set of four arrays representing
!!                        the second layer of the halo
!!
!!    MODULE MODD_ARGSLIST
!!
!!         HALO2LIST_ll : type for a list of "HALO2_lls"
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine NUM_DIFF )
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   15/10/94
!!      J. Stein   08/12/94 change the variable to be diffused rho*J*phi => phi
!!      J. Stein   20/03/95 remove R from the historical variables + switch 2D
!!                          + switch for TKE unused
!!                 01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                 16/10/95 (J. Stein)     change the budget calls
!!                 19/12/96 (J.-P. Pinty)  update the budget calls
!!                 09/06/98 (P. Kloos)     restructure to run in parallel
!!                 06/11/02 (V. Masson)    update the budget calls
!!                 03/11/04 (G. Zängl)     add calls for truly horizontal diffusion
!!                 05/06    (C.Lac)        Remove EPS
!!                 05/07    (C.Lac)        Separation between variables
!!                 07/09    (C.Lac)        Correction on budget calls
!!     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!     J.Escobar : 05/12/2017 : Pb SegFault , correct IF(ONUMDIFTH/OZDIFFU) nesting
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
USE MODD_BUDGET
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
USE MODI_SHUMAN
USE MODI_BUDGET
!
USE MODE_TYPE_ZDIFFU
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                  INTENT(IN)    :: KRR   ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV   ! Number of Scalar Variables
INTEGER,                  INTENT(IN)    :: KMI   ! Model index
!
REAL,                     INTENT(IN)    :: PDK2U  ! 2nd order dif. coef. /dx2
REAL,                     INTENT(IN)    :: PDK4U  ! 4th order dif. coef. /dx4
                                                  ! for momentum
REAL,                     INTENT(IN)    :: PDK2TH ! for theta and mixing ratios
REAL,                     INTENT(IN)    :: PDK4TH ! 
REAL,                     INTENT(IN)    :: PDK2SV ! for theta and mixing ratios
REAL,                     INTENT(IN)    :: PDK4SV ! 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSUM,PLSVM      ! Larger
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSWM,PLSTHM     ! Scales
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSRVM           ! Variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ            ! dry density of
!                                 ! anelastic reference state * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUM, PVM, PWM        ! Variables
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHM, PTKEM          !    at
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRM , PSVM           !   t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS     ! Source terms
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES        !
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS , PRSVS         !
!
! Fields for z-diffusion
LOGICAL, INTENT(IN) :: OZDIFFU
LOGICAL, INTENT(IN) :: ONUMDIFU
LOGICAL, INTENT(IN) :: ONUMDIFTH
LOGICAL, INTENT(IN) :: ONUMDIFSV
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST, TPHALO2LSLIST ! list for diffusion
TYPE(TYPE_ZDIFFU_HALO2) :: PZDIFFU_HALO2
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR           ! Loop index for  moist variables
INTEGER :: JSV           ! Loop index for Scalar Variables
INTEGER:: IIB,IJB        ! Begining useful area  in x,y directions
INTEGER:: IIE,IJE        ! End useful area in x,y directions
INTEGER :: IKU
!
LOGICAL     :: GTKEALLOC                 ! true if TKE arrays are not zero-sized
!
TYPE(HALO2LIST_ll), POINTER :: TZHALO2LIST, TZHALO2LSLIST
!
INTEGER :: IGRID ! localisation on the model grid
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTES THE DOMAIN DIMENSIONS
!               ------------------------------
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKU=SIZE(PUM,3)
!
GTKEALLOC = SIZE(PTKEM,1) /= 0
!
!-------------------------------------------------------------------------------
!
!*       2.     CALL THE NUM_DIFF_ALGO ROUTINE FOR EACH FIELD
!               ---------------------------------------------
!
IF (ONUMDIFU) THEN
 IGRID = 2
!!$ IF(NHALO == 1) THEN
  TZHALO2LIST => TPHALO2LIST
  TZHALO2LSLIST => TPHALO2LSLIST
  CALL NUM_DIFF_ALGO(PRUS, PUM, IGRID, MXM(PRHODJ), PDK2U, PDK4U, &
                     PLSUM,TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
!!$ ELSE
!!$  CALL NUM_DIFF_ALGO(PRUS, PUM, IGRID, MXM(PRHODJ), PDK2U, PDK4U, PLSUM )
!!$ ENDIF
!
 IGRID = 3
!!$ IF(NHALO == 1) THEN
  TZHALO2LIST => TZHALO2LIST%NEXT
  TZHALO2LSLIST => TZHALO2LSLIST%NEXT
  CALL NUM_DIFF_ALGO(PRVS, PVM, IGRID, MYM(PRHODJ), PDK2U, PDK4U, &
                     PLSVM, TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
!!$ ELSE
!!$  CALL NUM_DIFF_ALGO(PRVS, PVM, IGRID, MYM(PRHODJ), PDK2U, PDK4U, PLSVM )
!!$ ENDIF
!
 IGRID = 4
!
!!$ IF(NHALO == 1) THEN
  TZHALO2LIST => TZHALO2LIST%NEXT
  TZHALO2LSLIST => TZHALO2LSLIST%NEXT
  CALL NUM_DIFF_ALGO(PRWS, PWM, IGRID, MZM(1,IKU,1,PRHODJ), PDK2U, PDK4U, &
                     PLSWM, TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
!!$ ELSE
!!$  CALL NUM_DIFF_ALGO(PRWS, PWM, IGRID, MZM(1,IKU,1,PRHODJ), PDK2U, PDK4U, PLSWM )
!!$ ENDIF
ENDIF
!
IF (ONUMDIFTH) THEN
 IGRID = 1
!
!!$ IF(NHALO == 1) THEN
  TZHALO2LIST => TZHALO2LIST%NEXT
  TZHALO2LSLIST => TZHALO2LSLIST%NEXT
  IF (OZDIFFU) THEN   ! call z-diffusion for potential temperature
   CALL NUM_DIFF_ALGO_Z(PRTHS, PTHM, IGRID, PRHODJ,            &
                         PDK2TH, PDK4TH, PLSTHM,               &
                         TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
  ELSE
   CALL NUM_DIFF_ALGO(PRTHS, PTHM, IGRID, PRHODJ,              &
                      PDK2TH, PDK4TH, PLSTHM,                  &
                     TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
  ENDIF
!!$ ELSE
!!$  IF (OZDIFFU) THEN   ! call z-diffusion for potential temperature
!!$    CALL NUM_DIFF_ALGO_Z(PRTHS, PTHM, IGRID, PRHODJ, &
!!$                         PDK2TH, PDK4TH, PLSTHM )
!!$  ELSE
!!$    CALL NUM_DIFF_ALGO(PRTHS, PTHM, IGRID, PRHODJ, &
!!$                         PDK2TH, PDK4TH, PLSTHM )
!!$  ENDIF
!!$ ENDIF
!
 IF ( GTKEALLOC ) THEN
!!$  IF(NHALO == 1) THEN
    TZHALO2LIST => TZHALO2LIST%NEXT
    CALL NUM_DIFF_ALGO(PRTKES, PTKEM, IGRID, PRHODJ, &
                       PDK2TH, PDK4TH, TPHALO2=TZHALO2LIST%HALO2)
!!$  ELSE
!!$    CALL NUM_DIFF_ALGO(PRTKES, PTKEM, IGRID, PRHODJ, PDK2TH, PDK4TH)
!!$  ENDIF
 ENDIF
!

! Case with KRR moist variables
!
 IF(KRR >= 1) THEN
!!$  IF(NHALO == 1) THEN
    TZHALO2LIST => TZHALO2LIST%NEXT
    TZHALO2LSLIST => TZHALO2LSLIST%NEXT
    IF (OZDIFFU) THEN   ! call z-diffusion for wv mixing ratio
      CALL NUM_DIFF_ALGO_Z(PRRS(:,:,:,1), PRM(:,:,:,1), IGRID, PRHODJ, &
                           PDK2TH, PDK4TH,                             &
                           PLSRVM, TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
    ELSE
      CALL NUM_DIFF_ALGO(PRRS(:,:,:,1), PRM(:,:,:,1), IGRID, PRHODJ, & 
                         PDK2TH, PDK4TH, PLSRVM,                     &
                         TZHALO2LIST%HALO2, TZHALO2LSLIST%HALO2)
    ENDIF
!!$  ELSE
!!$    IF (OZDIFFU) THEN   ! call z-diffusion for wv mixing ratio
!!$      CALL NUM_DIFF_ALGO_Z(PRRS(:,:,:,1), PRM(:,:,:,1), IGRID, PRHODJ, &
!!$                           PDK2TH, PDK4TH, PLSRVM)
!!$    ELSE
!!$     CALL NUM_DIFF_ALGO(PRRS(:,:,:,1), PRM(:,:,:,1), IGRID, PRHODJ, &
!!$                        PDK2TH, PDK4TH, PLSRVM)
!!$    ENDIF
!!$  ENDIF
 ENDIF
!
! In some situations, it makes sense to use the z-diffusion also for cloud water, but it is
! not as necessary as for the water vapor mixing ratio
! This might be added later (is CLW stored in JRR = 2?)
!
 DO JRR=2, KRR
!!$  IF(NHALO == 1) THEN 
    TZHALO2LIST => TZHALO2LIST%NEXT
    CALL NUM_DIFF_ALGO(PRRS(:,:,:,JRR), PRM(:,:,:,JRR), IGRID, PRHODJ, &
                       PDK2TH, PDK4TH, TPHALO2=TZHALO2LIST%HALO2)
!!$  ELSE
!!$    CALL NUM_DIFF_ALGO(PRRS(:,:,:,JRR), PRM(:,:,:,JRR), IGRID, PRHODJ, &
!!$                       PDK2TH, PDK4TH )
!!$  ENDIF
 ENDDO
!
ENDIF
! Case with KSV tracers
!
IF (ONUMDIFSV) THEN
 DO JSV=1,KSV
!!$  IF(NHALO == 1) THEN
    TZHALO2LIST => TZHALO2LIST%NEXT
    CALL NUM_DIFF_ALGO(PRSVS(:,:,:,JSV), PSVM(:,:,:,JSV), IGRID, PRHODJ,&
                       PDK2SV, PDK4SV, TPHALO2=TZHALO2LIST%HALO2)
!!$  ELSE
!!$    CALL NUM_DIFF_ALGO(PRSVS(:,:,:,JSV), PSVM(:,:,:,JSV), IGRID, PRHODJ, &
!!$                       PDK2SV, PDK4SV )
!!$  ENDIF
 ENDDO
END IF
!
!*       3.     STORES FIELDS IN BUDGET ARRAYS
!           ------------------------------
!
IF (LBUDGET_U .AND. ONUMDIFU )   CALL BUDGET (PRUS,1,'DIF_BU_RU')
IF (LBUDGET_V .AND. ONUMDIFU )   CALL BUDGET (PRVS,2,'DIF_BU_RV')
IF (LBUDGET_W .AND. ONUMDIFU )   CALL BUDGET (PRWS,3,'DIF_BU_RW')
IF (LBUDGET_TH .AND. ONUMDIFTH )  CALL BUDGET (PRTHS,4,'DIF_BU_RTH')
IF (LBUDGET_TKE .AND. ONUMDIFTH ) CALL BUDGET (PRTKES,5,'DIF_BU_RTKE')
IF (LBUDGET_RV .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,1),6,'DIF_BU_RRV')
IF (LBUDGET_RC .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,2),7,'DIF_BU_RRC')
IF (LBUDGET_RR .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,3),8,'DIF_BU_RRR')
IF (LBUDGET_RI .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,4),9,'DIF_BU_RRI')
IF (LBUDGET_RS .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,5),10,'DIF_BU_RRS')
IF (LBUDGET_RG .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,6),11,'DIF_BU_RRG')
IF (LBUDGET_RH .AND. ONUMDIFTH )  CALL BUDGET (PRRS(:,:,:,7),12,'DIF_BU_RRH')
IF (LBUDGET_SV .AND. ONUMDIFSV ) THEN
  DO JSV=1,KSV
    CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'DIF_BU_RSV')
  END DO
END IF
!-------------------------------------------------------------------------------
!
!
CONTAINS
!
!     ############################################################
      SUBROUTINE NUM_DIFF_ALGO_Z(PRFIELDS, PFIELDM, KGRID, PRHODJ, &
                                 PDK2, PDK4, PLSFIELD, TPHALO2, TPHALO2LS )
!     ############################################################
!!
!!****  *NUM_DIFF_ALGO * - routine to apply truly horizontal diffusion to the
!!                    prognostic variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to add a 2sd or 4th order horizontal
!!    diffusion term to an advected variable.
!!
!!**  METHOD
!!    ------
!!      In case of cyclic LBCs, the routine adds a numerical diffusion
!!    term obtained by applying an horizontal nabla4 operator to the
!!    prognostic variable on each grid level. In case of open LBCs,
!!    the nabla4 operator degenerates to a nabla2 operator on the first
!!    ring inside the computationnal domain.
!!      In the current version no budget computation is performed.
!!      The "halo2" (or the second layer of the halo) of the prognostic
!!    variable is passed as argument. The "halo2" of the LS field
!!    is also passed as argument, if necessary.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    Module MODD_PARALLEL2
!!
!!      HALO2_ll : type for a set of arrays representing 
!!                 the second layer of the halo
!!
!!    REFERENCE
!!    ---------
!!
!!      Zängl, G., 2002: An improved method for computing horizontal diffusion in a
!!                       sigma-coordinate model and its application to simulations
!!                       over mountainous topography. Mon. Wea. Rev. 130, 1423-1432.
!!
!!    AUTHOR
!!    ------
!!
!!      G. Zängl       * University of Munich*
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!

USE MODE_ll
!
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
USE MODI_NABLA4
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRFIELDS ! source term
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFIELDM  ! variable at t-dt
INTEGER,                INTENT(IN)    :: KGRID    ! C grid localisation
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ   ! RHODJ at KGRID
REAL,                   INTENT(IN)    :: PDK2  ! 2nd order dif. coef. /dx2
REAL,                   INTENT(IN)    :: PDK4  ! 4th order dif. coef. /dx4
!
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PLSFIELD ! Larger Scale variable
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t-dt
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2LS    ! halo2 for the LS field

!!  *** Additional fields for z-diffusion ***
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IW,IE,IS,IN   ! Coordinate of forth order diffusion area

!!  *** Additional fields for z-diffusion ***

REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZPTBFIELD   ! local field to store the perturbation
                                                      ! field; contains PRFIELDS of no LS field is present
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZZMASS      ! local field to store the height at mass points
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZZDIFFU     ! field for truly horizontal diffusion
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZMDIFFU     ! field for diffusion on model levels
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZWGTFAC     ! weighting factor between ZZDIFFU and ZMDIFFU


! Indices for vertical interpolation
INTEGER:: IKIP1,IKIP2,IKIM1,IKIM2,IKJP1,IKJP2,IKJM1,IKJM2, &
          IKIP1P1,IKIP2P1,IKIM1P1,IKIM2P1,IKJP1P1,IKJP2P1,IKJM1P1,IKJM2P1

INTEGER:: IWZ,IEZ,ISZ,INZ   ! Field size of diffusion area including halos
INTEGER:: JI,JJ,JK          ! Do loop indices
INTEGER:: IKP1,IKM1         ! Auxiliary variables for do loops
INTEGER:: IIU,IJU,IKU,IKB,IKE    ! size parameters
INTEGER:: IERROR                 ! DUMMY VARIABLE FOR ERROR MESSAGES

! Real variables to store weighting coefficients
REAL   :: ZWFIP1,ZWFIP2,ZWFIM1,ZWFIM2,ZWFJP1,ZWFJP2,ZWFJM1,ZWFJM2

!
!-------------------------------------------------------------------------------
!

IIU = SIZE(PRFIELDS,1)
IJU = SIZE(PRFIELDS,2)
IKU = SIZE(PRFIELDS,3)
IKE = IKU - JPVEXT
IKB = 1 + JPVEXT

ALLOCATE (ZPTBFIELD(IIB-2:IIE+2,IJB-2:IJE+2,IKU),ZZMASS(IIU,IJU,IKU),ZZDIFFU(IIB-2:IIE+2,IJB-2:IJE+2,IKU),&
          ZMDIFFU(IIB-2:IIE+2,IJB-2:IJE+2,IKB:PZDIFFU_HALO2%NZDLB-1),&
          ZWGTFAC(IIB-2:IIE+2,IJB-2:IJE+2,IKB:PZDIFFU_HALO2%NZDLB-1))

ZPTBFIELD = 0
ZZMASS = 0
ZZDIFFU = 0
ZMDIFFU = 0
ZWGTFAC = 0



!*       1.     CALCULATE THE NUMERICAL DIFFUSION TERM IN THE X DIRECTION
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
!!$  IF(NHALO == 1) THEN
    IW=IIB+1
    IE=IIE-1
!!$  ELSE
!!$    IW=IIB
!!$    IE=IIE
!!$  END IF  
!


IF (PRESENT(PLSFIELD)) THEN
  ZPTBFIELD(IW-2:IE+2,IJB-1:IJE+1,:) = PFIELDM(IW-2:IE+2,IJB-1:IJE+1,:) - PLSFIELD(IW-2:IE+2,IJB-1:IJE+1,:)
!!$  IF(NHALO == 1) THEN
    ZPTBFIELD(IW-3,IJB-1:IJE+1,:) = TPHALO2%WEST(IJB-1:IJE+1,:) - TPHALO2LS%WEST(IJB-1:IJE+1,:)
    ZPTBFIELD(IE+3,IJB-1:IJE+1,:) = TPHALO2%EAST(IJB-1:IJE+1,:) - TPHALO2LS%EAST(IJB-1:IJE+1,:)
!!$  ENDIF
ELSE
  ZPTBFIELD(IW-2:IE+2,IJB-1:IJE+1,:) = PFIELDM(IW-2:IE+2,IJB-1:IJE+1,:)
!!$  IF(NHALO == 1) THEN
    ZPTBFIELD(IW-3,IJB-1:IJE+1,:) = TPHALO2%WEST(IJB-1:IJE+1,:)
    ZPTBFIELD(IE+3,IJB-1:IJE+1,:) = TPHALO2%EAST(IJB-1:IJE+1,:)
!!$  ENDIF
ENDIF

!
!
!*       1.2    NON CYCLIC CASE IN THE X DIRECTION 
!
CASE ('OPEN','WALL','NEST')
!
  IF (LWEST_ll()) THEN
    IF(KGRID == 2) THEN
      IW=IIB+2          ! special case of C grid
    ELSE
      IW=IIB+1
    END IF
  ELSE
!!$    IF(NHALO == 1) THEN
      IW=IIB+1
!!$    ELSE
!!$      IW=IIB
!!$    ENDIF
  ENDIF
!!$  IF (LEAST_ll() .OR. NHALO == 1) THEN
    IE=IIE-1
!!$  ELSE
!!$    IE=IIE
!!$  END IF  


IF (PRESENT(PLSFIELD)) THEN
  ZPTBFIELD(IW-2:IE+2,IJB-1:IJE+1,:) = PFIELDM(IW-2:IE+2,IJB-1:IJE+1,:) - PLSFIELD(IW-2:IE+2,IJB-1:IJE+1,:)
!!$  IF((NHALO == 1).AND.(.NOT.LWEST_ll())) THEN
  IF(.NOT.LWEST_ll()) THEN
    ZPTBFIELD(IW-3,IJB-1:IJE+1,:) = TPHALO2%WEST(IJB-1:IJE+1,:) - TPHALO2LS%WEST(IJB-1:IJE+1,:)
  ENDIF
!!$  IF((NHALO == 1).AND.(.NOT.LEAST_ll())) THEN
  IF(.NOT.LEAST_ll()) THEN
    ZPTBFIELD(IE+3,IJB-1:IJE+1,:) = TPHALO2%EAST(IJB-1:IJE+1,:) - TPHALO2LS%EAST(IJB-1:IJE+1,:)
  ENDIF
ELSE
  ZPTBFIELD(IW-2:IE+2,IJB-1:IJE+1,:) = PFIELDM(IW-2:IE+2,IJB-1:IJE+1,:)
!!$  IF((NHALO == 1).AND.(.NOT.LWEST_ll())) THEN
  IF(.NOT.LWEST_ll()) THEN
    ZPTBFIELD(IW-3,IJB-1:IJE+1,:) = TPHALO2%WEST(IJB-1:IJE+1,:)
  ENDIF
!!$  IF((NHALO == 1).AND.(.NOT.LEAST_ll())) THEN
  IF(.NOT.LEAST_ll()) THEN
    ZPTBFIELD(IE+3,IJB-1:IJE+1,:) = TPHALO2%EAST(IJB-1:IJE+1,:)
  ENDIF
ENDIF


    IF (LWEST_ll()) THEN
!
      PRFIELDS(IW-1,IJB:IJE,:) = PRFIELDS(IW-1,IJB:IJE,:) + PRHODJ(IW-1,IJB:IJE,:) *       &
        PDK2*(ZPTBFIELD(IW-2,IJB:IJE,:)-2.*ZPTBFIELD(IW-1,IJB:IJE,:)+ ZPTBFIELD(IW,IJB:IJE,:) )

    ENDIF
!
    IF (LEAST_ll()) THEN
!
      PRFIELDS(IE+1,IJB:IJE,:) = PRFIELDS(IE+1,IJB:IJE,:) + PRHODJ(IE+1,IJB:IJE,:) *       &
        PDK2*(ZPTBFIELD(IE+2,IJB:IJE,:)-2.*ZPTBFIELD(IE+1,IJB:IJE,:)+ ZPTBFIELD(IE,IJB:IJE,:) )

    ENDIF

!
END SELECT


! Compute diffusion kernel in a truly horizontally way (x direction)

! a) Determine E/W boundaries

!!$IF ((NHALO == 1).AND.(HLBCX(1) == 'CYCL').OR.((.NOT.LWEST_ll()).AND.(NHALO == 1)) ) THEN
IF ( (HLBCX(1) == 'CYCL') .OR. (.NOT.LWEST_ll()) ) THEN
  IWZ = IW-1
ELSE
  IWZ = IW
ENDIF

!!$IF ((NHALO == 1).AND.(HLBCX(1) == 'CYCL').OR.((.NOT.LEAST_ll()).AND.(NHALO == 1)) ) THEN
IF ((HLBCX(1) == 'CYCL').OR.(.NOT.LEAST_ll()) ) THEN
  IEZ = IE+1
ELSE
  IEZ = IE
ENDIF


!!! ATTENTION: HALOS MIGHT BE NEEDED FOR THE INTERPOLATION COEFFICIENTS IN THE CASE
!!! OF IWZ = IW-1 / IEZ = IE+1 (??)

! b) Compute truly horizontal diffusion


DO JI = IWZ,IEZ
  DO JJ = IJB-1,IJE+1   ! are these boundaries right or is ijb,ije sufficient?
    DO JK = IKB,IKE
! Calculate vertical indices and weighting factors for remote points
      IKIP2 = INT(PZDIFFU_HALO2%XRKIP2(JI,JJ,JK))
      IKIP2P1 = MIN(IKE,IKIP2+1)
      ZWFIP2 = MOD(PZDIFFU_HALO2%XRKIP2(JI,JJ,JK),1.)
      IKIP1 = INT(PZDIFFU_HALO2%XRKIP1(JI,JJ,JK))
      IKIP1P1 = MIN(IKE,IKIP1+1)
      ZWFIP1 = MOD(PZDIFFU_HALO2%XRKIP1(JI,JJ,JK),1.)
      IKIM2 = INT(PZDIFFU_HALO2%XRKIM2(JI,JJ,JK))
      IKIM2P1 = MIN(IKE,IKIM2+1)
      ZWFIM2 = MOD(PZDIFFU_HALO2%XRKIM2(JI,JJ,JK),1.)
      IKIM1 = INT(PZDIFFU_HALO2%XRKIM1(JI,JJ,JK))
      IKIM1P1 = MIN(IKE,IKIM1+1)
      ZWFIM1 = MOD(PZDIFFU_HALO2%XRKIM1(JI,JJ,JK),1.)
      ZZDIFFU(JI,JJ,JK) = ZWFIP2*    ZPTBFIELD(JI+2,JJ,IKIP2P1) + &
                         (1.-ZWFIP2)*ZPTBFIELD(JI+2,JJ,IKIP2) + &
                          ZWFIM2*    ZPTBFIELD(JI-2,JJ,IKIM2P1) + &
                         (1.-ZWFIM2)*ZPTBFIELD(JI-2,JJ,IKIM2) - 4*( &
                         ZWFIP1*    ZPTBFIELD(JI+1,JJ,IKIP1P1) + &
                         (1.-ZWFIP1)*ZPTBFIELD(JI+1,JJ,IKIP1) + &
                          ZWFIM1*    ZPTBFIELD(JI-1,JJ,IKIM1P1) + &
                         (1.-ZWFIM1)*ZPTBFIELD(JI-1,JJ,IKIM1) ) + &
                         6*ZPTBFIELD(JI,JJ,JK)

    ENDDO
  ENDDO
ENDDO

! c) Compute diffusion along model levels including the leading-order metric term
! It is sufficient to do this for vertical indices smaller than KZDLB because
! for higher levels no blending between z-diffusion and model-level diffusion is needed.

DO JI = IWZ,IEZ
  DO JJ = IJB-1,IJE+1
    DO JK = IKB,PZDIFFU_HALO2%NZDLB-1
      IKP1 = JK+1
      IKM1 = MAX(IKB,JK-1)
      ZMDIFFU(JI,JJ,JK) = ZPTBFIELD(JI+2,JJ,JK) + ZPTBFIELD(JI-2,JJ,JK) -4*( &
                          ZPTBFIELD(JI+1,JJ,JK) + ZPTBFIELD(JI-1,JJ,JK) ) +&
                          6*ZPTBFIELD(JI,JJ,JK) -  &
                         min(0.025,(ZPTBFIELD(JI,JJ,IKP1) - ZPTBFIELD(JI,JJ,IKM1))/ &
                         (PZDIFFU_HALO2%XZZ(JI,JJ,IKP1) - PZDIFFU_HALO2%XZZ(JI,JJ,IKM1)))* &
                         (PZDIFFU_HALO2%XZZ(JI+2,JJ,JK) + PZDIFFU_HALO2%XZZ(JI-2,JJ,JK) -4*( &
                          PZDIFFU_HALO2%XZZ(JI+1,JJ,JK) + PZDIFFU_HALO2%XZZ(JI-1,JJ,JK) ) +&
                          6*PZDIFFU_HALO2%XZZ(JI,JJ,JK) )
! Weighting factor for z-diffusion
      ZWGTFAC(JI,JJ,JK) = MAX(0,JK-PZDIFFU_HALO2%NZDI(JI,JJ)+1)/&
                        FLOAT(PZDIFFU_HALO2%NZDLB-PZDIFFU_HALO2%NZDI(JI,JJ)+1)
    ENDDO
  ENDDO
ENDDO

PRFIELDS(IWZ:IEZ,IJB:IJE,PZDIFFU_HALO2%NZDLB:IKE) = &
            PRFIELDS(IWZ:IEZ,IJB:IJE,PZDIFFU_HALO2%NZDLB:IKE) - &
            PRHODJ(IWZ:IEZ,IJB:IJE,PZDIFFU_HALO2%NZDLB:IKE)   * &
            PDK4*ZZDIFFU(IWZ:IEZ,IJB:IJE,PZDIFFU_HALO2%NZDLB:IKE)
DO JK = IKB,PZDIFFU_HALO2%NZDLB-1
  PRFIELDS(IWZ:IEZ,IJB:IJE,JK) = PRFIELDS(IWZ:IEZ,IJB:IJE,JK)-PRHODJ(IWZ:IEZ,IJB:IJE,JK) *   &
          PDK4*(ZWGTFAC(IWZ:IEZ,IJB:IJE,JK)*ZZDIFFU(IWZ:IEZ,IJB:IJE,JK) + &
           (1.-ZWGTFAC(IWZ:IEZ,IJB:IJE,JK))*PZDIFFU_HALO2%XREDFACI(IWZ:IEZ,IJB:IJE) * &
                                            ZMDIFFU(IWZ:IEZ,IJB:IJE,JK) )
ENDDO

!*       2.     COMPUTE THE 4TH ORDER DIFFUSION TERMS IN THE Y DIRECTION:
!               ---------------------------------------------------------
!
IF ( .NOT. L2D ) THEN
  SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
!
!*       2.1    CYCLIC CASE IN THE Y DIRECTION:
!
  CASE ('CYCL')          ! In that case one must have HLBCY(1) == HLBCY(2)
!
!
!!$    IF(NHALO == 1) THEN
      IS=IJB+1
      IN=IJE-1
!!$    ELSE
!!$      IS=IJB
!!$      IN=IJE
!!$    END IF


IF (PRESENT(PLSFIELD)) THEN
  ZPTBFIELD(IIB-1:IIE+1,IS-2:IN+2,:) = PFIELDM(IIB-1:IIE+1,IS-2:IN+2,:) - PLSFIELD(IIB-1:IIE+1,IS-2:IN+2,:)
!!$  IF(NHALO == 1) THEN
    ZPTBFIELD(IIB-1:IIE+1,IS-3,:) = TPHALO2%SOUTH(IIB-1:IIE+1,:) - TPHALO2LS%SOUTH(IIB-1:IIE+1,:)
    ZPTBFIELD(IIB-1:IIE+1,IN+3,:) = TPHALO2%NORTH(IIB-1:IIE+1,:) - TPHALO2LS%NORTH(IIB-1:IIE+1,:)
!!$  ENDIF
ELSE
  ZPTBFIELD(IIB-1:IIE+1,IS-2:IN+2,:) = PFIELDM(IIB-1:IIE+1,IS-2:IN+2,:)
!!$  IF(NHALO == 1) THEN
    ZPTBFIELD(IIB-1:IIE+1,IS-3,:) = TPHALO2%SOUTH(IIB-1:IIE+1,:)
    ZPTBFIELD(IIB-1:IIE+1,IN+3,:) = TPHALO2%NORTH(IIB-1:IIE+1,:)
!!$  ENDIF
ENDIF

!!! HALOS ARE PROBABLY ALSO NEEDED FOR THE INTERPOLATION COEFFICIENTS??!!


!
!
!*       2.2    NON CYCLIC CASE IN THE Y DIRECTION
!
  CASE ('OPEN','WALL','NEST')
!
    IF (LSOUTH_ll()) THEN
      IF(KGRID == 3) THEN
        IS=IJB+2          ! special case of C grid
      ELSE
        IS=IJB+1
      END IF
    ELSE
!!$      IF(NHALO == 1) THEN
        IS=IJB+1
!!$      ELSE
!!$        IS=IJB
!!$      ENDIF
    ENDIF
!!$    IF (LNORTH_ll() .OR. NHALO == 1) THEN
      IN=IJE-1
!!$    ELSE
!!$      IN=IJE
!!$    END IF

IF (PRESENT(PLSFIELD)) THEN
  ZPTBFIELD(IIB-1:IIE+1,IS-2:IN+2,:) = PFIELDM(IIB-1:IIE+1,IS-2:IN+2,:) - PLSFIELD(IIB-1:IIE+1,IS-2:IN+2,:)
!!$  IF((NHALO == 1).AND.(.NOT.LSOUTH_ll())) THEN
  IF(.NOT.LSOUTH_ll()) THEN
    ZPTBFIELD(IIB-1:IIE+1,IS-3,:) = TPHALO2%SOUTH(IIB-1:IIE+1,:) - TPHALO2LS%SOUTH(IIB-1:IIE+1,:)
  ENDIF
!!$  IF((NHALO == 1).AND.(.NOT.LNORTH_ll())) THEN
  IF(.NOT.LNORTH_ll()) THEN
    ZPTBFIELD(IIB-1:IIE+1,IN+3,:) = TPHALO2%NORTH(IIB-1:IIE+1,:) - TPHALO2LS%NORTH(IIB-1:IIE+1,:)
  ENDIF
ELSE
  ZPTBFIELD(IIB-1:IIE+1,IS-2:IN+2,:) = PFIELDM(IIB-1:IIE+1,IS-2:IN+2,:)
!!$  IF((NHALO == 1).AND.(.NOT.LSOUTH_ll())) THEN
  IF(.NOT.LSOUTH_ll()) THEN
    ZPTBFIELD(IIB-1:IIE+1,IS-3,:) = TPHALO2%SOUTH(IIB-1:IIE+1,:)
  ENDIF
!!$  IF((NHALO == 1).AND.(.NOT.LNORTH_ll())) THEN
  IF(.NOT.LNORTH_ll()) THEN
    ZPTBFIELD(IIB-1:IIE+1,IN+3,:) = TPHALO2%NORTH(IIB-1:IIE+1,:)
  ENDIF
ENDIF


      IF (LSOUTH_ll()) THEN
!
        PRFIELDS(IIB:IIE,IS-1,:) = PRFIELDS(IIB:IIE,IS-1,:) + PRHODJ(IIB:IIE,IS-1,:) *       &
          PDK2*(ZPTBFIELD(IIB:IIE,IS-2,:)-2.*ZPTBFIELD(IIB:IIE,IS-1,:)+ ZPTBFIELD(IIB:IIE,IS,:) )

      ENDIF
!
      IF (LNORTH_ll()) THEN
!
        PRFIELDS(IIB:IIE,IN+1,:) = PRFIELDS(IIB:IIE,IN+1,:) + PRHODJ(IIB:IIE,IN+1,:) *       &
          PDK2*(ZPTBFIELD(IIB:IIE,IN,:)-2.*ZPTBFIELD(IIB:IIE,IN+1,:)+ ZPTBFIELD(IIB:IIE,IN+2,:) )

      ENDIF

!
  END SELECT


! Compute diffusion kernel in a truly horizontally way (y direction)

! a) Determine E/W boundaries

!!$IF ((NHALO == 1).AND.(HLBCY(1) == 'CYCL').OR.((.NOT.LSOUTH_ll()).AND.(NHALO == 1)) ) THEN
IF ((HLBCY(1) == 'CYCL').OR.(.NOT.LSOUTH_ll()) ) THEN
  ISZ = IS-1
ELSE
  ISZ = IS
ENDIF

!!$IF ((NHALO == 1).AND.(HLBCY(1) == 'CYCL').OR.((.NOT.LNORTH_ll()).AND.(NHALO == 1)) ) THEN
IF ((HLBCY(1) == 'CYCL').OR.( .NOT.LNORTH_ll() ) ) THEN
  INZ = IN+1
ELSE
  INZ = IN
ENDIF

!!! ATTENTION: HALOS MIGHT BE NEEDED FOR THE INTERPOLATION COEFFICIENTS IN THE CASE
!!! OF ISZ = IS-1 / INZ = IN+1 (??)

! b) Compute truly horizontal diffusion


DO JI = IIB-1,IIE+1
  DO JJ = ISZ,INZ
    DO JK = IKB,IKE
! Calculate vertical indices and weighting factors for remote points

      IKJP2 = INT(PZDIFFU_HALO2%XRKJP2(JI,JJ,JK))
      IKJP2P1 = MIN(IKE,IKJP2+1)
      ZWFJP2 = MOD(PZDIFFU_HALO2%XRKJP2(JI,JJ,JK),1.)
      IKJP1 = INT(PZDIFFU_HALO2%XRKJP1(JI,JJ,JK))
      IKJP1P1 = MIN(IKE,IKJP1+1)
      ZWFJP1 = MOD(PZDIFFU_HALO2%XRKJP1(JI,JJ,JK),1.)
      IKJM2 = INT(PZDIFFU_HALO2%XRKJM2(JI,JJ,JK))
      IKJM2P1 = MIN(IKE,IKJM2+1)
      ZWFJM2 = MOD(PZDIFFU_HALO2%XRKJM2(JI,JJ,JK),1.)
      IKJM1 = INT(PZDIFFU_HALO2%XRKJM1(JI,JJ,JK))
      IKJM1P1 = MIN(IKE,IKJM1+1)
      ZWFJM1 = MOD(PZDIFFU_HALO2%XRKJM1(JI,JJ,JK),1.)
      ZZDIFFU(JI,JJ,JK) = ZWFJP2*    ZPTBFIELD(JI,JJ+2,IKJP2P1) + &
                         (1.-ZWFJP2)*ZPTBFIELD(JI,JJ+2,IKJP2) + &
                          ZWFJM2*    ZPTBFIELD(JI,JJ-2,IKJM2P1) + &
                         (1.-ZWFJM2)*ZPTBFIELD(JI,JJ-2,IKJM2) - 4*( &
                          ZWFJP1*    ZPTBFIELD(JI,JJ+1,IKJP1P1) + &
                         (1.-ZWFJP1)*ZPTBFIELD(JI,JJ+1,IKJP1) + &
                          ZWFJM1*    ZPTBFIELD(JI,JJ-1,IKJM1P1) + &
                         (1.-ZWFJM1)*ZPTBFIELD(JI,JJ-1,IKJM1) ) + &
                         6*ZPTBFIELD(JI,JJ,JK)

    ENDDO
  ENDDO
ENDDO

! c) Compute diffusion along model levels including the leading-order metric term
! It is sufficient to do this for vertical indices smaller than KZDLB because
! for higher levels no blending between z-diffusion and model-level diffusion is needed.

DO JI = IIB-1,IIE+1
  DO JJ = ISZ,INZ
    DO JK = IKB,PZDIFFU_HALO2%NZDLB-1
      IKP1 = JK+1
      IKM1 = MAX(IKB,JK-1)
      ZMDIFFU(JI,JJ,JK) = ZPTBFIELD(JI,JJ+2,JK) + ZPTBFIELD(JI,JJ-2,JK) -4*( &
                          ZPTBFIELD(JI,JJ+1,JK) + ZPTBFIELD(JI,JJ-1,JK) ) +&
                          6*ZPTBFIELD(JI,JJ,JK) -  &
                         min(0.025,(ZPTBFIELD(JI,JJ,IKP1) - ZPTBFIELD(JI,JJ,IKM1))/ &
                         (PZDIFFU_HALO2%XZZ(JI,JJ,IKP1) - PZDIFFU_HALO2%XZZ(JI,JJ,IKM1)))* &
                         (PZDIFFU_HALO2%XZZ(JI,JJ+2,JK) + PZDIFFU_HALO2%XZZ(JI,JJ-2,JK) -4*( &
                          PZDIFFU_HALO2%XZZ(JI,JJ+1,JK) + PZDIFFU_HALO2%XZZ(JI,JJ-1,JK) ) +&
                          6*PZDIFFU_HALO2%XZZ(JI,JJ,JK) )
! Weighting factor for z-diffusion
      ZWGTFAC(JI,JJ,JK) = MAX(0,JK-PZDIFFU_HALO2%NZDJ(JI,JJ)+1)/ & 
                        FLOAT(PZDIFFU_HALO2%NZDLB-PZDIFFU_HALO2%NZDJ(JI,JJ)+1)
    ENDDO
  ENDDO
ENDDO

PRFIELDS(IIB:IIE,ISZ:INZ,PZDIFFU_HALO2%NZDLB:IKE) = &
           PRFIELDS(IIB:IIE,ISZ:INZ,PZDIFFU_HALO2%NZDLB:IKE) - &
           PRHODJ(IIB:IIE,ISZ:INZ,PZDIFFU_HALO2%NZDLB:IKE)   * &
           PDK4*ZZDIFFU(IIB:IIE,ISZ:INZ,PZDIFFU_HALO2%NZDLB:IKE)
DO JK = IKB,PZDIFFU_HALO2%NZDLB-1
  PRFIELDS(IIB:IIE,ISZ:INZ,JK) = PRFIELDS(IIB:IIE,ISZ:INZ,JK)-PRHODJ(IIB:IIE,ISZ:INZ,JK) *   &
          PDK4*(ZWGTFAC(IIB:IIE,ISZ:INZ,JK)*ZZDIFFU(IIB:IIE,ISZ:INZ,JK) + &
           (1.-ZWGTFAC(IIB:IIE,ISZ:INZ,JK))*PZDIFFU_HALO2%XREDFACJ(IIB:IIE,ISZ:INZ) * &
                                            ZMDIFFU(IIB:IIE,ISZ:INZ,JK) )
ENDDO

!
ENDIF  ! 2D

DEALLOCATE (ZPTBFIELD,ZZDIFFU,ZMDIFFU,ZWGTFAC,ZZMASS)


!
!-------------------------------------------------------------------------------
END SUBROUTINE NUM_DIFF_ALGO_Z

!!! Original diffusion routine
!
!
!     ############################################################
      SUBROUTINE NUM_DIFF_ALGO(PRFIELDS, PFIELDM, KGRID, PRHODJ, &
                               PDK2, PDK4, PLSFIELD, TPHALO2, TPHALO2LS) 
!     ############################################################
!!
!!****  *NUM_DIFF_ALGO * - routine to apply an horizontal diffusion to the
!!                    prognostic variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to add a 2sd or 4th order horizontal
!!    diffusion term to an advected variable.
!!
!!**  METHOD
!!    ------
!!      In case of cyclic LBCs, the routine adds a numerical diffusion
!!    term obtained by applying an horizontal nabla4 operator to the
!!    prognostic variable on each grid level. In case of open LBCs,
!!    the nabla4 operator degenerates to a nabla2 operator on the first
!!    ring inside the computationnal domain.
!!      In the current version no budget computation is performed.
!!      The "halo2" (or the second layer of the halo) of the prognostic
!!    variable is passed as argument. The "halo2" of the LS field
!!    is also passed as argument, if necessary.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    Module MODD_PARALLEL2
!!
!!      HALO2_ll : type for a set of arrays representing 
!!                 the second layer of the halo
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine NUM_DIFF )
!!      User Interface for the MesoNH Parallel Package
!!
!!    AUTHOR
!!    ------
!!
!!      E. Richard       * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   15/10/94
!!      J. Stein   08/12/94 change the variable to be diffused rho*J*phi => phi
!!      J. Stein   20/03/95 remove R from the historical variables + switch 2D
!!                          + switch for TKE unused
!!                 01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                 16/10/95 (J. Stein)     change the budget calls
!!                 19/12/96 (J.-P. Pinty)  update the budget calls
!!                 19/06/98 (P. Kloos)     restructure to run in parallel
!!                 12/10    (J.Escobar) replace DX4 & DY4 by inline code for reproductibility 
!!                 12/10    (J.Escobar) reorder all operator from - to + index for reproductibility
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!

USE MODE_ll
!
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODI_NABLA4
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRFIELDS ! source term
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFIELDM  ! variable at t-dt
INTEGER,                INTENT(IN)    :: KGRID    ! C grid localisation
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ   ! RHODJ at KGRID
REAL,                   INTENT(IN)    :: PDK2  ! 2nd order dif. coef. /dx2
REAL,                   INTENT(IN)    :: PDK4  ! 4th order dif. coef. /dx4
!
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PLSFIELD ! Larger Scale variable
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t-dt
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2LS    ! halo2 for the LS field
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IW,IE,IS,IN   ! Coordinate of forth order diffusion area
!
!-------------------------------------------------------------------------------
!
!*       1.     CALCULATE THE NUMERICAL DIFFUSION TERM IN THE X DIRECTION
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
!!$  IF(NHALO == 1) THEN
    IW=IIB+1
    IE=IIE-1
!!$  ELSE
!!$    IW=IIB
!!$    IE=IIE
!!$  END IF  
!
  IF (PRESENT(PLSFIELD)) THEN
!
!*       1.1.1  Case with LS fields
!
!!$    IF(NHALO == 1) THEN
!
      PRFIELDS(IW-1,:,:) = PRFIELDS(IW-1,:,:) - PRHODJ(IW-1,:,:) *       &
      PDK4*(                                                             &
                TPHALO2%WEST(:,:)   +   PFIELDM(IW+1,:,:)                &
          -4.*( PFIELDM(IW-2,:,:)   +   PFIELDM(IW,:,:)     )            &
          +6.*  PFIELDM(IW-1,:,:)                                        &
            -   TPHALO2LS%WEST(:,:) -   PLSFIELD(IW+1,:,:)               &
          +4.*( PLSFIELD(IW-2,:,:)  +   PLSFIELD(IW,:,:)   )             &
          -6.*  PLSFIELD(IW-1,:,:)  )
!
      PRFIELDS(IE+1,:,:) = PRFIELDS(IE+1,:,:) - PRHODJ(IE+1,:,:) *       &
      PDK4*(                                                             &
                PFIELDM(IE-1,:,:)   +  TPHALO2%EAST(:,:)                 &
          -4.*( PFIELDM(IE,:,:)     +  PFIELDM(IE+2,:,:)   )             &
          +6.*  PFIELDM(IE+1,:,:)                                        &
             -  PLSFIELD(IE-1,:,:)  -  TPHALO2LS%EAST(:,:)               &
          +4.*( PLSFIELD(IE,:,:)    +  PLSFIELD(IE+2,:,:)  )             &
          -6.*  PLSFIELD(IE+1,:,:)  )
!
!!$    ENDIF
!
!!$    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:)-PRHODJ(IW:IE,:,:) *   &
!!$          PDK4*DX4(PFIELDM(IW-2:IE+2,:,:)-PLSFIELD(IW-2:IE+2,:,:))

    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:) - PRHODJ(IW:IE,:,:) *       &
         PDK4*(                                                           &
                  PFIELDM(IW-2:IE-2,:,:)  + PFIELDM(IW+2:IE+2,:,:)        &
            -4.*( PFIELDM(IW-1:IE-1,:,:)  + PFIELDM(IW+1:IE+1,:,:)   )    &
            +6.*  PFIELDM(IW:IE,:,:)                                      &                                                                         
              -   PLSFIELD(IW-2:IE-2,:,:) - PLSFIELD(IW+2:IE+2,:,:)       &
            +4.*( PLSFIELD(IW-1:IE-1,:,:) + PLSFIELD(IW+1:IE+1,:,:)  )    &
            -6.*  PLSFIELD(IW:IE,:,:))
!
  ELSE
!
!*       1.1.2  Case without LS fields
!
!!$    IF(NHALO == 1) THEN
!
      PRFIELDS(IW-1,:,:) = PRFIELDS(IW-1,:,:) - PRHODJ(IW-1,:,:) *       &
      PDK4*(                                                             &
                TPHALO2%WEST(:,:)  +  PFIELDM(IW+1,:,:)                  &
         -4.*(  PFIELDM(IW-2,:,:)  +  PFIELDM(IW,:,:)     )              &
         +6.*   PFIELDM(IW-1,:,:)  )
!
      PRFIELDS(IE+1,:,:) = PRFIELDS(IE+1,:,:) - PRHODJ(IE+1,:,:) *       &
      PDK4*(                                                             &
                PFIELDM(IE-1,:,:)  +  TPHALO2%EAST(:,:)                  &
         -4.*(  PFIELDM(IE,:,:)    +  PFIELDM(IE+2,:,:)    )             &
         +6.*   PFIELDM(IE+1,:,:)  )
!
!!$    ENDIF
!
!!$    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:)-PRHODJ(IW:IE,:,:) *   &
!!$          PDK4*DX4(PFIELDM(IW-2:IE+2,:,:))

    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:) - PRHODJ(IW:IE,:,:) *       &
         PDK4*(                                                           &
                  PFIELDM(IW-2:IE-2,:,:)  + PFIELDM(IW+2:IE+2,:,:)        &
            -4.*( PFIELDM(IW-1:IE-1,:,:)  + PFIELDM(IW+1:IE+1,:,:)   )    &
            +6.*  PFIELDM(IW:IE,:,:)     )                                 
                                                                          
    
!
  ENDIF
!
!
!*       1.2    NON CYCLIC CASE IN THE X DIRECTION 
!
CASE ('OPEN','WALL','NEST')
!
  IF (LWEST_ll()) THEN
    IF(KGRID == 2) THEN
      IW=IIB+2          ! special case of C grid
    ELSE
      IW=IIB+1
    END IF
  ELSE
!!$    IF(NHALO == 1) THEN
      IW=IIB+1
!!$    ELSE
!!$      IW=IIB
!!$    ENDIF
  ENDIF
!!$  IF (LEAST_ll() .OR. NHALO == 1) THEN
    IE=IIE-1
!!$  ELSE
!!$    IE=IIE
!!$  END IF  
!
  IF (PRESENT(PLSFIELD)) THEN
!
!*       1.2.1  Case with LS fields
!
!*             Use a second order scheme at the physical border
!
    IF (LWEST_ll()) THEN
!
      PRFIELDS(IW-1,:,:) = PRFIELDS(IW-1,:,:) + PRHODJ(IW-1,:,:) *              &
        PDK2*(                                                                  &
                 PFIELDM(IW-2,:,:)  -2.*PFIELDM(IW-1,:,:)  + PFIELDM(IW,:,:)    &
                -PLSFIELD(IW-2,:,:) +2.*PLSFIELD(IW-1,:,:) - PLSFIELD(IW,:,:)   )
!
!!$    ELSEIF (NHALO == 1) THEN
    ELSE
!
      PRFIELDS(IW-1,:,:) = PRFIELDS(IW-1,:,:) - PRHODJ(IW-1,:,:) *        &
      PDK4*(                                                              &
                 TPHALO2%WEST(:,:)   +  PFIELDM(IW+1,:,:)                 &
          -4.*(  PFIELDM(IW-2,:,:)   +  PFIELDM(IW,:,:)     )             &
          +6.*   PFIELDM(IW-1,:,:)                                        &
            -    TPHALO2LS%WEST(:,:) -  PLSFIELD(IW+1,:,:)                &
          +4.*(  PLSFIELD(IW-2,:,:)  +  PLSFIELD(IW,:,:)    )             &
          -6.*   PLSFIELD(IW-1,:,:)  ) 
!     
    ENDIF
!
    IF (LEAST_ll()) THEN
!
      PRFIELDS(IE+1,:,:) = PRFIELDS(IE+1,:,:) + PRHODJ(IE+1,:,:) *            &
        PDK2*(                                                                & 
                PFIELDM(IE,:,:)  -2.*PFIELDM(IE+1,:,:)  + PFIELDM(IE+2,:,:)   &
              - PLSFIELD(IE,:,:) +2.*PLSFIELD(IE+1,:,:) - PLSFIELD(IE+2,:,:)  )
!
!!$    ELSEIF (NHALO == 1) THEN
    ELSE
!
      PRFIELDS(IE+1,:,:) = PRFIELDS(IE+1,:,:) - PRHODJ(IE+1,:,:) *      &
      PDK4*(                                                            &
                 PFIELDM(IE-1,:,:)  +   TPHALO2%EAST(:,:)               &
          -4.*(  PFIELDM(IE  ,:,:)  +   PFIELDM(IE+2,:,:)    )          &
          +6.*   PFIELDM(IE+1,:,:)                                      &
            -    PLSFIELD(IE-1,:,:) -   TPHALO2LS%EAST(:,:)             &
          +4.*(  PLSFIELD(IE  ,:,:) +   PLSFIELD(IE+2,:,:))             &
          -6.*   PLSFIELD(IE+1,:,:))
!
    ENDIF

!*             Use a fourth order scheme 
!
!!$    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:)-PRHODJ(IW:IE,:,:) *   &
!!$          PDK4*DX4(PFIELDM(IW-2:IE+2,:,:)-PLSFIELD(IW-2:IE+2,:,:))

    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:) - PRHODJ(IW:IE,:,:) *       &
         PDK4*(                                                           &
                  PFIELDM(IW-2:IE-2,:,:)  + PFIELDM(IW+2:IE+2,:,:)        &
            -4.*( PFIELDM(IW-1:IE-1,:,:)  + PFIELDM(IW+1:IE+1,:,:)   )    &
            +6.*  PFIELDM(IW:IE,:,:)                                      & 
              -   PLSFIELD(IW-2:IE-2,:,:) - PLSFIELD(IW+2:IE+2,:,:)       &
            +4.*( PLSFIELD(IW-1:IE-1,:,:) + PLSFIELD(IW+1:IE+1,:,:)  )    &
            -6.*  PLSFIELD(IW:IE,:,:)) 
!
  ELSE
!
!*       1.2.2  Case without LS fields
!
!*             Use a second order scheme at the physical border
!
    IF (LWEST_ll()) THEN
!
      PRFIELDS(IW-1,:,:) = PRFIELDS(IW-1,:,:) + PRHODJ(IW-1,:,:) *       &
        PDK2*( PFIELDM(IW-2,:,:) -2.*PFIELDM(IW-1,:,:) + PFIELDM(IW,:,:) )
!
!!$    ELSEIF (NHALO == 1) THEN
    ELSE
!
      PRFIELDS(IW-1,:,:) = PRFIELDS(IW-1,:,:) - PRHODJ(IW-1,:,:) *       &
      PDK4*(                                                             & 
                 TPHALO2%WEST(:,:)  +  PFIELDM(IW+1,:,:)                 &
          -4.*(  PFIELDM(IW-2,:,:)  +  PFIELDM(IW,:,:)    )              &
          +6.*   PFIELDM(IW-1,:,:)   )
!
    ENDIF
!
    IF (LEAST_ll()) THEN
!
      PRFIELDS(IE+1,:,:) = PRFIELDS(IE+1,:,:) + PRHODJ(IE+1,:,:) *       &
        PDK2*( PFIELDM(IE,:,:) -2.*PFIELDM(IE+1,:,:) + PFIELDM(IE+2,:,:) )
!
!!$    ELSEIF (NHALO == 1) THEN
    ELSE
!
      PRFIELDS(IE+1,:,:) = PRFIELDS(IE+1,:,:) - PRHODJ(IE+1,:,:) *       &
      PDK4*(                                                             &
                PFIELDM(IE-1,:,:) + TPHALO2%EAST(:,:)                    &
          -4.*( PFIELDM(IE,:,:)   + PFIELDM(IE+2,:,:)  )                 &
          +6.*  PFIELDM(IE+1,:,:)  )
!
    ENDIF

!*             Use a fourth order scheme 
!
!!$    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:)-PRHODJ(IW:IE,:,:) * &
!!$                           PDK4*DX4(PFIELDM(IW-2:IE+2,:,:))

    PRFIELDS(IW:IE,:,:) = PRFIELDS(IW:IE,:,:) - PRHODJ(IW:IE,:,:) *       &
         PDK4*(                                                           &
                  PFIELDM(IW-2:IE-2,:,:)  + PFIELDM(IW+2:IE+2,:,:)        &
            -4.*( PFIELDM(IW-1:IE-1,:,:)  + PFIELDM(IW+1:IE+1,:,:)   )    &
            +6.*  PFIELDM(IW:IE,:,:)  )

!
!
  ENDIF
!
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTES THE 4TH ORDER DIFFUSION TERMS IN THE Y DIRECTION:
!               ---------------------------------------------------------
!
IF ( .NOT. L2D ) THEN
  SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
!
!*       2.1    CYCLIC CASE IN THE Y DIRECTION:
!
  CASE ('CYCL')          ! In that case one must have HLBCY(1) == HLBCY(2)
!
!
!!$    IF(NHALO == 1) THEN
      IS=IJB+1
      IN=IJE-1
!!$    ELSE
!!$      IS=IJB
!!$      IN=IJE
!!$    END IF
!
    IF (PRESENT(PLSFIELD)) THEN
!
!*       2.1.1  Case with LS fields
!
!!$      IF(NHALO == 1) THEN
!
        PRFIELDS(:,IS-1,:) = PRFIELDS(:,IS-1,:) - PRHODJ(:,IS-1,:) *      &
        PDK4*(                                                            &
                  TPHALO2%SOUTH(:,:)   +  PFIELDM(:,IS+1,:)               &
            -4.*( PFIELDM(:,IS-2,:)    +  PFIELDM(:,IS,:)    )            &
            +6.*  PFIELDM(:,IS-1,:)                                       &
              -   TPHALO2LS%SOUTH(:,:) -  PLSFIELD(:,IS+1,:)              &
            +4.*( PLSFIELD(:,IS-2,:)   +  PLSFIELD(:,IS,:)   )            &
            -6.*  PLSFIELD(:,IS-1,:)   )
!
        PRFIELDS(:,IN+1,:) = PRFIELDS(:,IN+1,:) - PRHODJ(:,IN+1,:) *      &
        PDK4*(                                                            &
                  PFIELDM(:,IN-1,:)    +  TPHALO2%NORTH(:,:)              &
            -4.*( PFIELDM(:,IN,:)      +  PFIELDM(:,IN+2,:)  )            &
            +6.*  PFIELDM(:,IN+1,:)                                       &
               -  PLSFIELD(:,IN-1,:)   -  TPHALO2LS%NORTH(:,:)            &
            +4.*( PLSFIELD(:,IN,:)     +  PLSFIELD(:,IN+2,:) )            &
            -6.*  PLSFIELD(:,IN+1,:)   )
!
!!$      ENDIF
!
!!$      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:)-PRHODJ(:,IS:IN,:) *   &
!!$            PDK4*DY4(PFIELDM(:,IS-2:IN+2,:)-PLSFIELD(:,IS-2:IN+2,:))

      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:) - PRHODJ(:,IS:IN,:) *      &
           PDK4*(                                                          &
                     PFIELDM(:,IS-2:IN-2,:)  +  PFIELDM(:,IS+2:IN+2,:)     &
             -4.*(   PFIELDM(:,IS-1:IN-1,:)  +  PFIELDM(:,IS+1:IN+1,:)  )  &
             +6.*    PFIELDM(:,IS:IN,:)                                    &
                   - PLSFIELD(:,IS-2:IN-2,:) -  PLSFIELD(:,IS+2:IN+2,:)    &
             +4.*(   PLSFIELD(:,IS-1:IN-1,:) +  PLSFIELD(:,IS+1:IN+1,:) )  &
             -6.*    PLSFIELD(:,IS:IN,:) )
!
    ELSE
!
!*       2.1.2  Case without LS fields
!
!
!!$      IF(NHALO == 1) THEN
!
        PRFIELDS(:,IS-1,:) = PRFIELDS(:,IS-1,:) - PRHODJ(:,IS-1,:) *      &
        PDK4*(                                                            &
                  TPHALO2%SOUTH(:,:)   +  PFIELDM(:,IS+1,:)               &
            -4.*( PFIELDM(:,IS-2,:)    +  PFIELDM(:,IS,:)    )            &
            +6.*  PFIELDM(:,IS-1,:)    )
!
        PRFIELDS(:,IN+1,:) = PRFIELDS(:,IN+1,:) - PRHODJ(:,IN+1,:) *      &
        PDK4*(                                                            &
                  PFIELDM(:,IN-1,:)    +  TPHALO2%NORTH(:,:)              &
            -4.*( PFIELDM(:,IN,:)      +  PFIELDM(:,IN+2,:)  )            &
            +6.*  PFIELDM(:,IN+1,:)    )
!
!!$      ENDIF
!
!!$      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:)-PRHODJ(:,IS:IN,:) *   &
!!$            PDK4*DY4(PFIELDM(:,IS-2:IN+2,:))

      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:) - PRHODJ(:,IS:IN,:) *      &
           PDK4*(                                                          &
                     PFIELDM(:,IS-2:IN-2,:)  +  PFIELDM(:,IS+2:IN+2,:)     &
             -4.*(   PFIELDM(:,IS-1:IN-1,:)  +  PFIELDM(:,IS+1:IN+1,:)  )  &
             +6.*    PFIELDM(:,IS:IN,:)  )

!
    ENDIF
!
!
!*       2.2    NON CYCLIC CASE IN THE Y DIRECTION
!
  CASE ('OPEN','WALL','NEST')
!
    IF (LSOUTH_ll()) THEN
      IF(KGRID == 3) THEN
        IS=IJB+2          ! special case of C grid
      ELSE
        IS=IJB+1
      END IF
    ELSE
!!$      IF(NHALO == 1) THEN
        IS=IJB+1
!!$      ELSE
!!$        IS=IJB
!!$      ENDIF
    ENDIF
!!$    IF (LNORTH_ll() .OR. NHALO == 1) THEN
      IN=IJE-1
!!$    ELSE
!!$      IN=IJE
!!$    END IF  
!*       2.2.1  Case with LS fields
!
    IF (PRESENT(PLSFIELD)) THEN
!
!*       2.2.1  Case with LS fields
!
!*             Use a second order scheme at the physical border
!
      IF (LSOUTH_ll()) THEN
!
        PRFIELDS(:,IS-1,:) = PRFIELDS(:,IS-1,:) + PRHODJ(:,IS-1,:) *            &
          PDK2*(                                                                &
                 PFIELDM(:,IS-2,:)  -2.*PFIELDM(:,IS-1,:)  + PFIELDM(:,IS,:)    &
                -PLSFIELD(:,IS-2,:) +2.*PLSFIELD(:,IS-1,:) - PLSFIELD(:,IS,:)   )
!
!!$      ELSEIF (NHALO == 1) THEN
      ELSE
!
        PRFIELDS(:,IS-1,:) = PRFIELDS(:,IS-1,:) - PRHODJ(:,IS-1,:) *            &
        PDK4*(                                                                  &
                  TPHALO2%SOUTH(:,:)   +  PFIELDM(:,IS+1,:)                     &
            -4.*( PFIELDM(:,IS-2,:)    +  PFIELDM(:,IS,:)    )                  &
            +6.*  PFIELDM(:,IS-1,:)                                             &
              -   TPHALO2LS%SOUTH(:,:) -  PLSFIELD(:,IS+1,:)                    &
            +4.*( PLSFIELD(:,IS-2,:)   +  PLSFIELD(:,IS,:)   )                  &
            -6.*  PLSFIELD(:,IS-1,:)   )
!
      ENDIF
!
      IF (LNORTH_ll()) THEN
!
        PRFIELDS(:,IN+1,:) = PRFIELDS(:,IN+1,:) + PRHODJ(:,IN+1,:) *             &
          PDK2*(                                                                 &
                   PFIELDM(:,IN,:)  -2.*PFIELDM(:,IN+1,:)   + PFIELDM(:,IN+2,:)  &
                  -PLSFIELD(:,IN,:) +2.*PLSFIELD(:,IN+1,:)  - PLSFIELD(:,IN+2,:) )
!
!!$      ELSEIF (NHALO == 1) THEN
      ELSE
!
        PRFIELDS(:,IN+1,:) = PRFIELDS(:,IN+1,:) - PRHODJ(:,IN+1,:) *         &
        PDK4*(                                                               &
                  PFIELDM(:,IN-1,:)    +  TPHALO2%NORTH(:,:)                 &
          -4.*(   PFIELDM(:,IN,:)      +  PFIELDM(:,IN+2,:)   )              &
          +6.*    PFIELDM(:,IN+1,:)                                          &
                - PLSFIELD(:,IN-1,:)   - TPHALO2LS%NORTH(:,:)                &
          +4.*(   PLSFIELD(:,IN,:)     + PLSFIELD(:,IN+2,:)   )              &
          -6.*    PLSFIELD(:,IN+1,:)   )
!
      ENDIF
!
!*             Use a fourth order scheme 
!
!!$      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:)-PRHODJ(:,IS:IN,:) *   &
!!$            PDK4*DY4(PFIELDM(:,IS-2:IN+2,:)-PLSFIELD(:,IS-2:IN+2,:))

      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:) - PRHODJ(:,IS:IN,:) *      &
           PDK4*(                                                          &
                     PFIELDM(:,IS-2:IN-2,:)  +  PFIELDM(:,IS+2:IN+2,:)     &
             -4.*(   PFIELDM(:,IS-1:IN-1,:)  +  PFIELDM(:,IS+1:IN+1,:)  )  &
             +6.*    PFIELDM(:,IS:IN,:)                                    &
                   - PLSFIELD(:,IS-2:IN-2,:) -  PLSFIELD(:,IS+2:IN+2,:)    &
             +4.*(   PLSFIELD(:,IS-1:IN-1,:) +  PLSFIELD(:,IS+1:IN+1,:) )  &
             -6.*    PLSFIELD(:,IS:IN,:)   )
      
!
    ELSE
!
!*       2.2.2  Case without LS fields 
!
!*             Use a second order scheme at the physical border
!
      IF (LSOUTH_ll()) THEN
!
        PRFIELDS(:,IS-1,:) = PRFIELDS(:,IS-1,:) + PRHODJ(:,IS-1,:) *       &
          PDK2*( PFIELDM(:,IS-2,:) -2.*PFIELDM(:,IS-1,:) + PFIELDM(:,IS,:) )
!
!!$      ELSEIF (NHALO == 1) THEN
      ELSE
!
        PRFIELDS(:,IS-1,:) = PRFIELDS(:,IS-1,:) - PRHODJ(:,IS-1,:) *      &
        PDK4*(                                                            &
                  TPHALO2%SOUTH(:,:)  +  PFIELDM(:,IS+1,:)                &
            -4.*( PFIELDM(:,IS-2,:)   +  PFIELDM(:,IS,:)    )             &
            +6.*  PFIELDM(:,IS-1,:)   )
!
      ENDIF
!
      IF (LNORTH_ll()) THEN
!
        PRFIELDS(:,IN+1,:) = PRFIELDS(:,IN+1,:) + PRHODJ(:,IN+1,:) *       &
          PDK2*( PFIELDM(:,IN,:) -2.*PFIELDM(:,IN+1,:) + PFIELDM(:,IN+2,:) )
!
!!$      ELSEIF (NHALO == 1) THEN
      ELSE
!
        PRFIELDS(:,IN+1,:) = PRFIELDS(:,IN+1,:) - PRHODJ(:,IN+1,:) *      &
        PDK4*(                                                            &
                  PFIELDM(:,IN-1,:)   +  TPHALO2%NORTH(:,:)               &
            -4.*( PFIELDM(:,IN,:)     +  PFIELDM(:,IN+2,:)   )            &
            +6.*  PFIELDM(:,IN+1,:)   )
!
      ENDIF
!
!*             Use a fourth order scheme 
!
!!$      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:)-PRHODJ(:,IS:IN,:) *   &
!!$            PDK4*DY4(PFIELDM(:,IS-2:IN+2,:))

      PRFIELDS(:,IS:IN,:) = PRFIELDS(:,IS:IN,:) - PRHODJ(:,IS:IN,:) *   &
           PDK4*(                                                          &
                     PFIELDM(:,IS-2:IN-2,:)  +  PFIELDM(:,IS+2:IN+2,:)     &
             -4.*(   PFIELDM(:,IS-1:IN-1,:)  +  PFIELDM(:,IS+1:IN+1,:)  )  &
             +6.*    PFIELDM(:,IS:IN,:)   )
!
    ENDIF
!
  END SELECT
!
ENDIF
!
!-------------------------------------------------------------------------------
END SUBROUTINE NUM_DIFF_ALGO
!
END SUBROUTINE NUM_DIFF

