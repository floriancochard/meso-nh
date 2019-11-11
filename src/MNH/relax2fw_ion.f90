!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ########################
      MODULE MODI_RELAX2FW_ION
!     ########################
!
INTERFACE
!
      SUBROUTINE RELAX2FW_ION( KTCOUNT, KMI, PTSTEP, PRHODJ, PSVM, KALBOT,  &
                               PALK, OMASK_RELAX, PKWRELAX, PRSVS )
!
INTEGER,                  INTENT(IN)    :: KTCOUNT! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KMI    ! model number
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! effective dry rho * Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM   ! Mixing ratio at t-dt         
INTEGER,                  INTENT(IN)    :: KALBOT !Vertical index corresponding
                                      ! to the absorbing layer base
REAL, DIMENSION(:),       INTENT(IN)    :: PALK   ! Function of the absorbing 
                                      ! layer damping coefficient defined for 
                                      ! u,v,and theta
LOGICAL, DIMENSION(:,:),  INTENT(IN)    :: OMASK_RELAX ! Mask for the locations
                                      ! where lateral relax. must be performed
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKWRELAX ! u, v and mass locations
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS    ! Source term
!
END SUBROUTINE RELAX2FW_ION
!
END INTERFACE
!
END MODULE MODI_RELAX2FW_ION
!
!     ##################################################################
      SUBROUTINE RELAX2FW_ION( KTCOUNT, KMI, PTSTEP, PRHODJ, PSVM, KALBOT,  &
                               PALK, OMASK_RELAX, PKWRELAX, PRSVS )
!     ##################################################################
!
!!****  *RELAXATION * - routine to apply a Rayleigh damping at the top  
!!                       and in the outermost verticals of the domain
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to damp the u, v, w, and theta fields  
!!    near the top of the domain and in the outermost vertical planes
!!    (vapor field added) to prevent spurious reflections on the top and on
!!    the side boundaries.
!!
!!
!!**  METHOD
!!    ------
!!      Implicit Rayleigh damping for the absorbing layer 
!!    upper part:
!!                                     -Kal(zhat)
!!    -Kal(zhat)(phi(t+dt)-LSphi)= ------------------ (phi(t-dt)-LSphi)
!!                                  1 + 2dt Kal(zhat)
!!    lateral part:
!!                                  -Krelax
!!    -Krelax(phi(t+dt)-LSphi)= ------------------ (phi(t-dt)-LSphi)
!!                                1 + 2dt Krelax
!!
!!     phi being u, v, w, theta or qv
!!     LSphi are the corresponding Larger Scale values or to the lateral 
!!     boundaries fields.
!!     The different sources terms are stored for the budget computations.
!!
!!
!!    EXTERNAL
!!    --------
!!      BUDGET              : Stores the different budget components
!!      GET_INDICE_ll       : get physical sub-domain bounds
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!       Module MODD_PARAMETERS : JPVEXT
!!
!!       Module MODD_CONF       : CCONF 
!!
!!       Module MODD_BUDGET     : NBUMOD,LBU_ENABLE,
!!                                LBU_RU,LBU_RV,LBU_RW,LBU_RTH,LBU_RRV,
!!                                LBU_RRC,LBU_RRR,LBU_RRI,LBU_RRH,LBU_RRG,
!!                                LBU_RTKE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RELAXATION )
!!
!!    AUTHOR
!!    ------
!! 
!!    M. Chong , LA,   Adapted from relaxation.f90    Dec. 2010
!!
!!    MODIFICATIONS
!!    -------------
!!      C.Lac, 07/11 : Avoid the horizontal relaxation if not father model
!!      C.Lac, 11/11 : Adaptation to FIT temporal scheme
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CONF 
USE MODD_BUDGET 
USE MODD_NSV, ONLY: NSV_ELECBEG, NSV_ELECEND
USE MODD_ELEC_n, ONLY: XCION_POS_FW, XCION_NEG_FW
!
USE MODE_ll
!
USE MODI_SHUMAN     
USE MODI_BUDGET     
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
INTEGER,                  INTENT(IN)    :: KTCOUNT! Temporal loop counter       
INTEGER,                  INTENT(IN)    :: KMI    ! model number
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! effective dry rho * Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM   ! Mixing ratio at t-dt         
INTEGER,                  INTENT(IN)    :: KALBOT !Vertical index corresponding
                                      ! to the absorbing layer base
REAL, DIMENSION(:),       INTENT(IN)    :: PALK   ! Function of the absorbing 
                                      ! layer damping coefficient defined for 
                                      ! u,v,and theta
LOGICAL, DIMENSION(:,:),  INTENT(IN)    :: OMASK_RELAX ! Mask for the locations
                                      ! where lateral relax. must be performed
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKWRELAX ! u, v and mass locations
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS    ! Source term
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK,JSV            ! Index for loops
INTEGER :: IIB,IJB,IIE,IJE   ! Bounds of the physical sub-domain             
INTEGER :: IKU,IKE           ! size of the array in the z direction and index
                             ! of the last inner mass point
INTEGER    :: IIBINT,IJBINT  ! Corner of relaxation domain 
INTEGER    :: IIEINT,IJEINT  ! global indices translated to local indices
!
REAL, DIMENSION(SIZE(PALK))          :: ZKV  ! Function of the upper absorbing 
                             ! layer damping coefficient for u,v,theta and qv
!
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2)) :: ZKH
                             ! Function of the lateral absorbing layer damping 
                             ! for u,v and mass points respectively
!  
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!	        -------------
IKU = SIZE(PSVM,3)
IKE = IKU - JPVEXT
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!
!-------------------------------------------------------------------------------
!
!*       2.     RELAXATION IN THE UPPER LAYERS
!	        ------------------------------
!
!*       2.1    set the top-level damping coef. (upstream or leapfrog)
!
  ZKV(:)  = PALK(:) / (1. - PTSTEP * PALK(:))
!
!
!*       2.2     applies the damping in the uppermost levels
!
  DO JK = KALBOT, IKE+1
    PRSVS(:,:,JK,NSV_ELECBEG) = PRSVS(:,:,JK,NSV_ELECBEG) - ZKV(JK) *  &
       (PSVM(:,:,JK,NSV_ELECBEG) - XCION_POS_FW(:,:,JK)) * PRHODJ(:,:,JK)
    PRSVS(:,:,JK,NSV_ELECEND) = PRSVS(:,:,JK,NSV_ELECEND) - ZKV(JK) *  &
       (PSVM(:,:,JK,NSV_ELECEND) - XCION_NEG_FW(:,:,JK)) * PRHODJ(:,:,JK)
  END DO  
!
!  
!-------------------------------------------------------------------------------
!
!*       3.     RELAXATION IN THE OUTERMOST VERTICAL PLANES
!	        -------------------------------------------
!
!*       3.1    set the rim zone damping coef. (upstream or leapfrog)
!                  only for the father model
!
IF (KMI == 1) THEN
!
  ZKH(:,:) = PKWRELAX(:,:) / (1. - PTSTEP * PKWRELAX(:,:))
!
!
!*       3.2    applies the damping near the lateral boundaries
!
DO JK = 1, IKU
  WHERE (OMASK_RELAX)
    PRSVS(:,:,JK,NSV_ELECBEG) = PRSVS(:,:,JK,NSV_ELECBEG) - ZKH(:,:) *  &
       (PSVM(:,:,JK,NSV_ELECBEG) - XCION_POS_FW(:,:,JK)) * PRHODJ(:,:,JK)
    PRSVS(:,:,JK,NSV_ELECEND) = PRSVS(:,:,JK,NSV_ELECEND) - ZKH(:,:) *  &
       (PSVM(:,:,JK,NSV_ELECEND) - XCION_NEG_FW(:,:,JK)) * PRHODJ(:,:,JK)
  END WHERE
ENDDO 
!
END IF 
!
!
!-------------------------------------------------------------------------------
!
!*       4.     STORES FIELDS IN BUDGET ARRAYS
!	        ------------------------------
!
IF (LBUDGET_SV) THEN
  JSV = NSV_ELECBEG
  CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'REL_BU_RSV')
  JSV = NSV_ELECEND
  CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'REL_BU_RSV')
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RELAX2FW_ION
