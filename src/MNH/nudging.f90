!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/24 18:05:52
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_NUDGING
!     ###################
!
INTERFACE
!
      SUBROUTINE NUDGING ( OUSERV, PRHODJ, PTNUDGING,                      &
                           PUM,  PVM,  PWM, PTHM, PRM,                     &
                           PLSUM,  PLSVM,  PLSWM,  PLSTHM,  PLSRVM,        &
                           PRUS, PRVS, PRWS, PRTHS, PRRS                   )
!
!
LOGICAL               , INTENT(IN) :: OUSERV  ! Logical to use rv
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! ( rhod J ) = dry density
              ! for reference state * Jacobian of the GCS transformation.
REAL                  , INTENT(IN) :: PTNUDGING  ! Nudging time scale
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUM,PVM,PWM,PTHM 
                                          ! wind, potential temperature and
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM !  moist variables at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM
                                          ! Large Scale fields
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS,PRVS,PRWS,PRTHS
                                          ! wind, potential temperature and
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS ! moist variables at time t+1
!
END SUBROUTINE NUDGING
!
END INTERFACE
!
END MODULE MODI_NUDGING
!
!     ######################################################################
      SUBROUTINE NUDGING ( OUSERV, PRHODJ, PTNUDGING,                      &
                           PUM,  PVM,  PWM, PTHM, PRM,                     &
                           PLSUM,  PLSVM,  PLSWM,  PLSTHM,  PLSRVM,        &
                           PRUS, PRVS, PRWS, PRTHS, PRRS                   )
!     ######################################################################
!
!!***  *NUDGING* - routine to compute the nudging terms
!!
!!    PURPOSE
!!    -------
!!      The purpuse is to compute a source term (nudging towards large scale
!!     field) for U,V,W,TH and Rv.
!!   
!!**  METHOD
!!    ------
!!      The nudging term writes:
!!
!!       
!!               d (RHODJ * A)                          t-1   
!!               -------------- = -(1/Tau)  * RHODJ * (A    - Als )
!!                    dt
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson      * Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/05/06
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
!
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL               , INTENT(IN) :: OUSERV  ! Logical to use rv
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ  ! ( rhod J ) = dry density
              ! for reference state * Jacobian of the GCS transformation.
REAL                  , INTENT(IN) :: PTNUDGING  ! Nudging time scale
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUM,PVM,PWM,PTHM 
                                          ! wind, potential temperature and
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM !  moist variables at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM
                                          ! Large Scale fields
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS,PRVS,PRWS,PRTHS
                                          ! wind, potential temperature and
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS ! moist variables at time t+1
!
!*       0.2   Declarations of local variables
!
REAL :: ZINVTAU ! inverse of nudging time scale
!
!----------------------------------------------------------------------------
!
!
ZINVTAU=1./PTNUDGING
!
!*        1.   NUGDGING TOWARDS LS FIELDS
!              --------------------------
!
PRTHS(:,:,:) = PRTHS(:,:,:) - ZINVTAU* PRHODJ(:,:,:)* (PTHM(:,:,:)-PLSTHM(:,:,:))
PRUS(:,:,:) = PRUS(:,:,:) - ZINVTAU* PRHODJ(:,:,:)* (PUM(:,:,:)-PLSUM(:,:,:))
PRVS(:,:,:) = PRVS(:,:,:) - ZINVTAU* PRHODJ(:,:,:)* (PVM(:,:,:)-PLSVM(:,:,:))
PRWS(:,:,:) = PRWS(:,:,:) - ZINVTAU* PRHODJ(:,:,:)* (PWM(:,:,:)-PLSWM(:,:,:))
IF (OUSERV) &
  PRRS(:,:,:,1) = PRRS(:,:,:,1) - ZINVTAU* PRHODJ(:,:,:)* (PRM(:,:,:,1)-PLSRVM(:,:,:))
!
!
!*       2.     BUDGET CALLS
!   	        ------------
!
IF (LBUDGET_U)   CALL BUDGET (PRUS,1,'NUD_BU_RU')
IF (LBUDGET_V)   CALL BUDGET (PRVS,2,'NUD_BU_RV')
IF (LBUDGET_W)   CALL BUDGET (PRWS,3,'NUD_BU_RW')
IF (LBUDGET_TH)  CALL BUDGET (PRTHS,4,'NUD_BU_RTH')
IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'NUD_BU_RRV')
!
END SUBROUTINE NUDGING
