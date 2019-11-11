!MNH_LIC Copyright 1999-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_TWO_WAY
!     ###################
!
INTERFACE 
!
      SUBROUTINE TWO_WAY   (     KRR,KSV,KTCOUNT,PRHODJ,KMI,PTSTEP,                &
                            PUM ,PVM, PWM, PTHM, PRM, PTKEM, PSVM,                 &
                            PRUS,PRVS,PRWS,PRTHS,PRRS,PRTKES,PRSVS,                &
                            PINPRC,PINPRR,PINPRS,PINPRG,PINPRH,PPRCONV,PPRSCONV,   &
                            PDIRFLASWD,PSCAFLASWD,PDIRSRFSWD,OMASKkids             )
! 
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of Scalar Variables
INTEGER,                  INTENT(IN)  :: KTCOUNT ! Temporal loop COUNTer
                                                 ! (=1 at the segment beginning)
INTEGER,                  INTENT(IN)  :: KMI     ! Model index     
!
REAL,                     INTENT(IN)  :: PTSTEP  ! Timestep duration
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PUM, PVM, PWM  ! Variables at t-dt 
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRM, PSVM
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS, PRSVS              !  terms
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC,PINPRR,PINPRS,PINPRG,PINPRH &
                                          ,PPRCONV,PPRSCONV     !  precipitating variables
LOGICAL, DIMENSION(:,:), INTENT(INOUT)  :: OMASKkids ! true where kids exist
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRFLASWD,PSCAFLASWD   ! Short wave radiation
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRSRFSWD

!
END SUBROUTINE TWO_WAY 
!
END INTERFACE
!
END MODULE MODI_TWO_WAY 
!
!     ########################################################################
      SUBROUTINE TWO_WAY   (     KRR,KSV,KTCOUNT,PRHODJ,KMI,PTSTEP,                &
                            PUM ,PVM, PWM, PTHM, PRM, PTKEM, PSVM,                 &
                            PRUS,PRVS,PRWS,PRTHS,PRRS,PRTKES,PRSVS,                &
                            PINPRC,PINPRR,PINPRS,PINPRG,PINPRH,PPRCONV,PPRSCONV,   &
                            PDIRFLASWD,PSCAFLASWD,PDIRSRFSWD,OMASKkids             )
!     ########################################################################
!
!!****  *TWO_WAY* - relaxation toward the fine-mesh model result
!!
!!    PURPOSE
!!    -------
!!      The purpose of TWO_WAY  is to select the KID model and to call the
!!    subroutine TWO_WAY$n which effectively performs the relaxation toward the 
!!    fine-mesh model result. 
!
!
!!**  METHOD
!!    ------
!!      
!!       We choose the right KID model and CALL the appropriate TWO_WAY$n
!!
!!    EXTERNAL
!!    --------
!!       TWO_WAYY1,...,8 : performs the relaxation
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODULE MODD_:   
!!
!!    REFERENCE
!!    ---------
!!    
!!
!!    AUTHOR
!!    ------
!!    J. Stein  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     8/4/99 
!!      N. Asencio   18/07/05  Add the surface parameters : precipitating
!!                             hydrometeors, the Short and Long Wave 
!!                              + MASKkids array
!!                   20/05/06 Remove EPS
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                 
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODD_CONF
USE MODD_NESTING
USE MODD_BUDGET

USE MODI_BUDGET
!
USE MODI_TWO_WAY_n
USE MODE_MODELN_HANDLER
!
IMPLICIT NONE
!
!
!
!*       0.1   declarations of arguments 
! 
!
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of Scalar Variables
INTEGER,                  INTENT(IN)  :: KTCOUNT ! Temporal loop COUNTer
                                                 ! (=1 at the segment beginning)
INTEGER,                  INTENT(IN)  :: KMI     ! Model index     
!
REAL,                     INTENT(IN)  :: PTSTEP  ! Timestep duration
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PUM, PVM, PWM  ! Variables at t-dt 
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRM, PSVM
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS, PRSVS              !  terms
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC,PINPRR,PINPRS,PINPRG,PINPRH &
                                          ,PPRCONV,PPRSCONV     !  precipitating variables
LOGICAL, DIMENSION(:,:), INTENT(INOUT)  :: OMASKkids ! true where kids exist
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRFLASWD,PSCAFLASWD   ! Short wave radiation
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRSRFSWD
!
!*       0.2   declarations of local variables
!
INTEGER :: JKID        ! loop index to look for the KID models
INTEGER :: JSV,JRR     ! Loop index for scalar and moist variables
!
!-------------------------------------------------------------------------------
!
!*       1.    CALL THE RIGHT TWO_WAY$n
!              ------------------------
!
DO JKID = KMI+1,NMODEL  ! min value of the possible kids
  IF (KMI == NDAD(JKID) .AND. (XWAY(JKID) == 2. ) &
   .AND. (CCONF == 'RESTA' .OR. (CCONF == 'START' .AND. KTCOUNT /= 1))) THEN
    CALL GOTO_MODEL(JKID)
    CALL TWO_WAY_n (KRR,KSV,KTCOUNT,PRHODJ,KMI,PTSTEP,                     &
                    PUM ,PVM, PWM, PTHM, PRM, PTKEM, PSVM,                 &
                    PRUS,PRVS,PRWS,PRTHS,PRRS,PRTKES,PRSVS,                &
                    PINPRC,PINPRR,PINPRS,PINPRG,PINPRH,PPRCONV,PPRSCONV,   &
                    PDIRFLASWD,PSCAFLASWD,PDIRSRFSWD,OMASKkids             )
  END IF
END DO
CALL GOTO_MODEL(KMI)
!
!*       2.    BUDGET COMPUTATION
!              ------------------
!
IF (LBUDGET_U)  CALL BUDGET (PRUS,1,'NEST_BU_RU')
IF (LBUDGET_V)  CALL BUDGET (PRVS,2,'NEST_BU_RV')
IF (LBUDGET_W)  CALL BUDGET (PRWS,3,'NEST_BU_RW')
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'NEST_BU_RTH')
DO JRR=1,KRR
  IF (JRR==1 .AND. LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,JRR),6,'NEST_BU_RRV')
  IF (JRR==2 .AND. LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,JRR),7,'NEST_BU_RRC')
  IF (JRR==3 .AND. LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,JRR),8,'NEST_BU_RRR')
  IF (JRR==4 .AND. LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,JRR),9,'NEST_BU_RRI')
  IF (JRR==5 .AND. LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,JRR),10,'NEST_BU_RRS')
  IF (JRR==6 .AND. LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,JRR),11,'NEST_BU_RRG')
  IF (JRR==7 .AND. LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,JRR),12,'NEST_BU_RRH')
ENDDO
DO JSV=1,KSV
  IF (LBUDGET_SV)              CALL BUDGET (PRSVS(:,:,:,JSV),12+JSV,'NEST_BU_RSV')
END DO
!------------------------------------------------------------------------------
!
END SUBROUTINE TWO_WAY
