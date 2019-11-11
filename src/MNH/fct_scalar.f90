!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/06/06 12:01:31
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_FCT_SCALAR
!     ######################
!
INTERFACE
  SUBROUTINE FCT_SCALAR  (HLBCX, HLBCY, KSV,             &
                          PTSTEP, PRHODJ, PSVM, PSVT,    &
                          PRUCT, PRVCT, PRWCT, PRSVS     )
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Double Time step 
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM
                                                  ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT
                                                  ! Contravariant component of
                                                  ! momentum 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                  ! Sources terms 
!
END SUBROUTINE FCT_SCALAR 
!
END INTERFACE
!
END MODULE MODI_FCT_SCALAR 
!
!
!
! #####################################################################
  SUBROUTINE FCT_SCALAR  (HLBCX, HLBCY, KSV,                          &
                          PTSTEP, PRHODJ, PSVM, PSVT,                 &
                          PRUCT, PRVCT, PRWCT, PRSVS                  )              
! #####################################################################
!
!!****  *FCT_SCALAR * - routine to call the Flux-Corrected Transport Scalar 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to call the FLUX-CORRected routine
!!    for scalar variables (tracers).
!!
!!**  METHOD
!!    ------
!!     The Flux-Corrected Transport method correct the fluxes using a limiting
!!    factor. This method ensures that the advection scheme is definite 
!!    positive. A minimum value of the scalar (MIN) equal to 0 is used for
!!    the positiveness of the scheme.
!!
!!    EXTERNAL
!!    --------
!!     Functions MXM,MYM,MZM : computes the averages along three directins
!!     Functions DXF,DYF,DZF  : computes the finite differences  
!!     Subroutine FLUX_CORR : corrects the advective fluxes 
!!     Subroutine BUDGET    : stores the sources for budget purposes
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_BUDGET : LBU_R* ( individual budget switches) 
!!                    CBUTYPE, NBUMOD
!!    REFERENCE
!!    ---------
!!      Book1 and book2 ( routine ADVECTION )
!!
!!    AUTHOR
!!    ------
!!      J. Vila & JP. Lafore    *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/10/94 
!!      Stein       27/06/96  add the budgets
!!      Pinty       20/12/96  update the budgets
!!      Lafore      27/03/98  call to DFLUX_CORR
!!      Lafore      01/04/98  FCT only on total water (rv+rc+ri) and
!!                            precipitating hydrometeors,
!!                            remove 4D flux local arrays
!!      Stein       20/04/99  remove KMI from the list of argument of DFLUX_CORR
!!      Jabouille   22/06/01  use XSVMIN
!!      Masson      07/11/02  update the budgets
!!      Lac         24/04/06  Split scalar and passive tracer routines
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_NSV, ONLY : XSVMIN
USE MODD_GRID_n
USE MODI_SHUMAN
USE MODI_DFLUX_CORR
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step 
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM
                                                  ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT
                                                  ! Contravariant component of
                                                  ! momentum 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                  ! Sources terms 
!
!
!*       0.2   declarations of local variables
!
INTEGER                                  :: JSV
                                                  ! Loop index 
!
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))  &
                     :: ZFX  ,ZFY  ,ZFZ    ! Advective flux components for each
INTEGER :: IKU                     
!
!-------------------------------------------------------------------------------
!
IKU=SIZE(XZHAT)
!*       1.   FLUX-CORRECTED TRANSPORT ADVECTION SCHEME for the HSV group
!
!
!
!
  DO JSV = 1, KSV
    CALL DFLUX_CORR (HLBCX,HLBCY,PTSTEP,XSVMIN(JSV),PRHODJ,                &
                     PSVM(:,:,:,JSV),PSVT(:,:,:,JSV),PRUCT,PRVCT,PRWCT,    &
                      ZFX(:,:,:),ZFY(:,:,:),ZFZ(:,:,:)                     )
!
    PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) - DXF(ZFX(:,:,:))        
    IF (LBUDGET_SV)                               &
                            CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVX_BU_RSV')
!
    PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) - DYF(ZFY(:,:,:))        
    IF (LBUDGET_SV)                               &
                            CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVY_BU_RSV')
!
    PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) - DZF(1,IKU,1,ZFZ(:,:,:)) 
    IF (LBUDGET_SV)                               &
                            CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'ADVZ_BU_RSV')
  END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FCT_SCALAR
