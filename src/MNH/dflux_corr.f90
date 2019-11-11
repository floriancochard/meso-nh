!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_DFLUX_CORR 
!     ##########################
!
INTERFACE
!
!     ######################################################################
      SUBROUTINE DFLUX_CORR   ( HLBCX, HLBCY, PTSTEP, PMIN,                &
                                PRHODJ, PAM, PAT, PRUCT, PRVCT, PRWCT,     &
                                PFX, PFY, PFZ)
!     ######################################################################
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
REAL,                     INTENT(IN)    :: PTSTEP ! Double time step
REAL,                     INTENT(IN)    :: PMIN
                                                  ! Absolute minimum variable
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
                                                  ! (Rho) dry *jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PAM 
                                                  ! Variable at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PAT
                                                  ! Variable at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT  ! Contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT  ! components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT  ! of momentum
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PFX, PFY, PFZ
                                                  ! Flux components
!
END SUBROUTINE DFLUX_CORR 
!
END INTERFACE
!
END MODULE MODI_DFLUX_CORR  
!
!
!
!     ######################################################################
      SUBROUTINE DFLUX_CORR   ( HLBCX, HLBCY,  PTSTEP, PMIN,               &
                                PRHODJ, PAM, PAT, PRUCT, PRVCT, PRWCT,     &
                                PFX, PFY, PFZ)
!     ######################################################################
!
!!****  *DFLUX_CORR* - calculates the advective tendencies fluxes by means of
!!                      the Directional Flux-Corrected Transport and  
!!                      the Flux-Corrected Transport advection schemes
!!
!!    PURPOSE
!!    -------
!!
!!    The purpose of the routine is to calculate the advection of a scalar.
!!    A centred advection scheme is used (leapfrog). Two corrections of
!!    the fluxes are is applied (DFCT and FCT) to insure that the total
!!    resulting scheme is positive definite.
!!    The advection scheme is second-order on time and on space.
!
!!**  METHOD
!!    ------
!!    
!!    First, the advective flux is calculated. Second the flux is corrected
!!    using the Directional-Flux-Corrected Transport method.
!!    Second, the advective flux is calculated. Eventually the flux is corrected
!!    using the Flux-Corrected Transport method. This method implies the
!!    calculation of one limiting factors: BETAOUT. The
!!    first factor insures that the calculated flux is less than its respective
!!    analytical value (nonoscillatory condition). 
!!
!!    EXTERNAL
!!    --------
!!    GET_DIM_EXT_ll : get extended sub-domain sizes
!!    ADD3DFIELD_ll  : add a field to 3D-list
!!    UPDATE_HALO_ll : update internal halos 
!!    UPDATE_BOUNDARIES_ll : update external boundaries
!!    LWEAST_ll,LEAST_ll,LNORTH_ll,LSOUTH_ll : position functions
!!    MXM,MYM,MZM   : Shuman functions (mean operators)
!!    DXF,DYF,DZF   : Shuman functions (finite difference operators) 
!!    CLEANLIST_ll  : deaalocate a list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    JPHEXT, JPVEXT
!!    MODD CONF:CCONF
!!    
!!    REFERENCE
!!    ---------
!!    Book1 of documentation (FCT scheme) 
!!
!!    AUTHOR
!!    ------
!!    J.-P. Lafore     *Meteo-France* 
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    original     27/03/98 
!!    V. Masson    24/11/97   removes the DO loops
!!    P. Jabouille 24/09/98   parallelize the code
!!    J. Stein     05/04/99 : bug for the case PMIN /= 0 + lbc
!!    JP Pinty &   12/10/98 : Vectorization of the first loops
!!      J Escobar 
!!    J. Stein &    20/03/01 : bug for the open case at the boundary
!!      P. Jabouille 
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_CONF
USE MODD_PARAMETERS
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*      0.1  DECLARATIONS OF ARGUMENTS
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
!
REAL,                     INTENT(IN)    :: PTSTEP
                                                  ! Double Time step
REAL,                     INTENT(IN)    :: PMIN
                                                  ! Absolute minimum variable
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
                                                  ! (Rho) dry *jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PAM 
                                                  ! Variable at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PAT
                                                  ! Variable at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT     ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT     !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT     ! of momentum
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PFX, PFY, PFZ
                                                  ! Flux component 
!
!*      0.2  DECLARATIONS OF LOCAL VARIABLES
!
!
INTEGER:: IIU,IJU,IKU                           ! Size array in the x, y,
                                                ! and z directions 
INTEGER:: JI,JJ,JK                              ! Loop index in the x, y,
                                                ! and z directions 
REAL   :: ZEPSILON                              ! Variable to ensure that the
                                                ! limiting factor is zero  
!
!
REAL,DIMENSION(SIZE(PAT,1),SIZE(PAT,2),SIZE(PAT,3)):: ZFOUT
                                                ! The outgoing flux of the grid cell 
                                                ! (located at mass point) 
REAL,DIMENSION(SIZE(PAT,1),SIZE(PAT,2),SIZE(PAT,3)):: ZBETAOUT
                                                ! The outgoing limiting factor 
INTEGER                :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!
!*     0.3   PROLOGUE 
!
NULLIFY(TZFIELDS_ll)
!
CALL GET_DIM_EXT_ll    ('B',IIU,IJU)
IKU=SIZE(PRUCT,3)
!
ZEPSILON=1.0E-15
!
!------------------------------------------------------------------------------
!
!
!*       1. First limitation on a directional base
!           --------------------------------------
!
ZBETAOUT(:,:,:) = -PRHODJ(:,:,:)*(PAM(:,:,:)-PMIN)/PTSTEP ! First limiter
!
!*      1.1 X-direction
!
ZFOUT(2:IIU,:,:) = -ZBETAOUT(1:IIU-1,:,:) ! Second limiter
ZFOUT(1,:,:)     = 0.0
!
PFX(:,:,:) = PRUCT(:,:,:) * MXM (PAT(:,:,:)) 
PFX(:,:,:) = (0.5+SIGN(0.5,PRUCT(:,:,:)))*MIN( PFX(:,:,:),ZFOUT(:,:,:) )  &
            +(0.5-SIGN(0.5,PRUCT(:,:,:)))*MAX( PFX(:,:,:),ZBETAOUT(:,:,:) )
!
!*      1.2 Y-direction
!
ZFOUT(:,2:IJU,:) = -ZBETAOUT(:,1:IJU-1,:) ! Second limiter
ZFOUT(:,1,:)     = 0.0
!
PFY(:,:,:) = PRVCT(:,:,:) * MYM (PAT(:,:,:))  
PFY(:,:,:) = (0.5+SIGN(0.5,PRVCT(:,:,:)))*MIN( PFY(:,:,:),ZFOUT(:,:,:) )  &
            +(0.5-SIGN(0.5,PRVCT(:,:,:)))*MAX( PFY(:,:,:),ZBETAOUT(:,:,:) )
!
!*      1.3 Z-direction
!
ZFOUT(:,:,2:IKU) = -ZBETAOUT(:,:,1:IKU-1) ! Second limiter
ZFOUT(:,:,1)     = 0.0
!
PFZ(:,:,:) = PRWCT(:,:,:) * MZM (1,IKU,1,PAT(:,:,:))  
PFZ(:,:,:) = (0.5+SIGN(0.5,PRWCT(:,:,:)))*MIN( PFZ(:,:,:),ZFOUT(:,:,:) )  &
            +(0.5-SIGN(0.5,PRWCT(:,:,:)))*MAX( PFZ(:,:,:),ZBETAOUT(:,:,:) )
!
!
!------------------------------------------------------------------------------
!*      3.  Flux-OUT calculation 
!           ---------------------
!
DO JK=2,IKU-1
  DO JJ=2,IJU-1
    DO JI=2,IIU-1
      ZFOUT(JI,JJ,JK) =  MAX(0.,PFX(JI+1,JJ,JK))       &
                       - MIN(0.,PFX(JI,  JJ,JK))       &
                       + MAX(0.,PFY(JI,JJ+1,JK))       &
                       - MIN(0.,PFY(JI,JJ  ,JK))       &
                       + MAX(0.,PFZ(JI,JJ,JK+1))       &
                       - MIN(0.,PFZ(JI,JJ,JK  ))
    END DO
  END DO
END DO 
!
!           
!------------------------------------------------------------------------------
!*      4.  BETAOUT calculation 
!           -------------------
ZBETAOUT(:,:,:) =(PAM(:,:,:)-PMIN)/                            &
                      (PTSTEP*ZFOUT(:,:,:)/PRHODJ(:,:,:)+ZEPSILON) 
!           
ZBETAOUT(:,:,1)   = 1.      ! no limitation outside the physical domain
ZBETAOUT(:,:,IKU) = 1.      ! because no velocity is available                  
!           
! 
! Update halo and apply possible cyclic boundary conditions
!
!!$IF(NHALO == 1 .OR. HLBCX(1)=='CYCL' .OR. HLBCY(1)=='CYCL') THEN
IF(HLBCX(1)=='CYCL' .OR. HLBCY(1)=='CYCL') THEN
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZBETAOUT)
!!$  IF(NHALO == 1) THEN
    CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
!!$  ELSE
!!$    IF(HLBCX(1)=='CYCL') CALL UPDATE_BOUNDARIES_ll('XX',TZFIELDS_ll,IINFO_ll)
!!$    IF(HLBCY(1)=='CYCL') CALL UPDATE_BOUNDARIES_ll('YY',TZFIELDS_ll,IINFO_ll)
!!$  END IF
  CALL CLEANLIST_ll(TZFIELDS_ll)
ENDIF
!
!
IF (HLBCX(1)/='CYCL') THEN
  IF (LWEST_ll( ))  ZBETAOUT(1,:,:) = 1.      ! no limitation outside the physical domain
  IF (LEAST_ll( ))  ZBETAOUT(IIU,:,:) = 1.    ! because no velocity is available
END IF
!
IF (HLBCY(1)/='CYCL') THEN
  IF (LSOUTH_ll( )) ZBETAOUT(:,1,:) = 1.      ! no limitation outside the physical domain
  IF (LNORTH_ll( )) ZBETAOUT(:,IJU,:) = 1.    ! because no velocity is available 
END IF
!
!            
!------------------------------------------------------------------------------
!*      4.  Second Flux limitation 
!           ----------------------
!
!
                                          ! x-component
                                          ! 
ZFOUT(2:IIU,:,:) = ZBETAOUT(1:IIU-1,:,:)
PFX(:,:,:) =  MIN(1., ZFOUT(:,:,:))    * MAX(0.,PFX(:,:,:)) &
            + MIN(1., ZBETAOUT(:,:,:)) * MIN(0.,PFX(:,:,:))
                                          ! y-component
                                          ! 
ZFOUT(:,2:IJU,:) = ZBETAOUT(:,1:IJU-1,:)
PFY(:,:,:) =  MIN(1., ZFOUT(:,:,:))    * MAX(0.,PFY(:,:,:)) &
            + MIN(1., ZBETAOUT(:,:,:)) * MIN(0.,PFY(:,:,:))
                                          ! z-component
                                          ! 
ZFOUT(:,:,2:IKU) = ZBETAOUT(:,:,1:IKU-1)
PFZ(:,:,:) =  MIN(1., ZFOUT(:,:,:))    * MAX(0.,PFZ(:,:,:)) &
            + MIN(1., ZBETAOUT(:,:,:)) * MIN(0.,PFZ(:,:,:))
!
!------------------------------------------------------------------------------
!
!*      5.  Boundary conditions for the flux cyclic case
!           --------------------------------------------
!
                                                      ! x-direction
IF (HLBCX(1)=='CYCL') THEN
  CALL ADD3DFIELD_ll(TZFIELDS_ll, PFX)
  CALL UPDATE_BOUNDARIES_ll('XX',TZFIELDS_ll, IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll) 
ENDIF
                                                      ! y-direction
IF (HLBCY(1)=='CYCL') THEN
  CALL ADD3DFIELD_ll(TZFIELDS_ll, PFY)
  CALL UPDATE_BOUNDARIES_ll('YY',TZFIELDS_ll, IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
ENDIF
!
!
!------------------------------------------------------------------------------
END SUBROUTINE DFLUX_CORR 
