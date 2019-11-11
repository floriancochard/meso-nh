!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_SLOW_TERMS
!     ###################### 
INTERFACE
      SUBROUTINE SLOW_TERMS ( KSPLITR, PTSTEP, KMI, HSUBG_AUCV,            &
                              PZZ, PRHODJ, PRHODREF, PCLDFR,               &
                              PTHT, PRVT, PRCT, PRRT, PPABST,              &
                              PTHS, PRVS, PRCS, PRRS, PINPRR,              &
                              PINPRR3D, PEVAP3D                            )
!
!
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI     ! Model index
CHARACTER(LEN=4),         INTENT(IN)   :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),     INTENT(IN)  :: PCLDFR  ! Cloud fraction
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRR  ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PEVAP3D ! Rain evap profil
!
END SUBROUTINE SLOW_TERMS
END INTERFACE
END MODULE MODI_SLOW_TERMS 
!
!
!     ######################################################################
!OPTION! -Ni
      SUBROUTINE SLOW_TERMS ( KSPLITR, PTSTEP, KMI, HSUBG_AUCV,            &
                              PZZ, PRHODJ, PRHODREF, PCLDFR,               &
                              PTHT, PRVT, PRCT, PRRT, PPABST,              &
                              PTHS, PRVS, PRCS, PRRS, PINPRR,              &
                              PINPRR3D, PEVAP3D                            )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the slow microphysical sources
!!    which can be computed explicitly 
!!
!!
!!**  METHOD
!!    ------
!!      The autoconversion computation follows Kessler (1969). 
!!      The sedimentation rate is computed with a time spliting technique and 
!!    an upstream scheme, written as a difference of non-advective fluxes. This
!!    source term is added to the future instant ( split-implicit process ).
!!      The others microphysical processes are evaluated at the central instant 
!!    (split-explicit process ): autoconversion, accretion and rain evaporation.
!!      These last 3 terms are bounded in order not to create negative values 
!!    for the water species at the future instant.  
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!      Module MODD_CLOUDPAR
!!          XCEXVT               ! constant in the rain drop fall velocity
!!          XC1RC, XC2RC         ! constants for autoconversion
!!          XCEXRA, XCRA         ! constants for accretion     
!!          XCEXRE, XC1RE, XC2RE ! constants for rain evaporation 
!!          XCEXRS, XCRS         ! constants for rain sedimentation
!!          XDIVA                ! vapor diffusivity in air
!!          XTHCO                ! thermal conductivity 
!!      Module MODD_CST     
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor 
!!                             !  pressure  function 
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH        
!!                        .FALSE. = no budget of RTH
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV 
!!                        .FALSE. = no budget of RRV 
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC 
!!                        .FALSE. = no budget of RRC 
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR 
!!                        .FALSE. = no budget of RRR 
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 and Book2 of documentation ( routine SLOW_TERMS )
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    22/12/94 
!!      Modifications: June 8, 1995 ( J.Stein )
!!                                   Cleaning to improve efficienty and clarity
!!                                  in agreement with the MESO-NH coding norm
!!                     16/10/95 (J. Stein)     change the budget calls
!!                     16/10/96 (J. Stein)     use PEXNT instead of PEXNREF 
!!                     19/12/96 (J.-P. Pinty)  update the budget calls
!!                     04/02/97 (J.Viviand)    convert precipitation rate in m/s
!!                                             & debug accumulated precipitation
!!                     01/09/97 (V. Masson & J. Stein) optimization: reduction
!!                                             of the number of exponentiations
!!                     14/09/97 (V. Masson) removes low rr non-physical values
!!                     06/11/02 (V. Masson) update the budget calls
!!     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CLOUDPAR  
USE MODD_CST
USE MODD_CONF
USE MODD_BUDGET
!
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV
                                        ! Kind of Subgrid autoconversion method
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),     INTENT(IN)  :: PCLDFR  ! Cloud fraction
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRR  ! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PEVAP3D ! Rain evap profil
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation 
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: IKB           !  Define the domain where is 
INTEGER :: IKE           !  the microphysical sources have to be computed
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) &
                               :: ZT,ZW,ZW1,ZW2,ZW3,ZEXNT,ZDZZ  ! Work arrays 
LOGICAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: G3D
!
INTEGER , DIMENSION(:), ALLOCATABLE  :: I1,I2       ! index for packed array
INTEGER                              :: JI,JJ,IC,JL ! loop control for packed array
!  
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS AND EXNER FUNCTION
!   	        ------------------------------------------
!
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
! compute Delta z at the mass point
DO JK = 1,IKE
  ZDZZ(:,:,JK) =  PZZ(:,:,JK+1)-PZZ(:,:,JK) 
END DO

!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
!
!*       2.1    time splitting loop initialization        
!
ZTSPLITR = PTSTEP / FLOAT(KSPLITR)       ! Small time step
! 
ZW1(:,:,:) = PRRS(:,:,:) * PTSTEP
ZW2(:,:,:) = 0.
ZW3(:,:,:) = 0.
!
G3D(:,:,:)=.FALSE.
DO JK=IKE,1,-1
  WHERE( ZW1(:,:,JK)>0. .OR. G3D(:,:,JK+1) )
    G3D(:,:,JK)=.TRUE.
  END WHERE
END DO
!
WHERE (G3D(:,:,:))
  ZW3(:,:,:) = PRHODREF(:,:,:) ** (XCEXRS-XCEXVT)
END WHERE
!
WHERE (ZW1(:,:,IKE+1)>0.)
  ZW2(:,:,IKE+1) =   XCRS                         &
                  * ZW1(:,:,IKE+1) ** XCEXRS      &
                  * PRHODREF(:,:,IKE+1) ** (XCEXRS-XCEXVT)
END WHERE
!
!
!*       2.2    small time step integration
!               ---------------------------
!
ALLOCATE (I1(SIZE(ZW1)))
ALLOCATE (I2(SIZE(ZW1)))
!
DO JN=1,KSPLITR
!
!*       2.2.1    test where computations should be made and pack arrays
!                 ------------------------------------------------------
!

  DO JK=IKE,IKB,-1
    !
    IC = 0
    DO JJ = 1,SIZE( ZW1,2)
      DO JI = 1,SIZE( ZW1,1)
        IF ( ( ZW1(JI,JJ,JK+1)>0. ) .AND. ( ZW1(JI,JJ,JK)>0. ) )  THEN
!!$       IF ( (ZW1(JI,JJ,JK)+ZW1(JI,JJ,JK+1)>0.) )  THEN
          IC = IC +1
          I1(IC) = JI
          I2(IC) = JJ
        ENDIF
      ENDDO
    ENDDO
!
!*       2.2.2    compute the rain flux
!                 ---------------------
!
!
    DO JL=1,IC
      ZW2(I1(JL),I2(JL),JK) =   XCRS                         &
                   * ZW1(I1(JL),I2(JL),JK) ** XCEXRS         &
                   * ZW3(I1(JL),I2(JL),JK)
!
!
!*       2.2.3    update the rain tendency with the small timestep
!                 ------------------------------------------------
! 
      ZW1(I1(JL),I2(JL),JK)  = ZW1(I1(JL),I2(JL),JK) + ZTSPLITR *       &
                 ( ZW2(I1(JL),I2(JL),JK+1)-ZW2(I1(JL),I2(JL),JK) ) /    &
                 ( ZDZZ(I1(JL),I2(JL),JK) * PRHODREF(I1(JL),I2(JL),JK) )
    ENDDO
!
  ENDDO
!
!*       2.2.4    compute the explicit accumulated precipitations
!                 -----------------------------------------------
! 
  IF( JN==1 ) THEN
    PINPRR(:,:) = ZW2(:,:,IKB)/XRHOLW                           ! in m/s
    PINPRR3D(:,:,:) = ZW2(:,:,:)/XRHOLW                           ! in m/s
  END IF
!
END DO
!
DEALLOCATE (I1)
DEALLOCATE (I2)
!
!*       2.4     update the rain tendency
!
PRRS(:,:,:) = ZW1(:,:,:) / PTSTEP
! 
!*       2.5     budget storage
!
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:),8,'SEDI_BU_RRR')
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTES THE ACCRETION SOURCE
!   	        -----------------------------
!
!
!*       3.1     compute the accretion and update the tendencies
!
WHERE ( (PRCT(:,:,:)>0.0) .AND. (PRRT(:,:,:)>0.0) &
                          .AND. (PRCS(:,:,:)>0.0) )
  ZW(:,:,:) = XCRA * PRCT(:,:,:)                        &
             * PRRT(:,:,:) ** XCEXRA                    &
             * PRHODREF(:,:,:) ** (XCEXRA - XCEXVT)
  ZW(:,:,:) = MIN ( ZW(:,:,:) , PRCS(:,:,:) )
!
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
END WHERE
!
!*       3.2     budget storage
!
IF (LBUDGET_RC) CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:),7,'ACCR_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:),8,'ACCR_BU_RRR')
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE AUTOCONVERSION SOURCE
!               ----------------------------------
!
!
!*       4.1     compute the autoconversion and update the tendencies
!
IF ( HSUBG_AUCV == 'CLFR' ) THEN
 WHERE ( (PRCT(:,:,:)>0.0) .AND. (PRCS(:,:,:)>0.0) .AND. (PCLDFR(:,:,:)>0.0) )
  ZW(:,:,:) = XC1RC * MAX(PRCT(:,:,:) /(PCLDFR(:,:,:)) - XC2RC / PRHODREF(:,:,:),0.)
  ZW(:,:,:) = MIN ( (ZW(:,:,:)* PCLDFR(:,:,:)), PRCS(:,:,:) )
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
 END WHERE
ELSE
 WHERE ( (PRCT(:,:,:)>0.0) .AND. (PRCS(:,:,:)>0.0) )
  ZW(:,:,:) = XC1RC * MAX ( PRCT(:,:,:) - XC2RC / PRHODREF(:,:,:), 0. )
  ZW(:,:,:) = MIN ( ZW(:,:,:) , PRCS(:,:,:) )
!
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
 END WHERE
END IF
!
!*       4.2     budget storage
!
IF (LBUDGET_RC) CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:),7,'AUTO_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:),8,'AUTO_BU_RRR')
!
!-------------------------------------------------------------------------------
!
!*       5.     COMPUTES THE RAIN EVAPORATION (RE) SOURCE
!   	        -----------------------------------------
!
PEVAP3D(:,:,:)=0.
WHERE ( (PRRT(:,:,:)>0.0) .AND. (PRCT(:,:,:)==0.0) ) 
!
!*       5.1    compute the Exner function
!
  ZEXNT(:,:,:) = (PPABST(:,:,:)/XP00)**(XRD/XCPD)
!
!*       5.2    compute the temperature
!
  ZT(:,:,:) = PTHT(:,:,:) * ZEXNT(:,:,:)
!
!*       5.3    compute the saturation vapor pressure
!
  ZW1(:,:,:) = EXP( XALPW - XBETAW/ZT(:,:,:) - XGAMW*ALOG(ZT(:,:,:) ) )
!
!*       5.4    compute the saturation mixing ratio
!
  ZW2(:,:,:) = (XMV/XMD) * ZW1(:,:,:)                &
                         / ( PPABST(:,:,:) - ZW1(:,:,:) )  
!
!*       5.5    compute the latent heat of vaporization
!
  ZW3(:,:,:) = XLVTT + ( XCPV - XCL ) * ( ZT(:,:,:) -XTT )
!
!*       5.6    compute the source
!
  ZW(:,:,:) = MAX( 1. - PRVT(:,:,:)/ZW2(:,:,:) , 0.0 ) *         & 
    ( XC1RE * SQRT( PRRT(:,:,:)/PRHODREF(:,:,:) )                & 
     +XC2RE * PRRT(:,:,:)**XCEXRE                                & 
            * PRHODREF(:,:,:)**(XCEXRE-0.5*XCEXVT-1.)            &
    ) /                                                          &
    ( XRV*ZT(:,:,:)/(ZW1(:,:,:)*XDIVA)                           & 
     +ZW3(:,:,:) ** 2 / (XTHCO*XRV*ZT(:,:,:)**2)                 &
    )     
  ZW(:,:,:) = MIN ( ZW(:,:,:) , PRRS(:,:,:) )
!
!*       5.7     update the total tendencies
!
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) + ZW(:,:,:)
!
  PTHS (:,:,:) = PTHS (:,:,:)  - ZW(:,:,:) * ZW3(:,:,:) /  &
    ( ZEXNT(:,:,:) * (  XCPD + XCPV * PRVT(:,:,:)          &
                      + XCL * (PRCT(:,:,:) + PRRT(:,:,:)) ) )
  PEVAP3D(:,:,:)=ZW(:,:,:)
END WHERE
!
!*       5.8     budget storage
!
IF (LBUDGET_RV) CALL BUDGET (PRVS(:,:,:)*PRHODJ(:,:,:),6,'REVA_BU_RRV')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:),8,'REVA_BU_RRR')
IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:)*PRHODJ(:,:,:),4,'REVA_BU_RTH')
!
!-------------------------------------------------------------------------------
!
!*       6.     REMOVES NON-PHYSICAL LOW VALUES
!   	        -------------------------------
!
WHERE (PRRS(:,:,:)<1.E-16)
  PRRS(:,:,:)=0.
END WHERE
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE SLOW_TERMS 
