!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_GRAVITY    
!     ###################
!
INTERFACE
!
      SUBROUTINE GRAVITY    ( KRR,KRRL, KRRI, PTHT, PRT,                   &
                              PRHODJ, PTHVREF, PRWS                        )
!
INTEGER,                  INTENT(IN)    :: KRR  ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT            !     at
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !  time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRWS ! Sources of Momentum
!
END SUBROUTINE GRAVITY
!
END INTERFACE
!
END MODULE MODI_GRAVITY 
!     ######################################################################
      SUBROUTINE GRAVITY( KRR,KRRL, KRRI, PTHT, PRT,          &
                          PRHODJ, PTHVREF, PRWS               )
!     ######################################################################
!
!!****  *GRAVITY * - routine to compute the curvature, coriolis and gravity terms
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the gravity on the vertical 
!!    component of the momentum.
!!
!!     
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      MXM,MYM,MZM : Shuman functions (mean operators)
!!      MXF,MYF,MZF : Shuman functions (mean operators)
!!      GZ_M_W      : projection along the vertical direction of the gradient
!!                    vector. It acts on a field localized in mass point and
!!                    the result is localized in w point.
!!      BUDGET      : Stores the different budget components
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONF     : contains configuration variables for all models
!!           LTHINSHELL  = .TRUE.  if the THINSHELL approximation is made
!!           LCARTESIAN  = .TRUE.  if the CARTESIAN approximation is made
!!      Module MODD_CST      : 
!!           XG           gravity acceleration
!!           XRADIUS      Earth radius
!!           XRD,XRV      Gas constant for dry air and wator vapor
!!           XCPD,XCPV    Cp for dry air and wator vapor
!!           XCL,XCI      C (calorific capacity) for liquid and solid water
!!      Module MODD_DYN      : contains the parameters for the dynamics
!!           LCORIO      = .FALSE. if the earth rotation is neglected
!!    
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
!!         LBU_RU       : logical for budget of RU (wind component along x)
!!                        .TRUE. = budget of RU         
!!                        .FALSE. = no budget of RU 
!!         LBU_RV       : logical for budget of RV (wind component along y)
!!                        .TRUE. = budget of RV         
!!                        .FALSE. = no budget of RV 
!!         LBU_RW        : logical for budget of RW (wind component along z)
!!                        .TRUE. = budget of RW         
!!                        .FALSE. = no budget of RW 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine DYN_SOURCE )
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      C.Lac - March 2011 - Splitted  from dyn_sources
!!      Q.Rodier 06/15 correction on budget
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
!
USE MODI_SHUMAN
USE MODI_BUDGET
USE MODI_GET_HALO
!  
IMPLICIT NONE
!  
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KRR  ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT            !     at
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !  time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRWS ! Sources of Momentum
!
!
!*       0.2   Declarations of local variables :
!
REAL       ::  ZRV_OV_RD    ! = RV / RD
INTEGER    ::  JWATER       ! loop index on the different types of water
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) ::           &
                              ZWORK1, ZWORK2
INTEGER :: IKU
!
!-------------------------------------------------------------------------------
!
!
!*       1.     COMPUTES THE GRAVITY TERM
!	        -------------------------
!
IKU=SIZE(PTHT,3)
!
IF( .NOT.L1D ) THEN     ! no buoyancy for 1D case
!
  IF(KRR > 0) THEN
!
!   compute the ratio : 1 + total water mass / dry air mass
!
    ZRV_OV_RD = XRV / XRD
    ZWORK1(:,:,:) = 1.
    DO JWATER = 1 , 1+KRRL+KRRI                
      CALL GET_HALO(PRT(:,:,:,JWATER))
      ZWORK1(:,:,:) = ZWORK1(:,:,:) + PRT(:,:,:,JWATER)
    END DO
!
!   compute the virtual potential temperature when water is present in any form
!
    CALL GET_HALO(PTHT)
!
    
    ZWORK2(:,:,:) = PTHT(:,:,:) * (1. + PRT(:,:,:,1)*ZRV_OV_RD) / ZWORK1(:,:,:)
  ELSE
!
!   compute the virtual potential temperature when water is absent
!
    ZWORK2(:,:,:) = PTHT(:,:,:)
  END IF
!
!   compute the gravity term
!
  PRWS(:,:,:) = PRWS + XG * MZM(1,IKU,1, ( (ZWORK2/PTHVREF) - 1. ) * PRHODJ )
!
!    the extrapolation for the PTHT and the THVREF must be the same at the
!    ground
!
!
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GRAVITY 
