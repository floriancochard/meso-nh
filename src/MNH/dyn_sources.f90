!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_DYN_SOURCES
!     #######################
!
INTERFACE
!
      SUBROUTINE DYN_SOURCES( KRR,KRRL, KRRI,                              &
                              PUT, PVT, PWT, PTHT, PRT,                    &
                              PCORIOX, PCORIOY, PCORIOZ, PCURVX, PCURVY,   &
                              PRHODJ, PZZ, PTHVREF, PEXNREF,               &
                              PRUS, PRVS, PRWS, PRTHS                      )
!
INTEGER,                  INTENT(IN)    :: KRR  ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT, PWT   ! variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT            !     at
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !  time t
!
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCORIOX   !    2D horizontal
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCORIOY   !      Coriolis 
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCORIOZ   !       factors
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCURVX    !    2D horizontal
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCURVY    !  curvature factors
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ       ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF   ! Exner function
                                          ! of the reference state
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS ! Sources of Momentum
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS            ! Sources of theta
!
END SUBROUTINE DYN_SOURCES
!
END INTERFACE
!
END MODULE MODI_DYN_SOURCES 
!     ######################################################################
      SUBROUTINE DYN_SOURCES( KRR,KRRL, KRRI,                              &
                              PUT, PVT, PWT, PTHT, PRT,                    &
                              PCORIOX, PCORIOY, PCORIOZ, PCURVX, PCURVY,   &
                              PRHODJ, PZZ, PTHVREF, PEXNREF,               &
                              PRUS, PRVS, PRWS, PRTHS                      )
!     ######################################################################
!
!!****  *DYN_SOURCES * - routine to compute the curvature, coriolis and gravity terms
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the dynamical sources
!!    (curvature, coriolis and gravity terms) for each component of the
!!    momentum. The curvature terms arise in case of non-cartesian geometry
!!    only. For both curvature and coriolis terms, the "thin shell" 
!!    approximation can be assumed or not. Gravity affects only the vertical
!!    component of the momentum.
!!      This routine also adds the Lagrangien derivative of the pressure as a
!!    source for the theta evolution.
!!
!!     
!!**  METHOD
!!    ------
!!      The horizontal curvature and coriolis field arrays are expanded in
!!    the vertical direction with the SPREAD intrinsics to match with the
!!    Shuman operators. The source term for theta due to the advection of the
!!    Exner function is treated in an advective form (or non-conservative form).
!!    Values of the calorific capacity of the melting and of the potential
!!    virtual temperature are recovered taking care of the number of water
!!    variables. The different sources terms are stored for the budget
!!    computations.
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
!!      Original    06/06/94 
!!	Corrections 06/08/94 (J.-P. Lafore) 
!!	Corrections 17/10/94 (Stein) For LCORIO
!!	Corrections 22/12/94 (Stein) add the pressure term in the theta evolution
!!  Corrections 30/12/94 (J.P. Lafore) bug corrections for the pressure term 
!!  Corrections 16/03/95 (Stein) remove R from the historical variables and
!!                                   correction of the pressure term   
!!  Corrections 14/04/95 (Masson Stein) bugs in the curvatureand Coriolis
!!                                   terms for PRUS,PRVS, PRWS         
!!  Corrections 01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!  Corrections 16/10/95 (J. Stein)     change the budget calls 
!!  Corrections 16/02/96 (J.-P. Pinty)  Introduce L1D switches
!!  Corrections 19/12/96 (J.-P. Pinty)  Update the CALL BUDGET
!!  Corrections 03/12/02  (P. Jabouille)  add no thinshell condition
!!  Correction     06/10  (C.Lac) Exclude L1D for Coriolis term 
!!  Modification   03/11  (C.Lac) Split the gravity term due to buoyancy
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_CONF
USE MODD_CST
USE MODD_DYN
!
USE MODI_SHUMAN
USE MODI_GRADIENT_M
USE MODI_BUDGET
!
USE MODE_MPPDB
!  
IMPLICIT NONE
!  
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KRR  ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT, PWT   ! variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT            !     at
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !  time t
!
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCORIOX   !    2D horizontal
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCORIOY   !      Coriolis 
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCORIOZ   !       factors
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCURVX    !    2D horizontal
REAL, DIMENSION(:,:),     INTENT(IN)    :: PCURVY    !  curvature factors
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ       ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF   ! Virtual Temperature
                                          ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF   ! Exner function
                                          ! of the reference state
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS ! Sources of Momentum
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS            ! Sources of theta
!
!
!*       0.2   Declarations of local variables :
!
REAL       ::  ZCPD_OV_RD   ! = CPD / RD
REAL       ::  ZG_OV_CPD    ! =-XG / XCPD
INTEGER    ::  JWATER       ! loop index on the different types of water
INTEGER    ::  IKU          ! array size along the k direction 
REAL       ::  ZD1          ! DELTA1 (switch 0/1) for thinshell approximation
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) ::           &
                              ZWORK1, ZWORK2, ZWORK3, ZRUT, ZRVT
!
!-------------------------------------------------------------------------------
!
!
!*       1.     COMPUTES THE TRUE VELOCITY COMPONENTS
!	        -------------------------------------
!
ZRUT(:,:,:) = PUT(:,:,:) * MXM(PRHODJ(:,:,:))
ZRVT(:,:,:) = PVT(:,:,:) * MYM(PRHODJ(:,:,:))
!
IKU = SIZE(PUT,3)
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTES THE CURVATURE TERMS
!    	        ----------------------------
!
! Only when earth rotation is considered but not in 1D and CARTESIAN cases
!
IF ((.NOT.L1D).AND.(.NOT.LCARTESIAN) )  THEN 
  IF ( LTHINSHELL ) THEN           !  THINSHELL approximation  
!
    ZWORK1(:,:,:) = SPREAD( PCURVX(:,:),DIM=3,NCOPIES=IKU ) / XRADIUS
    ZWORK2(:,:,:) = SPREAD( PCURVY(:,:),DIM=3,NCOPIES=IKU ) / XRADIUS
!
    PRUS(:,:,:) = PRUS                                   &
    + ZRUT * MXM(  MYF(PVT) * ZWORK1 )                   &
    + MXM( MYF(ZRVT*PVT) * ZWORK2 )
!
    PRVS(:,:,:) = PRVS                                   &
    - MYM( MXF(ZRUT*PUT) * ZWORK1 )                      &
    - ZRVT * MYM( MXF(PUT) * ZWORK2 ) 
!
  ELSE                           !  NO THINSHELL approximation
!
    ZWORK3(:,:,:) = 1.0 / ( XRADIUS + MZF(1,IKU,1,PZZ(:,:,:)) )
    ZWORK1(:,:,:) = SPREAD( PCURVX(:,:),DIM=3,NCOPIES=IKU )
    ZWORK2(:,:,:) = SPREAD( PCURVY(:,:),DIM=3,NCOPIES=IKU )
    CALL MPPDB_CHECK3DM("DYN_SOOURCES:ZWORK3,ZWORK1,ZWORK2",PRECISION,&
                      & ZWORK3,ZWORK1,ZWORK2,&
                      & MXM( MYF(ZRVT*PVT) * ZWORK2 * ZWORK3 ) , &
                      & MXM( ( MYF(PVT) * ZWORK1 - MZF(1,IKU,1,PWT) ) * ZWORK3 ) ,&
                      & MYF(PVT) * ZWORK1 - MZF(1,IKU,1,PWT) , &
                      & MYF(PVT) , MZF(1,IKU,1,PWT) , MXM(PWT) , MYM(PWT) )
    CALL MPPDB_CHECK3DM("DYN_SOOURCES:SUITE",PRECISION,&
         &  MXM(ZRVT),MXM(PVT),MXM(PWT),MXM(ZWORK1),MXM(ZWORK2),MXM(ZWORK3)  )
!
    PRUS(:,:,:) = PRUS                                              &
    + MXM( MYF(ZRVT*PVT) * ZWORK2 * ZWORK3 )                        &
    + ZRUT * MXM( ( MYF(PVT) * ZWORK1 - MZF(1,IKU,1,PWT) ) * ZWORK3 ) 
!
    PRVS(:,:,:) = PRVS                                              &
    - MYM( MXF(ZRUT*PUT) * ZWORK1 * ZWORK3 )                        &
    - ZRVT * MYM( (MXF(PUT) * ZWORK2 + MZF(1,IKU,1,PWT) ) * ZWORK3 )
!
    PRWS(:,:,:) = PRWS                                              &
    +MZM(1,IKU,1, ( MXF(ZRUT*PUT) + MYF(ZRVT*PVT) ) * ZWORK3 ) 
!
  END IF
!
END IF
!
IF (LBUDGET_U) CALL BUDGET (PRUS,1,'CURV_BU_RU')
IF (LBUDGET_V) CALL BUDGET (PRVS,2,'CURV_BU_RV')
IF (LBUDGET_W) CALL BUDGET (PRWS,3,'CURV_BU_RW')
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTES THE CORIOLIS TERMS
!	        ---------------------------
!
IF (LCORIO)   THEN 
!
  ZWORK3(:,:,:) = SPREAD( PCORIOZ(:,:),DIM=3,NCOPIES=IKU ) * PRHODJ(:,:,:)
!
  PRUS(:,:,:) = PRUS + MXM( ZWORK3 * MYF(PVT) )
  PRVS(:,:,:) = PRVS - MYM( ZWORK3 * MXF(PUT) )
!
  IF ((.NOT.LTHINSHELL) .AND. (.NOT.L1D)) THEN
!
    ZWORK1(:,:,:) = SPREAD( PCORIOX(:,:),DIM=3,NCOPIES=IKU) * PRHODJ(:,:,:) 
    ZWORK2(:,:,:) = SPREAD( PCORIOY(:,:),DIM=3,NCOPIES=IKU) * PRHODJ(:,:,:)
!
    PRUS(:,:,:) = PRUS - MXM( ZWORK2 * MZF(1,IKU,1,PWT) ) 
!
    PRVS(:,:,:) = PRVS - MYM( ZWORK1 * MZF(1,IKU,1,PWT) )
!
    PRWS(:,:,:) = PRWS + MZM( 1,IKU,1,ZWORK2 * MXF(PUT) + ZWORK1 * MYF(PVT) )
!
  END IF
!
END IF                      
!
IF (LBUDGET_U) CALL BUDGET (PRUS,1,'COR_BU_RU')
IF (LBUDGET_V) CALL BUDGET (PRVS,2,'COR_BU_RV')
IF (LBUDGET_W) CALL BUDGET (PRWS,3,'COR_BU_RW')
!
!-------------------------------------------------------------------------------
!
!*       4.     COMPUTES THE THETA SOURCE TERM DUE TO THE REFERENCE PRESSURE
!	        ------------------------------------------------------------
!
IF (LCARTESIAN .OR. LTHINSHELL) THEN
  ZD1=0.
ELSE
  ZD1=1.
ENDIF
!
IF( .NOT.L1D ) THEN
!
  IF (KRR > 0) THEN
    ZCPD_OV_RD = XCPD / XRD 
    ZG_OV_CPD  = -XG / XCPD
!
! stores the specific heat capacity (Cph) in ZWORK1
!                                        
    ZWORK1(:,:,:) = XCPD + XCPV * PRT(:,:,:,1)  ! gas mixing
    DO JWATER = 2,1+KRRL             !  loop on the liquid components  
      ZWORK1(:,:,:) = ZWORK1(:,:,:) + XCL * PRT(:,:,:,JWATER)
    END DO
!
    DO JWATER = 2+KRRL,1+KRRL+KRRI   !  loop on the solid components   
      ZWORK1(:,:,:) = ZWORK1(:,:,:) + XCI * PRT(:,:,:,JWATER)
    END DO
!
! computes the source term
!
    PRTHS(:,:,:) = PRTHS(:,:,:) +  PRHODJ(:,:,:)                             &
      * ( ( XRD + XRV * PRT(:,:,:,1) ) * ZCPD_OV_RD / ZWORK1(:,:,:)  - 1. )  &
      * PTHT(:,:,:)/PEXNREF(:,:,:)*MZF(1,IKU,1,PWT(:,:,:))*(ZG_OV_CPD/PTHVREF(:,:,:) &
      -ZD1*4./7.*PEXNREF(:,:,:)/( XRADIUS+MZF(1,IKU,1,PZZ(:,:,:)) ))
!
  END IF
!
  IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'PREF_BU_RTH')
!
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYN_SOURCES 
